{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}
module Taiji.Pipeline.ATACSeq.Functions.QC
    ( plotFragDistr
    , plotAlignQC
    , plotDupRate
    , plotTE
    , plotPeakCor

    , computeTE
    , fragDistr
    , alignStat
    , peakSignal
    ) where

import Language.Javascript.JMacro
import           Bio.Pipeline.Report
import Bio.Data.Bed.Utils
import Bio.Data.Bed.Types
import Bio.RealWorld.GENCODE
import Bio.Data.Bed
import Bio.Utils.Functions (slideAverage)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Data.Binary
import           Bio.Pipeline.Utils
import qualified Data.Map.Strict               as M
import qualified Data.Text as T
import Bio.HTS
import Statistics.Correlation (pearsonMatByRow)
import Statistics.Matrix (fromRows, toRowLists)

import           Taiji.Pipeline.ATACSeq.Types
import Taiji.Utils.Plot.ECharts
import           Taiji.Prelude
import qualified Taiji.Utils.DataFrame as DF

plotAlignQC :: [(String, Double, Double)]
            -> Maybe EChart
plotAlignQC xs
    | null xs = Nothing
    | otherwise = Just $ yAxisLabel "Percent" >+> options >+> toolbox >+> plt 
  where
    plt = (bar df){_height=480,_width=w}
    w = max 480 $ fromIntegral (length xs) * 30 + 40
    options = [jmacroE| { grid: {top: 60} } |]
    df = DF.mkDataFrame ["mappable", "chrM"] names $
        [p, chrM]
    (names, p, chrM) = unzip3 $
        map (\(a,b,c) -> (T.pack a, 100*b, 100*c)) xs

-- | Return name, percent_mapped_reads and percent_chrM_reads.
alignStat :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
          -> IO (String, Double, Double)
alignStat e = do
    stats <- either bamStat bamStat $ e^.replicates._2.files
    let p = fromIntegral (_mapped_tags stats) / fromIntegral (_total_tags stats)
        chrM = case lookup "chrM" (_mapped_tags_by_chrom stats) of
            Nothing -> case lookup "M" (_mapped_tags_by_chrom stats) of
                Nothing -> 0
                Just x -> fromIntegral x / fromIntegral (_mapped_tags stats)
            Just x -> fromIntegral x / fromIntegral (_mapped_tags stats)
    return (nm, p, chrM)
  where
    nm = printf "%s_rep%d" (T.unpack $ e^.eid) (e^.replicates._1)

plotDupRate :: [ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))]
            -> [EChart]
plotDupRate [] = []
plotDupRate input =
    [ yAxisLabel "Usable reads (Million)" >+> options >+> toolbox >+> (bar totalReads){_height=480,_width=w}
    , yAxisLabel "Percent duplicated" >+> options >+> toolbox >+> (bar dupRate){_height=480,_width=w}
    ]
  where
    w = max 480 $ fromIntegral (length input) * 20 + 40
    options = [jmacroE| { grid: {top: 60} }|]
    (dupRate, totalReads) = DF.unzip $ DF.mkDataFrame [""] names [dat]
    (names, dat) = unzip $ mapMaybe getDupRate input
    getDupRate e = case either getResult getResult (e^.replicates._2.files) of
        Nothing -> Nothing
        Just r ->  Just
            (T.pack $ printf "%s_rep%d" (T.unpack $ e^.eid) (e^.replicates._1), r)
    getResult fl = do
        total <- fmap (read . T.unpack) $ M.lookup "READ" $ fl^.info
        dup <- fmap (read . T.unpack) $ M.lookup "DUPLICATE TOTAL" $ fl^.info
        return (100 * dup / total, total / 1000000)

peakSignal :: ATACSeqConfig config
           => ( ATACSeq S (Either (File '[Gzip] 'Bed) (File '[PairedEnd, Gzip] 'Bed))
              , File '[] 'Bed )   -- ^ Reference peak list
           -> ReaderT config IO (ATACSeq S (File '[] 'Other))
peakSignal (atac, peakFl) = do 
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    let output = printf "%s/%s_rep%d.signal.bin" dir (T.unpack $ atac^.eid)
            (atac^.replicates._1)
    atac & replicates.traverse.files %%~ liftIO . ( \fl -> do
        regions <- runResourceT $ runConduit $ streamBed (peakFl^.location) .|
            concatMapC (splitBedBySizeLeft 500) .| sinkList :: IO [BED3]
        let readInput = either (streamBedGzip . (^.location))
                (\x -> streamBedGzip (x^.location) .| concatMapC f) 
        res <- runResourceT $ runConduit $ readInput fl .|
            rpkmBed regions :: IO (U.Vector Double)
        encodeFile output res
        return $ location .~ output $ emptyFile )
  where
    f :: BED -> [BED]
    f x = [ BED (x^.chrom) (x^.chromStart) (x^.chromStart + 1) (x^.name) Nothing (Just True)
          , BED (x^.chrom) (x^.chromEnd - 1) (x^.chromEnd) (x^.name) Nothing (Just False) ]

-- | Compute correlation between experiments.
plotPeakCor :: [ATACSeq S (File '[] 'Other)] -> IO (Maybe EChart)
plotPeakCor inputs
    | length inputs == 0 || length inputs > 100 = return Nothing
    | otherwise = do
        let names = flip map inputs $ \x ->
                (x^.eid) <> "_rep" <> T.pack (show $ x^.replicates._1)
        cor <- fmap (pearsonMatByRow . fromRows) $ forM inputs $ \input ->
            decodeFile $ input^.replicates._2.files.location
        let plt = heatmap $ DF.orderDataFrame id $
                DF.mkDataFrame names names $ toRowLists cor
        return $ Just $ addAttr (title "Pearson correlation") $ addAttr toolbox plt

plotTE :: [(T.Text, U.Vector Double)] -> [EChart]
plotTE [] = []
plotTE xs = [p1{_width=w,_height=480}, p2{_width=480,_height=480}]
  where
    p1 = addAttr toolbox $ addAttr (yAxisLabel "TSS Enrichment") $ stackBar $
        DF.mkDataFrame [""] names [map U.maximum vals]
    w = max 480 $ fromIntegral (length xs) * 30 + 40
    p2 = addAttr toolbox $ addAttr (yAxisLabel "TSS Enrichment") $ 
        line $ DF.mkDataFrame names (map (T.pack . show) [-2000..2000::Int]) $
        map U.toList vals
    (names, vals) = unzip xs

computeTE :: ATACSeqConfig config
          => ATACSeq S (Either (File '[Gzip] 'Bed) (File '[PairedEnd, Gzip] 'Bed))
          -> ReaderT config IO (T.Text, U.Vector Double)
computeTE input = do
    genes <- getAnnotation >>= liftIO . readGenes
    let tss = flip map genes $ \g -> 
            let chr = geneChrom g
                str = geneStrand g
                x = if str then geneLeft g else geneRight g
            in (chr, x, str)
        readInput = either (streamBedGzip . (^.location))
            (\x -> streamBedGzip (x^.location) .| concatMapC f)
    res <- liftIO $ runResourceT $ runConduit $
        readInput (input^.replicates._2.files) .| tssEnrichment tss
    return (input^.eid <> "_rep" <> T.pack (show $ input^.replicates._1), res)
  where
    f :: BED -> [BED]
    f x = [ BED (x^.chrom) (x^.chromStart) (x^.chromStart + 1) (x^.name) Nothing (Just True)
          , BED (x^.chrom) (x^.chromEnd - 1) (x^.chromEnd) (x^.name) Nothing (Just False) ]

-- | Transcription Start Site (TSS) Enrichment Score
-- The TSS enrichment calculation is a signal to noise calculation.
-- The reads around a reference set of TSSs are collected to
-- form an aggregate distribution of reads centered on the TSSs
-- and extending to 1000 bp in either direction (for a total of 2000bp).
-- This distribution is then normalized by taking the average
-- read depth in the 100 bps at each of the end flanks of
-- the distribution (for a total of 200bp of averaged data) and
-- calculating a fold change at each position over that average
-- read depth. This means that the flanks should start at 1,
-- and if there is high read signal at transcription start
-- sites (highly open regions of the genome) there should be an
-- increase in signal up to a peak in the middle. We take the signal
-- value at the center of the distribution after this normalization
-- as our TSS enrichment metric. Used to evaluate ATAC-seq. 
--
-- 1. For each TSS, get per base coverage for the 1000 bp flanking region.
--    We will get a matrix of #TSS x 2000 bp dimension.
-- 2. Do a column sum of the matrix.
-- 3. Sum of the coverage of the endFlank (100bp) at both ends and divide
--    by 200 bp to get a normalization factor.
-- 4. divide the the normalization factor for -1900 to + 1900 bp to get per base normalized coverage.
-- 5. do a smoothing with a defined window (50bp by default) using zoo::rollmean.
-- 6. select the highest value within a window (highest_tss_flank, 50 bp by default)
--    around the TSS because the highest peak is not necessary at exactly the TSS
--    site (position 0)
tssEnrichment :: [(B.ByteString, Int, Bool)]   -- ^ A list of TSSs
              -> ConduitT BED o (ResourceT IO) (U.Vector Double)
tssEnrichment tss = mapC getCutSite .| intersectBedWith f regions .| sink
  where
    sink = do
        vec <- liftIO $ UM.replicate 4000 0
        mapM_C $ mapM_ $ UM.unsafeModify vec (+1)
        liftIO $ normalize <$> U.unsafeFreeze vec
    normalize vec = slideAverage 5 $ U.map (/bk) vec
      where
        bk = (U.sum (U.take 100 vec) + U.sum (U.drop 3900 vec)) / 200 + 0.1
    getCutSite x = BED3 (x^.chrom) i $ i + 1
      where
        i = case x^.strand of
            Just False -> x^.chromEnd - 76
            _ -> x^.chromStart + 75
    f x ys = flip map ys $ \y -> case y^.strand of
        Just False -> 3999 - (x^.chromStart - y^.chromStart)
        _ -> x^.chromStart - y^.chromStart 
    regions = map ( \(chr, x, str) ->
        BED chr (x - 2000) (x + 2000) Nothing Nothing $ Just str ) tss

-- | Plot QC for fragment size distribution.
plotFragDistr :: [Maybe (T.Text, U.Vector Double)] -> [EChart]
plotFragDistr xs = flip map (catMaybes xs) $ \(nm, dat) -> 
    let df = DF.mkDataFrame ["fraction"] labels [U.toList dat]
    in addAttr (yAxisLabel "Fraction") $
        addAttr (title (T.unpack nm)) (stackLine df){_width=400, _height=480}
  where
    labels = map (T.pack . show) [0..1000 :: Int]

-- | Compute fragment size distribution.
fragDistr :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
          -> IO (Maybe (T.Text, U.Vector Double))
fragDistr input = case input^.replicates._2.files of
    Left _ -> return Nothing
    Right fl -> do
        r <- runResourceT $ runConduit $
            streamBam (fl^.location) .| fragmentSizeDistr 1001
        return $ Just (nm, r)
  where
    nm = (input^.eid) <> "_rep" <> (T.pack $ show $ input^.replicates._1)
