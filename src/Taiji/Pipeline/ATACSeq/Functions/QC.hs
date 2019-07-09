{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}
module Taiji.Pipeline.ATACSeq.Functions.QC where

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
import           Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts
import           Taiji.Prelude
import qualified Taiji.Utils.DataFrame as DF

alignQC :: ATACSeqConfig config
        => [(String, Double, Double)]
        -> ReaderT config IO ()
alignQC xs
    | null xs = return ()
    | otherwise = do
        dir <- qcDir
        let output = dir <> "mapping_qc.html"
        liftIO $ savePlots output [] [stackBar df]
  where
    df = DF.mkDataFrame ["percent mapped reads", "percent chrM reads"] names $
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

dupRate :: ATACSeqConfig config
        => [ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))]
        -> ReaderT config IO ()
dupRate input
    | null input = return ()
    | otherwise = do
        dir <- qcDir
        let output = dir <> "/qc_duplication_rate.html"
        liftIO $ savePlots output [] [stackBar df]
  where
    df = DF.mkDataFrame ["duplication_rate"] names [dat]
    (names, dat) = unzip $ mapMaybe getDupRate input
    getDupRate e = case either getResult getResult (e^.replicates._2.files) of
        Nothing -> Nothing
        Just r ->  Just
            (T.pack $ printf "%s_rep%d" (T.unpack $ e^.eid) (e^.replicates._1), r)
    getResult fl = case M.lookup "QC" (fl^.info) of
        Nothing -> Nothing
        Just txt -> let xs = T.lines txt
                        dup = read $ T.unpack $ last $ T.words $ last xs
                        total = read $ T.unpack $ (T.words $ head xs) !! 1
                    in Just $ 100 * dup / total

peakSignal :: ATACSeqConfig config
           => ( ATACSeq S (Either (File '[] 'Bed) (File '[Gzip] 'Bed))
              , File '[] 'Bed )   -- ^ Reference peak list
           -> ReaderT config IO (ATACSeq S (File '[] 'Other))
peakSignal (atac, peakFl) = do 
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    let output = printf "%s/%s_rep%d.signal.bin" dir (T.unpack $ atac^.eid)
            (atac^.replicates._1)
    atac & replicates.traverse.files %%~ liftIO . ( \fl -> do
        regions <- runResourceT $ runConduit $ streamBed (peakFl^.location) .|
            concatMapC (splitBedBySizeLeft 500) .| sinkList :: IO [BED3]
        let readInput = either (streamBed . (^.location))
                (streamBedGzip . (^.location)) 
        res <- runResourceT $ runConduit $ readInput fl .|
            rpkmBed regions :: IO (U.Vector Double)
        encodeFile output res
        return $ location .~ output $ emptyFile )

-- | Compute correlation between experiments.
peakCor :: ATACSeqConfig config
        => [ATACSeq S (File '[] 'Other)]
        -> ReaderT config IO ()
peakCor inputs = do
    dir <- qcDir
    let output = dir <> "/qc_correlation.html"
        names = flip map inputs $ \x ->
            (x^.eid) <> "_rep" <> T.pack (show $ x^.replicates._1)
    liftIO $ do
        cor <- fmap (pearsonMatByRow . fromRows) $ forM inputs $ \input ->
            decodeFile $ input^.replicates._2.files.location
        savePlots output [] $ [heatmap $ DF.mkDataFrame names names $ toRowLists cor]

teQC :: ATACSeqConfig config => [(T.Text, U.Vector Double)] -> ReaderT config IO ()
teQC xs = do
    dir <- qcDir
    let output = dir <> "/qc_tss_enrichment.html"
    liftIO $ savePlots output [] [p1,p2]
  where
    p1 = stackBar $ DF.mkDataFrame ["TSS Enrichment"] names [map U.maximum vals]
    p2 = stackLine $ DF.mkDataFrame names (map (T.pack . show) [-1000..999::Int]) $
        map U.toList vals
    (names, vals) = unzip xs

computeTE :: ATACSeqConfig config
          => ATACSeq S (Either (File '[] 'Bed) (File '[Gzip] 'Bed))
          -> ReaderT config IO (T.Text, U.Vector Double)
computeTE input = do
    genes <- asks _atacseq_annotation >>= liftIO . readGenes . fromJust
    let tss = flip map genes $ \g -> (geneChrom g, geneLeft g, geneStrand g)
        readInput = either (streamBed . (^.location))
            (streamBedGzip . (^.location))
    res <- liftIO $ runResourceT $ runConduit $
        readInput (input^.replicates._2.files) .| tssEnrichment tss
    return (input^.eid <> "_rep" <> T.pack (show $ input^.replicates._1), res)

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
        vec <- liftIO $ UM.replicate 2000 (0 :: Int)
        mapM_C $ mapM_ $ UM.unsafeModify vec (+1)
        liftIO $ normalize <$> U.unsafeFreeze vec
    normalize vec = slideAverage 25 $ U.map ((/bk) . fromIntegral) vec
      where
        bk = fromIntegral (U.sum (U.take 100 vec) + U.sum (U.drop 1900 vec)) / 200
    getCutSite x = BED3 (x^.chrom) i $ i + 1
      where
        i = case x^.strand of
            Just False -> x^.chromEnd - 1
            _ -> x^.chromStart
    f x ys = flip map ys $ \y -> case y^.strand of
        Just False -> 1999 - (x^.chromStart - y^.chromStart)
        _ -> x^.chromStart - y^.chromStart 
    regions = map ( \(chr, x, str) ->
        BED chr (x - 1000) (x + 1000) Nothing Nothing $ Just str ) tss

-- | Plot QC for fragment size distribution.
fragDistrQC :: ATACSeqConfig config 
            => [Maybe (T.Text, U.Vector Double)]
            -> ReaderT config IO ()
fragDistrQC xs = do
    dir <- qcDir
    let output = dir <> "fragment_size_distr.html"
    liftIO $ savePlots output [] plts
  where
    plts = flip map (catMaybes xs) $ \(nm, dat) -> 
        let df = DF.mkDataFrame ["fragment_count"] labels [U.toList dat]
        in stackLine df <> title (T.unpack nm)
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