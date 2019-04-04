{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}
module Taiji.Pipeline.ATACSeq.Functions.QC where

import           Bio.Data.Bed.Utils (rpkmSortedBed)
import           Bio.Pipeline.Report
import           Bio.Data.Experiment
import qualified Data.Vector.Unboxed as U
import Data.Maybe
import           Bio.Pipeline.Utils
import qualified Data.Map.Strict               as M
import           Control.Lens
import           Control.Monad.Reader          (asks)
import Data.Serialize (encode)
import qualified Data.Text as T
import qualified Data.ByteString as B
import Conduit
import           Scientific.Workflow
import Bio.HTS
import Text.Printf (printf)
import qualified Data.Vector.Unboxed.Mutable as UM

import           Taiji.Types
import           Taiji.Pipeline.ATACSeq.Config

saveQC :: ATACSeqConfig config
         => ((QC, QC), QC, [Maybe QC])
         -> WorkflowConfig config ()
saveQC ((q1,q2), q3, q4) = do
    dir <- asks _atacseq_output_dir >>= getPath
    let output = dir ++ "/atac_seq.qc"
    liftIO $ B.writeFile output $ encode $ [q1, q2, q3] ++ catMaybes q4

combineMappingQC :: [(String, Double, Double)] -> (QC, QC)
combineMappingQC xs = ( QC "percent_mapped_reads" (QCVector name p) Bar
    , QC "percent_chrM_reads" (QCVector name chrM) Bar )
  where
    (name, p, chrM) = unzip3 xs

-- | Return name, percent_mapped_reads and percent_chrM_reads.
getMappingQC :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
               -> IO (String, Double, Double)
getMappingQC e = do
    stats <- either bamStat bamStat $ e^.replicates._2.files
    let p = fromIntegral (_mapped_tags stats) / fromIntegral (_total_tags stats)
        chrM = case lookup "chrM" (_mapped_tags_by_chrom stats) of
            Nothing -> case lookup "M" (_mapped_tags_by_chrom stats) of
                Nothing -> 0
                Just x -> fromIntegral x / fromIntegral (_mapped_tags stats)
            Just x -> fromIntegral x / fromIntegral (_mapped_tags stats)
    return (name, p, chrM)
  where
    name = printf "%s_rep%d" (T.unpack $ e^.eid) (e^.replicates._1)

getDupQC :: [ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))]
         -> QC
getDupQC input = QC "duplication_rate" (QCVector names p) Bar
  where
    (names, p) = unzip $ mapMaybe getDupRate input
    getDupRate e = case either getResult getResult (e^.replicates._2.files) of
        Nothing -> Nothing
        Just r ->  Just
            (printf "%s_rep%d" (T.unpack $ e^.eid) (e^.replicates._1), r)
    getResult fl = case M.lookup "QC" (fl^.info) of
        Nothing -> Nothing
        Just txt -> let xs = T.lines txt
                        dup = read $ T.unpack $ last $ T.words $ last xs
                        total = read $ T.unpack $ (T.words $ head xs) !! 1
                    in Just $ dup / total

{-
getPeakQC :: (Elem 'Gzip tags1 ~ 'False, Elem 'Gzip tags2 ~ 'True)
          => ( ATACSeq N (Either (File tags1 'Bed) (File tags2 'Bed))
             , File tags3 'NarrowPeak )
          -> IO [QC]
getPeakQC (e, peakFl)
    | IM.size (e^.replicates) < 2 = return []
    | otherwise = do
        regions <- fmap sortBed $ runConduit $
            (readBed (peakFl^.location) :: ConduitT () BED IO ()) .|
            concatMapC (splitBedBySizeLeft 250) .| sinkList
        readcounts <- fmap IM.toList $ forM (e^.replicates) $ \r ->
            case r^.files of
                Left fl -> runConduit $ readBed (fl^.location) .|
                    rpkmSortedBed regions
                Right fl -> runResourceT $ runConduit $
                    sourceFile (fl^.location) .| ungzip .|
                    linesUnboundedAsciiC .| mapC (fromLine :: _ -> BED) .|
                    rpkmSortedBed regions
        return $ flip map (comb readcounts) $ \((i1, r1), (i2,r2)) ->
            let name = Pair 
                    (printf "%s_rep%d" (T.unpack $ e^.eid) i1)
                    (printf "%s_rep%d" (T.unpack $ e^.eid) i2)
            in QC "signal_correlation" name (pearson $ U.zip r1 r2) Nothing
  where
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb [] = []
    -}

{-
-- | Correlation of ATAC-seq signals between samples.
atacCorrelation :: ( (ATACSeq S (File '[Gzip] 'Bed), ATACSeq S (File '[Gzip] 'Bed))
                   , File '[] 'Bed )
                -> IO QC
atacCorrelation ((e1, e2), peakFl) = do
    regions <- fmap sortBed $ runConduit $
        (readBed (peakFl^.location) :: ConduitT () BED IO ()) .|
        concatMapC (splitBedBySizeLeft 250) .| sinkList
    let f1 = e1^.replicates._2.files.location
        f2 = e2^.replicates._2.files.location
    r1 <- runResourceT $ runConduit $
        sourceFile f1 .| ungzip .|
        linesUnboundedAsciiC .| mapC (fromLine :: _ -> BED) .|
        rpkmSortedBed regions
    r2 <- runResourceT $ runConduit $
        sourceFile f2 .| ungzip .|
        linesUnboundedAsciiC .| mapC (fromLine :: _ -> BED) .|
        rpkmSortedBed regions
    let name = Pair (T.unpack $ e1^.eid) (T.unpack $ e2^.eid)
    return $ QC "signal_correlation" name (pearson $ U.zip r1 r2) Nothing
-}

{-
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
atacTSSEnrichment :: -> WorkflowConfig config QC
atacTSSEnrichment = do
    anno <- 
-}

getFragmentQC :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
              -> IO (Maybe QC)
getFragmentQC e = case e^.replicates._2.files of
    Left _ -> return Nothing  
    Right x -> do
        r <- runResourceT $ runConduit $
            streamBam (x^.location) .| fragmentSizeDistr
        return $ Just $ QC ("fragment_size_" ++ T.unpack (e^.eid))
            (QCVector (map show [0..1000]) $ U.toList r) Density

fragmentSizeDistr :: PrimMonad m => ConduitT BAM o m (U.Vector Double)
fragmentSizeDistr = do
    vec <- lift $ UM.replicate n 0
    mapM_C $ f vec
    vec' <- lift $ U.unsafeFreeze vec
    return $ U.map (/ (U.sum vec')) vec'
  where
    f v x | s >= n = return ()
          | otherwise = UM.modify v (+1) s
      where
        s = abs $ tLen x
    n = 1001
{-# INLINE fragmentSizeDistr #-}