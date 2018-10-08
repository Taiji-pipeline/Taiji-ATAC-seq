{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE TemplateHaskell       #-}
module Taiji.Pipeline.ATACSeq.Core.Functions
    ( atacMkIndex
    , atacDownloadData
    , atacGetFastq
    , atacAlign
    , atacGetBam
    , atacFilterBamSort
    , atacGetBed
    , atacBamToBed
    , atacConcatBed
    , atacCallPeak

    -- * QC
    , alignQC
    , dupQC
    , reportQC
    , peakQC
    ) where

import           Bio.ChIPSeq                   (rpkmSortedBed)
import           Bio.Data.Bed                  (BED, chrom, fromLine, readBed,
                                                sortBed, splitBedBySizeLeft)
import           Bio.Data.Experiment
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Conduit
import           Control.Lens
import           Control.Monad                 (forM, mzero)
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Conduit.Zlib             (ungzip)
import           Data.Either                   (lefts)
import qualified Data.IntMap.Strict            as IM
import           Data.List                     (intercalate)
import           Data.List.Ordered             (nubSort)
import qualified Data.Map.Strict               as M
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Monoid                   ((<>))
import           Data.Promotion.Prelude.List   (Elem)
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Statistics.Correlation        (pearsonMatByRow,
                                                spearmanMatByRow)
import           Statistics.Matrix             (fromRows, toRowLists)
import           System.IO.Temp                (withTempFile)
import           Text.Printf                   (printf)

import           Taiji.Pipeline.ATACSeq.Config

type ATACSeqWithSomeFile = ATACSeq N [Either SomeFile (SomeFile, SomeFile)]

atacMkIndex :: ATACSeqConfig config => [a] -> WorkflowConfig config [a]
atacMkIndex input
    | null input = return input
    | otherwise = do
        genome <- asks (fromJust . _atacseq_genome_fasta)
        -- Generate BWA index
        dir <- asks (fromJust . _atacseq_bwa_index)
        _ <- liftIO $ bwaMkIndex genome dir
        return input

atacDownloadData :: ATACSeqConfig config
                 => [ATACSeqWithSomeFile]
                 -> WorkflowConfig config [ATACSeqWithSomeFile]
atacDownloadData dat = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Download"))
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ downloadFiles dir

atacGetFastq :: [ATACSeqWithSomeFile]
             -> [ ATACSeq S ( Either
                    (SomeTags 'Fastq) (SomeTags 'Fastq, SomeTags 'Fastq) )]
atacGetFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile (bimap castFile castFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

atacAlign :: ATACSeqConfig config
          => ATACSeq S ( Either
                (SomeTags 'Fastq) (SomeTags 'Fastq, SomeTags 'Fastq) )
          -> WorkflowConfig config ( ATACSeq S
                (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam)) )
atacAlign input = do
    dir <- asks ((<> "/Bam") . _atacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _atacseq_bwa_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> case fl of
        Left f ->
            let f' = fromSomeTags f :: File '[] 'Fastq
            in bwaAlign output idx (Left f') $ defaultBWAOpts & bwaCores .~ 2
        Right (f1,f2) ->
            let f1' = fromSomeTags f1 :: File '[] 'Fastq
                f2' = fromSomeTags f2 :: File '[] 'Fastq
            in bwaAlign output idx (Right (f1', f2')) $ defaultBWAOpts & bwaCores .~ 2
        )

-- | Grab bam files
atacGetBam :: [ATACSeqWithSomeFile]
           -> [ ATACSeq S (
                    Either (File '[] 'Bam) (File '[PairedEnd] 'Bam) )]
atacGetBam inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bam) $ lefts fls) $ \fl ->
        if fl `hasTag` PairedEnd
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl

atacFilterBamSort :: ATACSeqConfig config
    => ATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
    -> WorkflowConfig config ( ATACSeq S
        (Either (File '[CoordinateSorted] 'Bam) (File '[CoordinateSorted, PairedEnd] 'Bam)) )
atacFilterBamSort input = do
    dir <- asks ((<> "/Bam") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . fun output)
        (fmap Right . fun output)
  where
    fun output x = withTempFile "./" "tmp_file." $ \f _ ->
        filterBam "./" f x >>= sortBam "./" output

atacGetBed :: [ATACSeqWithSomeFile]
           -> [ATACSeq S (Either (File '[] 'Bed) (File '[Gzip] 'Bed))]
atacGetBed input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bed) $ lefts fls) $ \fl ->
        if fl `hasTag` Gzip
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl

atacBamToBed :: ATACSeqConfig config
             => ATACSeq S ( Either (File tags1 'Bam) (File tags2 'Bam) )
             -> WorkflowConfig config (ATACSeq S (File '[Gzip] 'Bed))
atacBamToBed input = do
    dir <- asks ((<> "/Bed") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        let fl' = either coerce coerce fl :: File '[] 'Bam
        bam2Bed output (\x -> not $ (x^.chrom) `elem` ["chrM", "M"]) fl' )

atacConcatBed :: ( ATACSeqConfig config
                 , Elem 'Gzip tags1 ~ 'False
                 , Elem 'Gzip tags2 ~ 'True )
              => ATACSeq N (Either (File tags1 'Bed) (File tags2 'Bed))
              -> WorkflowConfig config (ATACSeq S (File '[Gzip] 'Bed))
atacConcatBed input = do
    dir <- asks ((<> "/Bed") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep0_merged.bed.gz" dir (T.unpack $ input^.eid)
    fl <- liftIO $ concatBed output fls
    return $ input & replicates .~ (0, Replicate fl M.empty)
  where
    fls = input^..replicates.folded.files

atacCallPeak :: (ATACSeqConfig config, SingI tags)
             => ATACSeq S (File tags 'Bed)
             -> WorkflowConfig config (ATACSeq S (File '[] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    opts <- asks _atacseq_callpeak_opts
    let output = printf "%s/%s_rep%d.narrowPeak" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO .
        (\fl -> callPeaks output fl Nothing opts)

reportQC :: ATACSeqConfig config
         => ( [((T.Text, Int), (Int, Int))], [((T.Text, Int), Maybe Double)])
         -> WorkflowConfig config ()
reportQC (align, dup) = do
    dir <- asks _atacseq_output_dir >>= getPath
    let output = dir ++ "/QC.tsv"
        header = ["Sample_ID", "Replicate", "Total_Reads", "Mapped_Reads"
            , "Percent_Mapped", "Duplication_Rate"]
        content = flip map samples $ \x -> [T.unpack $ fst x, show $ snd x] ++
            ( case lookup x align of
                Just (mapped', total) -> [show total, show mapped'
                    , show $ fromIntegral mapped' / fromIntegral total]
                Nothing -> ["NA", "NA", "NA"] ) ++
            ( case fromMaybe Nothing (lookup x dup) of
                Nothing -> ["NA"]
                Just d  -> [show d] )
    liftIO $ writeFile output $ unlines $ map (intercalate "\t") $
        header : content
  where
    samples = nubSort $ map fst align ++ map fst dup

alignQC :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
        -> IO ((T.Text, Int), (Int, Int))
alignQC e = do
    result <- either bamStat bamStat $ e^.replicates._2.files
    return ((e^.eid, e^.replicates._1), result)

dupQC :: ATACSeq S (Either (File tags1 'Bam) (File tags2 'Bam))
      -> ((T.Text, Int), Maybe Double)
dupQC e = ( (e^.eid, e^.replicates._1)
          , either getResult getResult $ e^.replicates._2.files )
  where
    getResult fl = do
        txt <- M.lookup "QC" (fl^.info)
        let txt' = filter (\x -> not (T.null x || "#" `T.isPrefixOf` x)) $
                T.lines txt
        if length txt' > 1
            then do
                let [_,unpaired,paired,secondary,unmapped,unpaired_dup,paired_dup,optical_dup,percent_dup,_] = T.splitOn "\t" $ txt' !! 1
                return $ read $ T.unpack $ percent_dup
            else mzero

peakQC :: (Elem 'Gzip tags1 ~ 'False, Elem 'Gzip tags2 ~ 'True)
       => ( ATACSeq N (Either (File tags1 'Bed) (File tags2 'Bed))
          , File tags3 'NarrowPeak )
       -> IO (Maybe ([Int], [[Double]], [[Double]]))
peakQC (e, peakFl)
    | IM.size (e^.replicates) < 2 = return Nothing
    | otherwise = do
        regions <- fmap sortBed $ runConduit $
            (readBed (peakFl^.location) :: ConduitT () BED IO ()) .|
            concatMapC (splitBedBySizeLeft 250) .| sinkList

        (labels, values) <- fmap (unzip . IM.toList) $ forM (e^.replicates) $ \r -> do
            readcounts <- case r^.files of
                Left fl -> runConduit $ readBed (fl^.location) .|
                    rpkmSortedBed regions
                Right fl -> runResourceT $ runConduit $
                    sourceFile (fl^.location) .| ungzip .|
                    linesUnboundedAsciiC .| mapC (fromLine :: _ -> BED) .|
                    rpkmSortedBed regions
            return readcounts
        if null values
            then return Nothing
            else do
                let mat = fromRows values
                    cor1 = pearsonMatByRow mat
                    cor2 = spearmanMatByRow mat
                return $ Just (labels, toRowLists cor1, toRowLists cor2)
