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
    , atacGetBed
    , atacBamToBed
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
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Bio.Seq.IO                    (mkIndex)
import           Conduit
import           Control.Lens
import           Control.Monad                 (forM, mzero)
import           Control.Monad.Base            (liftBase)
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Morph           (hoist)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Conduit.Zlib             (ungzip)
import           Data.Either                   (lefts)
import           Data.List                     (intercalate)
import           Data.List.Ordered             (nubSort)
import qualified Data.Map.Strict               as M
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Monoid                   ((<>))
import           Data.Promotion.Prelude.List   (Elem)
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Shelly                        (fromText, mkdir_p, shelly)
import           Statistics.Correlation        (pearsonMatByRow,
                                                spearmanMatByRow)
import           Statistics.Matrix             (fromRows, toRowLists)

import           Taiji.Pipeline.ATACSeq.Config

type ATACSeqWithSomeFile = ATACSeq N [Either SomeFile (SomeFile, SomeFile)]

type ATACSeqMaybePair tag1 tag2 filetype =
    Either (ATACSeq S (File tag1 filetype))
           (ATACSeq S (File tag2 filetype, File tag2 filetype))

type ATACSeqEitherTag tag1 tag2 filetype = Either (ATACSeq S (File tag1 filetype))
                                                  (ATACSeq S (File tag2 filetype))

atacMkIndex :: ATACSeqConfig config => [a] -> WorkflowConfig config [a]
atacMkIndex input
    | null input = return input
    | otherwise = do
        genome <- asks (fromJust . _atacseq_genome_fasta)
        -- Generate BWA index
        dir <- asks (fromJust . _atacseq_bwa_index)
        liftIO $ bwaMkIndex genome dir
        return input

atacDownloadData :: ATACSeqConfig config
                 => [ATACSeqWithSomeFile]
                 -> WorkflowConfig config [ATACSeqWithSomeFile]
atacDownloadData dat = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Download"))
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ downloadFiles dir

atacGetFastq :: [ATACSeqWithSomeFile]
             -> [ATACSeqMaybePair '[] '[Pairend] 'Fastq]
atacGetFastq inputs = concatMap split $ concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap fromSomeFile (bimap fromSomeFile fromSomeFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

atacAlign :: ATACSeqConfig config
          => ATACSeqMaybePair '[] '[Pairend] 'Fastq
          -> WorkflowConfig config (ATACSeqEitherTag '[] '[Pairend] 'Bam)
atacAlign input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
    idx <- asks (fromJust . _atacseq_bwa_index)
    liftIO $ bwaAlign (dir, ".bam") idx (bwaCores .= 2) input

atacGetBam :: [ATACSeqWithSomeFile]
           -> [ATACSeqEitherTag '[] '[Pairend] 'Bam]
atacGetBam inputs = concatMap split $ concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bam) $ lefts fls) $ \fl ->
        if fl `hasTag` Pairend
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl

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
             => ATACSeqEitherTag tags1 tags2 'Bam
             -> WorkflowConfig config (ATACSeq S (File '[Gzip] 'Bed))
atacBamToBed input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    liftIO $ case input of
        Left x -> mapFileWithDefName (dir++"/") ".bed.gz" (\output fl ->
            coerce $ fun output fl) x
        Right x -> mapFileWithDefName (dir++"/") ".bed.gz" (\output fl ->
            coerce $ fun output fl) x
  where
    fun output fl = bam2Bed_ output (\x -> not $ chrom x `elem` ["chrM", "M"]) fl

atacCallPeak :: (ATACSeqConfig config, SingI tags)
             => ATACSeq S (File tags 'Bed)
             -> WorkflowConfig config (ATACSeq S (File '[] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    opts <- asks _atacseq_callpeak_opts
    let fn output fl = callPeaks output fl Nothing opts
    liftIO $ mapFileWithDefName (dir++"/") ".narrowPeak" fn input

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

alignQC :: ATACSeq S (File tags 'Bam)
        -> IO ((T.Text, Int), (Int, Int))
alignQC e = do
    result <- bamStat $ head $ e^..replicates.folded.files
    return ((e^.eid, runIdentity (e^.replicates) ^. number), result)

dupQC :: ATACSeq S (File tags 'Bam)
      -> ((T.Text, Int), Maybe Double)
dupQC e = ((e^.eid, runIdentity (e^.replicates) ^. number), result)
  where
    fl = runIdentity (e^.replicates) ^. files
    result = do
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
    | length (e^.replicates) < 2 = return Nothing
    | otherwise = do
        regions <- fmap sortBed $ (readBed (peakFl^.location) :: Source IO BED) =$=
            concatMapC (splitBedBySizeLeft 250) $$ sinkList

        (labels, values) <- fmap unzip $ forM (e^.replicates) $ \r -> do
            readcounts <- case r^.files of
                Left fl -> readBed (fl^.location) $$
                    hoist liftBase (rpkmSortedBed regions)
                Right fl -> runResourceT $ sourceFile (fl^.location) =$=
                    ungzip =$= linesUnboundedAsciiC =$=
                    mapC (fromLine :: _ -> BED) $$
                    hoist liftBase (rpkmSortedBed regions)
            return (r^.number, readcounts)
        if null values
            then return Nothing
            else do
                let mat = fromRows values
                    cor1 = pearsonMatByRow mat
                    cor2 = spearmanMatByRow mat
                return $ Just (labels, toRowLists cor1, toRowLists cor2)
