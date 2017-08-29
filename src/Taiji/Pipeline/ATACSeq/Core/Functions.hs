{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq.Core.Functions
    ( atacMkIndex
    , atacDownloadData
    , atacGetFastq
    , atacGetBam
    , atacBamToBed
    , atacCallPeak

    -- * QC
    , alignQC
    , dupQC
    , reportQC
    ) where

import           Bio.Data.Bed                  (chrom)
import           Bio.Data.Experiment
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Bio.Seq.IO                    (mkIndex)
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Either                   (lefts)
import           Data.List                     (intercalate)
import           Data.List.Ordered             (nubSort)
import qualified Data.Map.Strict               as M
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Monoid                   ((<>))
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)
import           System.IO

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
        -- Generate sequence index
        seqIndex <- asks (fromJust . _atacseq_genome_index)
        fileExist <- liftIO $ shelly $ test_f $ fromText $ T.pack seqIndex
        liftIO $ if fileExist
            then hPutStrLn stderr "Sequence index exists. Skipped."
            else do
                shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
                hPutStrLn stderr "Generating sequence index"
                mkIndex [genome] seqIndex
        -- Generate BWA index
        dir <- asks (fromJust . _atacseq_bwa_index)
        liftIO $ bwaMkIndex genome dir

        return input

atacDownloadData :: ATACSeqConfig config
                 => [ATACSeqWithSomeFile]
                 -> WorkflowConfig config [ATACSeqWithSomeFile]
atacDownloadData dat = do
    dir <- asks _atacseq_output_dir >>= getPath
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ download dir
  where
    download dir input@(Left (SomeFile fl)) = if getFileType fl == SRA
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
                sraToFastq dir (coerce fl :: File '[] 'SRA)
        else return input
    download _ x = return x

atacGetFastq :: [ATACSeqWithSomeFile]
             -> [ATACSeqMaybePair '[] '[Pairend] 'Fastq]
atacGetFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap fromSomeFile (bimap fromSomeFile fromSomeFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

atacGetBam :: [ATACSeqWithSomeFile]
           -> [ATACSeqEitherTag '[] '[Pairend] 'Bam]
atacGetBam inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bam) $ lefts fls) $ \fl ->
        if fl `hasTag` Pairend
            then Left $ fromSomeFile fl
            else Right $ fromSomeFile fl

atacBamToBed :: ATACSeqConfig config
             => ATACSeqEitherTag tags1 tags2 'Bam
             -> WorkflowConfig config (ATACSeq S (File '[Gzip] 'Bed))
atacBamToBed input = do
    dir <- asks _atacseq_output_dir >>= getPath
    liftIO $ case input of
        Left x -> mapFileWithDefName dir ".bed.gz" (\output fl ->
            coerce $ fun output fl) x
        Right x -> mapFileWithDefName dir ".bed.gz" (\output fl ->
            coerce $ fun output fl) x
  where
    fun output fl = bam2Bed_ output (\x -> not $ chrom x `elem` ["chrM", "M"]) fl

atacCallPeak :: (ATACSeqConfig config, SingI tags)
             => ATACSeq S (File tags 'Bed)
             -> WorkflowConfig config (ATACSeq S (File '[] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    let fn output fl = callPeaks output fl Nothing $ do
            callSummits .= False
            mode .= NoModel (-100) 200
    liftIO $ mapFileWithDefName dir ".narrowPeak" fn input

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
    result = case M.lookup "QC" (fl^.info) of
        Just txt -> let [l] = filter (\x -> not (T.null x || "#" `T.isPrefixOf` x)) $ T.lines txt
                        [_,unpaired,paired,secondary,unmapped,unpaired_dup,paired_dup,optical_dup,percent_dup,_] = T.splitOn "\t" l
                    in Just $ read $ T.unpack $ percent_dup
        Nothing -> Nothing
