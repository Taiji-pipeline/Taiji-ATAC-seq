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
    , alignQC
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
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)
import           System.IO

import           Taiji.Pipeline.ATACSeq.Config

type ATACSeqWithSomeFile = ATACSeq [Either SomeFile (SomeFile, SomeFile)]

type ATACSeqMaybePair tag1 tag2 filetype =
    Either (ATACSeq (File tag1 filetype))
           (ATACSeq (File tag2 filetype, File tag2 filetype))

type ATACSeqEitherTag tag1 tag2 filetype = Either (ATACSeq (File tag1 filetype))
                                                  (ATACSeq (File tag2 filetype))

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

atacGetFastq :: [ATACSeq [Either SomeFile (SomeFile, SomeFile)]]
         -> [ATACSeqMaybePair '[] '[Pairend] 'Fastq]
atacGetFastq inputs = flip concatMap inputs $ \input ->
    fromMaybe (error "A mix of single and pairend fastq was found") $
        splitExpByFileEither $ input & replicates.mapped.files %~ f
  where
    f fls = map (bimap fromSomeFile (bimap fromSomeFile fromSomeFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

atacGetBam :: [ATACSeq [Either SomeFile (SomeFile, SomeFile)]]
           -> [ATACSeqEitherTag '[] '[Pairend] 'Bam]
atacGetBam inputs = flip concatMap inputs $ \input ->
    fromMaybe (error "A mix of single and pairend fastq was found") $
        splitExpByFileEither $ input & replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bam) $ lefts fls) $ \fl ->
        if fl `hasTag` Pairend
            then Left $ fromSomeFile fl
            else Right $ fromSomeFile fl

atacBamToBed :: ATACSeqConfig config
             => Either (ATACSeq (File tags1 'Bam)) (ATACSeq (File tags2 'Bam))
             -> WorkflowConfig config (ATACSeq (File '[Gzip] 'Bed))
atacBamToBed input = do
    dir <- asks _atacseq_output_dir >>= getPath
    liftIO $ case input of
        Left x -> nameWith dir "bed.gz" (\output fl ->
            coerce $ fun output fl) x
        Right x -> nameWith dir "bed.gz" (\output fl ->
            coerce $ fun output fl) x
  where
    fun output fl = bam2Bed_ output (\x -> not $ chrom x `elem` ["chrM", "M"]) fl

atacCallPeak :: (ATACSeqConfig config, SingI tags)
             => ATACSeq (File tags 'Bed)
             -> WorkflowConfig config (ATACSeq (File '[] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath
    let fn output fl = callPeaks output fl Nothing $ do
            callSummits .= False
            mode .= NoModel (-100) 200
    liftIO $ nameWith dir "narrowPeak" fn input

alignQC :: Either (ATACSeq (File tags1 'Bam)) (ATACSeq (File tags2 'Bam))
        -> IO (T.Text, [(Int, Int)])
alignQC = either fun fun
  where
    fun x = do
        result <- mapM bamStat $ x^..replicates.folded.files
        return (x^.eid, result)
