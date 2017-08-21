{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq.Core.Functions
    ( atacMkIndex
    , atacDownloadData
    , getFastq
    , atacBamToBed
    , atacCallPeak
    , alignQC
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config

type ATACSeqWithSomeFile = ATACSeq [Either SomeFile (SomeFile, SomeFile)]

atacMkIndex :: ATACSeqConfig config => () -> WorkflowConfig config FilePath
atacMkIndex _ = do
    genome <- asks (fromJust . _atacseq_genome_fasta)
    dir <- asks (fromJust . _atacseq_bwa_index)
    liftIO $ bwaMkIndex genome dir

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

getFastq :: (FilePath, [ATACSeq [Either SomeFile (SomeFile, SomeFile)]])
         -> ContextData FilePath [ MaybePairExp ATACSeq '[] '[Pairend] 'Fastq ]
getFastq (idx, inputs) = ContextData idx $ flip concatMap inputs $ \input ->
    fromMaybe (error "A mix of single and pairend fastq was found") $
        splitExpByFileEither $ input & replicates.mapped.files %~ f
  where
    f :: [Either SomeFile (SomeFile, SomeFile)] -> [MaybePair '[] '[Pairend] 'Fastq]
    f fls = map (bimap fromSomeFile (bimap fromSomeFile fromSomeFile)) $
        filter (either (`someFileIs` Fastq) g) fls
      where
        g (x,y) = x `someFileIs` Fastq && y `someFileIs` Fastq

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
    fun output fl = bam2Bed_ output (const True) fl

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
