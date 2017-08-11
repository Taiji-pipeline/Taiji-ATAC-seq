{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Bitraversable            (bitraverse)
import           Data.Coerce                   (coerce)
import           Data.Either                   (lefts, rights)
import           Data.List                     (nub)
import           Data.Maybe                    (fromJust, fromMaybe, mapMaybe)
import           Data.Singletons               (SingI)
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config

builder :: Builder ()
builder = do
    nodeS "Make_Index" 'atacMkIndex $ return ()
    nodeS "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ readATACSeq input "ATAC-seq"
        |] $ submitToRemote .= Just False
    nodeS "Download_Data" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath
        liftIO $ downloadData dir input
        |] $ submitToRemote .= Just False
    node' "Get_Fastq" 'getFastq $ submitToRemote .= Just False
    nodeSharedPS 1 "Align" [| \(ContextData idx input) -> do
        dir <- asks _atacseq_output_dir >>= getPath
        liftIO $ bitraverse
            (bwaAlign1 (dir, "bam") idx (return ()))
            (bwaAlign2 (dir, "bam") idx (return ()))
            input
        |] $ return ()

    nodePS 1 "Filter_Bam" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath
        liftIO $ bitraverse
            (filterBam (dir, "filt.bam"))
            (filterBam (dir, "filt.bam"))
            input
        |] $ return ()
    nodePS 1 "Remove_Duplicates" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath
        picard <- fromJust <$> asks _atacseq_picard
        liftIO $ bitraverse
            (removeDuplicates picard (dir, "filt.dedup.bam"))
            (removeDuplicates picard (dir, "filt.dedup.bam"))
            input
        |] $ return ()
    nodePS 1 "Bam_To_Bed" 'atacBamToBed $ return ()
    node' "Merge_Bed_Prep" [| \input ->
        let f [x] = x
            f _   = error "Must contain exactly 1 file"
        in mapped.replicates.mapped.files %~ f $ mergeExps input
        |] $ submitToRemote .= Just False
    nodePS 1 "Merge_Bed" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath
        liftIO $ concatBed (dir, "merged.bed.gz") input
        |] $ return ()
    nodePS 1 "Call_Peak" 'atacCallPeak $ return ()

    path ["Make_Index", "Read_Input", "Download_Data"]
    ["Make_Index", "Download_Data"] ~> "Get_Fastq"
    path ["Get_Fastq", "Align", "Filter_Bam", "Remove_Duplicates"
        , "Bam_To_Bed", "Merge_Bed_Prep", "Merge_Bed", "Call_Peak"]

atacMkIndex :: ATACSeqConfig config => () -> WorkflowConfig config FilePath
atacMkIndex _ = do
    genome <- asks (fromJust . _atacseq_genome_fasta)
    dir <- asks (fromJust . _atacseq_bwa_index)
    liftIO $ bwaMkIndex genome dir

downloadData :: FilePath
             -> [ATACSeq [Either SomeFile (SomeFile, SomeFile)]]
             -> IO [ATACSeq [Either SomeFile (SomeFile, SomeFile)]]
downloadData dir = traverse.replicates.traverse.files.traverse %%~ download
  where
    download :: Either SomeFile (SomeFile, SomeFile)
             -> IO (Either SomeFile (SomeFile, SomeFile))
    download input@(Left (SomeFile fl)) = if getFileType fl == SRA
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
                sraToFastq dir (coerce fl :: File '[] 'SRA)
        else return input
    download x = return x

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
            callSummits .= True
            mode .= NoModel (-100) 200
    liftIO $ nameWith dir "narrowPeak" fn input
