{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.ATACSeq.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils
import           Control.Lens
import           Control.Monad.IO.Class                (liftIO)
import           Control.Monad.Reader                  (asks)
import           Data.Bitraversable                    (bitraverse)
import           Data.Maybe                            (fromJust)
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config
import           Taiji.Pipeline.ATACSeq.Core.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ readATACSeq input "ATAC-seq"
        |] $ submitToRemote .= Just False
    nodeS "Download_Data" 'atacDownloadData $ submitToRemote .= Just False
    node' "Get_Fastq" 'atacGetFastq $ submitToRemote .= Just False
    node' "Get_Bam" [| \(x,y) -> atacGetBam x ++ y |] $ submitToRemote .= Just False
    nodeS "Make_Index" 'atacMkIndex $ return ()
    nodePS 1 "Align" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath
        idx <- asks (fromJust . _atacseq_bwa_index)
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

    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index", "Align"]
    ["Download_Data", "Align"] ~> "Get_Bam"
    path ["Get_Bam", "Filter_Bam", "Remove_Duplicates", "Bam_To_Bed"
        , "Merge_Bed_Prep", "Merge_Bed", "Call_Peak"]

    nodeP 1 "Align_QC" 'alignQC $ return ()
    ["Align"] ~> "Align_QC"
