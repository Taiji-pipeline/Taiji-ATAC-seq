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
import           Data.Either                           (either)
import           Data.Maybe                            (fromJust)
import           Data.Monoid                           ((<>))
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config
import           Taiji.Pipeline.ATACSeq.Core.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ if ".tsv" == reverse $ take 4 $ reverse input
            then readATACSeqTSV input "ATAC-seq"
            else readATACSeq input "ATAC-seq"
        |] $ do
            submitToRemote .= Just False
            note .= "Read ATAC-seq data information from input file."

    nodeS "Download_Data" 'atacDownloadData $ do
        submitToRemote .= Just False
        note .= "Download data."

    node' "Get_Fastq" 'atacGetFastq $ submitToRemote .= Just False
    node' "Get_Bam" [| \(x,y) -> atacGetBam x ++ y |] $ submitToRemote .= Just False

    nodeS "Make_Index" 'atacMkIndex $ do
        note .= "Generate the BWA index."

    nodePS 1 "Align" 'atacAlign $ do
        remoteParam .= "--ntasks-per-node=2"  -- slurm
        --remoteParam .= "-pe smp 2"  -- sge
        note .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."

    nodePS 1 "Filter_Bam" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
        liftIO $ bitraverse
            (filterBam (dir, ".filt.bam"))
            (filterBam (dir, ".filt.bam"))
            input
        |] $ do
            note .= "Remove low quality tags using: samtools -F 0x70c -q 30"

    nodePS 1 "Remove_Duplicates" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
        picard <- fromJust <$> asks _atacseq_picard
        liftIO $ bitraverse
            (removeDuplicates picard (dir, ".filt.dedup.bam"))
            (removeDuplicates picard (dir, ".filt.dedup.bam"))
            input
        |] $ do
            note .= "Remove duplicated reads using picard."

    nodePS 1 "Bam_To_Bed" 'atacBamToBed $ do
        note .= "Convert Bam file to Bed file."

    node' "Get_Bed" [| \(input1, input2) ->
        let f [x] = x
            f _   = error "Must contain exactly 1 file"
        in mapped.replicates.mapped.files %~ f $ merge $ (atacGetBed input1) ++
            (input2 & mapped.replicates.mapped.files %~ Right)
        |] $ submitToRemote .= Just False

    nodePS 1 "Merge_Bed" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bed")
        liftIO $ concatBed (dir, ".merged.bed.gz") input
        |] $ return ()
    nodePS 1 "Call_Peak" 'atacCallPeak $ return ()

    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index", "Align"]
    ["Download_Data", "Align"] ~> "Get_Bam"
    path ["Get_Bam", "Filter_Bam", "Remove_Duplicates", "Bam_To_Bed"]
    ["Download_Data", "Bam_To_Bed"] ~> "Get_Bed"
    path ["Get_Bed", "Merge_Bed", "Call_Peak"]

    nodeP 1 "Align_QC" [| either alignQC alignQC |] $ return ()
    ["Align"] ~> "Align_QC"

    node' "Dup_QC" [| map (either dupQC dupQC) |] $ submitToRemote .= Just False
    ["Remove_Duplicates"] ~> "Dup_QC"

    node' "Correlation_QC_Prep" [| \(beds, peaks) ->
        let peaks' = map (\x -> (x^.eid, runIdentity (x^.replicates) ^. files)) peaks
        in flip map beds $ \bed -> ( bed, fromJust $ lookup (bed^.eid) peaks' )
        |] $ submitToRemote .= Just False
    nodeP 1 "Correlation_QC" 'peakQC $ return ()
    ["Get_Bed", "Call_Peak"] ~> "Correlation_QC_Prep"
    path ["Correlation_QC_Prep", "Correlation_QC"]

    nodeS "Report_QC" 'reportQC $ submitToRemote .= Just False
    ["Align_QC", "Dup_QC"] ~> "Report_QC"
