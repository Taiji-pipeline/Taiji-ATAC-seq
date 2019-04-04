{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.ATACSeq.Bulk (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Control.Lens
import           Control.Monad.IO.Class                (liftIO)
import           Control.Monad.Reader                  (asks)
import           Data.Either                           (either)
import           Data.Maybe                            (fromJust)
import           Data.Monoid                           ((<>))
import qualified Data.Text                             as T
import           Scientific.Workflow
import           Text.Printf                           (printf)

import           Taiji.Pipeline.ATACSeq.Config
import           Taiji.Pipeline.ATACSeq.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readATACSeqTSV input "ATAC-seq"
            else readATACSeq input "ATAC-seq"
        |] $ do
            submitToRemote .= Just False
            note .= "Read ATAC-seq data information from input file."
    nodeS "Download_Data" 'atacDownloadData $ do
        submitToRemote .= Just False
        note .= "Download data."
    node' "Get_Fastq" 'atacGetFastq $ submitToRemote .= Just False

    path ["Read_Input", "Download_Data", "Get_Fastq"]

    nodeS "Make_Index" 'atacMkIndex $ note .= "Generate the BWA index."
    node' "Align_Prep" [| fst |] $ submitToRemote .= Just False
    nodePS 1 "Align" 'atacAlign $ do
        remoteParam .= "--ntasks-per-node=2"  -- slurm
        --remoteParam .= "-pe smp 2"  -- sge
        note .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."

    ["Get_Fastq", "SC_Get_Fastq"] ~> "Make_Index"
    ["Get_Fastq", "Make_Index"] ~> "Align_Prep"
    ["Align_Prep"] ~> "Align"

    node' "Get_Bam" [| \(x,y) -> atacGetBam x ++ y |] $ submitToRemote .= Just False

    ["Download_Data", "Align"] ~> "Get_Bam"

    nodePS 1 "Filter_Bam" 'atacFilterBamSort $ do
        note .= "Remove low quality tags using: samtools -F 0x70c -q 30"

    nodePS 1 "Remove_Duplicates" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
        picard <- fromJust <$> asks _atacseq_picard
        let output = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
                (input^.replicates._1)
        input & replicates.traverse.files %%~ liftIO . either
            (fmap Left . removeDuplicates picard output)
            (fmap Right . removeDuplicates picard output)
        |] $ note .= "Remove duplicated reads using picard."

    nodePS 1 "Bam_To_Bed" 'atacBamToBed $ do
        note .= "Convert Bam file to Bed file."

    node' "Get_Bed" [| \(input1, input2) ->
        let f [x] = x
            f _   = error "Must contain exactly 1 file"
        in mapped.replicates.mapped.files %~ f $ mergeExp $ atacGetBed input1 ++
            (input2 & mapped.replicates.mapped.files %~ Right)
        |] $ submitToRemote .= Just False

    nodePS 1 "Merge_Bed" 'atacConcatBed $ return ()
    nodePS 1 "Call_Peak" 'atacCallPeak $ return ()

    node' "Get_Peak" [| \(input1, input2) -> atacGetNarrowPeak input1 ++ input2 
        |] $ submitToRemote .= Just False

    path ["Get_Bam", "Filter_Bam", "Remove_Duplicates", "Bam_To_Bed"]
    ["Download_Data", "Bam_To_Bed"] ~> "Get_Bed"
    path ["Get_Bed", "Merge_Bed", "Call_Peak"]
    ["Download_Data", "Call_Peak"] ~> "Get_Peak"

    nodeP 1 "Align_QC'" 'getMappingQC $ return ()
    node' "Align_QC" 'combineMappingQC $ submitToRemote .= Just False
    path ["Align", "Align_QC'", "Align_QC"]

    node' "Dup_QC" 'getDupQC $ submitToRemote .= Just False
    ["Remove_Duplicates"] ~> "Dup_QC"

    nodeP 1 "Fragment_Size_QC" 'getFragmentQC $ return ()
    ["Remove_Duplicates"] ~> "Fragment_Size_QC"

    {-
    node' "Correlation_QC_Prep" [| \(beds, peaks) ->
        let peaks' = map (\x -> (x^.eid, x^.replicates._2.files)) peaks
        in flip map beds $ \bed -> ( bed, fromJust $ lookup (bed^.eid) peaks' )
        |] $ submitToRemote .= Just False
    nodeP 1 "Correlation_QC" 'peakQC $ return ()
    ["Get_Bed", "Call_Peak"] ~> "Correlation_QC_Prep"
    path ["Correlation_QC_Prep", "Correlation_QC"]
    -}

    {-
    -- Calculate correlation between experiments
    node' "Correlation_Prep" [| \(beds, peak) ->
        let comb (x:xs) = zip (repeat x) xs ++ comb xs
            comb [] = []
        in zip (comb beds) $ repeat peak
        |] $ submitToRemote .= Just False
    ["Merge_Bed", "Merge_Peaks"] ~> "Correlation_Prep"
    nodeP 1 "Correlation" 'atacCorrelation $ return ()
    path ["Correlation_Prep", "Correlation"]
    -}

    nodeS "Report_QC" 'saveQC $ submitToRemote .= Just False
    ["Align_QC", "Dup_QC", "Fragment_Size_QC"] ~> "Report_QC"
