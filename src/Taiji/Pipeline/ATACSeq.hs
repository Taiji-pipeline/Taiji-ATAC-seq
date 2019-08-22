{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq (builder) where

import           Bio.Motif                              (readMEME)
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Data.List.Split                        (chunksOf)
import qualified Data.HashSet as S
import           Control.Workflow
import           Data.Either                           (either)
import qualified Data.Text                             as T

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Pipeline.ATACSeq.Functions
import           Taiji.Prelude

builder :: Builder ()
builder = do
    node "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readATACSeqTSV input "ATAC-seq"
            else readATACSeq input "ATAC-seq"
        |] $ doc .= "Read ATAC-seq data information from input file."
    node "Download_Data" 'atacDownloadData $ doc .= "Download data."
    node "Get_Fastq" [| return . atacGetFastq |] $ return ()

    path ["Read_Input", "Download_Data", "Get_Fastq"]

    node "Make_Index" 'atacMkIndex $ doc .= "Generate the BWA index."
    node "Align_Prep" [| return . fst |] $ return ()
    nodePar "Align" 'atacAlign $ do
        nCore .= 4
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."

    ["Get_Fastq"] ~> "Make_Index"
    ["Get_Fastq", "Make_Index"] ~> "Align_Prep"
    ["Align_Prep"] ~> "Align"

    node "Get_Bam" [| \(x,y) -> return $ atacGetBam x ++ y |] $ return ()

    ["Download_Data", "Align"] ~> "Get_Bam"

    nodePar "Filter_Bam" 'atacFilterBamSort $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"

    nodePar "Remove_Duplicates" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
        let output = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
                (input^.replicates._1)
        input & replicates.traverse.files %%~ liftIO . either
            (fmap Left . removeDuplicates output)
            (fmap Right . removeDuplicates output)
        |] $ doc .= "Remove duplicated reads using picard."

    nodePar "Bam_To_Bed" 'atacBamToBed $ do
        doc .= "Convert Bam file to Bed file."

    node "Get_Bed" [| \(input1, input2) ->
        let f [x] = x
            f _   = error "Must contain exactly 1 file"
        in return $ mapped.replicates.mapped.files %~ f $ mergeExp $ atacGetBed input1 ++
            (input2 & mapped.replicates.mapped.files %~ Left)
        |] $ return ()
    path ["Get_Bam", "Filter_Bam", "Remove_Duplicates", "Bam_To_Bed"]
    ["Download_Data", "Bam_To_Bed"] ~> "Get_Bed"
    nodePar "Merge_Bed" 'atacConcatBed $ return ()
    path ["Get_Bed", "Merge_Bed"]

    node "Call_Peak_Prep" [| \(beds, inputs) -> do
        let ids = S.fromList $ map (^.eid) $ atacGetNarrowPeak inputs
        return $ filter (\x -> not $ S.member (x^.eid) ids) beds
        |] $ return ()
    ["Merge_Bed", "Download_Data"] ~> "Call_Peak_Prep"
    nodePar "Call_Peak" 'atacCallPeak $ return ()
    path ["Call_Peak_Prep", "Call_Peak"]

    node "Get_Peak" [| \(input1, input2) -> return $
        atacGetNarrowPeak input1 ++ input2 
        |] $ return ()
    ["Download_Data", "Call_Peak"] ~> "Get_Peak"

    nodePar "Gene_Count" 'estimateExpr $ return ()
    node "Make_Expr_Table" 'mkTable $ return ()
    path ["Merge_Bed", "Gene_Count", "Make_Expr_Table"]


-------------------------------------------------------------------------------
-- Quality control
-------------------------------------------------------------------------------
    nodePar "Align_Stat" [| liftIO . alignStat |] $ return ()
    node "QC_Align" 'alignQC $ return ()
    path ["Align", "Align_Stat", "QC_Align"]

    node "QC_Duplication" 'dupRate $ return ()
    ["Remove_Duplicates"] ~> "QC_Duplication"

    nodePar "Fragment_Size_Distr" [| liftIO . fragDistr |] $ return ()
    node "QC_Fragment_Size" 'fragDistrQC $ return ()
    path ["Remove_Duplicates", "Fragment_Size_Distr", "QC_Fragment_Size"]

    node "Compute_TE_Prep" [| return . concatMap split |] $ return ()
    nodePar "Compute_TE" 'computeTE $ return ()
    node "QC_TE" 'teQC $ return ()
    path ["Get_Bed", "Compute_TE_Prep", "Compute_TE", "QC_TE"]


    node "Merge_Peaks" 'atacMergePeaks $ do
        doc .= "Merge peaks called from different samples together to form " <>
            "a non-overlapping set of open chromatin regions."

    node "Compute_Peak_Signal_Prep" [| \(xs, y) -> return $
        zip (concatMap split xs) $ repeat $ fromJust y|] $ return ()
    nodePar "Compute_Peak_Signal" 'peakSignal $ return ()
    node "Compute_Peak_Cor" 'peakCor $ return ()
    path ["Get_Peak", "Merge_Peaks"]
    ["Get_Bed", "Merge_Peaks"] ~> "Compute_Peak_Signal_Prep"
    path ["Compute_Peak_Signal_Prep", "Compute_Peak_Signal", "Compute_Peak_Cor"]

    node "Find_TFBS_Prep" [| \region -> do
        motifFile <- fromMaybe (error "Motif file is not specified!") <$>
            asks _atacseq_motif_file
        motifs <- liftIO $ readMEME motifFile
        return $ zip (repeat region) $ chunksOf 100 motifs
        |] $ doc .= "Prepare for parallel execution."
    nodePar "Find_TFBS_Union" [| \x -> atacFindMotifSiteAll 5e-5 x |] $ do
        doc .= "Identify TF binding sites in open chromatin regions using " <>
            "the FIMO's motif scanning algorithm. " <>
            "Use 1e-5 as the default p-value cutoff."

    node "Get_TFBS_Prep" [| \(x,y) -> return $ zip (repeat x) y|] $ return ()

    nodePar "Get_TFBS" [| atacGetMotifSite 50 |] $ do
        doc .= "Retrieve motif binding sites for each sample."
    path ["Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS_Union"]
    ["Find_TFBS_Union", "Get_Peak"] ~> "Get_TFBS_Prep"
    ["Get_TFBS_Prep"] ~> "Get_TFBS"
