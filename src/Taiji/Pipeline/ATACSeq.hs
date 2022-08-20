{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq (builder) where

import           Bio.Motif                              (readMEME)
import           Bio.Data.Experiment.Parser
import Bio.Data.Experiment.Types
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import Bio.Seq.IO (withGenome, getChrSizes)
import           Data.List.Split                        (chunksOf)
import qualified Data.HashSet as S
import           Control.Workflow
import qualified Data.Vector.Unboxed as U
import qualified Data.Text                             as T
import qualified Data.Text.IO                             as T

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Pipeline.ATACSeq.Functions
import           Taiji.Prelude
import           Taiji.Utils.Plot

builder :: Builder ()
builder = do
-------------------------------------------------------------------------------
-- Reads mapping
-------------------------------------------------------------------------------
    node "Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ mkInputReader input "ATAC-seq" (\_ x -> ATACSeq x)
        |] $ doc .= "Read input data information."
    nodePar "Download_Data" [| atacDownloadData |] $ doc .= "Download data."
    node "Make_Index" [| atacMkIndex |] $ doc .= "Generate genome indices."
    uNode "Get_Fastq" [| return . atacGetFastq |]
    nodePar "Align" [| atacAlign |] $ do
        nCore .= 4
        doc .= "Align reads using BWA: bwa mem -M -k 32."
    path ["Read_Input", "Download_Data", "Make_Index", "Get_Fastq", "Align"]

-------------------------------------------------------------------------------
-- Bam filtering
-------------------------------------------------------------------------------
    uNode "Get_Bam" [| \(x,y) -> return $ atacGetBam x ++ y |]
    ["Make_Index", "Align"] ~> "Get_Bam"
    nodePar "Filter_Bam" [| atacFilterBamSort |] $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
        nCore .= 2
    nodePar "Remove_Duplicates" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> asDir "/Bam")
        let output = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
                (input^.replicates._1)
        input & replicates.traverse.files %%~ liftIO . either
            (fmap Left . removeDuplicates output)
            (fmap Right . removeDuplicates output)
        |] $ doc .= "Remove PCR duplicates."
    path ["Get_Bam", "Filter_Bam", "Remove_Duplicates"]


-------------------------------------------------------------------------------
-- Bed
-------------------------------------------------------------------------------
    uNode "Bam_To_Bed_Prep" [| \(input, x) -> return $ atacGetFilteredBam input ++ x |]
    ["Make_Index", "Remove_Duplicates"] ~> "Bam_To_Bed_Prep"
    nodePar "Bam_To_Bed" [| atacBamToBed |] $ doc .= "Convert Bam file to Bed file."
    path ["Bam_To_Bed_Prep", "Bam_To_Bed"]

    uNode "Get_Bed" [| \(input1, input2) ->
        let f [x] = x
            f _   = error "Must contain exactly 1 file"
        in return $ mapped.replicates.mapped.files %~ f $
            mergeExp $ atacGetBed input1 ++
            (input2 & mapped.replicates.mapped.files %~ Left)
        |]
    ["Make_Index", "Bam_To_Bed"] ~> "Get_Bed"
    nodePar "Merge_Bed" [| atacConcatBed |] $
        doc .= "Merge Bed files from different replicates"
    nodePar "Make_BigWig" [| \input -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> "/BigWig/")
        tmpdir <- fromMaybe "./" <$> asks _atacseq_tmp_dir
        seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
            _atacseq_genome_index )
        let output = printf "%s/%s_rep%d.bw" dir (T.unpack $ input^.eid)
                (input^.replicates._1)
        liftIO $ do
            chrSize <- withGenome seqIndex $ return . getChrSizes
            bedToBigWig output chrSize [] tmpdir $ input^.replicates._2.files
        |] $ doc .= "Generate Bigwig files."
    path ["Get_Bed", "Merge_Bed", "Make_BigWig"]

-------------------------------------------------------------------------------
-- Call peaks
-------------------------------------------------------------------------------
    uNode "Call_Peak_Prep" [| \(beds, inputs) ->
        let ids = S.fromList $ map (^.eid) $ atacGetNarrowPeak inputs
        in return $ filter (\x -> not $ S.member (x^.eid) ids) beds
        |]
    ["Merge_Bed", "Make_Index"] ~> "Call_Peak_Prep"
    nodePar "Call_Peak" [| atacCallPeak |] $ doc .= "Call peaks using MACS2."
    path ["Call_Peak_Prep", "Call_Peak"]

    uNode "Get_Peak" [| \(input1, input2) ->
        let input2' = input2 & mapped.replicates.mapped.files %~ Right
        in return $ atacGetNarrowPeak input1 ++ input2' |]
    ["Make_Index", "Call_Peak"] ~> "Get_Peak"

    node "Merge_Peaks" [| atacMergePeaks |] $ do
        doc .= "Merge peaks called from different samples together to form " <>
            "a non-overlapping set of open chromatin regions."
    path ["Get_Peak", "Merge_Peaks"]

    nodePar "Gene_Count" [| estimateExpr |] $
        doc .= "Estimate gene expression using promoter accessibility."
    node "Make_Expr_Table" [| mkTable |] $ return ()
    path ["Merge_Bed", "Gene_Count", "Make_Expr_Table"]


-------------------------------------------------------------------------------
-- Quality control
-------------------------------------------------------------------------------
    nodePar "Align_Stat" [| liftIO . alignStat |] $ return ()
    path ["Align", "Align_Stat"]

    nodePar "Fragment_Size_Distr" [| liftIO . fragDistr |] $ return ()
    path ["Remove_Duplicates", "Fragment_Size_Distr"]

    uNode "Compute_TE_Prep" [| return . concatMap split |]
    nodePar "Compute_TE" [| computeTE |] $ doc .= "Compute TSS enrichment."
    path ["Get_Bed", "Compute_TE_Prep", "Compute_TE"]

    uNode "Compute_Peak_Signal_Prep" [| \(xs, y) -> 
        return $ zip (concatMap split xs) $ repeat $ fromJust y|]
    nodePar "Compute_Peak_Signal" [| peakSignal |] $ return ()
    ["Get_Bed", "Merge_Peaks"] ~> "Compute_Peak_Signal_Prep"
    path ["Compute_Peak_Signal_Prep", "Compute_Peak_Signal"]

    node "QC" [| \(ali, dup, frag, te, cor) -> do
        cor' <- liftIO $ plotPeakCor cor
        let plts = plotDupRate dup ++ catMaybes [plotAlignQC ali, cor'] ++
                plotFragDistr frag ++ plotTE te
        if null plts
            then return ()
            else do
                dir <- qcDir
                let output = dir <> "/TSSe.tsv"
                liftIO $ T.writeFile output $ T.unlines $
                    map (\(a,b) -> a <> "\t" <> T.pack (show $ U.maximum b)) te
                liftIO $ savePlots (dir <> "/qc.html") [] plts
        |] $ doc .= "Generating QC plots."
    [ "Align_Stat", "Remove_Duplicates", "Fragment_Size_Distr"
        , "Compute_TE", "Compute_Peak_Signal" ] ~> "QC"

    node "Find_TFBS_Prep" [| \case
        Nothing -> return []
        Just pk -> getMotif >>= liftIO . readMEME >>=
            return . zip (repeat pk) . chunksOf 100
        |] $ return ()
    nodePar "Find_TFBS_Union" [| \x -> atacFindMotifSiteAll 5e-5 x |] $ do
        doc .= "Identify TF binding sites in open chromatin regions using " <>
            "the FIMO's motif scanning algorithm. " <>
            "Use 5e-5 as the default p-value cutoff."

    uNode "Get_TFBS_Prep" [| \(x,y) -> return $ zip (repeat x) y|]

    nodePar "Get_TFBS" [| atacGetMotifSite 50 |] $ do
        doc .= "Retrieve motif binding sites for each sample."
    path ["Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS_Union"]
    ["Find_TFBS_Union", "Get_Peak"] ~> "Get_TFBS_Prep"
    ["Get_TFBS_Prep"] ~> "Get_TFBS"