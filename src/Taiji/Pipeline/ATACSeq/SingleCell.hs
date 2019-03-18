{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq.SingleCell (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Control.Lens
import           Scientific.Workflow
import           Control.Monad.IO.Class                (liftIO)
import           Control.Monad.Reader                  (asks)

import           Taiji.Pipeline.ATACSeq.Functions
import           Taiji.Pipeline.ATACSeq.Config

builder :: Builder ()
builder = do
    nodeS "SC_Read_Input" [| \_ -> do
        input <- asks _atacseq_input
        liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readATACSeqTSV input "scATAC-seq"
            else readATACSeq input "scATAC-seq"
        |] $ do
            submitToRemote .= Just False
            note .= "Read ATAC-seq data information from input file."

    nodeS "SC_Download_Data" 'atacDownloadData $ do
        submitToRemote .= Just False
        note .= "Download data."

    node' "SC_Get_Fastq" 'atacGetFastq $ submitToRemote .= Just False

    nodeS "SC_Make_Index" 'atacMkIndex $ do
        note .= "Generate the BWA index."

    nodePS 1 "SC_Align" 'atacAlign $ do
        remoteParam .= "--ntasks-per-node=2"  -- slurm
        --remoteParam .= "-pe smp 2"  -- sge
        note .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."

    path ["SC_Read_Input", "SC_Download_Data", "SC_Get_Fastq", "SC_Make_Index", "SC_Align"]