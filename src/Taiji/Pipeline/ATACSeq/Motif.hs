{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq.Motif (builder) where

import           Bio.Motif                              (readMEME)
import           Control.Lens
import           Control.Monad.IO.Class                 (liftIO)
import           Control.Monad.Reader                   (asks)
import           Data.List.Split                        (chunksOf)
import           Data.Maybe                             (fromMaybe)
import           Data.Monoid                            ((<>))
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config
import           Taiji.Pipeline.ATACSeq.Motif.Functions

builder :: Builder ()
builder = do
    nodeS "Merge_Peaks" 'atacMergePeaks $ do
        note .= "Merge peaks called from different samples together to form " <>
            "a non-overlapping set of open chromatin regions."
    nodeS "Find_TFBS_Prep" [| \region -> do
        motifFile <- fromMaybe (error "Motif file is not specified!") <$>
            asks _atacseq_motif_file
        motifs <- liftIO $ readMEME motifFile
        return $ ContextData region $ chunksOf 100 motifs
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for parallel execution."
    nodeSharedPS 1 "Find_TFBS_Union" [| \x -> atacFindMotifSiteAll 1e-4 x |] $ do
        note .= "Identify TF binding sites in open chromatin regions using " <>
            "the FIMO's motif scanning algorithm. " <>
            "Use 1e-5 as the default p-value cutoff."

    node' "Get_TFBS_Prep" [| uncurry ContextData |] $ do
        submitToRemote .= Just False

    nodeSharedPS 1 "Get_TFBS" [| atacGetMotifSite 50 |] $ do
        note .= "Retrieve motif binding sites for each sample."
    path ["Call_Peak", "Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS_Union"]
    ["Find_TFBS_Union", "Call_Peak"] ~> "Get_TFBS_Prep"
    ["Get_TFBS_Prep"] ~> "Get_TFBS"
