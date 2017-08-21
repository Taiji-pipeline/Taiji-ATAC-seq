{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq.Motif (builder) where

import           Bio.Motif                              (readMEME)
import           Control.Lens
import           Control.Monad.IO.Class                 (liftIO)
import           Control.Monad.Reader                   (asks)
import           Data.List.Split                        (chunksOf)
import           Data.Maybe                             (fromJust)
import           Scientific.Workflow

import           Taiji.Pipeline.ATACSeq.Config
import           Taiji.Pipeline.ATACSeq.Motif.Functions

builder :: Builder ()
builder = do
    nodeS "Merge_Peaks" 'atacMergePeaks $ return ()
    nodeS "Find_TFBS_prep" [| \region -> do
        motifFile <- fromJust <$> asks _atacseq_motif_file
        motifs <- liftIO $ readMEME motifFile
        return $ ContextData region $ chunksOf 100 motifs
        |] $ submitToRemote .= Just False
    nodeSharedPS 1 "Find_TFBS_Union" 'atacFindMotifSiteAll $ return ()
    nodeS "Get_TFBS" [| atacGetMotifSite 50 |] $ return ()
    path ["Call_Peak", "Merge_Peaks", "Find_TFBS_prep", "Find_TFBS_Union"]
    ["Find_TFBS_Union", "Call_Peak"] ~> "Get_TFBS"
