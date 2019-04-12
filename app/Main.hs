{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Utils
import           Control.Lens                  ((&), (.~))
import           Data.Aeson                    (FromJSON, ToJSON)
import           Data.Default
import           GHC.Generics                  (Generic)
import           Scientific.Workflow
import           Scientific.Workflow.Main      (MainOpts (..), defaultMainOpts,
                                                mainWith)
import           Taiji.Pipeline.ATACSeq        (builder)
import           Taiji.Pipeline.ATACSeq.Types

data ATACSeqOpts = ATACSeqOpts
    { outputDir :: Directory
    , bwaIndex  :: Maybe FilePath
    , genome    :: Maybe FilePath
    , input     :: FilePath
    , picard    :: Maybe FilePath
    , motifFile :: Maybe FilePath
    , genomeIndex :: Maybe FilePath
    } deriving (Generic)

instance ATACSeqConfig ATACSeqOpts where
    _atacseq_output_dir = outputDir
    _atacseq_bwa_index = bwaIndex
    _atacseq_genome_fasta = genome
    _atacseq_input = input
    _atacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
    _atacseq_genome_index = genomeIndex
    _atacseq_motif_file = motifFile

instance Default ATACSeqOpts where
    def = ATACSeqOpts
        { outputDir = asDir "output"
        , bwaIndex = Nothing
        , genome = Nothing
        , input = "input.yml"
        , picard = Nothing
        , motifFile = Nothing
        , genomeIndex = Nothing
        }

instance FromJSON ATACSeqOpts
instance ToJSON ATACSeqOpts

mainWith defaultMainOpts
    { programHeader = "Taiji-ATAC-Seq"
    , workflowConfigType = Just ''ATACSeqOpts
    } builder
