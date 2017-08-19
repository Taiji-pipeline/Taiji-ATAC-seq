{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                  ((.=))
import           Data.Aeson                    (FromJSON, ToJSON)
import           Data.Default
import           GHC.Generics                  (Generic)
import           Scientific.Workflow
import           Scientific.Workflow.Main      (MainOpts (..), defaultMainOpts,
                                                mainWith)
import           Taiji.Pipeline.ATACSeq        (builder)
import qualified Taiji.Pipeline.ATACSeq.Config as C

data ATACSeqOpts = ATACSeqOpts
    { outputDir :: Directory
    , bwaIndex  :: Maybe FilePath
    , genome    :: Maybe FilePath
    , input     :: FilePath
    , picard    :: Maybe FilePath
    } deriving (Generic)

instance C.ATACSeqConfig ATACSeqOpts where
    _atacseq_output_dir = outputDir
    _atacseq_bwa_index = bwaIndex
    _atacseq_genome_fasta = genome
    _atacseq_input = input
    _atacseq_picard = picard

instance Default ATACSeqOpts where
    def = ATACSeqOpts
        { outputDir = asDir "output"
        , bwaIndex = Nothing
        , genome = Nothing
        , input = "input.yml"
        , picard = Nothing
        }

instance FromJSON ATACSeqOpts
instance ToJSON ATACSeqOpts

-- | Instantiate the "ATACSeqConfig".
initialization :: () -> WorkflowConfig ATACSeqOpts ()
initialization _ = return ()

mainWith defaultMainOpts { programHeader = "Taiji-ATAC-Seq" } $ do
    nodeS "Initialization" 'initialization $ submitToRemote .= Just False
    ["Initialization"] ~> "Make_Index"
    builder
