{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Utils
import           Control.Lens                  ((&), (.~))
import           Data.Aeson                    (FromJSON)
import           Data.Default
import           GHC.Generics                  (Generic)
import Data.Binary (Binary)

import           Control.Workflow
import qualified Control.Workflow.Coordinator.Drmaa as D
import Control.Workflow.Main
import Data.Proxy (Proxy(..))

import           Taiji.Pipeline.ATACSeq        (builder)
import           Taiji.Pipeline.ATACSeq.Types

data ATACSeqOpts = ATACSeqOpts
    { outputDir :: Directory
    , bwaIndex  :: Maybe FilePath
    , genome    :: Maybe FilePath
    , input     :: FilePath
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

instance FromJSON ATACSeqOpts
instance Binary ATACSeqOpts

decodeDrmaa :: String -> Int -> FilePath -> IO D.DrmaaConfig
decodeDrmaa ip port _ = D.getDefaultDrmaaConfig
    ["remote", "--ip", ip, "--port", show port]

build "wf" [t| SciFlow ATACSeqOpts |] builder

main :: IO ()
main = defaultMain "" cmd wf
  where
    cmd = [ runParser decodeDrmaa
          , viewParser
          , remoteParser (Proxy :: Proxy D.Drmaa) ]

