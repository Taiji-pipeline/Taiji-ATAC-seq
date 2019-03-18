{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq (builder) where

import           Scientific.Workflow

import qualified Taiji.Pipeline.ATACSeq.Bulk as Bulk
import qualified Taiji.Pipeline.ATACSeq.Motif as Motif
import qualified Taiji.Pipeline.ATACSeq.SingleCell as SingleCell

builder :: Builder ()
builder = Bulk.builder >> Motif.builder >> SingleCell.builder
