{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.ATACSeq (builder) where

import           Scientific.Workflow

import qualified Taiji.Pipeline.ATACSeq.Core as Core

builder :: Builder ()
builder = Core.builder
