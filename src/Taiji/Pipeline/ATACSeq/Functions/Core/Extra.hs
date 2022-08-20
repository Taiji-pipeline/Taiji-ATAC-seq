-- https://www.well-typed.com/blog/2016/09/sharing-conduit/
{-# OPTIONS_GHC -fno-full-laziness #-}
{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedStrings     #-}
module Taiji.Pipeline.ATACSeq.Functions.Core.Extra
    ( atacConcatBed
    ) where

import           Bio.Data.Bed
import           Bio.Pipeline
import           Data.Either                   (lefts, rights)
import qualified Data.Map.Strict               as M
import           Data.List.Singletons (Elem)
import qualified Data.Text                     as T

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Prelude

atacConcatBed :: ( ATACSeqConfig config
                 , Elem 'Gzip tags1 ~ 'True
                 , Elem 'Gzip tags2 ~ 'True
                 , Elem 'PairedEnd tags2 ~ 'True )
              => ATACSeq N (Either (File tags1 'Bed) (File tags2 'Bed))
              -> ReaderT config IO (ATACSeq S (File '[Gzip] 'Bed))
atacConcatBed input = do
    dir <- asks ((<> "/Bed") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep0_merged.bed.gz" dir (T.unpack $ input^.eid)
    fl <- case input^..replicates.folded.files of
        [Left fl] -> return $ location .~ (fl^.location) $ emptyFile
        fls -> do
            let source = do
                    mapM_ (streamBedGzip . (^.location)) (lefts fls)
                    mapM_ (streamBedGzip . (^.location)) (rights fls) .| concatMapC f
            runResourceT $ runConduit $ source .| sinkFileBedGzip output
            return $ location .~ output $ emptyFile
    return $ input & replicates .~ (0, Replicate fl M.empty)
  where
    f :: BED -> [BED]
    f x = [ BED (x^.chrom) (x^.chromStart) (x^.chromStart + 1) (x^.name) Nothing (Just True)
          , BED (x^.chrom) (x^.chromEnd - 1) (x^.chromEnd) (x^.name) Nothing (Just False) ]