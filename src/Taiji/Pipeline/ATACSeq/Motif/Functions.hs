{-# LANGUAGE DataKinds        #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.ATACSeq.Motif.Functions
    ( atacMergePeaks
    , atacFindMotifSiteAll
    , atacGetMotifSite
    ) where

import           Bio.Data.Bed                  (BED, BED3 (..), BEDLike (..),
                                                getMotifPValue, getMotifScore,
                                                intersectBed, mergeBed,
                                                motifScan, readBed, readBed',
                                                writeBed, _npPeak)
import           Bio.Data.Experiment
import           Bio.Motif
import           Bio.Pipeline.Instances        ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils
import           Bio.Seq.IO
import           Conduit
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Default
import           Data.Maybe                    (fromJust)
import           Scientific.Workflow
import           System.IO.Temp                (emptyTempFile)

import           Taiji.Pipeline.ATACSeq.Config

atacMergePeaks :: ATACSeqConfig config
               => [ATACSeq S (File '[] 'NarrowPeak)]
               -> WorkflowConfig config (File '[] 'Bed)
atacMergePeaks input = do
    dir <- asks _atacseq_output_dir >>= getPath
    let fls = input^..folded.replicates.folded.files
        openChromatin = dir ++ "/openChromatin.bed"
    liftIO $ do
        peaks <- mapM (readBed' . (^.location)) fls :: IO [[BED3]]
        mergeBed (concat peaks) $$ writeBed openChromatin
        return $ location .~ openChromatin $ emptyFile

atacFindMotifSiteAll :: ATACSeqConfig config
                     => ContextData (File '[] 'Bed) [Motif]
                     -> WorkflowConfig config (File '[] 'Bed)
atacFindMotifSiteAll (ContextData openChromatin motifs) = do
    dir <- asks _atacseq_output_dir >>= getPath
    genome <- fromJust <$> asks _atacseq_genome_index
    liftIO $ withGenome genome $ \g -> do
        output <- emptyTempFile dir "motif_sites_part."
        (readBed (openChromatin^.location) :: Source IO BED3) =$=
            motifScan g motifs def p =$= getMotifScore g motifs def =$=
            getMotifPValue (Just (1 - p * 10)) motifs def $$ writeBed output
        return $ location .~ output $ emptyFile
  where
    p = 1e-5

atacGetMotifSite :: ATACSeqConfig config
                 => Int -- ^ region around summit
                 -> ([File '[] 'Bed], [ATACSeq S (File '[] 'NarrowPeak)])
                 -> WorkflowConfig config [ATACSeq S (File '[] 'Bed)]
atacGetMotifSite window (tfbs, experiment) = do
    dir <- asks _atacseq_output_dir >>= getPath
    mapM (mapFileWithDefName dir "" fun) experiment
  where
    fun output fl = liftIO $ do
        peaks <- readBed (fl^.location) =$= mapC getSummit $$ sinkList
        (mapM_ (readBed . (^.location)) tfbs :: Source IO BED) =$=
            intersectBed peaks $$ writeBed output
        return $ location .~ output $ emptyFile
    getSummit pk = let c = chromStart pk + (fromJust . _npPeak) pk
                   in BED3 (chrom pk) (c - window) (c + window)
