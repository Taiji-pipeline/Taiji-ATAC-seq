{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Taiji.Pipeline.ATACSeq.Functions.Motif
    ( atacMergePeaks
    , atacFindMotifSiteAll
    , atacGetMotifSite
    ) where

import           Bio.Data.Bed                  (BED, BED3, BEDLike (..),
                                                intersectBed, mergeBed,
                                                npPeak, streamBed,
                                                readBed, sinkFileBed)
import           Bio.Data.Bed.Utils  (scanMotif, mkCutoffMotif)
import           Bio.Data.Experiment
import           Bio.Motif                     hiding (score)
import           Bio.Pipeline.Instances        ()
import           Bio.Pipeline.Utils
import           Bio.Seq.IO
import           Conduit
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Default
import           Data.Maybe                    (fromJust, fromMaybe)
import           Data.Monoid                   ((<>))
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)
import           System.IO
import           System.IO.Temp                (emptyTempFile)
import           Text.Printf                   (printf)

import           Taiji.Pipeline.ATACSeq.Types

atacMergePeaks :: ATACSeqConfig config
               => [ATACSeq S (File '[] 'NarrowPeak)]
               -> WorkflowConfig config (File '[] 'Bed)
atacMergePeaks input = do
    dir <- asks _atacseq_output_dir >>= getPath
    let fls = input^..folded.replicates.folded.files
        openChromatin = dir ++ "/openChromatin.bed"
    liftIO $ do
        peaks <- mapM (readBed . (^.location)) fls :: IO [[BED3]]
        runResourceT $ runConduit $
            mergeBed (concat peaks) .| sinkFileBed openChromatin
        return $ location .~ openChromatin $ emptyFile

atacFindMotifSiteAll :: ATACSeqConfig config
                     => Double     -- ^ p value
                     -> ContextData (File '[] 'Bed) [Motif]
                     -> WorkflowConfig config (File '[] 'Bed)
atacFindMotifSiteAll p (ContextData openChromatin motifs) = do
    -- Generate sequence index
    genome <- asks ( fromMaybe (error "Genome fasta file was not specified!") .
        _atacseq_genome_fasta )
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _atacseq_genome_index )
    fileExist <- liftIO $ shelly $ test_f $ fromText $ T.pack seqIndex
    liftIO $ if fileExist
        then hPutStrLn stderr "Sequence index exists. Skipped."
        else do
            shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
            hPutStrLn stderr "Generating sequence index"
            mkIndex [genome] seqIndex

    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/TFBS/"))
    liftIO $ withGenome seqIndex $ \g -> do
        output <- emptyTempFile dir "motif_sites_part.bed"
        let motifs' = map (mkCutoffMotif def p) motifs
        runResourceT $ runConduit $
            (streamBed (openChromatin^.location) :: _ _ BED3 _ _) .|
            scanMotif g motifs' .| sinkFileBed output
        return $ location .~ output $ emptyFile

-- | Retrieve TFBS for each experiment
atacGetMotifSite :: ATACSeqConfig config
                 => Int -- ^ region around summit
                 -> ContextData [File '[] 'Bed] (ATACSeq S (File '[] 'NarrowPeak))
                 -> WorkflowConfig config (ATACSeq S (File '[] 'Bed))
atacGetMotifSite window (ContextData tfbs e) = do
    dir <- asks ((<> "/TFBS") . _atacseq_output_dir) >>= getPath
    e & replicates.traversed.files %%~ ( \fl -> liftIO $ do
        let output = printf "%s/%s_rep%d.bed" dir (T.unpack $ e^.eid)
                (e^.replicates._1)
        peaks <- runResourceT $ runConduit $
            streamBed (fl^.location) .| mapC getSummit .| sinkList
        runResourceT $ runConduit $
            (mapM_ (streamBed . (^.location)) tfbs :: _ _ BED _ _) .|
            intersectBed peaks .| sinkFileBed output
        return $ location .~ output $ emptyFile
        )
  where
    getSummit pk = let c = pk^.chromStart + fromJust (pk^.npPeak)
                   in pk & chromStart .~ c - window
                         & chromEnd .~ c + window
