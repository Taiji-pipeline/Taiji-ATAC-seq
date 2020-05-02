{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Taiji.Pipeline.ATACSeq.Functions.Motif
    ( atacMergePeaks
    , atacFindMotifSiteAll
    , atacGetMotifSite
    ) where

import           Bio.Data.Bed
import           Bio.Data.Bed.Utils  (scanMotif, mkCutoffMotif)
import           Bio.Data.Experiment
import           Bio.Motif                     hiding (score)
import           Bio.Pipeline.Instances        ()
import           Bio.Pipeline.Utils
import           Bio.Seq.IO
import           Data.Default
import qualified Data.Text                     as T
import           System.IO.Temp                (emptyTempFile)

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Prelude

-- | Merge overlapping peaks.
atacMergePeaks :: ATACSeqConfig config
               => [ATACSeq S (Either (File '[] 'NarrowPeak) (File '[Gzip] 'NarrowPeak))]
               -> ReaderT config IO (Maybe (File '[] 'Bed))
atacMergePeaks [] = return Nothing
atacMergePeaks input = do
    dir <- asks _atacseq_output_dir >>= getPath
    let output = dir ++ "/openChromatin.bed"
    liftIO $ do
        openChromatin <- foldM f [] $
            input^..folded.replicates.folded.files
        writeBed output openChromatin
        return $ Just $ location .~ output $ emptyFile
  where
    f acc (Left fl) = do
        peaks <- readBed (fl^.location) :: IO [BED3]
        runResourceT $ runConduit $ mergeBed (acc ++ peaks) .| sinkList
    f acc (Right fl) = do
        peaks <- runResourceT $ runConduit $ streamBedGzip (fl^.location) .| sinkList :: IO [BED3]
        runResourceT $ runConduit $ mergeBed (acc ++ peaks) .| sinkList

atacFindMotifSiteAll :: ATACSeqConfig config
                     => Double     -- ^ p value
                     -> (File '[] 'Bed, [Motif])
                     -> ReaderT config IO (File '[] 'Bed)
atacFindMotifSiteAll p (openChromatin, motifs) = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _atacseq_genome_index )
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
                 -> ( [File '[] 'Bed]
                    , ATACSeq S (Either (File '[] 'NarrowPeak) (File '[Gzip] 'NarrowPeak)) )
                 -> ReaderT config IO (ATACSeq S (File '[] 'Bed))
atacGetMotifSite window (tfbs, e) = do
    dir <- asks ((<> "/TFBS") . _atacseq_output_dir) >>= getPath
    e & replicates.traversed.files %%~ ( \fl -> liftIO $ do
        let output = printf "%s/%s_rep%d.bed" dir (T.unpack $ e^.eid)
                (e^.replicates._1)
        peaks <- runResourceT $ runConduit $ getBed fl .| mapC getSummit .| sinkList
        runResourceT $ runConduit $
            (mapM_ (streamBed . (^.location)) tfbs :: _ _ BED _ _) .|
            intersectBed peaks .| sinkFileBed output
        return $ location .~ output $ emptyFile
        )
  where
    getBed (Left fl) = streamBed $ fl^.location
    getBed (Right fl) = streamBedGzip $ fl^.location
    getSummit pk = pk & chromStart .~ center - window & chromEnd .~ center + window
      where
        center = case pk^.npPeak of
            Nothing -> (pk^.chromStart + pk^.chromEnd) `div` 2
            Just c -> pk^.chromStart + c

