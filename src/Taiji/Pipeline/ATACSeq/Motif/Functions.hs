{-# LANGUAGE DataKinds        #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.ATACSeq.Motif.Functions
    ( atacMergePeaks
    , atacFindMotifSiteAll
    , atacGetMotifSite
    ) where

import           Bio.Data.Bed                  (BED (..), BED3 (..),
                                                BEDLike (..), getMotifPValue,
                                                getMotifScore, intersectBedWith,
                                                mergeBed, motifScan, readBed,
                                                readBed', writeBed, _npEnd,
                                                _npPeak, _npPvalue, _npStart)
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
import           Data.Maybe                    (fromJust, isJust)
import           Data.Monoid                   ((<>))
import qualified Data.Text                     as T
import           Scientific.Workflow
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)
import           System.IO
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
    -- Generate sequence index
    genome <- asks (fromJust . _atacseq_genome_fasta)
    seqIndex <- asks (fromJust . _atacseq_genome_index)
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
        (readBed (openChromatin^.location) :: Source IO BED3) =$=
            motifScan g motifs def p =$= getMotifScore g motifs def =$=
            getMotifPValue (Just (1 - p * 10)) motifs def $$ writeBed output
        return $ location .~ output $ emptyFile
  where
    p = 1e-5

-- | Retrieve TFBS for each experiment
atacGetMotifSite :: ATACSeqConfig config
                 => Int -- ^ region around summit
                 -> ([File '[] 'Bed], [ATACSeq S (File '[] 'NarrowPeak)])
                 -> WorkflowConfig config [ATACSeq S (File '[] 'Bed)]
atacGetMotifSite window (tfbs, experiment) = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/TFBS/"))
    mapM (mapFileWithDefName dir ".bed" fun) experiment
  where
    fun output fl = liftIO $ do
        peaks <- readBed (fl^.location) =$= mapC getSummit $$ sinkList
        (mapM_ (readBed . (^.location)) tfbs :: Source IO BED) =$=
            intersectBedWith getPvalue peaks =$= filterC (isJust . snd) =$=
            mapC (\(bed, p) -> bed{_score = p}) $$ writeBed output
        return $ location .~ output $ emptyFile
    getSummit pk = let c = chromStart pk + (fromJust . _npPeak) pk
                   in pk { _npStart = c - window
                         , _npEnd = c + window }
    getPvalue [] = Nothing
    getPvalue xs = Just $ maximum $ map (fromJust . _npPvalue) xs
