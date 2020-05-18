{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.ATACSeq.Types where

import Bio.Pipeline.Utils (Directory)
import Bio.Pipeline.CallPeaks (CallPeakOpts)
import qualified Data.Text as T
import Shelly hiding (FilePath)
import           System.FilePath               (takeDirectory)
import           Bio.Seq.IO
import Control.Exception (catch, SomeException(..))

import Taiji.Prelude

class ATACSeqConfig config where
    _atacseq_output_dir :: config -> Directory
    _atacseq_input :: config -> FilePath
    _atacseq_assembly :: config -> Maybe String
    _atacseq_bwa_index :: config -> Maybe FilePath
    _atacseq_bwa_seed_length :: config -> Int
    _atacseq_genome_fasta :: config -> Maybe FilePath
    _atacseq_genome_index :: config -> Maybe FilePath
    _atacseq_motif_file :: config -> Maybe FilePath
    _atacseq_callpeak_opts :: config -> CallPeakOpts
    _atacseq_annotation :: config -> Maybe FilePath
    _atacseq_tmp_dir :: config -> Maybe FilePath

qcDir :: ATACSeqConfig config => ReaderT config IO FilePath
qcDir = asks _atacseq_output_dir >>= getPath . (<> "/QC/")

getGenomeFasta :: ATACSeqConfig config => ReaderT config IO FilePath
getGenomeFasta = asks _atacseq_genome_fasta >>= \case
    Nothing -> error "genome fasta is missing"
    Just fasta -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack fasta 
       if exist
           then return fasta
           else asks _atacseq_assembly >>= \case
               Nothing -> error "genome fasta is missing"
               Just assembly -> do
                   liftIO $ fetchGenome fasta assembly
                   return fasta

getAnnotation :: ATACSeqConfig config => ReaderT config IO FilePath
getAnnotation = asks _atacseq_annotation >>= \case
    Nothing -> error "annotation is missing"
    Just anno -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack anno
       if exist
           then return anno
           else asks _atacseq_assembly >>= \case
               Nothing -> error "annotation is missing"
               Just assembly -> do
                   liftIO $ fetchAnnotation anno assembly
                   return anno

getMotif :: ATACSeqConfig config => ReaderT config IO FilePath
getMotif = asks _atacseq_motif_file >>= \case
    Nothing -> error "motif file is missing"
    Just motif -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack motif
       if exist
           then return motif
           else asks _atacseq_assembly >>= \case
               Nothing -> error "motif file is missing"
               Just assembly -> do
                   liftIO $ fetchMotif motif assembly
                   return motif

mkGenomeIndex :: ATACSeqConfig config => ReaderT config IO ()
mkGenomeIndex = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _atacseq_genome_index )
    exist <- liftIO $ catch (withGenome seqIndex $ const $ return True) $
        \(SomeException _) -> return False
    unless exist $ do
        liftIO $ putStrLn "Generating genome index ..."
        genome <- getGenomeFasta
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
        liftIO $ mkIndex [genome] seqIndex
{-# INLINE mkGenomeIndex #-}

getCallPeakOpt :: ATACSeqConfig config => ReaderT config IO CallPeakOpts
getCallPeakOpt = do
    opt <- asks _atacseq_callpeak_opts 
    s <- case opt^.gSize of
        Nothing -> do
            idx <- asks ( fromMaybe (error "Genome index file was not specified!") .
                _atacseq_genome_index )
            s <- liftIO $ fmap fromIntegral $ withGenome idx $
                return . foldl1' (+) . map snd . getChrSizes
            return $ Just $ show (truncate $ 0.9 * (s :: Double) :: Int)
        x -> return x
    return $ gSize .~ s $ opt