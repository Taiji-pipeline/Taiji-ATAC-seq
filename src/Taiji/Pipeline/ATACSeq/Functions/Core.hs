{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE TemplateHaskell       #-}
module Taiji.Pipeline.ATACSeq.Functions.Core
    ( atacMkIndex
    , atacDownloadData
    , atacGetFastq
    , atacAlign
    , atacGetBam
    , atacGetFilteredBam
    , atacFilterBamSort
    , atacGetBed
    , atacBamToBed
    , atacCallPeak
    , atacGetNarrowPeak
    ) where

import           Bio.Data.Bed hiding (NarrowPeak)
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Either                   (lefts)
import Bio.Seq.IO (withGenome, getChrSizes)
import           Data.List.Singletons (Elem)
import qualified Data.Text                     as T
import Shelly hiding (FilePath)

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Prelude

type ATACSeqWithSomeFile = ATACSeq N [Either SomeFile (SomeFile, SomeFile)]

atacDownloadData :: ATACSeqConfig config
                 => ATACSeqWithSomeFile
                 -> ReaderT config IO ATACSeqWithSomeFile
atacDownloadData dat = do
    tmp <- fromMaybe "./" <$> asks _atacseq_tmp_dir
    dat & replicates.traverse.files.traverse %%~ (\fl -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Download"))
        liftIO $ downloadFiles dir tmp fl )

atacMkIndex :: ATACSeqConfig config
            => [ATACSeqWithSomeFile]
            -> ReaderT config IO [ATACSeqWithSomeFile]
atacMkIndex [] = return []
atacMkIndex input = do
    mkGenomeIndex 
    let fq = filter (isFq . either id fst) $ concat $
            input^..folded.replicates.folded.files
    unless (null fq) $ do
        genome <- getGenomeFasta
        -- Generate BWA index
        dir <- asks (fromJust . _atacseq_bwa_index)
        liftIO (bwaMkIndex genome dir) >> return ()
    return input
  where
    isFq x = getFileType x == Fastq

atacGetFastq :: [ATACSeqWithSomeFile]
             -> [ ATACSeq S ( Either
                    (SomeTags 'Fastq) (SomeTags 'Fastq, SomeTags 'Fastq) )]
atacGetFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile (bimap castFile castFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

atacAlign :: ATACSeqConfig config
          => ATACSeq S ( Either
                (SomeTags 'Fastq) (SomeTags 'Fastq, SomeTags 'Fastq) )
          -> ReaderT config IO ( ATACSeq S
                (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam)) )
atacAlign input = do
    dir <- asks ((<> "/Bam") . _atacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _atacseq_bwa_index)
    seedLen <- asks _atacseq_bwa_seed_length
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        opt = defaultBWAOpts & bwaCores .~ 4 & bwaSeedLen .~ seedLen
    input & replicates.traverse.files %%~ liftIO . ( \fl -> case fl of
        Left f ->
            let f' = fromSomeTags f :: File '[] 'Fastq
            in bwaAlign output idx (Left f') opt
        Right (f1,f2) ->
            let f1' = fromSomeTags f1 :: File '[] 'Fastq
                f2' = fromSomeTags f2 :: File '[] 'Fastq
            in bwaAlign output idx (Right (f1', f2')) opt
        )

-- | Grab bam files
atacGetBam :: [ATACSeqWithSomeFile]
           -> [ ATACSeq S (
                    Either (File '[] 'Bam) (File '[PairedEnd] 'Bam) )]
atacGetBam inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filterData $ lefts fls) $ \fl ->
        if fl `hasTag` PairedEnd
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl
    filterData = filter $ \x -> getFileType x == Bam && not (x `hasTag` Filtered)

atacGetFilteredBam :: Elem 'PairedEnd b ~ 'True
                   => [ATACSeqWithSomeFile]
                   -> [ ATACSeq S (
                       Either (File a 'Bam) (File b 'Bam) )]
atacGetFilteredBam inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filterData $ lefts fls) $ \fl ->
        if fl `hasTag` PairedEnd
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl
    filterData = filter $ \x -> getFileType x == Bam && x `hasTag` Filtered

atacFilterBamSort :: ATACSeqConfig config
    => ATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
    -> ReaderT config IO ( ATACSeq S
        (Either (File '[CoordinateSorted] 'Bam) (File '[CoordinateSorted, PairedEnd] 'Bam)) )
atacFilterBamSort input = do
    dir <- asks ((<> "/Bam") . _atacseq_output_dir) >>= getPath
    tmp <- fromMaybe "./" <$> asks _atacseq_tmp_dir
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . filterBamSort tmp output)
        (fmap Right . filterBamSort tmp output)

atacGetBed :: [ATACSeqWithSomeFile]
           -> [ATACSeq S (Either (File '[Gzip] 'Bed) (File '[PairedEnd, Gzip] 'Bed))]
atacGetBed input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~
        map f . filter (\x -> getFileType x == Bed && x `hasTag` Gzip) . lefts
  where
    f fl = if fl `hasTag` PairedEnd
        then Right $ fromSomeFile fl
        else Left $ fromSomeFile fl

atacBamToBed :: ATACSeqConfig config
             => ATACSeq S ( Either (File tags1 'Bam) (File tags2 'Bam) )
             -> ReaderT config IO (ATACSeq S (File '[Gzip] 'Bed))
atacBamToBed input = do
    dir <- asks ((<> "/Bed") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        let fl' = either coerce coerce fl :: File '[] 'Bam
        bam2Bed output (\x -> not $ (x^.chrom) `elem` ["chrM", "M"]) fl' )

atacCallPeak :: ATACSeqConfig config
             => ATACSeq S (File '[Gzip] 'Bed)
             -> ReaderT config IO (ATACSeq S (File '[Gzip] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    tmpdir <- fromMaybe "./" <$> asks _atacseq_tmp_dir
    opts <- getCallPeakOpt
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _atacseq_genome_index )
    let output = printf "%s/%s_rep%d.narrowPeak.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        --outputBW = printf "%s/%s_rep%d_signal.bw" dir (T.unpack $ input^.eid)
        --    (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        peak <- callPeaks output fl Nothing $ genSignal.~ True $ opts
        {-
        chrSize <- withGenome seqIndex $ return . getChrSizes
        bedGraphToBigWig outputBW chrSize [] tmpdir $ output <> ".bdg"
        shelly $ rm_f $ fromText $ T.pack $ output <> ".bdg"
        -}
        return peak
        )

-- | Fetch narrowpeaks in input data.
atacGetNarrowPeak :: [ATACSeqWithSomeFile]
                  -> [ATACSeq S (Either (File '[] 'NarrowPeak) (File '[Gzip] 'NarrowPeak))]
atacGetNarrowPeak input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ map f .
        filter (\x -> getFileType x == NarrowPeak) . lefts
  where
    f fl = if fl `hasTag` Gzip
        then Right $ fromSomeFile fl
        else Left $ fromSomeFile fl