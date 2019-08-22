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
    , atacFilterBamSort
    , atacGetBed
    , atacBamToBed
    , atacConcatBed
    , atacCallPeak
    , atacGetNarrowPeak
    ) where

import           Bio.Data.Bed
import           Bio.Data.Bed.Types (BED(..))
import           Bio.Pipeline
import           Data.Bifunctor                (bimap)
import           Data.Coerce                   (coerce)
import           Data.Either                   (lefts, rights)
import qualified Data.Map.Strict               as M
import           Data.Singletons.Prelude.List   (Elem)
import           Data.Singletons               (SingI)
import qualified Data.Text                     as T
import           System.IO.Temp                (withTempFile)

import           Taiji.Pipeline.ATACSeq.Types
import           Taiji.Prelude

type ATACSeqWithSomeFile = ATACSeq N [Either SomeFile (SomeFile, SomeFile)]

atacMkIndex :: ATACSeqConfig config => [a] -> ReaderT config IO ()
atacMkIndex input
    | null input = return ()
    | otherwise = do
        genome <- asks (fromJust . _atacseq_genome_fasta)
        -- Generate BWA index
        dir <- asks (fromJust . _atacseq_bwa_index)
        _ <- liftIO $ bwaMkIndex genome dir
        return ()

atacDownloadData :: ATACSeqConfig config
                 => [ATACSeqWithSomeFile]
                 -> ReaderT config IO [ATACSeqWithSomeFile]
atacDownloadData dat = dat & traverse.replicates.traverse.files.traverse %%~
    (\fl -> do
        dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Download"))
        liftIO $ downloadFiles dir fl )

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
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> case fl of
        Left f ->
            let f' = fromSomeTags f :: File '[] 'Fastq
            in bwaAlign output idx (Left f') $ defaultBWAOpts & bwaCores .~ 4
        Right (f1,f2) ->
            let f1' = fromSomeTags f1 :: File '[] 'Fastq
                f2' = fromSomeTags f2 :: File '[] 'Fastq
            in bwaAlign output idx (Right (f1', f2')) $ defaultBWAOpts & bwaCores .~ 4
        )

-- | Grab bam files
atacGetBam :: [ATACSeqWithSomeFile]
           -> [ ATACSeq S (
                    Either (File '[] 'Bam) (File '[PairedEnd] 'Bam) )]
atacGetBam inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = flip map (filter (\x -> getFileType x == Bam) $ lefts fls) $ \fl ->
        if fl `hasTag` PairedEnd
            then Right $ fromSomeFile fl
            else Left $ fromSomeFile fl

atacFilterBamSort :: ATACSeqConfig config
    => ATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
    -> ReaderT config IO ( ATACSeq S
        (Either (File '[CoordinateSorted] 'Bam) (File '[CoordinateSorted, PairedEnd] 'Bam)) )
atacFilterBamSort input = do
    dir <- asks ((<> "/Bam") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . fun output)
        (fmap Right . fun output)
  where
    fun output x = withTempFile "./" "tmp_file." $ \f _ ->
        filterBam "./" f x >>= sortBam "./" output

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

atacCallPeak :: (ATACSeqConfig config, SingI tags)
             => ATACSeq S (File tags 'Bed)
             -> ReaderT config IO (ATACSeq S (File '[] 'NarrowPeak))
atacCallPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    opts <- asks _atacseq_callpeak_opts
    let output = printf "%s/%s_rep%d.narrowPeak" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO .
        (\fl -> callPeaks output fl Nothing opts)

-- | Fetch narrowpeaks in input data.
atacGetNarrowPeak :: [ATACSeqWithSomeFile]
                  -> [ATACSeq S (File '[] 'NarrowPeak)]
atacGetNarrowPeak input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ map fromSomeFile .
        filter (\x -> getFileType x == NarrowPeak) . lefts