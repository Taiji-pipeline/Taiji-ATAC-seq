{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts      #-}
module Taiji.Pipeline.ATACSeq.Functions.GeneQuant
    ( estimateExpr
    , mkTable
    ) where

import qualified Data.ByteString.Char8 as B
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Bio.RealWorld.GENCODE (Gene(..))
import qualified Data.HashMap.Strict as M
import qualified Data.HashSet as S
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import           Bio.Pipeline.Utils
import qualified Data.Text as T
import Bio.Data.Bed.Utils

import Taiji.Prelude
import Taiji.Utils (readGenesValidated)
import           Taiji.Pipeline.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

-- | Estimate gene expression level by atac-seq signal.
estimateExpr :: ATACSeqConfig config
             => ATACSeq S (File '[Gzip] 'Bed)
             -> ReaderT config IO (ATACSeq S (File '[GeneQuant] 'Tsv))
estimateExpr input = do
    genes <- getAnnotation >>= liftIO . readGenesValidated
    dir <- asks ((<> "/GeneQuant") . _atacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_gene_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        counts <- runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            mkGeneCount genes
        B.writeFile output $ B.unlines $
            map (\(n, c) -> n <> "\t" <> toShortest c) counts
        return $ location .~ output $ emptyFile )

-- | Count number of tags in genes' promoters, normalized by total read counts.
-- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6385419/
mkGeneCount :: PrimMonad m
            => [Gene]
            -> ConduitT BED o m [(B.ByteString, Double)]
mkGeneCount genes = zip (map (original . geneName) genes) . U.toList <$>
    rpkmBed (map getTss genes)
  where
    getTss Gene{..} | geneStrand = BED3 geneChrom (geneLeft - 1000) (geneLeft + 1000)
                    | otherwise = BED3 geneChrom (geneRight - 1000) (geneRight + 1000)
{-# INLINE mkGeneCount #-}

-- | Combine RNA expression data into a table and output
mkTable :: ATACSeqConfig config
        => [ATACSeq S (File '[GeneQuant] 'Tsv)]
        -> ReaderT config IO (Maybe (File '[] 'Tsv))
mkTable [] = return Nothing
mkTable input = do
    outdir <- asks _atacseq_output_dir >>= getPath
    let output = outdir ++ "/GeneQuant/expression_profile.tsv"
    liftIO $ do
        geneNames <- getGeneName input
        (sampleNames, vals) <- unzip <$> mapM (f geneNames) input
        DF.writeTable output (T.pack . show) $ DF.mkDataFrame
            (map (T.pack . B.unpack . original) geneNames) sampleNames $
            transpose vals
        return $ Just $ location .~ output $ emptyFile
  where
    f genes fl = do
        (nm, vals) <- readExpr fl
        return $! (nm, map (\g -> M.lookupDefault 0.01 g vals) genes)

--------------------------------------------------------------------------------
-- Auxiliary functions
--------------------------------------------------------------------------------

readExpr :: ATACSeq S (File '[GeneQuant] 'Tsv)
         -> IO (T.Text, M.HashMap (CI B.ByteString) Double)
readExpr input = do
    c <- B.readFile $ input^.replicates._2.files.location
    return ( fromJust $ input^.groupName
           , M.fromListWith max $ map f $ B.lines c )
  where
    f xs = let fs = B.split '\t' xs in (mk $ head fs, readDouble $ fs!!1)
{-# INLINE readExpr #-}

getGeneName :: [ATACSeq S (File '[GeneQuant] 'Tsv)]
            -> IO [CI B.ByteString]
getGeneName = fmap S.toList . foldM f S.empty
  where
    f names x = foldl' (flip S.insert) names . map (mk . head . B.split '\t') .
        tail . B.lines <$> B.readFile (x^.replicates._2.files.location)
{-# INLINE getGeneName #-}