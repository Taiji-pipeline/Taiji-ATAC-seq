{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.ATACSeq.Functions.GeneQuant
    ( estimateExpr
    , mkTable
    ) where

import qualified Data.ByteString.Char8 as B
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.Double.Conversion.ByteString (toShortest)
import qualified Data.HashMap.Strict as M
import Bio.Utils.Misc (readDouble)
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)
import qualified Data.Text as T
import Bio.Data.Bed.Utils

import Taiji.Prelude
import           Taiji.Pipeline.ATACSeq.Types

-- | Estimate gene expression level by atac-seq signal.
estimateExpr :: ATACSeqConfig config
             => ATACSeq S (File '[Gzip] 'Bed)
             -> ReaderT config IO (ATACSeq S (File '[GeneQuant] 'Tsv))
estimateExpr input = do
    genes <- asks _atacseq_annotation >>= liftIO . readGenes . fromJust
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

mkTable :: ATACSeqConfig config
        => [ATACSeq S (File '[GeneQuant] 'Tsv)]
        -> ReaderT config IO (Maybe (File '[] 'Tsv))
mkTable es = do
    results <- liftIO $ readExpr es
    if null results
        then return Nothing
        else do
            outdir <- asks _atacseq_output_dir >>= getPath
            let output = outdir ++ "/GeneQuant/expression_profile.tsv"
                (expNames, values) = unzip $ M.toList $
                    fmap (map (second average) . combine) $ M.fromListWith (++) $ results
            liftIO $ B.writeFile output $ B.unlines $
                B.pack (T.unpack $ T.intercalate "\t" $ "Name" : expNames) :
                map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
                    (combine values)
            return $ Just $ location .~ output $ emptyFile

--------------------------------------------------------------------------------
-- Auxiliary functions
--------------------------------------------------------------------------------

readExpr :: [ATACSeq S (File '[GeneQuant] 'Tsv)]
         -> IO [(T.Text, [[(CI B.ByteString, Double)]])]
readExpr quantifications = forM quantifications $ \e -> do
    c <- B.readFile $ e^.replicates._2.files.location
    let result = map (\xs ->
            let fs = B.split '\t' xs in (mk $ head fs, readDouble $ fs!!1)) $
            tail $ B.lines c
    return (fromJust $ e^.groupName, [result])
{-# INLINE readExpr #-}

combine :: [[(CI B.ByteString, Double)]] -> [(CI B.ByteString, [Double])]
combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
  where
    names = nubSort $ concatMap (fst . unzip) xs
    xs' = map (fmap average . M.fromListWith (++) . map (second return)) xs
{-# INLINE combine #-}

average :: [Double] -> Double
average [a,b]   = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs      = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}