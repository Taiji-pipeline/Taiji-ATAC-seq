module Taiji.Pipeline.ATACSeq.Config where

import Bio.Pipeline.Utils (Directory)

class ATACSeqConfig config where
    _atacseq_output_dir :: config -> Directory
    _atacseq_input :: config -> FilePath
    _atacseq_picard :: config -> Maybe FilePath
    _atacseq_bwa_index :: config -> Maybe FilePath
    _atacseq_genome_fasta :: config -> Maybe FilePath
    _atacseq_genome_index :: config -> Maybe FilePath
    _atacseq_motif_file :: config -> Maybe FilePath
