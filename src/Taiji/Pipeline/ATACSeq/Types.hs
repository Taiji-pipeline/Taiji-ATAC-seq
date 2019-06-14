module Taiji.Pipeline.ATACSeq.Types where

import Bio.Pipeline.Utils (Directory)
import Bio.Pipeline.CallPeaks (CallPeakOpts)

class ATACSeqConfig config where
    _atacseq_output_dir :: config -> Directory
    _atacseq_input :: config -> FilePath
    _atacseq_bwa_index :: config -> Maybe FilePath
    _atacseq_genome_fasta :: config -> Maybe FilePath
    _atacseq_genome_index :: config -> Maybe FilePath
    _atacseq_motif_file :: config -> Maybe FilePath
    _atacseq_callpeak_opts :: config -> CallPeakOpts
    _atacseq_annotation :: config -> Maybe FilePath
