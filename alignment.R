library(msa)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(rentrez)
library(ggmsa)

# pre-2008

## Load in pre-2008 aligned sequences
pre2008_mafft_alm <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_alm_pre-2008.fasta'

## Convert to MSA object
seqs_pre2008 <- readAAStringSet(pre2008_mafft_alm)
msaobj_pre2008 <- new("MsaAAMultipleAlignment",
                      unmasked = seqs_pre2008)

## Create consensus sequence
cons <- msaConsensusSequence(msaobj_pre2008)
cons

ggmsa(aligned_fasta,
      color = "Clustal",   
      consensus_views = TRUE, 
      seq_name = TRUE)

ggsave("/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/figures/alignment-pre-2008:08-all-species-north_america.pdf",
       width = 12, height = 8)
