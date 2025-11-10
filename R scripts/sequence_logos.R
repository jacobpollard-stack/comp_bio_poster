library(Biostrings)
library(tidyverse)
library(ggseqlogo)
library(patchwork)
library(ggplot2)
library(msa)

# Load aligned MAFFT seqs

pre <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alignments/mafft_pre-2008_swine_NA.fasta')
post <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alignments/mafft_post-2008_human_NA.fasta')

# Define sialic acid binding region

start_pos <- 137
end_pos <- 220

# Function to crop sequences

crop_pre <- subseq(pre, start = start_pos, end = end_pos)
crop_post <- subseq(post, start = start_pos, end = end_pos)

#


