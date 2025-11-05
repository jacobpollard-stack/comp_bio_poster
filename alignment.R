library(msa)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(rentrez)

alignment <- msa('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/all_HA_world.fasta', method = "ClustalW", type = 'protein')
aligned <- msaConvert(alignment, type = "seqinr::alignment")
