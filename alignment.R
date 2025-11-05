library(msa)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(rentrez)

alignment <- msa("/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/all_HA.fasta", method = "ClustalW", type = 'protein')  # or method="Muscle", "ClustalOmega", "MAFFT"
aligned <- msaConvert(alignment, type = "seqinr::alignment")
