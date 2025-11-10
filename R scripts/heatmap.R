library(Biostrings)
library(tidyverse)
library(ggseqlogo)
library(patchwork)
library(ggplot2)
library(msa)

# Load MAFFT alignments

pre2008_mafft_alm <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_pre-2008_swine_NA.fasta'
post2008_mafft_alm <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_post-2008_human_NA.fasta'

# Make them msa objects

seqs_pre2008 <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_pre-2008_swine_NA.fasta')
seqs_post2008 <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_post-2008_human_NA.fasta')

# Define cropping positions

start_pos <- 137
end_pos <- 220
start_pos1 <- 137
end_pos1 <- 220
start_pos2 <- 149
end_pos2 <- 232

# Defining a function to crop and split sequences into lists

crop_and_split <- function(seqs, start, end) {
  seqs_vec <- as.character(seqs)
  alignment_list <- lapply(seqs_vec, function(s) {
    chars <- strsplit(s, "")[[1]]
    if(length(chars) < end) chars <- c(chars, rep("-", end - length(chars)))
    chars[start:end]
  })
  alignment_list
}

# Crop sequences and split

alignment_pre_list <- crop_and_split(seqs_pre2008, start_pos1, end_pos1)
alignment_post_list <- crop_and_split(seqs_post2008, start_pos2, end_pos2)

# Convert lists to matrices

alignment_pre_mat <- do.call(rbind, alignment_pre_list)
alignment_post_mat <- do.call(rbind, alignment_post_list)

# Make consensus seqs

cons_pre <- apply(alignment_pre_mat, 2, function(col) names(sort(table(col), decreasing=TRUE))[1])
cons_post <- apply(alignment_post_mat, 2, function(col) names(sort(table(col), decreasing=TRUE))[1])

# Convert character vectors to single strings

cons_pre_seq <- paste(cons_pre, collapse = "")
cons_post_seq <- paste(cons_post, collapse = "")

# Write to FASTA file

write_fasta <- function(seq, name, file) {
  cat(paste0(">", name, "\n", seq, "\n"), file = file)
}

# Save each consensus sequence

write_fasta(cons_pre_seq, "consensus_pre", "consensus_pre.fasta")
write_fasta(cons_post_seq, "consensus_post", "consensus_post.fasta")

# Define data frame for plotting

seq_df <- data.frame(
  Position = rep(start_pos:end_pos, 2),
  Sequence = rep(c("Pre-2009", "Post-2009"), each = end_pos-start_pos+1),
  AA = c(cons_pre, cons_post)
)

# Amino acid colours

aa_type <- c(
  A = "Hydrophobic", V = "Hydrophobic", L = "Hydrophobic", I = "Hydrophobic",
  M = "Hydrophobic", F = "Hydrophobic", W = "Hydrophobic", P = "Hydrophobic",
  G = "Special",
  S = "Polar", T = "Polar", C = "Polar", N = "Polar", Q = "Polar",
  D = "Acidic", E = "Acidic",
  K = "Basic", R = "Basic", H = "Basic",
  "-" = "Gap"
)

seq_df <- seq_df %>%
  mutate(Type = aa_type[AA])

# Define colours by type

type_colours <- c(
  Hydrophobic = "#1f77b4",
  Polar       = "#2ca02c",
  Acidic      = "#d62728",
  Basic       = "#ff7f0e",
  Special     = "#7f7BBB",
  Gap         = "#ffffff"
)

# Heatmap with letters and legend for type

heatmap_plot <- ggplot(seq_df, aes(x = Position, y = Sequence, fill = Type)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = AA), size = 3, fontface = "bold") +
  scale_fill_manual(values = type_colours) +
  labs(title = "HA1 Consensus Alignment (137–220) Pre- and Post-2009 Pandemic",
       x = "Amino Acid Position", y = "Sequence",
       fill = "AA Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

# Print plot

heatmap_plot

# Save heatmap

ggsave("Figures/HA1_consensus_heatmap.png", heatmap_plot, width = 15, height = 1.5, dpi = 300)

