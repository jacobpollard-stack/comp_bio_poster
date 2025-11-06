library(Biostrings)
library(tidyverse)
library(ggseqlogo)
library(patchwork)
library(ggplot2)

# --- Paths to MAFFT alignments ---

pre2008_mafft_alm <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_pre-2008:08_human_north-america.fasta'
post2008_mafft_alm <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alm/mafft_alm_post-2008.fasta'

# --- Load aligned sequences ---

seqs_pre2008 <- readAAStringSet(pre2008_mafft_alm)
seqs_post2008 <- readAAStringSet(post2008_mafft_alm)

# --- Define crop region ---

start_pos <- 137
end_pos <- 220

# --- Function to crop sequences and split into characters ---

crop_and_split <- function(seqs, start, end) {
  seqs_vec <- as.character(seqs)
  alignment_list <- lapply(seqs_vec, function(s) {
    chars <- strsplit(s, "")[[1]]
    if(length(chars) < end) chars <- c(chars, rep("-", end - length(chars)))
    chars[start:end]
  })
  alignment_list
}

# --- Crop sequences ---

alignment_pre_list <- crop_and_split(seqs_pre2008, start_pos, end_pos)
alignment_post_list <- crop_and_split(seqs_post2008, start_pos, end_pos)

# --- Filter out sequences that are all gaps ---

filter_nonempty <- function(alignment_list) {
  alignment_list[sapply(alignment_list, function(x) any(x != "-"))]
}

alignment_pre_list <- filter_nonempty(alignment_pre_list)
alignment_post_list <- filter_nonempty(alignment_post_list)

# --- Convert to matrices for consensus ---

alignment_pre_mat <- do.call(rbind, alignment_pre_list)
alignment_post_mat <- do.call(rbind, alignment_post_list)

# --- Consensus sequences ---

cons_pre <- apply(alignment_pre_mat, 2, function(col) names(sort(table(col), decreasing=TRUE))[1])
cons_post <- apply(alignment_post_mat, 2, function(col) names(sort(table(col), decreasing=TRUE))[1])

# --- Sequence logo plots ---

pre_logo <- ggseqlogo(alignment_pre_list, method = "prob", seq_type = "aa") +
  ggtitle("H1N1 HA1: Sialic acid binding region pre-2009 pandemic") +
  scale_x_continuous(breaks = 1:(end_pos-start_pos+1), labels = start_pos:end_pos) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        plot.title = element_text(face = "bold", size = 14))

post_logo <- ggseqlogo(alignment_post_list, method = "prob", seq_type = "aa") +
  ggtitle("H1N1 HA1: Sialic acid binding region post-2009 pandemic") +
  scale_x_continuous(breaks = 1:(end_pos-start_pos+1), labels = start_pos:end_pos) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        plot.title = element_text(face = "bold", size = 14))

combined_logos <- pre_logo / post_logo +
  plot_annotation(title = "Comparison of H1N1 HA1 Sialic Acid Binding Region Pre- and Post-2009 Pandemic",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
combined_logos

# --- Alignment heatmap with letters ---

seq_df <- data.frame(
  Position = rep(start_pos:end_pos, 2),
  Sequence = rep(c("Pre-2009", "Post-2009"), each = end_pos-start_pos+1),
  AA = c(cons_pre, cons_post)
)

# Define amino acid colours

aa_colors <- c(
  A = "#1f77b4", V = "#1f77b4", L = "#1f77b4", I = "#1f77b4", M = "#1f77b4", F = "#1f77b4", W = "#1f77b4", P = "#1f77b4",
  G = "#7f7f7f",
  S = "#2ca02c", T = "#2ca02c", C = "#2ca02c", N = "#2ca02c", Q = "#2ca02c",
  D = "#d62728", E = "#d62728",
  K = "#ff7f0e", R = "#ff7f0e", H = "#ff7f0e",
  "-" = "#ffffff"
)

heatmap_plot <- ggplot(seq_df, aes(x = Position, y = Sequence, fill = AA)) +
  geom_tile(color = "white") +
  geom_text(aes(label = AA), size = 3, fontface = "bold") +
  scale_fill_manual(values = aa_colors) +
  labs(title = "HA1 Consensus Alignment (137–220)",
       x = "Amino Acid Position", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

# --- Display plots ---

combined_logos
heatmap_plot
