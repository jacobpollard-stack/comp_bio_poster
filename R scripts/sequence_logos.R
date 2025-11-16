library(Biostrings)
library(tidyverse)
library(ggseqlogo)
library(patchwork)
library(ggplot2)
library(msa)
library(reshape2)

# Load aligned MAFFT seqs

alm_pre <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alignments/mafft_pre-2008_swine_NA.fasta')
alm_post <- readAAStringSet('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 1/BABS/Poster/comp_bio_poster/FASTAs/alignments/mafft_post-2008_human_NA.fasta')

# Define sialic acid binding region

start_pos <- 137
end_pos <- 220
start_pos2 <- 148
end_pos2 <- 231

# Function to crop sequences
crop_pre <- subseq(alm_pre, start = start_pos, end = end_pos)
crop_post <- subseq(alm_post, start = start_pos2, end = end_pos2)

# Create position weight matrix

pwm_pre <- consensusMatrix(crop_pre, as.prob=TRUE)
pwm_post <- consensusMatrix(crop_post, as.prob=TRUE)

# Making the x axis

vec <- seq(start_pos, end_pos)

# Plot sequence logos

## Pre-pandemic logo

logo_pre <- ggseqlogo(pwm_pre, method = 'bits') +
  ggtitle('Pre-pandemic') +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 16),
        axis.text.x = element_text(size = 10,
                                   angle = 90),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Bits') +
  scale_x_continuous(
    breaks = seq(1, length(vec), by = 1),
    labels = seq(start_pos, end_pos, by = 1)
  )

## Post-pandemic logo

logo_post <- ggseqlogo(pwm_post, method = 'bits', col_scheme = 'chemistry') +
  ggtitle('Post-pandemic') +
  theme(
    axis.text.y = element_text(hjust = 1, margin = margin(r = 2)),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(size = 10, angle = 90),
    legend.position = 'none') +
  ylab('Bits') +
  scale_x_continuous(
    breaks = seq(1, length(vec), by = 1),
    labels = seq(start_pos2, end_pos2, by = 1)
  )

# Combine plots

combined_logos <- logo_pre / logo_post +
  plot_annotation(
    title = 'HA1 Sequence Logos of Sialic Acid Binding Region Pre- and Post-2009 H1N1pdm09 Pandemic',
    caption = 'Data Source: NCBI Influenza Virus Resource',
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.1),
      plot.title.position = "plot",
      plot.caption = element_text(size = 10)))

combined_logos

# Save combined plot

ggsave('Figures/sequence_logos.png', combined_logos, width = 15, height = 5, dpi = 300)
