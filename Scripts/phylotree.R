library(tidyverse)
library(lubridate)
library(Biostrings)
library(ape)
library(phangorn)
library(ggtree)
library(dplyr)

# Load the data
ha_raw <- readAAStringSet("Files/combined.fa")

# Make this data have accession-only IDs
names(ha_raw) <- sapply(strsplit(names(ha_raw), " "), `[`, 1)

# Check uniqueness
stopifnot(!any(duplicated(names(ha_raw))))

# Load metadata
meta <- read_csv("Files/flu.txt")

# Filter metadata to only include sequences present in ha_raw
meta <- meta %>%
  filter(accession %in% names(ha_raw))

# Ensure the order of metadata matches the order of sequences in ha_raw
meta <- meta %>%
  arrange(match(accession, names(ha_raw)))

# Filter to full-length HA
meta_clean <- meta %>%
  filter(length >= 560)

# Apply same filter to sequences
ha_clean <- ha_raw[names(ha_raw) %in% meta_clean$accession]

# Clean host names and extract years
meta_clean <- meta_clean %>%
  mutate(
    host = tolower(host),
    host = if_else(host == "human", "human", "swine"),
    year = case_when(
      str_length(date) == 4 ~ as.integer(date),
      TRUE ~ year(ymd(date))
    )
  ) %>%
  filter(!is.na(year))

# Temporal thinning
set.seed(42)

meta_thin <- meta_clean %>%
  group_by(year, host) %>%
  group_modify(~ slice_sample(.x, n = min(2, nrow(.x)))) %>%
  ungroup()

# Subset sequences using thinned metadata
ha_thin <- ha_clean[names(ha_clean) %in% meta_thin$accession]

# Reorder
ha_thin <- ha_thin[match(meta_thin$accession, names(ha_thin))]

# Write FASTA
writeXStringSet(
  ha_thin,
  filepath = "Files/ha_thinned.fa",
  format = "fasta"
)


# Load aligned sequences
ha_aligned <- readAAStringSet("Files/alignments/thinned_alm.fasta")

# Convert to phyDat object
ha_phyDat <- phyDat(as.list(ha_aligned), type = "AA")

# NJ starting tree
dm <- dist.ml(ha_phyDat)
treeNJ <- NJ(dm)

# ML fit
fit <- pml(treeNJ, data = ha_phyDat)

fitML <- optim.pml(
  fit,
  model = "JTT",
  optGamma = TRUE,
  rearrangement = "stochastic"
)

# Plot tree coloured by host
meta_tree <- meta_thin %>%
  dplyr::select(accession, host, year) %>%
  dplyr::rename(label = accession)

# Update metadata to match pruned tree
meta_tree_pruned <- meta_tree %>%
  filter(label %in% tree_pruned$tip.label)

# Add bootstrap values
set.seed(42)

bs <- bootstrap.pml(
  fitML,
  bs = 100,
  optNni = TRUE,
  control = pml.control(trace = 0)
)

# Bootstrap
tree_bs <- plotBS(fitML$tree, bs, p = 0, type = "none")

# Prune the bootstrapped tree
tree_pruned <- drop.tip(tree_bs, c("QCT08133", "QCT08134"))

# Metadata must match pruned tips
meta_tree_pruned <- meta_thin %>%
  transmute(label = accession, host, year, name, date) %>%
  filter(label %in% tree_pruned$tip.label)

meta_tree_pruned <- meta_tree_pruned %>%
  mutate(strain = str_extract(name, "A/.+?\\)"))

# Plot with bootstrap labels on internal nodes
p_pruned <- ggtree(tree_pruned) %<+% meta_tree_pruned +
  geom_tree(aes(colour = host), size = 1) +
  geom_tiplab(aes(label = strain, colour = host), size = 3) +
  scale_color_manual(values = c(human = "#1f78b4", swine = "#e31a1c")) +
  theme_tree2()

p_pruned

# Save the plot
ggsave("Figures/ha_phylo_tree_pruned.png", p_pruned, width = 10, height = 12, dpi = 300)
