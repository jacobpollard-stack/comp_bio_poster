library(ggplot2)
library(grid)

# Data frame with start and end positions of gene features

df <- data.frame(
 start = c(11,64,115,133,153,155,183,185,189,192,219,225,263,276,278,324,329,489,509,514,545),
 end  = c(64,115,263,138,153,156,184,186,190,194,220,228,276,329,281,329,489,551,514,542,551),
 label = c("Structural domain and glycosylation site", 'Esterase domain', 'Receptor binding domain', 'Sialic acid binding site', 'Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Esterase domain', 'Structural domain and glycosylation site', 'Glycosylation site only found in H1N1pdm09','Cleavage site', 'Structural domain and glycosylation site', 'Structural domain and glycosylation site','Sialic acid binding site','Transmembrane','Sialic acid binding site')
)

# Filtered data

df_filtered <- df[!df$start %in% c(155, 185, 278), ]

# Define colour palette

cb_palette <- c(
 "#E69F00", "#56B4E9", "#009E73", "#F0E442",
 "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)
labels_unique <- unique(df_filtered$label)
fill_colors <- setNames(cb_palette[1:length(labels_unique)], labels_unique)

# Define feature ranges

feature_ranges <- data.frame(
 start = c(11, 64, 115, 263, 276, 489),
 end  = c(64, 115, 263, 276, 489, 551),
 label = c("Fusion peptide", "Esterase domain", "Receptor binding domain", 
      "Esterase domain", "Fusion peptide", "Structural domain")
)

# Define HA1/HA2 ranges

ha_ranges <- data.frame(
 start = c(11, 329),
 end  = c(329, 551),
 label = c("HA1", "HA2")
)

# Create plot

plot <- ggplot(df_filtered) +
 geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = label),
      colour = "black") +
 geom_text(aes(x = start, y = 1.25, label = start),
      size = 2, angle = 90, hjust = 1) +
 
 scale_fill_manual(name = NULL, values = fill_colors) +
 geom_segment(
  data = feature_ranges,
  aes(x = start, xend = end, y = -0.30, yend = -0.30),
  arrow = arrow(length = unit(0.18, "cm"), ends = "both", type = "closed"),
  size = 1, color = "black"
 ) +
 geom_text(
  data = feature_ranges,
  aes(x = (start+end)/2, y = -0.55, label = label),
  size = 3
 ) +
 geom_segment(
  data = ha_ranges,
  aes(x = start, xend = end, y = -0.85, yend = -0.85),
  arrow = arrow(length = unit(0.18, "cm"), ends = "both", type = "closed"),
  size = 1, color = "black"
 ) +
 geom_text(
  data = ha_ranges,
  aes(x = (start+end)/2, y = -1.10, label = label),
  size = 3
 ) +
 theme_void() +
 theme(
  plot.title = element_text(
   hjust = 0.065, vjust = 2, size = 16, face = "bold"
  ),
  legend.position = "bottom",
  legend.box.margin = margin(t = 10, b = 10),
  plot.margin = margin(10, 10, 10, 10)
 ) +
 ggtitle("HA0 Protein Map of H1N1 consensus sequence") +
 guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Save plot
ggsave(
 filename = "HA0_protein_map_H1N1_consensus_auto_spacing.jpg",
 plot = plot,
 width = 15,
 height = 3,
 dpi = 300
)
