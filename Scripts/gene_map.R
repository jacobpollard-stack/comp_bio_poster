library(ggplot2)
library(grid)

# Data frame with start and end positions of gene features

df <- data.frame(
  start = c(11,64,115,133,153,155,183,185,189,192,219,225,263,276,278,324,329,489,509,514,545),
  end   = c(64,115,263,138,153,156,184,186,190,194,220,228,276,329,281,329,489,551,514,542,551),
  label = c("Structural domain and glycosylation site", 'Esterase domain', 'Receptor binding domain', 'Sialic acid binding site', 'Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Sialic acid binding site','Esterase domain', 'Structural domain and glycosylation site', 'Glycosylation site only found in H1N1pdm09','Cleavage site', 'Structural domain and glycosylation site', 'Structural domain and glycosylation site','Sialic acid binding site','Transmembrane','Sialic acid binding site')
)

# Filtered data
df_filtered <- df[!df$start %in% c(155, 185, 278), ]

# Colorblind-friendly palette
cb_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)
labels_unique <- unique(df_filtered$label)
fill_colors <- setNames(cb_palette[1:length(labels_unique)], labels_unique)

# Define feature ranges for arrows
feature_ranges <- data.frame(
  start = c(11, 64, 115, 263, 276, 489),
  end   = c(64, 115, 263, 276, 489, 551),
  label = c("Fusion peptide", "Esterase domain", "Receptor binding domain", 
            "Esterase domain", "Fusion peptide", "Structural domain")
)

# Define HA1/HA2 ranges
ha_ranges <- data.frame(
  start = c(11, 329),
  end   = c(329, 551),
  label = c("HA1", "HA2")
)

# Plot
plot <- ggplot(df_filtered) +
  # Protein bars
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=label), colour="black") +
  geom_text(aes(x=start, y=1.5, label=start), size=2, angle=90, hjust=0.9) +
  scale_fill_manual(values = fill_colors) +
  
  # Feature barred arrows (double-headed)
  geom_segment(data=feature_ranges, aes(x=start, xend=end, y=-0.45, yend=-0.45),
               arrow=arrow(length=unit(0.2, "cm"), ends="both", type="closed"),
               color="black", size=1) +
  geom_text(data=feature_ranges, aes(x=(start+end)/2, y=-0.9, label=label), 
            size=3, vjust=1, hjust=0.5) +
  
  # HA1/HA2 barred arrows
  geom_segment(data=ha_ranges, aes(x=start, xend=end, y=-1.5, yend=-1.5),
               arrow=arrow(length=unit(0.2, "cm"), ends="both", type="closed"),
               color="black", size=1) +
  geom_text(data=ha_ranges, aes(x=(start+end)/2, y=-1.95, label=label), 
            size=3, vjust=1, hjust=0.5) +
  
  # Theme adjustments
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = margin(1,0.5,1,0.5, "cm")) +
  ggtitle("HA0 Protein Map of H1N1 consensus sequence") +
  theme(plot.title = element_text(hjust = 0.065, vjust = 2, size=16, face="bold"))

# Decrease y coordinate of legend

plot <- plot + guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))


# Save plot
ggsave("HA0_protein_map_H1N1_consensus.jpg", width=15, height=3, dpi=300)

