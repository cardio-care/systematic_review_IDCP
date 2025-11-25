library(ggplotify)
library(patchwork)

trait_dend <- png::readPNG("trait_dendrogram.png")
combined   <- png::readPNG("combined_trait_jaccard_and_gene_trait_heatmap_right_xaxis_clustered.pdf")
gene_dend  <- png::readPNG("gene_dendrogram.png")

# Convert to grobs
g_trait <- rasterGrob(trait_dend, interpolate = TRUE)
g_comb  <- rasterGrob(combined, interpolate = TRUE)
g_gene  <- rasterGrob(gene_dend, interpolate = TRUE)

# Arrange with gridExtra
grid.arrange(
  gene = g_gene,
  bottom = arrangeGrob(g_trait, g_comb, ncol = 2, widths = c(1, 4)),
  nrow = 2,
  heights = c(1, 4)
)
