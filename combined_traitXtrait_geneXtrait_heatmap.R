############################################
#####        REQUIRED LIBRARIES        #####
############################################

library(tidyverse)
library(forcats)
library(patchwork) # side-by-side alignment
library(ggdendro)



############################################
#####             STEP 1               #####
############################################

load("gene_lists.RData")
traits <- names(gene_lists)
n <- length(traits)

jaccard_matrix <- matrix(0, nrow=n, ncol=n, dimnames=list(traits, traits))

for (i in 1:n) {
  for (j in i:n) {
    inter <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
    uni <- length(union(gene_lists[[i]], gene_lists[[j]]))
    jaccard <- ifelse(uni == 0, 0, inter/uni)

    jaccard_matrix[i,j] <- jaccard
    jaccard_matrix[j,i] <- jaccard
  }
}

jaccard_dist <- as.dist(1 - jaccard_matrix)
hc <- hclust(jaccard_dist, method = "average")

trait_order <- hc$labels[hc$order]

ordered_matrix <- jaccard_matrix[trait_order, trait_order]
trans_matrix  <- sqrt(ordered_matrix)


############################################
#####             STEP 2               #####
############################################

trait_data <- read.csv("gene_lists.csv", stringsAsFactors = FALSE)

long_df <- trait_data %>%
  separate_rows(Genes, sep = ",") %>%
  rename(Gene = Genes)

gene_trait_matrix <- long_df %>%
  mutate(Present = 1) %>%
  group_by(Trait, Gene) %>%
  summarise(Present = max(Present), .groups = "drop") %>%
  pivot_wider(names_from = Trait, values_from = Present, values_fill = list(Present = 0))

gene_trait_matrix <- gene_trait_matrix %>%
  mutate(TraitCount = rowSums(select(., -Gene)))

df_long <- gene_trait_matrix %>%
  pivot_longer(cols = -c(Gene, TraitCount),
               names_to = "Trait",
               values_to = "Significant") %>%
  group_by(Gene) %>%
  mutate(TraitCount = first(TraitCount)) %>%
  ungroup() %>%
  mutate(Significant = factor(Significant))

df_pleio <- df_long %>% filter(TraitCount > 1)

# Get only genes present in df_pleio
pleio_genes <- unique(df_pleio$Gene)

# Subset the gene-trait matrix
gene_matrix_subset <- gene_trait_matrix %>%
  filter(Gene %in% pleio_genes) %>%
  select(-TraitCount)   # remove TraitCount before clustering

# Convert to matrix with genes as rownames
gene_matrix <- as.matrix(gene_matrix_subset[,-1])  # remove Gene column
rownames(gene_matrix) <- gene_matrix_subset$Gene

# Compute distance (binary for presence/absence)
gene_dist <- dist(gene_matrix, method = "binary")

# Hierarchical clustering
gene_hc <- hclust(gene_dist, method = "average")

# Order of genes for heatmap
gene_order <- gene_hc$labels[gene_hc$order]

## Step 2b: Cluster genes for x-axis ordering
#gene_matrix <- as.matrix(select(gene_trait_matrix, -Gene, -TraitCount))
#rownames(gene_matrix) <- gene_trait_matrix$Gene
#
## Compute distance and hierarchical clustering
#gene_dist <- dist(gene_matrix, method = "binary")   # or method = "euclidean" if you prefer
#gene_hc <- hclust(gene_dist, method = "average")
#gene_order <- gene_hc$labels[gene_hc$order]

############################################
#####             STEP 3               #####
############################################

left_df <- as.data.frame(trans_matrix)
left_df$Trait1 <- rownames(left_df)

left_df <- left_df %>%
  pivot_longer(cols = -Trait1,
               names_to = "Trait2",
               values_to = "Value")

# apply identical trait order
left_df$Trait1 <- factor(left_df$Trait1, levels = trait_order)
left_df$Trait2 <- factor(left_df$Trait2, levels = trait_order)

# Create dendrogram plot for traits (will go above p_left)
dend_data_trait <- dendro_data(as.dendrogram(hc), type = "rectangle")

p_trait_dend <- ggplot(segment(dend_data_trait)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 0, 5))

p_left <- ggplot(left_df, aes(x = Trait2, y = Trait1, fill = Value)) +
  geom_tile(color = "black", linewidth = 0.2) +        # <-- border
  geom_text(aes(label = sprintf("%.2f", Value)),       # <-- show numbers
            size = 2.2, color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 10, margin = margin(r = 10)), # margin for y-axis title
    axis.title.x = element_text(size = 10, margin = margin(t = 10)), # margin for x-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 10, 5)  # Add bottom margin
  ) +
  labs(x = "Cardiovascular trait categories",
       y = "Cardiovascular trait categories")

# Create dendrogram plot for traits (will go above p_left)
dend_data_trait <- dendro_data(as.dendrogram(hc), type = "rectangle")


############################################
#####             STEP 4               #####
############################################

# Apply trait clustering from left plot
df_pleio$Trait <- factor(df_pleio$Trait, levels = trait_order)

# Apply gene clustering order
df_pleio$Gene <- factor(df_pleio$Gene, levels = gene_order)

# Create dendrogram plot for genes (will go above p_right)
dend_data_gene <- dendro_data(as.dendrogram(gene_hc), type = "rectangle")

p_gene_dend <- ggplot(segment(dend_data_gene)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 0, 5))

p_right <- ggplot(df_pleio, aes(x = Gene, y = Trait, fill = Significant)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = "white", "1" = "#2C7BB6"),
                    guide = "none") +
  scale_y_discrete(limits = trait_order) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10, margin = margin(t = 10)), # margin for x-axis title
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 10, 5)  # Add bottom margin
  ) +
  labs(x = "Genes identified from association analyses")

# Create dendrogram plot for genes (will go above p_right)
dend_data_gene <- dendro_data(as.dendrogram(gene_hc), type = "rectangle")

############################################
#####             STEP 5               #####
############################################

combined_plot <- p_left | p_right

# First, combine each dendrogram with its heatmap
left_panel <- p_trait_dend / p_left + plot_layout(heights = c(1, 4))
right_panel <- p_gene_dend / p_right + plot_layout(heights = c(1, 4))

# Now combine side by side with equal heights
final_plot <- left_panel | right_panel + 
  plot_layout(widths = c(1, 2))

# Display
print(final_plot)

# Optional: Save with good dimensions
ggsave("heatmap_with_dendrograms.pdf", final_plot, 
       width = 14, height = 8, units = "in")
