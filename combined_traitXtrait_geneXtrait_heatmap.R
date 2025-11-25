############################################
#####        REQUIRED LIBRARIES        #####
############################################

library(tidyverse)
library(forcats)
library(patchwork) # side-by-side alignment




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

p_left <- ggplot(left_df, aes(x = Trait2, y = Trait1, fill = Value)) +
  geom_tile(color = "black", linewidth = 0.2) +        # <-- border
  geom_text(aes(label = sprintf("%.2f", Value)),       # <-- show numbers
            size = 2.2, color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Cardiovascular trait categories")


############################################
#####             STEP 4               #####
############################################

# Apply trait clustering from left plot
df_pleio$Trait <- factor(df_pleio$Trait, levels = trait_order)

# Apply gene clustering order
df_pleio$Gene <- factor(df_pleio$Gene, levels = gene_order)

p_right <- ggplot(df_pleio, aes(x = Gene, y = Trait, fill = Significant)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = "white", "1" = "#2C7BB6"),
                    guide = "none") +
  scale_y_discrete(limits = trait_order) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Genes identified from association analyses")

############################################
#####             STEP 5               #####
############################################

combined_plot <- p_left | p_right

tiff("combined_trait_jaccard_and_gene_trait_heatmap_right_xaxis_clustered.tiff",
     width = 9000, height = 6000, res = 600)
print(combined_plot)
dev.off()

# Save trait dendrogram
tiff("trait_dendrogram.tiff", width = 2000, height = 1500, res = 300)
plot(hc, main = "Trait Clustering Dendrogram",
     xlab = "", sub = "", cex = 0.8, hang = -1)
dev.off()

# Save gene dendrogram
tiff("gene_dendrogram.tiff", width = 10000, height = 6000, res = 600)
plot(gene_hc, main = "Gene Clustering Dendrogram",
     xlab = "", sub = "", cex = 0.6, las = 2, hang = -1)
dev.off()
