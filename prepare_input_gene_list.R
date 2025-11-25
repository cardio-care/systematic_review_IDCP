# Read file
df <- read.delim("gene_lists.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Split comma-separated genes into vectors
gene_lists <- setNames(
  lapply(df$Genes, function(x) unlist(strsplit(x, ","))),
  df$Trait
)

save(gene_lists, file = "gene_lists.RData")
