#!/usr/bin/env Rscript
# Script 3: Interactive demo of FindSimilarity
# Run from the seurat root directory:
#   Rscript scripts/03_demo.R

cat("=== Loading Seurat ===\n")
devtools::load_all(".")

# ─────────────────────────────────────────────────────────
# Setup: load PBMC test data and run standard pipeline
# ─────────────────────────────────────────────────────────
cat("\n=== Setting up PBMC test data ===\n")
pbmc.file <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")
counts <- as.sparse(as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)))
pbmc <- CreateSeuratObject(counts)
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 20, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 1, verbose = FALSE)

n.clusters <- length(levels(pbmc))
cat(sprintf("  Created %d clusters from %d cells\n", n.clusters, ncol(pbmc)))
cat(sprintf("  Cluster sizes: %s\n",
            paste(table(Idents(pbmc)), collapse = ", ")))

# ─────────────────────────────────────────────────────────
# Demo 1: FindSimilarity — single pairwise comparison
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 1: FindSimilarity (single pair) ===\n")
ids <- levels(pbmc)
result <- FindSimilarity(
  object = pbmc,
  ident.1 = ids[1],
  ident.2 = ids[2],
  reduction = "pca",
  dims = 1:20,
  method = "euclidean",
  verbose = FALSE
)
cat("Euclidean distance between cluster", ids[1], "and cluster", ids[2], ":\n")
print(result)

# ─────────────────────────────────────────────────────────
# Demo 2: FindSimilarity — one vs all others
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 2: FindSimilarity (one vs all) ===\n")
result.all <- FindSimilarity(
  object = pbmc,
  ident.1 = ids[1],
  ident.2 = NULL,
  reduction = "pca",
  dims = 1:20,
  method = "euclidean",
  min.cells = 3,
  verbose = FALSE
)
cat("Distances from cluster", ids[1], "to all others (sorted by similarity):\n")
print(result.all)

# ─────────────────────────────────────────────────────────
# Demo 3: Compare all three distance methods
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 3: Method comparison ===\n")
methods <- c("euclidean", "cosine", "correlation")
for (m in methods) {
  res <- FindSimilarity(
    pbmc, ident.1 = ids[1], ident.2 = ids[2],
    dims = 1:20, method = m, verbose = FALSE
  )
  cat(sprintf("  %-12s distance: %.6f\n", m, res$distance))
}

# ─────────────────────────────────────────────────────────
# Demo 4: FindAllSimilarity — full distance matrix
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 4: FindAllSimilarity (distance matrix) ===\n")
dist.matrix <- FindAllSimilarity(
  object = pbmc,
  reduction = "pca",
  dims = 1:20,
  method = "euclidean",
  min.cells = 3,
  return.matrix = TRUE,
  verbose = TRUE
)
cat("\nPairwise Euclidean distance matrix:\n")
print(round(dist.matrix, 3))

# ─────────────────────────────────────────────────────────
# Demo 5: FindAllSimilarity — long format
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 5: FindAllSimilarity (long format, top 10 closest pairs) ===\n")
dist.long <- FindAllSimilarity(
  object = pbmc,
  reduction = "pca",
  dims = 1:20,
  method = "euclidean",
  min.cells = 3,
  return.matrix = FALSE,
  verbose = FALSE
)
dist.long <- dist.long[order(dist.long$distance), ]
cat("Top 10 most similar cluster pairs:\n")
print(head(dist.long, 10))

# ─────────────────────────────────────────────────────────
# Demo 6: Variance weighting comparison
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 6: Effect of variance weighting ===\n")
dist.unweighted <- FindAllSimilarity(
  pbmc, dims = 1:20, method = "euclidean",
  weight.by.var = FALSE, min.cells = 3,
  return.matrix = TRUE, verbose = FALSE
)
dist.weighted <- FindAllSimilarity(
  pbmc, dims = 1:20, method = "euclidean",
  weight.by.var = TRUE, min.cells = 3,
  return.matrix = TRUE, verbose = FALSE
)
# Compare rank orders
unw.vals <- dist.unweighted[upper.tri(dist.unweighted)]
w.vals <- dist.weighted[upper.tri(dist.weighted)]
rank.cor <- cor(rank(unw.vals), rank(w.vals), method = "spearman")
cat(sprintf("  Spearman rank correlation between weighted and unweighted: %.4f\n", rank.cor))
cat("  (1.0 = identical ranking, <1.0 = weighting changes relative distances)\n")

# ─────────────────────────────────────────────────────────
# Demo 7: Permutation test
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 7: Permutation test (100 permutations) ===\n")
perm.result <- FindSimilarity(
  object = pbmc,
  ident.1 = ids[1],
  ident.2 = NULL,
  reduction = "pca",
  dims = 1:20,
  method = "euclidean",
  n.perm = 100,
  min.cells = 3,
  seed = 42,
  verbose = FALSE
)
cat("Distances with permutation p-values:\n")
print(perm.result)

# ─────────────────────────────────────────────────────────
# Demo 8: Custom distance function
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 8: Custom distance function (Manhattan) ===\n")
manhattan <- function(emb.1, emb.2) {
  sum(abs(colMeans(emb.1) - colMeans(emb.2)))
}
custom.result <- FindSimilarity(
  pbmc, ident.1 = ids[1], ident.2 = ids[2],
  dims = 1:20, method = manhattan, verbose = FALSE
)
cat("Manhattan distance:\n")
print(custom.result)

# ─────────────────────────────────────────────────────────
# Demo 9: Compare with BuildClusterTree
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 9: Consistency check vs BuildClusterTree ===\n")
if (requireNamespace("ape", quietly = TRUE)) {
  pbmc.tree <- BuildClusterTree(pbmc, dims = 1:20, verbose = FALSE)
  tree <- Tool(pbmc.tree, slot = "BuildClusterTree")
  cat("BuildClusterTree dendrogram tip labels: ", paste(tree$tip.label, collapse = ", "), "\n")
  cat("FindAllSimilarity distance matrix:\n")
  print(round(dist.matrix, 3))
  cat("\n(The dendrogram is built from the same Euclidean centroid distances.\n")
  cat(" You can verify by running: hclust(as.dist(dist.matrix)) )\n")
} else {
  cat("  Skipped: 'ape' package not installed\n")
}

# ─────────────────────────────────────────────────────────
# Demo 10: Hierarchical clustering from distance matrix
# ─────────────────────────────────────────────────────────
cat("\n=== Demo 10: Hierarchical clustering from FindAllSimilarity ===\n")
hc <- hclust(as.dist(dist.matrix), method = "average")
cat("Dendrogram merge order:\n")
print(hc$merge)
cat("\nCluster labels: ", paste(hc$labels, collapse = ", "), "\n")
cat("Merge heights: ", paste(round(hc$height, 3), collapse = ", "), "\n")

cat("\n=== All demos complete ===\n")
