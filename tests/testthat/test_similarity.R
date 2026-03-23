PATH_TO_DATA <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")

TOLERANCE <- 1.0e-4

setup_similarity_data <- function() {
  counts <- as.sparse(
    as.matrix(read.table(PATH_TO_DATA, sep = "\t", row.names = 1))
  )
  test.data <- CreateSeuratObject(counts)
  test.data <- NormalizeData(test.data, verbose = FALSE)
  test.data <- FindVariableFeatures(test.data, verbose = FALSE)
  test.data <- ScaleData(test.data, verbose = FALSE)
  test.data <- RunPCA(test.data, npcs = 20, verbose = FALSE)
  test.data <- FindNeighbors(test.data, dims = 1:20, verbose = FALSE)
  test.data <- FindClusters(test.data, resolution = 1, verbose = FALSE)
  return(test.data)
}

context("FindSimilarity")

test_that("FindSimilarity computes pairwise distance", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  result <- FindSimilarity(
    object = test.data,
    ident.1 = idents.all[1],
    ident.2 = idents.all[2],
    reduction = "pca",
    dims = 1:10,
    method = "euclidean",
    verbose = FALSE
  )
  expect_s3_class(object = result, class = "data.frame")
  expect_true(object = "distance" %in% colnames(x = result))
  expect_true(object = "ident.1" %in% colnames(x = result))
  expect_true(object = "ident.2" %in% colnames(x = result))
  expect_true(object = "method" %in% colnames(x = result))
  expect_true(object = "n.cells.1" %in% colnames(x = result))
  expect_true(object = "n.cells.2" %in% colnames(x = result))
  expect_equal(object = nrow(x = result), expected = 1)
  expect_true(object = result$distance > 0)
  expect_equal(object = result$method, expected = "euclidean")
})

test_that("FindSimilarity with ident.2=NULL returns all comparisons", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  result <- FindSimilarity(
    object = test.data,
    ident.1 = idents.all[1],
    ident.2 = NULL,
    reduction = "pca",
    dims = 1:10,
    min.cells = 3,
    verbose = FALSE
  )
  expect_s3_class(object = result, class = "data.frame")
  # Should have one row per other identity (minus any below min.cells)
  expect_true(object = nrow(x = result) >= 1)
  # Results sorted by distance ascending
  expect_true(object = all(diff(x = result$distance) >= 0))
})

test_that("FindAllSimilarity returns distance matrix", {
  test.data <- setup_similarity_data()
  idents.all <- sort(x = levels(x = test.data))
  n.idents <- length(x = idents.all)
  result <- FindAllSimilarity(
    object = test.data,
    reduction = "pca",
    dims = 1:10,
    method = "euclidean",
    min.cells = 3,
    return.matrix = TRUE,
    verbose = FALSE
  )
  expect_true(object = is.matrix(x = result))
  expect_equal(object = nrow(x = result), expected = n.idents)
  expect_equal(object = ncol(x = result), expected = n.idents)
  # Matrix should be symmetric
  expect_equal(object = result, expected = t(x = result), tolerance = TOLERANCE)
  # Diagonal should be zero
  expect_true(object = all(diag(x = result) == 0))
  # Off-diagonal should be positive
  off.diag <- result[row(x = result) != col(x = result)]
  expect_true(object = all(off.diag > 0))
})

test_that("FindAllSimilarity returns long-format data frame", {
  test.data <- setup_similarity_data()
  result <- FindAllSimilarity(
    object = test.data,
    reduction = "pca",
    dims = 1:10,
    min.cells = 3,
    return.matrix = FALSE,
    verbose = FALSE
  )
  expect_s3_class(object = result, class = "data.frame")
  expect_true(object = "ident.1" %in% colnames(x = result))
  expect_true(object = "ident.2" %in% colnames(x = result))
  expect_true(object = "distance" %in% colnames(x = result))
  expect_true(object = nrow(x = result) >= 1)
})

test_that("Different methods produce different results", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  euc <- FindSimilarity(
    test.data, ident.1 = idents.all[1], ident.2 = idents.all[2],
    dims = 1:10, method = "euclidean", verbose = FALSE
  )
  cos <- FindSimilarity(
    test.data, ident.1 = idents.all[1], ident.2 = idents.all[2],
    dims = 1:10, method = "cosine", verbose = FALSE
  )
  cor.result <- FindSimilarity(
    test.data, ident.1 = idents.all[1], ident.2 = idents.all[2],
    dims = 1:10, method = "correlation", verbose = FALSE
  )
  # Methods should produce different distances
  expect_false(object = euc$distance == cos$distance)
  expect_false(object = euc$distance == cor.result$distance)
  # All distances should be positive
  expect_true(object = euc$distance > 0)
  expect_true(object = cos$distance > 0)
  expect_true(object = cor.result$distance > 0)
})

test_that("Permutation test returns p-value", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  result <- FindSimilarity(
    object = test.data,
    ident.1 = idents.all[1],
    ident.2 = idents.all[2],
    dims = 1:10,
    n.perm = 50,
    seed = 42,
    verbose = FALSE
  )
  expect_true(object = "p.value" %in% colnames(x = result))
  expect_true(object = result$p.value >= 0 && result$p.value <= 1)
  # Reproducibility: same seed should give same p-value
  result2 <- FindSimilarity(
    object = test.data,
    ident.1 = idents.all[1],
    ident.2 = idents.all[2],
    dims = 1:10,
    n.perm = 50,
    seed = 42,
    verbose = FALSE
  )
  expect_equal(object = result$p.value, expected = result2$p.value)
})

test_that("min.cells validation works", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  expect_error(
    object = FindSimilarity(
      test.data,
      ident.1 = idents.all[1],
      ident.2 = idents.all[2],
      dims = 1:10,
      min.cells = 10000,
      verbose = FALSE
    )
  )
})

test_that("Invalid method raises error", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  expect_error(
    object = FindSimilarity(
      test.data,
      ident.1 = idents.all[1],
      ident.2 = idents.all[2],
      dims = 1:10,
      method = "nonexistent_method",
      verbose = FALSE
    ),
    regexp = "Unknown similarity method"
  )
})

test_that("Invalid identity raises error", {
  test.data <- setup_similarity_data()
  expect_error(
    object = FindSimilarity(
      test.data,
      ident.1 = "nonexistent_cluster",
      dims = 1:10,
      verbose = FALSE
    ),
    regexp = "not found"
  )
})

test_that("Missing ident.1 raises error", {
  test.data <- setup_similarity_data()
  expect_error(
    object = FindSimilarity(
      test.data,
      dims = 1:10,
      verbose = FALSE
    ),
    regexp = "ident.1"
  )
})

test_that("Custom distance function works", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  # Custom function: Manhattan distance between centroids
  manhattan.dist <- function(emb.1, emb.2) {
    sum(abs(colMeans(x = emb.1) - colMeans(x = emb.2)))
  }
  result <- FindSimilarity(
    object = test.data,
    ident.1 = idents.all[1],
    ident.2 = idents.all[2],
    dims = 1:10,
    method = manhattan.dist,
    verbose = FALSE
  )
  expect_equal(object = nrow(x = result), expected = 1)
  expect_true(object = result$distance > 0)
  expect_equal(object = result$method, expected = "custom")
})

test_that("Euclidean distances match manual centroid computation", {
  test.data <- setup_similarity_data()
  idents.all <- sort(x = levels(x = test.data))
  dims <- 1:10
  # Compute via FindAllSimilarity
  dist.matrix <- FindAllSimilarity(
    object = test.data,
    reduction = "pca",
    dims = dims,
    method = "euclidean",
    min.cells = 3,
    return.matrix = TRUE,
    verbose = FALSE
  )
  # Compute manually (same approach as BuildClusterTree)
  embeddings <- Embeddings(object = test.data[["pca"]])[, dims]
  centroids <- lapply(
    X = idents.all,
    FUN = function(x) {
      cells <- WhichCells(object = test.data, idents = x)
      colMeans(x = embeddings[cells, , drop = FALSE])
    }
  )
  centroid.matrix <- do.call(what = 'rbind', args = centroids)
  manual.dist <- as.matrix(x = dist(x = centroid.matrix))
  rownames(x = manual.dist) <- idents.all
  colnames(x = manual.dist) <- idents.all
  # Should match
  expect_equal(object = dist.matrix, expected = manual.dist, tolerance = TOLERANCE)
})

test_that("weight.by.var changes results", {
  test.data <- setup_similarity_data()
  idents.all <- levels(x = test.data)
  unweighted <- FindSimilarity(
    test.data, ident.1 = idents.all[1], ident.2 = idents.all[2],
    dims = 1:10, weight.by.var = FALSE, verbose = FALSE
  )
  weighted <- FindSimilarity(
    test.data, ident.1 = idents.all[1], ident.2 = idents.all[2],
    dims = 1:10, weight.by.var = TRUE, verbose = FALSE
  )
  expect_false(object = unweighted$distance == weighted$distance)
})
