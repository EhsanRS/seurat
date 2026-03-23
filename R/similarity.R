#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param ident.1 Identity class to compare (required)
#' @param ident.2 A second identity class for comparison. If NULL (default),
#' computes distance from \code{ident.1} to all other identities
#' @param reduction Reduction to use for distance computation (default: "pca")
#' @param dims Dimensions of the reduction to use. If NULL, uses all available
#' dimensions
#' @param method Distance method to use. One of "euclidean", "cosine", or
#' "correlation". Can also pass a custom function that takes two matrices
#' (cells x dims) and returns a single numeric distance
#' @param group.by Group cells by a metadata column instead of active identities
#' @param features Compute distance in gene space using these features instead
#' of a reduction. Mutually exclusive with \code{reduction}
#' @param assay Assay to use for gene-space mode (default: active assay)
#' @param slot Slot to use for gene-space mode (default: "data")
#' @param weight.by.var Weight each dimension by the variance explained
#' (requires a reduction with stored standard deviations)
#' @param n.perm Number of permutations for significance testing. Set to 0
#' (default) to skip permutation testing
#' @param min.cells Minimum number of cells required in each group (default: 10)
#' @param seed Random seed for reproducibility of permutation tests
#' @param verbose Print progress messages
#'
#' @importFrom stats cor p.adjust setNames
#' @importFrom pbapply pbsapply
#'
#' @rdname FindSimilarity
#' @method FindSimilarity Seurat
#' @export
#'
FindSimilarity.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  reduction = "pca",
  dims = NULL,
  method = "euclidean",
  group.by = NULL,
  features = NULL,
  assay = NULL,
  slot = "data",
  weight.by.var = FALSE,
  n.perm = 0,
  min.cells = 10,
  seed = 42,
  verbose = TRUE,
  ...
) {
  if (is.null(x = ident.1)) {
    stop("Please provide 'ident.1'", call. = FALSE)
  }
  # Handle group.by
  if (!is.null(x = group.by)) {
    if (!group.by %in% colnames(x = object[[]])) {
      stop(
        "'group.by' column '", group.by, "' not found in object metadata",
        call. = FALSE
      )
    }
    Idents(object = object) <- group.by
  }
  # Validate ident.1 exists
  idents.all <- levels(x = object)
  if (!ident.1 %in% idents.all) {
    stop(
      "Identity '", ident.1, "' not found. Available identities: ",
      paste(idents.all, collapse = ", "),
      call. = FALSE
    )
  }
  # Extract embeddings
  if (!is.null(x = features)) {
    # Gene-space mode
    if (!is.null(x = reduction) && reduction != "pca") {
      warning(
        "Both 'features' and non-default 'reduction' specified. ",
        "Using gene-space mode with 'features'.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- GetAssayData(object = object, assay = assay, slot = slot)
    features <- intersect(x = features, y = rownames(x = data.use))
    if (length(x = features) == 0) {
      stop("None of the requested features found in assay data", call. = FALSE)
    }
    embeddings <- t(x = as.matrix(x = data.use[features, , drop = FALSE]))
  } else {
    # Reduction mode
    if (is.null(x = object[[reduction]])) {
      stop(
        "Reduction '", reduction, "' not found in object. ",
        "Run a dimensional reduction first.",
        call. = FALSE
      )
    }
    embeddings <- Embeddings(object = object[[reduction]])
    if (!is.null(x = dims)) {
      if (max(dims) > ncol(x = embeddings)) {
        stop(
          "Requested dims exceed available dimensions in '", reduction,
          "' (", ncol(x = embeddings), " available)",
          call. = FALSE
        )
      }
      embeddings <- embeddings[, dims, drop = FALSE]
    }
    # Variance weighting
    if (isTRUE(x = weight.by.var)) {
      sdev <- Stdev(object = object[[reduction]])
      if (length(x = sdev) > 0) {
        if (!is.null(x = dims)) {
          sdev <- sdev[dims]
        } else {
          sdev <- sdev[1:ncol(x = embeddings)]
        }
        embeddings <- embeddings %*% diag(x = sdev)
      } else {
        warning(
          "No standard deviations found in reduction '", reduction,
          "'. Skipping variance weighting.",
          call. = FALSE,
          immediate. = TRUE
        )
      }
    }
  }
  # Determine comparison targets
  if (is.null(x = ident.2)) {
    targets <- setdiff(x = idents.all, y = ident.1)
  } else {
    if (!all(ident.2 %in% idents.all)) {
      missing <- setdiff(x = ident.2, y = idents.all)
      stop(
        "Identity(ies) not found: ", paste(missing, collapse = ", "),
        ". Available identities: ", paste(idents.all, collapse = ", "),
        call. = FALSE
      )
    }
    targets <- ident.2
  }
  # Get cells for ident.1
  cells.1 <- WhichCells(object = object, idents = ident.1)
  if (length(x = cells.1) < min.cells) {
    stop(
      "Identity '", ident.1, "' has ", length(x = cells.1),
      " cells, fewer than min.cells (", min.cells, ")",
      call. = FALSE
    )
  }
  # Compute distances for each target
  results <- list()
  for (target in targets) {
    cells.2 <- WhichCells(object = object, idents = target)
    if (length(x = cells.2) < min.cells) {
      if (verbose) {
        message(
          "Skipping identity '", target, "': only ",
          length(x = cells.2), " cells (min.cells = ", min.cells, ")"
        )
      }
      next
    }
    emb.1 <- embeddings[cells.1, , drop = FALSE]
    emb.2 <- embeddings[cells.2, , drop = FALSE]
    dist.val <- ComputeClusterDistance(
      embeddings.1 = emb.1,
      embeddings.2 = emb.2,
      method = method
    )
    result <- data.frame(
      ident.1 = ident.1,
      ident.2 = target,
      distance = dist.val,
      method = if (is.function(x = method)) "custom" else method,
      n.cells.1 = length(x = cells.1),
      n.cells.2 = length(x = cells.2),
      stringsAsFactors = FALSE
    )
    # Permutation testing
    if (n.perm > 0) {
      pval <- PermuteSimilarity(
        embeddings = embeddings,
        cells.1 = cells.1,
        cells.2 = cells.2,
        method = method,
        n.perm = n.perm,
        observed.dist = dist.val,
        seed = seed,
        verbose = verbose
      )
      result$p.value <- pval
    }
    results[[length(x = results) + 1]] <- result
  }
  if (length(x = results) == 0) {
    warning("No valid comparisons found", call. = FALSE, immediate. = TRUE)
    return(data.frame())
  }
  output <- do.call(what = 'rbind', args = results)
  rownames(x = output) <- NULL
  # Multiple testing correction for p-values
  if (n.perm > 0 && nrow(x = output) > 1) {
    output$p.value.adj.bonf <- p.adjust(p = output$p.value, method = "bonferroni")
    output$p.value.adj.bh <- p.adjust(p = output$p.value, method = "BH")
  }
  # Sort by distance ascending (most similar first)
  output <- output[order(output$distance), ]
  rownames(x = output) <- NULL
  return(output)
}

#' Compute all pairwise inter-cluster distances
#'
#' Iterates \code{\link{FindSimilarity}} over all pairs of clusters, producing
#' a full pairwise distance matrix or long-format data frame. Analogous to
#' \code{\link{FindAllMarkers}}.
#'
#' @param object A Seurat object
#' @param reduction Reduction to use (default: "pca")
#' @param dims Dimensions of the reduction to use
#' @param method Distance method (default: "euclidean")
#' @param group.by Group cells by a metadata column instead of active identities
#' @param features Compute distance in gene space using these features
#' @param assay Assay to use for gene-space mode
#' @param slot Slot to use for gene-space mode
#' @param weight.by.var Weight dimensions by variance explained
#' @param min.cells Minimum cells per group
#' @param return.matrix If TRUE (default), return a symmetric distance matrix.
#' If FALSE, return a long-format data frame
#' @param verbose Print progress messages
#' @param ... Additional arguments passed to \code{\link{FindSimilarity}}
#'
#' @return When \code{return.matrix = TRUE}, a symmetric numeric matrix with
#' cluster identities as row and column names and distances as values.
#' When \code{return.matrix = FALSE}, a data.frame with columns
#' \code{ident.1}, \code{ident.2}, \code{distance}, \code{method},
#' \code{n.cells.1}, \code{n.cells.2}
#'
#' @export
#'
#' @concept similarity
#'
FindAllSimilarity <- function(
  object,
  reduction = "pca",
  dims = NULL,
  method = "euclidean",
  group.by = NULL,
  features = NULL,
  assay = NULL,
  slot = "data",
  weight.by.var = FALSE,
  min.cells = 10,
  return.matrix = TRUE,
  verbose = TRUE,
  ...
) {
  # Handle group.by
  if (!is.null(x = group.by)) {
    if (!group.by %in% colnames(x = object[[]])) {
      stop(
        "'group.by' column '", group.by, "' not found in object metadata",
        call. = FALSE
      )
    }
    Idents(object = object) <- group.by
  }
  idents.all <- sort(x = levels(x = object))
  n.idents <- length(x = idents.all)
  if (n.idents < 2) {
    stop("Need at least 2 identity classes to compute distances", call. = FALSE)
  }
  # Generate all unique pairs
  pairs <- combn(x = idents.all, m = 2, simplify = FALSE)
  if (verbose) {
    message(
      "Computing ", length(x = pairs), " pairwise distances across ",
      n.idents, " identities"
    )
  }
  # Compute distance for each pair
  my.lapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
  results <- my.lapply(
    X = pairs,
    FUN = function(pair) {
      tryCatch(
        expr = {
          FindSimilarity(
            object = object,
            ident.1 = pair[1],
            ident.2 = pair[2],
            reduction = reduction,
            dims = dims,
            method = method,
            features = features,
            assay = assay,
            slot = slot,
            weight.by.var = weight.by.var,
            min.cells = min.cells,
            verbose = FALSE,
            ...
          )
        },
        error = function(e) {
          if (verbose) {
            message("Skipping pair (", pair[1], ", ", pair[2], "): ", e$message)
          }
          return(NULL)
        }
      )
    }
  )
  # Remove NULL results (failed pairs)
  results <- Filter(f = Negate(f = is.null), x = results)
  if (length(x = results) == 0) {
    stop("No valid pairwise distances computed", call. = FALSE)
  }
  all.results <- do.call(what = 'rbind', args = results)
  rownames(x = all.results) <- NULL
  if (isTRUE(x = return.matrix)) {
    # Build symmetric distance matrix
    dist.matrix <- matrix(
      data = 0,
      nrow = n.idents,
      ncol = n.idents,
      dimnames = list(idents.all, idents.all)
    )
    for (i in seq_len(length.out = nrow(x = all.results))) {
      id1 <- as.character(x = all.results$ident.1[i])
      id2 <- as.character(x = all.results$ident.2[i])
      dist.matrix[id1, id2] <- all.results$distance[i]
      dist.matrix[id2, id1] <- all.results$distance[i]
    }
    return(dist.matrix)
  } else {
    return(all.results)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Distance method registry
#'
#' Returns a named list of available distance methods and their implementing
#' functions. Each function takes two matrices (cells x dims) and returns a
#' single numeric distance value.
#'
#' @return A named list of distance functions
#'
#' @keywords internal
#' @noRd
#'
SimilarityMethods <- function() {
  list(
    euclidean   = EuclideanCentroidDist,
    cosine      = CosineCentroidDist,
    correlation = CorrelationCentroidDist
  )
}

#' Compute distance between two groups of cells
#'
#' Dispatches to the appropriate distance function based on the method name
#' or custom function.
#'
#' @param embeddings.1 Matrix of embeddings for group 1 (cells x dims)
#' @param embeddings.2 Matrix of embeddings for group 2 (cells x dims)
#' @param method Distance method name or custom function
#'
#' @return A single numeric distance value
#'
#' @keywords internal
#' @noRd
#'
ComputeClusterDistance <- function(embeddings.1, embeddings.2, method = "euclidean") {
  if (is.function(x = method)) {
    return(method(embeddings.1, embeddings.2))
  }
  methods.available <- SimilarityMethods()
  if (!method %in% names(x = methods.available)) {
    stop(
      "Unknown similarity method: '", method, "'. ",
      "Available methods: ", paste(names(x = methods.available), collapse = ", "),
      ". You can also pass a custom function.",
      call. = FALSE
    )
  }
  return(methods.available[[method]](embeddings.1, embeddings.2))
}

#' Permutation test for cluster distance significance
#'
#' Tests whether the observed distance between two groups is greater than
#' expected by chance. Pools all cells, randomly reassigns to two groups of
#' the original sizes, and recomputes the distance.
#'
#' @param embeddings Full embedding matrix (all cells x dims)
#' @param cells.1 Cell names for group 1
#' @param cells.2 Cell names for group 2
#' @param method Distance method
#' @param n.perm Number of permutations
#' @param observed.dist The observed distance to compare against
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return Empirical p-value
#'
#' @keywords internal
#' @noRd
#'
PermuteSimilarity <- function(
  embeddings,
  cells.1,
  cells.2,
  method,
  n.perm,
  observed.dist,
  seed = 42,
  verbose = TRUE
) {
  set.seed(seed = seed)
  all.cells <- c(cells.1, cells.2)
  n1 <- length(x = cells.1)
  null.dists <- vapply(
    X = seq_len(length.out = n.perm),
    FUN = function(i) {
      shuffled <- sample(x = all.cells, size = length(x = all.cells))
      perm.cells.1 <- shuffled[1:n1]
      perm.cells.2 <- shuffled[(n1 + 1):length(x = shuffled)]
      ComputeClusterDistance(
        embeddings.1 = embeddings[perm.cells.1, , drop = FALSE],
        embeddings.2 = embeddings[perm.cells.2, , drop = FALSE],
        method = method
      )
    },
    FUN.VALUE = numeric(length = 1)
  )
  # Empirical p-value with pseudocount
  p.value <- (sum(null.dists >= observed.dist) + 1) / (n.perm + 1)
  return(p.value)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Distance method implementations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Euclidean centroid distance
#'
#' Computes the L2 norm between cluster centroids (column means).
#'
#' @param embeddings.1 Matrix of embeddings for group 1 (cells x dims)
#' @param embeddings.2 Matrix of embeddings for group 2 (cells x dims)
#'
#' @return Euclidean distance between centroids
#'
#' @keywords internal
#' @noRd
#'
EuclideanCentroidDist <- function(embeddings.1, embeddings.2) {
  centroid.1 <- colMeans(x = embeddings.1)
  centroid.2 <- colMeans(x = embeddings.2)
  return(sqrt(x = sum((centroid.1 - centroid.2)^2)))
}

#' Cosine centroid distance
#'
#' Computes 1 minus the cosine similarity between cluster centroids.
#'
#' @param embeddings.1 Matrix of embeddings for group 1 (cells x dims)
#' @param embeddings.2 Matrix of embeddings for group 2 (cells x dims)
#'
#' @return Cosine distance between centroids
#'
#' @keywords internal
#' @noRd
#'
CosineCentroidDist <- function(embeddings.1, embeddings.2) {
  centroid.1 <- colMeans(x = embeddings.1)
  centroid.2 <- colMeans(x = embeddings.2)
  dot.product <- sum(centroid.1 * centroid.2)
  norm.1 <- sqrt(x = sum(centroid.1^2))
  norm.2 <- sqrt(x = sum(centroid.2^2))
  if (norm.1 == 0 || norm.2 == 0) {
    warning("Zero-norm centroid detected; cosine distance undefined", call. = FALSE)
    return(NaN)
  }
  return(1 - dot.product / (norm.1 * norm.2))
}

#' Correlation centroid distance
#'
#' Computes 1 minus the Pearson correlation between cluster centroids.
#'
#' @param embeddings.1 Matrix of embeddings for group 1 (cells x dims)
#' @param embeddings.2 Matrix of embeddings for group 2 (cells x dims)
#'
#' @return Correlation distance between centroids
#'
#' @keywords internal
#' @noRd
#'
CorrelationCentroidDist <- function(embeddings.1, embeddings.2) {
  centroid.1 <- colMeans(x = embeddings.1)
  centroid.2 <- colMeans(x = embeddings.2)
  return(1 - cor(x = centroid.1, y = centroid.2))
}
