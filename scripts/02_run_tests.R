#!/usr/bin/env Rscript
# Script 2: Run the FindSimilarity test suite
# Run from the seurat root directory:
#   Rscript scripts/02_run_tests.R

cat("=== Loading package ===\n")
devtools::load_all(".")

cat("\n=== Running FindSimilarity tests ===\n\n")
result <- devtools::test(filter = "similarity")

cat("\n=== Test Summary ===\n")
print(result)
