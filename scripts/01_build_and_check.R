#!/usr/bin/env Rscript
# Script 1: Rebuild docs, install, and check the package
# Run from the seurat root directory:
#   Rscript scripts/01_build_and_check.R

cat("=== Step 1: Regenerate documentation and NAMESPACE ===\n")
tryCatch({
  devtools::document()
  cat("OK: devtools::document() completed\n\n")
}, error = function(e) {
  cat("ERROR in devtools::document():", conditionMessage(e), "\n\n")
})

cat("=== Step 2: Verify NAMESPACE exports ===\n")
ns <- readLines("NAMESPACE")
checks <- c("export(FindSimilarity)", "export(FindAllSimilarity)",
             "S3method(FindSimilarity,Seurat)")
for (chk in checks) {
  found <- any(grepl(chk, ns, fixed = TRUE))
  cat(sprintf("  %s: %s\n", chk, if (found) "FOUND" else "MISSING"))
}
cat("\n")

cat("=== Step 3: Load package and check functions exist ===\n")
tryCatch({
  devtools::load_all(".")
  stopifnot(is.function(FindSimilarity))
  stopifnot(is.function(FindAllSimilarity))
  cat("OK: FindSimilarity and FindAllSimilarity are loadable\n\n")
}, error = function(e) {
  cat("ERROR loading package:", conditionMessage(e), "\n\n")
})

cat("=== Step 4: Quick R CMD check (docs only) ===\n")
cat("  Running: devtools::check(args = '--no-tests --no-vignettes --no-examples')\n")
cat("  (This checks for documentation/NAMESPACE issues without running tests)\n\n")
tryCatch({
  result <- devtools::check(args = c("--no-tests", "--no-vignettes", "--no-examples"),
                            quiet = TRUE)
  cat(sprintf("  Errors: %d, Warnings: %d, Notes: %d\n",
              length(result$errors), length(result$warnings), length(result$notes)))
  if (length(result$errors) > 0) {
    cat("\n  ERRORS:\n")
    for (e in result$errors) cat("  ", e, "\n")
  }
  if (length(result$warnings) > 0) {
    cat("\n  WARNINGS:\n")
    for (w in result$warnings) cat("  ", w, "\n")
  }
}, error = function(e) {
  cat("ERROR in check:", conditionMessage(e), "\n")
  cat("(This is OK if dependencies are missing -- focus on script 02 and 03)\n\n")
})

cat("\n=== Done ===\n")
