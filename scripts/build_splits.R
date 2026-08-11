# =============================================================================
# Build the blocked train / val / test split artefact
# =============================================================================
# Run this ONCE, then read the result everywhere:
#
#   Rscript scripts/build_splits.R                 # blocks on BLOCK_LEVEL
#   Rscript scripts/build_splits.R --level loose   # sensitivity check
#
# Downstream, never re-derive the assignment — read it:
#
#   df <- attach_splits(build_dataset("human"))
#
# Re-deriving it inside a modelling script means a rebuilt family.tsv, a
# different R version, or a changed row order can silently move genes between
# train and test. See FAMILY_CLUSTERING.md §2.2a.
# =============================================================================

source("R/load_all.R")

args  <- commandArgs(trailingOnly = TRUE)
i     <- match("--level", args)
level <- if (!is.na(i) && length(args) > i) args[i + 1L] else BLOCK_LEVEL

cat("Blocking level: ", level,
    if (level == BLOCK_LEVEL) "  (BLOCK_LEVEL)" else "  (override)", "\n",
    sep = "")

splits <- build_splits(level = level)

cat("\nWrote:\n  ", splits_path(level, "rds"),
    "\n  ", splits_path(level, "tsv"), "\n", sep = "")

prov <- attr(splits, "family_provenance")
if (!is.null(prov)) {
  cat("\nBuilt from ", prov$path, "\n  md5:   ", prov$md5,
      "\n  mtime: ", format(prov$mtime), "\n  genes: ", prov$n_rows, "\n",
      sep = "")
}
