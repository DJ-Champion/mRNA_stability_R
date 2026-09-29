# =============================================================================
# Compare the migrated correlation scripts to the step-1 baseline
# =============================================================================
# Runs the ranked, bands, dotplot and scatter functions with their new
# include / exclude / top_n arguments, exactly as their runners call them, and
# diffs the (stem, region) sets they plot against scripts/selection_baseline.csv.
#
# Expected differences (SELECTION_PLAN.md):
#   * probing (the four icSHAPE Gini columns) is gone from every "core" job;
#   * dotplot and scatter "core" jobs used to pin two hand-named codons and two
#     amino acids (bundle sequence_select); they now show the true top two, as
#     ranked and bands always did.
# Anything else is a regression.
#
# Usage (project root):  Rscript scripts/diff_selection_baseline.R
# =============================================================================

suppressMessages(source("R/load_all.R"))

load_fns <- function(path) {
  e <- new.env(parent = globalenv())
  suppressMessages(sys.source(path, envir = e, keep.source = FALSE))
  e
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

df_full <- suppressMessages(build_dataset("human"))
df_utr  <- suppressMessages(build_dataset("human", min_utr = MIN_UTR_LENGTH))
base    <- read.csv("scripts/selection_baseline.csv", stringsAsFactors = FALSE)

ranked  <- load_fns("analysis/correlations/feature_correlation_ranked.R")
dotplot <- load_fns("analysis/correlations/feature_correlation_dotplot.R")
scatter <- load_fns("analysis/correlations/feature_response_scatter.R")

top2 <- list(codon_freqs = 2, aa_freqs = 2)
now  <- list()   # key "script | job" -> character vector "feature_id | stem | region"
key3 <- function(g, stem, region) paste(g, stem, region, sep = " | ")
record <- function(script, job, tab) {
  stem <- if ("metric_stem" %in% names(tab)) tab$metric_stem else tab$variable
  now[[paste(script, job, sep = " | ")]] <<-
    unique(key3(tab$group, stem, tab$region))
}

for (j in list(
  list("halflife", "halflife", "core", top2),
  list("translation_efficiency", "translation_efficiency", "core", top2),
  list("halflife", "halflife_aa_full", "aa_freqs", NULL),
  list("halflife", "halflife_codon_full", "codon_freqs", NULL))) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = j[[1]], include = j[[3]], orientation = "horizontal",
    sig_threshold = "auto", top_n = j[[4]]))
  record("feature_correlation_ranked.R", j[[2]], out$table)
}
for (resp in c("halflife", "translation_efficiency")) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = resp, include = "core", orientation = "horizontal",
    top_n = top2))
  record("feature_correlation_bands.R", resp, out$table)
}
for (j in list(
  list("halflife", "halflife", "core", top2),
  list("translation_efficiency", "translation_efficiency", "core", top2),
  list("halflife", "halflife_codon_aa", c("codon_freqs", "aa_freqs"),
       list(codon_freqs = 15, aa_freqs = 20)),
  list("halflife", "halflife_nuc_ratios", c("nuc_ratios", "gc"), NULL))) {
  out <- quiet(dotplot$feature_correlation_dotplot(
    df_full, response = j[[1]], sig_threshold = 0.02, include = j[[3]],
    absolute = TRUE, top_n = j[[4]]))
  record("feature_correlation_dotplot.R", j[[2]], out$table)
}
for (j in list(
  list("included", list(top_n = top2, noise_filter = 0, label_quantile = 0.8)),
  list("top_3",    list(max_features = 3, label_quantile = 0.3)),
  list("top_1",    list(max_features = 1, label_quantile = 0)))) {
  out <- quiet(do.call(scatter$feature_response_scatter, c(list(df_full), j[[2]])))
  # the scatter table has region "none" for region-less columns; the baseline
  # recorded the same table, so compare like with like
  record("feature_response_scatter.R", j[[1]], out$table)
}

# Baseline sets, same key shape.
bp <- base[base$stage == "plotted", ]
was <- split(key3(bp$feature_id, bp$metric_stem, bp$region),
             paste(bp$script, bp$job, sep = " | "))

n_unexpected <- 0L
for (k in names(now)) {
  gone  <- setdiff(was[[k]], now[[k]])
  added <- setdiff(now[[k]], was[[k]])
  probing <- grepl("^probing \\|", gone)
  pinned  <- grepl("^(codon_freqs|aa_freqs) \\|", c(gone, added))
  unexpected <- c(gone[!probing & !grepl("^(codon_freqs|aa_freqs) \\|", gone)],
                  added[!grepl("^(codon_freqs|aa_freqs) \\|", added)])
  cat(sprintf("%-58s was %3d, now %3d | -probing %d | codon/aa changed %d | UNEXPECTED %d\n",
              k, length(was[[k]]), length(now[[k]]), sum(probing),
              sum(pinned), length(unexpected)))
  for (u in unexpected) cat("    ", u, "\n")
  n_unexpected <- n_unexpected + length(unexpected)
}
if (n_unexpected) { cat("\n", n_unexpected, " unexpected difference(s)\n", sep = ""); quit(status = 1) }
cat("\nOnly the expected differences.\n")
