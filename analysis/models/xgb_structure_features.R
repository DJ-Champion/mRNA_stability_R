# =============================================================================
# Shared feature-block and eligible-set definition for the structure comparison
# =============================================================================
# The one place the baseline block, the structure block and the eligible row
# set are defined. Both xgb_structure_comparison.R (the main held-out
# comparison) and xgb_structure_gini_subset.R (the icSHAPE secondary run)
# source this file, so the two runs cannot drift apart in what "baseline" or
# "structure" means.
#
# Validation checklist items 2, 3, 4 and 6 of the brief are properties of THIS
# file: Model B's predictors are constructed as BASELINE_COLS + STRUCTURE_COLS
# by definition rather than by two parallel edits, and the eligible set is
# computed once from the union so both models see identical rows.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


# --- Response ----------------------------------------------------------------
# `halflife` is PC1 of the Agarwal & Kelley (2022) consensus half-life measure,
# not a duration in hours: it is already a signed, roughly symmetric score
# (range -17.2 to +18.2, sd 4.81 on the human v9 cache). No project code
# transforms it, and a log is undefined on a variable that takes negative
# values, so it is modelled RAW. RMSE and MAE below are therefore in PC1 units.
TARGET_COL <- "halflife"


# --- Feature-block construction ----------------------------------------------

#' Build the baseline (non-structure) column list for a dataset
#'
#' Section 4 of the brief, resolved against this project's FEATURE_PATTERNS.
#' Deliberately EXCLUDED, and why:
#'
#'   nuc_ratios    (frac_*, 28 cols)  the brief's rule: GC%, AT-skew and GC-skew
#'   compositional (purine_/amino_, 14)  carry the same information without the
#'                                   simplex redundancy of four fractions that
#'                                   sum to 1 per region.
#'   translation_efficiency          measured phenotype, excluded by the brief;
#'                                   also 25.8% missing.
#'   every `structure` supergroup member  that is the experimental variable.
#'
#' Kept beyond the brief's explicit list: `exon_length_last_mrna`, the only
#' surviving member of the `exons` group. It is non-structure transcript
#' architecture (last-exon length is the classical NMD 50-nt-rule covariate),
#' so it belongs with the junction-distance features the brief does name.
#'
#' @param df A dataset from build_dataset() after drop_excluded().
#' @return Character vector of column names present in `df`.
#' @export
baseline_columns <- function(df) {
  cols <- c(
    fg_columns(df, "lengths"),       # 4   regional sequence length
    fg_columns(df, "gc"),            # 7   regional GC content
    fg_columns(df, "skews"),         # 14  AT-skew and GC-skew
    fg_columns(df, "codon_freqs"),   # 65  coding composition
    fg_columns(df, "aa_freqs"),      # 20
    intersect("cai", names(df)),     # 1   CAI (standalone; TE deliberately not)
    fg_columns(df, "stopfree"),      # 4   stop-free length
    fg_columns(df, "uorfs"),         # 1   uorf_present_mrna (numeric 0/1)
    fg_columns(df, "exon_density"),  # 4   CDS-exon density
    fg_columns(df, "eej_dist"),      # 2   junction distance
    fg_columns(df, "nmd"),           # 5   fragile-codon / alternative-stop
    fg_columns(df, "exons")          # 1   exon_length_last_mrna
  )
  unique(cols)
}


#' Build the structure column list for a dataset
#'
#' Length-normalised folding metrics only. Raw MFE scores (`rnafold_score_*`,
#' `rnalfold_score_*`) and per-nucleotide MFE (`rnafold_per_nt_*`) are excluded
#' on DJ's instruction: raw MFE is close to deterministic in length and GC,
#' both of which are already in the baseline, so an apparent structure win
#' could be length re-entering the model under a structure label.
#'
#' icSHAPE Gini is NOT here — see `gini_columns()` and the secondary run. It is
#' 80-91% missing, and a common complete analysis set containing it costs 93%
#' of the corpus.
#'
#' @param df A dataset from build_dataset() after drop_excluded().
#' @return Character vector of column names present in `df`.
#' @export
structure_columns <- function(df) {
  unique(c(
    fg_columns(df, "rnafold_zscores"),   # 8  MFE z-score, per region
    fg_columns(df, "rnalfold_zscores"),  # 7  local MFE z-score
    fg_columns(df, "mfe_deltas")         # 7  observed - expected MFE
  ))
}


#' The icSHAPE structural-Gini block (secondary analysis only)
#' @export
gini_columns <- function(df) fg_columns(df, "probing")


# --- Eligible set ------------------------------------------------------------

#' Load the human dataset and cut it to the common eligible analysis set
#'
#' Complete-case on the UNION of baseline and structure, computed once. Two
#' consequences the brief requires and this guarantees:
#'   - Model A and Model B are handed the same rows and the same gene ids;
#'   - Model B cannot lose rows to the extra missingness in its own block.
#'
#' Complete-case rather than XGBoost's native NA handling, deliberately. The
#' structure columns' missingness is INFORMATIVE — a missing 5'UTR MFE means a
#' 5'UTR too short to fold, not a failed computation — so letting XGBoost learn
#' a default direction for those splits would let Model B profit from an
#' annotation artefact and report it as a structure effect. Cost: 1,800 of
#' 13,601 genes (13%).
#'
#' Zero-variance baseline columns are removed here, on train+val only, so the
#' predictor sets are fixed before any model sees them and the two models
#' cannot end up with different baseline blocks via a recipe filter.
#'
#' @param species Character, passed to build_dataset().
#' @return list(data, baseline, structure, gini, dropped_zv)
#' @export
eligible_dataset <- function(species = "human") {

  df <- build_dataset(species) |>
    drop_excluded(verbose = FALSE) |>
    attach_splits()

  base_cols  <- baseline_columns(df)
  str_cols   <- structure_columns(df)
  gini_cols  <- gini_columns(df)

  stopifnot(length(intersect(base_cols, str_cols)) == 0,
            length(intersect(base_cols, gini_cols)) == 0,
            !TARGET_COL %in% c(base_cols, str_cols, gini_cols),
            length(intersect(META_COLS, c(base_cols, str_cols))) == 0)

  keep <- df |>
    filter(!is.na(.data[[TARGET_COL]]), !is.na(split))
  ok   <- stats::complete.cases(keep[, c(base_cols, str_cols), drop = FALSE])
  keep <- keep[ok, , drop = FALSE]

  # Constant on the data the model may learn from. Checked on train+val rather
  # than on everything, because inspecting test to decide the predictor set is
  # a (mild) use of held-out data.
  learnable <- keep[keep$split != "test", , drop = FALSE]
  zv <- names(which(vapply(learnable[base_cols], function(x)
    length(unique(x)) < 2L, logical(1))))
  base_cols <- setdiff(base_cols, zv)

  list(
    data       = keep,
    baseline   = base_cols,
    structure  = str_cols,
    gini       = gini_cols,
    dropped_zv = zv
  )
}


#' Print the feature lists and sample sizes (checklist item 11)
#' @export
report_feature_sets <- function(el) {
  d <- el$data
  cat("\n=== Eligible analysis set ===\n")
  cat(sprintf("Response          : %s (Agarwal & Kelley consensus PC1, untransformed)\n",
              TARGET_COL))
  cat(sprintf("Genes             : %d (1 row per gene, %d distinct gene_id)\n",
              nrow(d), dplyr::n_distinct(d$gene_id)))
  cat(sprintf("Split (blocked on family_id_%s):\n", BLOCK_LEVEL))
  print(table(d$split))
  cat(sprintf("\nBaseline features : %d\n", length(el$baseline)))
  cat(sprintf("Structure features: %d\n", length(el$structure)))
  cat(sprintf("Model A predictors: %d   Model B predictors: %d\n",
              length(el$baseline), length(el$baseline) + length(el$structure)))
  if (length(el$dropped_zv)) {
    cat("Dropped (zero variance on train+val): ",
        paste(el$dropped_zv, collapse = ", "), "\n")
  }
  cat("\n--- Structure block (the experimental difference) ---\n")
  print(el$structure)
  cat("\n--- Baseline block by family ---\n")
  for (g in names(FEATURE_PATTERNS)) {
    cc <- intersect(fg_columns(d, g), el$baseline)
    if (length(cc)) cat(sprintf("  %-16s %3d\n", g, length(cc)))
  }
  if ("cai" %in% el$baseline) cat("  (standalone: cai; translation_efficiency excluded)\n")
  invisible(el)
}
