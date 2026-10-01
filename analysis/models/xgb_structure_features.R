# =============================================================================
# Shared feature-block and eligible-set definition for the two XGBoost models
# =============================================================================
# The one place the two models are defined:
#
#   Without Structure   non-structure transcript features
#   With Structure      those + every computed secondary-structure feature
#
# Structure is the ONLY difference between them. Rows, gene ids, preprocessing,
# tuning resamples, tuning grid, budget, seeds and evaluation are shared by
# construction — `Structure` is built as BASELINE + the structure block, not as
# a second hand-maintained list that has to be kept in step.
#
# WHAT COUNTS AS STRUCTURE. Every member of the `structure` supergroup
# (config.R) EXCEPT `probing`. That is the computed, sequence-derived folding
# block: RNAfold and RNALfold MFE, their z-scores against shuffled sequence,
# the per-nucleotide normalisation, and MFE delta (observed - expected).
# icSHAPE structural Gini (`probing`) is experimental readout, not computed
# from sequence, and is deliberately outside both models — see PROBING_GROUP.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


# --- Response ----------------------------------------------------------------
# Chosen with the XGB_TARGET environment variable (default "halflife"):
#
#   XGB_TARGET=translation_efficiency Rscript analysis/models/xgb_structure_comparison.R
#
# Everything else — predictors, split, tuning, bootstrap, figures — is shared.
# Each response writes to its own folder (run_dir()), and the fit cache is keyed
# on the target, so responses never overwrite or reuse each other's results.
#
# ROW POLICY per response: every gene with a non-missing target and a split.
# For translation efficiency (missing for ~26% of human genes) that shrinks the
# eligible set; the committed family-blocked split is reused as it stands,
# because dropping genes cannot put a family on both sides of it.
#
# Both responses are modelled RAW. `halflife` is PC1 of the Agarwal & Kelley
# (2022) consensus half-life measure, not a duration in hours: a signed, roughly
# symmetric score (range -17.2 to +18.2, sd 4.81 on the human v10 cache), and a
# log is undefined on negative values. RMSE and MAE are in the response's units (for TE, a unitless residual).
# The response columns are never predictors: halflife is in the "Response /
# evaluation" group, which cannot be selected, and translation_efficiency is
# flagged out of the model in R/feature_table.csv.
RESPONSES <- list(
  halflife = list(
    label = "half-life",
    axis  = "PC1 score",
    note  = paste("Agarwal & Kelley 2022 consensus half-life PC1,",
                  "untransformed (signed score, not hours)")),
  translation_efficiency = list(
    label = "translation efficiency",
    axis  = "Translation Efficiency",
    note  = paste("Ribo-seq translation efficiency (mean_te), a unitless",
                  "compositional-regression residual (Liu et al. 2025;",
                  "Zheng et al. 2025), untransformed"))
)

TARGET_COL <- Sys.getenv("XGB_TARGET", "halflife")
if (!TARGET_COL %in% names(RESPONSES)) {
  stop("XGB_TARGET '", TARGET_COL, "' is not a known response. Known: ",
       paste(names(RESPONSES), collapse = ", "), call. = FALSE)
}

#' Labels for a response (label, axis, note).
#' @export
response_info <- function(target = TARGET_COL) RESPONSES[[target]]


# --- The two models ----------------------------------------------------------

REFERENCE_MODEL <- "Without Structure"
STRUCTURE_MODEL <- "With Structure"

# Order matters: reference first. Every table, figure and factor level in the
# downstream scripts reads this vector rather than restating the order.
MODELS <- c(REFERENCE_MODEL, STRUCTURE_MODEL)

# The `structure` supergroup member that is NOT part of the structure block.
# icSHAPE Gini is 80-91% missing on the human cache, and — more importantly —
# it is MEASURED rather than computed, so a model containing it cannot score a
# transcript nobody has probed. Anything using it is a different kind of model
# and belongs in its own, later, supplementary analysis.
PROBING_GROUP <- "probing"


#' The feature ids that make up the structure block.
#'
#' The model features (R/feature_table.csv, "Included in model") that sit in
#' the Structure supergroup. Nothing is hand-listed: a folding feature flagged
#' into the model joins this block, and one flagged out leaves it. probing and
#' mfe_expected are Structure rows the table keeps out of the model. A
#' function, not a constant, so this file can be sourced in any order.
#' @export
structure_groups <- function() {
  selected_features("model", exclude = setdiff(names(SUPERGROUPS), "structure"))
}


#' Where this analysis's artefacts live: one subfolder per response, so runs
#' for different responses (halflife, translation_efficiency, ...) never
#' overwrite each other.
#' @param what "root", "tables" or "plots".
#' @export
run_dir <- function(what = c("root", "tables", "plots")) {
  what <- match.arg(what)
  root <- file.path(OUTPUT_DIR, "xgb_structure", TARGET_COL)
  if (what == "root") root else file.path(root, what)
}


# --- Feature-block construction ----------------------------------------------

#' Build the baseline (non-structure) column list for a dataset
#'
#' Every model feature outside the Structure supergroup. The table's model flag
#' decides; its Notes column records each reason. For the record, the
#' non-structure features deliberately left out, and why:
#'
#'   aa_freqs      (aa_*, 20 cols)  an EXACT deterministic function of the
#'                                  codon columns that are kept. Verified on
#'                                  the human cache to 5.6e-17:
#'                                    aa_x = sum(codons encoding x)
#'                                           / (1 - stop_fraction - codon_other)
#'                                  So they add no information, and 20 exact
#'                                  substitutes for retained columns splinter
#'                                  gain importance and shrink the share of any
#'                                  column subsample the structure block can
#'                                  occupy. Dropped from the MODEL baseline
#'                                  only — they stay in the cache, because the
#'                                  correlation plots and the sequence_select
#'                                  bundle legitimately use them.
#'
#'   nuc_ratios    (frac_*, 28)     likewise exact. GC content and the two
#'   compositional (purine_/amino_, 14)  skews already in the baseline recover
#'                                  all four nucleotide fractions exactly:
#'                                    g = gc(1+gc_skew)/2   c = gc(1-gc_skew)/2
#'                                    a = (1-gc)(1+at_skew)/2
#'                                    u = (1-gc)(1-at_skew)/2
#'                                  and purine = a+g, amino = a+c. Verified to
#'                                  1e-16 on the same cache.
#'
#'   translation_efficiency         measured phenotype, not a sequence feature;
#'                                  also 25.8% missing.
#'
#'   exons         (exon_length_last_mrna, 1)  a 3'UTR-length proxy. The last
#'                                  exon carries the stop codon plus almost all
#'                                  of the 3'UTR, so on the v10 human cache it
#'                                  is Spearman 0.949 with length_3utr
#'                                  (rho^2 = 0.90) and exceeds the 3'UTR length
#'                                  in 96.3% of transcripts. length_3utr is
#'                                  already in the baseline via `lengths`, so
#'                                  this is a second, noisier copy of a
#'                                  retained column. The 50-nt-rule geometry it
#'                                  was once kept for is carried directly by
#'                                  the retained eej_dist_closest_* columns.
#'
#'   every `structure` supergroup member   that is the experimental variable.
#'
#' @param df A dataset from build_dataset() after drop_excluded().
#' @return Character vector of column names present in `df`.
#' @export
baseline_columns <- function(df) {
  selected_columns(select_features(df, "model", exclude = "structure"))
}


#' Build the structure column list
#'
#' Every computed folding feature: RNAfold and RNALfold MFE and z-scores, the
#' per-nucleotide normalisation, and MFE delta — the model features of the
#' Structure supergroup (structure_groups()). `mfe_expected` is flagged out of
#' the model as scaffolding for mfe_delta_*, so it is not among them.
#'
#' CONFOUNDING, TO REPORT WITH ANY RESULT. This block is NOT length- and
#' GC-neutral. Raw MFE scales almost linearly with sequence length and shifts
#' with GC content, both of which are already in the baseline; MFE delta is
#' observed minus an expected value that is itself a deterministic function of
#' GC and length. The z-scores are the only members normalised against
#' shuffled sequence. So a win for `Structure` is incremental predictive
#' information carried by the folding block AS A WHOLE, and cannot on its own
#' be attributed to secondary structure rather than to length and GC
#' re-entering the model under a structure label. Section 13e of the comparison
#' script quantifies exactly how much of each column the baseline already
#' explains, which is the check that separates the two readings.
#'
#' @param df A dataset from build_dataset() after drop_excluded().
#' @return Character vector of column names present in `df`.
#' @export
structure_columns <- function(df) {
  selected_columns(select_features(
    df, "model", exclude = setdiff(names(SUPERGROUPS), "structure")))
}


#' The icSHAPE structural-Gini block — excluded from both models
#'
#' Retained so the validation checklist can ASSERT its absence rather than
#' assume it. See PROBING_GROUP.
#' @export
probing_columns <- function(df) fg_columns(df, PROBING_GROUP)


#' The predictor list for one model
#' @param el Result of eligible_dataset().
#' @param model Character, one of MODELS.
#' @export
predictors_for <- function(el, model) {
  if (!model %in% names(el$models)) {
    stop("unknown model '", model, "'. Known: ",
         paste(names(el$models), collapse = ", "), call. = FALSE)
  }
  c(el$baseline, el$models[[model]])
}


# --- Eligible set ------------------------------------------------------------

#' Load the human dataset and cut it to the common eligible analysis set
#'
#' The guarantee provided here rather than downstream: both models are handed
#' the SAME rows and the same gene ids, and neither can lose rows to the extra
#' missingness in the structure block.
#'
#' ROW POLICY: every gene with a target and a split. XGBoost learns a default
#' split direction for NAs, so nothing is imputed and no row is discarded —
#' 13,601 genes on the human v10 cache. Eligibility does not depend on either
#' model's columns, which is what keeps the rows identical.
#'
#' The cost of that policy, to state rather than bury: structure missingness is
#' INFORMATIVE. A missing 5'UTR MFE means a 5'UTR too short to fold, not a
#' failed computation. `Structure` can therefore split on an annotation
#' artefact and bank it as a structure effect, and this design cannot rule that
#' out. report_feature_sets() prints the share of genes carrying at least one
#' missing predictor so the size of the channel is visible on every run.
#'
#' Zero-variance baseline columns are removed here, on train+val only, so the
#' predictor sets are fixed before any model sees them and the two models
#' cannot end up with different baseline blocks via a recipe filter.
#'
#' @param species Character, passed to build_dataset().
#' @return list(data, baseline, models, structure, probing, dropped_zv)
#' @export
eligible_dataset <- function(species = "human") {

  df <- build_dataset(species) |>
    drop_excluded(verbose = FALSE) |>
    attach_splits()

  base_cols <- baseline_columns(df)
  str_cols  <- structure_columns(df)
  gini_cols <- probing_columns(df)

  stopifnot(length(intersect(base_cols, str_cols)) == 0,
            length(intersect(base_cols, gini_cols)) == 0,
            length(intersect(str_cols, gini_cols)) == 0,
            !TARGET_COL %in% c(base_cols, str_cols),
            length(intersect(META_COLS, c(base_cols, str_cols))) == 0)

  keep <- df |>
    filter(!is.na(.data[[TARGET_COL]]), !is.na(split))

  # Constant on the data the model may learn from. Checked on train+val rather
  # than on everything, because inspecting test to decide the predictor set is
  # a (mild) use of held-out data.
  # NAs dropped before counting: a column with one observed value and the rest
  # missing has two distinct values (v and NA) and would sneak past a naive
  # uniqueness test while carrying no information.
  learnable <- keep[keep$split != "test", , drop = FALSE]
  zv <- names(which(vapply(learnable[base_cols], function(x)
    length(unique(x[!is.na(x)])) < 2L, logical(1))))
  base_cols <- setdiff(base_cols, zv)

  models <- list(character(), str_cols)
  names(models) <- MODELS

  list(
    data       = keep,
    baseline   = base_cols,
    models     = models,
    structure  = str_cols,
    probing    = gini_cols,
    dropped_zv = zv
  )
}


#' Print the feature lists and the sample sizes
#' @export
report_feature_sets <- function(el) {
  d <- el$data
  cat("\n=== Eligible analysis set ===\n")
  cat(sprintf("Response          : %s (%s)\n",
              TARGET_COL, response_info()$note))
  cat(sprintf("Genes             : %d (1 row per gene, %d distinct gene_id)\n",
              nrow(d), dplyr::n_distinct(d$gene_id)))
  miss <- mean(!stats::complete.cases(
    d[, c(el$baseline, el$structure), drop = FALSE]))
  cat(sprintf("                    %.1f%% of them have at least one missing predictor\n",
              100 * miss))
  cat(sprintf("Split (blocked on family_id_%s):\n", BLOCK_LEVEL))
  print(table(d$split))

  cat("\n--- The two models ---\n")
  cat(sprintf("  %-18s %-11s %-11s %s\n",
              "model", "structure", "predictors", "structure groups"))
  for (m in MODELS) {
    cat(sprintf("  %-18s %-11d %-11d %s\n", m,
                length(el$models[[m]]),
                length(el$baseline) + length(el$models[[m]]),
                if (length(el$models[[m]])) paste(structure_groups(), collapse = ", ")
                else "—"))
  }

  cat(sprintf("\nBaseline features : %d\n", length(el$baseline)))
  if (length(el$dropped_zv)) {
    cat("Dropped (zero variance on train+val): ",
        paste(el$dropped_zv, collapse = ", "), "\n")
  }

  cat("\n--- Baseline block by family ---\n")
  for (g in names(FEATURE_PATTERNS)) {
    cc <- intersect(fg_columns(d, g), el$baseline)
    if (length(cc)) cat(sprintf("  %-16s %3d\n", g, length(cc)))
  }
  cat("  (non-model features and why: the Notes column of R/feature_table.csv)\n")

  cat("\n--- Structure block by family ---\n")
  for (g in structure_groups()) {
    cc <- intersect(fg_columns(d, g), el$structure)
    # A model feature that contributes nothing has lost its columns upstream;
    # flag it rather than show a bare 0, so it is visible on the first run.
    cat(sprintf("  %-16s %3d%s\n", g, length(cc),
                if (length(cc) == 0) "   (empty — no columns in this build)" else ""))
  }
  cat(sprintf("  (%s excluded from both models: %d columns, measured rather\n",
              PROBING_GROUP, length(el$probing)))
  cat("   than computed from sequence — a later, supplementary model)\n")
  invisible(el)
}
