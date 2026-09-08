# =============================================================================
# Shared feature-block and eligible-set definition for the two XGBoost models
# =============================================================================
# The one place the two models are defined:
#
#   Baseline    non-structure transcript features
#   Structure   Baseline + every computed secondary-structure feature
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
# `halflife` is PC1 of the Agarwal & Kelley (2022) consensus half-life measure,
# not a duration in hours: it is already a signed, roughly symmetric score
# (range -17.2 to +18.2, sd 4.81 on the human v10 cache). No project code
# transforms it, and a log is undefined on a variable that takes negative
# values, so it is modelled RAW. RMSE and MAE are therefore in PC1 units.
TARGET_COL <- "halflife"


# --- The two models ----------------------------------------------------------

REFERENCE_MODEL <- "Baseline"
STRUCTURE_MODEL <- "Structure"

# Order matters: reference first. Every table, figure and factor level in the
# downstream scripts reads this vector rather than restating the order.
MODELS <- c(REFERENCE_MODEL, STRUCTURE_MODEL)

# The `structure` supergroup member that is NOT part of the structure block.
# icSHAPE Gini is 80-91% missing on the human cache, and — more importantly —
# it is MEASURED rather than computed, so a model containing it cannot score a
# transcript nobody has probed. Anything using it is a different kind of model
# and belongs in its own, later, supplementary analysis.
PROBING_GROUP <- "probing"


#' The FEATURE_PATTERNS keys that make up the structure block.
#'
#' Derived from SUPERGROUPS rather than hand-listed, so a folding family added
#' to the schema joins the structure block automatically instead of silently
#' sitting in neither model. A function, not a constant, so it resolves
#' SUPERGROUPS at call time and this file can be sourced in any order.
#' @export
structure_groups <- function() setdiff(SUPERGROUPS$structure, PROBING_GROUP)


#' Where this analysis's artefacts live.
#' @param what "root", "tables" or "plots".
#' @export
run_dir <- function(what = c("root", "tables", "plots")) {
  what <- match.arg(what)
  root <- file.path(OUTPUT_DIR, "xgb_structure")
  if (what == "root") root else file.path(root, what)
}


# --- Feature-block construction ----------------------------------------------

#' Build the baseline (non-structure) column list for a dataset
#'
#' Deliberately EXCLUDED, and why:
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
#'                                  correlation plots and the intrinsic_select
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
  c(
    fg_columns(df, "lengths"),       # 4   regional sequence length
    fg_columns(df, "gc"),            # 7   regional GC content
    fg_columns(df, "skews"),         # 14  AT-skew and GC-skew
    fg_columns(df, "codon_freqs"),   # 65  coding composition
    intersect("cai", names(df)),     # 1   CAI (standalone; TE deliberately not)
    fg_columns(df, "stopfree"),      # 4   stop-free length
    fg_columns(df, "uorfs"),         # 1   uorf_present_mrna (numeric 0/1)
    fg_columns(df, "exon_density"),  # 4   CDS-exon density
    fg_columns(df, "eej_dist"),      # 2   junction distance
    fg_columns(df, "nmd")            # 5   fragile-codon / alternative-stop
    # `exons` (exon_length_last_mrna) deliberately omitted — 3'UTR-length
    # proxy, see the exclusion list above.
  ) |> unique()
}


#' Build the structure column list
#'
#' Every computed folding feature: RNAfold and RNALfold MFE and z-scores, the
#' per-nucleotide normalisation, and MFE delta. `mfe_expected` is a member of
#' the supergroup but contributes nothing here — drop_excluded() removes it as
#' engineering scaffolding for mfe_delta_*, so fg_columns() finds none of it.
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
  unique(unlist(lapply(structure_groups(), function(g) fg_columns(df, g))))
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
  cat(sprintf("Response          : %s (Agarwal & Kelley consensus PC1, untransformed)\n",
              TARGET_COL))
  cat(sprintf("Genes             : %d (1 row per gene, %d distinct gene_id)\n",
              nrow(d), dplyr::n_distinct(d$gene_id)))
  miss <- mean(!stats::complete.cases(
    d[, c(el$baseline, el$structure), drop = FALSE]))
  cat(sprintf("                    %.1f%% of them have at least one missing predictor\n",
              100 * miss))
  cat(sprintf("Split (blocked on family_id_%s):\n", BLOCK_LEVEL))
  print(table(d$split))

  cat("\n--- The two models ---\n")
  cat(sprintf("  %-10s %-11s %-11s %s\n",
              "model", "structure", "predictors", "structure groups"))
  for (m in MODELS) {
    cat(sprintf("  %-10s %-11d %-11d %s\n", m,
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
  if ("cai" %in% el$baseline) {
    cat("  (standalone: cai; translation_efficiency excluded)\n")
  }
  cat("  (aa_freqs, frac_*, purine_/amino_* excluded as exact functions of\n")
  cat("   retained columns; exon_length_last_mrna excluded as a 3'UTR-length\n")
  cat("   proxy at rho 0.949 — see baseline_columns())\n")

  cat("\n--- Structure block by family ---\n")
  for (g in structure_groups()) {
    cc <- intersect(fg_columns(d, g), el$structure)
    # A registered folding family that contributes nothing is expected for
    # mfe_expected — drop_excluded() removes it as scaffolding for mfe_delta_*
    # — but is worth flagging rather than showing as a bare 0, so a family that
    # empties for any OTHER reason is visible on the run that first does it.
    cat(sprintf("  %-16s %3d%s\n", g, length(cc),
                if (length(cc) == 0) "   (empty — removed by drop_excluded)" else ""))
  }
  cat(sprintf("  (%s excluded from both models: %d columns, measured rather\n",
              PROBING_GROUP, length(el$probing)))
  cat("   than computed from sequence — a later, supplementary model)\n")
  invisible(el)
}
