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


# --- Variant registry --------------------------------------------------------
# A variant is a named, declarative deviation from the committed default. It
# may change TWO things and nothing else: which rows are eligible, and which
# columns go in each block. Resampling, the split artefact, the tuning grid and
# budget, and the seed are deliberately NOT variant-controlled, so any two
# variants differ in exactly the stated way and a difference between them is
# attributable.
#
# WHY THIS IS A SENSITIVITY ANALYSIS, NOT A SET OF TESTS. Running many
# specifications and reporting the one that clears zero is the garden of
# forking paths, and it would destroy the value of the committed result — which
# is defensible precisely because it was specified before anyone saw an answer.
# Every variant here must be reported, every time, in the table produced by
# xgb_structure_variant_summary.R. The claim to make is "the conclusion holds
# across every reasonable specification", which is stronger than a single run,
# not weaker. If one variant disagrees with the rest, that is a finding about
# the specification, to be explained rather than promoted.
#
# CROSS-VARIANT METRICS ARE NOT ALWAYS COMPARABLE. Variants that change row
# eligibility are scored on different genes, so their absolute R² values are
# not on the same footing. The paired DELTA within a variant is always
# meaningful; the absolute numbers between variants are not. The summary table
# prints n for every variant for exactly this reason.
#
# Each entry lists only what differs; the rest is inherited from `default`.

VARIANTS <- list(

  default = list(
    label = "complete cases; length-normalised structure only",
    # "complete"  = complete cases on baseline + structure (identical rows,
    #               no imputation, informative missingness cannot leak)
    # "native_na" = keep every gene with a target and a split, let XGBoost
    #               learn a default split direction for NAs
    rows             = "complete",
    structure_groups = c("rnafold_zscores", "rnalfold_zscores", "mfe_deltas"),
    baseline_add     = character(),
    baseline_drop    = character()
  ),

  keep_missing = list(
    label = "all genes with a target; XGBoost handles NA natively",
    # Tests whether the default's 13% complete-case loss changed the answer.
    # Read the result with care in ONE direction: structure missingness is
    # informative (a missing 5'UTR MFE means a 5'UTR too short to fold), so
    # this variant can flatter Model B by letting it split on an annotation
    # artefact. If it is the only variant that turns positive, that is the
    # first explanation to rule out, not a discovery.
    rows = "native_na"
  ),

  with_raw_mfe = list(
    label = "complete cases; adds raw and per-nucleotide MFE to structure",
    # Tests whether excluding length-confounded MFE hid a real signal. Same
    # caution: raw MFE is near-deterministic in length and GC, both already in
    # the baseline, so if this variant alone turns positive the likely
    # explanation is length re-entering the model, not structure.
    structure_groups = c("rnafold_zscores", "rnalfold_zscores", "mfe_deltas",
                         "rnafold_scores", "rnalfold_scores", "rnafold_per_nt")
  )
)


#' Resolve a variant name into a complete specification
#'
#' @param name Character, a key of VARIANTS.
#' @return The variant list, with every field filled in from `default`.
#' @export
resolve_variant <- function(name = "default") {
  if (!name %in% names(VARIANTS)) {
    stop("unknown variant '", name, "'. Known: ",
         paste(names(VARIANTS), collapse = ", "), call. = FALSE)
  }
  v <- utils::modifyList(VARIANTS$default, VARIANTS[[name]])
  if (!v$rows %in% c("complete", "native_na")) {
    stop("variant '", name, "' has an unknown rows policy: ", v$rows,
         call. = FALSE)
  }
  v$name <- name
  v
}


#' Where a variant's artefacts live. One self-contained directory per variant.
#' @export
variant_dir <- function(variant, what = c("root", "tables", "plots")) {
  what <- match.arg(what)
  root <- file.path(OUTPUT_DIR, "xgb_structure",
                    if (is.list(variant)) variant$name else variant)
  if (what == "root") root else file.path(root, what)
}


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
#' @param variant A resolved variant; `baseline_add` / `baseline_drop` adjust
#'   the list below without editing it.
#' @return Character vector of column names present in `df`.
#' @export
baseline_columns <- function(df, variant = resolve_variant()) {
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
  cols <- unique(c(cols, intersect(variant$baseline_add, names(df))))
  setdiff(cols, variant$baseline_drop)
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
#' @param variant A resolved variant; `structure_groups` names the
#'   FEATURE_PATTERNS keys that make up the block. The default is
#'   rnafold_zscores (8, MFE z-score per region), rnalfold_zscores (7, local
#'   MFE z-score) and mfe_deltas (7, observed minus expected MFE).
#' @return Character vector of column names present in `df`.
#' @export
structure_columns <- function(df, variant = resolve_variant()) {
  unique(unlist(lapply(variant$structure_groups, function(g) fg_columns(df, g))))
}


#' The icSHAPE structural-Gini block (secondary analysis only)
#' @export
gini_columns <- function(df) fg_columns(df, "probing")


# --- Eligible set ------------------------------------------------------------

#' Load the human dataset and cut it to the common eligible analysis set
#'
#' Whatever the row policy, the guarantee the brief requires is the same and is
#' provided here rather than downstream: Model A and Model B are handed the
#' SAME rows and the same gene ids, and Model B cannot lose rows to the extra
#' missingness in its own block. Eligibility is decided once, from the union of
#' the two blocks, before either model exists.
#'
#' Two policies, selected by the variant:
#'
#'   "complete"  (default) complete cases on baseline + structure. Chosen
#'               because the structure columns' missingness is INFORMATIVE — a
#'               missing 5'UTR MFE means a 5'UTR too short to fold, not a failed
#'               computation — so letting XGBoost learn a default split
#'               direction there would let Model B profit from an annotation
#'               artefact and report it as a structure effect. Costs 1,800 of
#'               13,601 genes (13%).
#'
#'   "native_na" every gene with a target and a split; XGBoost handles NA
#'               internally. Keeps all 13,601. Both models still see identical
#'               rows, because eligibility no longer depends on either block.
#'               The informative-missingness concern above is precisely what
#'               this variant exists to quantify — so read a gain here as a
#'               question, not an answer.
#'
#' Zero-variance baseline columns are removed here, on train+val only, so the
#' predictor sets are fixed before any model sees them and the two models
#' cannot end up with different baseline blocks via a recipe filter.
#'
#' @param species Character, passed to build_dataset().
#' @param variant A resolved variant (see resolve_variant()).
#' @return list(data, baseline, structure, gini, dropped_zv, variant)
#' @export
eligible_dataset <- function(species = "human", variant = resolve_variant()) {

  df <- build_dataset(species) |>
    drop_excluded(verbose = FALSE) |>
    attach_splits()

  base_cols  <- baseline_columns(df, variant)
  str_cols   <- structure_columns(df, variant)
  gini_cols  <- gini_columns(df)

  stopifnot(length(intersect(base_cols, str_cols)) == 0,
            length(intersect(base_cols, gini_cols)) == 0,
            !TARGET_COL %in% c(base_cols, str_cols, gini_cols),
            length(intersect(META_COLS, c(base_cols, str_cols))) == 0)

  keep <- df |>
    filter(!is.na(.data[[TARGET_COL]]), !is.na(split))

  if (identical(variant$rows, "complete")) {
    ok   <- stats::complete.cases(keep[, c(base_cols, str_cols), drop = FALSE])
    keep <- keep[ok, , drop = FALSE]
  }

  # Constant on the data the model may learn from. Checked on train+val rather
  # than on everything, because inspecting test to decide the predictor set is
  # a (mild) use of held-out data.
  # NAs dropped before counting: under the "native_na" policy a column with one
  # observed value and the rest missing has two distinct values (v and NA) and
  # would sneak past a naive uniqueness test while carrying no information.
  learnable <- keep[keep$split != "test", , drop = FALSE]
  zv <- names(which(vapply(learnable[base_cols], function(x)
    length(unique(x[!is.na(x)])) < 2L, logical(1))))
  base_cols <- setdiff(base_cols, zv)

  list(
    data       = keep,
    baseline   = base_cols,
    structure  = str_cols,
    gini       = gini_cols,
    dropped_zv = zv,
    variant    = variant
  )
}


#' Print the feature lists and sample sizes (checklist item 11)
#' @export
report_feature_sets <- function(el) {
  d <- el$data
  cat("\n=== Eligible analysis set ===\n")
  cat(sprintf("Variant           : %s — %s\n", el$variant$name, el$variant$label))
  cat(sprintf("Row policy        : %s\n", el$variant$rows))
  cat(sprintf("Response          : %s (Agarwal & Kelley consensus PC1, untransformed)\n",
              TARGET_COL))
  cat(sprintf("Genes             : %d (1 row per gene, %d distinct gene_id)\n",
              nrow(d), dplyr::n_distinct(d$gene_id)))
  if (identical(el$variant$rows, "native_na")) {
    miss <- mean(!stats::complete.cases(
      d[, c(el$baseline, el$structure), drop = FALSE]))
    cat(sprintf("                    %.1f%% of them have at least one missing predictor\n",
                100 * miss))
  }
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
