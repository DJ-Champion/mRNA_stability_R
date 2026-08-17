# =============================================================================
# Shared feature-block, model-ladder and eligible-set definition
# =============================================================================
# The one place the baseline block, the structure ladder and the eligible row
# set are defined. xgb_structure_comparison.R (the main held-out comparison)
# and xgb_structure_gini_subset.R (the icSHAPE secondary run) both source this
# file, so the two runs cannot drift apart in what "baseline" or "structure"
# means.
#
# The nesting the ladder argument depends on — Baseline < S-core < S-select <
# S-full — is a property of THIS file: each rung's predictors are constructed
# as BASELINE_COLS + that rung's structure groups, by definition rather than by
# four parallel edits kept in step by hand. resolve_variant() asserts the
# nesting rather than trusting it.
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


# --- The model ladder --------------------------------------------------------
# Four nested models. Each rung is the previous rung plus one more family of
# folding metrics, so a difference between adjacent rungs is attributable to
# the family that was added and nothing else.
#
#   Baseline   non-structure transcript features only
#   S-core     + MFE z-scores (8) and local MFE z-scores (7)          = 15
#   S-select   + MFE delta, observed minus expected (7)               = 22
#   S-full     + raw MFE (8+7) and per-nucleotide MFE (7)             = 44
#
# WHY THIS ORDER. It is not arbitrary: the rungs are ordered from least to most
# confounded with the baseline, so the ladder walks from the cleanest test of
# structure to the dirtiest.
#
#   S-core     the z-scores are the only folding metrics normalised against
#              shuffled sequence. The earlier redundancy regression put them at
#              R^2 0.015-0.45 from baseline (median ~0.13) — genuinely new
#              numbers. This is the PRIMARY test.
#   S-select   adds mfe_delta_*, which that same regression reconstructed from
#              baseline at R^2 0.68-0.92. Delta is observed minus EXPECTED MFE,
#              and expected MFE is a deterministic function of GC and length,
#              both already in the baseline. Expect it to add little.
#   S-full     adds raw and per-nucleotide MFE, which are near-deterministic in
#              length and GC. If the ladder turns positive only here, the
#              likely explanation is length re-entering the model under a
#              structure label, not structure.
#
# So the ladder has a PREDICTED SHAPE under each hypothesis, which is what
# makes it evidence rather than three more chances to clear zero. If structure
# carries information, the gain should appear at S-core and persist. A gain
# that appears only at S-full is the confounding signature, not a discovery.

LADDER <- list(
  "Baseline" = character(),
  "S-core"   = c("rnafold_zscores", "rnalfold_zscores"),
  "S-select" = c("rnafold_zscores", "rnalfold_zscores", "mfe_deltas"),
  "S-full"   = c("rnafold_zscores", "rnalfold_zscores", "mfe_deltas",
                 "rnafold_scores", "rnalfold_scores", "rnafold_per_nt")
)

# The rung every other rung is compared against.
REFERENCE_MODEL <- "Baseline"

# PRE-SPECIFIED before any of these models were fitted. The ladder produces
# three contrasts against Baseline and three rung-to-rung increments on one
# held-out set; naming the primary in the code, in the registry that the
# summary table reads, is what stops the most publishable of the six becoming
# the headline after the fact. Everything else is secondary and unadjusted, and
# is reported that way.
PRIMARY_MODEL <- "S-core"


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
# eligibility are scored on different genes, so their absolute R^2 values are
# not on the same footing. The paired contrasts WITHIN a variant are always
# meaningful; the absolute numbers between variants are not. The summary table
# prints n for every variant for exactly this reason.
#
# Each entry lists only what differs; the rest is inherited from `default`.

VARIANTS <- list(

  default = list(
    label = "all genes with a target; XGBoost handles NA natively",
    # "native_na" = keep every gene with a target and a split, let XGBoost
    #               learn a default split direction for NAs
    # "complete"  = complete cases on baseline + the widest structure block
    #               (no imputation, informative missingness cannot leak)
    rows          = "native_na",
    ladder        = LADDER,
    baseline_add  = character(),
    baseline_drop = character()
  ),

  complete_case = list(
    label = "complete cases on baseline + all structure columns",
    # The sensitivity arm for the default's row policy, and the ONLY run in
    # which the informative-missingness channel is closed by construction.
    #
    # Structure missingness is informative: a missing 5'UTR MFE means a 5'UTR
    # too short to fold, not a failed computation. Under `native_na` a
    # structure rung can therefore split on an annotation artefact and bank it
    # as a structure effect. Under complete cases it cannot, because the rows
    # carrying the artefact are gone from every rung equally.
    #
    # Costs ~13% of the corpus (1,800 of 13,601 genes on the v9 human cache).
    # If the default and this variant disagree, the missingness channel is the
    # first explanation to rule out, not a discovery.
    rows = "complete"
  )
)


#' Resolve a variant name into a complete specification
#'
#' Also validates the ladder, because every downstream claim rests on its
#' shape: the first rung must be the pure reference, and each rung must be a
#' superset of the one before it. A ladder that is not nested would still run
#' and still produce intervals, but the increments would no longer mean "what
#' this family of features added" — so this is checked, not assumed.
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

  lad <- v$ladder
  if (!is.list(lad) || length(lad) < 2 || is.null(names(lad))) {
    stop("variant '", name, "' has a malformed ladder", call. = FALSE)
  }
  if (!identical(names(lad)[[1]], REFERENCE_MODEL) ||
      length(lad[[1]]) != 0) {
    stop("variant '", name, "': the first rung must be '", REFERENCE_MODEL,
         "' with no structure groups", call. = FALSE)
  }
  for (i in seq_len(length(lad) - 1)) {
    if (!all(lad[[i]] %in% lad[[i + 1]])) {
      stop("variant '", name, "': ladder is not nested — '", names(lad)[i],
           "' is not a subset of '", names(lad)[i + 1], "'", call. = FALSE)
    }
  }
  if (!PRIMARY_MODEL %in% names(lad)) {
    stop("variant '", name, "': primary model '", PRIMARY_MODEL,
         "' is not a rung of the ladder", call. = FALSE)
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


#' The ordered contrasts the ladder implies
#'
#' Two families, both paired and both computed on the same bootstrap draws:
#'
#'   vs_baseline  each structure rung against the reference. Answers "does this
#'                much structure beat no structure at all".
#'   increment    each rung against the rung below it. Answers "what did THIS
#'                family of features add, given everything below it". This is
#'                the family that makes the ladder an argument rather than
#'                three separate comparisons.
#'
#' @param variant A resolved variant.
#' @return A tibble with columns kind, contrast, lhs, rhs, primary.
#' @export
ladder_contrasts <- function(variant = resolve_variant()) {
  nm   <- names(variant$ladder)
  rest <- setdiff(nm, REFERENCE_MODEL)

  vs_base <- tibble::tibble(
    kind     = "vs_baseline",
    lhs      = rest,
    rhs      = REFERENCE_MODEL
  )
  incr <- tibble::tibble(
    kind     = "increment",
    lhs      = nm[-1],
    rhs      = nm[-length(nm)]
  )

  # The first increment IS the first vs-baseline contrast — rung two against
  # the reference — so it would otherwise appear twice, once flagged primary
  # and once not. The same number turning up in both the headline row and the
  # "secondary, unadjusted" list is exactly the kind of double-count a reviewer
  # is right to be suspicious of, so it is removed here rather than filtered
  # in each consumer. vs_baseline wins because that is the framing the primary
  # contrast is specified in.
  incr <- dplyr::anti_join(incr, vs_base, by = c("lhs", "rhs"))

  dplyr::bind_rows(vs_base, incr) |>
    dplyr::mutate(
      contrast = paste(lhs, "vs", rhs),
      # The single pre-specified comparison; see PRIMARY_MODEL above.
      primary  = kind == "vs_baseline" & lhs == PRIMARY_MODEL
    ) |>
    dplyr::select(kind, contrast, lhs, rhs, primary)
}


# --- Feature-block construction ----------------------------------------------

#' Build the baseline (non-structure) column list for a dataset
#'
#' Deliberately EXCLUDED, and why:
#'
#'   aa_freqs      (aa_*, 20 cols)  an EXACT deterministic function of the
#'                                  codon columns that are kept. Verified on
#'                                  the v9 human cache to 5.6e-17:
#'                                    aa_x = sum(codons encoding x)
#'                                           / (1 - stop_fraction - codon_other)
#'                                  So they add no information, and 20 exact
#'                                  substitutes for retained columns splinter
#'                                  gain importance and shrink the share of any
#'                                  column subsample that the structure block
#'                                  can occupy. Dropped from the MODEL baseline
#'                                  only — they stay in the cache and in
#'                                  EXCLUDED_FEATURES' complement, because the
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
#'   every `structure` supergroup member   that is the experimental variable.
#'
#' Kept beyond the obvious list: `exon_length_last_mrna`, the only surviving
#' member of the `exons` group. It is non-structure transcript architecture
#' (last-exon length is the classical NMD 50-nt-rule covariate), so it belongs
#' with the junction-distance features.
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


#' Build the structure column list for one rung of the ladder
#'
#' icSHAPE Gini is NOT here — see `gini_columns()` and the secondary run. It is
#' 80-91% missing, and a common complete analysis set containing it costs 93%
#' of the corpus.
#'
#' @param df A dataset from build_dataset() after drop_excluded().
#' @param variant A resolved variant.
#' @param model Character, a rung name (a key of `variant$ladder`).
#' @return Character vector of column names present in `df`.
#' @export
structure_columns <- function(df, variant = resolve_variant(),
                              model = PRIMARY_MODEL) {
  if (!model %in% names(variant$ladder)) {
    stop("unknown model '", model, "'. Ladder: ",
         paste(names(variant$ladder), collapse = ", "), call. = FALSE)
  }
  unique(unlist(lapply(variant$ladder[[model]],
                       function(g) fg_columns(df, g))))
}


#' Every structure column any rung uses — the widest block
#'
#' Because the ladder is nested this is the top rung's block, but deriving it
#' from the union rather than from `ladder[[length(ladder)]]` keeps it correct
#' if a variant ever defines a ladder that widens non-monotonically in a way
#' the nesting check still permits.
#' @export
all_structure_columns <- function(df, variant = resolve_variant()) {
  unique(unlist(lapply(names(variant$ladder),
                       function(m) structure_columns(df, variant, m))))
}


#' The icSHAPE structural-Gini block (secondary analysis only)
#' @export
gini_columns <- function(df) fg_columns(df, "probing")


#' The predictor list for one rung
#' @param el Result of eligible_dataset().
#' @param model Character, a rung name.
#' @export
predictors_for <- function(el, model) {
  if (!model %in% names(el$models)) {
    stop("unknown model '", model, "'. Ladder: ",
         paste(names(el$models), collapse = ", "), call. = FALSE)
  }
  c(el$baseline, el$models[[model]])
}


# --- Eligible set ------------------------------------------------------------

#' Load the human dataset and cut it to the common eligible analysis set
#'
#' Whatever the row policy, the guarantee is the same and is provided here
#' rather than downstream: every rung of the ladder is handed the SAME rows and
#' the same gene ids, and no rung can lose rows to the extra missingness in its
#' own block. Eligibility is decided once, from the union of every block, before
#' any model exists — which is why `complete` screens on
#' all_structure_columns() and not on the primary rung's block.
#'
#' Two policies, selected by the variant:
#'
#'   "native_na" (default) every gene with a target and a split; XGBoost
#'               handles NA internally. Keeps all 13,601. Every rung still sees
#'               identical rows, because eligibility no longer depends on any
#'               block. The informative-missingness concern is real here and is
#'               what the `complete_case` variant exists to quantify.
#'
#'   "complete"  complete cases on baseline + every structure column. Closes
#'               the missingness channel by construction. Costs 1,800 of 13,601
#'               genes (13%).
#'
#' Zero-variance baseline columns are removed here, on train+val only, so the
#' predictor sets are fixed before any model sees them and the rungs cannot end
#' up with different baseline blocks via a recipe filter.
#'
#' @param species Character, passed to build_dataset().
#' @param variant A resolved variant (see resolve_variant()).
#' @return list(data, baseline, models, structure_all, gini, dropped_zv, variant)
#' @export
eligible_dataset <- function(species = "human", variant = resolve_variant()) {

  df <- build_dataset(species) |>
    drop_excluded(verbose = FALSE) |>
    attach_splits()

  base_cols <- baseline_columns(df, variant)
  all_str   <- all_structure_columns(df, variant)
  gini_cols <- gini_columns(df)

  stopifnot(length(intersect(base_cols, all_str)) == 0,
            length(intersect(base_cols, gini_cols)) == 0,
            !TARGET_COL %in% c(base_cols, all_str, gini_cols),
            length(intersect(META_COLS, c(base_cols, all_str))) == 0)

  keep <- df |>
    filter(!is.na(.data[[TARGET_COL]]), !is.na(split))

  if (identical(variant$rows, "complete")) {
    ok   <- stats::complete.cases(keep[, c(base_cols, all_str), drop = FALSE])
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

  models <- lapply(names(variant$ladder),
                   function(m) structure_columns(keep, variant, m))
  names(models) <- names(variant$ladder)

  list(
    data          = keep,
    baseline      = base_cols,
    models        = models,
    structure_all = all_str,
    gini          = gini_cols,
    dropped_zv    = zv,
    variant       = variant
  )
}


#' Print the feature lists, the ladder and the sample sizes
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
      d[, c(el$baseline, el$structure_all), drop = FALSE]))
    cat(sprintf("                    %.1f%% of them have at least one missing predictor\n",
                100 * miss))
  }
  cat(sprintf("Split (blocked on family_id_%s):\n", BLOCK_LEVEL))
  print(table(d$split))

  cat(sprintf("\nBaseline features : %d\n", length(el$baseline)))
  if (length(el$dropped_zv)) {
    cat("Dropped (zero variance on train+val): ",
        paste(el$dropped_zv, collapse = ", "), "\n")
  }

  cat("\n--- The ladder ---\n")
  cat(sprintf("  %-10s %-11s %-11s %s\n",
              "rung", "structure", "predictors", "structure groups"))
  for (m in names(el$models)) {
    star <- if (identical(m, PRIMARY_MODEL)) " *" else "  "
    cat(sprintf("%s%-10s %-11d %-11d %s\n", star, m,
                length(el$models[[m]]),
                length(el$baseline) + length(el$models[[m]]),
                paste(el$variant$ladder[[m]], collapse = ", ")))
  }
  cat("  * primary, pre-specified\n")

  cat("\n--- Structure columns by rung ---\n")
  for (m in setdiff(names(el$models), REFERENCE_MODEL)) {
    added <- setdiff(el$models[[m]], el$models[[which(names(el$models) == m) - 1]])
    cat(sprintf("  %-10s adds %2d: %s\n", m, length(added),
                paste(added, collapse = ", ")))
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
  cat("   retained columns — see baseline_columns())\n")
  invisible(el)
}
