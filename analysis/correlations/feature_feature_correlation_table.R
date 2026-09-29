# =============================================================================
# Feature-feature pairwise Spearman correlation table
# =============================================================================
# Pairwise Spearman correlations between all candidate predictor features
# (feature-vs-feature, not feature-vs-response). Intended as a preliminary
# redundancy triage pass before modelling.
#
# Output: one long-form TSV per species, sorted by descending |rho|.
# Columns: species, feature_a, display_a, feature_b, display_b,
#          group_a, group_b, region_a, region_b, n, rho.
#
# Excluded by default: codon_freqs and aa_freqs (within-group correlations
# are a mathematical artefact — each family is row-normalised to sum to
# exactly 1, so the columns are compositional and share a fixed total), plus
# the response variable halflife and derived prediction columns (R10).
#
# That closure held only accidentally before v8: a DNA-only regex in
# add_codon_aa_fractions() normalised just the 27 codons spelled without a U,
# which summed to 1 among themselves while the other 38 stayed raw counts.
# The exclusion was right for the wrong reason; it is now right for the
# stated one.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_feature_correlation_table.R")
#   df  <- build_dataset("human")
#   out <- compute_feature_correlation_table(df)
#   head(out$table)
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})


# Columns excluded regardless of group membership.
# Mirrors R10: drop the response and any derived predictions.
.FF_EXCLUDE_PATTERNS <- c(
  "^halflife$",
  "^saluki_prediction$",
  "^prediction_difference$"
)

# Pipeline identifier columns — never treated as features. From META_COLS
# (config.R), which also covers the family / split blocking columns.
.FF_ID_COLS <- META_COLS


# Internal: build column → group and column → region lookup maps from a
# select_features_v2() result. group = the feature id the column was selected
# under. region = the last underscore-delimited token when it is a legal
# REGIONS member; NA_character_ otherwise (e.g. standalone scalars like `cai`
# carry no region suffix).
.ff_build_col_maps <- function(sel) {
  tokens_last <- sub("^.*_", "", sel$column)
  has_region  <- tokens_last %in% REGIONS & grepl("_", sel$column, fixed = TRUE)
  list(
    group  = stats::setNames(sel$feature_id, sel$column),
    region = stats::setNames(ifelse(has_region, tokens_last, NA_character_),
                             sel$column)
  )
}


#' Compute pairwise Spearman correlations between all candidate predictor features.
#'
#' Feature-vs-feature only — no correlations against the response variable.
#' Uses stats::cor() with use="pairwise.complete.obs" for speed (vectorised,
#' not looped). Per-pair n derived via crossprod on a 0/1 non-NA indicator
#' matrix (O(n*p) to build + O(n*p^2) for crossprod — no O(p^2) loop needed).
#'
#' @param df      A dataframe from build_dataset(). Must contain exactly one
#'                species (pipeline invariant: pass build_dataset(sp), not
#'                build_all(), to this function — the runner loops species).
#' @param include Selection tokens (see select_features_v2()); default
#'                "exploratory".
#' @param exclude Tokens to subtract from `include`. Default codon_freqs and
#'                aa_freqs — their within-group correlations are a
#'                mathematical artefact (their values sum to ~1 within a
#'                region). Pass exclude = NULL to opt in.
#' @param method  Correlation method passed to stats::cor(). Default "spearman".
#'                Exposed as a parameter for consistency with
#'                feature_correlation_dotplot.R's signature; no CI machinery is
#'                used here — this is a cheaper triage pass.
#' @param pool    Logical. FALSE (default) — do not pool species. If df contains
#'                more than one species, an error is raised: call the function
#'                once per species and let the runner loop. TRUE pools all rows
#'                regardless of species and labels the output species = "pooled".
#'                Must be explicit opt-in; never the default.
#' @return list(plot = NULL, table = <tibble>). R9: no meaningful visual for a
#'   pairwise table of hundreds of pairs; plot = NULL signals this. Table has
#'   one row per unordered feature pair (upper triangle, no diagonal, no A-B
#'   and B-A duplicates), sorted descending by abs(rho). Columns:
#'   species, feature_a, display_a, feature_b, display_b, group_a, group_b,
#'   region_a, region_b, n, rho. display_a/display_b are the format_col_name()
#'   display strings (R4) matching what appears on plot axes.
#' @export
compute_feature_correlation_table <- function(df,
                                              include = "exploratory",
                                              exclude = c("codon_freqs", "aa_freqs"),
                                              method = "spearman",
                                              pool   = FALSE) {

  # --- R5: guard pipeline invariants ----------------------------------------
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  if (!pool && length(unique(df$species)) > 1) {
    stop(
      "df contains multiple species (",
      paste(unique(df$species), collapse = ", "), "). ",
      "Call compute_feature_correlation_table() separately per species ",
      "(the runner loops names(SPECIES_CONFIG)), or set pool = TRUE to ",
      "pool all rows (output labelled species = 'pooled')."
    )
  }

  sp_label <- if (pool) "pooled" else unique(df$species)[1]

  # --- Select candidate columns (R3) ----------------------------------------
  sel <- select_features_v2(df, include, exclude)
  if (nrow(sel) == 0) {
    stop("`include` / `exclude` resolved to zero columns")
  }
  maps           <- .ff_build_col_maps(sel)
  candidate_cols <- names(maps$group)

  # R5: keep only columns that are actually in df (all-NA cols are dropped by
  # engineer.R; a column from a sibling species may simply not exist here).
  candidate_cols <- intersect(candidate_cols, names(df))

  # Keep only numeric columns (binary / factor predictors are out of scope).
  is_numeric <- vapply(candidate_cols,
                       function(co) is.numeric(df[[co]]),
                       logical(1))
  candidate_cols <- candidate_cols[is_numeric]

  # R10 + identifier exclusions.
  is_excluded <- vapply(candidate_cols, function(co) {
    co %in% .FF_ID_COLS ||
      any(vapply(.FF_EXCLUDE_PATTERNS,
                 function(rgx) grepl(rgx, co),
                 logical(1)))
  }, logical(1))
  candidate_cols <- candidate_cols[!is_excluded]

  if (length(candidate_cols) < 2) {
    stop("Fewer than 2 numeric feature columns resolved for species '", sp_label,
         "'. Check include / exclude or data availability for this species.")
  }

  n_features <- length(candidate_cols)
  n_pairs    <- choose(n_features, 2L)
  message("  [", sp_label, "] ", n_features, " features → ", n_pairs, " pairs")

  # --- Build feature matrix -------------------------------------------------
  feature_mat <- as.matrix(df[, candidate_cols, drop = FALSE])

  # --- Per-pair n via crossprod on 0/1 indicator matrix --------------------
  # ind_mat[i, j] = 1 when df[i, col_j] is non-NA, 0 otherwise.
  # crossprod(ind_mat)[a, b] = sum of row-wise AND of columns a and b =
  # number of complete-observation pairs for (col_a, col_b).
  # This avoids an O(p^2) loop over pairs; O(n*p) to build + one BLAS call.
  ind_mat <- (!is.na(feature_mat)) + 0L   # n × p integer 0/1
  n_mat   <- crossprod(ind_mat)            # p × p pairwise non-NA counts

  # --- Vectorised pairwise correlation --------------------------------------
  rho_mat <- stats::cor(feature_mat,
                        method = method,
                        use    = "pairwise.complete.obs")

  # --- Extract upper triangle: each unordered pair exactly once -------------
  idx <- which(upper.tri(rho_mat), arr.ind = TRUE)

  col_a <- candidate_cols[idx[, 1]]
  col_b <- candidate_cols[idx[, 2]]

  result <- tibble::tibble(
    species   = sp_label,
    feature_a = col_a,
    display_a = format_col_name(col_a),
    feature_b = col_b,
    display_b = format_col_name(col_b),
    group_a   = unname(maps$group[col_a]),
    group_b   = unname(maps$group[col_b]),
    region_a  = unname(maps$region[col_a]),
    region_b  = unname(maps$region[col_b]),
    n         = as.integer(n_mat[idx]),
    rho       = rho_mat[idx]
  )

  # Drop constant-column pairs (cor() returns NaN when a column is invariant
  # within the pairwise-complete subset).
  n_nan <- sum(is.nan(result$rho))
  if (n_nan > 0) {
    message("  [", sp_label, "] Dropping ", n_nan,
            " pair(s) where rho = NaN (constant column in pairwise subset)")
    result <- result[!is.nan(result$rho), ]
  }

  result <- dplyr::arrange(result, dplyr::desc(abs(rho)))

  # R9: no meaningful plot for a p*(p-1)/2 row correlation table.
  list(plot = NULL, table = result)
}


# --- Top-to-bottom runner (Rule §6.1 step 5) ---------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  for (sp in names(SPECIES_CONFIG)) {
    message("\n--- ", sp, " ---")

    df <- build_dataset(sp)

    out <- compute_feature_correlation_table(df)
    tbl <- out$table

    out_path <- file.path(
      OUTPUT_DIR, "tables",
      paste0("feature_feature_correlation_", sp, ".tsv")
    )

    # R8: output under OUTPUT_DIR/tables.
    # write.table with sep="\t" (TSV) — no row names, NA rendered as "NA".
    write.table(tbl, out_path,
                sep       = "\t",
                row.names = FALSE,
                quote     = FALSE,
                na        = "NA")

    message("  Wrote ", nrow(tbl), " pairs → ", out_path)

    top5 <- head(tbl, 5)
    message("  Top 5 pairs by |rho|:")
    for (i in seq_len(nrow(top5))) {
      message(sprintf("    %+.4f  %-40s  vs  %s",
                      top5$rho[i], top5$feature_a[i], top5$feature_b[i]))
    }
  }

  message("\nDone. Tables: ",
          file.path(OUTPUT_DIR, "tables"),
          "/feature_feature_correlation_<species>.tsv")
}
