# =============================================================================
# Diagnostic feature-correlation heatmap workflow
# =============================================================================
#
# PURPOSE:
#   Generates reduced feature-feature correlation heatmaps from the mRNA
#   stability feature table for diagnostic exploration and candidate feature
#   selection. Outputs are labelled as diagnostic/candidate sets — NOT final
#   model feature sets.
#
# PLOT TYPES:
#   B. Target-filtered heatmap
#      Selects features by |Spearman rho| with half-life (top-N and/or
#      threshold). Shows redundancy among the features most associated with
#      the response.
#      When to use: diagnose whether top half-life features span biological
#      groups or cluster within one family (e.g. all length features).
#
#   C. Cluster-representative heatmap
#      Clusters all features by |rho| (distance = 1 - |rho|; average linkage).
#      One representative per cluster, chosen by highest |rho| with half-life
#      (ties broken by biological group priority, then alphabetically).
#      Shows a reduced non-redundant candidate feature set.
#      When to use: create a defensible baseline feature set; diagnose
#      redundancy before modelling.
#
#   D. Group-balanced heatmap
#      Top-K features per biological group (ranked by |rho| with half-life).
#      Prevents any one large feature family from dominating the visible set.
#      Preserves biological categories relevant to the paper thesis.
#      When to use: compare feature families and ensure all biological groups
#      are represented.
#
# FEATURE SELECTION:
#   Starts from INCLUDED_GROUPS (defined in config.R). All selection is
#   diagnostic/candidate. Do not automatically discard biologically central
#   corrected structure features (rnafold_zscores, mfe_deltas) purely
#   because their marginal correlation with half-life is weak.
#
# OUTPUTS:
#   Plots:  OUTPUT_DIR/plots/heatmaps/{prefix}_{type}_{species}.pdf
#   Tables: OUTPUT_DIR/tables/{prefix}_{table_name}_{species}.csv
#   Manifest: OUTPUT_DIR/tables/{prefix}_manifest.csv
#
# CHANGING THRESHOLDS:
#   Pass new values to run_correlation_heatmap_workflow():
#     top_n_target_features    <- c(30, 50, 100)   # sizes for Plot B
#     cluster_abs_rho_cutoff   <- 0.85             # redundancy for Plot C
#     top_k_per_group          <- 5                # features per group, Plot D
#     label_threshold          <- 0.3              # suppress labels < |rho|
#
# USAGE:
#   source("R/load_all.R")
#   source("analysis/correlations/correlation_heatmap_workflow.R")
#   df  <- build_dataset("human")
#   run_correlation_heatmap_workflow(df)
#
# CAUTION:
#   Feature selection decisions should integrate redundancy clusters,
#   biological category preservation, model performance, feature importance,
#   and known confounders (length, GC content). Do not use the reduced sets
#   from this workflow as automatic final model inputs.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
})


# =============================================================================
# Constants
# =============================================================================

.CW_ID_COLS <- c("transcript_id", "gene_id", "gene_name", "species")

# Priority for cluster representative tie-breaking (lower integer = preferred).
# Groups not listed get priority 999 (lowest priority).
.CW_REP_PRIORITY <- c(
  rnafold_zscores  = 1L,   # robust z-score / corrected structure
  mfe_deltas       = 2L,   # delta MFE
  mfe_expected     = 3L,   # expected MFE residual
  rnafold_per_nt   = 4L,   # MFE per nt (corrected structure family)
  rnaup            = 5L,   # RNAup / accessibility
  rnalfold_scores  = 6L,   # RNAlfold local structure
  rnalfold_zscores = 7L,   # RNAlfold z-score
  nmd              = 8L,   # NMD fragility
  junctions        = 9L,   # junctions
  eej_dist         = 10L,  # EEJ distances
  stopfree         = 11L,  # stop-free length
  gc               = 12L,  # GC content
  nuc_ratios       = 13L,  # nucleotide ratios
  nuc_combos       = 14L,  # purine / amino ratios
  skews            = 15L,  # AT/GC skew
  lengths          = 16L,  # length / architecture
  introns          = 17L,
  exons            = 18L,
  noncoding        = 19L
)


# =============================================================================
# Feature metadata helpers
# =============================================================================

#' Extract the region token from a feature column name.
#'
#' @param features Character vector of column names.
#' @return Character vector; NA_character_ where no recognised region suffix.
#' @export
parse_feature_region <- function(features) {
  vapply(features, function(f) {
    tokens <- strsplit(f, "_", fixed = TRUE)[[1]]
    last   <- tokens[length(tokens)]
    if (length(tokens) > 1L && last %in% REGIONS) last else NA_character_
  }, character(1), USE.NAMES = FALSE)
}


#' Assign each feature to a FEATURE_PATTERNS group key.
#'
#' First-group-wins when a column could match multiple patterns (order of
#' FEATURE_PATTERNS in config.R is authoritative).
#'
#' @param features         Character vector of column names.
#' @param feature_patterns Named list of regex patterns (default FEATURE_PATTERNS
#'                         from config.R).
#' @return Character vector of FEATURE_PATTERNS keys; NA_character_ for columns
#'   that match no pattern.
#' @export
assign_feature_group <- function(features, feature_patterns = FEATURE_PATTERNS) {
  vapply(features, function(f) {
    for (g in names(feature_patterns)) {
      if (grepl(feature_patterns[[g]], f)) return(g)
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}


#' Strip the region suffix from a feature column name to get the metric type.
#'
#' Delegates to format_metric_name() (R/utils/palettes.R) which applies
#' format_col_name() and then strips the formatted region display string.
#'
#' @param features Character vector of column names.
#' @return Character vector of human-readable metric type strings.
#' @export
parse_metric_type <- function(features) {
  format_metric_name(features)
}


#' Identify numeric feature columns, excluding identifiers and reserved columns.
#'
#' Optionally drops zero-variance columns and columns with fewer than
#' min_complete_n non-NA observations.
#'
#' @param data               Dataframe from build_dataset().
#' @param target_col         Response column to exclude from features.
#' @param id_cols            Identifier columns to exclude (default: pipeline
#'                           id columns).
#' @param drop_zero_variance Drop columns with zero variance (default TRUE).
#' @param min_complete_n     Minimum non-NA observations per column (default
#'                           100L). Columns with fewer are dropped.
#' @return Character vector of candidate feature column names.
#' @export
get_numeric_feature_cols <- function(data,
                                     target_col,
                                     id_cols            = .CW_ID_COLS,
                                     drop_zero_variance = TRUE,
                                     min_complete_n     = 100L) {
  # R10: also exclude derived predictions
  exclude_exact <- unique(c(id_cols, target_col,
                             "saluki_prediction", "prediction_difference"))

  cols <- names(data)
  cols <- cols[!cols %in% exclude_exact]
  cols <- cols[vapply(cols, function(co) is.numeric(data[[co]]), logical(1))]

  if (drop_zero_variance) {
    cols <- cols[vapply(cols, function(co) {
      x  <- data[[co]]
      ok <- !is.na(x)
      if (sum(ok) < 2L) return(FALSE)
      stats::var(x[ok]) > 0
    }, logical(1))]
  }

  if (min_complete_n > 0L) {
    cols <- cols[vapply(cols, function(co) {
      sum(!is.na(data[[co]])) >= as.integer(min_complete_n)
    }, logical(1))]
  }

  cols
}


# =============================================================================
# Correlation computation
# =============================================================================

#' Compute Spearman correlation of each feature with the target variable.
#'
#' @param data         Dataframe from build_dataset().
#' @param feature_cols Character vector of feature column names.
#' @param target_col   Response column name.
#' @param method       Correlation method (default "spearman").
#' @return Tibble: feature, target_col, spearman_rho, p_value, n_complete,
#'   abs_spearman_rho, feature_group, region, metric_type.
#' @export
compute_feature_target_correlations <- function(data,
                                                feature_cols,
                                                target_col,
                                                method = "spearman") {
  stopifnot(target_col %in% names(data))
  feature_cols <- intersect(feature_cols, names(data))

  y      <- data[[target_col]]
  grps   <- assign_feature_group(feature_cols)
  regs   <- parse_feature_region(feature_cols)
  mets   <- parse_metric_type(feature_cols)

  purrr::map_dfr(seq_along(feature_cols), function(i) {
    co   <- feature_cols[i]
    x    <- data[[co]]
    ok   <- !is.na(x) & !is.na(y)
    n    <- sum(ok)
    enough <- n >= 5L && length(unique(x[ok])) >= 2L

    if (!enough) {
      return(tibble::tibble(
        feature          = co,
        target_col       = target_col,
        spearman_rho     = NA_real_,
        p_value          = NA_real_,
        n_complete       = n,
        abs_spearman_rho = NA_real_,
        feature_group    = grps[i],
        region           = regs[i],
        metric_type      = mets[i]
      ))
    }

    ct <- suppressWarnings(
      stats::cor.test(x[ok], y[ok], method = method, exact = FALSE)
    )
    tibble::tibble(
      feature          = co,
      target_col       = target_col,
      spearman_rho     = unname(ct$estimate),
      p_value          = ct$p.value,
      n_complete       = n,
      abs_spearman_rho = abs(unname(ct$estimate)),
      feature_group    = grps[i],
      region           = regs[i],
      metric_type      = mets[i]
    )
  })
}


#' Compute the pairwise Spearman correlation matrix for a set of features.
#'
#' Uses vectorised stats::cor() for speed. NA values are handled with
#' pairwise.complete.obs.
#'
#' @param data         Dataframe from build_dataset().
#' @param feature_cols Character vector of feature column names.
#' @param method       Correlation method (default "spearman").
#' @return Named square numeric matrix (p × p) with row/col names = feature
#'   names. NaN entries (constant columns in pairwise subset) become NA.
#' @export
compute_feature_feature_correlations <- function(data,
                                                 feature_cols,
                                                 method = "spearman") {
  feature_cols <- intersect(feature_cols, names(data))
  feature_mat  <- as.matrix(data[, feature_cols, drop = FALSE])
  rho_mat      <- stats::cor(feature_mat, method = method,
                              use = "pairwise.complete.obs")
  rho_mat[is.nan(rho_mat)] <- NA_real_
  rho_mat
}


# =============================================================================
# Feature selection
# =============================================================================

#' Select top-N features by absolute Spearman correlation with the target.
#'
#' If both top_n and abs_rho_threshold are supplied, the threshold is applied
#' first, then the top-N cap is applied to the surviving features.
#' To use threshold only, pass top_n = NULL.
#' To use top-N only,   pass abs_rho_threshold = NULL.
#'
#' @param feature_target_tbl Tibble from compute_feature_target_correlations().
#' @param top_n              Integer or NULL. Number of top features to select.
#' @param abs_rho_threshold  Numeric or NULL. Minimum |rho| to retain.
#' @return Character vector of selected feature names, sorted descending |rho|.
#' @export
select_top_target_features <- function(feature_target_tbl,
                                       top_n             = 50L,
                                       abs_rho_threshold = NULL) {
  tbl <- feature_target_tbl |>
    dplyr::filter(!is.na(abs_spearman_rho)) |>
    dplyr::arrange(dplyr::desc(abs_spearman_rho))

  if (!is.null(abs_rho_threshold)) {
    tbl <- dplyr::filter(tbl, abs_spearman_rho >= abs_rho_threshold)
  }

  if (!is.null(top_n) && nrow(tbl) > as.integer(top_n)) {
    tbl <- dplyr::slice_head(tbl, n = as.integer(top_n))
  }

  tbl$feature
}


#' Select top-K features per biological group, balanced across groups.
#'
#' Prevents any one large feature family from dominating the feature set.
#' Features with NA feature_group are pooled under "other".
#'
#' @param feature_target_tbl    Tibble from compute_feature_target_correlations().
#' @param top_k_per_group       Max features per group (default 5L).
#' @param min_features_per_group Groups with fewer features than this are
#'                              excluded (default 1L).
#' @param max_features_total    Hard cap after group-balanced selection;
#'                              excess features dropped by descending |rho|
#'                              (default 100L, NULL = no cap).
#' @return Character vector of selected feature names.
#' @export
select_group_balanced_features <- function(feature_target_tbl,
                                           top_k_per_group        = 5L,
                                           min_features_per_group = 1L,
                                           max_features_total     = 100L) {
  tbl <- feature_target_tbl |>
    dplyr::filter(!is.na(abs_spearman_rho)) |>
    dplyr::mutate(feature_group = dplyr::coalesce(feature_group, "other"))

  selected <- tbl |>
    dplyr::group_by(feature_group) |>
    dplyr::arrange(dplyr::desc(abs_spearman_rho), .by_group = TRUE) |>
    dplyr::filter(dplyr::n() >= as.integer(min_features_per_group)) |>
    dplyr::slice_head(n = as.integer(top_k_per_group)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(abs_spearman_rho))

  feats <- selected$feature

  if (!is.null(max_features_total) && length(feats) > as.integer(max_features_total)) {
    feats <- feats[seq_len(as.integer(max_features_total))]
  }

  feats
}


#' Cluster features by pairwise |rho|, using hierarchical clustering.
#'
#' Distance metric: 1 - |rho|. Linkage: average. NA rho values are treated as
#' zero correlation (features assumed uncorrelated where data is sparse).
#'
#' @param feature_cor_mat Square numeric matrix from
#'                        compute_feature_feature_correlations().
#' @param abs_rho_cutoff  |rho| threshold for cluster membership. Features
#'                        with |rho| >= cutoff are considered redundant and
#'                        may end up in the same cluster (distance cutoff =
#'                        1 - abs_rho_cutoff; default 0.85).
#' @return Tibble: feature, cluster_id, cluster_size.
#' @export
cluster_features_by_correlation <- function(feature_cor_mat,
                                             abs_rho_cutoff = 0.85) {
  n_f  <- nrow(feature_cor_mat)
  nms  <- rownames(feature_cor_mat)

  if (n_f < 2L) {
    return(tibble::tibble(feature = nms, cluster_id = 1L, cluster_size = 1L))
  }

  mat           <- feature_cor_mat
  mat[is.na(mat)] <- 0
  diag(mat)       <- 1

  dist_mat <- stats::as.dist(1 - abs(mat))
  hc       <- stats::hclust(dist_mat, method = "average")
  clusters <- stats::cutree(hc, h = 1 - abs_rho_cutoff)

  tibble::tibble(feature = names(clusters), cluster_id = as.integer(clusters)) |>
    dplyr::group_by(cluster_id) |>
    dplyr::mutate(cluster_size = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::arrange(cluster_id, feature)
}


#' Select one representative feature per cluster.
#'
#' Selection rule (in order):
#'   1. Highest absolute Spearman rho with the target (primary).
#'   2. Biological group priority (.CW_REP_PRIORITY) — lower value = preferred.
#'      Biologically important corrected structure features rank first.
#'   3. Alphabetical feature name (deterministic tie-breaker).
#'
#' @param cluster_tbl        Tibble from cluster_features_by_correlation().
#' @param feature_target_tbl Tibble from compute_feature_target_correlations().
#' @return Tibble with cluster_tbl columns plus: representative_feature,
#'   is_representative, spearman_rho_with_target, abs_spearman_rho_with_target,
#'   feature_group, region, metric_type, cluster_size.
#' @export
select_cluster_representatives <- function(cluster_tbl, feature_target_tbl) {
  tbl <- cluster_tbl |>
    dplyr::left_join(
      dplyr::select(feature_target_tbl,
                    feature,
                    spearman_rho_with_target     = spearman_rho,
                    abs_spearman_rho_with_target = abs_spearman_rho,
                    feature_group, region, metric_type),
      by = "feature"
    ) |>
    dplyr::mutate(
      .group_pri     = dplyr::coalesce(
        as.integer(.CW_REP_PRIORITY[feature_group]), 999L
      ),
      .rho_for_sort  = dplyr::coalesce(abs_spearman_rho_with_target, -1)
    )

  reps <- tbl |>
    dplyr::group_by(cluster_id) |>
    dplyr::arrange(dplyr::desc(.rho_for_sort), .group_pri, feature,
                   .by_group = TRUE) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::select(cluster_id, representative_feature = feature)

  tbl |>
    dplyr::left_join(reps, by = "cluster_id") |>
    dplyr::mutate(is_representative = feature == representative_feature) |>
    dplyr::select(-.group_pri, -.rho_for_sort) |>
    dplyr::arrange(cluster_id, dplyr::desc(is_representative), feature)
}


# =============================================================================
# Heatmap plotting
# =============================================================================

#' Plot a feature × feature Spearman correlation heatmap.
#'
#' Visual style matches region_feature_heatmap.R: diverging blue-white-red
#' fill (RdBu colours), square tiles, white cell borders, in-cell rho labels,
#' group separator lines, heavy response separator.
#'
#' For all-region heatmaps, axis labels use format_col_name() (includes region
#' suffix) rather than format_metric_name(), because features from different
#' regions appear on the same axes.
#'
#' @param feature_cor_mat   Named square numeric matrix from
#'                          compute_feature_feature_correlations(). Row/column
#'                          names are feature column names.
#' @param feature_metadata  Tibble with columns: feature, feature_group.
#'                          Optional columns: region, metric_type.
#' @param target_cor_tbl    Optional tibble from
#'                          compute_feature_target_correlations(). When
#'                          non-NULL, a target row/column is appended to the
#'                          matrix separated by a heavy line.
#' @param target_col        Response column name for display (default
#'                          "halflife").
#' @param cluster           Logical. Reorder features by hierarchical
#'                          clustering on (1 - abs(rho)) / 2 distance,
#'                          average linkage (default TRUE).
#' @param triangle          "full" (default), "lower", or "upper".
#' @param color_limits      Numeric length-2 for fill scale (default c(-1, 1)).
#' @param cell_text_size    geom_text size for in-cell rho labels (default 3.2).
#' @param label_threshold   Only label cells where |rho| >= this value
#'                          (default 0 = label all cells).
#' @param base_size         Base font size in points (default 14).
#' @param title             Plot title (NULL = auto-generated).
#' @param subtitle          Plot subtitle (NULL = auto-generated).
#' @param output_path       If non-NULL, save PDF to this path.
#' @param size_mm           Square side length in mm for saved PDF.
#'                          NULL = auto-scaled from number of features.
#' @return list(plot = ggplot, table = tibble). Table columns: feature_x,
#'   feature_y, rho, feature_group_1, feature_group_2, region_1, region_2.
#'   Returned invisibly.
#' @export
plot_feature_correlation_heatmap <- function(feature_cor_mat,
                                              feature_metadata,
                                              target_cor_tbl   = NULL,
                                              target_col       = "halflife",
                                              cluster          = TRUE,
                                              triangle         = c("full", "lower", "upper"),
                                              color_limits     = c(-1, 1),
                                              cell_text_size   = 3.2,
                                              label_threshold  = 0,
                                              base_size        = 14,
                                              title            = NULL,
                                              subtitle         = NULL,
                                              output_path      = NULL,
                                              size_mm          = NULL) {

  triangle      <- match.arg(triangle)
  feature_names <- rownames(feature_cor_mat)
  n_features    <- length(feature_names)

  if (n_features < 2L) {
    warning("plot_feature_correlation_heatmap: fewer than 2 features — skipping")
    return(invisible(NULL))
  }

  # --- Build metadata lookup maps (R5: guard with %in% names) -----------------
  meta_grp <- if ("feature_group" %in% names(feature_metadata)) {
    setNames(as.character(feature_metadata$feature_group),
             feature_metadata$feature)
  } else {
    setNames(rep(NA_character_, n_features), feature_names)
  }
  meta_reg <- if ("region" %in% names(feature_metadata)) {
    setNames(as.character(feature_metadata$region),
             feature_metadata$feature)
  } else {
    setNames(rep(NA_character_, n_features), feature_names)
  }

  # --- Optionally augment matrix with target row/column -----------------------
  has_target <- !is.null(target_cor_tbl)

  if (has_target) {
    tbl_sub <- dplyr::filter(target_cor_tbl, feature %in% feature_names)
    t_rho   <- setNames(tbl_sub$spearman_rho, tbl_sub$feature)
    t_vals  <- vapply(feature_names,
                      function(f) if (f %in% names(t_rho)) t_rho[[f]] else NA_real_,
                      numeric(1))

    n_aug   <- n_features + 1L
    aug_mat <- matrix(NA_real_, nrow = n_aug, ncol = n_aug)
    aug_nms <- c(target_col, feature_names)
    dimnames(aug_mat) <- list(aug_nms, aug_nms)
    diag(aug_mat)                            <- 1
    aug_mat[target_col, feature_names]       <- t_vals
    aug_mat[feature_names, target_col]       <- t_vals
    aug_mat[feature_names, feature_names]    <- feature_cor_mat
    working_mat <- aug_mat
  } else {
    working_mat <- feature_cor_mat
  }

  # --- Clustering / ordering --------------------------------------------------
  feat_only <- feature_names

  if (cluster && length(feat_only) >= 2L) {
    sub_mat              <- working_mat[feat_only, feat_only, drop = FALSE]
    sub_mat[is.na(sub_mat)] <- 0
    diag(sub_mat)        <- 1
    dist_mat             <- stats::as.dist((1 - sub_mat) / 2)
    hc                   <- stats::hclust(dist_mat, method = "average")
    feat_order           <- feat_only[hc$order]
  } else {
    feat_order <- feat_only
  }

  all_names   <- if (has_target) c(target_col, feat_order) else feat_order
  working_mat <- working_mat[all_names, all_names, drop = FALSE]
  n_all       <- length(all_names)

  # --- Build tidy long-form data ----------------------------------------------
  idx <- which(upper.tri(working_mat), arr.ind = TRUE)
  pairs_df <- tibble::tibble(
    feature_x = all_names[idx[, 1]],
    feature_y = all_names[idx[, 2]],
    rho       = working_mat[idx],
    is_diag   = FALSE
  )

  sym_df <- dplyr::bind_rows(
    pairs_df,
    dplyr::rename(pairs_df, feature_x = feature_y, feature_y = feature_x),
    tibble::tibble(
      feature_x = all_names,
      feature_y = all_names,
      rho       = NA_real_,
      is_diag   = TRUE
    )
  )

  sym_df <- sym_df |>
    dplyr::mutate(
      cell_label = dplyr::if_else(
        is_diag | is.na(rho) | abs(rho) < label_threshold, "",
        sprintf("%.2f", rho)
      ),
      text_color = dplyr::case_when(
        is_diag | is.na(rho) ~ "grey55",
        abs(rho) > 0.62      ~ "white",
        TRUE                 ~ "grey10"
      )
    )

  sym_df$feature_x <- factor(sym_df$feature_x, levels = all_names)
  sym_df$feature_y <- factor(sym_df$feature_y, levels = all_names)

  # Triangle masking
  sym_df <- sym_df |>
    dplyr::mutate(
      .xp      = as.integer(feature_x),
      .yp      = as.integer(feature_y),
      is_shown = is_diag | dplyr::case_when(
        triangle == "lower" ~ .xp >= .yp,
        triangle == "upper" ~ .xp <= .yp,
        TRUE                ~ TRUE
      )
    ) |>
    dplyr::select(-.xp, -.yp)

  # --- Group separator geometry -----------------------------------------------
  feat_grp_vec <- vapply(all_names, function(f) {
    if (has_target && identical(f, target_col)) "response"
    else {
      g <- meta_grp[f]
      if (is.null(g) || is.na(g)) "other" else g
    }
  }, character(1))

  rle_grp      <- rle(feat_grp_vec)
  grp_ends     <- cumsum(rle_grp$lengths)
  boundary_pos <- if (length(grp_ends) > 1L) grp_ends[-length(grp_ends)] + 0.5 else numeric(0)

  is_resp_idx  <- which(
    rle_grp$values == "response" |
    c(rle_grp$values[-1L] == "response", FALSE)
  )
  resp_boundary <- if (has_target && length(is_resp_idx) > 0)
    boundary_pos[is_resp_idx] else numeric(0)
  feat_boundary <- if (cluster) numeric(0) else setdiff(boundary_pos, resp_boundary)

  # --- Axis labels (full name including region for all-region heatmaps) -------
  axis_labels <- setNames(
    vapply(all_names, function(f) format_col_name(f), character(1)),
    all_names
  )

  # --- PDF dimensions ---------------------------------------------------------
  tile_mm <- max(12, min(22, 280 / n_all))
  pdf_sq  <- if (!is.null(size_mm)) size_mm else n_all * tile_mm + 80

  # --- Titles -----------------------------------------------------------------
  n_feat_disp <- if (has_target) n_all - 1L else n_all
  auto_title  <- sprintf("Feature × feature correlation — %d features", n_feat_disp)
  auto_sub    <- sprintf(
    "Spearman rho%s%s%s",
    if (has_target) paste0("; + ", format_col_name(target_col)) else "",
    if (cluster) "; clustered by |rho|" else "",
    if (label_threshold > 0) sprintf("; labels |rho| ≥ %.2f", label_threshold) else ""
  )
  plot_title <- if (!is.null(title))    title    else auto_title
  plot_sub   <- if (!is.null(subtitle)) subtitle else auto_sub

  # --- Build ggplot2 ----------------------------------------------------------
  shown_df <- dplyr::filter(sym_df, is_shown & !is_diag)

  p <- ggplot2::ggplot(sym_df,
                       ggplot2::aes(x = feature_x, y = feature_y)) +

    ggplot2::geom_tile(
      data = shown_df, ggplot2::aes(fill = rho),
      colour = "white", linewidth = 0.25, na.rm = TRUE
    ) +

    ggplot2::geom_tile(
      data = dplyr::filter(sym_df, is_diag),
      fill = "grey88", colour = "white", linewidth = 0.25
    ) +

    ggplot2::geom_text(
      data = shown_df,
      ggplot2::aes(label = cell_label, colour = text_color),
      size = cell_text_size, fontface = "bold", na.rm = TRUE
    ) +

    {if (length(feat_boundary) > 0)
      list(
        ggplot2::geom_hline(yintercept = feat_boundary,
                            colour = "grey40", linewidth = 0.7),
        ggplot2::geom_vline(xintercept = feat_boundary,
                            colour = "grey40", linewidth = 0.7)
      )
    } +

    {if (length(resp_boundary) > 0)
      list(
        ggplot2::geom_hline(yintercept = resp_boundary,
                            colour = "grey10", linewidth = 1.4),
        ggplot2::geom_vline(xintercept = resp_boundary,
                            colour = "grey10", linewidth = 1.4)
      )
    } +

    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = color_limits, na.value = "grey88",
      name = "Spearman rho",
      guide = ggplot2::guide_colorbar(
        title.position  = "top",
        barwidth        = ggplot2::unit(0.8, "cm"),
        barheight       = ggplot2::unit(5.5, "cm"),
        ticks.linewidth = 0.5
      )
    ) +

    ggplot2::scale_colour_identity(guide = "none") +

    ggplot2::scale_x_discrete(
      labels = axis_labels,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_discrete(
      labels = axis_labels,
      expand = ggplot2::expansion(add = 0.5)
    ) +

    ggplot2::labs(title = plot_title, subtitle = plot_sub, x = NULL, y = NULL) +

    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      aspect.ratio   = 1,
      axis.text.x    = ggplot2::element_text(
        angle = 45, vjust = 1, hjust = 1, size = base_size - 2
      ),
      axis.text.y    = ggplot2::element_text(size = base_size - 2),
      plot.title     = ggplot2::element_text(size = base_size + 2, face = "bold"),
      plot.subtitle  = ggplot2::element_text(size = base_size - 2, colour = "grey30"),
      panel.grid     = ggplot2::element_blank(),
      panel.border   = ggplot2::element_rect(colour = "grey20", fill = NA, linewidth = 0.6),
      legend.position = "right",
      legend.title   = ggplot2::element_text(size = base_size - 2, face = "bold"),
      plot.margin    = ggplot2::margin(t = 8, r = 15, b = 8, l = 8, unit = "mm")
    )

  # --- Save PDF ---------------------------------------------------------------
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(
      output_path, plot = p,
      width = pdf_sq, height = pdf_sq,
      units = "mm", device = "pdf"
    )
    message("  Saved: ", output_path)
  }

  # --- Return table -----------------------------------------------------------
  # Unique pairs only (upper triangle, no diagonal); metadata attached.
  table_out <- sym_df |>
    dplyr::filter(!is_diag) |>
    dplyr::filter(as.integer(feature_x) < as.integer(feature_y)) |>
    dplyr::mutate(
      feature_group_1 = unname(meta_grp[as.character(feature_x)]),
      feature_group_2 = unname(meta_grp[as.character(feature_y)]),
      region_1        = unname(meta_reg[as.character(feature_x)]),
      region_2        = unname(meta_reg[as.character(feature_y)])
    ) |>
    dplyr::select(feature_x, feature_y, rho,
                  feature_group_1, feature_group_2, region_1, region_2)

  invisible(list(plot = p, table = table_out))
}


# =============================================================================
# Output helpers
# =============================================================================

#' Convert a feature-feature correlation matrix to a filtered long-form table.
#'
#' @param ff_mat            Square numeric matrix.
#' @param feat_meta         Tibble: feature, feature_group, region.
#' @param abs_rho_threshold Only retain pairs with |rho| >= threshold (default
#'                          0.7). NULL = keep all pairs.
#' @return Tibble: feature_1, feature_2, spearman_rho, abs_spearman_rho,
#'   feature_group_1, feature_group_2, region_1, region_2.
.cw_ff_mat_to_long <- function(ff_mat, feat_meta, abs_rho_threshold = 0.7) {
  nms <- rownames(ff_mat)
  idx <- which(upper.tri(ff_mat), arr.ind = TRUE)

  tbl <- tibble::tibble(
    feature_1        = nms[idx[, 1]],
    feature_2        = nms[idx[, 2]],
    spearman_rho     = ff_mat[idx],
    abs_spearman_rho = abs(ff_mat[idx])
  )

  meta_1 <- dplyr::select(feat_meta, feature,
                           feature_group_1 = feature_group, region_1 = region)
  meta_2 <- dplyr::select(feat_meta, feature,
                           feature_group_2 = feature_group, region_2 = region)

  tbl <- tbl |>
    dplyr::left_join(meta_1, by = c("feature_1" = "feature")) |>
    dplyr::left_join(meta_2, by = c("feature_2" = "feature"))

  if (!is.null(abs_rho_threshold)) {
    tbl <- dplyr::filter(tbl, !is.nan(abs_spearman_rho) &
                                abs_spearman_rho >= abs_rho_threshold)
  }

  dplyr::arrange(tbl, dplyr::desc(abs_spearman_rho))
}


#' Write correlation workflow tables to disk (R8: under OUTPUT_DIR).
#'
#' @param feature_target_tbl From compute_feature_target_correlations().
#' @param ff_cor_long        Long-form feature-feature table.
#' @param cluster_rep_tbl    From select_cluster_representatives().
#' @param table_dir          Directory for outputs.
#' @param file_prefix        Filename prefix.
#' @param species_label      Species label appended to filenames (or NULL).
#' @return Named list of written file paths (invisibly).
#' @export
write_correlation_outputs <- function(feature_target_tbl = NULL,
                                      ff_cor_long        = NULL,
                                      cluster_rep_tbl    = NULL,
                                      table_dir,
                                      file_prefix    = "correlation_workflow",
                                      species_label  = NULL) {
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  sp_suf <- if (!is.null(species_label)) paste0("_", species_label) else ""
  paths  <- list()

  .write_csv <- function(tbl, name) {
    p <- file.path(table_dir, sprintf("%s_%s%s.csv", file_prefix, name, sp_suf))
    write.csv(tbl, p, row.names = FALSE)
    message("  Wrote: ", p)
    p
  }

  if (!is.null(feature_target_tbl))
    paths$feature_target  <- .write_csv(feature_target_tbl, "feature_target_correlations")
  if (!is.null(ff_cor_long))
    paths$feature_feature <- .write_csv(ff_cor_long, "feature_feature_correlations")
  if (!is.null(cluster_rep_tbl))
    paths$cluster_reps    <- .write_csv(cluster_rep_tbl, "cluster_representatives")

  invisible(paths)
}


# =============================================================================
# Main workflow entry point
# =============================================================================

#' Run the complete feature-correlation heatmap diagnostic workflow.
#'
#' Generates Plot B (target-filtered), C (cluster-representative), and D
#' (group-balanced) heatmaps alongside supporting tables and a manifest.
#' Processes each species in df separately.
#'
#' @param df                           Dataframe from build_dataset() or
#'                                     build_all(). Must have a species column.
#' @param species                      Character vector of species to process,
#'                                     or NULL (default) to use all unique
#'                                     values in df$species.
#' @param groups                       Feature groups/bundles to include
#'                                     (default INCLUDED_GROUPS from config.R).
#' @param target_col                   Response column name (default
#'                                     "halflife").
#' @param correlation_method           Correlation method (default "spearman").
#' @param top_n_target_features        Integer vector of top-N sizes for Plot B
#'                                     (default c(30L, 50L, 100L)).
#' @param target_abs_rho_threshold     Minimum |rho| for Plot B, or NULL
#'                                     (default NULL; see select_top_target_features()).
#' @param cluster_abs_rho_cutoff       |rho| cutoff for Plot C clustering
#'                                     (default 0.85).
#' @param top_k_per_group              Features per group for Plot D
#'                                     (default 5L).
#' @param min_features_per_group       Min group size for Plot D (default 1L).
#' @param max_features_total           Hard cap for Plot D total features
#'                                     (default 100L).
#' @param feature_pair_abs_rho_threshold |rho| threshold for the feature-feature
#'                                     long table output (default 0.7). NULL =
#'                                     write all pairs (may be large).
#' @param min_complete_n               Minimum non-NA observations per feature
#'                                     column (default 100L).
#' @param drop_zero_variance           Drop constant feature columns (default
#'                                     TRUE).
#' @param output_dir                   Directory for plot PDFs (default
#'                                     OUTPUT_DIR/plots/heatmaps/).
#' @param table_dir                    Directory for CSV tables (default
#'                                     OUTPUT_DIR/tables/).
#' @param file_prefix                  Filename prefix (default
#'                                     "correlation_workflow").
#' @param label_threshold              Suppress in-cell labels below this |rho|
#'                                     (default 0.3).
#' @param cluster_plots                Apply hierarchical clustering to feature
#'                                     order in all plots (default TRUE).
#' @return Named list of all outputs (plots and tables), invisibly.
#' @export
run_correlation_heatmap_workflow <- function(df,
                                             species                        = NULL,
                                             groups                         = INCLUDED_GROUPS,
                                             target_col                     = "halflife",
                                             correlation_method             = "spearman",
                                             top_n_target_features          = c(30L, 50L, 100L),
                                             target_abs_rho_threshold       = NULL,
                                             cluster_abs_rho_cutoff         = 0.85,
                                             top_k_per_group                = 5L,
                                             min_features_per_group         = 1L,
                                             max_features_total             = 100L,
                                             feature_pair_abs_rho_threshold = 0.7,
                                             min_complete_n                 = 100L,
                                             drop_zero_variance             = TRUE,
                                             output_dir     = file.path(OUTPUT_DIR, "plots", "heatmaps"),
                                             table_dir      = file.path(OUTPUT_DIR, "tables"),
                                             file_prefix    = "correlation_workflow",
                                             label_threshold = 0.3,
                                             cluster_plots   = TRUE) {

  # --- Guards (R5, R6) --------------------------------------------------------
  if (!"species" %in% names(df)) stop("species column missing — pipeline invariant (R6)")
  if (!target_col %in% names(df)) stop("target_col '", target_col, "' not in df")

  species_vec <- if (!is.null(species)) species else unique(df$species)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)

  all_outputs   <- list()
  manifest_rows <- list()

  for (sp in species_vec) {
    message("\n=== Species: ", sp, " ===")

    df_sp <- if (length(unique(df$species)) > 1L) {
      dplyr::filter(df, species == sp)
    } else {
      df
    }
    df_sp <- df_sp[!is.na(df_sp[[target_col]]), , drop = FALSE]

    # ---- 1. Feature selection (R3: use select_features / fg) -----------------
    message("Selecting features from groups: ",
            paste(groups, collapse = ", "))
    candidate_features <- select_features(df_sp, groups = groups)

    feature_cols <- get_numeric_feature_cols(
      df_sp,
      target_col         = target_col,
      drop_zero_variance = drop_zero_variance,
      min_complete_n     = min_complete_n
    )
    feature_cols <- intersect(feature_cols, candidate_features)
    n_feats      <- length(feature_cols)
    message("  ", n_feats, " numeric feature columns after filtering")

    if (n_feats < 2L) {
      warning("Fewer than 2 features for '", sp, "' — skipping species")
      next
    }

    # ---- 2. Feature metadata -------------------------------------------------
    feat_meta <- tibble::tibble(
      feature       = feature_cols,
      feature_group = assign_feature_group(feature_cols),
      region        = parse_feature_region(feature_cols),
      metric_type   = parse_metric_type(feature_cols)
    )

    # ---- 3. Feature-target correlations --------------------------------------
    message("Computing feature-target (", target_col, ") correlations...")
    ft_tbl <- compute_feature_target_correlations(
      df_sp, feature_cols, target_col, method = correlation_method
    )

    # ---- 4. Feature-feature correlation matrix (full set) -------------------
    message("Computing feature-feature correlation matrix (",
            n_feats, " x ", n_feats, ")...")
    ff_mat <- compute_feature_feature_correlations(
      df_sp, feature_cols, method = correlation_method
    )

    # ---- 5. Write tables -----------------------------------------------------
    ff_long <- .cw_ff_mat_to_long(ff_mat, feat_meta,
                                   abs_rho_threshold = feature_pair_abs_rho_threshold)
    write_correlation_outputs(
      feature_target_tbl = ft_tbl,
      ff_cor_long        = ff_long,
      table_dir          = table_dir,
      file_prefix        = file_prefix,
      species_label      = sp
    )

    all_outputs[[paste0("feature_target_tbl_", sp)]] <- ft_tbl
    all_outputs[[paste0("ff_cor_long_", sp)]]         <- ff_long

    # ---- 6. Plot B: target-filtered heatmaps --------------------------------
    message("\n--- Plot B: target-filtered heatmaps ---")

    for (top_n in top_n_target_features) {
      sel_feats <- select_top_target_features(
        ft_tbl,
        top_n             = top_n,
        abs_rho_threshold = target_abs_rho_threshold
      )
      n_sel <- length(sel_feats)

      if (n_sel < 2L) {
        message("  top_n = ", top_n, ": fewer than 2 features selected — skipping")
        next
      }
      message("  top_n = ", top_n, ": ", n_sel, " features")

      sub_mat  <- ff_mat[sel_feats, sel_feats, drop = FALSE]
      sub_meta <- dplyr::filter(feat_meta, feature %in% sel_feats)

      pfx   <- sprintf("%s_B_top%d_%s", file_prefix, top_n, sp)
      fpath <- file.path(output_dir, paste0(pfx, ".pdf"))

      res <- plot_feature_correlation_heatmap(
        feature_cor_mat  = sub_mat,
        feature_metadata = sub_meta,
        target_cor_tbl   = ft_tbl,
        target_col       = target_col,
        cluster          = cluster_plots,
        label_threshold  = label_threshold,
        title            = sprintf(
          "Target-filtered feature correlation — %s [%s]",
          format_col_name(target_col), sp),
        subtitle         = sprintf(
          "Top %d features by |rho| with %s%s; Spearman rho",
          n_sel, format_col_name(target_col),
          if (!is.null(target_abs_rho_threshold))
            sprintf("; |rho| >= %.2f", target_abs_rho_threshold) else ""),
        output_path      = fpath
      )

      key <- paste0("B_top", top_n, "_", sp)
      all_outputs[[key]] <- res

      if (!is.null(res)) {
        tab_path <- file.path(table_dir, paste0(pfx, "_pairs.csv"))
        write.csv(res$table, tab_path, row.names = FALSE)
        manifest_rows[[length(manifest_rows) + 1L]] <- tibble::tibble(
          plot_name        = basename(fpath),
          plot_type        = "B_target_filtered",
          n_features       = n_sel,
          selection_method = sprintf("top_%d_by_abs_rho_with_%s", top_n, target_col),
          output_path      = fpath,
          table_path       = tab_path,
          parameters       = sprintf(
            "top_n=%d; threshold=%s; species=%s; cluster=%s",
            top_n,
            if (is.null(target_abs_rho_threshold)) "NULL"
            else as.character(target_abs_rho_threshold),
            sp, cluster_plots
          )
        )
      }
    }

    # ---- 7. Plot C: cluster-representative heatmap --------------------------
    message("\n--- Plot C: cluster-representative heatmap ---")

    clust_tbl <- cluster_features_by_correlation(ff_mat,
                                                  abs_rho_cutoff = cluster_abs_rho_cutoff)
    rep_tbl   <- select_cluster_representatives(clust_tbl, ft_tbl)

    n_clust   <- length(unique(rep_tbl$cluster_id))
    rep_feats <- rep_tbl$feature[rep_tbl$is_representative]
    n_reps    <- length(rep_feats)
    message("  ", n_feats, " features -> ", n_clust, " clusters -> ", n_reps,
            " representatives (cutoff |rho| = ", cluster_abs_rho_cutoff, ")")

    # Write cluster table
    write_correlation_outputs(
      cluster_rep_tbl = rep_tbl,
      table_dir       = table_dir,
      file_prefix     = file_prefix,
      species_label   = sp
    )
    all_outputs[[paste0("cluster_rep_tbl_", sp)]] <- rep_tbl

    if (n_reps >= 2L) {
      sub_mat  <- ff_mat[rep_feats, rep_feats, drop = FALSE]
      sub_meta <- dplyr::filter(feat_meta, feature %in% rep_feats)

      pfx   <- sprintf("%s_C_cluster_reps_%s", file_prefix, sp)
      fpath <- file.path(output_dir, paste0(pfx, ".pdf"))

      res <- plot_feature_correlation_heatmap(
        feature_cor_mat  = sub_mat,
        feature_metadata = sub_meta,
        target_cor_tbl   = ft_tbl,
        target_col       = target_col,
        cluster          = cluster_plots,
        label_threshold  = label_threshold,
        title            = sprintf(
          "Cluster-representative feature correlation — %s [%s]",
          format_col_name(target_col), sp),
        subtitle         = sprintf(
          "%d clusters (|rho| cutoff = %.2f) -> %d representatives; Spearman rho",
          n_clust, cluster_abs_rho_cutoff, n_reps),
        output_path      = fpath
      )

      key <- paste0("C_cluster_reps_", sp)
      all_outputs[[key]] <- res

      if (!is.null(res)) {
        tab_path <- file.path(table_dir, paste0(pfx, "_pairs.csv"))
        write.csv(res$table, tab_path, row.names = FALSE)
        manifest_rows[[length(manifest_rows) + 1L]] <- tibble::tibble(
          plot_name        = basename(fpath),
          plot_type        = "C_cluster_representative",
          n_features       = n_reps,
          selection_method = sprintf("cluster_rep_cutoff_%.2f", cluster_abs_rho_cutoff),
          output_path      = fpath,
          table_path       = tab_path,
          parameters       = sprintf(
            "cutoff=%.2f; n_clusters=%d; species=%s; cluster_plots=%s",
            cluster_abs_rho_cutoff, n_clust, sp, cluster_plots
          )
        )
      }
    } else {
      message("  Fewer than 2 cluster representatives — skipping Plot C")
    }

    # ---- 8. Plot D: group-balanced heatmap -----------------------------------
    message("\n--- Plot D: group-balanced heatmap ---")

    gbal_feats <- select_group_balanced_features(
      ft_tbl,
      top_k_per_group        = top_k_per_group,
      min_features_per_group = min_features_per_group,
      max_features_total     = max_features_total
    )
    n_gbal <- length(gbal_feats)
    message("  ", n_gbal, " features selected (top ", top_k_per_group, " per group)")

    if (n_gbal >= 2L) {
      sub_mat  <- ff_mat[gbal_feats, gbal_feats, drop = FALSE]
      sub_meta <- dplyr::filter(feat_meta, feature %in% gbal_feats)

      pfx   <- sprintf("%s_D_group_balanced_%s", file_prefix, sp)
      fpath <- file.path(output_dir, paste0(pfx, ".pdf"))

      res <- plot_feature_correlation_heatmap(
        feature_cor_mat  = sub_mat,
        feature_metadata = sub_meta,
        target_cor_tbl   = ft_tbl,
        target_col       = target_col,
        cluster          = cluster_plots,
        label_threshold  = label_threshold,
        title            = sprintf(
          "Group-balanced feature correlation — %s [%s]",
          format_col_name(target_col), sp),
        subtitle         = sprintf(
          "Top %d per biological group (max %d total); Spearman rho",
          top_k_per_group, max_features_total),
        output_path      = fpath
      )

      key <- paste0("D_group_balanced_", sp)
      all_outputs[[key]] <- res

      if (!is.null(res)) {
        tab_path <- file.path(table_dir, paste0(pfx, "_pairs.csv"))
        write.csv(res$table, tab_path, row.names = FALSE)
        manifest_rows[[length(manifest_rows) + 1L]] <- tibble::tibble(
          plot_name        = basename(fpath),
          plot_type        = "D_group_balanced",
          n_features       = n_gbal,
          selection_method = sprintf("top_%d_per_group", top_k_per_group),
          output_path      = fpath,
          table_path       = tab_path,
          parameters       = sprintf(
            "top_k=%d; min_per_group=%d; max_total=%d; species=%s; cluster_plots=%s",
            top_k_per_group, min_features_per_group,
            max_features_total, sp, cluster_plots
          )
        )
      }
    } else {
      message("  Fewer than 2 group-balanced features — skipping Plot D")
    }
  }  # end species loop

  # ---- 9. Write manifest ----------------------------------------------------
  if (length(manifest_rows) > 0L) {
    manifest <- dplyr::bind_rows(manifest_rows)
    man_path <- file.path(table_dir, sprintf("%s_manifest.csv", file_prefix))
    write.csv(manifest, man_path, row.names = FALSE)
    message("\nManifest (", nrow(manifest), " plots): ", man_path)
    all_outputs$manifest <- manifest
  }

  message("\nWorkflow complete. Plots: ", output_dir)
  invisible(all_outputs)
}


# =============================================================================
# Runner (executed when sourced directly, not via load_all.R)
# =============================================================================

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  df <- build_dataset("human")

  out <- run_correlation_heatmap_workflow(
    df,
    groups                         = INCLUDED_GROUPS,
    target_col                     = "halflife",
    correlation_method             = "spearman",
    top_n_target_features          = c(30L, 50L, 100L),
    target_abs_rho_threshold       = NULL,
    cluster_abs_rho_cutoff         = 0.85,
    top_k_per_group                = 5L,
    min_features_per_group         = 1L,
    max_features_total             = 100L,
    feature_pair_abs_rho_threshold = 0.7,
    min_complete_n                 = 100L,
    drop_zero_variance             = TRUE,
    label_threshold                = 0.3,
    cluster_plots                  = TRUE,
    output_dir  = file.path(OUTPUT_DIR, "plots", "heatmaps"),
    table_dir   = file.path(OUTPUT_DIR, "tables"),
    file_prefix = "correlation_workflow"
  )

  message("Done.")
}
