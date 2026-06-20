# =============================================================================
# Region feature heatmap: pairwise Spearman correlation matrix per region
# =============================================================================
# Builds a full feature × feature Spearman correlation matrix for each requested
# region. The response variable (default: halflife) is included as a column so
# you can read off every feature's relationship with the predictor, while also
# seeing how features co-vary with each other.
#
# Layout:
#   * One PDF per region.
#   * Feature order: user-supplied vector first (in order), then any remaining
#     features appended at the end sorted by |ρ with response| descending.
#   * Response column sits at the bottom-right, separated by a heavier line.
#   * Diverging blue (−1) → white (0) → red (+1) fill; limits fixed at ±1.
#   * Spearman ρ printed in every off-diagonal cell (2 decimal places).
#   * Square tiles via theme(aspect.ratio = 1).
#   * Region suffix stripped from axis labels (format_metric_name()).
#   * Group separator lines mark boundaries where the group changes.
#
# top_n mode:
#   Pass top_n = list(structure = 3, intrinsic = 2) (or a single integer to
#   apply the same cap to every group) to produce a focused heatmap showing
#   only the top-N features per group (by |ρ with response|), all combined.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/region_feature_heatmap.R")
#   df  <- build_dataset("human")
#
#   # Structure features, CDS, custom order
#   out <- region_feature_heatmap(
#     df, groups = "structure", regions = "cds",
#     feature_order = c("rnafold_zscore_cds", "mfe_delta_cds", "rnalfold_score_cds")
#   )
#   print(out[["cds"]]$plot)
#
#   # Top-3 per group across core regions
#   out <- region_feature_heatmap(
#     df,
#     groups    = c("structure", "intrinsic"),
#     regions   = c("5utr", "cds", "3utr"),
#     top_n     = 3,
#     output_dir = "data/outputs/plots/heatmaps"
#   )
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
})


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Compute Spearman ρ and p-value for every unique pair in `cols`.
#' Returns a long tibble: feature_x, feature_y, rho, p_value, n.
.compute_corr_pairs <- function(df, cols, min_n = 30) {
  pairs <- utils::combn(cols, 2, simplify = FALSE)

  purrr::map_dfr(pairs, function(p) {
    x  <- df[[p[1]]]
    y  <- df[[p[2]]]
    ok <- !is.na(x) & !is.na(y)
    n  <- sum(ok)

    if (n < min_n ||
        length(unique(x[ok])) < 2 ||
        length(unique(y[ok])) < 2) {
      return(tibble::tibble(
        feature_x = p[1], feature_y = p[2],
        rho = NA_real_, p_value = NA_real_, n = n
      ))
    }

    ct <- suppressWarnings(
      stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    )
    tibble::tibble(
      feature_x = p[1], feature_y = p[2],
      rho       = unname(ct$estimate),
      p_value   = ct$p.value,
      n         = n
    )
  })
}


#' Absolute Spearman ρ of a single column with the response.
.abs_rho_with_response <- function(df, col, response) {
  x  <- df[[col]]
  y  <- df[[response]]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 5) return(0)
  abs(suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman")))
}


#' Apply top-N per-group filter to `reg_features`.
#'
#' @param reg_features Character vector of feature column names for this region.
#' @param col_to_group Named list: column → group key.
#' @param top_n        Single integer, or named list (group → integer).
#'                     Groups absent from the list are not filtered.
#' @param response     Response column name (used for |ρ| ranking).
#' @param df_sub       Data subset (used for |ρ| computation).
#' @return Filtered character vector preserving group structure.
.apply_top_n <- function(reg_features, col_to_group, top_n, response, df_sub) {
  if (is.null(top_n)) return(reg_features)

  # Normalise to a named list: group → n (NULL entry = no cap for that group)
  all_groups <- unique(vapply(reg_features, function(f) {
    g <- col_to_group[[f]]
    if (is.null(g)) "other" else g
  }, character(1)))

  if (is.numeric(top_n) && length(top_n) == 1 && is.null(names(top_n))) {
    n_map <- setNames(
      rep(list(as.integer(top_n)), length(all_groups)),
      all_groups
    )
  } else if (is.list(top_n) || !is.null(names(top_n))) {
    n_map <- as.list(top_n)
  } else {
    warning("top_n must be a single integer or a named list — ignored")
    return(reg_features)
  }

  kept <- character(0)
  for (g in all_groups) {
    g_features <- reg_features[vapply(reg_features, function(f) {
      identical(
        {grp <- col_to_group[[f]]; if (is.null(grp)) "other" else grp},
        g
      )
    }, logical(1))]

    n_keep <- n_map[[g]]
    if (is.null(n_keep) || length(g_features) <= n_keep) {
      kept <- c(kept, g_features)
      next
    }

    rho_vals <- vapply(g_features,
                       .abs_rho_with_response, numeric(1), df = df_sub, response = response)
    kept <- c(kept, g_features[order(-rho_vals)][seq_len(n_keep)])
  }

  kept
}


#' Order features for display in the heatmap.
#'
#' Response is placed first (bottom-left corner). When cluster = TRUE, features
#' are ordered by hierarchical clustering on their pairwise Spearman ρ matrix
#' (distance = (1 - ρ) / 2, average linkage); feature_order is ignored.
#' When cluster = FALSE, features listed in `feature_order` come first (in
#' order); any remaining features are appended sorted by |ρ with response|.
#'
#' @param reg_features   Character vector of feature columns for this region.
#' @param response       Response column name.
#' @param feature_order  Character vector of preferred column names (or NULL).
#' @param df_sub         Data subset for computing |ρ|.
#' @param cluster        Logical. If TRUE, order by hierarchical clustering.
#' @return Ordered character vector with response first.
.arrange_features <- function(reg_features, response, feature_order, df_sub,
                               cluster = FALSE) {
  if (cluster && length(reg_features) >= 2) {
    n_f     <- length(reg_features)
    rho_mat <- diag(1, n_f, n_f)
    dimnames(rho_mat) <- list(reg_features, reg_features)

    for (i in seq_len(n_f - 1L)) {
      for (j in seq(i + 1L, n_f)) {
        x  <- df_sub[[reg_features[i]]]
        y  <- df_sub[[reg_features[j]]]
        ok <- !is.na(x) & !is.na(y)
        r  <- if (sum(ok) >= 5 &&
                  length(unique(x[ok])) >= 2 &&
                  length(unique(y[ok])) >= 2) {
          suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
        } else {
          0   # treat uncomputable pairs as uncorrelated for clustering
        }
        rho_mat[i, j] <- r
        rho_mat[j, i] <- r
      }
    }
    rho_mat[is.na(rho_mat)] <- 0
    dist_mat <- stats::as.dist((1 - rho_mat) / 2)
    hc       <- stats::hclust(dist_mat, method = "average")
    return(c(response, reg_features[hc$order]))
  }

  # --- Original ordering (feature_order then |ρ| sort) ----------------------
  if (!is.null(feature_order)) {
    in_order  <- feature_order[feature_order %in% reg_features]
    remaining <- setdiff(reg_features, in_order)
  } else {
    in_order  <- character(0)
    remaining <- reg_features
  }

  if (length(remaining) > 0) {
    rho_vals <- vapply(remaining,
                       .abs_rho_with_response, numeric(1), df = df_sub, response = response)
    remaining <- remaining[order(-rho_vals)]
  }

  c(response, in_order, remaining)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Pairwise Spearman correlation heatmaps, one per region.
#'
#' Each heatmap shows the full feature × feature correlation matrix for the
#' columns belonging to one region, with the response variable included as a
#' special row/column. The response sits at the bottom-left corner.
#'
#' @param df               Dataframe from build_dataset() or build_all().
#' @param response         Response column to include in the matrix
#'                         (default "halflife").
#' @param groups           Groups / supergroups / bundles to include
#'                         (NULL = all groups in FEATURE_PATTERNS).
#' @param pick             Named list: group key → column allow-list.
#' @param drop             Named list: group key → columns to remove.
#' @param regions          Character vector of region tokens to plot
#'                         (NULL = all regions present in the resolved
#'                         features, in REGIONS canonical order).
#' @param feature_order    Character vector of column names specifying the
#'                         display order. These are placed first (in order);
#'                         any remaining selected features are appended at
#'                         the end sorted by |ρ with response| descending.
#'                         Ignored when cluster = TRUE.
#'                         NULL (default) = all features sorted by |ρ|.
#' @param top_n            Retain only the top-N features per group before
#'                         plotting. Can be a single integer (same cap for
#'                         every group) or a named list (group key → integer).
#'                         Groups absent from a named list are not filtered.
#'                         NULL (default) = no filtering.
#' @param cluster          Logical. If TRUE (default FALSE), order features by
#'                         hierarchical clustering on their pairwise Spearman ρ
#'                         matrix (distance = (1 - ρ) / 2, average linkage).
#'                         Groups are interleaved so group separator lines are
#'                         suppressed; only the response separator is drawn.
#' @param triangle         Which part of the symmetric matrix to display.
#'                         "full" (default) = complete matrix; "lower" = lower-
#'                         left triangle + diagonal (response appears as a row
#'                         at the bottom); "upper" = upper-right triangle +
#'                         diagonal. Excluded cells are left blank (panel
#'                         background shows through).
#' @param min_n            Minimum non-NA pairs needed to compute a correlation
#'                         (default 30).
#' @param color_limits     Numeric length-2: fill scale limits (default c(-1,1)).
#' @param cell_text_size   geom_text size for in-cell labels (default 3.2).
#' @param label_threshold  Numeric in [0, 1]. Only print the ρ value in a cell
#'                         when |ρ| >= this threshold. Cell fill colours are
#'                         always shown. 0 (default) = label every cell.
#' @param base_size        Base font size in points (default 14).
#' @param output_dir       If non-NULL, write one PDF per region here.
#' @param file_prefix      Filename prefix for saved PDFs
#'                         (default "region_heatmap").
#' @param size_mm          Square side length for saved PDFs in mm.
#'                         NULL = auto (scales with number of features).
#' @return Named list keyed by region (or "{region}_{species}" for multi-species
#'   datasets). Each value is list(plot = ggplot, table = tibble).
#'   Table columns: feature_x, feature_y, rho, p_value, q_value, n.
#'   Returned invisibly.
#' @export
region_feature_heatmap <- function(df,
                                   response       = "halflife",
                                   groups         = NULL,
                                   pick           = list(),
                                   drop           = list(),
                                   regions        = NULL,
                                   feature_order  = NULL,
                                   top_n          = NULL,
                                   cluster        = FALSE,
                                   triangle       = c("full", "lower", "upper"),
                                   min_n          = 30,
                                   color_limits   = c(-1, 1),
                                   cell_text_size  = 3.2,
                                   label_threshold = 0,
                                   base_size       = 14,
                                   output_dir     = NULL,
                                   file_prefix    = "region_heatmap",
                                   size_mm        = NULL) {

  triangle <- match.arg(triangle)

  # --- Guards ---------------------------------------------------------------
  if (!response %in% names(df)) {
    stop("response '", response, "' not in df")
  }
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  # --- Feature selection ----------------------------------------------------
  all_features <- select_features(df, groups, pick, drop)
  all_features <- setdiff(all_features, response)

  if (length(all_features) == 0) {
    stop("No feature columns after selection — check `groups`")
  }

  # Map each feature column → its group key and region token
  sel           <- resolve_selection(groups, pick, drop)
  col_to_group  <- list()
  col_to_region <- list()

  for (g in sel$groups) {
    cols <- refine_group_columns(fg_columns(df, g), sel$pick[[g]], sel$drop[[g]])
    for (co in intersect(cols, all_features)) {
      if (!is.null(col_to_group[[co]])) next
      tokens <- strsplit(co, "_", fixed = TRUE)[[1]]
      last   <- tokens[length(tokens)]
      if (last %in% REGIONS) {
        col_to_group[[co]]  <- g
        col_to_region[[co]] <- last
      }
    }
  }

  # Determine which regions to process
  available_regions <- unique(unlist(col_to_region))
  if (is.null(regions)) {
    regions_to_plot <- intersect(REGIONS, available_regions)
  } else {
    unknown <- setdiff(regions, REGIONS)
    if (length(unknown) > 0) {
      warning("Unknown region tokens ignored: ", paste(unknown, collapse = ", "))
    }
    regions_to_plot <- intersect(regions, available_regions)
    missing_reg <- setdiff(intersect(regions, REGIONS), available_regions)
    if (length(missing_reg) > 0) {
      message("Regions absent from selected features (skipped): ",
              paste(missing_reg, collapse = ", "))
    }
  }

  if (length(regions_to_plot) == 0) {
    stop("No regions to plot — check `regions` and feature group selection")
  }

  # --- Multi-species split --------------------------------------------------
  species_vec <- unique(df$species)
  has_species <- length(species_vec) > 1

  .process_one_df <- function(sub_df, reg, sp_label = NULL) {
    reg_features <- names(col_to_region)[
      vapply(col_to_region, identical, logical(1), reg)
    ]
    reg_features <- intersect(all_features, reg_features)

    if (length(reg_features) == 0) {
      message("  No features for region '", reg, "' — skipping")
      return(NULL)
    }

    sub_df <- sub_df[!is.na(sub_df[[response]]), , drop = FALSE]

    # Apply top-N per group filter (uses |ρ| to rank within each group)
    reg_features <- .apply_top_n(
      reg_features, col_to_group, top_n, response, sub_df
    )

    if (length(reg_features) == 0) {
      message("  No features remain after top_n filter for region '", reg, "'")
      return(NULL)
    }

    # Arrange features in display order
    feat_order <- .arrange_features(
      reg_features, response, feature_order, sub_df, cluster = cluster
    )
    all_cols <- feat_order  # response is last element

    n_cols  <- length(all_cols)
    n_pairs <- choose(n_cols, 2)
    sp_str  <- if (!is.null(sp_label)) paste0(" [", sp_label, "]") else ""
    message("  Region ", reg, sp_str, ": ",
            n_cols - 1, " features + response → ", n_pairs, " pairs")

    pairs_df <- .compute_corr_pairs(sub_df, all_cols, min_n = min_n)

    # Symmetrise + add diagonal (NA fill → neutral grey via na.value)
    sym_df <- dplyr::bind_rows(
      dplyr::mutate(pairs_df, is_diag = FALSE),
      dplyr::mutate(
        dplyr::rename(pairs_df,
                      feature_x = feature_y,
                      feature_y = feature_x),
        is_diag = FALSE
      ),
      tibble::tibble(
        feature_x = all_cols,
        feature_y = all_cols,
        rho       = NA_real_,
        p_value   = NA_real_,
        n         = nrow(sub_df),
        is_diag   = TRUE
      )
    )

    # BH correction (kept in table; not displayed in cells)
    sym_df <- sym_df |>
      dplyr::mutate(
        q_value = dplyr::if_else(
          is_diag | is.na(p_value),
          NA_real_,
          stats::p.adjust(p_value, method = "BH")
        )
      )

    # Cell label: ρ only; blank on diagonal and when below label_threshold
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

    # Factor levels → axis order
    sym_df$feature_x <- factor(sym_df$feature_x, levels = all_cols)
    sym_df$feature_y <- factor(sym_df$feature_y, levels = all_cols)

    # --- Triangle masking ---------------------------------------------------
    # Mark cells that belong to the excluded half so their geoms can be skipped.
    # With response at position 1 (bottom-left), the diagonal runs bottom-left
    # to top-right. "lower" keeps x_pos >= y_pos; "upper" keeps x_pos <= y_pos.
    sym_df <- sym_df |>
      dplyr::mutate(
        .xp = as.integer(feature_x),
        .yp = as.integer(feature_y),
        is_shown = is_diag | dplyr::case_when(
          triangle == "lower" ~ .xp >= .yp,
          triangle == "upper" ~ .xp <= .yp,
          TRUE                ~ TRUE
        )
      ) |>
      dplyr::select(-.xp, -.yp)

    # --- Group separator geometry -------------------------------------------
    feat_grp <- vapply(all_cols, function(f) {
      if (f == response) "response" else {
        g <- col_to_group[[f]]
        if (is.null(g)) "other" else g
      }
    }, character(1))

    rle_grp  <- rle(feat_grp)
    grp_ends <- cumsum(rle_grp$lengths)

    boundary_pos  <- grp_ends[-length(grp_ends)] + 0.5
    is_resp_bound <- which(
      rle_grp$values == "response" |
      c(rle_grp$values[-1] == "response", FALSE)
    )
    resp_boundary <- boundary_pos[is_resp_bound]
    # When clustering, groups are interleaved so suppress group lines —
    # only keep the response separator.
    feat_boundary <- if (cluster) character(0) else setdiff(boundary_pos, resp_boundary)

    # --- Axis labels: region suffix stripped for features -------------------
    # Response is first (position 1); features follow.
    axis_labels <- setNames(
      c(format_col_name(response),
        format_metric_name(all_cols[-1])),
      all_cols
    )

    # --- PDF dimensions (square) -------------------------------------------
    tile_mm <- max(12, min(22, 280 / n_cols))
    pdf_sq  <- if (!is.null(size_mm)) size_mm else n_cols * tile_mm + 80

    # --- Titles -------------------------------------------------------------
    region_disp <- REGION_DISPLAYS[[reg]]
    sp_title    <- if (!is.null(sp_label)) paste0(" — ", sp_label) else ""
    tn_str      <- if (!is.null(top_n)) {
      if (is.numeric(top_n) && length(top_n) == 1 && is.null(names(top_n))) {
        sprintf("; top %d per group", as.integer(top_n))
      } else {
        "; top-N per group"
      }
    } else ""

    tri_str     <- if (triangle != "full") paste0("; ", triangle, " triangle") else ""
    clust_str   <- if (cluster) "; clustered by ρ" else ""
    lbl_str     <- if (label_threshold > 0)
      sprintf("; labels |ρ| ≥ %.2f", label_threshold) else ""
    hm_title    <- sprintf("Feature correlation matrix — %s%s",
                           region_disp, sp_title)
    hm_subtitle <- sprintf(
      "Spearman ρ; %d features + %s%s%s%s%s",
      n_cols - 1, format_col_name(response), tn_str, tri_str, clust_str, lbl_str
    )

    # --- Build ggplot -------------------------------------------------------
    # For triangle modes, only the shown cells are passed to fill/text geoms;
    # the unused half is left blank (panel background shows through).
    shown_df <- dplyr::filter(sym_df, is_shown & !is_diag)

    p <- ggplot2::ggplot(sym_df,
                         ggplot2::aes(x = feature_x, y = feature_y)) +

      # Correlation tiles (shown cells only, excluding diagonal)
      ggplot2::geom_tile(
        data      = shown_df,
        ggplot2::aes(fill = rho),
        colour    = "white",
        linewidth = 0.25,
        na.rm     = TRUE
      ) +

      # Diagonal tiles (neutral grey; always shown)
      ggplot2::geom_tile(
        data      = dplyr::filter(sym_df, is_diag),
        fill      = "grey88",
        colour    = "white",
        linewidth = 0.25
      ) +

      # In-cell ρ labels (shown cells only)
      ggplot2::geom_text(
        data     = shown_df,
        ggplot2::aes(label = cell_label, colour = text_color),
        size     = cell_text_size,
        fontface = "bold",
        na.rm    = TRUE
      ) +

      # Feature group separator lines
      {if (length(feat_boundary) > 0)
        list(
          ggplot2::geom_hline(yintercept = feat_boundary,
                              colour = "grey40", linewidth = 0.7),
          ggplot2::geom_vline(xintercept = feat_boundary,
                              colour = "grey40", linewidth = 0.7)
        )
      } +

      # Response separator (heavier)
      {if (length(resp_boundary) > 0)
        list(
          ggplot2::geom_hline(yintercept = resp_boundary,
                              colour = "grey10", linewidth = 1.4),
          ggplot2::geom_vline(xintercept = resp_boundary,
                              colour = "grey10", linewidth = 1.4)
        )
      } +

      # Fill scale
      ggplot2::scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = color_limits,
        na.value = "grey88",
        name     = "Spearman ρ",
        guide    = ggplot2::guide_colorbar(
          title.position  = "top",
          barwidth        = ggplot2::unit(0.8, "cm"),
          barheight       = ggplot2::unit(5.5, "cm"),
          ticks.linewidth = 0.5
        )
      ) +

      # Text colour: identity (column holds literal colour strings)
      ggplot2::scale_colour_identity(guide = "none") +

      # Axis labels (region suffix stripped)
      ggplot2::scale_x_discrete(
        labels = axis_labels,
        expand = expansion(add = 0.5)
      ) +
      ggplot2::scale_y_discrete(
        labels = axis_labels,
        expand = expansion(add = 0.5)
      ) +

      ggplot2::labs(
        title    = hm_title,
        subtitle = hm_subtitle,
        x = NULL, y = NULL
      ) +

      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        aspect.ratio     = 1,
        axis.text.x      = ggplot2::element_text(
          angle = 45, vjust = 1, hjust = 1,
          size  = base_size - 2
        ),
        axis.text.y      = ggplot2::element_text(size = base_size - 2),
        plot.title       = ggplot2::element_text(
          size = base_size + 2, face = "bold"
        ),
        plot.subtitle    = ggplot2::element_text(
          size = base_size - 2, colour = "grey30"
        ),
        panel.grid       = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(
          colour = "grey20", fill = NA, linewidth = 0.6
        ),
        legend.position  = "right",
        legend.title     = ggplot2::element_text(
          size = base_size - 2, face = "bold"
        ),
        plot.margin      = ggplot2::margin(
          t = 8, r = 15, b = 8, l = 8, unit = "mm"
        )
      )

    # --- Save PDF -----------------------------------------------------------
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      fname <- if (!is.null(sp_label)) {
        sprintf("%s_%s_%s.pdf", file_prefix, reg, sp_label)
      } else {
        sprintf("%s_%s.pdf", file_prefix, reg)
      }
      pdf_path <- file.path(output_dir, fname)
      ggplot2::ggsave(
        pdf_path, plot = p,
        width = pdf_sq, height = pdf_sq,
        units = "mm", device = "pdf"
      )
      message("  Saved: ", pdf_path)
    }

    # --- Return table (unique pairs only, no diagonal) ----------------------
    # Response is at position 1 so x_pos < y_pos captures all unique pairs
    # with the response in the feature_x column (consistent key direction).
    table_out <- sym_df |>
      dplyr::filter(!is_diag) |>
      dplyr::filter(as.integer(feature_x) < as.integer(feature_y)) |>
      dplyr::select(feature_x, feature_y, rho, p_value, q_value, n)

    list(plot = p, table = table_out)
  }

  # --- Loop over regions (and species if multi-species) --------------------
  results <- list()

  for (reg in regions_to_plot) {
    message("Building heatmap(s) for region: ", reg)

    if (has_species) {
      for (sp in species_vec) {
        key    <- paste(reg, sp, sep = "_")
        sub_df <- dplyr::filter(df, species == sp)
        res    <- .process_one_df(sub_df, reg, sp_label = sp)
        if (!is.null(res)) results[[key]] <- res
      }
    } else {
      res <- .process_one_df(df, reg, sp_label = NULL)
      if (!is.null(res)) results[[reg]] <- res
    }
  }

  if (length(results) == 0) stop("No heatmaps produced.")
  invisible(results)
}


# =============================================================================
# Runner (executed when sourced directly, not when sourced by load_all.R)
# =============================================================================

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  df <- build_dataset("human")

  out_dir <- file.path(OUTPUT_DIR, "plots", "heatmaps")
  tab_dir <- file.path(OUTPUT_DIR, "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Job 1: Structure features, core regions, clustered ------------------
  message("\n=== Structure × core regions ===")
  out_struct <- region_feature_heatmap(
    df,
    response    = "halflife",
    groups      = "structure",
    regions     = c("5utr", "cds", "3utr"),
    cluster     = TRUE,
    output_dir  = out_dir,
    file_prefix = "region_heatmap_structure"
  )

  # ---- Job 2: Top-3 per group, core regions, clustered --------------------
  message("\n=== Top-3 per group × core regions ===")
  out_top3 <- region_feature_heatmap(
    df,
    response    = "halflife",
    groups      = c("structure", "intrinsic", "splicing", "regulatory",
                    "lengths_core", "nmd_reported"),
    regions     = c("5utr", "cds", "3utr"),
    top_n       = 3,
    cluster     = TRUE,
    output_dir  = out_dir,
    file_prefix = "region_heatmap_top3"
  )

  # Export tables
  for (key in names(out_top3)) {
    csv_path <- file.path(tab_dir, sprintf("region_heatmap_top3_%s.csv", key))
    write.csv(out_top3[[key]]$table, csv_path, row.names = FALSE)
    message("Table: ", csv_path)
  }

  message("\nDone. Plots: ", out_dir)
}
