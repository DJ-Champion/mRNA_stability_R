# =============================================================================
# Region feature heatmap: pairwise Spearman correlation matrix per region
# =============================================================================
# Builds a full feature × feature Spearman correlation matrix for each requested
# region. The response variable (default: halflife) is included as a column so
# you can read off every feature's relationship with the predictor, while also
# seeing how features co-vary with each other.
#
# Layout choices (all controllable via arguments):
#   * One PDF per region.
#   * Features ordered by group (FEATURE_PATTERNS/SUPERGROUPS canonical order)
#     then by |ρ with response| descending within each group.
#   * Response column sits at the bottom-right, visually separated.
#   * Diverging blue (−1) → white (0) → red (+1) fill; limits fixed at ±1.
#   * Spearman ρ printed in every off-diagonal cell; BH significance stars
#     (* q<0.05, ** q<0.01, *** q<0.001) appended when show_significance=TRUE.
#   * Group separator lines and left-margin group labels with group colours.
#   * Diagonal cells rendered in neutral grey (self-correlation is trivial).
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/region_feature_heatmap.R")
#   df  <- build_dataset("human")
#
#   # Structure features, CDS only
#   out <- region_feature_heatmap(df, groups = "structure", regions = "cds")
#   print(out[["cds"]]$plot)
#
#   # Core regions, explicit group + bundle selection, save PDFs
#   out <- region_feature_heatmap(
#     df,
#     groups     = c("structure", "intrinsic", "lengths_core"),
#     regions    = c("5utr", "cds", "3utr"),
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

.sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q)  ~ "",
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    TRUE      ~ ""
  )
}


#' Compute Spearman ρ and p-value for every unique pair in `cols`.
#' Returns a long data frame with columns:
#'   feature_x, feature_y, rho, p_value, n
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


#' Determine feature order for one region's heatmap.
#' Primary sort: group position in FEATURE_PATTERNS canonical order.
#' Secondary sort within group: |ρ with response| descending.
#' Returns the ordered character vector; response is always last.
.order_region_features <- function(reg_features, response, col_to_group, df_sub) {
  group_order <- names(FEATURE_PATTERNS)
  resp_vals   <- df_sub[[response]]

  group_rank <- vapply(reg_features, function(f) {
    g <- col_to_group[[f]]
    if (is.null(g)) Inf else {
      idx <- match(g, group_order)
      if (is.na(idx)) Inf else idx
    }
  }, numeric(1))

  abs_rho_with_resp <- vapply(reg_features, function(f) {
    x  <- df_sub[[f]]
    ok <- !is.na(x) & !is.na(resp_vals)
    if (sum(ok) < 5) return(0)
    abs(suppressWarnings(
      stats::cor(x[ok], resp_vals[ok], method = "spearman")
    ))
  }, numeric(1))

  ord <- order(group_rank, -abs_rho_with_resp)
  c(reg_features[ord], response)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Pairwise Spearman correlation heatmaps, one per region.
#'
#' Each heatmap shows the full feature × feature correlation matrix for the
#' columns belonging to one region, with the response variable included as a
#' special row/column. Features are ordered by group then |ρ with response|;
#' the response sits at the bottom-right.
#'
#' @param df               Dataframe from build_dataset() or build_all().
#' @param response         Response column to include in the matrix
#'                         (default "halflife").
#' @param groups           Groups / supergroups / bundles to include
#'                         (NULL = all groups in FEATURE_PATTERNS).
#' @param pick             Named list: group key → column allow-list.
#' @param drop             Named list: group key → columns to remove.
#' @param regions          Character vector of region tokens to plot
#'                         (NULL = every region present in the resolved features,
#'                         in REGIONS canonical order).
#' @param min_n            Minimum non-NA pairs needed to compute a correlation
#'                         (default 30).
#' @param color_limits     Numeric length-2: fill scale limits (default c(-1,1)).
#' @param show_significance Logical. Append BH significance stars to the ρ
#'                         label in significant cells (default TRUE).
#' @param sig_threshold    Q-value threshold used for the subtitle annotation
#'                         (default 0.05; stars always use 0.05/0.01/0.001).
#' @param cell_text_size   Size passed to geom_text for in-cell labels
#'                         (default 3). Reduce for large feature sets.
#' @param base_size        Base font size in points (default 14).
#' @param output_dir       If non-NULL, write one PDF per region here.
#' @param file_prefix      Filename prefix for saved PDFs
#'                         (default "region_heatmap").
#' @param width_mm         PDF width in mm. NULL = auto (scales with N features).
#' @param height_mm        PDF height in mm. NULL = auto.
#' @return Named list keyed by region (or "{region}_{species}" for multi-species
#'   datasets). Each value is list(plot = ggplot, table = tibble).
#'   Table columns: feature_x, feature_y, rho, p_value, q_value, n.
#'   Returned invisibly; call print() on $plot to display.
#' @export
region_feature_heatmap <- function(df,
                                   response          = "halflife",
                                   groups            = NULL,
                                   pick              = list(),
                                   drop              = list(),
                                   regions           = NULL,
                                   min_n             = 30,
                                   color_limits      = c(-1, 1),
                                   show_significance = TRUE,
                                   sig_threshold     = 0.05,
                                   cell_text_size    = 3,
                                   base_size         = 14,
                                   output_dir        = NULL,
                                   file_prefix       = "region_heatmap",
                                   width_mm          = NULL,
                                   height_mm         = NULL) {

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
  sel         <- resolve_selection(groups, pick, drop)
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

    # Drop rows where response is NA
    sub_df <- sub_df[!is.na(sub_df[[response]]), , drop = FALSE]

    feat_order <- .order_region_features(
      reg_features, response, col_to_group, sub_df
    )
    all_cols <- feat_order  # includes response at end

    n_cols <- length(all_cols)
    n_pairs <- choose(n_cols, 2)
    sp_str  <- if (!is.null(sp_label)) paste0(" [", sp_label, "]") else ""
    message("  Region ", reg, sp_str, ": ",
            n_cols - 1, " features + response → ",
            n_pairs, " pairs")

    pairs_df <- .compute_corr_pairs(sub_df, all_cols, min_n = min_n)

    # Symmetrise: add reverse pairs + diagonal (NA fill → grey)
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

    # BH correction across all off-diagonal non-NA p-values
    sym_df <- sym_df |>
      dplyr::mutate(
        q_value = dplyr::if_else(
          is_diag | is.na(p_value),
          NA_real_,
          stats::p.adjust(p_value, method = "BH")
        )
      )

    # Cell label: ρ + optional significance stars
    sym_df <- sym_df |>
      dplyr::mutate(
        rho_str    = dplyr::if_else(
          is_diag | is.na(rho), "",
          sprintf("%.2f", rho)
        ),
        stars      = dplyr::if_else(
          !is_diag & show_significance,
          .sig_stars(q_value),
          ""
        ),
        cell_label = paste0(rho_str, stars),
        # White text on strongly coloured tiles; black on pale tiles
        text_color = dplyr::case_when(
          is_diag             ~ "grey60",
          is.na(rho)          ~ "grey60",
          abs(rho) > 0.62     ~ "white",
          TRUE                ~ "grey10"
        )
      )

    # Factor levels control axis order
    sym_df$feature_x <- factor(sym_df$feature_x, levels = all_cols)
    sym_df$feature_y <- factor(sym_df$feature_y, levels = all_cols)

    # --- Group annotation geometry -----------------------------------------
    feat_grp <- vapply(all_cols, function(f) {
      if (f == response) "response" else {
        g <- col_to_group[[f]]
        if (is.null(g)) "other" else g
      }
    }, character(1))

    rle_grp   <- rle(feat_grp)
    grp_ends  <- cumsum(rle_grp$lengths)
    grp_starts <- c(1L, head(grp_ends, -1L) + 1L)
    grp_mids  <- (grp_starts + grp_ends) / 2

    # Separator line positions: between groups (not at response boundary —
    # the response gets its own thicker separator added separately)
    is_resp_boundary <- which(
      rle_grp$values == "response" |
      c(rle_grp$values[-1] == "response", FALSE)
    )
    boundary_pos <- grp_ends[-length(grp_ends)] + 0.5
    resp_boundary <- boundary_pos[is_resp_boundary]
    feat_boundary <- setdiff(boundary_pos, resp_boundary)

    grp_colors <- ifelse(
      rle_grp$values == "response",
      "grey20",
      feature_colour(rle_grp$values)
    )
    grp_labels <- ifelse(
      rle_grp$values == "response",
      format_col_name(response),
      format_group_name(rle_grp$values)
    )

    bar_df <- data.frame(
      ymin  = grp_starts - 0.5,
      ymax  = grp_ends   + 0.5,
      color = grp_colors,
      label = grp_labels,
      ymid  = grp_mids,
      stringsAsFactors = FALSE
    )

    # Axis display labels
    axis_labels <- setNames(format_col_name(all_cols), all_cols)

    # Auto-size PDF dimensions (mm) from feature count
    tile_mm <- max(12, min(20, 260 / n_cols))
    pdf_w   <- if (!is.null(width_mm))  width_mm  else n_cols * tile_mm + 90
    pdf_h   <- if (!is.null(height_mm)) height_mm else n_cols * tile_mm + 90

    # Left margin: determined by longest group label (rough heuristic)
    max_label_nchar <- max(nchar(grp_labels))
    left_margin_mm  <- max(40, max_label_nchar * 2.2 + 12)

    # Positions for left-margin annotations (discrete axis: tile i sits at x = i)
    # With coord_cartesian(clip = "off"), negative x positions draw into the margin.
    bar_xmax <- -0.65
    bar_xmin <- bar_xmax - 0.35
    lbl_x    <- bar_xmin - 0.15

    # Build title
    region_disp <- REGION_DISPLAYS[[reg]]
    sp_title    <- if (!is.null(sp_label)) paste0(" — ", sp_label) else ""
    hm_title    <- sprintf("Feature correlation matrix — %s%s", region_disp, sp_title)
    hm_subtitle <- sprintf(
      "Spearman ρ; %d features + %s; group order → |ρ| with %s; BH correction",
      n_cols - 1, format_col_name(response), format_col_name(response)
    )

    # --- Build ggplot -------------------------------------------------------
    p <- ggplot2::ggplot(sym_df,
                         ggplot2::aes(x = feature_x, y = feature_y)) +

      # Correlation tiles
      ggplot2::geom_tile(
        ggplot2::aes(fill = rho),
        colour = "white", linewidth = 0.25,
        na.rm  = TRUE
      ) +

      # Diagonal grey tiles (rho = NA → na.value; this just makes it explicit)
      ggplot2::geom_tile(
        data    = dplyr::filter(sym_df, is_diag),
        fill    = "grey88", colour = "white", linewidth = 0.25
      ) +

      # In-cell labels
      ggplot2::geom_text(
        ggplot2::aes(label = cell_label, colour = text_color),
        size     = cell_text_size,
        fontface = "bold",
        na.rm    = TRUE
      ) +

      # Feature group separator lines (between groups, moderate weight)
      {if (length(feat_boundary) > 0)
        list(
          ggplot2::geom_hline(yintercept = feat_boundary,
                              colour = "grey30", linewidth = 0.7),
          ggplot2::geom_vline(xintercept = feat_boundary,
                              colour = "grey30", linewidth = 0.7)
        )
      } +

      # Response separator lines (heavier)
      {if (length(resp_boundary) > 0)
        list(
          ggplot2::geom_hline(yintercept = resp_boundary,
                              colour = "grey10", linewidth = 1.4),
          ggplot2::geom_vline(xintercept = resp_boundary,
                              colour = "grey10", linewidth = 1.4)
        )
      } +

      # Left margin: coloured group bars (drawn as thick line segments at x < 1)
      ggplot2::annotate(
        "rect",
        xmin  = rep(bar_xmin, nrow(bar_df)),
        xmax  = rep(bar_xmax, nrow(bar_df)),
        ymin  = bar_df$ymin,
        ymax  = bar_df$ymax,
        fill  = bar_df$color,
        colour = NA
      ) +

      # Left margin: group name labels
      ggplot2::annotate(
        "text",
        x        = rep(lbl_x, nrow(bar_df)),
        y        = bar_df$ymid,
        label    = bar_df$label,
        hjust    = 1, vjust = 0.5,
        size     = base_size * 0.28,
        fontface = "bold",
        colour   = bar_df$color
      ) +

      # Fill scale
      ggplot2::scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = color_limits,
        oob      = scales::squish,
        na.value = "grey88",
        name     = "Spearman ρ",
        guide    = ggplot2::guide_colorbar(
          title.position = "top",
          barwidth       = ggplot2::unit(0.8, "cm"),
          barheight      = ggplot2::unit(5.5, "cm"),
          ticks.linewidth = 0.5
        )
      ) +

      # Text colour scale (identity: the column contains literal colour strings)
      ggplot2::scale_colour_identity(guide = "none") +

      ggplot2::scale_x_discrete(labels = axis_labels) +
      ggplot2::scale_y_discrete(labels = axis_labels) +

      ggplot2::coord_cartesian(clip = "off") +

      ggplot2::labs(
        title    = hm_title,
        subtitle = hm_subtitle,
        x = NULL, y = NULL
      ) +

      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        axis.text.x       = ggplot2::element_text(
          angle = 45, vjust = 1, hjust = 1,
          size  = base_size - 2
        ),
        axis.text.y       = ggplot2::element_text(size = base_size - 2),
        plot.title        = ggplot2::element_text(
          size = base_size + 2, face = "bold"
        ),
        plot.subtitle     = ggplot2::element_text(
          size = base_size - 2, colour = "grey30"
        ),
        panel.grid        = ggplot2::element_blank(),
        panel.border      = ggplot2::element_rect(
          colour = "grey30", fill = NA, linewidth = 0.6
        ),
        legend.position   = "right",
        legend.title      = ggplot2::element_text(
          size = base_size - 2, face = "bold"
        ),
        plot.margin       = ggplot2::margin(
          t = 8, r = 15, b = 8, l = left_margin_mm, unit = "mm"
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
        width = pdf_w, height = pdf_h,
        units = "mm", device = "pdf"
      )
      message("  Saved: ", pdf_path)
    }

    # --- Return table (off-diagonal pairs only) ----------------------------
    table_out <- sym_df |>
      dplyr::filter(!is_diag) |>
      dplyr::select(feature_x, feature_y, rho, p_value, q_value, n) |>
      dplyr::filter(as.integer(feature_x) < as.integer(feature_y))  # upper triangle

    list(plot = p, table = table_out)
  }

  # --- Loop over regions (and species if multi-species) -------------------
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

  # ---- Job 1: Structure features, core regions ----------------------------
  message("\n=== Structure × core regions ===")
  out_struct <- region_feature_heatmap(
    df,
    response   = "halflife",
    groups     = "structure",
    regions    = c("5utr", "cds", "3utr"),
    output_dir = out_dir,
    file_prefix = "region_heatmap_structure"
  )

  # ---- Job 2: All supergroups, each core region ---------------------------
  message("\n=== All groups × core regions ===")
  out_all <- region_feature_heatmap(
    df,
    response   = "halflife",
    groups     = c("structure", "intrinsic", "splicing", "regulatory",
                   "lengths_core", "nmd_reported"),
    regions    = c("5utr", "cds", "3utr"),
    output_dir = out_dir,
    file_prefix = "region_heatmap_all"
  )

  # Export tables
  for (key in names(out_all)) {
    csv_path <- file.path(tab_dir, sprintf("region_heatmap_%s.csv", key))
    write.csv(out_all[[key]]$table, csv_path, row.names = FALSE)
    message("Table: ", csv_path)
  }

  message("\nDone. Plots: ", out_dir)
}
