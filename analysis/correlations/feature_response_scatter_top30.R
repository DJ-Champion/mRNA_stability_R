# =============================================================================
# Feature × two-response correlation scatter — top 30% labelled
# =============================================================================
# Single-plot variant of feature_response_scatter.R (the "included" plot; here
# halflife on x, TE on y). Differences from the original:
#   - label_quantile 0.7: the top 30% of points by distance from origin are
#     labelled (was 20%)
#   - larger axis titles and tick numbers
#   - both axes tick every 0.1, with a major gridline at each tick (no minor)
#
# The scatter function below is a copy of feature_response_scatter(); the
# original file runs its plots when sourced, so it is not sourced here. See it
# for full parameter documentation.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_response_scatter_top30.R")
#   out <- feature_response_scatter(build_dataset("human"), label_quantile = 0.7)
# =============================================================================

source("R/load_all.R")
#   source("analysis/correlations/feature_response_scatter.R")
#   df  <- build_dataset("human")
#
#   # Default: TE vs halflife, all features, per-column granularity
#   out <- feature_response_scatter(df)
#   print(out$plot)
#
#   # Cleaner: collapse to one point per (group, region), filter noise:
#   out <- feature_response_scatter(
#     df, collapse = "region", noise_filter = 0.1
#   )
#
#   # Restrict to structure features only:
#   out <- feature_response_scatter(df, include = "structure")
#
#   # Saluki diagnostic: which features explain Saluki's residuals?
#   out <- feature_response_scatter(
#     df,
#     response_x = "saluki_prediction",
#     response_y = "prediction_difference",
#     exclude_columns = c("^halflife$")
#   )
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(tibble)
})


#' Scatter of feature correlations against two responses.
#'
#' @param df             Dataframe from build_dataset() or build_all().
#' @param response_x     Character. Column for the x-axis correlation.
#' @param response_y     Character. Column for the y-axis correlation.
#' @param method         Correlation method (default "spearman").
#' @param include        Selection tokens: "core" (default), "exploratory",
#'                       "model", or supergroup / group / feature ids. See
#'                       select_features().
#' @param exclude        Tokens to subtract from `include`. NULL = none.
#' @param regions        Region tokens to keep ("5utr", "cds", ...). NULL
#'                       (default) = all.
#' @param collapse       "none" (one point per column, default), "region" (one
#'                       point per group × region — median r), or "group" (one
#'                       point per group — median r across all members).
#' @param max_features   Integer or NULL (default). When set, keep at most this
#'                       many points per (species, group) by distance from
#'                       origin: max_features = 3 shows only the three most
#'                       informative members per group. Applied after
#'                       noise_filter and `top_n`, before label selection.
#' @param top_n          Named list: feature id -> integer. Trims a family to
#'                       its N strongest members, so you can cap codon_freqs
#'                       at 2 while leaving other features unconstrained.
#'                       Ranked here by distance from origin (both responses),
#'                       max across species so panel selections stay
#'                       consistent; applied after noise_filter.
#' @param noise_filter   Numeric. Drop points with distance-from-origin below
#'                       this. 0 (default) = no filter. Try 0.1 to declutter.
#' @param label_quantile Numeric in [0, 1]. Label the top (1 - q) fraction by
#'                       distance from origin. Default 0.9 = label top 10%.
#' @param exclude_columns Regex patterns of columns to drop. Default = R10
#'                       derived-prediction set, unless response_x or
#'                       response_y is one of them (in which case it's
#'                       auto-removed from the list).
#' @param min_n          Minimum non-NA pairs to compute a correlation.
#' @param formatter      Display formatter (default format_col_name).
#' @param palette        Named colour vector keyed by feature id.
#'                       Default FEATURE_GROUP_COLOURS.
#' @param shapes         Named shape vector keyed by region token.
#'                       Default REGION_SHAPES.
#' @param axis_title_size Axis title font size.
#' @param axis_text_size  Axis tick-number font size.
#' @param axis_step       Spacing of tick labels and major gridlines on both
#'                       axes (default 0.1).
#' @param repel_args      Named list overriding ggrepel::geom_text_repel()
#'                       settings (see `repel_defaults` in the body).
#' @return list(plot, table). Table columns:
#'   species, variable, group, supergroup, region, n_x, correlation_x,
#'   p_value_x, q_x, n_y, correlation_y, p_value_y, q_y,
#'   distance_from_origin, labelled.
#' @export
feature_response_scatter <- function(df,
                                     response_x     = "halflife",
                                     response_y     = "translation_efficiency",
                                     method         = c("spearman", "pearson",
                                                        "kendall"),
                                     include        = "core",
                                     exclude        = NULL,
                                     regions        = NULL,
                                     collapse       = c("none", "region",
                                                        "group"),
                                     max_features   = NULL,
                                     top_n          = NULL,
                                     noise_filter   = 0,
                                     label_quantile = 0.9,
                                     exclude_columns = c("^saluki_prediction$",
                                                         "^prediction_difference$"),
                                     min_n          = 30,
                                     formatter      = format_col_name,
                                     palette        = NULL,
                                     shapes         = NULL,
                                     axis_title_size = 16,
                                     axis_text_size  = 14,
                                     axis_step       = 0.1,
                                     repel_args      = list()) {

  method   <- match.arg(method)
  collapse <- match.arg(collapse)
  if (is.null(palette)) palette <- FEATURE_GROUP_COLOURS
  if (is.null(shapes))  shapes  <- REGION_SHAPES

  # --- R5: guard ----------------------------------------------------------
  if (!response_x %in% names(df)) stop("response_x '", response_x, "' not in df")
  if (!response_y %in% names(df)) stop("response_y '", response_y, "' not in df")
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }
  if (response_x == response_y) {
    stop("response_x and response_y must differ")
  }
  if (label_quantile < 0 || label_quantile >= 1) {
    stop("label_quantile must be in [0, 1)")
  }
  if (!is.null(max_features) && (!is.numeric(max_features) ||
                                 length(max_features) != 1 ||
                                 max_features < 1 ||
                                 max_features != floor(max_features))) {
    stop("max_features must be a positive integer or NULL")
  }

  # Don't exclude a response if the user is correlating against it
  exclude_columns <- exclude_columns[!vapply(exclude_columns, function(rgx) {
    grepl(rgx, response_x) || grepl(rgx, response_y)
  }, logical(1))]

  # --- Enumerate candidate columns ----------------------------------------
  # Family trimming (`top_n`) is done below, on distance from origin, because
  # this figure has two responses; select_features() ranks against one.
  sel <- select_features(df, include, exclude, regions = regions)
  col_to_group <- stats::setNames(as.list(sel$feature_id), sel$column)

  candidates <- names(col_to_group)

  # Apply column exclusions
  for (rgx in exclude_columns) {
    candidates <- candidates[!vapply(candidates, function(c) grepl(rgx, c),
                                     logical(1))]
  }

  # Never include the responses themselves
  candidates <- setdiff(candidates, c(response_x, response_y))

  if (length(candidates) == 0) {
    stop("No candidate features after filtering — check `include`, `exclude` and `exclude_columns`")
  }

  # --- Per-species correlation computation --------------------------------
  has_species <- length(unique(df$species)) > 1

  compute_one <- function(sub, sp_label) {
    purrr::map_dfr(candidates, function(co) {
      v  <- sub[[co]]
      vx <- sub[[response_x]]
      vy <- sub[[response_y]]

      ok_x <- !is.na(v) & !is.na(vx)
      ok_y <- !is.na(v) & !is.na(vy)

      if (sum(ok_x) < min_n || sum(ok_y) < min_n) return(tibble::tibble())
      if (length(unique(v[ok_x])) < 2 ||
          length(unique(v[ok_y])) < 2) return(tibble::tibble())

      ct_x <- suppressWarnings(stats::cor.test(
        v[ok_x], vx[ok_x], method = method, exact = FALSE
      ))
      ct_y <- suppressWarnings(stats::cor.test(
        v[ok_y], vy[ok_y], method = method, exact = FALSE
      ))

      tibble::tibble(
        species       = sp_label,
        variable      = co,
        group         = col_to_group[[co]],
        n_x           = sum(ok_x),
        correlation_x = unname(ct_x$estimate),
        p_value_x     = ct_x$p.value,
        n_y           = sum(ok_y),
        correlation_y = unname(ct_y$estimate),
        p_value_y     = ct_y$p.value
      )
    })
  }

  result <- if (has_species) {
    purrr::map_dfr(unique(df$species), function(sp) {
      compute_one(df |> dplyr::filter(species == sp), sp)
    })
  } else {
    compute_one(df, unique(df$species)[1])
  }

  if (nrow(result) == 0) {
    stop("No correlations computed — try lowering min_n or check coverage")
  }

  # --- Region extraction (REGIONS-aware) ----------------------------------
  # A region-less feature (te) sits in the mrna slot, as in the other
  # region-aware plots; anything else with no region token stays "none".
  result <- result |>
    dplyr::mutate(
      region     = dplyr::coalesce(column_regions(variable, group), "none"),
      supergroup = dplyr::coalesce(supergroup_of(group), "other")
    )

  # --- BH q-values, per species, per response axis ------------------------
  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      q_x = stats::p.adjust(p_value_x, method = "BH"),
      q_y = stats::p.adjust(p_value_y, method = "BH")
    ) |>
    dplyr::ungroup()

  # --- Collapse if requested ----------------------------------------------
  if (collapse != "none") {
    grp_cols <- switch(collapse,
                       region = c("species", "group", "region"),
                       group  = c("species", "group"))

    result <- result |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
      dplyr::summarise(
        n_variables       = dplyr::n(),
        representative    = variable[which.max(abs(correlation_x) +
                                               abs(correlation_y))],
        correlation_x     = stats::median(correlation_x, na.rm = TRUE),
        correlation_y     = stats::median(correlation_y, na.rm = TRUE),
        n_x               = stats::median(n_x, na.rm = TRUE),
        n_y               = stats::median(n_y, na.rm = TRUE),
        q_x               = stats::median(q_x, na.rm = TRUE),
        q_y               = stats::median(q_y, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        # Variable name: synthetic, useful for the table identifier
        variable   = if (collapse == "region") {
          paste(group, region, sep = "__")
        } else {
          group
        },
        # Recompute supergroup (collapse may have dropped it)
        supergroup = dplyr::coalesce(supergroup_of(group), "other"),
        # Drop region in group-only collapse
        region     = if (collapse == "group") "none" else region
      )
  }

  # --- Distance + noise filter --------------------------------------------
  result <- result |>
    dplyr::mutate(
      distance_from_origin = sqrt(correlation_x^2 + correlation_y^2)
    ) |>
    dplyr::filter(distance_from_origin >= noise_filter)

  if (nrow(result) == 0) {
    stop("All points filtered out — noise_filter too high?")
  }

  # --- Family trim (top_n) --------------------------------------------------
  # Rank by max distance across species so panel selections stay consistent.
  if (length(top_n) > 0) {
    for (g in names(top_n)) {
      n_keep <- top_n[[g]]
      if (!any(result$group == g)) next

      keep_vars <- result |>
        dplyr::filter(group == g) |>
        dplyr::group_by(variable) |>
        dplyr::summarise(max_dist = max(distance_from_origin, na.rm = TRUE),
                         .groups = "drop") |>
        dplyr::arrange(dplyr::desc(max_dist)) |>
        dplyr::slice_head(n = n_keep) |>
        dplyr::pull(variable)

      result <- result |>
        dplyr::filter(group != g | variable %in% keep_vars)
    }
  }

  # --- max_features: keep top N per (species, group) by distance -----------
  if (!is.null(max_features)) {
    result <- result |>
      dplyr::group_by(species, group) |>
      dplyr::slice_max(order_by = distance_from_origin,
                       n = max_features, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  if (nrow(result) == 0) {
    stop("All points filtered out — max_features too restrictive?")
  }

  # --- Label selection by quantile ----------------------------------------
  # Per species: threshold at the requested quantile of distance.
  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      .thr     = stats::quantile(distance_from_origin,
                                 probs = label_quantile, na.rm = TRUE),
      labelled = distance_from_origin >= .thr
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.thr)

  # --- Display labels (R4) -------------------------------------------------
  # group_label() dispatches per element: selection keys (feature / group /
  # supergroup) via format_group_name(); anything else via format_col_name().
  # The key sets are read from the registries rather than hardcoded, so adding
  # a supergroup cannot silently mislabel it.
  group_label <- function(g) {
    selection_keys <- c(names(FEATURE_PATTERNS), names(FEATURE_GROUPS),
                        names(SUPERGROUPS), "other")
    ifelse(g %in% selection_keys,
           format_group_name(g, kind = "auto"),
           formatter(g))
  }

  result <- result |>
    dplyr::mutate(
      display_label = dplyr::case_when(
        collapse == "group"  ~ group_label(group),
        collapse == "region" ~ paste0(
          group_label(group),
          ifelse(region == "none", "",
                 paste0(" — ", format_col_name(region)))
        ),
        TRUE                 ~ format_metric_name(variable)
      ),
      label_text = ifelse(labelled, display_label, NA_character_)
    )

  # --- Build the plot ------------------------------------------------------
  # Colour aesthetic is group; legend label reuses the same per-element
  # dispatch as the in-plot labels.
  # The legend uses each feature's short name; collapsed group / supergroup
  # keys have none, so they keep the in-plot label.
  # Placeholders in a short name (nt.<x>%) are stripped: the legend names the
  # feature, not one of its members.
  legend_labeller <- function(g) {
    short <- gsub("\\.?<[^>]*>", "", FEATURE_SHORT_NAMES[g])
    ifelse(g %in% names(FEATURE_SHORT_NAMES), short, group_label(g))
  }

  # Make the categorical axes factors so legend ordering is consistent
  group_order <- intersect(names(palette), unique(result$group))
  result$group_f <- factor(result$group, levels = group_order)
  shape_order <- intersect(names(shapes), unique(result$region))
  result$region_f <- factor(result$region, levels = shape_order)

  axis_x_lab <- sprintf("Correlation with %s", formatter(response_x))
  axis_y_lab <- sprintf("Correlation with %s", formatter(response_y))

  title <- sprintf("Feature correlations: %s vs %s",
                   formatter(response_x), formatter(response_y))
  subtitle_bits <- character()
  if (collapse != "none") {
    subtitle_bits <- c(subtitle_bits,
                       paste0("collapsed by ", collapse,
                              " (median ", method, ")"))
  } else {
    subtitle_bits <- c(subtitle_bits, paste0("per-column ", method))
  }
  if (noise_filter > 0) {
    subtitle_bits <- c(subtitle_bits,
                       sprintf("noise filter |r| >= %.2f", noise_filter))
  }
  if (!is.null(max_features)) {
    subtitle_bits <- c(subtitle_bits,
                       sprintf("top %d per group", as.integer(max_features)))
  }
  subtitle_bits <- c(subtitle_bits,
                     sprintf("top %d%% labelled",
                             round((1 - label_quantile) * 100)))

  # Ticks every `axis_step` on both axes, covering the data range
  tick_seq <- function(v) {
    seq(floor(min(v, 0) / axis_step + 1e-9) * axis_step,
        ceiling(max(v, 0) / axis_step - 1e-9) * axis_step,
        by = axis_step)
  }
  fmt_tick <- function(x) sprintf("%.1f", round(x, 10))

  # Label placement. force_pull < 1 weakens the pull back towards each point,
  # so labels can sit further away (longer connecting segments); higher
  # max.iter / max.time let the solver keep searching for free positions.
  # Fixed seed = reproducible layout.
  repel_defaults <- list(
    size = 4, max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.2,
    force = 1, force_pull = 0.15,
    max.iter = 50000, max.time = 10,
    min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6,
    seed = 2, show.legend = FALSE, na.rm = TRUE
  )
  repel_args_final <- utils::modifyList(repel_defaults, repel_args)

  p <- ggplot2::ggplot(
    result,
    ggplot2::aes(x = correlation_x, y = correlation_y)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(colour = group_f, shape = region_f),
      size = 3, alpha = 0.85, stroke = 0.7
    ) +
    do.call(ggrepel::geom_text_repel, c(
      list(mapping = ggplot2::aes(label = label_text, colour = group_f)),
      repel_args_final
    )) +
    ggplot2::scale_colour_manual(
      values = palette,
      labels = legend_labeller,
      name   = "Feature group",
      drop   = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = shapes,
      labels = function(r) ifelse(r %in% names(REGION_DISPLAYS) &
                                  nzchar(REGION_DISPLAYS[r]),
                                  REGION_DISPLAYS[r], formatter(r)),
      name   = "Region",
      drop   = TRUE
    ) +
    ggplot2::scale_x_continuous(breaks = tick_seq(result$correlation_x),
                                labels = fmt_tick) +
    ggplot2::scale_y_continuous(breaks = tick_seq(result$correlation_y),
                                labels = fmt_tick) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title    = title,
      #subtitle = paste(subtitle_bits, collapse = " · "),
      x        = axis_x_lab,
      y        = axis_y_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30"),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(face = "bold"),
      legend.text      = ggplot2::element_text(size = 12),
      legend.key.size  = ggplot2::unit(0.8, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey82",
                                               linewidth = 0.3),
      axis.title       = ggplot2::element_text(size = axis_title_size),
      axis.text        = ggplot2::element_text(size = axis_text_size,
                                               colour = "black"),
      panel.background  = ggplot2::element_rect(fill = "grey95")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, ncol = 1,
                                     override.aes = list(size = 3, alpha = 1)),
      shape  = ggplot2::guide_legend(order = 2,
                                     override.aes = list(size = 3))
    )

  if (has_species) {
    p <- p + ggplot2::facet_wrap(~ species)
  }

  # --- R9: return table without the plot-only helper columns -------------
  table_out <- result |>
    dplyr::select(dplyr::any_of(c(
      "species", "variable", "group", "supergroup", "region",
      "n_variables", "representative",
      "n_x", "correlation_x", "p_value_x", "q_x",
      "n_y", "correlation_y", "p_value_y", "q_y",
      "distance_from_origin", "labelled", "display_label"
    )))

  list(plot = p, table = table_out)
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  df <- build_dataset("human")

  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)

  out <- feature_response_scatter(
    df,
    top_n          = list(codon_freqs = 2, aa_freqs = 2),
    noise_filter   = 0,
    label_quantile = 0.7
  )

  print(out$plot)

  base <- "feature_response_scatter_halflife_vs_te_top30"
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots", paste0(base, ".jpg")),
    plot = out$plot, width = 300, height = 220, units = "mm", dpi = 300
  )
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots", paste0(base, ".pdf")),
    plot = out$plot, width = 300, height = 220, units = "mm",
    device = grDevices::cairo_pdf
  )
  write.csv(out$table,
            file.path(OUTPUT_DIR, "tables", paste0(base, ".csv")),
            row.names = FALSE)

  message("\nFeature response scatter (top 30% labelled) complete: ",
          file.path(OUTPUT_DIR, "plots", base), ".{jpg,pdf}")
}
