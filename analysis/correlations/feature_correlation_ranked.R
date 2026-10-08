# =============================================================================
# Ranked feature correlation plot
# =============================================================================
# A sibling of feature_correlation_dotplot.R, not a replacement — that script
# is unchanged and still produces the figures it always did. This one exists
# to try three things the original does not do:
#
#   1. ORIENTATION. `orientation = "horizontal"` puts feature names on the
#      y axis and the correlation on the x axis. Feature labels are long and
#      of wildly uneven length ("MFE.z 5' UTR" vs "NMD alt-stop mRNA"), which
#      is exactly the case where a rotated axis beats a 45-degree one: the
#      labels are set horizontally and flush right, so they read at a glance
#      and the figure grows DOWNWARD as features are added rather than
#      sideways. That is what makes the full 64-codon and 20-amino-acid
#      panels below practical at all.
#
#   2. SUPERGROUP FACETS. Every supergroup in the selection gets its own
#      facet, ordered by the strongest correlation it contains. There is no
#      collapse into "other": what a facet holds is what the feature table
#      says the supergroup holds.
#
#   3. A DERIVED SIGNIFICANCE LINE. The original draws a reference line at a
#      hardcoded 0.02. Here `sig_threshold = "auto"` computes the |r| at
#      which p = alpha for the actual n, using the SAME standard error that
#      builds the confidence intervals. The line and the error bars therefore
#      agree by construction: a point whose CI excludes zero is exactly a
#      point beyond the line. See critical_correlation() for the caveat about
#      which n is used.
#
# COHORT. Nothing here filters transcripts. build_dataset() applies
# the UTR length filters (config.R) to the frame it returns, so this script
# sees the reduced set automatically — 12,302 human rows against 13,660
# built. Deliberately, none of that is written onto the figures: the subtitle
# stays to the method and the CI level, and every cohort number, per-feature
# n and dropped column goes into the run report (see write_run_report()) for
# the figure legend to be written from later.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_correlation_ranked.R")
#   df  <- build_dataset("human")
#
#   # The core set, one facet per supergroup, features on the y axis; the two
#   # strongest codons and amino acids stand for their families
#   out <- feature_correlation_ranked(df, orientation = "horizontal",
#                                     top_n = list(codon_freqs = 2, aa_freqs = 2))
#   print(out$plot)
#
#   # All 64 codons, one region, ranked top to bottom
#   out <- feature_correlation_ranked(df, include = "codon_freqs",
#                                     orientation = "horizontal")
#
#   # Core minus the Sequence supergroup, coding sequence only
#   out <- feature_correlation_ranked(df, exclude = "sequence", regions = "cds")
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(tibble)
})


# -----------------------------------------------------------------------------
# Helper: correlation with CI for one (x, y) pair
# -----------------------------------------------------------------------------

#' Compute a correlation and a Fisher's-Z confidence interval.
#' Returns NA estimate/CI if too few non-NA pairs.
#'
#' Identical in behaviour to the function of the same name in
#' feature_correlation_dotplot.R. Kept local rather than shared so the two
#' scripts stay independently sourceable in either order — whichever is
#' sourced last defines it, and both definitions agree.
#'
#' @param x,y      Numeric vectors of equal length.
#' @param method   "spearman", "pearson", or "kendall".
#' @param conf     Confidence level (default 0.95).
#' @param min_n    Minimum non-NA pairs to attempt.
#' @return list(estimate, conf.low, conf.high, p.value, n).
correlation_with_ci <- function(x, y,
                                method = c("spearman", "pearson", "kendall"),
                                conf   = 0.95,
                                min_n  = 30) {
  method <- match.arg(method)
  ok <- !is.na(x) & !is.na(y)
  n  <- sum(ok)

  if (n < min_n ||
      length(unique(x[ok])) < 2 ||
      length(unique(y[ok])) < 2) {
    return(list(estimate = NA_real_, conf.low = NA_real_,
                conf.high = NA_real_, p.value = NA_real_, n = n))
  }

  ct <- suppressWarnings(
    stats::cor.test(x[ok], y[ok], method = method, exact = FALSE,
                    conf.level = conf)
  )

  r <- unname(ct$estimate)

  if (method == "pearson" && !is.null(ct$conf.int)) {
    return(list(
      estimate  = r,
      conf.low  = ct$conf.int[1],
      conf.high = ct$conf.int[2],
      p.value   = ct$p.value,
      n         = n
    ))
  }

  z   <- atanh(r)
  se  <- switch(
    method,
    spearman = 1.06 / sqrt(n - 3),
    kendall  = sqrt(0.437 / (n - 4))
  )
  zc  <- stats::qnorm(1 - (1 - conf) / 2)
  list(
    estimate  = r,
    conf.low  = tanh(z - zc * se),
    conf.high = tanh(z + zc * se),
    p.value   = ct$p.value,
    n         = n
  )
}


# -----------------------------------------------------------------------------
# Helper: the |r| that reaches significance at a given n
# -----------------------------------------------------------------------------

#' Smallest |correlation| distinguishable from zero at level `alpha`.
#'
#' Inverts the interval machinery in correlation_with_ci(), so the reference
#' line this produces and the error bars on the plot are two views of one
#' calculation: |r| exceeds this value exactly when the (1 - alpha) CI for
#' that point excludes zero. Pearson uses the exact t-distribution result;
#' Spearman and Kendall invert the same Fisher's-Z standard errors used
#' above (Bonett-Wright for Spearman).
#'
#' WHICH n TO PASS. n varies per feature — a probing column is measured on a
#' fraction of the transcripts a length column is. A single line cannot be
#' correct for all of them, and the two defensible choices say different
#' things. Passing the LARGEST n across the plotted features (what
#' feature_correlation_ranked() does) makes the line a floor: it is the most
#' permissive threshold on the figure, so anything falling short of it is not
#' significant for ANY feature plotted, while clearing it only guarantees
#' significance for the best-powered ones. That one-directional reading is
#' the honest one for a single line, and the run report carries the n range
#' and the implied threshold at each end so the figure legend can state it.
#'
#' @param n      Integer sample size.
#' @param method "spearman", "pearson", or "kendall".
#' @param alpha  Two-sided significance level (default 0.05).
#' @return Numeric critical |r|, or NA if n is too small for the method.
critical_correlation <- function(n,
                                 method = c("spearman", "pearson", "kendall"),
                                 alpha  = 0.05) {
  method <- match.arg(method)
  if (!is.finite(n)) return(NA_real_)

  if (method == "pearson") {
    if (n < 3) return(NA_real_)
    tc <- stats::qt(1 - alpha / 2, df = n - 2)
    return(tc / sqrt(tc^2 + n - 2))
  }

  se <- switch(
    method,
    spearman = if (n > 3) 1.06 / sqrt(n - 3)     else NA_real_,
    kendall  = if (n > 4) sqrt(0.437 / (n - 4))  else NA_real_
  )
  if (!is.finite(se)) return(NA_real_)
  tanh(stats::qnorm(1 - alpha / 2) * se)
}


# -----------------------------------------------------------------------------
# Main plot function
# -----------------------------------------------------------------------------

#' Feature correlations with a chosen response, ranked, with confidence
#' intervals, faceted by supergroup.
#'
#' @param df                 Dataframe from build_dataset() or build_all().
#' @param response           Character. Response column (default "halflife").
#' @param method             Correlation method (default "spearman").
#' @param include            Selection tokens: "core" (default), "exploratory",
#'                           "model", or supergroup / group / feature ids. See
#'                           select_features().
#' @param exclude            Tokens to subtract from `include`. NULL = none.
#' @param regions            Region tokens to draw ("5utr", "cds", ...). NULL
#'                           (default) = all. Applied before `top_n`. A named
#'                           list restricts features individually, e.g.
#'                           list(stopfree = "mrna"); see select_features().
#' @param orientation        "vertical" (default) puts features on the x axis,
#'                           as feature_correlation_dotplot does. "horizontal"
#'                           puts features on the y axis and correlations on
#'                           the x axis, and stacks the supergroup facets as
#'                           rows.
#' @param absolute           Logical. If TRUE (default) plot |correlation|
#'                           and transform the CI accordingly.
#' @param min_abs_correlation Numeric. Drop points where |r| is below this.
#' @param sig_threshold      "auto" (default) draws the reference line at the
#'                           critical |r| for `sig_alpha` at the largest n on
#'                           the figure — see critical_correlation(). A number
#'                           draws the line there instead. NULL = no line.
#' @param sig_alpha          Significance level for the "auto" line
#'                           (default 0.05).
#' @param conf               Confidence level for the CI (default 0.95).
#' @param min_n              Minimum non-NA pairs to compute a correlation.
#' @param formatter          Display formatter for the response label.
#' @param region_colours     Named colour vector. Default REGION_COLOURS.
#' @param region_shapes      Named shape vector. Default REGION_SHAPES.
#' @param top_n              Named list: feature id -> keep only the N metric
#'                           stems with the largest max |r| against `response`.
#'                           NULL = whole families. See select_features().
#' @return list(plot, table, report). `report` is the diagnostic list consumed
#'   by write_run_report(); nothing in it is drawn on the figure.
#' @export
feature_correlation_ranked <- function(df,
                                       response             = "halflife",
                                       method               = c("spearman",
                                                                "pearson",
                                                                "kendall"),
                                       include              = "core",
                                       exclude              = NULL,
                                       regions              = NULL,
                                       orientation          = c("vertical",
                                                                "horizontal"),
                                       absolute             = TRUE,
                                       min_abs_correlation  = 0,
                                       sig_threshold        = "auto",
                                       sig_alpha            = 0.05,
                                       conf                 = 0.95,
                                       min_n                = 30,
                                       formatter            = format_col_name,
                                       region_colours       = NULL,
                                       region_shapes        = NULL,
                                       top_n                = NULL) {

  method      <- match.arg(method)
  orientation <- match.arg(orientation)
  horizontal  <- identical(orientation, "horizontal")

  if (is.null(region_colours)) region_colours <- REGION_COLOURS
  if (is.null(region_shapes))  region_shapes  <- REGION_SHAPES

  if (!response  %in% names(df)) stop("response '", response, "' not in df")
  if (!"species" %in% names(df)) {
    stop("species column missing - pipeline invariant violated")
  }

  # --- Enumerate region-bearing columns -----------------------------------
  # Same contract as feature_correlation_dotplot: a column earns a place on
  # the feature axis by ending in a real REGIONS token, or by being a
  # registered region-less scalar (mapped to `mrna`). Anything else cannot be
  # positioned against the region legend and is reported, not silently lost.
  sel <- select_features(df, include, exclude, top_n = top_n,
                            response = response, regions = regions,
                            method = method, min_n = min_n)
  sel <- sel[sel$column != response, , drop = FALSE]
  ranking_applied <- attr(sel, "top_n")
  sel$region <- column_regions(sel$column, sel$feature_id)
  sel$stem   <- ifelse(is.na(sel$region), NA_character_,
                       ifelse(sub("^.*_", "", sel$column) %in% REGIONS &
                                grepl("_", sel$column, fixed = TRUE),
                              sub("_[^_]+$", "", sel$column), sel$column))
  plottable    <- sel[!is.na(sel$region), , drop = FALSE]
  col_to_group <- stats::setNames(as.list(plottable$feature_id), plottable$column)
  col_region   <- stats::setNames(as.list(plottable$region),     plottable$column)
  col_stem     <- stats::setNames(as.list(plottable$stem),       plottable$column)
  dropped      <- sel$column[is.na(sel$region)]

  candidates <- setdiff(names(col_to_group), response)
  if (length(candidates) == 0) {
    stop("No plottable columns after filtering - check `include` / `exclude` / `regions`")
  }

  if (length(dropped) > 0) {
    message("feature_correlation_ranked: dropped ", length(dropped),
            " column(s) with no region token (not plottable here)")
  }

  # --- Per-species correlation + CI ---------------------------------------
  has_species <- length(unique(df$species)) > 1

  compute_one <- function(sub, sp_label) {
    purrr::map_dfr(candidates, function(co) {
      r <- correlation_with_ci(
        sub[[co]], sub[[response]],
        method = method, conf = conf, min_n = min_n
      )
      if (is.na(r$estimate)) return(tibble::tibble())

      tibble::tibble(
        species     = sp_label,
        variable    = co,
        group       = col_to_group[[co]],
        region      = col_region[[co]],
        metric_stem = col_stem[[co]],
        n           = r$n,
        correlation = r$estimate,
        conf.low    = r$conf.low,
        conf.high   = r$conf.high,
        p_value     = r$p.value
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
    stop("No correlations computed - try lowering min_n or check coverage")
  }

  # Columns that resolved to a plottable name but produced no correlation:
  # too few non-NA pairs against this response, or no variance left. Tracked
  # for the report — under a row filter this is exactly the failure mode
  # worth watching, since a feature can quietly leave the figure.
  no_estimate <- setdiff(candidates, unique(result$variable))

  # --- Attach supergroup, BH q ---------------------------------------------
  result <- result |>
    dplyr::mutate(supergroup = dplyr::coalesce(supergroup_of(group), "other"))

  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(q_value = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()

  # --- |r| transform (preserving CI semantics) ----------------------------
  result <- result |>
    dplyr::mutate(
      correlation_abs = abs(correlation),
      conf.low_abs = dplyr::case_when(
        conf.low < 0 & conf.high > 0 ~ 0,
        conf.low < 0 & conf.high < 0 ~ abs(conf.high),
        TRUE                         ~ conf.low
      ),
      conf.high_abs = dplyr::case_when(
        conf.low < 0 & conf.high > 0 ~ pmax(abs(conf.low), abs(conf.high)),
        conf.low < 0 & conf.high < 0 ~ abs(conf.low),
        TRUE                         ~ conf.high
      )
    )

  # --- Filter by |r| threshold --------------------------------------------
  if (min_abs_correlation > 0) {
    result <- result |>
      dplyr::filter(correlation_abs >= min_abs_correlation)
    if (nrow(result) == 0) {
      stop("All points filtered out - min_abs_correlation too high?")
    }
  }

  # --- Pick the value aesthetic + CI columns ------------------------------
  if (absolute) {
    result$.value <- result$correlation_abs
    result$.lo    <- result$conf.low_abs
    result$.hi    <- result$conf.high_abs
    value_lab <- sprintf("Absolute %s correlation with %s",
                         tools::toTitleCase(method), formatter(response))
    value_breaks <- seq(0, 1, by = 0.1)
  } else {
    result$.value <- result$correlation
    result$.lo    <- result$conf.low
    result$.hi    <- result$conf.high
    value_lab <- sprintf("%s correlation with %s",
                         tools::toTitleCase(method), formatter(response))
    value_breaks <- ggplot2::waiver()
  }

  # --- Resolve the significance line --------------------------------------
  # "auto" reads the largest n actually on the figure, so the line tracks the
  # cohort rather than asserting a number chosen against an older one.
  n_max <- max(result$n, na.rm = TRUE)
  n_min <- min(result$n, na.rm = TRUE)
  sig_value <- NULL
  sig_auto  <- FALSE
  if (!is.null(sig_threshold)) {
    if (identical(sig_threshold, "auto")) {
      sig_value <- critical_correlation(n_max, method = method,
                                        alpha = sig_alpha)
      sig_auto  <- TRUE
    } else if (is.numeric(sig_threshold)) {
      sig_value <- sig_threshold
    } else {
      stop("sig_threshold must be \"auto\", a number, or NULL")
    }
  }

  # --- Display labels ------------------------------------------------------
  result$metric_display <- format_metric_name(result$variable)

  # Ordering: by max |correlation| within each supergroup facet. The
  # supergroup prefix keeps the factor level unique when the same display
  # name occurs in two facets.
  result <- result |>
    dplyr::mutate(metric_key = paste(supergroup, metric_display, sep = "::"))

  order_tbl <- result |>
    dplyr::group_by(metric_key) |>
    dplyr::summarise(max_r = max(correlation_abs, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_r))

  # Facet order: by the strongest correlation the panel contains, mirroring
  # the within-panel ordering.
  supergroup_order <- result |>
    dplyr::group_by(supergroup) |>
    dplyr::summarise(max_r = max(correlation_abs, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_r)) |>
    dplyr::pull(supergroup)

  # A discrete y axis is drawn bottom-up, so the level order that puts the
  # strongest feature on the LEFT when vertical puts it at the BOTTOM when
  # horizontal. Reverse it so "strongest first" reads the same way in both
  # orientations: leftmost when vertical, topmost when horizontal.
  key_levels <- order_tbl$metric_key
  if (horizontal) key_levels <- rev(key_levels)

  result <- result |>
    dplyr::mutate(
      metric_key = factor(metric_key, levels = key_levels),
      supergroup = factor(supergroup, levels = supergroup_order),
      region_f   = factor(region,
                          levels = intersect(names(region_colours),
                                             unique(region)))
    )

  # --- Plot ---------------------------------------------------------------
  dodge_width <- 0.5
  dodge       <- ggplot2::position_dodge(width = dodge_width)

  # Tick labels: strip the supergroup prefix carried only to keep levels
  # unique.
  label_lookup <- setNames(
    sub("^[^:]+::", "", levels(result$metric_key)),
    levels(result$metric_key)
  )

  if (horizontal) {
    p <- ggplot2::ggplot(
      result,
      ggplot2::aes(x = .value, y = metric_key,
                   colour = region_f, shape = region_f, fill = region_f)
    ) +
      # orientation = "y" dodges along the discrete axis and draws the caps
      # perpendicular to it; geom_errorbarh is deprecated as of ggplot2 3.4.
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .lo, xmax = .hi), orientation = "y",
        width = 0.2, position = dodge, linewidth = 0.4
      ) +
      ggplot2::geom_point(size = 2.5, position = dodge, stroke = 0.5) +
      ggplot2::scale_y_discrete(labels = label_lookup) +
      ggplot2::scale_x_continuous(breaks = value_breaks) +
      ggplot2::labs(x = value_lab, y = NULL)
  } else {
    p <- ggplot2::ggplot(
      result,
      ggplot2::aes(x = metric_key, y = .value,
                   colour = region_f, shape = region_f, fill = region_f)
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .lo, ymax = .hi),
        width = 0.2, position = dodge, linewidth = 0.4
      ) +
      ggplot2::geom_point(size = 2.5, position = dodge, stroke = 0.5) +
      ggplot2::scale_x_discrete(labels = label_lookup) +
      ggplot2::scale_y_continuous(breaks = value_breaks) +
      ggplot2::labs(x = NULL, y = value_lab)
  }

  p <- p +
    ggplot2::scale_colour_manual(
      values = region_colours, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r], name = "Region"
    ) +
    ggplot2::scale_fill_manual(
      values = region_colours, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r], name = "Region"
    ) +
    ggplot2::scale_shape_manual(
      values = region_shapes, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r], name = "Region"
    ) +
    # Deliberately succinct: method and interval only. The cohort, the UTR
    # filter and the per-feature n live in the run report and belong in the
    # figure legend, not baked into the image.
    ggplot2::labs(
      title    = sprintf("%s correlation with %s",
                         tools::toTitleCase(method), formatter(response)),
      subtitle = sprintf("%d%% CI; ordered by |r|", round(conf * 100))
    )

  # --- Theme ---------------------------------------------------------------
  # The grid follows the CONTINUOUS axis in both orientations: major lines
  # run across the value scale so a reader can trace a point back to a
  # correlation, and the discrete axis gets the light dashed separators.
  base_theme <- ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 16, colour = "grey30"),
      strip.text       = ggplot2::element_text(size = 16, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90",
                                               colour = "black"),
      legend.title     = ggplot2::element_text(size = 18, face = "bold"),
      legend.text      = ggplot2::element_text(size = 16),
      panel.background = ggplot2::element_rect(fill = "grey95")
    )

  if (horizontal) {
    p <- p + base_theme + ggplot2::theme(
      axis.title.x       = ggplot2::element_text(
        size = 18, face = "bold", margin = ggplot2::margin(t = 10)),
      # Set horizontally and flush right against the panel — the whole point
      # of the rotation.
      axis.text.y        = ggplot2::element_text(size = 14, hjust = 1),
      axis.text.x        = ggplot2::element_text(size = 16),
      panel.grid.major.x = ggplot2::element_line(colour = "grey75",
                                                 linetype = "dotted"),
      panel.grid.minor.x = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dotted"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dashed"),
      panel.grid.minor.y = ggplot2::element_blank()
    )
  } else {
    p <- p + base_theme + ggplot2::theme(
      axis.title.y       = ggplot2::element_text(
        size = 18, face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text.x        = ggplot2::element_text(angle = 315, vjust = 0.5,
                                                 hjust = 0, size = 16),
      axis.text.y        = ggplot2::element_text(size = 16),
      panel.grid.major.y = ggplot2::element_line(colour = "grey75",
                                                 linetype = "dotted"),
      panel.grid.minor.y = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dashed"),
      panel.grid.minor.x = ggplot2::element_blank()
    )
  }

  # --- Reference lines -----------------------------------------------------
  # Perpendicular to the value axis, so they flip with the orientation.
  ref_line <- function(v, ...) {
    if (horizontal) ggplot2::geom_vline(xintercept = v, ...)
    else            ggplot2::geom_hline(yintercept = v, ...)
  }

  if (!is.null(sig_value) && is.finite(sig_value)) {
    p <- p + ref_line(sig_value, linetype = "dotted", colour = "red",
                      linewidth = 0.6)
  }
  if (!absolute) {
    p <- p + ref_line(0, linetype = "solid", colour = "grey40",
                      linewidth = 0.5)
  }

  # --- Faceting ------------------------------------------------------------
  # Vertical: facets side by side, each as wide as its feature count.
  # Horizontal: facets stacked, each as tall as its feature count. `space`
  # is what keeps the row pitch identical across panels of unequal size.
  #
  # A single supergroup and a single species means the strip would restate
  # what the title already says once per panel — true of every single-family
  # figure (all 64 codons, all 20 amino acids: one group, one region), so
  # drop the faceting entirely rather than ship a decorative strip.
  n_supergroups <- nlevels(droplevels(result$supergroup))
  if (n_supergroups > 1 || has_species) {
    # Horizontal strips sit at the right of short, stacked panels, so a long
    # name ("Transcript architecture") cannot be set along the strip: it is
    # wrapped and written horizontally instead.
    strip_labels <- ggplot2::labeller(
      supergroup = function(s) {
        nm <- format_group_name(s, kind = "supergroup")
        if (horizontal) vapply(strwrap(nm, width = 13, simplify = FALSE),
                               paste, character(1), collapse = "\n")
        else nm
      },
      species    = ggplot2::label_value
    )
    if (horizontal) {
      p <- p + ggplot2::theme(strip.text.y = ggplot2::element_text(
        angle = 0, size = 16, face = "bold"))
      p <- p + if (has_species) {
        ggplot2::facet_grid(supergroup ~ species, scales = "free_y",
                            space = "free_y", labeller = strip_labels)
      } else {
        ggplot2::facet_grid(supergroup ~ ., scales = "free_y",
                            space = "free_y", labeller = strip_labels)
      }
    } else {
      p <- p + if (has_species) {
        ggplot2::facet_grid(species ~ supergroup, scales = "free_x",
                            space = "free_x", labeller = strip_labels)
      } else {
        ggplot2::facet_grid(~ supergroup, scales = "free_x",
                            space = "free_x", labeller = strip_labels)
      }
    }
  }

  # Likewise the region legend: with one region there is nothing to look up,
  # and a one-entry key invites the reader to hunt for a contrast that the
  # figure does not draw.
  if (nlevels(droplevels(result$region_f)) < 2) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # --- Return table --------------------------------------------------------
  table_out <- result |>
    dplyr::select(dplyr::any_of(c(
      "species", "variable", "group", "supergroup",
      "metric_stem", "metric_display", "region",
      "n", "correlation", "conf.low", "conf.high",
      "p_value", "q_value",
      "correlation_abs", "conf.low_abs", "conf.high_abs"
    )))

  # --- Diagnostics for the run report --------------------------------------
  # Everything a figure legend might need to state, computed once here rather
  # than re-derived from the CSV later. None of it is drawn on the plot.
  sparse <- table_out |>
    dplyr::filter(n < n_max) |>
    dplyr::distinct(variable, n) |>
    dplyr::arrange(n)

  below <- if (!is.null(sig_value) && is.finite(sig_value)) {
    table_out |>
      dplyr::filter(correlation_abs < sig_value) |>
      dplyr::arrange(correlation_abs) |>
      dplyr::select(variable, region, n, correlation, q_value)
  } else {
    table_out[0, c("variable", "region", "n", "correlation", "q_value")]
  }

  report <- list(
    response         = response,
    method           = method,
    orientation      = orientation,
    absolute         = absolute,
    conf             = conf,
    supergroups      = levels(droplevels(result$supergroup)),
    n_rows_df        = nrow(df),
    n_response       = sum(!is.na(df[[response]])),
    n_points         = nrow(table_out),
    n_features       = dplyr::n_distinct(table_out$variable),
    n_rows_plot      = nlevels(droplevels(result$metric_key)),
    n_min            = n_min,
    n_max            = n_max,
    n_median         = stats::median(table_out$n),
    sig_value        = sig_value,
    sig_auto         = sig_auto,
    sig_alpha        = sig_alpha,
    sig_at_n_min     = critical_correlation(n_min, method, sig_alpha),
    sig_at_n_max     = critical_correlation(n_max, method, sig_alpha),
    sparse           = sparse,
    ranking_applied  = ranking_applied,
    below_threshold  = below,
    n_sig_q          = sum(table_out$q_value < sig_alpha, na.rm = TRUE),
    dropped_columns  = unique(dropped),
    no_estimate      = no_estimate,
    top              = table_out |>
      dplyr::arrange(dplyr::desc(correlation_abs)) |>
      dplyr::slice_head(n = 10) |>
      dplyr::select(variable, region, n, correlation, q_value)
  )

  list(plot = p, table = table_out, report = report)
}


# -----------------------------------------------------------------------------
# Run report
# -----------------------------------------------------------------------------

#' Render one job's diagnostics as a markdown section.
#'
#' The counterpart to keeping the figures clean. Everything the subtitle no
#' longer says — cohort size, the UTR threshold behind it, per-feature n
#' range, what the significance line actually is and which features fall
#' short of it, which columns produced no estimate at all — lands here, in
#' one file, ready to be turned into figure legends.
#'
#' @param report  The `report` element of a feature_correlation_ranked() result.
#' @param title   Section heading.
#' @param species Species label for the cohort line.
#' @return Character vector of markdown lines.
format_run_report <- function(report, title, species = "human") {
  fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
  fmt_r <- function(x) if (is.null(x) || !is.finite(x)) "n/a" else sprintf("%.4f", x)

  lines <- c(
    sprintf("## %s", title),
    "",
    sprintf("- **Response**: %s (%s correlation, %d%% CI, %s)",
            report$response, report$method, round(report$conf * 100),
            if (report$absolute) "absolute |r|" else "signed r"),
    sprintf("- **Orientation**: %s", report$orientation),
    sprintf("- **Facets**: %s", paste(report$supergroups, collapse = ", ")),
    "",
    "### Cohort",
    "",
    sprintf("- Transcripts in the analysed frame (%s): **%s**",
            species, fmt_n(report$n_rows_df)),
    sprintf("- Non-missing `%s`: **%s**",
            report$response, fmt_n(report$n_response)),
    sprintf(paste0("- Cohort definition: transcripts with a 5'UTR under %d nt or a ",
                   "3'UTR under %d nt are excluded from every analysis ",
                   "(`MIN_5UTR_LENGTH` / `MIN_3UTR_LENGTH`, applied in ",
                   "`build_dataset()`)"),
            MIN_5UTR_LENGTH, MIN_3UTR_LENGTH),
    "",
    "### Features plotted",
    "",
    sprintf("- Points: **%s** across **%s** columns, **%s** axis rows",
            fmt_n(report$n_points), fmt_n(report$n_features),
            fmt_n(report$n_rows_plot)),
    sprintf("- Per-feature n: min **%s**, median **%s**, max **%s**",
            fmt_n(report$n_min), fmt_n(report$n_median),
            fmt_n(report$n_max)),
    ""
  )

  # Per-feature n range is the number most likely to be misread off the
  # figure, so name the offenders rather than just quoting a range.
  if (nrow(report$sparse) > 0) {
    lines <- c(lines,
      sprintf(paste0("**%d column(s) measured on fewer than the full %s ",
                     "transcripts.** These carry wider intervals and a ",
                     "higher significance threshold than the figure's ",
                     "reference line implies:"),
              nrow(report$sparse), fmt_n(report$n_max)),
      "",
      "| Column | n |",
      "| --- | ---: |",
      sprintf("| `%s` | %s |", report$sparse$variable,
              fmt_n(report$sparse$n)),
      "")
  } else {
    lines <- c(lines,
      sprintf("All plotted columns are measured on the same %s transcripts.",
              fmt_n(report$n_max)), "")
  }

  # Which high-cardinality families were trimmed to their top N. The stems come
  # from the data (select_features), so the same response gives the same
  # stems in every figure; a different response may give different ones.
  if (length(report$ranking_applied) > 0) {
    lines <- c(lines, "### Family selection", "")
    for (g in names(report$ranking_applied)) {
      info <- report$ranking_applied[[g]]
      lines <- c(lines, if (info$bound) {
        sprintf(paste0("- `%s`: top **%d** of %d stems, ranked by max |r| ",
                       "on this cohort."),
                g, info$requested, info$available)
      } else {
        sprintf(paste0("- `%s`: all **%d** stem(s) shown; the requested top-%d ",
                       "did not bind."),
                g, info$available, info$requested)
      })
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, "### Significance", "")

  if (is.null(report$sig_value)) {
    lines <- c(lines, "- No reference line drawn.", "")
  } else {
    lines <- c(lines,
      sprintf("- Reference line at |r| = **%s**%s",
              fmt_r(report$sig_value),
              if (!report$sig_auto) {
                " (fixed value supplied by the caller)"
              } else if (report$n_min < report$n_max) {
                # Only a floor when the features differ in n; say so, and give
                # the worst case so a legend can quote both ends.
                sprintf(paste0(", the critical value at alpha = %.2f for the ",
                               "largest n on the figure (%s). It is a FLOOR: ",
                               "a feature falling short of it is not ",
                               "significant at any n present. At the ",
                               "smallest n (%s) the threshold rises to |r| ",
                               "= %s."),
                        report$sig_alpha, fmt_n(report$n_max),
                        fmt_n(report$n_min), fmt_r(report$sig_at_n_min))
              } else {
                sprintf(paste0(", the critical value at alpha = %.2f. Every ",
                               "feature here is measured on the same %s ",
                               "transcripts, so this threshold is exact for ",
                               "all of them."),
                        report$sig_alpha, fmt_n(report$n_max))
              }),
      sprintf("- Points with BH q < %.2f: **%s** of %s",
              report$sig_alpha, fmt_n(report$n_sig_q),
              fmt_n(report$n_points)),
      "")

    if (nrow(report$below_threshold) > 0) {
      lines <- c(lines,
        sprintf("**%d point(s) below the reference line:**",
                nrow(report$below_threshold)),
        "",
        "| Column | Region | n | r | q |",
        "| --- | --- | ---: | ---: | ---: |",
        sprintf("| `%s` | %s | %s | %.4f | %s |",
                report$below_threshold$variable,
                report$below_threshold$region,
                fmt_n(report$below_threshold$n),
                report$below_threshold$correlation,
                ifelse(is.na(report$below_threshold$q_value), "n/a",
                       sprintf("%.3g", report$below_threshold$q_value))),
        "")
    } else {
      lines <- c(lines, "Every plotted point clears the reference line.", "")
    }
  }

  # A column that resolved to a valid feature name but yielded no estimate is
  # the quiet failure mode under a row filter — it simply is not on the
  # figure, with nothing to see. Always state it, including when it is zero.
  lines <- c(lines, "### Exclusions", "")
  if (length(report$no_estimate) > 0) {
    lines <- c(lines,
      sprintf(paste0("**%d plottable column(s) produced no correlation** ",
                     "(fewer than the minimum non-missing pairs against ",
                     "`%s`, or no remaining variance):"),
              length(report$no_estimate), report$response),
      "",
      paste0("- `", report$no_estimate, "`"),
      "")
  } else {
    lines <- c(lines,
      "- Every plottable column produced a correlation.", "")
  }
  if (length(report$dropped_columns) > 0) {
    lines <- c(lines,
      sprintf(paste0("%d column(s) carried no region token and cannot sit ",
                     "on a region-dodged axis: %s"),
              length(report$dropped_columns),
              paste0("`", report$dropped_columns, "`", collapse = ", ")),
      "")
  }

  lines <- c(lines,
    "### Strongest correlations", "",
    "| Column | Region | n | r | q |",
    "| --- | --- | ---: | ---: | ---: |",
    sprintf("| `%s` | %s | %s | %.4f | %s |",
            report$top$variable, report$top$region, fmt_n(report$top$n),
            report$top$correlation,
            ifelse(is.na(report$top$q_value), "n/a",
                   sprintf("%.3g", report$top$q_value))),
    "")

  lines
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  species <- "human"

  # Filtered cohort, explicitly. build_dataset() applies the UTR length filters to the
  # frame it returns by default, so this is what you get anyway — naming it
  # and asserting it means a future change to that default cannot silently
  # move these figures onto the unfiltered table.
  df <- build_dataset(species, min_5utr = MIN_5UTR_LENGTH, min_3utr = MIN_3UTR_LENGTH)
  stopifnot(
    all(df$length_5utr >= MIN_5UTR_LENGTH, na.rm = TRUE),
    all(df$length_3utr >= MIN_3UTR_LENGTH, na.rm = TRUE),
    !anyNA(df$length_5utr), !anyNA(df$length_3utr)
  )
  message("Cohort: ", nrow(df), " transcripts (5'UTR >= ",
          MIN_5UTR_LENGTH, " nt, 3'UTR >= ", MIN_3UTR_LENGTH, " nt)")

  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "reports"),
             showWarnings = FALSE, recursive = TRUE)

  save_plot_jpg_pdf <- function(plot, filename_base, width, height,
                                units = "mm", dpi = 300, limitsize = TRUE) {
    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots", paste0(filename_base, ".jpg")),
      plot = plot, width = width, height = height, units = units, dpi = dpi,
      limitsize = limitsize
    )
    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots", paste0(filename_base, ".pdf")),
      plot = plot, width = width, height = height, units = units,
      limitsize = limitsize, device = grDevices::cairo_pdf
    )
  }

  # Horizontal figures grow downward, so height follows the feature count
  # instead of being fixed per job. `rows` is the number of discrete axis
  # ticks; `pitch` is millimetres per tick, larger when several regions dodge
  # within one tick.
  auto_height <- function(rows, pitch = 8, base = 70, min_h = 140) {
    max(min_h, base + rows * pitch)
  }

  # Four figures.
  #
  # The two broad ones use the core set, one facet per supergroup. Codons and
  # amino acids are trimmed by `top_n`, which ranks the whole family against
  # the figure's own response, so the stems shown are the genuine top two and
  # agree with the full-family figures below. Derived per response: half-life
  # and translation efficiency may show different codons, by design.
  #
  # The codon and amino-acid figures include one family each and take no top-N
  # cut: the whole family is the point. Both families are CDS-only, so there
  # is one region per row and no dodging, and a single supergroup, so no
  # facet strip.
  jobs <- list(
    list(response = "halflife",
         suffix   = "halflife",
         title    = "Half-life, core features",
         include  = "core",
         top_n    = list(codon_freqs = 2, aa_freqs = 2),
         width    = 260,
         pitch    = 9),
    list(response = "translation_efficiency",
         suffix   = "translation_efficiency",
         title    = "Translation efficiency, core features",
         include  = "core",
         top_n    = list(codon_freqs = 2, aa_freqs = 2),
         width    = 260,
         pitch    = 9),
    list(response = "halflife",
         suffix   = "halflife_aa_full",
         title    = "Half-life, all amino acids",
         include  = "aa_freqs",
         top_n    = NULL,
         width    = 240,
         pitch    = 8),
    list(response = "halflife",
         suffix   = "halflife_codon_full",
         title    = "Half-life, all codons",
         include  = "codon_freqs",
         top_n    = NULL,
         width    = 240,
         pitch    = 6)
  )

  report_lines <- c(
    "# Ranked feature correlation run report",
    "",
    sprintf("Generated %s from `analysis/correlations/feature_correlation_ranked.R`.",
            format(Sys.Date())),
    "",
    paste0("Figure subtitles are deliberately minimal. Every cohort number, ",
           "per-feature n and significance threshold behind the figures in ",
           "`", file.path(OUTPUT_DIR, "plots"), "/feature_correlation_ranked_*` ",
           "is recorded below, for the figure legends to be written from."),
    ""
  )

  for (job in jobs) {
    if (!job$response %in% names(df)) {
      message("Skipping: ", job$response, " not in dataset")
      next
    }

    message("\nRanked plot: ", job$title)

    out <- feature_correlation_ranked(
      df,
      response         = job$response,
      include          = job$include,
      orientation      = "horizontal",
      sig_threshold    = "auto",
      top_n            = job$top_n
    )

    height <- auto_height(out$report$n_rows_plot, pitch = job$pitch)
    message("  ", out$report$n_rows_plot, " rows -> ", height, " mm tall")

    save_plot_jpg_pdf(
      plot          = out$plot,
      filename_base = paste0("feature_correlation_ranked_", job$suffix),
      width         = job$width,
      height        = height,
      limitsize     = FALSE
    )

    write.csv(
      out$table,
      file.path(OUTPUT_DIR, "tables",
                paste0("feature_correlation_ranked_", job$suffix, ".csv")),
      row.names = FALSE
    )

    report_lines <- c(report_lines,
                      format_run_report(out$report, job$title, species))
  }

  report_path <- file.path(OUTPUT_DIR, "reports",
                           "feature_correlation_ranked_report.md")
  writeLines(report_lines, report_path)

  message("\nRanked plots complete:")
  message("  ", file.path(OUTPUT_DIR, "plots"),
          "/feature_correlation_ranked_*.{jpg,pdf}")
  message("  ", file.path(OUTPUT_DIR, "tables"),
          "/feature_correlation_ranked_*.csv")
  message("  ", report_path)
}
