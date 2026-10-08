# =============================================================================
# Banded feature correlation plot
# =============================================================================
# One panel, feature names down the y axis, correlation along the x axis. The
# supergroups are not separate facets: they are contiguous coloured bands
# behind the rows, each labelled on a bar at the right edge. That keeps a
# single shared x axis (rows are directly comparable across supergroups) while
# still showing which family a feature belongs to.
#
# Sibling of feature_correlation_ranked.R, not a replacement. The correlations,
# confidence intervals, BH q-values and significance line are computed by
# feature_correlation_ranked() itself, so the two figures cannot disagree on a
# number; this file only re-draws its returned table. The supergroups, their
# order, their display names and their colours all come from
# R/feature_table.csv (SUPERGROUPS, SUPERGROUP_DISPLAY_NAMES, FEATURE_TABLE) —
# nothing about the feature hierarchy is listed here.
#
# ORDER. Bands are ranked by their strongest member: the band holding the
# feature with the largest |r| (over its regions) is on top. Within a band the
# strongest feature is on top too. Band order therefore depends on the
# response, and ties fall back to the feature table's supergroup order.
#
# MATCHING ANOTHER PLOT. To flip between two figures (halflife vs TE), pass the
# first result as `order_from` and pick how much of its order to copy:
#   order_follow = "supergroup"  - bands in the reference's order, features
#                                  still ranked by this response's own |r|;
#   order_follow = "feature"     - bands AND the features inside them in the
#                                  reference's order.
# Anything the reference does not have (e.g. a codon that made the top_n cut
# for one response only) goes after the reference's rows in its band, by this
# response's own |r|; a band the reference lacks goes after its bands.
#
# BAND COLOURS. Set in R/colour_config.R (SUPERGROUP_COLOURS, and the tint
# strengths SUPERGROUP_BAND_ALPHA / SUPERGROUP_LABEL_ALPHA). The tints are pale
# on purpose — the points carry the region colours (the supergroup colours are
# vivid; the mixes with white make them subtle).
#
# LEGEND. Drawn inside the panel, over the emptiest stretch of the right-hand
# side (see choose_legend_position()), so the figure has no dead margin.
#
# COMPARABLE FIGURES. To flip between two figures, give them the same x range
# (`x_max`), the same drawn panel width (`panel_width_mm`) and the same legend
# spot. `legend_from` takes the other figures' results (same rows, same order)
# and picks the one window where the legend clears the points of all of them.
# The runner below does this in two passes: draw, take the common x maximum,
# draw again.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# feature_correlation_ranked.R runs its own jobs when sourced at top level.
# Loading it into a private environment (where the runner's guard is false)
# gives us its functions without regenerating its figures.
.ranked_env <- new.env(parent = globalenv())
sys.source("analysis/correlations/feature_correlation_ranked.R",
           envir = .ranked_env)


# -----------------------------------------------------------------------------
# Colour helpers
# -----------------------------------------------------------------------------

#' Mix colours with white. `amount` is the share of the original colour kept.
mix_with_white <- function(col, amount) {
  rgb_mat <- grDevices::col2rgb(col) / 255
  mixed   <- amount * rgb_mat + (1 - amount)
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

#' Where to put the legend inside the panel.
#'
#' Slides a window as tall as the legend down the rows and keeps the one whose
#' rows reach least far to the right, preferring windows that sit inside a
#' single band (a legend straddling two bands reads as belonging to neither).
#' Ties go to the lowest window. Returns the vertical position in npc units
#' of the panel, plus the rightmost data extent under the chosen window so the
#' caller can warn when the legend is likely to cover points.
choose_legend_position <- function(rows, tab, bands, legend_rows,
                                   other_reach = NULL) {
  ys    <- sort(rows$y)                           # bottom -> top
  n_win <- min(ceiling(legend_rows), length(ys))
  band_of <- stats::setNames(as.character(rows$supergroup), rows$y)
  reach <- tapply(tab$.hi, tab$y, max, na.rm = TRUE)
  # Other figures' reach for the same rows (named by row_key): the legend must
  # clear their points too.
  if (!is.null(other_reach)) {
    key   <- rows$row_key[match(names(reach), as.character(rows$y))]
    extra <- other_reach[key]
    reach <- pmax(reach, ifelse(is.na(extra), -Inf, extra))
  }

  starts <- seq_len(length(ys) - n_win + 1)
  cand <- do.call(rbind, lapply(starts, function(i) {
    w <- ys[i:(i + n_win - 1)]
    data.frame(lo = w[1], hi = w[n_win],
               one_band = length(unique(band_of[as.character(w)])) == 1,
               reach = max(reach[as.character(w)]))
  }))
  pool <- if (any(cand$one_band)) cand[cand$one_band, ] else cand
  best <- pool[order(pool$reach, pool$lo), ][1, ]
  ymin <- min(bands$ymin); ymax <- max(bands$ymax)
  list(y_npc = ((best$lo + best$hi) / 2 - ymin) / (ymax - ymin),
       reach = best$reach)
}


# -----------------------------------------------------------------------------
# Main plot function
# -----------------------------------------------------------------------------

#' Feature correlations with a response, one panel, supergroups as bands.
#'
#' Arguments other than the styling ones are passed straight to
#' feature_correlation_ranked(), which owns the selection, the correlation
#' method and the significance line.
#'
#' @param df           Dataframe from build_dataset().
#' @param response     Response column (default "halflife").
#' @param absolute     Plot |r| (default) or signed r.
#' @param band_alpha   Share of the supergroup colour kept in the band tint
#'                     (0 = white, 1 = full colour). Default from
#'                     R/colour_config.R.
#' @param label_alpha  The same for the label bar.
#' @param legend_position "auto" (default) puts the legend inside the panel over
#'                     the emptiest right-hand stretch; "right" puts it outside;
#'                     c(x, y) places it at those panel (npc) coordinates.
#' @param legend_from  Other feature_correlation_bands() results drawn with the
#'                     same rows in the same order (a result, or a list of
#'                     them). The legend goes where it clears the points of
#'                     this figure and of all of them, so every figure gets the
#'                     same spot. Default NULL: this figure alone decides.
#' @param x_max        Upper end of the x axis. Default NULL: just past the
#'                     largest confidence limit. Give the larger of several
#'                     figures' `report$x_max` to put them on one scale.
#' @param panel_width_mm Drawn width of the panel in mm (needs ggplot2 >= 4.0).
#'                     Default NULL: whatever the y labels leave over. Set the
#'                     same value on figures to be compared, so equal data
#'                     lengths are equal on paper.
#' @param label_size   Text size (mm) of the supergroup label bars.
#' @param row_fill     Fraction of a row's height the dodged regions span.
#' @param row_mm       Approximate millimetres per row on the saved figure;
#'                     used only to fit the supergroup labels into their bands.
#' @param order_from   A previous feature_correlation_bands() result whose band
#'                     (and optionally feature) order to reuse. Default NULL:
#'                     order by this response's own |r|.
#' @param order_follow With `order_from`: "supergroup" copies only the band
#'                     order, "feature" copies the feature order as well.
#' @param ...          Passed to feature_correlation_ranked(): method,
#'                     include, exclude, regions, top_n, min_abs_correlation,
#'                     sig_threshold, sig_alpha, conf, min_n.
#' @return list(plot, table, report), as feature_correlation_ranked() but with
#'   `table` gaining `row_order` and the report gaining `bands`.
#' @export
feature_correlation_bands <- function(df,
                                      response   = "halflife",
                                      absolute   = TRUE,
                                      band_alpha = SUPERGROUP_BAND_ALPHA,
                                      label_alpha = SUPERGROUP_LABEL_ALPHA,
                                      label_size = 4.2,
                                      legend_position = "auto",
                                      legend_from = NULL,
                                      x_max      = NULL,
                                      panel_width_mm = NULL,
                                      row_fill   = 0.75,
                                      row_mm     = 9,
                                      order_from = NULL,
                                      order_follow = c("supergroup", "feature"),
                                      ...) {
  order_follow <- match.arg(order_follow)
  ref_rows <- if (!is.null(order_from)) order_from$report$row_order
  if (!is.null(order_from) && is.null(ref_rows)) {
    stop("order_from must be a feature_correlation_bands() result")
  }

  out <- .ranked_env$feature_correlation_ranked(
    df,
    response         = response,
    absolute         = absolute,
    orientation      = "horizontal",
    ...
  )
  tab <- out$table
  if (length(unique(tab$species)) > 1) {
    stop("feature_correlation_bands draws one species at a time")
  }

  tab$.value <- if (absolute) tab$correlation_abs else tab$correlation
  tab$.lo    <- if (absolute) tab$conf.low_abs    else tab$conf.low
  tab$.hi    <- if (absolute) tab$conf.high_abs   else tab$conf.high

  # --- Row order: bands by their top member's |r|, then rows by max |r| -----
  # order() is stable, so equal maxima keep the feature table's order.
  sg_table  <- intersect(names(SUPERGROUPS), unique(tab$supergroup))
  sg_max    <- vapply(sg_table, function(g)
    max(tab$correlation_abs[tab$supergroup == g], na.rm = TRUE), numeric(1))
  sg_levels <- sg_table[order(-sg_max)]
  if (!is.null(ref_rows)) {
    ref_sg <- unique(as.character(ref_rows$supergroup))
    # Reference bands first, in its order; any extra bands after, by own |r|.
    sg_levels <- c(intersect(ref_sg, sg_levels), setdiff(sg_levels, ref_sg))
  }
  tab$supergroup <- factor(tab$supergroup, levels = sg_levels)
  tab$row_key    <- paste(tab$supergroup, tab$metric_display, sep = "::")

  rows <- tab |>
    dplyr::group_by(supergroup, row_key, metric_display) |>
    dplyr::summarise(max_r = max(correlation_abs, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(supergroup, dplyr::desc(max_r))
  if (!is.null(ref_rows) && order_follow == "feature") {
    ref_key <- paste(ref_rows$supergroup, ref_rows$metric_display, sep = "::")
    ref_pos <- match(rows$row_key, ref_key)
    rows <- rows[order(rows$supergroup,
                       is.na(ref_pos), ref_pos, -rows$max_r), ]
  }

  # A discrete y axis draws its first level at the bottom; reverse so the
  # first row of the first band is at the top.
  rows$y <- rev(seq_len(nrow(rows)))
  tab$y  <- rows$y[match(tab$row_key, rows$row_key)]
  tab$row_order <- match(tab$row_key, rows$row_key)

  # Manual dodge. The y axis is numeric (so bands and labels share it), and
  # position_dodge() needs a discrete axis, so offset each region within its
  # row by hand: regions keep a fixed order, centred on the row.
  region_order <- intersect(names(REGION_COLOURS), unique(tab$region))
  tab <- tab |>
    dplyr::group_by(row_key) |>
    dplyr::mutate(
      .k    = dplyr::n(),
      .rank = match(region, region_order),
      .slot = rank(.rank, ties.method = "first"),
      y_dodged = y + (((.k + 1) / 2 - .slot)) * (row_fill / pmax(.k, 1))
    ) |>
    dplyr::ungroup()

  # Position-based y (numeric) so bands, dodging and labels share one scale.
  y_labels <- stats::setNames(rows$metric_display, rows$y)

  # --- Bands ----------------------------------------------------------------
  bands <- rows |>
    dplyr::group_by(supergroup) |>
    dplyr::summarise(n_rows = dplyr::n(),
                     ymin = min(y) - 0.5, ymax = max(y) + 0.5,
                     .groups = "drop") |>
    dplyr::mutate(
      base  = supergroup_colour(as.character(supergroup)),
      tint  = mix_with_white(base, band_alpha),
      strong = mix_with_white(base, label_alpha),
      label = format_group_name(as.character(supergroup), kind = "supergroup")
    )

  # --- x extent and the label bar to the right of the panel -----------------
  xlo <- if (absolute) 0 else min(tab$.lo, na.rm = TRUE)
  xhi <- max(tab$.hi, na.rm = TRUE)
  pad <- 0.03 * (xhi - xlo)
  xlo <- if (absolute) 0 else xlo - pad
  xhi <- xhi + pad
  if (!is.null(x_max)) xhi <- max(xhi, x_max)
  bar_w <- 0.09 * (xhi - xlo)
  x_breaks <- round(seq(ceiling(xlo * 10 - 1e-9) / 10,
                        floor(xhi * 10 + 1e-9) / 10, by = 0.1), 1)
  x_guides <- if (absolute) setdiff(x_breaks, 0) else x_breaks  # 0 is the edge

  # Fit each label into its band once rotated: wrap to the band's length,
  # then shrink the text if the longest unbreakable word still overruns.
  char_mm <- 0.5                                  # bold char width per size unit
  band_mm <- bands$n_rows * row_mm
  bands$label_wrapped <- mapply(function(lab, mm) {
    paste(strwrap(lab, width = max(6, floor(mm / (label_size * char_mm)))),
          collapse = "\n")
  }, bands$label, band_mm)
  longest <- vapply(strsplit(bands$label_wrapped, "\n", fixed = TRUE),
                    function(l) max(nchar(l)), numeric(1))
  # One size for every label (the smallest that fits) so the bars look alike.
  bands$label_size <- min(label_size, 0.85 * band_mm / (longest * char_mm))
  # Dark bars get white text, light bars black.
  lum <- colSums(grDevices::col2rgb(bands$strong) * c(0.299, 0.587, 0.114))
  bands$text_col <- ifelse(lum < 140, "white", "black")

  sig_value <- out$report$sig_value
  region_colours <- REGION_COLOURS
  region_shapes  <- REGION_SHAPES
  tab$region_f <- factor(tab$region,
                         levels = intersect(names(region_colours),
                                            unique(tab$region)))

  value_lab <- sprintf("%s correlation with %s",
                       if (absolute) "Absolute" else "Signed",
                       format_col_name(response))

  p <- ggplot2::ggplot(tab) +
    # Bands first so they sit under everything.
    ggplot2::geom_rect(
      data = bands,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
                   fill = I(tint)),
      inherit.aes = FALSE
    ) +
    # White seam at every internal band boundary (all but the panel's bottom
    # edge, which is the lowest band's ymin).
    ggplot2::geom_hline(yintercept = setdiff(bands$ymin, min(bands$ymin)),
                        colour = "white", linewidth = 0.8) +
    # Vertical dashed guides every 0.1, over the bands.
    ggplot2::geom_vline(xintercept = x_guides, linetype = "dashed",
                        colour = "grey45", linewidth = 0.3, alpha = 0.6) +
    # Dashed separators between neighbouring feature rows, including across
    # band edges (over the white seams). The panel's bottom edge is skipped.
    ggplot2::geom_hline(yintercept = setdiff(rows$y - 0.5, min(bands$ymin)),
                        linetype = "dashed",
                        colour = "grey45", linewidth = 0.3, alpha = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(y = y_dodged, xmin = .lo, xmax = .hi, colour = region_f,
                   ),
      orientation = "y", width = 0.2, linewidth = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .value, y = y_dodged, colour = region_f, shape = region_f,
                   group = region_f),
      size = 2.5, stroke = 0.5
    ) +
    # Label bars sit just outside the panel; clip = "off" lets them draw there.
    ggplot2::geom_rect(
      data = bands,
      ggplot2::aes(xmin = xhi, xmax = xhi + bar_w, ymin = ymin, ymax = ymax,
                   fill = I(strong)),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = bands,
      ggplot2::aes(x = xhi + bar_w / 2, y = (ymin + ymax) / 2,
                   label = label_wrapped, size = label_size,
                   colour = I(text_col)),
      angle = 270, fontface = "bold", lineheight = 0.85,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_y_continuous(
      breaks = rows$y, labels = y_labels,
      expand = ggplot2::expansion(add = 0)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      expand = ggplot2::expansion(0)
    ) +
    ggplot2::coord_cartesian(xlim = c(xlo, xhi),
                             ylim = c(min(bands$ymin), max(bands$ymax)),
                             clip = "off") +
    ggplot2::scale_colour_manual(
      values = region_colours, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r], name = "Region"
    ) +
    ggplot2::scale_shape_manual(
      values = region_shapes, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r], name = "Region"
    ) +
    ggplot2::labs(
      x = value_lab, y = NULL,
      title    = sprintf("Correlation with %s",
                         format_col_name(response)),
      subtitle = sprintf("%d%% CI; within each supergroup ordered by |r|",
                         round(out$report$conf * 100))
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 16, colour = "grey30"),
      axis.title.x     = ggplot2::element_text(size = 18, face = "bold",
                                               margin = ggplot2::margin(t = 10)),
      axis.text.y      = ggplot2::element_text(size = 14, hjust = 1,
                                               colour = "black"),
      axis.text.x      = ggplot2::element_text(size = 16),
      legend.title     = ggplot2::element_text(size = 18, face = "bold"),
      legend.text      = ggplot2::element_text(size = 16),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "grey40", fill = NA,
                                               linewidth = 0.5),
      # Room on the right for the label bars, which draw outside the panel.
      plot.margin      = ggplot2::margin(10, 24, 10, 10, unit = "mm")
    )

  if (!is.null(sig_value) && is.finite(sig_value)) {
    p <- p + ggplot2::geom_vline(xintercept = sig_value, linetype = "dotted",
                                 colour = "red", linewidth = 0.6)
  }
  if (!absolute) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = "grey40",
                                 linewidth = 0.5)
  }
  n_regions <- nlevels(droplevels(tab$region_f))
  if (n_regions < 2) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else if (identical(legend_position, "right")) {
    p <- p + ggplot2::theme(legend.position = "right")
  } else {
    if (is.numeric(legend_position)) {
      pos <- list(x = legend_position[1], y = legend_position[2])
    } else {
      # Legend height in rows: title (~10 mm) plus ~7.5 mm per entry.
      other <- NULL
      if (!is.null(legend_from)) {
        if (!is.null(legend_from$report)) legend_from <- list(legend_from)
        other <- Reduce(function(a, b) {
          k <- union(names(a), names(b))
          pmax(a[k], b[k], na.rm = TRUE)
        }, lapply(legend_from, function(r) r$report$row_reach))
      }
      lg <- choose_legend_position(rows, tab, bands,
                                   (10 + n_regions * 7.5) / row_mm,
                                   other_reach = other)
      pos <- list(x = 0.995, y = lg$y_npc)
      # The legend is roughly a quarter of the panel wide.
      if (lg$reach > xhi - 0.26 * (xhi - xlo)) {
        message("feature_correlation_bands: no clear stretch for the legend; ",
                "it may cover points. Pass legend_position = c(x, y) or ",
                "\"right\".")
      }
    }
    p <- p + ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(pos$x, pos$y),
      legend.justification.inside = c(1, 0.5),
      legend.background = ggplot2::element_rect(
        fill = grDevices::adjustcolor("white", 0.9), colour = "grey40",
        linewidth = 0.4),
      legend.margin = ggplot2::margin(6, 8, 6, 8),
      legend.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.text  = ggplot2::element_text(size = 14),
      legend.key.height = ggplot2::unit(6, "mm")
    )
  }

  if (!is.null(panel_width_mm)) {
    p <- p + ggplot2::theme(panel.widths = ggplot2::unit(panel_width_mm, "mm"))
  }

  out$plot   <- p
  out$table  <- tab |>
    dplyr::select(-dplyr::any_of(c(".value", ".lo", ".hi", "row_key", "y",
                                   "region_f")))
  out$report$bands  <- as.character(bands$supergroup)
  out$report$row_order <- data.frame(
    supergroup     = as.character(rows$supergroup),
    metric_display = rows$metric_display)
  out$report$n_rows_plot <- nrow(rows)
  out$report$x_max <- xhi
  # Furthest right any point or CI reaches in each row, for legend_from.
  out$report$row_reach <- tapply(tab$.hi, tab$row_key, max, na.rm = TRUE)
  out
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  species <- "human"

  df <- build_dataset(species, min_5utr = MIN_5UTR_LENGTH, min_3utr = MIN_3UTR_LENGTH)
  message("Cohort: ", nrow(df), " transcripts (5'UTR >= ",
          MIN_5UTR_LENGTH, " nt, 3'UTR >= ", MIN_3UTR_LENGTH, " nt)")

  dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  # The core set, minus experimental probing (out of "core"), with codons and
  # amino acids trimmed to their top two against each figure's own response,
  # exactly as in the ranked figures. The Structure band therefore holds the
  # predicted folding features only.
  top_n <- list(codon_freqs = 2, aa_freqs = 2)

  # Millimetres per row. Also tells the label fitter how tall a band is.
  row_mm <- 14

  jobs <- list(
    list(response = "halflife",               suffix = "halflife"),
    list(response = "translation_efficiency", suffix = "translation_efficiency")
  )

  # Draw one job and save it under `suffix`.
  save_banded <- function(res, suffix) {
    height <- max(140, 70 + res$report$n_rows_plot * row_mm)
    base   <- file.path(OUTPUT_DIR, "plots",
                        paste0("feature_correlation_bands_", suffix))
    ggplot2::ggsave(paste0(base, ".jpg"), res$plot, width = 260,
                    height = height, units = "mm", dpi = 300,
                    limitsize = FALSE)
    ggplot2::ggsave(paste0(base, ".pdf"), res$plot, width = 260,
                    height = height, units = "mm", limitsize = FALSE,
                    device = grDevices::cairo_pdf)
    write.csv(res$table,
              file.path(OUTPUT_DIR, "tables",
                        paste0("feature_correlation_bands_", suffix, ".csv")),
              row.names = FALSE)
    message("  ", res$report$n_rows_plot, " rows, bands: ",
            paste(res$report$bands, collapse = ", "))
  }

  # When TRUE, the codons and amino acids shown are the union of each
  # response's own top-N, so a figure ordered like the reference carries
  # exactly the reference's rows. The reference figure itself uses the union
  # too; a later job's standalone figure keeps its own top-N. FALSE gives the
  # strict top-N everywhere (the feature-order figure then appends rows the
  # reference lacks at the end of their band).
  match_features <- TRUE

  # The trimmed families, as columns kept for a response under `top_n`.
  trimmed_cols <- function(response) {
    sel <- select_features(df, "core", top_n = top_n, response = response)
    sel$column[sel$feature_id %in% names(top_n)]
  }
  # `df` without the trimmed-family columns that no response's top-N kept.
  df_union <- function(responses) {
    fam  <- select_features(df, "core")
    fam  <- fam$column[fam$feature_id %in% names(top_n)]
    keep <- unique(unlist(lapply(responses, trimmed_cols)))
    df[, setdiff(names(df), setdiff(fam, keep)), drop = FALSE]
  }

  jobs <- Filter(function(j) j$response %in% names(df), jobs)
  df_match <- if (match_features) df_union(vapply(jobs, `[[`, "", "response"))

  # Draw one job. `d` is the data; `trim` the top_n it applies.
  draw <- function(job, d = df, trim = top_n, ...) {
    feature_correlation_bands(d, response = job$response, include = "core",
                              row_mm = row_mm, top_n = trim, ...)
  }

  # The first job is the reference. Every later job is drawn in its own |r|
  # order and again ordered like the reference (bands and features), so the
  # two figures can be flipped between. Returns list(figs, linked): `figs` is
  # every figure with its file suffix, `linked` those sharing the reference's
  # rows (the reference and the feature-ordered ones).
  #
  # Two passes make the figures comparable: the first finds the common x
  # maximum and the rows' reach; the second redraws on that x range with one
  # legend spot that clears every linked figure's points. Every figure gets
  # the same drawn panel width.
  draw_all <- function(x_max = NULL, legend_from = NULL) {
    reference <- NULL
    figs <- list(); linked <- list()
    common <- list(x_max = x_max, panel_width_mm = panel_width_mm)
    for (job in jobs) {
      message("\nBanded plot: ", job$response)
      if (is.null(reference)) {
        reference <- do.call(draw, c(list(job, if (match_features) df_match else df,
                                          if (match_features) NULL else top_n,
                                          legend_from = legend_from), common))
        ref_name  <- job$suffix
        figs[[job$suffix]] <- reference
        linked <- c(linked, list(reference))
      } else {
        figs[[job$suffix]] <- do.call(draw, c(list(job), common))
        message("  ordered like ", ref_name)
        matched <- if (match_features)
          do.call(draw, c(list(job, df_match, NULL, order_from = reference,
                               order_follow = "feature",
                               legend_from = legend_from), common))
        else
          do.call(draw, c(list(job, order_from = reference,
                               order_follow = "feature",
                               legend_from = legend_from), common))
        nm <- paste0(job$suffix, "_", ref_name, "_feature_order")
        figs[[nm]] <- matched
        linked <- c(linked, list(matched))
      }
    }
    list(figs = figs, linked = linked)
  }

  # Drawn panel width (mm) for every figure. Must leave room, on the 260 mm
  # page, for the y labels, the label bar and the margins.
  panel_width_mm <- 150

  first <- draw_all()
  x_max <- max(vapply(first$figs, function(r) r$report$x_max, numeric(1)))
  message("\nCommon x maximum: ", signif(x_max, 3))
  final <- draw_all(x_max = x_max, legend_from = first$linked)
  for (nm in names(final$figs)) save_banded(final$figs[[nm]], nm)
}
