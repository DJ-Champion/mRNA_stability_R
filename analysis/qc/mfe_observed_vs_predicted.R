# =============================================================================
# Observed vs predicted MFE, by region
# =============================================================================
# Predicted MFE (mfe_expected_*, "MFE.pred") is RNAexpected's average MFE of a
# random sequence with the same GC and length. Three figures, one facet per
# region (5' UTR, CDS, 3' UTR, mRNA):
#
#   1. Relative difference (observed - predicted) / |predicted| against length
#      (log x), hex-binned.
#   2. Predicted against observed MFE, zoomed to the central 95% of each region.
#   3. Predicted MFE against the per-transcript dinucleotide-shuffle median
#      (rnafold_median_*, 1,000 esl-shuffles), same layout as 2: how far a
#      GC-and-length-only baseline is from one that also keeps dinucleotides.
#
# Figures carry no subtitles or in-plot notes. Everything a legend needs
# (definitions, counts off scale, zoom windows, model limits, provenance) is
# written to reports/mfe_observed_vs_predicted_report.md.
#
# Known limits of the model to read them against:
#   < 100 nt     the model under-predicts structure (about -13% at 50-100 nt,
#                worse below) and returns 0 below ~30-60 nt
#   > 4,000 nt   outside the 20-4,000 nt fitting range: extrapolated
# Figure 1 marks both with dashed lines.
#
# Relative difference, not raw: random-sequence SD grows from ~2 kcal/mol at
# 50 nt to ~14 at 4 kb, so raw differences are not comparable across lengths.
# Rows with predicted MFE = 0 have no relative difference and are counted, not
# plotted.
#
# COHORT. The analysed cohort (build_dataset() defaults: 5' UTR >= 30 nt,
# 3' UTR >= 100 nt), the same frame the models see.
#
# Run after build_dataset() has populated the cache.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

species     <- "human"
MFE_REGIONS <- c("5utr", "cds", "3utr", "mrna")
SHORT_NT    <- 100
LONG_NT     <- 4000
REL_YLIM    <- c(-0.75, 0.75)   # shared across regions
HEX_BINS    <- 60
ZOOM_KEEP   <- 0.95             # figures 2-3: central share of each region shown

df <- build_dataset(species, min_5utr = MIN_5UTR_LENGTH, min_3utr = MIN_3UTR_LENGTH)


# --- Long form: one row per transcript x region --------------------------------
mfe_long <- purrr::map_dfr(MFE_REGIONS, function(r) {
  tibble::tibble(
    region    = r,
    observed  = df[[paste0("rnafold_score_",  r)]],
    predicted = df[[paste0("mfe_expected_",   r)]],
    shuffled  = df[[paste0("rnafold_median_", r)]],
    length    = df[[paste0("length_",         r)]]
  )
}) |>
  filter(!is.na(observed), !is.na(predicted), !is.na(length)) |>
  mutate(
    region       = factor(REGION_DISPLAYS[region], levels = REGION_DISPLAYS[MFE_REGIONS]),
    rel_observed = ifelse(predicted < 0, (observed - predicted) / abs(predicted), NA_real_),
    rel_shuffled = ifelse(predicted < 0, (shuffled - predicted) / abs(predicted), NA_real_)
  )


theme_mfe <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 17, face = "bold"),
    strip.text       = element_text(size = 13, face = "bold", hjust = 0),
    panel.background = element_rect(fill = "grey96", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "grey45", linewidth = 0.6),
    panel.grid.major = element_line(colour = "white", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.2, "lines"),
    legend.position  = "right"
  )

# Hex fill for all three figures: one blue, light to dark, log counts.
density_fill <- scale_fill_gradient(low = "#C6DBEF", high = "#08306B", trans = "log10",
                                    labels = scales::label_comma(),
                                    name = "Transcript\ndensity")


# --- Figure 1: relative difference vs length -----------------------------------
rel_on_scale <- filter(mfe_long, !is.na(rel_observed),
                       rel_observed >= REL_YLIM[1], rel_observed <= REL_YLIM[2])

# Hexes are drawn from the rows inside the y range only, so the colour scale
# describes what is on screen; the report counts the rows left out.
plot_rel_vs_length <- ggplot(rel_on_scale, aes(x = length, y = rel_observed)) +
  geom_hex(bins = HEX_BINS) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = c(SHORT_NT, LONG_NT), colour = "grey20",
             linewidth = 0.4, linetype = "22") +
  facet_wrap(~ region, ncol = 2) +
  scale_x_log10(labels = scales::label_comma(),
                breaks = c(30, 100, 300, 1000, 4000, 30000)) +
  scale_y_continuous(labels = scales::label_percent()) +
  coord_cartesian(ylim = REL_YLIM) +
  density_fill +
  labs(
    title = "Observed MFE relative to predicted MFE, by region length",
    x     = "Region length (nt, log scale)",
    y     = "Relative difference from predicted MFE"
  ) +
  theme_mfe


# --- Figures 2-3: predicted against a measured MFE -----------------------------
#' Zoom window per region: the central ZOOM_KEEP of both columns, with equal
#' x and y limits so y = x stays the diagonal.
zoom_window <- function(data, x, y) {
  tail <- (1 - ZOOM_KEEP) / 2
  data |>
    filter(!is.na(.data[[x]]), !is.na(.data[[y]])) |>
    group_by(region) |>
    summarise(
      n  = n(),
      lo = min(quantile(.data[[x]], tail),     quantile(.data[[y]], tail)),
      hi = max(quantile(.data[[x]], 1 - tail), quantile(.data[[y]], 1 - tail)),
      n_shown = sum(between(.data[[x]], lo, hi) & between(.data[[y]], lo, hi)),
      .groups = "drop"
    )
}

#' Hex-binned `y` against `x`, one panel per region, zoomed by zoom_window().
#' Rows outside the window are not binned.
mfe_scatter <- function(data, x, y, title, xlab, ylab) {
  zoom   <- zoom_window(data, x, y)
  shown  <- data |>
    inner_join(zoom, by = "region") |>
    filter(between(.data[[x]], lo, hi), between(.data[[y]], lo, hi))
  frame  <- tidyr::pivot_longer(zoom, c(lo, hi), values_to = "v")

  ggplot(shown, aes(x = .data[[x]], y = .data[[y]])) +
    geom_blank(data = frame, aes(x = v, y = v)) +
    geom_hex(bins = HEX_BINS) +
    density_fill +
    geom_abline(slope = 1, intercept = 0, colour = "#D55E00", linewidth = 0.5) +
    facet_wrap(~ region, ncol = 2, scales = "free") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(title = title, x = xlab, y = ylab) +
    theme_mfe
}

plot_obs_vs_pred <- mfe_scatter(
  mfe_long, x = "observed", y = "predicted",
  title = "Predicted vs observed MFE",
  xlab  = "Observed MFE (kcal/mol)", ylab = "Predicted MFE (kcal/mol)")

plot_shuf_vs_pred <- mfe_scatter(
  mfe_long, x = "shuffled", y = "predicted",
  title = "Predicted MFE vs dinucleotide-shuffle median",
  xlab  = "Dinucleotide-shuffle median MFE (kcal/mol)", ylab = "Predicted MFE (kcal/mol)")


# --- Report: the numbers for the figure legends --------------------------------
#' Markdown run report. The figures carry no subtitles or notes; this is where
#' their definitions and counts live, for writing captions.
mfe_report <- function(path) {
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
  pct <- function(x) sprintf("%.1f%%", 100 * x)
  md  <- function(d) c(paste("|", paste(names(d), collapse = " | "), "|"),
                       paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
                       apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
  k    <- attr(df, "mfe_model")
  if (is.null(k)) k <- mfe_model_provenance()

  by_region <- mfe_long |>
    group_by(region) |>
    summarise(
      n               = fmt(n()),
      `< 100 nt`      = pct(mean(length < SHORT_NT)),
      `> 4,000 nt`    = pct(mean(length > LONG_NT)),
      `predicted = 0` = fmt(sum(predicted == 0)),
      `off scale (fig. 1)` = fmt(sum(rel_observed < REL_YLIM[1] |
                                     rel_observed > REL_YLIM[2], na.rm = TRUE)),
      `median rel. diff, observed` = pct(median(rel_observed, na.rm = TRUE)),
      `median rel. diff, shuffle`  = pct(median(rel_shuffled, na.rm = TRUE)),
      .groups = "drop"
    )

  scatter_stats <- function(x, y) {
    zw <- zoom_window(mfe_long, x, y)
    st <- mfe_long |>
      filter(!is.na(.data[[x]]), !is.na(.data[[y]])) |>
      group_by(region) |>
      summarise(
        rho   = sprintf("%.3f", stats::cor(.data[[x]], .data[[y]], method = "spearman")),
        slope = sprintf("%.3f", unname(stats::coef(stats::lm(.data[[y]] ~ .data[[x]]))[2])),
        `above y = x` = pct(mean(.data[[y]] > .data[[x]])),
        .groups = "drop")
    zw |>
      inner_join(st, by = "region") |>
      transmute(region,
                `points shown` = sprintf("%s (%s)", fmt(n_shown), pct(n_shown / n)),
                n = fmt(n), `axis range (kcal/mol)` = sprintf("%.0f to %.0f", lo, hi),
                `Spearman rho` = rho, `slope (y on x)` = slope, `above y = x`) |>
      relocate(`points shown`, .after = `axis range (kcal/mol)`)
  }

  lines <- c(
    "# Observed vs predicted MFE: run report", "",
    sprintf("Generated %s from `analysis/qc/mfe_observed_vs_predicted.R`.", Sys.Date()), "",
    "The figures carry no subtitles or in-plot notes; the definitions and numbers for their captions are here.", "",
    "## Cohort and model", "",
    sprintf("- Human, analysed cohort: **%s** transcripts (`build_dataset()` defaults: 5' UTR >= %d nt, 3' UTR >= %d nt; the same frame the models use)",
            fmt(nrow(df)), MIN_5UTR_LENGTH, MIN_3UTR_LENGTH),
    sprintf("- Predicted MFE: RNAexpected model %s, (a GC^b + c) x length + d, capped at 0, GC as a fraction; a = %.7f, b = %.7f, c = %.7f, d = %.7f (RNAexpected %s, commit %s)",
            k$model, k$coefficients[["a"]], k$coefficients[["b"]], k$coefficients[["c"]],
            k$coefficients[["d"]], k$version, substr(k$sha %||% "unknown", 1, 7)),
    "- Predicted MFE is the average RNAfold MFE of a random sequence with the same GC content and length: a composition-matched baseline that does not see codon usage or dinucleotide bias",
    "- Fitted on 20-4,000 nt. Below ~100 nt it under-predicts structure (about -13% at 50-100 nt, worse below) and returns 0 below ~30-60 nt depending on GC; above 4,000 nt it is extrapolated",
    "- Observed MFE: RNAfold MFE of the region (`rnafold_score_<region>`). Dinucleotide-shuffle median: median MFE of 1,000 dinucleotide-preserving esl-shuffles of the same region (`rnafold_median_<region>`)",
    "- More negative MFE = more structured", "",
    md(by_region), "",
    "## Figure 1: relative difference vs length", "",
    "- y = (observed - predicted) / |predicted|; below 0 = more structured than predicted",
    "- x = region length, log scale; one panel per region, shared axes",
    sprintf("- Hexagonal bins (%d across each panel); fill = transcripts per hex, log scale", HEX_BINS),
    sprintf("- Dashed vertical lines at %s and %s nt: the model under-predicts structure below the first and extrapolates above the second",
            fmt(SHORT_NT), fmt(LONG_NT)),
    sprintf("- y axis clipped to %s to %s; rows outside are not drawn (\"off scale\" above, almost all short 5' UTRs where predicted MFE is near 0)",
            pct(REL_YLIM[1]), pct(REL_YLIM[2])),
    "- Rows with predicted MFE = 0 have no relative difference and are not drawn (\"predicted = 0\" above)",
    "- Relative rather than raw difference because random-sequence SD grows with length (~2 kcal/mol at 50 nt to ~14 at 4 kb)", "",
    "## Figure 2: predicted vs observed MFE", "",
    "- x = observed MFE, y = predicted MFE (kcal/mol); orange line y = x",
    "- Above the line = observed more negative than predicted = more structured than predicted",
    sprintf("- Each panel shows the central %s of its region: axis range from the lower of the two %s percentiles to the higher of the two %s percentiles, equal on x and y; points outside are not drawn. Scales differ by panel",
            pct(ZOOM_KEEP), pct((1 - ZOOM_KEEP) / 2), pct(1 - (1 - ZOOM_KEEP) / 2)),
    sprintf("- Hexagonal bins (%d across each panel); fill = transcripts per hex, log scale", HEX_BINS),
    "- Statistics are on all rows, not only those shown", "",
    md(scatter_stats("observed", "predicted")), "",
    "## Figure 3: predicted MFE vs dinucleotide-shuffle median", "",
    "- x = dinucleotide-shuffle median MFE, y = predicted MFE (kcal/mol); orange line y = x",
    "- Above the line = shuffles more structured than predicted: the dinucleotide composition the GC-only model does not see (plus model error)",
    "- Same zoom rule and hex binning as figure 2", "",
    md(scatter_stats("shuffled", "predicted")), "")
  writeLines(lines, path)
}


# --- Save ----------------------------------------------------------------------
for (d in c("plots", "reports"))
  dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  file.path(OUTPUT_DIR, "plots", sprintf("mfe_rel_vs_length_%s.jpg", species)),
  plot = plot_rel_vs_length, width = 250, height = 190, units = "mm", dpi = 300)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "plots", sprintf("mfe_obs_vs_pred_%s.jpg", species)),
  plot = plot_obs_vs_pred, width = 230, height = 200, units = "mm", dpi = 300)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "plots", sprintf("mfe_shuffle_vs_pred_%s.jpg", species)),
  plot = plot_shuf_vs_pred, width = 230, height = 200, units = "mm", dpi = 300)

report_path <- file.path(OUTPUT_DIR, "reports",
                         sprintf("mfe_observed_vs_predicted_report_%s.md", species))
mfe_report(report_path)

message("\nObserved vs predicted MFE:")
message("  ", file.path(OUTPUT_DIR, "plots"),
        "/mfe_{rel_vs_length,obs_vs_pred,shuffle_vs_pred}_", species, ".jpg")
message("  ", report_path)
