# =============================================================================
# Figures for the Baseline / Structure comparison — rendering only, no modelling
# =============================================================================
# Every figure here is built from tables already written to disk by
# xgb_structure_comparison.R. Nothing in this file fits, tunes or predicts
# anything, and it never loads the fitted models.
#
# WHY IT IS SEPARATE. Tuning is the great majority of the comparison script's
# runtime and every figure takes about a second. Restyling a plot should not
# cost a re-tune, so:
#
#   Rscript analysis/models/xgb_structure_plots.R      # figures only, ~5 s
#   Rscript analysis/models/xgb_structure_comparison.R # models, then figures
#
# The comparison script sources this file at the end, so a full run still
# produces everything in one command AND this file's inputs are exercised on
# every run — it cannot silently rot against a changed table format.
#
# THE CONTRACT. This script reads; it does not compute. If a number belongs on
# a figure it must first exist in a table written by the modelling script, so
# that what is on the slide can always be traced to a file. The one thing
# computed here is arithmetic already implied by those tables (a mean, an axis
# range) — never a metric, an interval or a model quantity.
#
# Anything genuinely derived from the fitted models (gain importance, the
# redundancy regression) is computed upstream and lands here as a CSV.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

if (!exists("OUTPUT_DIR")) source("R/load_all.R")
if (!exists("run_dir")) source("analysis/models/xgb_structure_features.R")


# ----------------------------- Style knobs ----------------------------------
# The things worth adjusting when a figure needs to fit a slide. Edit freely;
# nothing below depends on a model.

COL_BASELINE  <- "#4C6EF5"   # the reference model
COL_STRUCTURE <- "#E8590C"   # the structure model / a conclusive result
COL_NULL      <- "grey70"    # a CI that spans zero
BASE_SIZE     <- 12
PLOT_DPI      <- 300
PLOT_FORMATS  <- c("jpg", "pdf")

RUN_DIR   <- run_dir("root")
PLOT_DIR  <- run_dir("plots")
TABLE_DIR <- run_dir("tables")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# The metrics that appear ON A FIGURE, and their panel order. MAE is
# deliberately absent: it is still computed, bootstrapped, sign-flip tested and
# written to xgb_structure_delta_bootstrap.csv, so a reviewer asking about
# absolute error gets an answer without a re-run — it simply is not on a
# figure. Adding "mae" to this vector is the only change needed to put it back.
PLOT_METRICS <- c("rsq_trad", "rmse")

METRIC_LABEL <- c(rsq_trad = "R²", rmse = "RMSE", mae = "MAE",
                  rsq = "r² (correlation)")

#' Order a metric column as an axis / facet factor, dropping what is not drawn.
#'
#' One place decides both which metrics appear and in what order, so a figure
#' cannot quietly disagree with PLOT_METRICS.
metric_factor <- function(x) {
  keep <- intersect(PLOT_METRICS, unique(x))
  factor(METRIC_LABEL[x], levels = unname(METRIC_LABEL[keep]))
}

#' Restrict a table to the metrics that go on figures.
plot_metrics_only <- function(df) dplyr::filter(df, metric %in% PLOT_METRICS)

# Annotation policy: a short title, and a subtitle carrying SAMPLE SIZES only.
# Interpretation belongs in the paper's figure legend, where it can be as long
# as it needs to be and is read alongside the text. A figure that argues with
# itself in 9pt grey is harder to read and impossible to typeset.
theme_xgb <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9, colour = "grey35"),
          strip.text    = element_text(face = "bold"),
          plot.margin   = margin(6, 10, 6, 6))
}

#' Write one plot in every configured format.
#'
#' PDF goes through cairo_pdf, not R's default pdf device. The default device
#' has no UTF-8 support: it silently rendered "Spearman ρ" as "Spearman ." in
#' the PDF while the JPEG was fine, with nothing but a buried mbcsToSbcs
#' warning to show for it. Since the PDFs are the presentation assets, a
#' silently corrupted glyph is the worst possible failure mode. cairo_pdf
#' handles ρ, ² and the em dash correctly.
save_plot <- function(p, name, w = 210, h = 148, dir = PLOT_DIR) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in PLOT_FORMATS) {
    args <- list(filename = file.path(dir, paste0(name, ".", ext)),
                 plot = p, width = w, height = h, units = "mm", dpi = PLOT_DPI)
    if (ext == "pdf") args$device <- grDevices::cairo_pdf
    do.call(ggsave, args)
  }
  message("  wrote ", file.path(dir, name), " (",
          paste(PLOT_FORMATS, collapse = ", "), ")")
}


# ----------------------------- Artefact loading -----------------------------
# Named so a missing file names the script that produces it, rather than
# failing somewhere downstream with an unhelpful "object not found".

MAIN <- "analysis/models/xgb_structure_comparison.R"

need <- function(path) {
  if (!file.exists(path)) {
    stop("missing artefact: ", path,
         "\n  run `Rscript ", MAIN, "` first", call. = FALSE)
  }
  path
}

manifest <- readRDS(need(file.path(RUN_DIR, "run_manifest.rds")))

PLOT_MODELS <- manifest$models
REFMOD      <- manifest$reference_model
STRMOD      <- manifest$structure_model

as_model <- function(x) factor(x, levels = PLOT_MODELS)

preds  <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_test_predictions.csv")),
                   show_col_types = FALSE) |> mutate(model = as_model(model))
deltas <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_delta_bootstrap.csv")),
                   show_col_types = FALSE)
chunks <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_chunk_metrics.csv")),
                   show_col_types = FALSE) |> mutate(model = as_model(model))
imp    <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_gain_importance.csv")),
                   show_col_types = FALSE) |> mutate(model = as_model(model))
redun  <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_redundancy.csv")),
                   show_col_types = FALSE)

N_BOOT <- manifest$n_boot
N_TEST <- manifest$n_test

MODEL_COLS <- setNames(c(COL_BASELINE, COL_STRUCTURE), PLOT_MODELS)

# Predictor counts per model are no longer on any figure — they belong in the
# paper's methods, and manifest$model_metrics / the feature-list CSV remain the
# authority. Sample sizes are the one thing the subtitles still carry, because
# a reader cannot judge an interval without n.

message("Rendering figures from ", RUN_DIR, " (", N_TEST, " held-out genes)")


# ----------------------------- 1. Delta metrics -----------------------------
# The primary figure. Fill encodes CONCLUSIVE (does the CI clear zero), not the
# sign of the point estimate — colouring by sign would paint a point orange for
# landing a hair on the favourable side of zero with an interval ten times its
# own width, which is the exact misreading this figure exists to prevent.

# Fill still encodes CONCLUSIVE rather than the sign of the point estimate,
# but the explanation moves to the paper's legend. Both points are grey here,
# so nothing is ambiguous on the figure as it stands; the key exists so that a
# future run in which an interval DOES clear zero cannot be read as a win
# without checking which side of zero it cleared.
p_delta <- deltas |>
  plot_metrics_only() |>
  mutate(metric_lab = metric_factor(metric)) |>
  ggplot(aes(x = delta, y = "")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, linewidth = 0.8, colour = "grey25") +
  geom_point(aes(fill = conclusive), size = 4, shape = 21, colour = "grey20",
             show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = COL_STRUCTURE, `FALSE` = COL_NULL)) +
  facet_wrap(~ metric_lab, scales = "free_x") +
  labs(
    title    = "Held-out performance: Structure vs Baseline",
    subtitle = sprintf("%s held-out genes, %s bootstrap replicates, 95%% CI",
                       format(N_TEST, big.mark = ","),
                       format(N_BOOT, big.mark = ",")),
    x = "Difference (Structure − Baseline)",
    y = NULL
  ) +
  theme_xgb() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Short: one row of data, so height beyond the title block is dead space.
save_plot(p_delta, "xgb_structure_delta_metrics", w = 200, h = 58)


# ----------------------------- 2. Paired slices -----------------------------
# Is a pooled delta broad-based, or carried by one corner of the held-out genes?

# Only the TILT of each line is informative — panel heights differ because some
# slices contain intrinsically harder genes. That belongs in the legend, not on
# the figure.
p_paired <- chunks |>
  plot_metrics_only() |>
  mutate(metric_lab = metric_factor(metric)) |>
  ggplot(aes(x = model, y = value, group = chunk)) +
  geom_line(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = model), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = MODEL_COLS) +
  facet_wrap(~ metric_lab, scales = "free_y") +
  labs(
    title    = "Performance across held-out slices",
    subtitle = sprintf("%d family-blocked slices of the held-out set, one line each",
                       n_distinct(chunks$chunk)),
    x = NULL, y = NULL
  ) +
  theme_xgb()

save_plot(p_paired, "xgb_structure_paired_slices", w = 190, h = 110)


# ----------------------------- 3. Observed vs predicted ---------------------

# Three numbers, because they answer three different questions.
#   Pearson r    how tight the cloud is about SOME straight line
#   Spearman ρ   whether the RANKING is right, ignoring scale entirely
#   R²           the coefficient of determination, 1 - SSres/SStot
# R² is unqualified because it already means the coefficient of determination.
# The term needing care is the other one: yardstick's `rsq` is squared Pearson
# correlation, which many people also call R². They are identical for in-sample
# OLS with an intercept and diverge out of sample by exactly a calibration
# penalty, so reporting r beside R² is a free calibration check.
obspred_stats <- preds |>
  group_by(model) |>
  summarise(
    pearson   = cor(observed, predicted),
    spearman  = cor(observed, predicted, method = "spearman"),
    rsq_trad  = 1 - sum((observed - predicted)^2) /
                    sum((observed - mean(observed))^2),
    .groups   = "drop"
  ) |>
  mutate(label = sprintf("r = %.3f\n%s = %.3f\nR%s = %.3f",
                         pearson, "ρ", spearman, "²", rsq_trad))

lims <- range(c(preds$observed, preds$predicted))

# The stats block stays: it is data, not commentary. Kept compact — the symbols
# are defined in the legend rather than spelled out in the panel.
p_obspred <- preds |>
  ggplot(aes(observed, predicted)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(aes(colour = model), alpha = 0.25, size = 0.7, show.legend = FALSE) +
  geom_text(data = obspred_stats, inherit.aes = FALSE,
            aes(x = lims[1], y = lims[2], label = label),
            hjust = 0, vjust = 1, size = 2.9, lineheight = 1.15,
            colour = "grey15") +
  scale_colour_manual(values = MODEL_COLS) +
  coord_equal(xlim = lims, ylim = lims) +
  facet_wrap(~ model) +
  labs(title = "Observed vs predicted half-life",
       subtitle = sprintf("%s held-out genes, plotted in both panels",
                          format(N_TEST, big.mark = ",")),
       x = "Observed (PC1 score)", y = "Predicted (PC1 score)") +
  theme_xgb()

save_plot(p_obspred, "xgb_structure_observed_vs_predicted", w = 190, h = 115)


# ----------------------------- 4. Gain importance ---------------------------
# Exploratory context only, drawn for the Structure model. The gain values are
# computed upstream from the fitted models; this only draws them.

# Correlated predictors share and redistribute gain, so this is not an effect
# estimate and does not rank independent contributions — a caveat for the
# legend, not the panel. The block legend stays: which bars are structure is
# the whole point of the figure.
p_imp <- imp |>
  filter(model == STRMOD) |>
  slice_max(Gain, n = 25) |>
  mutate(block = factor(block, levels = c("baseline", "structure"),
                        labels = c("Baseline", "Structure"))) |>
  ggplot(aes(Gain, fct_reorder(Feature, Gain), fill = block)) +
  geom_col() +
  scale_fill_manual(values = setNames(c("grey65", COL_STRUCTURE),
                                      c("Baseline", "Structure"))) +
  scale_y_discrete(labels = function(x) format_col_name(x)) +
  labs(title = "Gain importance, Structure model (top 25)",
       subtitle = "Exploratory: gain is shared among correlated predictors",
       x = "Gain", y = NULL, fill = NULL) +
  theme_xgb() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

save_plot(p_imp, "xgb_structure_gain_importance", w = 190, h = 160)


# ----------------------------- 5. Structure redundancy ----------------------
# The companion to figure 1, and the reason it is a figure rather than a
# buried CSV: the structure block is not length- and GC-neutral, so a positive
# delta has two readings and this is what tells them apart. Each point is one
# structure column; the x-axis is the share of its variance a linear model on
# the whole baseline block already reproduces.

if (nrow(redun) > 0) {

  # The 0.6 reference line is gone. It was a reading aid I chose, it has no
  # statistical standing, and a dashed line invites being read as a threshold
  # or a test. Replaced by each family's MEDIAN R² as a crossbar — the same
  # at-a-glance ordering, but a property of the data rather than of an
  # arbitrary cutoff. Colour now encodes the family, not which side of a line
  # a point fell on.
  redun_fam <- redun |>
    group_by(family) |>
    summarise(median_r2 = median(r2_from_baseline), .groups = "drop")

  p_redun <- redun |>
    left_join(redun_fam, by = "family") |>
    mutate(family_lab = fct_reorder(format_group_name(family, "group"),
                                    median_r2)) |>
    ggplot(aes(r2_from_baseline, family_lab)) +
    geom_crossbar(data = ~ distinct(.x, family_lab, median_r2),
                  aes(x = median_r2, xmin = median_r2, xmax = median_r2),
                  width = 0.55, linewidth = 0.5, colour = "grey35",
                  middle.linewidth = 0.5) +
    geom_point(alpha = 0.75, size = 2.4, colour = COL_STRUCTURE,
               position = position_jitter(height = 0.10, width = 0, seed = 42)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = "Structure features explained by the baseline",
      subtitle = sprintf("One point per structure column (n = %d); bar marks the family median",
                         nrow(redun)),
      x = "R² from the baseline block", y = NULL
    ) +
    theme_xgb() +
    theme(panel.grid.major.y = element_blank())

  save_plot(p_redun, "xgb_structure_redundancy", w = 190, h = 115)
} else {
  message("  skipping the redundancy figure — the regression was not run")
}

message("Figures written to ", PLOT_DIR)
