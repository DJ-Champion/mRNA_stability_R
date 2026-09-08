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

METRIC_LABEL <- c(rsq_trad = "R² (held-out)", rmse = "RMSE", mae = "MAE",
                  rsq = "r² (correlation)")

theme_xgb <- function(subtitle_size = 10) {
  theme_minimal(base_size = BASE_SIZE) +
    theme(plot.title       = element_text(face = "bold", size = 15),
          plot.subtitle    = element_text(size = subtitle_size, colour = "grey30"),
          strip.text       = element_text(face = "bold"))
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

# Predictor counts per model, read from the manifest — the one place they are
# recorded — so a caption cannot drift from what was actually fitted.
n_pred <- manifest$model_metrics |>
  mutate(txt = sprintf("%s %d", model, n_predictors))

message("Rendering figures from ", RUN_DIR, " (", N_TEST, " held-out genes)")


# ----------------------------- 1. Delta metrics -----------------------------
# The primary figure. Fill encodes CONCLUSIVE (does the CI clear zero), not the
# sign of the point estimate — colouring by sign would paint a point orange for
# landing a hair on the favourable side of zero with an interval ten times its
# own width, which is the exact misreading this figure exists to prevent.

delta_plot_data <- deltas |>
  mutate(metric_lab = factor(METRIC_LABEL[metric],
                             levels = METRIC_LABEL[c("rsq_trad", "rmse", "mae")]))

p_delta <- delta_plot_data |>
  ggplot(aes(x = delta, y = "")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, linewidth = 0.8, colour = "grey25") +
  geom_point(aes(fill = conclusive), size = 4, shape = 21, colour = "grey20",
             show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = COL_STRUCTURE, `FALSE` = COL_NULL)) +
  facet_wrap(~ metric_lab, scales = "free_x") +
  labs(
    title    = "Change in held-out performance: Structure vs Baseline",
    # The fill key is not decoration. Fill encodes "the CI clears zero", NOT
    # "Structure won" — a filled point can sit on the side favouring Baseline,
    # and unlabelled it reads as a win.
    subtitle = sprintf(paste0("n = %s genes trained, %s held out, %s bootstrap ",
                              "replicates.\nPredictors: %s.\n",
                              "Filled points: CI excludes zero — read which ",
                              "side. Grey points: CI spans zero."),
                       format(manifest$n_trainval, big.mark = ","),
                       format(N_TEST, big.mark = ","),
                       format(N_BOOT, big.mark = ","),
                       paste(n_pred$txt, collapse = ", ")),
    x = "Structure - Baseline (R² unitless; RMSE and MAE in half-life score units)",
    y = NULL
  ) +
  theme_xgb(subtitle_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

save_plot(p_delta, "xgb_structure_delta_metrics", w = 250, h = 95)


# ----------------------------- 2. Paired slices -----------------------------
# Is a pooled delta broad-based, or carried by one corner of the held-out genes?

p_paired <- chunks |>
  mutate(metric_lab = factor(METRIC_LABEL[metric],
                             levels = METRIC_LABEL[c("rsq_trad", "rmse", "mae")])) |>
  ggplot(aes(x = model, y = value, group = chunk)) +
  geom_line(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = model), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = MODEL_COLS) +
  facet_wrap(~ metric_lab, scales = "free_y") +
  labs(
    title    = "Paired performance across family-blocked slices of the held-out set",
    subtitle = paste0("One line per slice, Baseline on the left and Structure ",
                      "on the right. A consistent effect tilts every line the\n",
                      "same way; lines that cross mean the pooled delta is ",
                      "carried by one corner of the test set. Heights differ\n",
                      "because some slices contain intrinsically harder genes ",
                      "— only the TILT is informative."),
    x = NULL, y = NULL
  ) +
  theme_xgb(subtitle_size = 9)

save_plot(p_paired, "xgb_structure_paired_slices", w = 290, h = 145)


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
  mutate(label = sprintf("Pearson r = %.3f\nSpearman %s = %.3f\nR%s = %.3f",
                         pearson, "ρ", spearman, "²", rsq_trad))

lims <- range(c(preds$observed, preds$predicted))

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
  labs(title = "Observed vs held-out predicted half-life",
       subtitle = sprintf(paste0("n = %s genes trained, %s held out and plotted ",
                                 "in both panels.\nPredictors: %s."),
                          format(manifest$n_trainval, big.mark = ","),
                          format(N_TEST, big.mark = ","),
                          paste(n_pred$txt, collapse = ", ")),
       x = "Observed half-life score", y = "Predicted half-life score") +
  theme_xgb(subtitle_size = 9)

save_plot(p_obspred, "xgb_structure_observed_vs_predicted", w = 230, h = 135)


# ----------------------------- 4. Gain importance ---------------------------
# Exploratory context only, drawn for the Structure model. The gain values are
# computed upstream from the fitted models; this only draws them.

p_imp <- imp |>
  filter(model == STRMOD) |>
  slice_max(Gain, n = 25) |>
  ggplot(aes(Gain, fct_reorder(Feature, Gain), fill = block)) +
  geom_col() +
  scale_fill_manual(values = c(baseline = "grey65", structure = COL_STRUCTURE)) +
  scale_y_discrete(labels = function(x) vapply(x, format_col_name, character(1))) +
  labs(title = paste0(STRMOD, " gain importance (top 25): exploratory only"),
       subtitle = paste0("Correlated predictors share and redistribute gain, so ",
                         "this is not an effect estimate and does not rank\n",
                         "independent contributions. The statistical evidence is ",
                         "the held-out improvement, not this figure."),
       x = "Gain", y = NULL, fill = NULL) +
  theme_xgb() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

save_plot(p_imp, "xgb_structure_gain_importance", h = 180)


# ----------------------------- 5. Structure redundancy ----------------------
# The companion to figure 1, and the reason it is a figure rather than a
# buried CSV: the structure block is not length- and GC-neutral, so a positive
# delta has two readings and this is what tells them apart. Each point is one
# structure column; the x-axis is the share of its variance a linear model on
# the whole baseline block already reproduces.

if (nrow(redun) > 0) {
  p_redun <- redun |>
    mutate(family = fct_reorder(family, r2_from_baseline, .fun = median)) |>
    ggplot(aes(r2_from_baseline, family)) +
    geom_vline(xintercept = 0.6, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = r2_from_baseline > 0.6), alpha = 0.75, size = 2.4,
               show.legend = FALSE,
               position = position_jitter(height = 0.12, width = 0, seed = 42)) +
    scale_colour_manual(values = c(`TRUE` = COL_STRUCTURE, `FALSE` = "grey45")) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title    = "How much of each structure feature the baseline already explains",
      subtitle = paste0("One point per structure column: R² from a linear ",
                        "regression on the whole baseline block, training\n",
                        "genes only. Right of the dashed line the column is ",
                        "largely a restatement of length and GC, both\n",
                        "already in the baseline. Linear R², so this is a LOWER ",
                        "bound on redundancy."),
      x = "R² reconstructible from the baseline block", y = NULL
    ) +
    theme_xgb(subtitle_size = 9) +
    theme(panel.grid.major.y = element_blank())

  save_plot(p_redun, "xgb_structure_redundancy", w = 230, h = 130)
} else {
  message("  skipping the redundancy figure — the regression was not run")
}

message("Figures written to ", PLOT_DIR)
