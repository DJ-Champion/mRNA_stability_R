# =============================================================================
# Figures for the structure ladder — rendering only, no modelling
# =============================================================================
# Every figure here is built from tables already written to disk by
# xgb_structure_comparison.R and xgb_structure_gini_subset.R. Nothing in this
# file fits, tunes or predicts anything, and it never loads the fitted models.
#
# WHY IT IS SEPARATE. Tuning is the great majority of the comparison script's
# runtime and every figure takes about a second. Restyling a plot should not
# cost a re-tune, and for a while it did. Now:
#
#   Rscript analysis/models/xgb_structure_plots.R      # figures only, ~5 s
#   Rscript analysis/models/xgb_structure_comparison.R # models, then figures
#
# The comparison script sources this file at the end, so a full run still
# produces everything in one command AND this file's inputs are exercised on
# every run — it cannot silently rot against a changed table format.
#
# THE CONTRACT. This script reads; it does not compute. If a number belongs on
# a figure it must first exist in a table written by the modelling scripts, so
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
if (!exists("resolve_variant")) source("analysis/models/xgb_structure_features.R")

# Which variant's artefacts to draw. Inherited from the calling script when
# sourced; otherwise taken from XGB_VARIANT, defaulting to `default`:
#
#   Rscript analysis/models/xgb_structure_plots.R                     # default
#   XGB_VARIANT=complete_case Rscript analysis/models/xgb_structure_plots.R
PLOT_VARIANT <- if (exists("VARIANT")) VARIANT else
  resolve_variant(Sys.getenv("XGB_VARIANT", "default"))


# ----------------------------- Style knobs ----------------------------------
# The things worth adjusting when a figure needs to fit a slide. Edit freely;
# nothing below depends on a model.

COL_REFERENCE <- "#4C6EF5"   # the baseline rung
COL_PRIMARY   <- "#E8590C"   # the pre-specified primary rung / conclusive result
COL_NULL      <- "grey70"    # a CI that spans zero
BASE_SIZE     <- 12
PLOT_DPI      <- 300
PLOT_FORMATS  <- c("jpg", "pdf")

RUN_DIR   <- variant_dir(PLOT_VARIANT, "root")
PLOT_DIR  <- variant_dir(PLOT_VARIANT, "plots")
TABLE_DIR <- variant_dir(PLOT_VARIANT, "tables")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

METRIC_LABEL <- c(rsq_trad = "R² (held-out)", rmse = "RMSE", mae = "MAE",
                  rsq = "r² (correlation)")

KIND_LABEL <- c(vs_baseline = "Each rung vs the baseline",
                increment   = "Rung-to-rung increments")

# Any figure from a non-default specification carries the variant in its title.
# Without it, a sensitivity figure is visually indistinguishable from the
# committed result and could be presented as the headline by accident — which
# is precisely the failure mode the sensitivity framing exists to prevent.
title_for <- function(txt) {
  if (identical(PLOT_VARIANT$name, "default")) txt
  else paste0("[", PLOT_VARIANT$name, "] ", txt)
}
variant_note <- function() {
  if (identical(PLOT_VARIANT$name, "default")) ""
  else sprintf("\nSENSITIVITY VARIANT '%s': %s. Not the committed result.",
               PLOT_VARIANT$name, PLOT_VARIANT$label)
}

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

need <- function(path, produced_by) {
  if (!file.exists(path)) {
    stop("missing artefact: ", path,
         "\n  run `Rscript ", produced_by, "` first", call. = FALSE)
  }
  path
}

optional <- function(path) if (file.exists(path)) path else NULL

MAIN <- "analysis/models/xgb_structure_comparison.R"
GINI <- "analysis/models/xgb_structure_gini_subset.R"

manifest <- readRDS(need(file.path(RUN_DIR, "run_manifest.rds"), MAIN))

MODELS  <- manifest$models
PRIMARY <- manifest$primary_model
REFMOD  <- manifest$reference_model

as_model <- function(x) factor(x, levels = MODELS)

preds   <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_test_predictions.csv"), MAIN),
                    show_col_types = FALSE) |> mutate(model = as_model(model))
deltas  <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_delta_bootstrap.csv"), MAIN),
                    show_col_types = FALSE)
chunks  <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_chunk_metrics.csv"), MAIN),
                    show_col_types = FALSE) |> mutate(model = as_model(model))
imp     <- read_csv(need(file.path(TABLE_DIR, "xgb_structure_gain_importance.csv"), MAIN),
                    show_col_types = FALSE) |> mutate(model = as_model(model))

N_BOOT <- manifest$n_boot
N_TEST <- manifest$n_test

# Reference blue, then warming toward the primary orange as structure is added.
# Deliberately NOT a rainbow: the ladder is ordered, so the palette is ordered.
MODEL_COLS <- setNames(
  c(COL_REFERENCE,
    colorRampPalette(c("#FFA94D", COL_PRIMARY))(length(MODELS) - 1)),
  MODELS)

message("Rendering figures from ", RUN_DIR, " (", N_TEST, " held-out genes, ",
        length(MODELS), " rungs)")


# ----------------------------- 1. Delta metrics -----------------------------
# The primary figure. Fill encodes CONCLUSIVE (does the CI clear zero), not the
# sign of the point estimate — colouring by sign would paint a point orange for
# landing a hair on the favourable side of zero with an interval ten times its
# own width, which is the exact misreading this figure exists to prevent.
#
# The pre-specified contrast is marked with a bold label and a caret. Every
# other row is secondary and unadjusted, and the subtitle says so on the
# figure, not only in the report — a plot travels further than its caption.

delta_plot_data <- deltas |>
  mutate(
    metric_lab = factor(METRIC_LABEL[metric],
                        levels = METRIC_LABEL[c("rsq_trad", "rmse", "mae")]),
    kind_lab   = factor(KIND_LABEL[kind], levels = KIND_LABEL),
    contrast_lab = fct_rev(factor(contrast, levels = unique(deltas$contrast))),
    face       = if_else(primary, "bold", "plain")
  )

p_delta <- delta_plot_data |>
  ggplot(aes(x = delta, y = contrast_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.2, linewidth = 0.8, colour = "grey25") +
  geom_point(aes(fill = conclusive), size = 4, shape = 21, colour = "grey20",
             show.legend = FALSE) +
  geom_text(data = filter(delta_plot_data, primary),
            aes(x = delta, y = contrast_lab, label = "▲"),
            vjust = 2.1, size = 3, colour = COL_PRIMARY) +
  scale_fill_manual(values = c(`TRUE` = COL_PRIMARY, `FALSE` = COL_NULL)) +
  # Free x per metric column, free y per contrast row, row heights proportional
  # to the number of contrasts. Sharing the x scale down a column is deliberate:
  # it puts the vs-baseline and increment families on one axis per metric, so
  # their magnitudes are read against each other rather than separately rescaled.
  facet_grid(kind_lab ~ metric_lab, scales = "free", space = "free_y") +
  labs(
    title    = title_for("Change in held-out performance along the structure ladder"),
    subtitle = sprintf(paste0("%d held-out genes; point estimate and 95%% paired-",
                              "bootstrap CI (%d reps, every rung scored on the same draw).\n",
                              "R²: higher is better, so positive favours the larger model. ",
                              "RMSE / MAE: lower is better, so negative does.\n",
                              "Filled points clear zero; grey points do not. ",
                              "▲ marks the ONE pre-specified contrast (%s vs %s);\n",
                              "every other row is secondary and unadjusted — read them ",
                              "for the SHAPE of the ladder, not individually.%s"),
                       N_TEST, N_BOOT, PRIMARY, REFMOD, variant_note()),
    x = "Delta (larger model - smaller model), PC1 units", y = NULL
  ) +
  theme_xgb(subtitle_size = 9) +
  theme(panel.grid.major.y = element_blank())

save_plot(p_delta, "xgb_structure_delta_metrics", w = 280, h = 175)


# ----------------------------- 2. Paired slices -----------------------------
# The brief asked for a paired plot across outer CV folds. Under a single
# holdout there are no outer folds, so this slices the test set instead and
# answers the same question: is a pooled delta broad-based, or carried by one
# corner of the held-out genes?

p_paired <- chunks |>
  mutate(metric_lab = factor(METRIC_LABEL[metric],
                             levels = METRIC_LABEL[c("rsq_trad", "rmse", "mae")])) |>
  ggplot(aes(x = model, y = value, group = chunk)) +
  geom_line(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = model), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = MODEL_COLS) +
  facet_wrap(~ metric_lab, scales = "free_y") +
  labs(
    title    = title_for("Paired performance across family-blocked slices of the held-out set"),
    subtitle = paste0("One line per slice, walking the ladder left to right. ",
                      "A consistent effect tilts every line the same way;\n",
                      "lines that cross mean the pooled delta is carried by one ",
                      "corner of the test set. Heights differ because some\n",
                      "slices contain intrinsically harder genes — only the TILT ",
                      "is informative.", variant_note()),
    x = NULL, y = NULL
  ) +
  theme_xgb(subtitle_size = 9) +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))

save_plot(p_paired, "xgb_structure_paired_slices", w = 290, h = 155)


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
    cal_slope = stats::coef(stats::lm(observed ~ predicted))[2],
    sd_pred   = sd(predicted),
    .groups   = "drop"
  ) |>
  mutate(label = sprintf("Pearson r = %.3f\nSpearman %s = %.3f\nR%s = %.3f",
                         pearson, "ρ", spearman, "²", rsq_trad))

lims <- range(c(preds$observed, preds$predicted))
ref_stats <- obspred_stats |> filter(model == REFMOD)

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
  facet_wrap(~ model, nrow = 1) +
  labs(title = title_for("Observed vs held-out predicted half-life"),
       subtitle = sprintf(paste0("Identical axis limits. Dashed line is y = x, not a fit. ",
                                 "Agarwal & Kelley consensus PC1, untransformed.\n",
                                 "R² is the coefficient of determination (1 - SSres/SStot), ",
                                 "so it measures scatter about the dashed line.\n",
                                 "Predictions span less than observations (sd %.1f vs %.1f). ",
                                 "That is calibrated shrinkage, not a defect: a well-\n",
                                 "calibrated model shrinks by about factor r, and r x %.1f = ",
                                 "%.1f. Regressing observed on predicted gives slope %.2f.\n",
                                 "Consequence: the model never predicts the most extreme ",
                                 "half-lives. Applies equally to every rung.\n",
                                 "Secondary to the paired comparison."),
                          ref_stats$sd_pred, sd(preds$observed[preds$model == REFMOD]),
                          sd(preds$observed[preds$model == REFMOD]),
                          ref_stats$pearson * sd(preds$observed[preds$model == REFMOD]),
                          ref_stats$cal_slope),
       x = "Observed half-life (PC1)", y = "Predicted half-life (PC1)") +
  theme_xgb(subtitle_size = 8.5)

save_plot(p_obspred, "xgb_structure_observed_vs_predicted", w = 300, h = 155)


# ----------------------------- 4. Gain importance ---------------------------
# Exploratory context only, drawn for the PRIMARY rung. The gain values are
# computed upstream from the fitted models; this only draws them.
#
# The primary rung rather than the widest: S-full's block is the most
# confounded with the baseline, so its importance figure would show structure
# collecting gain for length and GC and invite exactly the misreading this
# caption warns against. The primary rung is where "genuinely new numbers, used
# by the model, buying nothing" is visible cleanly.

imp_primary <- imp |> filter(model == PRIMARY)

p_imp <- imp_primary |>
  slice_max(Gain, n = 25) |>
  ggplot(aes(Gain, fct_reorder(Feature, Gain), fill = block)) +
  geom_col() +
  scale_fill_manual(values = c(baseline = "grey65", structure = COL_PRIMARY)) +
  scale_y_discrete(labels = function(x) vapply(x, format_col_name, character(1))) +
  labs(title = title_for(paste0(PRIMARY, " gain importance (top 25): exploratory only")),
       subtitle = paste0("Correlated predictors share and redistribute gain, so ",
                         "this is not an effect estimate and does not rank\n",
                         "independent contributions. The statistical evidence is ",
                         "the held-out improvement, not this figure."),
       x = "Gain", y = NULL, fill = NULL) +
  theme_xgb() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

save_plot(p_imp, "xgb_structure_primary_importance", h = 180)


# ----------------------------- 5. icSHAPE secondary -------------------------
# Rendered only if the secondary run has been done. Absent is not an error —
# the main comparison stands on its own.

# The icSHAPE run has its own directory, not a variant's — it is a different
# design on a different eligible set, always built from the primary rung's
# feature blocks. So this section is independent of PLOT_VARIANT and renders
# the same figure whichever variant is being drawn.
GINI_DIR  <- file.path(OUTPUT_DIR, "xgb_structure", "gini")
gini_path <- optional(file.path(GINI_DIR, "tables", "xgb_structure_gini_deltas.csv"))
gini_man  <- optional(file.path(GINI_DIR, "gini_run_manifest.rds"))

if (is.null(gini_path) || is.null(gini_man)) {
  message("  skipping the icSHAPE figure — run `Rscript ", GINI, "` to produce it")
} else {
  gd  <- read_csv(gini_path, show_col_types = FALSE)
  gm  <- readRDS(gini_man)
  lab <- c(rsq_trad = "R² (out-of-fold)", rmse = "RMSE", mae = "MAE")

  p_gini <- gd |>
    mutate(metric_lab = factor(lab[metric], levels = lab),
           comparison = factor(comparison, levels = rev(unique(comparison)))) |>
    ggplot(aes(delta, comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                  width = 0.15, linewidth = 0.8, colour = "grey25") +
    geom_point(aes(fill = conclusive), size = 4, shape = 21, colour = "grey20",
               show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = COL_PRIMARY, `FALSE` = "grey75")) +
    facet_wrap(~ metric_lab, scales = "free_x") +
    labs(
      title = "SECONDARY: icSHAPE Gini on the genes that have it",
      subtitle = sprintf(paste0("Out-of-fold, %d genes, %s. ",
                                "Grey points span zero. NOT the headline comparison.\n",
                                "Caveats: (1) icSHAPE coverage tracks expression, so ",
                                "this subset is biased toward\nabundant transcripts. ",
                                "(2) Gini is MEASURED, not computed from sequence, so ",
                                "a model\nusing it cannot score an unprobed transcript. ",
                                "(3) Probing read depth tracks abundance,\nwhich this ",
                                "design cannot separate from a structure effect."),
                         gm$n_genes, gm$design),
      x = "Delta (candidate - reference), PC1 units", y = NULL
    ) +
    theme_xgb(subtitle_size = 9) +
    theme(panel.grid.major.y = element_blank())

  save_plot(p_gini, "xgb_structure_gini_deltas", w = 300, h = 155,
            dir = file.path(GINI_DIR, "plots"))
}

message("Figures written to ", PLOT_DIR)
