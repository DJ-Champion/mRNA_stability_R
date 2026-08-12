# =============================================================================
# Sensitivity table: every variant that has been run, side by side
# =============================================================================
# This exists to make selective reporting difficult.
#
# The committed result (variant `default`) is defensible because it was
# specified before anyone saw an answer. Run enough alternative specifications
# and one will clear zero by chance — and it will be the one that feels most
# publishable. Reporting that one as the finding is the garden of forking
# paths, and it would be indefensible the moment anyone asks how many
# specifications were tried.
#
# So this script does not take arguments and has no filters. It discovers every
# variant directory on disk and reports all of them, including any that
# disagree with the rest. Adding a variant to the registry and running it
# automatically adds a row here; there is no supported way to run a variant and
# leave it out of the table. If the answer to "how many did you try?" is not
# "all of the ones in this table", something has gone wrong.
#
# THE CLAIM TO MAKE is "the conclusion holds across every specification we
# tried", which is stronger than any single run. If one variant dissents, that
# is a finding about the specification and needs explaining — not promoting.
#
# HOW TO READ IT. Within a variant, the paired delta and its CI are meaningful:
# both models saw the same genes. BETWEEN variants, the absolute R² values are
# NOT comparable when the row policy differs, because the models are scored on
# different genes. That is why n is printed on every row and every point.
#
# Usage:
#   Rscript analysis/models/xgb_structure_variant_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("R/load_all.R")
source("analysis/models/xgb_structure_features.R")

ROOT      <- file.path(OUTPUT_DIR, "xgb_structure")
PLOT_DIR  <- file.path(ROOT, "_summary", "plots")
TABLE_DIR <- file.path(ROOT, "_summary", "tables")
for (d in c(PLOT_DIR, TABLE_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)


# ----------------------------- Discover runs --------------------------------
# Driven by the registry, not by a glob of the filesystem: a variant that is
# defined but not yet run must appear as MISSING rather than silently vanish.
# A glob would let "delete the directory" work as a way to drop an inconvenient
# result, which is exactly what this file is meant to prevent.

rows <- map(names(VARIANTS), function(nm) {
  v      <- resolve_variant(nm)
  deltas <- file.path(variant_dir(v, "tables"), "xgb_structure_delta_bootstrap.csv")
  man    <- file.path(variant_dir(v, "root"), "run_manifest.rds")

  if (!file.exists(deltas) || !file.exists(man)) {
    return(tibble(variant = nm, label = v$label, status = "NOT RUN",
                  n_test = NA_integer_, n_eligible = NA_integer_,
                  n_structure = NA_integer_, metric = NA_character_,
                  baseline = NA_real_, structure = NA_real_, delta = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_, conclusive = NA))
  }

  m <- readRDS(man)
  read_csv(deltas, show_col_types = FALSE) |>
    transmute(variant = nm, label = v$label, status = "run",
              n_test = m$n_test, n_eligible = m$n_eligible,
              n_structure = length(m$structure_cols),
              metric, baseline, structure, delta, ci_low, ci_high, conclusive)
}) |> list_rbind()

write_csv(rows, file.path(TABLE_DIR, "xgb_structure_variant_summary.csv"))

not_run <- rows |> filter(status == "NOT RUN") |> distinct(variant)
if (nrow(not_run) > 0) {
  message("\nDefined but not yet run: ", paste(not_run$variant, collapse = ", "))
  message("  Rscript analysis/models/xgb_structure_comparison.R <variant>")
}

done <- rows |> filter(status == "run")
if (nrow(done) == 0) {
  stop("no variant has been run yet — run the default first", call. = FALSE)
}


# ----------------------------- Report ---------------------------------------

cat("\n=== Sensitivity across specifications ===\n")
cat("Within a variant the delta is paired and meaningful.\n")
cat("BETWEEN variants, absolute R² is not comparable when n differs.\n\n")

print(as.data.frame(
  done |>
    filter(metric == "rsq_trad") |>
    select(variant, n_eligible, n_test, n_structure,
           baseline, structure, delta, ci_low, ci_high, conclusive)),
  row.names = FALSE, digits = 4)

n_conclusive <- done |> filter(metric == "rsq_trad", conclusive) |> nrow()
n_total      <- done |> filter(metric == "rsq_trad") |> nrow()

cat(sprintf("\n%d of %d specifications show a delta R² whose 95%% CI excludes zero.\n",
            n_conclusive, n_total))
if (n_conclusive == 0) {
  cat("The conclusion is unchanged across every specification tried.\n")
} else {
  cat("At least one specification dissents. Before treating it as a finding:\n")
  cat("  - how many specifications were run in total? (that is this table)\n")
  cat("  - does it change the eligible genes, so it is a different sample?\n")
  cat("  - is there a mechanism specific to it (informative missingness for\n")
  cat("    keep_missing, length confounding for with_raw_mfe) that explains it\n")
  cat("    without structure carrying information?\n")
}


# ----------------------------- Figure ---------------------------------------

metric_label <- c(rsq_trad = "R² (held-out)", rmse = "RMSE", mae = "MAE")

p <- done |>
  filter(metric %in% names(metric_label)) |>
  mutate(metric_lab = factor(metric_label[metric], levels = metric_label),
         variant_lab = fct_rev(factor(sprintf("%s\n(n = %s)", variant,
                                              format(n_test, big.mark = ","))))) |>
  ggplot(aes(delta, variant_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, linewidth = 0.8, colour = "grey25") +
  geom_point(aes(fill = conclusive), size = 4, shape = 21, colour = "grey20",
             show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#E8590C", `FALSE` = "grey70")) +
  facet_wrap(~ metric_lab, scales = "free_x") +
  labs(
    title = "Sensitivity: does the conclusion survive every specification?",
    subtitle = paste0("Every variant in the registry, reported together. ",
                      "Grey points span zero.\n",
                      "n is the held-out gene count, which differs when a variant ",
                      "changes row eligibility — so the\ndeltas are comparable in ",
                      "DIRECTION but the absolute scores behind them are not."),
    x = "Delta (Model B - Model A), PC1 units", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

for (ext in c("jpg", "pdf")) {
  args <- list(filename = file.path(PLOT_DIR, paste0("xgb_structure_variant_summary.", ext)),
               plot = p, width = 280, height = 60 + 30 * n_total, units = "mm", dpi = 300)
  if (ext == "pdf") args$device <- grDevices::cairo_pdf
  do.call(ggsave, args)
}

message("\nSummary written to ", TABLE_DIR, " and ", PLOT_DIR)
