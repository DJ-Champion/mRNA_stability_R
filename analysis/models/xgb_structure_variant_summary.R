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
# So this script does not take arguments and has no filters. It reports every
# variant in the registry, including any that disagree with the rest, and it
# additionally reports any run directory it finds on disk that the registry
# does NOT define — see "Orphans" below. Adding a variant to the registry and
# running it automatically adds rows here; there is no supported way to run a
# variant and leave it out of the table. If the answer to "how many did you
# try?" is not "all of the ones in this table", something has gone wrong.
#
# THE CLAIM TO MAKE is "the conclusion holds across every specification we
# tried", which is stronger than any single run. If one variant dissents, that
# is a finding about the specification and needs explaining — not promoting.
#
# HOW TO READ IT. Within a variant, the paired contrasts and their CIs are
# meaningful: every rung saw the same genes. BETWEEN variants, the absolute R²
# values are NOT comparable when the row policy differs, because the models are
# scored on different genes. That is why n is printed on every row.
#
# MULTIPLICITY. Each variant contributes several contrasts, only one of which
# is pre-specified. The headline table shows the primary contrast per variant;
# the secondary contrasts are written to the CSV and drawn on the figure, and
# are labelled as secondary in both. Counting "how many cleared zero" across
# secondary contrasts is not a meaningful summary and this script does not
# print one.
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
# defined but not yet run must appear as NOT RUN rather than silently vanish.
# A glob would let "delete the directory" work as a way to drop an inconvenient
# result, which is exactly what this file is meant to prevent.

rows <- map(names(VARIANTS), function(nm) {
  v      <- resolve_variant(nm)
  deltas <- file.path(variant_dir(v, "tables"), "xgb_structure_delta_bootstrap.csv")
  man    <- file.path(variant_dir(v, "root"), "run_manifest.rds")

  if (!file.exists(deltas) || !file.exists(man)) {
    return(tibble(variant = nm, label = v$label, status = "NOT RUN"))
  }

  m <- readRDS(man)
  read_csv(deltas, show_col_types = FALSE) |>
    mutate(variant = nm, label = v$label, status = "run",
           n_test = m$n_test, n_eligible = m$n_eligible,
           rows_policy = m$variant_rows,
           .before = 1)
}) |> list_rbind()

write_csv(rows, file.path(TABLE_DIR, "xgb_structure_variant_summary.csv"))

not_run <- rows |> filter(status == "NOT RUN") |> distinct(variant)
if (nrow(not_run) > 0) {
  message("\nDefined but not yet run: ", paste(not_run$variant, collapse = ", "))
  message("  Rscript analysis/models/xgb_structure_comparison.R <variant>")
}


# ----------------------------- Orphaned run directories ---------------------
# The registry drives the table, which closes the "delete a directory to hide a
# result" hole. It opens a smaller one in the other direction: RENAMING or
# removing a variant from the registry would leave its results on disk and
# quietly drop them from the table, which is the same failure wearing a
# different hat.
#
# So directories that exist but are not in the registry are reported loudly.
# They are usually legitimate — a superseded design whose numbers no longer
# apply — but "superseded" is a claim someone has to make out loud, not a thing
# that happens by a directory going unmentioned.

known   <- c(names(VARIANTS), "_summary", "gini")
on_disk <- list.dirs(ROOT, full.names = FALSE, recursive = FALSE)
orphans <- setdiff(on_disk, known)

if (length(orphans) > 0) {
  message("\n", strrep("!", 70))
  message("ORPHANED RUN DIRECTORIES — on disk but not in the variant registry:")
  for (o in orphans) {
    man <- file.path(ROOT, o, "run_manifest.rds")
    if (file.exists(man)) {
      m <- readRDS(man)
      message(sprintf("  %-16s run %s, n_test = %s", o,
                      format(m$fitted_at, "%Y-%m-%d"), m$n_test))
    } else {
      message(sprintf("  %-16s (no manifest)", o))
    }
  }
  message("These are NOT in the table above. Either re-register them or delete")
  message("them deliberately — do not leave results on disk unaccounted for.")
  message(strrep("!", 70))
}

writeLines(orphans, file.path(TABLE_DIR, "xgb_structure_orphan_dirs.txt"))

done <- rows |> filter(status == "run")
if (nrow(done) == 0) {
  stop("no variant has been run yet — run the default first", call. = FALSE)
}


# ----------------------------- Report ---------------------------------------

cat("\n=== Sensitivity across specifications: the PRIMARY contrast ===\n")
cat("Within a variant the contrast is paired and meaningful.\n")
cat("BETWEEN variants, absolute R² is not comparable when n differs.\n\n")

primary_rows <- done |> filter(primary, metric == "rsq_trad")

print(as.data.frame(
  primary_rows |>
    select(variant, rows_policy, n_eligible, n_test, contrast,
           rhs_value, lhs_value, delta, ci_low, ci_high, conclusive)),
  row.names = FALSE, digits = 4)

n_conclusive <- sum(primary_rows$conclusive)
n_total      <- nrow(primary_rows)

cat(sprintf("\n%d of %d specifications show a PRIMARY delta R² whose 95%% CI excludes zero.\n",
            n_conclusive, n_total))
if (n_conclusive == 0) {
  cat("The conclusion is unchanged across every specification tried.\n")
} else {
  cat("At least one specification dissents. Before treating it as a finding:\n")
  cat("  - how many specifications were run in total? (that is this table)\n")
  cat("  - does it change the eligible genes, so it is a different sample?\n")
  cat("  - is there a mechanism specific to it — informative missingness under\n")
  cat("    `default`, closed by construction under `complete_case` — that\n")
  cat("    explains it without structure carrying information?\n")
}

cat("\n=== Secondary contrasts (unadjusted; reported, not headlined) ===\n")
print(as.data.frame(
  done |>
    filter(!primary, metric == "rsq_trad") |>
    select(variant, kind, contrast, delta, ci_low, ci_high, conclusive)),
  row.names = FALSE, digits = 4)


# ----------------------------- Figure ---------------------------------------

metric_label <- c(rsq_trad = "R² (held-out)", rmse = "RMSE", mae = "MAE")

plot_data <- done |>
  filter(metric %in% names(metric_label)) |>
  mutate(metric_lab = factor(metric_label[metric], levels = metric_label),
         row_lab = fct_rev(factor(sprintf("%s%s\n%s (n = %s)",
                                          if_else(primary, "* ", "  "),
                                          contrast, variant,
                                          format(n_test, big.mark = ",")))))

p <- plot_data |>
  ggplot(aes(delta, row_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, linewidth = 0.8, colour = "grey25") +
  geom_point(aes(fill = conclusive, shape = primary), size = 4, colour = "grey20",
             show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#E8590C", `FALSE` = "grey70")) +
  scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21)) +
  facet_wrap(~ metric_lab, scales = "free_x") +
  labs(
    title = "Sensitivity: does the conclusion survive every specification?",
    subtitle = paste0("Every variant in the registry and every contrast in the ",
                      "ladder, reported together. Grey points span zero.\n",
                      "Diamonds marked * are the pre-specified primary contrast; ",
                      "circles are secondary and unadjusted.\n",
                      "n is the held-out gene count, which differs when a variant ",
                      "changes row eligibility — so the deltas\nare comparable in ",
                      "DIRECTION but the absolute scores behind them are not."),
    x = "Delta (larger model - smaller model), PC1 units", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_blank())

n_rows <- n_distinct(plot_data$row_lab)
for (ext in c("jpg", "pdf")) {
  args <- list(filename = file.path(PLOT_DIR, paste0("xgb_structure_variant_summary.", ext)),
               plot = p, width = 300, height = 70 + 16 * n_rows, units = "mm", dpi = 300)
  if (ext == "pdf") args$device <- grDevices::cairo_pdf
  do.call(ggsave, args)
}

message("\nSummary written to ", TABLE_DIR, " and ", PLOT_DIR)
