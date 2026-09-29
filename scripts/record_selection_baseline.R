# =============================================================================
# Record which columns every figure / job selects TODAY (SELECTION_PLAN.md, step 1)
# =============================================================================
# The reference for the selection consolidation: after each migration step,
# re-run this (or its equivalent for the new API) and diff against the CSV it
# wrote. The only intended difference is probing leaving the "core" figures.
#
# Two layers, one long CSV (scripts/selection_baseline.csv):
#
#   stage = "selected"  Columns the job's selection arguments resolve to, before
#                       any ranking. Deterministic: a function of the arguments
#                       and the columns of the dataset, not of the data values.
#                       Region filters (region heatmap `regions`) are NOT applied.
#   stage = "plotted"   What the job's plot function returns after top_n_per_group
#                       / top_n / min_n filters, one row per (feature, stem,
#                       region). Depends on the data. Only for jobs with a
#                       ranking step.
#
# Job definitions are copied from each script's runner, because the runners
# are not callable. If a runner changes, change it here and say so in the
# commit; the whole point is that this file records the PRE-migration state.
#
# NOTE (step 4): this script calls the PRE-migration API (groups / pick /
# keep_supergroups / top_n_per_group) and so runs only at commit 53f1a3e's
# parent state of the four migrated scripts. selection_baseline.csv is the
# record; scripts/diff_selection_baseline.R compares the migrated scripts to it.
#
# Usage (project root; needs data/cache/human_dataset_v*.rds):
#   Rscript scripts/record_selection_baseline.R
# =============================================================================

suppressMessages(source("R/load_all.R"))

OUT_FILE <- "scripts/selection_baseline.csv"

# Load a script's functions without running its runner (the guard is
# `sys.nframe() == 0 || identical(environment(), globalenv())`).
load_fns <- function(path) {
  e <- new.env(parent = globalenv())
  suppressMessages(sys.source(path, envir = e, keep.source = FALSE))
  e
}

df_full <- suppressMessages(build_dataset("human"))                          # dotplot, scatter, heatmaps, ...
df_utr  <- suppressMessages(build_dataset("human", min_utr = MIN_UTR_LENGTH))  # same call ranked/bands make
message("human: ", nrow(df_full), " rows, ", ncol(df_full), " columns")

rows <- list()
add <- function(script, job, stage, feature_id, column, stem = NA, region = NA) {
  rows[[length(rows) + 1L]] <<- data.frame(
    script = script, job = job, stage = stage, feature_id = feature_id,
    column = column, metric_stem = stem, region = region,
    stringsAsFactors = FALSE)
}

# --- Layer 1: selection as resolve_selection() + refine_group_columns() ------
record_selected <- function(script, job, df, groups, pick = list(), drop = list()) {
  sel <- resolve_selection(groups, pick, drop)
  n <- 0L
  for (g in sel$groups) {
    cols <- refine_group_columns(fg_columns(df, g), sel$pick[[g]], sel$drop[[g]])
    for (co in cols) add(script, job, "selected", g, co)
    n <- n + length(cols)
  }
  message(sprintf("  %-34s %-26s %4d columns", script, job, n))
}

# Selection with plain fg_columns() over a vector of feature ids.
record_features <- function(script, job, df, features) {
  n <- 0L
  for (g in features) {
    for (co in fg_columns(df, g)) add(script, job, "selected", g, co)
    n <- n + length(fg_columns(df, g))
  }
  message(sprintf("  %-34s %-26s %4d columns", script, job, n))
}

release <- list(codon_freqs = NULL, aa_freqs = NULL)
top2    <- list(codon_freqs = 2, aa_freqs = 2)

message("\nLayer 1: selected columns")

# ranked (4 jobs)
S <- "feature_correlation_ranked.R"
record_selected(S, "halflife",            df_utr, DEFAULT_PLOT_GROUPS, pick = release)
record_selected(S, "translation_efficiency", df_utr, DEFAULT_PLOT_GROUPS, pick = release)
record_selected(S, "halflife_aa_full",    df_utr, "aa_freqs")
record_selected(S, "halflife_codon_full", df_utr, "codon_freqs")

# bands (2 jobs, same selection)
S <- "feature_correlation_bands.R"
bands_groups <- c("nmd_core", "junction_core", "global_folding",
                  "local_folding", "sequence_select", "translation_core")
record_selected(S, "halflife",               df_utr, bands_groups, pick = release)
record_selected(S, "translation_efficiency", df_utr, bands_groups, pick = release)

# dotplot (4 jobs)
S <- "feature_correlation_dotplot.R"
record_selected(S, "halflife",               df_full, DEFAULT_PLOT_GROUPS)
record_selected(S, "translation_efficiency", df_full, DEFAULT_PLOT_GROUPS)
record_selected(S, "halflife_codon_aa",      df_full, c("codon_freqs", "aa_freqs"))
record_selected(S, "halflife_nuc_ratios",    df_full, c("nuc_ratios", "gc"))

# scatter (3 jobs)
S <- "feature_response_scatter.R"
for (j in c("included", "top_3", "top_1"))
  record_selected(S, j, df_full, DEFAULT_PLOT_GROUPS)

# region heatmap (2 jobs; `regions` filter not applied here)
S <- "region_feature_heatmap.R"
record_selected(S, "structure", df_full, "structure")
record_selected(S, "included",  df_full, DEFAULT_PLOT_GROUPS)

# correlation heatmap workflow
record_selected("correlation_heatmap_workflow.R", "workflow", df_full, DEFAULT_PLOT_GROUPS)

# hex panels
S <- "feature_response_hex_panels.R"
record_selected(S, "sequence",  df_full, "sequence")
record_selected(S, "structure", df_full, "structure")

# group panel sweep: groups = NULL -> expand_groups(NULL) minus DEFAULT_SWEEP_SKIP
sweep_skip <- c("codon_freqs", "aa_freqs")   # DEFAULT_SWEEP_SKIP in group_panel_sweep.R
record_features("group_panel_sweep.R", "sweep", df_full,
                setdiff(expand_groups(NULL), sweep_skip))

# feature-feature correlation table: .FF_DEFAULT_GROUPS
ff_default <- setdiff(EXPLORATORY_FEATURES, c("codon_freqs", "aa_freqs"))
record_features("feature_feature_correlation_table.R", "default", df_full, ff_default)

# QC overview: groups = NULL becomes names(FEATURE_PATTERNS), i.e. EVERY
# selectable feature, not just the exploratory ones (unlike the other NULLs).
record_selected("dataset_overview.R", "overview", df_full, names(FEATURE_PATTERNS))

# models: baseline / structure blocks (model flag), and probing (excluded from both)
record_features("xgb_structure_features.R", "baseline_columns", df_full,
                setdiff(MODEL_FEATURES, SUPERGROUPS$structure))
record_features("xgb_structure_features.R", "structure_columns", df_full,
                intersect(MODEL_FEATURES, SUPERGROUPS$structure))

# cross-species concordance: hard-wired probing
record_features("cross_species_probing_concordance.R", "probing", df_full, "probing")

# example_analysis.R uses direct fg("...") calls; not selection, not recorded.

# --- Layer 2: what the plot functions keep after ranking ----------------------
message("\nLayer 2: plotted stems (runs the correlations)")
ranked  <- load_fns("analysis/correlations/feature_correlation_ranked.R")
dotplot <- load_fns("analysis/correlations/feature_correlation_dotplot.R")
scatter <- load_fns("analysis/correlations/feature_response_scatter.R")

record_plotted <- function(script, job, tab) {
  if (!"group" %in% names(tab)) tab$group <- NA_character_
  if (!"region" %in% names(tab)) tab$region <- NA_character_
  stem <- if ("metric_stem" %in% names(tab)) tab$metric_stem else tab$variable
  u <- unique(data.frame(group = tab$group, stem = stem, region = tab$region,
                         stringsAsFactors = FALSE))
  for (i in seq_len(nrow(u)))
    add(script, job, "plotted", u$group[i], NA_character_, u$stem[i], u$region[i])
  message(sprintf("  %-34s %-26s %4d (stem, region) rows", script, job, nrow(u)))
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# ranked
for (j in list(
  list("halflife", "halflife", DEFAULT_PLOT_GROUPS, release, top2),
  list("translation_efficiency", "translation_efficiency", DEFAULT_PLOT_GROUPS, release, top2),
  list("halflife", "halflife_aa_full", "aa_freqs", list(), list()),
  list("halflife", "halflife_codon_full", "codon_freqs", list(), list()))) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = j[[1]], groups = j[[3]], pick = j[[4]],
    orientation = "horizontal", sig_threshold = "auto", top_n_per_group = j[[5]],
    keep_supergroups = if (j[[2]] %in% c("halflife", "translation_efficiency"))
      c("structure", "sequence") else NULL))
  record_plotted("feature_correlation_ranked.R", j[[2]], out$table)
}

# bands: ranked with keep_supergroups = NULL, horizontal, its own groups
for (resp in c("halflife", "translation_efficiency")) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = resp, groups = bands_groups, pick = release,
    keep_supergroups = NULL, orientation = "horizontal", top_n_per_group = top2))
  record_plotted("feature_correlation_bands.R", resp, out$table)
}

# dotplot
for (j in list(
  list("halflife", "halflife", DEFAULT_PLOT_GROUPS, top2),
  list("translation_efficiency", "translation_efficiency", DEFAULT_PLOT_GROUPS, top2),
  list("halflife", "halflife_codon_aa", c("codon_freqs", "aa_freqs"),
       list(codon_freqs = 15, aa_freqs = 20)),
  list("halflife", "halflife_nuc_ratios", c("nuc_ratios", "gc"), list()))) {
  out <- quiet(dotplot$feature_correlation_dotplot(
    df_full, response = j[[1]], sig_threshold = 0.02, groups = j[[3]],
    absolute = TRUE, top_n_per_group = j[[4]]))
  record_plotted("feature_correlation_dotplot.R", j[[2]], out$table)
}

# scatter
for (j in list(
  list("included", list(top_n_per_group = top2, noise_filter = 0, label_quantile = 0.8)),
  list("top_3",    list(top_n = 3, label_quantile = 0.3)),
  list("top_1",    list(top_n = 1, label_quantile = 0)))) {
  out <- quiet(do.call(scatter$feature_response_scatter,
                       c(list(df_full, groups = DEFAULT_PLOT_GROUPS), j[[2]])))
  record_plotted("feature_response_scatter.R", j[[1]], out$table)
}

# --- Write --------------------------------------------------------------------
res <- do.call(rbind, rows)
res <- res[order(res$script, res$job, res$stage, res$feature_id, res$column,
                 res$metric_stem, res$region), ]
write.csv(res, OUT_FILE, row.names = FALSE, na = "")
message("\nWrote ", nrow(res), " rows to ", OUT_FILE)
