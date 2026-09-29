# =============================================================================
# Compare the migrated scripts to the step-1 baseline
# =============================================================================
# Part 1 (step 4): runs the ranked, bands, dotplot and scatter functions with their new
# include / exclude / top_n arguments, exactly as their runners call them, and
# diffs the (stem, region) sets they plot against scripts/selection_baseline.csv.
#
# Part 2 (step 5): compares the columns the other migrated scripts select
# (region heatmap, correlation workflow, hex panels, sweep, feature-feature
# table, QC overview, xgb blocks) with the baseline "selected" rows.
#
# Expected differences (SELECTION_PLAN.md):
#   * probing (the four icSHAPE Gini columns) is gone from every "core" job;
#   * dotplot and scatter "core" jobs used to pin two hand-named codons and two
#     amino acids (bundle sequence_select); they now show the true top two, as
#     ranked and bands always did.
#   * the QC overview covers exploratory features only (owner's decision), so
#     it loses exactly the columns of features the table excludes.
# Anything else is a regression.
#
# Usage (project root):  Rscript scripts/diff_selection_baseline.R
# =============================================================================

suppressMessages(source("R/load_all.R"))

load_fns <- function(path) {
  e <- new.env(parent = globalenv())
  suppressMessages(sys.source(path, envir = e, keep.source = FALSE))
  e
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

df_full <- suppressMessages(build_dataset("human"))
df_utr  <- suppressMessages(build_dataset("human", min_utr = MIN_UTR_LENGTH))
base    <- read.csv("scripts/selection_baseline.csv", stringsAsFactors = FALSE)

ranked  <- load_fns("analysis/correlations/feature_correlation_ranked.R")
dotplot <- load_fns("analysis/correlations/feature_correlation_dotplot.R")
scatter <- load_fns("analysis/correlations/feature_response_scatter.R")

top2 <- list(codon_freqs = 2, aa_freqs = 2)
now  <- list()   # key "script | job" -> character vector "feature_id | stem | region"
key3 <- function(g, stem, region) paste(g, stem, region, sep = " | ")
record <- function(script, job, tab) {
  stem <- if ("metric_stem" %in% names(tab)) tab$metric_stem else tab$variable
  now[[paste(script, job, sep = " | ")]] <<-
    unique(key3(tab$group, stem, tab$region))
}

for (j in list(
  list("halflife", "halflife", "core", top2),
  list("translation_efficiency", "translation_efficiency", "core", top2),
  list("halflife", "halflife_aa_full", "aa_freqs", NULL),
  list("halflife", "halflife_codon_full", "codon_freqs", NULL))) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = j[[1]], include = j[[3]], orientation = "horizontal",
    sig_threshold = "auto", top_n = j[[4]]))
  record("feature_correlation_ranked.R", j[[2]], out$table)
}
for (resp in c("halflife", "translation_efficiency")) {
  out <- quiet(ranked$feature_correlation_ranked(
    df_utr, response = resp, include = "core", orientation = "horizontal",
    top_n = top2))
  record("feature_correlation_bands.R", resp, out$table)
}
for (j in list(
  list("halflife", "halflife", "core", top2),
  list("translation_efficiency", "translation_efficiency", "core", top2),
  list("halflife", "halflife_codon_aa", c("codon_freqs", "aa_freqs"),
       list(codon_freqs = 15, aa_freqs = 20)),
  list("halflife", "halflife_nuc_ratios", c("nuc_ratios", "gc"), NULL))) {
  out <- quiet(dotplot$feature_correlation_dotplot(
    df_full, response = j[[1]], sig_threshold = 0.02, include = j[[3]],
    absolute = TRUE, top_n = j[[4]]))
  record("feature_correlation_dotplot.R", j[[2]], out$table)
}
for (j in list(
  list("included", list(top_n = top2, noise_filter = 0, label_quantile = 0.8)),
  list("top_3",    list(max_features = 3, label_quantile = 0.3)),
  list("top_1",    list(max_features = 1, label_quantile = 0)))) {
  out <- quiet(do.call(scatter$feature_response_scatter, c(list(df_full), j[[2]])))
  # the scatter table has region "none" for region-less columns; the baseline
  # recorded the same table, so compare like with like
  record("feature_response_scatter.R", j[[1]], out$table)
}

# Baseline sets, same key shape.
bp <- base[base$stage == "plotted", ]
was <- split(key3(bp$feature_id, bp$metric_stem, bp$region),
             paste(bp$script, bp$job, sep = " | "))

n_unexpected <- 0L
for (k in names(now)) {
  gone  <- setdiff(was[[k]], now[[k]])
  added <- setdiff(now[[k]], was[[k]])
  probing <- grepl("^probing \\|", gone)
  pinned  <- grepl("^(codon_freqs|aa_freqs) \\|", c(gone, added))
  unexpected <- c(gone[!probing & !grepl("^(codon_freqs|aa_freqs) \\|", gone)],
                  added[!grepl("^(codon_freqs|aa_freqs) \\|", added)])
  cat(sprintf("%-58s was %3d, now %3d | -probing %d | codon/aa changed %d | UNEXPECTED %d\n",
              k, length(was[[k]]), length(now[[k]]), sum(probing),
              sum(pinned), length(unexpected)))
  for (u in unexpected) cat("    ", u, "\n")
  n_unexpected <- n_unexpected + length(unexpected)
}
# =============================================================================
# Part 2: selected columns of the step-5 scripts
# =============================================================================
cat("\n--- step 5: selected columns ---\n")
bs <- base[base$stage == "selected", ]
was_cols <- split(bs$column, paste(bs$script, bs$job, sep = " | "))
was_feat <- split(bs$feature_id, paste(bs$script, bs$job, sep = " | "))
cols_of <- function(x) selected_columns(x)
now_cols <- list()

hm  <- load_fns("analysis/correlations/region_feature_heatmap.R")
wf  <- load_fns("analysis/correlations/correlation_heatmap_workflow.R")
hx  <- load_fns("analysis/correlations/feature_response_hex_panels.R")
ff  <- load_fns("analysis/correlations/feature_feature_correlation_table.R")
qc  <- load_fns("analysis/qc/dataset_overview.R")
xg  <- load_fns("analysis/models/xgb_structure_features.R")

# region heatmap: every feature that appears in any region's matrix
hm_core <- c("5utr", "cds", "3utr")   # the runner's regions for each job
hm_all  <- c("5utr", "cds", "3utr", "mrna", "start", "stop", "last100")
for (j in list(list("structure", "structure", NULL, hm_core),
               list("included", "core", top2, hm_all))) {
  out <- quiet(hm$region_feature_heatmap(
    df_full, response = "halflife", include = j[[2]], top_n = j[[3]],
    regions = j[[4]],
    output_dir = NULL))
  feats <- unique(unlist(lapply(out, function(o) c(o$table$feature_x, o$table$feature_y))))
  now_cols[[paste("region_feature_heatmap.R", j[[1]], sep = " | ")]] <- setdiff(feats, "halflife")
}

# workflow: the candidate set it builds (same call it makes)
now_cols[["correlation_heatmap_workflow.R | workflow"]] <-
  cols_of(select_features_v2(df_full, "core", top_n = top2, response = "halflife"))

# hex panels: the table keeps every selected feature
for (g in c("sequence", "structure")) {
  out <- quiet(hx$feature_response_hex_panels(df_full, include = g))
  now_cols[[paste("feature_response_hex_panels.R", g, sep = " | ")]] <- out$table$variable
}

# sweep: features -> columns
now_cols[["group_panel_sweep.R | sweep"]] <-
  unlist(lapply(selected_features("exploratory", c("codon_freqs", "aa_freqs")),
                fg_columns, df = df_full))

# feature-feature table
out <- quiet(ff$compute_feature_correlation_table(df_full))
now_cols[["feature_feature_correlation_table.R | default"]] <-
  unique(c(out$table$feature_a, out$table$feature_b))

# QC overview: n_columns per feature from the returned table
out <- quiet(qc$missingness_by_group_plot(df_full))
now_cols[["dataset_overview.R | overview"]] <-
  unlist(lapply(unique(as.character(out$table$group)), fg_columns, df = df_full))

# xgb blocks
now_cols[["xgb_structure_features.R | baseline_columns"]]  <- xg$baseline_columns(df_full)
now_cols[["xgb_structure_features.R | structure_columns"]] <- xg$structure_columns(df_full)

# The ff table drops non-numeric / response / id columns and the heatmaps drop
# columns without a region, so those compare against the baseline restricted
# to what they can show.
non_expl <- unlist(lapply(setdiff(names(FEATURE_PATTERNS), EXPLORATORY_FEATURES),
                          fg_columns, df = df_full))
for (k in names(now_cols)) {
  w <- unique(was_cols[[k]]); n <- unique(now_cols[[k]])
  if (grepl("^region_feature_heatmap", k)) {
    # the heatmaps draw only columns ending in one of the runner's regions
    hm_regions <- if (grepl("structure$", k)) hm_core else hm_all
    w <- w[sub("^.*_", "", w) %in% hm_regions]
  }
  gone <- setdiff(w, n); added <- setdiff(n, w)
  is_probing <- gone %in% fg_columns(df_full, "probing")
  is_pin     <- c(gone, added) %in% unlist(lapply(c("codon_freqs", "aa_freqs"),
                                                   fg_columns, df = df_full))
  # allowed: probing leaving core (region heatmap "included" only), pinned
  # codons/aa becoming the top two (heatmap "included", workflow), and the QC
  # overview shedding non-exploratory columns
  allow_probing <- k == "region_feature_heatmap.R | included" ||
                   k == "correlation_heatmap_workflow.R | workflow"
  allow_pin     <- allow_probing
  allow_qc      <- k == "dataset_overview.R | overview"
  bad_gone <- gone[!(allow_probing & is_probing) &
                   !(allow_pin & gone %in% unlist(lapply(c("codon_freqs", "aa_freqs"), fg_columns, df = df_full))) &
                   !(allow_qc & gone %in% non_expl)]
  bad_add  <- added[!(allow_pin & added %in% unlist(lapply(c("codon_freqs", "aa_freqs"), fg_columns, df = df_full)))]
  cat(sprintf("%-58s was %3d, now %3d | gone %3d, added %3d | UNEXPECTED %d\n",
              k, length(w), length(n), length(gone), length(added),
              length(bad_gone) + length(bad_add)))
  for (u in c(bad_gone, bad_add)) cat("    ", u, "\n")
  n_unexpected <- n_unexpected + length(bad_gone) + length(bad_add)
}

if (n_unexpected) { cat("\n", n_unexpected, " unexpected difference(s)\n", sep = ""); quit(status = 1) }
cat("\nOnly the expected differences.\n")
