# =============================================================================
# SECONDARY: does icSHAPE structural Gini add anything, on the genes that have it?
# =============================================================================
# This is NOT the headline result. Read xgb_structure_comparison.R first.
#
# WHY IT IS SEPARATE. The icSHAPE Gini columns are 80-91% missing on the human
# v9 cache. Only 861 of 13,601 modellable genes are complete on baseline +
# structure + Gini. Putting Gini into the main comparison would either cost 93%
# of the corpus, or require XGBoost to learn a default split direction for a
# feature absent in seven of every eight genes. Neither is a good trade for the
# primary question, so Gini gets its own run on its own genes, clearly labelled.
#
# WHY CROSS-VALIDATION HERE AND A HOLDOUT THERE. The committed `test` split
# contains only 98 Gini-complete genes — too few to put a usable interval on a
# paired difference. So this run uses 5-fold cross-validation, blocked on
# family_id_medium, over all 861, giving every gene one out-of-fold prediction
# from each model. Hyperparameters are re-tuned inside each outer training set,
# so no fold's assessment genes influence its own model.
#
# THREE CAVEATS THAT MUST TRAVEL WITH ANY NUMBER FROM THIS FILE.
#
# 1. SELECTION. icSHAPE coverage is not a random sample of the transcriptome —
#    probing depth tracks expression, so these 861 genes are biased toward
#    abundant transcripts. A result here describes well-probed genes, not the
#    corpus. n = 861 is also small, so intervals are wide.
#
# 2. GINI IS MEASURED, NOT PREDICTED. Every other feature in this project is
#    computed from sequence. icSHAPE Gini is an experimental readout, so a
#    model containing it is not a sequence-to-half-life predictor: you cannot
#    evaluate it on a transcript nobody has probed. That is a different kind of
#    model from Models A and B, and the comparison should be described as
#    "does measured structure carry information", not "does structure improve
#    prediction from sequence".
#
# 3. CONFOUNDING WITH ABUNDANCE. Gini is computed from probing reads, read
#    depth tracks expression, and expression is associated with half-life. So
#    an apparent Gini effect has a live alternative explanation — that it is
#    partly an abundance proxy — which this design cannot rule out. Raise it
#    before someone else does. Testing it needs an expression covariate, which
#    the current cache does not carry.
#
# Three models, so the run answers two questions at once:
#   A  baseline                            (main comparison, replicated here)
#   B  baseline + structure                (main comparison, replicated here)
#   C  baseline + structure + icSHAPE Gini (the increment this file exists for)
#
# Usage:
#   Rscript analysis/models/xgb_structure_gini_subset.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(finetune)
  library(xgboost)
})

source("R/load_all.R")
source("analysis/models/xgb_structure_features.R")

tidymodels_prefer(quiet = TRUE)


# ----------------------------- Configuration --------------------------------

SEED      <- 42L
OUTER_V   <- 5L
INNER_V   <- 4L
N_BOOT    <- 2000L
N_WORKERS <- 4L
N_THREADS <- 3L
GRID_SIZE <- as.integer(Sys.getenv("XGB_GINI_GRID_SIZE", "20"))

RUN_DIR   <- file.path(OUTPUT_DIR, "xgb_structure")
PLOT_DIR  <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
for (d in c(RUN_DIR, PLOT_DIR, TABLE_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

METRICS <- metric_set(rsq_trad, rmse, mae)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  utils::flush.console()
}

future::plan(future::multisession, workers = N_WORKERS)


# ----------------------------- 1. Gini-complete subset ----------------------
# Built from the SAME block definitions as the main run (both source
# xgb_structure_features.R), then cut down further to complete cases on Gini.

log_msg("Building the Gini-complete subset")
el <- eligible_dataset("human")

BASELINE  <- el$baseline
STRUCTURE <- el$structure
GINI      <- el$gini

sub <- el$data
sub <- sub[stats::complete.cases(sub[, GINI, drop = FALSE]), , drop = FALSE]

PREDS <- list(
  `A: baseline`                   = BASELINE,
  `B: baseline + structure`       = c(BASELINE, STRUCTURE),
  `C: + icSHAPE Gini`             = c(BASELINE, STRUCTURE, GINI)
)

cat("\n=== icSHAPE secondary subset ===\n")
cat(sprintf("Genes with complete baseline + structure + Gini: %d of %d eligible\n",
            nrow(sub), nrow(el$data)))
cat(sprintf("  (%.1f%% of the main analysis set)\n", 100 * nrow(sub) / nrow(el$data)))
cat("Their distribution across the committed splits (unused here, shown for context):\n")
print(table(sub$split))
cat(sprintf("\nPredictor counts:  A %d   B %d   C %d (+%d Gini)\n",
            length(PREDS[[1]]), length(PREDS[[2]]), length(PREDS[[3]]), length(GINI)))
cat("Gini block: ", paste(GINI, collapse = ", "), "\n\n")

# Zero-variance guard on the smaller frame: a codon that is constant across 861
# genes but varied across 11,801 is entirely possible.
zv <- names(which(vapply(sub[unique(unlist(PREDS))],
                         function(x) length(unique(x)) < 2L, logical(1))))
if (length(zv)) {
  cat("Dropped (zero variance on this subset): ", paste(zv, collapse = ", "), "\n")
  PREDS <- lapply(PREDS, setdiff, y = zv)
}


# ----------------------------- 2. Outer folds -------------------------------
# One shared set, created once, reused by all three models — the same guarantee
# the main run gets from the committed split artefact.

set.seed(SEED)
outer_folds <- group_vfold_cv(sub,
                              group   = paste0("family_id_", BLOCK_LEVEL),
                              v       = OUTER_V,
                              balance = "observations")
saveRDS(outer_folds, file.path(RUN_DIR, "gini_outer_folds.rds"))


# ----------------------------- 3. Model machinery ---------------------------
# Same spec and same ranges as the main run, so the two are comparable.

xgb_spec <- boost_tree(
  trees = tune(), tree_depth = tune(), min_n = tune(), learn_rate = tune(),
  loss_reduction = tune(), sample_size = tune(), mtry = tune()
) |>
  # !! so the value, not the symbol, reaches the future workers — see the note
  # in xgb_structure_comparison.R.
  set_engine("xgboost", nthread = !!N_THREADS, counts = FALSE,
             tree_method = "hist", lambda = tune("lambda")) |>
  set_mode("regression")

build_workflow <- function(preds) {
  rec <- recipe(x = sub[0, c(TARGET_COL, "gene_id", preds)]) |>
    update_role(all_of(preds),      new_role = "predictor") |>
    update_role(all_of(TARGET_COL), new_role = "outcome") |>
    update_role(gene_id,            new_role = "id") |>
    step_zv(all_predictors())
  workflow() |> add_recipe(rec) |> add_model(xgb_spec)
}

xgb_params <- build_workflow(PREDS[[1]]) |>
  extract_parameter_set_dials() |>
  update(
    mtry           = mtry_prop(c(0.2, 1.0)),
    trees          = trees(c(300L, 1200L)),
    tree_depth     = tree_depth(c(2L, 8L)),
    min_n          = min_n(c(5L, 100L)),
    learn_rate     = learn_rate(c(-2.5, -1)),
    loss_reduction = loss_reduction(c(-4, 1)),
    sample_size    = sample_prop(c(0.4, 1.0)),
    lambda         = dials::new_quant_param(
      type = "double", range = c(-3, 2), inclusive = c(TRUE, TRUE),
      trans = scales::log10_trans(), label = c(lambda = "L2 regularisation"))
  )

set.seed(SEED)
xgb_grid <- grid_space_filling(xgb_params, size = GRID_SIZE)


# ----------------------------- 4. Nested CV ---------------------------------
# For each outer fold: split the training genes into family-blocked inner
# folds, tune there, refit the winner on the whole outer training set, predict
# the outer assessment genes. The assessment genes take no part in choosing
# their own model's hyperparameters.

run_outer_fold <- function(i) {
  spl   <- outer_folds$splits[[i]]
  tr    <- analysis(spl)
  te    <- assessment(spl)

  set.seed(SEED + i)
  inner <- group_vfold_cv(tr, group = paste0("family_id_", BLOCK_LEVEL),
                          v = INNER_V, balance = "observations")

  map_dfr(names(PREDS), function(nm) {
    wf <- build_workflow(PREDS[[nm]])
    set.seed(SEED + i)
    res <- tune_race_anova(
      wf, resamples = inner, grid = xgb_grid,
      metrics = metric_set(rmse),
      control = control_race(verbose_elim = FALSE, save_pred = FALSE,
                             burn_in = 3L, parallel_over = "everything")
    )
    fitted <- fit(finalize_workflow(wf, select_best(res, metric = "rmse")), tr)
    tibble(fold = i, model = nm, gene_id = te$gene_id,
           observed = te[[TARGET_COL]],
           predicted = predict(fitted, te)$.pred)
  })
}

log_msg("Nested CV: ", OUTER_V, " outer folds x ", length(PREDS), " models x ",
        GRID_SIZE, " configurations x ", INNER_V, " inner folds")

oof <- map_dfr(seq_len(nrow(outer_folds)), function(i) {
  t0 <- Sys.time()
  out <- run_outer_fold(i)
  log_msg("  outer fold ", i, "/", OUTER_V, " done in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  out
})

future::plan(future::sequential)

stopifnot(nrow(oof) == nrow(sub) * length(PREDS),
          !anyNA(oof$predicted))

oof_wide <- oof |>
  select(fold, gene_id, observed, model, predicted) |>
  pivot_wider(names_from = model, values_from = predicted)

write_csv(oof_wide, file.path(TABLE_DIR, "xgb_structure_gini_oof_predictions.csv"))


# ----------------------------- 5. Metrics and paired bootstrap --------------

pooled <- oof |>
  group_by(model) |>
  group_modify(~ METRICS(.x, truth = observed, estimate = predicted)) |>
  ungroup() |>
  select(model, metric = .metric, value = .estimate) |>
  pivot_wider(names_from = metric, values_from = value)

cat("\n=== Out-of-fold metrics on the", nrow(sub), "Gini-complete genes ===\n")
print(as.data.frame(pooled), row.names = FALSE, digits = 4)

score <- function(o, p, idx) {
  o <- o[idx]; p <- p[idx]
  c(rsq_trad = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
    rmse     = sqrt(mean((o - p)^2)),
    mae      = mean(abs(o - p)))
}

#' Paired bootstrap of one model against another over the same genes.
paired_boot <- function(lhs, rhs, label) {
  o  <- oof_wide$observed
  a  <- oof_wide[[rhs]]      # reference
  b  <- oof_wide[[lhs]]      # candidate
  n  <- length(o)
  set.seed(SEED)
  bt <- vapply(seq_len(N_BOOT), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    score(o, b, idx) - score(o, a, idx)
  }, numeric(3))
  pt <- score(o, b, seq_len(n)) - score(o, a, seq_len(n))
  tibble(comparison = label,
         metric     = names(pt),
         delta      = as.numeric(pt),
         ci_low     = apply(bt, 1, quantile, 0.025),
         ci_high    = apply(bt, 1, quantile, 0.975)) |>
    mutate(conclusive = !(ci_low <= 0 & ci_high >= 0))
}

delta_table <- bind_rows(
  paired_boot("B: baseline + structure", "A: baseline",
              "structure vs baseline (replicated on subset)"),
  paired_boot("C: + icSHAPE Gini", "B: baseline + structure",
              "icSHAPE Gini increment"),
  paired_boot("C: + icSHAPE Gini", "A: baseline",
              "structure + Gini vs baseline")
) |>
  mutate(favours_candidate = if_else(metric == "rsq_trad",
                                     "delta > 0", "delta < 0"))

cat("\n=== Paired bootstrap,", N_BOOT, "replicates, 95% percentile CI ===\n")
print(as.data.frame(delta_table), row.names = FALSE, digits = 4)
write_csv(delta_table, file.path(TABLE_DIR, "xgb_structure_gini_deltas.csv"))


# ----------------------------- 6. Figure ------------------------------------
# Drawn by xgb_structure_plots.R from the CSV written above, like every other
# figure in this analysis. Rendered at the end of section 7.


# ----------------------------- 7. Manifest ----------------------------------

saveRDS(list(
  role          = "SECONDARY analysis; see xgb_structure_comparison.R for the primary",
  n_genes       = nrow(sub),
  n_eligible_main = nrow(el$data),
  design        = sprintf("%d-fold CV blocked on family_id_%s, %d-fold inner tuning",
                          OUTER_V, BLOCK_LEVEL, INNER_V),
  design_note   = paste("cross-validation rather than the committed holdout,",
                        "because `test` holds only ~98 Gini-complete genes"),
  caveats       = c(
    selection   = paste("icSHAPE coverage tracks expression; these genes are",
                        "biased toward abundant transcripts and do not represent",
                        "the corpus"),
    measured    = paste("Gini is an experimental readout, not computed from",
                        "sequence, so model C cannot score an unprobed transcript"),
    confounding = paste("Gini is derived from probing read depth, which tracks",
                        "abundance, which is associated with half-life; this design",
                        "cannot separate a structure effect from an abundance proxy")),
  predictors    = PREDS,
  gini_block    = GINI,
  seed          = SEED,
  grid_size     = GRID_SIZE,
  pooled        = pooled,
  delta_table   = delta_table,
  fitted_at     = Sys.time()
), file.path(RUN_DIR, "gini_run_manifest.rds"))

log_msg("Secondary run complete. n = ", nrow(sub), " genes.")


# ----------------------------- 8. Figure ------------------------------------
# Sourced after the manifest, because the figure's subtitle reads its design
# string and gene count from it. See section 16 of the comparison script for
# why rendering lives in a separate file.

log_msg("Rendering figures")
source("analysis/models/xgb_structure_plots.R", local = new.env())
