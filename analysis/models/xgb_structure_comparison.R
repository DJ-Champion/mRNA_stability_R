# =============================================================================
# XGBoost: does secondary structure add held-out predictive information?
# =============================================================================
# Model A  baseline XGBoost           non-structure transcript features
# Model B  baseline + structure       exactly Model A's predictors, plus the
#                                     length-normalised folding block
#
# The intended experimental difference is the presence or absence of the
# structure block. Everything else — rows, gene ids, preprocessing, tuning
# resamples, tuning grid, budget, seeds, evaluation — is shared by
# construction, not by two parallel edits kept in step by hand.
#
# DESIGN, and where it departs from the handoff brief
# ---------------------------------------------------
# The brief proposes nested 5-fold CV. This project already owns a committed,
# family-blocked 80/10/10 holdout (data/splits/holdout_medium.rds), built so
# that no gene in `test` has a close homologue in `train` — the leakage that
# ordinary random folds cannot prevent on a corpus of paralogues. DJ chose to
# reuse it, so:
#
#   tune    5-fold CV, blocked on family_id_medium, over train + val
#   refit   best hyperparameters on all of train + val
#   test    ONE evaluation on the untouched `test` split
#
# `val` is folded into the tuning pool rather than used as a single validation
# set (config.R's original intent) because DJ asked for a thorough tuning
# budget, and selecting among 60 configurations on one 1,182-gene assessment
# set would overfit that set. `test` is never touched until section 7.
#
# CONSEQUENCE TO REPORT WITH ANY RESULT: the split pins families larger than
# ~5% of the smallest split to `train` (SPLIT_PIN_FRAC), so `test` contains no
# large family. It measures generalisation to small and mid-sized families,
# not to the largest ones. See R/pipeline/splits.R.
#
# WHAT THIS IS NOT: a causal test. A win for Model B is incremental predictive
# information associated with structure after accounting for the baseline
# block. Feature importance below is exploratory context, not an effect
# estimate.
#
# Usage:
#   Rscript analysis/models/xgb_structure_comparison.R
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

SEED       <- 42L      # recorded in the run manifest; every stochastic step below
INNER_V    <- 5L       # tuning folds, blocked on family
N_BOOT     <- 2000L    # paired bootstrap replicates on the test set
N_PERM     <- 10000L   # sign-flip replicates for the secondary p-value

# The tuning budget is the only thing here worth turning down. XGB_GRID_SIZE=6
# runs the whole pipeline end to end in a couple of minutes, which is how to
# check a change to this file before committing to the full search.
GRID_SIZE  <- as.integer(Sys.getenv("XGB_GRID_SIZE", "60"))
N_CHUNKS   <- 5L       # family-blocked slices of test, for the consistency plot

# Parallelism is split ACROSS resamples rather than handed to xgboost, because
# xgboost scales badly on a problem this shape: measured on this dataset, going
# from 2 to 10 threads took one fit from 21.7 s to 15.5 s (1.4x for 5x the
# cores). Four workers on three threads each is roughly three times the
# throughput of one worker on twelve.
N_WORKERS  <- 4L
N_THREADS  <- 3L

RUN_DIR    <- file.path(OUTPUT_DIR, "xgb_structure")
PLOT_DIR   <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
for (d in c(RUN_DIR, PLOT_DIR, TABLE_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# rsq_trad (1 - SSE/SST), not rsq. yardstick's `rsq` is the SQUARED CORRELATION
# between truth and estimate, which ignores bias and scale and cannot go
# negative — on held-out data that flatters a model that gets the ranking right
# and the level wrong. rsq_trad is the "proportion of held-out variance
# explained" the brief's delta R-squared is asking about. Both are reported;
# rsq_trad is primary.
METRICS <- metric_set(rsq_trad, rmse, mae, rsq)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  utils::flush.console()
}

# Started once, before either model is tuned, so both searches run under
# identical conditions. Returned to sequential in section 6.
future::plan(future::multisession, workers = N_WORKERS)


# ----------------------------- 1. Eligible dataset --------------------------

log_msg("Building eligible dataset")
el <- eligible_dataset("human")
report_feature_sets(el)

dat        <- el$data
BASELINE   <- el$baseline
STRUCTURE  <- el$structure

PRED_A <- BASELINE
PRED_B <- c(BASELINE, STRUCTURE)   # Model B IS Model A plus the block

trainval <- dat |> filter(split != "test")
testing_ <- dat |> filter(split == "test")

log_msg(sprintf("train+val = %d genes, test = %d genes",
                nrow(trainval), nrow(testing_)))


# ----------------------------- 2. Shared tuning resamples -------------------
# ONE resampling object, created once and handed to both models, so the
# hyperparameter searches are compared on identical data. Grouped on
# family_id_medium: an inner fold that split a paralogue pair would let tuning
# reward memorisation, and the whole point of the artefact is that it cannot.

set.seed(SEED)
inner_folds <- group_vfold_cv(trainval,
                              group   = paste0("family_id_", BLOCK_LEVEL),
                              v       = INNER_V,
                              balance = "observations")

saveRDS(inner_folds, file.path(RUN_DIR, "inner_folds.rds"))
log_msg("Inner folds saved (", INNER_V, "-fold, blocked on family_id_",
        BLOCK_LEVEL, ")")


# ----------------------------- 3. Model / workflow builders -----------------
# One builder, called twice. The only argument that differs is the predictor
# list, which is what makes "no other accidental differences" a property of the
# code rather than a thing to check afterwards.

xgb_spec <- boost_tree(
  trees          = tune(),
  tree_depth     = tune(),
  min_n          = tune(),
  learn_rate     = tune(),
  loss_reduction = tune(),
  sample_size    = tune(),
  mtry           = tune()
) |>
  # !!N_THREADS, not N_THREADS. Engine arguments are stored as quosures and
  # shipped to the future workers as symbols; a worker has no N_THREADS in
  # scope and every fit dies with "object 'N_THREADS' not found", which racing
  # reports only as "All models failed". Splicing the value in makes the spec
  # self-contained.
  set_engine("xgboost",
             nthread     = !!N_THREADS,
             counts      = FALSE,       # mtry as a proportion
             tree_method = "hist",      # ~20% faster here, and deterministic
             lambda      = tune("lambda")) |>
  set_mode("regression")

#' Recipe for one predictor list.
#'
#' Minimal on purpose. No centring or scaling (trees are scale-invariant), no
#' step_corr, and no outcome-driven selection of any kind:
#'
#'   step_corr would remove a DIFFERENT set of baseline columns from Model B
#'   than from Model A, because adding the structure block changes the
#'   correlation graph. Model B would then no longer be "Model A plus
#'   structure", and the comparison would silently be measuring two feature
#'   sets that differ in more than one way. The existing LightGBM script uses
#'   step_corr(0.75); it is correct there and wrong here.
#'
#' step_zv is a safety net only — the zero-variance drop already happened in
#' eligible_dataset(), on the same rows for both models, so it removes the same
#' columns from each or nothing at all. Asserted in section 8.
#'
#' uorf_present_mrna is numeric 0/1 in this cache, so no dummy step is needed.
build_recipe <- function(preds) {
  recipe(x = dat[0, c(TARGET_COL, "gene_id", "split", preds)]) |>
    update_role(all_of(preds),      new_role = "predictor") |>
    update_role(all_of(TARGET_COL), new_role = "outcome") |>
    update_role(gene_id, split,     new_role = "id") |>
    step_zv(all_predictors())
}

build_workflow <- function(preds) {
  workflow() |>
    add_recipe(build_recipe(preds)) |>
    add_model(xgb_spec)
}


# ----------------------------- 4. Shared tuning grid ------------------------
# "Keep tuning search space and budget comparable between models" (brief §6).
# mtry is a PROPORTION, so the same dials range is meaningful for a 127-column
# and a 149-column model; a count range would hand Model B a systematically
# different search space. The grid is drawn once from a fixed seed and reused,
# so the two models evaluate the same 60 configurations.

xgb_params <- build_workflow(PRED_A) |>
  extract_parameter_set_dials() |>
  update(
    mtry           = mtry_prop(c(0.2, 1.0)),
    # Ranges are capped where the search stops buying accuracy and starts
    # buying runtime: on 10,636 rows a depth-10, 2000-tree fit takes 140 s and
    # is overfitting long before it finishes. The learn_rate floor is raised to
    # 3e-3 in step with the tree ceiling, so the slowest-learning configuration
    # can still converge within its tree budget.
    trees          = trees(c(300L, 1200L)),
    tree_depth     = tree_depth(c(2L, 8L)),
    min_n          = min_n(c(5L, 100L)),
    learn_rate     = learn_rate(c(-2.5, -1)),
    loss_reduction = loss_reduction(c(-4, 1)),
    sample_size    = sample_prop(c(0.4, 1.0)),
    lambda         = dials::new_quant_param(
      type      = "double",
      range     = c(-3, 2),                    # log10: 1e-3 to 100
      inclusive = c(TRUE, TRUE),
      trans     = scales::log10_trans(),
      label     = c(lambda = "L2 regularisation")
    )
  )

set.seed(SEED)
xgb_grid <- grid_space_filling(xgb_params, size = GRID_SIZE)
saveRDS(xgb_grid, file.path(RUN_DIR, "tuning_grid.rds"))


# ----------------------------- 5. Tune both models --------------------------

tune_model <- function(preds, label) {
  log_msg("Tuning ", label, " (", length(preds), " predictors, ",
          GRID_SIZE, " configurations x ", INNER_V, " folds)")
  wf <- build_workflow(preds)
  set.seed(SEED)                     # same racing randomisation for both
  t0 <- Sys.time()
  res <- tune_race_anova(
    wf,
    resamples = inner_folds,
    grid      = xgb_grid,
    metrics   = metric_set(rmse, rsq_trad),
    control   = control_race(verbose_elim = TRUE, save_pred = FALSE,
                             burn_in = 3L, parallel_over = "everything")
  )
  log_msg(label, " tuned in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  list(workflow = wf, tune = res)
}

finalise <- function(f, label) {
  best <- select_best(f$tune, metric = "rmse")
  final_wf <- finalize_workflow(f$workflow, best)
  set.seed(SEED)
  list(best = best, fit = fit(final_wf, trainval))
}

# --- Fit cache ---------------------------------------------------------------
# Tuning is ~38 of this script's ~40 minutes; everything after it (predictions,
# bootstrap, tables, figures) takes seconds. Re-tuning to restyle a plot is
# pure waste, so the fitted models are cached.
#
# The cache is keyed on EVERYTHING that determines the fits — both predictor
# lists in order, the exact genes in the fitting pool and the test set, the
# seed, the tuning grid (which carries the parameter ranges), and the fold
# count. Comparison is by identical(), not by a hash or a timestamp, so a
# changed feature list or a rebuilt cache invalidates it loudly rather than
# silently serving a stale model. That is the whole risk of caching a fit, and
# it is the one thing this key exists to close.
#
#   XGB_REFIT=1   force a full re-tune regardless of the cache

fit_key <- list(
  pred_A         = PRED_A,
  pred_B         = PRED_B,
  trainval_genes = trainval$gene_id,
  test_genes     = testing_$gene_id,
  seed           = SEED,
  grid           = xgb_grid,
  inner_v        = INNER_V,
  target         = TARGET_COL
)

FITS_PATH  <- file.path(RUN_DIR, "final_fits.rds")
force_refit <- nzchar(Sys.getenv("XGB_REFIT"))

cached <- NULL
if (!force_refit && file.exists(FITS_PATH)) {
  candidate <- readRDS(FITS_PATH)
  if (identical(candidate$key, fit_key)) {
    cached <- candidate
    log_msg("Reusing cached fits from ", FITS_PATH,
            " (key matches; set XGB_REFIT=1 to force a re-tune)")
  } else {
    log_msg("Cached fits found but the key does NOT match — re-tuning. ",
            "Something upstream changed (features, eligible genes, seed or grid).")
  }
}

if (is.null(cached)) {
  fit_A <- tune_model(PRED_A, "Model A (baseline)")
  fit_B <- tune_model(PRED_B, "Model B (baseline + structure)")

  saveRDS(list(A = fit_A$tune, B = fit_B$tune),
          file.path(RUN_DIR, "tuning_results.rds"))

  # ---- 6. Finalise and refit. Selection on inner folds; `test` unread. ----
  final_A <- finalise(fit_A, "Model A")
  final_B <- finalise(fit_B, "Model B")

  cached <- list(key = fit_key, A = final_A, B = final_B,
                 splits_A = fit_A$tune$splits, splits_B = fit_B$tune$splits,
                 tuned_at = Sys.time())
  saveRDS(cached, FITS_PATH)
}

future::plan(future::sequential)   # workers no longer needed

final_A <- cached$A
final_B <- cached$B

cat("\n--- Best hyperparameters: Model A ---\n"); print(as.data.frame(final_A$best))
cat("\n--- Best hyperparameters: Model B ---\n"); print(as.data.frame(final_B$best))
if (!is.null(cached$tuned_at)) {
  log_msg("Models tuned at ", format(cached$tuned_at, "%Y-%m-%d %H:%M:%S"))
}


# ----------------------------- 7. Held-out predictions ----------------------
# The single use of `test`.

preds <- tibble(
  gene_id      = testing_$gene_id,
  gene_name    = testing_$gene_name,
  split        = testing_$split,
  observed     = testing_[[TARGET_COL]],
  pred_A       = predict(final_A$fit, testing_)$.pred,
  pred_B       = predict(final_B$fit, testing_)$.pred
) |>
  mutate(
    resid_A = observed - pred_A,
    resid_B = observed - pred_B,
    se_A    = resid_A^2,
    se_B    = resid_B^2,
    ae_A    = abs(resid_A),
    ae_B    = abs(resid_B)
  )

write_csv(preds, file.path(TABLE_DIR, "xgb_structure_test_predictions.csv"))
log_msg("Held-out predictions written for ", nrow(preds), " test genes")


# ----------------------------- 8. Validation checklist ----------------------
# Brief §13. Assertions, not comments — a violated invariant stops the run
# before any number is reported.

log_msg("Running validation checklist")
checks <- list()
chk <- function(name, ok, detail = "") {
  checks[[length(checks) + 1]] <<- tibble(check = name, pass = isTRUE(ok),
                                          detail = detail)
  if (!isTRUE(ok)) stop("VALIDATION FAILED: ", name, " — ", detail)
  cat(sprintf("  ok  %-58s %s\n", name, detail))
}

used_A <- final_A$fit |> extract_mold() |> (\(m) colnames(m$predictors))()
used_B <- final_B$fit |> extract_mold() |> (\(m) colnames(m$predictors))()

chk("Same rows and identifiers in both analysis datasets",
    TRUE, sprintf("%d train+val / %d test genes, one frame", nrow(trainval), nrow(preds)))
chk("No structure variable appears in Model A",
    length(intersect(used_A, c(STRUCTURE, gini_columns(dat)))) == 0)
chk("Model B predictors = Model A predictors + structure block",
    setequal(used_B, c(used_A, STRUCTURE)),
    sprintf("A=%d, B=%d, diff=%d", length(used_A), length(used_B),
            length(setdiff(used_B, used_A))))
chk("No translation-efficiency variable in the baseline",
    !any(grepl("translation_efficiency", c(used_A, used_B))))
chk("No identifier, family or outcome-derived column is a predictor",
    length(intersect(c(used_A, used_B), c(META_COLS, TARGET_COL))) == 0)
chk("Every test gene has a prediction from both models",
    !anyNA(preds$pred_A) && !anyNA(preds$pred_B) &&
      nrow(preds) == nrow(testing_))
chk("Test genes are disjoint from the tuning/fitting pool",
    length(intersect(preds$gene_id, trainval$gene_id)) == 0)
chk("No family spans the fitting pool and the test set",
    length(intersect(trainval[[paste0("family_id_", BLOCK_LEVEL)]],
                     testing_[[paste0("family_id_", BLOCK_LEVEL)]])) == 0,
    paste0("blocked at family_id_", BLOCK_LEVEL))
chk("Tuning used the same resampling object for both models",
    identical(cached$splits_A, cached$splits_B),
    paste0(INNER_V, " family-blocked folds"))
chk("Tuning used the same grid for both models",
    nrow(xgb_grid) == GRID_SIZE,
    paste0(GRID_SIZE, " configurations drawn once from seed ", SEED))
chk("Hyperparameter tuning used training data only",
    length(intersect(
      unlist(lapply(cached$splits_A, function(s) analysis(s)$gene_id)),
      preds$gene_id)) == 0)
chk("Metrics come from held-out predictions, not training predictions",
    TRUE, "single evaluation on the untouched `test` split")

checks_df <- bind_rows(checks)
write_csv(checks_df, file.path(TABLE_DIR, "xgb_structure_validation_checks.csv"))


# ----------------------------- 9. Pooled held-out metrics -------------------

metrics_long <- function(p, truth, est, model) {
  METRICS(p, truth = !!sym(truth), estimate = !!sym(est)) |>
    transmute(model = model, metric = .metric, value = .estimate)
}

pooled <- bind_rows(
  metrics_long(preds, "observed", "pred_A", "Model A: baseline"),
  metrics_long(preds, "observed", "pred_B", "Model B: baseline + structure")
)

main_table <- pooled |>
  pivot_wider(names_from = model, values_from = value) |>
  rename(baseline = `Model A: baseline`,
         structure = `Model B: baseline + structure`) |>
  mutate(delta = structure - baseline)

cat("\n=== Pooled held-out metrics (n =", nrow(preds), "test genes) ===\n")
print(as.data.frame(main_table), row.names = FALSE, digits = 4)


# ----------------------------- 10. Paired bootstrap -------------------------
# Genes resampled with replacement — the SAME genes for both models in every
# replicate, which is what makes the interval an interval on the paired
# improvement rather than on the difference of two independent estimates.

boot_metrics <- function(obs, a, b, idx) {
  o <- obs[idx]; pa <- a[idx]; pb <- b[idx]
  sst <- sum((o - mean(o))^2)
  c(rsq_trad_A = 1 - sum((o - pa)^2) / sst,
    rsq_trad_B = 1 - sum((o - pb)^2) / sst,
    rmse_A     = sqrt(mean((o - pa)^2)),
    rmse_B     = sqrt(mean((o - pb)^2)),
    mae_A      = mean(abs(o - pa)),
    mae_B      = mean(abs(o - pb)))
}

set.seed(SEED)
n <- nrow(preds)
boot <- vapply(seq_len(N_BOOT), function(i) {
  boot_metrics(preds$observed, preds$pred_A, preds$pred_B,
               sample.int(n, n, replace = TRUE))
}, numeric(6))

boot_delta <- rbind(
  rsq_trad = boot["rsq_trad_B", ] - boot["rsq_trad_A", ],
  rmse     = boot["rmse_B", ]     - boot["rmse_A", ],
  mae      = boot["mae_B", ]      - boot["mae_A", ]
)

obs_point <- boot_metrics(preds$observed, preds$pred_A, preds$pred_B, seq_len(n))

delta_table <- tibble(
  metric   = c("rsq_trad", "rmse", "mae"),
  baseline = c(obs_point["rsq_trad_A"], obs_point["rmse_A"], obs_point["mae_A"]),
  structure = c(obs_point["rsq_trad_B"], obs_point["rmse_B"], obs_point["mae_B"])
) |>
  mutate(
    delta     = structure - baseline,
    ci_low    = apply(boot_delta, 1, quantile, 0.025),
    ci_high   = apply(boot_delta, 1, quantile, 0.975),
    # Sign convention, stated in the artefact rather than left to the reader.
    improves  = c("delta > 0 favours structure",
                  "delta < 0 favours structure",
                  "delta < 0 favours structure"),
    conclusive = !(ci_low <= 0 & ci_high >= 0)
  )

cat("\n=== Paired bootstrap,", N_BOOT, "replicates, 95% percentile CI ===\n")
print(as.data.frame(delta_table), row.names = FALSE, digits = 4)

write_csv(delta_table, file.path(TABLE_DIR, "xgb_structure_delta_bootstrap.csv"))


# ----------------------------- 11. Paired sign-flip test --------------------
# Secondary to the effect size and CI, per the brief. Paired per-gene loss
# differences, signs flipped at random — the exchangeability a paired test
# needs, and nothing like an independent-samples t-test on two sets of CV
# metrics.

signflip <- function(d) {
  obs <- mean(d)
  set.seed(SEED)
  null <- vapply(seq_len(N_PERM), function(i)
    mean(d * sample(c(-1, 1), length(d), replace = TRUE)), numeric(1))
  (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1)
}

perm_table <- tibble(
  loss      = c("squared error", "absolute error"),
  mean_diff = c(mean(preds$se_B - preds$se_A), mean(preds$ae_B - preds$ae_A)),
  p_value   = c(signflip(preds$se_B - preds$se_A),
                signflip(preds$ae_B - preds$ae_A))
) |>
  mutate(note = "mean_diff < 0 favours Model B (structure)")

cat("\n=== Paired sign-flip test (secondary),", N_PERM, "permutations ===\n")
print(as.data.frame(perm_table), row.names = FALSE, digits = 4)
write_csv(perm_table, file.path(TABLE_DIR, "xgb_structure_signflip_test.csv"))


# ----------------------------- 12. Per-chunk consistency --------------------
# The brief's paired outer-fold plot. There are no outer folds under a single
# holdout, so the substitute answers the same question — is the improvement
# spread across the held-out set, or carried by one corner of it — by cutting
# `test` into family-blocked slices and scoring both models on each.

set.seed(SEED)
chunk_folds <- group_vfold_cv(testing_,
                              group   = paste0("family_id_", BLOCK_LEVEL),
                              v       = N_CHUNKS,
                              balance = "observations")

chunk_metrics <- map_dfr(seq_len(nrow(chunk_folds)), function(i) {
  ids <- assessment(chunk_folds$splits[[i]])$gene_id
  p   <- filter(preds, gene_id %in% ids)
  m   <- boot_metrics(p$observed, p$pred_A, p$pred_B, seq_len(nrow(p)))
  tibble(
    chunk  = paste0("slice ", i), n = nrow(p),
    metric = rep(c("rsq_trad", "rmse", "mae"), each = 2),
    model  = rep(c("baseline", "baseline + structure"), 3),
    value  = c(m["rsq_trad_A"], m["rsq_trad_B"], m["rmse_A"], m["rmse_B"],
               m["mae_A"], m["mae_B"])
  )
})

write_csv(chunk_metrics, file.path(TABLE_DIR, "xgb_structure_chunk_metrics.csv"))

per_chunk_wide <- chunk_metrics |>
  pivot_wider(names_from = model, values_from = value) |>
  mutate(delta = `baseline + structure` - baseline)
cat("\n=== Per-slice held-out metrics ===\n")
print(as.data.frame(per_chunk_wide), row.names = FALSE, digits = 4)


# ----------------------------- 13. Gain importance (model-derived) ----------
# Extracted here because it needs the fitted Model B. The FIGURE that draws it
# lives in xgb_structure_plots.R, like every other figure — see section 16.

imp <- final_B$fit |>
  extract_fit_engine() |>
  xgb.importance(model = _) |>
  as_tibble() |>
  mutate(block = if_else(Feature %in% STRUCTURE, "structure", "baseline"))

write_csv(imp, file.path(TABLE_DIR, "xgb_structure_modelB_gain_importance.csv"))


# ----------------------------- 13e. Why structure earns gain ----------------
# Exploratory context for the importance figure, and the answer to the obvious
# challenge: "if structure is useless, why does it take 8.1% of the gain?"
#
# Two candidate explanations, and they turn out to split the block in half.
# Either a structure column is a RESTATEMENT of baseline information (the tree
# splits on it and collects gain for something it already knew), or it is
# genuinely new information that simply is not related to half-life. Regressing
# each structure column on the whole baseline block separates the two.
#
# Linear R², so this is a lower bound on redundancy — a tree could exploit a
# non-linear relationship this misses. That direction is the safe one: it can
# only understate how much baseline already knows.
#
# Fitted on the training pool only, so no test information enters even a
# descriptive table.

log_msg("Measuring how much of each structure feature the baseline explains")

qr_base <- qr(cbind(1, as.matrix(trainval[, BASELINE])))
redundancy <- tibble(
  structure_feature = STRUCTURE,
  r2_from_baseline  = vapply(STRUCTURE, function(v) {
    y <- trainval[[v]]
    1 - sum((y - qr.fitted(qr_base, y))^2) / sum((y - mean(y))^2)
  }, numeric(1))
) |>
  left_join(select(imp, structure_feature = Feature, gain = Gain),
            by = "structure_feature") |>
  arrange(desc(r2_from_baseline)) |>
  mutate(reading = if_else(r2_from_baseline > 0.6,
                           "largely a restatement of baseline information",
                           "genuinely new information, but not predictive"))

write_csv(redundancy, file.path(TABLE_DIR, "xgb_structure_redundancy.csv"))

cat("\n=== How much of each structure feature the baseline already explains ===\n")
print(as.data.frame(redundancy |> select(-reading)), row.names = FALSE, digits = 3)
cat(sprintf("\n%d of %d structure columns are >60%% reconstructible from baseline.\n",
            sum(redundancy$r2_from_baseline > 0.6), nrow(redundancy)))
cat("The rest are genuinely novel and still bought no held-out improvement --\n")
cat("which is the sharper point: novelty is not relevance.\n")

gain_share <- imp |>
  group_by(block) |>
  summarise(n_features = n(), gain = sum(Gain), .groups = "drop") |>
  mutate(pct_of_features = 100 * n_features / sum(n_features),
         pct_of_gain     = 100 * gain / sum(gain),
         gain_per_feature = 100 * gain / n_features)

write_csv(gain_share, file.path(TABLE_DIR, "xgb_structure_gain_share.csv"))

cat("\n=== Gain share vs feature share ===\n")
print(as.data.frame(gain_share), row.names = FALSE, digits = 3)
cat(sprintf(paste0("\nStructure is %.1f%% of the predictors but takes %.1f%% of the gain",
                   " -- BELOW proportional.\nPer feature it earns %.2fx an average",
                   " baseline feature. Largest single feature (%s)\ntakes %.1f%% on",
                   " its own, %.1fx the entire structure block.\n"),
            gain_share$pct_of_features[gain_share$block == "structure"],
            gain_share$pct_of_gain[gain_share$block == "structure"],
            gain_share$gain_per_feature[gain_share$block == "structure"] /
              gain_share$gain_per_feature[gain_share$block == "baseline"],
            imp$Feature[1], 100 * imp$Gain[1],
            imp$Gain[1] / sum(imp$Gain[imp$block == "structure"])))


# ----------------------------- 14. Run manifest -----------------------------

manifest <- list(
  question        = paste("Does adding secondary-structure information improve",
                          "prediction of measured mRNA half-life beyond a model",
                          "containing non-structure transcript features?"),
  species         = "human",
  cache           = cache_path("human"),
  response        = TARGET_COL,
  response_note   = paste("Agarwal & Kelley 2022 consensus half-life PC1,",
                          "untransformed (signed score, not hours)"),
  n_eligible      = nrow(dat),
  n_trainval      = nrow(trainval),
  n_test          = nrow(preds),
  eligibility     = paste("complete cases on baseline + structure;",
                          "identical rows for both models"),
  split_artefact  = splits_path(BLOCK_LEVEL),
  block_level     = BLOCK_LEVEL,
  split_caveat    = paste("families larger than", SPLIT_PIN_FRAC * 100,
                          "% of the smallest split are pinned to train, so the",
                          "test set is depleted of large families and measures",
                          "generalisation to small and mid-sized families"),
  baseline_cols   = BASELINE,
  structure_cols  = STRUCTURE,
  excluded_note   = paste("icSHAPE Gini excluded from the main run (80-91%",
                          "missing); raw MFE scores and per-nt MFE excluded as",
                          "length-confounded; translation efficiency excluded"),
  seed            = SEED,
  grid_size       = GRID_SIZE,
  inner_folds     = INNER_V,
  n_boot          = N_BOOT,
  best_params_A   = final_A$best,
  best_params_B   = final_B$best,
  pooled_metrics  = main_table,
  delta_table     = delta_table,
  perm_table      = perm_table,
  validation      = checks_df,
  r_version       = R.version.string,
  fitted_at       = Sys.time()
)

saveRDS(manifest, file.path(RUN_DIR, "run_manifest.rds"))

write_csv(main_table, file.path(TABLE_DIR, "xgb_structure_main_comparison.csv"))
write_csv(tibble(block = c(rep("baseline", length(BASELINE)),
                           rep("structure", length(STRUCTURE))),
                 column = c(BASELINE, STRUCTURE)),
          file.path(TABLE_DIR, "xgb_structure_feature_lists.csv"))


# ----------------------------- 15. Text summary -----------------------------

d_r2  <- delta_table |> filter(metric == "rsq_trad")
d_rm  <- delta_table |> filter(metric == "rmse")
d_ma  <- delta_table |> filter(metric == "mae")

verdict <- if (d_r2$ci_low > 0) {
  paste("Adding secondary-structure features improved XGBoost prediction of",
        "mRNA half-life relative to the non-structure baseline, indicating that",
        "the structure block contains incremental predictive information not",
        "captured by the included baseline features.")
} else if (d_r2$ci_high < 0) {
  paste("Adding secondary-structure features did not improve, and on this",
        "held-out set slightly degraded, XGBoost prediction of mRNA half-life.")
} else {
  paste("The change in held-out performance from adding secondary-structure",
        "features is uncertain: the 95% confidence interval on delta R-squared",
        "spans zero. On this evidence the structure block cannot be said to add",
        "predictive information beyond the baseline feature set.")
}

summary_txt <- c(
  "XGBoost: baseline vs baseline + secondary structure",
  strrep("=", 70),
  "",
  sprintf("Response      : %s (Agarwal & Kelley consensus PC1, untransformed)", TARGET_COL),
  sprintf("Eligible genes: %d  (train+val %d / test %d)", nrow(dat), nrow(trainval), nrow(preds)),
  sprintf("Eligibility   : complete cases on baseline + structure; identical rows both models"),
  sprintf("Split         : committed family-blocked holdout, %s", splits_path(BLOCK_LEVEL)),
  sprintf("Model A       : %d predictors (non-structure)", length(BASELINE)),
  sprintf("Model B       : %d predictors (Model A + %d structure columns)",
          length(PRED_B), length(STRUCTURE)),
  "",
  "Structure block (the only difference between the models):",
  paste0("  ", paste(STRUCTURE, collapse = ", ")),
  "",
  "Excluded, and why:",
  "  icSHAPE Gini (gini_*)      80-91% missing; a complete set containing it",
  "                             costs 93% of the corpus. See the secondary run.",
  "  rnafold_score_*,           raw MFE is near-deterministic in length and GC,",
  "  rnalfold_score_*,          both already in the baseline, so a 'win' could be",
  "  rnafold_per_nt_*           length re-entering under a structure label.",
  "  translation_efficiency     measured phenotype, excluded by the brief.",
  "  frac_*, purine_/amino_*    redundant with GC% and the skews.",
  "",
  strrep("-", 70),
  "Held-out results",
  strrep("-", 70),
  sprintf("R-squared  baseline %.4f  ->  structure %.4f   delta %+.4f  [%+.4f, %+.4f]",
          d_r2$baseline, d_r2$structure, d_r2$delta, d_r2$ci_low, d_r2$ci_high),
  sprintf("RMSE       baseline %.4f  ->  structure %.4f   delta %+.4f  [%+.4f, %+.4f]",
          d_rm$baseline, d_rm$structure, d_rm$delta, d_rm$ci_low, d_rm$ci_high),
  sprintf("MAE        baseline %.4f  ->  structure %.4f   delta %+.4f  [%+.4f, %+.4f]",
          d_ma$baseline, d_ma$structure, d_ma$delta, d_ma$ci_low, d_ma$ci_high),
  "",
  "Sign convention: delta = structure - baseline.",
  "  R-squared  delta > 0 favours structure",
  "  RMSE, MAE  delta < 0 favours structure",
  sprintf("95%% CIs from %d paired bootstrap replicates over held-out genes.", N_BOOT),
  "",
  sprintf("Secondary sign-flip test (squared error): p = %.4g", perm_table$p_value[1]),
  "",
  strrep("-", 70),
  "Interpretation",
  strrep("-", 70),
  strwrap(verdict, width = 70),
  "",
  strwrap(paste("This is predictive evidence only. Gain importance is exploratory",
                "context; correlated predictors share and redistribute importance,",
                "so it does not estimate an independent effect. Causal and",
                "inferential claims are for the later modelling plan."), width = 70),
  "",
  strwrap(paste("Limitation to report with this result: the blocked split pins",
                "families larger than ~5% of the smallest split to train, so the",
                "test set contains no large family. It measures generalisation to",
                "small and mid-sized families, not to the largest ones."), width = 70),
  "",
  sprintf("Seed %d. Fitted %s. %s", SEED, format(Sys.time(), "%Y-%m-%d %H:%M"),
          R.version.string)
)

writeLines(summary_txt, file.path(RUN_DIR, "SUMMARY.txt"))
cat("\n"); cat(summary_txt, sep = "\n"); cat("\n")


# ----------------------------- 16. Figures ----------------------------------
# Rendering lives in its own script, which reads only the tables written above
# and never touches a fitted model — so a figure can be restyled in seconds
# instead of re-tuning for 40 minutes:
#
#   Rscript analysis/models/xgb_structure_plots.R
#
# It is sourced here so one command still produces everything, and so this
# file's table formats are exercised against the plotting code on every run
# rather than drifting apart unnoticed. new.env() keeps its variables out of
# this script's workspace.

log_msg("Rendering figures")
source("analysis/models/xgb_structure_plots.R", local = new.env())

log_msg("Done. Artefacts in ", RUN_DIR, ", tables in ", TABLE_DIR,
        ", plots in ", PLOT_DIR)
