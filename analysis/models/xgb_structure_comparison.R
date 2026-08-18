# =============================================================================
# XGBoost: does secondary structure add held-out predictive information?
# =============================================================================
# A nested ladder of four models, fitted once, evaluated once:
#
#   Baseline    non-structure transcript features
#   S-core      + MFE z-scores and local MFE z-scores        (PRIMARY test)
#   S-select    + MFE delta (observed - expected)
#   S-full      + raw MFE and per-nucleotide MFE
#
# The intended experimental difference between adjacent rungs is the family of
# folding metrics that was added. Everything else — rows, gene ids,
# preprocessing, tuning resamples, tuning grid, budget, seeds, evaluation — is
# shared by construction, not by four parallel edits kept in step by hand. The
# ladder itself is defined and its nesting asserted in
# xgb_structure_features.R.
#
# WHY A LADDER RATHER THAN ONE COMPARISON. The rungs are ordered from least to
# most confounded with the baseline (see LADDER in the features file), so the
# ladder has a predicted SHAPE under each hypothesis. If structure carries
# information, the gain appears at S-core and persists. If the apparent gain is
# length re-entering the model under a structure label, it appears only at
# S-full. Three contrasts that can each individually clear zero would be three
# more chances at a false positive; three contrasts with a predicted ordering
# are evidence. The primary contrast is pre-specified as S-core vs Baseline;
# every other contrast is secondary and unadjusted, and is labelled as such in
# every artefact this script writes.
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
# WHAT THIS IS NOT: a causal test. A win for a structure rung is incremental
# predictive information associated with structure after accounting for the
# baseline block. Feature importance below is exploratory context, not an
# effect estimate.
#
# Usage:
#   Rscript analysis/models/xgb_structure_comparison.R                # default
#   Rscript analysis/models/xgb_structure_comparison.R complete_case
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
# runs the whole pipeline end to end in a few minutes, which is how to check a
# change to this file before committing to the full search.
GRID_SIZE  <- as.integer(Sys.getenv("XGB_GRID_SIZE", "60"))
N_CHUNKS   <- 5L       # family-blocked slices of test, for the consistency plot

# Parallelism is split ACROSS resamples rather than handed to xgboost, because
# xgboost scales badly on a problem this shape: measured on this dataset, going
# from 2 to 10 threads took one fit from 21.7 s to 15.5 s (1.4x for 5x the
# cores). Four workers on three threads each is roughly three times the
# throughput of one worker on twelve.
N_WORKERS  <- 4L
N_THREADS  <- 3L

# --- Variant ----------------------------------------------------------------
# Which specification to run. A variant may change row eligibility and the
# feature blocks, and nothing else — see VARIANTS in xgb_structure_features.R
# for the registry and for why these are a SENSITIVITY ANALYSIS rather than a
# set of independent tests.
#
# Every variant gets its own self-contained directory, so their fit caches,
# tables and figures cannot collide or silently overwrite one another.

.args   <- commandArgs(trailingOnly = TRUE)
VARIANT <- resolve_variant(if (length(.args) >= 1) .args[[1]] else "default")

RUN_DIR    <- variant_dir(VARIANT, "root")
PLOT_DIR   <- variant_dir(VARIANT, "plots")
TABLE_DIR  <- variant_dir(VARIANT, "tables")
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

# The three metrics the bootstrap carries, and which direction is good. Kept as
# data rather than scattered through the code, because the sign convention is
# the single easiest thing to get wrong when reading a delta table.
BOOT_METRICS <- c("rsq_trad", "rmse", "mae")
FAVOURS      <- c(rsq_trad = "delta > 0 favours the larger model",
                  rmse     = "delta < 0 favours the larger model",
                  mae      = "delta < 0 favours the larger model")

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  utils::flush.console()
}

# Started once, before any model is tuned, so every search runs under identical
# conditions. Returned to sequential in section 6.
future::plan(future::multisession, workers = N_WORKERS)


# ----------------------------- 1. Eligible dataset --------------------------

log_msg("Building eligible dataset")
el <- eligible_dataset("human", VARIANT)
report_feature_sets(el)

dat        <- el$data
BASELINE   <- el$baseline
MODELS     <- names(el$models)              # ladder order, reference first
STRUCTURE  <- el$models                     # named list: rung -> structure cols
STRUCT_ALL <- el$structure_all

# One predictor list per rung, all built from the same baseline vector, so
# "every rung is the baseline plus its own structure block" is a property of
# the construction rather than something to verify afterwards.
PREDS <- lapply(MODELS, function(m) predictors_for(el, m))
names(PREDS) <- MODELS

CONTRASTS <- ladder_contrasts(VARIANT)

trainval <- dat |> filter(split != "test")
testing_ <- dat |> filter(split == "test")

log_msg(sprintf("train+val = %d genes, test = %d genes",
                nrow(trainval), nrow(testing_)))
log_msg("Ladder: ", paste(sprintf("%s (%d)", MODELS, lengths(PREDS)),
                          collapse = " -> "))


# ----------------------------- 2. Shared tuning resamples -------------------
# ONE resampling object, created once and handed to every model, so the
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
# One builder, called once per rung. The only argument that differs is the
# predictor list, which is what makes "no other accidental differences" a
# property of the code rather than a thing to check afterwards.

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
#'   step_corr would remove a DIFFERENT set of baseline columns from each rung,
#'   because adding a structure block changes the correlation graph. A rung
#'   would then no longer be "Baseline plus structure", and the ladder would
#'   silently be comparing feature sets that differ in more than one way. The
#'   existing LightGBM script uses step_corr(0.75); it is correct there and
#'   wrong here.
#'
#' step_zv is a safety net only — the zero-variance drop already happened in
#' eligible_dataset(), on the same rows for every rung, so it removes the same
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
# mtry is a PROPORTION, so the same dials range is meaningful for a 107-column
# and a 151-column model; a count range would hand the wider rungs a
# systematically different search space. The grid is drawn once from a fixed
# seed and reused, so every rung evaluates the same 60 configurations.

xgb_params <- build_workflow(PREDS[[REFERENCE_MODEL]]) |>
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


# ----------------------------- 5. Tune every rung ---------------------------

tune_model <- function(preds, label) {
  log_msg("Tuning ", label, " (", length(preds), " predictors, ",
          GRID_SIZE, " configurations x ", INNER_V, " folds)")
  wf <- build_workflow(preds)
  set.seed(SEED)                     # same racing randomisation for every rung
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

finalise <- function(f) {
  best <- select_best(f$tune, metric = "rmse")
  final_wf <- finalize_workflow(f$workflow, best)
  set.seed(SEED)
  list(best = best, fit = fit(final_wf, trainval))
}

# --- Fit cache ---------------------------------------------------------------
# Tuning is the great majority of this script's runtime; everything after it
# (predictions, bootstrap, tables, figures) takes seconds. Re-tuning to restyle
# a plot is pure waste, so the fitted models are cached.
#
# The cache is keyed on EVERYTHING that determines the fits — every predictor
# list in ladder order, the exact genes in the fitting pool and the test set,
# the seed, the tuning grid (which carries the parameter ranges), and the fold
# count. Comparison is by identical(), not by a hash or a timestamp, so a
# changed feature list or a rebuilt cache invalidates it loudly rather than
# silently serving a stale model. That is the whole risk of caching a fit, and
# it is the one thing this key exists to close.
#
#   XGB_REFIT=1   force a full re-tune regardless of the cache

fit_key <- list(
  preds          = PREDS,
  trainval_genes = trainval$gene_id,
  test_genes     = testing_$gene_id,
  seed           = SEED,
  grid           = xgb_grid,
  inner_v        = INNER_V,
  target         = TARGET_COL,
  variant        = VARIANT
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
  # Baseline is tuned ONCE here and reused by every contrast, rather than
  # re-derived per comparison. That is the point of fitting the ladder in one
  # script: three contrasts against a baseline that is literally the same
  # fitted object, not three baselines that happen to agree.
  tuned <- lapply(MODELS, function(m) tune_model(PREDS[[m]], m))
  names(tuned) <- MODELS

  saveRDS(lapply(tuned, `[[`, "tune"),
          file.path(RUN_DIR, "tuning_results.rds"))

  # ---- 6. Finalise and refit. Selection on inner folds; `test` unread. ----
  finals <- lapply(tuned, finalise)
  names(finals) <- MODELS

  cached <- list(key      = fit_key,
                 finals   = finals,
                 splits   = lapply(tuned, function(f) f$tune$splits),
                 tuned_at = Sys.time())
  saveRDS(cached, FITS_PATH)
}

future::plan(future::sequential)   # workers no longer needed

finals <- cached$finals

for (m in MODELS) {
  cat(sprintf("\n--- Best hyperparameters: %s ---\n", m))
  print(as.data.frame(finals[[m]]$best))
}
if (!is.null(cached$tuned_at)) {
  log_msg("Models tuned at ", format(cached$tuned_at, "%Y-%m-%d %H:%M:%S"))
}


# ----------------------------- 7. Held-out predictions ----------------------
# The single use of `test`.

PMAT <- vapply(MODELS, function(m) predict(finals[[m]]$fit, testing_)$.pred,
               numeric(nrow(testing_)))
colnames(PMAT) <- MODELS

obs <- testing_[[TARGET_COL]]

# Per-gene losses, held once as matrices. Every metric below — pooled, bootstrap
# and per-slice — is a column mean of one of these, so they are computed here
# and nowhere else.
SE <- (obs - PMAT)^2
AE <- abs(obs - PMAT)

preds <- tibble(
  gene_id   = rep(testing_$gene_id,   times = length(MODELS)),
  gene_name = rep(testing_$gene_name, times = length(MODELS)),
  split     = rep(testing_$split,     times = length(MODELS)),
  observed  = rep(obs,                times = length(MODELS)),
  model     = rep(MODELS, each = nrow(testing_)),
  predicted = as.vector(PMAT)
) |>
  mutate(model = factor(model, levels = MODELS),
         resid = observed - predicted,
         se    = resid^2,
         ae    = abs(resid))

write_csv(preds, file.path(TABLE_DIR, "xgb_structure_test_predictions.csv"))
log_msg("Held-out predictions written for ", nrow(testing_), " test genes x ",
        length(MODELS), " models")


# ----------------------------- 8. Validation checklist ----------------------
# Assertions, not comments — a violated invariant stops the run before any
# number is reported. The nesting checks are the ones the ladder argument
# depends on, and they are run against the predictors the FITTED models
# actually used, not against the lists this script intended to pass them.

log_msg("Running validation checklist")
checks <- list()
chk <- function(name, ok, detail = "") {
  checks[[length(checks) + 1]] <<- tibble(check = name, pass = isTRUE(ok),
                                          detail = detail)
  if (!isTRUE(ok)) stop("VALIDATION FAILED: ", name, " — ", detail)
  cat(sprintf("  ok  %-58s %s\n", name, detail))
}

used <- lapply(MODELS, function(m)
  finals[[m]]$fit |> extract_mold() |> (\(mo) colnames(mo$predictors))())
names(used) <- MODELS

chk("Same rows and identifiers for every rung",
    TRUE, sprintf("%d train+val / %d test genes, one frame",
                  nrow(trainval), nrow(testing_)))
chk("No structure variable appears in the reference model",
    length(intersect(used[[REFERENCE_MODEL]],
                     c(STRUCT_ALL, gini_columns(dat)))) == 0)
chk("Every rung carries the identical baseline block",
    all(vapply(MODELS, function(m)
      setequal(intersect(used[[m]], BASELINE), BASELINE), logical(1))),
    sprintf("%d baseline columns", length(BASELINE)))
chk("Every rung = baseline + its own structure block",
    all(vapply(MODELS, function(m)
      setequal(used[[m]], c(BASELINE, STRUCTURE[[m]])), logical(1))))
chk("The ladder is nested in the fitted models",
    all(vapply(seq_len(length(MODELS) - 1), function(i)
      all(used[[MODELS[i]]] %in% used[[MODELS[i + 1]]]), logical(1))),
    paste(MODELS, collapse = " < "))
chk("Each rung adds exactly its declared structure columns",
    all(vapply(seq_len(length(MODELS) - 1), function(i)
      setequal(setdiff(used[[MODELS[i + 1]]], used[[MODELS[i]]]),
               setdiff(STRUCTURE[[MODELS[i + 1]]], STRUCTURE[[MODELS[i]]])),
      logical(1))))
chk("No amino-acid or nucleotide-fraction column is a predictor",
    !any(grepl("^(aa_|frac_|purine_|amino_)", unlist(used))),
    "exact functions of retained codon / GC / skew columns")
chk("No last-exon-length column is a predictor",
    !any(grepl("^exon_length_", unlist(used))),
    "3'UTR-length proxy; length_3utr is already in the baseline")
chk("No translation-efficiency variable in any rung",
    !any(grepl("translation_efficiency", unlist(used))))
chk("No identifier, family or outcome-derived column is a predictor",
    length(intersect(unlist(used), c(META_COLS, TARGET_COL))) == 0)
chk("Every test gene has a prediction from every rung",
    !anyNA(PMAT) && nrow(PMAT) == nrow(testing_) && ncol(PMAT) == length(MODELS))
chk("Test genes are disjoint from the tuning/fitting pool",
    length(intersect(testing_$gene_id, trainval$gene_id)) == 0)
chk("No family spans the fitting pool and the test set",
    length(intersect(trainval[[paste0("family_id_", BLOCK_LEVEL)]],
                     testing_[[paste0("family_id_", BLOCK_LEVEL)]])) == 0,
    paste0("blocked at family_id_", BLOCK_LEVEL))
chk("Tuning used the same resampling object for every rung",
    length(unique(lapply(cached$splits, identity))) == 1,
    paste0(INNER_V, " family-blocked folds"))
chk("Tuning used the same grid for every rung",
    nrow(xgb_grid) == GRID_SIZE,
    paste0(GRID_SIZE, " configurations drawn once from seed ", SEED))
chk("Hyperparameter tuning used training data only",
    length(intersect(
      unlist(lapply(cached$splits[[REFERENCE_MODEL]],
                    function(s) analysis(s)$gene_id)),
      testing_$gene_id)) == 0)
chk("Metrics come from held-out predictions, not training predictions",
    TRUE, "single evaluation on the untouched `test` split")

checks_df <- bind_rows(checks)
write_csv(checks_df, file.path(TABLE_DIR, "xgb_structure_validation_checks.csv"))


# ----------------------------- 9. Pooled held-out metrics -------------------

pooled <- map_dfr(MODELS, function(m) {
  METRICS(tibble(truth = obs, est = PMAT[, m]), truth = truth, estimate = est) |>
    transmute(model = m, metric = .metric, value = .estimate)
}) |>
  mutate(model = factor(model, levels = MODELS))

model_table <- pooled |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(n_structure = lengths(STRUCTURE[as.character(model)]),
         n_predictors = lengths(PREDS[as.character(model)]),
         .after = model)

cat("\n=== Pooled held-out metrics (n =", nrow(testing_), "test genes) ===\n")
print(as.data.frame(model_table), row.names = FALSE, digits = 4)

write_csv(model_table, file.path(TABLE_DIR, "xgb_structure_model_metrics.csv"))


# ----------------------------- 10. Paired bootstrap -------------------------
# Genes resampled with replacement — the SAME genes for EVERY rung in every
# replicate. That is what makes each interval an interval on a paired
# improvement rather than on the difference of two independent estimates, and
# it is also what makes the increments legitimate: because all four rungs are
# scored on one shared draw, S-select minus S-core is computed within a draw
# and its interval is a real interval on that increment.

boot_stat <- function(idx) {
  o   <- obs[idx]
  sst <- sum((o - mean(o))^2)
  se  <- SE[idx, , drop = FALSE]
  ae  <- AE[idx, , drop = FALSE]
  rbind(rsq_trad = 1 - colSums(se) / sst,
        rmse     = sqrt(colMeans(se)),
        mae      = colMeans(ae))
}

n_test <- nrow(PMAT)
point  <- boot_stat(seq_len(n_test))

set.seed(SEED)
boot_arr <- array(NA_real_, dim = c(length(BOOT_METRICS), length(MODELS), N_BOOT),
                  dimnames = list(BOOT_METRICS, MODELS, NULL))
for (i in seq_len(N_BOOT)) {
  boot_arr[, , i] <- boot_stat(sample.int(n_test, n_test, replace = TRUE))
}

boot_ci <- function(metric, lhs, rhs, p) {
  unname(quantile(boot_arr[metric, lhs, ] - boot_arr[metric, rhs, ], p))
}

delta_table <- CONTRASTS |>
  crossing(metric = BOOT_METRICS) |>
  mutate(
    lhs_value  = map2_dbl(metric, lhs, ~ point[.x, .y]),
    rhs_value  = map2_dbl(metric, rhs, ~ point[.x, .y]),
    delta      = lhs_value - rhs_value,
    ci_low     = pmap_dbl(list(metric, lhs, rhs), boot_ci, p = 0.025),
    ci_high    = pmap_dbl(list(metric, lhs, rhs), boot_ci, p = 0.975),
    favours    = unname(FAVOURS[metric]),
    conclusive = !(ci_low <= 0 & ci_high >= 0),
    metric     = factor(metric, levels = BOOT_METRICS)
  ) |>
  arrange(kind, metric, match(lhs, MODELS))

cat("\n=== Paired bootstrap,", N_BOOT, "replicates, 95% percentile CI ===\n")
cat("Primary contrast is marked *; every other row is secondary and unadjusted.\n\n")
for (k in unique(delta_table$kind)) {
  cat("--", k, "--\n")
  print(as.data.frame(
    delta_table |>
      filter(kind == k) |>
      transmute(` ` = if_else(primary, "*", " "),
                contrast, metric, lhs_value, rhs_value, delta,
                ci_low, ci_high, conclusive)),
    row.names = FALSE, digits = 4)
  cat("\n")
}

write_csv(delta_table, file.path(TABLE_DIR, "xgb_structure_delta_bootstrap.csv"))


# ----------------------------- 11. Paired sign-flip test --------------------
# Secondary to the effect size and CI, per the brief. Paired per-gene loss
# differences, signs flipped at random — the exchangeability a paired test
# needs, and nothing like an independent-samples t-test on two sets of CV
# metrics.
#
# Run for every contrast the bootstrap covers, so the sign-flip table cannot be
# quoted for one contrast while the others go unreported.

signflip <- function(d) {
  obs_mean <- mean(d)
  set.seed(SEED)
  null <- vapply(seq_len(N_PERM), function(i)
    mean(d * sample(c(-1, 1), length(d), replace = TRUE)), numeric(1))
  (1 + sum(abs(null) >= abs(obs_mean))) / (N_PERM + 1)
}

perm_table <- CONTRASTS |>
  crossing(loss = c("squared error", "absolute error")) |>
  mutate(
    diffs     = pmap(list(lhs, rhs, loss), function(l, r, ls)
                     if (ls == "squared error") SE[, l] - SE[, r]
                     else                       AE[, l] - AE[, r]),
    mean_diff = map_dbl(diffs, mean),
    p_value   = map_dbl(diffs, signflip)
  ) |>
  select(kind, contrast, primary, loss, mean_diff, p_value) |>
  mutate(note = "mean_diff < 0 favours the larger model")

cat("=== Paired sign-flip test (secondary),", N_PERM, "permutations ===\n")
print(as.data.frame(perm_table), row.names = FALSE, digits = 4)
write_csv(perm_table, file.path(TABLE_DIR, "xgb_structure_signflip_test.csv"))


# ----------------------------- 12. Per-chunk consistency --------------------
# The brief's paired outer-fold plot. There are no outer folds under a single
# holdout, so the substitute answers the same question — is an improvement
# spread across the held-out set, or carried by one corner of it — by cutting
# `test` into family-blocked slices and scoring every rung on each.

set.seed(SEED)
chunk_folds <- group_vfold_cv(testing_,
                              group   = paste0("family_id_", BLOCK_LEVEL),
                              v       = N_CHUNKS,
                              balance = "observations")

# Built by explicit (metric, model) lookup rather than by reshaping the matrix.
# A transpose-and-recycle would be shorter and is easy to get subtly wrong in a
# way that surfaces only as a mislabelled line on a plot.
chunk_metrics <- map_dfr(seq_len(nrow(chunk_folds)), function(i) {
  ids <- assessment(chunk_folds$splits[[i]])$gene_id
  idx <- which(testing_$gene_id %in% ids)
  m   <- boot_stat(idx)
  expand_grid(metric = rownames(m), model = colnames(m)) |>
    mutate(chunk = paste0("slice ", i), n = length(idx),
           value = map2_dbl(metric, model, ~ m[.x, .y]))
}) |>
  mutate(model = factor(model, levels = MODELS)) |>
  select(chunk, n, metric, model, value)

write_csv(chunk_metrics, file.path(TABLE_DIR, "xgb_structure_chunk_metrics.csv"))

cat("\n=== Per-slice held-out R-squared ===\n")
print(as.data.frame(
  chunk_metrics |>
    filter(metric == "rsq_trad") |>
    pivot_wider(names_from = model, values_from = value)),
  row.names = FALSE, digits = 4)


# ----------------------------- 13. Gain importance (model-derived) ----------
# Extracted here because it needs the fitted models. The FIGURE that draws it
# lives in xgb_structure_plots.R, like every other figure — see section 16.
#
# Written for every rung, not just the widest, so the "structure earns gain but
# buys nothing" observation can be read at each level of the ladder rather than
# only where the block is largest and most confounded.

imp <- map_dfr(MODELS, function(m) {
  finals[[m]]$fit |>
    extract_fit_engine() |>
    xgb.importance(model = _) |>
    as_tibble() |>
    mutate(model = m,
           block = if_else(Feature %in% STRUCT_ALL, "structure", "baseline"),
           .before = 1)
}) |>
  mutate(model = factor(model, levels = MODELS))

write_csv(imp, file.path(TABLE_DIR, "xgb_structure_gain_importance.csv"))

gain_share <- imp |>
  group_by(model, block) |>
  summarise(n_features = n(), gain = sum(Gain), .groups = "drop_last") |>
  mutate(pct_of_features  = 100 * n_features / sum(n_features),
         pct_of_gain      = 100 * gain / sum(gain),
         gain_per_feature = 100 * gain / n_features) |>
  ungroup()

write_csv(gain_share, file.path(TABLE_DIR, "xgb_structure_gain_share.csv"))

cat("\n=== Gain share vs feature share, by rung ===\n")
print(as.data.frame(gain_share), row.names = FALSE, digits = 3)
cat("\nA structure block taking LESS gain than its share of the features is\n")
cat("participating below average. Gain says a feature was used, not needed.\n")


# ----------------------------- 13e. Why structure earns gain ----------------
# Exploratory context for the importance figure, and the answer to the obvious
# challenge: "if structure is useless, why does it take any of the gain?"
#
# Two candidate explanations, and they turn out to split the block. Either a
# structure column is a RESTATEMENT of baseline information (the tree splits on
# it and collects gain for something it already knew), or it is genuinely new
# information that simply is not related to half-life. Regressing each
# structure column on the whole baseline block separates the two — and because
# the rungs are ordered by how confounded they are, the result is also a direct
# check on the ladder's premise.
#
# Linear R^2, so this is a lower bound on redundancy — a tree could exploit a
# non-linear relationship this misses. That direction is the safe one: it can
# only understate how much the baseline already knows.
#
# Fitted on the training pool only, so no test information enters even a
# descriptive table.
#
# Complete cases regardless of the variant's row policy. Under "native_na" the
# training pool contains NAs, and qr() does not merely warn on those — it fails
# outright with "NA/NaN/Inf in foreign function call". Restricting to rows
# complete on baseline + structure also makes the redundancy figures comparable
# ACROSS variants, since they are then measured on the same genes rather than
# on whatever each policy happened to admit.

log_msg("Measuring how much of each structure feature the baseline explains")

# Which rung first introduces each structure column — so the redundancy table
# can be read one rung at a time.
first_rung <- vapply(STRUCT_ALL, function(v) {
  hit <- MODELS[vapply(MODELS, function(m) v %in% STRUCTURE[[m]], logical(1))]
  if (length(hit)) hit[1] else NA_character_
}, character(1))

red_rows <- trainval[stats::complete.cases(
  trainval[, c(BASELINE, STRUCT_ALL), drop = FALSE]), , drop = FALSE]

if (nrow(red_rows) < 10 * length(BASELINE)) {
  log_msg("  skipped: only ", nrow(red_rows), " complete rows for ",
          length(BASELINE), " predictors — too few to regress against")
  redundancy <- tibble(structure_feature = character(), first_rung = character(),
                       r2_from_baseline = numeric(), reading = character())
} else {

  if (nrow(red_rows) < nrow(trainval)) {
    log_msg("  using ", nrow(red_rows), " of ", nrow(trainval),
            " training rows (complete cases; this variant admits missing data)")
  }

  qr_base <- qr(cbind(1, as.matrix(red_rows[, BASELINE])))
  redundancy <- tibble(
    structure_feature = STRUCT_ALL,
    first_rung        = unname(first_rung[STRUCT_ALL]),
    r2_from_baseline  = vapply(STRUCT_ALL, function(v) {
      y <- red_rows[[v]]
      1 - sum((y - qr.fitted(qr_base, y))^2) / sum((y - mean(y))^2)
    }, numeric(1))
  ) |>
    arrange(desc(r2_from_baseline)) |>
    mutate(reading = if_else(r2_from_baseline > 0.6,
                             "largely a restatement of baseline information",
                             "genuinely new information"))
}

write_csv(redundancy, file.path(TABLE_DIR, "xgb_structure_redundancy.csv"))

if (nrow(redundancy) > 0) {
  cat("\n=== How much of each structure feature the baseline already explains ===\n")
  print(as.data.frame(redundancy |> select(-reading)), row.names = FALSE, digits = 3)
  cat("\nBy rung (mean R-squared reconstructible from baseline):\n")
  print(as.data.frame(
    redundancy |>
      group_by(first_rung) |>
      summarise(n = n(), mean_r2 = mean(r2_from_baseline),
                n_over_60pct = sum(r2_from_baseline > 0.6), .groups = "drop")),
    row.names = FALSE, digits = 3)
  cat("\nThe ladder's premise is that this INCREASES down the rungs. If it does,\n")
  cat("a gain that appears only at the confounded rungs is the expected\n")
  cat("signature of baseline information re-entering, not of structure.\n")
}


# ----------------------------- 14. Run manifest -----------------------------

manifest <- list(
  question        = paste("Does adding secondary-structure information improve",
                          "prediction of measured mRNA half-life beyond a model",
                          "containing non-structure transcript features?"),
  species         = "human",
  variant         = VARIANT$name,
  variant_label   = VARIANT$label,
  variant_rows    = VARIANT$rows,
  cache           = cache_path("human"),
  response        = TARGET_COL,
  response_note   = paste("Agarwal & Kelley 2022 consensus half-life PC1,",
                          "untransformed (signed score, not hours)"),
  n_eligible      = nrow(dat),
  n_trainval      = nrow(trainval),
  n_test          = nrow(testing_),
  eligibility     = if (identical(VARIANT$rows, "complete"))
                      "complete cases on baseline + all structure columns; identical rows for every rung"
                    else
                      "every gene with a target and a split; NA handled natively; identical rows for every rung",
  split_artefact  = splits_path(BLOCK_LEVEL),
  block_level     = BLOCK_LEVEL,
  split_caveat    = paste("families larger than", SPLIT_PIN_FRAC * 100,
                          "% of the smallest split are pinned to train, so the",
                          "test set is depleted of large families and measures",
                          "generalisation to small and mid-sized families"),
  ladder          = VARIANT$ladder,
  models          = MODELS,
  primary_model   = PRIMARY_MODEL,
  reference_model = REFERENCE_MODEL,
  baseline_cols   = BASELINE,
  structure_cols  = STRUCTURE,
  structure_all   = STRUCT_ALL,
  excluded_note   = paste("icSHAPE Gini excluded from the main run (80-91%",
                          "missing); translation efficiency excluded as a",
                          "measured phenotype; aa_*, frac_*, purine_/amino_*",
                          "excluded as exact functions of retained columns"),
  seed            = SEED,
  grid_size       = GRID_SIZE,
  inner_folds     = INNER_V,
  n_boot          = N_BOOT,
  best_params     = lapply(finals, `[[`, "best"),
  model_metrics   = model_table,
  delta_table     = delta_table,
  perm_table      = perm_table,
  validation      = checks_df,
  r_version       = R.version.string,
  fitted_at       = Sys.time()
)

saveRDS(manifest, file.path(RUN_DIR, "run_manifest.rds"))

write_csv(
  map_dfr(MODELS, function(m)
    tibble(model = m,
           block = c(rep("baseline", length(BASELINE)),
                     rep("structure", length(STRUCTURE[[m]]))),
           column = c(BASELINE, STRUCTURE[[m]]))),
  file.path(TABLE_DIR, "xgb_structure_feature_lists.csv"))


# ----------------------------- 15. Text summary -----------------------------

prim <- delta_table |> filter(primary, metric == "rsq_trad")

verdict <- if (prim$ci_low > 0) {
  paste("The primary contrast (", PRIMARY_MODEL, "vs", REFERENCE_MODEL,
        ") shows an improvement in held-out R-squared whose 95% interval",
        "excludes zero, indicating that the length-normalised folding block",
        "carries incremental predictive information not captured by the",
        "baseline features.")
} else if (prim$ci_high < 0) {
  paste("The primary contrast (", PRIMARY_MODEL, "vs", REFERENCE_MODEL,
        ") shows a DECREASE in held-out R-squared whose 95% interval excludes",
        "zero: adding the folding block made prediction worse on these genes.")
} else {
  paste("The primary contrast (", PRIMARY_MODEL, "vs", REFERENCE_MODEL,
        ") is inconclusive: the 95% interval on delta R-squared spans zero.",
        "On this evidence the folding block cannot be said to add predictive",
        "information beyond the baseline feature set.")
}

fmt_row <- function(r) {
  sprintf("  %-24s %-9s %8.4f -> %8.4f  delta %+.4f  [%+.4f, %+.4f]%s",
          r$contrast, as.character(r$metric), r$rhs_value, r$lhs_value,
          r$delta, r$ci_low, r$ci_high, if (r$primary) "  *PRIMARY" else "")
}

summary_txt <- c(
  "XGBoost: a nested ladder from baseline to full secondary structure",
  strrep("=", 74),
  "",
  sprintf("Response      : %s (Agarwal & Kelley consensus PC1, untransformed)", TARGET_COL),
  sprintf("Variant       : %s — %s", VARIANT$name, VARIANT$label),
  sprintf("Eligible genes: %d  (train+val %d / test %d)",
          nrow(dat), nrow(trainval), nrow(testing_)),
  sprintf("Eligibility   : %s", manifest$eligibility),
  sprintf("Split         : committed family-blocked holdout, %s", splits_path(BLOCK_LEVEL)),
  "",
  "The ladder (every rung is the baseline block plus its own structure block):",
  paste0("  ", sprintf("%-10s %3d structure cols, %3d predictors%s",
                       MODELS, lengths(STRUCTURE[MODELS]), lengths(PREDS[MODELS]),
                       if_else(MODELS == PRIMARY_MODEL, "   *PRIMARY", ""))),
  "",
  "Excluded from the baseline, and why:",
  "  aa_*                       exact function of the retained codon columns:",
  "                             aa_x = sum(codons for x)/(1-stops-other).",
  "  frac_*, purine_/amino_*    exact function of GC content and the two skews.",
  "  exon_length_last_mrna      3'UTR-length proxy: Spearman 0.949 with",
  "                             length_3utr, which is itself in the baseline.",
  "  translation_efficiency     measured phenotype, not a sequence feature.",
  "  gini_* (icSHAPE)           80-91% missing. See the secondary run.",
  "",
  strrep("-", 74),
  "Held-out results — each rung against the baseline",
  strrep("-", 74),
  vapply(which(delta_table$kind == "vs_baseline"),
         function(i) fmt_row(delta_table[i, ]), character(1)),
  "",
  strrep("-", 74),
  "Held-out results — rung-to-rung increments",
  strrep("-", 74),
  vapply(which(delta_table$kind == "increment"),
         function(i) fmt_row(delta_table[i, ]), character(1)),
  "",
  "Sign convention: delta = larger model - smaller model.",
  "  R-squared  delta > 0 favours the larger model",
  "  RMSE, MAE  delta < 0 favours the larger model",
  sprintf("95%% CIs from %d paired bootstrap replicates over held-out genes.", N_BOOT),
  "Every rung is scored on the SAME draw in every replicate, so the increments",
  "are within-draw differences and their intervals are intervals on increments.",
  "",
  sprintf("Secondary sign-flip test, primary contrast (squared error): p = %.4g",
          perm_table$p_value[perm_table$primary &
                             perm_table$loss == "squared error"][1]),
  "",
  strrep("-", 74),
  "Interpretation",
  strrep("-", 74),
  strwrap(verdict, width = 74),
  "",
  strwrap(paste("MULTIPLICITY. The ladder produces", nrow(CONTRASTS),
                "contrasts on one held-out set. Only", PRIMARY_MODEL, "vs",
                REFERENCE_MODEL, "is pre-specified; the rest are secondary and",
                "unadjusted, and none of them should be promoted to the",
                "headline after the fact. What the secondary rows are FOR is",
                "the shape of the ladder, not their individual intervals."),
          width = 74),
  "",
  strwrap(paste("READING THE SHAPE. The rungs are ordered from least to most",
                "confounded with the baseline. A real structure effect should",
                "appear at", PRIMARY_MODEL, "and persist up the ladder. A gain",
                "that appears only at S-full, whose block is near-deterministic",
                "in length and GC, is the expected signature of baseline",
                "information re-entering under a structure label."), width = 74),
  "",
  strwrap(paste("This is predictive evidence only. Gain importance is exploratory",
                "context; correlated predictors share and redistribute importance,",
                "so it does not estimate an independent effect. Causal and",
                "inferential claims are for the later modelling plan."), width = 74),
  "",
  strwrap(paste("Limitation to report with this result: the blocked split pins",
                "families larger than ~5% of the smallest split to train, so the",
                "test set contains no large family. It measures generalisation to",
                "small and mid-sized families, not to the largest ones."), width = 74),
  "",
  sprintf("Seed %d. Fitted %s. %s", SEED, format(Sys.time(), "%Y-%m-%d %H:%M"),
          R.version.string)
)

writeLines(summary_txt, file.path(RUN_DIR, "SUMMARY.txt"))
cat("\n"); cat(summary_txt, sep = "\n"); cat("\n")


# ----------------------------- 16. Figures ----------------------------------
# Rendering lives in its own script, which reads only the tables written above
# and never touches a fitted model — so a figure can be restyled in seconds
# instead of re-tuning:
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
