# =============================================================================
# Diagnostic: is the complete_case S-select deficit features, or tuning?
# =============================================================================
# In the `complete_case` run, S-select vs Baseline is conclusive on all three
# metrics IN THE DIRECTION FAVOURING THE BASELINE — the only contrast anywhere
# in the study whose interval excludes zero. Three things argue it is noise
# rather than a finding (it does not replicate in `default`; it does not persist
# at S-full, which contains every mfe_delta column; it is 1 of 8 unadjusted
# secondary contrasts). This script tests the fourth explanation.
#
# THE HYPOTHESIS. Each rung tunes independently, so a rung can land on a poor
# configuration by luck of the racing draw. In complete_case, S-core and S-full
# both selected `mod16` (498 trees, lambda 0.16, mtry 0.40) while S-select
# selected `mod54` (345 trees, lambda 45.8, mtry 0.92) — a far more heavily
# regularised model. If the deficit follows the CONFIGURATION rather than the
# FEATURES, it is tuning variance.
#
# THE TEST. Fit all four rungs under each configuration in turn and score them
# on the same held-out genes. Within one configuration block the only thing that
# differs between rungs is the feature set, so the ladder is read feature-only.
#   - deficit present in BOTH blocks  -> the mfe_delta columns really do hurt
#   - deficit present in NEITHER      -> it was the tuning draw
#
# WHAT THIS IS NOT. It is not a replacement estimate and must never be quoted as
# one. Both configurations were chosen after seeing the result, and transplanting
# one rung's winner onto another rung is not the pre-specified design — under
# which every rung is entitled to its own search. This answers "what explains
# that row", not "what is the effect of structure". The headline stays with the
# tuned runs.
#
#   Rscript analysis/models/xgb_structure_tuning_diagnostic.R
#
# Writes complete_case/tables/xgb_structure_tuning_diagnostic.csv and prints the
# comparison. No figure: this is a diagnostic, not a result.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(xgboost)
})

source("R/load_all.R")
source("analysis/models/xgb_structure_features.R")

tidymodels_prefer(quiet = TRUE)

SEED       <- 42L
N_BOOT     <- 2000L
N_THREADS  <- 12L     # sequential fits here, so xgboost gets every core
BOOT_METRICS <- c("rsq_trad", "rmse", "mae")

# The rung under investigation, and the rung whose configuration is borrowed.
SUSPECT   <- "S-select"
REFERENCE <- "Baseline"

VARIANT   <- resolve_variant("complete_case")
RUN_DIR   <- variant_dir(VARIANT, "root")
TABLE_DIR <- variant_dir(VARIANT, "tables")

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  utils::flush.console()
}


# ----------------------------- Inputs ---------------------------------------
# The configurations come from the run manifest rather than being retyped, so
# this cannot drift from the models that produced the result being explained.

man <- readRDS(file.path(RUN_DIR, "run_manifest.rds"))
if (is.null(man$best_params)) stop("manifest has no best_params", call. = FALSE)

el       <- eligible_dataset("human", VARIANT)
dat      <- el$data
MODELS   <- names(el$models)
PREDS    <- lapply(MODELS, function(m) predictors_for(el, m))
names(PREDS) <- MODELS

trainval <- dat |> filter(split != "test")
testing_ <- dat |> filter(split == "test")

# Which configurations to try: the suspect's own, and the one its neighbours
# chose. Named by the rung they came from, deduplicated by .config so that
# neighbours agreeing (as S-core and S-full do here) costs nothing.
cfg_source <- c(SUSPECT, setdiff(MODELS, c(SUSPECT, REFERENCE)))
CONFIGS <- man$best_params[cfg_source] |>
  bind_rows(.id = "from") |>
  distinct(.config, .keep_all = TRUE)

log_msg("Configurations under test:")
print(as.data.frame(CONFIGS), row.names = FALSE)

log_msg(nrow(CONFIGS), " configurations x ", length(MODELS), " rungs = ",
        nrow(CONFIGS) * length(MODELS), " fits on ", nrow(trainval),
        " genes, scored on ", nrow(testing_))


# ----------------------------- Fitting --------------------------------------
# Same recipe as the main script: no centring, no scaling, no step_corr, and no
# outcome-driven selection. step_zv is the same safety net. Anything else here
# would make the diagnostic incomparable with the run it is diagnosing.

build_recipe <- function(preds) {
  recipe(x = dat[0, c(TARGET_COL, "gene_id", "split", preds)]) |>
    update_role(all_of(preds),      new_role = "predictor") |>
    update_role(all_of(TARGET_COL), new_role = "outcome") |>
    update_role(gene_id, split,     new_role = "id") |>
    step_zv(all_predictors())
}

# Values spliced in, not tuned. lambda is an engine argument and must be !!'d
# for the same reason the main script splices N_THREADS: engine args are stored
# as quosures and evaluated later, in a scope that has no `p`.
spec_from <- function(p) {
  boost_tree(trees          = p$trees,
             tree_depth     = p$tree_depth,
             min_n          = p$min_n,
             learn_rate     = p$learn_rate,
             loss_reduction = p$loss_reduction,
             sample_size    = p$sample_size,
             mtry           = p$mtry) |>
    set_engine("xgboost",
               nthread     = !!N_THREADS,
               counts      = FALSE,
               tree_method = "hist",
               lambda      = !!p$lambda) |>
    set_mode("regression")
}

grid <- crossing(config = CONFIGS$.config, model = MODELS)

preds_mat <- matrix(NA_real_, nrow = nrow(testing_), ncol = nrow(grid),
                    dimnames = list(NULL, paste(grid$model, grid$config)))

for (i in seq_len(nrow(grid))) {
  p  <- CONFIGS |> filter(.config == grid$config[i])
  wf <- workflow() |>
    add_recipe(build_recipe(PREDS[[grid$model[i]]])) |>
    add_model(spec_from(p))

  set.seed(SEED)
  fitted <- fit(wf, trainval)
  preds_mat[, i] <- predict(fitted, testing_)$.pred
  log_msg("  fitted ", grid$model[i], " under ", grid$config[i],
          " (from ", p$from, ")")
}


# ----------------------------- Paired bootstrap -----------------------------
# Identical machinery and seed to the main script, so an interval here is
# directly comparable with an interval there. Every column is scored on the same
# resampled gene list in every replicate, which is what makes the within-block
# rung differences paired.

obs <- testing_[[TARGET_COL]]
SE  <- (preds_mat - obs)^2
AE  <- abs(preds_mat - obs)

boot_stat <- function(idx) {
  o   <- obs[idx]
  sst <- sum((o - mean(o))^2)
  rbind(rsq_trad = 1 - colSums(SE[idx, , drop = FALSE]) / sst,
        rmse     = sqrt(colMeans(SE[idx, , drop = FALSE])),
        mae      = colMeans(AE[idx, , drop = FALSE]))
}

n_test <- nrow(preds_mat)
point  <- boot_stat(seq_len(n_test))

set.seed(SEED)
boot_arr <- array(NA_real_, dim = c(length(BOOT_METRICS), ncol(preds_mat), N_BOOT),
                  dimnames = list(BOOT_METRICS, colnames(preds_mat), NULL))
for (i in seq_len(N_BOOT)) {
  boot_arr[, , i] <- boot_stat(sample.int(n_test, n_test, replace = TRUE))
}

boot_ci <- function(metric, lhs, rhs, p) {
  unname(quantile(boot_arr[metric, lhs, ] - boot_arr[metric, rhs, ], p))
}

# The contrast that produced the finding, recomputed inside each configuration
# block. Same rungs, same genes, same bootstrap — the only thing held fixed that
# was not fixed before is the hyperparameter configuration.
contrasts <- CONFIGS |>
  transmute(config = .config,
            from,
            lhs = paste(SUSPECT, .config),
            rhs = paste(REFERENCE, .config))

delta_table <- contrasts |>
  crossing(metric = BOOT_METRICS) |>
  mutate(
    contrast   = sprintf("%s vs %s (both under %s, %s's config)",
                         SUSPECT, REFERENCE, config, from),
    lhs_value  = map2_dbl(metric, lhs, ~ point[.x, .y]),
    rhs_value  = map2_dbl(metric, rhs, ~ point[.x, .y]),
    delta      = lhs_value - rhs_value,
    ci_low     = pmap_dbl(list(metric, lhs, rhs), boot_ci, p = 0.025),
    ci_high    = pmap_dbl(list(metric, lhs, rhs), boot_ci, p = 0.975),
    conclusive = (ci_low > 0) | (ci_high < 0),
    metric     = factor(metric, levels = BOOT_METRICS)
  ) |>
  arrange(metric, config)

write_csv(delta_table |> select(config, config_from = from, contrast, metric,
                                lhs_value, rhs_value, delta, ci_low, ci_high,
                                conclusive),
          file.path(TABLE_DIR, "xgb_structure_tuning_diagnostic.csv"))

cat("\n=== Held-out R-squared, every rung under every configuration ===\n")
as_tibble(t(point["rsq_trad", , drop = FALSE]), rownames = "fit") |>
  separate_wider_regex(fit, c(model = ".*", " ", config = "pre.*")) |>
  pivot_wider(names_from = config, values_from = rsq_trad) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 4)

cat("\n=== ", SUSPECT, " vs ", REFERENCE,
    " within each configuration, 95% paired-bootstrap CI ===\n", sep = "")
delta_table |>
  select(metric, config, from, delta, ci_low, ci_high, conclusive) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\nAs tuned (each rung its own winner), for reference:\n")
read_csv(file.path(TABLE_DIR, "xgb_structure_delta_bootstrap.csv"),
         show_col_types = FALSE) |>
  filter(contrast == paste(SUSPECT, "vs", REFERENCE)) |>
  select(metric, delta, ci_low, ci_high, conclusive) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

n_conc <- sum(delta_table$conclusive[delta_table$metric == "rsq_trad"])
cat("\n", n_conc, " of ", nrow(CONFIGS),
    " configurations still show a conclusive R-squared deficit.\n", sep = "")
cat(if (n_conc == 0) {
  paste0("Reading: the deficit does not survive holding the configuration\n",
         "fixed, so it is tuning variance, not the mfe_delta columns.\n")
} else if (n_conc == nrow(CONFIGS)) {
  paste0("Reading: the deficit survives under every configuration tried, so\n",
         "it is not explained by the tuning draw. The mfe_delta block needs\n",
         "a proper look.\n")
} else {
  paste0("Reading: mixed — the deficit depends on the configuration, which is\n",
         "itself evidence that tuning variance is large relative to it.\n")
})

log_msg("Written to ", file.path(TABLE_DIR, "xgb_structure_tuning_diagnostic.csv"))
