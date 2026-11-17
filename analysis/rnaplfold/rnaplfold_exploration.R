# =============================================================================
# RNAplfold start-codon accessibility: is it worth modelling?  (step 1)
# =============================================================================
# Exploratory only. Asks, before any model is fitted:
#   1. coverage: how many cohort rows carry a value
#   2. redundancy among the five accessibility columns
#   3. how much of each column length / GC / uORF already explain (R^2)
#   4. Spearman vs each response, raw and partial (adjusting for those)
#   5. overlap with the start-region folding features already in the table
#
# Usage: Rscript analysis/rnaplfold/rnaplfold_exploration.R
# Writes CSVs to data/outputs/rnaplfold_exploration/.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("R/load_all.R")

OUT <- file.path(OUTPUT_DIR, "rnaplfold_exploration")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

df <- build_dataset("human") |> drop_excluded(verbose = FALSE)

acc       <- grep("^rnaplfold_", names(df), value = TRUE)
responses <- c("halflife", "translation_efficiency")
# What the baseline already knows about a 5' UTR / start region.
adjust    <- c("length_5utr", "length_cds", "length_3utr", "gc_content_5utr",
               "gc_content_cds", "gc_content_3utr", "gc_content_start",
               "uorf_present_mrna")
folding   <- c("rnafold_zscore_start", "rnafold_score_start",
               "rnalfold_zscore_start", "mfe_delta_start",
               "rnafold_zscore_5utr", "rnafold_zscore_mrna")
adjust  <- intersect(adjust, names(df)); folding <- intersect(folding, names(df))
df[adjust] <- lapply(df[adjust], as.numeric)

# --- 1. Coverage --------------------------------------------------------------
cov <- data.frame(n_rows = nrow(df), n_value = sum(!is.na(df[[acc[1]]])))
cov$share_na <- 1 - cov$n_value / cov$n_rows
cat("1. coverage\n"); print(cov)
na_by <- df |> group_by(has_value = !is.na(.data[[acc[1]]])) |>
  summarise(n = n(), len5 = median(length_5utr, na.rm = TRUE),
            halflife = mean(halflife, na.rm = TRUE))
print(na_by)

# --- 2. Redundancy ------------------------------------------------------------
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman",
                                          use = "pairwise.complete.obs"))
red <- cor(df[acc], method = "spearman", use = "pairwise.complete.obs")
cat("\n2. Spearman among accessibility columns\n"); print(round(red, 2))

# --- 3. Variance explained by what the baseline already carries ---------------
r2 <- sapply(acc, function(v) {
  d <- df[c(v, adjust)]; d <- d[complete.cases(d), ]
  d[] <- lapply(d, rank)
  summary(lm(reformulate(adjust, v), d))$r.squared
})
cat("\n3. rank-R^2 of each column on length / GC / uORF\n"); print(round(r2, 3))

# --- 4. Association with responses: raw and partial ---------------------------
partial <- function(v, y) {
  d <- df[c(v, y, adjust)]; d <- d[complete.cases(d), ]
  d[] <- lapply(d, rank)
  rv <- resid(lm(reformulate(adjust, v), d))
  ry <- resid(lm(reformulate(adjust, y), d))
  c(n = nrow(d), raw = sp(d[[v]], d[[y]]), partial = cor(rv, ry))
}
assoc <- do.call(rbind, lapply(responses, function(y)
  do.call(rbind, lapply(acc, function(v)
    data.frame(response = y, feature = v, t(partial(v, y)))))))
cat("\n4. Spearman vs response (partial = adjusted for length/GC/uORF)\n")
print(assoc |> mutate(across(c(raw, partial), ~ round(., 3))), row.names = FALSE)

# --- 5. Overlap with folding features ----------------------------------------
fold <- outer(acc, folding, Vectorize(function(a, f) sp(df[[a]], df[[f]])))
dimnames(fold) <- list(acc, folding)
cat("\n5. Spearman with existing folding features\n"); print(round(fold, 2))

# Same comparison, for scale: the folding features' own association.
ref <- do.call(rbind, lapply(responses, function(y)
  do.call(rbind, lapply(folding, function(v)
    data.frame(response = y, feature = v, t(partial(v, y)))))))
cat("\nfor scale: folding features vs response\n")
print(ref |> mutate(across(c(raw, partial), ~ round(., 3))), row.names = FALSE)

write.csv(assoc, file.path(OUT, "association.csv"), row.names = FALSE)
write.csv(red,   file.path(OUT, "redundancy.csv"))
write.csv(data.frame(feature = acc, rank_r2 = r2), file.path(OUT, "baseline_r2.csv"),
          row.names = FALSE)
write.csv(fold,  file.path(OUT, "folding_overlap.csv"))
