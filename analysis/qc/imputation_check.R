# =============================================================================
# QC: mRNA MFE imputation diagnostics
# =============================================================================
# Checks how well the length-weighted region-based imputation of mRNA z-score
# (and additive region-based imputation of mRNA score) agrees with the
# directly-computed mRNA-level values, on transcripts that have both.
#
# Fixes vs. old pipeline:
#   - The old impute_Zscores.R filtered out rows with any NA region z-score
#     before imputing, which defeats the purpose of imputation. The corrected
#     pipeline imputes from whatever regions are present.
#   - Labels here say "Spearman \u03c1" not "Pearson" — the old code used
#     Spearman but labelled plots as Pearson.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

species <- "human"
# min_5utr = NULL, min_3utr = NULL: QC describes the whole built table, including the short-UTR
# transcripts the cohort filter removes from every analysis. Their missingness
# is exactly what a coverage diagnostic should be able to see. See
# MIN_5UTR_LENGTH / MIN_3UTR_LENGTH in config.R.
df <- build_dataset(species, min_5utr = NULL, min_3utr = NULL)

# impute_mrna_mfe() is switched off in engineer_features(), so the cache holds
# only the directly-folded mRNA values and no *_imputed columns. Run the
# imputation here, on a copy, to see how well it would agree if switched on.
#
# The direct values are captured BEFORE imputing: impute_mrna_mfe() coalesces
# the imputed value into rnafold_*_mrna wherever the direct one is missing,
# and comparing after that would count every filled row as a perfect match.

need <- c(paste0("rnafold_", c("zscore", "score"), "_mrna"),
          paste0("rnafold_", rep(c("zscore", "score"), each = 3), "_",
                 c("5utr", "cds", "3utr")),
          paste0("length_", c("5utr", "cds", "3utr")))
missing <- setdiff(need, names(df))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "),
       ". Rebuild the dataset with build_dataset('", species, "', rebuild = TRUE).")
}

imp <- impute_mrna_mfe(df)

cmp <- tibble(
  rnafold_zscore_mrna         = df$rnafold_zscore_mrna,
  rnafold_zscore_mrna_imputed = imp$rnafold_zscore_mrna_imputed,
  rnafold_score_mrna          = df$rnafold_score_mrna,
  rnafold_score_mrna_imputed  = imp$rnafold_score_mrna_imputed
) |>
  filter(!is.na(rnafold_zscore_mrna), !is.na(rnafold_zscore_mrna_imputed),
         !is.na(rnafold_score_mrna),  !is.na(rnafold_score_mrna_imputed))



# --- Plots 1-2: direct against estimated, hex-binned ------------------------
# Estimate on x, the directly folded value on y (Pineiro et al. 2008): the
# direct value is what the estimate is judged against.
# Same look as the figures in mfe_observed_vs_predicted.R: log-count hexes,
# y = x in orange, Spearman rho top left. Axes are equal so y = x stays the
# diagonal; the view is zoomed past the long tails (ZOOM_KEEP of rows on each
# axis), not filtered, so every row is still binned.
HEX_BINS  <- 60
ZOOM_KEEP <- 0.99

theme_mfe <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 17, face = "bold"),
    panel.background = element_rect(fill = "grey96", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "grey45", linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey88", linewidth = 0.3),
    legend.position  = "right"
  )

hex_vs_direct <- function(x, y, title, xlab, ylab) {
  tail <- (1 - ZOOM_KEEP) / 2
  lo   <- min(quantile(x, tail), quantile(y, tail))
  hi   <- max(quantile(x, 1 - tail), quantile(y, 1 - tail))
  rho  <- stats::cor(x, y, method = "spearman")
  ggplot(tibble(x = x, y = y), aes(x = x, y = y)) +
    geom_hex(binwidth = rep((hi - lo) / HEX_BINS, 2)) +
    geom_abline(slope = 1, intercept = 0, colour = "#D55E00", linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("Spearman \u03c1 = %.3f", rho),
             hjust = -0.08, vjust = 1.6, size = 4.2) +
    coord_cartesian(xlim = c(lo, hi), ylim = c(lo, hi)) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_gradient(low = "#C6DBEF", high = "#08306B", trans = "log10",
                        labels = scales::label_comma(), name = "Transcript\ndensity") +
    labs(title = title, x = xlab, y = ylab) +
    theme_mfe
}

p_score <- hex_vs_direct(
  cmp$rnafold_score_mrna_imputed, cmp$rnafold_score_mrna,
  title = "Observed vs estimated mRNA MFE",
  xlab  = "Sum of MFE components (kcal/mol)",
  ylab  = "Observed mRNA MFE (kcal/mol)")

# plotmath labels, so the z can be italic (bold italic in the bold title)
p_zscore <- hex_vs_direct(
  cmp$rnafold_zscore_mrna_imputed, cmp$rnafold_zscore_mrna,
  title = expression(bold("Observed vs estimated mRNA MFE " * bolditalic(z) * "-scores")),
  xlab  = expression("Length-weighted region " * italic(z) * "-score"),
  ylab  = expression("Observed mRNA MFE " * italic(z) * "-score"))


# --- Plot 3: region MFE delta distributions ---------------------------------
# Distribution of (observed - expected) by region, per the thermodynamic model.
delta_cols <- grep("^mfe_delta_", names(df), value = TRUE)
delta_cols <- setdiff(delta_cols, "mfe_delta_mrna")       # cleaner contrast

if (length(delta_cols) >= 1 && "mfe_delta_mrna" %in% names(df)) {
  long_delta <- df |>
    select(all_of(c(delta_cols, "mfe_delta_mrna"))) |>
    pivot_longer(everything(),
                 names_to  = "region",
                 values_to = "delta") |>
    mutate(region = sub("^mfe_delta_", "", region)) |>
    filter(!is.na(delta))

  p_delta <- ggplot(long_delta, aes(x = delta, fill = region)) +
    geom_density(alpha = 0.5, colour = "black") +
    coord_cartesian(xlim = c(-60, 35)) +
    labs(
      title = "MFE delta (observed - expected) by region",
      x = "MFE delta",
      y = "Density",
      fill = "Region"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12)
    )
  print(p_delta)
}

print(p_score)
print(p_zscore)

dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "plots", sprintf("mfe_mrna_direct_vs_region_sum_%s.jpg", species)),
  plot = p_score, width = 180, height = 150, units = "mm", dpi = 300)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "plots", sprintf("mfe_zscore_mrna_direct_vs_weighted_%s.jpg", species)),
  plot = p_zscore, width = 180, height = 150, units = "mm", dpi = 300)
