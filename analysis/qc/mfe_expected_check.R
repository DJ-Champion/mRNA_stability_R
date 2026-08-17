# =============================================================================
# QC: Observed vs expected MFE
# =============================================================================
# Diagnostic for the in-house thermodynamic model (R/features/mfe_model.R):
#   Expected MFE = (a · GC^b + c) · length + d,  GC as a fraction in [0, 1]
#
# If the model is well-calibrated, observed and expected should lie near the
# y = x line. These plots also colour by GC content and length to expose any
# remaining residual structure.
#
# The numeric summary below is the part that actually catches things. A units
# or scale error in the model is invisible to the eye — the scatter still
# looks like a tight band — but it shows up immediately as an observed/expected
# slope away from 1. It was a silent GC-units error (gc/100 applied to a column
# already stored as a fraction) that flattened the GC term to nothing and left
# mfe_delta_* encoding little but length; the slope check is what would have
# surfaced it, so run the summary, not just the plots.
#
# Run this *after* build_dataset() has populated the cache.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

species <- "human"                              # change as needed
region  <- "mrna"                               # 5utr / cds / 3utr / mrna
df <- build_dataset(species)

score_col <- paste0("rnafold_score_", region)
exp_col   <- paste0("mfe_expected_",  region)
gc_col    <- paste0("gc_content_",    region)
len_col   <- paste0("length_",        region)
delta_col <- paste0("mfe_delta_",     region)

req <- c(score_col, exp_col, gc_col, len_col)
missing <- setdiff(req, names(df))
if (length(missing)) {
  stop("missing columns: ", paste(missing, collapse = ", "),
       "\n  (gc is `gc_content_<region>`, not `gc_<region>`)")
}


# --- Numeric calibration summary, all regions --------------------------------
#' Report observed/expected slope and residual length-dependence per region.
#'
#' slope ~ 1 means the model tracks observed MFE at the right magnitude.
#' rho(delta, length) ~ 0 means delta is a genuine residual rather than a
#' length proxy wearing a structure label.
mfe_calibration <- function(df, regions = c("5utr", "cds", "3utr", "mrna")) {
  purrr::map_dfr(regions, function(r) {
    sc <- paste0("rnafold_score_", r); ex <- paste0("mfe_expected_", r)
    dl <- paste0("mfe_delta_", r);     ln <- paste0("length_", r)
    if (!all(c(sc, ex, dl, ln) %in% names(df))) return(tibble::tibble())

    ok <- stats::complete.cases(df[, c(sc, ex, dl, ln)])
    if (sum(ok) < 30) return(tibble::tibble())

    tibble::tibble(
      region      = r,
      n           = sum(ok),
      slope       = unname(stats::coef(stats::lm(df[[sc]][ok] ~ df[[ex]][ok]))[2]),
      obs_per_nt  = stats::median(df[[sc]][ok] / df[[ln]][ok]),
      exp_per_nt  = stats::median(df[[ex]][ok] / df[[ln]][ok]),
      rho_delta_length = stats::cor(df[[dl]][ok], df[[ln]][ok],
                                    method = "spearman")
    )
  })
}

calib <- mfe_calibration(df)

message("\nMFE model calibration — ", species)
print(as.data.frame(calib), digits = 3, row.names = FALSE)

# Slope far from 1 is the units/scale signature; near-|1| delta-vs-length means
# delta is still mostly length. Warn rather than stop so the plots still render.
bad_slope <- calib$region[abs(calib$slope - 1) > 0.25]
if (length(bad_slope)) {
  warning("observed/expected slope is off by >25% in: ",
          paste(bad_slope, collapse = ", "),
          " — check the GC units reaching calculate_mfe_expected()",
          call. = FALSE)
}
leaky <- calib$region[abs(calib$rho_delta_length) > 0.5]
if (length(leaky)) {
  message("  note: mfe_delta still tracks length in: ",
          paste(leaky, collapse = ", "),
          " — the excess over the random null is multiplicative there, so a",
          " ratio or z-score form is a better length-free metric than delta.")
}

complete <- df |>
  filter(!is.na(.data[[score_col]]), !is.na(.data[[exp_col]]))

rho <- cor(complete[[score_col]], complete[[exp_col]], method = "spearman")
rho_label <- sprintf("Spearman \u03c1 = %.3f", rho)

x_lims <- quantile(complete[[score_col]], probs = c(0.01, 0.99), na.rm = TRUE)
y_lims <- quantile(complete[[exp_col]],   probs = c(0.01, 0.99), na.rm = TRUE)

# --- Plot 1: coloured by GC ---------------------------------------------------
plot_by_gc <- ggplot(complete,
                     aes(x = .data[[score_col]], y = .data[[exp_col]])) +
  geom_point(aes(colour = .data[[gc_col]]), alpha = 0.8, size = 1) +
  scale_colour_viridis_c(option = "plasma") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "red", linewidth = 1) +
  annotate("text", x = -Inf, y = Inf, label = rho_label,
           hjust = -0.2, vjust = 2, size = 4) +
  coord_cartesian(xlim = x_lims, ylim = y_lims) +
  labs(
    title    = paste0("Observed vs modelled MFE — ", format_col_name(region)),
    subtitle = bquote(Modelled ~ MFE == (-0.84 ~ GC^2.35 - 0.15) %.% length + 13.56 ~
                      ~ ~ "(GC as a fraction)"),
    x = format_col_name(score_col),
    y = format_col_name(exp_col),
    colour = format_col_name(gc_col)
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14)
  )

# --- Plot 2: coloured by length -----------------------------------------------
plot_by_length <- plot_by_gc %+%
  complete +
  aes(colour = .data[[len_col]])
plot_by_length <- plot_by_length +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = format_col_name(len_col))

# --- Plot 3: delta vs length, coloured by GC ---------------------------------
if (delta_col %in% names(complete)) {
  plot_delta_vs_length <- ggplot(complete,
                                 aes(x = .data[[len_col]],
                                     y = .data[[delta_col]])) +
    geom_point(aes(colour = .data[[gc_col]]), alpha = 0.7, size = 2) +
    scale_colour_viridis_c(option = "viridis") +
    labs(
      title  = paste0("Length vs MFE \u0394 — ", format_col_name(region)),
      x      = format_col_name(len_col),
      y      = format_col_name(delta_col),
      colour = format_col_name(gc_col)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14)
    )
}


# --- Save (R8: outputs go under OUTPUT_DIR) ----------------------------------
# Written to disk rather than only printed, so running this under Rscript does
# not drop a stray Rplots.pdf in the repo root.
dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

qc_plots <- list(
  obs_vs_exp_by_gc     = plot_by_gc,
  obs_vs_exp_by_length = plot_by_length
)
if (exists("plot_delta_vs_length")) {
  qc_plots$delta_vs_length <- plot_delta_vs_length
}

for (nm in names(qc_plots)) {
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots",
              sprintf("qc_mfe_%s_%s_%s.jpg", nm, species, region)),
    plot = qc_plots[[nm]], width = 250, height = 180, units = "mm", dpi = 300
  )
}

write.csv(calib,
          file.path(OUTPUT_DIR, "tables",
                    sprintf("qc_mfe_calibration_%s.csv", species)),
          row.names = FALSE)

message("\nMFE QC complete:")
message("  ", file.path(OUTPUT_DIR, "plots"),  "/qc_mfe_*_", species, "_", region, ".jpg")
message("  ", file.path(OUTPUT_DIR, "tables"), "/qc_mfe_calibration_", species, ".csv")
