# =============================================================================
# Figure 5: experimental probing vs thermodynamic folding predictions
# =============================================================================
# A: icSHAPE Gini vs robust MFE z-score (full mRNA, cytoplasm). The z-score is
#    (MFE - null median) / (1.4826 * MAD) over esl-shuffle nulls (worker.sh).
# B: probing-MFE correlations across regions and thermodynamic metrics.
# C: associations of probing vs predicted structure with half-life and TE
#    (unadjusted; adjusted estimates are in the table, see f5_plot_c()).
#
# ORIENTATION. Predicted values are plotted as they are. Gini is high when a
# region is more structured; every MFE-type metric is more NEGATIVE when it is.
# So agreement between probing and prediction is a NEGATIVE correlation (panel
# B), and "more structured goes with lower half-life/TE" is a negative rho for
# Gini but a positive rho for the MFE metrics (panel C).
#
# COVERAGE. Gini depends on icSHAPE coverage, and coverage tracks expression,
# which also tracks half-life and TE. icshape_nvalid_* (usable positions) and
# `expression` (Xpresso median RPKM, 57 Roadmap epigenomes) come from the
# dataset; covariates are rank-transformed in f5_partial(), so no log is taken.
#
# COHORT. build_dataset() applies the UTR length filters. Gini exists for only
# the genes with enough icSHAPE coverage (roughly 10-20 % depending on region
# and compartment), so every number here is for that subset. For each region
# the cohort is the genes with Gini AND the predicted metric(s) in that region
# (and, in C, the response and the adjustment covariates), so probing and
# prediction are always compared on identical genes. One transcript per gene
# is kept: Gini is gene-level and is broadcast to a gene's transcripts.
#
# Usage:
#   source("R/load_all.R"); source("analysis/structure/figure5.R")
#   df  <- build_dataset("human")
#   tab <- figure5_tables(df)
# =============================================================================

if (!exists("build_dataset")) source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# correlation_with_ci() lives in the ranked script, which runs its own jobs
# when sourced at top level; load it into a private environment.
.ranked_env <- new.env(parent = globalenv())
sys.source("analysis/correlations/feature_correlation_ranked.R",
           envir = .ranked_env)
correlation_with_ci <- .ranked_env$correlation_with_ci

F5_REGIONS      <- c("mrna", "5utr", "cds", "3utr")
F5_COMPARTMENTS <- c("cytoplasm", "nucleoplasm")
F5_METRICS      <- c("rnafold_score", "rnafold_zscore", "rnafold_per_nt",
                     "mfe_delta", "rnalfold_score", "rnalfold_zscore")
F5_RESPONSES    <- c(halflife = "halflife", te = "translation_efficiency")

#' One row per gene (first transcript), so gene-level Gini and responses are
#' not counted once per isoform.
f5_gene_level <- function(df) df[!duplicated(df$gene_id), , drop = FALSE]

#' Rank-based partial correlation of x and y given covariates `z` (data frame),
#' on complete cases. Ranks are residualised on the ranked covariates; the CI
#' uses the Spearman standard error with n reduced by the number of covariates.
f5_partial <- function(x, y, z, conf = 0.95, min_n = 30) {
  ok <- stats::complete.cases(x, y, z)
  n  <- sum(ok); k <- ncol(z)
  if (n < min_n) return(list(estimate = NA_real_, conf.low = NA_real_,
                             conf.high = NA_real_, p.value = NA_real_, n = n))
  rk <- function(v) rank(v[ok])
  Z  <- as.data.frame(lapply(z, rk))
  rx <- stats::resid(stats::lm(rk(x) ~ ., data = Z))
  ry <- stats::resid(stats::lm(rk(y) ~ ., data = Z))
  r  <- stats::cor(rx, ry)
  se <- 1.06 / sqrt(n - 3 - k)
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  p  <- 2 * stats::pnorm(-abs(atanh(r)) / se)
  list(estimate = r, conf.low = tanh(atanh(r) - zc * se),
       conf.high = tanh(atanh(r) + zc * se), p.value = p, n = n)
}

#' Panel A data: Gini vs MFE z-score, full mRNA.
f5_panel_a <- function(g, compartment = "cytoplasm", region = "mrna") {
  gc <- paste0("gini_", compartment, "_", region)
  zc <- paste0("rnafold_zscore_", region)
  d  <- data.frame(gini = g[[gc]], mfe_z = g[[zc]])
  d  <- d[stats::complete.cases(d), ]
  list(data = d, cor = correlation_with_ci(d$gini, d$mfe_z))
}

#' Panel B table: Gini vs every predicted metric, per region and compartment.
f5_panel_b <- function(g) {
  rows <- list()
  for (comp in F5_COMPARTMENTS) for (r in F5_REGIONS) for (m in F5_METRICS) {
    gc <- paste0("gini_", comp, "_", r); mc <- paste0(m, "_", r)
    if (!all(c(gc, mc) %in% names(g))) next
    ct <- correlation_with_ci(g[[gc]], g[[mc]])
    rows[[length(rows) + 1]] <- data.frame(
      compartment = comp, region = r, metric = m,
      rho = ct$estimate, conf.low = ct$conf.low, conf.high = ct$conf.high,
      p.value = ct$p.value, n = ct$n)
  }
  out <- dplyr::bind_rows(rows)
  out$q.value <- stats::p.adjust(out$p.value, "BH")
  out
}

#' Panel C table: Gini and MFE z-score vs half-life and TE, per region. For
#' each, the plain correlation and the one adjusted for the other measure,
#' log10 region length and C+G ("adjusted"), and additionally for Gini
#' coverage (log10 n_valid) and log expression ("adjusted_coverage"). All on
#' the covariate-complete subset of the plain cohort (n differs by row).
f5_panel_c <- function(g, compartment = "cytoplasm") {
  rows <- list()
  for (r in F5_REGIONS) {
    gi <- g[[paste0("gini_", compartment, "_", r)]]
    st <- g[[paste0("rnafold_zscore_", r)]]
    cov_len <- log10(g[[paste0("length_", r)]])
    cov_gc  <- g[[paste0("gc_content_", r)]]
    cov_nv  <- g[[paste0("icshape_nvalid_", compartment, "_", r)]]
    cov_ex  <- g$expression
    for (rn in names(F5_RESPONSES)) {
      y  <- g[[F5_RESPONSES[[rn]]]]
      # Plain: genes with both measures and the response. Adjusted levels need
      # the covariates too, so they run on that (smaller) subset; n is
      # reported per row.
      ok     <- stats::complete.cases(gi, st, y)
      ok_cov <- ok & stats::complete.cases(cov_len, cov_gc, cov_nv, cov_ex)
      for (meas in c("probing", "prediction")) {
        x     <- if (meas == "probing") gi else st
        other <- if (meas == "probing") st else gi
        plain <- correlation_with_ci(x[ok], y[ok])
        adj   <- f5_partial(x[ok_cov], y[ok_cov],
                            data.frame(other = other[ok_cov],
                                       len = cov_len[ok_cov], gc = cov_gc[ok_cov]))
        adj_cov <- f5_partial(x[ok_cov], y[ok_cov],
                              data.frame(other = other[ok_cov], len = cov_len[ok_cov],
                                         gc = cov_gc[ok_cov], nvalid = cov_nv[ok_cov],
                                         expr = cov_ex[ok_cov]))
        for (kind in c("plain", "adjusted", "adjusted_coverage")) {
          v <- switch(kind, plain = plain, adjusted = adj,
                      adjusted_coverage = adj_cov)
          rows[[length(rows) + 1]] <- data.frame(
            compartment = compartment, region = r, response = rn,
            measure = meas, estimate_type = kind, rho = v$estimate,
            conf.low = v$conf.low, conf.high = v$conf.high,
            p.value = v$p.value, n = v$n)
        }
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  out$q.value <- stats::p.adjust(out$p.value, "BH")
  out
}

figure5_tables <- function(df, compartment = "cytoplasm") {
  g <- f5_gene_level(df)
  list(n_genes = nrow(g),
       a = f5_panel_a(g, compartment),
       b = f5_panel_b(g),
       c = f5_panel_c(g, compartment))
}


# =============================================================================
# Panel A
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

#' Gini vs MFE z-score: hex density (about 1.4k genes is too many points for a
#' scatter to read), a loess trend and the Spearman rho, n and 95 % CI.
#' Gini rises with structure; the MFE z-score falls with it (x axis as is).
#'
#' DISPLAY RANGE. The z-scores are heavy-tailed (median 0.9, 99th percentile
#' 21, maximum 56), which squeezes the bulk into a sliver. The view is cropped
#' to the central `view_q` of each axis; the correlation and the trend use all
#' the genes. The number outside the view is returned as attr(, "n_outside").
f5_plot_a <- function(a, compartment = "cytoplasm", view_q = 0.99,
                      base_size = 16) {
  ct <- a$cor; d <- a$data
  lo <- (1 - view_q) / 2
  xr <- stats::quantile(d$mfe_z, c(lo, 1 - lo))
  yr <- stats::quantile(d$gini, c(lo, 1 - lo))
  inside  <- d$mfe_z >= xr[1] & d$mfe_z <= xr[2] &
             d$gini >= yr[1] & d$gini <= yr[2]
  outside <- sum(!inside)
  lab <- sprintf("\u03c1 = %.2f, n = %s", ct$estimate,
                 format(ct$n, big.mark = ","))
  p <- ggplot(d, aes(mfe_z, gini)) +
    geom_hex(data = d[inside, ], bins = 35) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = "firebrick", linewidth = 0.9) +
    scale_fill_viridis_c(trans = "log10", name = "Genes") +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(title = "Probing and predicted structure", subtitle = lab, x = "Robust MFE z-score (\u2190 more structured)",
         y = sprintf("icSHAPE Gini, %s\n(more structured \u2192)", compartment)) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base_size),
          plot.subtitle = element_text(size = base_size * 0.75))
  attr(p, "n_outside") <- outside
  p
}


# =============================================================================
# Panel B
# =============================================================================
F5_METRIC_LABELS <- c(
  rnafold_score  = "MFE",
  rnafold_zscore = "Robust MFE z-score",
  rnafold_per_nt = "MFE / nt",
  mfe_delta      = "MFE \u0394 (vs C+G, length)",
  rnalfold_score = "Local MFE",
  rnalfold_zscore = "Local robust MFE z-score")
F5_REGION_LABELS <- c(mrna = "mRNA", `5utr` = "5\u2032 UTR", cds = "CDS",
                      `3utr` = "3\u2032 UTR")

#' Heatmap of Spearman rho between Gini and each predicted metric (not
#' sign-changed, so negative = both say "more structured"): metrics down, regions across, one
#' facet per compartment. * marks BH q < 0.05 across the whole table.
f5_plot_b <- function(b, base_size = 15) {
  b$metric_f <- factor(F5_METRIC_LABELS[b$metric],
                       levels = rev(F5_METRIC_LABELS))
  b$region_f <- factor(F5_REGION_LABELS[b$region], levels = F5_REGION_LABELS)
  b$comp_f   <- factor(tools::toTitleCase(b$compartment),
                       levels = c("Cytoplasm", "Nucleoplasm"))
  b$label    <- sprintf("%.2f%s", b$rho, ifelse(b$q.value < 0.05, "*", ""))
  lim <- max(abs(b$rho), na.rm = TRUE)
  ggplot(b, aes(region_f, metric_f, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = base_size / 15 * 4.2) +
    facet_wrap(~comp_f, nrow = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-lim, lim), name = "Spearman \u03c1") +
    labs(x = NULL, y = NULL,
         title = "Probing\u2013prediction correlation by region and metric",
         caption = "* BH q < 0.05") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size))
}


# =============================================================================
# Panel C
# =============================================================================
F5_MEASURE_COLOURS <- c(`icSHAPE Gini` = "#D55E00",
                        `Robust MFE z-score` = "#0072B2")

#' Spearman rho of each measure with half-life and TE, per region, with 95 %
#' CIs. Unadjusted only for now: whether and how to control for length, C+G,
#' coverage and expression is undecided. The adjusted estimates are still
#' computed (f5_panel_c(), figure5_C_associations.csv) and can be plotted by
#' passing `types`. Not sign-changed: more structure is a higher Gini but a lower
#' MFE z-score, so the two measures' rho with a response have opposite signs
#' when they agree.
f5_plot_c <- function(cc, base_size = 15, types = "plain") {
  d <- cc[cc$estimate_type %in% types, ]
  d$measure_f <- factor(ifelse(d$measure == "probing", "icSHAPE Gini",
                               "Robust MFE z-score"),
                        levels = names(F5_MEASURE_COLOURS))
  d$type_f    <- factor(d$estimate_type, levels = c("plain", "adjusted_coverage"))
  d$region_f  <- factor(F5_REGION_LABELS[d$region],
                        levels = rev(F5_REGION_LABELS))
  d$resp_f    <- factor(c(halflife = "Half-life", te = "Translation efficiency")[d$response],
                        levels = c("Half-life", "Translation efficiency"))
  dodge <- position_dodge(width = 0.6)
  p <- ggplot(d, aes(rho, region_f, colour = measure_f, shape = type_f,
                     group = interaction(measure_f, type_f))) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    # separators between regions (the default y grid runs through each
    # region's centre, which reads as belonging to one of its points)
    geom_hline(yintercept = seq_len(nlevels(d$region_f) - 1) + 0.5,
               colour = "grey80", linewidth = 0.4) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.25,
                  position = dodge, orientation = "y", linewidth = 0.6) +
    geom_point(aes(fill = after_scale(ifelse(shape == 21, "white", colour))),
               position = dodge, size = base_size / 15 * 3.2, stroke = 1) +
    facet_wrap(~resp_f, nrow = 1) +
    scale_colour_manual(values = F5_MEASURE_COLOURS, name = NULL) +
    scale_shape_manual(values = c(plain = 16, adjusted_coverage = 21),
                       labels = c(plain = "Unadjusted", adjusted_coverage = "Adjusted"),
                       name = NULL) +
    labs(x = "Spearman \u03c1 with response", y = NULL,
         title = "Association with half-life and translation efficiency") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "grey30", fill = NA,
                                      linewidth = 0.8),
          panel.spacing = unit(5, "mm"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size),
          legend.position = "bottom")
  if (length(types) == 1) p <- p + guides(shape = "none")
  p
}

# =============================================================================
# Run: save the panels, the stitched figure, the tables and the run report
# =============================================================================
f5_save <- function(p, name, width, height) {
  base <- file.path(OUTPUT_DIR, "plots", paste0("figure5_", name))
  ggsave(paste0(base, ".jpg"), p, width = width, height = height,
         units = "mm", dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height,
         units = "mm", device = grDevices::cairo_pdf)
}

f5_run_report <- function(df, tab, pa, path) {
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
  g   <- f5_gene_level(df)
  b   <- tab$b; cc <- tab$c
  md  <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"),
                       paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
                       apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  cohort <- do.call(rbind, lapply(F5_REGIONS, function(r) {
    z  <- g[[paste0("rnafold_zscore_", r)]]
    gi <- g[[paste0("gini_cytoplasm_", r)]]; gn <- g[[paste0("gini_nucleoplasm_", r)]]
    data.frame(region = F5_REGION_LABELS[[r]],
               `cytoplasm Gini` = fmt(sum(!is.na(gi))),
               `cytoplasm Gini + MFE z` = fmt(sum(!is.na(gi) & !is.na(z))),
               `nucleoplasm Gini` = fmt(sum(!is.na(gn))),
               `nucleoplasm Gini + MFE z` = fmt(sum(!is.na(gn) & !is.na(z))),
               check.names = FALSE)
  }))
  cn <- cc |> dplyr::filter(estimate_type == "plain") |>
    dplyr::distinct(region, response, n) |>
    dplyr::mutate(region = F5_REGION_LABELS[region]) |>
    tidyr::pivot_wider(names_from = response, values_from = n) |>
    dplyr::mutate(dplyr::across(-region, fmt))
  cw <- cc |> dplyr::mutate(
      cell = sprintf("%.3f [%.3f, %.3f]", rho, conf.low, conf.high),
      region = F5_REGION_LABELS[region],
      sig = ifelse(q.value < 0.05, "", " (ns)"), cell = paste0(cell, sig)) |>
    dplyr::select(region, response, measure, estimate_type, cell) |>
    tidyr::pivot_wider(names_from = estimate_type, values_from = cell)
  lines <- c(
    "# Figure 5 run report", "",
    sprintf("Generated %s from `analysis/structure/figure5.R`.", Sys.Date()), "",
    "Figure subtitles are deliberately minimal; the numbers behind them are here, for the legends.", "",
    "## Cohort", "",
    sprintf("- Transcripts in the analysed frame (human): **%s** (`build_dataset()`, UTR length filters applied)", fmt(nrow(df))),
    sprintf("- Genes after keeping one transcript per gene: **%s**", fmt(nrow(g))),
    "- icSHAPE Gini exists only for genes with enough coverage; every statistic below is on genes with Gini AND the predicted metric(s) in that region.", "",
    md(cohort), "",
    "## Panel A", "",
    sprintf("- Gini (cytoplasm, mRNA) vs robust MFE z-score: Spearman rho **%.3f** (95%% CI %.3f to %.3f), n **%s**",
            pa$cor$estimate, pa$cor$conf.low, pa$cor$conf.high, fmt(pa$cor$n)),
    sprintf("- Display cropped to the central 99%% of each axis; **%d** genes lie outside the view (all genes enter the correlation and the loess trend)", attr(pa$plot, "n_outside")),
    "- MFE z-score = (MFE - null median) / (1.4826 x MAD) over esl-shuffle nulls, i.e. a robust z-score; more negative = more structured (plotted as is)", "",
    "## Panel B", "",
    sprintf("- **%d** Spearman correlations (6 metrics x 4 regions x 2 compartments), n **%s** to **%s** each; BH-adjusted over all %d; **%d** at q < 0.05",
            nrow(b), fmt(min(b$n)), fmt(max(b$n)), nrow(b), sum(b$q.value < 0.05)),
    "- Predicted metrics are not sign-changed: negative rho = probing and prediction agree (high Gini, low/negative MFE)", "",
    "## Panel C", "",
    "- Cytoplasm Gini and robust MFE z-score vs half-life and translation efficiency, Spearman; the plotted (unadjusted) estimates use the genes with Gini, MFE z-score and the response in that region (n below, same for both measures). Adjusted levels use the subset that also has all covariates (n in the table)",
    "- The figure shows the unadjusted estimates only (confounder control is undecided). The table also carries two adjusted levels: `adjusted` (rank partial correlation given the other measure, log10 region length, C+G) and `adjusted_coverage` (additionally icSHAPE usable positions `icshape_nvalid_*` and expression); neither is plotted",
    "- Partial-correlation CIs use the Spearman standard error 1.06/sqrt(n - 3 - k) (k covariates): an approximation",
    "- BH-adjusted over all 96 tests of the table", "",
    "### n per region and response", "", md(cn), "",
    "### Estimates, rho [95% CI]", "(ns) = BH q >= 0.05", "", md(as.data.frame(cw)), "",
    "## Covariate provenance", "",
    "- `expression`: Xpresso label, median RPKM across 57 Roadmap epigenomes, linear scale (agarwal_features.csv; verified against the Xpresso authors' 57epigenomes.median_expr.txt). A cross-tissue level, not specific to the icSHAPE cell type, so a proxy for coverage; adjusting for it may over-correct if highly expressed genes are truly less structured",
    "- `icshape_nvalid_<compartment>_<region>`: positions with usable icSHAPE reactivity (gini_long.tsv `n_valid`)",
    "- Covariates are rank-transformed before partialling, so no log transform is applied", "")
  writeLines(lines, path)
}

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  for (d in c("plots", "tables", "reports"))
    dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
  df  <- build_dataset("human")
  tab <- figure5_tables(df)
  pa  <- list(cor = tab$a$cor, plot = f5_plot_a(tab$a))
  pb  <- f5_plot_b(tab$b)
  pc  <- f5_plot_c(tab$c)

  suppressPackageStartupMessages(library(patchwork))
  # Stitched figure: smaller text, A and B on top, C across the bottom.
  fs  <- 11
  tag <- theme(plot.tag = element_text(face = "bold", size = 18))
  top  <- (f5_plot_a(tab$a, base_size = fs) + tag |
           f5_plot_b(tab$b, base_size = fs) + tag) +
    plot_layout(widths = c(0.8, 1.5))
  full <- (top / (f5_plot_c(tab$c, base_size = fs) + tag)) +
    plot_layout(heights = c(1, 0.9)) +
    plot_annotation(tag_levels = "A")
  f5_save(pa$plot, "A", 150, 130)
  f5_save(pb,      "B", 260, 110)
  f5_save(pc,      "C", 200, 125)
  f5_save(full,    "full", 280, 230)

  write.csv(tab$a$data, file.path(OUTPUT_DIR, "tables", "figure5_A_points.csv"), row.names = FALSE)
  write.csv(tab$b, file.path(OUTPUT_DIR, "tables", "figure5_B_correlations.csv"), row.names = FALSE)
  write.csv(tab$c, file.path(OUTPUT_DIR, "tables", "figure5_C_associations.csv"), row.names = FALSE)
  f5_run_report(df, tab, pa, file.path(OUTPUT_DIR, "reports", "figure5_report.md"))
  message("Figure 5 written to ", OUTPUT_DIR)
}
