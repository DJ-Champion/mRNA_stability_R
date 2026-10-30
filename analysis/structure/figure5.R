# =============================================================================
# Figure 5: experimental probing vs thermodynamic folding predictions
# =============================================================================
# A: icSHAPE Gini vs robust MFE z-score (full mRNA, cytoplasm). The z-score is
#    (MFE - null median) / (1.4826 * MAD) over esl-shuffle nulls (worker.sh).
# B: probing-MFE correlations across regions and thermodynamic metrics, with
#    C+G%, length and (full mRNA only) expression alongside as the
#    confounders the paper is about.
#
# NO ADJUSTMENT. The point of the figure is to show what might be confounding,
# so nothing is adjusted for length, C+G, coverage or expression: every number
# is a plain Spearman correlation.
#
# ORIENTATION. Predicted values are used as they are. Gini is high when a
# region is more structured; every MFE-type metric is more NEGATIVE when it is.
# So agreement between probing and prediction is a NEGATIVE correlation (panel
# B).
#
# COHORT. build_dataset() applies the UTR length filters. Gini exists for only
# the genes with enough icSHAPE coverage (roughly 10-20 % depending on region
# and compartment), so every number here is for that subset. "Enough" is set
# upstream (Gini-from-icSHAPE, lib/gini_core.py, VALID_CUTOFF = 15): a region
# gets a Gini only if more than 15 of its nucleotides have a valid icSHAPE
# score (0-1) and the scores are not all zero. There is no read-depth filter. For each region
# the cohort is the genes with Gini AND the predicted metric(s) in that region
# (and, in C, the response), so probing and prediction are always compared on
# identical genes. One transcript per gene is kept: Gini is gene-level and is
# broadcast to a gene's transcripts.
#
# NAMES. Every label comes from the feature table through the project's naming
# functions (REGION_DISPLAYS for regions, format_metric_name() for metrics,
# measures and responses), so the figure follows the table if names change.
# Nothing is spelled out here.
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
F5_REGION_LABELS <- REGION_DISPLAYS[F5_REGIONS]
F5_COMPARTMENTS <- c("cytoplasm", "nucleoplasm")
F5_METRICS      <- c("rnafold_score", "rnafold_zscore", "rnafold_per_nt",
                     "mfe_delta", "rnalfold_score", "rnalfold_zscore")
F5_SEQ_METRICS  <- c("gc_content", "length", "expression")

# Labels, all from the feature table.
F5_METRIC_LABELS <- stats::setNames(
  c(format_metric_name(paste0(c(F5_METRICS, "gc_content", "length"), "_mrna")),
    format_metric_name("expression")),
  c(F5_METRICS, F5_SEQ_METRICS))
F5_COMPARTMENT_LABELS <- stats::setNames(
  format_metric_name(paste0("gini_", F5_COMPARTMENTS, "_mrna")), F5_COMPARTMENTS)
# Panel A axes: probing (Gini, per compartment) and prediction (MFE z).
f5_measure_labels <- function(compartment = "cytoplasm") c(
  probing    = unname(F5_COMPARTMENT_LABELS[compartment]),
  prediction = format_metric_name("rnafold_zscore_mrna"))

#' One row per gene (first transcript), so gene-level Gini and responses are
#' not counted once per isoform.
f5_gene_level <- function(df) df[!duplicated(df$gene_id), , drop = FALSE]

#' Panel A data: Gini vs MFE z-score, full mRNA.
f5_panel_a <- function(g, compartment = "cytoplasm", region = "mrna") {
  gc <- paste0("gini_", compartment, "_", region)
  zc <- paste0("rnafold_zscore_", region)
  d  <- data.frame(gini = g[[gc]], mfe_z = g[[zc]])
  d  <- d[stats::complete.cases(d), ]
  list(data = d, cor = correlation_with_ci(d$gini, d$mfe_z))
}

#' Panel B table: Gini vs every predicted metric, and vs C+G%, length and
#' expression, per region and compartment. Expression is a property of the
#' gene, so it is only correlated for the full mRNA.
f5_panel_b <- function(g) {
  rows <- list()
  for (comp in F5_COMPARTMENTS) for (r in F5_REGIONS)
    for (m in c(F5_METRICS, F5_SEQ_METRICS)) {
      if (m == "expression" && r != "mrna") next
      gc <- paste0("gini_", comp, "_", r)
      mc <- if (m == "expression") m else paste0(m, "_", r)
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

figure5_tables <- function(df, compartment = "cytoplasm") {
  g <- f5_gene_level(df)
  list(n_genes = nrow(g),
       a = f5_panel_a(g, compartment),
       b = f5_panel_b(g))
}


# =============================================================================
# Panel A
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

#' Gini (x) vs robust MFE z-score (y): hex density (about 1.4k genes is too
#' many points for a scatter to read), with the Spearman rho and n.
#' Gini rises with structure; the MFE z-score falls with it.
#'
#' The axes show the full range of both measures (no cropping), so the heavy
#' right tail of the MFE z-scores is visible.
f5_plot_a <- function(a, compartment = "cytoplasm", base_size = 16) {
  ct <- a$cor; d <- a$data
  lab <- sprintf("ρ = %.2f, n = %s", ct$estimate,
                 format(ct$n, big.mark = ","))
  p <- ggplot(d, aes(gini, mfe_z)) +
    geom_hex(bins = 35) +
    scale_fill_viridis_c(trans = "log10", name = "Genes") +
    labs(title = "Probing and predicted structure", subtitle = lab,
         x = f5_measure_labels(compartment)[["probing"]],
         y = f5_measure_labels(compartment)[["prediction"]]) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base_size),
          plot.subtitle = element_text(size = base_size * 0.75))
  p
}


# =============================================================================
# Panel B
# =============================================================================
#' Heatmap of Spearman rho between Gini and each predicted metric (not
#' sign-changed, so negative = both say "more structured"), then C+G%, length
#' and expression in a second block below (expression: full mRNA only, other
#' regions left blank): metrics down, regions across, one facet
#' per compartment. * marks BH q < 0.05 across the whole table.
f5_plot_b <- function(b, base_size = 15) {
  b$metric_f <- factor(F5_METRIC_LABELS[b$metric],
                       levels = rev(F5_METRIC_LABELS))
  b$block_f  <- factor(ifelse(b$metric %in% F5_SEQ_METRICS, "Covariates",
                              "Predicted structure"),
                       levels = c("Predicted structure", "Covariates"))
  b$region_f <- factor(F5_REGION_LABELS[b$region], levels = F5_REGION_LABELS)
  b$comp_f   <- factor(F5_COMPARTMENT_LABELS[b$compartment],
                       levels = F5_COMPARTMENT_LABELS)
  b$label    <- sprintf("%.2f%s", b$rho, ifelse(b$q.value < 0.05, "*", ""))
  lim <- max(abs(b$rho), na.rm = TRUE)
  ggplot(b, aes(region_f, metric_f, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = base_size / 15 * 4.2) +
    facet_grid(block_f ~ comp_f, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-lim, lim), name = "Spearman ρ") +
    labs(x = NULL, y = NULL,
         title = "Probing and predicted structure correlation",
         caption = "* BH q < 0.05") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid = element_blank(),
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_blank(),
          plot.title = element_text(face = "bold", size = base_size))
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
  b   <- tab$b
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
  lines <- c(
    "# Figure 5 run report", "",
    sprintf("Generated %s from `analysis/structure/figure5.R`.", Sys.Date()), "",
    "Figure subtitles are deliberately minimal; the numbers behind them are here, for the legends.", "",
    "No correlation is adjusted for length, C+G, coverage or expression: the figure is about what might confound.", "",
    "## Cohort", "",
    sprintf("- Transcripts in the analysed frame (human): **%s** (`build_dataset()`, UTR length filters applied)", fmt(nrow(df))),
    sprintf("- Genes after keeping one transcript per gene: **%s**", fmt(nrow(g))),
    "- icSHAPE Gini exists only for regions with more than 15 nucleotides carrying a valid icSHAPE score (and non-zero total), per Gini-from-icSHAPE; every statistic below is on genes with Gini AND the predicted metric(s) in that region.", "",
    md(cohort), "",
    "## Panel A", "",
    sprintf("- Gini (cytoplasm, mRNA) vs robust MFE z-score: Spearman rho **%.3f** (95%% CI %.3f to %.3f), n **%s**",
            pa$cor$estimate, pa$cor$conf.low, pa$cor$conf.high, fmt(pa$cor$n)),
    "- Axes show the full data range (no cropping)",
    "- MFE z-score = (MFE - null median) / (1.4826 x MAD) over esl-shuffle nulls, i.e. a robust z-score; more negative = more structured (plotted as is)", "",
    "## Panel B", "",
    sprintf("- **%d** Spearman correlations (6 predicted metrics, C+G%% and length x 4 regions x 2 compartments, plus expression for the full mRNA), n **%s** to **%s** each; BH-adjusted over all %d; **%d** at q < 0.05",
            nrow(b), fmt(min(b$n)), fmt(max(b$n)), nrow(b), sum(b$q.value < 0.05)),
    "- Not sign-changed: negative rho = probing and prediction agree (high Gini, low/negative MFE)",
    "- C+G% and length are the region's own (`gc_content_<region>`, `length_<region>`); length is not log-transformed (Spearman is rank-based)",
    "- Expression (Xpresso median RPKM) is gene-level, so it is correlated for the full mRNA only; other regions are blank", "",
    "## Provenance", "",
    "- `expression`: Xpresso label, median RPKM across 57 Roadmap epigenomes, linear scale (agarwal_features.csv; verified against the Xpresso authors' 57epigenomes.median_expr.txt). A cross-tissue level, not specific to the icSHAPE cell type", "")
  writeLines(lines, path)
}

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  for (d in c("plots", "tables", "reports"))
    dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
  df  <- build_dataset("human")
  tab <- figure5_tables(df)
  pa  <- list(cor = tab$a$cor, plot = f5_plot_a(tab$a))
  pb  <- f5_plot_b(tab$b)

  suppressPackageStartupMessages(library(patchwork))
  # Stitched figure: smaller text, A and B side by side.
  fs  <- 11
  tag <- theme(plot.tag = element_text(face = "bold", size = 18))
  full <- (f5_plot_a(tab$a, base_size = fs) + tag |
           f5_plot_b(tab$b, base_size = fs) + tag) +
    plot_layout(widths = c(0.8, 1.5)) +
    plot_annotation(tag_levels = "A")
  f5_save(pa$plot, "A", 150, 130)
  f5_save(pb,      "B", 260, 130)
  f5_save(full,    "full", 280, 130)

  write.csv(tab$a$data, file.path(OUTPUT_DIR, "tables", "figure5_A_points.csv"), row.names = FALSE)
  write.csv(tab$b, file.path(OUTPUT_DIR, "tables", "figure5_B_correlations.csv"), row.names = FALSE)
  f5_run_report(df, tab, pa, file.path(OUTPUT_DIR, "reports", "figure5_report.md"))
  message("Figure 5 written to ", OUTPUT_DIR)
}
