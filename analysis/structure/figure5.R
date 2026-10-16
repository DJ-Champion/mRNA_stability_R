# =============================================================================
# Figure 5: experimental probing vs thermodynamic folding predictions
# =============================================================================
# A: icSHAPE Gini vs robust MFE z-score (full mRNA, cytoplasm). The z-score is
#    (MFE - null median) / (1.4826 * MAD) over esl-shuffle nulls (worker.sh).
# B: probing-MFE correlations across regions and thermodynamic metrics, with
#    C+G% and length alongside as the confounders the paper is about.
# C: absolute association of probing, predicted structure (and, for the full
#    mRNA, expression) with half-life and TE.
#
# NO ADJUSTMENT. The point of the figure is to show what might be confounding,
# so nothing is adjusted for length, C+G, coverage or expression: every number
# is a plain Spearman correlation.
#
# ORIENTATION. Predicted values are used as they are. Gini is high when a
# region is more structured; every MFE-type metric is more NEGATIVE when it is.
# So agreement between probing and prediction is a NEGATIVE correlation (panel
# B), and "more structured goes with lower half-life/TE" is a negative rho for
# Gini but a positive rho for the MFE z-score (panel C plots |rho|; the signed
# values are in the table).
#
# COHORT. build_dataset() applies the UTR length filters. Gini exists for only
# the genes with enough icSHAPE coverage (roughly 10-20 % depending on region
# and compartment), so every number here is for that subset. For each region
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
F5_SEQ_METRICS  <- c("gc_content", "length")
F5_RESPONSES    <- c(halflife = "halflife", te = "translation_efficiency")

# Labels, all from the feature table.
F5_METRIC_LABELS <- stats::setNames(
  format_metric_name(paste0(c(F5_METRICS, F5_SEQ_METRICS), "_mrna")),
  c(F5_METRICS, F5_SEQ_METRICS))
F5_RESPONSE_LABELS <- stats::setNames(format_metric_name(unname(F5_RESPONSES)),
                                      names(F5_RESPONSES))
F5_COMPARTMENT_LABELS <- stats::setNames(
  format_metric_name(paste0("gini_", F5_COMPARTMENTS, "_mrna")), F5_COMPARTMENTS)
# Panel A / C measures: probing (Gini, per compartment), prediction (MFE z),
# expression. Colours are keyed by role, labels by the table.
f5_measure_labels <- function(compartment = "cytoplasm") c(
  probing    = unname(F5_COMPARTMENT_LABELS[compartment]),
  prediction = format_metric_name("rnafold_zscore_mrna"),
  expression = format_metric_name("expression"))

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

#' Panel B table: Gini vs every predicted metric, and vs C+G% and length, per
#' region and compartment.
f5_panel_b <- function(g) {
  rows <- list()
  for (comp in F5_COMPARTMENTS) for (r in F5_REGIONS)
    for (m in c(F5_METRICS, F5_SEQ_METRICS)) {
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

#' Panel C table: Gini and MFE z-score vs half-life and TE, per region, on the
#' genes that have both measures and the response. Gene expression is added for
#' the full mRNA only, on those of the same genes that also have a value. Also
#' carries |rho| and its interval folded as in feature_correlation_ranked():
#' an interval that spans zero starts at 0.
f5_panel_c <- function(g, compartment = "cytoplasm") {
  rows <- list()
  add <- function(r, rn, meas, ct) {
    rows[[length(rows) + 1]] <<- data.frame(
      compartment = compartment, region = r, response = rn, measure = meas,
      rho = ct$estimate, conf.low = ct$conf.low, conf.high = ct$conf.high,
      p.value = ct$p.value, n = ct$n)
  }
  for (r in F5_REGIONS) {
    gi <- g[[paste0("gini_", compartment, "_", r)]]
    st <- g[[paste0("rnafold_zscore_", r)]]
    for (rn in names(F5_RESPONSES)) {
      y  <- g[[F5_RESPONSES[[rn]]]]
      ok <- stats::complete.cases(gi, st, y)
      add(r, rn, "probing",    correlation_with_ci(gi[ok], y[ok]))
      add(r, rn, "prediction", correlation_with_ci(st[ok], y[ok]))
      if (r == "mrna") {
        ok_ex <- ok & !is.na(g$expression)
        add(r, rn, "expression", correlation_with_ci(g$expression[ok_ex], y[ok_ex]))
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  out$q.value <- stats::p.adjust(out$p.value, "BH")
  out$rho_abs <- abs(out$rho)
  out$conf.low_abs <- dplyr::case_when(
    out$conf.low < 0 & out$conf.high > 0 ~ 0,
    out$conf.low < 0 & out$conf.high < 0 ~ abs(out$conf.high),
    TRUE ~ out$conf.low)
  out$conf.high_abs <- dplyr::case_when(
    out$conf.low < 0 & out$conf.high > 0 ~ pmax(abs(out$conf.low), abs(out$conf.high)),
    out$conf.low < 0 & out$conf.high < 0 ~ abs(out$conf.low),
    TRUE ~ out$conf.high)
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

#' Gini (x) vs robust MFE z-score (y): hex density (about 1.4k genes is too
#' many points for a scatter to read), with the Spearman rho and n.
#' Gini rises with structure; the MFE z-score falls with it.
#'
#' DISPLAY RANGE. The z-scores are heavy-tailed (median 0.9, 99th percentile
#' 21, maximum 56), which squeezes the bulk into a sliver. The view is cropped
#' to the central `view_q` of each axis; the correlation uses all the genes.
#' The number outside the view is returned as attr(, "n_outside").
f5_plot_a <- function(a, compartment = "cytoplasm", view_q = 0.99,
                      base_size = 16) {
  ct <- a$cor; d <- a$data
  lo <- (1 - view_q) / 2
  xr <- stats::quantile(d$gini,  c(lo, 1 - lo))
  yr <- stats::quantile(d$mfe_z, c(lo, 1 - lo))
  inside  <- d$gini >= xr[1] & d$gini <= xr[2] &
             d$mfe_z >= yr[1] & d$mfe_z <= yr[2]
  outside <- sum(!inside)
  lab <- sprintf("ρ = %.2f, n = %s", ct$estimate,
                 format(ct$n, big.mark = ","))
  p <- ggplot(d[inside, ], aes(gini, mfe_z)) +
    geom_hex(bins = 35) +
    scale_fill_viridis_c(trans = "log10", name = "Genes") +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(title = "Probing and predicted structure", subtitle = lab,
         x = sprintf("%s (more structured →)",
                     f5_measure_labels(compartment)[["probing"]]),
         # the y label is rotated, so "down" on the page is "left" in its text
         y = sprintf("%s\n(← more structured)",
                     f5_measure_labels(compartment)[["prediction"]])) +
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
#' Heatmap of Spearman rho between Gini and each predicted metric (not
#' sign-changed, so negative = both say "more structured"), then C+G% and
#' length in a second block below: metrics down, regions across, one facet
#' per compartment. * marks BH q < 0.05 across the whole table.
f5_plot_b <- function(b, base_size = 15) {
  b$metric_f <- factor(F5_METRIC_LABELS[b$metric],
                       levels = rev(F5_METRIC_LABELS))
  b$block_f  <- factor(ifelse(b$metric %in% F5_SEQ_METRICS, "Sequence",
                              "Predicted structure"),
                       levels = c("Predicted structure", "Sequence"))
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
         title = "Probing–prediction correlation by region and metric",
         caption = "* BH q < 0.05") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid = element_blank(),
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_blank(),
          plot.title = element_text(face = "bold", size = base_size))
}


# =============================================================================
# Panel C
# =============================================================================
F5_MEASURE_COLOURS <- c(probing = "#D55E00", prediction = "#0072B2",
                        expression = "#009E73")

#' Absolute Spearman correlation of each measure with half-life and TE, per
#' region, with 95 % CIs folded to |rho| (an interval spanning zero starts at
#' 0). The sign is in the table (figure5_C_associations.csv): Gini is negative
#' with both responses, the MFE z-score and expression positive. Gene
#' expression is shown for the full mRNA only.
f5_plot_c <- function(cc, base_size = 15, compartment = "cytoplasm") {
  d <- cc
  d$measure_f <- factor(d$measure, levels = names(F5_MEASURE_COLOURS))
  d$region_f  <- factor(F5_REGION_LABELS[d$region],
                        levels = rev(F5_REGION_LABELS))
  d$resp_f    <- factor(F5_RESPONSE_LABELS[d$response],
                        levels = F5_RESPONSE_LABELS)
  dodge <- position_dodge(width = 0.6, preserve = "single")
  ggplot(d, aes(rho_abs, region_f, colour = measure_f, group = measure_f)) +
    # separators between regions (the default y grid runs through each
    # region's centre, which reads as belonging to one of its points)
    geom_hline(yintercept = seq_len(nlevels(d$region_f) - 1) + 0.5,
               colour = "grey80", linewidth = 0.4) +
    geom_errorbar(aes(xmin = conf.low_abs, xmax = conf.high_abs), width = 0.25,
                  position = dodge, orientation = "y", linewidth = 0.6) +
    geom_point(position = dodge, size = base_size / 15 * 3.2) +
    facet_wrap(~resp_f, nrow = 1) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_colour_manual(values = F5_MEASURE_COLOURS, name = NULL,
                        labels = f5_measure_labels(compartment)) +
    labs(x = "Absolute Spearman correlation with response", y = NULL,
         title = "Association of probing and predicted structure with half-life and translation efficiency") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "grey30", fill = NA,
                                      linewidth = 0.8),
          panel.spacing = unit(5, "mm"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size),
          legend.position = "bottom")
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
  cn <- cc |>
    dplyr::mutate(region = F5_REGION_LABELS[region],
                  measure = f5_measure_labels()[measure]) |>
    dplyr::select(region, response, measure, n) |>
    tidyr::pivot_wider(names_from = c(response, measure), values_from = n,
                       names_sep = " / ") |>
    dplyr::mutate(dplyr::across(-region, ~ ifelse(is.na(.x), "-", fmt(.x))))
  cw <- cc |>
    dplyr::mutate(region = F5_REGION_LABELS[region],
                  rho_ci = sprintf("%.3f [%.3f, %.3f]%s", rho, conf.low,
                                   conf.high, ifelse(q.value < 0.05, "", " (ns)"))) |>
    dplyr::select(region, response, measure, rho_ci, n) |>
    dplyr::mutate(n = fmt(n))
  lines <- c(
    "# Figure 5 run report", "",
    sprintf("Generated %s from `analysis/structure/figure5.R`.", Sys.Date()), "",
    "Figure subtitles are deliberately minimal; the numbers behind them are here, for the legends.", "",
    "No correlation is adjusted for length, C+G, coverage or expression: the figure is about what might confound.", "",
    "## Cohort", "",
    sprintf("- Transcripts in the analysed frame (human): **%s** (`build_dataset()`, UTR length filters applied)", fmt(nrow(df))),
    sprintf("- Genes after keeping one transcript per gene: **%s**", fmt(nrow(g))),
    "- icSHAPE Gini exists only for genes with enough coverage; every statistic below is on genes with Gini AND the predicted metric(s) in that region.", "",
    md(cohort), "",
    "## Panel A", "",
    sprintf("- Gini (cytoplasm, mRNA) vs robust MFE z-score: Spearman rho **%.3f** (95%% CI %.3f to %.3f), n **%s**",
            pa$cor$estimate, pa$cor$conf.low, pa$cor$conf.high, fmt(pa$cor$n)),
    sprintf("- Display cropped to the central 99%% of each axis; **%d** genes lie outside the view (all genes enter the correlation)", attr(pa$plot, "n_outside")),
    "- MFE z-score = (MFE - null median) / (1.4826 x MAD) over esl-shuffle nulls, i.e. a robust z-score; more negative = more structured (plotted as is)", "",
    "## Panel B", "",
    sprintf("- **%d** Spearman correlations (%d metrics x 4 regions x 2 compartments: 6 predicted metrics plus C+G%% and length), n **%s** to **%s** each; BH-adjusted over all %d; **%d** at q < 0.05",
            nrow(b), length(c(F5_METRICS, F5_SEQ_METRICS)), fmt(min(b$n)), fmt(max(b$n)), nrow(b), sum(b$q.value < 0.05)),
    "- Not sign-changed: negative rho = probing and prediction agree (high Gini, low/negative MFE)",
    "- C+G% and length are the region's own (`gc_content_<region>`, `length_<region>`); length is not log-transformed (Spearman is rank-based)", "",
    "## Panel C", "",
    "- Cytoplasm Gini and robust MFE z-score vs half-life and translation efficiency, Spearman, on the genes with Gini, MFE z-score and the response in that region (same genes for both measures)",
    "- Gene expression (Xpresso median RPKM) is added for the full mRNA only, on those of the same genes that also have an expression value, so its n is slightly smaller",
    "- The figure plots |rho| with each interval folded to |rho| (one spanning zero starts at 0); signed values are below and in the table",
    sprintf("- BH-adjusted over all %d tests of the table", nrow(cc)), "",
    "### n per region, response and measure", "", md(cn), "",
    "### Estimates, rho [95% CI]", "(ns) = BH q >= 0.05", "", md(as.data.frame(cw)), "",
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
  f5_save(pb,      "B", 260, 130)
  f5_save(pc,      "C", 330, 125)
  f5_save(full,    "full", 280, 250)

  write.csv(tab$a$data, file.path(OUTPUT_DIR, "tables", "figure5_A_points.csv"), row.names = FALSE)
  write.csv(tab$b, file.path(OUTPUT_DIR, "tables", "figure5_B_correlations.csv"), row.names = FALSE)
  write.csv(tab$c, file.path(OUTPUT_DIR, "tables", "figure5_C_associations.csv"), row.names = FALSE)
  f5_run_report(df, tab, pa, file.path(OUTPUT_DIR, "reports", "figure5_report.md"))
  message("Figure 5 written to ", OUTPUT_DIR)
}
