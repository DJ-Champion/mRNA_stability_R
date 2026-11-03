# =============================================================================
# Figure 6: shuffle z-score distributions vs native z-scores (8 transcripts)
# =============================================================================
# A: frequency histogram of each transcript's shuffle z-scores (shaded, overlaid)
#    with its native z marked by a vertical line in the same colour.
#
# NOTHING IS RECOMPUTED HERE. Every z-score, MFE and shuffle comes from the
# mRNA_stability repo's special small run; this script only reads, checks and
# plots:
#   shuffle_distributions.tsv  one row per shuffled sequence (shuffle_z)
#   shuffle_summary.tsv        one row per transcript (native_z, stored_z, ...)
#   candidates.tsv             gene symbols, joined on transcript_id
# This run is not part of the database; candidates.tsv is used only for symbols.
#
# WHAT THE FIGURE CAN AND CANNOT SAY.
#  - shuffle_z = (MFE - median) / (1.4826 * MAD) over that transcript's own
#    shuffles, so every curve is centred near 0 with scale near 1 by
#    construction. The overlap of the curves is NOT evidence on its own.
#  - n = 100 shuffles per transcript, so the empirical p-value
#    (shuffles <= native + 1) / (n + 1) cannot go below 1/101.
#  - native_z is recomputed from this run's 100 shuffles; stored_z is the
#    earlier transcriptome-wide value. They should be close; any transcript
#    where they differ by more than STORED_Z_TOLERANCE is flagged, not hidden.
#  - The x-axis is one continuous axis (no break, no rescaling), so the shuffle
#    curves are narrow next to the natives. That is the data, not a defect.
#
# Colours: Okabe-Ito (colour-blind-safe, 8 categories) with the yellow swapped
# for purple, see F6_COLOURS. The repo has no
# 8-category palette; REGION_COLOURS and the supergroup colours mean something
# else, so they are not reused.
#
# Usage:
#   source("analysis/structure/figure6.R")   # runs when sourced at top level
# Inputs (edit the paths below; data/raw/ is only ever read):
#   data/raw/figure6/shuffle_distributions.tsv
#   data/raw/figure6/shuffle_summary.tsv
#   ../mRNA_stability/runs/human_all/ultrafold/candidates.tsv
# =============================================================================

if (!exists("OUTPUT_DIR")) source("R/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

F6_DIR            <- file.path(RAW_DIR, "figure6")
F6_DISTRIBUTIONS  <- file.path(F6_DIR, "shuffle_distributions.tsv")
F6_SUMMARY        <- file.path(F6_DIR, "shuffle_summary.tsv")
F6_CANDIDATES     <- Sys.getenv("F6_CANDIDATES", "../mRNA_stability/runs/human_all/ultrafold/candidates.tsv")

F6_GENES <- c("GEMIN6", "ZNF713", "GATD1", "VHL", "ZNF101", "TYRO3", "TRUB2", "AP1M1")
F6_N_SHUFFLES        <- 100L
STORED_Z_TOLERANCE   <- 3

# Okabe-Ito with its yellow (#F0E442) swapped for a dark purple: yellow is
# unreadable on the grey panel and as label text.
F6_COLOURS <- stats::setNames(
  sub("#F0E442", "#6A3D9A", unname(grDevices::palette.colors(8, "Okabe-Ito"))),
  F6_GENES)

F6_BINWIDTH     <- 0.25              # histogram bin width, in z units
F6_ARROW_LEVELS <- c(0.15, 0.35, 0.55, 0.75)  # arrow heights, as a share of the tallest bin, cycled


#' Read the three inputs and join gene symbols. Stops on missing transcripts,
#' a symbol that is not found, or n_shuffles != F6_N_SHUFFLES.
f6_load <- function(distributions = F6_DISTRIBUTIONS, summary = F6_SUMMARY,
                    candidates = F6_CANDIDATES, genes = F6_GENES,
                    n_expected = F6_N_SHUFFLES) {
  for (f in c(distributions, summary, candidates))
    if (!file.exists(f)) stop("Input not found: ", f, call. = FALSE)
  rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  dist <- rd(distributions); summ <- rd(summary); cand <- rd(candidates)

  cand <- cand[!duplicated(cand$transcript_id), c("transcript_id", "gene_name")]
  summ <- merge(summ, cand, by = "transcript_id", all.x = TRUE, sort = FALSE)

  no_symbol <- summ$transcript_id[is.na(summ$gene_name)]
  if (length(no_symbol))
    stop("No gene symbol in candidates.tsv for: ",
         paste(no_symbol, collapse = ", "), call. = FALSE)
  missing <- setdiff(genes, summ$gene_name)
  if (length(missing))
    stop("Transcripts missing from shuffle_summary.tsv: ",
         paste(missing, collapse = ", "), call. = FALSE)
  extra <- setdiff(summ$gene_name, genes)
  if (length(extra))
    stop("Unexpected transcripts in shuffle_summary.tsv: ",
         paste(extra, collapse = ", "), call. = FALSE)
  bad_n <- summ$gene_name[summ$n_shuffles != n_expected]
  if (length(bad_n))
    stop("n_shuffles is not ", n_expected, " for: ",
         paste(bad_n, collapse = ", "), call. = FALSE)

  dist <- merge(dist, summ[, c("transcript_id", "gene_name")],
                by = "transcript_id", sort = FALSE)
  rows <- table(dist$gene_name)
  bad_rows <- names(rows)[rows != n_expected]
  if (length(bad_rows) || length(rows) != length(genes))
    stop("shuffle_distributions.tsv does not have ", n_expected,
         " rows for every transcript (", paste(names(rows), rows, sep = "=",
                                               collapse = ", "), ")", call. = FALSE)

  summ$gene_name <- factor(summ$gene_name, levels = genes)
  dist$gene_name <- factor(dist$gene_name, levels = genes)
  list(dist = dist, summary = summ[order(summ$gene_name), ])
}


#' Console checks table; flags native_z vs stored_z differences.
f6_checks <- function(summ, tolerance = STORED_Z_TOLERANCE) {
  chk <- data.frame(
    gene           = as.character(summ$gene_name),
    native_z       = round(summ$native_z, 2),
    stored_z       = round(summ$stored_z, 2),
    `native-stored` = round(summ$native_z - summ$stored_z, 2),
    n_shuffles     = summ$n_shuffles,
    skew           = round(summ$shuffle_skew, 2),
    kurtosis       = round(summ$shuffle_excess_kurtosis, 2),
    empirical_p    = signif(summ$empirical_p, 3),
    check.names = FALSE)
  print(chk, row.names = FALSE)
  off <- !is.na(chk$`native-stored`) & abs(chk$`native-stored`) > tolerance
  if (any(off))
    warning("native_z and stored_z differ by more than ", tolerance,
            " z units for: ", paste(chk$gene[off], collapse = ", "), call. = FALSE)
  if (any(is.na(summ$stored_z)))
    warning("stored_z is missing for: ",
            paste(summ$gene_name[is.na(summ$stored_z)], collapse = ", "), call. = FALSE)
  invisible(chk)
}


#' Panel A: shuffle z histograms (frequency = fraction of the transcript's
#' shuffles in each bin, i.e. count / n_shuffles; one per
#' transcript, shaded and overlaid) with each native z marked by a short arrow pointing at the x-axis.
#' Arrows are staggered over F6_ARROW_LEVELS (in rank order of native z) so the
#' gene labels above them do not collide.
f6_plot_a <- function(dist, summ, base_size = 14, binwidth = F6_BINWIDTH) {
  n <- unique(summ$n_shuffles)
  summ <- summ[order(summ$native_z), ]
  # Tallest bin of any transcript sets the arrow heights (frequency scale).
  peak <- max(table(dist$gene_name, floor(dist$shuffle_z / binwidth))) / n
  nat  <- transform(summ, y = peak * rep_len(F6_ARROW_LEVELS, nrow(summ)))

  ggplot() +
    geom_histogram(data = dist,
                   aes(shuffle_z, after_stat(count / n), fill = gene_name),
                   binwidth = binwidth, position = "identity",
                   alpha = 0.3, colour = NA) +
    geom_segment(data = nat,
                 aes(x = native_z, xend = native_z, y = y, yend = 0,
                     colour = gene_name),
                 linewidth = 1, show.legend = FALSE,
                 arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
    geom_text(data = nat, aes(native_z, y, label = gene_name, colour = gene_name),
              vjust = -0.4, size = base_size * 0.28, show.legend = FALSE) +
    scale_fill_manual(values = F6_COLOURS, name = NULL) +
    scale_colour_manual(values = F6_COLOURS, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "Native z-scores vs per-transcript shuffle distributions",
      x = "z-score",
      y = "frequency") +
    theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey93"),
          plot.title = element_text(face = "bold", size = base_size),
          plot.title.position = "plot",
          legend.key.height = unit(4, "mm"))
}


#' Journal-style legend text for panel A, built from the data so the numbers
#' cannot drift from the figure. Makes only claims the data support: the
#' empirical-p sentence is added only if every transcript sits at the floor.
f6_caption <- function(summ, binwidth = F6_BINWIDTH) {
  n    <- unique(summ$n_shuffles)
  nat  <- range(summ$native_z)
  floor_p <- 1 / (n + 1)
  at_floor <- all(abs(summ$empirical_p - floor_p) < 1e-3 * floor_p)
  paste0(
    "Figure 6A. Native predicted folding z-scores against shuffled-sequence ",
    "distributions for ", nrow(summ), " transcripts. Histograms show the ",
    "z-scores of n = ", n, " shuffled sequences per transcript (bin width ",
    binwidth, "); frequency is the fraction of that transcript's shuffles in ",
    "each bin. For each transcript the shuffle z-score is (MFE \u2212 median) / ",
    "(1.4826 \u00d7 MAD), with the median and median absolute deviation (MAD) ",
    "taken over its own shuffles, so every histogram is centred near 0 with a ",
    "spread near 1 by construction. Arrows mark the native z-score, computed ",
    "in the same way from the native MFE (", sprintf("%.1f", nat[2]), " to ",
    sprintf("%.1f", nat[1]), "). ",
    if (at_floor) paste0(
      "The empirical p-value was ", signif(floor_p, 2), " (1/", n + 1L,
      ", the minimum possible with ", n, " shuffles) for every transcript."),
    "\n")
}


f6_save <- function(p, name, width, height) {
  base <- file.path(OUTPUT_DIR, "plots", paste0("figure6_", name))
  ggsave(paste0(base, ".jpg"), p, width = width, height = height,
         units = "mm", dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height,
         units = "mm", device = grDevices::cairo_pdf)
}


if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
  d <- f6_load()
  f6_checks(d$summary)
  f6_save(f6_plot_a(d$dist, d$summary), "A", 180, 64)
  dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE, recursive = TRUE)
  writeLines(f6_caption(d$summary), file.path(OUTPUT_DIR, "reports", "figure6_A_caption.txt"),
             useBytes = TRUE)
  message("Figure 6A written to ", file.path(OUTPUT_DIR, "plots"))
}
