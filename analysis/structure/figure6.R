# =============================================================================
# Figure 6: shuffle z-score distributions vs native z-scores (6 transcripts)
# =============================================================================
# A: frequency histogram of each transcript's shuffle z-scores (shaded, overlaid)
#    with its native z marked by a vertical line in the same colour. Small grey
#    arrows mark every other mRNA z-score below F6_BACKGROUND_Z in the current
#    database, for context (no names, no shuffle curves).
#
# NOTHING IS RECOMPUTED HERE. Every z-score, MFE and shuffle comes from the
# mRNA_stability repo's special small run; this script only reads, checks and
# plots:
#   shuffle_distributions.tsv  one row per shuffled sequence (shuffle_z)
#   shuffle_summary.tsv        one row per transcript (native_z, stored_z, ...)
#   candidates.tsv             gene symbols, joined on transcript_id
# This run is not part of the database; candidates.tsv is used only for symbols.
# The background arrows are rnafold_zscore_mrna from the current database,
# build_dataset("human") with its default UTR length filters, and every named
# gene must pass them too (f6_load_background stops otherwise). The run also
# holds ZNF713 (no longer named; a background arrow) and GATD1 (5'UTR 25 nt,
# fails the filter; not drawn at all).
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
#  - The named transcripts were rerun only to keep their shuffles for the
#    histograms (the transcriptome-wide run did not keep them). native_z and
#    stored_z differ only by sampling variation between two shuffle sets
#    (TYRO3 most: -31.8 vs -50.2); the caption says so. The database is not
#    patched with native_z: replacing hand-picked transcripts would bias it.
#  - The background arrows are database z-scores (stored_z: same definition,
#    the earlier shuffles), so they are on stored_z's footing, not native_z's.
#    The named transcripts are drawn only once, at native_z, with their names.
#  - The x-axis is one continuous axis (no break, no rescaling), so the shuffle
#    curves are narrow next to the natives. That is the data, not a defect.
#
# Colours: Okabe-Ito (colour-blind-safe, 8 categories) with the yellow swapped
# for purple; the orange and sky blue (ZNF713's and GATD1's, before they were
# dropped) are unused, so each gene keeps its colour; see F6_COLOURS. The repo has no
# 8-category palette; REGION_COLOURS and the supergroup colours mean something
# else, so they are not reused.
#
# Usage:
#   source("analysis/structure/figure6.R")   # runs when sourced at top level
# Inputs (edit the paths below; data/raw/ is only ever read):
#   data/raw/figure6/shuffle_distributions.tsv
#   data/raw/figure6/shuffle_summary.tsv
#   ../mRNA_stability/runs/human_all/ultrafold/candidates.tsv
#   data/cache/human_dataset_v<CACHE_VERSION>.rds, via build_dataset()
# =============================================================================

if (!exists("build_dataset")) source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

F6_DIR            <- file.path(RAW_DIR, "figure6")
F6_DISTRIBUTIONS  <- file.path(F6_DIR, "shuffle_distributions.tsv")
F6_SUMMARY        <- file.path(F6_DIR, "shuffle_summary.tsv")
F6_CANDIDATES     <- Sys.getenv("F6_CANDIDATES", "../mRNA_stability/runs/human_all/ultrafold/candidates.tsv")

F6_GENES <- c("GEMIN6", "VHL", "ZNF101", "TYRO3", "TRUB2", "AP1M1")
F6_N_SHUFFLES        <- 100L
STORED_Z_TOLERANCE   <- 3

# Okabe-Ito with its yellow (#F0E442) swapped for a dark purple: yellow is
# unreadable on the grey panel and as label text. The orange (#E69F00) and sky
# blue (#56B4E9) were ZNF713's and GATD1's; they are left out so the other
# genes keep their colours.
F6_COLOURS <- stats::setNames(
  sub("#F0E442", "#6A3D9A", unname(grDevices::palette.colors(8, "Okabe-Ito")))[-(2:3)],
  F6_GENES)

F6_BINWIDTH     <- 0.25              # histogram bin width, in z units
F6_ARROW_LEVEL  <- 0.35              # named arrow height, as a share of the tallest bin
F6_BACKGROUND_Z <- -5                # mRNA z-scores below this are drawn small
F6_BACKGROUND_LEVEL <- 0.08          # their arrow height, as a share of the tallest bin
F6_BACKGROUND_ALPHA <- 0.25          # their opacity; overlaps add up, so dense runs read as density


#' Read the three inputs and join gene symbols. Keeps only `genes` (the run
#' has more). Stops on missing transcripts, a symbol that is not found, or
#' n_shuffles != F6_N_SHUFFLES.
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
    message("Not named, dropped from the shuffle tables: ",
            paste(extra, collapse = ", "))
  summ <- summ[summ$gene_name %in% genes, ]
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


#' Database mRNA z-scores below `threshold` (NA dropped), minus the transcripts
#' in `exclude` (the named ones, which are drawn at native_z). `df` is
#' build_dataset("human"). Stops if a named transcript is not in it.
f6_load_background <- function(df, threshold = F6_BACKGROUND_Z,
                               exclude = character()) {
  missing <- setdiff(exclude, df$transcript_id)
  if (length(missing))
    stop("Named transcripts not in the database (UTR length filters?): ",
            paste(missing, collapse = ", "), call. = FALSE)
  keep <- !is.na(df$rnafold_zscore_mrna) & df$rnafold_zscore_mrna < threshold &
          !df$transcript_id %in% exclude
  data.frame(transcript_id = df$transcript_id[keep],
             z = df$rnafold_zscore_mrna[keep])
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
#' Named arrows share one height (F6_ARROW_LEVEL); their labels read bottom to
#' top and are repelled sideways only, since some natives are < 0.2 z apart.
#' `background` (from f6_load_background) adds small grey unlabelled arrows.
f6_plot_a <- function(dist, summ, background = NULL, base_size = 14,
                      binwidth = F6_BINWIDTH) {
  n <- unique(summ$n_shuffles)
  summ <- summ[order(summ$native_z), ]
  # Tallest bin of any transcript sets the arrow heights (frequency scale).
  peak <- max(table(dist$gene_name, floor(dist$shuffle_z / binwidth))) / n
  nat  <- transform(summ, y = peak * F6_ARROW_LEVEL)
  bg   <- if (is.null(background)) data.frame(z = numeric()) else background

  ggplot() +
    geom_segment(data = bg,
                 aes(x = z, xend = z, y = peak * F6_BACKGROUND_LEVEL, yend = 0),
                 colour = "grey35", linewidth = 0.25, alpha = F6_BACKGROUND_ALPHA,
                 arrow = arrow(length = unit(0.8, "mm"), type = "closed")) +
    geom_histogram(data = dist,
                   aes(shuffle_z, after_stat(count / n), fill = gene_name),
                   binwidth = binwidth, position = "identity",
                   alpha = 0.3, colour = NA) +
    geom_segment(data = nat,
                 aes(x = native_z, xend = native_z, y = y, yend = 0,
                     colour = gene_name),
                 linewidth = 1, show.legend = FALSE,
                 arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
    geom_text_repel(data = nat, aes(native_z, y, label = gene_name,
                                    colour = gene_name),
                    angle = 90, hjust = 0, nudge_y = 0.12 * peak, direction = "x",
                    size = base_size * 0.28, segment.size = 0.3,
                    min.segment.length = 0, box.padding = 0.1, seed = 1,
                    show.legend = FALSE) +
    scale_fill_manual(values = F6_COLOURS, name = NULL) +
    scale_colour_manual(values = F6_COLOURS, name = NULL) +
    scale_x_continuous(breaks = scales::breaks_width(5),
                       expand = expansion(mult = c(0.06, 0.04))) +
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
f6_caption <- function(summ, background = NULL, binwidth = F6_BINWIDTH) {
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
      ", the minimum possible with ", n, " shuffles) for every transcript. "),
    if (!is.null(background) && nrow(background)) paste0(
      "Small grey arrows mark the other ", format(nrow(background), big.mark = ","),
      " transcripts with an mRNA z-score below ", F6_BACKGROUND_Z, " (",
      sprintf("%.1f", max(background$z)), " to ", sprintf("%.1f", min(background$z)),
      "), taken from the transcriptome-wide run (100 or 1,000 shuffles per ",
      "transcript). "),
    "The named transcripts were refolded with a fresh set of ", n, " shuffles ",
    "so that their shuffle distributions could be shown; their native z-scores ",
    "therefore differ from their transcriptome-wide values (by ",
    sprintf("%.1f", min(abs(summ$native_z - summ$stored_z))), " to ",
    sprintf("%.1f", max(abs(summ$native_z - summ$stored_z))), " z units) ",
    "through sampling variation between the two shuffle sets.",
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
  df <- build_dataset("human")
  bg <- f6_load_background(df, exclude = d$summary$transcript_id)
  f6_save(f6_plot_a(d$dist, d$summary, bg), "A", 180, 64)
  dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE, recursive = TRUE)
  writeLines(f6_caption(d$summary, bg), file.path(OUTPUT_DIR, "reports", "figure6_A_caption.txt"),
             useBytes = TRUE)
  message("Figure 6A written to ", file.path(OUTPUT_DIR, "plots"))
}
