# =============================================================================
# Figure 6: shuffle z-score distributions vs database z-scores (6 transcripts)
# =============================================================================
# A: frequency histogram of each transcript's shuffle z-scores (shaded, overlaid)
#    with its database mRNA z-score marked by an arrow in the same colour. Small
#    grey arrows mark every other mRNA z-score below F6_BACKGROUND_Z in the
#    current database, for context (no names, no shuffle curves).
#
# The named transcripts are the six lowest rnafold_zscore_mrna values in the
# current database, build_dataset("human") with its default UTR length filters
# (f6_caption checks this before saying so).
#
# NOTHING IS RECOMPUTED HERE. This script only reads, checks and plots:
#   shuffle_distributions.tsv  one row per shuffled sequence (shuffle_z); the
#                              histograms only
#   shuffle_summary.tsv        one row per transcript (n_shuffles, stored_z, ...)
#   structure/model.tsv        gene symbols, joined on transcript_id
#   rnafold_zscore_mrna, rnafold_pval_mrna from build_dataset("human"); every
#                              arrow, named or grey, and every number in the
#                              caption apart from the histogram n
# f6_add_database stops if stored_z does not match the database value.
#
# WHAT THE FIGURE CAN AND CANNOT SAY.
#  - shuffle_z = (MFE - median) / (1.4826 * MAD) over that transcript's own
#    shuffles, so every curve is centred near 0 with scale near 1 by
#    construction. The overlap of the curves is NOT evidence on its own.
#  - The histograms are a simplified picture of the shuffle null (100 shuffles
#    per transcript). The arrows are database z-scores (same definition, 100 or
#    1,000 shuffles per transcript); the histogram run's own native z is not
#    used anywhere in the figure or caption.
#  - The database empirical p-value (shuffles <= native + 1) / (n + 1) cannot go
#    below 1/(n + 1); the caption says so only if every named transcript is there.
#  - The x-axis is one continuous axis (no break, no rescaling), so the shuffle
#    curves are narrow next to the natives. That is the data, not a defect.
#
# Colours: Okabe-Ito (colour-blind-safe, 8 categories) with the yellow swapped
# for purple and the orange and sky blue left out, so the figure keeps the
# colours of its earlier versions; see F6_COLOURS. The repo has no 8-category
# palette; REGION_COLOURS and the supergroup colours mean something else, so
# they are not reused.
#
# Usage:
#   source("analysis/structure/figure6.R")   # runs when sourced at top level
# Inputs (edit the paths below; data/raw/ is only ever read):
#   data/raw/ultrafolders/shuffle_distributions.tsv
#   data/raw/ultrafolders/shuffle_summary.tsv
#   data/raw/ultrafolders/structure/model.tsv
#   data/cache/human_dataset_v<CACHE_VERSION>.rds, via build_dataset()
# =============================================================================

if (!exists("build_dataset")) source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

F6_DIR            <- file.path(RAW_DIR, "ultrafolders")
F6_DISTRIBUTIONS  <- file.path(F6_DIR, "shuffle_distributions.tsv")
F6_SUMMARY        <- file.path(F6_DIR, "shuffle_summary.tsv")
F6_SYMBOLS        <- file.path(F6_DIR, "structure", "model.tsv")

F6_GENES <- c("AP1M1", "ORAI2", "TYRO3", "LRRC57", "CYP20A1", "ZBTB25")
F6_N_SHUFFLES        <- 100L
F6_DB_SHUFFLE_TIERS  <- c(100L, 1000L)   # shuffles per transcript in the database
STORED_Z_TOLERANCE   <- 1e-3             # stored_z vs database z (same number, rounded)

# Okabe-Ito with its yellow (#F0E442) swapped for a dark purple: yellow is
# unreadable on the grey panel and as label text. The orange (#E69F00) and sky
# blue (#56B4E9) are left out so the figure keeps its earlier colours.
F6_COLOURS <- stats::setNames(
  sub("#F0E442", "#6A3D9A", unname(grDevices::palette.colors(8, "Okabe-Ito")))[-(2:3)],
  F6_GENES)

F6_BINWIDTH     <- 0.25              # histogram bin width, in z units
F6_ARROW_LEVEL  <- 0.35              # named arrow height, as a share of the tallest bin
F6_BACKGROUND_Z <- -5                # mRNA z-scores below this are drawn small
F6_BACKGROUND_LEVEL <- 0.08          # their arrow height, as a share of the tallest bin
F6_BACKGROUND_ALPHA <- 0.25          # their opacity; overlaps add up, so dense runs read as density


#' Read the three inputs and join gene symbols. Keeps only `genes` (the run
#' could have more). Stops on missing transcripts, a symbol that is not found,
#' or n_shuffles != F6_N_SHUFFLES.
f6_load <- function(distributions = F6_DISTRIBUTIONS, summary = F6_SUMMARY,
                    symbols = F6_SYMBOLS, genes = F6_GENES,
                    n_expected = F6_N_SHUFFLES) {
  for (f in c(distributions, summary, symbols))
    if (!file.exists(f)) stop("Input not found: ", f, call. = FALSE)
  rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  dist <- rd(distributions); summ <- rd(summary); sym <- rd(symbols)

  sym <- sym[!duplicated(sym$transcript_id), c("transcript_id", "gene_name")]
  summ <- merge(summ, sym, by = "transcript_id", all.x = TRUE, sort = FALSE)

  no_symbol <- summ$transcript_id[is.na(summ$gene_name)]
  if (length(no_symbol))
    stop("No gene symbol in ", basename(symbols), " for: ",
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


#' Add the database mRNA z-score (z), empirical p-value (p) and rank among all
#' database z-scores to `summ`. `df` is build_dataset("human"). Stops if a
#' transcript is not in it or its stored_z is not the database value.
f6_add_database <- function(summ, df, tolerance = STORED_Z_TOLERANCE) {
  m <- match(summ$transcript_id, df$transcript_id)
  if (anyNA(m))
    stop("Named transcripts not in the database (UTR length filters?): ",
         paste(summ$gene_name[is.na(m)], collapse = ", "), call. = FALSE)
  summ$z <- df$rnafold_zscore_mrna[m]
  summ$p <- df$rnafold_pval_mrna[m]
  summ$rank <- rank(df$rnafold_zscore_mrna, na.last = "keep")[m]
  off <- is.na(summ$z) | abs(summ$z - summ$stored_z) > tolerance
  if (any(off))
    stop("stored_z is not the database z-score for: ",
         paste(summ$gene_name[off], collapse = ", "), call. = FALSE)
  summ
}


#' Database mRNA z-scores below `threshold` (NA dropped), minus the transcripts
#' in `exclude` (the named ones, which are drawn in colour). `df` is
#' build_dataset("human").
f6_load_background <- function(df, threshold = F6_BACKGROUND_Z,
                               exclude = character()) {
  keep <- !is.na(df$rnafold_zscore_mrna) & df$rnafold_zscore_mrna < threshold &
          !df$transcript_id %in% exclude
  data.frame(transcript_id = df$transcript_id[keep],
             z = df$rnafold_zscore_mrna[keep])
}


#' Database p-value floor tier for each p: the shuffle count n in
#' F6_DB_SHUFFLE_TIERS with p == 1/(n + 1) (to the stored rounding), else NA.
f6_floor_tier <- function(p, tiers = F6_DB_SHUFFLE_TIERS) {
  vapply(p, function(x) {
    hit <- tiers[abs(x * (tiers + 1) - 1) < 0.02]
    if (length(hit) == 1L) hit else NA_integer_
  }, integer(1))
}


#' Console checks table (database values, plus the shape of each histogram).
f6_checks <- function(summ) {
  chk <- data.frame(
    gene           = as.character(summ$gene_name),
    database_z     = round(summ$z, 2),
    database_rank  = summ$rank,
    database_p     = signif(summ$p, 3),
    p_floor_tier   = f6_floor_tier(summ$p),
    hist_n         = summ$n_shuffles,
    hist_skew      = round(summ$shuffle_skew, 2),
    hist_kurtosis  = round(summ$shuffle_excess_kurtosis, 2))
  print(chk, row.names = FALSE)
  invisible(chk)
}


#' Panel A: shuffle z histograms (frequency = fraction of the transcript's
#' shuffles in each bin, i.e. count / n_shuffles; one per
#' transcript, shaded and overlaid) with each database z (summ$z, from
#' f6_add_database) marked by a short arrow pointing at the x-axis.
#' Named arrows share one height (F6_ARROW_LEVEL); their labels read bottom to
#' top and are repelled sideways only, in case two z-scores are close.
#' `background` (from f6_load_background) adds small grey unlabelled arrows.
f6_plot_a <- function(dist, summ, background = NULL, base_size = 14,
                      binwidth = F6_BINWIDTH) {
  n <- unique(summ$n_shuffles)
  summ <- summ[order(summ$z), ]
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
                 aes(x = z, xend = z, y = y, yend = 0,
                     colour = gene_name),
                 linewidth = 1, show.legend = FALSE,
                 arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
    geom_text_repel(data = nat, aes(z, y, label = gene_name,
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
#' cannot drift from the figure. `summ` comes from f6_add_database and
#' `n_database` is the number of transcripts with a database mRNA z-score. Makes
#' only claims the data support: "the k lowest" only if no other transcript is
#' lower, and the empirical-p sentence only if every transcript sits at a floor.
f6_caption <- function(summ, n_database, background = NULL, binwidth = F6_BINWIDTH) {
  n    <- unique(summ$n_shuffles)
  zr   <- range(summ$z)
  lowest <- all(summ$rank <= nrow(summ))
  tier <- f6_floor_tier(summ$p)
  p_txt <- if (!anyNA(tier)) {
    parts <- vapply(sort(unique(tier)), function(t) {
      g <- as.character(summ$gene_name[tier == t])
      paste0(signif(1 / (t + 1), 2), " with ", format(t, big.mark = ","), " shuffles",
             if (length(unique(tier)) > 1 && length(g) < nrow(summ) / 2)
               paste0(" for ", paste(g, collapse = ", ")) else "")
    }, character(1))
    paste0("For every named transcript no shuffle folded as stably as the native ",
           "sequence, so its empirical p-value is the minimum possible, 1/(n + 1) for ",
           "n shuffles (", paste(parts, collapse = "; "), "). ")
  }
  paste0(
    "Figure 6A. Native predicted folding z-scores against shuffled-sequence ",
    "distributions for ", nrow(summ), " transcripts. Histograms show the ",
    "z-scores of n = ", n, " shuffled sequences per transcript (bin width ",
    binwidth, "); frequency is the fraction of that transcript's shuffles in ",
    "each bin. Each z-score is (MFE \u2212 median) / (1.4826 \u00d7 MAD), with ",
    "the median and median absolute deviation (MAD) taken over the transcript's ",
    "shuffles, so every histogram is centred near 0 with a spread near 1 by ",
    "construction. Coloured arrows mark each transcript's mRNA z-score (",
    sprintf("%.1f", zr[2]), " to ", sprintf("%.1f", zr[1]), ")",
    if (lowest) paste0(", the ", nrow(summ), " lowest of the ",
                       format(n_database, big.mark = ","), " transcripts analysed"),
    ". ", p_txt,
    if (!is.null(background) && nrow(background)) paste0(
      "Small grey arrows mark the other ", format(nrow(background), big.mark = ","),
      " transcripts with an mRNA z-score below ", F6_BACKGROUND_Z, " (",
      sprintf("%.1f", max(background$z)), " to ", sprintf("%.1f", min(background$z)),
      "). "),
    "All mRNA z-scores are from the transcriptome-wide analysis (",
    paste(format(F6_DB_SHUFFLE_TIERS, big.mark = ",", trim = TRUE), collapse = " or "),
    " shuffles per transcript).",
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
  df <- build_dataset("human")
  d$summary <- f6_add_database(d$summary, df)
  f6_checks(d$summary)
  bg <- f6_load_background(df, exclude = d$summary$transcript_id)
  f6_save(f6_plot_a(d$dist, d$summary, bg), "A", 180, 64)
  dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE, recursive = TRUE)
  writeLines(f6_caption(d$summary, sum(!is.na(df$rnafold_zscore_mrna)), bg), file.path(OUTPUT_DIR, "reports", "figure6_A_caption.txt"),
             useBytes = TRUE)
  message("Figure 6A written to ", file.path(OUTPUT_DIR, "plots"))
}
