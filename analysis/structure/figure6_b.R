# =============================================================================
# Figure 6B: ultrafolder repeat & structure architecture (MFE figure, genes stacked)
# =============================================================================
# One figure: the genes in F6B_GENES stacked under one main title, one shared
# x-axis label ("position"), each gene with its own panel pair (own axis in nt),
# ordered by native MFE z-score, lowest on top. Within a panel, top to bottom:
#   1. arcs        one arc per collapsed MFE stem (arcs_collapsed.tsv), spanning
#                  the midpoints of its two arms, width proportional to n_bp,
#                  coloured by class (Alu-Alu opposite / Alu-Alu same / other).
#   2. model       5'UTR thin, CDS thick, 3'UTR thin, to scale, with one arrow
#                  per Alu floating just above it: at the Alu's overall span
#                  (position and length), pointing right for sense and left for
#                  antisense; black if paired with an opposite Alu at the gene's
#                  cutoff, grey if not.
# There is no legend; the caption explains the colours, arrows and mark.
# Panel title: gene symbol. Transcript ID, region lengths and everything else (stem
# counts, widths, cutoff sensitivity) is in the caption files, not the figure.
#
# NOTHING IS RECOMPUTED HERE. Folding, Alu annotation, partner calls and
# cutoffs all come from the mRNA_stability repo; this script reads, checks and
# plots. Coordinates are 0-based half-open, 5'->3' on the spliced transcript.
# Every input row is selected by transcript_id (the input files hold all 8
# transcripts); the gene symbol and mRNA length come from model.tsv.
#
# CAVEATS (also in each caption file).
#  - The arcs are ONE predicted minimum-free-energy (MFE) structure (RNAfold).
#    The pairings were not tested for robustness (no base-pair probabilities or
#    suboptimal structures are used), so a missing or present arc is a property
#    of this one prediction.
#  - "Marked" means the Alu is paired with an oppositely oriented Alu in this
#    structure at the gene's cutoff, not merely oppositely oriented. Sequence
#    orientation and structure differ a lot (transcript_summary.tsv has both:
#    n_opposite_pairs_orientation vs n_opposite_pairs_structure); the caption
#    quotes the gene's own numbers.
#  - Marking is not selective: for most of these genes most Alus are marked, so
#    it does not single out a few special Alus. The caption says how many.
#  - Alu orientation is the RepeatMasker strand converted to sense/antisense
#    relative to the transcript. Coordinates are transcript-space. A split-exon
#    Alu is one row with its overall span (alu_partners.tsv).
#  - Per-Alu calls use the Alu's total opposite base pairs (all partners and
#    stems together), compared with >=.
#  - Arc width is scaled per gene panel (to that gene's drawn stems), so widths are
#    not comparable between genes; the caption gives each gene's n_bp range.
#
# Colours: region colours (5'UTR/CDS/3'UTR) come from REGION_COLOURS. The repo
# has no palette for arc class, so F6B_ARC_COLOURS defines it.
#
# Usage:
#   source("analysis/structure/figure6_b.R")   # runs when sourced at top level
# Inputs (data/raw/ is only ever read): data/raw/figure6/structure/*.tsv
# Outputs: plots/figure6_B_<gene>.{pdf,png,jpg}, reports/figure6_B_<gene>_caption.txt
# =============================================================================

if (!exists("OUTPUT_DIR")) source("R/config.R")
if (!exists("REGION_COLOURS")) {      # palettes.R needs the colour config and its sibling utils
  source("R/colour_config.R")
  for (f in list.files("R/utils", pattern = "\\.R$", full.names = TRUE)) source(f)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

F6B_DIR <- file.path(RAW_DIR, "figure6", "structure")

# One row per gene.
#   cutoff_type / cutoff_value: "Alu paired with an opposite Alu" is
#       "nbp" (>= N base pairs) or "frac" (>= N % of the Alu's length); >= is used.
#   min_bp: stems with fewer base pairs are not drawn (the caption says how many).
F6B_GENES <- data.frame(
  gene         = c("TRUB2",             "AP1M1",             "ZNF101"),
  transcript_id = c("ENST00000372890.6", "ENST00000291439.8", "ENST00000592502.2"),
  cutoff_type  = c("nbp",               "nbp",               "nbp"),
  cutoff_value = c(20,                  20,                  20),
  min_bp       = c(5,                   5,                   5),
  stringsAsFactors = FALSE)

F6B_WIDTH_MM    <- 240
F6B_HEIGHT_MM   <- 190     # three stacked panels (~58 mm each) + title + shared x label
F6B_ARROW_INSET <- 10       # nt trimmed from each end of an Alu arrow so abutting Alus stay distinct
F6B_ARC_HEIGHT  <- 0.25     # arc height / half-span (shape only; y axis is blank)

F6B_ARC_COLOURS <- c(`Alu-Alu opposite` = "#CC3311",
                     `Alu-Alu same`     = "#009988",
                     other              = "grey65")


f6b_read <- function(name, transcript_id) {
  f <- file.path(F6B_DIR, name)
  if (!file.exists(f)) stop("Input not found: ", f, call. = FALSE)
  d <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  d[d$transcript_id == transcript_id, , drop = FALSE]
}


#' Read every input for one row of F6B_GENES (filtered by transcript_id).
f6b_load <- function(cfg) {
  tx <- cfg$transcript_id
  model    <- f6b_read("model.tsv", tx)
  partners <- f6b_read("alu_partners.tsv", tx)
  bycut    <- f6b_read("alu_paired_by_cutoff.tsv", tx)
  summ     <- f6b_read("transcript_summary.tsv", tx)
  arcs     <- f6b_read("arcs_collapsed.tsv", tx)
  if (nrow(model) != 1L)
    stop("model.tsv: expected one row for ", cfg$gene, " (", tx, "), found ", nrow(model),
         call. = FALSE)

  alus <- bycut[bycut$cutoff_type == cfg$cutoff_type & bycut$cutoff_value == cfg$cutoff_value, ]
  if (nrow(alus) == 0L || nrow(alus) != nrow(partners))
    stop(cfg$gene, ": cutoff ", cfg$cutoff_type, " ", cfg$cutoff_value,
         " is not available for every Alu", call. = FALSE)
  alus$marked <- alus$paired_with_opposite %in% c(TRUE, "true", "TRUE")
  alus <- merge(alus, partners[, c("alu", "alu_len")], by = "alu")
  alus$orientation <- factor(alus$orientation, levels = c("sense", "antisense"))
  list(cfg = cfg, model = model, alus = alus[order(alus$tx_start), ], summ = summ, arcs = arcs)
}


#' One checks table for all genes; stops if any check fails.
f6b_checks <- function(ds) {
  rows <- lapply(ds, function(d) {
    m <- d$model; cfg <- d$cfg
    s <- d$summ[d$summ$cutoff_type == cfg$cutoff_type & d$summ$cutoff_value == cfg$cutoff_value, ]
    reg_sum <- m$utr5_len + m$cds_len + m$utr3_len
    data.frame(
      gene          = m$gene_name,
      transcript_id = m$transcript_id,
      mrna_len      = m$mrna_len,
      n_alu         = nrow(d$alus),
      cutoff        = paste0(cfg$cutoff_type, " >= ", cfg$cutoff_value),
      n_marked      = sum(d$alus$marked),
      mfe_match     = m$mfe_match %in% c(TRUE, "true", "TRUE"),
      regions_sum   = reg_sum,
      regions_ok    = reg_sum == m$mrna_len,
      symbol_ok     = m$gene_name == cfg$gene,
      # the summary table must agree with the per-Alu table at this cutoff
      summary_ok    = nrow(s) == 1L && s$n_alu == nrow(d$alus) &&
                      s$n_paired_opposite == sum(d$alus$marked))
  })
  chk <- do.call(rbind, rows)
  print(chk, row.names = FALSE)
  bad <- c("mfe_match", "regions_ok", "symbol_ok", "summary_ok")
  failed <- vapply(bad, function(k) !all(chk[[k]]), logical(1))
  if (any(failed)) stop("Check failed: ", paste(bad[failed], collapse = ", "), call. = FALSE)
  invisible(chk)
}


#' Half-ellipse path for one arc between x = a and x = b.
f6b_arc_path <- function(a, b, id, n = 40) {
  th <- seq(0, pi, length.out = n)
  data.frame(arc = id, x = (a + b) / 2 + (b - a) / 2 * cos(th),
             y = (b - a) / 2 * F6B_ARC_HEIGHT * sin(th))
}

f6b_arc_paths <- function(arcs) {
  do.call(rbind, lapply(seq_len(nrow(arcs)), function(k)
    data.frame(f6b_arc_path(arcs$a[k], arcs$b[k], k),
               arcs[k, setdiff(names(arcs), c("a", "b")), drop = FALSE], row.names = NULL)))
}


#' Arc layer data: stems with n_bp >= min_bp, positioned at arm midpoints.
#' Also returns what was hidden (count, n_bp range, Alu-Alu opposite stems).
f6b_mfe_arcs <- function(arcs, min_bp) {
  arcs$class <- ifelse(arcs$arc_class == "alu_opposite", "Alu-Alu opposite",
                ifelse(arcs$arc_class == "alu_same", "Alu-Alu same", "other"))
  arcs$class <- factor(arcs$class, levels = names(F6B_ARC_COLOURS))
  arcs$a <- (arcs$i_start + arcs$i_end) / 2
  arcs$b <- (arcs$j_start + arcs$j_end) / 2
  keep <- arcs[arcs$n_bp >= min_bp, ]
  hid  <- arcs[arcs$n_bp <  min_bp, ]
  keep <- keep[order(keep$class != "other", keep$n_bp), ]   # other at the back
  list(paths = f6b_arc_paths(keep[, c("a", "b", "class", "n_bp")]),
       n_shown = nrow(keep), n_hidden = nrow(hid),
       hidden_range = if (nrow(hid)) range(hid$n_bp) else c(NA_real_, NA_real_),
       hidden_opposite = sum(hid$class == "Alu-Alu opposite"))
}


#' Console table of what each figure hides.
f6b_hidden_table <- function(ds) {
  rows <- lapply(ds, function(d) {
    a <- f6b_mfe_arcs(d$arcs, d$cfg$min_bp)
    data.frame(gene = d$cfg$gene, min_bp = d$cfg$min_bp, stems = a$n_shown + a$n_hidden,
               shown = a$n_shown, hidden = a$n_hidden,
               hidden_bp = paste(a$hidden_range, collapse = "-"),
               hidden_alu_opposite = a$hidden_opposite)
  })
  print(do.call(rbind, rows), row.names = FALSE)
}


f6b_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey96"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.title.position = "plot",
          legend.position = "bottom")
}

f6b_blank_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank())
f6b_blank_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


#' Panel for one gene: arcs over the transcript model with Alu arrows. The gene
#' symbol (plus native z-score, if d$z is set) is its title; the main title and x-axis label are added by f6b_stack.
f6b_main <- function(d, base_size = 11) {
  cfg <- d$cfg
  L <- d$model$mrna_len
  xs <- scale_x_continuous(limits = c(0, L), expand = c(0, 0), labels = scales::comma,
                           breaks = function(lim) { b <- pretty(lim, n = 6); b[b <= lim[2]] })
  th <- f6b_theme(base_size)

  a <- f6b_mfe_arcs(d$arcs, cfg$min_bp)
  p_arc <- ggplot() +
    geom_path(data = a$paths, aes(x, y, group = arc, colour = class, linewidth = n_bp),
              alpha = 0.75, lineend = "round") +
    scale_linewidth_continuous(range = c(0.15, 2), guide = "none") +
    scale_colour_manual(values = F6B_ARC_COLOURS, name = "Pair class", drop = TRUE) +
    xs + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = if (is.null(d$z)) d$model$gene_name else
           parse(text = sprintf('"%s"~~"("*italic(z)~"="~"%s"*")"', d$model$gene_name,
                                sub("-", "\u2212", sprintf("%.1f", d$z)))), y = NULL) +
    th + f6b_blank_x + f6b_blank_y +
    theme(plot.title = element_text(face = "plain", size = base_size))

  cut_lab <- sprintf("Paired with an opposite Alu (≥ %s%s)", cfg$cutoff_value,
                     if (cfg$cutoff_type == "frac") "% of Alu length" else " bp")
  alus <- d$alus
  alus$mark <- factor(ifelse(alus$marked, cut_lab, "not marked"), levels = c(cut_lab, "not marked"))
  lo <- alus$tx_start + F6B_ARROW_INSET; hi <- alus$tx_end - F6B_ARROW_INSET
  alus$from <- ifelse(alus$orientation == "sense", lo, hi)   # sense points right
  alus$to   <- ifelse(alus$orientation == "sense", hi, lo)
  alus$y    <- 0.72
  m <- d$model
  reg <- data.frame(
    region = c("5′UTR", "CDS", "3′UTR"),
    xmin = c(0, m$utr5_len, m$utr5_len + m$cds_len),
    xmax = c(m$utr5_len, m$utr5_len + m$cds_len, m$mrna_len),
    h    = c(0.15, 0.45, 0.15),
    col  = unname(REGION_COLOURS[c("5utr", "cds", "3utr")]))
  reg$label <- sprintf("%s (%s nt)", reg$region, format(reg$xmax - reg$xmin, big.mark = ",", trim = TRUE))
  reg$lx <- c(reg$xmin[1], (reg$xmin[2:3] + reg$xmax[2:3]) / 2)
  reg$ly <- c(0.55, -0.55, -0.55)
  reg$vj <- c(0, 1, 1)
  p_mod <- ggplot(reg) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -h, ymax = h, fill = I(col)), colour = "grey20",
              linewidth = 0.25) +
    geom_text(aes(lx, ly, label = label, hjust = c(0, 0.5, 0.5), vjust = vj), size = base_size * 0.3) +
    geom_segment(data = alus, aes(x = from, xend = to, y = y, yend = y, colour = mark),
                 linewidth = 0.7, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
    scale_colour_manual(values = stats::setNames(c("black", "grey60"), levels(alus$mark)),
                        name = NULL, breaks = cut_lab) +   # unmarked is the default look: no key
    xs + scale_y_continuous(limits = c(-1.1, 1.0)) +
    labs(x = NULL, y = NULL) +
    th + f6b_blank_y + theme(panel.background = element_blank(), panel.border = element_blank(),
                             axis.line.x = element_line(colour = "grey20", linewidth = 0.3))

  # No legend: arc colours, arrow direction and the black/grey mark are explained
  # in the caption (f6b_caption), which keeps the figure clean.
  (p_arc / p_mod) + plot_layout(heights = c(3.6, 1.8)) & theme(legend.position = "none")
}


#' All genes stacked under one main title with one shared x-axis label.
f6b_stack <- function(ds, base_size = 11) {
  # lowest (most negative) native MFE z-score on top, from Figure 6A's table
  zs <- utils::read.delim(file.path(RAW_DIR, "figure6", "shuffle_summary.tsv"),
                          stringsAsFactors = FALSE)
  z <- stats::setNames(zs$native_z, zs$transcript_id)
  tid <- vapply(ds, function(d) d$cfg$transcript_id, "")
  if (anyNA(z[tid])) stop("native_z missing for: ", paste(names(ds)[is.na(z[tid])], collapse = ", "))
  ds <- ds[order(z[tid])]
  for (g in names(ds)) ds[[g]]$z <- unname(z[ds[[g]]$cfg$transcript_id])
  message("Panel order (native z): ", paste(sprintf("%s %.1f", names(ds), sort(z[tid])), collapse = ", "))
  panels <- lapply(ds, f6b_main, base_size = base_size)
  xlab <- wrap_elements(full = grid::textGrob("position", gp = grid::gpar(fontsize = base_size)))
  wrap_plots(c(panels, list(xlab)), ncol = 1,
             heights = c(rep(1, length(panels)), 0.06)) +
    plot_annotation(title = "Ultrafolder repeat & structure architecture",
                    theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
}


#' Journal-style caption, built from the gene's own data so the numbers cannot
#' drift from the figure.
f6b_caption <- function(d) {
  cfg <- d$cfg; m <- d$model; a <- d$alus; ss <- d$summ
  counts <- function(type) {
    x <- ss[ss$cutoff_type == type, ]
    v <- x$cutoff_value; n <- x$n_paired_opposite
    paste0(paste(n, collapse = ", "), " for cutoffs of ",
           paste(v[-length(v)], collapse = ", "), " and ", v[length(v)],
           if (type == "frac") "% of Alu length" else " base pairs")
  }
  s0 <- ss[ss$cutoff_type == cfg$cutoff_type & ss$cutoff_value == cfg$cutoff_value, ]
  cut_txt <- paste0("≥ ", cfg$cutoff_value,
                    if (cfg$cutoff_type == "frac") "% of the Alu's length" else " base pairs")
  n_ant <- sum(a$orientation == "antisense"); n_sen <- sum(a$orientation == "sense")
  n_mk <- sum(a$marked); n_al <- nrow(a)
  nr <- range(ss$n_paired_opposite)
  kind <- sub(".*:", "", a$alu)                    # repeat name, e.g. AluSx, FAM
  n_fam <- sum(kind == "FAM"); n_oth <- sum(!grepl("^Alu", kind) & kind != "FAM")
  extra <- c(if (n_fam) sprintf("%d %s a fossil Alu monomer, FAM", n_fam, if (n_fam == 1) "is" else "are"),
             if (n_oth) sprintf("%d %s not an Alu (%s)", n_oth, if (n_oth == 1) "is" else "are",
                     paste(unique(kind[!grepl("^Alu", kind) & kind != "FAM"]), collapse = ", ")))
  extra <- if (length(extra)) paste0("; ", paste(extra, collapse = "; ")) else ""

  ar <- f6b_mfe_arcs(d$arcs, cfg$min_bp)
  bp <- range(ar$paths$n_bp)
  hid_txt <- if (ar$n_hidden) paste0(
    "Stems with fewer than ", cfg$min_bp, " base pairs are omitted (", ar$n_shown, " of ",
    ar$n_shown + ar$n_hidden, " stems shown; the ", ar$n_hidden, " omitted stems have ",
    ar$hidden_range[1], "–", ar$hidden_range[2], " base pairs, ", ar$hidden_opposite,
    " of them between Alus of opposite orientation). ") else
    paste0("All ", ar$n_shown, " stems are shown. ")

  structure_txt <- if (!is.na(s0$n_opposite_pairs_structure[1]))
    paste0(" opposite Alu pairs, but only ", s0$n_opposite_pairs_structure[1],
           " pairs of Alus of opposite orientation are joined by ≥ ", cfg$cutoff_value,
           " base pairs in total (summed over stems).") else
    " opposite Alu pairs."
  sens_txt <- if (nr[1] == nr[2])
    paste0("The number of marked Alus is the same (", nr[1], ") at every cutoff tested. ") else
    paste0("The number of marked Alus ranges from ", nr[1], " to ", nr[2], " across cutoffs (",
           counts("nbp"), "; ", counts("frac"), "). ")
  share_txt <- if (n_mk / n_al >= 0.5)
    paste0("Because ", n_mk, " of the ", n_al, " elements are marked, the mark identifies most of ",
           "the repeats in this transcript, not a small subset. ") else ""

  paste0(
    "Figure 6B (", m$gene_name, "). Ultrafolder repeat & structure architecture of ", m$gene_name,
    " (", m$transcript_id, "), a ", format(m$mrna_len, big.mark = ","), "-nt spliced mRNA (5′UTR ",
    format(m$utr5_len, big.mark = ","), " nt, CDS ", format(m$cds_len, big.mark = ","), " nt, 3′UTR ", format(m$utr3_len, big.mark = ","),
    " nt). Top: base pairing in a single predicted minimum-free-energy (MFE) structure ",
    "(RNAfold). Each arc represents one stem (a run of consecutive base pairs) and joins the ",
    "centres of the stem’s two strands; arc width is proportional to the number of base pairs ",
    "in the stem (", bp[1], "–", bp[2], " in this figure). ", hid_txt,
    "Red arcs connect Alus of opposite orientation, teal arcs Alus of the same orientation, and ",
    "grey arcs all other pairings (within one Alu, between an Alu and non-Alu sequence, or ",
    "between non-Alu sequences). ",
    "Bottom: the transcript regions, drawn to scale, with the ", n_al, " annotated repeat ",
    "elements (", n_ant, " antisense, ", n_sen, " sense", extra, ") shown as arrows above them. ",
    "Each arrow spans the element’s position and length on the spliced transcript (each ",
    "trimmed by ", F6B_ARROW_INSET, " nt at both ends so that adjoining elements remain distinct) ",
    "and points right for sense and left for antisense orientation (the RepeatMasker strand ",
    "relative to the transcript; coordinates are transcript-space). Black arrows mark the ", n_mk,
    " elements whose total base pairing with Alus of opposite orientation in this structure, ",
    "summed over all of their stems and partners, is ", cut_txt, "; grey arrows are the others. ",
    "‘Marked’ therefore means paired in this structure, not merely oppositely oriented: ",
    "sequence orientation alone allows ", s0$n_opposite_pairs_orientation[1], structure_txt, " ",
    share_txt, sens_txt,
    "The structure is one prediction; the pairings were not tested for robustness.\n")
}


#' Combined caption for the whole set (shared definitions plus per-gene numbers),
#' built from the same data. Settings that differ between genes are not merged:
#' the caption then points to the per-gene captions.
f6b_caption_group <- function(ds) {
  cf  <- do.call(rbind, lapply(ds, function(d) d$cfg))
  same <- function(x) length(unique(x)) == 1L
  join <- function(x) if (length(x) <= 2) paste(x, collapse = " and ") else
    paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
  per <- function(f) vapply(ds, f, character(1))      # one string per gene, in order
  genes <- cf$gene

  cut_txt <- if (same(paste(cf$cutoff_type, cf$cutoff_value)))
    paste0("\u2265 ", cf$cutoff_value[1],
           if (cf$cutoff_type[1] == "frac") "% of the Alu's length" else " base pairs") else
    "the gene-specific cutoff given in the individual captions"
  min_txt <- if (same(cf$min_bp)) paste0("fewer than ", cf$min_bp[1], " base pairs") else
    "fewer base pairs than the gene-specific threshold given in the individual captions"

  info <- lapply(ds, function(d) {
    a <- f6b_mfe_arcs(d$arcs, d$cfg$min_bp); s0 <- d$summ[d$summ$cutoff_type == d$cfg$cutoff_type &
                                                           d$summ$cutoff_value == d$cfg$cutoff_value, ]
    kind <- sub(".*:", "", d$alus$alu)
    list(len = d$model$mrna_len, n_al = nrow(d$alus), n_mk = sum(d$alus$marked),
         shown = a$n_shown, stems = a$n_shown + a$n_hidden,
         n_or = s0$n_opposite_pairs_orientation[1], n_st = s0$n_opposite_pairs_structure[1],
         other = unique(kind[!grepl("^Alu", kind)]))
  })
  gtxt <- join(vapply(seq_along(ds), function(i) with(info[[i]], sprintf(
    "%s (%s nt; %d repeat elements, %d marked; %d of %d stems drawn)",
    genes[i], format(len, big.mark = ","), n_al, n_mk, shown, stems)), character(1)))
  oth <- vapply(seq_along(ds), function(i) if (length(info[[i]]$other))
    paste0(genes[i], ": ", paste(info[[i]]$other, collapse = ", ")) else NA_character_, character(1))
  oth <- oth[!is.na(oth)]
  all_major <- all(vapply(info, function(x) x$n_mk / x$n_al >= 0.5, logical(1)))
  have_st <- all(vapply(info, function(x) !is.na(x$n_st), logical(1)))

  paste0(
    "Figure 6B. Ultrafolder repeat & structure architecture of ", join(genes), ". ",
    "Transcripts are ordered by native MFE z-score, lowest at the top. ",
    "For each transcript, the top panel shows base pairing in a single predicted ",
    "minimum-free-energy (MFE) structure (RNAfold). Each arc represents one stem (a run of ",
    "consecutive base pairs) and joins the centres of the stem\u2019s two strands; arc width is ",
    "proportional to the number of base pairs in the stem, and stems with ", min_txt,
    " are omitted. Red arcs connect Alus of opposite orientation, teal arcs Alus of the same ",
    "orientation, and grey arcs all other pairings (within one Alu, between an Alu and non-Alu ",
    "sequence, or between non-Alu sequences). The bottom panel shows the transcript regions ",
    "(5\u2032UTR, CDS, 3\u2032UTR) to scale, with the annotated repeat elements as arrows above ",
    "them. Each arrow spans the element\u2019s position and length on the spliced transcript (each ",
    "trimmed by ", F6B_ARROW_INSET, " nt at both ends so that adjoining elements remain distinct) ",
    "and points right for sense and left for antisense orientation (the RepeatMasker strand ",
    "relative to the transcript; coordinates are transcript-space). Black arrows mark elements ",
    "whose total base pairing with Alus of opposite orientation in the structure, summed over all ",
    "of their stems and partners, is ", cut_txt, "; grey arrows are the others. ",
    "Transcripts: ", gtxt, ". ",
    if (length(oth)) paste0("The repeat elements include non-Alu members (", paste(oth, collapse = "; "),
                            "). ") else "",
    if (all_major) paste0("In every transcript most elements are marked, so the mark does not ",
                          "single out a small subset of repeats. ") else "",
    "\u2018Marked\u2019 means paired in the structure, not merely oppositely oriented",
    if (have_st) paste0(": sequence orientation alone allows ",
                        join(vapply(info, function(x) as.character(x$n_or), character(1))),
                        " opposite Alu pairs (", join(genes), ", respectively), of which only ",
                        join(vapply(info, function(x) as.character(x$n_st), character(1))),
                        " are joined by ", cut_txt, " in total") else "",
    ". The horizontal axis of each panel is scaled to that transcript\u2019s own length, and arc ",
    "width to that panel\u2019s own stems, so neither is comparable between transcripts. Each ",
    "structure is one prediction; the pairings were not tested for robustness. Per-gene ",
    "captions give each transcript\u2019s identifiers, region lengths, stem counts and cutoff ",
    "sensitivity.\n")
}


# Minimal save: PDF (vector), JPG (the repo's usual raster) and PNG.
f6b_save <- function(p, name, width, height) {
  base <- file.path(OUTPUT_DIR, "plots", paste0("figure6_", name))
  ggsave(paste0(base, ".pdf"), p, width = width, height = height, units = "mm",
         device = grDevices::cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = width, height = height, units = "mm", dpi = 300)
  ggsave(paste0(base, ".jpg"), p, width = width, height = height, units = "mm", dpi = 300)
}


if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE, recursive = TRUE)
  ds <- stats::setNames(lapply(seq_len(nrow(F6B_GENES)), function(i) f6b_load(F6B_GENES[i, ])),
                        F6B_GENES$gene)
  f6b_checks(ds)
  f6b_hidden_table(ds)
  f6b_save(f6b_stack(ds), "B", F6B_WIDTH_MM, F6B_HEIGHT_MM)
  for (d in ds) {
    writeLines(f6b_caption(d), file.path(OUTPUT_DIR, "reports",
                                         paste0("figure6_B_", d$cfg$gene, "_caption.txt")),
               useBytes = TRUE)
  }
  writeLines(f6b_caption_group(ds), file.path(OUTPUT_DIR, "reports", "figure6_B_caption.txt"),
             useBytes = TRUE)
  message("Figure 6B written to ", file.path(OUTPUT_DIR, "plots"))
}
