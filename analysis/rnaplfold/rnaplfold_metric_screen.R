# =============================================================================
# RNAplfold metric screen: does any summary beat the plain start-window mean?
# =============================================================================
# Every rnaplfold_* column from rnaplfold_combined.tsv against half-life and
# translation efficiency, as raw Spearman rho with a 95% CI and as a PARTIAL rho
# adjusted for 5' UTR / CDS length and 5' UTR / start-window GC. Meant for a
# pilot run on a few hundred genes, to decide whether a full run is worth it.
#
# The RNAplfold columns are joined onto the cached dataset IN MEMORY from the
# loader, replacing whatever the cache holds, so a pilot file never touches the
# cache. Rows are restricted to genes that carry a value, so n is the pilot's.
#
# Two views of every column, because they answer different questions:
#   unadjusted  rho: what the metric carries, GC and length included. This is
#               what a banded correlation plot shows.
#   adjusted    partial rho, after regressing both sides on 5'/CDS/3' UTR length
#               and 5'UTR/CDS/3'UTR/start-window GC (ranks): what is left once
#               composition and length are removed.
# The benchmark is the plain start mean u<l> (the feature already in the model).
#
# Usage: Rscript analysis/rnaplfold/rnaplfold_metric_screen.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
source("R/load_all.R")

OUT <- file.path(OUTPUT_DIR, "rnaplfold_exploration")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

responses <- c("halflife", "translation_efficiency")
adjust    <- c("length_5utr", "length_cds", "length_3utr", "gc_content_5utr",
               "gc_content_cds", "gc_content_3utr", "gc_content_start")

rp <- load_rnaplfold("human", all = TRUE)
df <- build_dataset("human") |> drop_excluded(verbose = FALSE)
df <- df |> select(-starts_with("rnaplfold_")) |>
  inner_join(rp, by = "transcript_id")
acc <- setdiff(names(rp), "transcript_id")
cat("genes in pilot after cohort filter:", nrow(df), "\n")

rk <- function(x) rank(x, na.last = "keep")
screen_one <- function(v, y) {
  d <- df[c(v, y, adjust)]; d <- d[complete.cases(d), ]
  n <- nrow(d); d[] <- lapply(d, rank)
  raw <- cor(d[[v]], d[[y]])
  res <- function(z) resid(lm(reformulate(adjust, z), d))
  part <- cor(res(v), res(y))
  se <- 1.06 / sqrt(n - 3)
  data.frame(feature = v, response = y, n = n, rho = raw,
             lo = tanh(atanh(raw) - 1.96 * se), hi = tanh(atanh(raw) + 1.96 * se),
             partial = part,
             gc_start = cor(d[[v]], d$gc_content_start))
}
tab <- do.call(rbind, lapply(responses, function(y)
  do.call(rbind, lapply(acc, screen_one, y = y))))

# Parse the column name into stat / u length / region for grouping.
p <- regmatches(tab$feature, regexec(
  "^rnaplfold_(?:(.*)_)?unpaired_u([0-9]+)_([a-z0-9]+)$", tab$feature, perl = TRUE))
tab$stat <- vapply(p, function(m) if (length(m) == 4 && nzchar(m[2])) m[2] else "mean", "")
tab$u    <- vapply(p, function(m) if (length(m) == 4) as.integer(m[3]) else NA_integer_, 1L)
tab$region <- vapply(p, function(m) if (length(m) == 4) m[4] else "start", "")
tab$stat[tab$feature == "rnaplfold_aug_unpaired_start"] <- "aug"

# Family names as RNAPLFOLD_FEATURES.md calls them. "Mean log" is the mean of
# log10 stretch probability (mean opening energy up to a constant): closer to
# 0 = more accessible. "Flank" windows are f nt each side of the AUG; they are
# not RNAplfold's -W.
family_label <- function(stat, region) {
  base <- c(mean = "Mean", meanlog = "Mean log", min = "Min (bottleneck)",
            max = "Max (best landing pad)", median = "Median",
            openfrac = "Open fraction", up_meanlog = "Upstream mean log",
            down_meanlog = "Downstream mean log", aug = "AUG")
  ifelse(region == "start",
    ifelse(grepl("^flank", stat),
           paste0("Flank ", sub("^flank([0-9]+)_mean(log)?$", "\\1 nt, mean\\2", stat) |> sub(pattern = "meanlog", replacement = "mean log")),
           paste0("Start window: ", base[stat])),
    ifelse(stat == "cap50_meanlog", "Cap-proximal (nt 1-50) mean log",
           paste0(c(`5utr` = "5' UTR", cds = "CDS", `3utr` = "3' UTR",
                    mrna = "mRNA")[region], " mean log")))
}
tab$family <- family_label(tab$stat, tab$region)

# The short list RNAPLFOLD_FEATURES.md suggests deciding on up front, so the
# headline read is not a search across 65 correlated columns. CDS is not on
# it; its mean log is shown with the rest in the full table.
STARTING_SET <- c("rnaplfold_meanlog_unpaired_u4_start",
                  "rnaplfold_min_unpaired_u8_start",
                  "rnaplfold_max_unpaired_u8_start",
                  "rnaplfold_up_meanlog_unpaired_u4_start",
                  "rnaplfold_meanlog_unpaired_u4_5utr",
                  "rnaplfold_meanlog_unpaired_u4_3utr",
                  "rnaplfold_cap50_meanlog_unpaired_u4_mrna")

mde <- 2.8 * 1.06 / sqrt(min(tab$n) - 3)   # ~80% power, two-sided 0.05
tab$q <- p.adjust(2 * pnorm(-abs(atanh(tab$rho)) * sqrt(tab$n - 3) / 1.06), "BH")


bench <- tab |> filter(stat == "mean", region == "start") |>
  select(response, u, bench_rho = rho, bench_partial = partial)
tab <- tab |> left_join(bench, by = c("response", "u"))
write.csv(tab, file.path(OUT, "metric_screen.csv"), row.names = FALSE)

show <- function(y, by) {  # feature = canonical column name
  cat("\n=== ", y, ", top 12 by |", by, "|\n", sep = "")
  tab |> filter(response == y) |> arrange(desc(abs(.data[[by]]))) |>
    transmute(feature = sub("^rnaplfold_", "", feature), n, rho = round(rho, 3),
              partial = round(partial, 3), q = signif(q, 2),
              vs_gc = round(gc_start, 2)) |>
    head(12) |> as.data.frame() |> print(row.names = FALSE)
}
for (y in responses) for (by in c("rho", "partial")) show(y, by)
cat("\nbenchmark (plain start mean u<l>):\n")
print(as.data.frame(bench |> mutate(across(c(bench_rho, bench_partial), ~ round(., 3)))),
      row.names = FALSE)

# --- the pre-specified starting set (q over these 14 tests only) --------------
ss <- tab |> filter(feature %in% STARTING_SET) |>
  mutate(q_set = p.adjust(2 * pnorm(-abs(atanh(rho)) * sqrt(n - 3) / 1.06), "BH"))
cat("\nStarting set from RNAPLFOLD_FEATURES.md:\n")
ss |> transmute(response, feature = sub("^rnaplfold_", "", feature), n,
                rho = round(rho, 3), partial = round(partial, 3),
                q_set = signif(q_set, 2)) |>
  as.data.frame() |> print(row.names = FALSE)
write.csv(ss, file.path(OUT, "metric_screen_starting_set.csv"), row.names = FALSE)

# --- figure: unadjusted rho (95% CI), by family x response ---------------------
resp_lab <- c(halflife = "Half-life", translation_efficiency = "Translation efficiency")
fam <- tab |> filter(!is.na(u)) |>
  mutate(u = factor(paste0("u", u), paste0("u", sort(unique(u)))),
         response = resp_lab[response])
fam$family <- factor(fam$family, rev(unique(fam$family[order(-abs(fam$rho))])))
pl <- ggplot(fam, aes(rho, family, colour = u)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7),
                  size = 0.25) +
  facet_wrap(~ response) +
  labs(x = "Spearman rho (95% CI)", y = NULL, colour = "stretch",
       title = "RNAplfold metrics vs response (unadjusted)",
       subtitle = sprintf("n = %d..%d genes.\nPositive rho = more accessible goes with a higher response\n(mean log: closer to 0 = more accessible)", min(tab$n), max(tab$n))) +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "metric_screen.jpg"), pl, width = 230, height = 110,
       units = "mm", dpi = 200)

# --- figure: mean log in the four regions (start window, 5' UTR, CDS, 3' UTR) --
reg_lab <- c(start = "Start codon", `5utr` = "5' UTR", cds = "CDS", `3utr` = "3' UTR")
reg <- tab |> filter(stat == "meanlog", region %in% names(reg_lab)) |>
  mutate(region = factor(reg_lab[region], rev(reg_lab)),
         u = factor(paste0("u", u), paste0("u", sort(unique(u)))),
         response = resp_lab[response])
pr <- ggplot(reg, aes(rho, region, colour = u)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6),
                  size = 0.35) +
  facet_wrap(~ response) +
  labs(x = "Spearman rho (95% CI)", y = NULL, colour = "stretch",
       title = "RNAplfold mean log accessibility by region",
       subtitle = sprintf("n = %d..%d genes", min(reg$n), max(reg$n))) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "metric_screen_regions.jpg"), pr, width = 200, height = 90,
       units = "mm", dpi = 200)
