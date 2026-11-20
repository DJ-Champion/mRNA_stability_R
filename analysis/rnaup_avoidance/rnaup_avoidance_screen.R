# =============================================================================
# RNAup ncRNA avoidance at the start codon vs half-life and TE
# =============================================================================
# Every rnaup_g* column from rnaup_avoidance_combined.tsv (load_rnaup_avoidance(
# all = TRUE)) against half-life and translation efficiency, as raw Spearman
# rho with a 95% CI.
#
# The columns are joined onto the cached dataset IN MEMORY from the loader,
# replacing the two the cache holds, so the per-family and G_total columns can
# be screened without entering the cache.
#
# Unadjusted: the z-scores are composition-controlled by construction; the raw
# sums are not and track window GC (the `vs_gc` column reports how much).
#
# Sign: energies are kcal/mol, more negative = stronger binding; z > 0 =
# avoidance (binds the ncRNAs less than its dinucleotide shuffles). The
# pre-specified pair is the G_int total z and sum (the two cache columns).
#
# Usage: Rscript analysis/rnaup_avoidance/rnaup_avoidance_screen.R
# Writes CSVs and figures to data/outputs/rnaup_avoidance/.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/load_all.R")

OUT <- file.path(OUTPUT_DIR, "rnaup_avoidance")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

responses <- c("halflife", "translation_efficiency")
# Start-region features already in the table, for overlap.
neighbours <- c("gc_content_start", "rnafold_score_start", "rnafold_zscore_start",
                "mfe_delta_start", "rnaplfold_meanlog_unpaired_u1_start")
PRIMARY <- RNAUP_AVOIDANCE_COLUMNS

ru <- load_rnaup_avoidance("human", all = TRUE)
df <- build_dataset("human") |> drop_excluded(verbose = FALSE)
df <- df |> select(-starts_with("rnaup_g")) |> inner_join(ru, by = "transcript_id")
feat <- setdiff(names(ru), "transcript_id")
cat("cohort transcripts with a value:", nrow(df), "\n")

screen_one <- function(v, y) {
  d <- df[c(v, y, "gc_content_start")]; d <- d[complete.cases(d), ]
  n <- nrow(d); d[] <- lapply(d, rank)
  raw <- cor(d[[v]], d[[y]])
  se <- 1.06 / sqrt(n - 3)
  data.frame(feature = v, response = y, n = n, rho = raw,
             lo = tanh(atanh(raw) - 1.96 * se), hi = tanh(atanh(raw) + 1.96 * se),
             gc_start = cor(d[[v]], d$gc_content_start))
}
tab <- do.call(rbind, lapply(responses, function(y)
  do.call(rbind, lapply(feat, screen_one, y = y))))
tab$q <- p.adjust(2 * pnorm(-abs(atanh(tab$rho)) * sqrt(tab$n - 3) / 1.06), "BH")

# Parse rnaup_<energy>[_<family>]_<stat>_start.
p <- regmatches(tab$feature, regexec(
  "^rnaup_(gint|gtotal)_(?:(trna|rrna|srp|snrna|yrna|snorna)_)?(.+)_start$",
  tab$feature, perl = TRUE))
tab$energy <- vapply(p, `[`, "", 2)
tab$family <- vapply(p, function(m) if (nzchar(m[3])) m[3] else "all 34", "")
tab$stat   <- vapply(p, `[`, "", 4)
write.csv(tab, file.path(OUT, "screen.csv"), row.names = FALSE)

fmt <- function(t) t |>
  transmute(response, feature = sub("^rnaup_", "", feature), n,
            rho = round(rho, 3), lo = round(lo, 3), hi = round(hi, 3),
            q = signif(q, 2), vs_gc = round(gc_start, 2)) |>
  as.data.frame()

cat("\n=== Pre-specified (cache columns) ===\n")
print(fmt(tab |> filter(feature %in% PRIMARY)), row.names = FALSE)
for (y in responses) {
  cat("\n=== ", y, ": all columns by |rho| ===\n", sep = "")
  print(fmt(tab |> filter(response == y) |> arrange(desc(abs(rho)))),
        row.names = FALSE)
}

# --- Overlap with start-region features already in the table ------------------
nb <- intersect(neighbours, names(df))
ov <- cor(df[c(PRIMARY, "rnaup_gtotal_zscore_start", "rnaup_gint_trna_zscore_start")],
          df[nb], method = "spearman", use = "pairwise.complete.obs")
cat("\n=== Spearman with existing start-region features ===\n"); print(round(ov, 2))
write.csv(round(ov, 3), file.path(OUT, "overlap_start_features.csv"))

# Z distribution: is there avoidance on average?
z <- df$rnaup_gint_zscore_start
cat(sprintf("\nG_int total z: median %.3f, mean %.3f, share > 0 %.3f, |z| > 2 %.3f, n %d\n",
            median(z, na.rm = TRUE), mean(z, na.rm = TRUE), mean(z > 0, na.rm = TRUE),
            mean(abs(z) > 2, na.rm = TRUE), sum(!is.na(z))))

# --- Figure: rho by family, for the z-scores and sums --------------------------
resp_lab <- c(halflife = "Half-life", translation_efficiency = "Translation efficiency")
fam_lev  <- c("all 34", "trna", "rrna", "srp", "snrna", "yrna", "snorna")
fam_lab  <- c(`all 34` = "All 34", trna = "tRNA (20)", rrna = "rRNA (1)", srp = "SRP (1)",
              snrna = "snRNA (5)", yrna = "Y RNA (3)", snorna = "snoRNA (4)")
pd <- tab |> filter(stat %in% c("zscore", "sum")) |>
  mutate(family = factor(fam_lab[family], rev(fam_lab[fam_lev])),
         panel = paste(ifelse(energy == "gint", "G_int", "G_total"),
                       ifelse(stat == "zscore", "z-score", "sum")),
         response = resp_lab[response])
pl <- ggplot(pd, aes(rho, family)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, colour = "grey25") +
  facet_grid(panel ~ response) +
  labs(x = "Spearman rho (95% CI)", y = NULL,
       title = "RNAup ncRNA pairing at the start codon vs response",
       subtitle = sprintf("n = %d..%d transcripts", min(tab$n), max(tab$n))) +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "screen.jpg"), pl, width = 210, height = 210,
       units = "mm", dpi = 200)
message("wrote ", OUT)
