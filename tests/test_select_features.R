# Tests for select_features(). Run from the project root:
#   Rscript tests/test_select_features.R
# Uses the cached human dataset (data/cache/), like scripts/check_feature_table.R.

suppressMessages(source("R/load_all.R"))
df <- suppressMessages(build_dataset("human", min_utr = MIN_UTR_LENGTH))

n_fail <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  if (!ok) n_fail <<- n_fail + 1L
  cat(if (ok) "  ok    " else "  FAIL  ", desc, "\n", sep = "")
}
errors <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
cols_of <- function(...) selected_columns(select_features(df, ...))
fcols   <- function(ids) unlist(lapply(ids, fg_columns, df = df), use.names = FALSE)

cat("expansion\n")
core <- select_features(df)
check("default is core", identical(core, select_features(df, "core")))
check("core = columns of CORE_FEATURES, table order",
      identical(core$column, fcols(CORE_FEATURES)))
check("feature_id carries the source feature",
      all(mapply(function(co, g) co %in% fg_columns(df, g), core$column, core$feature_id)))
check("plain view is a character vector", is.character(selected_columns(core)))
check("exploratory", identical(cols_of("exploratory"), fcols(EXPLORATORY_FEATURES)))
check("model", identical(cols_of("model"), fcols(MODEL_FEATURES)))
check("supergroup expands to exploratory members only",
      identical(cols_of("structure"),
                fcols(intersect(SUPERGROUPS$structure, EXPLORATORY_FEATURES))))
check("group expands to exploratory members only",
      identical(cols_of("global_folding"),
                fcols(intersect(FEATURE_GROUPS$global_folding, EXPLORATORY_FEATURES))))
check("a non-exploratory feature id is reachable directly",
      length(cols_of("rnafold_median")) > 0 && !"rnafold_median" %in% EXPLORATORY_FEATURES)
check("mixed levels, no duplicates",
      !anyDuplicated(cols_of(c("structure", "global_folding", "rnafold_scores", "core"))))
check("include order does not change output order",
      identical(cols_of(c("gc", "lengths")), cols_of(c("lengths", "gc"))))
check("unknown include errors", errors(select_features(df, "nonsense")))

cat("exclude\n")
check("exclude a feature", !any(cols_of(exclude = "probing") %in% fg_columns(df, "probing")))
check("exclude removes exactly those columns",
      setequal(cols_of("exploratory", exclude = c("codon_freqs", "aa_freqs")),
               setdiff(cols_of("exploratory"), c(fg_columns(df, "codon_freqs"),
                                                 fg_columns(df, "aa_freqs")))))
check("exclude a supergroup removes all of it",
      !any(cols_of("exploratory", exclude = "sequence") %in% fcols(SUPERGROUPS$sequence)))
check("exclude a flag", length(cols_of("exploratory", exclude = "core")) ==
      length(setdiff(cols_of("exploratory"), cols_of("core"))))
check("excluding everything gives an empty result",
      nrow(select_features(df, "core", exclude = "core")) == 0)
check("unknown exclude errors", errors(select_features(df, exclude = "nonsense")))
check("probing is not in core", !any(fg_columns(df, "probing") %in% core$column))

cat("top_n\n")
r <- function(co) abs(cor(df[[co]], df$halflife, method = "spearman", use = "complete.obs"))
codons <- fg_columns(df, "codon_freqs")
top2 <- select_features(df, "codon_freqs", top_n = list(codon_freqs = 2),
                           response = "halflife")
best2 <- names(sort(sapply(codons, r), decreasing = TRUE))[1:2]
check("keeps exactly the two strongest codon stems", setequal(top2$column, best2))
check("matches the ranked figure's halflife top two (AGU, AUC)",
      setequal(top2$column, c("codon_agu_cds", "codon_auc_cds")))
check("independent of what else is selected",
      identical(
        subset(select_features(df, c("codon_freqs", "gc"),
                                  top_n = list(codon_freqs = 2), response = "halflife"),
               feature_id == "codon_freqs")$column, top2$column))
check("differs by response",
      !setequal(top2$column,
                selected_columns(select_features(df, "codon_freqs",
                  top_n = list(codon_freqs = 2), response = "translation_efficiency"))))
untrimmed <- function(x) { x <- subset(x, !feature_id %in% c("codon_freqs", "aa_freqs"))
                           rownames(x) <- NULL; x }
check("leaves other features alone",
      identical(untrimmed(select_features(df, top_n = list(codon_freqs = 2, aa_freqs = 2),
                                             response = "halflife")),
                untrimmed(core)))
rg <- select_features(df, "lengths", top_n = list(lengths = 1), response = "halflife")
check("keeps every region of a kept stem",
      length(unique(.column_stem(rg$column))) == 1 && nrow(rg) >= 1)
check("N larger than the family keeps it all",
      nrow(select_features(df, "aa_freqs", top_n = list(aa_freqs = 999),
                              response = "halflife")) == length(fg_columns(df, "aa_freqs")))
check("top_n reports what it trimmed",
      identical(attr(top2, "top_n")$codon_freqs,
                list(requested = 2, available = 64L, bound = TRUE)))

cat("regions\n")
cds <- select_features(df, "core", regions = "cds")
check("regions keeps only that region",
      all(column_regions(cds$column, cds$feature_id) == "cds") && nrow(cds) > 0)
check("regions = NULL keeps everything", nrow(select_features(df, "core", regions = NULL)) == nrow(core))
check("region-less features count as mrna",
      "cai" %in% select_features(df, "core", regions = "mrna")$feature_id)
check("unknown region errors", errors(select_features(df, regions = "nowhere")))
utr <- select_features(df, "codon_freqs", regions = c("cds", "5utr"),
                          top_n = list(codon_freqs = 2), response = "halflife")
check("top_n ranks after the region filter", nrow(utr) == 2 && all(grepl("_cds$", utr$column)))

cat("top_n (errors)\n")
check("top_n without response errors",
      errors(select_features(df, top_n = list(codon_freqs = 2))))
check("top_n with an unknown feature errors",
      errors(select_features(df, top_n = list(nope = 2), response = "halflife")))

cat("against the step-1 baseline (ranked figure, core = old default less probing)\n")
base <- read.csv("scripts/selection_baseline.csv", stringsAsFactors = FALSE)
for (resp in c("halflife", "translation_efficiency")) {
  b <- subset(base, script == "feature_correlation_ranked.R" & job == resp &
                stage == "plotted" & feature_id != "probing")
  sel <- select_features(df, top_n = list(codon_freqs = 2, aa_freqs = 2), response = resp)
  key <- function(stem, region) paste(stem, region)
  region <- sub("^.*_", "", sel$column); stem <- .column_stem(sel$column)
  # region-less columns (cai, ...) are plotted in the mrna slot
  region <- ifelse(region %in% REGIONS, region, "mrna")
  check(paste(resp, ": same (stem, region) set as the baseline"),
        setequal(key(stem, region), key(b$metric_stem, b$region)))
}

if (n_fail) { cat("\n", n_fail, " FAILED\n", sep = ""); quit(status = 1) }
cat("\nAll passed.\n")
