# =============================================================================
# Display-name formatter for plot labels
# =============================================================================
# Converts canonical column names (lowercase, underscore-separated) into
# human-readable strings suitable for plot axes, legends and titles.
#
# The pipeline stores data with a canonical schema — e.g. "rnafold_zscore_5utr",
# "length_cds". This function is the bridge between that schema and the
# pretty labels end-users see on figures.
#
# Labels come from R/feature_table.csv; REPLACEMENTS covers only what one
# short name cannot (see below).
#
# --- v4 pseudo-region retirement (CACHE_VERSION 4L) --------------------------
# v3 introduced four pseudo-region tokens — `transcript`, `window`, `core`,
# `full`. v4 removes all of them. Every region-bearing column ends in one of
# the eight real region tokens (see `REGIONS` in R/config.R).
#
# The display-formatting consequences:
#   * Whole-transcript scalars (architecture, uORF, probing) carry the real
#     `mrna` suffix instead of `transcript`. They now render with a trailing
#     "mRNA" like every other mRNA-region column — e.g.
#     `intron_length_mean_mrna` -> "Mean intron length mRNA". The old
#     `transcript` rule (which rendered to "" and was stripped) is gone.
#   * The three NMD fragility windows are collapsed to one model carrying
#     `mrna`. The NMD metric-prefix rules are unchanged; combined with the
#     `mrna` region rule they yield e.g.
#     `nmd_fragile_codon_density_mrna` -> "NMD fragile codon density mRNA".
#     The `distal`/`nearest`/`any` region rules are gone.
# =============================================================================


# --- Labels from the feature table ---------------------------------------------
# The table's "Display Short name" is the label for every column of a row whose
# regex is one literal stem followed by a region token (`^mfe_delta_` +
# `cds` -> "MFE.Δ CDS") or is one literal column (`^cai$` -> "CAI"). That covers
# most rows; edit the label in R/feature_table.csv.
#
# A row whose short name has a <placeholder>, or whose columns differ by more
# than their region (uORF count vs n overlapping uORFs, intron mean/median/sd,
# icSHAPE nucleoplasm/cytoplasm, nucleotide fractions), cannot be labelled by
# one string: those fall through to the codon/amino-acid special case and the
# REPLACEMENTS rules below. scripts/check_feature_table.R checks their output
# against the table's short name too.

.TABLE_LABELS <- local({
  rows <- FEATURE_TABLE[nzchar(FEATURE_TABLE$columns) &
                        !grepl("<", FEATURE_TABLE$`Display Short name`, fixed = TRUE), ]
  lit <- regmatches(rows$columns, regexpr("^\\^[a-z0-9_]*", rows$columns))
  lit <- sub("^\\^", "", lit)
  data.frame(pattern = rows$columns, literal = lit,
             whole   = rows$columns == paste0("^", lit, "$"),
             label   = rows$`Display Short name`, stringsAsFactors = FALSE)
})

.table_label <- function(name) {
  hit <- which(vapply(.TABLE_LABELS$pattern, grepl, logical(1), x = name))
  if (length(hit) != 1) return(NULL)
  r <- .TABLE_LABELS[hit, ]
  toks <- strsplit(name, "_", fixed = TRUE)[[1]]
  region <- if (r$whole) {
    if (length(toks) > 1 && toks[length(toks)] %in% REGIONS) toks[length(toks)] else ""
  } else substring(name, nchar(r$literal) + 1)
  if (nzchar(region) && !region %in% REGIONS) return(NULL)   # not one stem
  if (nzchar(region)) paste(r$label, REGION_DISPLAYS[[region]]) else r$label
}


# --- Fallback rules ----------------------------------------------------------------
# Substitutions applied in order, for columns the table cannot label (see
# above) and for identifier columns outside the table. Region-suffix rules are
# last so they see the bare region token after the prefix rules.
REPLACEMENTS <- list(
  # --- Multi-column rows ---
  list("^intron_length_mean_",         "Mean intron length "),
  list("^uorf_count_",                 "uORF count "),
  list("^codon_",                      "codon."),     # codon_other_cds only
  list("^frac_a",                      "nt.A% "),
  list("^frac_c",                      "nt.C% "),
  list("^frac_g",                      "nt.G% "),
  list("^frac_u",                      "nt.U% "),
  list("^gini_nucleoplasm_",           "icSHAPE.nuc "),
  list("^gini_cytoplasm_",             "icSHAPE.cyto "),

  # --- Identifiers / metadata ---
  list("^gene_name$",                  "Gene name"),
  list("^gene_id$",                    "Ensembl gene ID"),
  list("^transcript_id$",              "Transcript ID"),
  list("^prediction_difference$",      "Prediction difference"),
  list("^species$",                    "Species"),

  # --- Region suffixes ---
  # Each rule matches a single leading separator — space OR underscore —
  # plus the region token, end-anchored (`[ _]<region>$`), so a residual `_`
  # before the region cannot leak a lowercase token, and `[ _]stop$` cannot
  # fire on the `stop` inside `alt-stop`.
  list("[ _]5utr$",                        " 5' UTR"),
  list("[ _]3utr$",                        " 3' UTR"),
  list("[ _]cds$",                         " CDS"),
  list("[ _]mrna$",                        " mRNA"),
  list("[ _]utrpair$",                     " UTR interactions"),
  list("[ _]last100$",                     " last 100 nt"),
  list("[ _]start$",                       " start codon"),
  list("[ _]stop$",                        " stop codon")
)


#' Format a canonical column name into a display string
#'
#' @param col_name Character (length-1 or vector). A canonical column name.
#' @return Character vector of the same length.
#' @examples
#' format_col_name("rnafold_zscore_5utr")   # "MFE.z 5' UTR"
#' format_col_name("length_cds")            # "Length CDS"
#' format_col_name("halflife")              # "Half-life"
#' format_col_name("gc_content_5utr")       # "C+G% 5' UTR"
#' format_col_name("eej_dist_closest_start")            # "EEJ.closest start codon"
#' format_col_name("intron_length_mean_mrna")           # "Mean intron length mRNA"
#' format_col_name("nmd_snv_fragile_codon_density_mrna")# "NMD.frag. mRNA"
#' format_col_name("codon_aaa_cds")         # "codon.AAA%"
#' format_col_name("aa_l_cds")              # "aa.L%"
#' @export
format_col_name <- function(col_name) {
  vapply(col_name, format_single_name, character(1), USE.NAMES = FALSE)
}

format_single_name <- function(name) {
  # Codon / amino-acid composition: token uppercased, region suffix dropped.
  # Cheap special-case avoids 84 exact-match rules in REPLACEMENTS and
  # avoids needing PCRE2 case-folding (which R's sub() doesn't support).
  # No trailing space on either: these return EARLY, so they never reach the
  # trimws() at the foot of this function the way every REPLACEMENTS-driven
  # name does. They previously carried one, which the documented examples above
  # ("Codon.AAA%", "aa.L%") already said they should not — visible as a ragged
  # right edge wherever codon labels are set flush, i.e. most of the gain
  # importance figure, since codons are 64 of the 106 baseline columns.
  m <- regmatches(name, regexec("^codon_([acgtu]{3})_cds$", name))[[1]]
  if (length(m) == 2) return(paste0("codon.", toupper(m[2]), "%"))

  m <- regmatches(name, regexec("^aa_([a-z])_cds$", name))[[1]]
  if (length(m) == 2) return(paste0("aa.", toupper(m[2]), "%"))

  lab <- .table_label(name)
  if (!is.null(lab)) return(lab)

  for (rule in REPLACEMENTS) {
    name <- sub(rule[[1]], rule[[2]], name)
  }
  name <- gsub("_", " ", name, fixed = TRUE)
  name <- gsub("  +", " ", name)
  trimws(name)
}