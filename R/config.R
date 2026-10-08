# =============================================================================
# Pipeline configuration
# =============================================================================
# Central configuration for the RNA half-life analysis pipeline.
# Edit values here to add species or change paths. Features are defined in
# R/feature_table.csv, which the feature-table section below reads.
# =============================================================================

# --- Paths -------------------------------------------------------------------

DATA_ROOT  <- "data"
RAW_DIR    <- file.path(DATA_ROOT, "raw")
SHARED_DIR <- file.path(RAW_DIR, "shared")
CACHE_DIR  <- file.path(DATA_ROOT, "cache")
OUTPUT_DIR <- file.path(DATA_ROOT, "outputs")
SPLITS_DIR <- file.path(DATA_ROOT, "splits")

# Bump this integer when feature-engineering logic changes so stale caches
# are regenerated instead of silently reused.
#
# v11: Saluki predictions reach the cache. Up to v10 the .rds keyed on
# `ensembl_gene_id`, join_gene_level() skipped it silently, and
# `saluki_prediction` never existed. See load_saluki_predictions().
#
# v12: the CAI column is `cai_cds` (was `cai`), so it carries its real region.
# v13: eej_dist_closest_start dropped; eej_dist_closest_stop renamed
#      eej_dist_closest_stop_mrna (whole-mRNA region token).
CACHE_VERSION <- 13L


# --- Region vocabulary -------------------------------------------------------
# Canonical internal region names are lowercase. Display names for plots are
# handled by `format_col_name()` in R/utils/naming.R.
#
# Regions are an ordered 5' -> 3' traversal of the transcript. `mrna` is the
# whole-mRNA member of regional families (length_mrna = length_5utr+cds+3utr)
# and, as of v4, is also the suffix used for genuinely whole-transcript
# scalar metrics (architecture, uORF, probing, NMD) — the v3 `transcript`
# pseudo-region and the `window`/`core`/`full` NMD pseudo-regions have been
# retired. There are no pseudo-regions: every token here is a real region.

REGIONS <- c("5utr", "cds", "3utr", "mrna", "utrpair",
             "last100", "start", "stop")

#' Raw → canonical region-token aliases.
#'
#' Some upstream pipelines emit verbose / suffixed region names. This is the
#' single mapping `normalise_region()` consults to bring them in line with
#' REGIONS. Add an entry when a new variant appears upstream; do NOT scatter
#' equivalent renames through individual loaders.
REGION_ALIASES <- c(
  tail_region        = "last100",
  start_codon_region = "start",
  stop_codon_region  = "stop"
)


# --- Nucleotide alphabet -----------------------------------------------------
# The schema is RNA-canonical: uracil, never thymine. Upstream `nuc_U_ratio_*`
# becomes `frac_u_*`, and format_col_name() labels it "nt.U%". Species differ
# upstream on codon spelling (human `codon_AAU`, mouse `codon_AAT`);
# normalise_codon_alphabet() in R/io/load_raw.R folds them onto this alphabet
# at load time, the same way normalise_region() applies REGION_ALIASES.
#
# Consequence for anyone writing a regex over composition columns: the triplet
# character class is `[acgtu]`, not `[acgt]`. Matching DNA-only silently drops
# the 37 U-containing codons rather than erroring.


# --- Species registry --------------------------------------------------------
# Add a new species by appending an entry. Each entry defines where to find
# its raw data and which columns to pull from shared files.

SPECIES_CONFIG <- list(
  
  human = list(
    dir        = "human",
    saluki_rds = "saluki_predictions.rds"    # relative to species dir; NULL if none
  ),
  
  mouse = list(
    dir        = "mouse",
    saluki_rds = "saluki_predictions.rds"
  )
)

# Species the analyses use by default — build_all() with no `species`
# argument stacks exactly these. The project is human-only as of 2026-09-23:
# mouse has no RNAfold/RNALfold features (never run through the folding
# pipeline) and carries mouse-only stop-free/stop-codon fraction columns, so a
# cross-species comparison would silently be sequence-only. Mouse stays
# registered above and still builds; add it back here to re-enable it
# everywhere at once.
ANALYSIS_SPECIES <- c("human")


# --- Feature table: the single source of truth -------------------------------
# R/feature_table.csv defines every feature: its columns (a regex), its place
# in the Supergroup > Group > Feature hierarchy, whether it is used in the
# exploratory analysis and in the model, its display names, its colour, and
# why. Everything below is DERIVED from it — edit the table, not these
# objects. scripts/check_feature_table.R verifies the table against a built
# cache (every column claimed by exactly one row, labels, regions, flags).
#
# Three levels, one flat selection namespace (ids must not collide):
#   feature     a table row; `feature_id`. The unit of FEATURE_PATTERNS.
#   group       the table's Group column, snake-cased ("Global folding" ->
#               global_folding). Every group sits in exactly one supergroup.
#   supergroup  the table's Supergroup column, snake-cased.
#
# The three flags:
#   Included in exploratory analysis   Selecting a group or supergroup returns
#       only its exploratory features. Naming a feature id directly always
#       works, so a plot can still reach an excluded feature on purpose.
#   Included in model                  The model's predictor list
#       (model_columns()). Nothing else reads it.
#   Included in core plots             The default set for the correlation
#       figures (CORE_FEATURES). A subset of the exploratory features.
# A row Excluded from both is never used anywhere; drop_excluded() removes it.
#
# Rows in the "Response / evaluation" supergroup (halflife, saluki) and rows
# with no `columns` (features not yet built) are documented but are not
# features: they get no FEATURE_PATTERNS entry and cannot be selected.

FEATURE_TABLE_PATH <- file.path(
  if (exists(".PIPELINE_ROOT")) .PIPELINE_ROOT else "R", "feature_table.csv")

.snake <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(x)))

.read_feature_table <- function(path) {
  t <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8", na.strings = character())
  flag <- function(col) {
    bad <- setdiff(t[[col]], c("Included", "Excluded"))
    if (length(bad)) stop(path, ": `", col, "` must be Included/Excluded; found ",
                          paste(bad, collapse = ", "), call. = FALSE)
    t[[col]] == "Included"
  }
  t$exploratory   <- flag("Included in exploratory analysis")
  t$model         <- flag("Included in model")
  t$core          <- flag("Included in core plots")
  t$group_id      <- .snake(t$Group)
  t$supergroup_id <- .snake(t$Supergroup)
  t$is_feature    <- nzchar(t$columns) & t$supergroup_id != "response_evaluation"

  # Invariants that need no data. The data-dependent ones are in
  # scripts/check_feature_table.R.
  dup <- t$feature_id[duplicated(t$feature_id)]
  if (length(dup)) stop(path, ": duplicate feature_id ", paste(dup, collapse = ", "),
                        call. = FALSE)
  split_groups <- names(which(tapply(t$supergroup_id, t$group_id,
                                     function(s) length(unique(s))) > 1))
  if (length(split_groups)) stop(path, ": group(s) in more than one supergroup: ",
                                 paste(split_groups, collapse = ", "), call. = FALSE)
  ids <- c(t$feature_id, unique(t$group_id), unique(t$supergroup_id))
  clash <- unique(ids[duplicated(ids)])
  if (length(clash)) stop(path, ": id used at more than one level: ",
                          paste(clash, collapse = ", "), call. = FALSE)
  t
}

FEATURE_TABLE <- .read_feature_table(FEATURE_TABLE_PATH)

.feature_rows <- FEATURE_TABLE[FEATURE_TABLE$is_feature, ]

#' Regex per feature. Mutually exclusive by construction — the checker
#' asserts no column matches two rows.
FEATURE_PATTERNS <- stats::setNames(as.list(.feature_rows$columns),
                                    .feature_rows$feature_id)

#' Group id -> feature ids; supergroup id -> feature ids. Table order.
FEATURE_GROUPS <- split(.feature_rows$feature_id,
                        factor(.feature_rows$group_id, unique(.feature_rows$group_id)))
SUPERGROUPS    <- split(.feature_rows$feature_id,
                        factor(.feature_rows$supergroup_id,
                               unique(.feature_rows$supergroup_id)))

EXPLORATORY_FEATURES <- .feature_rows$feature_id[.feature_rows$exploratory]
MODEL_FEATURES       <- .feature_rows$feature_id[.feature_rows$model]
CORE_FEATURES        <- .feature_rows$feature_id[.feature_rows$core]

#' Built columns that are used nowhere: rows Excluded from both flags. Regexes,
#' because the table is data-independent; drop_excluded() resolves them.
NEVER_USED_PATTERNS <- FEATURE_TABLE$columns[nzchar(FEATURE_TABLE$columns) &
                                             !FEATURE_TABLE$exploratory &
                                             !FEATURE_TABLE$model]
rm(.feature_rows)


# --- Cohort definition -------------------------------------------------------
# Two minimum UTR lengths, in nucleotides, that a transcript must reach to
# enter the analysis. A transcript failing either is dropped entirely — this
# is a ROW filter, and the counterpart to drop_excluded(), which is a column
# filter.
#
#   MIN_5UTR_LENGTH   30   5'UTR must be at least this long
#   MIN_3UTR_LENGTH  100   3'UTR must be at least this long
#
# WHY 5'UTR. A UTR of a few nucleotides is not a short UTR so much as an
# absent or mis-annotated one, and it poisons the regional features rather
# than merely weakening them. Folding energy over 12 nt is not comparable to
# folding energy over 1,200, GC content over 12 nt takes a handful of distinct
# values, and the length-normalised z-scores divide by a shuffled-sequence
# distribution that is itself near-degenerate.
#
# WHY 3'UTR 100. The "tail" region is the last 100 nt of the transcript. The
# sequence extraction extended a short 3'UTR's tail UPSTREAM into the CDS
# (wrong: it should extend downstream) and never extended a 3'UTR shorter
# than the tail. Rather than re-extract, transcripts with a 3'UTR under 100 nt
# are removed, which guarantees the tail lies wholly inside the 3'UTR. About
# 200 transcripts are affected by the extraction issue; the filter removes
# every one of them.
#
# NA COUNTS AS FAILING for both. A missing UTR length cannot be shown to clear
# the threshold, and reads as "no annotated UTR" rather than a failed
# measurement.
#
# WHERE IT IS APPLIED. build_dataset() applies it to the frame it RETURNS,
# after the cache is read or written — so the cache on disk stays complete and
# this needs no CACHE_VERSION bump. It is selection intent, like the feature
# table's flags, not a schema change. Pass `min_5utr = NULL, min_3utr = NULL`
# to build_dataset() / build_all() for the unfiltered table; the QC scripts do
# exactly that, because a coverage and missingness diagnostic should describe
# the whole built table including what this removes.
#
# THE SPLIT ARTEFACT DOES NOT NEED REBUILDING. Blocking is preserved under any
# subsetting (removing genes cannot make a family span two splits). Split
# proportions shift slightly with the larger removal; re-check with
# validate_splits() if the tolerance matters.

MIN_5UTR_LENGTH <- 30L
MIN_3UTR_LENGTH <- 100L


# --- Identity, family and split columns --------------------------------------
# Columns that identify a row rather than describe it. Distinct from
# never-used features in what happens to them: drop_excluded() REMOVES a
# never-used feature, whereas these must SURVIVE into a modelling frame — the
# family label is the grouping variable for blocked CV and the cluster for
# robust standard errors, so a script that dropped it could not do its job.
#
# The distinction that matters: never a predictor, always carried.

ID_COLS <- c("species", "transcript_id", "gene_id", "gene_name")

# Ingested from family.tsv by load_family() (R/io/load_raw.R). See
# FAMILY_CLUSTERING.md §1.7 for the seam's full column list; three of its
# columns are deliberately not ingested — `transcript_id` (supplied by
# load_transcripts(); keeping it would collide on join), `dataset` (constant
# per species, recorded in the cache's family provenance attribute instead)
# and `protein_len` (an exact restatement of length_cds/3 - 1; verified to
# correlate 1.000 with length_cds on the v8 human cache).
FAMILY_COLS <- c(
  "family_id_strict",   "family_size_strict",
  "family_id_medium",   "family_size_medium",
  "family_id_loose",    "family_size_loose",
  "family_searched",    "family_had_internal_stop"
)

# Everything a modelling script must exclude from its predictor matrix. The
# scripts build features as setdiff(names(df), c(META_COLS, TARGET_COL)), so
# any column NOT named here becomes a predictor by default — which is why
# family_size_* has to be listed. It is numeric, plausible-looking, and a
# property of the corpus rather than of the transcript.
# External model outputs carried for benchmarking. Never predictors: using a
# published half-life model's prediction to predict half-life would be
# leakage, and v11 is the first cache in which this column actually exists —
# so without this entry it would enter every setdiff()-built predictor matrix
# the moment the cache rebuilt.
BENCHMARK_COLS <- c("saluki_prediction")

META_COLS <- c(ID_COLS, FAMILY_COLS, BENCHMARK_COLS, "split")

# NOTE: family columns are deliberately absent from FEATURE_PATTERNS. Adding a
# `family` key there would make them reachable through select_features() and
# fg(), i.e. selectable AS FEATURES, which is the opposite of the intent.


# --- Blocking and split configuration ----------------------------------------
# See FAMILY_CLUSTERING.md §1.7 and §2.2a.

# Which clustering level blocks the splits. `medium` is the measured choice on
# human MANE: max family 282 (2.07% of the corpus), 10,605 families, and it
# splits the Ras superfamily along known subfamily lines where `loose` merges
# the whole superfamily. `strict` is unusable despite a reassuring 0.18% — its
# three largest families are all ZNF fragments, i.e. it shreds the largest real
# gene family in the genome.
#
# Every level is a column in the cache, so refitting at "loose" is a one-line
# sensitivity check rather than a rebuild.
BLOCK_LEVEL <- "medium"

# 80-10-10 holdout. `val` serves the purpose a nested inner CV loop would.
SPLIT_PROPS <- c(train = 0.8, val = 0.1, test = 0.1)

# Families larger than this fraction of the SMALLEST split are pinned to
# `train`. At 13,601 genes the smallest split is ~1,360, so the ceiling is ~68
# genes and the 282-member family cannot land in test and dominate 21% of it.
# The consequence must be reported with any result: the test set is depleted
# of large families, so it measures generalisation to small and singleton
# families rather than to all of them.
SPLIT_PIN_FRAC <- 0.05

# Fixed so the artefact is reproducible from the same family.tsv.
SPLIT_SEED <- 42L


# --- Helpers -----------------------------------------------------------------

#' Return the absolute path to a raw file for a given species.
#' @param species Character, one of names(SPECIES_CONFIG).
#' @param filename Character, filename within the species folder.
species_path <- function(species, filename) {
  stopifnot(species %in% names(SPECIES_CONFIG))
  file.path(RAW_DIR, SPECIES_CONFIG[[species]]$dir, filename)
}

#' Return the absolute path to a shared raw file.
shared_path <- function(filename) file.path(SHARED_DIR, filename)

#' Return the cache file path for a species.
cache_path <- function(species) {
  file.path(CACHE_DIR, sprintf("%s_dataset_v%d.rds", species, CACHE_VERSION))
}

#' Return the split-artefact path for a blocking level.
#'
#' Unversioned by CACHE_VERSION on purpose: the split is a property of
#' family.tsv and the seed, not of feature-engineering logic. Its traceability
#' comes from the family.tsv checksum stored inside the artefact, not from the
#' filename.
#'
#' @param level Character, a clustering level (strict / medium / loose).
#' @param ext Character, "rds" (the artefact) or "tsv" (the readable copy).
splits_path <- function(level = BLOCK_LEVEL, ext = "rds") {
  file.path(SPLITS_DIR, sprintf("holdout_%s.%s", level, ext))
}

#' Prefix payload columns with a tool or other prefix, leaving keys untouched.
#' Retained for backward compatibility; new code should prefer affix_payload().
prefix_payload <- function(df, prefix, keys = c("transcript_id", "region")) {
  dplyr::rename_with(df, ~ paste0(prefix, "_", .x), -dplyr::all_of(keys))
}

#' Affix payload columns with a prefix and/or suffix, leaving key columns
#' untouched. Generalises prefix_payload — used where a loader needs to push
#' a token to the END of the column name (e.g. a trailing region suffix)
#' rather than the front.
#'
#' @param df     A dataframe / tibble.
#' @param prefix Character prepended to every non-key column (default "").
#' @param suffix Character appended to every non-key column (default "").
#' @param keys   Character vector of key columns to leave untouched. Keys not
#'   present in `df` are ignored (unlike prefix_payload, which errors).
#' @return df with non-key columns renamed.
affix_payload <- function(df, prefix = "", suffix = "",
                          keys = c("transcript_id", "region")) {
  keep <- intersect(keys, names(df))
  dplyr::rename_with(df, ~ paste0(prefix, .x, suffix), -dplyr::all_of(keep))
}


# --- Hierarchy lookups -------------------------------------------------------

#' Reverse lookups: the supergroup / group a feature id belongs to.
#'
#' Works for every table row, including rows that are not selectable features
#' (never-used rows, the response rows).
#' @param feature Character vector of feature ids.
#' @return Character vector of ids; NA_character_ for unknown feature ids.
#' @examples
#' supergroup_of("rnafold_zscores")   # "structure"
#' group_of(c("gc", "at_skew"))       # c("nucleotide_composition", "nucleotide_asymmetry")
supergroup_of <- function(feature) {
  unname(stats::setNames(FEATURE_TABLE$supergroup_id, FEATURE_TABLE$feature_id)[feature])
}

group_of <- function(feature) {
  unname(stats::setNames(FEATURE_TABLE$group_id, FEATURE_TABLE$feature_id)[feature])
}
