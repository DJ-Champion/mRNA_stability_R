# =============================================================================
# Pipeline configuration
# =============================================================================
# Central configuration for the RNA half-life analysis pipeline.
# Edit values here to add species, change paths, or register new feature groups.
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
CACHE_VERSION <- 10L


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


# --- Feature groups ----------------------------------------------------------
# Regex patterns defining named groups of columns. Use downstream with
# `fg()` / `fg_columns()`. Add a new group by appending a named entry.
#
# Display labels for the selection-key namespace below — FEATURE_PATTERNS keys,
# SUPERGROUPS names, and GROUP_BUNDLES names — live in R/utils/palettes.R as
# FEATURE_GROUP_DISPLAY_NAMES / SUPERGROUP_DISPLAY_NAMES / BUNDLE_DISPLAY_NAMES,
# formatted via format_group_name(). Add a label there, not here, when a key
# needs a nicer display string than its title-cased name. The individual
# columns inside the `standalone` group (cai, translation_efficiency,
# orfexondensity) are labelled as COLUMNS via format_col_name()
# (R/utils/naming.R), not as group keys — see the note on `standalone` below.
#
# INVARIANT: the patterns are mutually exclusive. No column may match two
# groups — a plot that builds a column -> group map would double-assign it.
# When one family's prefix is a prefix of another's, anchor the broader one
# more tightly (e.g. `gc` is `^gc_content_`, not `^gc_`, so it does not
# swallow the `skews` columns; `exons` excludes `^exon_density_`).

FEATURE_PATTERNS <- list(
  lengths        = "^length_",
  gc             = "^gc_content_",
  nmd            = "^nmd_",
  introns        = "^intron_",
  exons          = "^exon_(count|length)_",
  noncoding      = "^noncoding_",
  rnafold_scores = "^rnafold_score_",
  rnafold_zscores = "^rnafold_zscore_",
  rnafold_per_nt = "^rnafold_per_nt_",
  mfe_deltas     = "^mfe_delta_",
  mfe_expected   = "^mfe_expected_",
  rnalfold_scores    = "^rnalfold_score_",
  rnalfold_zscores   = "^rnalfold_zscore_",
  junctions      = "^junctions_",
  eej_dist       = "^eej_dist_",
  uorfs          = "^(uorf_|dist_cap_)",
  exon_density   = "^exon_density_",
  stopfree       = "^stopfree_",
  skews          = "^(gc|at)_skew_",
  codon_freqs    = "^codon_",
  aa_freqs       = "^aa_",
  nuc_ratios     = "^frac_",
  compositional  = "^(purine_|amino_)",
  probing        = "^gini_",

  # Genuinely region-less whole-transcript scalars. One group rather than one
  # group per column: they behave identically to every consumer (no region
  # suffix, mapped to the `mrna` slot in region-aware plots) and differ only
  # in their display label, which comes from format_col_name().
  standalone     = "^(cai|translation_efficiency|orfexondensity)$"
)


# --- Feature supergroups -----------------------------------------------------
# Coarse-grained categorisation of FEATURE_PATTERNS keys for plotting, palette
# assignment, and any analysis that wants to colour or facet by category
# family.
#
# INVARIANT: every FEATURE_PATTERNS key belongs to exactly one supergroup —
# no key omitted, no key listed twice. `standalone` lives in `other`.
#
# Update this when adding new groups to FEATURE_PATTERNS.

SUPERGROUPS <- list(
  structure  = c("rnafold_scores", "rnafold_zscores", "rnafold_per_nt",
                 "mfe_deltas", "mfe_expected",
                 "rnalfold_scores", "rnalfold_zscores",
                 "probing"),

  intrinsic  = c("lengths", "gc", "stopfree", "skews", "codon_freqs", "aa_freqs",
                 "nuc_ratios", "compositional"),

  splicing   = c("junctions", "introns", "exons", "noncoding", "eej_dist"),

  translation = c("uorfs", "exon_density"),
  decay       = c("nmd"),
  other       = c("standalone")
)

GROUP_BUNDLES <- list(
  nmd_core = list(
    groups = "nmd",
    pick = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                        "nmd_alt_stop_codon_density_mrna"))
  ),
  lengths_core = list(
    groups = "lengths",
    pick = list(lengths = c("length_5utr", "length_cds",
                            "length_3utr", "length_mrna"))
  ),
  splicing_core = list(
    groups = "eej_dist",
    pick = list(eej_dist = c("eej_dist_closest_start", "eej_dist_closest_stop"))
  ),
  structure_core = list(
    groups = "structure",
    pick = list(probing = c("gini_cytoplasm_mrna", "gini_cytoplasm_5utr",
                            "gini_cytoplasm_cds", "gini_cytoplasm_3utr"))
  ),
  intrinsic_core = list(
    groups = c("intrinsic", "standalone"),
    pick = list(lengths = c("length_5utr", "length_cds",
                            "length_3utr", "length_mrna"),
                stopfree = c("stopfree_length_5utr", "stopfree_length_3utr",
                             "stopfree_length_cds", "stopfree_length_mrna"),
                standalone = "cai")
  ),
  intrinsic_select = list(
    groups = c("intrinsic", "standalone"),
    pick = list(lengths = c("length_5utr", "length_cds",
                            "length_3utr", "length_mrna"),
                stopfree = c("stopfree_length_5utr", "stopfree_length_3utr",
                             "stopfree_length_cds", "stopfree_length_mrna"),
                codon_freqs = c("codon_agu_cds", "codon_uca_cds"),
                aa_freqs = c("aa_s_cds", "aa_v_cds"),
                standalone = "cai")
  ),
  translation_core = list(
    groups = c("uorfs", "exon_density"),
    pick = list(uorfs = c("uorf_present_mrna"))
  )
)

# The main set of groups we've decided to focus on for almost everything.
# Pass anywhere a `groups =` argument is accepted; it is the default for the
# correlation dotplot, the response scatter, the region heatmap and the
# correlation-heatmap workflow. Editing this changes what those plots show by
# default — it is selection intent, so it never needs a CACHE_VERSION bump.
INCLUDED_GROUPS <- c("nmd_core", "splicing_core", "structure_core",
                     "intrinsic_select", "translation_core")


# --- Cohort definition -------------------------------------------------------
# The minimum length, in nucleotides, that BOTH UTRs must reach for a
# transcript to enter the analysis. A transcript failing it is dropped
# entirely — this is a ROW filter, and the counterpart to EXCLUDED_FEATURES
# below, which is a column filter.
#
# WHY. A UTR of a few nucleotides is not a short UTR so much as an absent or
# mis-annotated one, and it poisons the regional features rather than merely
# weakening them. Every per-region metric becomes degenerate at that scale:
# folding energy over 12 nt is not comparable to folding energy over 1,200,
# GC content over 12 nt takes a handful of distinct values, and the
# length-normalised z-scores divide by a shuffled-sequence distribution that
# is itself near-degenerate. The models cannot tell "this UTR is unstructured"
# from "this UTR is barely there".
#
# 30 nt is the project's agreed threshold. Measured on the v10 caches:
#
#   human  13,660 -> 12,302 built rows      (1,358 removed, 9.9%)
#          13,601 -> 12,277 modellable      (1,324 removed, 9.7%)
#   mouse  14,197 -> 13,215 built rows        (982 removed, 6.9%)
#
# The cut falls overwhelmingly on the 5' side: 1,244 human transcripts have a
# 5'UTR under 30 nt against 50 with a short 3'UTR.
#
# NA COUNTS AS FAILING. A missing UTR length cannot be shown to clear the
# threshold, and no transcript in either cache has a recorded length of 0
# while the minimum observed is 1 nt — so NA reads as "no annotated UTR", not
# as a failed measurement. 72 human and 302 mouse rows are dropped on this
# branch, and they are included in the totals above.
#
# WHERE IT IS APPLIED. build_dataset() applies it to the frame it RETURNS,
# after the cache is read or written — so the cache on disk stays complete and
# this needs no CACHE_VERSION bump. It is selection intent, like
# INCLUDED_GROUPS and EXCLUDED_FEATURES, not a schema change. Pass
# `min_utr = NULL` to build_dataset() / build_all() for the unfiltered table;
# the QC scripts do exactly that, because a coverage and missingness diagnostic
# should describe the whole built table including what this removes.
#
# THE SPLIT ARTEFACT DOES NOT NEED REBUILDING. Blocking is preserved under any
# subsetting (removing genes cannot make a family span two splits), and the
# proportions barely move: 80.12 / 9.99 / 9.90 against a target of 80/10/10,
# inside validate_splits()'s tolerance. holdout_medium.rds therefore remains
# valid, and results stay traceable to the clustering run behind it.

MIN_UTR_LENGTH <- 30L


# --- Model-pipeline exclusions -----------------------------------------------
# Columns that are built and cached, but are NOT admissible as covariates.
# Applied by drop_excluded() (R/utils/feature_groups.R) at the head of a
# modelling / screening script — never inside build_dataset(), for two reasons:
#
#   1. Several entries here are ENGINEERING SCAFFOLDING: engineer.R derives
#      kept columns from them. mfe_expected_* feeds mfe_delta_*;
#      junctions_count_* feeds exon_density_*; eej_dist_{up,down}stream_*
#      feeds eej_dist_closest_*. Remove them before engineering and the
#      retained derivatives disappear too.
#   2. QC scripts legitimately want them — analysis/qc/mfe_expected_check.R
#      reads the mfe_expected_* family this list removes wholesale.
#
# The cache therefore stays complete and this needs no CACHE_VERSION bump: it
# is selection intent, like INCLUDED_GROUPS above, not a schema change.
#
# Why a flat denylist rather than GROUP_BUNDLES pick/drop: 43 of these match no
# FEATURE_PATTERNS key at all (the rnafold/rnalfold median+pval families, and
# the legacy scalars). select_features() only ever returns columns belonging to
# a registered group, so the pick/drop machinery cannot address them. A growing
# table also fails safer with a denylist — a newly ingested feature arrives
# INCLUDED and surfaces in screening, rather than being silently excluded.
#
# Species differ: mouse lacks 37 of these. drop_excluded() uses any_of().

EXCLUDED_FEATURES <- c(
  # Vienna auxiliary distribution statistics. The score/zscore/per-nt families
  # are the modelled encoding of the same folding runs; median and p-value are
  # retained upstream for provenance only. Human-only (absent in mouse).
  "rnafold_median_3utr", "rnafold_median_5utr", "rnafold_median_cds",
  "rnafold_median_last100", "rnafold_median_mrna", "rnafold_median_start",
  "rnafold_median_stop", "rnafold_median_utrpair",
  "rnafold_pval_3utr", "rnafold_pval_5utr", "rnafold_pval_cds",
  "rnafold_pval_last100", "rnafold_pval_mrna", "rnafold_pval_start",
  "rnafold_pval_stop", "rnafold_pval_utrpair",
  "rnalfold_median_3utr", "rnalfold_median_5utr", "rnalfold_median_cds",
  "rnalfold_median_last100", "rnalfold_median_mrna", "rnalfold_median_start",
  "rnalfold_median_stop",
  "rnalfold_pval_3utr", "rnalfold_pval_5utr", "rnalfold_pval_cds",
  "rnalfold_pval_last100", "rnalfold_pval_mrna", "rnalfold_pval_start",
  "rnalfold_pval_stop",

  # Scaffolding for mfe_delta_* (= rnafold_score - mfe_expected). Deterministic
  # in gc_content_{region} and length_{region}, so it carries no information the
  # delta and its two inputs do not already hold. Empties the `mfe_expected`
  # group key; `structure` keeps its other seven members.
  "mfe_expected_5utr", "mfe_expected_cds", "mfe_expected_3utr",
  "mfe_expected_mrna", "mfe_expected_last100", "mfe_expected_start",
  "mfe_expected_stop",

  # Fixed-width analysis windows, not measured transcript properties.
  # length_last100 takes 2 distinct values across the human table.
  "length_last100", "length_start", "length_stop",

  # Exon / intron architecture summaries superseded by the retained
  # exon_density_* encoding.
  "exon_count_internal_mrna", "internal_exon_mean", "internal_exon_median",
  "internal_exon_sd", "intron_length_mean_mrna", "intron_median",
  "intron_sd", "n_exons",

  # Junction counts and their per-kb densities. engineer.R turns the counts
  # into exon_density_* ((junctions + 1) / kb), which is what the models use;
  # this drops the whole `junctions` group and keeps `exon_density`.
  "junctions_count_5utr", "junctions_count_cds", "junctions_count_3utr",
  "junctions_count_mrna",
  "junctions_density_5utr", "junctions_density_cds", "junctions_density_3utr",
  "junctions_density_mrna",

  # Directional exon-exon junction distances. add_eej_min_distance() collapses
  # each up/downstream pair into the retained eej_dist_closest_{start,stop}.
  "eej_dist_upstream_stop", "eej_dist_downstream_stop",
  "stop_dist_last_downstream",
  "eej_dist_upstream_start", "eej_dist_downstream_start",

  # Exact duplicate of length_5utr (verified identical on the v8 human cache).
  "utr5_length",

  # uORF counts and distances. Retains uorf_present_mrna, which is what
  # translation_core picks.
  "uorf_count_mrna", "n_overlapping_uorfs", "total_classical_uorf_codons",
  "max_classical_uorf_codons", "dist_cap_to_first_uatg_mrna",
  "dist_last_uorf_stop_to_main_atg",

  # CDS size / denominator columns. The codon and aa families are row-normalised
  # fractions (add_codon_aa_fractions), so these are the length proxies that
  # normalisation exists to remove.
  "cds_length_codons_cds", "n_codons_scored_cds", "n_stops_cds",

  # Leaves `standalone` as cai + translation_efficiency, both retained as
  # candidate covariates.
  "orfexondensity"
)


# --- Identity, family and split columns --------------------------------------
# Columns that identify a row rather than describe it. Distinct from
# EXCLUDED_FEATURES in what happens to them: drop_excluded() REMOVES an
# excluded feature, whereas these must SURVIVE into a modelling frame — the
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
META_COLS <- c(ID_COLS, FAMILY_COLS, "split")

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


# --- Supergroup helpers ------------------------------------------------------

#' Reverse lookup: which supergroup does a FEATURE_PATTERNS key belong to?
#' @param group Character vector of FEATURE_PATTERNS keys.
#' @return Character vector of supergroup names; NA_character_ for unknown keys.
#' @examples
#' supergroup_of("rnafold_zscores")   # "structure"
#' supergroup_of(c("gc", "nmd"))      # c("intrinsic", "splicing")
supergroup_of <- function(group) {
  vapply(group, function(g) {
    hits <- names(SUPERGROUPS)[vapply(
      SUPERGROUPS, function(members) g %in% members, logical(1)
    )]
    if (length(hits) == 0) NA_character_ else hits[1]
  }, character(1), USE.NAMES = FALSE)
}


