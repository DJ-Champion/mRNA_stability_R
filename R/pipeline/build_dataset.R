# =============================================================================
# Build dataset — main pipeline entry point
# =============================================================================
# `build_dataset(species)` is what analysis scripts call. It:
#   1. Checks the cache and returns it unless rebuild = TRUE.
#   2. Loads every raw source (tolerant of missing files).
#   3. Assembles into a wide tibble keyed by transcript_id.
#   4. Joins gene-level attributes (halflife, TE, probing, etc.).
#   5. Runs feature engineering.
#   6. Saves to cache.
#
# The return is a single tibble with a `species` column, ready to stack with
# other species via `bind_rows()`.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


#' Drop transcripts whose UTRs are too short to carry meaningful features
#'
#' The cohort filter. A transcript is kept only if `length_5utr` is present
#' and at least `min_5utr` AND `length_3utr` is present and at least
#' `min_3utr`; NA fails, because a missing length cannot be shown to clear the
#' threshold. See the MIN_5UTR_LENGTH / MIN_3UTR_LENGTH block in config.R for
#' the rationale (the 3' threshold guarantees the last-100-nt tail lies wholly
#' inside the 3'UTR).
#'
#' This is the row-filter counterpart to drop_excluded(). Unlike that function
#' it is applied for you, by build_dataset(), rather than called at the head of
#' each analysis script.
#'
#' Safe to apply at any point after the join: engineer_features() derives every
#' column row-wise, with no statistic pooled across transcripts, so removing
#' rows cannot change the feature values of the rows that remain.
#'
#' @param df A dataset from build_dataset() / build_all().
#' @param min_5utr Integer. Minimum 5'UTR length. NULL or NA skips this filter.
#' @param min_3utr Integer. Minimum 3'UTR length. NULL or NA skips this filter.
#' @param verbose Logical. If TRUE (default) report how many rows went.
#' @return `df` without the short-UTR transcripts.
#' @export
drop_short_utr <- function(df, min_5utr = MIN_5UTR_LENGTH,
                           min_3utr = MIN_3UTR_LENGTH, verbose = TRUE) {
  skip5 <- is.null(min_5utr) || is.na(min_5utr)
  skip3 <- is.null(min_3utr) || is.na(min_3utr)
  if (skip5 && skip3) return(df)

  cols <- c(if (!skip5) "length_5utr", if (!skip3) "length_3utr")
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    warning("drop_short_utr: no ", paste(missing_cols, collapse = " / "),
            " column — cohort filter NOT applied")
    return(df)
  }

  pass5 <- if (skip5) rep(TRUE, nrow(df))
           else !is.na(df$length_5utr) & df$length_5utr >= min_5utr
  pass3 <- if (skip3) rep(TRUE, nrow(df))
           else !is.na(df$length_3utr) & df$length_3utr >= min_3utr
  keep <- pass5 & pass3

  if (verbose) {
    message("drop_short_utr: removed ", sum(!keep), " of ", nrow(df),
            " transcripts (", sum(!pass5), " with 5'UTR < ", min_5utr,
            " nt or missing, ", sum(!pass3), " with 3'UTR < ", min_3utr,
            " nt or missing); ", sum(keep), " remain")
  }

  # Belt and braces, and currently redundant: base `[` on a tibble preserves
  # custom attributes, as do filter/mutate/left_join/bind_rows (measured, not
  # assumed). It is kept because attribute preservation for arbitrary
  # attributes is not part of dplyr's API contract, and because losing this one
  # fails SILENTLY rather than loudly — attach_splits() skips its provenance
  # comparison when either md5 is NULL, so a dataset built from a stale
  # family.tsv would simply stop being flagged. Two lines against a silent
  # failure is a trade worth making; delete them if the attribute ever becomes
  # a documented guarantee.
  prov <- attr(df, "family_provenance")
  out  <- df[keep, , drop = FALSE]
  if (!is.null(prov)) attr(out, "family_provenance") <- prov
  out
}


#' Build (or load) the full dataset for one species
#'
#' @param species Character. Must be a key of SPECIES_CONFIG.
#' @param rebuild Logical. If FALSE (default) and a cache exists for this
#'   species at the current CACHE_VERSION, returns the cache. If TRUE, rebuilds
#'   from raw files.
#' @param min_5utr,min_3utr Integer. Transcripts with a 5'UTR / 3'UTR shorter
#'   than these are dropped from the RETURNED frame — see drop_short_utr() and
#'   the MIN_5UTR_LENGTH / MIN_3UTR_LENGTH block in config.R. Pass NULL for
#'   either to skip that filter (the QC scripts pass NULL for both). The cache written to disk is always complete, so this
#'   never invalidates it and never needs a CACHE_VERSION bump.
#' @return A wide-form tibble with a `species` column.
#' @export
build_dataset <- function(species, rebuild = FALSE,
                          min_5utr = MIN_5UTR_LENGTH,
                          min_3utr = MIN_3UTR_LENGTH) {
  if (!species %in% names(SPECIES_CONFIG)) {
    stop("Unknown species '", species, "'. Known: ",
         paste(names(SPECIES_CONFIG), collapse = ", "))
  }

  if (!rebuild) {
    cached <- load_snapshot(species)
    if (!is.null(cached)) {
      message("Loaded cache for ", species, " (", nrow(cached), " rows, ",
              ncol(cached), " columns)")
      return(drop_short_utr(cached, min_5utr, min_3utr))
    }
  }

  message("Building dataset for ", species, " (no cache, or rebuild forced)...")
  
  # --- 1. Load raw sources --------------------------------------------------
  message("Loading raw sources...")
  transcripts <- load_transcripts(species)
  if (is.null(transcripts)) {
    stop("transcripts.csv is required but missing for species '", species, "'")
  }
  
  # Regional (Long form)
  regional <- list(
    rnafold        = load_rnafold(species),
    rnalfold       = load_rnalfold(species),
    rnaup          = load_rnaup(species),
    sequence_basic = load_sequence_basic(species), # REPLACEMENT
    stopfree       = load_stopfree(species)
  )
  
  # Transcript Level (Wide form)
  transcript_level <- list(
    nmd           = load_nmd_fragility(species),                  # v4: single model
    architecture  = load_architecture(species),
    junctions     = load_junctions_wide(species),                 # REPLACEMENT
    uorfs         = load_uorfs(species),                          # UPDATED
    codon_aa      = load_codon_aa_counts(species),                # REPLACEMENT
    cai           = load_cai(species),
    nte           = load_nte(species)
  )
  
  halflife_df  <- load_halflife(species)
  te_df        <- load_translation_efficiency(species)
  agarwal_df   <- load_agarwal_features(species)
  saluki_df    <- load_saluki_predictions(species)
  gini_df      <- load_gini_probing(species)        # gene-level, joins on gene_id
  family_df    <- load_family(species)              # gene-level, joins on gene_id
  
  purrr::imap(regional, function(df, nm) {
    if (is.null(df) || !all(c("transcript_id", "region") %in% names(df))) return(NULL)
    dups <- df |>
      dplyr::count(transcript_id, region) |>
      dplyr::filter(n > 1)
    if (nrow(dups) > 0) {
      message(nm, ": ", nrow(dups), " duplicate (transcript_id, region) combos")
      print(dplyr::slice_head(dups, n = 5))
    }
  })
  
  # --- 2. Pivot regional data to wide --------------------------------------
  message("Pivoting regional data to wide form...")
  wide <- pivot_regional_to_wide(regional)
  if (is.null(wide)) {
    # Fall back to transcripts metadata if no regional data exists yet.
    wide <- tibble::tibble(transcript_id = transcripts$transcript_id)
  }
  
  # --- 3. Attach metadata and gene-level features --------------------------
  wide <- left_join(wide, transcripts, by = "transcript_id")
  
  # Transcript-level features
  wide <- join_transcript_level(wide, transcript_level)
  
  # Some external sources key on gene_name rather than transcript_id/gene_id.
  # Apply halflife first (it's the response variable — loudest if missing).
  if (!is.null(halflife_df)) {
    join_key <- if ("gene_id" %in% names(halflife_df)) "gene_id"
    else if ("gene_name"   %in% names(halflife_df)) "gene_name"
    else NA_character_
    if (is.na(join_key)) {
      warning("halflife data has no gene_id or gene_name — skipped")
    } else {
      wide <- left_join(wide, halflife_df, by = join_key)
    }
  }
  
  if (!is.null(te_df)) {
    join_key <- if ("gene_id" %in% names(te_df)) "gene_id"
    else if ("gene_name"   %in% names(te_df)) "gene_name"
    else NA_character_
    if (!is.na(join_key)) wide <- left_join(wide, te_df, by = join_key)
  }
  
  if (!is.null(agarwal_df))  wide <- join_gene_level(wide, list(agarwal_df))
  if (!is.null(saluki_df))   wide <- join_gene_level(wide, list(saluki_df))
  if (!is.null(gini_df))     wide <- join_gene_level(wide, list(gini_df))

  # Family labels. Not features — blocking metadata (see FAMILY_COLS in
  # config.R). The row count is asserted either side because this join is the
  # one place a malformed seam could silently multiply the corpus, and a
  # dataset that has quietly doubled still looks entirely plausible.
  if (!is.null(family_df)) {
    n_before <- nrow(wide)
    wide <- join_gene_level(wide, list(family_df))
    if (nrow(wide) != n_before) {
      stop("joining family.tsv changed the row count (", n_before, " -> ",
           nrow(wide), ") — gene_id is not unique on one side of the join")
    }

    # Genes the clustering never saw. Expected to be rows that arrived via the
    # half-life join without sequence features, so they carry no family and
    # cannot be modelled anyway. Report rather than assert: the count is
    # informative, and a rise in it is the signal that family.tsv has drifted
    # out of step with the rest of the raw inputs.
    n_unfamilied <- sum(is.na(wide$family_id_medium) & !is.na(wide$gene_id))
    if (n_unfamilied > 0) {
      message("  ", n_unfamilied, " gene(s) present in the dataset but absent ",
              "from family.tsv (no family label; excluded from blocked splits)")
    }
  }


  # --- 4. Add the species marker -------------------------------------------
  wide$species <- species
  
  # Put identifier columns first for readability
  id_cols <- intersect(ID_COLS, names(wide))
  wide <- wide[, c(id_cols, setdiff(names(wide), id_cols))]

  # --- 5. Feature engineering ----------------------------------------------
  message("Running feature engineering...")
  wide <- engineer_features(wide)

  # --- 6. Save cache -------------------------------------------------------
  # Family provenance is attached AFTER engineering, not at load time: dplyr
  # verbs are not reliable about carrying custom attributes through, and a
  # left_join keeps the attributes of x only. Setting it here, on the finished
  # frame, is the one point where survival to disk is guaranteed.
  if (!is.null(family_df)) {
    attr(wide, "family_provenance") <- family_provenance()
  }

  save_snapshot(wide, species)

  # AFTER save_snapshot, deliberately. The cache is the complete built table;
  # the cohort filter is selection intent applied to what callers receive, so
  # the two never have to be kept in step and changing the UTR thresholds does not
  # invalidate a single cache.
  drop_short_utr(wide, min_5utr, min_3utr)
}


#' Build datasets for multiple species and stack them
#'
#' @param species Character vector. Defaults to ANALYSIS_SPECIES (config.R);
#'   pass names(SPECIES_CONFIG) for every registered species.
#' @param rebuild Logical, passed through to build_dataset.
#' @param min_5utr,min_3utr Integer or NULL, passed through to build_dataset. Applied
#'   per species before stacking, which is the same result as filtering after —
#'   the threshold is a property of one transcript.
#' @return A single tibble with a `species` column. Columns absent from a
#'   given species are NA for that species' rows.
#' @export
build_all <- function(species = ANALYSIS_SPECIES, rebuild = FALSE,
                      min_5utr = MIN_5UTR_LENGTH, min_3utr = MIN_3UTR_LENGTH) {
  dfs <- lapply(species, build_dataset, rebuild = rebuild,
                min_5utr = min_5utr, min_3utr = min_3utr)
  dplyr::bind_rows(dfs)
}