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


#' Build (or load) the full dataset for one species
#'
#' @param species Character. Must be a key of SPECIES_CONFIG.
#' @param rebuild Logical. If FALSE (default) and a cache exists for this
#'   species at the current CACHE_VERSION, returns the cache. If TRUE, rebuilds
#'   from raw files.
#' @return A wide-form tibble with a `species` column.
#' @export
build_dataset <- function(species, rebuild = FALSE) {
  if (!species %in% names(SPECIES_CONFIG)) {
    stop("Unknown species '", species, "'. Known: ",
         paste(names(SPECIES_CONFIG), collapse = ", "))
  }
  
  if (!rebuild) {
    cached <- load_snapshot(species)
    if (!is.null(cached)) {
      message("Loaded cache for ", species, " (", nrow(cached), " rows, ",
              ncol(cached), " columns)")
      return(cached)
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

  wide
}


#' Build datasets for multiple species and stack them
#'
#' @param species Character vector. Defaults to all species in SPECIES_CONFIG.
#' @param rebuild Logical, passed through to build_dataset.
#' @return A single tibble with a `species` column. Columns absent from a
#'   given species are NA for that species' rows.
#' @export
build_all <- function(species = names(SPECIES_CONFIG), rebuild = FALSE) {
  dfs <- lapply(species, build_dataset, rebuild = rebuild)
  dplyr::bind_rows(dfs)
}