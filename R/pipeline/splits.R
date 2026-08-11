# =============================================================================
# Blocked train / val / test splits
# =============================================================================
# Assigns every gene to a split such that no sequence family is ever divided
# across two splits. The guarantee this exists to provide:
#
#   no sequence in the test split has a close homologue in the train split
#
# Training on mouse `Rpl13a` and testing on human `RPL13A` is leakage, and so
# is training on one member of a paralogue pair and testing on the other. The
# family labels ingested by load_family() are connected components over a
# pairwise homology graph, so holding a whole family in one split is exactly
# the guarantee above — leakage is a pairwise property, and connected
# components are the transitive closure of pairwise edges.
#
# THE ARTEFACT IS THE POINT. build_splits() is run ONCE and writes
# gene_id -> split to disk; every downstream script reads that file via
# load_splits() / attach_splits(). Regenerating the assignment inside
# modelling code means a rebuilt family.tsv, a different R version, or even a
# changed row order can silently reshuffle which genes are in test — results
# stop being reproducible with nothing looking wrong. The checksum of the
# family.tsv that produced the split travels inside the artefact so the
# assignment is always traceable to the clustering run behind it.
#
# See FAMILY_CLUSTERING.md §2.2a for the design and §2.4 for the assertions.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


# --- Packers -----------------------------------------------------------------

#' Assign families to k equal folds, balancing on gene count
#'
#' Greedy longest-processing-time (LPT) bin-packing: sort families by size
#' descending, assign each to the currently-emptiest fold. The same algorithm
#' `sklearn.model_selection.GroupKFold` uses; it lands within a few percent of
#' optimal balance.
#'
#' Balancing is not optional. Family sizes are wildly unequal — one family of
#' 282 alongside 8,797 singletons on human MANE — so assigning families to
#' folds at random gives folds that differ several-fold in size.
#' `caret::groupKFold` partitions groups without balancing on group size and
#' is not adequate here.
#'
#' Not used by this project's 80-10-10 design (see assign_holdout, which
#' generalises it), and retained for a k-fold one.
#'
#' @param fam A data frame with a `family_id` column, one row per gene.
#' @param k Integer, number of folds.
#' @param seed Integer, RNG seed for tie-breaking.
#' @return A tibble(family_id, fold).
#' @export
assign_folds <- function(fam, k, seed = SPLIT_SEED) {
  stopifnot("family_id" %in% names(fam), k >= 2L)

  sizes <- fam |>
    dplyr::filter(!is.na(family_id)) |>
    dplyr::count(family_id, name = "n")

  # Randomise BEFORE ordering by size. Without it, equal-sized families break
  # ties by whatever order the table happened to arrive in, and every seed
  # produces an identical partition.
  set.seed(seed)
  sizes <- sizes[sample(nrow(sizes)), ]
  sizes <- dplyr::arrange(sizes, dplyr::desc(n))

  load <- numeric(k)
  out  <- integer(nrow(sizes))
  for (i in seq_len(nrow(sizes))) {
    j       <- which.min(load)
    out[i]  <- j
    load[j] <- load[j] + sizes$n[i]
  }

  dplyr::tibble(family_id = sizes$family_id, fold = out)
}


#' Assign families to named splits with unequal targets
#'
#' Two changes from assign_folds.
#'
#' **Unequal bins.** The equal-fold packer assigns each family to the emptiest
#' bin, which is only correct when bins are the same size. For unequal targets
#' the family goes to whichever bin is furthest below its quota — measured as
#' a FRACTION of that bin's quota, not in absolute genes. The distinction is
#' invisible at equal bin sizes (the quotas divide out) and decisive at
#' 80-10-10; see the comment on the packing loop below.
#'
#' **Large families are pinned to the first split.** With 13,601 genes a 10%
#' test split is ~1,360 genes, and blocking keeps families intact — so the
#' 282-member family would be ~21% of the test set wherever it landed. A
#' single holdout gives one draw and no error bar, so rather than hope, the
#' biggest families go to training deliberately, by a stated criterion
#' ("larger than `pin_frac` of the smallest split") rather than by eye.
#'
#' The consequence must be reported with any result that uses this split: the
#' held-out splits contain no family above the pin threshold, so their
#' performance is not an unbiased estimate for a randomly chosen gene. They
#' measure generalisation to small and mid-sized families, not to the largest
#' ones. (Only families above the threshold are affected — the packing rule
#' below keeps the rest of the size distribution stratified across all three
#' splits.)
#'
#' @param fam A data frame with a `family_id` column, one row per gene.
#' @param props Named numeric vector of split proportions, summing to 1. The
#'   FIRST name receives the pinned families.
#' @param pin_frac Numeric. Families larger than `pin_frac * min(props) *
#'   n_genes` are pinned to the first split. NULL or 0 disables pinning.
#' @param seed Integer, RNG seed for tie-breaking.
#' @return A tibble(family_id, split).
#' @export
assign_holdout <- function(fam,
                           props    = SPLIT_PROPS,
                           pin_frac = SPLIT_PIN_FRAC,
                           seed     = SPLIT_SEED) {
  stopifnot("family_id" %in% names(fam),
            !is.null(names(props)),
            abs(sum(props) - 1) < 1e-9)

  sizes <- fam |>
    dplyr::filter(!is.na(family_id)) |>
    dplyr::count(family_id, name = "n")
  n_tot <- sum(sizes$n)

  # Families too large to sit in the smallest split without dominating it.
  pin_above <- if (is.null(pin_frac) || pin_frac <= 0) Inf else {
    pin_frac * min(props) * n_tot
  }
  pinned <- dplyr::filter(sizes, n >  pin_above)
  rest   <- dplyr::filter(sizes, n <= pin_above)

  target <- props * n_tot
  load   <- stats::setNames(numeric(length(props)), names(props))
  load[1] <- sum(pinned$n)                 # pinned genes already in split 1

  set.seed(seed)
  rest <- rest[sample(nrow(rest)), ]
  rest <- dplyr::arrange(rest, dplyr::desc(n))

  # Furthest below quota AS A FRACTION OF THAT QUOTA — not in absolute genes.
  #
  # Note this rule only ever chooses a DESTINATION for a whole family; no
  # branch of it can divide one. What follows is about which families each
  # split ends up with, never about their integrity.
  #
  # The absolute form `which.max(target - load)` balances gene counts
  # perfectly and silently destroys the family structure of the small splits.
  # On human MANE at `medium`: train's deficit starts at ~10,600 genes against
  # val/test's 1,360, so train wins every comparison until it has absorbed
  # ~9,500 genes. Families are placed largest-first, and all 4,804 genes in
  # multi-member families fit inside that head start — so every last one lands
  # in train, and val and test receive nothing but singletons. Measured:
  #
  #   absolute:  test = 1,360 families / 1,360 genes, max family 1,  0% in
  #              multi-member families
  #   relative:  test = 1,061 families / 1,360 genes, max family 20, 35% in
  #              multi-member families (train 35%, val 35%)
  #
  # Both give an 80.00/10.00/10.00 gene split, which is why the failure is
  # invisible in the obvious summary. But a test set of pure singletons does
  # not measure what the model will be asked to do, and singletons are not a
  # random sample of genes — their median 3'UTR is 1,181 nt against 1,576 for
  # genes with paralogues. Dividing by `target` makes every split give up the
  # same PROPORTION of itself, so family structure is stratified across the
  # three rather than pooled into the largest.
  #
  # Equal bins make the two forms identical (the quotas divide out), so this
  # is a correction to the unequal-bin generalisation specifically, and
  # assign_folds() above is right as it stands.
  out <- character(nrow(rest))
  for (i in seq_len(nrow(rest))) {
    j       <- which.max((target - load) / target)
    out[i]  <- names(load)[j]
    load[j] <- load[j] + rest$n[i]
  }

  dplyr::bind_rows(
    dplyr::tibble(family_id = pinned$family_id, split = names(props)[1]),
    dplyr::tibble(family_id = rest$family_id,   split = out)
  )
}


# --- Validation --------------------------------------------------------------

#' Assert the properties a blocked split must have
#'
#' Cheap, and they catch the failure this whole exercise exists to prevent.
#' Run before modelling, not after.
#'
#' @param assigned A data frame with `gene_id`, `family_id` and a split column.
#' @param split_col Character, the column holding the assignment.
#' @param props Named numeric vector of intended proportions, or NULL to skip
#'   the balance check.
#' @param tol Numeric, allowed absolute deviation from each intended
#'   proportion. Blocking makes exact quotas unreachable — a family is
#'   indivisible, so the packer can only get close.
#' @return `assigned`, invisibly.
#' @export
validate_splits <- function(assigned, split_col = "split",
                            props = SPLIT_PROPS, tol = 0.01) {
  stopifnot(all(c("gene_id", "family_id", split_col) %in% names(assigned)))

  sp <- assigned[[split_col]]

  # 1. No family spans two splits — the core guarantee.
  spanning <- assigned |>
    dplyr::group_by(family_id) |>
    dplyr::summarise(n_splits = dplyr::n_distinct(.data[[split_col]]),
                     .groups = "drop") |>
    dplyr::filter(n_splits > 1)
  if (nrow(spanning) > 0) {
    stop(nrow(spanning), " family/families span more than one ", split_col,
         " — the blocking guarantee is broken (e.g. ",
         paste(utils::head(spanning$family_id, 3), collapse = ", "), ")")
  }

  # 2. Every gene assigned.
  if (anyNA(sp)) {
    stop(sum(is.na(sp)), " gene(s) have no ", split_col)
  }

  # 3. No gene in more than one split. Distinct from (1): this catches a
  #    duplicated gene row rather than a divided family.
  dup <- assigned |>
    dplyr::count(species, gene_id) |>
    dplyr::filter(n > 1)
  if (nrow(dup) > 0) {
    stop(nrow(dup), " gene(s) appear in more than one row of the split table")
  }

  # 4. Splits are close to their intended proportions. A large miss means the
  #    blocking level produced a family too big to pack — go back to
  #    family_qc.tsv and pick a stricter level.
  if (!is.null(props)) {
    actual <- prop.table(table(factor(sp, levels = names(props))))
    off    <- abs(as.numeric(actual[names(props)]) - as.numeric(props))
    if (any(off > tol)) {
      warning("split proportions deviate from target by more than ", tol,
              ": ", paste(sprintf("%s %.4f (want %.4f)", names(props),
                                  as.numeric(actual[names(props)]),
                                  as.numeric(props)), collapse = "; "),
              ". A large miss means a family too big to pack at this ",
              "blocking level.")
    }
  }

  # 5. Family STRUCTURE is comparable across splits, not just gene counts.
  #    Checks (4) can pass perfectly while a split holds only singletons —
  #    that was a real failure of the absolute-deficit packer, and nothing
  #    else here detects it. A held-out split of pure singletons measures
  #    generalisation to genes with no paralogues, which is not the question
  #    being asked, and singletons differ systematically from the rest.
  if ("family_size" %in% names(assigned)) {
    struct <- assigned |>
      dplyr::group_by(.data[[split_col]]) |>
      dplyr::summarise(pct_multi = mean(family_size > 1), .groups = "drop")
    spread <- diff(range(struct$pct_multi))
    if (spread > 0.10) {
      warning("family structure differs markedly across ", split_col, "s: ",
              paste(sprintf("%s %.1f%% of genes in multi-member families",
                            struct[[split_col]], 100 * struct$pct_multi),
                    collapse = "; "),
              ". A split of mostly singletons does not measure generalisation ",
              "to the corpus. Check the packer's deficit rule.")
    }
  }

  invisible(assigned)
}


# --- The artefact ------------------------------------------------------------

#' Build the split artefact and write it to disk
#'
#' Reads family.tsv directly rather than the built cache. That is deliberate:
#' the split depends only on the family labels, so building it this way needs
#' no cache to exist first and cannot be perturbed by a feature-engineering
#' change.
#'
#' @param level Character, clustering level to block on (a suffix of the
#'   `family_id_*` columns).
#' @param species Character vector. Defaults to every registered species;
#'   those absent from the cohort are skipped, since load_family() returns
#'   NULL for them.
#' @param props,pin_frac,seed Passed to assign_holdout().
#' @param write Logical. FALSE returns the table without touching disk.
#' @return A tibble(species, gene_id, family_id, family_size, split) with a
#'   `family_provenance` attribute, invisibly.
#' @export
build_splits <- function(level    = BLOCK_LEVEL,
                         species  = names(SPECIES_CONFIG),
                         props    = SPLIT_PROPS,
                         pin_frac = SPLIT_PIN_FRAC,
                         seed     = SPLIT_SEED,
                         write    = TRUE) {

  id_col   <- paste0("family_id_",   level)
  size_col <- paste0("family_size_", level)

  # Pool every species in the cohort into ONE packing run. Families span
  # species by construction (that is what makes them block cross-species
  # leakage), so packing each species separately would put two halves of an
  # orthologue group in different splits — precisely the leakage the families
  # exist to prevent.
  fam <- lapply(species, function(sp) {
    df <- load_family(sp)
    if (is.null(df)) return(NULL)
    if (!id_col %in% names(df)) {
      stop("family.tsv has no column '", id_col, "' — is '", level,
           "' a level the clustering run emitted?")
    }
    dplyr::tibble(species     = sp,
                  gene_id     = df$gene_id,
                  family_id   = df[[id_col]],
                  family_size = df[[size_col]])
  })
  fam <- dplyr::bind_rows(fam)

  if (nrow(fam) == 0) {
    stop("no family assignments found for species: ",
         paste(species, collapse = ", "),
         ". Is data/raw/shared/family.tsv present?")
  }
  if (anyNA(fam$family_id)) {
    stop(sum(is.na(fam$family_id)), " gene(s) in family.tsv have no ", id_col,
         " — the seam guarantees every gene a family at every level, ",
         "singletons included, so this is an upstream bug")
  }

  assigned <- assign_holdout(fam, props = props, pin_frac = pin_frac,
                             seed = seed) |>
    dplyr::right_join(fam, by = "family_id") |>
    dplyr::select(species, gene_id, family_id, family_size, split) |>
    dplyr::arrange(species, gene_id)

  validate_splits(assigned, props = props)

  attr(assigned, "family_provenance") <- family_provenance()
  attr(assigned, "split_params") <- list(
    level = level, props = props, pin_frac = pin_frac, seed = seed,
    built_at = Sys.time(), r_version = R.version.string
  )

  summarise_splits(assigned, props = props, pin_frac = pin_frac)

  if (write) {
    dir.create(SPLITS_DIR, showWarnings = FALSE, recursive = TRUE)
    saveRDS(assigned, splits_path(level, "rds"))
    # A plain TSV alongside, for eyeballing and for non-R consumers. The RDS
    # is the artefact of record — it is the one that carries the provenance.
    readr::write_tsv(assigned, splits_path(level, "tsv"))
    message("Split artefact written: ", splits_path(level, "rds"),
            " (", nrow(assigned), " genes)")
  }

  invisible(assigned)
}


#' Print what the split actually did, including its known bias
#'
#' @param assigned Output of build_splits().
#' @param props,pin_frac As passed to assign_holdout().
#' @return `assigned`, invisibly.
#' @export
summarise_splits <- function(assigned, props = SPLIT_PROPS,
                             pin_frac = SPLIT_PIN_FRAC) {
  n_tot <- nrow(assigned)
  tab <- assigned |>
    dplyr::group_by(split) |>
    dplyr::summarise(n_genes    = dplyr::n(),
                     pct        = round(100 * dplyr::n() / n_tot, 2),
                     n_families = dplyr::n_distinct(family_id),
                     max_family = max(family_size),
                     # The column that exposes a singleton-only split. Should
                     # be near-identical across splits; see validate_splits().
                     pct_multi  = round(100 * mean(family_size > 1), 1),
                     .groups    = "drop")

  message("\nSplit summary (", n_tot, " genes, ",
          dplyr::n_distinct(assigned$family_id), " families):")
  print(as.data.frame(tab), row.names = FALSE)
  message(
    "  max_family = largest family that split received, ENTIRELY ",
    "(families are never divided).\n",
    "  pct_multi  = % of the split's genes having a relative, which is ",
    "always in the same split.\n",
    "               These should be close across splits — that is what makes ",
    "the held-out\n",
    "               splits resemble train.")

  # The pinning bias, stated rather than left implicit.
  if (!is.null(pin_frac) && pin_frac > 0) {
    pin_above <- pin_frac * min(props) * n_tot
    n_pinned  <- assigned |>
      dplyr::filter(family_size > pin_above) |>
      dplyr::summarise(f = dplyr::n_distinct(family_id), g = dplyr::n())
    if (n_pinned$f > 0) {
      message(sprintf(
        paste0("\n%d famil(y/ies) larger than %.0f genes (%.0f%% of the ",
               "smallest split) pinned to '%s', covering %d genes.\n",
               "  => the held-out splits contain NO family above %.0f genes. ",
               "They measure\n     generalisation to small and mid-sized ",
               "families, not to the largest ones.\n     Report it as such."),
        n_pinned$f, pin_above, 100 * pin_frac, names(props)[1], n_pinned$g,
        pin_above))
    }
  }

  invisible(assigned)
}


#' Read the split artefact from disk
#'
#' @param level Character, clustering level.
#' @return The split tibble, or NULL if no artefact exists.
#' @export
load_splits <- function(level = BLOCK_LEVEL) {
  path <- splits_path(level, "rds")
  if (!file.exists(path)) {
    message("  skip (missing): ", path,
            " — run scripts/build_splits.R to create it")
    return(NULL)
  }
  readRDS(path)
}


#' Attach the split assignment to a dataset
#'
#' Joins on (species, gene_id) so a stacked multi-species frame is handled
#' correctly. Genes with no assignment get `NA` and are reported — they are
#' outside the clustered cohort and must not be quietly treated as training
#' data.
#'
#' @param df A dataset from build_dataset() / build_all().
#' @param level Character, clustering level.
#' @param check_provenance Logical. If TRUE (default), warn when the
#'   family.tsv behind the split differs from the one behind the dataset.
#' @return `df` with a `split` column.
#' @export
attach_splits <- function(df, level = BLOCK_LEVEL, check_provenance = TRUE) {
  sp <- load_splits(level)
  if (is.null(sp)) return(df)

  if (check_provenance) {
    a <- attr(df, "family_provenance")$md5
    b <- attr(sp, "family_provenance")$md5
    if (!is.null(a) && !is.null(b) && !identical(a, b)) {
      warning("the split artefact was built from a DIFFERENT family.tsv than ",
              "this dataset (", substr(b, 1, 8), " vs ", substr(a, 1, 8),
              "). The clustering has been re-run since the split was made, so ",
              "family membership and the split no longer correspond. Rebuild ",
              "the split with scripts/build_splits.R.")
    }
  }

  out <- dplyr::left_join(df, dplyr::select(sp, species, gene_id, split),
                          by = c("species", "gene_id"))

  n_missing <- sum(is.na(out$split))
  if (n_missing > 0) {
    message(n_missing, " row(s) have no split assignment (absent from the ",
            "clustered cohort) — filter them out before modelling rather ",
            "than letting them default into training")
  }
  out
}
