# =============================================================================
# Feature-group tidyselect helpers + group resolution / selection
# =============================================================================
# Three layers, kept deliberately separate:
#
#   1. SCHEMA      R/feature_table.csv, read by config.R into FEATURE_PATTERNS
#                  (regex per feature), FEATURE_GROUPS and SUPERGROUPS (the two
#                  coarser levels), and the exploratory / model flags.
#
#   2. INTENT      Reusable named selections live in GROUP_BUNDLES (config.R).
#                  A bundle is the reusable form of plotting/modelling intent —
#                  the proper home for what the old `_some` groups were trying
#                  to be. One-off intent is the per-call pick/drop arguments.
#                  Neither touches the schema.
#
#   3. RESOLUTION  resolve_selection() normalises (groups, pick, drop) +
#                  bundles into a (group_keys, pick, drop) triple.
#                  expand_groups() is the group-key-only view of that.
#                  select_features() flattens the triple to actual COLUMN
#                  names present in `df`. Plots that need group identity per
#                  column (e.g. the correlation dotplot) consume the triple
#                  directly. Column selection is thus defined in one place.
#
# --- Bundles -----------------------------------------------------------------
# A GROUP_BUNDLES entry is a list with any of:
#   groups : character vector of group / supergroup / OTHER bundle names
#   pick   : named list  group_key -> columns to KEEP from that group
#   drop   : named list  group_key -> columns to REMOVE from that group
# A bare character vector is accepted as shorthand for list(groups = <vec>).
#
# Merge policy when a caller passes pick/drop AND names a bundle carrying
# pick/drop for the SAME group key: the CALLER wins for that key (per-group
# replacement); the bundle's entry applies only where the caller is silent.
# Resolution order is always pick-then-drop, so a caller drop trims whatever
# a bundle pick produced.
#
# Selection keys are feature ids, group ids, supergroup ids and bundle names.
# Expanding a group or supergroup yields only its EXPLORATORY features (the
# table's "Included in exploratory analysis" flag); naming a feature id
# directly always yields it. pick/drop are keyed by feature id.
#
# Usage:
#   df %>% select(fg("rnafold_zscores"))
#   df %>% select(all_of(select_features(df, groups = "structure")))
#   select_features(df, groups = "global_folding")        # a group
#   select_features(df, groups = "lengths_core")          # a bundle
#   select_features(df, groups = "nmd_susceptibility",
#                   pick = list(nmd_snv_fragile = "nmd_snv_fragile_codon_density_mrna"))
# =============================================================================


#' Return a tidyselect spec for a named feature group
#'
#' @param group Character, a feature id (a key of FEATURE_PATTERNS).
#' @return A tidyselect spec usable inside dplyr::select().
#' @export
fg <- function(group) {
  if (!group %in% names(FEATURE_PATTERNS)) {
    stop("Unknown feature '", group, "'. Known features: ",
         paste(names(FEATURE_PATTERNS), collapse = ", "))
  }
  tidyselect::matches(FEATURE_PATTERNS[[group]])
}


#' List the columns that match a named feature group in a given dataframe
#'
#' @param df A dataframe.
#' @param group Character, a feature id (a key of FEATURE_PATTERNS).
#' @return Character vector of column names, in `df` column order.
#' @export
fg_columns <- function(df, group) {
  if (!group %in% names(FEATURE_PATTERNS)) {
    stop("Unknown feature '", group, "'.")
  }
  grep(FEATURE_PATTERNS[[group]], names(df), value = TRUE)
}


# --- internal: bundle registry accessor --------------------------------------
.group_bundles <- function() {
  if (exists("GROUP_BUNDLES", inherits = TRUE)) GROUP_BUNDLES else list()
}

# --- internal: normalise a bundle entry to list(groups, pick, drop) ----------
.as_bundle <- function(entry) {
  if (is.character(entry)) entry <- list(groups = entry)
  list(
    groups = if (is.null(entry$groups)) character() else entry$groups,
    pick   = if (is.null(entry$pick))   list()       else entry$pick,
    drop   = if (is.null(entry$drop))   list()       else entry$drop
  )
}


#' Resolve a selection (groups + caller pick/drop) into a normalised triple
#'
#' Expands supergroups and bundles, merges any pick/drop the bundles carry
#' with the caller's (caller wins per group key), and returns the pieces every
#' downstream consumer needs. This is the single source of truth for what a
#' selection means; expand_groups() and select_features() are thin views over
#' it.
#'
#' @param groups Character vector of feature / group / supergroup / bundle
#'   ids, or NULL for every exploratory feature.
#' @param pick   Named list: caller's per-feature keep-lists.
#' @param drop   Named list: caller's per-feature drop-lists.
#' @return list(groups = <feature ids>, pick = <named list>,
#'   drop = <named list>). The element keeps the name `groups` so existing
#'   consumers read it unchanged.
#' @export
resolve_selection <- function(groups = NULL, pick = list(), drop = list()) {
  bundles <- .group_bundles()

  if (is.null(groups)) {
    return(list(groups = EXPLORATORY_FEATURES, pick = pick, drop = drop))
  }

  out_groups <- character()
  bun_pick   <- list()
  bun_drop   <- list()

  # `seen` guards against bundle self / mutual reference.
  walk <- function(tokens, seen) {
    for (g in tokens) {
      if (g %in% names(SUPERGROUPS)) {
        out_groups <<- c(out_groups, intersect(SUPERGROUPS[[g]], EXPLORATORY_FEATURES))
      } else if (g %in% names(FEATURE_GROUPS)) {
        out_groups <<- c(out_groups, intersect(FEATURE_GROUPS[[g]], EXPLORATORY_FEATURES))
      } else if (g %in% names(bundles)) {
        if (g %in% seen) {
          warning("Bundle '", g, "' is self-referential — cycle broken")
          next
        }
        b <- .as_bundle(bundles[[g]])
        # Bundle pick/drop accumulate; later bundles override earlier ones for
        # the same group key (caller still overrides all of them, below).
        bun_pick[names(b$pick)] <<- b$pick
        bun_drop[names(b$drop)] <<- b$drop
        walk(b$groups, c(seen, g))
      } else if (g %in% names(FEATURE_PATTERNS)) {
        out_groups <<- c(out_groups, g)
      } else {
        warning("Unknown feature, group, supergroup or bundle: '", g, "' — skipped")
      }
    }
  }
  walk(groups, character())

  # Caller wins per group key; bundle entries fill the gaps.
  merged_pick <- utils::modifyList(bun_pick, pick)
  merged_drop <- utils::modifyList(bun_drop, drop)

  list(groups = unique(out_groups), pick = merged_pick, drop = merged_drop)
}


#' Expand selection keys into feature ids
#'
#' The group-key-only view of resolve_selection(); pick/drop carried by any
#' named bundles are resolved but not returned (use resolve_selection() or
#' select_features() if you need them).
#'
#' @param groups Character vector, or NULL for every exploratory feature.
#' @return Character vector of feature ids.
#' @examples
#' expand_groups()                                      # every exploratory feature
#' expand_groups("structure")                           # exploratory structure features
#' expand_groups(c("structure", "junction_abundance"))  # mixed
#' @export
expand_groups <- function(groups = NULL) {
  resolve_selection(groups)$groups
}


#' Apply pick/drop refinement to one group's columns (shared semantics)
#'
#' pick = allow-list, honouring caller order; drop = remove from the whole
#' group. Order: pick then drop. Used by both select_features() and the
#' correlation dotplot so the two never diverge.
#'
#' @param members Character vector of columns for a group, in df order.
#' @param pick_g  Character vector or NULL.
#' @param drop_g  Character vector or NULL.
#' @return Refined character vector.
#' @export
refine_group_columns <- function(members, pick_g = NULL, drop_g = NULL) {
  if (!is.null(pick_g)) members <- pick_g[pick_g %in% members]   # caller order
  if (!is.null(drop_g)) members <- setdiff(members, drop_g)      # members order
  members
}


# --- Discovery helpers -------------------------------------------------------

#' Identify which namespace a selection key belongs to.
#'
#' Returns "supergroup", "group", "bundle", "feature", or "unknown". Useful
#' interactively when you have a string and aren't sure which it is.
#'
#' @param key Character scalar — a token you want to use as a selection key.
#' @return Character scalar.
#' @examples
#' lookup_key("structure")       # "supergroup"
#' lookup_key("global_folding")  # "group"
#' lookup_key("nmd_core")        # "bundle"
#' lookup_key("rnafold_zscores") # "feature"
#' lookup_key("typo")            # "unknown"
#' @export
lookup_key <- function(key) {
  stopifnot(is.character(key), length(key) == 1)
  bundles <- .group_bundles()
  if (exists("SUPERGROUPS", inherits = TRUE) && key %in% names(SUPERGROUPS)) {
    return("supergroup")
  }
  if (key %in% names(FEATURE_GROUPS))   return("group")
  if (key %in% names(bundles))          return("bundle")
  if (key %in% names(FEATURE_PATTERNS)) return("feature")
  "unknown"
}


#' List every known selection key with its namespace and display name.
#'
#' Prints (and invisibly returns) a data.frame of all supergroups, groups,
#' bundles and features. Call this interactively to browse what you can pass
#' to `select_features()`, `fg()`, or any plot's `groups =` argument.
#'
#' @param kind Character vector. Which namespace(s) to show. Any combination
#'   of "supergroup", "group", "bundle", "feature". Default shows all four.
#' @param verbose Logical. If TRUE (default) print a formatted table.
#' @return Invisibly, a data.frame with columns `key`, `kind`, `display`.
#' @examples
#' list_selection_keys()                       # everything
#' list_selection_keys(kind = "feature")       # only feature ids
#' list_selection_keys(kind = c("supergroup", "bundle"))
#' @export
list_selection_keys <- function(kind = c("supergroup", "group", "bundle", "feature"),
                                verbose = TRUE) {
  kind <- match.arg(kind, several.ok = TRUE)
  rows <- list()

  if ("supergroup" %in% kind) {
    sgs <- if (exists("SUPERGROUPS", inherits = TRUE)) names(SUPERGROUPS) else character()
    for (k in sgs) {
      members <- paste(SUPERGROUPS[[k]], collapse = ", ")
      display <- format_group_name(k, "supergroup")
      rows[[length(rows) + 1]] <- data.frame(
        key = k, kind = "supergroup", display = display,
        members = members, stringsAsFactors = FALSE
      )
    }
  }

  if ("bundle" %in% kind) {
    bundles <- .group_bundles()
    for (k in names(bundles)) {
      b       <- .as_bundle(bundles[[k]])
      display <- format_group_name(k, "bundle")
      members <- paste(b$groups, collapse = ", ")
      if (length(b$pick) > 0) members <- paste0(members, " [pick]")
      if (length(b$drop) > 0) members <- paste0(members, " [drop]")
      rows[[length(rows) + 1]] <- data.frame(
        key = k, kind = "bundle", display = display,
        members = members, stringsAsFactors = FALSE
      )
    }
  }

  if ("group" %in% kind) {
    for (k in names(FEATURE_GROUPS)) {
      rows[[length(rows) + 1]] <- data.frame(
        key = k, kind = "group", display = format_group_name(k, "group"),
        members = paste(FEATURE_GROUPS[[k]], collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  }

  if ("feature" %in% kind) {
    for (k in names(FEATURE_PATTERNS)) {
      rows[[length(rows) + 1]] <- data.frame(
        key = k, kind = "feature", display = format_group_name(k, "feature"),
        members = paste0(group_of(k),
                         if (k %in% EXPLORATORY_FEATURES) "" else " [not exploratory]",
                         if (k %in% MODEL_FEATURES) " [model]" else ""),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  if (verbose && nrow(out) > 0) {
    # Print grouped by kind, with aligned columns.
    for (k in intersect(c("supergroup", "group", "bundle", "feature"), unique(out$kind))) {
      sub <- out[out$kind == k, , drop = FALSE]
      cat(sprintf("\n--- %ss ---\n", k))
      fmt <- paste0("  %-25s  %-28s  %s\n")
      cat(sprintf(fmt, "key", "display", if (k == "feature") "group" else "members"))
      cat(sprintf(fmt,
                  strrep("-", 25), strrep("-", 28), strrep("-", 20)))
      for (i in seq_len(nrow(sub))) {
        cat(sprintf(fmt, sub$key[i], sub$display[i], sub$members[i]))
      }
    }
    cat("\n")
  }

  invisible(out)
}


#' Columns of a dataset that the feature table marks as used nowhere
#'
#' The columns matching a row Excluded from both the exploratory analysis and
#' the model (NEVER_USED_PATTERNS, derived from R/feature_table.csv): Vienna
#' auxiliary statistics, engineering scaffolding, duplicates and superseded
#' encodings. The table's Notes column says why for each row.
#'
#' @param df A dataframe.
#' @return Character vector of column names, in `df` column order.
#' @export
excluded_columns <- function(df) {
  hit <- Reduce(`|`, lapply(NEVER_USED_PATTERNS, grepl, x = names(df)),
                logical(ncol(df)))
  names(df)[hit]
}


#' Remove the never-used columns from a built dataset
#'
#' Call it once, immediately after build_dataset() / build_all(), before
#' screening or modelling. Never inside build_dataset(): several of these
#' columns are engineering scaffolding (mfe_expected_* feeds mfe_delta_*,
#' eej_dist_{up,down}stream_* feed eej_dist_closest_*), and the QC scripts
#' legitimately read them. The cache stays complete, so changing the table's
#' flags never needs a CACHE_VERSION bump.
#'
#' @param df      Dataframe from build_dataset() / build_all().
#' @param exclude Character vector of column names. Defaults to
#'   excluded_columns(df); pass your own to screen a different pool.
#' @param verbose Logical. If TRUE (default) report how many columns went.
#' @return `df` without the excluded columns.
#' @examples
#' df <- build_dataset("human") |> drop_excluded()
#' # keep the Vienna auxiliary stats for one investigation:
#' df <- build_dataset("human")
#' df <- drop_excluded(df, exclude = grep("_(median|pval)_", excluded_columns(df),
#'                                        value = TRUE, invert = TRUE))
#' @export
drop_excluded <- function(df, exclude = excluded_columns(df), verbose = TRUE) {
  if (verbose) {
    message("drop_excluded: removed ", length(intersect(exclude, names(df))),
            " never-used columns; ",
            ncol(df) - length(intersect(exclude, names(df))), " remain")
  }
  dplyr::select(df, -dplyr::any_of(exclude))
}


#' The model's predictor columns, from the table's "Included in model" flag
#'
#' @param df       A dataframe.
#' @param features Feature ids to draw from; default every model feature.
#'   Intersected with MODEL_FEATURES, so passing a supergroup's members yields
#'   that supergroup's model block.
#' @return Character vector of column names, ordered by table row then `df`.
#' @examples
#' model_columns(df)                                           # every predictor
#' model_columns(df, SUPERGROUPS$structure)                    # structure block
#' @export
model_columns <- function(df, features = MODEL_FEATURES) {
  ids <- intersect(MODEL_FEATURES, features)
  unique(unlist(lapply(ids, function(g) fg_columns(df, g)), use.names = FALSE))
}


#' Is a feature a single region-less column (e.g. cai, translation efficiency)?
#'
#' Region-aware plots place such columns in the `mrna` slot. True when the
#' feature's regex is one literal column name whose last token is not a region.
#' @param feature Character vector of feature ids.
#' @return Logical vector.
#' @export
is_regionless_feature <- function(feature) {
  vapply(feature, function(g) {
    p <- FEATURE_PATTERNS[[g]]
    if (is.null(p) || !grepl("^\\^[a-z0-9_]+\\$$", p)) return(FALSE)
    toks <- strsplit(gsub("[$^]", "", p), "_", fixed = TRUE)[[1]]
    !toks[length(toks)] %in% REGIONS
  }, logical(1), USE.NAMES = FALSE)
}


#' Resolve a feature selection into an ordered vector of column names
#'
#' The single entry point for "which columns does this analysis use". Accepts
#' groups / supergroups / bundles plus optional per-group pick (allow-list) and
#' drop (remove). See file header for the bundle data model and merge policy.
#'
#' Missing columns are silently skipped (loaders drop missing data — R5).
#' pick/drop names that match no column are reported via message().
#'
#' @param df     Dataframe from build_dataset() / build_all().
#' @param groups Character vector of group / supergroup / bundle names, or NULL.
#' @param pick   Named list: group key -> columns to keep.
#' @param drop   Named list: group key -> columns to drop.
#' @return Ordered, de-duplicated character vector of column names in `df`.
#' @export
select_features <- function(df, groups = NULL,
                            pick = list(), drop = list()) {
  sel      <- resolve_selection(groups, pick, drop)
  expanded <- sel$groups

  unknown_pick <- setdiff(names(sel$pick), expanded)
  unknown_drop <- setdiff(names(sel$drop), expanded)
  if (length(unknown_pick) > 0) {
    message("select_features: pick names not in resolved groups (ignored): ",
            paste(unknown_pick, collapse = ", "))
  }
  if (length(unknown_drop) > 0) {
    message("select_features: drop names not in resolved groups (ignored): ",
            paste(unknown_drop, collapse = ", "))
  }

  cols <- character()
  for (g in expanded) {
    members <- fg_columns(df, g)

    if (!is.null(sel$pick[[g]])) {
      missing <- setdiff(sel$pick[[g]], members)
      if (length(missing) > 0) {
        message("select_features: pick[['", g, "']] columns not found: ",
                paste(missing, collapse = ", "))
      }
    }
    if (!is.null(sel$drop[[g]])) {
      missing <- setdiff(sel$drop[[g]], members)
      if (length(missing) > 0) {
        message("select_features: drop[['", g, "']] columns not found: ",
                paste(missing, collapse = ", "))
      }
    }

    cols <- c(cols, refine_group_columns(members, sel$pick[[g]], sel$drop[[g]]))
  }

  unique(cols)
}


# =============================================================================
# select_features_v2(): include / exclude / top_n  (SELECTION_PLAN.md, step 3)
# =============================================================================
# The replacement for select_features() + resolve_selection() + bundles.
# Lives beside the old API until the scripts are migrated; step 6 deletes the
# old one and renames this select_features().
#
#   result = expand(include) minus expand(exclude), then top_n trims families.
#
# Nothing else applies silently: no bundle pick/drop, no default skips.
#
# Tokens (include and exclude alike):
#   "core", "exploratory", "model"   the table's flags (CORE_FEATURES, ...)
#   a supergroup or group id         its members; in `include` only the
#                                    EXPLORATORY ones (the table's eligibility
#                                    rule), in `exclude` all of them
#   a feature id                     that feature, always, whatever its flags
# =============================================================================

.flag_sets <- function() {
  list(core = CORE_FEATURES, exploratory = EXPLORATORY_FEATURES,
       model = MODEL_FEATURES)
}

# Expand selection tokens to feature ids, in table order.
.expand_tokens <- function(tokens, eligible_only, arg) {
  flags <- .flag_sets()
  ids <- character()
  for (tk in tokens) {
    members <- if (tk %in% names(flags)) {
      flags[[tk]]
    } else if (tk %in% names(SUPERGROUPS)) {
      if (eligible_only) intersect(SUPERGROUPS[[tk]], EXPLORATORY_FEATURES) else SUPERGROUPS[[tk]]
    } else if (tk %in% names(FEATURE_GROUPS)) {
      if (eligible_only) intersect(FEATURE_GROUPS[[tk]], EXPLORATORY_FEATURES) else FEATURE_GROUPS[[tk]]
    } else if (tk %in% names(FEATURE_PATTERNS)) {
      tk
    } else {
      stop("select_features: unknown ", arg, " '", tk, "'. Use \"core\", ",
           "\"exploratory\", \"model\", or a supergroup / group / feature id ",
           "(see list_selection_keys()).", call. = FALSE)
    }
    ids <- c(ids, members)
  }
  intersect(names(FEATURE_PATTERNS), ids)   # table order, unique
}

# Metric stem of a column: the name without its region token. Columns with no
# region token are their own stem.
.column_stem <- function(cols) {
  # (no feature context here, so a region-less column is its own stem)
  last <- sub("^.*_", "", cols)
  has_region <- last %in% REGIONS & grepl("_", cols, fixed = TRUE)
  ifelse(has_region, sub("_[^_]+$", "", cols), cols)
}


#' Region token of each column, as the region-dodged plots place it
#'
#' A column ending in a real REGIONS token has that region. A column of a
#' region-less feature (cai, translation efficiency) sits in the "mrna" slot.
#' Any other column has no region (NA) and cannot be drawn on a region axis.
#'
#' @param cols       Character vector of column names.
#' @param feature_id Feature id of each column (same length).
#' @return Character vector of regions, NA where there is none.
#' @export
column_regions <- function(cols, feature_id) {
  last <- sub("^.*_", "", cols)
  has_region <- last %in% REGIONS & grepl("_", cols, fixed = TRUE)
  regionless <- vapply(feature_id, is_regionless_feature, logical(1),
                       USE.NAMES = FALSE)
  ifelse(has_region, last, ifelse(regionless, "mrna", NA_character_))
}


#' Feature ids selected by include / exclude
#'
#' The id-level half of select_features_v2(): the same expansion, before any
#' column is looked up. For consumers that show features with no columns in the
#' data too (a coverage tile reading "no cols").
#'
#' @inheritParams select_features_v2
#' @return Character vector of feature ids, in table order.
#' @export
selected_features <- function(include = "core", exclude = NULL) {
  if (length(include) == 0) stop("select_features: `include` is empty", call. = FALSE)
  ids <- .expand_tokens(include, eligible_only = TRUE, arg = "include")
  if (length(exclude)) {
    ids <- setdiff(ids, .expand_tokens(exclude, eligible_only = FALSE, arg = "exclude"))
  }
  ids
}


#' Select feature columns: include, exclude, trim families
#'
#' @param df       Dataframe from build_dataset() / build_all().
#' @param include  Tokens to start from (see above). Default "core".
#' @param exclude  Tokens to subtract. NULL for none.
#' @param top_n    Named list, feature id -> N: keep only the N metric stems of
#'   that feature with the largest |correlation| with `response` (the maximum
#'   over regions, and over species when `df` holds several), with all their
#'   regions. Ranking is over what `include` minus `exclude` left, so a
#'   feature trimmed here is ranked on the whole family unless you exclude
#'   part of it. Needs `response`.
#' @param regions  Character vector of region tokens to keep ("5utr", "cds", ...),
#'   or NULL for all. Applied before `top_n`, so a family is ranked only on the
#'   regions asked for. Region-less features count as "mrna"; columns with no
#'   region at all are dropped whenever `regions` is given.
#' @param response Column to rank against for `top_n`.
#' @param method   Correlation for `top_n`; default "spearman".
#' @param min_n    Minimum complete pairs for a column to be ranked; columns
#'   below it rank last. Default 30, as in the correlation figures.
#' @return A data.frame with columns `column` and `feature_id`, in table order
#'   then `df` column order. Columns absent from `df` are skipped silently.
#'   Use selected_columns() for the plain character vector. When `top_n` is
#'   given, attribute "top_n" is a list, one entry per trimmed feature, with
#'   `requested`, `available` (stems before trimming) and `bound` (did the
#'   trim remove anything).
#' @examples
#' select_features_v2(df)                                    # the core set
#' select_features_v2(df, exclude = "sequence")              # core minus Sequence
#' select_features_v2(df, "exploratory", exclude = c("codon_freqs", "aa_freqs"))
#' select_features_v2(df, top_n = list(codon_freqs = 2), response = "halflife")
#' @export
select_features_v2 <- function(df, include = "core", exclude = NULL,
                               top_n = NULL, response = NULL, regions = NULL,
                               method = "spearman", min_n = 30) {
  ids <- selected_features(include, exclude)

  out <- do.call(rbind, lapply(ids, function(g) {
    cols <- fg_columns(df, g)
    if (length(cols)) data.frame(column = cols, feature_id = g, stringsAsFactors = FALSE)
  }))
  if (is.null(out)) out <- data.frame(column = character(), feature_id = character(),
                                      stringsAsFactors = FALSE)

  if (!is.null(regions) && nrow(out)) {
    bad <- setdiff(regions, REGIONS)
    if (length(bad)) stop("select_features: unknown region(s): ",
                          paste(bad, collapse = ", "), call. = FALSE)
    out <- out[column_regions(out$column, out$feature_id) %in% regions, , drop = FALSE]
  }

  trimmed <- list()
  if (length(top_n)) {
    if (is.null(response) || !response %in% names(df))
      stop("select_features: `top_n` ranks against `response`, which must be a ",
           "column of df", call. = FALSE)
    bad <- setdiff(names(top_n), names(FEATURE_PATTERNS))
    if (length(bad)) stop("select_features: `top_n` names unknown feature(s): ",
                          paste(bad, collapse = ", "), call. = FALSE)
    groups <- if ("species" %in% names(df)) split(seq_len(nrow(df)), df$species)
              else list(seq_len(nrow(df)))
    abs_r <- function(co) {
      max(vapply(groups, function(ix) {
        x <- df[[co]][ix]; y <- df[[response]][ix]
        ok <- stats::complete.cases(x, y)
        if (sum(ok) < min_n || stats::sd(x[ok]) == 0) return(NA_real_)
        abs(stats::cor(x[ok], y[ok], method = method))
      }, numeric(1)), na.rm = TRUE)   # -Inf (with a warning) if all NA
    }
    keep <- rep(TRUE, nrow(out))
    for (g in names(top_n)) {
      rows <- which(out$feature_id == g)
      if (!length(rows)) next
      stem <- .column_stem(out$column[rows])
      r <- suppressWarnings(vapply(out$column[rows], abs_r, numeric(1)))
      by_stem <- tapply(r, stem, max)
      by_stem <- by_stem[order(-by_stem, names(by_stem))]   # ties: name order
      keep_stems <- names(by_stem)[seq_len(min(top_n[[g]], length(by_stem)))]
      keep[rows] <- stem %in% keep_stems
      trimmed[[g]] <- list(requested = top_n[[g]], available = length(by_stem),
                           bound = length(by_stem) > top_n[[g]])
    }
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  if (length(top_n)) attr(out, "top_n") <- trimmed
  out
}


#' The plain column-name view of a select_features_v2() result
#' @param sel Result of select_features_v2().
#' @return Character vector of column names.
#' @export
selected_columns <- function(sel) sel$column
