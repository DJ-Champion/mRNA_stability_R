# =============================================================================
# Feature selection: fg(), select_features(), and the helpers around them
# =============================================================================
# Two layers, kept deliberately separate:
#
#   1. SCHEMA      R/feature_table.csv, read by config.R into FEATURE_PATTERNS
#                  (regex per feature), FEATURE_GROUPS and SUPERGROUPS (the two
#                  coarser levels), and the three flags: exploratory, core,
#                  model (EXPLORATORY_FEATURES, CORE_FEATURES, MODEL_FEATURES).
#
#   2. INTENT      select_features(df, include, exclude, top_n) says which
#                  columns one analysis wants, in one vocabulary, identically
#                  everywhere:  expand(include) minus expand(exclude), then
#                  top_n trims families to their strongest stems. Nothing else
#                  applies silently: no named bundles, no pick/drop lists, no
#                  per-script default skips.
#
# Tokens (include and exclude alike):
#   "core", "exploratory", "model"   the table's flags
#   a supergroup or group id         its members; in `include` only the
#                                    EXPLORATORY ones (the table's eligibility
#                                    rule), in `exclude` all of them
#   a feature id                     that feature, always, whatever its flags
#
# Usage:
#   df %>% select(fg("rnafold_zscores"))                     # one feature
#   select_features(df)                                      # the core set
#   select_features(df, "structure")                         # a supergroup
#   select_features(df, "core", exclude = "sequence")        # core minus one
#   select_features(df, "exploratory", exclude = c("codon_freqs", "aa_freqs"))
#   select_features(df, top_n = list(codon_freqs = 2), response = "halflife")
#   selected_columns(select_features(df, "structure"))       # plain vector
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


# --- Discovery helpers -------------------------------------------------------

#' Identify which namespace a selection key belongs to.
#'
#' Returns "flag", "supergroup", "group", "feature", or "unknown". Useful
#' interactively when you have a string and aren't sure which it is.
#'
#' @param key Character scalar — a token you want to use as a selection key.
#' @return Character scalar.
#' @examples
#' lookup_key("core")            # "flag"
#' lookup_key("structure")       # "supergroup"
#' lookup_key("global_folding")  # "group"
#' lookup_key("rnafold_zscores") # "feature"
#' lookup_key("typo")            # "unknown"
#' @export
lookup_key <- function(key) {
  stopifnot(is.character(key), length(key) == 1)
  if (key %in% names(.flag_sets()))     return("flag")
  if (key %in% names(SUPERGROUPS))      return("supergroup")
  if (key %in% names(FEATURE_GROUPS))   return("group")
  if (key %in% names(FEATURE_PATTERNS)) return("feature")
  "unknown"
}


#' List every known selection key with its namespace and display name.
#'
#' Prints (and invisibly returns) a data.frame of all flags, supergroups,
#' groups and features. Call this interactively to browse what you can pass to
#' `include` / `exclude` in `select_features()` or any plot.
#'
#' @param kind Character vector. Which namespace(s) to show. Any combination
#'   of "flag", "supergroup", "group", "feature". Default shows all four.
#' @param verbose Logical. If TRUE (default) print a formatted table.
#' @return Invisibly, a data.frame with columns `key`, `kind`, `display`,
#'   `members`.
#' @examples
#' list_selection_keys()                       # everything
#' list_selection_keys(kind = "feature")       # only feature ids
#' list_selection_keys(kind = c("flag", "supergroup"))
#' @export
list_selection_keys <- function(kind = c("flag", "supergroup", "group", "feature"),
                                verbose = TRUE) {
  kind <- match.arg(kind, several.ok = TRUE)
  row <- function(key, kind, display, members)
    data.frame(key = key, kind = kind, display = display, members = members,
               stringsAsFactors = FALSE)
  rows <- list()

  if ("flag" %in% kind) {
    for (k in names(.flag_sets())) {
      rows[[length(rows) + 1]] <- row(k, "flag", k,
                                      paste(.flag_sets()[[k]], collapse = ", "))
    }
  }
  if ("supergroup" %in% kind) {
    for (k in names(SUPERGROUPS)) {
      rows[[length(rows) + 1]] <- row(k, "supergroup",
                                      format_group_name(k, "supergroup"),
                                      paste(SUPERGROUPS[[k]], collapse = ", "))
    }
  }
  if ("group" %in% kind) {
    for (k in names(FEATURE_GROUPS)) {
      rows[[length(rows) + 1]] <- row(k, "group", format_group_name(k, "group"),
                                      paste(FEATURE_GROUPS[[k]], collapse = ", "))
    }
  }
  if ("feature" %in% kind) {
    for (k in names(FEATURE_PATTERNS)) {
      flags <- c(if (k %in% CORE_FEATURES) "core",
                 if (k %in% EXPLORATORY_FEATURES) "exploratory",
                 if (k %in% MODEL_FEATURES) "model")
      rows[[length(rows) + 1]] <- row(
        k, "feature", format_group_name(k, "feature"),
        paste0(group_of(k), if (length(flags)) paste0(" [", paste(flags, collapse = ", "), "]")))
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  if (verbose && nrow(out) > 0) {
    for (k in intersect(c("flag", "supergroup", "group", "feature"), unique(out$kind))) {
      sub <- out[out$kind == k, , drop = FALSE]
      cat(sprintf("\n--- %ss ---\n", k))
      fmt <- paste0("  %-25s  %-28s  %s\n")
      cat(sprintf(fmt, "key", "display", if (k == "feature") "group [flags]" else "members"))
      cat(sprintf(fmt, strrep("-", 25), strrep("-", 28), strrep("-", 20)))
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


#' Is a feature a single region-less column (e.g. translation efficiency)?
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


# =============================================================================
# select_features(): include / exclude / top_n
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
#' region-less feature (translation efficiency) sits in the "mrna" slot.
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
#' The id-level half of select_features(): the same expansion, before any
#' column is looked up. For consumers that show features with no columns in the
#' data too (a coverage tile reading "no cols").
#'
#' @inheritParams select_features
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
#' select_features(df)                                    # the core set
#' select_features(df, exclude = "sequence")              # core minus Sequence
#' select_features(df, "exploratory", exclude = c("codon_freqs", "aa_freqs"))
#' select_features(df, top_n = list(codon_freqs = 2), response = "halflife")
#' @export
select_features <- function(df, include = "core", exclude = NULL,
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


#' The plain column-name view of a select_features() result
#' @param sel Result of select_features().
#' @return Character vector of column names.
#' @export
selected_columns <- function(sel) sel$column
