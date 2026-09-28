# =============================================================================
# Feature group palettes and shape mappings
# =============================================================================
# Feature colours and every group/supergroup display label come from
# R/feature_table.csv (the Colour, Feature, Group and Supergroup columns) —
# edit them there. Colours follow supergroup families: structure purples,
# magentas and oranges; sequence blues, greens and teals; translation
# lavenders; transcript architecture browns and reds; RNA decay rusts.
#
# Plus Okabe-Ito colourblind-friendly palette for REGIONS, used when region
# is the visual variable rather than the grouping (see feature_correlation_
# dotplot vs feature_response_scatter for the two patterns).
#
# v4: the v3 pseudo-region tokens (`transcript`, `window`, `core`, `full`)
# have been retired. Whole-transcript scalars and the single NMD model now
# use the real `mrna` region, so REGION_COLOURS/SHAPES/DISPLAYS cover only
# the eight real regions plus the `none` fallback.
#
# Adapted from the legacy plotting scripts. To override: pass `palette = ...`
# or `region_colours = ...` to plot functions. To extend: add keys here.
# =============================================================================


# -----------------------------------------------------------------------------
# Feature-group colours (used when GROUP is the visual variable)
# -----------------------------------------------------------------------------

#' Colour for each feature id (the table's Colour column), plus the `other`
#' catch-all. Features without a colour (never-used rows) fall back to it via
#' feature_colour(); scripts/check_feature_table.R requires a colour on every
#' exploratory feature.
#' @export
FEATURE_GROUP_COLOURS <- c(
  with(FEATURE_TABLE[nzchar(FEATURE_TABLE$Colour), ], stats::setNames(Colour, feature_id)),
  other = "#C7C7C7"
)


# -----------------------------------------------------------------------------
# Region colours and shapes (used when REGION is the visual variable)
# -----------------------------------------------------------------------------
# Okabe-Ito colourblind-friendly palette for the eight real regions. Order
# intentionally follows the 5' -> 3' traversal of the transcript so the
# legend reads naturally and intersect()-based ordering preserves it.

#' Colours for each REGIONS token plus a fallback.
#' @export
REGION_COLOURS <- c(
  `5utr`     = "#009E73",  # bluish green
  start      = "#006442",  # green
  cds        = "#56B4E9",  # sky blue
  mrna       = "#0072B2",  # blue
  stop       = "#E69F00",  # orange
  `3utr`     = "#D55E00",  # vermillion
  last100    = "#CC79A7",  # reddish purple
  utrpair    = "#000000",  # black
  none       = "#999999"   # grey (unrecognised / no region)
)


#' Shape codes for each REGIONS token plus a fallback.
#' @export
REGION_SHAPES <- c(
  `5utr`     = 17,   # filled triangle (up)
  start      = 2,    # open triangle (up)
  cds        = 15,   # filled square
  mrna       = 16,   # filled circle
  stop       = 6,    # open triangle (down)
  `3utr`     = 18,   # filled diamond
  last100    = 8,    # asterisk
  utrpair    = 11,   # crossed square
  none       = 4     # x (unrecognised / no region)
)


#' Display strings for each REGIONS token.
#' MUST stay in sync with the Region tokens table in PIPELINE_GUIDE §2.4 and
#' with the region rules in R/utils/naming.R. Used by format_metric_name()
#' to strip the region from a formatted column name.
#' @export
REGION_DISPLAYS <- c(
  `5utr`     = "5' UTR",
  cds        = "CDS",
  `3utr`     = "3' UTR",
  mrna       = "mRNA",
  utrpair    = "UTR interactions",
  last100    = "last 100 nt",
  start      = "start codon",
  stop       = "stop codon",
  none       = ""
)


# -----------------------------------------------------------------------------
# Display formatters
# -----------------------------------------------------------------------------

#' Display labels for selection keys, from the table: a feature is labelled by
#' its Feature column, a group by Group, a supergroup by Supergroup. `other` is
#' the display bucket plots use for anything outside a named supergroup.
FEATURE_GROUP_DISPLAY_NAMES <- c(
  stats::setNames(FEATURE_TABLE$Feature, FEATURE_TABLE$feature_id), other = "Other")
GROUP_DISPLAY_NAMES <- c(
  stats::setNames(FEATURE_TABLE$Group, FEATURE_TABLE$group_id)[!duplicated(FEATURE_TABLE$group_id)],
  other = "Other")
SUPERGROUP_DISPLAY_NAMES <- c(
  stats::setNames(FEATURE_TABLE$Supergroup, FEATURE_TABLE$supergroup_id)[!duplicated(FEATURE_TABLE$supergroup_id)],
  other = "Other")

#' Bundle name -> display label. MUST cover every GROUP_BUNDLES key, or the
#' toTitleCase fallback produces things like "Nmd Core".
BUNDLE_DISPLAY_NAMES <- c(
  nmd_core         = "NMD (core)",
  lengths_core     = "Core lengths",
  junction_core    = "Core junction proximity",
  structure_core   = "Core structure",
  sequence_core    = "Core sequence",
  sequence_select  = "Sequence (selected)",
  translation_core = "Core translation"
)


#' Display name for a selection key (feature / group / supergroup / bundle).
#'
#' Looks the key up in the table for its `kind`; on a miss, falls back to
#' title-cased underscore replacement. Vectorised.
#'
#' @param g    Character vector of selection keys.
#' @param kind One of "feature", "group", "supergroup", "bundle", or "auto".
#'   "auto" resolves each token with the same precedence as
#'   resolve_selection(): supergroup -> group -> bundle -> feature.
#' @return Character vector of display strings, same length as `g`.
#' @export
format_group_name <- function(g, kind = c("auto", "feature", "group",
                                          "supergroup", "bundle")) {
  kind <- match.arg(kind)
  bundles <- if (exists("GROUP_BUNDLES", inherits = TRUE)) names(GROUP_BUNDLES) else character()

  pick_table <- function(k) switch(k,
    feature    = FEATURE_GROUP_DISPLAY_NAMES,
    group      = GROUP_DISPLAY_NAMES,
    supergroup = SUPERGROUP_DISPLAY_NAMES,
    bundle     = BUNDLE_DISPLAY_NAMES
  )

  vapply(g, function(key) {
    this_kind <- kind
    if (this_kind == "auto") {
      this_kind <-
        if      (key %in% names(SUPERGROUPS))    "supergroup"
        else if (key %in% names(FEATURE_GROUPS)) "group"
        else if (key %in% bundles)               "bundle"
        else                                     "feature"
    }
    tbl <- pick_table(this_kind)
    if (key %in% names(tbl)) unname(tbl[[key]])
    else tools::toTitleCase(gsub("_", " ", key))
  }, character(1), USE.NAMES = FALSE)
}


#' Display name for a column with its region suffix stripped.
#'
#' Useful when region is encoded as a separate visual (colour, shape) and
#' shouldn't be repeated in the axis label. Strategy: format the full column
#' via format_col_name(), then strip the formatted region from the end.
#' Falls back to format_col_name() unchanged for columns without a region
#' suffix.
#'
#' @param col Character vector of column names.
#' @return Character vector of display strings.
#' @examples
#' format_metric_name("length_cds")                    # "Length"
#' format_metric_name("rnafold_zscore_5utr")            # "MFE z-score"
#' format_metric_name("cai")                            # "CAI"
#' format_metric_name("intron_length_mean_mrna")        # "Mean intron length"
#' format_metric_name("nmd_fragile_codon_count_mrna")   # "NMD fragile codon count"
#' @export
format_metric_name <- function(col) {
  vapply(col, function(co) {
    tokens <- strsplit(co, "_", fixed = TRUE)[[1]]
    if (length(tokens) <= 1) return(format_col_name(co))

    last <- tokens[length(tokens)]
    if (!last %in% REGIONS) return(format_col_name(co))

    full   <- format_col_name(co)
    suffix <- paste0(" ", REGION_DISPLAYS[[last]])

    # Use endsWith + substr to avoid regex escaping. If the formatted name
    # doesn't actually end in the region display string, return it unchanged.
    if (nchar(suffix) > 1 && endsWith(full, suffix)) {
      trimws(substr(full, 1, nchar(full) - nchar(suffix)))
    } else {
      full
    }
  }, character(1), USE.NAMES = FALSE)
}


# -----------------------------------------------------------------------------
# Convenience accessors
# -----------------------------------------------------------------------------

#' Get colours for a vector of feature ids, falling back to "other" grey.
#' @export
feature_colour <- function(group) {
  out <- FEATURE_GROUP_COLOURS[group]
  out[is.na(out)] <- FEATURE_GROUP_COLOURS[["other"]]
  unname(out)
}

#' Get colours for a vector of region tokens, falling back to "none" grey.
#' @export
region_colour <- function(region) {
  region <- ifelse(is.na(region) | region == "", "none", region)
  out <- REGION_COLOURS[region]
  out[is.na(out)] <- REGION_COLOURS[["none"]]
  unname(out)
}

#' Get shapes for a vector of region tokens, falling back to "none".
#' @export
region_shape <- function(region) {
  region <- ifelse(is.na(region) | region == "", "none", region)
  out <- REGION_SHAPES[region]
  out[is.na(out)] <- REGION_SHAPES[["none"]]
  unname(out)
}
