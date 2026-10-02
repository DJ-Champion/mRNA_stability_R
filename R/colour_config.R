# =============================================================================
# Colour configuration
# =============================================================================
# Project-wide colours that are a matter of taste rather than of the feature
# schema. Edit the values here; every plot that reads them follows.
#
# Region colours and shapes live in R/utils/palettes.R. This file holds the
# SUPERGROUP colours: one hue per supergroup, used wherever a whole supergroup
# needs a colour (the background bands in feature_correlation_bands.R, and
# anything that wants to say "which family is this" without colouring every
# feature), and the FEATURE colours, which are derived from them: each feature
# is a saturated colour in its supergroup's hue family (see feature_colours()
# below), so the two levels always read as one system. Change a supergroup colour and its features
# follow.
#
# Keys are supergroup ids (the table's Supergroup column, snake-cased). Every
# supergroup with built features (the names of SUPERGROUPS) must have an entry,
# or loading fails, so adding a supergroup to the table cannot silently leave
# it uncoloured. Supergroups with no columns yet (expression) may have one
# ready for when they do.
# =============================================================================

SUPERGROUP_COLOURS <- c(
  structure               = "#FF0000",  
  sequence                = "#00FFFF",  
  translation             = "#FF00FF",  
  transcript_architecture = "#00FF00",  
  rna_decay               = "#0000FF",  
  expression              = "#FFFF00"   
)

# Colour for anything outside a named supergroup (the "other" bucket).
SUPERGROUP_OTHER_COLOUR <- "#999999"

# Background bands: the share of the supergroup colour kept when it is mixed
# with white (0 = white, 1 = full colour), and the same for the label bar. The
# supergroup colours are vivid (they are the base of the feature colours, which
# have to stand out on white), so these mixes are what keep the bands subtle.
SUPERGROUP_BAND_ALPHA  <- 0.10
SUPERGROUP_LABEL_ALPHA <- 0.40

# Feature colours: each feature stays in its supergroup's hue family, at a
# higher saturation than the pale band tints. Within a supergroup the features
# (table order) are spread evenly in lightness (FEATURE_LUMINANCE, dark ->
# light) and in hue, over SUPERGROUP_HUE_SPREAD degrees centred on the
# supergroup hue (0 = one hue). Bigger supergroups get bigger spreads so their
# features stay apart. Chroma is the colourfulness, shared by every feature.
FEATURE_CHROMA    <- 100
FEATURE_LUMINANCE <- c(30, 68)
SUPERGROUP_HUE_SPREAD <- c(
  structure               = 40,
  sequence                = 40,
  translation             = 14,
  transcript_architecture = 16,
  rna_decay               = 22,
  expression              = 0
)


# --- Validation --------------------------------------------------------------
local({
  is_hex <- function(x) grepl("^#[0-9A-Fa-f]{6}$", x)
  missing <- setdiff(names(SUPERGROUPS), names(SUPERGROUP_COLOURS))
  if (length(missing)) {
    stop("R/colour_config.R: no SUPERGROUP_COLOURS entry for: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(names(SUPERGROUP_COLOURS), FEATURE_TABLE$supergroup_id)
  if (length(unknown)) {
    stop("R/colour_config.R: SUPERGROUP_COLOURS names not in the feature table: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  bad <- c(names(SUPERGROUP_COLOURS)[!is_hex(SUPERGROUP_COLOURS)],
           if (!is_hex(SUPERGROUP_OTHER_COLOUR)) "SUPERGROUP_OTHER_COLOUR")
  if (length(bad)) {
    stop("R/colour_config.R: not a #RRGGBB colour: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
})


#' Colour for each supergroup id; `other` and unknown ids get the fallback.
#' @param supergroups Character vector of supergroup ids.
#' @export
supergroup_colour <- function(supergroups) {
  out <- unname(SUPERGROUP_COLOURS[supergroups])
  out[is.na(out)] <- SUPERGROUP_OTHER_COLOUR
  out
}


# --- Feature colours (vivid versions of the supergroup hue) ------------------

#' Colour for every exploratory feature
#'
#' Within a supergroup the features, in table order, are spaced evenly in
#' lightness (FEATURE_LUMINANCE) and in hue (over SUPERGROUP_HUE_SPREAD degrees
#' centred on the supergroup hue) at a shared chroma (FEATURE_CHROMA). A
#' supergroup with one feature gets the middle of both ranges.
#' @return Named character vector keyed by feature id.
#' @export
feature_colours <- function() {
  ids <- EXPLORATORY_FEATURES
  sg  <- FEATURE_TABLE$supergroup_id[match(ids, FEATURE_TABLE$feature_id)]
  out <- stats::setNames(rep(SUPERGROUP_OTHER_COLOUR, length(ids)), ids)
  for (g in unique(sg[sg %in% names(SUPERGROUP_COLOURS)])) {
    in_g <- ids[sg == g]
    n    <- length(in_g)
    t    <- if (n == 1) 0.5 else seq(0, 1, length.out = n)
    hue0 <- farver::decode_colour(SUPERGROUP_COLOURS[[g]], to = "hcl")[1, "h"]
    spread <- if (g %in% names(SUPERGROUP_HUE_SPREAD)) SUPERGROUP_HUE_SPREAD[[g]] else 0
    out[in_g] <- grDevices::hcl(
      h = (hue0 + spread * (t - 0.5)) %% 360,
      c = FEATURE_CHROMA,
      l = FEATURE_LUMINANCE[1] + diff(FEATURE_LUMINANCE) * t,
      fixup = TRUE
    )
  }
  out
}

FEATURE_COLOURS <- feature_colours()
