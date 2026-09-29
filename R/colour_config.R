# =============================================================================
# Colour configuration
# =============================================================================
# Project-wide colours that are a matter of taste rather than of the feature
# schema. Edit the values here; every plot that reads them follows.
#
# Feature colours live in R/feature_table.csv (Colour column). Region colours
# and shapes live in R/utils/palettes.R. This file holds the SUPERGROUP colours:
# one hue per supergroup, used wherever a whole supergroup needs a colour (the
# background bands in feature_correlation_bands.R, and anything that wants to
# say "which family is this" without colouring every feature).
#
# Keys are supergroup ids (the table's Supergroup column, snake-cased). Every
# supergroup with built features (the names of SUPERGROUPS) must have an entry,
# or loading fails, so adding a supergroup to the table cannot silently leave
# it uncoloured. Supergroups with no columns yet (expression) may have one
# ready for when they do.
# =============================================================================

SUPERGROUP_COLOURS <- c(
  structure               = "#2D004B",  # deep purple
  sequence                = "#002642",  # deep navy
  translation             = "#9E9AC8",  # lavender
  transcript_architecture = "#DFC27D",  # sand
  rna_decay               = "#B15928",  # rust
  expression              = "#1B7837"   # green
)

# Colour for anything outside a named supergroup (the "other" bucket).
SUPERGROUP_OTHER_COLOUR <- "#999999"

# Background bands: the share of the supergroup colour kept when it is mixed
# with white (0 = white, 1 = full colour), and the same for the label bar.
SUPERGROUP_BAND_ALPHA  <- 0.18
SUPERGROUP_LABEL_ALPHA <- 0.75


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
