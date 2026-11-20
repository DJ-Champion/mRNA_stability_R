# =============================================================================
# MFE thermodynamic model
# =============================================================================
# Expected MFE comes from the RNAexpected package
# (github.com/DJ-Champion/RNAexpectedMFE), the single source of truth for the
# coefficients. This file only wraps it so callers keep one stable name.
#
# Expected MFE = (a * GC^b + c) * length + d,  capped at 0,  GC as a fraction.
# The average MFE of a random sequence of the same GC and length: a
# composition-matched baseline, not a dinucleotide-matched one (that is
# rnafold_median_*, from the per-transcript shuffles).
#
# Known limits of the model: it still under-predicts structure below ~100 nt
# and returns 0 below ~30-60 nt (GC-dependent), and it was fitted on
# 20-4000 nt, so longer regions are extrapolations.
# =============================================================================

#' RNAexpected model used for every mfe_expected_* column.
MFE_EXPECTED_MODEL <- "C"


#' Calculate expected MFE for given GC content and sequence length
#'
#' Thin wrapper around RNAexpected::expected_mfe(). Vectorised; NA in either
#' argument gives NA. GC must be a FRACTION (0-1): RNAexpected stops on values
#' outside [0, 1], since passing percent would drive gc^b to ~0 and silently
#' collapse the model to a GC-free constant.
#'
#' @param gc Numeric vector of GC content as a fraction in [0, 1].
#' @param length Numeric vector of sequence lengths (nt).
#' @return Numeric vector of expected MFE values (kcal/mol), capped at 0.
#' @export
calculate_mfe_expected <- function(gc, length) {
  if (!requireNamespace("RNAexpected", quietly = TRUE)) {
    stop("calculate_mfe_expected() needs the RNAexpected package:\n",
         "  remotes::install_github(\"DJ-Champion/RNAexpectedMFE\")",
         call. = FALSE)
  }
  RNAexpected::expected_mfe(gc, length, model = MFE_EXPECTED_MODEL)
}


#' Provenance of the expected-MFE model, for stamping on the cache
#'
#' There is no lockfile, so this is what records which coefficients a cached
#' dataset's mfe_expected_* / mfe_delta_* columns were built with.
#' @return List: model, coefficients, package version and (if installed from
#'   GitHub) commit SHA.
mfe_model_provenance <- function() {
  calculate_mfe_expected(0.5, 100)  # fails early with the install hint
  desc <- utils::packageDescription("RNAexpected")
  list(model        = MFE_EXPECTED_MODEL,
       coefficients = RNAexpected::expected_mfe_coefficients(MFE_EXPECTED_MODEL),
       version      = desc$Version,
       sha          = desc$RemoteSha)
}


#' Calculate the thermodynamic delta between observed and expected MFE
#'
#' @param mfe Observed MFE score.
#' @param mfe_expected Expected MFE for the same (GC, length).
#' @return mfe - mfe_expected. Proportional to the log-odds of the folded
#'   state given sequence composition.
#' @export
calculate_mfe_delta <- function(mfe, mfe_expected) mfe - mfe_expected
