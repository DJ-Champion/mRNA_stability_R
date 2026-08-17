# =============================================================================
# MFE thermodynamic model
# =============================================================================
# Parameters from in-house model
#
# Expected MFE = (a * GC^b + c) * length + d,  capped at 0,  GC as a fraction.
# =============================================================================


#' Calculate expected MFE for given GC content and sequence length
#'
#' Vectorised. GC is expected as a FRACTION (0-1), which is the scale the
#' underlying nls was fit on (see RNAexpected/A-DataInput.R, eq3) and the
#' scale the gc_content_* columns are stored in. Passing percent (0-100)
#' drives gc^b to ~0 and silently collapses the model to a constant
#' c = -0.1534 kcal/mol/nt with no GC dependence at all, so it is rejected.
#'
#' @param gc Numeric vector of GC content as a fraction in [0, 1].
#' @param length Numeric vector of sequence lengths (nt).
#' @return Numeric vector of expected MFE values. Capped at 0 (folding free
#'   energy cannot be positive).
#' @export
calculate_mfe_expected <- function(gc, length) {
  a <- -0.8403211
  b <-  2.3521348
  c <- -0.1534111
  d <- 13.5600933
  if (any(gc > 1, na.rm = TRUE)) {
    stop("gc must be a fraction in [0, 1], not percent — got max ",
         format(max(gc, na.rm = TRUE)))
  }
  expected <- (a * gc^b + c) * length + d
  pmin(expected, 0)
}


#' Calculate the thermodynamic delta between observed and expected MFE
#'
#' @param mfe Observed MFE score.
#' @param mfe_expected Expected MFE for the same (GC, length).
#' @return mfe - mfe_expected. Proportional to the log-odds of the folded
#'   state given sequence composition.
#' @export
calculate_mfe_delta <- function(mfe, mfe_expected) mfe - mfe_expected
