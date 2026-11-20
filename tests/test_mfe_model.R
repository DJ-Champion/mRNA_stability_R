# Tests for the expected-MFE wrapper (R/features/mfe_model.R). Run from the
# project root:
#   Rscript tests/test_mfe_model.R
# Needs RNAexpected but not the cache.

suppressMessages(source("R/load_all.R"))

n_fail <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  if (!ok) n_fail <<- n_fail + 1L
  cat(if (ok) "  ok    " else "  FAIL  ", desc, "\n", sep = "")
}
errors <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

k  <- RNAexpected::expected_mfe_coefficients(MFE_EXPECTED_MODEL)
gc <- c(0.3, 0.5, 0.7, 0.5)
L  <- c(80, 500, 3000, 25)   # in the fitted range; 25 nt still caps at 0

cat("calculate_mfe_expected\n")
check("matches the package formula",
      isTRUE(all.equal(calculate_mfe_expected(gc, L),
                       pmin((k[["a"]] * gc^k[["b"]] + k[["c"]]) * L + k[["d"]], 0))))
check("identical to RNAexpected::expected_mfe",
      identical(calculate_mfe_expected(gc, L),
                RNAexpected::expected_mfe(gc, L, model = MFE_EXPECTED_MODEL)))
check("capped at 0 for short sequences", calculate_mfe_expected(0.5, 25) == 0)
check("NA in, NA out",
      identical(is.na(calculate_mfe_expected(c(0.5, NA, 0.5), c(100, 100, NA))),
                c(FALSE, TRUE, TRUE)))
check("percent GC errors", errors(calculate_mfe_expected(50, 100)))

cat("provenance\n")
prov <- mfe_model_provenance()
check("records the model", identical(prov$model, MFE_EXPECTED_MODEL))
check("records the coefficients", identical(prov$coefficients, k))
check("check_mfe_model warns on a cache with no stamp",
      inherits(tryCatch(check_mfe_model(data.frame()), warning = function(w) w), "warning"))
stamped <- structure(data.frame(), mfe_model = prov)
check("check_mfe_model is quiet on a current stamp",
      is.null(tryCatch(check_mfe_model(stamped), warning = function(w) w)))

cat(if (n_fail) sprintf("\n%d FAILED\n", n_fail) else "\nall passed\n")
quit(status = as.integer(n_fail > 0))
