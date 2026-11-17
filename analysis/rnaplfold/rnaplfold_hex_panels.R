# =============================================================================
# RNAplfold accessibility vs both responses (hex panels)
# =============================================================================
# Thin runner over feature_response_hex_panels() (analysis/correlations/):
# every RNAplfold column against half-life and against TE, with the existing
# start-region MFE and 5' UTR GC alongside for scale.
#
# Usage: Rscript analysis/rnaplfold/rnaplfold_hex_panels.R
# =============================================================================
source("R/load_all.R")
suppressPackageStartupMessages(library(ggplot2))

# The hex-panel file runs its own jobs when sourced at top level; a private
# environment (where the runner's guard is false) gives us the function only.
.hex_env <- new.env(parent = globalenv())
sys.source("analysis/correlations/feature_response_hex_panels.R", envir = .hex_env)

OUT <- file.path(OUTPUT_DIR, "rnaplfold_exploration")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

df <- build_dataset("human") |> drop_excluded(verbose = FALSE)

for (response in c("halflife", "translation_efficiency")) {
  out <- .hex_env$feature_response_hex_panels(
    df, response = response, include = "local_accessibility",
    columns = c("rnafold_score_start", "gc_content_5utr"))
  n <- sum(out$table$plotted)
  ggsave(file.path(OUT, paste0("hex_panels_", response, ".jpg")), out$plot,
         width = 270, height = 190, units = "mm", dpi = 200, limitsize = FALSE)
  write.csv(out$table, file.path(OUT, paste0("hex_panels_", response, ".csv")),
            row.names = FALSE)
  message(response, ": ", n, " panels")
}
