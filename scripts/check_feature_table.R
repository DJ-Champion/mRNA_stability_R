# =============================================================================
# Check R/feature_table.csv against a built cache
# =============================================================================
# The feature table is the single source of truth for every feature. This
# script verifies it against the built data, so the table and the pipeline
# cannot drift apart silently. Exits non-zero if any check fails.
#
# Usage (from the project root):
#   Rscript scripts/check_feature_table.R            # ANALYSIS_SPECIES
#   Rscript scripts/check_feature_table.R human
#
# The checks config.R already makes when it reads the table (flag values,
# unique feature ids, each group in one supergroup, no id at two levels) are
# not repeated here — a table that fails them does not load.
# =============================================================================

suppressMessages(source("R/load_all.R"))

args    <- commandArgs(trailingOnly = TRUE)
species <- if (length(args)) args else ANALYSIS_SPECIES

failures <- character()
fail <- function(...) failures <<- c(failures, paste0(...))
note <- function(...) cat("  note: ", ..., "\n", sep = "")

tbl  <- FEATURE_TABLE
rows <- tbl[nzchar(tbl$columns), ]

region_of <- function(col) {
  toks <- strsplit(col, "_", fixed = TRUE)[[1]]
  last <- toks[length(toks)]
  if (length(toks) > 1 && last %in% REGIONS) last else NA_character_
}
region_token <- stats::setNames(names(REGION_DISPLAYS), REGION_DISPLAYS)


# --- Data-independent checks --------------------------------------------------
cat("Table: ", nrow(tbl), " rows, ", length(FEATURE_PATTERNS), " selectable features, ",
    length(MODEL_FEATURES), " model features\n", sep = "")

dup_names <- tbl$Feature[duplicated(tbl$Feature)]
if (length(dup_names)) fail("Feature names not unique (they label legends): ",
                            paste(unique(dup_names), collapse = ", "))

for (id in EXPLORATORY_FEATURES) {
  col <- tbl$Colour[tbl$feature_id == id]
  if (!grepl("^#[0-9A-Fa-f]{6}$", col)) fail("exploratory feature '", id,
                                              "' has no valid Colour (got '", col, "')")
}
cols_used <- tbl$Colour[tbl$feature_id %in% EXPLORATORY_FEATURES]
if (anyDuplicated(cols_used)) fail("Colour reused across exploratory features: ",
                                   paste(unique(cols_used[duplicated(cols_used)]), collapse = ", "))

not_expl <- setdiff(MODEL_FEATURES, EXPLORATORY_FEATURES)
if (length(not_expl)) note("in the model but not the exploratory analysis: ",
                           paste(not_expl, collapse = ", "))

for (b in names(GROUP_BUNDLES)) {
  bb <- .as_bundle(GROUP_BUNDLES[[b]])
  bad <- bb$groups[vapply(bb$groups, lookup_key, character(1)) == "unknown"]
  if (length(bad)) fail("bundle '", b, "' names unknown key(s): ", paste(bad, collapse = ", "))
  bad <- setdiff(c(names(bb$pick), names(bb$drop)), names(FEATURE_PATTERNS))
  if (length(bad)) fail("bundle '", b, "' pick/drop keyed by non-feature(s): ",
                        paste(bad, collapse = ", "))
  if (b %in% c(tbl$feature_id, tbl$group_id, tbl$supergroup_id))
    fail("bundle name '", b, "' collides with a table id")
}
bad <- DEFAULT_PLOT_GROUPS[vapply(DEFAULT_PLOT_GROUPS, lookup_key, character(1)) == "unknown"]
if (length(bad)) fail("DEFAULT_PLOT_GROUPS names unknown key(s): ", paste(bad, collapse = ", "))

unbuilt <- tbl$feature_id[!nzchar(tbl$columns)]
if (length(unbuilt)) note("not yet built (no columns): ", paste(unbuilt, collapse = ", "))


# --- Checks against each species' cache -----------------------------------------
for (sp in species) {
  cat("\n== ", sp, " ==\n", sep = "")
  df <- suppressMessages(build_dataset(sp, min_utr = NULL))
  cols <- setdiff(names(df), setdiff(META_COLS, BENCHMARK_COLS))   # benchmarks have rows

  # 1. Every column claimed by exactly one row.
  hits <- vapply(rows$columns, function(p) grepl(p, cols), logical(length(cols)))
  if (!is.matrix(hits)) hits <- matrix(hits, nrow = length(cols))
  n_hit <- rowSums(hits)
  if (any(n_hit == 0)) fail(sp, ": column(s) claimed by no table row: ",
                            paste(cols[n_hit == 0], collapse = ", "))
  if (any(n_hit > 1)) fail(sp, ": column(s) claimed by more than one row: ",
                           paste(cols[n_hit > 1], collapse = ", "))

  # 2. Every row with a regex matches something.
  empty <- rows$feature_id[colSums(hits) == 0]
  if (length(empty)) fail(sp, ": row(s) whose columns regex matches nothing: ",
                          paste(empty, collapse = ", "))

  # 3. Regions column agrees with the columns' region tokens.
  for (i in seq_len(nrow(rows))) {
    rc <- cols[hits[, i]]
    have <- unique(stats::na.omit(vapply(rc, region_of, character(1))))
    if (!length(have)) next                      # region-less columns (cai, te)
    listed <- region_token[trimws(strsplit(rows$Regions[i], ",")[[1]])]
    if (anyNA(listed) || !setequal(have, listed))
      fail(sp, ": '", rows$feature_id[i], "' Regions lists {", rows$Regions[i],
           "} but its columns carry {", paste(REGION_DISPLAYS[have], collapse = ", "), "}")
  }

  # 4. Labels match the table's short name. Rows whose columns differ by more
  #    than their region cannot share one label: they must at least be distinct.
  esc <- function(s) gsub("([][.+?()|^$\\\\{}*])", "\\\\\\1", s, perl = TRUE)
  for (i in seq_len(nrow(rows))) {
    rc <- cols[hits[, i]]
    if (!length(rc)) next
    short <- rows$`Display Short name`[i]
    labs  <- format_col_name(rc)
    stems <- unique(sub("_[^_]+$", "", rc[!is.na(vapply(rc, region_of, character(1)))]))
    stems <- c(stems, rc[is.na(vapply(rc, region_of, character(1)))])
    if (!grepl("<", short, fixed = TRUE) && length(unique(stems)) > 1) {
      if (anyDuplicated(labs)) fail(sp, ": '", rows$feature_id[i], "' gives two columns one label")
      next
    }
    parts <- strsplit(short, "<[^>]+>")[[1]]
    if (grepl("<[^>]+>$", short)) parts <- c(parts, "")
    rx  <- paste0("^", paste(esc(parts), collapse = "\\S+"), "( |$)")
    bad <- rc[!grepl(rx, labs, perl = TRUE)]
    if (length(bad)) fail(sp, ": '", rows$feature_id[i], "' short name \"", short,
                          "\" but format_col_name() gives \"", format_col_name(bad[1]), "\"")
  }

  # 5. The model never sees an identifier, blocking or benchmark column.
  leak <- intersect(model_columns(df), c(META_COLS, "halflife"))
  if (length(leak)) fail(sp, ": model columns include non-predictors: ",
                         paste(leak, collapse = ", "))
  cat("  ", length(cols), " columns checked; model uses ", length(model_columns(df)),
      "; exploratory features select ", length(select_features(df, EXPLORATORY_FEATURES)),
      "; never used: ", length(excluded_columns(df)), "\n", sep = "")

  # 6. Fallback label rules that no longer fire (dead code in naming.R).
  fallback <- Filter(function(r) !grepl("^\\[ _\\]", r[[1]]), REPLACEMENTS)
  dead <- vapply(fallback, function(r) !any(grepl(r[[1]], names(df))), logical(1))
  if (any(dead)) note("naming.R rule(s) matching no column: ",
                      paste(vapply(fallback[dead], `[[`, "", 1), collapse = ", "))
}


# --- Verdict -----------------------------------------------------------------
if (length(failures)) {
  cat("\nFAILED (", length(failures), "):\n", paste0("  - ", failures, collapse = "\n"),
      "\n", sep = "")
  quit(status = 1)
}
cat("\nOK: the feature table matches the pipeline.\n")
