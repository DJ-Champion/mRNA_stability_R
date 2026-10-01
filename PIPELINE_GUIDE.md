# RNAstab — Pipeline Extension Guide

**Audience:** a contributor (or AI coding instance) tasked with writing new plotting scripts, refactoring legacy scripts onto this pipeline, or extending the pipeline with new features, loaders, or species.

**Authority:** this document is a rulebook, not a tour. Where it says **MUST** / **MUST NOT** / **DO** / **DO NOT**, treat those as constraints, not suggestions. The README explains *what the pipeline is*. This document explains *how to interface with and extend it without breaking the invariants*.

**Companion reading (in this order):** `README.md` → this document → `R/config.R` → `R/io/load_raw.R` → `R/pipeline/build_dataset.R`. Do not write code before reading all four.

---

## 1. Mental model

### 1.1 Data flow

```
data/raw/<species>/         ──┐
data/raw/shared/            ──┤  load_raw.R (one loader per source file)
                              │   returns NULL on missing files (silent)
                              ▼
                      regional (long)  +  transcript-level (wide)
                              │
                              ▼  pipeline/assemble.R
                       pivot_regional_to_wide()
                              │
                              ▼  join_transcript_level / join_gene_level
                       wide tibble, one row per transcript
                              │
                              ▼  features/engineer.R
                       engineer_features()  (derived columns)
                              │
                              ▼  io/cache.R
                  data/cache/<species>_dataset_v<N>.rds
                              │
                              ▼
                  consumer code: build_dataset(species)
```

Blocked splits are a **second, independent artefact**, built from `family.tsv` alone and joined on at use:

```
data/raw/shared/family.tsv  ──▶  load_family()  ──┬──▶  build_dataset()   [family_* columns]
      (the clustering seam)                       │
                                                  └──▶  build_splits()    [pipeline/splits.R]
                                                             │  run ONCE
                                                             ▼
                                              data/splits/holdout_<level>.rds
                                                             │
                                                             ▼
                                    consumer code: attach_splits(build_dataset(...))
```

`build_splits()` reads the seam directly rather than the cache, so it needs no cache to exist and cannot be perturbed by a feature-engineering change. See §6.7.

### 1.2 Shape invariants

After `build_dataset(species)` returns, the dataframe satisfies **all** of these. Any consumer code MAY rely on them; any extender MUST preserve them:

- **One row per transcript.** `transcript_id` is unique within a species.
- **`species` column is present** and equals the species name on every row.
- **Identifier columns come first:** `species`, `transcript_id`, `gene_id`, `gene_name` (whichever exist).
- **Column names are canonical** (lowercase, snake_case, no species prefix). See §2.
- **Region suffix is the last token** of any regional column: `length_5utr`, `gc_content_cds`, `rnafold_zscore_mrna`. *Aspirational, not yet universal* — a handful of loader columns still violate this (`utr5_length`, `n_exons`, `internal_exon_mean`, …) and are consequently invisible to `fg()` and dropped by region-aware plots. See §10. New code MUST NOT add to that set.
- **All-NA columns are dropped** by `drop_all_na_columns()` at the end of `engineer_features()`. DO NOT rely on a specific column existing without guarding.
- **`halflife` is the canonical response variable.** No species prefix, no suffix.

### 1.3 Cross-species idiom

```r
combined <- bind_rows(build_dataset("human"), build_dataset("mouse"))
# or
combined <- build_all()
combined |> group_by(species) |> summarise(...)
```

Any plot or analysis that supports multi-species **MUST** facet, group, or filter on the `species` column. DO NOT bake species into column names anywhere.

---

## 2. Canonical schema reference

### 2.1 Region vocabulary

These are the only legal region suffixes. Defined in `REGIONS` in `R/config.R`:

| Suffix     | Meaning                                |
|------------|----------------------------------------|
| `5utr`     | 5' untranslated region                 |
| `cds`      | Coding sequence                        |
| `3utr`     | 3' untranslated region                 |
| `mrna`     | Whole mRNA                             |
| `utrpair`  | 5' UTR + 3' UTR (interaction features) |
| `last100`  | Last 100 nt of CDS                     |
| `start`    | Start codon region                     |
| `stop`     | Stop codon region                      |

DO NOT invent new suffixes. If a feature genuinely needs a new region, add it to `REGIONS` first and document it in this table.

### 2.2 The feature table — the single source of truth

**`R/feature_table.csv` defines every feature.** One row per feature. Each row gives:

- `feature_id` and `columns`, the regex that claims that feature's columns;
- its place in the **Supergroup > Group > Feature** hierarchy;
- the three flags, **Included in exploratory analysis**, **Included in model** and **Included in core plots**;
- its display names;
- its Regions, Description and Notes (why it is in or out).

`R/config.R` reads the table and derives `FEATURE_PATTERNS`, `FEATURE_GROUPS`, `SUPERGROUPS`, `EXPLORATORY_FEATURES`, `CORE_FEATURES`, `MODEL_FEATURES` and `NEVER_USED_PATTERNS`. `R/utils/palettes.R` takes display names from it, and `R/colour_config.R` derives each feature's colour from its supergroup hue, and `R/utils/naming.R` takes labels from it. **Edit the table, never those objects.** Then run:

```bash
Rscript scripts/check_feature_table.R
```

It fails if any of these is broken:

- a cache column is claimed by no row, or by two;
- a row's regex matches nothing;
- a row's Regions disagree with its columns;
- a label differs from the table's short name;
- two exploratory features share a colour;
- a core feature is not an exploratory one, or a flag name collides with a table id;
- a retired selection name (`GROUP_BUNDLES`, `DEFAULT_PLOT_GROUPS`, `resolve_selection`, `pick`/`drop` machinery, …) reappears in `R/`, `analysis/` or `scripts/`;
- the model would see an identifier or benchmark column.

**Three levels, one flat selection namespace.** A selection key is any of the following, and ids must not collide across levels (config.R refuses a table where they do):

- a **feature id** (a table row);
- a **group id**, the table's Group snake-cased: "Global folding" → `global_folding`;
- a **supergroup id**: "Transcript architecture" → `transcript_architecture`;
- a **flag** name: `core`, `exploratory` or `model`.

Browse them with `list_selection_keys()`.

**What the flags do.**

| Flag | Effect |
|---|---|
| Included in exploratory analysis | Expanding a group or supergroup yields only its exploratory features. Naming a feature id directly always yields it, so a plot can reach an excluded feature deliberately. |
| Included in core plots | `"core"`: the default `include` for the correlation figures. A subset of the exploratory features; probing is deliberately outside it. |
| Included in model | `model_columns(df)`: the model's predictor list. The structure model's blocks are its Structure and non-Structure parts. |
| Excluded from both | Used nowhere. `drop_excluded()` removes these columns (`excluded_columns(df)` lists them). |

Rows in the "Response / evaluation" supergroup (`halflife`, `saluki`) and rows with no `columns` (a feature not yet built) are documented but are not features: they cannot be selected.

**Two schema invariants,** both enforced:

1. **The patterns are mutually exclusive.** No column may match two rows. Plots that build a column → feature map (the dotplot, the response scatter, the heatmap workflow) would otherwise draw a column twice. The checker enforces this.
2. **Every group sits in exactly one supergroup.** Enforced when config.R reads the table.

> **Region-less features.** `translation_efficiency` is a single column with no region suffix (CAI is `cai_cds`: it is a CDS measure and carries the token). `is_regionless_feature()` recognises them from their literal regex, and region-aware plots place them in the `mrna` slot.

> **Whole-transcript scalars carry the real `mrna` suffix.** Architecture (`intron_length_mean_mrna`, …), uORF (`uorf_count_mrna`, …) and NMD fragility (`nmd_snv_fragile_codon_density_mrna`, …) all end in `mrna`. The old `transcript` pseudo-region and the `window`/`core`/`full` NMD pseudo-regions are gone; there is one NMD fragility model and "NMD" appears only as a metric-name prefix in the display string.

### 2.3 Display labels

Every column that appears on a plot **MUST** have a human-readable `format_col_name()` result.

**Where labels come from.** Most rows are either one literal stem followed by a region token (`^mfe_delta_` + `cds`) or one literal column (`^translation_efficiency$`). For these, the table's **Display Short name** is the label, plus the region's display string: `mfe_delta_cds` → `MFE.Δ CDS`. **To change a label, edit the table.**

**The exceptions.** A few rows cannot share one label across their columns. Their labels come from `format_single_name()` and a short `REPLACEMENTS` list in `R/utils/naming.R`:

| Row | Why one label cannot cover it |
|---|---|
| codon and amino-acid frequencies | The codon or residue is uppercased into the label: `codon.AAA%`, `aa.L%` |
| nucleotide fractions | One label per base: `nt.A%` |
| icSHAPE | One label per compartment: `icSHAPE.nuc`, `icSHAPE.cyto` |
| multi-column rows (uORF count, intron length, …) | The columns differ by more than their region |

The checker verifies these against the short name too. Placeholders such as `<codon>` match any token.

The house style is compact and dotted (`MFE.z`, `C+G%`, `nt.A%`, `icSHAPE.cyto`) rather than prose. The labels are sized for dense small-multiples panels.

> **PDF output.** The base `pdf()` device cannot draw non-Latin-1 characters. That affects `Δ` in `MFE.Δ` and `→` in two never-used labels, as well as the `ρ` and em dashes already in some titles. Save figures that carry them with `device = cairo_pdf`, or as PNG/JPG.

#### Region-suffix rules

For fallback labels only, these fire **last**, after every prefix rule has
consumed the leading portion of the column name. Each matches a single
leading separator (space *or* underscore) plus a region token, end-anchored —
`[ _]<region>$`. Table-driven labels append `REGION_DISPLAYS` directly.

#### Region tokens

The eight tokens are the complete `REGIONS` vocabulary (§2.1). The table's
Regions column uses these display strings, and the checker compares it with
each row's columns.

| Token       | Display string         |
|-------------|------------------------|
| `5utr`      | `5' UTR`               |
| `3utr`      | `3' UTR`               |
| `cds`       | `CDS`                  |
| `mrna`      | `mRNA`                 |
| `utrpair`   | `UTR interactions`     |
| `last100`   | `last 100 nt`          |
| `start`     | `start codon`          |
| `stop`      | `stop codon`           |

**If a column has no table row,** the checker fails. Add a row rather than working around it in the plot.

---

## 3. The public interface — the only things you should call

These are the functions every extender code should rely on. DO NOT reach into internals (`load_*`, `pivot_regional_to_wide`, `engineer_features` directly, etc.) from analysis scripts.

| Function                          | File                                       | Purpose                                          |
|-----------------------------------|--------------------------------------------|--------------------------------------------------|
| `build_dataset(species, rebuild)` | `R/pipeline/build_dataset.R`               | Get the wide dataframe for one species           |
| `build_all(species, rebuild)`     | `R/pipeline/build_dataset.R`               | Stack multiple species, one `species` column     |
| `attach_splits(df, level)`        | `R/pipeline/splits.R`                      | Add the blocked `split` column from the on-disk artefact |
| `load_splits(level)`              | `R/pipeline/splits.R`                      | Read the split artefact (`NULL` if not built yet) |
| `build_splits(level, ...)`        | `R/pipeline/splits.R`                      | **Run once**, via `scripts/build_splits.R` — writes the artefact |
| `validate_splits(assigned)`       | `R/pipeline/splits.R`                      | Assert the blocking guarantee before modelling   |
| `META_COLS` / `ID_COLS` / `FAMILY_COLS` | constants in `R/config.R`            | Columns carried on every row but never predictors |
| `BLOCK_LEVEL`                     | constant in `R/config.R`                   | Which clustering level blocks the splits (`medium`) |
| `R/feature_table.csv`             | —                                          | **The feature definitions.** Everything below that describes features is derived from it |
| `scripts/check_feature_table.R`   | —                                          | Verify the table against a built cache; non-zero exit on any mismatch |
| `fg(feature)`                     | `R/utils/feature_groups.R`                 | Tidyselect spec for one feature                  |
| `fg_columns(df, feature)`         | `R/utils/feature_groups.R`                 | Inspect what a feature resolves to in this df    |
| `select_features(df, include, exclude, top_n, response, regions)` | `R/utils/feature_groups.R` | **The** column selector: `include` minus `exclude` (flags / supergroups / groups / features), families trimmed by `top_n`. Returns `column` + `feature_id` |
| `selected_columns(sel)` / `selected_features(include, exclude)` | `R/utils/feature_groups.R` | The plain column vector of a result / the feature ids before any column lookup |
| `model_columns(df, features)`             | `R/utils/feature_groups.R` | The model's predictor columns (table flag)                     |
| `excluded_columns(df)` / `drop_excluded(df)` | `R/utils/feature_groups.R` | List / remove the never-used columns                        |
| `lookup_key(key)`                         | `R/utils/feature_groups.R` | "flag" / "supergroup" / "group" / "feature" / "unknown"        |
| `list_selection_keys(kind)`               | `R/utils/feature_groups.R` | Browse every selection key with display names                  |
| `supergroup_of(f)` / `group_of(f)`        | `R/config.R`               | Reverse lookup: feature id → its supergroup / group            |
| `FEATURE_GROUPS` / `SUPERGROUPS`          | constants in `R/config.R`  | Group / supergroup id → feature ids (from the table)           |
| `EXPLORATORY_FEATURES` / `CORE_FEATURES` / `MODEL_FEATURES` | constants in `R/config.R`  | Feature ids carrying each flag (from the table)                |
| `format_col_name(x)`              | `R/utils/naming.R`                         | Canonical column name → display string (vectorised) |
| `format_group_name(x, kind)`      | `R/utils/palettes.R`                       | Feature / group / supergroup key → display string (vectorised) |
| `format_metric_name(x)`           | `R/utils/palettes.R`                       | Column → display string with the region suffix stripped |
| `feature_colour(group)` / `region_colour(r)` / `region_shape(r)` | `R/utils/palettes.R` | Palette accessors with a documented fallback |
| `clear_snapshot(species)`         | `R/io/cache.R`                             | Force next build to rebuild                      |
| `REGIONS`                         | constant in `R/config.R`                   | The legal region suffix vocabulary               |
| `FEATURE_PATTERNS`                | constant in `R/config.R`                   | Feature id → regex (from the table)              |
| `REGION_COLOURS` / `REGION_SHAPES` / `REGION_DISPLAYS` | constants in `R/utils/palettes.R` | Per-region visual vocabulary        |
| `SPECIES_CONFIG`                  | constant in `R/config.R`                   | Species registry                                 |
| `OUTPUT_DIR`                      | constant in `R/config.R`                   | `data/outputs` — base path for all outputs       |
| `CACHE_VERSION`                   | constant in `R/config.R`                   | Bump when feature-engineering logic changes      |

Calling `source("R/load_all.R")` once at the top of any script loads all of these.

---

## 4. Hard rules

The following are absolute. Violating them will silently corrupt data, break caching, or make plots unrenderable.

### R1 — Always source `R/load_all.R` first

Every script under `analysis/`, `scripts/`, or anywhere else MUST begin with:

```r
source("R/load_all.R")
```

DO NOT `source()` individual pipeline files. DO NOT redefine pipeline functions locally.

**Corollary — nothing may sit under `R/` except live pipeline code.** `load_all.R` sources *every* `.R` file in `R/utils/`, `R/io/`, `R/features/` and `R/pipeline/` via `list.files()`. A backup, a scratch file, or a `foo_old.R` left in place is loaded too, and if it redefines a function it may win or lose depending on **locale collation order** — `list.files()` sorts differently under `LC_COLLATE=C` (typical for `Rscript` on a server or in CI) than under a UTF-8 desktop locale. This is not hypothetical: `R/utils/naming_old.R` shadowed `format_col_name()` this way and was only harmless by luck of sort order. Park work-in-progress outside `R/`, or on a branch.

### R2 — Always get data via `build_dataset()`

DO NOT read `.rds` cache files directly. DO NOT call loaders. DO NOT bypass the engineering step. The only legal entry points are `build_dataset(species)` and `build_all()`.

### R3 — Select features with `fg()` / `select_features()`, never hand-rolled regex

**DO:**
```r
df |> select(halflife, fg("rnafold_zscores"), fg("rnalfold_zscores"))
```

**DO NOT:**
```r
df |> select(halflife, matches("^rnafold_zscore_"), starts_with("rnal"))
```

If the feature you need has no row in `R/feature_table.csv`, **add one** (see §6.4). DO NOT inline a regex.

### R3a — Select *subsets* through the selection layer, never new schema groups

`fg()` (R3) selects a whole feature. When you need **less than a whole
feature** — a supergroup minus one family, the top few codons, everything but
probing — that is *selection intent*, and it MUST be expressed through
`select_features()`, not by adding a narrower row to the feature table.

The feature table is **schema**: one row per real column family, each in
exactly one group and supergroup. DO NOT add subset or alias rows (the deleted
`*_some`, `mfe_scores`, `mfe_zscores` keys were exactly this mistake — a subset
masquerading as a family). A subset row double-claims columns, which the
checker rejects.

There is exactly one selection mechanism, with one vocabulary:

```r
# select_features(df, include = "core", exclude = NULL, top_n = NULL, response = NULL, regions = NULL)
select_features(df)                                        # the core set
select_features(df, "structure")                           # a supergroup
select_features(df, "core", exclude = "sequence")          # core minus a supergroup
select_features(df, "exploratory", exclude = c("codon_freqs", "aa_freqs"))
select_features(df, top_n = list(codon_freqs = 2, aa_freqs = 2),
                response = "halflife")                     # trim families to their top stems
selected_columns(select_features(df, "global_folding"))    # plain character vector
```

- **`include`** / **`exclude`** take flags (`"core"`, `"exploratory"`, `"model"`),
  supergroup ids, group ids and feature ids. Result = expand(`include`) minus
  expand(`exclude`); nothing else applies silently. In `include`, a group or
  supergroup expands to its *exploratory* members; a feature id always works.
- **`top_n`** is the only family trimmer: a named list, feature id → N, keeping
  the N metric stems with the largest |r| against `response` (max over regions
  and species). Two figures with the same response can therefore never disagree
  about which codons are "top". `top_n` never means anything else; an argument
  limiting how many things are *drawn* is called `max_features`.
- **`regions`** restricts to region tokens before `top_n` ranks.
- Do **not** hand-name columns to keep or drop, and do not add default skip
  lists inside an analysis: if a default excludes something (the sweep and the
  feature-feature table skip `codon_freqs` and `aa_freqs`), it is a visible
  `exclude = ` default in the function signature.

**Plot functions that choose features SHOULD take `include`, `exclude` and,
where ranking makes sense, `top_n` (and `regions`)** with exactly these
meanings, and pass them to `select_features()` (see §6.1, step 3). That keeps
selection identical across every plot.

**Labelling a selection key.** `format_col_name()` is for *column* names and
produces wrong output on selection keys. When a plot renders a feature, group
or supergroup **key** as visible text — a facet strip, a legend — pass it
through `format_group_name(key, kind)`, where `kind` is one of `"feature"`,
`"group"`, `"supergroup"`, or `"auto"` (supergroup → group → feature
precedence). Display strings come from the table's Feature, Group and
Supergroup columns. Edit those rather than hardcoding a label in the plot.

**Region-less features.** `te` (translation efficiency) is a single
column with no region suffix; `is_regionless_feature()` identifies them and
`column_regions()` maps them to the `mrna` slot for region-aware plots.

**Exception — single-feature tools.** A tool whose entire premise is "one panel
per schema family" (e.g. `feature_group_panel_sweep()`) makes one panel per
*feature id* that `include`/`exclude` resolve to. Its default excludes the
high-cardinality families `codon_freqs` and `aa_freqs` in the signature.

### R4 — Format every plot label through `format_col_name()`

Axis labels, legend titles, plot titles, facet strip labels, summary table column headers, and CSV column labels intended for human eyes **MUST** be passed through `format_col_name()` (or use a scale labeller that calls it: `scale_y_discrete(labels = format_col_name)`).

DO NOT hardcode display strings inside plot functions. If a column needs a better label, edit its Display Short name in `R/feature_table.csv` (or, for the few rows §2.3 lists, `R/utils/naming.R`); do not work around it in the plot.

### R5 — Guard every column access

Loaders return `NULL` silently when raw files are missing, and `engineer_features()` drops all-NA columns. A column you expect may not be present. Always guard:

```r
if (!"halflife" %in% names(df)) stop("halflife missing; this species lacks the response variable")
if ("mfe_delta_cds" %in% names(df)) { ... }
```

The single exception is `transcript_id` and `species`, which are pipeline invariants (§1.2).

### R6 — Never bake species into column names

DO NOT create columns like `human_halflife`, `mouse_length_cds`. Species belongs in the `species` column. If you need per-species comparisons, pivot or facet.

### R7 — Bump `CACHE_VERSION` when feature engineering changes

Edit `R/config.R`. Increment `CACHE_VERSION` integer. This invalidates every species' cache on next call to `build_dataset()`.

**When to bump:** any change to `R/features/engineer.R`, any change to a loader's output schema, any change to assembly logic.

**When NOT to bump:** new plotting script, new analysis script, anything under `analysis/`.

For a one-off rebuild without bumping the global version, use `build_dataset(species, rebuild = TRUE)` or `clear_snapshot(species)` first.

### R8 — Outputs go under `OUTPUT_DIR`, never anywhere else

Plots: `data/outputs/plots/<name>.<ext>`. Tables: `data/outputs/tables/<name>.csv`. Use `file.path(OUTPUT_DIR, "plots", ...)`, never a relative path or `getwd()`.

DO NOT create new top-level output directories. DO NOT save to the working directory.

### R9 — Plot functions return `list(plot, table)`

Any function that produces a plot derived from summary statistics (correlations, group means, lasso coefficients, top-N rankings) MUST return:

```r
list(plot = <ggplot>, table = <tibble of the underlying data>)
```

The `table` element is the source of truth that the plot visualises. Consumers can write it to CSV without re-running the computation. Pure scatter plots that visualise raw rows of `df` are an exception — they may return just the ggplot.

### R10 — Exclude derived response variables from response-correlation analyses

When computing correlations against `halflife`, you MUST exclude columns that are derived predictions of half-life — they are circular. The standard exclusion list is:

```r
exclude = c("^saluki_prediction$", "^prediction_difference$")
```

This is the default in `correlate_with_response()`. Any new model-output column added to the pipeline MUST be added to this default exclusion.

### R11 — Long-form region data uses `(transcript_id, region)` keys

Any new long-form regional loader MUST return a tibble with both `transcript_id` and a lowercase `region` column. The `normalise_region()` helper in `R/io/load_raw.R` enforces lowercase. The assembly step pivots on these two keys.

### R12 — Transcript-level wide loaders MUST NOT contain a `gene_id` column

`gene_id` is supplied by `load_transcripts()` only. Other loaders that happen to carry it MUST drop it (`select(-any_of("gene_id"))`) before returning. Failing to do so produces `gene_id.x` / `gene_id.y` suffix collisions on join.

### R13 — Family and split columns are metadata, never features

The `family_*` columns and `split` are carried on every row but MUST NOT enter a predictor matrix. They are listed in `META_COLS` (`R/config.R`); build a feature set as `setdiff(names(df), c(META_COLS, TARGET_COL))` and they are handled.

The trap is `family_size_medium`: it is numeric, sits among genuine features, and describes the *corpus* rather than the transcript. Nothing about its dtype or name stops a model consuming it.

Note they are deliberately absent from `FEATURE_PATTERNS`. Adding a `family` key there would make them reachable through `fg()` and `select_features()` — i.e. selectable *as features*, the exact opposite of the intent. `META_COLS` is the right register; a feature group is not.

### R14 — Generate the split once; read it everywhere

`build_splits()` is run from `scripts/build_splits.R` and writes `data/splits/holdout_<level>.rds`. Every consumer reads it via `attach_splits(df)`. NEVER re-derive the assignment inside an analysis or modelling script.

A re-derived split is not reproducible even with a fixed seed: a rebuilt `family.tsv`, a different R version, or a changed row order can move genes between train and test, and results stop being comparable with nothing looking wrong. `family.tsv` carries no record of the clustering run behind it, so the artefact stores that file's md5; `attach_splits()` warns if the dataset and the split were built from different ones.

---

## 5. Anti-patterns

Things that look reasonable and will silently break the pipeline or downstream analyses. If you find yourself writing any of these, stop.

| Anti-pattern                                                    | Why it breaks                                       | Correct form                                  |
|-----------------------------------------------------------------|-----------------------------------------------------|-----------------------------------------------|
| `readRDS("data/cache/human_dataset_v2.rds")`                    | Bypasses cache version mgmt, breaks on next bump    | `build_dataset("human")`                      |
| `select(df, matches("^rnafold_zscore_"))`                       | Group definition fragmented across files            | `select(df, fg("rnafold_zscores"))`           |
| `labs(x = "MFE z-score (CDS)")`                                 | Display layer drifts from canonical schema          | `labs(x = format_col_name("rnafold_zscore_cds"))` |
| `df$mfe_delta_cds + df$mfe_delta_3utr`                          | Column may not exist; no guard                      | Guard with `%in% names(df)` or `coalesce`     |
| Saving to `"plots/foo.png"`                                     | Lands in working dir, not `data/outputs/`           | `file.path(OUTPUT_DIR, "plots", "foo.png")`   |
| Hard-coding `c("5utr","cds","3utr")` in a loop                  | Misses `mrna`, `utrpair`, etc.                      | Iterate over `REGIONS`                        |
| Adding `species` filter without `if(species %in% ...)` check    | Silent empty plot when species missing              | `if (!any(df$species == "human")) stop(...)`  |
| Calling `engineer_features()` from an analysis script           | Skips cache, may double-engineer                    | `build_dataset()` does this for you            |
| Mutating column names with `rename_with(toupper)` for display   | Breaks `format_col_name()` round-trip               | Format only at the moment of display          |
| Renaming a column produced by a loader inside `engineer.R`      | Downstream `fg()` patterns break                    | Either rename in the loader, or add new col   |
| Adding a `nmd_core` row with regex `^nmd_(snv\|alt)` to the feature table | Subset masquerading as a schema family; double-claims columns, which the checker rejects | Express it in the call: `select_features(df, include = ..., exclude = ..., top_n = ...)` |
| `assign_holdout(fam)` inside a modelling script                 | Split silently reshuffles on any upstream change; results stop being reproducible | `attach_splits(df)` — read the artefact (R14) |
| `group_vfold_cv(df, group = gene_id)`                           | Blocks on the gene, not the family; paralogues still split across folds | `group = family_id_medium`, or block on `split`  |
| Judging a split by its gene counts alone                        | An 80/10/10 split can be exact while a held-out split holds only genes with no relatives | Check `pct_multi` is similar across splits in the `build_splits()` summary |
---

## 6. Extension recipes

Pick the recipe that matches your task. Follow every numbered step.

### 6.1 Adding a new plot or analysis script

**Location:**
- Correlation-style plot → `analysis/correlations/<name>.R`
- Model fit / coefficient plot → `analysis/models/<name>.R`
- QC / diagnostic plot → `analysis/qc/<name>.R`
- Anything comparing species to each other → `analysis/cross_species/<name>.R`
- Anything else → propose a new subdirectory in your PR description

**Steps:**

1. Create the file. Begin with `source("R/load_all.R")`. Load only the additional packages you need (`ggplot2`, `forcats`, `viridis`, etc.).
2. Define a single primary function. Signature pattern:
   ```r
   <verb>_<noun>_plot <- function(df, ..., formatter = format_col_name) { ... }
   ```
   The first argument MUST be the dataframe; the formatter MUST default to `format_col_name`.
3. Inside the function:
   - Use `fg()` for whole-feature selection, or `select_features(df, include,
     exclude, top_n, response)` when the plot exposes feature selection to its
     caller. Take `include = "core"` and `exclude = NULL` (and `top_n = NULL`
     where you rank) so behaviour matches every other plot (R3a).
   - Filter NA rows on the variables you actually plot.
   - Compute the summary table.
   - Build the ggplot using `format_col_name()` (or `formatter`) for every label.
4. Return `list(plot = p, table = <tibble>)`. (See R9 for the exception.)
5. Append a "top-to-bottom run" block at the bottom of the file, guarded so it only runs when the script is executed directly:
   ```r
   if (sys.nframe() == 0 || identical(environment(), globalenv())) {
     df  <- build_dataset("human")
     out <- <your_function>(df)
     print(out$plot)
     dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
     dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
     ggsave(file.path(OUTPUT_DIR, "plots", "<name>.jpg"),
            plot = out$plot, width = 210, height = 148, units = "mm", dpi = 300)
     write.csv(out$table, file.path(OUTPUT_DIR, "tables", "<name>.csv"), row.names = FALSE)
   }
   ```
6. DO NOT bump `CACHE_VERSION`. Analysis scripts never invalidate the cache.

### 6.2 Refactoring a legacy plotting script

Most legacy scripts will have one or more of these problems. Walk this checklist:

| Symptom                                       | Fix                                                  |
|-----------------------------------------------|------------------------------------------------------|
| Hard-coded column names like `human_halflife` | Strip species prefix → `halflife`. Use `species` column for filtering / faceting. |
| Local `format_col_name_v2()` definition       | Delete it. Use the pipeline's `format_col_name()`. Fix labels in `R/feature_table.csv`. |
| Reads CSV/RDS from disk directly              | Replace with `build_dataset("<species>")`.           |
| Long hand-rolled `select(matches(...))`       | Replace with `fg()` calls.                           |
| Imputes data inline                           | Move the imputation into `R/features/engineer.R` and bump `CACHE_VERSION`. Plots consume the engineered column. |
| Bare global parameters (`a <- 0.8`, etc.)     | Scope inside the function that uses them. Move thermodynamic constants to `R/features/mfe_model.R`. |
| Writes outputs to working directory           | Use `file.path(OUTPUT_DIR, ...)`.                    |
| Returns just a ggplot from a stats-driven plot| Wrap as `list(plot, table)` per R9.                  |
| Uses `mfe_` legacy alias instead of `rnafold_`| Use the canonical column name in code (`rnafold_score_*`); display will render correctly via `format_col_name()`. |

### 6.3 Adding a new derived feature

A "derived" feature is one computed from columns already in the assembled dataframe.

1. Open `R/features/engineer.R`.
2. Add a function `add_<feature_name>(df) -> df`. The function MUST guard on the presence of every input column it reads (R5). It MUST add columns, not mutate existing ones.
3. Add the call to the pipe inside `engineer_features()`. Order matters: place it after any dependencies and before `drop_all_na_columns()`.
4. Add a row for the new column(s) to `R/feature_table.csv` (§6.4). Every built column needs one; the checker fails otherwise.
5. **Bump `CACHE_VERSION`** in `R/config.R`.
6. Verify: `build_dataset("human", rebuild = TRUE)`, then `Rscript scripts/check_feature_table.R`.

### 6.4 Adding a new feature

One row in `R/feature_table.csv`. Nothing else needs editing: patterns, groups, supergroups, flags and labels are all derived from it; the colour follows the supergroup hue and is set in `R/colour_config.R`.

| Column | What to put |
|---|---|
| `feature_id` | Short snake_case id; unique, and not equal to any group or supergroup id |
| `columns` | Regex claiming exactly this feature's columns, mutually exclusive with every other row. Prefer a literal stem ending in `_` (`^my_metric_`): then the label is the short name plus the region. Leave empty for a feature not yet built. |
| `Included in exploratory analysis` / `Included in model` / `Included in core plots` | `Included` or `Excluded`. Core is the default plotting set and must be a subset of exploratory |
| `Supergroup`, `Group` | An existing pair, or a new one; a group sits in exactly one supergroup |
| `Feature`, `Display Short name`, `Display Long name` | The Feature name labels the feature in legends and must be unique; the short name is the plot label |
| `Regions` | The region display strings its columns carry (`5' UTR, CDS, …`) |
| `Description`, `Notes` | What it is; why it is in or out |

No cache bump needed (the table is a query layer, not data). **Verify:**

```bash
Rscript scripts/check_feature_table.R
```

### 6.4a Choosing a different set of columns for one analysis

There is no registry to edit. Say what the analysis wants in the call:
`include = ` / `exclude = ` (and `top_n = ` to trim a family). If the set should
be the *default* for routine plots, change the feature's **Included in core
plots** flag in the table and run `Rscript scripts/check_feature_table.R`. No
`CACHE_VERSION` bump is needed either way: selection is a query, not data.

### 6.5 Adding a new raw input source

1. Add a loader to `R/io/load_raw.R`. Required behaviours:
   - Takes a single `species` argument.
   - Calls `read_if_exists(species_path(species, "<filename>"))`.
   - Returns `NULL` if missing (the helper does this).
   - Returns canonical lowercase snake_case column names.
   - If regional/long-form: includes `transcript_id` and lowercase `region` (use `normalise_region()`).
   - If wide-form: includes `transcript_id` and drops `gene_id` (R12).
2. Wire it into `R/pipeline/build_dataset.R`:
   - Long-form → add to the `regional` list.
   - Transcript-level wide-form → add to the `transcript_level` list.
   - Gene-level → add a `load_*` call and an explicit `left_join` block with the right key.
3. Add a row per feature to `R/feature_table.csv` (§6.4).
4. Bump `CACHE_VERSION`, rebuild, and run `Rscript scripts/check_feature_table.R`.

### 6.6 Adding a new species

1. Append a block to `SPECIES_CONFIG` in `R/config.R`. Copy the `human` block, change `dir`, adjust `dani.shape_pattern` / `dani.keth_pattern` / `dani.extra_cols` to match the genome suffixes in the shared probing file.
2. Place raw files under `data/raw/<dir>/` matching the names the loaders expect (`grep -h "species_path" R/io/load_raw.R` lists every expected filename).
3. Add a runner: `scripts/build_<species>.R`. Copy `scripts/build_human.R` verbatim, change one string.
4. `Rscript scripts/build_<species>.R`. Watch the `skip (missing): ...` messages — they tell you which files are absent.

No code changes to any other file should be needed. If they are, you've found a leak — fix the leak rather than patching around it.

### 6.7 Family blocking — using it, and refreshing it

**Using it.** Two lines:

```r
df <- attach_splits(build_dataset("human"))
train <- df |> filter(split == "train", !is.na(halflife))
```

**The invariant, first, because the summary output reads as if it might be violated:** a family is *never* divided across splits. The packer's unit is the whole family — it sorts families by size and places each one, entire, into one split. Its only decision is which split. `validate_splits()` asserts this, and `build_splits()` refuses to write an artefact that fails.

So in the summary table, `max_family = 20` for `test` means *the largest family test received has 20 members, all 20 of them in test* — not that a family was cut. `pct_multi` is the share of a split's genes that have at least one relative, and that relative is necessarily in the same split. The three `pct_multi` figures should be close to one another: that is what makes the held-out splits resemble the training data. If one is near zero, that split holds only genes with no relatives at all — see the anti-pattern below.

`split` is `train` / `val` / `test`. `val` serves the purpose a nested inner CV loop would, so there is no inner-fold machinery to set up. For inference rather than prediction, the family label is the **clustering unit for uncertainty**, not just the split key — a 282-member family is not 282 independent observations, so `family_id_medium` should enter as a random effect or as the cluster for robust standard errors, or the confidence interval comes out too narrow and a null result cannot be trusted either.

**Two things to state in any write-up** that uses this split, because neither is visible in the numbers:

1. It tests **generalisation to novel gene families** — a stronger and different claim than random-over-genes, and different again from leave-one-species-out.
2. Families above 5% of the smallest split (68 genes at present) are **pinned to `train`**, so the held-out splits are depleted of the largest families.

**Sensitivity check.** Every level is a column, so re-blocking at `loose` is one flag: `Rscript scripts/build_splits.R --level loose`. "Conclusions unchanged under a looser grouping" is worth more than the level choice itself.

**Refreshing after a re-clustering.** When the Python side emits a new `family.tsv`:

1. Drop it in `data/raw/shared/family.tsv`.
2. Bump `CACHE_VERSION` and rebuild — the family columns are cached, so a stale cache would keep the old labels.
3. Re-run `Rscript scripts/build_splits.R`. **Gene membership of train/val/test will change.** Any model fitted against the old split is no longer comparable.

Skipping step 3 is caught, not silent: `attach_splits()` compares the md5 of the `family.tsv` behind the dataset against the one behind the split artefact and warns when they diverge.

**Adding a species to the blocking.** Family assignment comes from a *cohort* — the set of datasets clustered together — and merging orthologues across species is the whole point, since training on mouse `Rpl13a` and testing on human `RPL13A` is leakage. A species not in the cohort therefore gets no family columns and `NA` splits; `load_family()` returns `NULL` for it, which is correct behaviour, not a degradation. Making mouse blockable means adding it to the cohort **on the Python side** and re-clustering both species together — clustering mouse on its own would produce labels that silently fail to block cross-species leakage. `build_splits()` already pools every species the seam contains into one packing run.

---

## 7. Worked examples

Two worked examples covering different plot shapes. The first is correlation-shaped (every column in a group vs the response); the second is distribution-shaped (response distribution stratified by bins of a predictor). Together they exercise most of the public interface and rules.

### 7.1 Correlation panel — every column in a feature group vs `halflife`

**Task:** for a chosen feature group, produce a faceted scatter panel showing every column in the group plotted against `halflife`, with Spearman correlation annotated on each panel.

This example exercises: `build_dataset()`, `fg_columns()`, vectorised `format_col_name()`, the `list(plot, table)` return contract, the runner-block pattern, and `OUTPUT_DIR`.

File: `analysis/correlations/feature_group_panel.R`

```r
# =============================================================================
# Feature-group correlation panel
# =============================================================================
# For a chosen feature group, plot every member column against halflife as a
# small-multiples scatter panel, annotated with Spearman correlation.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_group_panel.R")
#   df  <- build_dataset("human")
#   out <- feature_group_panel(df, group = "rnafold_zscores")
#   print(out$plot)
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})


#' Faceted scatter panel of every column in a feature group vs halflife
#'
#' @param df      A dataframe from build_dataset().
#' @param group   Character. A feature id (a key of FEATURE_PATTERNS).
#' @param response Character. Response column (default "halflife").
#' @param formatter Function. Display formatter (default format_col_name).
#' @return list(plot, table). `table` has columns variable, n, spearman, p_value.
#' @export
feature_group_panel <- function(df,
                                group,
                                response  = "halflife",
                                formatter = format_col_name) {

  # --- Rule R5: guard every column access -----------------------------------
  if (!response %in% names(df)) {
    stop("response '", response, "' not in df — wrong species or missing data")
  }

  # --- Rule R3: use fg_columns to enumerate ---------------------------------
  cols <- fg_columns(df, group)
  if (length(cols) == 0) {
    stop("Feature group '", group, "' resolved to zero columns in this df")
  }

  # --- Long-format for faceting ---------------------------------------------
  long <- df |>
    dplyr::select(dplyr::all_of(c(response, cols))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(cols),
                        names_to = "variable",
                        values_to = "value") |>
    dplyr::filter(!is.na(value), !is.na(.data[[response]]))

  # --- Summary table (Rule R9: this is what the plot visualises) ------------
  summary_tbl <- long |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n        = dplyr::n(),
      spearman = suppressWarnings(
        cor(value, .data[[response]], method = "spearman",
            use = "pairwise.complete.obs")),
      p_value  = suppressWarnings(
        cor.test(value, .data[[response]], method = "spearman",
                 exact = FALSE)$p.value),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p_adj_bh   = p.adjust(p_value, method = "BH"),
      annotation = sprintf("ρ = %.2f\nq = %.2g", spearman, p_adj_bh)
    )

  # --- Rule R4: every label goes through the formatter ----------------------
  facet_labels <- setNames(formatter(summary_tbl$variable),
                           summary_tbl$variable)

  p <- ggplot(long, aes(x = value, y = .data[[response]])) +
    geom_point(alpha = 0.3, size = 0.6, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#4B0082", linewidth = 0.6) +
    geom_text(
      data = summary_tbl,
      aes(x = -Inf, y = Inf, label = annotation),
      hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
    ) +
    facet_wrap(~ variable, scales = "free_x",
               labeller = labeller(variable = facet_labels)) +
    labs(
      title    = paste0("Feature group: ", group, " vs ", formatter(response)),
      subtitle = sprintf("%d columns, Spearman with BH-adjusted q",
                         nrow(summary_tbl)),
      x        = NULL,
      y        = formatter(response)
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text    = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  # --- Rule R9: return both halves ------------------------------------------
  list(plot = p, table = summary_tbl)
}


# --- Top-to-bottom run (Rule §6.1 step 5) ------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  df  <- build_dataset("human")
  out <- feature_group_panel(df, group = "rnafold_zscores")
  print(out$plot)

  # Rule R8: outputs under OUTPUT_DIR
  dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(OUTPUT_DIR, "plots", "feature_group_panel_rnafold_zscores.jpg"),
         plot = out$plot, width = 297, height = 210, units = "mm", dpi = 300)
  write.csv(out$table,
            file.path(OUTPUT_DIR, "tables", "feature_group_panel_rnafold_zscores.csv"),
            row.names = FALSE)
}
```

Things to notice in the example:

- **Every rule reference is inline as a comment** so the reader can trace which rule motivates which line.
- **`fg_columns()` is used to enumerate**, not `fg()` — because we need the names as strings for the `summary_tbl$variable` column. `fg()` is for tidyselect, `fg_columns()` is for string vectors.
- **The summary table is what gets written to CSV** — that's the data source of truth, not the plot.
- **`labeller = labeller(variable = facet_labels)`** is the idiomatic way to apply `format_col_name()` to facet strips.

### 7.2 Distribution by predictor bins — `halflife` stratified by a chosen column

**Task:** plot the distribution of `halflife` across bins of a chosen predictor. Auto-detect whether the predictor is numeric (→ quantile bins) or categorical (→ native levels). Run a Kruskal-Wallis test for differences across bins. Facet by `species` automatically if the input contains more than one.

This example demonstrates a different shape from §7.1: a single predictor (not a group), distributions rather than scatter, categorical-vs-numeric input handling, and multi-species facetting via the `species` column (R6).

File: `analysis/correlations/halflife_by_predictor.R`

```r
# =============================================================================
# Halflife distribution stratified by a predictor
# =============================================================================
# Plot the distribution of halflife across bins of a chosen predictor:
#   - numeric predictor → quantile bins (default 4 bins = quartiles)
#   - categorical / logical / factor predictor → its native levels
# Tests whether halflife differs across bins via Kruskal-Wallis.
# Multi-species inputs are faceted by species (Rule R6).
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/halflife_by_predictor.R")
#
#   # Categorical:
#   df  <- build_dataset("human")
#   out <- halflife_by_predictor(df, predictor = "uorf_present_mrna")
#
#   # Numeric → quartiles:
#   out <- halflife_by_predictor(df, predictor = "mfe_delta_cds", n_bins = 4)
#
#   # Multi-species facet:
#   all <- build_all()
#   out <- halflife_by_predictor(all, predictor = "uorf_present_mrna")
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})


#' Halflife distribution stratified by a predictor (numeric → quantile bins,
#' categorical → native levels). Multi-species inputs are faceted by species.
#'
#' @param df         A dataframe from build_dataset() or build_all().
#' @param predictor  Character. Column name of the predictor.
#' @param n_bins     Integer. Number of quantile bins for numeric predictors.
#'                   Ignored for categorical predictors.
#' @param response   Character. Response column (default "halflife").
#' @param formatter  Function. Display formatter (default format_col_name).
#' @return list(plot, table). table columns: species (if present), bin, n,
#'   median, q25, q75, kw_chi2, kw_p, kw_p_adj_bh.
#' @export
halflife_by_predictor <- function(df,
                                  predictor,
                                  n_bins    = 4,
                                  response  = "halflife",
                                  formatter = format_col_name) {

  # --- Rule R5: guard every column access ----------------------------------
  if (!response  %in% names(df)) stop("response '",  response,  "' not in df")
  if (!predictor %in% names(df)) stop("predictor '", predictor, "' not in df")

  has_species <- "species" %in% names(df) &&
                 length(unique(df$species)) > 1

  # --- Auto-detect predictor type & bin if numeric -------------------------
  is_numeric_pred <- is.numeric(df[[predictor]]) &&
                     length(unique(stats::na.omit(df[[predictor]]))) >= n_bins

  if (is_numeric_pred) {
    qs <- stats::quantile(df[[predictor]],
                          probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE)
    if (length(unique(qs)) < n_bins + 1) {
      message("Quantile breaks not unique for '", predictor,
              "' — collapsing to ", length(unique(qs)) - 1, " bins")
      qs <- unique(qs)
    }
    df$.bin <- cut(df[[predictor]], breaks = qs, include.lowest = TRUE,
                   labels = paste0("Q", seq_len(length(qs) - 1)))
  } else {
    df$.bin <- as.factor(df[[predictor]])
  }

  # --- Drop missing --------------------------------------------------------
  df <- df |> dplyr::filter(!is.na(.bin), !is.na(.data[[response]]))
  if (nrow(df) == 0) stop("No non-NA rows for predictor / response combination")

  # --- Per-bin summary (Rule R9: this is the data the plot visualises) -----
  group_vars <- if (has_species) c("species", ".bin") else ".bin"
  summary_tbl <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n      = dplyr::n(),
      median = stats::median(.data[[response]], na.rm = TRUE),
      q25    = stats::quantile(.data[[response]], 0.25, na.rm = TRUE),
      q75    = stats::quantile(.data[[response]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(bin = .bin)

  # --- Kruskal-Wallis: per species if faceting, else overall ---------------
  kw_one <- function(sub) {
    res <- suppressWarnings(stats::kruskal.test(sub[[response]] ~ sub$.bin))
    tibble::tibble(kw_chi2 = unname(res$statistic), kw_p = res$p.value)
  }
  if (has_species) {
    kw_tbl <- df |>
      dplyr::group_by(species) |>
      dplyr::group_modify(~ kw_one(.x)) |>
      dplyr::ungroup() |>
      dplyr::mutate(kw_p_adj_bh = stats::p.adjust(kw_p, method = "BH"))
    summary_tbl <- dplyr::left_join(summary_tbl, kw_tbl, by = "species")
  } else {
    kw <- kw_one(df)
    summary_tbl$kw_chi2     <- kw$kw_chi2
    summary_tbl$kw_p        <- kw$kw_p
    summary_tbl$kw_p_adj_bh <- kw$kw_p  # only one test, no correction
  }

  # --- Plot (Rule R4: every label through formatter) -----------------------
  subtitle <- if (has_species) {
    "Kruskal-Wallis per species (BH-adjusted q in facet annotation)"
  } else {
    sprintf("Kruskal-Wallis: chi^2 = %.2f, p = %s",
            summary_tbl$kw_chi2[1],
            format.pval(summary_tbl$kw_p[1], digits = 2, eps = 0.001))
  }

  title_suffix <- if (is_numeric_pred) " (quantile bins)" else ""

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .bin, y = .data[[response]])) +
    ggplot2::geom_violin(fill = "#c8d8e4", alpha = 0.7, scale = "width") +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA,
                          fill = "white", alpha = 0.9) +
    ggplot2::labs(
      title    = paste0(formatter(response), " by ",
                        formatter(predictor), title_suffix),
      subtitle = subtitle,
      x        = formatter(predictor),
      y        = formatter(response)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 13),
      axis.text     = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank()
    )

  # --- Multi-species facet (Rule R6) ---------------------------------------
  if (has_species) {
    p <- p + ggplot2::facet_wrap(~ species)
    kw_labels <- summary_tbl |>
      dplyr::distinct(species, kw_chi2, kw_p_adj_bh) |>
      dplyr::mutate(label = sprintf("chi^2 = %.2f\nq = %.2g",
                                    kw_chi2, kw_p_adj_bh))
    p <- p + ggplot2::geom_text(
      data = kw_labels,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
    )
  }

  list(plot = p, table = summary_tbl)
}


# --- Top-to-bottom run (Rule §6.1 step 5) ------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  df <- build_dataset("human")

  # Example 1: categorical predictor
  out_cat <- halflife_by_predictor(df, predictor = "uorf_present_mrna")
  print(out_cat$plot)

  # Example 2: numeric predictor → quartiles
  out_num <- halflife_by_predictor(df, predictor = "mfe_delta_cds", n_bins = 4)
  print(out_num$plot)

  dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(OUTPUT_DIR, "plots", "halflife_by_uorf_present_mrna.jpg"),
         plot = out_cat$plot, width = 210, height = 148, units = "mm", dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "plots", "halflife_by_mfe_delta_cds_quartiles.jpg"),
         plot = out_num$plot, width = 210, height = 148, units = "mm", dpi = 300)

  write.csv(out_cat$table,
            file.path(OUTPUT_DIR, "tables", "halflife_by_uorf_present_mrna.csv"),
            row.names = FALSE)
  write.csv(out_num$table,
            file.path(OUTPUT_DIR, "tables", "halflife_by_mfe_delta_cds_quartiles.csv"),
            row.names = FALSE)
}
```

Things to notice that differ from §7.1:

- **Predictor type auto-detection** — `is.numeric()` plus a uniqueness check decides whether to bin or factorise. Plot functions that accept "any column" should always branch on the column's actual type, not assume.
- **Quantile binning with collapse fallback** — when the predictor has heavy ties (e.g. lots of zeros), quantile breaks can be non-unique. The example detects this and reduces the bin count rather than failing.
- **Hypothesis test stored in the table, not just the subtitle.** Per R9, the table must contain the data the plot visualises. The Kruskal-Wallis result is per-facet, so it's joined onto the summary table as columns that repeat across bins within a facet. A downstream consumer reading the CSV can recover the test without re-running it.
- **`species` faceting is conditional** — only applied when the dataframe actually contains multiple species. Single-species inputs get a flat plot. This is how R6 plays out in practice: code branches on the *content* of the `species` column, never on column names.
- **Two `ggsave()` calls in the runner block.** A script may produce multiple outputs; just give each a distinctive filename under `OUTPUT_DIR/plots/`.
- **`stats::` and `ggplot2::` prefixes** on namespace-ambiguous functions (`quantile`, `median`, `aes`) — safe practice in any script that may run in environments where many packages are attached.

---

## 8. Pre-submit checklist

Before considering any extension complete, run through this list. Each item maps to a numbered rule.

- [ ] Script starts with `source("R/load_all.R")` (R1)
- [ ] Data acquired via `build_dataset()` only (R2)
- [ ] All column-group selection uses `fg()` or `fg_columns()` (R3)
- [ ] All plot labels go through `format_col_name()` (R4)
- [ ] Every column access is guarded (R5)
- [ ] No species names in column names (R6)
- [ ] `CACHE_VERSION` bumped iff pipeline logic changed (R7)
- [ ] All outputs under `OUTPUT_DIR` (R8)
- [ ] Stats-driven plot returns `list(plot, table)` (R9)
- [ ] Response-correlation analyses exclude derived predictions (R10)
- [ ] Long-form loaders use `(transcript_id, region)` (R11)
- [ ] Wide-form loaders drop `gene_id` (R12)
- [ ] Any new column has a working `format_col_name()` result (no leftover underscores)
- [ ] Any new column has a row in `R/feature_table.csv`, and `Rscript scripts/check_feature_table.R` passes (§6.4)
- [ ] Column *subsets* use `include` / `exclude` / `top_n` through `select_features()` — never a subset row in the feature table (R3a)
- [ ] Plot functions that choose features take `include`, `exclude` (and `top_n`) and call `select_features()` (R3a)
- [ ] No backup / scratch / `*_old.R` file left anywhere under `R/` (R1)
- [ ] Script runs cleanly from a fresh R session via `Rscript <file>`
- [ ] Outputs are written to the correct subdirectory under `data/outputs/`
- [ ] Display-name spot check: `format_col_name(c(<your new cols>))` returns clean strings

---

## 9. Reference: where things live

```
R/                                  pipeline core (DO NOT scatter analysis here)
├── feature_table.csv               THE feature definitions (source of truth)
├── config.R                        paths, REGIONS, reads feature_table.csv,
│                                   CACHE_VERSION
├── load_all.R                      sources everything in dependency order
├── utils/                          pure helpers, no pipeline state
│   ├── normalise.R                 z_score_normalize, min_max_normalize
│   ├── naming.R                    format_col_name (labels from the table)
│   ├── palettes.R                  FEATURE_GROUP_COLOURS, REGION_COLOURS/SHAPES,
│   │                               format_group_name, format_metric_name
│   └── feature_groups.R            fg, fg_columns, select_features,
│                                   selected_columns, lookup_key
├── io/
│   ├── load_raw.R                  one function per raw source file
│   └── cache.R                     save/load/clear snapshot
├── features/
│   ├── mfe_model.R                 thermodynamic constants + math
│   └── engineer.R                  derived features, engineer_features()
└── pipeline/
    ├── assemble.R                  long→wide pivot, join helpers
    └── build_dataset.R             build_dataset, build_all

analysis/                           consumer code (this is where new plots go)
├── qc/                             diagnostics, sanity checks
├── correlations/                   anything correlation- or scatter-shaped
├── cross_species/                  species-vs-species comparisons
└── models/                         fitted models, coefficient plots

scripts/                            CLI runners (one per build, one-off jobs)
data/
├── raw/                            inputs you place
├── cache/                          .rds snapshots, auto-generated
└── outputs/                        all generated artifacts
    ├── plots/
    └── tables/
```

---

## 10. Known gotchas (read once, remember forever)

- **`build_dataset()` is cached by integer version, not by file hash.** Edit a raw input file without bumping `CACHE_VERSION` or passing `rebuild = TRUE`, and the stale cache silently wins.
- **Loaders skip missing files silently.** Watch the `skip (missing): ...` messages on first build — they're the only signal that an expected input wasn't found.
- **Numeric colour variables with fewer than 10 unique values are auto-discretised** in `create_scatter_plot()`. Don't be surprised when an integer `species_index`-style column shows up as a factor.
- **`add_mfe_expected_and_delta()` keys on `gc_content_{region}` and `length_{region}`.** If your loader produces `gc_{region}` instead, no expected/delta columns are computed. (Current state: `sequence_basic.tsv` produces `gc_content_*`.)
- **`engineer_features()` drops all-NA columns at the end.** A column whose loader produced only `NA` for this species will simply not appear in the output. This is by design but easy to forget when debugging "where did my column go".
- **The `species` column is added before `engineer_features()` is called.** Any engineering step that operates row-wise has access to `species` if it needs species-specific behaviour. Use this sparingly — most logic should be species-agnostic.
- **The pseudo-region tokens are long gone.** Whole-transcript scalars (architecture, uORF, NMD) end in the real `mrna` suffix; the `transcript` token and the `window`/`core`/`full` NMD tokens no longer exist. There is one NMD fragility model. Any analysis script written against the old names (`*_transcript`, `nmd_*_window`, …) will silently select nothing — `fg()` returns an empty set rather than erroring.
- **The schema is RNA-canonical: `u`, never `t`.** Nucleotide columns are `frac_u_*` and codons are `codon_aau_cds`, in every species. Upstream disagrees — the human counts file spells codons with U, the mouse one with T — so `normalise_codon_alphabet()` folds them at load time, the same way `normalise_region()` applies `REGION_ALIASES`. The consequence for anyone writing a regex over composition columns: the triplet class is `[acgtu]`, not `[acgt]`. A DNA-only class matches nothing rather than erroring, which is exactly how the v7 fraction bug survived — `add_codon_aa_fractions()` normalised only the 27 codons spelled without a U, leaving the other 38 as raw counts that correlated with `length_cds` at |rho| up to 0.82 while the normalised 27 summed to 1 among themselves and looked fine.
- **Every built column has a table row now, including the never-used ones** (Vienna median/p-value, analysis-window lengths, `utr5_length`, the CDS codon denominators, `stop_dist_last_downstream`). A few of these names still break the region-suffix-last invariant of §1.2 (`internal_exon_mean`, `n_overlapping_uorfs`, …); they are never-used, so no plot meets them, but renaming them means changing the loader and bumping `CACHE_VERSION`.
- **Not sure which namespace a string belongs to?** Call `lookup_key("mytoken")` — it returns `"flag"`, `"supergroup"`, `"group"`, `"feature"`, or `"unknown"`. Call `list_selection_keys()` to print all of them in one table.

---

## 11. When in doubt

Call `lookup_key("mytoken")` to find out which namespace a string belongs to, or `list_selection_keys()` to browse everything. Run `fg_columns(df, "<group>")` and `format_col_name("<column>")` to inspect what the schema actually produces in the current dataframe. Read the source file for deeper context (each module is ~50–200 lines).
