# Plan: one way to choose features

Handoff document. Written 2026-09-29 after an audit of every script; the
decisions in "Decisions already made" are the owner's and are not open.
Start this work on a **new branch from `main`**, not on the banded-plot branch.

## Why

Choosing which features go into a figure has grown many mechanisms. The
symptom that started this: removing "probing" from a figure meant tracing
three files, because probing entered through the `structure_core` bundle
(the whole Structure supergroup, with a `pick` keeping four icSHAPE Gini
columns from probing), and the table's "Included in exploratory analysis"
flag said nothing about that.

Goal: `R/feature_table.csv` says what a feature is and what it is eligible
for. Each script says what it wants with one vocabulary (`include` /
`exclude`), with identical behaviour everywhere.

## What exists today (audit)

Layers, all in effect at once:

| Layer | Where | Role |
|---|---|---|
| Table flags | `R/feature_table.csv`: "Included in exploratory analysis", "Included in model" | Eligibility. Group/supergroup expansion returns only exploratory features. Naming a feature id directly always works. |
| Hierarchy ids | Table Feature / Group / Supergroup columns; `FEATURE_PATTERNS`, `FEATURE_GROUPS`, `SUPERGROUPS` in `R/config.R` | The selectable names (one flat namespace). |
| Bundles | `GROUP_BUNDLES` in `R/config.R` | Named selections that can nest and carry `pick`/`drop`. |
| Default plot set | `DEFAULT_PLOT_GROUPS` in `R/config.R` | `c("nmd_core","junction_core","structure_core","sequence_select","translation_core")`. |
| Per-call `pick` / `drop` | `resolve_selection()`, `select_features()`, `refine_group_columns()` in `R/utils/feature_groups.R` | Column-level keep/drop keyed by feature id; caller wins over bundle per key. |
| Per-script extras | inside scripts | `keep_supergroups`, `top_n_per_group`, `standalones`, `DEFAULT_SWEEP_SKIP`, `min_abs_correlation`, `release_families`, `top_n`. |

Per script:

| Script | Selection today |
|---|---|
| `analysis/correlations/feature_correlation_ranked.R` | `DEFAULT_PLOT_GROUPS` + `pick` with `NULL` values (`release_families`) + `top_n_per_group` + `keep_supergroups` (collapse to "other") + `standalones`. |
| `analysis/correlations/feature_correlation_bands.R` (new, from the banded-plot branch) | Explicit group list in its runner + `release_families` + `top_n_per_group`. Calls `feature_correlation_ranked()` for the numbers. |
| `analysis/correlations/feature_correlation_dotplot.R` | `DEFAULT_PLOT_GROUPS` or named groups per job + `top_n_per_group` + `standalones` + `min_abs_correlation`. |
| `analysis/correlations/feature_response_scatter.R` | `DEFAULT_PLOT_GROUPS` + `top_n_per_group` + `standalones` + a separate `top_n` (different meaning: how many features shown). |
| `analysis/correlations/region_feature_heatmap.R` | `"structure"` or `DEFAULT_PLOT_GROUPS` + a `regions` argument + `top_n`. |
| `analysis/correlations/correlation_heatmap_workflow.R` | `DEFAULT_PLOT_GROUPS` via `select_features` + `top_n = 50` (pairs shown). |
| `analysis/correlations/feature_response_hex_panels.R` | Whole `"sequence"` and `"structure"` supergroups. |
| `analysis/correlations/feature_feature_correlation_table.R` | `setdiff(EXPLORATORY_FEATURES, c("codon_freqs","aa_freqs"))` hard-coded. |
| `analysis/correlations/group_panel_sweep.R` | Everything minus `DEFAULT_SWEEP_SKIP = c("codon_freqs","aa_freqs")`. |
| `analysis/qc/dataset_overview.R` | `groups = NULL` (all exploratory). |
| `analysis/models/xgb_structure_features.R` | Table model flag via `MODEL_FEATURES` + `model_columns()`; hand-written structure vs non-structure split with `SUPERGROUPS$structure`. |
| `analysis/models/xgb_structure_comparison.R`, `xgb_structure_plots.R` | No selection of their own. |
| `analysis/cross_species/cross_species_probing_concordance.R` | Hard-wired `fg_columns(df, "probing")`. |
| `scripts/example_analysis.R` | Direct `fg("...")`. |
| `scripts/check_feature_table.R` | Validates flags and bundles against a built cache. |
| Everything else in `analysis/` and `scripts/` | No feature selection. |

Findings: only `DEFAULT_PLOT_GROUPS` is truly shared; "skip codons and amino
acids" is written four different ways; trimming a big family is done both by
naming columns (bundle `pick`) and by `top_n_per_group`, and the named
codons are what made the broad figure disagree with the full 64-codon figure;
`top_n` means two different things; `standalones` is never given a non-empty
value by any caller.

## Decisions already made

- **Probing is Excluded from "core plots"** (it stays reachable directly, and
  the concordance script keeps using it).
- **Cache versions:** no bump needed. The one behaviour change is the probing
  exclusion, which alters default figures; that is accepted.
- **Hex panels and the sweep keep using "exploratory"**, not "core".
- **The ranked plot drops the collapse to "other"** and gets real supergroup
  facets. `keep_supergroups` goes away.
- **`standalones` is removed.** It was for features the owner had trouble
  placing; no longer needed. No caller passes it, so no figure changes.
- **Regions become their own argument** (`regions = NULL` means all) on the
  ranked, bands and dotplot scripts.

## Target design

1. **Table.** Keep the hierarchy and the two existing flags. Add a column
   "Included in core plots" (Included/Excluded per feature row, same style and
   parsing as the other two flags in `.read_feature_table()`). Fill it so the
   resolution equals today's `DEFAULT_PLOT_GROUPS` resolution, except probing
   is Excluded. For codon and amino-acid families the flag is for the whole
   family; the top-N rule below trims within them.
2. **One resolver** in `R/utils/feature_groups.R`:
   ```r
   select_features(df, include = "core", exclude = NULL, top_n = NULL)
   ```
   - `include`: ids at any level (feature, group, supergroup), or `"core"`
     (default; the flag), `"exploratory"`, `"model"`.
   - Result = expand `include`, then subtract expanded `exclude`. Nothing else
     applies silently.
   - `top_n = list(codon_freqs = 2, aa_freqs = 2)` is the only family
     trimmer: rank by |r| against the current response, keep N. Two figures
     can then never disagree about which codons are top.
   - Return column names together with the feature id each came from (plots
     need it for grouping and colour). Keep a plain character-vector view for
     callers that only want columns.
   - Rename the "how many pairs/features to show" arguments in the heatmaps
     and scatter to `max_features` so `top_n` only ever means "trim a family".
3. **Retire:** `GROUP_BUNDLES`, `DEFAULT_PLOT_GROUPS`, `DEFAULT_SWEEP_SKIP`,
   bundle handling and `pick`/`drop` merging in `resolve_selection()`,
   `release_families`, `standalones`, `keep_supergroups`,
   `BUNDLE_DISPLAY_NAMES` if nothing else uses it. Keep `fg()` and
   `fg_columns()` (explicit single-feature access).
4. **Regions** argument as above.

## Migration order

Each step is its own commit and leaves the project working.

1. **Baseline.** Write a script (scratchpad or `scripts/`) that records, per
   figure/job in each script, the exact column list it selects today, using
   the cached human dataset (`data/` is a symlink to the main checkout's data
   in worktrees; it is git-ignored). This is the reference for step 4.
2. **Add the "core plots" column** and extend `scripts/check_feature_table.R`
   to assert the flag reproduces the old `DEFAULT_PLOT_GROUPS` resolution
   (minus probing). Note `structure_core` currently means the Structure
   supergroup with probing narrowed to four `gini_cytoplasm_*` columns; after
   this change probing is out entirely.
3. **New `select_features()`** alongside the old API, with tests for
   expansion, exclusion, `"core"`, and the top-N rule.
4. **Migrate the four scripts that share `DEFAULT_PLOT_GROUPS` and the codon
   exclusion:** ranked (with supergroup facets, `regions`), bands, dotplot,
   scatter. Diff column lists against the step-1 baseline; the only intended
   difference is probing.
5. **Migrate the rest:** region heatmap, correlation heatmap workflow, hex
   panels (exploratory), sweep (exploratory), feature-feature table, QC
   overview, `xgb_structure_features.R` structure split.
6. **Delete the old machinery.** Grep must show no remaining references to
   the retired names; make `check_feature_table.R` fail if any reappear.
7. **Docs.** README / PIPELINE_GUIDE feature-table sections, the header
   comments in `config.R` and `feature_groups.R`, and any memory notes.

## Things to keep in mind

- Memory/decisions from earlier work still hold: human-only scope; table names
  win; the model flag stays independent of plotting flags and "follows the
  code for now"; the exon/junction rows are still unsettled by the owner.
- `analysis/correlations/feature_correlation_bands.R` and `R/colour_config.R`
  arrive from the banded-plot branch. If this work starts before that branch is
  merged, expect to rebase; the bands script's runner selection
  (`plot_groups`, `row_mm`) is the part that changes.
- The ranked script's runner has a guard (`sys.nframe() == 0 ||
  identical(environment(), globalenv())`) so it can be loaded into a private
  environment without running its jobs; the bands script relies on that.
- Verify by rendering figures and comparing column lists, not by reading code.
