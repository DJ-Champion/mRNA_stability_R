# RNAstab

A modular R pipeline for analysing mRNA transcript half-life and the sequence,
structural, and translational features that predict it. Designed to run on any
species with an Ensembl-style transcript annotation and the appropriate raw
inputs in place.

Built around the half-life data from Agarwal & Kelley, *The genetic and
biochemical determinants of mRNA degradation rates in mammals*, Genome Biology
23:245 (2022), but agnostic to data provenance — drop in any `halflife.csv`
keyed on `ensembl_gene_id` or `gene_name`.


## Quick start

```r
# 1. Place raw inputs under data/raw/<species>/
# 2. Load the pipeline:
source("R/load_all.R")

# 3. Build (or load from cache) a single species:
human <- build_dataset("human")

# 4. Or build everything you have configured and stack:
all_species <- build_all()      # tibble with a `species` column

# 5. Force a rebuild after changing raw inputs or feature logic:
human <- build_dataset("human", rebuild = TRUE)

# 6. For modelling, attach the blocked train/val/test split:
human <- attach_splits(build_dataset("human"))
```


## Directory layout

```
RNAstab/
├── R/
│   ├── feature_table.csv       # THE feature definitions (source of truth)
│   ├── config.R                # paths, species registry, REGIONS,
│   │                           # reads feature_table.csv, CACHE_VERSION
│   ├── load_all.R              # sources every R/ file in dependency order
│   ├── utils/
│   │   ├── normalise.R         # z_score_normalize, min_max_normalize,
│   │   │                       # normalize_numeric
│   │   ├── naming.R            # format_col_name() — canonical → display
│   │   ├── palettes.R          # FEATURE_GROUP_COLOURS, REGION_COLOURS/SHAPES,
│   │   │                       # format_group_name(), format_metric_name()
│   │   └── feature_groups.R    # fg(), fg_columns(), select_features(),
│   │                           # selected_features(), selected_columns(),
│   │                           # lookup_key(), list_selection_keys()
│   ├── io/
│   │   ├── load_raw.R          # one loader per source file
│   │   └── cache.R             # save_snapshot, load_snapshot, clear_snapshot
│   ├── features/
│   │   ├── mfe_model.R         # calculate_mfe_expected, calculate_mfe_delta
│   │   └── engineer.R          # add_mfe_expected_and_delta, impute_mrna_mfe,
│   │                           # add_mfe_per_nt, add_junction_density,
│   │                           # add_eej_min_distance, add_codon_aa_fractions,
│   │                           # drop_all_na_columns, engineer_features
│   └── pipeline/
│       ├── assemble.R          # pivot_regional_to_wide, join helpers
│       ├── build_dataset.R     # build_dataset(), build_all()
│       └── splits.R            # assign_folds(), assign_holdout(),
│                               # build_splits(), load_splits(),
│                               # attach_splits(), validate_splits()
├── analysis/
│   ├── qc/
│   │   ├── dataset_overview.R          # coverage / missingness diagnostics
│   │   ├── mfe_expected_check.R
│   │   └── imputation_check.R
│   ├── correlations/
│   │   ├── scatter_plot.R              # create_scatter_plot()
│   │   ├── halflife_correlation.R      # correlate_with_response() + plot
│   │   ├── top_correlations.R          # top_n_response_correlations()
│   │   ├── feature_correlation_dotplot.R
│   │   ├── feature_response_scatter.R
│   │   ├── feature_response_hex_panels.R
│   │   ├── feature_feature_correlation_table.R
│   │   ├── correlation_heatmap_workflow.R
│   │   ├── region_feature_heatmap.R
│   │   ├── group_panel_sweep.R         # one panel per schema family
│   │   └── quadrant_export.R
│   ├── cross_species/
│   │   └── cross_species_probing_concordance.R
│   └── models/                         # two XGBoost models: Baseline, Structure
│       ├── xgb_structure_features.R    # feature blocks + eligible row set
│       ├── xgb_structure_comparison.R  # tune, fit, evaluate, write tables
│       └── xgb_structure_plots.R       # figures, from the tables only
├── scripts/
│   ├── preprocess_saluki.R     # one-off: HDF5 → .rds
│   ├── build_human.R           # CLI runner for human
│   ├── build_mouse.R           # CLI runner for mouse
│   ├── build_splits.R          # CLI runner: writes the blocked split (once)
│   └── example_analysis.R      # worked end-to-end example
├── data/
│   ├── raw/                    # populated by you
│   │   ├── human/
│   │   ├── mouse/
│   │   └── shared/             # incl. family.tsv, the clustering seam
│   ├── cache/                  # auto-generated .rds snapshots
│   ├── splits/                 # auto-generated blocked train/val/test
│   └── outputs/
│       ├── plots/
│       └── tables/
├── PIPELINE_GUIDE.md           # rulebook for extending the pipeline
└── README.md
```


## How it works

### One row per (transcript, species)

Every row of the built dataset is one transcript for one species. A `species`
column makes cross-species analysis trivial:

```r
combined <- bind_rows(build_dataset("human"), build_dataset("mouse"))
combined |> group_by(species) |> summarise(median_hl = median(halflife, na.rm = TRUE))
```

Or just call `build_all()`.

### No species prefix on column names

Columns are named `{metric}_{region}` in lowercase, e.g. `rnafold_zscore_5utr`,
`length_cds`, `gc_3utr`. Half-life is just `halflife`. The species lives in a
column, not in column names. This is what makes the cross-species stack work.

For pretty plot labels, `format_col_name()` (in `R/utils/naming.R`) turns
canonical column names into display strings (`"rnafold_zscore_5utr"` →
`"MFE.z 5' UTR"`, `"gc_content_cds"` → `"C+G% CDS"`), taking each label from
the feature table's short name. Its sibling `format_group_name()` (in
`R/utils/palettes.R`) does the same for *selection keys* — feature, group,
and supergroup ids (`"rnafold_zscores"` → `"MFE z-score"`,
`"structure"` → `"Structure"`) — for plots that label a facet strip or legend. A third, `format_metric_name()`, strips the region
suffix for plots where region is already encoded as colour or shape
(`"length_cds"` → `"Length"`).

### Region vocabulary

```
5utr   cds   3utr   mrna   utrpair   last100   start   stop
```

These appear as suffixes on metric columns. `utrpair` is the combined non-coding
regions (5' UTR + 3' UTR)

### The feature table

**`R/feature_table.csv` is the single source of truth for every feature.** Each row is one feature:

- its columns (a regex);
- its place in the **Supergroup > Group > Feature** hierarchy;
- whether it is **Included in exploratory analysis**, **Included in model** and **Included in core plots**;
- its display names, colour, regions, and the reason it is in or out.

The code derives everything from it, so to change a feature, edit the table. Then check it against the built cache:

```bash
Rscript scripts/check_feature_table.R
```

The three flags do different jobs:

- **Exploratory analysis:** selecting a group or supergroup returns only its exploratory features. Naming a feature id directly always works, so a plot can still show an excluded feature on purpose.
- **Core plots:** the default set for the correlation figures (`"core"`). A subset of the exploratory features; experimental probing (icSHAPE) is deliberately outside it, and reachable by naming `"probing"`.
- **Model:** `model_columns(df)` gives the model's predictors, and nothing else reads this flag.
- **Excluded from both:** the column is used nowhere, and `drop_excluded()` removes it.

Reach for `fg("rnafold_zscores")` instead of typing out region patterns:

```r
df |> select(halflife, fg("rnafold_zscores"), fg("rnalfold_zscores"))

# What does a feature resolve to?
fg_columns(df, "rnafold_zscores")
# [1] "rnafold_zscore_3utr"    "rnafold_zscore_5utr"    "rnafold_zscore_cds"
# [4] "rnafold_zscore_last100" "rnafold_zscore_mrna"    "rnafold_zscore_start"
# [7] "rnafold_zscore_stop"    "rnafold_zscore_utrpair"
```

`list_selection_keys()` prints every flag, supergroup, group and feature with its display name. `lookup_key("foo")` tells you which kind a single token is.

### Choosing which columns to plot

Every analysis chooses its columns the same way, with two arguments:
`include` and `exclude`. The feature table says what a feature *is* and what it
is eligible for; the call says what this analysis wants. Narrowing a plot never
means editing the table.

Both arguments take **tokens**:

- **A flag** — `"core"` (the default), `"exploratory"` or `"model"`: the
  features the table flags for that purpose.
- **A supergroup** — e.g. `"structure"`. Expands to its exploratory features in
  `include`, and to all of its features in `exclude`.
- **A group** — the table's Group column snake-cased, e.g. `"global_folding"`.
  Same rule as a supergroup.
- **A feature** — a table row, e.g. `"codon_freqs"`. Always works, whatever its
  flags.

The result is `include` expanded, minus `exclude` expanded. Nothing else applies
silently. `select_features()` returns the columns present in your dataframe,
each with the feature it came from; `selected_columns()` gives the plain
vector:

```r
select_features(df)                                   # the core set
select_features(df, "structure")                      # a whole supergroup
select_features(df, "core", exclude = "sequence")     # core minus a supergroup
select_features(df, "exploratory",
                exclude = c("codon_freqs", "aa_freqs"))
df |> select(all_of(selected_columns(select_features(df, "global_folding"))))
```

Big families are trimmed with `top_n`, the only family trimmer. It keeps the N
metric stems with the largest |r| against a `response` (the maximum over
regions), with all their regions, so every figure using the same response
agrees on which codons are "top":

```r
select_features(df, top_n = list(codon_freqs = 2, aa_freqs = 2),
                response = "halflife")
```

Half-life and translation efficiency may then show different codons; that is
intended. `regions = c("5utr", "cds")` restricts columns to those regions
before `top_n` ranks them.

The plotting functions (`feature_correlation_ranked()`, `_bands()`,
`_dotplot()`, `feature_response_scatter()`, `region_feature_heatmap()`,
`run_correlation_heatmap_workflow()`, `feature_response_hex_panels()`, the
sweep, the feature-feature table, the QC overview) all take `include`,
`exclude` and, where it makes sense, `top_n` and `regions`, with these meanings.
Arguments that limit *how many things are drawn* are called `max_features`;
`top_n` always means "trim a family". The default is `"core"` for the
correlation figures and `"exploratory"` for the hex panels, sweep,
feature-feature table and QC overview.

No `CACHE_VERSION` bump is ever needed for any of this — it is selection logic,
not feature engineering.

### Caching

After every successful build, the dataset is saved to
`data/cache/{species}_dataset_v{CACHE_VERSION}.rds`. Subsequent calls to
`build_dataset()` return the cache instantly. Three ways to invalidate:

```r
build_dataset("human", rebuild = TRUE)   # one-off rebuild
clear_snapshot("human")                  # delete cache file
# Bump CACHE_VERSION in R/config.R       # invalidates everyone's cache
```

Bump `CACHE_VERSION` whenever feature-engineering logic changes — that's the
mechanism for keeping caches honest after refactors.

### Sequence families and blocked splits

Genes are grouped into **sequence families** so that homologous genes never
land on opposite sides of a train/test boundary. Training on one member of a
paralogue pair and testing on the other inflates the score, and nothing in the
model output reveals it.

Families are computed upstream in Python (see `FAMILY_CLUSTERING.md`) and
arrive as `data/raw/shared/family.tsv`. `load_family()` ingests them, giving
every gene a `family_id_{strict,medium,loose}` label. `BLOCK_LEVEL` in
`R/config.R` picks which one blocks the splits — `medium`, a measured choice.

The split itself is a **separate artefact, generated once**:

```bash
Rscript scripts/build_splits.R                # writes data/splits/holdout_medium.rds
Rscript scripts/build_splits.R --level loose  # sensitivity check
```

then read everywhere via `attach_splits(build_dataset("human"))`. Never
re-derive it inside an analysis script: a rebuilt `family.tsv` or a changed row
order can move genes between train and test, and results stop being comparable
with nothing looking wrong. `attach_splits()` warns if the dataset and the
split were built from different `family.tsv` files.

**Whole families move together — a family is never divided across splits.**
The packer's only decision is *which* split each intact family goes to. What
it balances is threefold: gene counts (80/10/10), and family structure, and
the size of the largest family any split receives. The summary printed by
`build_splits()` reports all three:

```
 split n_genes pct n_families max_family pct_multi
  test    1360  10       1061         20      35.3
 train   10881  80       8484        282      35.3
   val    1360  10       1060         29      35.4
```

- `max_family` — the largest family that split received, *entirely*. Test's
  biggest family has 20 members and all 20 are in test.
- `pct_multi` — the share of that split's genes having at least one relative
  (which is necessarily in the same split). The three figures should be close
  to each other; that is what makes test resemble train.

A split whose `pct_multi` is near zero holds nothing but genes with no
relatives, which is an unrepresentative slice of the genome even though the
gene counts look perfect. `validate_splits()` checks for it.

Two caveats to state in any write-up using this split:

1. It measures **generalisation to novel gene families** — a stronger claim
   than random-over-genes, and a different one from cross-species transfer.
2. Families larger than 5% of the smallest split (68 genes) are pinned to
   `train` wholesale, so held-out splits are depleted of the very largest
   families. Only one family currently qualifies.

Species absent from the clustering cohort (currently mouse) get no family
columns and `NA` splits. Making mouse blockable means adding it to the cohort
**upstream** and re-clustering both species together — cross-species
orthologues must merge into one family, or the blocking does not prevent
cross-species leakage.

### Adding a new species

Three steps:

1. Add an entry to `SPECIES_CONFIG` in `R/config.R` (copy the human or mouse
   block; each entry is just `dir` plus an optional `saluki_rds`).
2. Place raw files under `data/raw/<species_dir>/` matching the filenames the
   loaders expect — `grep -h "species_path" R/io/load_raw.R` lists them all.
3. `build_dataset("rat")`.

No other code changes needed.

### Adding a new feature

Two patterns depending on what it is.

**A new derived feature** (computed from existing columns): add a function to
`R/features/engineer.R`, then call it in `engineer_features()`. Bump
`CACHE_VERSION`. Done.

**A new raw input source**: add a loader to `R/io/load_raw.R`, then add it to
the appropriate join block in `R/pipeline/build_dataset.R`. Bump
`CACHE_VERSION`.

Either way, give every new column a row in `R/feature_table.csv` (PIPELINE_GUIDE
§6.4) and run `Rscript scripts/check_feature_table.R`; it fails on any column
without a row.


## Common analysis recipes

```r
source("R/load_all.R")
df <- build_dataset("human")

# --- Top correlations with half-life ---
source("analysis/correlations/halflife_correlation.R")
out <- halflife_correlation_plot(df, top_n = 20)   # list(plot, table)
print(out$plot)
write.csv(out$table,
          file.path(OUTPUT_DIR, "tables", "halflife_spearman.csv"),
          row.names = FALSE)

# --- A custom scatter ---
source("analysis/correlations/scatter_plot.R")
create_scatter_plot(df,
                    x_var = "rnafold_zscore_mrna",
                    y_var = "halflife",
                    color_var = "length_mrna",
                    log_color = TRUE,
                    add_density_rings = TRUE)

# --- The default feature set, as a dotplot ---
source("analysis/correlations/feature_correlation_dotplot.R")
out <- feature_correlation_dotplot(df)          # include = "core"
print(out$plot)
```


## The analysed cohort

`build_dataset()` drops any transcript whose 5'UTR **or** 3'UTR is shorter than
`MIN_UTR_LENGTH` (30 nt), or whose UTR length is missing. A UTR of a few
nucleotides is an absent or mis-annotated one, and it makes every regional
feature degenerate — folding energy over 12 nt is not comparable to folding
energy over 1,200, and the length-normalised z-scores divide by a
shuffled-sequence distribution that is itself near-degenerate.

```
human  13,660 → 12,302 built rows  (9.9% removed; 12,277 modellable)
mouse  14,197 → 13,215 built rows  (6.9% removed)
```

The filter applies to the frame `build_dataset()` **returns**, not to what it
writes — so `data/cache/*.rds` stays complete, changing the threshold never
invalidates a cache, and no `CACHE_VERSION` bump is involved. It is selection
intent, like the feature table's flags.

```r
df  <- build_dataset("human")                   # filtered — the default
all <- build_dataset("human", min_utr = NULL)   # everything, for QC
```

The three `analysis/qc/` scripts pass `min_utr = NULL`: a coverage and
missingness diagnostic should describe the whole built table, including what
the filter removes.

**The split artefact does not need rebuilding.** Blocking survives any
subsetting — removing genes cannot make a family span two splits — and the
proportions barely move (80.12 / 9.99 / 9.90 against a target of 80/10/10,
inside `validate_splits()`'s tolerance). `holdout_medium.rds` stays valid, so
results remain traceable to the clustering run behind it.


## Models

Two, and only two:

| Model | Predictors |
|---|---|
| `Baseline` | non-structure transcript features (lengths, GC, skews, codon composition, CAI, stop-free lengths, uORF presence, exon density, junction distances, NMD) |
| `Structure` | `Baseline` **+** every computed secondary-structure feature (RNAfold / RNALfold MFE, their z-scores against shuffled sequence, per-nucleotide MFE, MFE delta) |

Structure is the only difference between them: both are handed the same rows,
the same gene ids, the same preprocessing, the same family-blocked tuning
resamples, the same 60-configuration grid and the same seed. The single
pre-specified contrast is `Structure` vs `Baseline` on held-out R², evaluated
once on the untouched `test` split of the committed family-blocked holdout.

```bash
Rscript analysis/models/xgb_structure_comparison.R      # tune, fit, evaluate, plot
XGB_GRID_SIZE=6 Rscript analysis/models/xgb_structure_comparison.R   # fast smoke test
Rscript analysis/models/xgb_structure_plots.R           # figures only, from the tables
```

Fits are cached in `data/outputs/xgb_structure/final_fits.rds`, keyed on
everything that determines them (both predictor lists, the exact genes, the
seed, the grid, the fold count). Change any of those and it re-tunes; set
`XGB_REFIT=1` to force it.

Two things to carry with any number from this comparison, both stated in
`SUMMARY.txt` on every run:

- **The structure block is not length- and GC-neutral.** Raw and
  per-nucleotide MFE scale with length and shift with GC, both already in the
  baseline. `xgb_structure_redundancy.csv` reports how much of each structure
  column the baseline already explains — read it alongside the headline delta.
- **Structure missingness is informative.** A missing 5'UTR MFE means a 5'UTR
  too short to fold. XGBoost handles NA natively, so `Structure` can split on
  an annotation artefact; this design cannot rule that out.

icSHAPE structural Gini (`probing`) is deliberately in neither model: it is
measured rather than computed from sequence, so a model using it cannot score
an unprobed transcript. If it is ever modelled, that is a separate,
supplementary analysis.


## Required R packages

Core pipeline: `dplyr`, `tidyr`, `readr`, `purrr`, `stringr`, `tibble`,
`tidyselect`, `rlang`.

Analysis layer: `ggplot2`, `forcats`, `viridis`, `scales`.

Modelling: `tidymodels`, `finetune`, `xgboost`, `future`.

Saluki preprocessing only: `rhdf5` (Bioconductor).


## Known limitations

- Cache invalidation is by integer version, not by hashing the source files.
  If you edit a raw file without bumping `CACHE_VERSION` or passing
  `rebuild = TRUE`, the stale cache wins. (For per-stage caching with automatic
  staleness detection, the {targets} package is the natural next step.)
- The Saluki preprocessing step is one-off and not integrated into the main
  build. If Saluki outputs change, run `scripts/preprocess_saluki.R` again
  before rebuilding.
- Loaders silently skip missing files. Good for incomplete species, bad if you
  expected a file to be picked up. Watch the `skip (missing): …` messages on
  the first build.
- **A few never-used columns break the region-suffix-last convention**
  (`internal_exon_mean`, `n_overlapping_uorfs`, `stop_dist_last_downstream`,
  …). Every one has a feature-table row marked excluded from both flags, so no
  plot meets them; renaming them means changing the loader and bumping
  `CACHE_VERSION`.
