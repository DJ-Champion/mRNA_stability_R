# Feature table audit — supplementary table vs. pipeline

Reconciliation of `mRNA_stability_supplementary_table - All Features.csv` (54 feature
rows) against the code and the built data. **Report only — no code, config, CSV or
analysis script was modified in producing this.**

## Sources compared

| Source | What was taken from it |
|---|---|
| The CSV | 54 rows: Included/Excluded flag, Supergroup, Group, Feature, display names, Calculated-using, Regions, Notes |
| `R/config.R` | `FEATURE_PATTERNS` (25 groups), `SUPERGROUPS` (6), `GROUP_BUNDLES` (7), `INCLUDED_GROUPS`, `EXCLUDED_FEATURES` (76 entries), `REGIONS`, `META_COLS` |
| `R/utils/naming.R` | `REPLACEMENTS` (the `format_col_name()` rule list) |
| `data/cache/human_dataset_v10.rds` | 308 real columns |
| `data/cache/mouse_dataset_v10.rds` | 233 real columns |
| `data/raw/{human,mouse,shared}/` | which raw inputs exist and which are actually loaded |

"Pool" below means the covariate pool a modelling script sees: every column surviving
`drop_excluded()` after removing `META_COLS` and the `halflife` response.
**Human pool = 222. Mouse pool = 193. Intersection = 176.**

### What is already correct

Worth stating up front, because most of this document is problems. Three things
`config.R` asserts in comments and nothing tests — all three hold on the real data:

- `FEATURE_PATTERNS` patterns are mutually exclusive. No column matches two groups.
- Every `FEATURE_PATTERNS` key sits in exactly one supergroup. None omitted, none doubled.
- `format_col_name()` leaves no underscore or lowercase region token in any of the 541
  distinct column names across both caches. The v4 region-suffix fix holds.

`EXCLUDED_FEATURES` is also clean: no duplicates, and every one of the 76 entries
matches a real column in at least one species.

---

# Part 1 — Scientific decisions

These change what the models see or what the paper claims. They are not typos, and I
have not guessed at them.

## S1. Junction density vs. exon density — the table and the code disagree about which encoding is the feature

The single largest conflict, and it runs in both directions.

| | CSV says | Code does |
|---|---|---|
| `junctions_density_{5utr,cds,3utr,mrna}` | **Included** (row 23, "Junction density") | **Excluded** — the whole `junctions` group is in `EXCLUDED_FEATURES`; 0 of 8 columns reach the pool |
| `exon_density_{5utr,cds,3utr,mrna}` | **Excluded** (row 39, "CDS-exon density") | **Included** — all 4 in the pool, and picked by the `translation_core` bundle, so they are in the `INCLUDED_GROUPS` default |

So the table marks the excluded encoding Included and the included encoding Excluded.
Whatever the resolution, one of the two documents is currently misdescribing every
figure that uses the default selection.

Three things make this more than a label swap:

1. **The two are not the same quantity.** `engineer.R` computes
   `exon_density = (junctions + 1) / kb`, against `junctions_density = junctions / kb`.
   They correlate ρ = 0.97 on human CDS, but they differ exactly where it matters: for
   the **353 single-exon human transcripts**, `junctions_density` is 0 and
   `exon_density` is `1/kb`. Whether an intronless transcript reads as "zero junctions"
   or "one exon" is a modelling choice with a biological argument on each side (EJC
   deposition vs. exon count).
2. **They sit in different supergroups.** The CSV files junction density under
   *Transcript architecture*; the code files `exon_density` under `translation`. A
   supergroup-coloured figure and the table will not agree on the category.
3. **The display label is actively misleading.** `naming.R` renders `exon_density_cds`
   as **"EEJ dens. CDS"**, and `junctions_density_cds` as "Junction density CDS". So a
   plot of `exon_density_*` is currently labelled as if it were the junction density —
   the very column the pipeline excluded. Anyone reading a figure would conclude the
   table's row 23 was right.

**Recommendation.** Decide the encoding on the single-exon argument, then make all four
of these agree: the `EXCLUDED_FEATURES` entry, the supergroup, the `naming.R` label, and
the CSV row. My inclination is to keep `exon_density` (the `+1` makes it a genuine
per-kb exon count, and it is what every fitted model has actually used) but **relabel it
"Exon dens."** and move it to `splicing`, since calling it EEJ density is what created
this confusion. That is a call for you.

## S2. Saluki predictions never reach the dataset — a silent join failure

The CSV marks Saluki prediction / prediction difference **Included** (row 55). Neither
`saluki_prediction` nor `prediction_difference` exists in either cache.

The cause, in [R/pipeline/assemble.R:60](R/pipeline/assemble.R:60):

```r
join_gene_level <- function(wide_df, gene_level_dfs) {
  for (df in dfs) {
    if (!"gene_id" %in% names(df)) next     # <- silently skips
    wide_df <- left_join(wide_df, df, by = "gene_id")
  }
}
```

`data/raw/human/saluki_predictions.rds` is a 12,968-row data frame whose key column is
named **`ensembl_gene_id`**, not `gene_id`. The docstring on `load_saluki_predictions()`
([R/io/load_raw.R:180](R/io/load_raw.R:180)) says "Expected columns: gene_id,
saluki_prediction" — the file does not meet that contract, the guard fires, and the
join is skipped with no message. The file loads fine; it is simply discarded.

Downstream consequences:
- The two `naming.R` rules for `^saluki_prediction$` and `^prediction_difference$` are dead.
- `prediction_difference` is never computed at all, since it needs the join to have happened.
- Any analysis benchmarking against Saluki has been silently running on nothing.

**Recommendation.** Rename the key in `preprocess_saluki.R` or in the loader, rebuild,
and confirm the column appears. Separately, `join_gene_level()`'s `next` should `warn`
rather than skip silently — that guard is what turned a one-word mismatch into an
invisible one. Both need a `CACHE_VERSION` bump.

## S3. Mouse has no structure features whatsoever

The CSV has no species column, so it reads as though it describes both species. It does
not. Seven whole groups are human-only:

| Group | Human | Mouse |
|---|---|---|
| `rnafold_scores` | 8 | **0** |
| `rnafold_zscores` | 8 | **0** |
| `rnafold_per_nt` | 7 | **0** |
| `mfe_deltas` | 7 | **0** |
| `mfe_expected` | 7 | **0** |
| `rnalfold_scores` | 7 | **0** |
| `rnalfold_zscores` | 7 | **0** |
| `standalone` (`cai`, `translation_efficiency`, `orfexondensity`) | 3 | **0** |

That is the entire `structure` supergroup apart from `probing`, plus CAI and TE. The
mouse pool of 193 features contains **zero folding energies**.

The comment on `drop_excluded()` in [R/utils/feature_groups.R:308](R/utils/feature_groups.R:308)
says mouse "legitimately lacks the human-only Vienna median/pval families". That is a
substantial understatement — those are the *excluded* families. Mouse lacks the modelled
ones too.

Two knock-on effects:
- `INCLUDED_GROUPS` resolves to 21 groups including all seven above. Run the default
  correlation dotplot on mouse and it silently returns a structure-free figure.
- `intrinsic_core` and `intrinsic_select` both carry `pick = list(standalone = "cai")`,
  which resolves to nothing on mouse. `select_features()` emits a `message()` about it,
  which in a scripted run scrolls past unnoticed.

**Recommendation.** This determines the shape of the tidied table. If mouse folding is
planned, the table needs a `Species` column with "human only (mouse pending)". If mouse
is deliberately a sequence-only replicate, that is a paper-level caveat and the table
should say so explicitly. Please tell me which.

## S4. The same regex means different features in different species

`stopfree = "^stopfree_"` resolves to **4 columns on human and 14 on mouse**. Mouse
carries ten `stopfree_*` columns human lacks, and none are excluded, so they enter the
mouse pool as covariates:

- `stopfree_length_{start,stop,last100}` — the fixed-width analysis windows, whose
  `length_*` counterparts are explicitly excluded as "not measured transcript properties"
- `stopfree_fraction_{3utr,5utr,cds,mrna,start,stop,last100}` — a fraction encoding that
  has no human equivalent

Worse, mouse also carries seven `stopcodon_fraction_*` columns that **match no
`FEATURE_PATTERNS` group at all**. They are therefore invisible to `select_features()`
and `fg()`, but fully visible to any modelling script that builds its matrix as
`setdiff(names(df), c(META_COLS, TARGET_COL))` — which, per the comment on `META_COLS`,
is exactly how the modelling scripts do it. Seven undocumented, ungrouped, mouse-only
covariates.

Together these 17 columns are the entire mouse-pool-only set.

**Recommendation.** Decide whether the fraction encodings are wanted. If yes, build them
for human and give `stopcodon_fraction_*` a group. If no, add all 17 to
`EXCLUDED_FEATURES`. What must not persist is one regex resolving to two different
feature sets by species — that breaks the cross-species comparison the pipeline exists
to support.

## S5. RNAplfold accessibility — raw data exists, nothing consumes it

CSV row 31 marks "Pairing probability / Accessibility (RNAplfold)" as **Excluded**.

`data/raw/human/rnaplfold_results.csv` exists and is well-formed long-form data
(`transcript_id, region, score`). There is **no loader, no `FEATURE_PATTERNS` entry, no
`naming.R` rule and no column** anywhere in the pipeline. `grep -r "plfold"` over `R/`,
`scripts/` and `analysis/` returns nothing.

So "Excluded" here means "never built", which is a different claim from the "Excluded"
on row 28 (MFE null median — built, cached, deliberately kept out of the pool). The
table currently uses one word for both.

**Recommendation.** Split the flag into three states — `Included` / `Built, excluded` /
`Not built` — and decide whether to ingest RNAplfold. Given the raw data is sitting
there and accessibility is a distinct structural axis from MFE, this looks like the
cheapest real feature addition available.

## S6. `exon_length_last_mrna` is in the covariate pool, and it is a near-duplicate of `length_3utr`

CSV row 45 marks "Last exon length" **Excluded**. Every sibling is in
`EXCLUDED_FEATURES` — `exon_count_internal_mrna`, `intron_length_mean_mrna`,
`internal_exon_{mean,median,sd}`, `n_exons`. `exon_length_last_mrna` is **not**, so it
survives `drop_excluded()` and sits in both pools. It is the sole surviving member of the
`exons` group (`keptH=1`, `keptM=1`).

This matters more than a bookkeeping slip because on the human cache:

```
cor(exon_length_last_mrna, length_3utr, method = "spearman") = 0.949
```

The last exon is mostly the 3' UTR, so this is a collinear restatement of a retained
core feature. Every model fitted to date has carried both, which will have split their
importance and destabilised any ranking between them.

**Recommendation.** Add it to `EXCLUDED_FEATURES` (the table already says it should be),
and re-check any feature-importance result that ranks `length_3utr`. Flagged as
scientific rather than clerical because it has already affected fitted results.

## S7. The three stop codons are in the pool, and they encode something real

CSV row 17 says "61 sense codons". `codon_freqs` matches **64** columns — the 61 sense
codons plus `codon_uaa_cds`, `codon_uag_cds`, `codon_uga_cds`, all three of which are in
the pool. (`codon_other_cds` is the excluded 65th bucket, correctly documented.)

These three are not noise. On the human cache they are close to a one-hot encoding of
which stop codon terminates the CDS:

| Column | Transcripts with non-zero value |
|---|---|
| `codon_uaa_cds` | 3,740 |
| `codon_uag_cds` | 3,036 |
| `codon_uga_cds` | 6,843 |
| (total) | 13,619 vs. 13,601 rows |

Stop-codon identity is a plausible NMD-relevant covariate, so including them may well be
right — but it should be a decision, and the table should say 64, not 61, and explain
what the three extra columns mean. As written, a reader would assume the stop codons
were excluded.

**Recommendation.** Keep them, correct the count to 64, and add a sentence that the three
stop-codon columns act as a stop-identity indicator rather than a usage frequency.

---

# Part 2 — Clerical corrections

Real errors, but none require a judgement call.

## C1. Group names in the CSV Notes that do not exist in the code

| CSV rows | Notes say | Actual `FEATURE_PATTERNS` key |
|---|---|---|
| 40–46 | `architecture` | split across `introns`, `exons`, `noncoding` — no `architecture` key exists |
| 33 | `orfs` | no such key |
| 24, 48, 49 | `junctions` | `eej_dist` |
| 13, 14 | `nuc_ratios` | `compositional` (`^(purine_\|amino_)`) |

## C2. Supergroup vocabulary does not match `SUPERGROUPS`

| CSV Supergroup | Code supergroup | Note |
|---|---|---|
| Structure | `structure` | ✓ |
| Sequence | `intrinsic` | name differs |
| Transcript architecture | `splicing` | name differs |
| RNA decay | `decay` | ✓ |
| Translation | `translation` | ✓ but membership differs — CSV files CAI and TE here; code puts them in `other` via `standalone` |
| Expression | — | no code equivalent, and no expression column exists |
| Response / evaluation | — | no code equivalent; `halflife` is the response, correctly outside `FEATURE_PATTERNS` |

## C3. Region lists that do not match the built columns

- **Row 21, uORF presence** — Regions given as "5' UTR, CDS". The column is
  `uorf_present_mrna`; region is `mrna`.
- **`utrpair` appears nowhere in the CSV.** It is a real token in `REGIONS` ("UTR
  interactions") with 8 human columns, two of which — `rnafold_score_utrpair` and
  `rnafold_zscore_utrpair` — are in the pool. Rows 2 and 3 list seven regions and should
  list eight.
- **Rows 7, 8 (RNALfold)** correctly omit `utrpair` — RNALfold has 7 regions, not 8. The
  asymmetry between RNAfold and RNALfold is real and worth a note.
- **Row 9, icSHAPE Gini** — regions are right, but the compartment axis is only hinted at
  in the long name. There are two compartments (`nucleoplasm`, `cytoplasm`) × 4 regions =
  8 columns, all in the pool.
- **Row 19, Stop-free length** — the only feature row with an entirely empty Notes cell;
  no `fg group` recorded. It is `stopfree`. See also S4.

## C4. Excluded columns with no CSV row

Built, cached, in `EXCLUDED_FEATURES`, absent from the table:

| Column(s) | Why excluded (per `config.R`) |
|---|---|
| `rnalfold_median_*` (7) | Vienna auxiliary stats. The CSV has the local-MFE *p-value* row (30) but no local-MFE *median* row, though it has the global median row (28). |
| `length_{last100,start,stop}` | fixed-width analysis windows |
| `utr5_length` | exact duplicate of `length_5utr` |
| `cds_length_codons_cds`, `n_codons_scored_cds`, `n_stops_cds` | CDS-size denominators that row normalisation exists to remove |
| `codon_other_cds` | the unresolvable-codon bucket; identically zero in both species |
| `stop_dist_last_downstream` | superseded by `eej_dist_closest_stop` |

## C5. CSV rows describing nothing that exists

| Row | Feature | Status |
|---|---|---|
| 31 | Accessibility (RNAplfold) | raw data exists, never ingested — see S5 |
| 53 | Expression | no column in either cache. Note also claims "in code 'standalone' group" — it is not; `standalone` is `cai`, `translation_efficiency`, `orfexondensity` |
| 33 | uORF length | no column; also cites the non-existent `orfs` group |
| 46 | Non-coding length fraction | no column in either species — see C6 |
| 44 | First exon length | no column in either species — see C6 |
| 32 | uORF count | the column is `uorf_count_mrna`; the Notes name it `n_uorfs`, which does not exist |

## C6. Dead code with no data behind it

- **`noncoding` is an empty group.** `FEATURE_PATTERNS$noncoding = "^noncoding_"` matches
  **0 columns in both species**, yet it holds a `FEATURE_PATTERNS` key, a `SUPERGROUPS`
  membership under `splicing`, a `naming.R` rule, and CSV row 46. It is the only group
  that is empty in both species.
- **Four dead `naming.R` rules** — patterns matching no column in either cache:
  `^exon_length_first_`, `^noncoding_length_fraction_`, `^saluki_prediction$`,
  `^prediction_difference$`. The last two are dead because of S2 and should be revived,
  not deleted; the first two have no data source at all.
- **`mfe_expected` is fully excluded but still in the default selection.** All 7 columns
  are in `EXCLUDED_FEATURES`, yet `structure_core` pulls in the whole `structure`
  supergroup, so `INCLUDED_GROUPS` resolves to a group with zero surviving members.
  Harmless, but the default selection advertises a group that can never return anything.
- **`introns` (3 cols) and `junctions` (8 cols) are likewise fully excluded** —
  `keptH=0` for both. Unlike `mfe_expected` they are not reachable from `INCLUDED_GROUPS`,
  so this is only a documentation point: three of the 25 groups can never contribute a
  covariate.

## C7. Display names that do not round-trip through `format_col_name()`

`PIPELINE_GUIDE.md` §2.3 makes the round-trip a rule, so these are violations. I ran
every pooled column through `format_col_name()` and diffed against the CSV.

| CSV "Display Short name" | `format_col_name()` actually returns | Column |
|---|---|---|
| `C+G%` | `G+C%` | `gc_content_*` |
| `MFE.pred` | `MFE.pred.` | `mfe_expected_*` |
| `MFE.Δ` | `MFE.delta` | `mfe_delta_*` |
| `codon.<x>%` | `Codon.<X>%` | `codon_*_cds` |
| `stop-free` | `Stop-free` | `stopfree_length_*` |
| `EEJ.dens.` | `EEJ dens.` (space, not dot) | `exon_density_*` — and see S1 |
| `icSHAPE` | `icSHAPE.nuc` / `icSHAPE.cyto` | `gini_*` — one CSV row covers two distinct labels |
| `NMD.frag.` | `NMD frag.` | `nmd_snv_fragile_codon_density_mrna` |
| `NMD.alt-stop` | `NMD alt-stop` | `nmd_alt_stop_codon_density_mrna` |
| `uORF` | `uORF mRNA` | `uorf_present_mrna` |
| `min.local.MFE` | `min.local.MFE` ✓ | `rnalfold_score_*` |

The pattern is consistent: the CSV uses `.` as a separator where `naming.R` uses a space,
and differs on capitalisation. Whichever convention wins, it should win everywhere.

## C8. Stale column names in analysis scripts

Independent of the CSV, found while tracing feature usage:

- [analysis/qc/dataset_overview.R:362](analysis/qc/dataset_overview.R:362) references three
  columns that do not exist in either cache:

  | Referenced | Actual |
  |---|---|
  | `stopfree_cds` | `stopfree_length_cds` |
  | `intron_mean` | `intron_length_mean_mrna` |
  | `aa_freq_leu` | `aa_l_cds` |

- [analysis/qc/imputation_check.R](analysis/qc/imputation_check.R) requires
  `rnafold_{score,zscore}_mrna_imputed`. No `*_imputed` column exists in v10. The script
  is dead; its `need <- c(...)` guard means it fails cleanly rather than silently, but it
  cannot run.

---

# Part 3 — Master reconciliation

All 25 groups. `keptH`/`keptM` = columns surviving `drop_excluded()`.

| Group | Supergroup | H | M | keptH | keptM | CSV row(s) | Status |
|---|---|---|---|---|---|---|---|
| `lengths` | intrinsic | 7 | 7 | 4 | 4 | 10 | OK (CSV lists included regions only) |
| `gc` | intrinsic | 7 | 7 | 7 | 7 | 11 | C7 label |
| `nmd` | decay | 5 | 5 | 5 | 5 | 25, 26, 50–52 | OK |
| `introns` | splicing | 3 | 3 | **0** | **0** | 43 | C1, C6 |
| `exons` | splicing | 2 | 2 | **1** | **1** | 40, 41, 44, 45 | **S6 leak** |
| `noncoding` | splicing | **0** | **0** | 0 | 0 | 46 | **C6 empty group** |
| `rnafold_scores` | structure | 8 | **0** | 8 | **0** | 2 | **S3**; C3 utrpair |
| `rnafold_zscores` | structure | 8 | **0** | 8 | **0** | 3 | **S3**; C3 utrpair |
| `rnafold_per_nt` | structure | 7 | **0** | 7 | **0** | 4 | **S3** |
| `mfe_deltas` | structure | 7 | **0** | 7 | **0** | 6 | **S3**; C7 label |
| `mfe_expected` | structure | 7 | **0** | **0** | **0** | 5 | **S3**; C6 in default selection |
| `rnalfold_scores` | structure | 7 | **0** | 7 | **0** | 7 | **S3** |
| `rnalfold_zscores` | structure | 7 | **0** | 7 | **0** | 8 | **S3** |
| `junctions` | splicing | 8 | 8 | **0** | **0** | 23, 47 | **S1** |
| `eej_dist` | splicing | 6 | 6 | 2 | 2 | 24, 48, 49 | C1 group name |
| `uorfs` | translation | 3 | 3 | 1 | 1 | 21, 32, 34–37 | C3 region, C5 |
| `exon_density` | translation | 4 | 4 | 4 | 4 | 39 | **S1** |
| `stopfree` | intrinsic | 4 | **14** | 4 | **14** | 19 | **S4** |
| `skews` | intrinsic | 14 | 14 | 14 | 14 | 15, 16 | OK |
| `codon_freqs` | intrinsic | 64 | 64 | 64 | 64 | 17 | **S7** count |
| `aa_freqs` | intrinsic | 20 | 20 | 20 | 20 | 18 | OK |
| `nuc_ratios` | intrinsic | 28 | 28 | 28 | 28 | 12 | OK |
| `compositional` | intrinsic | 14 | 14 | 14 | 14 | 13, 14 | C1 group name |
| `probing` | structure | 8 | 8 | 8 | 8 | 9 | C3 compartment axis |
| `standalone` | other | 3 | **0** | 2 | **0** | 20, 22, 38 | **S3** |

Plus, outside any group: `halflife` (response, correct), 45 excluded ungrouped columns
(C4), and the 7 mouse-only `stopcodon_fraction_*` (S4).

---

# Part 4 — Proposed schema for the checked-in table

For approval, not built. The goal is a CSV a script can verify, so this audit never has
to be done by hand again.

## Columns to add

| Column | Why |
|---|---|
| `column_pattern` | The regex or literal column name. **The join key** — without it, nothing can be checked mechanically. |
| `fg_group` | Promoted out of free-text Notes into its own validated column (C1 exists only because it was buried in prose). |
| `status` | `included` / `built_excluded` / `not_built` — resolves the S5 ambiguity. |
| `species` | `both` / `human` / `mouse` — makes S3 and S4 visible in the table itself. |
| `n_cols_human`, `n_cols_mouse` | Catches S7-style count errors and S4-style asymmetries automatically. |

`Supergroup` should be restricted to the six `SUPERGROUPS` keys (C2), and
`Display Short name` should be the literal `format_col_name()` output (C7).

## What the checking script would assert

Run against a built cache, exits non-zero on any failure:

1. Every `column_pattern` matches ≥ 1 real column, per the row's `species`.
2. Every real column is claimed by exactly one row — nothing undocumented (this is what
   would have caught `stopcodon_fraction_*`).
3. `n_cols_human` / `n_cols_mouse` match the cache.
4. `status` agrees with `EXCLUDED_FEATURES` and with actual presence in the cache.
5. `fg_group` is a real `FEATURE_PATTERNS` key; `Supergroup` is a real `SUPERGROUPS` key
   and matches `supergroup_of(fg_group)`.
6. `Display Short name` equals `format_col_name()` on a representative column.
7. The three `config.R` invariants (mutual exclusivity, supergroup coverage, no empty
   groups) — all currently pass except the last, which fails on `noncoding`.

Suggested home: `docs/feature_table.csv` + `scripts/check_feature_table.R`, with the
script also runnable as a test.

---

# Suggested order of work

1. **S2** (Saluki join) and **S6** (`exon_length_last_mrna`) — both are defects with
   results consequences, and both are small fixes. S2 needs a `CACHE_VERSION` bump.
2. **S1**, **S3**, **S4**, **S5**, **S7** — the decisions. S3 in particular gates the
   table's shape, so it is worth settling early.
3. **C1–C7** — apply as a batch to the CSV once the decisions above are made, since
   several of them (C1 group names, C7 labels) depend on the outcome of S1.
4. **C8** — independent of everything else; fix whenever.
5. **Part 4** — build the checked-in table and the checking script last, so it encodes
   the resolved state rather than the current one.
