# Does secondary structure help predict mRNA half-life? — `default` variant

**Answer: no.** The pre-specified contrast is −0.0001 in held-out R², with a 95%
interval of [−0.0100, +0.0097]. That is not a small positive effect we were
underpowered to detect; it is zero to four decimal places, sitting in the middle
of a symmetric interval.

---

## Scope of this report

This report covers **one run and one run only**: the `default` variant, fitted
2026-08-18 13:48 against `data/cache/human_dataset_v10.rds`.

It is deliberately separate from `XGB_STRUCTURE_REPORT.md`, which documents the
older 107-predictor baseline and the pre-fix `mfe_delta` columns. Two things
changed underneath that report, and both change numbers it quotes:

1. **The MFE GC-units fix** (PR #17). `calculate_mfe_expected()` applied
   `gc / 100` to a column already stored as a fraction. Because the GC term is
   raised to `b = 2.352`, this erased it rather than shrinking it, collapsing
   the model to a constant −0.1534 kcal/mol/nt for every transcript.
   `mfe_expected_*` became a pure function of length and `mfe_delta_*`
   inherited that. Cache rebuilt at `CACHE_VERSION = 10`.
2. **`exon_length_last_mrna` dropped from the baseline** (PR #18). Last-exon
   length is a 3'UTR-length proxy — Spearman 0.949 with `length_3utr`, which is
   itself in the baseline. Baseline 107 → 106 predictors.

**Nothing here is comparable to the older report line-for-line.** Where a number
has moved for a reason worth knowing, it is flagged inline.

This report does **not** cover the `complete_case` variant, the icSHAPE `gini`
secondary run, the cross-variant summary, or the tuning diagnostic. Those
artefacts under `data/outputs/xgb_structure/` still reflect the 107-predictor
baseline and must be regenerated before they can be read alongside this.

---

## The question

> Does adding secondary-structure information improve prediction of measured
> mRNA half-life beyond a model containing non-structure transcript features?

**Response:** `halflife` — Agarwal & Kelley (2022) consensus half-life PC1,
untransformed. This is a signed score, not a duration in hours, so RMSE and MAE
below are in PC1 units and are not interpretable as time.

---

## The ladder

Four nested models. Each rung is the baseline block plus one more family of
folding metrics, so a difference between adjacent rungs is attributable to the
family that was added and nothing else.

| Rung | Structure cols | Predictors | Adds |
|---|---|---|---|
| Baseline | 0 | 106 | — |
| **S-core** ★ | 15 | 121 | MFE z-scores, local MFE z-scores |
| S-select | 22 | 128 | + MFE delta (observed − expected) |
| S-full | 44 | 150 | + raw MFE, per-nucleotide MFE |

★ = primary, pre-specified before any model was fitted.

The order is not arbitrary: rungs run from least to most confounded with the
baseline. A real structure effect should appear at **S-core** and persist. A
gain appearing only at **S-full** — whose block is near-deterministic in length
and GC — is the signature of baseline information re-entering under a structure
label, not a discovery.

### Excluded from the baseline

| Excluded | Why |
|---|---|
| `aa_*` (20) | exact function of retained codon columns: `aa_x = sum(codons for x)/(1−stops−other)` |
| `frac_*`, `purine_/amino_*` (42) | exact function of GC content and the two skews |
| `exon_length_last_mrna` (1) | 3'UTR-length proxy, Spearman 0.949 with `length_3utr` |
| `translation_efficiency` | measured phenotype, not a sequence feature; also 25.8% missing |
| `gini_*` (icSHAPE) | 80–91% missing — see the secondary run |

---

## Method

- **Eligible genes:** 13,601 (12,241 train+val / 1,360 test), one row per gene.
  Every rung sees identical rows; NAs handled natively by XGBoost.
- **Split:** committed family-blocked holdout, `data/splits/holdout_medium.rds`,
  blocked at `family_id_medium`.
- **Tuning:** 60 configurations drawn once from seed 42, raced over 5
  family-blocked inner folds. Same grid and same resampling object for every
  rung.
- **Inference:** 2,000 paired bootstrap replicates over held-out genes. Every
  rung is scored on the **same draw** in every replicate, so increments are
  within-draw differences and the intervals are intervals on increments.
- Seed 42, R 4.5.2.

---

## Results

### Primary contrast — S-core vs Baseline, 1,360 held-out genes

| Metric | Baseline | S-core | Delta | 95% CI |
|---|---|---|---|---|
| **R² (trad)** | 0.4976 | 0.4975 | **−0.0001** | **[−0.0100, +0.0097]** |
| RMSE | 3.5335 | 3.5339 | +0.0004 | [−0.0342, +0.0362] |
| MAE | 2.7606 | 2.7562 | −0.0044 | [−0.0387, +0.0289] |

Secondary sign-flip test: **p = 0.9832** (squared error), **p = 0.7983**
(absolute error).

Inconclusive, and about as cleanly null as this design can produce. The point
estimate is essentially zero and the interval is near-symmetric about it.

### Every contrast

| Contrast | Δ R² | 95% CI |
|---|---|---|
| **S-core vs Baseline** ★ | **−0.0001** | [−0.0100, +0.0097] |
| S-select vs Baseline | +0.0016 | [−0.0065, +0.0094] |
| S-full vs Baseline | −0.0032 | [−0.0116, +0.0052] |
| S-select vs S-core | +0.0017 | [−0.0077, +0.0117] |
| S-full vs S-select | −0.0048 | [−0.0121, +0.0026] |

**Every interval spans zero.** No contrast is conclusive in either direction.

Only the starred row is pre-specified. The other four are secondary and
unadjusted; none should be promoted to the headline after the fact. What they
are *for* is the shape of the ladder, not their individual intervals.

### Reading the shape

The ladder is flat, then slightly negative. S-full — the rung most confounded
with length and GC — is the *worst* performer at −0.0032 against baseline and
−0.0048 against S-select. Adding 22 more near-deterministic folding columns made
held-out prediction slightly worse, which is what you expect when the added
columns are redundant with information already present and simply dilute the
column subsample.

There is no rung at which structure helps, and the confounding signature the
design was built to detect does not appear either. The honest summary is that
the folding block is inert here.

---

## Why the structure block is inert: redundancy

Each structure column regressed on the baseline block. R² near 1 means the
column is a restatement of what the baseline already knows.

| Group | n | Min | Median | Max | First rung |
|---|---|---|---|---|---|
| `rnafold_zscore` | 8 | 0.012 | **0.071** | 0.253 | S-core |
| `rnalfold_zscore` | 7 | 0.013 | **0.120** | 0.449 | S-core |
| `mfe_delta` | 7 | 0.293 | 0.384 | 0.605 | S-select |
| `rnalfold_score` | 7 | 0.570 | 0.641 | 0.765 | S-full |
| `rnafold_per_nt` | 7 | 0.676 | 0.782 | 0.971 | S-full |
| `rnafold_score` | 8 | 0.676 | 0.912 | 0.965 | S-full |

Overall: 23 columns read as "genuinely new information", 21 as "largely a
restatement of baseline information".

**This is the key result of the run.** The z-scores in S-core are *not*
redundant — median R² 0.071 and 0.120 means they carry genuinely new numbers.
They were handed to the model as new information and still bought nothing. That
is a much stronger null than "the structure features were collinear with the
baseline, so of course they added nothing."

> **Changed since the old report.** `mfe_delta` was previously reconstructible
> from baseline at R² 0.68–0.92. It is now **0.293–0.605**. That is the direct
> consequence of the units fix: `mfe_expected` is now a genuine function of GC
> *and* length rather than length alone, so the residual is a real residual.
> The old report's stated rationale for expecting little from S-select — "expected
> MFE is a deterministic function of GC and length, both already in the baseline"
> — rested on a premise that was only accidentally true. S-select still adds
> nothing (+0.0016 [−0.0065, +0.0094]), but now for a different and better reason.

---

## Gain importance — exploratory context only

Correlated predictors share and redistribute importance, so this does not
estimate an independent effect. It is included because it corroborates the
redundancy story, not because it tests anything.

| Rung | Block | Features | % of features | % of gain | Gain per feature |
|---|---|---|---|---|---|
| S-core | baseline | 106 | 87.6% | 92.6% | 0.873 |
| S-core | structure | 15 | 12.4% | **7.4%** | **0.497** |
| S-select | structure | 22 | 17.2% | 9.7% | 0.440 |
| S-full | structure | 44 | 29.3% | 14.3% | 0.325 |

Structure features earn roughly **half** the gain per feature that baseline
features do, and the ratio falls as more of them are added.

In the primary rung, the highest-ranked structure feature is
`rnafold_zscore_5utr` at **rank 48 of 121**. The top of the table is
architecture and composition:

| Rank | Feature | Gain |
|---|---|---|
| 1 | `exon_density_cds` | 0.1421 |
| 2 | `eej_dist_closest_stop` | 0.0586 |
| 3 | `exon_density_mrna` | 0.0542 |
| 4 | `length_mrna` | 0.0227 |
| 5 | `codon_agg_cds` | 0.0142 |

---

## Validation

All **17** checks pass, including the two that guard the exclusions:

```
Same rows and identifiers for every rung                          TRUE
No structure variable appears in the reference model              TRUE
Every rung carries the identical baseline block                   TRUE   106 baseline columns
Every rung = baseline + its own structure block                   TRUE
The ladder is nested in the fitted models                         TRUE   Baseline < S-core < S-select < S-full
Each rung adds exactly its declared structure columns             TRUE
No amino-acid or nucleotide-fraction column is a predictor        TRUE
No last-exon-length column is a predictor                         TRUE   3'UTR-length proxy
No translation-efficiency variable in any rung                    TRUE
No identifier, family or outcome-derived column is a predictor    TRUE
Every test gene has a prediction from every rung                  TRUE
Test genes are disjoint from the tuning/fitting pool              TRUE
No family spans the fitting pool and the test set                 TRUE
Tuning used the same resampling object for every rung             TRUE   5 family-blocked folds
Tuning used the same grid for every rung                          TRUE   60 configurations, seed 42
Hyperparameter tuning used training data only                     TRUE
Metrics come from held-out predictions, not training predictions  TRUE
```

---

## Caveats

1. **Selected configurations differ across rungs.** Baseline, S-select and
   S-full all landed on `pre0_mod57_post0`; S-core landed on
   `pre0_mod12_post0`. Each rung tunes independently by racing, so a rung can
   land on a different configuration and carry that difference into its held-out
   score. Here the primary delta is −0.0001, so no visible artefact — but the
   confound is present and would matter if any contrast were near the decision
   boundary. **No tuning diagnostic has been run against this baseline.** The
   existing diagnostic is hardcoded to `complete_case` and its output is stale.

2. **Tuning variance may exceed the effects measured.** Under the previous
   baseline, the same feature set moved 0.0185 in held-out R² between two
   legitimately selected configurations — larger than any delta in this run.
   This is the strongest argument for reading the paired contrasts rather than
   individual R² values.

3. **Multiplicity.** Five contrasts on one held-out set; one pre-specified.
   The rest are unadjusted.

4. **The blocked split excludes large families.** Families larger than ~5% of
   the smallest split are pinned to train, so the test set contains no large
   family. This measures generalisation to small and mid-sized families only.

5. **The response is a PC1 score, not hours.** RMSE and MAE are in score units.

6. **Predictive evidence only.** No causal or inferential claim is made here.

7. **Structure is represented by computed folding metrics, not measurement.**
   The icSHAPE probing data is excluded from this run at 80–91% missing.
   "Structure adds nothing" here means *predicted* structure.

---

## What this means

A model of mRNA half-life built from transcript architecture and composition —
exon density, junction distances, lengths, codon usage, GC and skews — is not
improved by adding computed secondary-structure metrics, on 1,360 held-out
genes.

The result is strengthened rather than weakened by the redundancy analysis: the
z-score block was genuinely new information (median R² 0.07–0.12 from baseline)
and still bought nothing. The null is not an artefact of handing the model
information it already had.

It is also strengthened by the units fix. Under the bug, `mfe_delta` was largely
a length proxy, and a sceptic could have argued the ladder never really tested
delta at all. It now is a real residual, and the answer is unchanged.

---

## Files

### Scripts
```
analysis/models/xgb_structure_features.R      ladder, baseline, eligible set
analysis/models/xgb_structure_comparison.R    the run
analysis/models/xgb_structure_plots.R         figures
```

### Outputs — `data/outputs/xgb_structure/default/`
```
SUMMARY.txt                              generated summary
run_manifest.rds                         full provenance
tables/xgb_structure_model_metrics.csv   per-rung held-out metrics
tables/xgb_structure_delta_bootstrap.csv contrasts with bootstrap CIs
tables/xgb_structure_redundancy.csv      R2 of each structure col from baseline
tables/xgb_structure_gain_share.csv      block-level gain shares
tables/xgb_structure_gain_importance.csv per-feature gain
tables/xgb_structure_signflip_test.csv   secondary sign-flip test
tables/xgb_structure_chunk_metrics.csv   per-slice stability
tables/xgb_structure_validation_checks.csv
plots/xgb_structure_delta_metrics.*
plots/xgb_structure_observed_vs_predicted.*
plots/xgb_structure_paired_slices.*
plots/xgb_structure_primary_importance.*
```

### Re-run
```
Rscript analysis/models/xgb_structure_comparison.R default
```

### Related
- `XGB_STRUCTURE_REPORT.md` — the older multi-variant report. **Documents the
  107-predictor baseline and pre-fix `mfe_delta`; numbers do not carry over.**
- `XGB_ICSHAPE_REPORT.md` — icSHAPE Gini secondary run. Also predates both fixes.

---

*Generated from the run fitted 2026-08-18 13:48 against
`data/cache/human_dataset_v10.rds`. Seed 42. R 4.5.2.*
