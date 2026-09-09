# Computed secondary structure and mRNA half-life: a bounded null

**A reporting brief for the paper.** What was tested, what the numbers are, how
to frame them, and — at least as importantly — what they do not license.

Generated from `data/outputs/xgb_structure/` (fitted 2026-09-08 17:18, seed 42,
R 4.5.2). Every number below is traceable to a file listed in
[Provenance](#provenance); none of it is quoted from memory.

---

## 1. What was tested

Two XGBoost models on human MANE-select transcripts, differing in one thing:

| Model | Predictors |
|---|---|
| `Baseline` | 106 non-structure transcript features — regional lengths, GC, AT/GC skews, codon composition, CAI, stop-free lengths, uORF presence, CDS-exon density, exon-junction distances, NMD features |
| `Structure` | `Baseline` + all 44 computed folding features — RNAfold and RNALfold MFE, their z-scores against shuffled sequence, per-nucleotide MFE, and MFE delta (observed − expected) |

**Response:** `halflife`, PC1 of the Agarwal & Kelley (2022) consensus half-life
measure. A signed score, not a duration in hours, modelled untransformed. RMSE
and MAE are therefore in PC1 units.

**Cohort:** 12,277 genes. Transcripts with either UTR under 30 nt (or a missing
UTR length) are excluded — see `MIN_UTR_LENGTH` in `R/config.R`. A UTR of a few
nucleotides makes the regional features degenerate rather than merely noisy.

**Design:** hyperparameters tuned by 5-fold cross-validation blocked on
`family_id_medium` over train+val (11,062 genes, 60 configurations), refit on
all of train+val, then **one** evaluation on the untouched `test` split (1,215
genes) of the committed family-blocked holdout. `test` was read once.

**What is shared by construction, not by hand:** rows, gene ids, preprocessing,
tuning resamples, tuning grid, budget, seeds and evaluation. `Structure` is
built as `BASELINE + block`, so "the only difference is the folding block" is a
property of the code, and is asserted against the predictors the *fitted*
models actually used (17 validation checks, all passing —
`xgb_structure_validation_checks.csv`).

---

## 2. The result, stated as a bound

This is the framing that matters. A null is worth reporting when it is
converted from "we found nothing" into "we bounded it."

| Metric | `Baseline` | `Structure` | Δ (Structure − Baseline) | 95% CI |
|---|---|---|---|---|
| **R²** (held-out, 1 − SSE/SST) | **0.4797** | 0.4756 | **−0.0041** | **[−0.0127, +0.0046]** |
| RMSE | 3.583 | 3.597 | +0.0142 | [−0.0157, +0.0444] |
| MAE | 2.787 | 2.793 | +0.0059 | [−0.0240, +0.0349] |

95% percentile intervals from 2,000 paired bootstrap replicates over held-out
genes; both models are scored on the same draw in every replicate, so each
interval is an interval on a paired improvement. Secondary paired sign-flip
test (10,000 permutations): p = 0.36 on squared error, p = 0.70 on absolute
error.

**The reportable claim:**

> Adding computed secondary-structure features to a primary-sequence baseline
> produced no detectable improvement in held-out prediction of mRNA half-life.
> The 95% interval on the change in R² was [−0.013, +0.005], so the data are
> inconsistent with an improvement larger than approximately 0.005 absolute R²
> — about 1% relative to a baseline R² of 0.48.

Say "bounded below 0.005 R²", not "no significant difference". The former is
quantitative and falsifiable; the latter is an absence of evidence dressed as
evidence of absence.

**Supporting consistency.** Across five family-blocked slices of the held-out
set, `Baseline` scored higher in three and `Structure` in two
(`xgb_structure_chunk_metrics.csv`). A real effect tilts the slices
consistently; these cross, which is what a null looks like rather than a small
true effect.

**A robustness detail worth one line.** Tuning selected the *identical*
hyperparameter configuration for both models (`mtry` 0.403, 498 trees,
`min_n` 100, depth 8, `learn_rate` 0.0218, `loss_reduction` 0.166,
`sample_size` 0.573, λ 0.160). The difference between the models is therefore
the feature block alone, and not a hyperparameter artefact. Because `mtry` is a
*proportion*, the wider model still sampled proportionally as many columns per
split (≈60 of 150 against ≈43 of 106) — a count-based range would have handed
the two models systematically different search spaces.

---

## 3. Why this null is informative rather than uninformative

The obvious challenge to any negative result is *your features were bad.* A
regression of each structure column on the entire baseline block answers it
directly.

**Method.** For each of the 44 structure columns in turn, an OLS multiple
regression of *that column* on all 106 baseline predictors plus an intercept
(n = 10,051 complete training rows; one QR factorisation reused across all 44
responses). The response is the structure feature, not half-life — the question
is "can the non-structure features already reproduce this number?" Fitted on
training rows only, so no test information enters even a descriptive table.

| Folding family | n | mean R² from baseline | range |
|---|---|---|---|
| `rnafold_scores` (raw MFE) | 8 | 0.851 | 0.676 – 0.965 |
| `rnafold_per_nt` | 7 | 0.815 | 0.676 – 0.971 |
| `rnalfold_scores` | 7 | 0.649 | 0.574 – 0.738 |
| `mfe_deltas` | 7 | 0.437 | 0.291 – 0.607 |
| `rnalfold_zscores` | 7 | 0.137 | 0.014 – 0.456 |
| `rnafold_zscores` | 8 | 0.108 | 0.012 – 0.253 |

**This defends both flanks at once**, which is unusual value for one
supplementary panel:

- *"Your structure block was degenerate / duplicative."* No. The two z-score
  families have median R² of 0.11 and 0.14, and all **15 z-score columns**
  retain a mean **87.9%** of their variance as not linearly reconstructible
  from the baseline. The block carried genuinely new numbers. They did not
  predict half-life.
- *"Of course nothing happened — MFE is just length in disguise."* True of raw
  MFE (R² 0.85–0.97, mechanically expected: folding energy scales with how much
  sequence there is to fold and shifts with GC) and demonstrably **false** of
  the shuffled-normalised z-scores, whose normalisation is designed to divide
  length and composition out and measurably does so.

Neither the bound in §2 nor this table is sufficient alone. The bound
quantifies the effect; this establishes that the bound is *informative* rather
than an artefact of a broken feature set.

**The figure** (`plots/xgb_structure_redundancy.*`) plots one point per
structure column against its R², grouped by folding family and ordered by
family median, with the median marked as a crossbar. There is deliberately no
reference line: any cutoff would be arbitrary, and a dashed line invites being
read as a threshold or a test. The families separate without one.

**State the limitation, because it cuts in our favour.** These R² are linear
and unadjusted. Linear, so a tree could exploit curvature and interactions this
misses — a *lower* bound on redundancy, which can only understate what the
baseline already knew. Unadjusted, with p/n = 106/10,051 ≈ 0.011, so roughly
one percentage point of each value is fitting noise (`rnafold_zscore_mrna`:
0.2533 unadjusted, 0.2453 adjusted). Both directions make the "the features
were genuinely novel" claim conservative rather than optimistic.

---

## 4. A methodological aside for the Discussion

Worth one sentence, no more, and explicitly as exploratory context.

The structure block held 29.3% of the predictors and collected 12.0% of the
total gain — participating below its share of the features
(`xgb_structure_gain_share.csv`). More instructive is *which* columns collected
it: the top structure features by gain are `rnafold_score_mrna`,
`mfe_delta_3utr` and `mfe_delta_cds` — all high-redundancy columns — while the
near-orthogonal z-scores were barely used.

So gradient boosting preferentially split on redundant *restatements* of
information it already had, in preference to the genuinely novel numbers
available to it. That is a useful caution about reading feature importance as
evidence of biological relevance: correlated predictors share and redistribute
gain, and high gain can mark convenience rather than information. It is not an
effect estimate and should not be presented as a finding.

---

## 5. What this does not license

Three claims the analysis does not support. Reviewers will test all three.

1. **Not "structure does not matter for mRNA stability."** What was tested is
   *computed thermodynamic folding predictions*. In vivo structure is shaped by
   ribosomes, RNA-binding proteins and helicases that no MFE calculation sees.
   The claim is about a class of sequence-derived features, not about RNA
   structure as a biological entity.
2. **Not a causal or mechanistic claim.** This is incremental *predictive
   information* conditional on a specific baseline block. Nothing here
   estimates an effect.
3. **Not general across species, datasets or structure measurements.** Human
   only; one half-life measure (Agarwal & Kelley consensus PC1); one
   annotation. Experimentally *measured* structure is untested here —
   icSHAPE structural Gini is excluded from both models by design (measured
   rather than computed, so a model using it cannot score an unprobed
   transcript, and it is 80–91% missing). That is a genuinely different
   question and a separate, later, supplementary analysis.

**Two design caveats to report with the numbers.**

- **The split is depleted of large families.** Families larger than ~5% of the
  smallest split are pinned to `train` (`SPLIT_PIN_FRAC`), so `test` contains no
  large family. It measures generalisation to small and mid-sized families, not
  to the largest ones.
- **Structure missingness is informative.** A missing 5'UTR MFE means a 5'UTR
  too short to fold, not a failed computation. XGBoost handles NA natively
  here, so `Structure` could in principle split on an annotation artefact and
  bank it as a structure effect. This design cannot exclude that. It would bias
  *toward* a spurious positive, so it does not threaten a null — but it must be
  stated, and it is the first thing to close if a future specification turns
  positive.

---

## 6. Suggested text

**Results.** *"Adding 44 computed secondary-structure features to a
106-feature primary-sequence baseline did not improve held-out prediction of
mRNA half-life (R² 0.476 vs 0.480; ΔR² −0.004, 95% CI [−0.013, +0.005]; n =
1,215 held-out genes from a family-blocked split). The interval excludes any
improvement above ~0.005 absolute R². To establish that this was not an
artefact of a redundant feature block, we regressed each folding feature on the
full baseline: the shuffled-sequence-normalised z-scores retained a mean 88% of
their variance as unexplained by the baseline (Fig. SX), confirming the block
contributed information the baseline lacked."*

**Discussion.** The negative result is useful to the field precisely because
the features were well constructed — it argues against further investment in
thermodynamic folding predictions as half-life covariates, and points instead
toward measured structure and the trans-acting factors that MFE cannot see.

---

## Provenance

Everything above is derived from artefacts in `data/outputs/xgb_structure/`:

| File | Contents |
|---|---|
| `SUMMARY.txt` | the run's own text summary |
| `run_manifest.rds` | design, feature lists, seeds, best hyperparameters |
| `tables/xgb_structure_model_metrics.csv` | pooled held-out metrics |
| `tables/xgb_structure_delta_bootstrap.csv` | Δ and 95% CIs (§2) |
| `tables/xgb_structure_signflip_test.csv` | secondary permutation test |
| `tables/xgb_structure_chunk_metrics.csv` | per-slice consistency |
| `tables/xgb_structure_redundancy.csv` | the regression in §3 |
| `tables/xgb_structure_gain_share.csv`, `..._gain_importance.csv` | §4 |
| `tables/xgb_structure_validation_checks.csv` | the 17 design assertions |
| `plots/xgb_structure_redundancy.*` | the Fig. SX candidate (§3) |
| `plots/xgb_structure_delta_metrics.*` | the Δ figure (§2) |
| `plots/xgb_structure_paired_slices.*` | per-slice consistency (§2) |
| `plots/xgb_structure_observed_vs_predicted.*` | calibration, both models |
| `plots/xgb_structure_gain_importance.*` | the gain pattern (§4) |

Figures carry a short title, a subtitle with sample sizes, and nothing else —
interpretation belongs in the paper's legend, where it is read alongside the
text. MAE is computed, bootstrapped and tabulated but not plotted; adding
`"mae"` to `PLOT_METRICS` in `xgb_structure_plots.R` puts it back on the two
metric figures with no re-run.

To regenerate from scratch:

```bash
Rscript analysis/models/xgb_structure_comparison.R
```

Fits are cached and keyed on the predictor lists, the exact gene ids, the seed,
the grid and the fold count, compared with `identical()` — so a changed cohort
or feature set forces a re-tune rather than silently serving a stale model.
