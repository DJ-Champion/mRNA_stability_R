# Does secondary structure help predict mRNA half-life?

A four-rung nested ladder, fitted with XGBoost. Human, run 2026-08-17, seed 42.

**Bottom line.** No. The pre-specified contrast — adding folding z-scores to a
107-feature non-structure baseline — moved held-out R² by **−0.0034**, with a 95%
interval of **−0.0140 to +0.0067**. Every rung of the ladder and every metric
tells the same story, under both row policies tested: **nothing favours
structure anywhere**. The method could have detected an R² gain of about 0.01;
it found nothing, so the honest claim is *no improvement larger than ~0.01 in
R²*, not merely *we couldn't tell*. (One secondary contrast in the sensitivity
arm does clear zero — in the direction favouring the baseline. It is unadjusted,
does not replicate, and is dissected in the sensitivity section.)

**The interesting exception.** On the 861 genes with experimental icSHAPE
probing, *measured* structural Gini does add predictive information — ΔR²
**+0.026**, CI **+0.007 to +0.045** — while computed folding on those same genes
adds nothing. Computed and measured structure are not interchangeable here. See
the secondary test, and the three caveats that come with it.

---

## The question

Does adding computationally predicted secondary-structure information improve
prediction of measured mRNA half-life, beyond what non-structure transcript
features already give us?

---

## The design: a nested ladder

Four models. Every rung is **the identical baseline block plus its own structure
block**, so a difference between two rungs is attributable to the features that
were added and to nothing else.

| Rung | Structure block added | Structure cols | Predictors | |
|---|---|---|---|---|
| **Baseline** | — | 0 | 107 | reference |
| **S-core** | MFE z-scores (8) + local MFE z-scores (7) | 15 | 122 | **PRIMARY** |
| **S-select** | + MFE delta, observed − expected (7) | 22 | 129 | secondary |
| **S-full** | + raw MFE (15) and per-nucleotide MFE (7) | 44 | 151 | secondary |

**Why this order, and why it is the whole point.** The rungs run from *least* to
*most* confounded with the baseline, so the ladder walks from the cleanest test
of structure to the dirtiest:

- **z-scores** are normalised against shuffled sequence, so they are the closest
  thing available to "structure, holding composition constant".
- **MFE deltas** are observed minus *expected* MFE, and expected MFE is largely a
  function of GC and length — both already in the baseline.
- **raw MFE** is close to a deterministic function of length and GC.

That ordering gives the shape a diagnostic meaning. A genuine structure effect
should appear at S-core and persist upward. An effect that appears *only* at
S-full is the expected signature of baseline information re-entering the model
wearing a structure label — not a discovery.

**Pre-specification.** S-core vs Baseline is the single pre-specified contrast.
The ladder generates 5 contrasts on one held-out set; the other 4 are secondary
and unadjusted, and exist to show the *shape*, not to be promoted to the
headline after the fact.

---

## The data

- **Species:** human, cache `data/cache/human_dataset_v9.rds`
- **Response:** `halflife` — PC1 of the Agarwal & Kelley (2022) consensus
  half-life. A signed score, not hours (range −17 to +18, sd 4.99). Used
  untransformed: no project code transforms it, and a log is undefined on a
  variable that takes negative values. **RMSE and MAE are therefore in
  half-life score units.**
- **Eligible:** **13,601** genes — every gene with a target and a split
- **Split:** the committed family-blocked holdout, `data/splits/holdout_medium.rds`
  - **12,241** train+validate / **1,360** held out
  - Genes are grouped into sequence families and a whole family stays in one
    split, so no test gene has a close paralogue in training
- **Missing data:** handled natively — XGBoost learns a default split direction
  per feature. No imputation, no dropped rows, and every rung sees **identical
  rows** so no rung can gain from a different gene set.

### Excluded from the baseline, and why

| Excluded | Why |
|---|---|
| `aa_*` (20 cols) | An exact deterministic function of the retained codon columns: `aa_x = sum(codons for x) / (1 − stops − other)`. Verified to 5.6e-17 on the v9 cache. |
| `frac_*`, `purine_*`, `amino_*` | Exact functions of GC content and the two skews, all retained. |
| `translation_efficiency` | A measured phenotype, not a sequence feature. |
| `gini_*` (icSHAPE) | 80–91% missing. See the secondary run. |

Removing the redundant blocks took the baseline from 127 predictors to **107**.
This matters for interpretation: a baseline padded with exact restatements of
its own columns makes the structure block look proportionally more important
than it is.

---

## The method

- **Tuning:** 60 hyperparameter settings, drawn once from seed 42, evaluated by
  5-fold cross-validation *inside the training data* (folds also family-blocked).
  Every rung got the same 60 settings and the same folds, and each picked its
  own winner — they chose differently, so no rung was handicapped.
- **Testing:** the 1,360 held-out genes were scored **once**, at the very end.
- **Uncertainty:** 2,000 **paired** bootstrap replicates over the held-out genes.
- **Secondary test:** a sign-flip (permutation) test on per-gene losses.

All 16 automated validation checks pass, including *no structure variable in the
reference model*, *the ladder is nested in the fitted models*, *no family spans
the fitting pool and the test set*, and *metrics come from held-out predictions*.
See `tables/xgb_structure_validation_checks.csv`.

### Sidebar: why "paired" bootstrap

**The problem.** S-core is behind Baseline by 0.0034. But our test set is one
particular sample of 1,360 genes. If a different 1,360 genes had landed there,
would the ordering hold?

**The bootstrap.** Simulate other test sets: draw 1,360 genes at random from our
1,360, with replacement. Score every rung. Repeat 2,000 times.

**The "paired" bit.** In each draw, *every rung is scored on the exact same gene
list*. This is the whole trick. In this run:

| | sd across draws | 95% range |
|---|---|---|
| Baseline R² alone | 0.019 | 0.461 to 0.537 |
| S-core R² alone | 0.020 | — |
| **Difference (S-core − Baseline)** | **0.005** | **−0.0140 to +0.0067** |

Either rung's R² swings across a range of about 0.076 depending on which genes
get drawn — some draws contain more easy genes and *both* rungs score well. That
wobble is shared: the two rungs' scores correlate at **0.96** across draws.
Subtracting within each draw cancels the shared luck, leaving a difference
roughly **3.7× less variable** than either rung alone.

**Why the obvious alternative fails.** You could put a CI on each rung
separately and check for overlap. Ours overlap almost entirely — but that proves
nothing, because they would overlap almost entirely *even if one rung were
genuinely better*. Both are dominated by the same shared noise. Pairing removes
it and asks the sharper question: *for the same genes, does the larger model
win?*

**Analogy.** Two runners. Timing them on different days means wind differences
swamp the result. Racing them side by side cancels the wind. The bootstrap draws
are the races; pairing is running them side by side.

**What this buys us.** Because pairing sharpened the measurement, the null is
*informative*. The interval is about ±0.010 wide on either side of zero, so an
R² improvement of that size would have shown up. It didn't.

*(The per-rung figures in the table above are reproducible from
`tables/xgb_structure_test_predictions.csv`; the delta row is read straight from
`tables/xgb_structure_delta_bootstrap.csv`.)*

---

## Results

### Primary contrast — S-core vs Baseline, 1,360 held-out genes

| Metric | Baseline | S-core | Change | 95% CI | |
|---|---|---|---|---|---|
| **R²** | 0.5003 | 0.4970 | **−0.0034** | −0.0140 to +0.0067 | no effect |
| RMSE | 3.5239 | 3.5357 | +0.0119 | −0.0235 to +0.0492 | no effect |
| MAE | 2.7553 | 2.7528 | −0.0025 | −0.0339 to +0.0306 | no effect |

Sign convention: delta = larger model − smaller model. For R², positive favours
the larger model; for RMSE and MAE, negative does. Here the R² point estimate is
*negative* and the RMSE point estimate is *positive* — both nominally favouring
the baseline — but all three intervals span zero, so the reading is "no
detectable difference", not "structure hurts".

Sign-flip test on the primary contrast: **p = 0.53** (squared error), **p = 0.88**
(absolute error). Same conclusion by a different route.

### Secondary rows — the shape of the ladder

| Contrast | ΔR² | 95% CI |
|---|---|---|
| S-core vs Baseline **(primary)** | −0.0034 | −0.0140 to +0.0067 |
| S-select vs Baseline | −0.0027 | −0.0146 to +0.0078 |
| S-full vs Baseline | −0.0022 | −0.0140 to +0.0086 |
| S-select vs S-core *(increment)* | +0.0006 | −0.0084 to +0.0101 |
| S-full vs S-select *(increment)* | +0.0005 | −0.0085 to +0.0088 |

**0 of 5 contrasts clear zero, on any of the three metrics.**

Read the shape, not the rows. Held-out R² is flat along the entire ladder —
0.5003, 0.4970, 0.4976, 0.4981 — drifting by less than 0.004 while the structure
block grows from 0 to 44 columns. There is no rung at which structure starts
paying, and the tiny upward creep from S-core to S-full is exactly what the
ordering predicts for baseline information re-entering under a structure label:
it arrives with the *most* confounded blocks, not the cleanest one.

### How to read the slice plot

`xgb_structure_paired_slices` cuts the 1,360 held-out genes into **5 slices of
272**, then scores every rung on each slice separately — one line per slice.
Slices are family-blocked, so a gene and its paralogue always land together.

| Slice | n | Baseline R² | S-core R² | Change | |
|---|---|---|---|---|---|
| 1 | 272 | 0.5406 | 0.5265 | −0.0141 | favours baseline |
| 2 | 272 | 0.4816 | 0.4781 | −0.0035 | favours baseline |
| 3 | 272 | 0.5063 | 0.5082 | +0.0020 | favours structure |
| 4 | 272 | 0.4831 | 0.4788 | −0.0044 | favours baseline |
| 5 | 272 | 0.4660 | 0.4714 | +0.0054 | favours structure |

**Only the tilt of each line matters, not its height.** Lines sit at different
heights because some slices contain intrinsically harder genes — baseline R²
ranges from 0.466 to 0.541, a spread of 0.075. That is about 22× larger than the
effect being tested, and it says nothing about structure.

- If structure genuinely helped, **all five lines would tilt the same way.**
- Ours tilt in both directions, 3 to 2. The pooled near-zero is not a small
  consistent gain being diluted — it is noise.

**What this plot is not.** It is a descriptive check on *where the pooled number
comes from*, not extra statistical evidence. All five slices come from the same
test set and the same fitted models, so they are a partition of data already
used. The intervals come from the bootstrap; this plot only shows whether the
pooled result is broad-based or driven by one corner of the test set.

### How to read the observed-vs-predicted plot

`xgb_structure_observed_vs_predicted` plots the held-out genes for every rung,
one panel each, on identical axes with a dashed y = x line (not a fit).

**This figure is secondary to the paired comparison, and cannot replace it.**
Each panel scores one model on its own, so reading the ladder off four panels is
exactly the unpaired comparison the sidebar above rejects: the panel-to-panel
differences are swamped by the shared luck of which genes landed in the test
set. The evidence about structure is the paired delta and its CI. This figure
answers a different question — *is the model any good at all, and how does it
fail* — and it applies equally to every rung.

**The clouds are squashed vertically, and that is correct.** Predictions span
noticeably less than observations, so the cloud is flatter than the dashed line.
This is calibrated shrinkage, not a defect: a well-calibrated model shrinks its
predictions toward the mean by roughly a factor of r, so sd(predicted) should be
about r × sd(observed). In this run r = 0.71 and the ratio is 0.67 on every
rung, and regressing observed on predicted gives a slope of 1.05 — the same fact
stated the other way round. A model that did *not* shrink would be over-confident
on the extremes and would score worse.

The consequence is worth saying out loud when the figure is on a slide: **the
model never predicts the most extreme half-lives**, and the shortest- and
longest-lived transcripts will always be pulled toward the middle. That is a
property of squared-error regression at this level of r, not of the structure
features, and it is unchanged along the ladder.

### The importance trap

At the primary rung, the 15 structure features take **6.7% of S-core's total
"importance"** while contributing zero predictive improvement. Their ranks run
**37th to 101st out of 122**.

**What gain importance actually measures.** Each time a tree splits on a
feature, XGBoost records how much that split reduced the loss *on the training
data, at the moment the split was made*. Those reductions are summed per feature
and normalised to 100%. So it is measured in-sample, it is a forced share (all
122 features carve up exactly 100% whether or not any of them help), and it
records that a feature was *used*, not that it was *needed*.

**6.7% is not even a large share.** The 15 structure features are 12.3% of the
122 predictors, so proportionally they should take about 12.3%.

| | Features | % of features | % of gain | Gain per feature |
|---|---|---|---|---|
| Baseline | 107 | 87.7% | 93.3% | 0.872 |
| Structure | 15 | 12.3% | **6.7%** | **0.445** |

A structure feature earns **0.51×** what an average baseline feature earns. For
scale, `exon_density_cds` alone takes **26.7%** — 4× the entire structure block.
The top 10 features hold 43.7% of all gain. The pattern holds up the ladder and
gets worse: at S-full the structure block is 29.1% of the features and takes
13.2% of the gain, for 0.37× per feature.

**Why they earn anything at all — and this is the useful bit.** We tested
whether structure is simply a restatement of baseline features, by regressing
each structure column on all 107 baseline columns:

| Block | First appears | R² from baseline | Reading |
|---|---|---|---|
| z-scores (15) | S-core | 0.012 – 0.449, median 0.087 | **genuinely new information** |
| MFE deltas (7) | S-select | 0.676 – 0.917 | largely a restatement |
| raw / per-nt MFE (22) | S-full | 0.570 – 0.971 | largely a restatement |

Overall, **27 of 44** structure columns are more than 60% reconstructible from
the baseline — and every one of those 27 sits in the S-select or S-full blocks.
The split runs exactly along the ladder:

- The **MFE deltas and raw MFE are substitutes.** The baseline reconstructs them
  at R² 0.57–0.97, which makes sense: expected MFE is a deterministic function
  of GC and length, both already present. Splitting on them collects gain for
  information the model already had. This is why they are the *upper*, dirtier
  rungs.
- The **15 z-scores are genuinely novel** (R² 0.012 to 0.449, median 0.087). The
  model really is seeing numbers it could not reconstruct. It uses them. They
  earn gain. **They still buy nothing on held-out genes.**

(These are linear R², so they are a *lower bound* on redundancy — a tree could
exploit a non-linear relationship the regression misses. That errs in the safe
direction: it can only understate how much the baseline already knows.)

That second case is the stronger lesson, and the ladder was built to isolate it.
The problem is not that structure duplicates the baseline — it is that **a
feature can be genuinely new, get used by the model, accumulate importance, and
simply not be related to half-life.** Novelty is not relevance.

**If challenged.** Gain importance says a feature got used; the held-out
comparison asks whether the model predicts new genes better. When those
disagree, held-out wins. And if 6.7% sounds substantial: 15 of 122 features
should take about 12.3% by proportion, so it is below-average participation, and
one baseline feature carries 4× the whole block.

*Numbers in this section come from `xgb_structure_gain_share.csv` and
`xgb_structure_redundancy.csv`.*

### Secondary test — icSHAPE Gini

Run 2026-08-17 against this same 107-feature baseline, on the **861 genes that
actually have icSHAPE data** (7.3% of the complete-case eligible set). Because
the committed test split holds only ~98 Gini-complete genes, this design uses
5-fold cross-validation blocked on `family_id_medium`, with 4-fold inner tuning
over a 20-point grid, and reports out-of-fold predictions.

Three models: Baseline (107) → S-core (122, the primary rung's 15 z-scores) →
S-core + Gini (130, adding 8 measured Gini columns).

| Out-of-fold, 861 genes | R² | RMSE | MAE |
|---|---|---|---|
| Baseline | 0.4785 | 3.568 | 2.734 |
| S-core | 0.4805 | 3.561 | 2.741 |
| **S-core + icSHAPE Gini** | **0.5066** | **3.470** | **2.653** |

| Comparison | ΔR² | 95% CI | |
|---|---|---|---|
| S-core vs Baseline *(the ladder's primary, replicated on this subset)* | +0.0020 | −0.0070 to +0.0111 | no effect |
| **icSHAPE Gini increment** | **+0.0261** | **+0.0074 to +0.0452** | **real effect** |
| **S-core + Gini vs Baseline** | **+0.0281** | **+0.0084 to +0.0485** | **real effect** |

All three metrics agree, and all three clear zero for both Gini rows. The
computed-structure row does not, reproducing the main null on this subset.

So: **computed folding energy adds nothing; experimentally measured structure
does.** This replicates the pre-ladder result (which gave +0.029 against the old
127-column baseline) now that the baseline has been cleaned of its redundant
blocks — the increment shrank slightly, from +0.029 to +0.026, and survived.

**Read the three caveats below before this goes on a slide.** The effect is
real, but this design cannot say it is a *structure* effect.

---

## Sensitivity — does the answer depend on our choices?

The row policy could reasonably have gone the other way, so it was re-run as a
full specification. Both are reported here, always, in one table.

| Variant | Rows | Eligible | Held out | Primary ΔR² | 95% CI | |
|---|---|---|---|---|---|---|
| `default` | native NA | 13,601 | 1,360 | −0.0034 | −0.0140 to +0.0067 | no effect |
| `complete_case` | complete cases | 11,801 | 1,165 | −0.0039 | −0.0123 to +0.0046 | no effect |

**0 of 2 specifications** produce a primary ΔR² whose interval excludes zero.

`complete_case` is more than a routine robustness check. Structure missingness
is *informative* — a missing 5'UTR MFE means a 5'UTR too short to fold, not a
failed computation — so under native-NA handling a structure rung can split on
an annotation artefact and bank it as a structure effect. Complete cases close
that channel by construction, at a cost of ~13% of the corpus. It was the one
run that could have revealed the null as an artefact. It didn't: the primary
delta is if anything slightly more negative.

Both arms use byte-identical 107-column baselines and identical structure
blocks; only the rows differ. Caveat when reading the table: *within* a variant
the delta is paired and meaningful, because every rung saw the same genes.
*Between* variants the absolute R² values are not comparable — `complete_case`
scores 1,165 genes and `default` 1,360. Compare directions, not levels.

> **Why this is reported as a block, not one at a time.** Run enough
> specifications and one will clear zero by chance — and it will be the one that
> looks most publishable. `xgb_structure_variant_summary.R` reads the variant
> registry rather than globbing the filesystem, so a defined-but-unrun variant
> shows as `NOT RUN` and deleting a directory cannot quietly remove an
> inconvenient row.

### The one contrast that did clear zero — and why it is not a finding

In `complete_case`, **S-select vs Baseline** is conclusive on all three metrics,
in the direction favouring the **baseline**:

| Metric | Delta | 95% CI | Sign-flip p |
|---|---|---|---|
| R² | −0.0106 | −0.0188 to −0.0025 | 0.010 |
| RMSE | +0.0348 | +0.0084 to +0.0610 | 0.037 |
| MAE | +0.0284 | +0.0004 to +0.0544 | — |

Stated plainly: on complete cases, adding the 7 `mfe_delta` columns made the
model measurably *worse* than the baseline. Before anyone builds a story on
that, four reasons it should not be treated as a result:

1. **It is 1 of 8 secondary, unadjusted contrasts** across the two variants. At
   the 5% level you expect about 0.4 such rows by chance; finding one is
   unremarkable. It is not the pre-specified contrast.
2. **It does not replicate in `default`**, where the same comparison gives
   −0.0027 with an interval spanning zero.
3. **It does not persist up the ladder.** S-full *contains* every `mfe_delta`
   column, and S-full vs Baseline is −0.0032, not conclusive. A block that
   genuinely degraded the model should keep degrading it when it stays in.
   The non-monotonicity is the tell.
4. **The hyperparameters point at tuning variance.** Each rung tunes
   independently by racing, so a rung can land on a poor configuration. In this
   run S-core and S-full both selected config `mod16` (498 trees, λ = 0.16,
   mtry = 0.40) while S-select selected `mod54` (345 trees, λ = 45.8,
   mtry = 0.92) — a substantially different, more heavily regularised model.
   The S-select dip tracks that choice more closely than it tracks the seven
   columns that were added.

The three metrics agreeing is *not* three independent confirmations — they are
computed from the same predictions on the same genes and are strongly
correlated. The sign-flip p-values are likewise on the same predictions.

**What would settle it:** refit S-select under S-core's winning configuration
and see whether the gap survives. If it vanishes, it was tuning variance; if it
holds, the `mfe_delta` block is genuinely harmful on complete cases and deserves
a proper look. That refit has not been run.

Note also what this row is *not* evidence for: it points the wrong way for the
structure hypothesis. Nothing in either variant favours structure at any rung.

---

## What this means

**Fair to say:**

> Adding computationally predicted secondary-structure features did not improve
> XGBoost prediction of mRNA half-life beyond a 107-feature non-structure
> baseline. This held at every rung of a nested ladder running from
> composition-normalised folding z-scores to the full raw-MFE block, on a
> family-blocked held-out set of 1,360 genes, with an interval tight enough to
> exclude improvements larger than about 0.01 in R², and it held under both
> row policies tested. On the subset of 861 genes
> with experimental icSHAPE probing, measured structural Gini *did* add
> incremental predictive information (ΔR² +0.026, CI +0.007 to +0.045) where
> computed folding on those same genes did not.

**Not fair to say:**

> XGBoost proves secondary structure does / doesn't cause changes in half-life.

This is a prediction test, not a causal one. It tells us whether these
representations of structure carry information the baseline lacks — nothing
about mechanism. Causal and inferential claims belong to the later modelling
plan.

> Measured structure beats computed structure, therefore folding algorithms are
> the weak link.

Tempting, and it may even be true, but this design cannot support it. The Gini
increment is measured on a non-random 7% of the corpus, and Gini is derived from
probing read depth, which tracks abundance, which is itself associated with
half-life. An abundance proxy would produce exactly this result. Caveats 5–7
below are the ones to have ready.

---

## Caveats to have ready

1. **A null result isn't proof of absence.** It means: at this sample size, with
   these features, no detectable improvement. A different structure
   representation could still work.
2. **The test set has no large gene families.** The blocked split pins families
   larger than ~5% of the smallest split to training, so this measures
   generalisation to small and mid-sized families, not to the largest ones.
3. **One conclusive secondary row exists**, and it favours the baseline
   (`complete_case`, S-select). It is unadjusted, does not replicate in
   `default`, and does not persist at S-full — see the sensitivity section. Have
   the explanation ready rather than being surprised by it.
4. **The response is a PC1 score, not hours.** RMSE and MAE are in score units
   and are not interpretable as a duration.
5. **The icSHAPE genes aren't a random sample.** Probing coverage tracks
   expression, so those 861 genes skew toward abundant transcripts and do not
   represent the corpus.
6. **Gini is measured, not computed.** A model using it cannot score a
   transcript nobody has probed — a fundamentally different kind of model from
   the ladder, which runs on sequence alone.
7. **Gini may be an abundance proxy.** It derives from probing read depth, read
   depth tracks abundance, and abundance relates to half-life. This design
   cannot separate those; it would need an expression covariate, which is not
   in the dataset.
8. **Predictions are shrunk toward the mean** by design, so the model never
   calls the most extreme half-lives. Applies to every rung equally.
9. **Gain importance is exploratory context**, not an effect estimate.
   Correlated predictors share and redistribute importance.

---

## Files

### Scripts

| What | Where |
|---|---|
| Main comparison (models + figures) | `analysis/models/xgb_structure_comparison.R [variant]` |
| Figures only (~5 s) | `analysis/models/xgb_structure_plots.R` |
| Sensitivity table across variants | `analysis/models/xgb_structure_variant_summary.R` |
| icSHAPE secondary | `analysis/models/xgb_structure_gini_subset.R` |
| Feature blocks, ladder + variant registry | `analysis/models/xgb_structure_features.R` |

### Outputs

One self-contained directory per run, under `data/outputs/xgb_structure/`:

```
xgb_structure/
├── default/        tables/  plots/  fits, folds, manifest, SUMMARY.txt
├── complete_case/  same, for the sensitivity arm
├── gini/           the icSHAPE secondary
└── _summary/       the cross-variant sensitivity table and figure
```

### How to re-run

```bash
Rscript analysis/models/xgb_structure_comparison.R              # default, ~40 min cold
```

```bash
Rscript analysis/models/xgb_structure_comparison.R complete_case
```

```bash
Rscript analysis/models/xgb_structure_plots.R                   # figures only, ~5 s
```

Reproducible: seed 42, re-runs give identical numbers. Fitted models are cached
per variant and keyed on the predictor lists, gene IDs, seed and tuning grid, so
a re-run takes about 70 seconds unless something that affects the models
changed — in which case it re-tunes on its own. `XGB_REFIT=1` forces it.

**Key figures for the talk** (in `data/outputs/xgb_structure/default/plots/`
unless noted):

1. `xgb_structure_delta_metrics` — each rung vs the baseline, with intervals: the null
2. `xgb_structure_paired_slices` — the near-zero is noise, not a diluted gain
3. `xgb_structure_primary_importance` — the importance trap
4. `gini/plots/xgb_structure_gini_deltas` — the icSHAPE secondary
5. `_summary/plots/xgb_structure_variant_summary` — the null survives both specifications
6. `xgb_structure_observed_vs_predicted` — model quality and calibration; secondary

In the delta figures, a **filled** point is one whose CI excludes zero — check
which side of the line it sits on. `complete_case` has one, and it favours the
baseline.

---

*Generated from `data/outputs/xgb_structure/` — `default/` (fitted 2026-08-17
17:47), `complete_case/` (22:55), `gini/` (22:12) and `_summary/`, via their
manifests, SUMMARY.txt files and tables. R 4.5.2, seed 42.*
