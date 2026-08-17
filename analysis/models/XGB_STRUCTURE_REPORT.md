> ## ⚠️ SUPERSEDED — do not quote the numbers below
>
> This report describes the **two-model design** (Model A vs Model B) that ran
> on 2026-08-12. That design has been replaced by a **four-rung nested ladder**
> — Baseline → S-core → S-select → S-full — and the baseline itself changed:
> the 20 `aa_*` amino-acid frequency columns were removed, because they are an
> exact deterministic function of the codon columns that are kept
> (`aa_x = sum(codons for x) / (1 - stop_fraction - codon_other)`, verified to
> 5.6e-17 on the v9 human cache). The baseline went from 127 predictors to 107.
>
> **Every number in this report was produced against the old baseline and is
> therefore not comparable to the current run.** The old variants map onto the
> new ladder as follows:
>
> | Old | New |
> |---|---|
> | `default` structure block (22 cols) | **S-select** rung |
> | `with_raw_mfe` (44 cols) | **S-full** rung |
> | `keep_missing` row policy | the new `default` variant's row policy |
> | — | **S-core** (15 cols) — new, and now the pre-specified primary |
>
> The prose about method, the paired-bootstrap sidebar, the importance trap and
> the caveats all still apply and carry over. The **result tables do not.**
> Rewrite this file from the new `SUMMARY.txt` and `_summary/` tables once the
> full-budget ladder has been run.

# Does secondary structure help predict mRNA half-life?

A short summary of the XGBoost test. Run 2026-08-12.

---

## The question

Does adding secondary-structure information improve prediction of mRNA half-life,
beyond what non-structure transcript features already give us?

Two models, identical in every way except one:

- **Model A** — baseline features only
- **Model B** — exactly Model A's features, **plus** the structure block

If Model B predicts held-out genes better, structure carries information the
baseline doesn't already have.

---

## The data

- **Species:** human
- **Response:** `halflife` — PC1 of the Agarwal & Kelley (2022) consensus half-life
  - This is a signed score, not hours. Range −17 to +18. Used untransformed.
- **Starting point:** 13,601 genes with a half-life value (one row per gene)
- **After dropping genes missing any feature:** **11,801 genes**
  - Same genes for both models. Neither model gets more data than the other.

---

## What went into each model

**Baseline — 127 features**

| Family | n | What it is |
|---|---|---|
| Codon frequencies | 64 | coding composition |
| Amino-acid frequencies | 20 | coding composition |
| GC/AT skews | 14 | strand asymmetry, per region |
| GC content | 7 | per region |
| NMD features | 5 | fragile codons, alternative stops |
| Lengths | 4 | 5'UTR, CDS, 3'UTR, whole mRNA |
| Stop-free lengths | 4 | per region |
| Exon density | 4 | per region |
| Junction distances | 2 | to start, to stop |
| uORF present | 1 | yes/no |
| Last-exon length | 1 | NMD-relevant |
| CAI | 1 | codon adaptation index |

**Structure block — 22 features** (the only difference)

- MFE z-scores, 8 regions
- Local MFE z-scores, 7 regions
- MFE delta (observed − expected), 7 regions

**Deliberately left out**

| Excluded | Why |
|---|---|
| Raw MFE scores, MFE per nucleotide | Almost entirely determined by length and GC, which are already in the baseline. A "win" could just be length sneaking back in. |
| icSHAPE Gini | 80–91% missing. Keeping it would cost 93% of our genes. Tested separately — see below. |
| Translation efficiency | A measured phenotype, not a sequence feature. Also 26% missing. |
| Individual nucleotide fractions, purine/amino ratios | Redundant with GC% and the skews. |

---

## The method

- **Train/test split:** the project's existing family-blocked 80/10/10 split
  - Genes are grouped into sequence families; a whole family stays in one split
  - So no test gene has a close paralogue in the training data
- **Tuning:** 60 hyperparameter settings, tested with 5-fold cross-validation
  inside the training data (folds also family-blocked)
  - Both models got the same 60 settings and the same folds
  - Each model picked its own winner — they chose differently, so neither was
    handicapped
- **Testing:** the held-out 1,165 genes were used **once**, at the very end
- **Uncertainty:** 2,000 paired bootstrap replicates
  - Same genes resampled for both models each time, so we get an interval on
    the *improvement*, not on two separate numbers
  - See the sidebar below for why the pairing matters

---

## Sidebar: why "paired" bootstrap

**The problem.** Model B is ahead by 0.0006. But our test set is one particular
sample of 1,165 genes. If a different 1,165 genes had landed there, would B still
be ahead?

**The bootstrap.** Simulate other test sets: draw 1,165 genes at random from our
1,165, with replacement. Score both models. Repeat 2,000 times.

**The "paired" bit.** In each draw, *both models are scored on the exact same gene
list*. This is the whole trick. Here is what our run actually produced:

| | mean | sd | 95% range |
|---|---|---|---|
| Model A R² | 0.4708 | 0.0224 | 0.4255 to 0.5140 |
| Model B R² | 0.4714 | 0.0226 | 0.4258 to 0.5148 |
| **Difference (B − A)** | **0.0006** | **0.0048** | **−0.0086 to 0.0101** |

Each model's R² swings across a range of about 0.09 depending on which genes get
drawn — some draws contain more easy genes, and *both* models score well. That
wobble is shared: the two models' scores correlate at **0.977** across draws.
Subtracting them within each draw cancels the shared luck, leaving a difference
5× less variable than either model alone.

**Why the obvious alternative fails.** You could put a CI on each model separately
and check for overlap. Ours overlap almost entirely — but that proves nothing,
because they would overlap almost entirely *even if B were genuinely better*.
Both are dominated by the same shared noise. Pairing removes it and asks the
sharper question: for the same genes, does B beat A?

**Analogy.** Two runners. Timing them on different days means wind differences
swamp the result. Racing them side by side cancels the wind. The bootstrap draws
are the races; pairing is running them side by side.

**What this buys us.** Because the pairing sharpened the measurement, our null is
informative. The method could have detected an R² improvement of about 0.01. It
found nothing. So the honest claim is "no improvement larger than ~0.01 in R²" —
not merely "we couldn't tell."

---

## Results

### Main comparison — 1,165 held-out genes

| Metric | Baseline | + Structure | Change | 95% CI |
|---|---|---|---|---|
| R² | 0.4717 | 0.4723 | +0.0006 | −0.0086 to +0.0101 |
| RMSE | 3.4992 | 3.4972 | −0.0020 | −0.0333 to +0.0285 |
| MAE | 2.7322 | 2.7495 | +0.0173 | −0.0131 to +0.0492 |

**Every interval crosses zero.** Structure made no measurable difference.

Two extra checks that say the same thing:

- Split the test set into 5 slices: 3 slices slightly favour structure, 2 favour
  baseline. No consistent direction. (See the next section.)
- Formal paired test: p = 0.90.

### How to read the slice plot

`xgb_structure_paired_slices` cuts the 1,165 held-out genes into **5 slices of 233
genes each**, then scores both models on each slice separately. One line per slice
— hence five lines. Slices are family-blocked, so a gene and its paralogue always
land in the same slice.

| Slice | n | Baseline R² | + Structure R² | Change | |
|---|---|---|---|---|---|
| 1 | 233 | 0.4960 | 0.4891 | −0.0069 | favours baseline |
| 2 | 233 | 0.5236 | 0.5096 | −0.0140 | favours baseline |
| 3 | 233 | 0.3735 | 0.3955 | +0.0220 | favours structure |
| 4 | 233 | 0.4497 | 0.4527 | +0.0030 | favours structure |
| 5 | 233 | 0.4728 | 0.4749 | +0.0020 | favours structure |

**Only the tilt of each line matters, not its height.** Lines sit at different
heights because some slices contain intrinsically harder genes — baseline R² ranges
from 0.37 to 0.52 across slices, a spread of 0.15. That is 250× larger than the
effect we are testing for, and it tells us nothing about structure.

- If structure genuinely helped, **all five lines would tilt the same way.**
- Ours tilt in both directions. The pooled near-zero is not a small consistent
  gain being diluted — it is genuine noise.

**What this plot is not.** It is a descriptive check on *where the pooled number
comes from*, not extra statistical evidence. All five slices come from the same
test set and the same two fitted models, so they are a partition of data we have
already used. The confidence intervals come from the bootstrap; this plot shows
whether the pooled result is broad-based or driven by one corner of the test set.

### The importance trap

Structure features take **8.1% of Model B's total "importance"** while contributing
**zero** predictive improvement. Their ranks run 47th to 125th out of 149.

**What gain importance actually measures.** Each time a tree splits on a feature,
XGBoost records how much that split reduced the loss *on the training data, at the
moment the split was made*. Those reductions are summed per feature and normalised
to 100%. So it is measured in-sample, it is a forced share (all 149 features carve
up exactly 100% whether or not any of them help), and it records that a feature was
*used*, not that it was *needed*.

**8.1% is not even a large share.** The 22 structure features are 14.8% of the 149
predictors, so proportionally they should take about 14.8%.

| | Features | % of features | % of gain | Gain per feature |
|---|---|---|---|---|
| Baseline | 127 | 85.2% | 91.9% | 0.724 |
| Structure | 22 | 14.8% | **8.1%** | **0.368** |

A structure feature earns **0.51×** what an average baseline feature earns. For
scale, `exon_density_cds` alone takes **23.6%** — 2.9× the entire structure block.
The top 10 features hold 41.3% of all gain.

**Why they earn anything at all — and this is the useful bit.** We tested whether
structure is simply a restatement of baseline features, by regressing each
structure feature on all 127 baseline features:

| Structure feature | R² from baseline | |
|---|---|---|
| `mfe_delta_cds` | 0.921 | baseline already knows this |
| `mfe_delta_5utr` | 0.917 | |
| `mfe_delta_mrna` | 0.908 | |
| … | | |
| `rnafold_zscore_5utr` | 0.058 | genuinely new information |
| `rnafold_zscore_stop` | 0.015 | |

The block splits almost exactly along family lines — **7 of 22 columns are more
than 60% reconstructible from baseline, and all 7 are the `mfe_delta_*` family:**

- The **`mfe_delta_*` features are substitutes.** Baseline reconstructs them at
  R² 0.68–0.92, which makes sense — delta is observed minus *expected* MFE, and
  expected MFE is a deterministic function of GC and length, both already in the
  baseline. Splitting on them collects gain for information the model already had.
- The **15 z-scores are genuinely novel** (R² 0.015 to 0.45, median about 0.13).
  The model really is seeing numbers it could not reconstruct. It uses them. They
  earn gain. They still buy nothing on held-out genes.

(These are linear R², so they are a *lower bound* on redundancy — a tree could
exploit a non-linear relationship the regression misses. That errs in the safe
direction: it can only understate how much the baseline already knows.)

That second case is the stronger lesson. The problem is not that structure
duplicates the baseline — it is that **a feature can be genuinely new, get used by
the model, accumulate importance, and simply not be related to half-life.** Novelty
is not relevance.

**If challenged.** Gain importance says a feature got used; the held-out comparison
asks whether the model predicts new genes better. When those disagree, held-out
wins. And if someone says 8.1% sounds substantial: 22 of 149 features should get
about 14.8% by proportion, so 8.1% is below-average participation, and one baseline
feature carries 2.9× the whole block.

*Numbers in this section come from `xgb_structure_gain_share.csv` and
`xgb_structure_redundancy.csv`, both written by the main script.*

### Sensitivity — does the answer depend on our choices?

Two of the design decisions above could reasonably have gone the other way, so
both were re-run as full specifications. Every variant is reported here, always,
in one table.

| Variant | Eligible genes | Held-out | Structure cols | ΔR² | 95% CI | |
|---|---|---|---|---|---|---|
| `default` | 11,801 | 1,165 | 22 | +0.0006 | −0.0086 to +0.0101 | no effect |
| `keep_missing` | 13,601 | 1,360 | 22 | +0.0008 | −0.0088 to +0.0107 | no effect |
| `with_raw_mfe` | 11,801 | 1,165 | 44 | −0.0050 | −0.0138 to +0.0044 | no effect |

- **`keep_missing`** drops the complete-case rule and lets XGBoost handle NA
  natively, recovering all 13,601 genes (13.2% of which have at least one missing
  predictor). This is the variant most likely to *favour* structure, because
  missingness in the structure block is informative. It doesn't: +0.0008.
- **`with_raw_mfe`** adds raw MFE and per-nucleotide MFE back, doubling the
  structure block to 44 columns. If excluding length-confounded MFE had hidden a
  signal, this is where it would appear. It doesn't — the point estimate goes
  slightly *negative*, consistent with 22 extra length-correlated columns adding
  noise rather than information.

**0 of 3 specifications** produce a ΔR² whose interval excludes zero.

> **Why this is reported as a block, not one at a time.** Run enough
> specifications and one will clear zero by chance, and it will be the one that
> looks most publishable. `xgb_structure_variant_summary.R` reads the variant
> registry rather than globbing the filesystem, so a defined-but-unrun variant
> shows as `NOT RUN` and deleting a directory cannot quietly remove an
> inconvenient row. The claim to make is "the conclusion holds across every
> specification tried", which is stronger than any single run.

Caveat when reading the table: within a variant the delta is paired and
meaningful, because both models saw the same genes. *Between* variants the
absolute R² values are not comparable when the row policy differs — `keep_missing`
scores on 1,360 genes, the others on 1,165. Compare directions, not levels.

### Secondary test — icSHAPE Gini

Run on the 861 genes that actually have icSHAPE data, using cross-validation.

| Comparison | Change in R² | 95% CI | Verdict |
|---|---|---|---|
| Computational structure vs baseline | −0.007 | −0.019 to +0.005 | No effect |
| **icSHAPE Gini, added on top** | **+0.029** | **+0.009 to +0.049** | **Real effect** |

All three metrics (R², RMSE, MAE) agree and all three clear zero.

So: **computed folding energy adds nothing; experimentally measured structure does.**

---

## What this means

**Fair to say:**

> Adding computationally predicted secondary-structure features did not improve
> XGBoost prediction of mRNA half-life beyond the non-structure baseline. On the
> subset with experimental probing data, icSHAPE structural Gini did add
> incremental predictive information.

**Not fair to say:**

> XGBoost proves secondary structure does / doesn't cause changes in half-life.

This is a prediction test, not a causal one. It tells us whether structure carries
information the baseline lacks — nothing about mechanism.

---

## Caveats to have ready

1. **The icSHAPE genes aren't a random sample.** Probing depth follows expression,
   so those 861 genes skew toward abundant transcripts.
2. **Gini is measured, not computed.** A model using it can't score a transcript
   nobody has probed. That's a different kind of model from A and B.
3. **Gini might be an abundance proxy.** It comes from read depth, read depth
   tracks abundance, and abundance relates to half-life. This design can't
   separate those. We'd need an expression covariate, which isn't in the dataset.
4. **The test set has no large gene families.** The split deliberately pins big
   families to training, so results measure generalisation to small and mid-sized
   families.
5. **13% of genes were dropped** to keep both models on identical rows — but the
   `keep_missing` variant puts them back and the answer is unchanged, so this
   caveat is now tested rather than merely acknowledged.
6. **A null result isn't proof of absence.** It means: at this sample size, with
   these features, no detectable improvement. A different structure
   representation could still work — though the two most obvious alternatives
   have now been tried and neither changed the answer.

---

## Files

### Scripts

| What | Where |
|---|---|
| Main comparison | `analysis/models/xgb_structure_comparison.R [variant]` |
| Figures only (~5 s) | `analysis/models/xgb_structure_plots.R` |
| Sensitivity table | `analysis/models/xgb_structure_variant_summary.R` |
| icSHAPE secondary | `analysis/models/xgb_structure_gini_subset.R` |
| Feature blocks + variant registry | `analysis/models/xgb_structure_features.R` |

### Outputs

One self-contained directory per run, under `data/outputs/xgb_structure/`:

```
xgb_structure/
├── default/       tables/  plots/  fits, folds, manifest, SUMMARY.txt
├── keep_missing/  same
├── with_raw_mfe/  same
├── gini/          the icSHAPE secondary run
└── _summary/      the cross-variant sensitivity table and figure
```

### How to re-run

```bash
Rscript analysis/models/xgb_structure_comparison.R              # default, ~40 min cold
Rscript analysis/models/xgb_structure_comparison.R keep_missing
Rscript analysis/models/xgb_structure_variant_summary.R         # seconds
Rscript analysis/models/xgb_structure_plots.R                   # figures only, ~5 s
```

Reproducible: seed 42, reruns give identical numbers. Fitted models are cached
per variant and keyed on the predictor lists, gene IDs, seed and tuning grid, so
a re-run takes about 70 seconds unless something that affects the models changed
— in which case it re-tunes on its own. `XGB_REFIT=1` forces it.

**Key figures for the talk** (in `data/outputs/xgb_structure/default/plots/`
unless noted)

1. `xgb_structure_delta_metrics` — the three changes with error bars (the null)
2. `xgb_structure_paired_slices` — improvement is not consistent across the test set
3. `_summary/plots/xgb_structure_variant_summary` — the null survives all three specifications
4. `gini/plots/xgb_structure_gini_deltas` — the icSHAPE result
5. `xgb_structure_modelB_importance` — the importance trap
