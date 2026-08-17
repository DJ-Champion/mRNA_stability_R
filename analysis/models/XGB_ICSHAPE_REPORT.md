# Does *measured* structure help predict mRNA half-life?

The icSHAPE Gini secondary analysis. Human, run 2026-08-17, seed 42.

> **This is a SECONDARY result.** The headline analysis is the four-rung
> structure ladder in [`XGB_STRUCTURE_REPORT.md`](XGB_STRUCTURE_REPORT.md),
> which asks whether *computed* folding features help and answers no. Read that
> first. This report exists because the answer changes when the structure is
> measured rather than predicted — and because that change comes with three
> caveats heavy enough to need their own document.

**Bottom line.** On the 861 genes with icSHAPE probing data, adding 8 measured
structural Gini columns raises out-of-fold R² from **0.4805 to 0.5066** — a
paired increment of **+0.026 (95% CI +0.007 to +0.045)**, conclusive on all
three metrics. Computed folding on the *same* genes adds nothing (+0.002, CI
spans zero). So measured and computed structure are not interchangeable.

**And the caveat that has to travel with it.** These 861 genes are not a random
sample of the transcriptome, and the way they are selected is aligned with the
outcome: they are **10× over-represented among the most stable transcripts**
(26.8% of the corpus's top half-life decile, against 2.7% of the bottom). Gini
is computed from probing read depth, read depth tracks abundance, and abundance
tracks stability. This design cannot separate "measured structure carries
information" from "probing depth is an abundance proxy". The effect is real;
its interpretation is not settled.

---

## What this run asks

Three models, so it answers two questions at once:

| | Predictors | |
|---|---|---|
| **Baseline** | 107 | non-structure transcript features |
| **S-core** | 122 | + the primary rung's 15 folding z-scores |
| **S-core + icSHAPE Gini** | 130 | + 8 measured Gini columns |

The Gini block is `gini_{nucleoplasm,cytoplasm}_{mrna,5utr,cds,3utr}` — a Gini
coefficient of icSHAPE reactivity, per compartment and per region, measuring how
unevenly structured a transcript is rather than how structured it is on average.

Comparing Baseline → S-core replicates the main analysis on this subset.
Comparing S-core → S-core + Gini is the increment this run exists for.

---

## Why it is separate from the main analysis

**The missingness.** The Gini columns are 80–91% missing on the human v9 cache.
Only 861 of 11,801 complete-case-eligible genes have all eight. Putting Gini
into the main ladder would either cost 93% of the corpus or ask XGBoost to learn
a default split direction for a feature absent in seven of every eight genes.
Neither is a good trade for the primary question.

**Cross-validation, not the committed holdout.** The committed `test` split
contains only ~98 Gini-complete genes — too few for a usable interval on a
paired difference. So this run uses **5-fold cross-validation blocked on
`family_id_medium`** over all 861 genes, with **4-fold inner tuning** over a
20-point grid inside each outer training set. Every gene gets one out-of-fold
prediction per model, and no fold's assessment genes influenced their own
model's hyperparameters.

This is a real design difference from the main run, and it means **absolute R²
here is not comparable to absolute R² there.** Compare within this report only.

---

## Results

### Out-of-fold, 861 genes

| Model | R² | RMSE | MAE |
|---|---|---|---|
| Baseline | 0.4785 | 3.568 | 2.734 |
| S-core | 0.4805 | 3.561 | 2.741 |
| **S-core + icSHAPE Gini** | **0.5066** | **3.470** | **2.653** |

### Paired bootstrap, 2,000 replicates

| Comparison | Metric | Delta | 95% CI | |
|---|---|---|---|---|
| S-core vs Baseline | R² | +0.0020 | −0.0070 to +0.0111 | no effect |
| *(computed structure)* | RMSE | −0.0067 | −0.0377 to +0.0244 | no effect |
| | MAE | +0.0072 | −0.0224 to +0.0349 | no effect |
| **icSHAPE Gini increment** | **R²** | **+0.0261** | **+0.0074 to +0.0452** | **real** |
| | **RMSE** | **−0.0907** | **−0.1518 to −0.0261** | **real** |
| | **MAE** | **−0.0879** | **−0.1421 to −0.0315** | **real** |
| **Gini + structure vs Baseline** | **R²** | **+0.0281** | **+0.0084 to +0.0485** | **real** |
| | **RMSE** | **−0.0974** | **−0.1646 to −0.0303** | **real** |
| | **MAE** | **−0.0808** | **−0.1406 to −0.0192** | **real** |

Both Gini rows clear zero on all three metrics; the computed-structure row does
not, reproducing the main null on this subset.

**Effect size in plain terms.** The increment is 2.6 percentage points of
explained variance, a ~5% relative improvement over the baseline's 0.479. RMSE
falls by 0.091 half-life score units against a within-subset outcome sd of 4.94
— about 2% of a standard deviation. It is a real effect and a modest one.

The three metrics agreeing is not three independent confirmations: they are
computed from the same out-of-fold predictions on the same genes and are
strongly correlated. What the agreement does rule out is an artefact of one
metric's sensitivity to outliers.

---

## Who these 861 genes are

This is the part that determines what the result means, and it is worth more
than the assertion "probing tracks expression". Comparing the 861 Gini-complete
genes against the other 10,940 eligible genes:

| | Gini-complete (861) | Rest (10,940) | Std. diff |
|---|---|---|---|
| **Half-life score, mean** | **4.18** | **0.01** | **+0.88** |
| Half-life score, median | 4.67 | 0.12 | |
| Half-life 10th–90th pct | −2.03 to 9.95 | −5.88 to 5.81 | |
| Transcript length (mean, nt) | 2,717 | 4,047 | −0.60 |
| CDS length (mean, nt) | 1,218 | 1,909 | −0.52 |
| Exon density (CDS) | 8.23 | 6.76 | +0.50 |
| CAI | 0.769 | 0.768 | +0.01 |

**The selection is aligned with the outcome.** The probed genes sit almost a
full standard deviation higher in half-life. Split the corpus into deciles by
half-life and the Gini-complete share runs from **2.7% in the least stable
decile to 26.8% in the most stable** — against a base rate of 7.3%. Probing
coverage is roughly 10× denser among stable transcripts than unstable ones.

They are also **shorter** and more **exon-dense** than the corpus. CAI is the
one thing that does not move, but CAI is a weak proxy for expression and its
flatness should not be read as reassurance.

Two consequences:

1. **Generalisation is limited.** The increment describes well-probed, stable,
   short transcripts. There is no basis here for claiming it holds on the
   unstable transcripts that make up most of the corpus and most of the
   biological interest in decay.
2. **The abundance confound is not hypothetical.** The mechanism that makes a
   gene probable — high abundance — is the same mechanism that correlates with
   the outcome. A Gini column can carry that signal without carrying anything
   about structure.

One thing the selection does *not* do is inflate R² by compressing the outcome:
the within-subset half-life sd is 4.94 against 4.54 in the rest, so the subset
is if anything slightly *harder* to explain, not easier.

---

## The three caveats, in the order they will be raised

**1. Selection.** Covered above. The honest framing is "on well-probed
transcripts", never "on human mRNAs". n = 861 is also small, so the intervals
are wide — the R² increment's interval spans a factor of six, +0.007 to +0.045.

**2. Gini is measured, not computed.** Every other feature in this project is
derived from sequence. Gini is an experimental readout, so a model containing it
**cannot score a transcript nobody has probed**. This is a different kind of
model from the ladder: it does not answer "can we predict half-life from
sequence", it answers "does measured structure carry information about
half-life". Those are different claims and the second is much weaker as a
prediction tool, even though it is the more interesting one biologically.

**3. Confounding with abundance.** Gini is derived from probing read depth;
depth tracks abundance; abundance is associated with half-life. An apparent
Gini effect therefore has a live alternative explanation — that it is partly an
abundance proxy — and this design cannot rule it out. The selection analysis
above shows the confound is present in the data, not merely conceivable.

**What would settle it:** an expression covariate. Put transcript abundance in
the baseline and re-run the increment. If Gini survives, the structure reading
is much stronger; if it collapses, it was a depth proxy. The current cache does
not carry an expression column, so this is a data question, not a modelling one.
A weaker version — regressing Gini on read depth, if depth is recoverable from
the icSHAPE source data — would at least bound the confound.

---

## What this means

**Fair to say:**

> On 861 human transcripts with icSHAPE probing data, a measured index of
> structural unevenness (Gini) added out-of-fold predictive information about
> half-life beyond a 107-feature sequence baseline (ΔR² +0.026, 95% CI +0.007 to
> +0.045), while computationally predicted folding features on the same genes
> added nothing. The probed transcripts are markedly more stable and shorter
> than the corpus, and Gini derives from read depth, so an abundance-proxy
> explanation cannot be excluded.

**Not fair to say:**

> Measured structure beats computed structure, so folding algorithms are the
> weak link.

Tempting and possibly true, but unsupported. Three things differ between the two
arms at once: what is measured, which genes it is measured on, and whether the
signal could be an abundance proxy. Only a design that closes the third can
attribute the difference to measurement quality.

> icSHAPE Gini predicts mRNA half-life.

Not as stated. It contributes ~2.6 points of R² on a non-random 7% of the
corpus, in a model that cannot be run on unprobed transcripts.

---

## Files

| What | Where |
|---|---|
| Script | `analysis/models/xgb_structure_gini_subset.R` |
| Deltas | `data/outputs/xgb_structure/gini/tables/xgb_structure_gini_deltas.csv` |
| Out-of-fold predictions | `.../tables/xgb_structure_gini_oof_predictions.csv` |
| Manifest (design, caveats, n) | `.../gini_run_manifest.rds` |
| Figure | `.../plots/xgb_structure_gini_deltas` |

```bash
Rscript analysis/models/xgb_structure_gini_subset.R
```

About 80 minutes: 5 outer folds × 3 models × 20 configurations × 4 inner folds.
`XGB_GINI_GRID_SIZE=6` runs it end to end in a few minutes for a code check.

The subset characterisation table above is computed from the v9 human cache and
is not written to disk by the script; it is reproducible from `eligible_dataset()`
with the `complete_case` variant, comparing rows complete on the Gini block
against the rest.

---

*Secondary to [`XGB_STRUCTURE_REPORT.md`](XGB_STRUCTURE_REPORT.md). Generated
from `data/outputs/xgb_structure/gini/`, fitted 2026-08-17 22:12. R 4.5.2,
seed 42.*
