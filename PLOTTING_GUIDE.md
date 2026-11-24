# Plotting style guide

How figures in this project should be drawn and saved, followed by an audit of
every plot that existed on 2026-11-25 against these rules. The guide is the
target; the audit is the to-do list. Plots are brought into line one at a time,
by choice, not by a project-wide rewrite (see "Why not one global change").

---

## 1. Axis assignment

### 1.1 Observed vs predicted: observed on y, predicted on x

Any plot comparing a model's output with the thing it models puts the
**prediction on x** and the **observation on y** (Piñeiro et al. 2008,
*Ecological Modelling* 216:316–322). The observed values are the response the
model is being judged against; regressing observed on predicted (O ~ P) gives
an unbiased slope and intercept for judging agreement, while P ~ O does not.

- Always draw the 1:1 line (`y = x`), with equal limits on both axes and a
  square panel (`coord_equal()` or `theme(aspect.ratio = 1)` with matching
  `xlim`/`ylim`), so the diagonal really sits at 45°.
- Reading rule, worth stating in the legend: points **above** the line are
  under-predicted, points **below** are over-predicted.
- If a fitted line is drawn, it is `lm(observed ~ predicted)`, and its slope and
  intercept are what get reported (ideal: 1 and 0).
- Goodness of fit is R² = 1 − SS_res/SS_tot on held-out data (as
  `xgb_structure_plots.R` already does), not the squared correlation; RMSE
  alongside it. A rank correlation is fine for QC plots.
- Title wording follows the axes: "Observed vs predicted <thing>" (y vs x).

What counts as "predicted": any model, fit or estimate (XGBoost, the
GC/length expected-MFE model, the sum-of-components imputation, Saluki).
What counts as "observed": the measured or directly computed value the
prediction is standing in for (measured half-life, RNAfold MFE on the full
sequence). A null/reference distribution such as dinucleotide-shuffled MFE is
the "observed" side when it is what the prediction is checked against.

Note: Spearman ρ and R² are symmetric or computed independently of the plot,
so flipping axes on the existing plots changes no reported number. The rule
matters for reading the plot and for any slope/intercept that gets added.

### 1.2 Feature vs response: feature on x, response on y

Where one variable is the outcome (half-life, translation efficiency), it goes
on y and the feature on x. All current feature/response plots already do this.

### 1.3 Symmetric comparisons

Two measurements with no direction (human vs mouse, probing vs probing) can go
either way; pick one and keep it across related figures. Default: the
reference species/condition on x.

### 1.4 Residuals

Residual (observed − predicted, or a relative version) on y, against the
predicted value or the covariate of interest on x, with a horizontal line at 0.

---

## 2. Output targets

Decide the target **first**: font sizes only mean something at the size the
figure will be printed or projected. Draw at final size; never shrink a 300 mm
figure into a column.

| Target | Width | `base_size` | Notes |
|---|---|---|---|
| Paper, single column | 85 mm | 8 | most journals 85–89 mm |
| Paper, 1.5 column | 120 mm | 8 | |
| Paper, full width | 180 mm | 8 | max height ~ 230 mm |
| Slide (16:9) | 254 × 143 mm | 14 | the xgb PDFs are slide assets |
| Exploratory / QC | any | 11 | goes in run folders, not the paper |

Printed text should be no smaller than 6 pt after any scaling.

### 2.1 Saving

- Units: `mm`, always (no bare inches).
- Vector: PDF via `grDevices::cairo_pdf` (the default `pdf()` device drops
  UTF-8 glyphs such as ρ and ²; see `xgb_structure_plots.R`).
- Raster: PNG at 300 dpi. JPEG blurs text and thin lines; keep it only where a
  downstream tool needs it.
- Paper and slide figures save both PDF and PNG; QC figures may save PNG only.

---

## 3. Theme

- Base: `theme_bw()` for faceted or gridded plots, `theme_minimal()` only for
  heatmaps and bar/dot charts where the panel border adds nothing. Within one
  figure, one theme.
- Title: bold, left-aligned (`hjust = 0`), size `base_size + 2`. No centred or
  right-aligned titles.
- Subtitle: sample sizes and definitions only; interpretation belongs in the
  figure legend (the policy already written in `xgb_structure_plots.R`).
- Facet strips: bold, `base_size`.
- Minor gridlines off.
- Axis labels in sentence case with units in parentheses:
  "Half-life (h)", "MFE (kcal/mol)". Use `format_col_name()` /
  `format_metric_name()` rather than hand-typed labels where possible.

---

## 4. Colour and marks

Project palettes live in `R/utils/palettes.R` and `R/colour_config.R`; use them
for regions, feature groups and supergroups rather than local hex codes.

| Role | Setting |
|---|---|
| 1:1 / identity line | `"grey40"`, dashed, `linewidth = 0.4` |
| Zero / reference line | `"grey20"`, solid, `linewidth = 0.3` |
| Linear trend | `"firebrick"`, `linewidth = 0.6`, `se = FALSE` |
| LOESS trend | `"steelblue"`, `linewidth = 0.6`, `span = 0.75` |
| Density (hex) fill | `scale_fill_viridis_c(trans = "log10")`, legend "Genes" or "Transcripts" |
| Correlation / signed fill | `scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0)`, symmetric limits |
| Proportion / coverage fill | sequential (viridis or a single-hue ramp), never diverging |
| Two-species contrast | human `"#1F78B4"`, mouse `"#E66101"` |

- Over ~2,000 points: hex-bin (`geom_hex`) instead of a scatter.
- Scatters: `alpha` 0.2–0.4, `size` 0.5–1.2, `shape = 16`.
- Colour must not be the only carrier of a distinction that matters; pair it
  with shape or facet when it does (regions already use `REGION_SHAPES`).

---

## 5. Statistics on the panel

- Correlation: `"Spearman ρ = 0.42"` (Unicode ρ, two decimals; three only in QC
  plots where small differences matter), then `n = 1,234` with thousands
  separators. Never the ASCII word "rho" on a figure.
- q-values: `"q = 3.1e-05"`, BH-adjusted, method named once in the subtitle
  or legend.
- CIs: bracketed interval ("95% CI 0.41–0.44"), not ±, since percentile
  intervals are asymmetric.
- Position: top-left inside the panel (`x = -Inf, y = Inf, hjust = -0.08,
  vjust = 1.6`), `size = base_size * 0.3` (i.e. in mm for `geom_text`).

---

## 6. Why not one global change

A single project theme applied everywhere would not be safe:

- Font sizes, widths and heights are tuned together per plot, and range from
  base 10 at 110 mm to base 18 at 300 mm. Changing one without the other breaks
  layouts (labels overlap, facets collide), so each plot has to be looked at.
- Some figures (Figure 5, Figure 6, the xgb set) were designed deliberately and
  recently; others are exploratory sweeps whose look matters little.
- Axis flips change titles, limits and annotation positions, which have to be
  checked by eye.

What *is* safe to do project-wide, because it changes no pixel, is to replace
the six near-identical save helpers (`save_plot_jpg_pdf` ×3, `f5_save`,
`f6_save`, `save_plot`) with one shared helper, and to add shared constants for
the line, fill and label settings above so plots adopt them by name as they are
revised.

---

## 7. Audit (2026-11-25)

Deviations from the guide, per script. **Bold** = axis-assignment issues
(section 1), the priority.

### 7.1 Observed vs predicted plots

Flipped on 2026-11-25 (prediction on x, observation on y); titles and labels
follow. What remains on each is styling only, listed in 7.2.

| Plot | Status |
|---|---|
| `analysis/models/xgb_structure_plots.R` observed vs predicted half-life | Flipped. Already `coord_equal`, dashed 1:1, R² with CI |
| `analysis/qc/mfe_observed_vs_predicted.R` observed vs predicted MFE (fig. 2) | Flipped. Report's "slope (y on x)" is now observed on predicted, and its "above y = x" share is now the complement (e.g. 3' UTR slope 0.83 -> 1.18, mRNA 0.89 -> 1.10) |
| `analysis/qc/mfe_observed_vs_predicted.R` shuffle median vs predicted MFE (fig. 3) | Flipped, same report changes; lower-case axis labels fixed |
| `analysis/qc/mfe_expected_check.R` observed vs predicted <region> MFE | Flipped; now matches the `observed ~ expected` calibration slope it already reported. Still: red dashed 1:1 at `linewidth = 1`, centred 20 pt title, x and y zoomed separately so 1:1 is not at 45°, x label is the raw "MFE.pred" name |
| `analysis/qc/imputation_check.R` observed vs estimated mRNA MFE (+ z-scores) | Flipped |
| `analysis/structure/figure5.R` Gini vs MFE z-score | Not changed. Judgement call: two structure measures, not a model fit. If MFE z is framed as predicting probing structure, swap; otherwise leave (section 1.3) |

### 7.2 Theme, size and export

| Script | Theme / base | Title | Export | Other deviations |
|---|---|---|---|---|
| `correlations/correlation_heatmap_workflow.R` | bw / 14 | +2, bold | PDF (cairo), mm | ok |
| `correlations/region_feature_heatmap.R` | bw / 14 | +2, bold | PDF (cairo), mm | ok |
| `correlations/feature_correlation_bands.R` | minimal / 18 | 22 | JPG 300 + PDF | large base for its width |
| `correlations/feature_correlation_dotplot.R` | bw / 18 | 22 | JPG 300 + PDF | large base |
| `correlations/feature_correlation_ranked.R` | bw / 18 | 22 | JPG 300 + PDF | large base |
| `correlations/feature_response_scatter.R` | bw / 11 | 14 | JPG 300 + PDF | ok |
| `correlations/feature_response_scatter_top30.R` | bw / 11 | 14 | JPG + PDF, 300×220 mm | wider than any print target |
| `correlations/feature_response_hex_panels.R` | bw / 11 | 14 | JPG only | density legend "count" |
| `correlations/group_panel_sweep.R` | bw / 11 | 14 | JPG only | density legend "count" |
| `correlations/scatter_plot.R` | bw / 11, text 16–19 | 19, centred | (caller) | lm `#4B0082`, loess `#FF4500` (off-palette) |
| `correlations/halflife_correlation.R` | minimal / 11 | 20, right-aligned | JPG only | |
| `correlations/top_correlations.R` | bw / 11 | 16 | JPG only, 210×240 | |
| `qc/dataset_overview.R` | bw / 11 | 16 | JPG only | coverage % on a diverging red→blue scale (reversed relative to every other diverging fill, and the quantity is not signed) |
| `qc/mfe_expected_check.R` | minimal / 11, text 14–20 | 20, centred | JPG only, 250×180 | `grey96` panel; three different viridis options |
| `qc/mfe_observed_vs_predicted.R` | `theme_mfe` (minimal / 13) | 17 | JPG only | blue gradient fill differs from viridis used elsewhere; `theme_mfe` duplicated in `imputation_check.R` |
| `qc/imputation_check.R` | `theme_mfe` (minimal / 13) | 17 | JPG only | 1:1 line `#D55E00` solid |
| `models/xgb_structure_plots.R` | `theme_xgb` (minimal / 12) | 14 | JPG + PDF, `save_plot()` | 1:1 line grey50 dashed; good reference implementation |
| `structure/figure5.R` | minimal / 15–16 | = base | JPG + PDF | ok |
| `structure/figure6.R` | bw / 14 | = base | JPG + PDF | `grey93` panel |
| `structure/figure6_b.R` | bw / 11 | +1 | PDF + PNG + JPG | ok |
| `cross_species/cross_species_probing_concordance.R` | bw / 11 | default | PDF only, **inches** (10×7, 7×6, 14×6) | ASCII "rho" in labels; out of scope (mouse) per current project scope |
| `rnaplfold/rnaplfold_metric_screen.R` | bw / 10 and 12 | default | JPG 200 dpi | two bases in one script |
| `rnaplfold/rnaplfold_hex_panels.R` | (from hex panels) | | JPG 200 dpi | |
| `rnaup_avoidance/rnaup_avoidance_screen.R` | bw / 10 | default | JPG 200 dpi | |
| `scripts/example_analysis.R` | minimal / 11 | default | JPG 300 | example only |

### 7.3 Cross-cutting inconsistencies

- **1:1 lines**: four styles (red dashed 1.0; `#D55E00` solid 0.5; grey50
  dashed; grey55 dashed 0.4).
- **Density fills**: viridis log10 ("count", "Genes") vs blue gradient log10
  ("Transcript density").
- **Correlation labels**: "Spearman ρ = %.3f", "ρ = %.2f, n = …",
  "rho = %.2f, p = …".
- **Title alignment**: left (most), centred (`scatter_plot.R`,
  `mfe_expected_check.R`, `imputation_check.R` delta plot), right
  (`halflife_correlation.R`).
- **Save helpers**: six local copies of the same JPG+PDF logic.
- **Supergroup colours** (`R/colour_config.R`) are pure RGB primaries
  (`#FF0000`, `#00FF00`, `#FFFF00`, `#00FFFF`); yellow and cyan have poor
  contrast on white and red/green is the commonest colour-vision confusion.
  Features derived from them are spread in lightness, which helps, but the
  bands and labels use the raw hues.
- **Stray `Rplots.pdf`**: `imputation_check.R` calls `print()` on its plots,
  which under `Rscript` writes the tracked `Rplots.pdf` in the repo root through
  the default `pdf()` device (and drops ρ with an `mbcsToSbcs` warning).
  `mfe_expected_check.R` already saves instead of printing for this reason.
