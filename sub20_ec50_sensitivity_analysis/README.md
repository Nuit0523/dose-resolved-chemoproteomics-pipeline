# Sub-20 uM EC50 Boundary and Detectability Sensitivity Analysis

## Important interpretation

This folder documents an exploratory sensitivity analysis performed during thesis revision. It should not be interpreted as a definitive biological conclusion or as evidence that true sub-20 uM probe-responsive sites are absent.

The TMT experiments used six non-zero measured probe concentrations: 20, 40, 80, 160, 320, and 640 uM. Because no untreated/0 uM channel and no concentrations below 20 uM were measured, responses whose transition occurred entirely below 20 uM are structurally difficult, and in some cases impossible, to identify from these data alone.

Therefore, this analysis asks a limited technical question:

> Did the existing 20-640 uM data produce fitted EC50 estimates below the measured lower boundary, and how sensitive is the current design to hypothetical sub-20 uM responses?

It does not assign real sites to confirmed sub-20 uM EC50 classes.

## Boundary definition

The lower concentration boundary is defined by the lowest measured concentration:

```text
lowest measured concentration = 20 uM = 20 x 10^-6 M
pEC50 boundary = -log10(20 x 10^-6) = 4.699
```

The upper measured concentration boundary is:

```text
highest measured concentration = 640 uM
pEC50 upper-range boundary = -log10(640 x 10^-6) = 3.194
```

Curve-level pEC50 values were annotated as:

| Status | Rule | Interpretation |
| --- | --- | --- |
| fitted_EC50_below_20uM | pEC50 > 4.699 | Fitted EC50 extrapolates below the lowest measured concentration |
| fitted_EC50_in_20_640uM | 3.194 <= pEC50 <= 4.699 | Fitted EC50 lies within the measured range |
| fitted_EC50_above_640uM_or_weak | pEC50 < 3.194 | Fitted EC50 is above the measured range or weak |
| no_pEC50_or_not_modeled | missing pEC50 | No interpretable fitted EC50 value |

Near-lower-boundary curves were also annotated using:

```text
4.5 <= pEC50 <= 4.699
```

This corresponds approximately to EC50 values between 20 and 31.6 uM. These curves are not treated as confirmed sub-20 uM responses.

## Pipeline

### 1. Boundary annotation of CurveCurator outputs

Script:

```text
scripts/annotate_pec50_boundary_status.py
```

Purpose:

- Read CurveCurator-derived pEC50 outputs for 8170 and 8171.
- Annotate fitted EC50 status relative to the measured 20-640 uM range.
- Summarize site-level and gene-level counts.

Key output:

```text
results/pec50_boundary_annotation_summary.csv
results/pec50_boundary_annotation_summary.md
```

### 2. Raw-row and strict-hit overlap analysis

Script:

```text
scripts/raw_boundary_hit_analysis_8170_8171.py
```

Purpose:

- Merge raw TMTMosaic site rows with CurveCurator fitted metrics.
- Count detected sub-20 uM fitted EC50 values.
- Check whether sub-20 uM fitted values overlap with strict-hit definitions.
- Annotate near-lower-boundary curves.

Strict-hit definitions used for sensitivity checks:

```text
Original strict-hit definition:
adjusted p <= 0.05
pEC50 >= 4
R2 >= 0.99
abs(log2 curve fold change) >= 2
```

```text
Alternative positive-response definition:
adjusted p <= 0.05
pEC50 >= 4
R2 >= 0.99
log2 curve fold change >= 1
```

Important note: these are alternative prioritization definitions, not a new validated hit-calling strategy.

Key output:

```text
results/raw_boundary_hit_summary.csv
results/raw_boundary_hit_summary.md
```

### 3. 8170-trained exploratory detectability simulation

Script:

```text
scripts/sub20_detectability_simulation_8170_trained.py
```

Purpose:

- Use 8170 as a pilot/training dataset to estimate realistic noise, baseline intensity, and fitted slope distributions.
- Simulate hypothetical dose-response curves under the current 20-640 uM design.
- Compare with an extended low-dose design that includes concentrations below 20 uM.

Simulation settings:

| Setting | Values |
| --- | --- |
| Current design | 20, 40, 80, 160, 320, 640 uM |
| Extended low-dose design | 1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640 uM |
| True EC50 values | 1, 5, 10, 15, 20, 40 uM |
| True log2 fold changes | 1, 2, 3 |
| Signal groups | low, medium, high |
| Pilot replicates per parameter cell | 20 |

Detection definitions:

```text
recovered_sub20_fit:
fitted pEC50 > 4.699
```

```text
detected_significant:
BH-adjusted p <= 0.05 and R2 >= 0.8
```

```text
strict_detected:
BH-adjusted p <= 0.05 and R2 >= 0.99 and observed-range log2FC >= 1
```

The simulation was intentionally treated as exploratory. The pilot replicate count is small and was used to illustrate the experimental-design blind spot, not to produce a formal detection-power estimate.

Key outputs:

```text
results/simulation_detection_summary_by_true_EC50.csv
results/simulation_detection_summary_by_cell.csv
results/simulation_detection_summary.md
```

## Key results

### Existing fitted data

- 8170: no fitted pEC50 values exceeded 4.699; no detected sub-20 uM fitted EC50 sites.
- 8171: no fitted pEC50 values exceeded 4.699; no detected sub-20 uM fitted EC50 sites.
- 8171 original strict-hit branch: 220 sites / 189 genes; sub-20 uM overlap = 0.
- 8171 alternative positive-response branch: 75 sites / 74 genes; sub-20 uM overlap = 0.
- 8171 near-lower-boundary annotation: 758 sites / 681 genes with 4.5 <= pEC50 <= 4.699.

The correct interpretation is:

> No sub-20 uM EC50 values were recovered by CurveCurator from residual dose-dependent variation in the measured 20-640 uM range.

This is not equivalent to:

> No true sub-20 uM biological responses exist.

### Simulation

The exploratory 8170-trained simulation suggested that a design beginning at 20 uM has limited ability to recover strong responses whose transition occurs entirely before the first measured concentration. Adding measured concentrations below 20 uM substantially improves detectability in simulated scenarios.

This simulation is a design-sensitivity analysis only. It should not be used to reclassify real sites as sub-20 uM hits.

## Thesis-revision positioning

Recommended wording:

> Because 20 uM was the lowest experimentally measured concentration, responses whose transition and front plateau occurred entirely below 20 uM were not identifiable from the available data. No artificial zero-dose or simulated low-dose observations were incorporated into experimental curve fitting.

> No fitted pEC50 values exceeded 4.699, the value corresponding to an EC50 of 20 uM. This result indicates that no below-range responses were identified by CurveCurator from the residual concentration-dependent variation observed between 20 and 640 uM.

> The absence of fitted sub-20 uM EC50 values should not be interpreted as evidence that such responses were biologically absent. Additional measurements below 20 uM, ideally including an untreated control, would be required to resolve these possibilities.

## Files intentionally not committed

Full per-site and per-gene derived annotation tables were not committed to keep the repository lightweight and avoid distributing large row-level data files. They can be regenerated from the original local raw data and CurveCurator outputs using the scripts in this folder.

## Relation to the main thesis

This folder supports a cautious revision of the thesis. The overall thesis should be framed as a large-scale exploratory workflow for prioritization and annotation of probe-responsive sites, not as definitive proof of target engagement or mechanism. Docking and AlphaFold-based structural inspection should likewise be described as hypothesis-generating computational prioritization rather than experimental validation.
