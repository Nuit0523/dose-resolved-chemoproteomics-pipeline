# Thesis Revision Notes: 20 uM Lower-Boundary Limitation and pEC50 Interpretation

## Purpose of This Revision

This document summarizes the final interpretation and thesis-ready wording for the 20 uM lower-boundary issue in the CurveCurator analysis.

The key issue is that the TMT experiment did not include a 0 uM or untreated channel. The lowest measured concentration was 20 uM. Therefore, the current data can support dose-response fitting within the measured 20-640 uM range, but cannot reliably identify responses whose transition and front plateau occurred entirely below 20 uM.

This is an experimental design / structural identifiability limitation, not a CurveCurator parameter problem.

## Final Conceptual Position

### What the data can say

The data can answer:

```text
Did CurveCurator produce any fitted EC50 estimates below the lowest measured concentration of 20 uM?
```

The answer is:

```text
No.
```

No fitted pEC50 values exceeded 4.699 in either 8170 or 8171.

### What the data cannot say

The data cannot answer:

```text
Were there true biological responses that reached plateau before 20 uM?
```

That is because a true early-saturating curve and a non-responsive or flat profile can look similar if all measured points fall after the transition region.

Therefore, the thesis should not claim that sub-20 uM responses were absent. It should claim only that no sub-20 uM fitted estimates were detected.

## Why pEC50 = 4.699 Matters

pEC50 is defined as:

```text
pEC50 = -log10(EC50 in M)
```

The lowest measured concentration was:

```text
20 uM = 20 x 10^-6 M
```

Therefore:

```text
pEC50 = -log10(20 x 10^-6)
      = 4.699
```

Thus:

```text
pEC50 = 4.699 corresponds to EC50 = 20 uM.
```

This value is not a biological activity cutoff. It is a lower-boundary marker for the measured concentration range.

## Boundary Classification

Measured concentration range:

```text
20, 40, 80, 160, 320, 640 uM
```

Corresponding pEC50 boundaries:

```text
20 uM  -> pEC50 = 4.699
640 uM -> pEC50 = 3.194
```

Curve classification:

| pEC50 range | EC50 interpretation | Recommended label |
|---|---|---|
| pEC50 > 4.699 | fitted EC50 < 20 uM | detected below-range estimate |
| 3.194 <= pEC50 <= 4.699 | fitted EC50 within 20-640 uM | in measured range |
| pEC50 < 3.194 | fitted EC50 > 640 uM | weak / above-range estimate |

## 8171 Boundary Results

8171 raw site-level rows:

```text
15,331
```

Successfully modeled CurveCurator profiles:

```text
13,615
```

8171 site-level fitted EC50 classification:

| Status | Sites | % of raw 15,331 rows | % of modeled 13,615 curves |
|---|---:|---:|---:|
| EC50 within 20-640 uM | 13,214 | 86.19% | 97.05% |
| EC50 above 640 uM / weak | 401 | 2.62% | 2.95% |
| No pEC50 / not modeled | 1,716 | 11.19% | not applicable |
| EC50 below 20 uM | 0 | 0% | 0% |

Recommended interpretation:

```text
Among the 13,615 successfully modeled 8171 site-level profiles, 13,214 sites (97.05%) had fitted EC50 values within the measured 20-640 uM range, whereas 401 sites (2.95%) had fitted EC50 values above 640 uM or weak/out-of-range estimates. No fitted EC50 values fell below the 20 uM lower boundary.
```

## 8171 Near-Lower-Boundary Curves

Near-lower-boundary was defined as:

```text
4.5 <= pEC50 <= 4.699
```

This corresponds approximately to:

```text
20 <= EC50 <= 31.6 uM
```

8171 near-lower-boundary curves:

| Dataset | Sites | Genes | % of raw rows | % of modeled curves |
|---|---:|---:|---:|---:|
| 8171 | 758 | 681 | 4.94% | 5.57% |

Recommended interpretation:

```text
A subset of 758 modeled sites (5.57%) had fitted pEC50 values between 4.5 and 4.699, corresponding to apparent EC50 values of approximately 20-31.6 uM. These profiles were annotated as near-lower-boundary curves but were not interpreted as sub-20 uM responses.
```

## 8171 Strict-Hit Results

Strict-hit counts should use the successfully modeled curves as the denominator:

```text
13,615 modeled curves
```

| Candidate definition | Sites | Genes | % of modeled curves | Sub-20 uM strict hits |
|---|---:|---:|---:|---:|
| Original strict definition: adjusted p <= 0.05, pEC50 >= 4, R2 >= 0.99, abs(log2FC) >= 2 | 220 | 189 | 1.62% | 0 |
| Alternative positive-response definition: adjusted p <= 0.05, pEC50 >= 4, R2 >= 0.99, log2FC >= 1 | 75 | 74 | 0.55% | 0 |

Important wording:

```text
The 75 positive-response candidates should be described as an alternative prioritization definition, not as a direct subset of the 220 original strict hits, because the effect-size criteria are different.
```

## 8171 Gene-Level Boundary Results

Gene-level results should use the gene-level best-pEC50 annotation, not the non-mutually-exclusive gene counts from site-level categories.

8171 gene-level best-pEC50 summary:

| Gene-level status | Genes |
|---|---:|
| Best fitted EC50 within 20-640 uM | 4,958 |
| Best fitted EC50 above 640 uM / weak | 84 |
| Best fitted EC50 below 20 uM | 0 |
| Total gene-level summaries | 5,042 |

Recommended interpretation:

```text
At the gene level, best-pEC50 summarization produced 5,042 genes. Among these, 4,958 genes had best fitted EC50 values within the measured range, 84 genes had best fitted EC50 values above 640 uM or weak/out-of-range estimates, and none had best fitted EC50 values below 20 uM.
```

## 8170 Boundary Results

8170 was used as a pilot/training dataset for assessing analysis behavior.

8170 site-level classification:

| Status | Sites | Genes |
|---|---:|---:|
| Fitted EC50 within 20-640 uM | 6,405 | 2,497 |
| Fitted EC50 above 640 uM / weak | 434 | 376 |
| No pEC50 / not modeled or not matched | 3,908 | 1,839 |
| Fitted EC50 below 20 uM | 0 | 0 |

8170 near-lower-boundary curves:

| Dataset | Sites | Genes |
|---|---:|---:|
| 8170 | 436 | 352 |

Recommended interpretation:

```text
The 8170 pilot dataset similarly did not contain any fitted EC50 estimates below 20 uM. A subset of profiles was close to the lower measured boundary but remained within the measured concentration range.
```

## Simulation Sensitivity Analysis

### Purpose

The simulation does not identify real missed sites. Its purpose is to evaluate whether the current 20-640 uM experimental design has limited power to detect true sub-20 uM responses.

Simulation answers:

```text
If a true response had EC50 < 20 uM, how often would the current experimental design detect it?
```

Simulation does not answer:

```text
Which real sites have EC50 < 20 uM?
```

### Simulation setup

Training dataset:

```text
8170 raw TMTMosaic data
```

Current design:

```text
20, 40, 80, 160, 320, 640 uM
```

Extended low-dose design:

```text
1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640 uM
```

Simulated true EC50 values:

```text
1, 5, 10, 15, 20, 40 uM
```

Simulated effect sizes:

```text
log2FC = 1, 2, 3
```

Signal/noise:

```text
Signal intensity, replicate CV, and slope distributions were estimated from the 8170 dataset.
```

### Pilot simulation result

This was a pilot-scale simulation and should be described as illustrative or sensitivity analysis, not as a definitive power calculation.

| True EC50 | Current design significant detection | Extended design significant detection |
|---:|---:|---:|
| 1 uM | 0% | 46.7% |
| 5 uM | 0% | 96.7% |
| 10 uM | 0% | 99.4% |
| 15 uM | 0% | 98.9% |
| 20 uM | 0% | 100% |
| 40 uM | 0% | 99.4% |

Recommended interpretation:

```text
A pilot simulation suggested that the current 20-640 uM design has limited power to detect responses whose transition occurs entirely below the lowest measured concentration. In contrast, inclusion of additional low-dose points substantially improved recovery of simulated below-range responses. This result supports the interpretation that the absence of detected sub-20 uM EC50 estimates should not be treated as evidence that such responses are biologically absent.
```

## Recommended Methods Revision

Use this in the CurveCurator / dose-response fitting Methods section:

```text
Curve fitting was performed using the experimentally measured TMT concentration series, which included six non-zero probe concentrations ranging from 20 to 640 uM. No untreated or 0 uM TMT channel was available, and no artificial zero-dose or simulated low-dose values were introduced during CurveCurator input generation or fitting. The lowest measured concentration of 20 uM corresponds to a pEC50 value of 4.699. Therefore, fitted pEC50 values above 4.699 would indicate model-estimated EC50 values below the measured concentration range. Such estimates were treated as below-range model outputs rather than directly observed potency measurements.
```

Optional additional Methods text:

```text
To assess the effect of the measured concentration range on detectability of potential below-range responses, a pilot simulation-based sensitivity analysis was performed using empirical signal intensity, replicate variability, and curve-shape distributions estimated from the 8170 pilot dataset. Synthetic dose-response profiles with known EC50 values were sampled under the current 20-640 uM design and under a hypothetical extended low-dose design. Simulated profiles were analyzed using the same general curve-fitting and prioritization logic to estimate the probability of recovering below-range responses under each design.
```

## Recommended Results Revision

Use this in the CurveCurator Results section:

```text
The 8171 dataset contained 15,331 site-level rows, of which 13,615 profiles were successfully modeled by CurveCurator. Among the modeled profiles, 13,214 sites (97.05%) had fitted EC50 values within the measured 20-640 uM concentration range, while 401 sites (2.95%) had fitted EC50 values above 640 uM or weak/out-of-range estimates. No fitted pEC50 values exceeded 4.699, the value corresponding to an EC50 of 20 uM. Therefore, no below-range EC50 estimates were detected from the measured 20-640 uM profiles.
```

Add this near strict hits:

```text
Using the original strict-hit definition, 220 sites corresponding to 189 genes passed the combined significance, potency, curve-quality, and effect-size criteria, representing 1.62% of the 13,615 modeled profiles. None of these strict-hit sites had fitted EC50 values below 20 uM.
```

Optional near-boundary sentence:

```text
A subset of 758 modeled sites (5.57%) had fitted pEC50 values between 4.5 and 4.699, corresponding to apparent EC50 values of approximately 20-31.6 uM. These sites were considered near the lower measured boundary but were not interpreted as sub-20 uM responses.
```

Optional simulation sentence:

```text
Pilot simulation analysis further suggested that responses with true EC50 values below 20 uM may have limited detectability under the current 20-640 uM design, particularly when the transition region occurs entirely below the lowest measured concentration.
```

## Recommended Discussion / Limitations Revision

Use this in the Discussion or Limitations section:

```text
A limitation of the present TMT design is that 20 uM was the lowest experimentally measured concentration and no untreated control channel was included. As a result, responses whose transition and front plateau occurred entirely below 20 uM were not identifiable from the available data. Although no fitted EC50 estimates below 20 uM were detected in either dataset, this should be interpreted as the absence of detected below-range estimates rather than evidence that below-range biological responses were absent. Additional measurements below 20 uM, ideally including an untreated control, would be required to distinguish early-saturating responses from flat or non-responsive profiles within the current measured range.
```

Optional simulation limitation sentence:

```text
Consistent with this limitation, pilot simulation analysis indicated that the current concentration design may have low power to recover simulated responses with true EC50 values below 20 uM, whereas inclusion of additional low-dose measurements would substantially improve detectability.
```

## What Not to Claim

Do not write:

```text
No sub-20 uM responses existed.
```

Do not write:

```text
The 758 near-boundary curves are sub-20 uM hits.
```

Do not write:

```text
Simulation recovered the missed real sites.
```

Do not write:

```text
Artificial 0 uM values were added to improve fitting.
```

## What to Claim Instead

Write:

```text
No fitted EC50 estimates below 20 uM were detected.
```

Write:

```text
Near-lower-boundary curves were annotated but not interpreted as sub-20 uM responses.
```

Write:

```text
The absence of detected below-range EC50 estimates does not exclude biological responses that occurred below the measured concentration range.
```

Write:

```text
The current dataset is best suited for prioritizing responses with measurable transitions within the 20-640 uM range.
```

## Final Thesis Position

The final thesis should frame this issue as follows:

```text
The CurveCurator workflow was retained as the primary analysis because it used only measured experimental data. Boundary analysis showed that no fitted EC50 estimates fell below the 20 uM lower measured boundary. However, because the experiment did not include concentrations below 20 uM or an untreated control, the analysis could not rule out early-saturating responses whose transition occurred entirely below the measured range. This limitation was therefore treated as an experimental design constraint and a direction for future optimization rather than corrected by artificial input values or reinterpretation of flat curves.
```

## Integration with Advisor Revision: Overall Thesis Positioning

The 20 uM boundary discussion should be integrated into the broader revised thesis position requested by the advisor:

```text
This thesis should be framed as a large-scale exploratory workflow for prioritization and annotation, not as definitive target discovery or mechanistic validation.
```

This means:

- pEC50 should be described as a curve-derived prioritization metric, not as a fully validated potency measurement for every site.
- ORA/GSEA should be described as exploratory pathway annotation, not proof of pathway regulation.
- AlphaFold3 and docking should be described as structure-guided hypothesis generation, not evidence of direct binding or mechanism.
- Candidate proteins should be described as prioritized candidates for future validation, not confirmed targets.

Recommended combined limitation sentence:

```text
Because the experiment lacked an untreated channel and concentrations below 20 uM, and because docking was used only as a computational prediction, the resulting candidate sites should be interpreted as prioritized hypotheses rather than validated targets or mechanisms.
```

Recommended final thesis-level statement:

```text
This project developed and applied a reproducible exploratory workflow for large-scale prioritization and annotation of candidate probe-responsive cysteine sites. The workflow integrates dose-response modeling, enrichment analysis, and structural prediction to generate biologically interpretable hypotheses, but experimental validation is required before drawing firm conclusions about direct binding, target identity, or mechanism.
```
