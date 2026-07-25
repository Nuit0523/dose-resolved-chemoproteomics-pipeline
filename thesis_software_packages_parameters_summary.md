# Thesis Software, Packages, Settings, and Analysis Parameters

This document summarizes the software, packages, major settings, and parameter cutoffs used in the thesis analysis. Items marked as "not recorded" were not found in the local project files and should be filled in manually before final submission if needed.

## 1. Raw Data Source and TMT/LC-MS Processing

### Input data

| Item | Value |
|---|---|
| Raw quantitative file | `E:/R/TMT_analysis/TMT4/Raw_sq_8171_TMTMosaic.csv` |
| Raw rows | 15,331 |
| Quantitative columns | TMT reporter-ion `sn_sum` columns |
| TMT channels | 18 channels |
| Experimental design | 6 concentrations x 3 replicates |
| Concentrations | 20, 40, 80, 160, 320, 640 uM |
| Probe/compound | WRX-035 / wxr35 |

### Preprocessing confirmed from project files

| Step | Setting / output |
|---|---|
| Reshaping | Raw TMTMosaic table was reshaped into long format |
| Long table | 245,214 site-by-channel rows |
| Replicate averaging | Three TMT channels per concentration were averaged |
| CurveCurator input | 13,615 site-level curves in `8171_conc.txt` |
| Input scale to CurveCurator | Mean linear TMT reporter S/N values, not log2-transformed |
| Baseline comparison in summary tables | `log2FC_vs_baseline` was calculated against the lowest dose/baseline group |

### Important interpretation

The CurveCurator input values are averaged raw reporter S/N values per dose. CurveCurator then performs its own internal normalization when `normalization = true`. Therefore, CurveCurator curve fold change reflects the normalized curve behavior, not simply the raw increase/decrease in TMT intensity.

## 2. CurveCurator

### Software version

| Software / package | Version |
|---|---|
| Python in CurveCurator environment | 3.12.13 |
| CurveCurator pipeline | 0.6.0 |
| `curve-curator` package | 0.6.0 |
| `curve_curator` package | 0.6.0 |
| numpy | 2.4.3 |
| pandas | 2.3.3 |
| scipy | 1.17.1 |
| statsmodels | 0.14.6 |
| bokeh | 3.7.3 |

### 8171 CurveCurator configuration

Source file: `E:/gradthesis/curvecurator/8171/8171_conc_config.toml`

| Category | Parameter | Value |
|---|---|---|
| Meta | id | `TMT4_direct_labeling_wxr35` |
| Meta | condition | `wxr35_treatment` |
| Meta | treatment_time | `2 h` |
| Experiment | experiments | `1, 2, 3, 4, 5, 6` |
| Experiment | doses | `20, 40, 80, 160, 320, 640` |
| Experiment | dose_scale | `1e-6` |
| Experiment | dose_unit | `M` |
| Experiment | control_experiment | `1` |
| Experiment | measurement_type | `TMT` |
| Experiment | data_type | `PROTEIN` |
| Processing | available_cores | `4` |
| Processing | imputation | `true` |
| Processing | max_missing | `2` |
| Processing | max_imputation | `2` |
| Processing | normalization | `true` |
| Curve fit | type | `OLS` |
| Curve fit | speed | `standard` |
| Curve fit | max_iterations | `1000` |
| Curve fit | control_fold_change | `true` |
| Curve fit | interpolation | `false` |
| F-statistic | alpha | `0.05` |
| F-statistic | fc_lim | `1.0` |
| F-statistic | optimized_dofs | `true` |
| F-statistic | mtc_method | `fdr_bh` |
| F-statistic | not_rmse_limit | `0.1` |
| Dashboard | backend | `webgl` |

Note: `imputation_pct = 0.005` appears in the TOML file, but the CurveCurator log reports this key as unknown. The effective imputation value was therefore determined by CurveCurator during processing.

### 8171 CurveCurator log-confirmed settings and results

Source log: `E:/gradthesis/curvecurator/doc/curveCurator.log`

| Item | Value |
|---|---|
| Curves loaded | 13,615 |
| Curves removed for >2 missing values | 0 |
| Imputation value | 2.95 |
| Curves removed for >2 imputed values | 0 |
| Normalization factors | Raw1 2.07; Raw2 1.08; Raw3 0.08; Raw4 -0.80; Raw5 -1.17; Raw6 -1.27 |
| Statistical fitting | Standard OLS |
| CPU cores | 4 |
| Multiple-testing correction | Benjamini-Hochberg FDR (`fdr_bh`) |
| Decoys simulated | 13,615 |

### 8170 CurveCurator configuration

Source file: `E:/gradthesis/curvecurator/8170/8170_conc_config.toml`

| Category | Parameter | Value |
|---|---|---|
| Experiments | experiments | `1, 2, 3, 4, 5, 6, 7` |
| Doses | doses | `20, 40, 80, 160, 320, 640, 0` |
| Control | control_experiment | `7` |
| Input | file | `8170_curvecurator_input_rowmedian_control.txt` |
| Output curves | file | `8170_output_curves_row.tsv` |
| Processing | imputation | `true` |
| Processing | max_missing | `2` |
| Processing | max_imputation | `2` |
| Processing | normalization | `true` |
| Curve fit | type | `OLS` |
| Multiple testing | mtc_method | `fdr_bh` |
| RMSE threshold | not_rmse_limit | `0.1` |

8170 log-confirmed values: imputation value 8.35; 189 curves removed for >2 imputed values; 5,695 decoys simulated.

## 3. Curve-Level QC Metrics and Hit Definitions

### Curve quality metrics

| Metric | Meaning | Thesis interpretation |
|---|---|---|
| R2 | Goodness of fit | Whether the fitted 4PL curve captures the overall dose-response trend |
| RMSE | Root mean square error | Point-wise deviation between observed and fitted values |
| Adjusted p-value | FDR-adjusted curve significance | Whether the dose-response curve is statistically significant after multiple testing |
| pEC50 | Potency metric | `pEC50 = -log10(EC50)`; higher pEC50 means lower EC50 and stronger apparent potency |
| Curve fold change | Response magnitude | Normalized curve-level response magnitude from CurveCurator |

### Original strict hit definition used in thesis slides

| Criterion | Threshold | Meaning |
|---|---|---|
| Adjusted p-value | `<= 0.05` | Statistically significant curve after BH-FDR correction |
| Adjusted log-p | `>= 1.30103` | Equivalent to adjusted p <= 0.05 |
| pEC50 | `>= 4` | EC50 <= 100 uM |
| R2 | `>= 0.99` | High-quality curve fit |
| Absolute log2 curve fold change | `>= 2` | Large response magnitude |

Original result: 220 hit sites corresponding to 189 genes.

### Revised positive-response strict hit sensitivity check

After inspecting the raw intensity trend, a stricter directional check was also tested:

| Criterion | Threshold |
|---|---|
| Adjusted p-value | `<= 0.05` |
| pEC50 | `>= 4` |
| R2 | `>= 0.99` |
| Curve fold change | `>= 1` |

Result: 75 hit sites corresponding to 74 genes.

## 4. Gene-Level Summarization

### Why gene-level summarization was required

CurveCurator fits curves at the site level, but ORA and GSEA require gene-level input. Therefore, cysteine-site measurements were collapsed into gene-level summaries.

### Gene summary metrics

| Gene-level metric | Meaning |
|---|---|
| `best_pEC50` | Highest pEC50 among sites assigned to the same gene |
| `mean_pEC50` | Average pEC50 across quantified sites |
| `median_pEC50` | Median pEC50 across quantified sites |
| `weighted_score` | Composite ranking score used for GSEA robustness analysis |
| `n_sites` | Number of quantified responsive sites assigned to the gene |

### Main gene-level counts

| Step | Count |
|---|---:|
| Gene-level pEC50 summaries | 5,042 genes |
| Genes mapped for GSEA | 4,994 genes |
| Unmapped genes for GSEA | 48 genes |
| Candidate genes with `n_sites >= 2` | 2,967 genes |
| Candidate genes with `best_pEC50 >= 4` | 3,075 genes |
| Intersection candidate set | 2,256 genes |

## 5. ORA: Over-Representation Analysis

### R software and packages

| Software / package | Version confirmed in current analysis environment |
|---|---|
| R | 4.5.2 |
| clusterProfiler | 4.18.4 |
| org.Hs.eg.db | 3.22.0 |
| ReactomePA | 1.54.0 |
| enrichplot | 1.30.5 |
| AnnotationDbi | 1.72.0 |
| dplyr | 1.2.1 |
| stringr | 1.6.0 |
| ggplot2 | 4.0.3 |
| BiocManager | 1.30.27 |
| readr | Used in original script; version not recorded in current project environment |

### Original 8171 ORA parameters

Source: `E:/gradthesis/curvecurator/8171/enrichment_from_raw_name_pEC50/pipeline_parameters.csv`

| Parameter | Value |
|---|---|
| Input file | `E:/python/curvecurator/8171/8171_EC50.csv` |
| Candidate criterion 1 | `n_sites >= 2` |
| Candidate criterion 2 | `best_pEC50 >= 4` |
| Candidate genes | 2,256 |
| TopN lists | Top100, Top200, Top400 |
| Ranking metric | `best_pEC50` descending |
| GO p-value cutoff | 0.05 |
| GO q-value/FDR cutoff | 0.2 |
| KEGG p-value cutoff | 0.05 |
| Reactome p-value cutoff | 0.05 |
| Multiple-testing correction | BH |
| Plot display | Top 15 categories |

### Original 8171 ORA results

| TopN | Input genes | Mapped genes | GO terms | KEGG terms | Reactome terms |
|---:|---:|---:|---:|---:|---:|
| 100 | 100 | 100 | 0 | 0 | 0 |
| 200 | 200 | 200 | 0 | 0 | 0 |
| 400 | 400 | 399 | 0 | 0 | 0 |

Mapping note: the full 2,256-gene candidate set contains 12 unmapped genes. In the Top400 ORA set, one gene was unmapped: `TTC37`.

### 8170 ORA parameters

Source: `E:/gradthesis/curvecurator/doc/enrichment.R`

| Parameter | Value |
|---|---|
| Input file | `E:/python/curvecurator/8170/gene_level_pEC50_summary.csv` |
| Candidate criterion 1 | `n_sites >= 2` |
| Candidate criterion 2 | `best_pEC50 >= 4` |
| TopN lists | Top30, Top50, Top80 |
| Ranking | `best_pEC50` descending, then `n_sites` descending |
| GO method | `clusterProfiler::enrichGO` |
| GO ontology | Biological Process (`ont = "BP"`) |
| KEGG method | `clusterProfiler::enrichKEGG` |
| Reactome method | `ReactomePA::enrichPathway` |
| p-value cutoff | 0.2 in original pilot script |
| GO q-value cutoff | 0.5 in original pilot script |
| Multiple-testing correction | BH |

## 6. GSEA: Gene Set Enrichment Analysis

### R packages

GSEA was performed using the Bioconductor enrichment workflow centered on `clusterProfiler`/`fgsea`.

| Package | Version confirmed in current analysis environment |
|---|---|
| fgsea | 1.36.2 |
| clusterProfiler | 4.18.4 |
| DOSE | 4.4.0 |
| org.Hs.eg.db | 3.22.0 |
| ReactomePA | 1.54.0 |
| enrichplot | 1.30.5 |

### 8171 GSEA parameters

Source: `E:/gradthesis/curvecurator/8171/GSEA_from_raw_name_pEC50/GSEA_parameters.csv`

| Parameter | Value |
|---|---|
| Input file | `E:/python/curvecurator/8171/8171_EC50.csv` |
| Gene-level input | Full ranked gene list |
| Ranked genes after mapping | 4,994 |
| `min_valid_pEC50` | 0 |
| `minGSSize` | 10 |
| `maxGSSize` | 500 |
| `pvalueCutoff` | 0.2 |
| `nPermSimple` | 10,000 |
| Plot display | Top 20 categories |
| Ranking strategies compared | `best_pEC50`, `mean_pEC50`, `weighted_score` |

### Key GSEA result summary

| Ranking strategy | Ranked genes | Main enriched signal |
|---|---:|---|
| `best_pEC50` | 4,994 | Intracellular protein transport; nucleocytoplasmic transport |
| `mean_pEC50` | 4,994 | Weaker enrichment signal |
| `weighted_score` | 4,994 | Protein localization to organelle; nucleocytoplasmic transport |

Important distinction: GSEA used the full mapped ranked list of 4,994 genes. It did not use the `n_sites >= 2` candidate cutoff, and it did not require Top100/Top200/Top400 lists.

## 7. AlphaFold3 Structural Modeling

### Confirmed local files

| Item | Value |
|---|---|
| AlphaFold3 input-format file | `E:/gradthesis/curvecurator/doc/AF3_FORMAT.csv` |
| AF3 model output folder | `E:/gradthesis/curvecurator/doc/Q5NVN0_KPYM_PONAB (2)` |
| Ranking-score CSV files found | 113 |
| CIF model files found | 1,243 |

### AF3 input fields observed

`AF3_FORMAT.csv` contains protein name/accession, full protein sequence, local motif sequence, and cysteine position. Example proteins include TXNDC5, PSME1, IPO5, and DHX57.

### Version status

The exact AlphaFold3 software/platform version was not recorded in the files inspected. If this was run through the AlphaFold Server or a local AlphaFold3 installation, the final thesis methods should state the platform and access/run date or software version.

## 8. Docking and Structural Interpretation

### Confirmed use in thesis workflow

Structural analysis was used after quantitative prioritization to provide mechanistic context for selected high-confidence cysteine sites. Representative examples discussed in the defense include:

| Protein | Site | Role in interpretation |
|---|---:|---|
| IPO5 | Cys682 | Transport-associated GSEA driver linked to a strict-hit cysteine site |
| IPO7 | Cys477 | Transport-associated GSEA driver linked to RNA localization/nucleocytoplasmic transport |
| PSME1 | Cys106 | High-confidence strict-hit docking example |
| DHX57 | Cys1369 | High-confidence strict-hit structural/docking example |

### Docking software/version status

The exact docking software name, version, scoring settings, grid size, covalent docking protocol, and non-covalent docking protocol were not clearly recorded in the project files inspected. These should be filled in from the actual docking run logs or the software used for the final thesis methods section.

## 9. Python Support Scripts Used During Thesis/Defense Preparation

These scripts were used for data auditing, plotting, and defense figure generation in the Codex working folder.

| Purpose | Packages/imports observed |
|---|---|
| CurveCurator input audit | pandas, numpy |
| Dose-response plotting | pandas, numpy, matplotlib |
| PNG schematic generation | Pillow (`PIL`) |
| DOCX generation | python-docx |
| XML/text extraction from DOCX/PPTX | Python standard library: zipfile, xml.etree.ElementTree, pathlib, re, csv, json |

### Python environment used for support scripts

| Software / package | Version |
|---|---|
| Python | 3.13.12 |
| pandas | 3.0.1 |
| numpy | 2.4.3 |
| matplotlib | 3.10.8 |
| scipy | 1.17.1 |

Pillow and python-docx were used by local support scripts, but their versions were not separately recorded in the current summary.

## 10. Parameters to Report Clearly in the Defense

### Candidate branch for enrichment

| Parameter | Value | Purpose |
|---|---|---|
| `best_pEC50 >= 4` | EC50 <= 100 uM | Retain genes with at least moderate potency |
| `n_sites >= 2` | At least two responsive sites | Increase gene-level evidence strength |
| TopN | 100, 200, 400 | Define ORA foreground lists |
| ORA p-value/FDR cutoffs | GO p < 0.05, GO FDR/q < 0.2; KEGG p < 0.05; Reactome p < 0.05 | Select enriched pathway terms |

### GSEA branch

| Parameter | Value | Purpose |
|---|---|---|
| Ranked gene list | 4,994 mapped genes | Avoid hard TopN cutoff |
| Ranking strategies | `best_pEC50`, `mean_pEC50`, `weighted_score` | Robustness check |
| `minGSSize` / `maxGSSize` | 10 / 500 | Gene-set size filter |
| `nPermSimple` | 10,000 | Permutation-based null distribution |
| `pvalueCutoff` | 0.2 | Exploratory pathway-level significance threshold |

### Strict hit branch

| Parameter | Original thesis setting | Revised directional sensitivity check |
|---|---|---|
| Adjusted p-value | <= 0.05 | <= 0.05 |
| pEC50 | >= 4 | >= 4 |
| R2 | >= 0.99 | >= 0.99 |
| Fold change | `|log2FC| >= 2` | `log2FC >= 1` |
| Result | 220 sites / 189 genes | 75 sites / 74 genes |

## 11. Missing Items to Verify Before Final Thesis Submission

The following items were not fully recoverable from the files inspected:

1. Exact LC-MS/MS instrument model, acquisition software, and raw search/export software version.
2. Exact TMTMosaic software version.
3. Exact AlphaFold3 platform/version and run date.
4. Exact docking software, version, grid/protocol settings, scoring settings, and whether covalent docking was constrained to the target cysteine.
5. `readr` package version from the original R environment.
6. Full `sessionInfo()` from the original R enrichment runs.

These should be added if required by the final thesis Methods or defense appendix.
