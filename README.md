# Dose-Resolved Chemoproteomics Analysis Scripts

This repository contains the analysis scripts and configuration files used for a thesis project on dose-resolved cysteine-focused chemoproteomic profiling.

## Repository structure

- SQL/: SQL scripts used for raw-data cleaning and quality checks.
- Python/: Python scripts used for preprocessing checks, data auditing, plotting, and defense figure generation.
- R/: R scripts used for ORA/GSEA-related enrichment analysis and sensitivity reruns.
- CurveCurator_configs/: CurveCurator TOML configuration files for 8170 and 8171 dose-response modeling.
- manifests/script_manifest.csv: file-level manifest linking each script to its source and purpose.
- README_pipeline_scripts.md: additional notes generated during script organization.

## Main analysis logic

The computational workflow follows these major steps:

1. Raw TMT/LC-MS output checking and cleaning.
2. Dose-level aggregation of reporter-ion signal.
3. CurveCurator 4PL dose-response modeling.
4. Curve-level quality control using pEC50, curve fold change, R2, RMSE, and adjusted p-value.
5. Site-to-gene summarization for pathway analysis.
6. ORA using TopN candidate gene lists.
7. GSEA using full ranked gene lists.
8. Structure-guided prioritization of selected high-confidence cysteine sites.

## Reproducibility notes

Several scripts contain project-local paths such as E:/gradthesis/... or C:/Users/.... These paths document the original analysis environment and should be edited before running the scripts on another machine.

Raw data files, large result tables, figures, manuscripts, and slide decks are intentionally excluded from Git tracking by .gitignore.

## Software summary

The main software summary is maintained separately in the thesis documentation. Key tools include:

- CurveCurator 0.6.0
- R 4.5.2 with clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, and gsea
- Python with pandas, 
umpy, matplotlib, and Pillow

## Status

This repository is intended as a code archive for thesis defense and reproducibility documentation.

## Sub-20 uM EC50 sensitivity analysis

The folder `sub20_ec50_sensitivity_analysis/` documents an exploratory revision analysis of the 20 uM lower concentration boundary. It includes scripts, summarized outputs, and thesis-ready wording for interpreting below-range EC50 detectability. This analysis is explicitly treated as a sensitivity/limitation assessment, not as evidence for definitive sub-20 uM hits or absence of such hits.
