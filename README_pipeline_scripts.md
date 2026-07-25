# Thesis Pipeline Scripts

This folder collects the R, Python, SQL, and CurveCurator configuration files used or audited for the thesis analysis.

## Folder structure

- SQL/: SQL scripts for raw-data cleaning and checks.
- Python/: Python scripts for preprocessing, data auditing, plotting, and figure generation.
- R/: R scripts for ORA/GSEA-related enrichment analysis and sensitivity reruns.
- CurveCurator_configs/: TOML configuration files used for CurveCurator 4PL fitting.
- manifests/script_manifest.csv: table linking each copied file to its original source and purpose.

## Notes

- Original files in the root script folder were not deleted.
- Codex-generated helper scripts are prefixed with codex_.
- CurveCurator output files, CSV results, figures, logs, and large package libraries were not copied here unless they are configuration/code files.
- The main parameter summary is available in 	hesis_software_packages_parameters_summary.md.
