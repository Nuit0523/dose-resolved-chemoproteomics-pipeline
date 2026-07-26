from pathlib import Path

import pandas as pd


P_EC50_LOWER_BOUNDARY_20UM = -pd.np.log10(20e-6) if hasattr(pd, "np") else 4.698970004336019
P_EC50_UPPER_BOUNDARY_640UM = 3.1938200260161127


def ec50_um_from_pec50(pec50):
    return 10 ** (-pec50) * 1e6


def classify_pec50(pec50):
    if pd.isna(pec50):
        return "missing_pEC50"
    if pec50 > P_EC50_LOWER_BOUNDARY_20UM:
        return "lower_boundary_limited_EC50_below_20uM"
    if pec50 >= P_EC50_UPPER_BOUNDARY_640UM:
        return "in_measured_range_20_640uM"
    return "upper_boundary_or_weak_EC50_above_640uM"


def add_boundary_columns(df, pec50_col, prefix=""):
    out = df.copy()
    status_col = f"{prefix}ec50_range_status" if prefix else "ec50_range_status"
    ec50_col = f"{prefix}estimated_EC50_uM" if prefix else "estimated_EC50_uM"
    boundary_col = f"{prefix}boundary_limited_flag" if prefix else "boundary_limited_flag"
    out[ec50_col] = out[pec50_col].apply(ec50_um_from_pec50)
    out[status_col] = out[pec50_col].apply(classify_pec50)
    out[boundary_col] = out[status_col].eq("lower_boundary_limited_EC50_below_20uM")
    return out


def summarize_status(df, status_col, dataset, level):
    counts = df[status_col].value_counts(dropna=False).reset_index()
    counts.columns = ["ec50_range_status", "count"]
    counts.insert(0, "level", level)
    counts.insert(0, "dataset", dataset)
    counts["percent"] = counts["count"] / len(df) * 100
    return counts


def process_site_level(dataset, ec50_file, curve_file=None):
    site = pd.read_csv(ec50_file)
    site["gene_symbol"] = site["Name"].astype(str).str.split("_").str[0]
    site = add_boundary_columns(site, "pEC50")

    if curve_file and Path(curve_file).exists():
        curve = pd.read_csv(curve_file, sep="\t")
        if "Name" in curve.columns and curve["Name"].astype(str).isin(site["Name"].astype(str)).mean() > 0.5:
            keep_cols = [
                col
                for col in [
                    "Name",
                    "Curve Fold Change",
                    "Curve RMSE",
                    "Curve R2",
                    "Curve P_Value adjusted",
                    "Curve Log P_Value adjusted",
                    "Curve Regulation",
                    "Signal Quality",
                ]
                if col in curve.columns
            ]
            site = site.merge(curve[keep_cols], on="Name", how="left")

    return site


def process_gene_level(dataset, gene_file):
    gene = pd.read_csv(gene_file)
    gene = add_boundary_columns(gene, "best_pEC50", prefix="best_")
    if "mean_pEC50" in gene.columns:
        gene = add_boundary_columns(gene, "mean_pEC50", prefix="mean_")
    if "median_pEC50" in gene.columns:
        gene = add_boundary_columns(gene, "median_pEC50", prefix="median_")
    return gene


def write_outputs():
    out_dir = Path("C:/Users/owner/Documents/Codex/2026-07-07/ni-ha/outputs/pec50_boundary_annotation_20uM")
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = {
        "8171": {
            "ec50": Path("E:/gradthesis/curvecurator/8171/8171_EC50.csv"),
            "curve": Path("E:/gradthesis/curvecurator/8171/8171_conc_curve.tsv"),
            "gene": Path("E:/gradthesis/curvecurator/8171/enrichment_from_raw_name_pEC50/gene_level_pEC50_summary_generated.csv"),
        },
        "8170": {
            "ec50": Path("E:/gradthesis/curvecurator/8170/8170_EC50.csv"),
            "curve": Path("E:/gradthesis/curvecurator/8170/row/8170_output_curves_row.tsv"),
            "gene": Path("E:/gradthesis/curvecurator/8170/gene_level_pEC50_summary.csv"),
        },
    }

    summaries = []
    for dataset, paths in datasets.items():
        site = process_site_level(dataset, paths["ec50"], paths["curve"])
        site_out = out_dir / f"{dataset}_site_pEC50_boundary_annotation.csv"
        site.to_csv(site_out, index=False)
        summaries.append(summarize_status(site, "ec50_range_status", dataset, "site"))

        gene = process_gene_level(dataset, paths["gene"])
        gene_out = out_dir / f"{dataset}_gene_pEC50_boundary_annotation.csv"
        gene.to_csv(gene_out, index=False)
        summaries.append(summarize_status(gene, "best_ec50_range_status", dataset, "gene_best_pEC50"))

    summary = pd.concat(summaries, ignore_index=True)
    summary_out = out_dir / "pec50_boundary_annotation_summary.csv"
    summary.to_csv(summary_out, index=False)

    md_lines = [
        "# pEC50 Boundary Annotation Summary",
        "",
        "Measured concentration range: 20-640 uM.",
        "",
        f"Lower measured boundary: 20 uM -> pEC50 = {P_EC50_LOWER_BOUNDARY_20UM:.3f}.",
        f"Upper measured boundary: 640 uM -> pEC50 = {P_EC50_UPPER_BOUNDARY_640UM:.3f}.",
        "",
        "Classification rules:",
        "",
        "- pEC50 > 4.699: lower-boundary-limited; fitted EC50 is below the lowest measured concentration.",
        "- 3.194 <= pEC50 <= 4.699: fitted EC50 lies within the measured 20-640 uM range.",
        "- pEC50 < 3.194: fitted EC50 is above the highest measured concentration or weak/out-of-range.",
        "",
        "## Counts",
        "",
        "```text",
        summary.to_string(index=False),
        "```",
        "",
        "Interpretation note: lower-boundary-limited curves should not be interpreted as precise sub-20 uM potency estimates.",
    ]
    md_out = out_dir / "pec50_boundary_annotation_summary.md"
    md_out.write_text("\n".join(md_lines), encoding="utf-8")

    print(f"Wrote outputs to: {out_dir}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    write_outputs()
