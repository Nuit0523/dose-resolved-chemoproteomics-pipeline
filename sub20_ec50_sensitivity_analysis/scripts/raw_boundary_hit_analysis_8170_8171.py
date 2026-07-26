from pathlib import Path

import numpy as np
import pandas as pd


LOWER_BOUNDARY_PEC50_20UM = -np.log10(20e-6)
UPPER_BOUNDARY_PEC50_640UM = -np.log10(640e-6)

TMT_CHANNELS = [
    "default~126_sn_sum",
    "default~127n_sn_sum",
    "default~127c_sn_sum",
    "default~128n_sn_sum",
    "default~128c_sn_sum",
    "default~129n_sn_sum",
    "default~129c_sn_sum",
    "default~130n_sn_sum",
    "default~130c_sn_sum",
    "default~131n_sn_sum",
    "default~131c_sn_sum",
    "default~132n_sn_sum",
    "default~132c_sn_sum",
    "default~133n_sn_sum",
    "default~133c_sn_sum",
    "default~134n_sn_sum",
    "default~134c_sn_sum",
    "default~135_sn_sum",
]

# Channel order in the raw TMTMosaic files is arranged as six dose groups with
# three channels per dose. This reproduces the previously audited CurveCurator
# input construction.
DOSES_UM = [20, 40, 80, 160, 320, 640]
DOSE_GROUPS = {
    dose: TMT_CHANNELS[i * 3 : (i + 1) * 3] for i, dose in enumerate(DOSES_UM)
}


def site_name_from_raw(df):
    return df["gene_symbol"].astype(str) + "_" + df["site_id"].astype(str) + "_" + df["Site Position"].astype(str)


def summarize_raw(raw_path):
    raw = pd.read_csv(raw_path)
    out = raw[["site_id", "Protein Id", "gene_symbol", "prot_description", "Site Position", "Motif"]].copy()
    out["Name"] = site_name_from_raw(raw)

    for dose, cols in DOSE_GROUPS.items():
        values = raw[cols].apply(pd.to_numeric, errors="coerce")
        out[f"raw_mean_{dose}uM"] = values.mean(axis=1, skipna=True)
        out[f"raw_sd_{dose}uM"] = values.std(axis=1, skipna=True)
        out[f"raw_cv_{dose}uM"] = out[f"raw_sd_{dose}uM"] / out[f"raw_mean_{dose}uM"]
        out[f"missing_n_{dose}uM"] = values.isna().sum(axis=1)

    mean_cols = [f"raw_mean_{dose}uM" for dose in DOSES_UM]
    cv_cols = [f"raw_cv_{dose}uM" for dose in DOSES_UM]
    out["raw_mean_signal_overall"] = out[mean_cols].mean(axis=1, skipna=True)
    out["raw_median_cv"] = out[cv_cols].median(axis=1, skipna=True)
    out["raw_log2FC_640_vs_20"] = np.log2(out["raw_mean_640uM"] / out["raw_mean_20uM"])
    out["raw_dynamic_range_log2_max_min"] = np.log2(out[mean_cols].max(axis=1) / out[mean_cols].min(axis=1))
    out["raw_missing_total"] = out[[f"missing_n_{dose}uM" for dose in DOSES_UM]].sum(axis=1)
    return out


def load_curve(curve_path):
    curve = pd.read_csv(curve_path, sep="\t")
    keep = [
        "Name",
        "pEC50",
        "Curve Fold Change",
        "Curve RMSE",
        "Curve R2",
        "Curve P_Value adjusted",
        "Curve Log P_Value adjusted",
        "Curve Regulation",
    ]
    keep = [col for col in keep if col in curve.columns]
    return curve[keep].copy()


def load_ec50(ec50_path):
    ec50 = pd.read_csv(ec50_path)
    return ec50[["Name", "pEC50"]].copy()


def classify_range(pec50):
    if pd.isna(pec50):
        return "missing"
    if pec50 > LOWER_BOUNDARY_PEC50_20UM:
        return "detected_EC50_below_20uM"
    if pec50 >= UPPER_BOUNDARY_PEC50_640UM:
        return "fitted_EC50_in_20_640uM"
    return "fitted_EC50_above_640uM_or_weak"


def add_hit_flags(df):
    out = df.copy()
    out["estimated_EC50_uM"] = 10 ** (-out["pEC50"]) * 1e6
    out["ec50_range_status"] = out["pEC50"].apply(classify_range)
    out["detected_sub20_EC50"] = out["pEC50"] > LOWER_BOUNDARY_PEC50_20UM

    adjp = out.get("Curve P_Value adjusted", pd.Series(np.nan, index=out.index))
    r2 = out.get("Curve R2", pd.Series(np.nan, index=out.index))
    fc = out.get("Curve Fold Change", pd.Series(np.nan, index=out.index))

    out["strict_hit_absFC2"] = (adjp <= 0.05) & (out["pEC50"] >= 4) & (r2 >= 0.99) & (fc.abs() >= 2)
    out["strict_hit_positiveFC1"] = (adjp <= 0.05) & (out["pEC50"] >= 4) & (r2 >= 0.99) & (fc >= 1)
    out["sub20_strict_hit_absFC2"] = out["detected_sub20_EC50"] & out["strict_hit_absFC2"]
    out["sub20_strict_hit_positiveFC1"] = out["detected_sub20_EC50"] & out["strict_hit_positiveFC1"]
    out["near_lower_boundary_pEC50_4p5_to_4p699"] = (out["pEC50"] >= 4.5) & (out["pEC50"] <= LOWER_BOUNDARY_PEC50_20UM)
    return out


def count_flags(df, dataset):
    rows = []
    flags = [
        "detected_sub20_EC50",
        "strict_hit_absFC2",
        "strict_hit_positiveFC1",
        "sub20_strict_hit_absFC2",
        "sub20_strict_hit_positiveFC1",
        "near_lower_boundary_pEC50_4p5_to_4p699",
    ]
    for flag in flags:
        n = int(df[flag].sum())
        rows.append(
            {
                "dataset": dataset,
                "metric": flag,
                "site_count": n,
                "gene_count": int(df.loc[df[flag], "gene_symbol"].nunique()),
                "percent_sites": n / len(df) * 100,
            }
        )
    for status, sub in df.groupby("ec50_range_status"):
        rows.append(
            {
                "dataset": dataset,
                "metric": f"range_status::{status}",
                "site_count": len(sub),
                "gene_count": int(sub["gene_symbol"].nunique()),
                "percent_sites": len(sub) / len(df) * 100,
            }
        )
    return pd.DataFrame(rows)


def process_dataset(dataset, raw_path, curve_path, ec50_path=None, curve_merge_mode="name"):
    raw_summary = summarize_raw(raw_path)
    if ec50_path is not None:
        ec50 = load_ec50(ec50_path)
        merged = raw_summary.merge(ec50, on="Name", how="left")
    else:
        merged = raw_summary.copy()

    curve = load_curve(curve_path)
    curve_no_pec50 = curve.drop(columns=["pEC50"], errors="ignore")
    if curve_merge_mode == "site_id":
        curve_no_pec50 = curve_no_pec50.rename(columns={"Name": "site_id"})
        curve_no_pec50["site_id"] = pd.to_numeric(curve_no_pec50["site_id"], errors="coerce")
        merged["site_id"] = pd.to_numeric(merged["site_id"], errors="coerce")
        merged = merged.merge(curve_no_pec50, on="site_id", how="left")
    else:
        if "pEC50" not in merged.columns:
            merged = raw_summary.merge(curve, on="Name", how="left")
        else:
            merged = merged.merge(curve_no_pec50, on="Name", how="left")
    merged = add_hit_flags(merged)
    return merged, count_flags(merged, dataset)


def main():
    out_dir = Path("C:/Users/owner/Documents/Codex/2026-07-07/ni-ha/outputs/raw_boundary_hit_analysis")
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = {
        "8170": {
            "raw": Path("E:/gradthesis/Raw_sq_8170_TMTMosaic.csv"),
            "curve": Path("E:/gradthesis/curvecurator/8170/row/8170_output_curves_row.tsv"),
            "ec50": Path("E:/gradthesis/curvecurator/8170/8170_EC50.csv"),
            "curve_merge_mode": "site_id",
        },
        "8171": {
            "raw": Path("E:/gradthesis/Raw_sq_8171_TMTMosaic.csv"),
            "curve": Path("E:/gradthesis/curvecurator/8171/8171_conc_curve.tsv"),
            "ec50": None,
            "curve_merge_mode": "name",
        },
    }

    all_counts = []
    for dataset, paths in datasets.items():
        merged, counts = process_dataset(
            dataset,
            paths["raw"],
            paths["curve"],
            ec50_path=paths.get("ec50"),
            curve_merge_mode=paths.get("curve_merge_mode", "name"),
        )
        merged.to_csv(out_dir / f"{dataset}_raw_curve_boundary_hit_annotation.csv", index=False)
        all_counts.append(counts)

    summary = pd.concat(all_counts, ignore_index=True)
    summary.to_csv(out_dir / "raw_boundary_hit_summary.csv", index=False)

    md = [
        "# Raw Data + CurveCurator Boundary Hit Analysis",
        "",
        f"20 uM lower boundary pEC50 = {LOWER_BOUNDARY_PEC50_20UM:.6f}.",
        f"640 uM upper boundary pEC50 = {UPPER_BOUNDARY_PEC50_640UM:.6f}.",
        "",
        "A real-data site was counted as detected below 20 uM only if fitted pEC50 > 4.699.",
        "",
        "```text",
        summary.to_string(index=False),
        "```",
        "",
        "Important: this analysis counts detected below-range fitted EC50 values. It cannot identify true missed sub-20 uM sites from raw data alone.",
    ]
    (out_dir / "raw_boundary_hit_summary.md").write_text("\n".join(md), encoding="utf-8")

    print(summary.to_string(index=False))
    print(f"\nWrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
