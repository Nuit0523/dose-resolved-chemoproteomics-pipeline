from pathlib import Path

import numpy as np
import pandas as pd


raw_file = Path(r"E:\R\TMT_analysis\TMT4\Raw_sq_8171_TMTMosaic.csv")
meta_file = Path(r"E:\R\TMT_analysis\TMT4\direct_labeling_desthiobiotin_pipeline_8171\01_long_and_metadata\tmt_column_metadata.csv")
mean_file = Path(r"E:\R\TMT_analysis\TMT4\direct_labeling_desthiobiotin_pipeline_8171\02_summary_tables\site_level_mean_8171.csv")
conc_full_file = Path(r"E:\gradthesis\curvecurator\8171\8171_conc_full.csv")
conc_file = Path(r"E:\gradthesis\curvecurator\8171\8171_conc.txt")
curve_file = Path(r"E:\gradthesis\curvecurator\8171\8171_conc_curve.tsv")
norm_file = Path(r"E:\gradthesis\curvecurator\8171\8171_normalization_factors.tsv")

out_dir = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\curvecurator_input_audit")
out_dir.mkdir(parents=True, exist_ok=True)

raw = pd.read_csv(raw_file)
meta = pd.read_csv(meta_file)
means = pd.read_csv(mean_file)
conc = pd.read_csv(conc_file, sep="\t")
conc_full = pd.read_csv(conc_full_file)
curves = pd.read_csv(curve_file, sep="\t")
norm = pd.read_csv(norm_file, sep="\t", header=None, names=["channel", "log2_factor"])

raw = raw.rename(columns={"Site Position": "Site.Position"})
raw["Name"] = raw["gene_symbol"].astype(str) + "_" + raw["site_id"].astype(str) + "_" + raw["Site.Position"].astype(str)

checks = []
for dose_idx, dose in enumerate([20, 40, 80, 160, 320, 640], start=1):
    dose_cols = meta.loc[meta["concentration"] == dose, "column"].tolist()
    raw[f"calc_mean_{dose}"] = raw[dose_cols].mean(axis=1, skipna=True)
    raw_subset = raw[["Name", f"calc_mean_{dose}"]].copy()
    merged = conc[["Name", f"Raw {dose_idx}"]].merge(raw_subset, on="Name", how="inner")
    diff = (merged[f"Raw {dose_idx}"] - merged[f"calc_mean_{dose}"]).abs()
    checks.append(
        {
            "dose_uM": dose,
            "raw_channels": ";".join(dose_cols),
            "n_matched_rows": len(merged),
            "max_abs_difference": float(diff.max()),
            "median_abs_difference": float(diff.median()),
        }
    )

check_df = pd.DataFrame(checks)
check_df.to_csv(out_dir / "raw_to_curvecurator_mean_check.csv", index=False)

example_name = "PKM_0_49"
example_raw = raw.loc[raw["Name"] == example_name].iloc[0]
example = []
for dose_idx, dose in enumerate([20, 40, 80, 160, 320, 640], start=1):
    cols = meta.loc[meta["concentration"] == dose, "column"].tolist()
    vals = [example_raw[c] for c in cols]
    cc_val = conc.loc[conc["Name"] == example_name, f"Raw {dose_idx}"].iloc[0]
    example.append(
        {
            "Name": example_name,
            "dose_uM": dose,
            "channels": ";".join(cols),
            "raw_channel_values": ";".join(f"{v:.4f}" for v in vals),
            "manual_mean": np.mean(vals),
            "curvecurator_input_Raw": cc_val,
        }
    )
pd.DataFrame(example).to_csv(out_dir / "example_PKM_0_49_raw_mean_to_input.csv", index=False)

site_counts = {
    "raw_rows": len(raw),
    "site_level_mean_rows": len(means),
    "site_level_mean_unique_sites": means["site_id"].nunique(),
    "curvecurator_input_rows": len(conc),
    "curvecurator_full_rows": len(conc_full),
    "curvecurator_curve_rows": len(curves),
}
pd.DataFrame([site_counts]).to_csv(out_dir / "row_count_summary.csv", index=False)

classified = pd.read_csv(
    Path(r"E:\R\TMT_analysis\TMT4\direct_labeling_desthiobiotin_pipeline_8171\TMT_protein_classification_round22_8171.csv")
)
classified = classified.rename(columns={"Site Position": "Site.Position"})
classified["Name"] = (
    classified["gene_symbol"].astype(str)
    + "_"
    + classified["site_id"].astype(str)
    + "_"
    + classified["Site.Position"].astype(str)
)
missing_from_input = classified.loc[~classified["Name"].isin(set(conc["Name"]))].copy()
missing_from_input.to_csv(out_dir / "classified_rows_missing_from_curvecurator_input.csv", index=False)

raw_channel_cols = meta["column"].tolist()
missing_reason_rows = []
for _, row in missing_from_input.iterrows():
    values = pd.to_numeric(row[raw_channel_cols], errors="coerce")
    missing_reason_rows.append(
        {
            "Name": row["Name"],
            "n_missing_channels": int(values.isna().sum()),
            "n_nonpositive_channels": int((values <= 0).sum()),
            "min_value": values.min(skipna=True),
            "max_value": values.max(skipna=True),
        }
    )
pd.DataFrame(missing_reason_rows).to_csv(out_dir / "classified_rows_missing_reason_summary.csv", index=False)

norm_example = curves.loc[curves["Name"] == example_name].iloc[0]
norm_rows = []
for dose_idx in range(1, 7):
    raw_val = norm_example[f"Raw {dose_idx}"]
    normalized_val = norm_example[f"Normalized {dose_idx}"]
    factor = norm.loc[norm["channel"] == f"Raw {dose_idx}", "log2_factor"].iloc[0]
    reconstructed = raw_val * (2 ** factor)
    norm_rows.append(
        {
            "Name": example_name,
            "channel": f"Raw {dose_idx}",
            "raw_input_value": raw_val,
            "normalization_log2_factor": factor,
            "raw_times_2pow_factor": reconstructed,
            "CurveCurator_normalized_value": normalized_val,
            "difference": abs(reconstructed - normalized_val),
        }
    )
pd.DataFrame(norm_rows).to_csv(out_dir / "example_PKM_0_49_curvecurator_normalization_check.csv", index=False)

print("Wrote", out_dir)
print(check_df.to_string(index=False))
print(pd.DataFrame([site_counts]).to_string(index=False))
