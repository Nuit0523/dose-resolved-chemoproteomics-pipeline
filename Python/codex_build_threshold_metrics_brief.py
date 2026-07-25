from __future__ import annotations

import csv
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha")
OUT = ROOT / "outputs" / "thesis_threshold_metrics_brief.md"
GRAD = Path(r"E:\gradthesis")
CURVE = GRAD / "curvecurator"
TMT = Path(r"E:\R\TMT_analysis\TMT4")


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def as_float(value: str | None) -> float | None:
    if value in (None, "", "NA", "NaN", "nan"):
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    return out if math.isfinite(out) else None


def raw_profile(path: Path) -> dict:
    rows = read_csv(path)
    cols = list(rows[0].keys()) if rows else []
    tmt_cols = [
        c for c in cols
        if c.startswith("default~") and c.endswith("_sn_sum")
        or re.match(r"^\d+(n|c)?_wxr35_", c)
    ]
    nulls = Counter()
    for row in rows:
        for col in cols:
            if row.get(col) in ("", "NA", "NaN", None):
                nulls[col] += 1
    return {
        "rows": len(rows),
        "columns": len(cols),
        "tmt_intensity_columns": len(tmt_cols),
        "missing_gene_symbol": nulls.get("gene_symbol", 0),
        "missing_best_scan": nulls.get("best_scan", 0),
        "missing_quantified_peptides": nulls.get("default~quantified_peptides", 0),
        "top_missing": nulls.most_common(8),
    }


def long_table_profile(path: Path) -> dict:
    rows = read_csv(path)
    concentrations = sorted({as_float(r.get("concentration")) for r in rows if as_float(r.get("concentration")) is not None})
    sites = {r.get("site_id") for r in rows}
    genes = {r.get("gene_symbol") for r in rows if r.get("gene_symbol")}
    reps = Counter(r.get("concentration") for r in rows)
    nulls = Counter()
    numeric_zero = 0
    numeric_negative = 0
    for row in rows:
        for col in ["intensity", "log2_intensity", "log2fc", "log2FC_vs_baseline"]:
            val = as_float(row.get(col))
            if val is None:
                if col in row:
                    nulls[col] += 1
                continue
            if val == 0:
                numeric_zero += 1
            elif val < 0:
                numeric_negative += 1
    return {
        "rows": len(rows),
        "unique_sites": len(sites),
        "unique_genes": len(genes),
        "concentrations": concentrations,
        "rows_per_concentration": dict(reps),
        "numeric_zero_values_in_checked_cols": numeric_zero,
        "numeric_negative_values_in_checked_cols": numeric_negative,
        "missing_checked_cols": dict(nulls),
    }


def curve_input_profile(path: Path) -> dict:
    rows = read_csv(path, delimiter="\t")
    raw_cols = [c for c in rows[0].keys() if c.startswith("Raw ")] if rows else []
    row_missing = Counter()
    negative_rows = 0
    zero_rows = 0
    all_values = []
    for row in rows:
        miss = 0
        has_negative = False
        has_zero = False
        for col in raw_cols:
            val = as_float(row.get(col))
            if val is None:
                miss += 1
            else:
                all_values.append(val)
                has_negative = has_negative or val < 0
                has_zero = has_zero or val == 0
        row_missing[miss] += 1
        negative_rows += int(has_negative)
        zero_rows += int(has_zero)
    return {
        "rows": len(rows),
        "raw_columns": raw_cols,
        "row_missing_distribution": dict(sorted(row_missing.items())),
        "rows_with_more_than_2_missing": sum(n for miss, n in row_missing.items() if miss > 2),
        "rows_with_negative_values": negative_rows,
        "rows_with_zero_values": zero_rows,
        "min_raw": round(min(all_values), 4) if all_values else None,
        "median_raw": round(sorted(all_values)[len(all_values) // 2], 4) if all_values else None,
        "max_raw": round(max(all_values), 4) if all_values else None,
    }


def curve_profile(path: Path) -> dict:
    rows = read_csv(path, delimiter="\t")
    regs = Counter((r.get("Curve Regulation") or "blank") for r in rows)
    nums: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        for col in ["pEC50", "Curve R2", "Curve Fold Change", "Signal Quality", "Curve RMSE", "Curve Slope", "Curve P_Value", "Curve Log P_Value adjusted"]:
            val = as_float(row.get(col))
            if val is not None:
                nums[col].append(val)

    def stats(col: str) -> dict:
        values = sorted(nums[col])
        return {
            "min": round(values[0], 5),
            "median": round(values[len(values) // 2], 5),
            "mean": round(sum(values) / len(values), 5),
            "max": round(values[-1], 5),
        } if values else {}

    hit_rows = [
        r for r in rows
        if (as_float(r.get("Curve Log P_Value adjusted")) or -999) >= 1.30103
        and (as_float(r.get("pEC50")) or -999) >= 4
        and (as_float(r.get("Curve R2")) or -999) >= 0.99
        and abs(as_float(r.get("Curve Fold Change")) or 0) >= 2
    ]
    genes = [r.get("Name", "").split("_")[0] for r in hit_rows]
    gene_counts = Counter(genes)
    return {
        "rows": len(rows),
        "regulation": dict(regs),
        "stats": {col: stats(col) for col in nums},
        "pEC50_ge_4": sum((as_float(r.get("pEC50")) or -999) >= 4 for r in rows),
        "r2_lt_0_8": sum((as_float(r.get("Curve R2")) or 999) < 0.8 for r in rows),
        "r2_lt_0_5": sum((as_float(r.get("Curve R2")) or 999) < 0.5 for r in rows),
        "rmse_gt_0_2": sum((as_float(r.get("Curve RMSE")) or -999) > 0.2 for r in rows),
        "slope_eq_10": sum((as_float(r.get("Curve Slope")) or -999) == 10 for r in rows),
        "padj_log_ge_1_30103_and_pEC50_ge_4": sum(
            (as_float(r.get("Curve Log P_Value adjusted")) or -999) >= 1.30103
            and (as_float(r.get("pEC50")) or -999) >= 4 for r in rows
        ),
        "strict_hit_sites": len(hit_rows),
        "strict_unique_genes": len(gene_counts),
        "strict_single_site_genes": sum(1 for n in gene_counts.values() if n == 1),
        "strict_multi_site_genes": sum(1 for n in gene_counts.values() if n > 1),
        "strict_top20": sorted(
            [
                {
                    "gene": r.get("Name", "").split("_")[0],
                    "name": r.get("Name"),
                    "pEC50": as_float(r.get("pEC50")),
                    "fold_change": as_float(r.get("Curve Fold Change")),
                    "r2": as_float(r.get("Curve R2")),
                    "log_p_adj": as_float(r.get("Curve Log P_Value adjusted")),
                    "regulation": r.get("Curve Regulation"),
                }
                for r in hit_rows
            ],
            key=lambda r: (r["pEC50"] or -999, r["log_p_adj"] or -999, r["r2"] or -999, abs(r["fold_change"] or 0)),
            reverse=True,
        )[:20],
    }


def config_thresholds(path: Path) -> dict:
    txt = path.read_text(encoding="utf-8-sig", errors="replace")
    keys = [
        "doses", "dose_scale", "dose_unit", "control_experiment", "imputation",
        "imputation_pct", "max_missing", "max_imputation", "normalization",
        "max_iterations", "alpha", "fc_lim", "mtc_method", "not_rmse_limit",
    ]
    out = {}
    for key in keys:
        m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.+)$", txt, re.M)
        if m:
            out[key] = m.group(1).strip()
    return out


def docx_threshold_mentions(path: Path) -> list[str]:
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    lines = []
    for para in root.findall(".//w:p", ns):
        line = "".join(t.text or "" for t in para.findall(".//w:t", ns)).strip()
        if re.search(r"threshold|cutoff|filter|missing|imputation|FDR|pEC|R2|R²|fold|MAD|normalization|quality", line, re.I):
            lines.append(line)
    return lines[:120]


def write_brief() -> None:
    raw8170 = raw_profile(TMT / "Raw_sq_8170_TMTMosaic.csv")
    raw8171 = raw_profile(TMT / "Raw_sq_8171_TMTMosaic.csv")
    long8170 = long_table_profile(TMT / "direct_labeling_desthiobiotin_pipeline_8170" / "01_long_and_metadata" / "site_level_long_table.csv")
    long8171 = long_table_profile(TMT / "direct_labeling_desthiobiotin_pipeline_8171" / "01_long_and_metadata" / "site_level_long_table.csv")
    inp8171 = curve_input_profile(CURVE / "8171" / "8171_conc.txt")
    curve8171 = curve_profile(CURVE / "8171" / "8171_conc_curve.tsv")
    conf8171 = config_thresholds(CURVE / "8171" / "8171_conc_config.toml")
    conf8170 = config_thresholds(CURVE / "8170" / "8170_conc_config.toml")
    mentions = docx_threshold_mentions(ROOT / "work" / "thesis.docx") if (ROOT / "work" / "thesis.docx").exists() else []

    lines = ["# Data Cleaning, Thresholds, and QC Metrics Brief\n"]
    lines.append("## Why This Belongs in the Defense\n")
    lines.append("This section explains how raw TMTMosaic outputs were converted into reliable site-level curves and how thresholds were used to move from all quantified curves to candidate hit sites and genes.\n")

    lines.append("## 1. Raw Input QC\n")
    lines.append("| Dataset | Raw rows | Columns | TMT intensity columns | Missing gene symbols | Missing best scans | Missing quantified peptide lists |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    lines.append(f"| 8170 | {raw8170['rows']:,} | {raw8170['columns']} | {raw8170['tmt_intensity_columns']} | {raw8170['missing_gene_symbol']} | {raw8170['missing_best_scan']} | {raw8170['missing_quantified_peptides']} |")
    lines.append(f"| 8171 | {raw8171['rows']:,} | {raw8171['columns']} | {raw8171['tmt_intensity_columns']} | {raw8171['missing_gene_symbol']} | {raw8171['missing_best_scan']} | {raw8171['missing_quantified_peptides']} |")
    lines.append("- These checks show that most rows have core annotation and scan metadata, while quantified peptide lists are the most commonly incomplete field.\n")

    lines.append("## 2. Long Table Processing Metrics\n")
    lines.append("| Dataset | Long-table rows | Unique sites | Unique genes | Concentrations |")
    lines.append("|---|---:|---:|---:|---|")
    lines.append(f"| 8170 | {long8170['rows']:,} | {long8170['unique_sites']:,} | {long8170['unique_genes']:,} | {long8170['concentrations']} |")
    lines.append(f"| 8171 | {long8171['rows']:,} | {long8171['unique_sites']:,} | {long8171['unique_genes']:,} | {long8171['concentrations']} |")
    lines.append("- The processing script groups each site by concentration and averages replicate intensities before pivoting into CurveCurator's six-column input format.")
    lines.append("- The six dose levels are 20, 40, 80, 160, 320, and 640 uM.\n")

    lines.append("## 3. CurveCurator Preprocessing Thresholds\n")
    lines.append("| Parameter | 8170 | 8171 | Meaning |")
    lines.append("|---|---|---|---|")
    for key in ["doses", "control_experiment", "imputation", "imputation_pct", "max_missing", "max_imputation", "normalization", "max_iterations", "alpha", "fc_lim", "mtc_method", "not_rmse_limit"]:
        lines.append(f"| `{key}` | `{conf8170.get(key, '')}` | `{conf8171.get(key, '')}` | |")
    lines.append("- Key cleaning rule for curve fitting: CurveCurator allowed imputation, but rows with more than 2 missing dose values were flagged as likely to be filtered.")
    lines.append("- `max_missing = 2` and `max_imputation = 2` were the explicit missing-data limits.")
    lines.append("- `alpha = 0.05`, `mtc_method = fdr_bh`, and `fc_lim = 1.0` define the CurveCurator statistical testing and fold-change framework.\n")

    lines.append("## 4. 8171 CurveCurator Input QC\n")
    lines.append(f"- CurveCurator input rows: {inp8171['rows']:,}.")
    lines.append(f"- Missing-dose distribution across Raw 1 to Raw 6: {inp8171['row_missing_distribution']}.")
    lines.append(f"- Rows with more than 2 missing dose values: {inp8171['rows_with_more_than_2_missing']}.")
    lines.append(f"- Rows with negative raw values: {inp8171['rows_with_negative_values']}; rows with zero raw values: {inp8171['rows_with_zero_values']}.")
    lines.append(f"- Raw intensity range after formatting: min {inp8171['min_raw']}, median {inp8171['median_raw']}, max {inp8171['max_raw']}.\n")

    lines.append("## 5. 8171 Curve-Level QC Metrics\n")
    lines.append(f"- Total fitted/curated curves: {curve8171['rows']:,}.")
    lines.append(f"- Regulation labels: {curve8171['regulation']}.")
    for col, stats in curve8171["stats"].items():
        lines.append(f"- {col}: min {stats['min']}, median {stats['median']}, mean {stats['mean']}, max {stats['max']}.")
    lines.append(f"- Curves with pEC50 >= 4: {curve8171['pEC50_ge_4']:,}.")
    lines.append(f"- Curves with R2 < 0.8: {curve8171['r2_lt_0_8']:,}; R2 < 0.5: {curve8171['r2_lt_0_5']:,}.")
    lines.append(f"- Curves with RMSE > 0.2: {curve8171['rmse_gt_0_2']:,}.")
    lines.append(f"- Curves where slope hit the upper bound of 10: {curve8171['slope_eq_10']:,}.")
    lines.append(f"- Curves passing adjusted log-p >= 1.30103 and pEC50 >= 4: {curve8171['padj_log_ge_1_30103_and_pEC50_ge_4']:,}.\n")

    lines.append("## 6. Final Strict Hit Threshold\n")
    lines.append("The strict hit-site SQL view used the following combined criteria:")
    lines.append("- `Curve Log P_Value adjusted >= 1.30103`, equivalent to adjusted p <= 0.05 on a -log10 scale.")
    lines.append("- `pEC50 >= 4`.")
    lines.append("- `Curve R2 >= 0.99`.")
    lines.append("- `abs(Curve Fold Change) >= 2`.")
    lines.append(f"Applying these criteria to 8171 gives {curve8171['strict_hit_sites']:,} strict hit sites from {curve8171['strict_unique_genes']:,} unique genes.")
    lines.append(f"- Single-site hit genes: {curve8171['strict_single_site_genes']:,}.")
    lines.append(f"- Multi-site hit genes: {curve8171['strict_multi_site_genes']:,}.")
    lines.append("- Top strict hit sites by pEC50:")
    lines.append("| Gene | Site name | pEC50 | Fold change | R2 | adj log-p | Regulation |")
    lines.append("|---|---|---:|---:|---:|---:|---|")
    for row in curve8171["strict_top20"][:10]:
        lines.append(f"| {row['gene']} | {row['name']} | {row['pEC50']:.4f} | {row['fold_change']:.4f} | {row['r2']:.4f} | {row['log_p_adj']:.4f} | {row['regulation'] or ''} |")
    lines.append("")

    lines.append("## 7. How to Present This in the PPT\n")
    lines.append("- Add a dedicated slide called `Data cleaning converted raw TMT channels into analysis-ready dose curves`.")
    lines.append("- Show the path: raw TMTMosaic table -> long table -> replicate mean by concentration -> CurveCurator input -> normalized fitted curves -> strict hit sites.")
    lines.append("- Put the final strict hit thresholds in one clean callout, because these are the numbers committee members are most likely to ask about.")
    lines.append("- Present 8170 as the test case for verifying this processing logic, then 8171 as the scaled production run.\n")

    if mentions:
        lines.append("## 8. Thesis Text Mentions to Reconcile\n")
        lines.append("These lines from the thesis text mention thresholds/QC and should be checked against the final slide wording:")
        for m in mentions[:25]:
            clean = m.replace("|", "\\|")
            lines.append(f"- {clean}")
        lines.append("")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(str(OUT))


if __name__ == "__main__":
    write_brief()
