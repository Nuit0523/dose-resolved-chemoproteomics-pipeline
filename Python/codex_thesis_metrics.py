from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path


TMT = Path(r"E:\R\TMT_analysis\TMT4")
GRAD = Path(r"E:\gradthesis")


def read_rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def as_float(value: str) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def fit_summary(path: Path) -> dict:
    rows = read_rows(path)
    statuses = Counter(row.get("fit_status", "") for row in rows)
    r2s = [as_float(row.get("r_squared", "")) for row in rows]
    r2s = [x for x in r2s if x is not None]
    ec50s = [as_float(row.get("ec50", "")) for row in rows]
    ec50s = [x for x in ec50s if x is not None]
    return {
        "rows": len(rows),
        "fit_status": dict(statuses),
        "success_rate_pct": round(statuses.get("success", 0) / len(rows) * 100, 2) if rows else 0,
        "r2_median": round(sorted(r2s)[len(r2s) // 2], 4) if r2s else None,
        "r2_ge_0_8": sum(x >= 0.8 for x in r2s),
        "r2_ge_0_9": sum(x >= 0.9 for x in r2s),
        "ec50_min": round(min(ec50s), 4) if ec50s else None,
        "ec50_median": round(sorted(ec50s)[len(ec50s) // 2], 4) if ec50s else None,
        "ec50_max": round(max(ec50s), 4) if ec50s else None,
    }


def class_summary(path: Path) -> dict:
    rows = read_rows(path)
    total = sum(int(row["Protein_Count"]) for row in rows if row.get("Protein_Count"))
    return {
        "total": total,
        "classes": len(rows),
        "top10": rows[:10],
    }


def curve_summary(path: Path) -> dict:
    rows = read_rows(path, delimiter="\t")
    p_ec50 = [as_float(row.get("pEC50", "")) for row in rows]
    p_ec50 = [x for x in p_ec50 if x is not None]
    r2 = [as_float(row.get("Curve R2", "")) for row in rows]
    r2 = [x for x in r2 if x is not None]
    reg = Counter(row.get("Curve Regulation", "") or "unlabeled" for row in rows)
    return {
        "rows": len(rows),
        "regulation_counts": dict(reg),
        "pEC50_min": round(min(p_ec50), 4) if p_ec50 else None,
        "pEC50_median": round(sorted(p_ec50)[len(p_ec50) // 2], 4) if p_ec50 else None,
        "pEC50_max": round(max(p_ec50), 4) if p_ec50 else None,
        "curve_r2_median": round(sorted(r2)[len(r2) // 2], 4) if r2 else None,
        "curve_r2_ge_0_8": sum(x >= 0.8 for x in r2),
        "curve_r2_ge_0_9": sum(x >= 0.9 for x in r2),
    }


def main() -> None:
    out = {
        "raw_rows": {
            "8170": sum(1 for _ in (TMT / "Raw_sq_8170_TMTMosaic.csv").open("r", encoding="utf-8-sig")) - 1,
            "8171": sum(1 for _ in (TMT / "Raw_sq_8171_TMTMosaic.csv").open("r", encoding="utf-8-sig")) - 1,
        },
        "site_mean_rows": {
            "8170": sum(1 for _ in (TMT / "direct_labeling_desthiobiotin_pipeline_8170" / "02_summary_tables" / "site_level_mean_8170.csv").open("r", encoding="utf-8-sig")) - 1,
            "8171": sum(1 for _ in (TMT / "direct_labeling_desthiobiotin_pipeline_8171" / "02_summary_tables" / "site_level_mean_8171.csv").open("r", encoding="utf-8-sig")) - 1,
        },
        "fit": {
            "8170": fit_summary(TMT / "direct_labeling_desthiobiotin_pipeline_8170" / "02_summary_tables" / "fit_8170.csv"),
            "8171": fit_summary(TMT / "direct_labeling_desthiobiotin_pipeline_8171" / "02_summary_tables" / "fit_8171.csv"),
        },
        "protein_classes": {
            "8170": class_summary(TMT / "TMT_protein_class_summary_8170.csv"),
            "8171": class_summary(TMT / "TMT_protein_class_summary_8171.csv"),
        },
        "curvecurator_8171": curve_summary(GRAD / "curvecurator" / "8171" / "8171_conc_curve.tsv"),
    }
    print(json.dumps(out, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
