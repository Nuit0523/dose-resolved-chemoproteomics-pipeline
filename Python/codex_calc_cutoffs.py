from __future__ import annotations

import csv
import json
from pathlib import Path


BASE = Path(r"E:\gradthesis\curvecurator")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def as_int(value: str) -> int:
    try:
        return int(float(value))
    except Exception:
        return 0


def summarize_8171() -> dict:
    gene_summary = read_csv(BASE / "8171" / "enrichment_from_raw_name_pEC50" / "gene_level_pEC50_summary_generated.csv")
    candidates = read_csv(BASE / "8171" / "enrichment_from_raw_name_pEC50" / "candidate_hits_after_cutoff.csv")
    n_total = len(gene_summary)
    n_sites_ge_2 = sum(as_int(r.get("n_sites", "0")) >= 2 for r in gene_summary)
    n_best_ge_4 = sum(as_float(r.get("best_pEC50", "nan")) >= 4 for r in gene_summary)
    n_both = sum(as_int(r.get("n_sites", "0")) >= 2 and as_float(r.get("best_pEC50", "nan")) >= 4 for r in gene_summary)
    return {
        "total_gene_level_rows": n_total,
        "n_sites_ge_2": n_sites_ge_2,
        "best_pEC50_ge_4": n_best_ge_4,
        "both_cutoffs": n_both,
        "candidate_file_rows": len(candidates),
        "top_candidates": candidates[:10],
    }


def summarize_topn(folder: Path, ns: list[int]) -> dict:
    out = {}
    for n in ns:
        d = folder / f"top{n}"
        out[str(n)] = {}
        for name in [f"hit_top{n}_table.csv", f"hit_top{n}_genes.csv", f"GO_BP_top{n}.csv", f"KEGG_top{n}.csv", f"Reactome_top{n}.csv"]:
            p = d / name
            if p.exists() and p.stat().st_size > 4:
                out[str(n)][name] = len(read_csv(p))
            else:
                out[str(n)][name] = 0
    return out


def main() -> None:
    out = {
        "8171_candidate_cutoff_funnel": summarize_8171(),
        "8171_topn_outputs": summarize_topn(BASE / "8171" / "enrichment_from_raw_name_pEC50", [100, 200, 400]),
        "8170_topn_outputs": summarize_topn(BASE / "8170" / "enrichment_compare_topN", [30, 50, 80]),
    }
    print(json.dumps(out, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
