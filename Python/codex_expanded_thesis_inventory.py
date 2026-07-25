from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


GRAD = Path(r"E:\gradthesis")
CURVE = GRAD / "curvecurator"
TMT = Path(r"E:\R\TMT_analysis\TMT4")


def rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def head_csv(path: Path, n: int = 10, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return [row for _, row in zip(range(n), reader)]


def text(path: Path, max_chars: int = 4000) -> str:
    return path.read_text(encoding="utf-8-sig", errors="replace")[:max_chars]


def compact_tree(root: Path) -> dict:
    files = [p for p in root.rglob("*") if p.is_file()]
    ext_counts = Counter(p.suffix.lower() or "[no_ext]" for p in files)
    top_level = Counter()
    for p in files:
        rel = p.relative_to(root)
        top_level[rel.parts[0] if rel.parts else "."] += 1
    return {
        "root": str(root),
        "file_count": len(files),
        "extension_counts": dict(ext_counts.most_common(20)),
        "top_level_file_counts": dict(top_level.most_common(30)),
    }


def gsea_summary(base: Path) -> dict:
    out: dict[str, object] = {}
    comparison = base / "GSEA_comparison_summary.csv"
    params = base / "GSEA_parameters.csv"
    if comparison.exists():
        out["comparison"] = head_csv(comparison, 20)
    if params.exists():
        out["parameters"] = head_csv(params, 20)
    for mode in ["best_pEC50", "mean_pEC50", "weighted_score"]:
        mode_dir = base / mode
        mode_out = {}
        for csv_file in sorted(mode_dir.glob("GSEA_*.csv")):
            if csv_file.stat().st_size <= 4:
                mode_out[csv_file.name] = {"rows": 0, "top": []}
                continue
            data = head_csv(csv_file, 5)
            mode_out[csv_file.name] = {"rows_at_least": len(rows(csv_file)), "top": data}
        summary_files = sorted(mode_dir.glob("summary_*.txt"))
        if summary_files:
            mode_out["summary_text_preview"] = text(summary_files[0], 1200)
        out[mode] = mode_out
    return out


def enrichment_summary(base: Path) -> dict:
    out: dict[str, object] = {}
    for name in ["pipeline_parameters.csv", "enrichment_comparison_summary.csv", "candidate_hits_after_cutoff.csv", "gene_level_pEC50_summary_generated.csv"]:
        p = base / name
        if p.exists():
            data = head_csv(p, 10)
            out[name] = {"rows_at_least": len(rows(p)), "top": data}
    for cutoff in ["top100", "top200", "top400"]:
        cutoff_dir = base / cutoff
        cutoff_out = {}
        for p in sorted(cutoff_dir.glob("*")):
            if p.suffix.lower() == ".csv":
                cutoff_out[p.name] = {"rows_at_least": len(rows(p)) if p.stat().st_size > 4 else 0, "top": head_csv(p, 5) if p.stat().st_size > 4 else []}
            elif p.suffix.lower() == ".txt":
                cutoff_out[p.name] = text(p, 800)
        out[cutoff] = cutoff_out
    return out


def alphafold_summary(root: Path) -> dict:
    doc = root / "doc"
    files = [p for p in doc.rglob("*") if p.is_file()] if doc.exists() else []
    sample_dirs = sorted({p.parent for p in files if p.name.endswith("_data.json")})
    targets = []
    for d in sample_dirs[:30]:
        target_files = list(d.glob("*"))
        summary = next((p for p in target_files if p.name.endswith("_summary_confidences.json")), None)
        ranking = next((p for p in target_files if p.name.endswith("_ranking_scores.csv")), None)
        model = next((p for p in target_files if p.name.endswith("_model.cif")), None)
        entry = {
            "target_dir": str(d),
            "target": d.name,
            "has_model_cif": model is not None,
            "has_summary_confidences": summary is not None,
            "has_ranking_scores": ranking is not None,
        }
        if summary:
            try:
                entry["summary_confidences"] = json.loads(summary.read_text(encoding="utf-8", errors="replace"))
            except Exception as exc:
                entry["summary_confidences_error"] = str(exc)
        if ranking:
            entry["ranking_scores"] = head_csv(ranking, 5)
        targets.append(entry)
    ext_counts = Counter(p.suffix.lower() or "[no_ext]" for p in files)
    protein_dirs = sorted({p.parent.name for p in files if p.name.endswith("_model.cif")})
    return {
        "doc_exists": doc.exists(),
        "file_count": len(files),
        "extension_counts": dict(ext_counts.most_common(20)),
        "protein_model_dir_count": len(protein_dirs),
        "protein_model_dirs_sample": protein_dirs[:50],
        "targets_sample": targets,
    }


def tmt_development_summary() -> dict:
    summaries = {}
    for p in sorted(TMT.glob("Round*_summary.csv")):
        summaries[p.name] = head_csv(p, 20)
    round_txt = {}
    for p in sorted(TMT.glob("TMT_protein_classification_round*_summary.txt")):
        round_txt[p.name] = text(p, 1600)
    return {
        "8170_files": [str(p) for p in (TMT / "direct_labeling_desthiobiotin_pipeline_8170").rglob("*") if p.is_file()],
        "8171_files": [str(p) for p in (TMT / "direct_labeling_desthiobiotin_pipeline_8171").rglob("*") if p.is_file()],
        "round_csv_summaries": summaries,
        "round_txt_summaries_preview": round_txt,
        "final_pipeline_outputs": [str(p) for p in (TMT / "direct_labeling_desthiobiotin_pipeline_output").rglob("*") if p.is_file()],
    }


def main() -> None:
    out = {
        "curvecurator_tree": compact_tree(CURVE),
        "curvecurator_8171_tree": compact_tree(CURVE / "8171"),
        "curvecurator_config_8171": text(CURVE / "8171" / "8171_conc_config.toml", 4000),
        "normalization_factors_8171": text(CURVE / "8171" / "8171_normalization_factors.tsv", 1000),
        "mad_analysis_8171": text(CURVE / "8171" / "8171_mad_analysis.tsv", 1000),
        "driver_genes_top_go": head_csv(CURVE / "8171" / "driver_genes_top_GO_term.csv", 20),
        "enrichment": enrichment_summary(CURVE / "8171" / "enrichment_from_raw_name_pEC50"),
        "gsea": gsea_summary(CURVE / "8171" / "GSEA_from_raw_name_pEC50"),
        "alphafold": alphafold_summary(CURVE),
        "tmt_development": tmt_development_summary(),
    }
    print(json.dumps(out, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
