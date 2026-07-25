from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path


GRAD = Path(r"E:\gradthesis")
CURVE = GRAD / "curvecurator"
TMT = Path(r"E:\R\TMT_analysis\TMT4")
OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\thesis_defense_expanded_brief.md")


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def file_count(root: Path, pattern: str = "*") -> int:
    return sum(1 for p in root.rglob(pattern) if p.is_file())


def ext_counts(root: Path) -> Counter:
    c = Counter()
    for p in root.rglob("*"):
        if p.is_file():
            c[p.suffix.lower() or "[no_ext]"] += 1
    return c


def parse_residuals(text: str) -> tuple[int | None, int | None, int | None]:
    before = re.search(r"Residual(?:-unclassified)? before(?: Round \d+)?:\s*(\d+)", text, re.I)
    new = re.search(r"Newly classified(?: in round \d+| in Round \d+)?:\s*(\d+)", text, re.I)
    after = re.search(r"Residual(?:-unclassified)? after(?: Round \d+)?:\s*(\d+)", text, re.I)
    return (
        int(before.group(1)) if before else None,
        int(new.group(1)) if new else None,
        int(after.group(1)) if after else None,
    )


def tmt_round_table() -> list[tuple[str, int | None, int | None, int | None]]:
    rows = []
    for p in sorted(TMT.glob("TMT_protein_classification_round*_summary.txt")):
        text = p.read_text(encoding="utf-8-sig", errors="replace")
        m = re.search(r"round(\d+)", p.name, re.I)
        rows.append((m.group(1) if m else p.name, *parse_residuals(text)))
    return rows


def curve_counts() -> dict:
    curve_rows = read_csv(CURVE / "8171" / "8171_conc_curve.tsv", delimiter="\t")
    regs = Counter(row.get("Curve Regulation", "") or "unlabeled" for row in curve_rows)
    ec50_rows = read_csv(CURVE / "8171" / "8171_EC50.csv")
    candidates = read_csv(CURVE / "8171" / "enrichment_from_raw_name_pEC50" / "candidate_hits_after_cutoff.csv")
    gene_summary = read_csv(CURVE / "8171" / "enrichment_from_raw_name_pEC50" / "gene_level_pEC50_summary_generated.csv")
    return {
        "curves": len(curve_rows),
        "regulation": regs,
        "ec50_rows": len(ec50_rows),
        "gene_summary_rows": len(gene_summary),
        "candidate_rows": len(candidates),
        "candidate_top10": candidates[:10],
    }


def enrichment_counts() -> dict:
    comp = read_csv(CURVE / "8171" / "enrichment_from_raw_name_pEC50" / "enrichment_comparison_summary.csv")
    params = read_csv(CURVE / "8171" / "enrichment_from_raw_name_pEC50" / "pipeline_parameters.csv")
    return {"comparison": comp, "parameters": params}


def gsea_counts() -> dict:
    base = CURVE / "8171" / "GSEA_from_raw_name_pEC50"
    out = {"comparison": read_csv(base / "GSEA_comparison_summary.csv"), "sets": {}}
    for mode in ["best_pEC50", "mean_pEC50", "weighted_score"]:
        mode_dir = base / mode
        mode_data = {}
        for p in sorted(mode_dir.glob("GSEA_*.csv")):
            if p.stat().st_size <= 4:
                mode_data[p.name] = {"count": 0, "top": []}
            else:
                data = read_csv(p)
                mode_data[p.name] = {"count": len(data), "top": data[:5]}
        out["sets"][mode] = mode_data
    return out


def alphafold_counts() -> dict:
    doc = CURVE / "doc"
    exts = ext_counts(doc)
    model_files = sorted(doc.rglob("*_model.cif"))
    target_dirs = sorted({p.parent for p in doc.rglob("*_data.json")})
    protein_names = sorted({p.parent.name for p in model_files})
    ranking_files = sorted(doc.rglob("*_ranking_scores.csv"))
    ranking_examples = []
    for p in ranking_files[:12]:
        ranking_examples.append({"target": p.parent.name, "rows": read_csv(p)[:5]})
    return {
        "files": file_count(doc),
        "ext_counts": exts,
        "model_files": len(model_files),
        "target_dirs": len(target_dirs),
        "protein_names_sample": protein_names[:40],
        "ranking_files": len(ranking_files),
        "ranking_examples": ranking_examples,
    }


def write_brief() -> None:
    curve = curve_counts()
    enrich = enrichment_counts()
    gsea = gsea_counts()
    af = alphafold_counts()
    rounds = tmt_round_table()
    curve_ext = ext_counts(CURVE)
    c8171_ext = ext_counts(CURVE / "8171")

    lines: list[str] = []
    lines.append("# Expanded Thesis Defense Brief\n")
    lines.append("## Corrected Scope\n")
    lines.append("The thesis is not only a 4PL fitting project. It has four connected layers: code/pipeline development on 8170, production-scale application to 8171, CurveCurator-based dose-response and enrichment analysis, and AlphaFold3/docking structural interpretation.\n")

    lines.append("## 1. Pipeline Development: 8170 as the Test Dataset\n")
    lines.append("- 8170 was used as the smaller development/test dataset before scaling the workflow to 8171.")
    lines.append("- 8170 outputs include raw TMTMosaic input, long-format site table, site-level mean table, fit table, and a final classification table.")
    lines.append("- 8170 raw TMTMosaic rows: 8,707; 8170 site-level mean rows: 35,304; 8170 fit rows: 5,884.")
    lines.append("- The 8170 pipeline folder contains `01_long_and_metadata`, `02_summary_tables`, and `TMT_protein_classification_round13_8170.csv`, which makes it useful for explaining how the code was validated before the larger 8171 run.\n")

    lines.append("## 2. Production Application: 8171\n")
    lines.append("- 8171 is the larger production dataset: 15,331 raw TMTMosaic rows and 81,738 site-level mean rows.")
    lines.append("- 8171 classification proceeded through later refinement rounds, ending at round 22.")
    lines.append("- Key residual-classification milestones:")
    lines.append("| Round | Residual before | Newly classified | Residual after |")
    lines.append("|---:|---:|---:|---:|")
    for round_id, before, new, after in rounds:
        if before is not None or new is not None or after is not None:
            lines.append(f"| {round_id} | {before if before is not None else ''} | {new if new is not None else ''} | {after if after is not None else ''} |")
    lines.append("- Round 22 reduced residual unclassified genes from 12 to 0, giving a fully assigned final 8171 classification table.\n")

    lines.append("## 3. CurveCurator Layer\n")
    lines.append(f"- `E:/gradthesis/curvecurator` contains {file_count(CURVE):,} files; the 8171 CurveCurator folder contains {file_count(CURVE / '8171'):,} files.")
    lines.append(f"- Major 8171 output types: {dict(c8171_ext.most_common(8))}.")
    lines.append("- 8171 CurveCurator config: 6 dose points, 20 to 640 uM, TMT measurement, control experiment 1, imputation enabled, normalization enabled, FDR Benjamini-Hochberg correction.")
    lines.append("- Normalization factors were generated for Raw 1 through Raw 6, and MAD analysis was recorded for Ratio 1 through Ratio 6.")
    lines.append(f"- Curve table rows: {curve['curves']:,}; EC50 rows: {curve['ec50_rows']:,}; gene-level pEC50 summaries: {curve['gene_summary_rows']:,}; candidate genes after cutoff: {curve['candidate_rows']:,}.")
    lines.append(f"- Curve regulation counts: {dict(curve['regulation'])}.")
    lines.append("- Top candidate genes by best pEC50 after cutoff include: " + ", ".join(row["gene_symbol"] for row in curve["candidate_top10"]) + ".\n")

    lines.append("## 4. Enrichment and GSEA Layer\n")
    lines.append("- Simple over-representation enrichment on top 100/200/400 pEC50-ranked genes found no significant GO, KEGG, or Reactome terms under the chosen cutoffs; this is useful as a negative-control style result.")
    lines.append("- Enrichment parameters included min_sites_cutoff = 2, min_pEC50_cutoff = 4, and top_n_list = 100/200/400.")
    lines.append("- GSEA was then run with ranked gene lists rather than only discrete top-N cutoffs.")
    lines.append("- GSEA comparison summary:")
    for row in gsea["comparison"]:
        lines.append(f"  - {row}")
    lines.append("- Recurrent driver genes for the top GO term include nuclear transport/RNA trafficking-related genes such as NUP205, TPR, RANBP2, XPO1, NUP155, and KPNB1.")
    lines.append("- This supports the thesis-level interpretation that nuclear transport and RNA trafficking emerged from ranked dose-response behavior, not just from a manually curated target list.\n")

    lines.append("## 5. AlphaFold3 and Structural Interpretation\n")
    lines.append(f"- The structural `doc` folder contains {af['files']:,} files, including {af['model_files']:,} model CIF files and {af['ranking_files']:,} ranking-score CSV files.")
    lines.append(f"- Structural file types include: {dict(af['ext_counts'].most_common(8))}.")
    lines.append("- Example protein/model folders include: " + ", ".join(af["protein_names_sample"][:20]) + ".")
    lines.append("- This layer provides structural context for selected probe-responsive cysteine sites: AlphaFold3 models, confidence files, ranking scores, and downstream docking/visualization assets can be used to argue whether a cysteine is accessible or located near a plausible binding pocket.")
    lines.append("- For the defense, this should be framed as mechanistic prioritization rather than definitive biochemical validation.\n")

    lines.append("## 6. Defense Narrative Update\n")
    lines.append("A stronger defense storyline is:\n")
    lines.append("1. Develop and debug a reproducible TMT direct-labeling analysis pipeline on 8170.")
    lines.append("2. Apply the validated workflow to the larger 8171 dataset.")
    lines.append("3. Use CurveCurator to transform concentration-resolved TMT signals into interpretable curve metrics.")
    lines.append("4. Move from site-level curves to gene-level candidate prioritization and pathway interpretation.")
    lines.append("5. Use AlphaFold3 and docking to add structural hypotheses for selected candidates.")
    lines.append("6. Conclude that the thesis contributes an integrated quantitative-to-structural workflow for prioritizing covalent-probe responsive targets.\n")

    lines.append("## 7. Slide Implications\n")
    lines.append("- Add one methods slide specifically titled `8170 was used to validate the pipeline before scaling to 8171`.")
    lines.append("- Add one CurveCurator slide showing config choices, curve outputs, dashboard/curve TSV, and regulation counts.")
    lines.append("- Add one enrichment/GSEA slide contrasting top-N ORA with ranked GSEA.")
    lines.append("- Add one AlphaFold3/docking slide focused on target prioritization and structural plausibility.")
    lines.append("- Keep 4PL as one part of the quantitative modeling section, not the whole thesis.\n")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(str(OUT))


if __name__ == "__main__":
    write_brief()
