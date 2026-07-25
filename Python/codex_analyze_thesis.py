from __future__ import annotations

import csv
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET


ROOTS = {
    "gradthesis": Path(r"E:\gradthesis"),
    "tmt4": Path(r"E:\R\TMT_analysis\TMT4"),
}


def csv_profile(path: Path, max_rows: int | None = None) -> dict:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        rows = 0
        nulls = Counter()
        samples = []
        for row in reader:
            rows += 1
            if len(samples) < 3:
                samples.append(row)
            for key, val in row.items():
                if val in ("", "NA", "NaN", "nan", None):
                    nulls[key] += 1
            if max_rows and rows >= max_rows:
                break
    return {
        "path": str(path),
        "rows_scanned": rows,
        "columns": headers,
        "sample_rows": samples,
        "top_nulls": nulls.most_common(10),
    }


def docx_text(path: Path) -> str:
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    parts: list[str] = []
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    for para in root.findall(".//w:p", ns):
        texts = [node.text or "" for node in para.findall(".//w:t", ns)]
        line = "".join(texts).strip()
        if line:
            parts.append(line)
    return "\n".join(parts)


def summarize_docx(path: Path) -> dict:
    text = docx_text(path)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    heading_like = [
        line
        for line in lines
        if re.match(r"^(\d+(\.\d+)*\.?\s+|Abstract|Introduction|Methods|Results|Discussion|Conclusion|References)", line, re.I)
    ]
    return {
        "path": str(path),
        "characters": len(text),
        "lines": len(lines),
        "heading_like": heading_like[:80],
        "first_3000_chars": text[:3000],
        "last_2000_chars": text[-2000:],
    }


def read_csv_table(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    profiles = {
        "raw_8170": csv_profile(ROOTS["tmt4"] / "Raw_sq_8170_TMTMosaic.csv", max_rows=None),
        "raw_8171": csv_profile(ROOTS["tmt4"] / "Raw_sq_8171_TMTMosaic.csv", max_rows=None),
        "fit_8170": csv_profile(ROOTS["tmt4"] / "direct_labeling_desthiobiotin_pipeline_8170" / "02_summary_tables" / "fit_8170.csv", max_rows=None),
        "fit_8171": csv_profile(ROOTS["tmt4"] / "direct_labeling_desthiobiotin_pipeline_8171" / "02_summary_tables" / "fit_8171.csv", max_rows=None),
        "site_mean_8170": csv_profile(ROOTS["tmt4"] / "direct_labeling_desthiobiotin_pipeline_8170" / "02_summary_tables" / "site_level_mean_8170.csv", max_rows=None),
        "site_mean_8171": csv_profile(ROOTS["tmt4"] / "direct_labeling_desthiobiotin_pipeline_8171" / "02_summary_tables" / "site_level_mean_8171.csv", max_rows=None),
    }

    class_8170 = read_csv_table(ROOTS["tmt4"] / "TMT_protein_class_summary_8170.csv")
    class_8171 = read_csv_table(ROOTS["tmt4"] / "TMT_protein_class_summary_8171.csv")

    ec50_8171 = csv_profile(ROOTS["gradthesis"] / "curvecurator" / "8171" / "8171_EC50.csv", max_rows=None)
    conc_curve_8171 = csv_profile(ROOTS["gradthesis"] / "curvecurator" / "8171" / "8171_conc_curve.tsv", max_rows=None)

    out = {
        "profiles": profiles,
        "class_summary_8170_top20": class_8170[:20],
        "class_summary_8171_top20": class_8171[:20],
        "ec50_8171": ec50_8171,
        "conc_curve_8171": conc_curve_8171,
        "thesis_docx": summarize_docx(Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\work\thesis.docx")),
    }

    print(json.dumps(out, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
