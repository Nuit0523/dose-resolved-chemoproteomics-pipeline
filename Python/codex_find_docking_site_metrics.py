from pathlib import Path
import csv

targets = {
    ("DHX57", "1369"),
    ("PSME1", "106"),
    ("IPO5", "682"),
    ("IPO7", "477"),
}

base = Path(r"E:\gradthesis\curvecurator\8171")
curve_file = base / "8171_conc_curve.tsv"
site_file = base / "enrichment_from_raw_name_pEC50" / "site_level_with_gene_symbol.csv"
out_file = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\work\docking_site_metrics.csv")

site_symbol = {}
with site_file.open(newline="", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        site_symbol[row["Name"]] = row.get("gene_symbol", "")

rows = []
with curve_file.open(newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name = row["Name"]
        gene = site_symbol.get(name, name.split("_", 1)[0])
        parts = name.split("_")
        cys = parts[-1] if parts else ""
        if (gene, cys) in targets:
            rows.append({
                "Gene": gene,
                "Cys": f"Cys{cys}",
                "Name": name,
                "pEC50": row.get("pEC50", ""),
                "Curve Fold Change": row.get("Curve Fold Change", ""),
                "Curve R2": row.get("Curve R2", ""),
                "Curve RMSE": row.get("Curve RMSE", ""),
                "Curve P_Value adjusted": row.get("Curve P_Value adjusted", ""),
                "Curve Log P_Value adjusted": row.get("Curve Log P_Value adjusted", ""),
                "Curve Regulation": row.get("Curve Regulation", ""),
            })

rows.sort(key=lambda r: ["DHX57", "PSME1", "IPO5", "IPO7"].index(r["Gene"]))

with out_file.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [
        "Gene", "Cys", "Name", "pEC50", "Curve Fold Change", "Curve R2", "Curve RMSE",
        "Curve P_Value adjusted", "Curve Log P_Value adjusted", "Curve Regulation"
    ])
    writer.writeheader()
    writer.writerows(rows)

print(f"matched {len(rows)} rows")
for row in rows:
    print(row)
