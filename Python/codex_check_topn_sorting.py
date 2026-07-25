from pathlib import Path
import csv

base = Path(r"E:\gradthesis\curvecurator\8171\enrichment_from_raw_name_pEC50")
candidate = base / "candidate_hits_after_cutoff.csv"
top100 = base / "top100" / "hit_top100_table.csv"
top200 = base / "top200" / "hit_top200_table.csv"
top400 = base / "top400" / "hit_top400_table.csv"

def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

cand = read_csv(candidate)
tops = {100: read_csv(top100), 200: read_csv(top200), 400: read_csv(top400)}

def is_sorted_desc(rows, col):
    vals = [float(r[col]) for r in rows]
    return all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))

print("candidate rows", len(cand))
print("columns", cand[0].keys())
for col in ["best_pEC50", "mean_pEC50", "median_pEC50", "n_sites"]:
    if col != "n_sites":
        print(f"candidate sorted desc by {col}:", is_sorted_desc(cand, col))

for n, rows in tops.items():
    same = [r["gene_symbol"] for r in rows] == [r["gene_symbol"] for r in cand[:n]]
    print(f"top{n} rows", len(rows), "matches first candidate rows:", same)
    print(f"top{n} first/last best_pEC50", rows[0]["best_pEC50"], rows[-1]["best_pEC50"])
