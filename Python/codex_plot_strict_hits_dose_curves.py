from pathlib import Path
import csv
import html
import math


ROOT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha")
INPUT = ROOT / "outputs" / "strict_hits_positive_log2fc_ge1_8171" / "strict_hits_positive_log2FC_ge1_sites.csv"
OUTDIR = ROOT / "outputs" / "strict_hits_positive_log2fc_ge1_8171" / "dose_response_plot"
OUTDIR.mkdir(parents=True, exist_ok=True)


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def site_label(name):
    parts = name.split("_")
    gene = parts[0]
    site = parts[-1] if len(parts) > 2 else ""
    return f"{gene}-Cys{site}" if site else gene


with INPUT.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

for row in rows:
    row["_pEC50"] = to_float(row.get("pEC50"))
    row["_fc"] = to_float(row.get("Curve Fold Change"))
    row["_r2"] = to_float(row.get("Curve R2"))
    row["_padj"] = to_float(row.get("Curve P_Value adjusted"))

ranked = sorted(rows, key=lambda r: (r["_pEC50"], r["_fc"]), reverse=True)
highlight = ranked[:12]
highlight_names = {r["Name"] for r in highlight}

doses_m = [20e-6, 40e-6, 80e-6, 160e-6, 320e-6, 640e-6]
log_doses = [math.log10(x) for x in doses_m]
ratio_cols = [f"Ratio {i}" for i in range(1, 7)]

width, height = 1320, 760
plot_x, plot_y = 92, 104
plot_w, plot_h = 695, 470
legend_x = 830
legend_y = 128

x_min = math.log10(min(doses_m) * 0.72)
x_max = math.log10(max(doses_m) * 1.65)


def sx(x):
    return plot_x + (math.log10(x) - x_min) / (x_max - x_min) * plot_w


def sy(y):
    return plot_y + (1.05 - y) / 1.12 * plot_h


def normalized(row):
    vals = [to_float(row.get(c)) for c in ratio_cols]
    vals = [v for v in vals if math.isfinite(v)]
    if len(vals) != 6:
        return None
    lo, hi = min(vals), max(vals)
    if hi == lo:
        return None
    return [(v - lo) / (hi - lo) for v in vals]


def interp_points(yvals, steps=180):
    pts = []
    for i in range(steps):
        t = i / (steps - 1)
        xlog = log_doses[0] + t * (log_doses[-1] - log_doses[0])
        j = 0
        while j < len(log_doses) - 2 and log_doses[j + 1] < xlog:
            j += 1
        denom = log_doses[j + 1] - log_doses[j]
        local = (xlog - log_doses[j]) / denom if denom else 0
        y = yvals[j] + local * (yvals[j + 1] - yvals[j])
        pts.append((sx(10**xlog), sy(y)))
    return pts


def path_from_points(points):
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in points)


palette = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#17becf",
    "#bcbd22", "#7f7f7f", "#0f766e", "#dc2626",
]

svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
svg.append('<rect width="100%" height="100%" fill="#ffffff"/>')
svg.append('<text x="62" y="52" font-family="Arial, Helvetica, sans-serif" font-size="27" font-weight="700" fill="#10233f">Dose-response curves for positive strict-hit sites</text>')
svg.append('<text x="62" y="82" font-family="Arial, Helvetica, sans-serif" font-size="15" fill="#4b5563">Criteria: BH-adjusted p ≤ 0.05, pEC50 ≥ 4, R² ≥ 0.99, log₂(curve FC) ≥ 1</text>')

for yt in [0, 0.25, 0.5, 0.75, 1.0]:
    y = sy(yt)
    svg.append(f'<line x1="{plot_x}" y1="{y:.1f}" x2="{plot_x + plot_w}" y2="{y:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
    svg.append(f'<text x="{plot_x - 14}" y="{y + 5:.1f}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#111827">{yt:g}</text>')

xticks = [20e-6, 40e-6, 80e-6, 160e-6, 320e-6, 640e-6]
for xval in xticks:
    x = sx(xval)
    svg.append(f'<line x1="{x:.1f}" y1="{plot_y}" x2="{x:.1f}" y2="{plot_y + plot_h}" stroke="#eef2f7" stroke-width="1"/>')
    label = f"{int(xval * 1e6)}"
    svg.append(f'<text x="{x:.1f}" y="{plot_y + plot_h + 25}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#111827">{label}</text>')

svg.append(f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y + plot_h}" stroke="#4b5563" stroke-width="1.5"/>')
svg.append(f'<line x1="{plot_x}" y1="{plot_y + plot_h}" x2="{plot_x + plot_w}" y2="{plot_y + plot_h}" stroke="#4b5563" stroke-width="1.5"/>')
svg.append(f'<text x="{plot_x + plot_w / 2}" y="{plot_y + plot_h + 60}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="15" fill="#111827">WRX-035 concentration (µM, log scale)</text>')
svg.append(f'<text x="25" y="{plot_y + plot_h / 2}" transform="rotate(-90 25 {plot_y + plot_h / 2})" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="15" fill="#111827">Relative response (normalized)</text>')

for row in ranked:
    if row["Name"] in highlight_names:
        continue
    yvals = normalized(row)
    if not yvals:
        continue
    pts = interp_points(yvals, 80)
    svg.append(f'<polyline points="{path_from_points(pts)}" fill="none" stroke="#c7ced8" stroke-width="1.4" opacity="0.34"/>')

legend_rows = []
for idx, row in enumerate(highlight, start=1):
    yvals = normalized(row)
    if not yvals:
        continue
    color = palette[idx - 1]
    pts = interp_points(yvals)
    svg.append(f'<polyline points="{path_from_points(pts)}" fill="none" stroke="{color}" stroke-width="3.2" stroke-linecap="round" stroke-linejoin="round"/>')
    for dose, yv in zip(doses_m, yvals):
        svg.append(f'<circle cx="{sx(dose):.1f}" cy="{sy(yv):.1f}" r="4.4" fill="{color}" stroke="#ffffff" stroke-width="1.4"/>')
    end_x, end_y = pts[-1]
    svg.append(f'<circle cx="{end_x + 20:.1f}" cy="{end_y:.1f}" r="10" fill="{color}"/>')
    svg.append(f'<text x="{end_x + 20:.1f}" y="{end_y + 4:.1f}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="11" font-weight="700" fill="#ffffff">{idx}</text>')
    legend_rows.append((idx, site_label(row["Name"]), row["_pEC50"], row["_fc"], color, row["Name"]))

svg.append(f'<text x="{legend_x}" y="{legend_y - 32}" font-family="Arial, Helvetica, sans-serif" font-size="18" font-weight="700" fill="#10233f">Highlighted sites</text>')
svg.append(f'<text x="{legend_x}" y="{legend_y - 10}" font-family="Arial, Helvetica, sans-serif" font-size="12.5" fill="#4b5563">Top 12 ranked by pEC50; grey lines show remaining 63 sites.</text>')
for i, label, pec50, fc, color, _name in legend_rows:
    y = legend_y + (i - 1) * 42
    svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 34}" y2="{y}" stroke="{color}" stroke-width="4" stroke-linecap="round"/>')
    svg.append(f'<circle cx="{legend_x - 15}" cy="{y}" r="10" fill="{color}"/>')
    svg.append(f'<text x="{legend_x - 15}" y="{y + 4}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="10.5" font-weight="700" fill="#ffffff">{i}</text>')
    svg.append(f'<text x="{legend_x + 46}" y="{y - 4}" font-family="Arial, Helvetica, sans-serif" font-size="13.2" font-weight="700" fill="#111827">{html.escape(label)}</text>')
    svg.append(f'<text x="{legend_x + 46}" y="{y + 14}" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#4b5563">pEC50 {pec50:.2f}  |  log₂FC {fc:.2f}</text>')

svg.append('<text x="62" y="708" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#4b5563">All curves were normalized per site using Ratio 1–6 to emphasize dose-response shape rather than absolute signal intensity.</text>')
svg.append("</svg>")

svg_path = OUTDIR / "positive_strict_hits_dose_response_top12.svg"
svg_path.write_text("\n".join(svg), encoding="utf-8")

table_path = OUTDIR / "highlighted_curve_legend_table.csv"
with table_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["plot_number", "plot_label", "Name", "pEC50", "Curve Fold Change", "Curve R2", "Curve P_Value adjusted", "Curve Log P_Value adjusted"])
    for i, label, pec50, fc, _color, name in legend_rows:
        row = next(r for r in highlight if r["Name"] == name)
        writer.writerow([i, label, name, row["pEC50"], row["Curve Fold Change"], row["Curve R2"], row["Curve P_Value adjusted"], row["Curve Log P_Value adjusted"]])

print(svg_path)
print(table_path)
