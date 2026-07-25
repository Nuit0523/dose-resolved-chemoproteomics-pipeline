from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha")
INPUT = ROOT / "outputs" / "strict_hits_positive_log2fc_ge1_8171" / "strict_hits_positive_log2FC_ge1_sites.csv"
OUTDIR = ROOT / "outputs" / "strict_hits_positive_log2fc_ge1_8171" / "dose_response_plot_python"
OUTDIR.mkdir(parents=True, exist_ok=True)


def site_label(name: str) -> str:
    parts = str(name).split("_")
    gene = parts[0]
    site = parts[-1] if len(parts) > 2 else ""
    return f"{gene}-Cys{site}" if site else gene


df = pd.read_csv(INPUT)
numeric_cols = [
    "pEC50",
    "Curve Slope",
    "Curve Front",
    "Curve Back",
    "Curve Fold Change",
    "Curve R2",
    "Curve P_Value adjusted",
    "Curve Log P_Value adjusted",
]
for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

doses_uM = np.array([20, 40, 80, 160, 320, 640], dtype=float)
doses_M = doses_uM * 1e-6
pconc = -np.log10(doses_M)
pconc_smooth = np.linspace(pconc.max() + 0.07, pconc.min() - 0.07, 320)
doses_smooth_M = 10 ** (-pconc_smooth)

ranked = (
    df.dropna(subset=["pEC50", "Curve Fold Change", "Curve R2"])
    .sort_values(["pEC50", "Curve Fold Change"], ascending=[False, False])
    .reset_index(drop=True)
)
highlight = ranked.head(12).copy()
highlight_names = set(highlight["Name"])


def normalized_response(row: pd.Series):
    y = four_pl(row, pconc_smooth)
    if y is None or np.isnan(y).any():
        return None
    y_min = np.min(y)
    y_max = np.max(y)
    if not np.isfinite(y_min) or not np.isfinite(y_max) or y_max == y_min:
        return None
    return (y - y_min) / (y_max - y_min)


def four_pl(row: pd.Series, p_conc_values):
    """CurveCurator-style 4PL using p-concentration: low dose -> front, high dose -> back."""
    pec50 = row["pEC50"]
    slope = row["Curve Slope"]
    front = row["Curve Front"]
    back = row["Curve Back"]
    if not np.all(np.isfinite([pec50, slope, front, back])):
        return None
    return front + (back - front) / (1 + np.power(10, slope * (p_conc_values - pec50)))


plt.rcParams.update(
    {
        "font.family": "Arial",
        "axes.titlesize": 18,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 8.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

fig, ax = plt.subplots(figsize=(12.8, 6.9), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

for _, row in ranked.iterrows():
    if row["Name"] in highlight_names:
        continue
    y_norm = normalized_response(row)
    if y_norm is None:
        continue
    ax.plot(
        doses_smooth_M,
        y_norm,
        color="#c8d0da",
        lw=1.0,
        alpha=0.30,
        zorder=1,
    )

palette = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#17becf",
    "#bcbd22",
    "#7f7f7f",
    "#0f766e",
    "#dc2626",
]

legend_handles = []
for idx, (_, row) in enumerate(highlight.iterrows(), start=1):
    y_norm = normalized_response(row)
    if y_norm is None:
        continue
    color = palette[idx - 1]
    ax.plot(
        doses_smooth_M,
        y_norm,
        color=color,
        lw=2.55,
        zorder=4,
    )
    label = site_label(row["Name"])
    ax.annotate(
        str(idx),
        xy=(doses_smooth_M[-1], y_norm[-1]),
        xytext=(8, 0),
        textcoords="offset points",
        ha="center",
        va="center",
        fontsize=7.5,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="circle,pad=0.24", fc=color, ec="none"),
        zorder=6,
    )
    legend_handles.append(
        plt.Line2D(
            [0],
            [0],
            color=color,
            lw=2.35,
            label=f"{idx}. {label}   pEC50={row['pEC50']:.2f}   log2FC={row['Curve Fold Change']:.2f}",
        )
    )

ax.set_xscale("log")
ax.set_xlim(doses_M.min() * 0.72, doses_M.max() * 1.9)
ax.set_ylim(-0.05, 1.09)

ax.set_title("Dose-response curves for positive strict-hit sites", loc="left", fontweight="bold", pad=16)
ax.text(
    0.0,
    1.02,
    "Criteria: BH-adjusted p <= 0.05, pEC50 >= 4, R2 >= 0.99, log2(curve FC) >= 1",
    transform=ax.transAxes,
    ha="left",
    va="bottom",
    fontsize=10,
    color="#4b5563",
)

ax.set_xlabel("WRX-035 concentration (M)")
ax.set_ylabel("Relative response (normalized)")
ax.set_xticks(doses_M)
ax.set_xticklabels(["2e-5", "4e-5", "8e-5", "1.6e-4", "3.2e-4", "6.4e-4"])
ax.set_yticks(np.linspace(0, 1, 5))
ax.grid(True, which="major", color="#e5e7eb", lw=0.75)
ax.grid(True, which="minor", axis="x", color="#eef2f7", lw=0.45)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#6b7280")
ax.spines["bottom"].set_color("#6b7280")

leg = ax.legend(
    handles=legend_handles,
    title="Highlighted sites",
    title_fontsize=9.5,
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=False,
    borderaxespad=0,
    handlelength=2.6,
)
for text in leg.get_texts():
    text.set_color("#111827")
leg.get_title().set_color("#111827")

ax.text(
    0.0,
    -0.2,
    "Grey lines: remaining positive strict-hit sites (n = 63). Colored lines: top 12 sites ranked by pEC50. "
    "Curves are fitted from CurveCurator 4PL parameters and normalized per site.",
    transform=ax.transAxes,
    fontsize=9.2,
    color="#4b5563",
    ha="left",
    va="top",
)

fig.subplots_adjust(left=0.08, right=0.62, top=0.84, bottom=0.23)

png_path = OUTDIR / "positive_strict_hits_dose_response_top12_python.png"
pdf_path = OUTDIR / "positive_strict_hits_dose_response_top12_python.pdf"
svg_path = OUTDIR / "positive_strict_hits_dose_response_top12_python.svg"
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
fig.savefig(svg_path, bbox_inches="tight", facecolor="white")

legend_table = highlight.copy()
legend_table.insert(0, "plot_number", range(1, len(legend_table) + 1))
legend_table["plot_label"] = legend_table["Name"].map(site_label)
legend_table[
    [
        "plot_number",
        "plot_label",
        "Name",
        "pEC50",
        "Curve Fold Change",
        "Curve R2",
        "Curve P_Value adjusted",
        "Curve Log P_Value adjusted",
    ]
].to_csv(OUTDIR / "highlighted_curve_legend_table_python.csv", index=False)

print(png_path)
print(pdf_path)
print(svg_path)
print(OUTDIR / "highlighted_curve_legend_table_python.csv")
