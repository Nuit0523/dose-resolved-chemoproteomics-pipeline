from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\qc_parameter_images")
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10233F"
MUTED = "#5D6675"
LINE = "#C8D0DA"
FAINT = "#EEF2F6"
BLUE = "#2E6E9E"
TEAL = "#238C82"
ORANGE = "#D36B2A"
GREEN = "#238C82"
WHITE = "#FFFFFF"

TOTAL = 13615
W, H = 1600, 560


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


F_TITLE = font(36, True)
F_SUB = font(22)
F_LABEL = font(21)
F_LABEL_B = font(22, True)
F_COUNT = font(21)
F_FOOT = font(18)


def text(draw, xy, value, fill, fnt, anchor=None):
    draw.text(xy, value, fill=fill, font=fnt, anchor=anchor)


def draw_binned_chart(filename, title, subtitle, bins, median_label=None, threshold_label=None):
    img = Image.new("RGB", (W, H), WHITE)
    d = ImageDraw.Draw(img)

    x0 = 70
    text(d, (x0, 55), title, NAVY, F_TITLE)
    text(d, (x0, 108), subtitle, MUTED, F_SUB)
    if median_label:
        text(d, (W - 70, 55), median_label, NAVY, F_LABEL_B, anchor="ra")
    if threshold_label:
        text(d, (W - 70, 108), threshold_label, MUTED, F_SUB, anchor="ra")

    max_count = max(item["count"] for item in bins)
    y0 = 245
    row_h = 78 if len(bins) <= 3 else 64
    label_x = 70
    bar_x = 520
    bar_w = 690
    count_x = 1265
    pct_x = 1485

    for idx, item in enumerate(bins):
        y = y0 + idx * row_h
        count = item["count"]
        pct = count / TOTAL * 100
        text(d, (label_x, y + 10), item["label"], MUTED, F_LABEL)
        d.rounded_rectangle([bar_x, y + 5, bar_x + bar_w, y + 29], radius=2, fill=FAINT)
        fill_w = max(18, int(bar_w * count / max_count))
        d.rounded_rectangle([bar_x, y + 5, bar_x + fill_w, y + 29], radius=2, fill=item["color"])
        text(d, (count_x, y + 7), f"{count:,}", NAVY, F_COUNT)
        text(d, (pct_x, y + 7), f"{pct:.1f}%", MUTED, F_COUNT, anchor="ra")

    d.line([70, H - 65, W - 70, H - 65], fill=LINE, width=2)
    text(d, (70, H - 38), "8171 CurveCurator output; binned from confirmed summary counts.", MUTED, F_FOOT)

    img.save(OUT / filename)


def make_panel():
    files = [
        "qc_pec50_distribution.png",
        "qc_r2_distribution.png",
        "qc_rmse_distribution.png",
        "qc_adjlogp_potency_screen.png",
    ]
    imgs = [Image.open(OUT / f).convert("RGB") for f in files]
    thumb_w = 1020
    thumb_h = int(H * thumb_w / W)
    thumbs = [img.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS) for img in imgs]
    gap = 42
    panel = Image.new("RGB", (thumb_w * 2 + gap * 3, thumb_h * 2 + gap * 3), WHITE)
    positions = [(gap, gap), (thumb_w + gap * 2, gap), (gap, thumb_h + gap * 2), (thumb_w + gap * 2, thumb_h + gap * 2)]
    for im, pos in zip(thumbs, positions):
        panel.paste(im, pos)
    panel.save(OUT / "qc_four_parameter_panel.png")


def main():
    draw_binned_chart(
        "qc_pec50_distribution.png",
        "pEC50 potency screen",
        "Higher pEC50 indicates stronger dose responsiveness.",
        [
            {"label": "pEC50 < 4", "count": TOTAL - 5423, "color": BLUE},
            {"label": "pEC50 >= 4", "count": 5423, "color": ORANGE},
        ],
        median_label="median = 3.914",
        threshold_label="cutoff: pEC50 >= 4",
    )

    draw_binned_chart(
        "qc_r2_distribution.png",
        "R2 distribution (binned)",
        "R2 describes overall curve fit quality.",
        [
            {"label": "R2 < 0.5", "count": 1911, "color": ORANGE},
            {"label": "0.5 <= R2 < 0.8", "count": 2713, "color": BLUE},
            {"label": "R2 >= 0.8", "count": 8991, "color": TEAL},
        ],
        median_label="median = 0.919",
        threshold_label="QC concern: R2 < 0.8",
    )

    draw_binned_chart(
        "qc_rmse_distribution.png",
        "RMSE error screen",
        "RMSE measures point-wise deviation from the fitted curve.",
        [
            {"label": "RMSE <= 0.2", "count": TOTAL - 1159, "color": TEAL},
            {"label": "RMSE > 0.2", "count": 1159, "color": ORANGE},
        ],
        median_label="median = 0.068",
        threshold_label="QC concern: RMSE > 0.2",
    )

    draw_binned_chart(
        "qc_adjlogp_potency_screen.png",
        "Adjusted log-p + potency screen",
        "Curves passing both statistical significance and potency screens.",
        [
            {"label": "pEC50 < 4", "count": TOTAL - 5423, "color": BLUE},
            {"label": "pEC50 >= 4 only", "count": 5423 - 3083, "color": ORANGE},
            {"label": "adj log-p >= 1.301 and pEC50 >= 4", "count": 3083, "color": GREEN},
        ],
        median_label="median adj log-p = 1.611",
        threshold_label="1.301 = adjusted p <= 0.05",
    )

    make_panel()


if __name__ == "__main__":
    main()
