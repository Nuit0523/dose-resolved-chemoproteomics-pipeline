from pathlib import Path
from math import exp

from PIL import Image, ImageDraw, ImageFont


OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\qc_metric_icons")
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10233F"
MUTED = "#5D6675"
LINE = "#C8D0DA"
FAINT = "#F3F6F9"
BLUE = "#2E6E9E"
TEAL = "#238C82"
ORANGE = "#D36B2A"
RED = "#B94A48"
WHITE = "#FFFFFF"


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


F_TITLE = font(42, True)
F_SMALL = font(20)


def new_icon(size=512):
    img = Image.new("RGBA", (size, size), (255, 255, 255, 0))
    return img, ImageDraw.Draw(img)


def rounded_card(d, size=512):
    pad = 20
    d.rounded_rectangle([pad, pad, size - pad, size - pad], radius=54, fill=WHITE, outline=LINE, width=6)


def draw_axes(d):
    d.line([104, 380, 418, 380], fill=MUTED, width=7)
    d.line([104, 380, 104, 104], fill=MUTED, width=7)


def draw_polyline(d, pts, fill, width=10):
    if len(pts) > 1:
        d.line(pts, fill=fill, width=width, joint="curve")


def r2_icon():
    img, d = new_icon()
    rounded_card(d)
    draw_axes(d)
    curve = [(112, 352), (160, 330), (210, 286), (262, 226), (318, 165), (388, 118)]
    draw_polyline(d, curve, TEAL, 12)
    for x, y in [(122, 350), (170, 318), (218, 278), (268, 228), (322, 168), (382, 124)]:
        d.ellipse([x - 11, y - 11, x + 11, y + 11], fill=NAVY)
    d.text((74, 54), "R2", fill=NAVY, font=F_TITLE)
    d.text((254, 405), "fit quality", fill=MUTED, font=F_SMALL, anchor="mm")
    return img


def rmse_icon():
    img, d = new_icon()
    rounded_card(d)
    draw_axes(d)
    curve = [(112, 330), (165, 292), (222, 250), (282, 210), (346, 175), (406, 150)]
    draw_polyline(d, curve, TEAL, 10)
    points = [(132, 340), (178, 252), (226, 292), (286, 165), (352, 222), (398, 122)]
    curve_y = [322, 284, 247, 207, 172, 152]
    for (x, y), cy in zip(points, curve_y):
        d.line([x, y, x, cy], fill=ORANGE, width=5)
        d.ellipse([x - 12, y - 12, x + 12, y + 12], fill=NAVY)
    d.text((70, 54), "RMSE", fill=NAVY, font=font(38, True))
    d.text((254, 405), "fitting error", fill=MUTED, font=F_SMALL, anchor="mm")
    return img


def adjp_icon():
    img, d = new_icon()
    rounded_card(d)
    base_y = 365
    left = 92
    scale_x = 320 / 100
    pts = []
    for i in range(101):
        x = left + i * scale_x
        z = (i - 38) / 17
        y = base_y - 205 * exp(-0.5 * z * z)
        pts.append((x, y))
    area = [(left, base_y)] + pts + [(left + 100 * scale_x, base_y)]
    d.polygon(area, fill="#EAF1F7")
    draw_polyline(d, pts, BLUE, 8)
    threshold_x = left + 66 * scale_x
    d.line([threshold_x, 130, threshold_x, base_y], fill=ORANGE, width=8)
    d.polygon([(threshold_x - 15, 128), (threshold_x + 15, 128), (threshold_x, 104)], fill=ORANGE)
    d.line([left, base_y, left + 320, base_y], fill=MUTED, width=7)
    d.text((70, 54), "adj p", fill=NAVY, font=font(38, True))
    d.text((254, 405), "significance", fill=MUTED, font=F_SMALL, anchor="mm")
    return img


def fc_icon():
    img, d = new_icon()
    rounded_card(d)
    center_x = 256
    d.line([center_x, 132, center_x, 380], fill=LINE, width=5)
    d.line([138, 256, 374, 256], fill=LINE, width=5)
    # Up arrow
    d.line([196, 332, 196, 178], fill=TEAL, width=22)
    d.polygon([(196, 118), (154, 188), (238, 188)], fill=TEAL)
    # Down arrow
    d.line([316, 178, 316, 332], fill=ORANGE, width=22)
    d.polygon([(316, 392), (274, 322), (358, 322)], fill=ORANGE)
    d.text((70, 54), "log2FC", fill=NAVY, font=font(36, True))
    d.text((256, 415), "response size", fill=MUTED, font=F_SMALL, anchor="mm")
    return img


def save_all():
    icons = {
        "icon_r2_fit_quality.png": r2_icon(),
        "icon_rmse_fitting_error.png": rmse_icon(),
        "icon_adjusted_p_significance.png": adjp_icon(),
        "icon_log2fc_response_size.png": fc_icon(),
    }
    for name, img in icons.items():
        img.save(OUT / name)

    # White-background overview for quick inspection.
    panel = Image.new("RGB", (1120, 620), WHITE)
    positions = [(40, 40), (310, 40), (580, 40), (850, 40)]
    for img, pos in zip(icons.values(), positions):
        panel.paste(img.resize((230, 230), Image.Resampling.LANCZOS), pos, img.resize((230, 230), Image.Resampling.LANCZOS))
    d = ImageDraw.Draw(panel)
    labels = ["R2", "RMSE", "Adjusted p", "log2FC"]
    for label, (x, y) in zip(labels, positions):
        d.text((x + 115, y + 260), label, fill=NAVY, font=font(24, True), anchor="mm")
    panel.save(OUT / "qc_metric_icons_overview.png")


if __name__ == "__main__":
    save_all()
