from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\gene_summary_diagram")
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10233F"
MUTED = "#5D6675"
LINE = "#C8D0DA"
BLUE = "#2E6E9E"
TEAL = "#238C82"
ORANGE = "#D36B2A"
FAINT = "#F5F7FA"
PEACH = "#FFF4EC"
WHITE = "#FFFFFF"
TRANSPARENT = (255, 255, 255, 0)


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


F_TITLE = font(34, True)
F_HEAD = font(25, True)
F_TEXT = font(22)
F_TEXT_B = font(22, True)
F_SMALL = font(18)


def draw_arrow(draw, start, end, color=ORANGE, width=7):
    x1, y1 = start
    x2, y2 = end
    draw.line([x1, y1, x2, y2], fill=color, width=width)
    if x2 > x1:
        pts = [(x2, y2), (x2 - 24, y2 - 15), (x2 - 24, y2 + 15)]
    elif x2 < x1:
        pts = [(x2, y2), (x2 + 24, y2 - 15), (x2 + 24, y2 + 15)]
    else:
        pts = [(x2, y2), (x2 - 15, y2 - 24), (x2 + 15, y2 - 24)]
    draw.polygon(pts, fill=color)


def rounded_rect(draw, box, radius, fill, outline=LINE, width=3):
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def make_diagram(background="transparent"):
    W, H = 1420, 520
    bg = WHITE if background == "white" else TRANSPARENT
    img = Image.new("RGBA", (W, H), bg)
    d = ImageDraw.Draw(img)

    # Main surface
    rounded_rect(d, [28, 28, W - 28, H - 28], 28, WHITE if background == "transparent" else "#FBFCFE", LINE, 3)

    d.text((70, 62), "Site-level curves -> gene-level summary", fill=NAVY, font=F_TITLE)
    d.text((70, 108), "One gene can contain multiple quantified cysteine sites; enrichment requires one gene-level score.", fill=MUTED, font=F_SMALL)

    # Left: protein/gene schematic
    rounded_rect(d, [70, 165, 520, 390], 20, FAINT, LINE, 2)
    d.text((102, 190), "Gene X / Protein X", fill=NAVY, font=F_HEAD)

    # Protein backbone
    d.rounded_rectangle([118, 276, 470, 306], radius=15, fill="#E7EDF3", outline=LINE, width=2)
    sites = [
        ("Cys45", 175, TEAL),
        ("Cys108", 282, ORANGE),
        ("Cys302", 398, TEAL),
    ]
    for label, x, color in sites:
        d.ellipse([x - 18, 258, x + 18, 294], fill=color, outline=WHITE, width=3)
        d.line([x, 294, x, 338], fill=color, width=3)
        d.text((x, 346), label, fill=NAVY, font=F_SMALL, anchor="mm")

    # Middle: site pEC50 values
    rounded_rect(d, [600, 158, 870, 397], 20, WHITE, LINE, 2)
    d.text((630, 188), "Site pEC50", fill=NAVY, font=F_HEAD)

    rows = [
        ("Cys45", "4.1", TEAL, False),
        ("Cys108", "5.2", ORANGE, True),
        ("Cys302", "3.8", TEAL, False),
    ]
    y = 244
    for site, val, color, is_best in rows:
        if is_best:
            d.rounded_rectangle([622, y - 10, 848, y + 30], radius=12, fill=PEACH)
        d.ellipse([638, y, 656, y + 18], fill=color)
        d.text((672, y - 3), site, fill=NAVY, font=F_TEXT_B if is_best else F_TEXT)
        d.text((820, y - 3), val, fill=color if is_best else MUTED, font=F_TEXT_B, anchor="ra")
        y += 52

    # Right: selected summary
    rounded_rect(d, [980, 190, 1340, 360], 20, PEACH, ORANGE, 3)
    d.text((1160, 222), "Gene-level score", fill=NAVY, font=F_HEAD, anchor="mm")
    d.text((1160, 282), "best_pEC50 = 5.2", fill=ORANGE, font=font(31, True), anchor="mm")
    d.text((1160, 326), "strongest responsive site", fill=MUTED, font=F_SMALL, anchor="mm")

    draw_arrow(d, (522, 278), (592, 278), ORANGE, 7)
    draw_arrow(d, (872, 278), (970, 278), ORANGE, 7)

    d.text((650, 438), "Rationale: averaging can dilute a meaningful reactive cysteine.", fill=NAVY, font=F_TEXT_B, anchor="mm")

    return img


def main():
    transparent = make_diagram("transparent")
    transparent.save(OUT / "best_pec50_gene_summary_diagram_transparent.png")

    white = make_diagram("white").convert("RGB")
    white.save(OUT / "best_pec50_gene_summary_diagram_white.png")


if __name__ == "__main__":
    main()
