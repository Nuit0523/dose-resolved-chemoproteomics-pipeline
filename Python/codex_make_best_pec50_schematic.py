from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\best_pec50_schematic")
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10233F"
MUTED = "#5D6675"
LINE = "#C8D0DA"
FAINT = "#F5F7FA"
BLUE = "#2E6E9E"
TEAL = "#238C82"
ORANGE = "#D36B2A"
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


F_TITLE = font(34, True)
F_SUB = font(20)
F_LABEL = font(22)
F_BOLD = font(24, True)
F_SMALL = font(18)
F_TINY = font(15)


def draw_arrow(d, start, end, color=ORANGE, width=7):
    x1, y1 = start
    x2, y2 = end
    d.line([x1, y1, x2, y2], fill=color, width=width)
    if x2 > x1:
        d.polygon([(x2, y2), (x2 - 18, y2 - 12), (x2 - 18, y2 + 12)], fill=color)
    else:
        d.polygon([(x2, y2), (x2 + 18, y2 - 12), (x2 + 18, y2 + 12)], fill=color)


def make_card(path, transparent=False):
    W, H = 1240, 520
    bg = (255, 255, 255, 0) if transparent else WHITE
    img = Image.new("RGBA", (W, H), bg)
    d = ImageDraw.Draw(img)

    if not transparent:
        d.rounded_rectangle([18, 18, W - 18, H - 18], radius=16, fill=WHITE, outline=LINE, width=3)

    # Header
    d.text((54, 48), "Site-level curves are collapsed to one gene-level score", fill=NAVY, font=F_TITLE)
    d.text((54, 92), "best_pEC50 keeps the strongest dose-responsive cysteine as the gene summary.", fill=MUTED, font=F_SUB)

    # Protein schematic
    px0, py = 92, 250
    px1 = 510
    d.line([px0, py, px1, py], fill=NAVY, width=16)
    d.rounded_rectangle([px0 - 18, py - 28, px1 + 18, py + 28], radius=28, outline=NAVY, width=4)
    d.text((px0, 174), "Gene X / Protein X", fill=NAVY, font=F_BOLD)

    sites = [
        ("Cys45", 150, "5.2", TEAL),
        ("Cys108", 288, "3.8", BLUE),
        ("Cys302", 440, "4.1", ORANGE),
    ]
    for label, x, pec50, color in sites:
        d.ellipse([x - 22, py - 22, x + 22, py + 22], fill=color, outline=WHITE, width=5)
        d.line([x, py + 24, x, py + 62], fill=color, width=4)
        d.text((x, py + 74), label, fill=NAVY, font=F_SMALL, anchor="mm")
        d.text((x, py + 100), f"pEC50={pec50}", fill=color, font=F_TINY, anchor="mm")

    # Table
    tx, ty = 630, 165
    tw, row_h = 270, 56
    d.rounded_rectangle([tx, ty, tx + tw, ty + row_h * 4], radius=8, fill=FAINT, outline=LINE, width=2)
    d.rectangle([tx, ty, tx + tw, ty + row_h], fill=NAVY)
    d.text((tx + 26, ty + 16), "Site", fill=WHITE, font=F_LABEL)
    d.text((tx + 156, ty + 16), "pEC50", fill=WHITE, font=F_LABEL)

    table_rows = [("Cys45", "5.2", TEAL), ("Cys108", "3.8", BLUE), ("Cys302", "4.1", ORANGE)]
    for i, (site, val, color) in enumerate(table_rows, start=1):
        y = ty + row_h * i
        if i > 1:
            d.line([tx, y, tx + tw, y], fill=LINE, width=2)
        d.text((tx + 26, y + 16), site, fill=NAVY, font=F_LABEL)
        d.text((tx + 166, y + 16), val, fill=color, font=F_LABEL if val != "5.2" else F_BOLD)
        if val == "5.2":
            d.rounded_rectangle([tx + 146, y + 10, tx + 230, y + 45], radius=6, outline=TEAL, width=3)

    draw_arrow(d, (920, 276), (1000, 276), ORANGE, 8)

    # Output box
    ox, oy = 1010, 200
    d.rounded_rectangle([ox, oy, ox + 170, oy + 150], radius=12, fill="#FFF4EC", outline="#E5B28A", width=3)
    d.text((ox + 85, oy + 38), "Gene-level", fill=NAVY, font=F_SMALL, anchor="mm")
    d.text((ox + 85, oy + 75), "best_pEC50", fill=NAVY, font=F_BOLD, anchor="mm")
    d.text((ox + 85, oy + 116), "5.2", fill=ORANGE, font=font(42, True), anchor="mm")

    # Bottom takeaway
    d.text(
        (W / 2, 462),
        "Rationale: one strongly reactive cysteine can define a targetable protein; averaging can dilute that signal.",
        fill=NAVY,
        font=F_SUB,
        anchor="mm",
    )

    img.save(path)


if __name__ == "__main__":
    make_card(OUT / "best_pec50_site_to_gene_schematic.png", transparent=False)
    make_card(OUT / "best_pec50_site_to_gene_schematic_transparent.png", transparent=True)
