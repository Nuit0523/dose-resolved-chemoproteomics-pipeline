from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\gsea_result_redesign")
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10233F"
INK = "#172033"
MUTED = "#5D6675"
LINE = "#C8D0DA"
FAINT = "#F5F7FA"
BLUE = "#2E6E9E"
TEAL = "#238C82"
ORANGE = "#D36B2A"
ORANGE_LIGHT = "#FFF4EC"
WHITE = "#FFFFFF"

W, H = 1600, 900


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for p in candidates:
        if Path(p).exists():
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


F_TITLE = font(39, True)
F_SUB = font(20)
F_HEAD = font(20, True)
F_BODY = font(17)
F_BODY_B = font(17, True)
F_SMALL = font(14)
F_TINY = font(12)
F_TAKE = font(24, True)


def txt(d, xy, value, fill=INK, fnt=F_BODY, anchor=None, align="left", width=None, line_gap=6):
    if width is None:
        d.text(xy, value, fill=fill, font=fnt, anchor=anchor)
        return
    words = value.split()
    lines, cur = [], ""
    for word in words:
        trial = (cur + " " + word).strip()
        if d.textbbox((0, 0), trial, font=fnt)[2] <= width:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    x, y = xy
    for line in lines:
        d.text((x, y), line, fill=fill, font=fnt, anchor=anchor, align=align)
        y += fnt.size + line_gap


def rect(d, box, fill, outline=LINE, width=1, radius=0):
    if radius:
        d.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)
    else:
        d.rectangle(box, fill=fill, outline=outline, width=width)


def draw_table(d, x, y):
    cols = [160, 330, 105, 130]
    headers = ["Ranking", "Pathway signal", "NES", "adjusted p"]
    rows = [
        ["best_pEC50", "intracellular protein transport", "2.923", "0.00097"],
        ["best_pEC50", "nucleocytoplasmic transport", "2.366", "0.0419"],
        ["weighted_score", "protein localization to organelle", "1.688", "3.42E-07"],
    ]
    row_h = 58
    total_w = sum(cols)
    rect(d, [x, y, x + total_w, y + row_h], NAVY, NAVY)
    cx = x
    for h, w in zip(headers, cols):
        txt(d, (cx + 18, y + 20), h, WHITE, font(15, True))
        cx += w
    for i, row in enumerate(rows):
        yy = y + row_h * (i + 1)
        rect(d, [x, yy, x + total_w, yy + row_h], "#F2F5F8" if i % 2 == 0 else WHITE, LINE)
        cx = x
        for j, (cell, w) in enumerate(zip(row, cols)):
            col = TEAL if j == 1 and i < 2 else (ORANGE if j == 1 else INK)
            txt(d, (cx + 18, yy + 19), cell, col, F_BODY_B if j == 1 else F_BODY)
            cx += w
    # vertical lines
    cx = x
    for w in cols[:-1]:
        cx += w
        d.line([cx, y, cx, y + row_h * 4], fill=LINE, width=1)


def draw_rank_schematic(d, x, y):
    rect(d, [x, y, x + 470, y + 345], WHITE, LINE, width=2)
    txt(d, (x + 36, y + 32), "Full ranked gene list", NAVY, F_HEAD)
    txt(d, (x + 36, y + 64), "best_pEC50 ranking, n = 4,994", MUTED, F_SMALL)

    axis_x, axis_y, axis_w = x + 95, y + 132, 300
    txt(d, (axis_x, y + 98), "Top", MUTED, F_TINY, anchor="mm")
    txt(d, (axis_x + axis_w, y + 98), "Bottom", MUTED, F_TINY, anchor="mm")
    rect(d, [axis_x, axis_y, axis_x + axis_w, axis_y + 12], "#DDE5EC", "#DDE5EC")
    txt(d, (axis_x + axis_w / 2, axis_y + 28), "rank 1 -> 4,994", MUTED, F_TINY, anchor="mm")

    genes = [
        ("RANBP2", 71),
        ("TNPO1", 167),
        ("IPO5", 599),
        ("NUP205", 856),
        ("XPO5", 1329),
    ]
    for name, rank in genes:
        gx = axis_x + axis_w * rank / 4994
        rect(d, [gx - 6, axis_y - 22, gx + 6, axis_y + 42], TEAL, TEAL)

    # small rank table
    table_y = y + 195
    txt(d, (x + 38, table_y), "Example genes", NAVY, F_SMALL if False else font(14, True))
    txt(d, (x + 198, table_y), "rank", NAVY, font(14, True))
    txt(d, (x + 282, table_y), "Example genes", NAVY, font(14, True))
    txt(d, (x + 442, table_y), "rank", NAVY, font(14, True), anchor="ra")
    positions = [
        (x + 38, table_y + 32, "RANBP2", "71"),
        (x + 282, table_y + 32, "TNPO1", "167"),
        (x + 38, table_y + 62, "IPO5", "599"),
        (x + 282, table_y + 62, "NUP205", "856"),
        (x + 38, table_y + 92, "XPO5", "1329"),
    ]
    for xx, yy, g, r in positions:
        txt(d, (xx, yy), g, TEAL, font(14, True))
        txt(d, (xx + 180, yy), r, MUTED, F_SMALL, anchor="ra")

    txt(d, (x + 36, y + 315), "Early-shifted, but not confined to one TopN block.", MUTED, F_SMALL)


def main():
    img = Image.new("RGB", (W, H), WHITE)
    d = ImageDraw.Draw(img)

    # top bar
    rect(d, [0, 0, W, 20], NAVY, NAVY)
    txt(d, (74, 62), "GSEA detected distributed transport-related pathway signals.", NAVY, F_TITLE)
    txt(d, (76, 118), "The full ranked gene list revealed coordinated transport signals that hard-cut ORA foreground lists missed.", MUTED, F_SUB)

    # left result table
    rect(d, [70, 180, 895, 545], WHITE, LINE, width=2)
    rect(d, [70, 180, 78, 545], BLUE, BLUE)
    txt(d, (110, 212), "Top GSEA results", NAVY, F_HEAD)
    txt(d, (110, 244), "The strongest signals point to transport and localization biology.", MUTED, F_BODY)
    draw_table(d, 110, 295)

    # right schematic
    rect(d, [955, 180, 1515, 610], ORANGE_LIGHT, LINE, width=2)
    rect(d, [955, 180, 963, 610], ORANGE, ORANGE)
    txt(d, (995, 212), "Why GSEA can detect it", NAVY, F_HEAD)
    draw_rank_schematic(d, 1010, 255)

    # bottom takeaway
    rect(d, [118, 650, 1482, 775], "#F8FAFC", LINE, width=2)
    txt(d, (170, 697), "Take-home", NAVY, F_HEAD)
    txt(
        d,
        (360, 681),
        "Transport-related genes were not confined to Top100/200/400,",
        NAVY,
        F_TAKE,
    )
    txt(
        d,
        (360, 715),
        "but were collectively shifted toward the upper part of the full ranked list.",
        NAVY,
        F_TAKE,
    )

    d.line([70, 830, 1450, 830], fill=LINE, width=1)
    txt(d, (74, 850), "GSEA result summary", MUTED, F_SMALL)
    txt(d, (1518, 850), "8171", MUTED, F_SMALL, anchor="ra")

    img.save(OUT / "gsea_result_redesigned_slide.png")


if __name__ == "__main__":
    main()
