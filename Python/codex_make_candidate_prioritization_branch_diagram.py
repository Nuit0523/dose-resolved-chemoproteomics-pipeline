from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs")
PNG = OUT_DIR / "candidate_prioritization_branch_diagram.png"
SVG = OUT_DIR / "candidate_prioritization_branch_diagram.svg"

W, H = 1500, 820
NAVY = "#10233F"
MUTED = "#5D6675"
LINE = "#C8D0DA"
TEAL = "#238C82"
TEAL_LIGHT = "#F3FAF8"
ORANGE = "#D36B2A"
ORANGE_LIGHT = "#FFF4EC"
BLUE = "#2E6E9E"
BLUE_LIGHT = "#F0F6FA"
WHITE = "#FFFFFF"

FONT = Path(r"C:\Windows\Fonts\arial.ttf")
BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")
font_title = ImageFont.truetype(str(BOLD), 42)
font_box = ImageFont.truetype(str(BOLD), 34)
font_small = ImageFont.truetype(str(FONT), 22)


def text_center(draw, xy, text, font, fill=NAVY):
    x, y, w, h = xy
    lines = text.split("\n")
    line_heights = []
    max_width = 0
    for line in lines:
        bbox = draw.textbbox((0, 0), line, font=font)
        line_heights.append(bbox[3] - bbox[1])
        max_width = max(max_width, bbox[2] - bbox[0])
    total_h = sum(line_heights) + (len(lines) - 1) * 8
    cy = y + (h - total_h) / 2
    for line, lh in zip(lines, line_heights):
        bbox = draw.textbbox((0, 0), line, font=font)
        tw = bbox[2] - bbox[0]
        draw.text((x + (w - tw) / 2, cy), line, font=font, fill=fill)
        cy += lh + 8


def box(draw, x, y, w, h, label, fill, accent):
    draw.rounded_rectangle((x, y, x + w, y + h), radius=8, fill=fill, outline=LINE, width=2)
    draw.rounded_rectangle((x, y, x + 12, y + h), radius=6, fill=accent)
    text_center(draw, (x + 22, y + 10, w - 44, h - 20), label, font_box, NAVY)


def line(draw, p1, p2, color=NAVY, width=6):
    draw.line((p1, p2), fill=color, width=width)


def arrow_down(draw, x, y1, y2, color=NAVY):
    line(draw, (x, y1), (x, y2 - 18), color)
    draw.polygon([(x, y2), (x - 14, y2 - 22), (x + 14, y2 - 22)], fill=color)


def arrow_to_box(draw, x1, y1, x2, y2, color=NAVY):
    line(draw, (x1, y1), (x2 - 22, y2), color)
    draw.polygon([(x2, y2), (x2 - 24, y2 - 14), (x2 - 24, y2 + 14)], fill=color)


def make_png():
    img = Image.new("RGB", (W, H), WHITE)
    d = ImageDraw.Draw(img)

    d.text((90, 64), "Candidate Prioritization Branch", font=font_title, fill=NAVY)
    d.text(
        (92, 120),
        "Candidate prioritization generates two parallel outputs: pathway-level biology and site-level structural follow-up.",
        font=font_small,
        fill=MUTED,
    )

    top = (500, 190, 500, 94)
    left = (170, 430, 390, 110)
    right = (940, 430, 390, 110)
    bottom = (940, 625, 390, 110)

    box(d, *top, "Candidate\nprioritization", ORANGE_LIGHT, ORANGE)
    box(d, *left, "Pathway\nenrichment", TEAL_LIGHT, TEAL)
    box(d, *right, "High-confidence\nsites", BLUE_LIGHT, BLUE)
    box(d, *bottom, "Structure-guided\nprioritization", ORANGE_LIGHT, ORANGE)

    # Branch from candidate prioritization.
    stem_x = top[0] + top[2] / 2
    stem_y1 = top[1] + top[3]
    junction_y = 355
    arrow_down(d, stem_x, stem_y1 + 4, junction_y, NAVY)
    line(d, (stem_x, junction_y), (left[0] + left[2] / 2, junction_y), NAVY)
    arrow_down(d, left[0] + left[2] / 2, junction_y, left[1], TEAL)
    line(d, (stem_x, junction_y), (right[0] + right[2] / 2, junction_y), NAVY)
    arrow_down(d, right[0] + right[2] / 2, junction_y, right[1], BLUE)

    # Structural branch.
    arrow_down(d, right[0] + right[2] / 2, right[1] + right[3] + 8, bottom[1], ORANGE)

    d.text((226, 564), "ORA / GSEA", font=font_small, fill=TEAL)
    d.text((1000, 564), "strict site evidence", font=font_small, fill=BLUE)
    d.text((980, 758), "AlphaFold3 + docking follow-up", font=font_small, fill=ORANGE)

    img.save(PNG, quality=95)


def make_svg():
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <rect width="100%" height="100%" fill="{WHITE}"/>
  <text x="90" y="96" font-family="Arial" font-size="42" font-weight="700" fill="{NAVY}">Candidate Prioritization Branch</text>
  <text x="92" y="132" font-family="Arial" font-size="22" fill="{MUTED}">Candidate prioritization generates two parallel outputs: pathway-level biology and site-level structural follow-up.</text>
  <defs>
    <marker id="arrow-navy" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth"><path d="M2,2 L10,6 L2,10 Z" fill="{NAVY}"/></marker>
    <marker id="arrow-teal" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth"><path d="M2,2 L10,6 L2,10 Z" fill="{TEAL}"/></marker>
    <marker id="arrow-blue" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth"><path d="M2,2 L10,6 L2,10 Z" fill="{BLUE}"/></marker>
    <marker id="arrow-orange" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth"><path d="M2,2 L10,6 L2,10 Z" fill="{ORANGE}"/></marker>
  </defs>
  <g font-family="Arial" font-weight="700" font-size="34" fill="{NAVY}" text-anchor="middle">
    <rect x="500" y="190" width="500" height="94" rx="8" fill="{ORANGE_LIGHT}" stroke="{LINE}" stroke-width="2"/><rect x="500" y="190" width="12" height="94" rx="6" fill="{ORANGE}"/>
    <text x="750" y="228"><tspan x="750">Candidate</tspan><tspan x="750" dy="40">prioritization</tspan></text>
    <rect x="170" y="430" width="390" height="110" rx="8" fill="{TEAL_LIGHT}" stroke="{LINE}" stroke-width="2"/><rect x="170" y="430" width="12" height="110" rx="6" fill="{TEAL}"/>
    <text x="365" y="475"><tspan x="365">Pathway</tspan><tspan x="365" dy="40">enrichment</tspan></text>
    <rect x="940" y="430" width="390" height="110" rx="8" fill="{BLUE_LIGHT}" stroke="{LINE}" stroke-width="2"/><rect x="940" y="430" width="12" height="110" rx="6" fill="{BLUE}"/>
    <text x="1135" y="475"><tspan x="1135">High-confidence</tspan><tspan x="1135" dy="40">sites</tspan></text>
    <rect x="940" y="625" width="390" height="110" rx="8" fill="{ORANGE_LIGHT}" stroke="{LINE}" stroke-width="2"/><rect x="940" y="625" width="12" height="110" rx="6" fill="{ORANGE}"/>
    <text x="1135" y="670"><tspan x="1135">Structure-guided</tspan><tspan x="1135" dy="40">prioritization</tspan></text>
  </g>
  <path d="M750 284 V355" stroke="{NAVY}" stroke-width="6" fill="none" marker-end="url(#arrow-navy)"/>
  <path d="M750 355 H365" stroke="{NAVY}" stroke-width="6" fill="none"/>
  <path d="M365 355 V430" stroke="{TEAL}" stroke-width="6" fill="none" marker-end="url(#arrow-teal)"/>
  <path d="M750 355 H1135" stroke="{NAVY}" stroke-width="6" fill="none"/>
  <path d="M1135 355 V430" stroke="{BLUE}" stroke-width="6" fill="none" marker-end="url(#arrow-blue)"/>
  <path d="M1135 540 V625" stroke="{ORANGE}" stroke-width="6" fill="none" marker-end="url(#arrow-orange)"/>
  <text x="226" y="586" font-family="Arial" font-size="22" fill="{TEAL}">ORA / GSEA</text>
  <text x="1000" y="586" font-family="Arial" font-size="22" fill="{BLUE}">strict site evidence</text>
  <text x="980" y="780" font-family="Arial" font-size="22" fill="{ORANGE}">AlphaFold3 + docking follow-up</text>
</svg>'''
    SVG.write_text(svg, encoding="utf-8")


OUT_DIR.mkdir(parents=True, exist_ok=True)
make_png()
make_svg()
print(PNG)
print(SVG)
