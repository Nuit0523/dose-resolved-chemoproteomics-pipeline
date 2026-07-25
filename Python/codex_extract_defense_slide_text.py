import json
import pathlib
import re
import zipfile
import xml.etree.ElementTree as ET


PPT = pathlib.Path("work/defense_zip_input/extracted/defense.pptx")
OUT_TXT = pathlib.Path("work/defense_zip_input/slide_text_dump.txt")
OUT_JSON = pathlib.Path("work/defense_zip_input/slide_text_dump.json")

NS = {
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
}


def natural_key(value):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", value)]


def get_text(xmlbytes):
    root = ET.fromstring(xmlbytes)
    blocks = []
    for txbody in root.findall(".//p:txBody", NS):
        paras = []
        for para in txbody.findall("./a:p", NS):
            runs = []
            for text_node in para.findall(".//a:t", NS):
                if text_node.text:
                    runs.append(text_node.text)
            if runs:
                paras.append("".join(runs))
        if paras:
            blocks.append("\n".join(paras))
    if not blocks:
        blocks = [node.text for node in root.findall(".//a:t", NS) if node.text]
    return blocks


def main():
    data = []
    lines = []
    with zipfile.ZipFile(PPT) as archive:
        slides = sorted(
            [
                name
                for name in archive.namelist()
                if re.match(r"ppt/slides/slide\d+\.xml$", name)
            ],
            key=natural_key,
        )
        for slide_number, slide_name in enumerate(slides, 1):
            texts = get_text(archive.read(slide_name))
            notes = []
            rels_name = slide_name.replace("slides/", "slides/_rels/") + ".rels"
            if rels_name in archive.namelist():
                relroot = ET.fromstring(archive.read(rels_name))
                for rel in relroot:
                    if rel.attrib.get("Type", "").endswith("/notesSlide"):
                        target = rel.attrib.get("Target", "")
                        notes_name = "ppt/notesSlides/" + target.split("/")[-1]
                        if notes_name in archive.namelist():
                            notes = get_text(archive.read(notes_name))
            lines.append(f"===== SLIDE {slide_number} =====")
            lines.extend(texts)
            if notes:
                lines.append("--- NOTES ---")
                lines.extend(notes)
            data.append({"slide": slide_number, "texts": texts, "notes": notes})
    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    OUT_JSON.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(OUT_TXT)
    print(f"slides {len(data)}")


if __name__ == "__main__":
    main()
