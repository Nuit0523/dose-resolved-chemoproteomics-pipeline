from pathlib import Path
from zipfile import ZipFile
import re
import xml.etree.ElementTree as ET

docx = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\work\thesis_for_docking_lookup.docx")
out = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\work\thesis_text.txt")
ctx_out = Path(r"C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\work\docking_contexts.txt")

ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
with ZipFile(docx) as z:
    xml = z.read("word/document.xml")

root = ET.fromstring(xml)
paras = []
for para in root.findall(".//w:p", ns):
    texts = [t.text or "" for t in para.findall(".//w:t", ns)]
    s = "".join(texts).strip()
    if s:
        paras.append(s)

text = "\n".join(paras)
out.write_text(text, encoding="utf-8")

patterns = [
    "docking",
    "Docking",
    "AlphaFold",
    "structure",
    "structural",
    "Cys",
    "NUP205",
    "IPO5",
    "IPO7",
    "PSME1",
    "DHX57",
    "TXNDC5",
]

lines = [f"paras\t{len(paras)}", f"chars\t{len(text)}"]
for pat in patterns:
    lines.append(f"{pat}\t{len(re.findall(re.escape(pat), text))}")

contexts = []
for m in re.finditer(r"(?i)docking|alphafold|structural|binding pose|pocket", text):
    start = max(0, m.start() - 700)
    end = min(len(text), m.end() + 1200)
    contexts.append("\n---CTX---\n" + text[start:end].replace("\n", " "))

ctx_out.write_text("\n".join(lines) + "\n" + "\n".join(contexts), encoding="utf-8")
print("\n".join(lines))
