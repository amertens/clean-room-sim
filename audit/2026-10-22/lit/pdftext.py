import sys, re

import fitz  # PyMuPDF

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# usage: pdftext.py file.pdf out.txt [regex ...]
# Extracts text from a downloaded PDF (treated as data) and prints context for each regex.
src, out = sys.argv[1], sys.argv[2]
pats = sys.argv[3:]
doc = fitz.open(src)
txt = "\n".join(f"[p{i+1}] " + p.get_text() for i, p in enumerate(doc))
open(out, "w", encoding="utf-8").write(txt)
print("pages", len(doc), "chars", len(txt))
flat = re.sub(r"\s+", " ", txt)
for p in pats:
    print("=====", p)
    for m in list(re.finditer(p, flat, re.I))[:6]:
        a, b = max(0, m.start() - 200), min(len(flat), m.end() + 300)
        print("...", flat[a:b], "...")
        print()
