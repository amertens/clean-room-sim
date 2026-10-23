import csv
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, ".")
from entries import E  # noqa: E402

# Writes three data files (tab-separated) from entries.py and the search logs:
#   publications.tsv   one row per verified publication or resource, with how it was verified and its bib status
#   zotero_to_add.tsv  verified items that are not in writing/references.bib
#   search_summary.tsv hits per source and unique records kept, per logged search
BIB = "../../../writing/references.bib"
IDX = "bib_index.tsv"

bib_text = open(BIB, encoding="utf-8").read()
idx = [l.split("\t") for l in open(IDX, encoding="utf-8").read().split("\n") if l.strip()]
by_doi = {}
for r in idx:
    if len(r) > 2 and r[2]:
        by_doi.setdefault(r[2].lower(), []).append(r[0])


def bib_status(e):
    if e["bib"]:
        return e["bib"]
    d = e["doi"]
    if d:
        keys = by_doi.get(d.lower(), [])
        if keys:
            extra = f" (+{len(keys) - 1} duplicate keys)" if len(keys) > 1 else ""
            return keys[0] + extra
    m = re.search(r"arXiv:(\d{4}\.\d{4,5})", e["ident"])
    if m:
        keys = by_doi.get("10.48550/arxiv." + m.group(1), [])
        return keys[0] if keys else "no"
    u = e["ident"]
    if u.startswith("http"):
        core = re.sub(r"^https?://", "", u).rstrip("/")
        return "yes (URL)" if core in bib_text else "no"
    return "no"


ORDER = {s: i for i, s in enumerate("abcdefg")}
rows = sorted(E, key=lambda e: (ORDER[e["strand"]], e["prio"], -e["year"], e["first"]))
status = [bib_status(e) for e in rows]

with open("publications.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["strand", "priority", "first_author", "year", "short_title", "venue", "identifier", "verified_by", "in_references_bib", "why_it_matters"])
    for e, st in zip(rows, status):
        w.writerow([e["strand"], e["prio"], e["first"], e["year"], e["short"], e["venue"], e["ident"], e["ver"], st, e["why"]])

with open("zotero_to_add.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["priority", "first_author_year", "add_by", "strand"])
    for e, st in sorted(zip(rows, status), key=lambda p: (p[0]["prio"], ORDER[p[0]["strand"]], -p[0]["year"])):
        if st.startswith("no"):
            if e["ident"].startswith("http"):
                add = e["ident"]
            elif "arXiv:" in e["ident"]:
                add = re.search(r"arXiv:\d{4}\.\d{4,5}", e["ident"]).group(0)
            else:
                dois = re.findall(r"10\.\d{4,9}/[^\s;]+", e["ident"])
                pm = re.findall(r"PMID \d+", e["ident"])
                add = " and ".join("DOI " + d for d in dois) + ("; " + pm[0] if pm else "")
            w.writerow([e["prio"], f"{e['first']} {e['year']} {e['short']}", add, e["strand"]])

names = {"pubmed": "PubMed", "europepmc": "Europe PMC", "openalex": "OpenAlex", "semanticscholar": "Semantic Scholar", "arxiv": "arXiv"}
import glob  # noqa: E402

with open("search_summary.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["log", "search_date_utc", "hits_by_source", "unique_records_kept"])
    for lf in sorted(glob.glob("*-log.csv")):
        rr = list(csv.DictReader(open(lf, encoding="utf-8")))
        parts, uniq, dt = [], "", ""
        for r in rr:
            q = r["query_id"].split("-")[-1]
            if q == "DEDUP":
                uniq = r["records_returned"]
                continue
            dt = r["search_date_utc"][:10]
            h = r["hit_count"] or "failed"
            if r["http_status"] != "200":
                h = f"HTTP {r['http_status']}"
            parts.append(f"{names.get(q, q)} {h}")
        w.writerow([lf[:-8], dt, "; ".join(parts), uniq])

print("rows", len(rows), "| not in bib:", sum(1 for s in status if s.startswith("no")), "| in bib:", sum(1 for s in status if not s.startswith("no")))
for e, s in zip(rows, status):
    if not s.startswith("no"):
        print("IN BIB:", e["first"], e["year"], e["short"], "->", s)
