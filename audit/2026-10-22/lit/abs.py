import sys, json, urllib.request, urllib.parse, re

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# usage: abs.py <doi> [regex] -- prints the Europe PMC abstract for a DOI, and context for the regex in the abstract.
UA = {"User-Agent": "cleanTMLE-audit/1.0 (mailto:amertens@berkeley.edu)"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


for doi in sys.argv[1:]:
    if doi.startswith("--"):
        continue
    q = urllib.parse.quote(f'DOI:"{doi}"')
    d = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core&pageSize=1")
    rs = d["resultList"]["result"]
    if not rs:
        print("NO EUROPEPMC RECORD", doi)
        continue
    r = rs[0]
    print("=====", doi, "|", r.get("title"), "|", r.get("journalTitle"), r.get("pubYear"), "| pmid", r.get("pmid"))
    ab = re.sub(r"<[^>]+>", "", r.get("abstractText") or "")
    print(ab[:1800] if ab else "(no abstract)")
    print()
