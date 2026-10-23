import sys, json, urllib.request, urllib.parse, time, datetime

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# usage: cites.py <doi> <label> [from_year]
# Lists works that cite a DOI, using the OpenAlex API (cited-by search). Writes cites-<label>.jsonl as a search log.
UA = {"User-Agent": "cleanTMLE-audit/1.0 (mailto:amertens@berkeley.edu)"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


doi, label = sys.argv[1], sys.argv[2]
from_year = int(sys.argv[3]) if len(sys.argv) > 3 else 2024
w = get("https://api.openalex.org/works/doi:" + urllib.parse.quote(doi))
wid = w["id"].split("/")[-1]
print("OpenAlex id", wid, "|", w.get("title"), "| cited_by_count", w.get("cited_by_count"))
out = []
cursor = "*"
while cursor and len(out) < 400:
    url = (f"https://api.openalex.org/works?filter=cites:{wid},from_publication_date:{from_year}-01-01"
           f"&per-page=100&cursor={urllib.parse.quote(cursor)}&select=id,doi,title,publication_year,primary_location,authorships,type")
    d = get(url)
    out += d["results"]
    cursor = d["meta"].get("next_cursor")
    if not d["results"]:
        break
    time.sleep(0.3)
with open(f"cites-{label}.jsonl", "w", encoding="utf-8") as f:
    meta = {"query": "cites " + doi, "openalex_id": wid, "from_year": from_year, "n": len(out),
            "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    f.write(json.dumps({"_meta": meta}) + "\n")
    for r in out:
        f.write(json.dumps(r) + "\n")
print("citing works since", from_year, ":", len(out))
for r in sorted(out, key=lambda r: -(r.get("publication_year") or 0)):
    au = (r.get("authorships") or [{}])[0].get("author", {}).get("display_name", "")
    ven = ((r.get("primary_location") or {}).get("source") or {}).get("display_name", "")
    print(r.get("publication_year"), "|", au, "|", (r.get("title") or "")[:130], "|", ven[:40], "|", r.get("doi"))
