import sys, json, urllib.request, urllib.parse, datetime, re, time

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# usage: verify.py doi:<doi> ... | pmid:<id> ... | arxiv:<id>
# Opens the record at Crossref, PubMed (NCBI esummary) or arXiv and prints the essentials.
# Every lookup is appended to verify-log.jsonl in the working directory.
UA = {"User-Agent": "cleanTMLE-audit/1.0 (mailto:amertens@berkeley.edu)"}
log = open("verify-log.jsonl", "a", encoding="utf-8")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=40) as r:
        return r.read().decode("utf-8")


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    try:
        d = json.loads(get(url))["message"]
    except Exception as e:
        print(f"CROSSREF FAIL {doi}: {e}")
        log.write(json.dumps({"id": doi, "src": "crossref", "status": "fail", "err": str(e), "t": now()}) + "\n")
        return
    au = d.get("author", [])
    a = ", ".join((x.get("family") or x.get("name") or "?") for x in au[:3]) + (" et al." if len(au) > 3 else "")
    yr = (d.get("issued", {}).get("date-parts", [[None]])[0] or [None])[0]
    ct = (d.get("container-title") or [""])[0]
    upd = d.get("update-to") or ""
    title = (d.get("title") or [""])[0]
    print(
        f"CROSSREF OK {doi} | {a} | {yr} | {title[:170]} | {ct} | vol={d.get('volume')} iss={d.get('issue')} "
        f"p={d.get('page') or d.get('article-number')} | type={d.get('type')} | updates={upd}"
    )
    log.write(json.dumps({"id": doi, "src": "crossref", "status": "ok", "authors": a, "year": yr, "title": title,
                          "container": ct, "type": d.get("type"), "t": now()}) + "\n")


def pubmed(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    try:
        d = json.loads(get(url))["result"][pmid]
    except Exception as e:
        print(f"PUBMED FAIL {pmid}: {e}")
        log.write(json.dumps({"id": pmid, "src": "pubmed", "status": "fail", "err": str(e), "t": now()}) + "\n")
        return
    a = ", ".join(x["name"] for x in d["authors"][:3]) + (" et al." if len(d["authors"]) > 3 else "")
    doi = [x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"]
    print(f"PUBMED OK {pmid} | {a} | {d.get('pubdate')} | {d.get('title')[:170]} | {d.get('fulljournalname')} | doi={doi} | types={d.get('pubtype')}")
    log.write(json.dumps({"id": pmid, "src": "pubmed", "status": "ok", "authors": a, "pubdate": d.get("pubdate"),
                          "title": d.get("title"), "journal": d.get("fulljournalname"), "doi": doi, "t": now()}) + "\n")


def arxiv(aid):
    url = f"http://export.arxiv.org/api/query?id_list={aid}"
    try:
        x = get(url)
    except Exception as e:
        print(f"ARXIV FAIL {aid}: {e}")
        log.write(json.dumps({"id": aid, "src": "arxiv", "status": "fail", "err": str(e), "t": now()}) + "\n")
        return
    t = re.search(r"<entry>.*?<title>(.*?)</title>", x, re.S)
    au = re.findall(r"<name>(.*?)</name>", x)
    pub = re.search(r"<published>(.*?)</published>", x)
    jr = re.search(r"<arxiv:journal_ref[^>]*>(.*?)</arxiv:journal_ref>", x, re.S)
    dd = re.search(r"<arxiv:doi[^>]*>(.*?)</arxiv:doi>", x, re.S)
    title = re.sub(r"\s+", " ", t.group(1)).strip() if t else ""
    print(
        f"ARXIV OK {aid} | {', '.join(au[:3])}{' et al.' if len(au) > 3 else ''} | {pub.group(1)[:10] if pub else ''} | "
        f"{title[:170]} | jref={jr.group(1).strip() if jr else ''} | doi={dd.group(1) if dd else ''}"
    )
    log.write(json.dumps({"id": aid, "src": "arxiv", "status": "ok", "authors": au[:6],
                          "published": pub.group(1) if pub else "", "title": title, "t": now()}) + "\n")


def title_search(q):
    """Crossref bibliographic search; prints the top three candidates so a DOI can be located, then verified."""
    url = "https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(q)
    try:
        items = json.loads(get(url))["message"]["items"]
    except Exception as e:
        print(f"TITLE-SEARCH FAIL {q}: {e}")
        return
    for d in items:
        au = d.get("author", [])
        a = ", ".join((x.get("family") or "?") for x in au[:2])
        yr = (d.get("issued", {}).get("date-parts", [[None]])[0] or [None])[0]
        print(f"CANDIDATE for '{q[:50]}': {d.get('DOI')} | {a} | {yr} | {(d.get('title') or [''])[0][:120]} | {(d.get('container-title') or [''])[0]}")
    log.write(json.dumps({"id": q, "src": "crossref-search", "status": "ok", "t": now()}) + "\n")


for a in sys.argv[1:]:
    k, v = a.split(":", 1)
    {"doi": crossref, "pmid": pubmed, "arxiv": arxiv, "title": title_search}[k](v)
    time.sleep(0.4)
