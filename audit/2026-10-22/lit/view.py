import sys
sys.stdout.reconfigure(encoding="utf-8",errors="replace")
import json,sys,re
# usage: view.py file.jsonl [regex] [--abs]
f=sys.argv[1]; rx=sys.argv[2] if len(sys.argv)>2 and not sys.argv[2].startswith('--') else None
ab='--abs' in sys.argv
r=[json.loads(l) for l in open(f,encoding='utf-8')]
r.sort(key=lambda x:(-(int(x['year']) if x.get('year') else 0)))
n=0
for x in r:
    t=(x['title'] or '')
    a=(x['abstract'] or '')
    if rx and not re.search(rx,t+' '+a,re.I): continue
    n+=1
    au=(x['authors'] or [''])[0] if x.get('authors') else ''
    print(f"{x['rec_id']} {x['year']} | {au} | {t[:150]} | {(x['venue'] or '')[:40]} | doi={x['doi']} pmid={x['pmid']} arx={x['arxiv_id']}")
    if ab: print('    ',a[:500].replace('\n',' '))
print('shown',n,'of',len(r))
