"""Compare all published reading pages and theme assets to the accepted local build."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import argparse
import json
import os
import time
import requests
import yaml
from bs4 import BeautifulSoup
from urllib.parse import urlsplit, unquote

ROOT=Path(__file__).resolve().parents[1]
HTML=ROOT/'book/_build/html'
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/refresh_20260913_exhibition')/'publication.json'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',default='https://phonchi.github.io/vizualization/')
    parser.add_argument('--attempts',type=int,default=12)
    args=parser.parse_args()
    toc=yaml.safe_load((ROOT/'book/_toc.yml').read_text())
    files=[toc['root']+'.html']+[c['file']+'.html' for p in toc['parts'] for c in p['chapters']]
    reading_count=len(files)
    reading_files=list(files)
    images=set()
    for name in reading_files:
        for image in BeautifulSoup((HTML/name).read_text(),'html.parser').select('img[src]'):
            url=urlsplit(image['src'])
            if url.scheme or url.netloc:continue
            target=(HTML/name).parent/unquote(url.path)
            if target.is_file():images.add(str(target.relative_to(HTML)))
    files+=['index.html',*[f'_static/{n}' for n in ['teaching.css','museum.css','exhibits.css','teaching_theme.js','museum.js','exhibits.js']], '_static/reading/biondini2023_fig8.png','_static/diagrams/data/italy_exhibit_data.json']
    files += [str(p.relative_to(HTML)) for p in sorted((HTML/'_static/diagrams/standalone').glob('d*.html')) if '.archify.' not in p.name and '.visual-check.' not in p.name]
    files=list(dict.fromkeys(files+sorted(images)))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    local={name:sha256((HTML/name).read_bytes()).hexdigest() for name in files}
    history=[]
    def fetch(name):
        try:
            response=requests.get(args.base.rstrip('/')+'/'+name,timeout=30)
            digest=sha256(response.content).hexdigest()
            return dict(file=name,status=response.status_code,local_sha256=local[name],
                        remote_sha256=digest,matched=response.status_code==200 and digest==local[name])
        except requests.RequestException as exc:
            return dict(file=name,matched=False,error=str(exc))
    for attempt in range(1,args.attempts+1):
        with ThreadPoolExecutor(max_workers=4) as pool:
            results=list(pool.map(fetch,files))
        passed=all(r['matched'] for r in results)
        history.append(dict(attempt=attempt,matched=sum(r['matched'] for r in results)))
        report=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),base=args.base,
                    passed=passed,reading_pages=reading_count,results=results,attempts=history)
        OUT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
        print(json.dumps(dict(attempt=attempt,matched=history[-1]['matched'],total=len(files),passed=passed)),flush=True)
        if passed:return
        if attempt<args.attempts:time.sleep(15)
    raise SystemExit('Published bytes do not match; see '+str(OUT))


if __name__=='__main__':main()
