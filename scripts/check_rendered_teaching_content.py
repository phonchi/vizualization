"""Fail a build if diagnostics or unparsed teaching markup reach a reading page."""
from pathlib import Path
import hashlib,json,os,re
from bs4 import BeautifulSoup
import yaml
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/refresh_20260913_exhibition')


def main():
    toc=yaml.safe_load((ROOT/'book/_toc.yml').read_text());names=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']];rows=[]
    for name in names:
        path=ROOT/'book/_build/html'/f'{name}.html'
        soup=BeautifulSoup(path.read_text(),'html.parser');article=soup.select_one('article.bd-article');problems=[]
        if article is None:problems.append('missing article')
        else:
            text=article.get_text()
            if re.search(r'Matplotlib created|not a\s+writable directory|Traceback \(most recent|RuntimeWarning:|UserWarning:',text):problems.append('execution diagnostic in content')
            for node in article.find_all(string=True):
                if node.find_parent(['code','pre','script','style']):continue
                if '**' in str(node):problems.append('unparsed emphasis: '+str(node).strip()[:130])
            if article.select('.diagram-missing'):problems.append('missing diagram placeholder')
        rows.append({'page':name,'passed':not problems,'problems':problems,'html_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    OUT.mkdir(parents=True,exist_ok=True);report={'passed':all(r['passed'] for r in rows),'pages':rows};(OUT/'rendered_content.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('RENDERED CONTENT',json.dumps({'passed':report['passed'],'pages':len(rows),'failed':[r for r in rows if not r['passed']]},ensure_ascii=False),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
