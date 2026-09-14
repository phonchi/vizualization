"""Check paired sources and local links of the currently built reading pages."""
from pathlib import Path
from urllib.parse import urlsplit,unquote
import ast,hashlib,json,os,re
import jupytext,nbformat,yaml
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[1];BOOK=ROOT/'book';HTML=BOOK/'_build/html'
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/reader_review_20260914')
toc=yaml.safe_load((BOOK/'_toc.yml').read_text());names=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']]
failures=[];pairs=[];pages={};documents={};references=0
for name in names:
    source=BOOK/f'{name}.py'
    if source.exists():
        text=source.read_text();ast.parse(text)
        a=jupytext.read(source);b=nbformat.read(source.with_suffix('.ipynb'),4)
        same=[(c.cell_type,c.source) for c in a.cells]==[(c.cell_type,c.source) for c in b.cells]
        if not same:failures.append({'page':name,'issue':'paired cell sources differ'})
        if re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]',text):failures.append({'page':name,'issue':'control character'})
        pairs.append({'page':name,'matched':same,'sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
    file=HTML/f'{name}.html';pages[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
    doc=BeautifulSoup(file.read_text(),'html.parser')
    for tag in doc.select('[href],[src]'):
        raw=tag.get('href',tag.get('src',''));url=urlsplit(raw)
        if url.scheme or url.netloc:continue
        path=(file.parent/unquote(url.path)).resolve() if url.path else file.resolve()
        references+=1
        if not path.exists():failures.append({'page':name,'issue':'missing local target','url':raw});continue
        if url.fragment and path.suffix=='.html':
            if path not in documents:
                target=BeautifulSoup(path.read_text(),'html.parser');documents[path]={e.get('id') for e in target.select('[id]')}|{e.get('name') for e in target.select('a[name]')}
            if unquote(url.fragment) not in documents[path]:failures.append({'page':name,'issue':'missing anchor','url':raw})
OUT.mkdir(parents=True,exist_ok=True)
report={'passed':not failures,'paired_sources':pairs,'reading_pages':len(pages),'html_sha256':pages,'internal_references':references,'failures':failures}
(OUT/'source_integrity.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'passed':not failures,'pairs':len(pairs),'pages':len(pages),'references':references,'failures':failures},ensure_ascii=False))
if failures:raise SystemExit(1)
