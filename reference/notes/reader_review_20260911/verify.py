from pathlib import Path
from urllib.parse import urlsplit,unquote
from bs4 import BeautifulSoup
import json,yaml,jupytext,nbformat,re,subprocess
out=Path('reference/notes/reader_review_20260911');toc=yaml.safe_load(Path('book/_toc.yml').read_text());names=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']]
old_labels={};new_labels={};refs=[];docs=[];pairs=[]
for stem in names:
 p=Path('book')/(stem+'.md')
 if not p.exists():p=p.with_suffix('.py')
 current=p.read_text();old=subprocess.check_output(['git','show','44d70a0:'+str(p)],text=True)
 if p.suffix=='.py':
  a=jupytext.read(p);b=nbformat.read(p.with_suffix('.ipynb'),4)
  assert [(c.cell_type,c.source,c.metadata.get('tags',[])) for c in a.cells]==[(c.cell_type,c.source,c.metadata.get('tags',[])) for c in b.cells],stem
  for cell in b.cells:
   if cell.cell_type=='code':
    compile(cell.source,str(p),'exec');assert 'remove-input' in cell.metadata.get('tags',[])
   assert not any(o.get('output_type')=='error' for o in cell.get('outputs',[]))
  pairs.append({'page':stem,'outputs':sum(len(c.get('outputs',[])) for c in b.cells)})
  text='\n'.join(c.source for c in a.cells if c.cell_type=='markdown')
  old='\n'.join(c.source for c in jupytext.reads(old,fmt='py:percent').cells if c.cell_type=='markdown')
 else:text=current
 for body,target in [(old,old_labels),(text,new_labels)]:
  for label in re.findall(r'\((eq:[\w-]+)\)',body):target.setdefault(label,[]).append(stem)
 assert text.count('## 參考資料與延伸閱讀')==1,stem
 refs.extend(re.findall(r'\{eq\}`([^`]+)`',text))
 docs.extend(d.split('<')[-1].rstrip('>') if '<' in d else d for d in re.findall(r'\{doc\}`([^`]+)`',text))
assert old_labels==new_labels
assert all(len(v)==1 for v in new_labels.values())
assert all(r in new_labels for r in refs)
assert all(d in names for d in docs)
html=Path('book/_build/html').resolve();ids={};count=0
for stem in names:
 p=html/(stem+'.html')
 for link in BeautifulSoup(p.read_text(),'html.parser').select('a[href]'):
  u=urlsplit(link['href'])
  if u.scheme or u.netloc or not u.path and not u.fragment:continue
  dest=(p.parent/unquote(u.path)).resolve() if u.path else p;count+=1
  assert dest.exists(),(stem,dest)
  if dest.suffix=='.html' and u.fragment:
   if dest not in ids:ids[dest]={e['id'] for e in BeautifulSoup(dest.read_text(),'html.parser').select('[id]')}
   assert unquote(u.fragment) in ids[dest],(stem,link['href'])
r={'pages':len(names),'paired_notebooks':len(pairs),'outputs':sum(x['outputs'] for x in pairs),'numbered_equations_preserved':len(new_labels),'source_doc_links':len(docs),'rendered_internal_links':count,'errors':[]}
(out/'structural_validation.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(r)
