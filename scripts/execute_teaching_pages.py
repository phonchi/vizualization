"""Execute explicitly selected teaching pages with external downloads disabled."""
from pathlib import Path
import argparse
import json
import time
import jupytext
import nbformat
from nbclient import NotebookClient

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reference/notes/rewrite_20260911'
GUARD='''import requests
def _offline_request(*args, **kwargs):
    raise RuntimeError("Teaching validation is offline; prepare the local cache first")
requests.sessions.Session.request = _offline_request
import gdms_toolkit as _gt
_gt.GDMSSession = _offline_request
'''

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('pages',nargs='+')
    args=parser.parse_args()
    results=[]
    OUT.mkdir(parents=True,exist_ok=True)
    for name in args.pages:
        start=time.monotonic()
        source=ROOT/'book'/f'{name}.py'
        notebook=jupytext.read(source)
        notebook.cells.insert(0,nbformat.v4.new_code_cell(GUARD,metadata={'tags':['remove-cell']}))
        print(f'START {name}',flush=True)
        try:
            client=NotebookClient(notebook,timeout=600,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}})
            client.execute()
            notebook.cells.pop(0)
            for c in notebook.cells:
                if c.cell_type=='code':
                    tags=[t for t in c.metadata.get('tags',[]) if t!='hide-input']
                    c.metadata['tags']=list(dict.fromkeys(tags+['remove-input']))
            nbformat.write(notebook,source.with_suffix('.ipynb'))
            # Keep a reusable executed copy whose code sources can be matched by sync.
            cache=ROOT/'book/_build/jupyter_execute'/f'{name}.ipynb'
            cache.parent.mkdir(parents=True,exist_ok=True)
            nbformat.write(notebook,cache)
            row={'page':name,'status':'passed','seconds':round(time.monotonic()-start,2),
                 'outputs':sum(len(c.get('outputs',[])) for c in notebook.cells)}
        except Exception as exc:
            row={'page':name,'status':'failed','seconds':round(time.monotonic()-start,2),'error':str(exc)}
            print(str(exc),flush=True)
        results.append(row)
        print(json.dumps(row,ensure_ascii=False),flush=True)
        (OUT/'execution_latest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    if any(r['status']!='passed' for r in results):
        raise SystemExit(1)

if __name__=='__main__':main()
