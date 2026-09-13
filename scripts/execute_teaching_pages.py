"""Execute explicitly selected teaching pages with external downloads disabled."""

from _teaching_runtime import prepare_plotting_environment
prepare_plotting_environment()

from pathlib import Path
import argparse
import ast
import hashlib
from datetime import datetime, timezone
import json
import os
import time
import jupytext
import nbformat
from nbclient import NotebookClient

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/refresh_20260913_exhibition/execution')
GUARD='''import requests
def _offline_request(*args, **kwargs):
    raise RuntimeError("Teaching validation is offline; prepare the local cache first")
requests.sessions.Session.request = _offline_request
import gdms_toolkit as _gt
_gt.GDMSSession = _offline_request
'''

def execute_in_process(notebook, diagnostics_path):
    """Socket-free execution with IPython rich outputs, one fresh namespace per page."""
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    from IPython.display import display
    shell = InteractiveShell.instance()
    shell.reset(new_session=True)
    timings=[]
    for i,cell in enumerate(notebook.cells):
        if cell.cell_type != 'code':
            continue
        started=time.monotonic()
        stamp=datetime.now(timezone.utc).isoformat()
        tree=ast.parse(cell.source)
        last=tree.body.pop() if tree.body and isinstance(tree.body[-1],ast.Expr) else None
        with capture_output() as captured:
            exec(compile(tree,f'<cell {i}>','exec'),shell.user_ns)
            if last is not None:
                value=eval(compile(ast.Expression(last.value),f'<cell {i}>','eval'),shell.user_ns)
                if value is not None:
                    display(value)
        seconds=time.monotonic()-started
        cell.outputs=[]
        if captured.stderr:
            with diagnostics_path.open('a') as diagnostic:
                diagnostic.write(f'CELL {i-1} STDERR\n'+captured.stderr+'\n')
            raise RuntimeError(f'Unexpected stderr in cell {i-1}; see {diagnostics_path}')
        if captured.stdout:
            cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        for output in captured.outputs:
            cell.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
        cell.execution_count=len(timings)+1
        cell.metadata['execution']={'iopub.execute_input':stamp,'iopub.status.idle':datetime.now(timezone.utc).isoformat()}
        timings.append({'cell':i-1,'seconds':round(seconds,4)})
        print(f'  cell {i-1}: {seconds:.3f}s',flush=True)
        if seconds>=30:
            raise RuntimeError(f'Cell {i-1} exceeded 30 seconds: {seconds:.2f}')
    return timings[1:]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('pages',nargs='+')
    parser.add_argument('--in-process',action='store_true',help='Use socket-free IPython execution')
    args=parser.parse_args()
    if os.name == "posix":
        import resource
        cap=3*1024**3
        _,hard=resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS,(min(cap,hard) if hard>0 else cap,hard))
    results=[]
    OUT.mkdir(parents=True,exist_ok=True)
    for name in args.pages:
        start=time.monotonic()
        source=ROOT/'book'/f'{name}.py'
        notebook=jupytext.read(source)
        notebook.cells.insert(0,nbformat.v4.new_code_cell(GUARD,metadata={'tags':['remove-cell']}))
        print(f'START {name}',flush=True)
        try:
            if args.in_process:
                timings=execute_in_process(notebook,OUT/f'{name}_diagnostics.log')
            else:
                client=NotebookClient(notebook,timeout=30,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}})
                client.execute()
                timings=[]
            notebook.cells.pop(0)
            stderr=[o.get('text','') for c in notebook.cells for o in c.get('outputs',[]) if o.get('output_type')=='stream' and o.get('name')=='stderr']
            if stderr:
                (OUT/f'{name}_diagnostics.log').write_text('\n'.join(stderr))
                raise RuntimeError(f'Unexpected notebook stderr: {name}')
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
                 'outputs':sum(len(c.get('outputs',[])) for c in notebook.cells),
                 'method':'in-process IPython' if args.in_process else 'nbclient',
                 'cells':timings,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                 'toolkit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'gdms_toolkit').glob('*.py'))}}
        except Exception as exc:
            row={'page':name,'status':'failed','seconds':round(time.monotonic()-start,2),'error':str(exc)}
            print(str(exc),flush=True)
        (OUT/f'{name}_execution.json').write_text(json.dumps(row,ensure_ascii=False,indent=2)+'\n')
        results.append(row)
        print(json.dumps(row,ensure_ascii=False),flush=True)
        (OUT/'execution_latest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    if any(r['status']!='passed' for r in results):
        raise SystemExit(1)

if __name__=='__main__':main()
