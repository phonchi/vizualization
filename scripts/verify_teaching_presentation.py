"""Regression checks for cached diagnostics and CommonMark emphasis failures."""
from pathlib import Path
import contextlib,io,json,os,subprocess,sys,tempfile
import nbformat,jupytext
import fix_teaching_format as formatting
import sync_teaching_notebooks as synchronizer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reference/notes/refresh_20260913_exhibition'
checks=[]

def check(name,condition):
    checks.append({'name':name,'passed':bool(condition)})
    assert condition,name

source='**規模（magnitude）**衡量地震大小。'
check('original Chinese punctuation case reproduces unparsed emphasis','**' in formatting.PARSER.renderInline(source))
fixed=formatting.normalize(source)
check('normalized definition has a real strong element','<strong>規模</strong>' in formatting.PARSER.renderInline(fixed) and '**' not in formatting.PARSER.renderInline(fixed))
code='```python\nprint("**規模（magnitude）**接字", 2**3)\n```\n'
check('code fence is preserved exactly',formatting.normalize(code)==code)
check('content after code fence is still normalized','<strong>規模</strong>' in formatting.PARSER.render(formatting.normalize(code+source)))
admonition='```{admonition} 定義\n'+source+'\n```\n'+source
check('admonition and following prose both normalize',formatting.normalize(admonition).count('**規模**（magnitude）')==2)
with tempfile.TemporaryDirectory(prefix='teaching-output-regression-') as temporary:
    root=Path(temporary);book=root/'book';cache=book/'_build/jupyter_execute';cache.mkdir(parents=True)
    sourcefile=book/'test.py';sourcefile.write_text('# %% tags=["remove-input"]\nimport math\n')
    old=nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell('import math',execution_count=1,outputs=[nbformat.v4.new_output('stream',name='stderr',text='old Matplotlib cache warning')])])
    new=nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell('import math',execution_count=2,outputs=[])])
    nbformat.write(old,cache/'test.ipynb');nbformat.write(new,book/'test.ipynb')
    synchronizer.ROOT=root;synchronizer.BOOK=book
    with contextlib.redirect_stdout(io.StringIO()):synchronizer.sync()
    saved=nbformat.read(book/'test.ipynb',as_version=4)
    check('fresh executed empty output overrides stale cached stderr',saved.cells[0].outputs==[] and saved.cells[0].execution_count==2)
    nbformat.write(old,book/'test.ipynb')
    try:
        with contextlib.redirect_stdout(io.StringIO()):synchronizer.sync()
    except RuntimeError:
        check('unresolved cached stderr stops synchronization',True)
    else:check('unresolved cached stderr stops synchronization',False)
    bad=root/'not-a-directory';bad.write_text('fixture')
    env={**os.environ,'MPLCONFIGDIR':str(bad),'PYTHONPATH':str(ROOT/'scripts'),'OPENBLAS_NUM_THREADS':'1'}
    code='from _teaching_runtime import prepare_plotting_environment; prepare_plotting_environment(); import matplotlib; print(matplotlib.get_configdir())'
    run=subprocess.run([sys.executable,'-c',code],env=env,capture_output=True,text=True)
    check('unusable configured path falls back without Matplotlib stderr',run.returncode==0 and not run.stderr and run.stdout.strip()==str(ROOT/'data/cache/matplotlib'))
report={'passed':all(r['passed'] for r in checks),'checks':checks,'corrected_example':fixed}
OUT.mkdir(parents=True,exist_ok=True);(OUT/'presentation_regressions.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps(report,ensure_ascii=False))
