"""Recompute only the final learning window after clipping at 2012-01-01.
Run from the repository root with PYTHONPATH=. and one BLAS thread.
"""
from pathlib import Path
import hashlib,json,time,resource
import numpy as np
from gdms_toolkit import italy,italy_models as m
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,resource.RLIM_INFINITY))
out=Path(__file__).parent
cat=italy.experiment_catalog()
window=m.learning_windows().tail(1)
rows=[]
for name in m.MODELS:
    path=italy.FORECAST_DIR/f'{name}_learning.npy'
    before=hashlib.sha256(path.read_bytes()).hexdigest()
    old=np.load(path); previous=old[-1].sum(); started=time.monotonic()
    arr=getattr(m,name.lower()+'_forecast')(cat,window)
    old[-1]=arr[0];np.save(path,old)
    row=dict(model=name,window_start=str(window.start.iloc[0]),window_end=str(window.end.iloc[0]),old_final_total=float(previous),new_final_total=float(arr.sum()),seconds=time.monotonic()-started,before_sha256=before,after_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    rows.append(row);print(json.dumps(row),flush=True)
    (out/'learning_boundary_repair.json').write_text(json.dumps(rows,indent=2)+'\n')
