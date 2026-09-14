from pathlib import Path
import json,re,nbformat,pandas as pd
ROOT=Path(__file__).resolve().parents[3];out=Path(__file__).parent;checks={}
n=nbformat.read(ROOT/'book/05_seismic.ipynb',4)
text=''.join(o.get('text','') for c in n.cells for o in c.get('outputs',[]) if o.get('output_type')=='stream')
rows=[x for x in text.splitlines() if '估計 N=' in x]
assert len(rows)==2 and 'N=59' in rows[0] and 'b = 0.88' in rows[0] and 'N=1608' in rows[1] and 'b = 0.73' in rows[1]
images=sum('image/png' in o.get('data',{}) for c in n.cells for o in c.get('outputs',[]));assert images==2
checks['seismic']={'b_outputs':rows,'raster_images':images}
n=nbformat.read(ROOT/'book/06_gnss.ipynb',4)
for c in n.cells:
    for o in c.get('outputs',[]):
        h=o.get('data',{}).get('text/html','');m=re.search(r'Plotly.newPlot\(\s*"[^"]+",\s*',h)
        if not m:continue
        data,_=json.JSONDecoder().raw_decode(h[m.end():]);x=data[0].get('x',[])
        if len(x)==2880:
            t=pd.DatetimeIndex(x)
            assert t.nunique()==2880 and t[-1]==pd.Timestamp('2024-04-02 23:59:30')
            assert (t.to_series().diff().dropna()==pd.Timedelta(seconds=30)).all()
            checks['gnss']={'epochs':len(x),'unique_epochs':t.nunique(),'first':x[0],'last':x[-1],'interval_seconds':30}
assert 'gnss' in checks
checks['passed']=True;(out/'updated_outputs.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n');print(json.dumps(checks,ensure_ascii=False))
