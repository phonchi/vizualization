"""Refresh three literal plot labels while preserving every saved data value."""
from pathlib import Path
import copy,hashlib,json,re
import jupytext,nbformat
ROOT=Path(__file__).resolve().parents[3];BOOK=ROOT/'book'
source=BOOK/'03_groundwater.py';paired=source.with_suffix('.ipynb')
old=nbformat.read(paired,4);new=jupytext.read(source)
labels={'壯圍（TUN）觀測井，2024/03–04（時間為 UTC）':'壯圍（TUN）觀測井，2024/03/01–05/01（UTC）',
        '一週的水位：接近半日尺度的起伏':'3 月 10–17 日水位：接近半日尺度的起伏',
        'PSD':'功率譜密度（cm²/Hz）'}
old_codes=[c for c in old.cells if c.cell_type=='code'];new_codes=[c for c in new.cells if c.cell_type=='code'];assert len(old_codes)==len(new_codes)
changes=[]
for before,after in zip(old_codes,new_codes):
    expected=before.source
    for a,b in labels.items():expected=expected.replace('"'+a+'"','"'+b+'"')
    assert expected==after.source,'A computational change cannot use label-only refresh'
    after.outputs=copy.deepcopy(before.outputs);after.execution_count=before.execution_count
    if before.source==after.source:continue
    replaced=[]
    for output in after.outputs:
        for mime,content in output.get('data',{}).items():
            if not isinstance(content,str):continue
            def replace(match):
                value=json.loads(match.group(1))
                if value not in labels:return match.group(0)
                replaced.append({'before':value,'after':labels[value]})
                return json.dumps(labels[value],ensure_ascii=True)
            output.data[mime]=re.sub(r'(?<="text":)("(?:[^"\\]|\\.)*")',replace,content)
    assert len(replaced)==1,replaced
    changes.extend(replaced)
assert len(changes)==3
new.metadata=copy.deepcopy(old.metadata)
nbformat.write(new,paired)
nbformat.write(new,BOOK/'_build/jupyter_execute'/paired.name)
report={'passed':True,'method':'Only three literal Plotly labels refreshed; every numeric array and all other outputs preserved',
        'reason':'Full raw archive read exceeded existing 30-second cell limit; no new computation required for these labels',
        'changed_labels':changes,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(Path(__file__).parent/'groundwater_label_refresh.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
