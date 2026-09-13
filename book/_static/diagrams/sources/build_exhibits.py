#!/usr/bin/env python3
"""Export fixed teaching data, then emit exhibit fragments with the shared renderer.
No model fitting, forecast generation, or network access. Requires the existing
Italy .mat files and cached HORUS CSV. Run with the project's Python environment.
"""
from pathlib import Path
import csv,json,hashlib
from datetime import datetime,timedelta
import numpy as np
from scipy.io import loadmat
from pyproj import Transformer
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
META=[
('d01_forecast_pipeline','pipeline','一份預報，怎樣成為證據？','沿著資料的時間順序，走完一次實驗。','process','步驟'),
('d02_regions_s_r','map','同一片義大利，兩種區域','疊上真實區域與歷史震央，辨認誰提供資料、誰接受評分。','real','歷史截止年'),
('d03_timeline_windows','timeline','每次發報，只能看見過去','移動發報窗，觀看歷史集合與資料延遲。','real','預報窗'),
('d03_thresholds','magnitude','一個規模，兩種刻度','沿規模軸移動，看原始數值如何取整、進入哪個箱。','synthetic','原始規模'),
('d04_density_to_counts','integration','把一片率，累積成一個數','依序積分小格，再把總期望數轉成計數機率。','synthetic','累積進度'),
('d06_quantile_pvalue','distribution','觀測值，在分布的哪裡？','一起閱讀柱狀分布、累積機率與含同值的尾端。','synthetic','觀測統計量'),
('d09_kernel_smoothing','kernel','震央是一個點，貢獻是一片','移動剖面線，比較等值圈與同一條截線上的貢獻。','synthetic','剖面位置'),
('d11_hawkes_process','hawkes','事件如何改寫下一刻的率？','讓時間前進，看每顆已知事件加入自己的衰減貢獻。','synthetic','時間'),
('d12_etas_branching','branching','一顆事件，可以接出一個家族','逐代揭露親子連線，分開看直接後代與整個家族。','synthetic','世代'),
('d15_eepas_kernels','product','三個核，指向同一份貢獻','看規模、時間與空間的分配如何相乘，再落入一個目標格箱。','synthetic','累積階段'),
('d16_csep_tests_flow','matrix','先固定什麼，再模擬什麼？','選一種檢驗，追蹤總数、位置、規模與時間的角色。','synthetic','模擬進度'),
('d19_convex_mix','mix','兩張預報，逐格相加','播放一份固定的混合示例，看每個成分留下多少期望數。','synthetic','組合階段'),
('d20_model_family','family','模型之間，交換了哪些資訊？','沿關係線探索固定基準、歷史更新與模型組合。','process','模型'),
('d03_forecast_protocols','protocol','回溯、擬前瞻、真前瞻','沿時間軸比較測試資料、模型與建模者的關係。','process','時間位置')]
manifest=[]
for slug,kind,title,deck,source,label in META:manifest.append(dict(slug=slug,kind=kind,title=title,deck=deck,source=source,control_label=label))
# Authoritative local region coordinates. Projection remains kilometres throughout.
mat=ROOT/'EEPAS/data/CELLE_ter.mat';poly=ROOT/'EEPAS/data/CPTI15.mat';cat=ROOT/'data/cache/italy/horus_clean.csv'
cells=loadmat(mat)['CELLESD'];polygon=loadmat(poly)['cpti15'];trans=Transformer.from_crs('EPSG:4326','EPSG:7794',always_xy=True)
px,py=trans.transform(polygon[:,0],polygon[:,1]);points=[]
with cat.open() as f:
 for r in csv.DictReader(f):
  if r['is_earthquake']!='True' or r['in_cpti15']!='True' or float(r['depth'])>40 or float(r['mb'])<5 or not ('1960'<=r['time']<'2022'):continue
  x,y=trans.transform(float(r['lon']),float(r['lat']));xx,yy=x/1000,y/1000
  ci=np.flatnonzero((cells[:,0]<=xx)&(xx<cells[:,1])&(cells[:,2]<=yy)&(yy<cells[:,3]))
  points.append([round(xx,3),round(yy,3),int(r['time'][:4]),float(r['mb']),r['time'].replace(' ','T')+'Z',int(ci[0]) if len(ci) else -1])
data=dict(cells=np.round(cells,3).tolist(),polygon=np.round(np.c_[px/1000,py/1000],3).tolist(),events=points,crs='EPSG:7794',time_standard='UTC',source='CELLE_ter.mat, CPTI15.mat; HORUS 2024-06, M≥5, S內、深度≤40 km、1960–2021')
data['windows']=[dict(index=i+1,start=(datetime(2012,1,1)+timedelta(days=91.31*i)).isoformat()+'Z',end=(datetime(2012,1,1)+timedelta(days=91.31*(i+1))).isoformat()+'Z') for i in range(40)]
data['targets']=[p for p in points if p[2]>=2012 and p[5]>=0 and p[3]<7.5]
(OUT.parent/'data').mkdir(exist_ok=True)
(OUT.parent/'data/italy_exhibit_data.json').write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n')
for m in manifest:
 if m['kind'] in ['map','timeline']:m['data']=data
(OUT/'exhibits.json').write_text(json.dumps(manifest,ensure_ascii=False,separators=(',',':'))+'\n')
receipt=dict(inputs=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [mat,poly,cat]],cells=len(cells),polygon_vertices=len(polygon),events=len(points),scope='lightweight export only; no model computation')
(OUT/'exhibit_data_provenance.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False))
