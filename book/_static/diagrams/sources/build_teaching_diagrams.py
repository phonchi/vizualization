#!/usr/bin/env python3
"""Reproduce twelve inline teaching figures and their offline standalone previews.

The two workflow JSON files preserve the Archify source/acceptance receipts.
Embedded views restate that topology at a readable mobile scale. Remaining
figures are conceptual examples, not the cached Italy experiment's results.
"""
from pathlib import Path
from html import escape
import json, math, sys
SELECTED=set(sys.argv[1:])
WRITTEN=0

ROOT=Path(__file__).resolve().parents[1]
TOKENS='''--dg-paper:#fbfcfd;--dg-paper-2:#eef2f7;--dg-ink:#243040;--dg-muted:#56657a;--dg-soft:#7d8797;--dg-rule:rgba(36,48,64,.12);--dg-accent:#2a78d6;--dg-accent-tint:rgba(42,120,214,.08);--dg-quake:#e34948;'''
DARK='''--dg-paper:#1e232b;--dg-paper-2:#2a313c;--dg-ink:#e9eef4;--dg-muted:#b6c0cd;--dg-soft:#8d98a8;--dg-rule:rgba(233,238,244,.18);--dg-accent:#6aa5ec;--dg-accent-tint:rgba(106,165,236,.12);--dg-quake:#f07a79;'''

def text(x,y,s,cls='',anchor='start'):
 return f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{escape(s)}</text>'
def box(x,y,w,h,label,sub='',focal=False):
 return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" class="node {"focal" if focal else ""}"/>'+text(x+w/2,y+32,label,'name','middle')+(text(x+w/2,y+60,sub,'sub','middle') if sub else '')
def path(d,cls='edge'):
 return f'<path d="{d}" class="{cls}"/>'
def group(key,body):
 return f'<g data-detail-mark="{key}">{body}</g>'
def chart_line(points,cls='curve'):
 return f'<polyline points="{" ".join(f"{x:.2f},{y:.2f}" for x,y in points)}" class="{cls}"/>'
def arrow(x,y):
 return path(f'M {x} {y} v 24 m -8 -8 l 8 8 l 8 -8')


def make(name,title,desc,height,body,choices,lead='',kind='concept',extra=''):
 global WRITTEN
 if SELECTED and name not in SELECTED:return
 WRITTEN+=1
 prefix='dg-'+name
 styles=f'''<style>
:where(:root){{{TOKENS}}}
.{prefix}{{max-width:100%;color:var(--dg-ink);font-family:'Geist','Noto Sans TC','PingFang TC','Microsoft JhengHei',system-ui,sans-serif}}
.{prefix} *{{box-sizing:border-box}}
.{prefix} svg{{display:block;width:100%;max-width:480px;height:auto;margin:0 auto;background:var(--dg-paper)}}
.{prefix} svg text{{fill:var(--dg-ink);font-size:20px}}
.{prefix} svg .name{{font-weight:600}}
.{prefix} svg .sub{{fill:var(--dg-muted);font-size:20px}}
.{prefix} .node{{fill:var(--dg-paper-2);stroke:var(--dg-muted);stroke-width:1.2}}
.{prefix} .focal{{fill:var(--dg-accent-tint);stroke:var(--dg-accent);stroke-width:2}}
.{prefix} .edge{{fill:none;stroke:var(--dg-muted);stroke-width:2;stroke-linejoin:round;stroke-linecap:round}}
.{prefix} .curve{{fill:none;stroke:var(--dg-accent);stroke-width:3;stroke-linejoin:round}}
.{prefix} .curve-other{{fill:none;stroke:var(--dg-muted);stroke-width:2;stroke-dasharray:8 4;stroke-linejoin:round}}
.{prefix} .curve-faint{{fill:none;stroke:var(--dg-muted);stroke-width:2;stroke-dasharray:2 4;stroke-linejoin:round}}
.{prefix} .fill{{fill:var(--dg-accent-tint);stroke:var(--dg-accent);stroke-width:2}}
.{prefix} .dot{{fill:var(--dg-ink);stroke:var(--dg-paper);stroke-width:2}}
.{prefix} .event{{fill:var(--dg-quake);stroke:var(--dg-paper);stroke-width:2}}
.{prefix} .grid{{fill:none;stroke:var(--dg-rule);stroke-width:1}}
.{prefix} .choices{{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}}
.{prefix} .choices label{{display:block;font:inherit;font-size:15px;line-height:1.4;min-height:44px;padding:8px 12px;border:1px solid var(--dg-muted);border-radius:6px;background:var(--dg-paper);color:var(--dg-ink);cursor:pointer}}
.{prefix} .choices input:checked + label{{border-color:var(--dg-accent);background:var(--dg-accent-tint);font-weight:600}}
.{prefix} .choices input:focus-visible + label{{outline:3px solid var(--dg-accent);outline-offset:3px}}
.{prefix} .details{{border-top:1px solid var(--dg-rule);padding-top:12px;margin-top:8px;min-height:100px}}
.{prefix} .details p,.{prefix} .intro{{font-size:16px;line-height:1.8;margin:0 0 8px}}
.{prefix} .details h4{{font-family:inherit;font-weight:600;font-size:16px;line-height:1.5;margin:0 0 4px}}
.{prefix} [hidden]{{display:none}}
.{prefix} .choices input{{position:absolute;width:1px;height:1px;opacity:0}}
.{prefix} .legend{{font-size:14px;line-height:1.7;color:var(--dg-muted);margin:8px 0}}
@media(prefers-reduced-motion:reduce){{.{prefix} svg [data-detail-mark]{{transition:none}}}}
@media print{{.{prefix} .choices{{display:none}}.{prefix} [hidden]{{display:block}}.{prefix}.enhanced svg [data-detail-mark]{{opacity:1}}.{prefix} .details{{min-height:0}}}}
{extra}
</style>'''
 svg=f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 480 {height}" role="img" aria-labelledby="{prefix}-title {prefix}-desc"><title id="{prefix}-title">{escape(title)}</title><desc id="{prefix}-desc">{escape(desc)}</desc>{body}</svg>'
 rules=[]
 hover_rules=[]
 for k,_,_ in choices:
  rules.append(f'.{prefix}:has(#{prefix}-choice-{k}:checked) [data-detail]:not([data-detail="{k}"]){{display:none}}')
  hover_rules.extend([f'.{prefix}:has(label:hover) [data-detail]{{display:none!important}}',f'.{prefix}:has(label[for="{prefix}-choice-{k}"]:hover) [data-detail="{k}"]{{display:block!important}}'])
 rules.append(f'.{prefix} svg [data-detail-mark]{{opacity:1}}')
 rules.append('@media (hover:hover) and (pointer:fine){'+'\n'.join(hover_rules)+'}')
 styles=styles.replace('</style>','\n'.join(rules)+f'\n@media print{{.{prefix} [data-detail]{{display:block!important}}.{prefix} [data-detail-mark]{{opacity:1!important}}}}\n</style>')
 buttons=''.join(f'<span><input type="radio" id="{prefix}-choice-{k}" name="{prefix}-choice" value="{k}" aria-controls="{prefix}-panel-{k}" {"checked" if i==0 else ""}><label for="{prefix}-choice-{k}">{escape(label)}</label></span>' for i,(k,label,_) in enumerate(choices))
 panels=''.join(f'<section id="{prefix}-panel-{k}" data-detail="{k}" aria-label="{escape(label)}"><h4>{escape(label)}</h4><p>{escape(detail)}</p></section>' for k,label,detail in choices)
 frag=styles+f'<div class="{prefix}" data-teaching-widget="{name}">'+(f'<p class="intro">{escape(lead)}</p>' if lead else '')+svg+f'<div class="choices" role="group" aria-label="選擇圖中概念；也可用 Tab 與方向鍵">{buttons}</div><div class="details" aria-live="polite">{panels}</div><p class="legend">選擇按鈕可查看說明；滑鼠停留、鍵盤聚焦與觸控皆可操作。</p></div>'+'\n'
 (ROOT/f'{name}.html').write_text(frag)
 standalone=f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(title)}</title><style>:root{{{TOKENS}}}html[data-theme=dark]{{{DARK}}}@media(prefers-color-scheme:dark){{:root:not([data-theme=light]){{{DARK}}}}}body{{margin:0;padding:24px 16px;background:var(--dg-paper);color:var(--dg-ink);font-family:system-ui,sans-serif}}main{{max-width:760px;margin:auto}}h1{{font-family:'Instrument Serif','Noto Serif TC',serif;font-weight:400;font-size:28px}}.eyebrow{{color:var(--dg-muted);font-size:14px}}footer{{margin-top:24px;font-size:14px;color:var(--dg-muted)}}</style></head><body><main><p class="eyebrow">義大利預報實驗 · {kind}</p><h1>{escape(title)}</h1>{frag}<footer>可離線閱讀；字型沿用系統已安裝的字型。</footer></main></body></html>'''
 (ROOT/'standalone'/f'{name}.html').write_text(standalone)

# 1: topology from the accepted Archify source. A narrow vertical reading path
# replaces lane chrome in the embedded version; ordered edges carry the handoff.
a=json.loads((ROOT/'sources/d01_forecast_pipeline.workflow.json').read_text())
body=''.join(arrow(240,96+i*104) for i in range(4))
choices=[]
detail=["讀取 HORUS 的時間、震央、深度與規模；先決定哪些事件可收進實驗。","把收集區 S、測試區 R、規模門檻、發報窗與評分方式寫清楚，再開始比較。","用學習期決定參數；發報後，模型只能依事先約定的規則加入新歷史。","輸出的是各格、各規模箱、各預報窗的期望數，不是某顆地震必定發生的宣告。","用測試期真正發生的目標地震檢驗預報；看完結果再改模型，便需要新的測試資料。"]
for i,n in enumerate(a['nodes']):
 body+=group(n['id'],box(64,16+i*104,352,80,n['label'],n['sublabel'].split('：')[-1] if i==0 else ['','先把題目寫清楚','只用當時已知資訊','各格箱的期望數','與之後的事件對照'][i],i==3))
 choices.append((n['id'],n['label'],detail[i]))
make('d01_forecast_pipeline',a['meta']['title'],'从 HORUS 目錄依序訂規格、定模型、發出格內期望數預報，最後對照測試期事件。'.replace('从','從'),536,body,choices,kind='workflow')

# 9: a spatial cross-section makes the additivity visible without claiming a
# particular quantitative PPE kernel or a numerical Italy result.
body=path('M 48 232 H 440 M 48 232 V 40')+text(48,28,'率面的一條剖面','name')+text(240,272,'位置（示意）','sub','middle')
xs=list(range(48,441,4)); centers=[144,304]
def bell(x,c):return math.exp(-.5*((x-c)/52)**2)
for i,c in enumerate(centers):
 points=[(x,232-104*bell(x,c)) for x in xs]
 body+=group('source'+str(i+1),chart_line(points,'curve-other' if i==0 else 'curve-faint')+f'<circle class="dot" cx="{c}" cy="232" r="8"/>'+text(c,304,'震央 '+str(i+1),'sub','middle'))
body+=group('sum',chart_line([(x,232-104*(bell(x,144)+bell(x,304))) for x in xs]))
body+=text(240,344,'同一位置的貢獻可以相加','name','middle')
make('d09_kernel_smoothing','從過去震央，疊成一張率面','兩個過去震央各貢獻平滑空間核；同一位置的兩份貢獻相加，形成總率面。',372,body,[('source1','震央 1','核把一顆過去事件的貢獻分配到附近位置。距離越遠通常越小，並不只填入震央所在那一格。'),('source2','震央 2','第二顆事件也有自己的核。帶寬決定貢獻攤得多開；比較不同帶寬時，要同時注意核的正規化。'),('sum','相加後的率面','實線是兩條虛線逐點相加。真正預報還要把空間率積分到測試格，乘上規模分配與預報窗長度。')],lead='概念剖面：虛線是單顆事件的核，實線是總和。這裡不代入義大利 PPE 的參數。',kind='line concept')

# 12: bounded illustrative branching tree; every child has a single parent.
body=path('M 216 92 V 112 Q 216 120 208 120 H 120 Q 112 120 112 128 V 152')+path('M 264 92 V 112 Q 264 120 272 120 H 360 Q 368 120 368 128 V 152')+path('M 112 216 V 248 Q 112 256 120 256 H 232 Q 240 256 240 264 V 288')
body+=group('background',box(128,16,224,76,'背景事件','第 0 代',True))
body+=group('first',box(24,152,176,64,'直接後代 A')+box(280,152,176,64,'直接後代 B'))
body+=group('second',box(128,288,224,76,'A 的直接後代','根事件的第 2 代'))
body+=text(240,412,'連線表示觸發親子關係','sub','middle')
make('d12_etas_branching','ETAS：事件也可以生出事件','一顆背景事件產生兩個第一代後代，其中一個再產生第二代；分支比是平均直接後代數。',440,body,[('background','背景事件','背景機制產生第 0 代事件。這些事件可以觸發後代，也可能沒有後代。'),('first','第 1 代','這張示意樹有兩顆直接後代；真實模型的數量是隨機的，不能把這張樹的 2 當作分支比。'),('second','第 2 代','後代還能觸發自己的後代。分支比 n 指平均直接後代數；多代總數是另一個量。')],kind='tree')

# 15: three independent visual directions, with the mathematical operation
# stated explicitly. These shapes are conceptual, not calibrated forecasts.
body=box(96,16,288,72,'一顆輸入事件','時間、位置、規模')
for k,(y,label,sub) in enumerate([(136,'規模核 g(m)','未來規模'),(248,'時間核 f(t)','距輸入事件多久'),(360,'空間核 h(x,y)','距輸入震央多遠')]):
 visual=text(24,y+24,label,'name')+text(24,y+56,sub,'sub')
 if k<2:
  visual+=path(f'M 252 {y+72} H 452')
  pts=[]
  for j in range(49):
   u=j/48
   val=math.exp(-.5*((u-.5)/.18)**2) if k==0 else (0 if j==0 else math.exp(-.5*((math.log(u)+1.2)/.5)**2)/u/3.8)
   pts.append((252+200*u,y+72-64*val))
  visual+=chart_line(pts)
 else:
  visual+=''.join(f'<circle cx="352" cy="{y+40}" r="{r}" class="curve-other"/>' for r in [16,32,48])+f'<circle cx="352" cy="{y+40}" r="8" class="dot"/>'
 body+=group(['magnitude','time','space'][k],visual)
body+=text(240,480,'單顆貢獻 ∝ g × f × h','name','middle')+text(240,516,'乘上權重，再加總所有事件','sub','middle')
make('d15_eepas_kernels','EEPAS：一顆事件的三個核','输入事件的規模決定規模、時間、空間三核的尺度；三核相乘、加權後，再與其他事件和背景相加。'.replace('输入','輸入'),544,body,[('magnitude','規模核','輸入事件的規模，決定它對哪些未來規模較有貢獻。這是密度分配，不是把輸入規模直接當作下一次地震的規模。'),('time','時間核','時間核在正的時間差上分配貢獻。對數時間的中心與散布決定高貢獻大致出現在哪段時間。'),('space','空間核','空間核以輸入震央為中心，依其規模決定尺度。將三核乘積積分到某個窗、格與規模箱，才得到格箱期望數的貢獻。')],lead='圖中的乘積是核結構。完整預報還含產能／門檻修正、事件權重與背景混合。',kind='kernel concept')

# 16: testing flow is invariant; the selector changes conditioning assumptions.
body=box(80,16,320,72,'一張格箱預報','先選檢驗問題')+arrow(240,88)+box(48,128,384,76,'固定條件 → 生成模擬','各檢驗固定的量不同',True)+arrow(240,204)+box(48,244,384,76,'觀測與模擬用同一統計量','比較尾機率或分位數')
make('d16_csep_tests_flow','每一種檢驗，都有自己的模擬問題','N 檢驗總數；S、M 與 cL 固定觀測總數後模擬位置或規模；L 同時模擬總數與格箱分配。',348,body,[('n','N：數量','固定預报的總期望數，模擬總事件數；以觀測總數計算含等號的上下尾機率。Poisson 與負二項版本有不同的計數假設。'.replace('預报','預報')),('s','S：位置','固定觀測總數 N，將預報對規模箱加總，正規化成各空間格機率，再用多項分布模擬 N 顆事件的位置。'),('m','M：規模','固定觀測總數 N，將預報對空間格加總，正規化成各規模箱機率，再用多項分布模擬規模分配。'),('cl','cL：條件聯合','固定觀測總數 N，依各空間格與規模箱的聯合機率模擬分配；總數不變，但位置與規模配置會變。'),('l','L：完整聯合','先依預報的計數分布模擬總數，再分配到格箱；等價的 Poisson 版本可逐格箱抽計數。總數、位置、規模都會變。')],lead='先指定問題與模擬規則，再把觀測統計量放進模擬分布。',kind='workflow')

# 19: exact two-bin example. Lambda_A=(2,6), Lambda_B=(6,2), pi=.5.
body=text(232,20,'格 1','sub','middle')+text(344,20,'格 2','sub','middle')
for key,y,label,vals in [('a',24,'預報 A',(2,6)),('b',152,'預報 B',(6,2)),('mix',280,'各占一半',(4,4))]:
 g=text(32,y+24,label,'name')
 for i,value in enumerate(vals):
  x=200+i*112
  g+=f'<rect x="{x}" y="{y+72-value*8}" width="64" height="{value*8}" class="{"fill" if key=="mix" else "node"}"/>'+text(x+32,y+100,str(value),'name','middle')
 g+=path(f'M 192 {y+72} H 416')
 body+=group(key,g)
body+=text(240,420,'每張總數都是 8','name','middle')+text(240,456,'加總守恆 ≠ 概似守恆','sub','middle')
make('d19_convex_mix','凸組合：逐格相加，總期望數如何變','兩格預報 A 為 2 與 6，B 為 6 與 2；各占一半得到 4 與 4，總期望數仍為 8。',484,body,[('a','預報 A','兩個格箱的期望數為 2、6，總和為 8。這是合成例子，不是某個義大利預報窗的結果。'),('b','預報 B','相同兩個格箱的期望數為 6、2，總和也是 8。兩張預報必須使用相同區域、窗與規模箱，才能逐格混合。'),('mix','各占一半','混合結果為 4、4。一般而言，總數等於 πΛ_A + (1−π)Λ_B；只有兩個成分總数相同時，才對所有 π 保持同一總數。對數概似不服從這個守恆式。'.replace('總数','總數'))],lead='合成的兩格例子。柱高表示格內期望數；權重先定為 1/2。',kind='bar concept')

# 20: distinguish a fixed SUP baseline from the models using event history.
# Time-scale tendencies belong to model settings, not universal family bounds.
body=box(112,16,256,64,'格箱率預報',focal=True)
body+=path('M 216 80 V 104 Q 216 112 208 112 H 96 Q 88 112 88 120 V 144')+path('M 264 80 V 104 Q 264 112 272 112 H 336 Q 344 112 344 120 V 144')
body+=group('fixed',box(16,144,144,88,'SUP','固定基準'))
body+=box(224,144,240,88,'使用歷史資訊','依發報規則更新')
body+=path('M 280 232 V 252 Q 280 260 272 260 H 96 Q 88 260 88 268 V 328')+path('M 344 232 V 276 Q 344 284 336 284 H 248 Q 240 284 240 292 V 328')+path('M 392 232 V 328')
body+=group('ppe',box(16,328,144,80,'PPE','緩慢更新'))+group('etas',box(168,328,144,80,'ETAS','觸發活動'))+group('eepas',box(320,328,144,80,'EEPAS','前兆尺度'))
body+=text(240,456,'成分可組合，但權重須先定','name','middle')
make('d20_model_family','四種預報：差別在於如何使用歷史','SUP 是固定空間基準；本站 PPE 隨已知來源緩慢更新並使用 50 天延遲；ETAS 與 EEPAS 分別由觸發與前兆尺度描述歷史貢獻。',488,body,[('fixed','SUP','本站 SUP 在等面積格中均勻配置空間率，使用學習期決定的固定率。它是後面三個模型共同的簡單比較基準。'),('ppe','PPE','本站 PPE 的來源目錄逐窗更新，來源須早於發報時刻 50 天。固定的是參數，並非來源或空間權重；時間項也隨目錄累積長度緩慢改變。'),('etas','ETAS','ETAS 以觸發核描述事件後的叢集。本站格網版本只納入背景與已知歷史的第一代貢獻。變化的時間尺度依核參數、歷史與發報窗而定。'),('eepas','EEPAS','EEPAS 依前兆尺度把事件貢獻分配到規模、時間與空間。本站使用 50 天資料延遲；較長時間尺度是本實驗的設定與解讀，不是所有 EEPAS 的固定範圍。')],kind='tree')



# Original five figures rebuilt for narrow reading surfaces. Long text lives in
# HTML, not a scaled SVG. Full static semantics remain visible on every choice.
def stable_marks(name):
 return f'.dg-{name} svg [data-detail-mark]{{opacity:1!important}}'

# 2: this is a set diagram, deliberately not a geographic polygon.
name='d02_regions_s_r'
body=group('collection','<rect x="24" y="24" width="432" height="352" rx="8" class="node"/>'+text(48,60,'收集區 S','name'))
body+=group('testing','<rect x="72" y="120" width="224" height="192" rx="8" class="focal"/>'+text(96,156,'測試區 R','name'))
for x in [128,184,240]:body+=path(f'M {x} 180 V 288','grid')
for y in [180,216,252,288]:body+=path(f'M 88 {y} H 280','grid')
body+=group('edge','<circle cx="376" cy="216" r="8" class="dot"/>'+path('M 364 216 H 304 m 12 -8 l -12 8 l 12 8')+text(376,252,'R 外事件','sub','middle'))
body+=text(184,344,'只在 R 內評分','name','middle')+text(240,420,'收集歷史的範圍，要比評分範圍大','sub','middle')
make(name,'為什麼需要收集區與測試區？','集合示意：測試區 R 位於收集區 S 內；S 內但 R 外的歷史事件仍可能影響 R 邊緣，評分只計 R 內目標。',448,body,[('collection','S：收集歷史','收集區 S（collection region）提供符合模型條件的歷史事件。真實實驗使用 CPTI15 的涵蓋範圍；這個矩形只表示集合關係，不是假畫義大利地圖。'),('testing','R：預報與評分','測試區 R（testing region）是預報及評分的範圍。實際義大利實驗有 177 格；圖中的少量格線只示範格內計數，沒有聲稱格數或形狀與實際相同。'),('edge','邊界效應','R 外、S 內的事件仍可能對 R 內貢獻活動。只收 R 內事件會漏掉這部分已知歷史；收集區放大能減少這類邊界效應，但不保證已知所有區外來源。')],lead='集合關係示意，不是實際 Italy polygon。',kind='nested',extra=stable_marks(name))

# 3a: period cards are not a scaled axis; the issue-time inset below is scaled.
name='d03_timeline_windows'
body=group('periods',box(40,16,400,72,'暖機期 1960–1989','提供較早的歷史')+box(40,104,400,72,'學習期 1990–2011','估計或選定參數')+box(40,192,400,72,'測試期 2012–2021','依序發報，再對照觀測',True))
body+=group('windows',text(240,308,'40 個窗 × 91.31 天','name','middle')+text(240,340,'每窗固定自己的已知歷史','sub','middle'))
# x = 240 + 2.2 * days from issue; 50-day gap and 91.31-day window are proportional.
body+=group('cutoff',text(24,384,'一個發報窗的時間軸','name')+'<rect x="130" y="424" width="110" height="32" class="node"/><rect x="240" y="424" width="200.882" height="32" class="focal"/>'+path('M 24 440 H 456 M 130 404 V 468 M 240 400 V 468 M 440.882 416 V 468')+text(130,496,'−50 天','sub','middle')+text(240,496,'發報 t₀','name','middle')+text(432,496,'窗終點','sub','middle')+text(240,540,'PPE／EEPAS：只收截止前的來源','name','middle'))
make(name,'三段期間，與每次發報的資料截止','暖機提供歷史、學習期定參數、測試期發40個91.31天窗；PPE與EEPAS的來源截止在發報前50天，ETAS使用發報前歷史。',568,body,[('periods','1．三段期間','上方是期間用途卡，卡片寬度不表示年數。暖機期提供早期歷史，但小規模事件不完整；學習期決定模型參數，測試期不拿未來事件重選參數。'),('windows','2．逐窗發報','測試期共有 40 個互不重疊的 91.31 天窗，不是精確日曆季。最後一窗到 2021-12-31 09:36，年底餘下 0.6 天不評分。新事件可以進入後續窗的歷史，不能改寫先前預報。'),('cutoff','3．50 天延遲','下方時間軸的灰帶是發報前 50 天，藍帶是接下來 91.31 天。PPE 與 EEPAS 在發報時排除這段近期來源，整窗固定來源集合；ETAS 使用發報前歷史。SUP 使用學習期估出的固定率。')],kind='timeline',extra=stable_marks(name))

# 3b: explicitly compare bin centres and edges; never subtract unlike units.
name='d03_thresholds'
body=group('complete',box(24,16,432,80,'Mc：目錄能否收錄完整？','採用的年代門檻，需要文獻依據'))
body+=group('input',box(24,120,432,80,'m₀：哪些事件進入歷史？','名目 2.5 ↔ 箱界 2.45'))
body+=group('target',box(24,224,432,80,'mT：哪些事件會被評分？','名目 5.0 ↔ 箱界 4.95',True))
body+=group('bins',text(240,344,'輸出：25 個規模箱','name','middle'))
for i in range(25):body+=f'<rect x="{88+i*12}" y="368" width="12" height="32" class="fill"/>'
body+=text(88,432,'5.0','sub','middle')+text(388,432,'7.5','sub','middle')+text(240,472,'每箱寬 0.1；右端不包含','sub','middle')
make(name,'三個規模門檻，各回答不同問題','Mc是目錄完整度，m0是輸入門檻，mT是評分門檻；名目規模2.5與5.0分別對應原始箱界2.45與4.95，输出有25個寬0.1的規模箱。'.replace('输出','輸出'),500,body,[('complete','Mc：觀測條件','完整度規模 Mc 是指定地區與年代的收錄能力。1990 年起採用的整體文獻門檻是名目 2.5；轉成箱左界為 2.45。暖機期與局部區域仍可能需要更高門檻。'),('input','m₀：模型輸入','本實驗 m₀ = 2.45 是原始規模的有效邊界，對應取整後的名目 2.5。與 Mc 比較前，兩者必須都用箱心或都用箱界；數值相容不保證每一格、每一時段都完整。'),('target','mT：預報目標','目標用取整規模至少 5.0；原始規模 4.95 會取整成 5.0，因此也可能是目標。同一把尺上，目標與輸入相差 2.5 個規模單位，而不是將 5.0 和 2.45 直接相減。'),('bins','25 箱與差距','取整規模按 [5.0,5.1)、[5.1,5.2) 等輸出箱指派，直到 [7.4,7.5)。輸入下限比目標低至少兩級，是相關模型文獻的設計建議，不是完整度或預報有效性的普遍保證。')],kind='threshold concept',extra=stable_marks(name))

# 4: preserve the integral -> expected count -> distribution distinction.
name='d04_density_to_counts'
body=''.join(arrow(240,96+i*112) for i in range(3))
body+=group('density',box(48,16,384,80,'率密度 λ(t,x,y,m)','每單位時間、面積與規模'))
body+=group('integral',box(48,128,384,80,'在窗、格、箱上積分','累積範圍內的全部貢獻'))
body+=group('expected',box(48,240,384,80,'得到期望數 Λ','合成例子：Λ = 0.03',True))
body+=group('count',box(48,352,384,80,'再指定計數分布','例如 N ∼ Poisson(Λ)'))
body+=text(240,476,'期望數不會自動指定全部機率','name','middle')
make(name,'從率密度，到可以比較的顆數','率密度經過時間、空間與規模積分得到期望數；再加入Poisson等計數假設，才能計算實際顆數的機率。',504,body,[('density','1．率密度','率密度 λ 不是顆數，也不是機率。若使用天、平方公里與規模單位，單位就是顆／天／平方公里／規模單位。'),('integral','2．積分','積分範圍是一個 91.31 天窗、一個邊長 30√2 km 的方格及寬 0.1 的規模箱。數值積分可將小塊的密度乘體積後加總；只有密度為常數時，才等於單一密度值乘全部範圍。'),('expected','3．期望數','Λ = 0.03 是合成示例，表示重複相同條件時的平均顆數可以是小數。實際觀測仍是整數；沒有指定計數模型前，不能直接將 0.03 當作至少一顆的機率。'),('count','4．計數分布','若另採 Poisson 假設，Λ = 0.03 時 P(N=0)≈0.97045、P(N=1)≈0.02911、P(N≥2)≈0.00044。至少一顆的機率是 1−exp(−Λ)≈0.02955。其他計數分布可能給不同機率。')],kind='density concept',extra=stable_marks(name))

# 6: exact finite illustrative sample, so ties and both tails can be checked.
name='d06_quantile_pvalue'
body=text(24,32,'1．固定統計量與判斷方向','name')+text(24,72,'2．模擬，也計算同一統計量','name')+text(24,112,'3．看觀測值在模擬分布的位置','name')
counts=[1,2,4,6,4,2,1];base=336
body+=path(f'M 48 {base} H 440')
for i,c in enumerate(counts):
 x=64+i*52;body+=f'<rect x="{x}" y="{base-c*24}" width="36" height="{c*24}" class="{"fill" if i<=2 else "node"}"/>'+text(x+18,base-c*24-12,str(c),'sub','middle')+text(x+18,368,str(i),'sub','middle')
body+=path('M 186 160 V 340','curve')+text(194,152,'觀測值 = 2','name')+text(240,408,'20 次合成模擬；7 次不大於 2','sub','middle')+text(240,452,'q = 7 / 20 = 0.35','name','middle')
make(name,'分位數分數，怎麼接回 p 值？','20次合成模擬中有7次統計量不大於觀測值2，因此q為0.35；其中4次同值，上尾需包含這些同值，為17除以20。',480,body,[('statistic','先指定問題','先選統計量及「極端」的方向。例如總數太多是上尾問題；對數概似太低是下尾問題。不能看完觀測位置，再挑讓結果最顯著的一端。'),('quantile','q：排序比例','每根柱表示合成模擬次數，從統計量 0 到 6 依序為 1、2、4、6、4、2、1，合計 20。觀測值為 2；模擬值 ≤ 2 有 7 次，所以 q = 0.35。'),('tails','同值也要算進尾端','此例與觀測值同為 2 的模擬有 4 次。下尾比例是 7/20；上尾比例是 1−7/20+4/20=17/20。直接用 1−q 會漏掉同值。'),('decision','q 不總是 p 值','q 是觀測值的排序位置；只有事先指定下尾檢定時，它才近似該下尾 p 值。上尾及雙尾需採各自規則。模擬比例有抽樣誤差，q 小也不是模型為假的機率。')],lead='下面是可逐顆數出的合成例子，不是義大利的檢驗分數。',kind='quantile concept',extra=stable_marks(name))

print(f'Generated {WRITTEN} embedded figures and {WRITTEN} offline previews.')
