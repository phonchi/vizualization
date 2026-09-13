"""Measure SVG text bounds in the real browser at narrow and wide exhibit sizes."""
from pathlib import Path
from bs4 import BeautifulSoup
import json
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reference/notes/refresh_20260913_exhibition'
configs=[]
for path in sorted((ROOT/'book/_static/diagrams').glob('d*.html')):
    soup=BeautifulSoup(path.read_text(),'html.parser');node=soup.select_one('.quake-exhibit')
    if node:configs.append(json.loads(node['data-config']))
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
    page=browser.new_page(viewport={'width':1440,'height':1000})
    page.goto((ROOT/'book/_build/html/03_experiment_spec.html').as_uri(),wait_until='load')
    result=page.evaluate('''async configs=>{
      await document.fonts.ready;
      const root=document.createElement('div');root.className='quake-exhibit';root.dataset.mounted='true';root.style.cssText='position:absolute;left:-5000px;top:0';document.body.append(root);
      const rows=[];
      for(const cfg of configs)for(const width of [354,1100])for(const progress of [0,.55,1]){
        const def=EarthquakeExhibits.kinds[cfg.kind];const state={title:cfg.title,variant:def.select?.at(-1),S:true,R:true,events:true,contours:true,sources:true};
        const r=EarthquakeExhibits.svgHTML(cfg.kind,progress,state,cfg.data||{},width,'geometry-probe');
        root.style.width=width+'px';root.innerHTML=r.svg;const svg=root.querySelector('svg');
        const outside=[...svg.querySelectorAll('text')].map(e=>{const b=e.getBBox();return{text:e.textContent,x:b.x,y:b.y,right:b.x+b.width,bottom:b.y+b.height}}).filter(b=>b.x < -1||b.y < -1||b.right>width+1||b.bottom>r.height+1);
        rows.push({diagram:cfg.slug,width,progress,outside});
      }
      root.remove();return rows;
    }''',configs)
    browser.close()
report={'passed':all(not r['outside'] for r in result),'states':len(result),'results':result};(OUT/'svg_text_geometry.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'passed':report['passed'],'states':len(result),'failures':[r for r in result if r['outside']]},ensure_ascii=False))
if not report['passed']:raise SystemExit(1)
