from playwright.sync_api import sync_playwright
from pathlib import Path
import json
with sync_playwright() as p:
 b=p.chromium.launch(headless=True,args=['--allow-file-access-from-files']);page=b.new_page(viewport={'width':390,'height':844})
 page.goto(Path('book/_build/html/11_conditional_intensity.html').resolve().as_uri(),wait_until='load')
 print(json.dumps(page.evaluate('''()=>['.bd-container','.bd-main','.bd-content','.bd-article-container','.bd-article','.museum-shell','.museum-chapter-header','.museum-section','.cell','.cell_output','.output.text_html','.teaching-diagram','.quake-exhibit','.ex-stage'].map(s=>{const e=document.querySelector(s),c=getComputedStyle(e),r=e.getBoundingClientRect();return{selector:s,width:r.width,x:r.x,min:c.minWidth,max:c.maxWidth,cssWidth:c.width,flex:c.flex,display:c.display,overflow:c.overflowX}})'''),indent=2))
 b.close()
