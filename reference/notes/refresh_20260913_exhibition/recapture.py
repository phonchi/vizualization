import json
from pathlib import Path
from playwright.sync_api import sync_playwright
root=Path('/home/phonchi/vizualization');out=root/'reference/notes/refresh_20260913_exhibition/browser_final'
rows=json.loads((out/'results.json').read_text());seen=set();checks=[]
with sync_playwright() as p:
 b=p.chromium.launch(args=['--allow-file-access-from-files'])
 for device in ['desktop-light','mobile-light','desktop-dark','mobile-dark']:
  w=390 if 'mobile' in device else 1440;c=b.new_context(viewport={'width':w,'height':844 if w==390 else 1000},color_scheme=device.split('-')[1]);page=c.new_page()
  for row in [r for r in rows if r['device']==device and r['exhibit_checks']]:
   page.goto((root/'book/_build/html'/f"{row['page']}.html").as_uri());page.evaluate('(t)=>{document.documentElement.dataset.theme=t;document.documentElement.dataset.mode=t}',device.split('-')[1]);page.wait_for_timeout(300)
   for check in row['exhibit_checks']:
    name=check['name']
    if (device,name) in seen:continue
    e=page.locator(f'[data-exhibit="{name}"]');e.evaluate('e=>{document.activeElement?.blur();scrollTo(0,e.getBoundingClientRect().top+scrollY-100)}');page.wait_for_timeout(180)
    clip=e.evaluate('e=>{const r=e.getBoundingClientRect();return {x:r.x+scrollX,y:r.y+scrollY,width:r.width,height:r.height}}')
    skip=page.locator('.skip-link').evaluate('e=>e.getBoundingClientRect().bottom');assert skip<=0
    page.screenshot(path=str(out/f'{device}-{name}.png'),clip=clip,full_page=True);seen.add((device,name));checks.append({'device':device,'name':name,'skip_link_hidden_without_focus':True,'capture':'document clip below fixed navigation'})
  c.close()
 b.close()
(out/'capture_checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(len(checks),'recaptured')
