"""Browser acceptance of exhibition content, rendering, navigation and real controls."""
from pathlib import Path
import argparse,hashlib,json,os,re
import yaml
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
HTML=ROOT/'book/_build/html'
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/refresh_20260913_exhibition/browser')
DEVICES={'desktop-light':(1440,1000,'light'),'mobile-light':(390,844,'light'),'desktop-dark':(1440,1000,'dark'),'mobile-dark':(390,844,'dark')}


def main():
    toc=yaml.safe_load((ROOT/'book/_toc.yml').read_text());allpages=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']]
    parser=argparse.ArgumentParser();parser.add_argument('--pages',nargs='+',choices=allpages);parser.add_argument('--devices',nargs='+',choices=DEVICES,default=list(DEVICES));parser.add_argument('--preview',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);results=[];failures=[];seen=set()
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
        for device in args.devices:
            width,height,theme=DEVICES[device];context=browser.new_context(viewport={'width':width,'height':height},has_touch=width<600,color_scheme=theme,device_scale_factor=1)
            page=context.new_page()
            for stem in args.pages or allpages:
                errors=[];capture=lambda error:errors.append(str(error));page.on('pageerror',capture)
                page.goto((HTML/f'{stem}.html').as_uri(),wait_until='load',timeout=60000)
                page.evaluate('theme=>{document.documentElement.dataset.theme=theme;document.documentElement.dataset.mode=theme}',theme)
                page.wait_for_function("[...document.querySelectorAll('.quake-exhibit')].every(e=>e.dataset.mounted==='true')",timeout=15000)
                page.evaluate('async()=>{if(window.MathJax?.startup?.promise)await MathJax.startup.promise;}')
                page.wait_for_function("()=>{const c=getComputedStyle(document.documentElement).getPropertyValue('--tc-paper').trim();return [...document.querySelectorAll('.js-plotly-plot')].every(e=>e._fullLayout?.paper_bgcolor===c)}",timeout=15000)
                page.wait_for_timeout(180)
                metrics=page.evaluate(r'''()=>{
                    const article=document.querySelector('.bd-article');const text=article.innerText;
                    const walker=document.createTreeWalker(article,NodeFilter.SHOW_TEXT);let n,bold=[];
                    while(n=walker.nextNode()){const e=n.parentElement;if(!e.closest('code,pre,script,style,.math,mjx-container')&&e.getClientRects().length&&n.textContent.includes('**'))bold.push(n.textContent.trim().slice(0,150));}
                    return {viewport:innerWidth,width:document.documentElement.scrollWidth,heading:article.querySelector('h1')?.innerText,
                      unparsedBold:bold,warnings:/Matplotlib created|not a\s+writable directory|Traceback \(most recent|RuntimeWarning:|UserWarning:/.test(text),
                      weakStrong:[...article.querySelectorAll('strong')].filter(e=>parseInt(getComputedStyle(e).fontWeight)<600).length,
                      mathErrors:document.querySelectorAll('mjx-merror').length,visibleCode:[...document.querySelectorAll('.cell_input')].filter(e=>e.getBoundingClientRect().height>0).length,
                      brokenImages:[...document.images].filter(e=>!e.complete||e.naturalWidth===0).map(e=>e.src),
                      plots:[...document.querySelectorAll('.js-plotly-plot')].length,exhibits:document.querySelectorAll('.quake-exhibit').length};}''')
                row={'page':stem,'device':device,**metrics,'errors':errors,'html_sha256':hashlib.sha256((HTML/f'{stem}.html').read_bytes()).hexdigest(),'exhibit_checks':[]}
                if width>=960:
                    sidebar=page.locator('.bd-sidebar-primary')
                    row['menu_opens']=sidebar.is_visible() and sidebar.locator('.bd-sidenav a').count()>=37
                    row['menu_closes']=not page.get_by_role('button',name='開啟章節目錄',exact=True).is_visible()
                    current=sidebar.locator('a[aria-current=page]')
                    row['sidebar_current_visible']=bool(current.count()) and current.evaluate('e=>{const a=e.getBoundingClientRect(),b=e.closest(".bd-sidebar-primary").getBoundingClientRect();return a.top>=b.top&&a.bottom<=b.bottom}')
                else:
                    page.get_by_role('button',name='開啟章節目錄',exact=True).click();row['menu_opens']=page.locator('.museum-menu').is_visible();page.keyboard.press('Escape');row['menu_closes']=not page.locator('.museum-menu').is_visible()
                exhibits=page.locator('.quake-exhibit')
                for i in range(exhibits.count()):
                    exhibit=exhibits.nth(i);name=exhibit.get_attribute('data-exhibit');kind=exhibit.get_attribute('data-kind');slider=exhibit.locator('[data-control=progress]')
                    slider.scroll_into_view_if_needed();slider.focus();page.keyboard.press('Home');page.keyboard.press('ArrowRight');keyboard=float(slider.input_value())>0
                    before=exhibit.locator('.ex-stage').inner_html();slider.evaluate("e=>{e.value='80';e.dispatchEvent(new Event('input',{bubbles:true}))}")
                    changed=exhibit.locator('.ex-stage').inner_html()!=before
                    select=exhibit.locator('select');compare=True
                    if select.count():
                        options=select.locator('option');last=options.last.get_attribute('value') or options.last.inner_text();chosen=options.first if select.input_value()==last else options.last;value=chosen.get_attribute('value') or chosen.inner_text();before_compare=exhibit.locator('.ex-stage').inner_html();select.select_option(value)
                        compare=select.input_value()==value and exhibit.locator('.ex-stage').inner_html()!=before_compare
                    layers=exhibit.locator('[data-layer]');layers_ok=True
                    for j in range(layers.count()):
                        control=layers.nth(j);prior=control.is_checked();before_layer=exhibit.locator('.ex-stage').inner_html();control.click();layers_ok&=control.is_checked()!=prior and exhibit.locator('.ex-stage').inner_html()!=before_layer;control.click()
                    supports_play=page.evaluate('kind=>EarthquakeExhibits.kinds[kind].play',kind)
                    playing=stopped=True
                    if supports_play:
                        play=exhibit.locator('[data-action=play]');play.click();page.wait_for_timeout(160);playing=exhibit.get_attribute('data-playing')=='true';play.click()
                        stopped=exhibit.get_attribute('data-playing')=='false'
                    start_value=page.evaluate('kind=>EarthquakeExhibits.kinds[kind].start*100',kind)
                    slider.evaluate('(e,v)=>{e.value=v;e.dispatchEvent(new Event("input",{bubbles:true}))}',start_value)
                    svg=exhibit.locator('.ex-stage svg');geometry=svg.evaluate('''e=>({width:e.getBoundingClientRect().width,minText:Math.min(...[...e.querySelectorAll('text')].map(t=>parseFloat(getComputedStyle(t).fontSize)*e.getBoundingClientRect().width/e.viewBox.baseVal.width))})''')
                    check={'name':name,'kind':kind,'keyboard':keyboard,'scrub_changes_scene':changed,'compare':compare,'layers':layers_ok,'play':playing,'pause':stopped,'supports_play':supports_play,**geometry}
                    row['exhibit_checks'].append(check)
                    if (device,name) not in seen:
                        exhibit.evaluate('e=>{document.activeElement?.blur();scrollTo(0,e.getBoundingClientRect().top+scrollY-100)}')
                        page.wait_for_timeout(180)
                        clip=exhibit.evaluate('e=>{const r=e.getBoundingClientRect();return {x:r.x+scrollX,y:r.y+scrollY,width:r.width,height:r.height}}')
                        page.screenshot(path=str(OUT/f'{device}-{name}.png'),clip=clip,full_page=True);seen.add((device,name))
                page.evaluate('scrollTo(0,0)');page.screenshot(path=str(OUT/f'{device}-{stem}-top.png'));page.screenshot(path=str(OUT/f'{device}-{stem}-full.png'),full_page=True)
                results.append(row);(OUT/'results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n');print(json.dumps(row,ensure_ascii=False),flush=True)
                failed=(row['width']>width+2 or errors or row['mathErrors'] or row['visibleCode'] or row['brokenImages'] or not row['menu_opens'] or not row['menu_closes'] or row.get('sidebar_current_visible') is False or any(not c[k] for c in row['exhibit_checks'] for k in ['keyboard','scrub_changes_scene','compare','layers','play','pause']))
                if not args.preview:failed|=bool(row['warnings'] or row['unparsedBold'] or row['weakStrong'])
                if failed:failures.append(row)
                page.remove_listener('pageerror',capture)
            context.close()
        browser.close()
    (OUT/'summary.json').write_text(json.dumps({'passed':not failures,'combinations':len(results),'failures':[{'page':r['page'],'device':r['device']} for r in failures]},indent=2)+'\n')
    print('FAILURES',json.dumps(failures,ensure_ascii=False),flush=True)
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
