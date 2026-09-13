"""Browser acceptance of the generated local teaching book (no publication)."""
from pathlib import Path
import json
import argparse
import hashlib
import os
import time
import yaml
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
HTML=ROOT/'book/_build/html'
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/rewrite_20260913')/'browser'
OUT.mkdir(parents=True,exist_ok=True)
toc=yaml.safe_load((ROOT/'book/_toc.yml').read_text())
pages=[toc['root']]+[c['file'] for part in toc['parts'] for c in part['chapters']]
parser=argparse.ArgumentParser()
parser.add_argument('--pages',nargs='+',choices=pages)
args=parser.parse_args()
if args.pages:pages=args.pages
representative=set(pages)
seen_diagrams=set()
results=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
    for size,width,height in [('desktop-light',1440,900),('mobile-light',390,844),('desktop-dark',1440,900),('mobile-dark',390,844)]:
        context=browser.new_context(viewport={'width':width,'height':height},device_scale_factor=1)
        page=context.new_page()
        for stem in pages:
            errors=[]
            def capture(error):errors.append(str(error))
            page.on('pageerror',capture)
            page.goto((HTML/f'{stem}.html').as_uri(),wait_until='load',timeout=60000)
            page.evaluate("theme => {document.documentElement.dataset.theme=theme; document.documentElement.dataset.mode=theme;}",size.split('-')[-1])
            page.wait_for_timeout(500)
            page.evaluate('''async () => {
                if (window.MathJax && MathJax.startup && MathJax.startup.promise) {
                    await MathJax.startup.promise;
                }
            }''')
            page.wait_for_function('''() => {const c=document.documentElement.dataset.theme==='dark'?'#1e232b':'#fbfcfd';return [...document.querySelectorAll('.js-plotly-plot')].every(e=>e._fullLayout?.paper_bgcolor===c);}''',timeout=15000)
            metrics=page.evaluate('''() => ({
                heading: document.querySelector('h1')?.textContent,
                viewport: innerWidth,
                width: document.documentElement.scrollWidth,
                visibleCode: [...document.querySelectorAll('.cell_input')].filter(e=>e.getBoundingClientRect().height>0).length,
                plots: document.querySelectorAll('.js-plotly-plot').length,
                plotsReady: [...document.querySelectorAll('.js-plotly-plot')].filter(e=>e._fullLayout).length,
                plotBackgrounds: [...document.querySelectorAll('.js-plotly-plot')].map(e=>e._fullLayout?.paper_bgcolor),
                missingDiagrams: document.querySelectorAll('.diagram-missing').length,
                brokenImages: [...document.images].filter(e=>!e.complete || e.naturalWidth===0).map(e=>e.src),
                mathErrors: document.querySelectorAll('mjx-merror').length,
                mathRendered: document.querySelectorAll('mjx-container, .MathJax').length,
                rawMath: document.querySelectorAll('.math').length,
                dropdowns: document.querySelectorAll('.admonition.dropdown,details.sd-dropdown,details.toggle-details').length,
                initiallyOpen: document.querySelectorAll('.admonition.dropdown:not(.toggle-hidden),details.sd-dropdown[open],details.toggle-details[open]').length,
            })''')
            row={'page':stem,'device':size,**metrics,'errors':errors,
                 'html_sha256':hashlib.sha256((HTML/f'{stem}.html').read_bytes()).hexdigest()}
            expected_background='#1e232b' if size.endswith('dark') else '#fbfcfd'
            row['plotThemeCorrect']=all(c==expected_background for c in row['plotBackgrounds'])
            row['diagramChecks']=[]
            diagrams=page.locator('.teaching-diagram')
            for index in range(diagrams.count()):
                diagram=diagrams.nth(index)
                name=diagram.get_attribute('data-diagram')
                check={'name':name,'radioWorks':True,'buttonsWork':True}
                check['minimumTextPixels']=diagram.evaluate('''el=>{const svg=el.querySelector('svg'); const scale=svg.getBoundingClientRect().width/svg.viewBox.baseVal.width; return Math.min(...[...svg.querySelectorAll('text')].map(t=>parseFloat(getComputedStyle(t).fontSize)*scale));}''')
                radios=diagram.locator('input[type=radio]')
                if radios.count()>1:
                    radio=radios.nth(1)
                    identifier=radio.get_attribute('id')
                    diagram.locator(f'label[for="{identifier}"]').click()
                    check['radioWorks']=radio.is_checked()
                    radio.focus()
                    page.keyboard.press('ArrowLeft')
                    check['keyboardWorks']=radios.first.is_checked()
                    # Restore the first explanation before capturing.
                    first_id=radios.first.get_attribute('id')
                    diagram.locator(f'label[for="{first_id}"]').click()
                buttons=diagram.locator('button')
                for j in range(buttons.count()):
                    button=buttons.nth(j)
                    if button.is_visible() and button.is_enabled():
                        button.click()
                key=(size,name)
                if key not in seen_diagrams:
                    diagram.screenshot(path=str(OUT/f'{size}-{name}.png'))
                    seen_diagrams.add(key)
                row['diagramChecks'].append(check)
            admonitions=page.locator('.admonition.dropdown')
            if admonitions.count():
                first=admonitions.first
                first.locator('.admonition-title').click()
                row['admonition_opens']=not first.evaluate('(e)=>e.classList.contains("toggle-hidden")')
                first.locator('.admonition-title').click()
                row['admonition_closes']=first.evaluate('(e)=>e.classList.contains("toggle-hidden")')
            # MapLibre attribution uses <details open>; it is not a lesson
            # supplement and its legally required attribution must remain intact.
            details=page.locator('details.sd-dropdown,details.toggle-details')
            if details.count():
                first=details.first
                first.locator('summary').click()
                row['dropdown_opens']=first.evaluate('(e)=>e.open')
                first.locator('summary').click()
                row['dropdown_closes']=not first.evaluate('(e)=>e.open')
            if stem in representative:
                page.evaluate('window.scrollTo(0,0)')
                page.screenshot(path=str(OUT/f'{size}-{stem}-top.png'))
                page.screenshot(path=str(OUT/f'{size}-{stem}-full.png'),full_page=True)
                plot=page.locator('.js-plotly-plot').first
                if plot.count():
                    plot.evaluate('(e)=>window.scrollTo(0,e.getBoundingClientRect().top+window.scrollY-80)')
                    page.wait_for_timeout(100)
                    page.screenshot(path=str(OUT/f'{size}-{stem}-plot.png'))
                    row['plot_trace_count']=plot.evaluate('(e)=>e.data?.length || 0')
                    if size.startswith('mobile'):
                        row['figure_scrolls']=plot.evaluate('''(e)=>{
                            // MyST-NB may put overflow on the inner output.
                            let container=e;
                            while(container && container.tagName!=='BODY') {
                                if(container.scrollWidth>container.clientWidth+1 &&
                                   ['auto','scroll'].includes(getComputedStyle(container).overflowX)) {
                                    container.scrollLeft=container.scrollWidth;
                                    return container.scrollLeft>0;
                                }
                                container=container.parentElement;
                            }
                            return false;
                        }''')
                        page.screenshot(path=str(OUT/f'{size}-{stem}-plot-right.png'))
                        plot.evaluate("(e)=>{while(e && e.tagName!=='BODY'){e.scrollLeft=0;e=e.parentElement;}}")
            results.append(row)
            print(json.dumps(row,ensure_ascii=False),flush=True)
            page.remove_listener('pageerror',capture)
            (OUT/'results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
        context.close()
    browser.close()
failures=[r for r in results if r['width']>r['viewport']+2 or not r['plotThemeCorrect'] or r['missingDiagrams'] or r['brokenImages'] or r['plotsReady']!=r['plots'] or any(not d.get('radioWorks',True) or not d.get('keyboardWorks',True) or d.get('minimumTextPixels',12)<12 for d in r['diagramChecks']) or r['visibleCode'] or r['mathErrors'] or r['errors'] or r.get('dropdown_opens') is False or r.get('dropdown_closes') is False or r.get('admonition_opens') is False or r.get('admonition_closes') is False or r.get('figure_scrolls') is False]
print('FAILURES',json.dumps(failures,ensure_ascii=False),flush=True)
if failures:raise SystemExit(1)
