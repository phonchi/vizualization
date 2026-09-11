"""Browser acceptance of the generated local teaching book (no publication)."""
from pathlib import Path
import json
import time
import yaml
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
HTML=ROOT/'book/_build/html'
OUT=ROOT/'reference/notes/rewrite_20260911/browser'
OUT.mkdir(parents=True,exist_ok=True)
toc=yaml.safe_load((ROOT/'book/_toc.yml').read_text())
pages=[toc['root']]+[c['file'] for part in toc['parts'] for c in part['chapters']]
representative={'00_intro','foundation_catalog','foundation_randomness','foundation_inference',
                '10_point_process','13_etas_structure','15_psi_phenomenon','18_testing_comparison',
                '01_overview','07_case_hualien2024','appendix_a_point_process','appendix_d_eepas'}
results=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
    for size,width,height in [('desktop',1440,900),('mobile',390,844)]:
        context=browser.new_context(viewport={'width':width,'height':height},device_scale_factor=1)
        page=context.new_page()
        for stem in pages:
            errors=[]
            def capture(error):errors.append(str(error))
            page.on('pageerror',capture)
            page.goto((HTML/f'{stem}.html').as_uri(),wait_until='load',timeout=60000)
            page.wait_for_timeout(250)
            page.evaluate('''async () => {
                if (window.MathJax && MathJax.startup && MathJax.startup.promise) {
                    await MathJax.startup.promise;
                }
            }''')
            metrics=page.evaluate('''() => ({
                heading: document.querySelector('h1')?.textContent,
                viewport: innerWidth,
                width: document.documentElement.scrollWidth,
                visibleCode: [...document.querySelectorAll('.cell_input')].filter(e=>e.getBoundingClientRect().height>0).length,
                plots: document.querySelectorAll('.js-plotly-plot').length,
                mathErrors: document.querySelectorAll('mjx-merror').length,
                mathRendered: document.querySelectorAll('mjx-container, .MathJax').length,
                rawMath: document.querySelectorAll('.math').length,
                dropdowns: document.querySelectorAll('.admonition.dropdown,details.sd-dropdown,details.toggle-details').length,
                initiallyOpen: document.querySelectorAll('.admonition.dropdown:not(.toggle-hidden),details.sd-dropdown[open],details.toggle-details[open]').length,
            })''')
            row={'page':stem,'device':size,**metrics,'errors':errors}
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
                plot=page.locator('.js-plotly-plot').first
                if plot.count():
                    plot.evaluate('(e)=>window.scrollTo(0,e.getBoundingClientRect().top+window.scrollY-80)')
                    page.wait_for_timeout(100)
                    page.screenshot(path=str(OUT/f'{size}-{stem}-plot.png'))
                    row['plot_trace_count']=plot.evaluate('(e)=>e.data?.length || 0')
                    if size=='mobile':
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
failures=[r for r in results if r['width']>r['viewport']+2 or r['visibleCode'] or r['mathErrors'] or r['errors'] or r.get('dropdown_opens') is False or r.get('dropdown_closes') is False or r.get('admonition_opens') is False or r.get('admonition_closes') is False or r.get('figure_scrolls') is False]
print('FAILURES',json.dumps(failures,ensure_ascii=False),flush=True)
if failures:raise SystemExit(1)
