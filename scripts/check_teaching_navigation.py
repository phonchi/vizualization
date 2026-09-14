"""Check persistent chapter navigation, resizing, and narrow reading layouts."""
from pathlib import Path
import json, os
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/os.environ.get('TEACHING_REPORT_DIR','reference/notes/reader_review_20260914/navigation')
OUT.mkdir(parents=True,exist_ok=True)
rows=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
    page=browser.new_page(viewport={'width':1024,'height':900})
    for width in [959,960,1024,1280,1440]:
        page.set_viewport_size({'width':width,'height':900})
        for stem in ['00_intro','11_conditional_intensity','appendix_g_observations']:
            page.goto((ROOT/'book/_build/html'/f'{stem}.html').as_uri(),wait_until='load')
            page.wait_for_timeout(300)
            menu=page.get_by_role('button',name='開啟章節目錄',exact=True)
            sidebar=page.locator('.bd-sidebar-primary')
            if width>=960:
                assert sidebar.is_visible() and not menu.is_visible()
                active=sidebar.locator('[aria-current=page]');assert active.count()==1
                assert active.evaluate('e=>{const r=e.getBoundingClientRect(),s=e.closest(".bd-sidebar-primary").getBoundingClientRect();return r.top>=s.top&&r.bottom<=s.bottom}')
            else:
                assert not sidebar.is_visible() and menu.is_visible()
                menu.click();assert page.locator('.museum-menu').is_visible()
                page.keyboard.press('Escape');assert not page.locator('.museum-menu').is_visible()
            assert page.evaluate('document.documentElement.scrollWidth')<=width+2
            page.screenshot(path=str(OUT/f'{width}-{stem}.png'))
            rows.append({'width':width,'page':stem,'passed':True})
    sidebar.locator('a[href="20_beyond_forecast.html"]').click()
    page.wait_for_url('**/20_beyond_forecast.html')
    assert sidebar.locator('[aria-current=page]').inner_text().startswith('20.')
    page.set_viewport_size({'width':390,'height':844})
    page.get_by_role('button',name='開啟章節目錄',exact=True).click()
    assert page.locator('.museum-menu').is_visible()
    page.set_viewport_size({'width':1024,'height':900});page.wait_for_timeout(200)
    assert not page.locator('.museum-menu').is_visible() and sidebar.is_visible()
    browser.close()
(OUT/'results.json').write_text(json.dumps({'passed':True,'layouts':rows,'chapter_link':True,'resize_closes_modal':True},indent=2)+'\n')
print('NAVIGATION PASSED',len(rows))
