"""Exercise reduced motion, offscreen pausing, and the real navigation controls."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reference/notes/refresh_20260913_exhibition'
PAGES={'d01_forecast_pipeline':'01_forecast_question','d02_regions_s_r':'02_reading_horus','d03_timeline_windows':'03_experiment_spec','d03_forecast_protocols':'03_experiment_spec','d03_thresholds':'03_experiment_spec','d04_density_to_counts':'04_poisson_and_sup','d06_quantile_pvalue':'06_simulation_tests_scores','d09_kernel_smoothing':'09_ppe_forecast','d11_hawkes_process':'11_conditional_intensity','d12_etas_branching':'12_etas_structure','d15_eepas_kernels':'15_eepas_italy_forecast','d16_csep_tests_flow':'16_test_number','d19_convex_mix':'19_ensembles','d20_model_family':'20_beyond_forecast'}
rows=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
    context=browser.new_context(viewport={'width':390,'height':844},has_touch=True,reduced_motion='reduce',color_scheme='dark');page=context.new_page()
    for name,stem in PAGES.items():
        page.goto((ROOT/'book/_build/html'/f'{stem}.html').as_uri(),wait_until='load')
        page.emulate_media(reduced_motion='reduce')
        exhibit=page.locator(f'.quake-exhibit[data-exhibit="{name}"]')
        page.wait_for_function('(name)=>document.querySelector(`[data-exhibit="${name}"]`).dataset.mounted==="true"',arg=name)
        slider=exhibit.locator('[data-control=progress]');slider.scroll_into_view_if_needed();slider.focus();page.keyboard.press('Home');page.keyboard.press('ArrowRight');assert float(slider.input_value())>0
        controls=exhibit.locator('[data-action=play]');offscreen=None
        if controls.count():
            controls.tap();before=exhibit.get_attribute('data-progress');page.wait_for_timeout(250)
            assert exhibit.get_attribute('data-playing')=='false' and exhibit.get_attribute('data-progress')==before,(name,'reduced motion must not run')
            page.emulate_media(reduced_motion='no-preference');controls.click();page.wait_for_timeout(80)
            assert exhibit.get_attribute('data-playing')=='true',(name,'play')
            page.evaluate('scrollTo(0,0)')
            visible=exhibit.evaluate('e=>{const r=e.getBoundingClientRect();return r.top<innerHeight&&r.bottom>0}')
            if visible:page.evaluate('scrollTo(0,document.documentElement.scrollHeight)')
            page.wait_for_timeout(180);offscreen=exhibit.get_attribute('data-playing')=='false';assert offscreen,(name,'offscreen')
        row={'diagram':name,'reduced_motion_keyboard':True,'animated':bool(controls.count()),'offscreen_pauses':offscreen};rows.append(row);print(json.dumps(row),flush=True)
    page.goto((ROOT/'book/_build/html/00_intro.html').as_uri(),wait_until='load')
    data=page.evaluate('JSON.parse(document.querySelector("#museum-atlas-data").textContent)')
    slider=page.locator('[data-atlas-range]');slider.evaluate('e=>{e.value="2012";e.dispatchEvent(new Event("input",{bubbles:true}))}')
    assert page.locator('.atlas-event').count()==sum(e[2]<=2012 for e in data['targets'])
    slider.evaluate('e=>{e.value="2021";e.dispatchEvent(new Event("input",{bubbles:true}))}')
    assert page.locator('.atlas-event').count()==25
    page.locator('.museum-search').click();assert page.locator('#pst-search-dialog').is_visible();page.keyboard.press('Escape')
    page.get_by_role('button',name='開啟章節目錄',exact=True).click();assert page.locator('.museum-menu a').count()==36;page.keyboard.press('Escape')
    seen=set()
    for _ in range(3):
        page.get_by_role('button',name='切換深淺主題',exact=True).click();page.wait_for_timeout(80);seen.add(page.evaluate('document.documentElement.dataset.theme'))
    assert seen=={'light','dark'},seen
    browser.close()
report={'passed':True,'exhibits':rows,'atlas_data_filter':True,'search_dialog':True,'chapter_menu':True,'theme_button':True};(OUT/'lifecycle_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print('LIFECYCLE PASSED')
