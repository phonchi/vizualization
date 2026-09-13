"""Touch, hover, keyboard and reduced-motion checks for all twelve embedded diagrams."""
from pathlib import Path
from hashlib import sha256
import json
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
HTML=ROOT/'book/_build/html'
OUT=ROOT/'reference/notes/rewrite_20260913/continuation/interactions'
PAGES={'d01_forecast_pipeline':'01_forecast_question','d02_regions_s_r':'02_reading_horus',
       'd03_timeline_windows':'03_experiment_spec','d03_thresholds':'03_experiment_spec',
       'd04_density_to_counts':'04_poisson_and_sup','d06_quantile_pvalue':'06_simulation_tests_scores',
       'd09_kernel_smoothing':'09_ppe_forecast','d12_etas_branching':'12_etas_structure',
       'd15_eepas_kernels':'15_eepas_italy_forecast','d16_csep_tests_flow':'16_test_number',
       'd19_convex_mix':'19_ensembles','d20_model_family':'20_beyond_forecast'}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,args=['--allow-file-access-from-files'])
        for theme in ['light','dark']:
            context=browser.new_context(viewport={'width':390,'height':844},has_touch=True,
                                        reduced_motion='reduce',color_scheme=theme)
            page=context.new_page()
            for name,stem in PAGES.items():
                page.goto((HTML/f'{stem}.html').as_uri(),wait_until='load')
                page.evaluate('t=>document.documentElement.dataset.theme=t',theme)
                diagram=page.locator(f'.teaching-diagram[data-diagram="{name}"]')
                radios=diagram.locator('input[type=radio]')
                for i in range(radios.count()):
                    radio=radios.nth(i)
                    label=diagram.locator('label[for="'+radio.get_attribute('id')+'"]')
                    label.tap()
                    assert radio.is_checked(),(name,'tap',i)
                    panel=diagram.locator('#'+radio.get_attribute('aria-controls'))
                    assert panel.is_visible(),(name,'tap panel',i)
                    label.hover()
                    assert panel.is_visible(),(name,'hover',i)
                page.mouse.move(0,0)
                radios.first.focus()
                page.keyboard.press('ArrowRight')
                assert radios.nth(1).is_checked(),(name,'keyboard')
                panel=diagram.locator('#'+radios.nth(1).get_attribute('aria-controls'))
                assert panel.is_visible(),(name,'keyboard panel')
                diagram.screenshot(path=str(OUT/f'{theme}-{name}.png'))
                row=dict(diagram=name,page=stem,theme=theme,touch=True,keyboard=True,
                         hover_media_supported=page.evaluate('matchMedia("(hover:hover) and (pointer:fine)").matches'),
                         reduced_motion=page.evaluate('matchMedia("(prefers-reduced-motion: reduce)").matches'),
                         options=radios.count(),html_sha256=sha256((HTML/f'{stem}.html').read_bytes()).hexdigest())
                rows.append(row)
                print(json.dumps(row),flush=True)
                (OUT/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
            context.close()
        context=browser.new_context(viewport={'width':1440,'height':900})
        page=context.new_page()
        hover_rows=[]
        for name,stem in PAGES.items():
            page.goto((HTML/f'{stem}.html').as_uri(),wait_until='load')
            diagram=page.locator(f'.teaching-diagram[data-diagram="{name}"]')
            radios=diagram.locator('input[type=radio]')
            for i in range(radios.count()):
                radio=radios.nth(i)
                diagram.locator('label[for="'+radio.get_attribute('id')+'"]').hover()
                assert diagram.locator('#'+radio.get_attribute('aria-controls')).is_visible(),(name,'desktop hover',i)
            hover_rows.append(dict(diagram=name,passed=True,options=radios.count()))
        (OUT/'hover_results.json').write_text(json.dumps(hover_rows,indent=2)+'\n')
        page.goto((HTML/'05_seismic.html').as_uri(),wait_until='load')
        themes=[]
        for _ in range(3):
            page.locator('button.theme-switch-button').click()
            page.wait_for_function('''()=>{const c=document.documentElement.dataset.theme==='dark'?'#1e232b':'#fbfcfd'; return [...document.querySelectorAll('.js-plotly-plot')].every(e=>e._fullLayout?.paper_bgcolor===c)}''',timeout=15000)
            themes.append(page.evaluate('document.documentElement.dataset.theme'))
        assert set(themes)=={'light','dark'},themes
        (OUT/'theme_button.json').write_text(json.dumps({'passed':True,'themes':themes},indent=2)+'\n')
        browser.close()
    print('PASSED',len(rows),'diagram/theme combinations and actual theme button',flush=True)


if __name__=='__main__':main()
