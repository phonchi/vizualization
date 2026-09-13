"""Server-render the science-exhibition shell without changing document routes."""
from pathlib import Path
from html import escape
import json
import re
import yaml
from bs4 import BeautifulSoup

KINDS = {
 '01_forecast_question':('question','從一個問題出發'), '02_reading_horus':('catalog','讓資料說話'),
 '03_experiment_spec':('protocol','先把實驗說清楚'), '04_poisson_and_sup':('statistics','把直覺寫成機率'),
 '05_likelihood_estimation':('statistics','從資料找參數'), '06_simulation_tests_scores':('evaluation','讓模型接受檢驗'),
 '07_gr_bvalue':('catalog','讀出大小地震的比例'), '08_completeness':('catalog','看見目錄的邊界'),
 '09_ppe_forecast':('model','從歷史位置到預報'), '10_clustering_laws':('pattern','找出叢集的規律'),
 '11_conditional_intensity':('model','每顆事件，都能改變下一刻'), '12_etas_structure':('model','一顆事件，展開一個家族'),
 '13_etas_italy_forecast':('case','讓模型走進義大利'), '14_psi_precursory_scale':('pattern','從現象走向可檢驗的假說'),
 '15_eepas_italy_forecast':('model','把一顆事件放進三種尺度'), '16_test_number':('evaluation','數量對不對？'),
 '17_test_space_magnitude':('evaluation','位置與規模對不對？'), '18_test_comparison':('evaluation','哪一張預報增加了資訊？'),
 '19_ensembles':('synthesis','讓不同觀點相遇'), '20_beyond_forecast':('outlook','把模型帶回論文與現實'),
 '01_overview':('observation','回到地表，回到觀測'), '02_download':('practice','資料從哪裡來'),
 '03_groundwater':('observation','地下水留下的線索'), '04_geomagnetic':('observation','讀懂看不見的磁場'),
 '05_seismic':('observation','把搖晃變成可以分析的資料'), '06_gnss':('observation','量出地表的微小位移'),
 '07_case_hualien2024':('case','讓不同儀器對準同一場地震'), '08_explore_ideas':('outlook','從觀察走向證據'),
 '23_taiwan_outlook':('case','把義大利的問題，帶回臺灣'),
}


def decorate(app,pagename,templatename,context,doctree):
    src=Path(app.srcdir)
    toc=yaml.safe_load((src/'_toc.yml').read_text())
    names=[toc['root']]+[c['file'] for part in toc['parts'] for c in part['chapters']]
    if pagename not in names or not context.get('body'):
        return
    link=context['pathto']
    titles={n:app.env.titles[n].astext() if n in app.env.titles else n for n in names}
    soup=BeautifulSoup(context['body'],'html.parser')
    heading=soup.find('h1')
    if heading is None:return
    home=pagename==toc['root']
    kind,kicker=KINDS.get(pagename,('appendix','深入推導，回查方法'))
    chapter_match=re.match(r'(\d+)\.\s*',heading.get_text())
    chapter=chapter_match.group(1) if chapter_match else (pagename.split('_')[1].upper() if pagename.startswith('appendix_') else '')
    section_items=[]
    for section in soup.find_all('section'):
        h=section.find(['h2','h3'],recursive=False)
        if h is None:continue
        classes=list(section.get('class',[]))+['museum-section']
        if section.select('.teaching-diagram'):classes+=['section-stage']
        elif section.select('.js-plotly-plot'):classes+=['section-data']
        elif section.find('table'):classes+=['section-comparison']
        elif section.select('.math'):classes+=['section-theory']
        section['class']=classes
        if h.name=='h2':section_items.append((section.get('id',''),h.get_text().replace('¶','').strip().rstrip('#').strip()))
    nav=''.join(f'<a href="{escape(link(n))}">{label}</a>' for n,label in [('01_forecast_question','義大利實驗'),('01_overview','臺灣觀測'),('appendix_a_point_process','技術附錄')])
    groups=[]
    for part in toc['parts']:
        links=''.join(f'<a {"aria-current=page" if c["file"]==pagename else ""} href="{escape(link(c["file"]))}">{escape(titles[c["file"]])}</a>' for c in part['chapters'])
        groups.append(f'<section><h3>{escape(part["caption"])}</h3>{links}</section>')
    chrome=f'''<header class="museum-topbar"><a class="museum-brand" href="{escape(link(toc['root']))}"><b>地震統計與預報</b><span>從義大利實驗到臺灣觀測</span></a><nav class="museum-primary" aria-label="主要入口">{nav}</nav><div class="museum-actions"><a class="museum-search" href="{escape(link('search'))}" aria-label="搜尋全書"><svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="10" cy="10" r="6.5"/><path d="m15 15 6 6"/></svg></a><span data-museum-theme-slot></span><button type="button" data-museum-menu aria-label="開啟章節目錄">目錄 <span aria-hidden="true">☰</span></button></div><div class="museum-reading-progress" aria-hidden="true"></div></header><dialog class="museum-menu" aria-label="全書章節目錄"><div class="museum-menu-heading"><h2>選一條路，開始探索。</h2><button type="button" data-museum-close aria-label="關閉目錄">關閉 ×</button></div><div class="museum-menu-columns">{''.join(groups)}</div></dialog>'''
    if home:
        data_path=src/'_static/diagrams/data/italy_exhibit_data.json'
        data=json.loads(data_path.read_text()) if data_path.exists() else {}
        payload=json.dumps({k:data.get(k,[]) for k in ['cells','polygon','targets']},ensure_ascii=False).replace('</','<\\/')
        hero=f'''<section class="museum-home-hero"><div class="museum-home-copy"><p class="museum-eyebrow">一份真實資料，兩條探索路線</p><h1>地震會來。<br><em>預報能說多少？</em></h1><p class="museum-home-lead">從義大利的一份預報實驗出發，學會讀資料、看模型、辨認證據。再回到臺灣，聽見儀器記下的地球。</p><div class="museum-home-links"><a class="museum-enter" href="{escape(link('01_forecast_question'))}">走進義大利實驗 <span aria-hidden="true">↗</span></a><a href="{escape(link('01_overview'))}">從臺灣觀測開始 →</a></div><div class="museum-facts"><span><b>177</b>測試方格</span><span><b>40</b>預報窗口</span><span><b>{len(data.get('targets',[]))}</b>目標地震</span></div></div><div class="museum-atlas"><div class="museum-atlas-label"><span>義大利實驗格網</span><b data-atlas-year>2021</b></div><svg data-museum-atlas viewBox="0 0 520 600" role="img" aria-label="義大利177格測試區與2012至2021年目標地震的事後回放"></svg><div class="museum-atlas-controls"><button type="button" data-atlas-play aria-label="播放目標地震回放">播放回放</button><input type="range" min="2012" max="2021" step="1" value="2021" aria-label="目標地震回放年份" data-atlas-range><output data-atlas-count></output></div><p class="museum-atlas-caption">紅點是事後觀測；格網沿用既有研究區。這是資料回放，不是即時預警。</p></div><script type="application/json" id="museum-atlas-data">{payload}</script></section>'''
        anchor=heading.parent.get('id','')
        heading.replace_with(BeautifulSoup(hero,'html.parser'))
    else:
        clean_title=re.sub(r'^\d+\.\s*','',heading.get_text().replace('¶','').strip().rstrip('#').strip())
        headid=heading.get('id')
        first=heading.find_next_sibling('p')
        lead=str(first.extract()) if first is not None else ''
        waypoints=''.join(f'<a href="#{escape(sid)}"><span>{i+1:02d}</span>{escape(title)}</a>' for i,(sid,title) in enumerate(section_items[:3]))
        hero=f'''<header class="museum-chapter-header" data-chapter-kind="{kind}"><div class="museum-chapter-label"><span>{'技術附錄' if kind=='appendix' else ('第一部 · 義大利實驗' if pagename in list(KINDS)[:20] else '第二部 · 臺灣觀測')}</span><b>{escape(chapter)}</b></div><div class="museum-chapter-heading"><p class="museum-kicker">{escape(kicker)}</p><h1{f' id="{escape(headid)}"' if headid else ''}>{escape(clean_title)}</h1><div class="museum-chapter-lead">{lead}</div></div><nav class="museum-waypoints" aria-label="本章起點">{waypoints}</nav></header>'''
        heading.replace_with(BeautifulSoup(hero,'html.parser'))
    context['body']=f'<div class="museum-shell museum-{kind}{" museum-home" if home else ""}" data-page="{escape(pagename)}">{chrome}{str(soup)}</div>'


def setup(app):
    app.connect('html-page-context',decorate,priority=800)
    return {'version':'1.0','parallel_read_safe':True,'parallel_write_safe':True}
