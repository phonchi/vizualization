# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: kernelspec,jupytext
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 3. 標準預報實驗規格：把題目一次說清楚
#
# 這個實驗的題目可以寫成一句話：從 2012 年 1 月 1 日起的十年，
# 義大利 177 格內，每三個月、每一格、每 0.1 個規模箱，預期有幾顆 Mw ≥ 5.0 的地震？
# 這一句裡藏了五個決定：哪幾年、多久發一次、多大的格、多寬的箱、多大的地震。
#
# 上一章填好了資料側的欄位。這一章把題目本身的欄位填完，
# 順便把論文方法段裡的用語，和另一套常見的統計軟體用語對起來。
# 讀完你能把一篇網格預報論文的方法段，逐項對回這張規格卡。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import HTML, display

from gdms_toolkit import italy, italy_models
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

spec = italy.SPEC
cat = italy.experiment_catalog()
cells = italy.testing_cells()
poly = italy.collection_polygon()
windows = italy.forecast_windows()
learn = italy.target_events(cat, "learning")
test = italy.target_events(cat, "testing")

# %% [markdown]
# ## 3.1 三段期間：暖機、學習、測試
#
# 工作目錄從 1960 年到 2021 年，共 62 年。實驗把它切成三段，各有各的用途。
#
# - **暖機期**（warm-up period，1960–1989）：只當歷史。這 30 年的地震會進入
#   模型的輸入，但不用來估參數，也不評分。上一章看過，這段期間只有 Mw ≥ 4.0
#   是完整的。
# - **學習期**（learning period，1990–2011）：估參數。模型在這 22 年的資料上
#   調整自己的參數，直到最能解釋這段期間的目標地震。
# - **測試期**（testing period，2012–2021）：只評分。參數在學習期結束時凍結，
#   這 10 年的目標地震只用來檢驗預報。
#
# **目標地震**（target earthquake）是被評分的那些地震：在測試區 R 內、
# 取整規模 $m_b \ge 5.0$、落在指定期間內。兩段期間各有幾顆，直接數出來。

# %% tags=["remove-input"]
print(f"學習期 {spec.learning[0]}–{spec.learning[1]}：目標地震 {len(learn)} 顆"
      f"（{learn.time.min():%Y-%m-%d} 至 {learn.time.max():%Y-%m-%d}）")
print(f"測試期 {spec.testing[0]}–{spec.testing[1]}：目標地震 {len(test)} 顆"
      f"（{test.time.min():%Y-%m-%d} 至 {test.time.max():%Y-%m-%d}）")

# %% [markdown]
# 學習期本站數到 27 顆，與論文報告一致。測試期本站數到 25 顆，論文為 27 顆。
# 本站使用 2024 年 6 月版 HORUS，與論文使用的版本不同。
# 版本修訂可能改變規模、位置與收錄情形，但尚未取得逐事件差異對照。
# 因此目前不能斷定是哪兩顆、也不能把差額全部歸因於某一欄位修訂。
# 比較數字前，須一起核對目錄版本、篩選條件與分箱規則。
#
# 為什麼要把學習和測試分開？因為模型有很多可以調的地方。
# 如果調參數時看得到測試期的地震，模型會不知不覺往那些地震靠。
# 這叫**資料窺探**（data snooping）：任何用測試期資料影響模型的決定。
# 注意界線在哪裡。用測試期資料**檢驗**預報不是窺探，那正是測試期的用途；
# 看了檢驗結果**再回頭改模型**，才是窺探。
#
# 上一章提過一個例外：177 格的選區用到了 2021 年前的資訊。
# 本站沿用既有研究區，這件事要在讀結果時記住。

# %% [markdown]
# ## 3.2 發報時間軸：每三個月發一張預報
#
# 測試期不是只發一張十年預報。論文用三個月的預報間隔；本站把十年切成 40 個
# 91.31 天的窗，近似一年的四分之一。發報步長同為 91.31 天。
# 相鄰窗口不重疊，每顆目標地震只在一窗計數。
# 最後一窗到 2021 年 12 月 31 日 09:36，年底剩餘 0.6 天不評分。
# 學習期另切 88 窗，最後一窗裁到 2012 年 1 月 1 日，避免跨入測試期。

# %% tags=["remove-input"]
show_diagram("d03_timeline_windows",
             caption="三段期間、40 個滾動窗、每窗的發報時刻與資料截止。PPE／EEPAS 另加 50 天 delay。")

# %% [markdown]
# 每一窗開始的那一刻叫**發報時刻**（issue time）：預報在這一刻算出。
# 發報時刻是資料可用性的最晚界線。模型還可以把截止時間往前移。
# **資料截止**（data cutoff）指模型最後可納入歷史的時刻。
# 窗內發生的地震，模型在發報時看不到。
#
# 這裡有一組容易混淆的概念。**參數凍結**（frozen parameters）指參數在學習期
# 結束後就不再改。**歷史不凍結**指下一窗仍可接入剛發生的事件。
# SUP 保留學習期估出的固定平均率。PPE、ETAS 與 EEPAS 則逐窗更新歷史。
# 這種作法叫**滾動更新**（rolling update）：每個窗口有自己的發報時刻與輸入截止。
# 已在測試期發生的事件可以成為下一窗歷史，但不能反過來改上一窗的預報。
#
# PPE 與 EEPAS 都另設 50 天延遲；ETAS 使用發報前的歷史，沒有這項延遲。
# 這 50 天叫 **delay**，Biondini 等人在附錄寫「通常取 50 天」。
# 目的是避開餘震群的短期叢集，讓模型只看長期的變化；第 15 章會再談。
#
# 下圖把 40 個窗攤在時間軸上，紅點是 25 顆目標地震。下半放大 2016 年第三窗，
# 標出發報時刻、資料截止與 50 天 delay 的位置。

# %% tags=["remove-input"]
w_idx = np.searchsorted(windows.t2.to_numpy(), test.t_days.to_numpy(), side="right") + 1
test = test.assign(window=w_idx)
event_windows = sorted(set(w_idx))
k = int(windows.index[windows.start <= pd.Timestamp("2016-08-24")].max())   # 含 Amatrice 主震的窗
w = windows.loc[k]
fig = make_subplots(rows=2, cols=1, vertical_spacing=0.16, row_heights=[0.55, 0.45],
                    subplot_titles=(f"測試期 {len(windows)} 個窗，每窗 {spec.window_days} 天",
                                    f"放大：第 {w.window} 窗（{w.start:%Y-%m-%d} 發報）"))
for _, r in windows.iterrows():
    if r.window % 2 == 0:
        fig.add_vrect(x0=r.start, x1=r.end, fillcolor="rgba(42,120,214,0.10)", line_width=0, row=1, col=1)
fig.add_trace(go.Scatter(x=test.time, y=test.mb, mode="markers", name="目標地震（25 顆）",
                         marker=dict(color=QUAKE_COLOR, size=7, line=dict(color="white", width=0.8)),
                         text=[f"{t:%Y-%m-%d} Mw {m:.1f}，第 {wi} 窗" for t, m, wi in zip(test.time, test.mb, test.window)],
                         hovertemplate="%{text}<extra></extra>"), row=1, col=1)
for wi in event_windows:
    r = windows.loc[wi - 1]
    fig.add_annotation(x=r.start + (r.end - r.start) / 2, y=6.95, text=str(wi), showarrow=False,
                       font=dict(size=10, color="rgba(60,70,90,0.9)"), row=1, col=1)
fig.update_xaxes(range=[windows.start.iloc[0], windows.end.iloc[-1]], row=1, col=1)
fig.update_yaxes(title_text="取整規模 Mw", range=[4.8, 7.1], row=1, col=1)
# 放大一窗
t_issue = w.start
t_delay = t_issue - pd.Timedelta(days=spec.delay_days)
zoom0, zoom1 = t_delay - pd.Timedelta(days=40), w.end + pd.Timedelta(days=15)
hist = cat[(cat.time >= zoom0) & (cat.time < zoom1) & (cat.mb >= 3.5)]
fig.add_vrect(x0=t_issue, x1=w.end, fillcolor="rgba(42,120,214,0.10)", line_width=0, row=2, col=1)
fig.add_vrect(x0=t_delay, x1=t_issue, fillcolor="rgba(235,104,52,0.15)", line_width=0, row=2, col=1)
fig.add_vline(x=t_issue, line_color="rgba(36,48,64,0.8)", line_width=2, row=2, col=1)
fig.add_trace(go.Scatter(x=hist.time, y=hist.mb, mode="markers", name="S 內 Mw ≥ 3.5（含 R 外）",
                         marker=dict(color=ACCENT, size=5, opacity=0.6)), row=2, col=1)
in_win = test[(test.time >= t_issue) & (test.time < w.end)]
fig.add_trace(go.Scatter(x=in_win.time, y=in_win.mb, mode="markers", showlegend=False,
                         marker=dict(color=QUAKE_COLOR, size=9, line=dict(color="white", width=0.8)),
                         hoverinfo="skip"), row=2, col=1)
for xx, txt, ax in [(t_issue, "發報時刻（ETAS 截止）", -10), (t_delay, "PPE／EEPAS 截止（50 天前）", -10)]:
    fig.add_annotation(x=xx, y=6.6, text=txt, showarrow=False, xanchor="right", xshift=ax,
                       font=dict(size=11), row=2, col=1)
fig.add_annotation(x=t_issue + (w.end - t_issue) / 2, y=6.6, text="這一窗被評分", showarrow=False,
                   font=dict(size=11), row=2, col=1)
fig.update_xaxes(range=[zoom0, zoom1], row=2, col=1)
fig.update_yaxes(title_text="取整規模 Mw", range=[3.4, 6.9], row=2, col=1)
gap_days = italy.year_start_days(spec.testing[1] + 1) - windows.t2.iloc[-1]
n_after = int((test.t_days >= windows.t2.iloc[-1]).sum())
apply_layout(fig, height=680, hovermode="closest", legend=dict(orientation="h", y=-0.06),
             title=f"{len(event_windows)} 個窗有目標地震（窗號標在上緣）；"
                   f"最後一窗到年底差 {gap_days:.1f} 天，其中 {n_after} 顆目標")
fig

# %% [markdown]
# 上半圖：25 顆目標地震集中在 12 個窗裡，其餘 28 個窗一顆都沒有。
# 這是預報要面對的常態：大多數時候什麼都沒發生，少數時候一次來一串。
#
# 下半圖：第 19 窗在 2016 年 7 月 1 日發報。灰線左邊的藍點是發報前 S 內的地震，
# ETAS 可以使用這些歷史；PPE 與 EEPAS 則排除橘色帶中的最近 50 天。
# 紅點是 8 月 24 日的 Amatrice 主震與它 1 小時後的 Mw 5.5，兩顆都在窗內，
# 模型發報時看不到它們。
#
# ### 檢驗方式：哪些資料已經被誰看過？
#
# **回溯檢驗**（retrospective testing）以已存在的資料測試預報。
# 測試資料可能用於擬合，也可能完全留出。
# 所以「回溯」不等於「用訓練資料評分」。
# 關鍵還包括建模者是否知道測試期的活動特徵。
#
# **擬前瞻檢驗**（pseudo-prospective testing）是回溯檢驗的子類。
# 它把模型估計資料放在前，把測試資料放在後。
# 重演中每次只使用指定截止前的事件。
# 這降低模型直接使用未來事件的風險。
# 但選區、門檻或模型選擇仍可能受已知結果影響。
#
# **前瞻檢驗**（prospective testing）使用模型及建模者尚未見過的測試結果。
# **即時前瞻**（real-time prospective）在當時產生並保存預報。
# **延遲前瞻**（delayed prospective）可在稍後由受控程序執行。
# 後者仍須事先固定流程及測試規格，嚴格隔離未來資料。
# 單純事後把資料表切兩半，不能取得這個資格。
#
# **部分前瞻**（partially prospective）同時評估已知與未知的時段。
# 例如主震後幾天才出預報，卻從主震時刻開始計分。
# 最初幾天已發生的餘震，須和真正未來的部分分開說明。
# Mizrahi 等（2024）的圖 2 用「資料／模型／建模者」呈現這些差異。
#
# | 檢驗安排 | 模型取得的資料 | 建模者與結果的關係 |
# |---|---|---|
# | 回溯、樣本內 | 擬合可使用被評估的資料 | 結果已存在 |
# | 回溯、樣本外 | 評估資料未用於擬合，未必依時間分離 | 結果仍可能已知 |
# | 回溯、擬前瞻 | 依過去截止製作後續預報 | 模仿未知，不能保證未受影響 |
# | 前瞻、即時或受控延遲 | 依預先鎖定的資料邊界 | 未來結果隔離，保留可稽核紀錄 |
# | 部分前瞻 | 預報範圍同時含過去與未來 | 部分結果可能已知 |
#
# ### 模型與發報時程是另外幾個問題
#
# **時間相依**（time-dependent）描述率如何隨時間或歷史改變。
# **時間獨立**（time-independent）表示指定期間內沒有這種時間變化。
# 這與模型有沒有在空間上均勻，是兩個問題。
# SUP 空間均勻；SVP 空間不均勻，兩者仍可時間獨立。
# SVP 的定義與 PPE 的區別見第 4、9 章。
#
# **滾動更新**（rolling update）描述重新發報及更新可用歷史的安排。
# **即時**（real-time）描述資料與發報的執行時機。
# **作業預報**（operational forecasting）描述持續提供權威資訊供決策的用途。
# 作業系統中的模型可以事後接受回溯評估。
# 時間相依模型也可以只做擬前瞻重演。
# 這些名稱不能排成同一列互斥選項。
#
# 本站使用 2024 年修訂 HORUS，依事件時間回放 2012–2021。
# 它沒有保存各次歷史發報當時可得的速報版本。
# 177 格又沿用含 2021 年前資訊的選區。
# 因此，本站稱「修訂目錄上的擬前瞻教學重演」。
# 它屬回溯檢驗，不是真前瞻、即時作業或完全事前凍結的實驗。

# %% tags=["remove-input"]
show_diagram("d03_forecast_protocols", caption="先看誰能接觸哪些結果，再區分檢驗方式；更新頻率與作業用途另列。")

# %% [markdown]
# ## 3.3 三個門檻：$M_c$、$m_0$、$m_T$
#
# 規模門檻有三個，各管一件事。上一章已經遇過第一個。

# %% tags=["remove-input"]
show_diagram("d03_thresholds",
             caption="規模取整與輸入／目標門檻的有效邊界；完整度 M_c 的角色見下方定義。")

# %% [markdown]
# ```{admonition} 定義：三個規模門檻
# :class: definition
# - **完整度規模** $M_c$（magnitude of completeness）：目錄的性質。從這個規模起，
#   目錄幾乎記到了全部的地震。由資料估出，第 8 章給數值。
# - **輸入門檻** $m_0$（input threshold）：模型的決定。取整規模 $\ge m_0$ 的地震
#   才進入模型當歷史。本實驗 $m_0 = 2.45$，名目寫 2.5。
# - **目標門檻** $m_T$（target threshold）：題目的決定。取整規模 $\ge m_T$ 的地震
#   才被評分。本實驗 $m_T = 5.0$，有效值 4.95。
# ```
#
# 三者的關係是 $M_c \le m_0 < m_T$。$m_0$ 不能低於 $M_c$：
# 若把不完整的小地震餵給模型，模型會把「沒記到」學成「沒發生」。
# 論文引 Rhoades 與 Evison（2004）：輸入目錄的完整度最好比目標門檻低兩級以上，
# 這是 $m_0$ 取 2.5、$m_T$ 取 5.0 的來由。
# 學習期從 1990 年開始，也正是 HORUS 對 Mw ≥ 2.5 完整的起點。
#
# 為什麼目標是 5.0？論文給的理由是義大利 Mw ≥ 5 的地震可能造成建築損害，
# 威脅居民安全。這是一個社會決定，不是統計決定。

# %% [markdown]
# ## 3.4 測試期的 25 顆目標地震
#
# 下圖沿用紅色標出目標地震，之後各章的預報圖也使用相同約定。

# %% tags=["remove-input"]
ring_lon, ring_lat = [], []
for _, c in cells.iterrows():
    ring_lon += list(c.corner_lon) + [None]
    ring_lat += list(c.corner_lat) + [None]
emilia = test[(test.time >= "2012-05-20") & (test.time < "2012-07-01")]
central = test[(test.time >= "2016-08-24") & (test.time < "2017-02-01")]
fig = go.Figure()
fig.add_trace(go.Scattergeo(lon=poly.lon, lat=poly.lat, mode="lines", name="收集區 S",
                            line=dict(color="rgba(36,48,64,0.7)", width=1.5), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(lon=ring_lon, lat=ring_lat, mode="lines", name="測試區 R（177 格）",
                            line=dict(color="rgba(36,48,64,0.45)", width=0.6), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(
    lon=test.lon, lat=test.lat, mode="markers", name="目標地震（測試期）",
    marker=dict(color=QUAKE_COLOR, size=4 + 6 * (test.mb - 5.0), line=dict(color="white", width=0.8),
                opacity=0.85),
    text=[f"{t:%Y-%m-%d} Mw {m:.1f}，格 {c}" for t, m, c in zip(test.time, test.mb, test.cell)],
    hovertemplate="%{text}<extra></extra>"))
for df, label, dy in [(emilia, f"Emilia 2012：{len(emilia)} 顆", 1.0),
                      (central, f"中義大利 2016–17：{len(central)} 顆", -1.2)]:
    fig.add_annotation(x=df.lon.mean(), y=df.lat.mean() + dy, text=label, showarrow=False,
                       font=dict(size=12, color=QUAKE_COLOR), xref="x", yref="y")
fig.add_trace(go.Scattergeo(lon=[emilia.lon.mean(), central.lon.mean()],
                            lat=[emilia.lat.mean() + 0.9, central.lat.mean() - 1.0],
                            mode="text", text=[f"Emilia 2012：{len(emilia)} 顆",
                                               f"中義大利 2016–17：{len(central)} 顆"],
                            textfont=dict(size=12, color=QUAKE_COLOR), showlegend=False,
                            hoverinfo="skip"))
fig.layout.annotations = ()
fig.update_geos(fitbounds="locations", visible=False, showcountries=True,
                countrycolor="rgba(120,130,150,0.5)", showcoastlines=True,
                coastlinecolor="rgba(120,130,150,0.5)", showland=True,
                landcolor="rgba(238,242,247,0.6)", projection_type="mercator")
apply_layout(fig, height=560, hovermode="closest", margin=dict(l=10, r=10, t=50, b=10),
             legend=dict(orientation="h", y=-0.02),
             title=f"測試期 {len(test)} 顆目標地震，落在 {test.cell.nunique()} 格；點的大小隨規模")
fig

# %% [markdown]
# 25 顆裡有兩串。2012 年 5 月 20 日到 29 日的 Emilia 序列貢獻 7 顆，
# 2016 年 8 月到 2017 年 1 月的中義大利序列（Amatrice、Norcia、Campotosto）
# 貢獻 10 顆。兩串加起來 17 顆，占了三分之二。其餘 8 顆散在各地。
# 這個結構會一路影響到檢驗：一張預報能不能拿到好分數，很大一部分取決於
# 它在序列開始後有沒有跟上。下表列出全部 25 顆。

# %% tags=["remove-input"]
tbl = pd.DataFrame({
    "發震時刻（UTC）": test.time.dt.strftime("%Y-%m-%d %H:%M"),
    "緯度、經度": [f"{la:.3f}、{lo:.3f}" for la, lo in zip(test.lat, test.lon)],
    "深度（km）": test.depth.round(1),
    "Mw（取整）": [f"{m:.2f}（{b:.1f}）" for m, b in zip(test.mw, test.mb)],
    "格、窗": [f"{c}、{wi}" for c, wi in zip(test.cell, test.window)],
})
tbl.index = np.arange(1, len(tbl) + 1)
tbl

# %% [markdown]
# 表裡有三顆原始 Mw 低於 5.0（4.98、4.98、4.95），取整後才成為目標。
# 這就是上一章說的：門檻 5.0 的有效值是 4.95。
#
# ## 3.5 格與箱：預報的解析度
#
# 一張預報不是一個數字，而是一個陣列。三個維度分別是時間、空間、規模：
#
# - 時間：40 個窗，每窗 91.31 天。
# - 空間：177 格，每格 1,800 km²。
# - 規模：25 箱，從 5.0 到 7.5，每箱寬 0.1；取整規模落在 $[m_1, m_2)$ 算該箱。
#
# 所以一個模型在測試期的完整輸出，是形狀 (40, 177, 25) 的陣列，
# 共 177,000 個數字。每個數字是一個**期望數**（expected number）：
# 那一窗、那一格、那一箱，預期發生幾顆地震。本站記作 $\Lambda_{jk}$，
# $j$ 是格、$k$ 是規模箱，窗的下標省略。期望數是什麼、怎麼從率算出來，
# 是下一章的主題。
#
# 格與箱的大小就是預報的**解析度**（resolution）。格越小，預報越精細，
# 但每格的期望數越小、越難檢驗。177 格的設計是承接 Gasperini 等人（2021）
# 的既有工作，讓結果能和先前研究比較。

# %% tags=["remove-input"]
sup = italy_models.get_forecast("SUP", "testing")     # 讀快取，不重算
print(f"SUP 測試期預報陣列形狀：{sup.shape}，共 {sup.size:,} 個期望數")
print(f"規模箱邊界：{italy.magnitude_edges()[:4]} … {italy.magnitude_edges()[-2:]}，共 {len(italy.magnitude_edges()) - 1} 箱")

# %% [markdown]
# ## 3.6 輸出格式：CSEP 十欄格式
#
# 陣列要交給別人檢驗，得有一個約定的檔案格式。CSEP（Collaboratory for the
# Study of Earthquake Predictability，地震可預測性研究合作組織）用的是
# **十欄格式**：每一列是一個空間–規模箱，前八欄是箱的邊界，第九欄是期望數，
# 第十欄是遮罩旗標。下表用 SUP 模型的十年合計示範前 10 列。
# SUP 是下一章要做的第一張預報，這裡只借它的格式。

# %% tags=["remove-input"]
ten_year = sup.sum(axis=0)                                  # (177, 25)：40 窗加總
out_path = Path(tempfile.mkdtemp()) / "SUP_10yr_177cells.dat"
italy.write_csep_10col(ten_year, out_path)
ten_col = pd.read_csv(out_path, sep="\t")
print(f"檔案 {out_path.name}：{len(ten_col):,} 列（177 格 × 25 箱），十年合計期望數 {ten_col.RATE.sum():.2f}")
ten_col.head(10).style.format({c: "{:.3f}" for c in ["LON_0", "LON_1", "LAT_0", "LAT_1"]} |
                              {"RATE": "{:.3e}"})

# %% [markdown]
# 欄位意義：LON_0／LON_1 與 LAT_0／LAT_1 是格的經緯度邊界，Z_0／Z_1 是深度範圍
# （0–40 km），MAG_0／MAG_1 是規模箱的左右邊界，RATE 是期望數，
# FLAG 為 1 表示這一箱參與評分。前 10 列是同一格的前 10 個規模箱，
# 期望數隨規模遞減，這是 GR 律的形狀，第 7 章會講。
#
# 這張表有一個簡化。177 格在公里座標上是正方形，在經緯度上是梯形；
# 上表寫的是每格的經緯度外接框，鄰格會略有重疊。正式交給 CSEP 檢驗時，
# 要先把 177 格重新分配到 0.1° 的經緯度格，第 9 章與附錄 E 會做。
#
# Savran 等人（2022）的 pyCSEP 套件把這件事拆成兩個物件。
# **Region** 是空間格加上規模箱的左邊界清單；格的下邊與左邊界包含在內。
# **GriddedForecast** 是一個二維陣列（格 × 箱）加上一個 Region，
# 而且一個物件只管一段時間。所以本站的 (40, 177, 25) 陣列，
# 在 pyCSEP 裡是 40 個 GriddedForecast。

# %% [markdown]
# ## 3.7 兩套術語：CSEP 用語與 ETAS R 套件用語
#
# 讀點過程軟體的文件時，會遇到另一套說法。Jalilian（2019）的 ETAS R 套件
# 是最常見的例子。下表把兩套用語對起來；本站正文用左欄，第一次出現時括號附右欄。
#
# | CSEP／EEPAS（本站用語） | ETAS R 套件（Jalilian 2019） |
# |---|---|
# | 測試區 R（testing region） | study region |
# | 收集區 S（collection region） | complementary events 所在區 |
# | 學習期／測試期 | study period |
# | 目標地震（target） | target events |
# | 提供歷史的事件（precursor／trigger） | complementary events |
# | decimal days 自 1960-01-01 | decimal days 自 `time.begin` |
#
# 有一個字母陷阱。Jalilian 把 study region 記作 $S$，那是我們的 **R**；
# 我們的 S 在他的架構裡沒有名字，只是 complementary events 所在的地方。
# 看到 ETAS R 的公式裡有 $S$，先確認它指的是哪一個。
#
# 兩套術語背後的集合位置是一樣的。**target events** 是同時滿足三件事的地震：
# 在 study region 內、在 study period 內、規模不低於門檻。
# **complementary events** 是其餘仍進入歷史的地震：在區域外，或在 study period
# 開始之前。它們不被評分，但會影響 target events 的率。這正是上一章的
# 收集區 S 與暖機期的角色。
#
# Jalilian（2019）第 16 頁的範例可以直接對回規格卡：`study.start` 與 `study.end`
# 是我們的期間邊界；`region.poly` 是我們的 R，頂點要反時針、首點不重複；
# `mag.threshold` 是我們的 $m_0$；時間軸是自 `time.begin` 起算的 decimal days，
# 我們的零點固定在 1960 年 1 月 1 日。

# %% [markdown]
# ## 3.8 本章填入的規格欄位
#
# 規格卡到這裡只剩兩個空格：$M_c$ 的數值（第 8 章）與評分方式（第 16–18 章）。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("目錄與規模尺度", "HORUS 2024 年 6 月版（Lolli et al. 2020），Mw 均一化，規模取整到 0.1"),
    ("時間範圍", "1960–2021 年"),
    ("深度", "≤ 40 km；排除非地震事件"),
    ("收集區 S", "CPTI15 多邊形（14 個頂點）"),
    ("測試區 R", "177 個 30√2 km 方格（EPSG:7794 公里座標）；選格屬既有研究區的重演"),
    ("三段期間", f"暖機 {spec.warmup[0]}–{spec.warmup[1]}、學習 {spec.learning[0]}–{spec.learning[1]}"
                f"（目標 {len(learn)} 顆）、測試 {spec.testing[0]}–{spec.testing[1]}（目標 {len(test)} 顆；論文 27）"),
    ("發報時間軸", f"測試期切 {len(windows)} 個 {spec.window_days} 天的窗；每窗只用發報時刻前的事件；"
                  f"參數凍結於 2011 年底；PPE／EEPAS 另加 {spec.delay_days:.0f} 天 delay"),
    ("門檻", f"m_0 = {spec.m0}（名目 2.5）、m_T = {spec.mT}（有效 4.95）；M_c 待填（第 8 章）"),
    ("規模箱", f"{spec.mT}–{spec.m_max}，寬 {spec.mag_bin}，共 {len(italy.magnitude_edges()) - 1} 箱"),
    ("輸出格式", f"每窗 × {spec.n_cells} 格 × 25 箱的期望數陣列；CSEP 十欄格式（外接框示範）；0.1° 重新分格待第 9 章"),
    ("評分方式", "待填（第 16–18 章）"),
])))

# %% [markdown]
# 題目說清楚了。下一步是做出第一張預報：假設地震在 177 格內均勻、
# 在時間上也均勻，那每窗每格該預期幾顆？這需要先弄懂「期望數」是什麼。
# 下一章：{doc}`統計工具箱 I：計數、Poisson 過程，與第一張預報 SUP <04_poisson_and_sup>`。
#
# 分類依據：Mizrahi 等（2024）§2.2.1.4–7、圖 2。
# 重演限制應隨每份預報結果一起交代。
#
# ## 參考資料與延伸閱讀
#
# - Mizrahi 等（2024），[預報發展、檢驗與溝通](https://doi.org/10.1029/2023RG000823)。免費綜述；§2.2.1.4–7 與圖 2 區分前瞻、回溯、擬前瞻及部分前瞻。
#
# - Savran, W. H. 等（2022），[pyCSEP: A Python Toolkit for Earthquake Forecast Developers](https://doi.org/10.1785/0220220033)，*Seismological Research Letters* 93, 2858–2870。從「Core classes」一節讀 Region 與 GriddedForecast，就能對上本章 3.5–3.6 節；套件與文件在 GitHub 免費取得。
# - Jalilian, A.（2019），[ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01)，*Journal of Statistical Software* 88, Code Snippet 1，免費全文。第 2.1 節定義 target／complementary events，第 4 節的 `catalog()` 範例對應本章 3.7 節；不需要會 R 也能讀懂。
# - Schorlemmer, D., Gerstenberger, M. C., Wiemer, S., Jackson, D. D. 與 Rhoades, D. A.（2007），Earthquake likelihood model testing，*Seismological Research Letters* 78(1), 17–29。RELM 實驗的規格文件：網格、規模箱、十欄格式與一致性檢驗的原始出處；出版社全文需訂閱或館藏權限。
# - Gasperini, P., Biondini, E., Lolli, B., Petruccelli, A. 與 Vannucci, G.（2021），[Retrospective short-term forecasting experiment in Italy based on the occurrence of strong (fore) shocks](https://doi.org/10.1093/gji/ggaa592)，*Geophysical Journal International* 225, 1192–1206。177 格研究區的來源；讀方法段可以看到選格規則的原始敘述。
# - Biondini, E., Rhoades, D. A. 與 Gasperini, P.（2023），[Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123)，*Geophysical Journal International* 234, 1681–1700。本章的規格全部來自它的「Application to Italy」一節與附錄 A；建議拿著規格卡逐句對照。
