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
# # 20. 預報之後：讀一篇論文、進入決策，與回到儀器
#
# 一張格網預報，已經走過三類檢驗。
# 第 19 章還把兩張預報加成一張。
# 本章先拿規格卡去讀原論文。
# 再追問：格內的地震期望數，如何連到地面搖動？
#
# 一列儀器記錄，是整條工作線的起點。
# 從波形到目錄，再到預報與決策，
# 每一步都會改變數字所代表的東西。
# 本章沿這條線回望第一部。
# 最後把待查欄位交給第二部。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats
from plotly.subplots import make_subplots
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, apply_layout, show_diagram

# %% [markdown]
# ## 20.1 在方法段旁邊，寫下規格欄位
#
# 1990–2011 與 2012–2021 是兩段明確的日期。
# 先在紙上圈出這些數字，再讀作者的理由。
# 旁註（annotation）是緊貼原句的閱讀說明。
# 每條旁註都指出設定、證據或尚待查明之處。
# 以下使用摘述，只保留必要的英文短語。
#
# [Biondini 等（2023）](https://doi.org/10.1093/gji/ggad123)
# 在 APPLICATION TO ITALY 節寫出學習與測試分工。
# 下面選讀該節第二段，共兩句。
# 右欄將每句接回本站的規格卡。
# 兩句必須連著讀，才能看見資料如何分工。
#
# | 方法句摘述 | 逐句旁註 |
# |---|---|
# | 作者用 1990–2011 擬合 EEPAS。理由是目錄品質改善。短語是「accuracy and completeness」。 | 先填學習期間。完整度指目錄可靠收錄的規模範圍。英文為 completeness。這個理由需連到門檻證據。 |
# | 作者用 2012–2021 作擬前瞻測試。 | 再填測試期間。擬前瞻（pseudo-prospective）是保留時間因果的回溯檢驗子類。核對模型與建模者可見的資料，不只看期間是否分開。 |
#
# 「獨立測試期」先描述資料用途的分離。
# 它不能單獨證明所有研究設計都事前固定。
# 同節的研究區選擇用到 2021 年前的資訊。
# 因此，本站稱為既有研究區的重演。
# 這個旁註把時間分工與選區證據放回同一頁。
#
# 177 格是已採用的研究區結果。
# 格子的選法會影響哪些事件進入評估。
# 讀者應另查深度、規模與目錄版本。
# 本站測試目錄為 25 顆，原論文為 27 顆。
# 日期相同，仍可能因版本而得到不同計數。
#
# ## 20.2 在結果圖旁邊，寫下分數的分母
#
# Fig. 8 比較不同發報窗長的模型表現。
# 下圖保留原論文的座標、模型標籤與區間。
#
# [![Biondini 2023 Fig. 8：五種預報窗長下各模型相對 SUP 的每活動箱資訊增益與信賴區間](_static/reading/biondini2023_fig8.png)](_static/reading/biondini2023_fig8.png)
#
# 點圖可開啟大圖，放大比較模型標籤與信賴區間。
#
# 圖：Biondini、Rhoades 與 Gasperini（2023），Fig. 8，p. 1692。
# 來源：[原論文](https://doi.org/10.1093/gji/ggad123)。黑線為 95%，灰線為 97.5% 信賴區間。
# 先從最左的三個月窗開始讀，再往右比較窗長。
# 原圖使用每活動箱資訊增益，簡稱 IGPA。
# 英文為 information gain per active bin。
# 活動箱（active bin）指至少有一顆事件的格箱。
#
# 一箱若有兩顆，IGPA 仍只算一個活動箱。
# 第 18 章的 IGPE 則按事件數平均。
# 原論文還使用二元分數作比較。
# 因此，兩種增益的分母與評分都須核對。
# 本站數字不能直接填回原圖。
#
# | 結果描述摘述 | 逐句旁註 |
# |---|---|
# | 正文指出，最短窗長「3 months」由 ETAS 系列表現較好。 | 先圈窗長，再圈模型版本。這是作者那份測試目錄的結果。本站 ETAS 採第一代近似。 |
# | 正文接著指出，較長窗長由 EEPAS 系列表現較好。 | 結論隨窗長改變。它沒有給出所有情境共用的冠軍。 |
# | 圖說以 SUP 為零點，並畫出兩組信賴界限。 | 先看零點，再看區間。基準、區間水準與多重判斷要一起讀。 |
#
# 一個高於零的中心點，先表示相對分數較高。
# 區間是否跨零，才影響差異的判讀。
# 原圖的黑、灰線對應不同信賴水準。
# 旁註對照上方原圖，沒有重畫或估讀數值。
# 讀者可以沿章末免費全文逐項核對。
#
# ## 20.3 在結論旁邊，留下適用範圍
#
# 結論段把 EEPAS 的優勢連到較長預報窗。
# 同一句也提醒目標樣本有限。
# 短語「small number of target shocks」保留了這層語氣。
# 摘述是：EEPAS 可能較適合較長期預報。
# 小樣本使這個比較仍需審慎看待。
#
# 這一句的旁註應保留「可能」與「較長期」。
# 前者對應證據強度，後者對應使用情境。
# 若只抄「EEPAS 較好」，兩項條件都消失了。
# 把句子變短，仍要保留作者真正比較的範圍。
# 這就是從原文走到可讀教材的關鍵。
#
# 一篇論文可以用三種記號來讀。
# 日期與門檻旁邊寫「設定」。
# 圖表與分數旁邊寫「證據」。
# 結論的情境與語氣旁邊寫「範圍」。
# 三者能互相對上，才形成完整的論證。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d20_model_family", caption="時間尺度與更新方式是閱讀模型論文的兩個線索。")

# %% [markdown]
# ## 20.4 同一顆地震，各場址搖動不同
#
# 一格預報「可能有地震」，還沒描述某棟樓的搖動。
# 地動（ground motion）是場址的地面運動。
# 地震到場址的距離會改變它。
# 規模與地下介質也會改變它。
# 所以事件發生率還需要一個地動轉換。
#
# 地動預估式稱為 GMPE。
# 英文是 ground-motion prediction equation。
# 它按規模、距離與場址條件給出地動分布。
# 同樣條件下，地動仍會有散布。
# 預估式因此提供分布，而非單一確定搖動值。
#
# 一個場址的地動門檻可記作 $a$。
# 地動危害（seismic hazard）描述超過此門檻的可能性。
# 機率式地震危害分析簡稱 PSHA。
# 英文是 probabilistic seismic hazard analysis。
# 它把可能地震與地動分布合計起來。
#
# Baker 的入門教材把這件事拆成五步。
# 先辨認能影響場址的震源。
# 再描述各規模的發生率。
# 接著描述震源到場址的距離分布。
# 第四步給出條件地動分布。
# 第五步把大小、位置與地動的不確定性積分。

# %% [markdown]
# $$\nu(a)=\int\!\!\int
# \lambda(m,r)P(Y>a\mid m,r)\,\mathrm{d}r\,\mathrm{d}m.$$
#
# $Y$ 是場址地動，$r$ 是震源距離。
# $\lambda(m,r)$ 是每年每單位規模與距離的事件率。
# $\nu(a)$ 是每年的地動超越率。
# 超越率（exceedance rate）表示每年超過門檻的期望次數。
# 這個年率式假定評估期間的發生率可視為固定。
#
# 下圖只用兩個合成震源，展示合計的樣子。
# 峰值地表加速度稱為 PGA。
# 英文是 peak ground acceleration。
# 橫軸用重力加速度 $g$ 作單位。
# 曲線愈往右，要求的搖動門檻愈高。
# 合成參數用於教學，沒有代入義大利場址。

# %% tags=["remove-input"]
pga = np.geomspace(0.01, 1.5, 150)
annual_rates = np.array([0.02, 0.005])
medians = np.array([0.12, 0.30])
log_sigma = 0.55
contributions = np.array([rate*stats.lognorm.sf(pga, s=log_sigma, scale=median)
                          for rate, median in zip(annual_rates, medians)])
nu = contributions.sum(axis=0)
assert np.all(np.diff(nu) <= 0)
fig = go.Figure()
for i, values in enumerate(contributions):
    fig.add_scatter(x=pga, y=values, name=f"合成震源 {i+1}", line=dict(dash="dot", color=["#008300", "#4a3aa7"][i]))
fig.add_scatter(x=pga, y=nu, name="合計", line=dict(color=ACCENT, width=3))
apply_layout(fig, title="合成地動危害曲線", xaxis_title="PGA 門檻（g）", yaxis_title="每年超越率", xaxis_type="log", yaxis_type="log")
display(Markdown(f"圖：{len(annual_rates)} 個合成震源；年率 {annual_rates[0]:.3f}、{annual_rates[1]:.3f}。地動中位數 {medians[0]:.2f}、{medians[1]:.2f} g；對數標準差 {log_sigma:.2f}。"))
fig

# %% [markdown]
# ## 20.5 等待時間的 hazard，與地動危害分開讀
#
# $T$ 表示下一次指定事件的等待時間。
# 存活函數（survival function）寫作 $S(t)=P(T>t)$。
# 它表示等到 $t$ 仍未發生的機率。
# 危害函數（hazard function）則是已等到此刻後的瞬時率。
# 兩個英文都有 hazard，指的量卻需看上下文。
#
# $$h(t)=\frac{f_T(t)}{S(t)},\qquad
# P(t<T\le t+\mathrm{d}t\mid T>t)\approx h(t)\mathrm{d}t.$$
#
# $f_T$ 是等待時間的機率密度。
# $h(t)$ 的單位是時間的倒數。
# 它描述指定事件在尚未發生條件下的發生率。
# 前節的地動超越率則針對場址門檻。
# 兩者不能只因中文都譯作危害就互換。
#
# 一條常數 $h(t)$ 對應指數等待時間。
# 一條上升的 $h(t)$ 則表示等愈久，條件率愈高。
# 下圖是兩種合成分布的對照。
# 它不替任何真實斷層選定複發模型。
# 完整推導可接著讀附錄 F。

# %% tags=["remove-input"]
time_years = np.linspace(0, 200, 150)
scale = 100.0
survival_exp = np.exp(-time_years/scale)
survival_increasing = np.exp(-(time_years/scale)**2)
fig = make_subplots(rows=1, cols=2, subplot_titles=["仍未發生的機率", "已等到此刻的瞬時率"])
for name, surv, hazard, color in [
    ("常數瞬時率", survival_exp, np.full_like(time_years, 1/scale), ACCENT),
    ("上升瞬時率", survival_increasing, 2*time_years/scale**2, "#4a3aa7")]:
    fig.add_trace(go.Scatter(x=time_years, y=surv, name=name, legendgroup=name, line=dict(color=color)), row=1, col=1)
    fig.add_trace(go.Scatter(x=time_years, y=hazard, name=name, legendgroup=name, showlegend=False, line=dict(color=color)), row=1, col=2)
fig.update_xaxes(title_text="已等待年數")
fig.update_yaxes(title_text="S(t)", row=1, col=1)
fig.update_yaxes(title_text="h(t)（每年）", row=1, col=2)
apply_layout(fig, title="合成等待時間：機率與瞬時率", legend=dict(orientation="h", y=-0.24), margin=dict(l=60, r=30, t=80, b=100))
display(Markdown(f"圖：兩條合成分布都使用尺度 {scale:.0f} 年；尺度參數在不同分布中未必等於平均等待時間。"))
fig

# %% [markdown]
# $1/\nu(a)$ 稱為回歸期（return period）。
# 它把長期超越率換成平均間隔的尺度。
# 例如年率越小，回歸期就越長。
# 它仍容許短時間內連續超越門檻。
# 若超越次數服從固定率 Poisson，
# 期間 $\Delta t$ 的機率為 $1-e^{-\nu(a)\Delta t}$。
#
# 一個場址的損失還需要建築與暴露資訊。
# 風險（risk）描述事件帶來的損失可能性。
# 同樣的地動，在不同建築可能造成不同損失。
# 因此，危害曲線還沒給出全部決策答案。
# 行動成本與能減少的損失也必須具體列出。
#
# ## 20.6 預報表如何走到使用者手上？
#
# USGS 的 OAF 頁面同時呈現窗長與規模門檻。
# 作業型餘震預報稱為 OAF。
# 英文是 operational aftershock forecasting。
# 它會隨序列觀測更新預報。
# 讀表時要一起保留發布時間、區域與有效期間。
# 細節見 [USGS 官方說明](https://earthquake.usgs.gov/data/oaf/overview.php)。
#
# 下表模仿「選定門檻、比較窗長」的讀法。
# 所有期望數都是合成值。
# 時間窗彼此包含，所以較長窗機率較高。
# 這張表僅示範格式，沒有提供官方預報。
# 它也沒有使用任何美國或第二部的資料。

# %% tags=["remove-input"]
durations = np.array([1, 7, 30, 365])
synthetic_expected = 0.08*np.log1p(durations/0.5)
synthetic_probability = -np.expm1(-synthetic_expected)
display(Markdown(f"OAF 讀表示意：固定合成規模門檻 Mw≥{5.0:.1f}；合成期望數採 0.08 ln(1+t/0.5)。"))
display(HTML(pd.DataFrame({"發布後天數": durations, "合成期望數": synthetic_expected,
    "至少一顆機率": [f"{x:.1%}" for x in synthetic_probability]}).to_html(index=False, float_format=lambda x:f"{x:.3f}")))

# %% [markdown]
# OEF-Italy 則提供另一種制度案例。
# 作業型地震預報簡稱 OEF。
# 英文是 operational earthquake forecasting。
# 它持續更新地震機率，供社會使用。
# Mizrahi 等人的 2024 年回顧記錄了義大利系統。
# 該系統把多個成分模型結合成預報。
# 這裡引用出版時的描述，未宣稱今日配置相同。
#
# 一個預報訊息應先讓使用者找到自己。
# 哪個地區、哪段時間、哪個規模門檻？
# 接著才讀機率、期望數與更新時間。
# 如果行動需要地動門檻，
# 資訊提供者就應說明如何從事件率轉換。
# 第 19 章的混合率本身還沒有完成這項工作。
#
# ## 20.7 本章填入的規格欄位
#
# 一張義大利卡換成另一地區的卡，
# 最先要重查的是儀器、目錄與尺度。
# 第二部會回到觀測來源，逐步補齊資料條件。
# 下卡將轉用台灣時的欄位全部標為待定。
# 章號標示首次接手位置，後續章節再補細節。
# 這些空欄也避免把義大利參數直接搬走。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("儀器與資料來源", "待定；待填（第 21 章）"),
    ("目錄與規模尺度", "待定；待填（第 29 章）"),
    ("深度、收集區 S 與測試區 R", "待定；待填（第 29 章）"),
    ("完整度與輸入、目標門檻", "待定；待填（第 29 章）"),
    ("三期、發報窗與資料延遲", "待定；待填（第 29 章）"),
    ("格箱、模型參數與比較基準", "待定；待填（第 29 章）"),
    ("地動、場址與使用情境", "待定；待填（第 29 章）"),
], title="義大利卡換台灣卡：待查欄位")))

# %% [markdown]
# 一個機率數字，最後仍要回到可查的觀測。
# 第二部從儀器與測網重新接起這條線。
# 接著讀：{doc}`21. 回到儀器與觀測資料 <01_overview>`。
#
# 回到第 3 章的分類，本站是修訂目錄上的擬前瞻教學重演。
# 2024 年 HORUS 與含事後資訊的研究區，不能還原成當時的即時前瞻紀錄。
# 時間相依模型、滾動發報與作業用途，都不會改變這項限制。
# 下一步若設計作業實驗，應先鎖定資料版本、更新規則與保存機制，
# 再明確選擇即時前瞻或受控延遲前瞻的安排。
#
# ## 參考資料與延伸閱讀
#
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。
#   [免費機構全文](https://www.earth-prints.org/handle/2122/17084)。先對照研究區與測試期，再讀檢驗。
# - Baker（2013），[PSHA 入門教材](https://scits.stanford.edu/sites/g/files/sbiybj22081/files/media/file/baker_2013_intro_psha_v2_0.pdf)。
#   免費全文。依五步流程閱讀危害積分與回歸期。
# - NIST，[存活函數與危害函數](https://itl.nist.gov/div898/handbook/eda/section3/eda362.htm)。
#   免費手冊。核對密度、存活機率與瞬時率的關係。
# - USGS，[餘震預報說明](https://earthquake.usgs.gov/data/oaf/overview.php)。
#   免費官方頁面。對照時間窗、機率與預報表。
# - Mizrahi 等（2024），[預報發展、檢驗與溝通](https://doi.org/10.1029/2023RG000823)。
#   免費回顧。從義大利制度追讀模型與使用者的連結。
