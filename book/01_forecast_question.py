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
# # 1. 地震可以「預報」嗎？一份預報長什麼樣
#
# 「未來三個月，這個地區會有幾次較大的地震？」
# 這是一個可以交出答案、事後核對的問題。
# 答案可以是一張地圖，每格放一個數字。
# 本書從這張圖出發，走完一份義大利實驗。
# 你不需要先懂地震學，也不需要寫程式。
# 需要的是基礎統計，以及把問題問清楚的習慣。
#
# ## 1.1 三種常被混在一起的問題
#
# **地震預測（earthquake prediction）**，在此指斷言某個時地會發生某種大小的地震。
# **地震預報（earthquake forecast）**，則為指定範圍提供機率或期望數。
# 兩者都談未來，但交付的答案不同。
# 一份機率預報允許事件發生，也允許不發生。
# 判斷它的品質，需要累積證據。
#
# **地震預警（earthquake early warning）**是在地震已開始後，利用先到的訊號發出警示。
# 它爭取的是較強搖晃抵達前的時間。
# 本書的三個月預報，回答的是另一個問題。
# 把預報圖當成即時警報，會誤讀圖上的數字。
#
# | 問題 | 答案的形式 | 何時取得主要訊號 |
# |---|---|---|
# | 某時某地會不會發生？ | 明確的事件斷言 | 所指事件發生前 |
# | 指定期間可能發生多少？ | 分布、機率或期望數 | 每次發報時已有的資料 |
# | 強烈搖晃快到哪裡？ | 已發生地震的即時警示 | 地震開始後 |
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d01_forecast_pipeline", caption="先讀目錄，再固定規格、製作預報與檢驗。")

# %% [markdown]
# ## 1.2 先認識地圖上的一個點
#
# 一筆地震至少需要時間、位置和大小。
# **規模（magnitude）**衡量地震本身的大小。
# 本書主要使用矩規模 $M_w$。
# 同一場地震在不同地方造成的搖晃可不同。
# 因此，規模不能直接當成你所在地的搖晃程度。
#
# **震源（hypocentre）**是地震破裂開始的位置。
# **震央（epicentre）**是震源在地表的投影。
# 地圖上的點標出震央；**深度（depth）**則描述震源在地下多深。
# 同一經緯度的兩顆點，仍可能有不同深度。
#
# **地震目錄（earthquake catalogue）**是一張逐事件的資料表。
# 本書使用 HORUS 義大利目錄。
# 下一章會解釋它的全名、欄位與篩選。
# 下面先把研究區內較大的地震畫出來。
# 紅點大小對應規模，位置對應震央。
# 密集的點可能重疊，不能只靠目視數事件。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout, show_diagram, plot_forecast_map

# %% tags=["remove-input"]
cat = italy.experiment_catalog()
large = cat.loc[cat.in_R & (cat.mb >= 5.0)]
fig = go.Figure(go.Scatter(x=large.lon, y=large.lat, mode="markers",
    marker=dict(size=4 + 3 * (large.mb - 5), color=QUAKE_COLOR, opacity=.65),
    customdata=np.column_stack([large.time.astype(str), large.mb]),
    hovertemplate="%{customdata[0]}<br>規模 %{customdata[1]:.1f}<extra></extra>"))
apply_layout(fig, title=f"HORUS：1960–2021 年研究區內規模 ≥ 5.0，共 {len(large)} 顆",
    xaxis_title="經度（度）", yaxis_title="緯度（度）", height=480)
fig

# %% [markdown]
# 這張經緯度散點圖用來找事件位置。
# 它沒有呈現地形，也沒有校正成等距地圖。
# 後面的預報會使用投影公里座標分格。
# 同一份目錄可以畫成地圖或時間圖。
# 兩種讀法分別顯示「在哪裡」與「何時」。
#
# ## 1.3 終點預覽：每格放的數字是什麼
#
# **期望數（expected count）**是重複相同條件時的平均次數。
# 某格寫 0.2，不是有五分之一顆地震。
# 它表示多次相同實驗的平均可以是 0.2 顆。
# **機率（probability）**則描述指定事件的可能性。
# 例如「至少一顆」的機率必須介於 0 與 1。
# 期望數可以大於 1；兩者不能直接互換。
#
# 下面是後面會完成的 EEPAS 預報。
# 現在只把 EEPAS 當成模型名字。
# 第 15 章才解釋它如何把歷史事件變成預報。
# 藍色越深，代表每格預報的期望數越高。
# 色階採對數刻度，方便看出小數值的差異。
# 紅點是測試期間實際發生的目標地震。
# **目標地震（target event）**就是規格要求預報、也要計分的事件。

# %% tags=["remove-input"]
from gdms_toolkit.italy_models import get_forecast
forecast = get_forecast("EEPAS", "testing")
targets = italy.target_events(cat)
fig = plot_forecast_map(forecast.sum(axis=(0, 2)), targets=targets,
    title="終點預覽：EEPAS 測試期逐窗期望數合計與目標事件")
fig

# %% [markdown]
# 這張圖加總測試期各窗的預報。
# 它不是在 2012 年初一次發布的十年預報。
# 每個窗可以使用發報前新增的歷史。
# 紅點也是發生之後才疊上去的。
# 本章只疊圖，不評分；評分留到第 16–18 章。
#
# 紅點落在深色格，是值得檢查的線索。
# 但高率若鋪滿全區，也容易包住紅點。
# 只看命中的位置，會忽略空白區的預報。
# **基準模型（baseline model）**是比較時的簡單參照。
# 第 4 章會先做出空間均勻的基準。
# 後面所有複雜模型，都要說明多提供了什麼。
#
# ## 1.4 同樣平均數，可以有不同時間排列
#
# **餘震（aftershock）**指較大地震後，在相關區域出現的後續地震。
# 一串相關的事件稱為**地震序列（earthquake sequence）**。
# 圖上時間接近，只能先描述為群聚。
# 是否存在觸發關係，需要模型與其他證據。
#
# 下面取 Amatrice 周圍 50 公里作示例。
# 起點固定為 2016 年 8 月 24 日零時。
# 範圍是後續 30 天，規模至少 2.5。
# 左圖呈現相鄰事件的時間差。
# 右圖把同一批事件按天計數。
# 這是地域與日期篩選，沒有替事件指定家族。

# %% tags=["remove-input"]
from plotly.subplots import make_subplots
x0, y0 = italy.lonlat_to_km(13.23, 42.70)
t0 = pd.Timestamp("2016-08-24")
seq = cat.loc[(cat.time >= t0) & (cat.time < t0 + pd.Timedelta(days=30))
    & (cat.mb >= 2.5) & (np.hypot(cat.x_km-x0, cat.y_km-y0) <= 50)].sort_values("time")
day = ((seq.time - t0).dt.total_seconds() / 86400).to_numpy()
n_day, edges = np.histogram(day, bins=np.arange(31))
fig = make_subplots(rows=1, cols=2, subplot_titles=("相鄰事件間隔", "每天事件數"))
fig.add_trace(go.Scatter(x=day[1:], y=np.diff(day)*24, mode="markers",
    marker=dict(color=ACCENT, size=4), showlegend=False), row=1, col=1)
fig.add_trace(go.Bar(x=edges[:-1], y=n_day, marker_color=ACCENT, showlegend=False), row=1, col=2)
fig.update_xaxes(title_text="自 8 月 24 日起的天數")
fig.update_yaxes(title_text="間隔（小時）", row=1, col=1)
fig.update_yaxes(title_text="顆／天", row=1, col=2)
apply_layout(fig, title=f"Amatrice 周圍 50 公里、30 天內規模 ≥ 2.5：{len(seq)} 顆", height=400)
fig

# %% [markdown]
# 早期的短間隔與高計數，可以互相對照。
# 但規模較小的地震可能在繁忙時漏掉。
# 所以，計數下降可能同時含有觀測影響。
# 下一章會先檢查目錄，再談地震規律。
#
# 下圖刻意讓三份合成資料有相同總數。
# 每條線都有相同長度，每個刻痕是一個事件。
# 這些刻痕不是 HORUS 資料。
# 中間的隨機撒點，也會自然形成空白與密集段。
# 因此，「不規則」本身還不足以證明群聚機制。

# %% tags=["remove-input"]
from gdms_toolkit.teaching import learning_catalog
kinds = [("regular", "接近等間隔"), ("poisson", "均勻隨機，固定總數"), ("cluster", "分組生成，固定總數")]
fig = make_subplots(rows=3, cols=1, shared_xaxes=True, subplot_titles=[v for _, v in kinds])
for row, (kind, label) in enumerate(kinds, 1):
    demo = learning_catalog(kind)
    fig.add_trace(go.Scatter(x=demo.day, y=np.ones(len(demo)), mode="markers",
        marker=dict(symbol="line-ns", size=14, color=ACCENT, line=dict(width=2, color=ACCENT)),
        showlegend=False), row=row, col=1)
    fig.update_yaxes(visible=False, row=row, col=1)
fig.update_xaxes(title_text="合成時間（天）", row=3, col=1)
apply_layout(fig, title=f"三種合成排列：每份 {len(demo)} 顆，均為 120 天", height=450)
fig

# %% [markdown]
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d20_model_family", caption="SUP 固定學習期率；其餘模型依各自規則更新歷史。")

# %% [markdown]
# ## 1.5 一份能核對的預報實驗
#
# **前瞻實驗（prospective experiment）**在結果尚未出現時，就固定規則並發布預報。
# **擬前瞻實驗（pseudo-prospective experiment）**則在已有目錄上重演發報。
# 每次只讓模型讀到當時之前的資料。
# 兩者都需要清楚記錄資料截止時間。
# 重演仍可能受事後選區或模型選擇影響。
#
# 本書重演義大利 2012–2021 年實驗。
# 研究區沿用已發表論文的選區。
# 該選區使用過 2021 年前的資訊。
# 所以不能把整個設計宣稱為事前凍結。
# 這個限制不妨礙學習，但會限制結論。
#
# **地震可預測性實驗協作平台（Collaboratory for the Study of Earthquake Predictability，CSEP）**推動共同規格與可比較的檢驗。
# 共同規格讓不同模型回答同一題。
# **作業地震預報（Operational Earthquake Forecasting，OEF）**則持續提供及更新預報，供實際決策參考。
# 從研究分數走到決策，還需要溝通與用途評估。
# 第 20 章會回到這個問題。
#
# ## 1.6 從論文帶走三種判斷
#
# 一張圖旁邊寫著「預報表現良好」，先找它回答的題目。
# 預報對象是規模至少 5 的地震，還是所有小地震？
# 期間是一天、三個月，還是十年？
# 這些設定改變，數字的意思也會改變。
# 圖名與圖說因此是閱讀的一部分。
#
# 接著找資料的時間順序。
# 學習用的資料若混入後來的事件，模型就得到額外資訊。
# 相反地，事件發生後拿來核對預報，是正常的評估。
# 兩者差在資料是否參與了當時的預報選擇。
# 第 3 章會把這條界線畫成時間軸。
#
# 最後找支持結論的證據種類。
# 位置看起來接近，只能先談空間配置。
# 總數接近，只能先談數量。
# 一個平均分數較高，也還需要知道比較對象。
# 讀完第一部後，你應能指出這些差別。
# 遇到陌生模型，也能先讀懂它的輸入與交付物。
#
# 第一部以同一份義大利實驗練習上述讀法。
# 第二部再回到儀器、訊號與資料取得。
# 目錄裡的一列，是許多觀測與處理工作的結果。
# 了解這條來源路徑，才能知道模型看不到什麼。
#
# ## 1.7 本章填入的規格欄位
#
# 先保留一張空白卡。下一章從目錄與區域填起。
# 每填一欄，就縮小「我們到底在預報什麼」的歧義。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("目錄、規模尺度、深度", "待填（第 2 章）"),
    ("收集區與測試區", "待填（第 2 章）"),
    ("三段期間、發報窗與門檻", "待填（第 3 章）"),
    ("預報模型與期望數", "待填（第 4 章起）"),
    ("評分方式", "待填（第 6、16–18 章）"),
])))

# %% [markdown]
# 帶著地圖上的一個點，進入 {doc}`02_reading_horus`。
# 先知道資料怎麼來，才知道哪些數字值得相信。
#
# ## 參考資料與延伸閱讀
#
# - CSEP，[官方入門網站](https://cseptesting.org/)。免費認識共同預報格式與測試的目的。
# - Lolli 等（2020），[HORUS 義大利均一化目錄](https://doi.org/10.1785/0220200148)。了解規模整合及目錄年代差異。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。本書的實驗來源；先看研究區與預報圖，再逐章返回方法。
