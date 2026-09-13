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
# # 4. 統計工具箱 I：計數、Poisson 過程，與第一張預報 SUP
#
# 一張預報的某格寫著 0.01，單位是「顆」。
# 這個小數如何對應未來的整數事件？
# 上一章定好格子、期間與規模箱。
# 本章把率轉成期望數，再接上計數分布。
# 最後做出第一張可供比較的預報。
#
# ## 4.1 先把每一年完整數出來
#
# 下圖只取測試區內、規模至少 5.0 的事件。
# 年份從 1960 排到 2021，沒有事件的年也保留。
# 若只留下有地震的年份，平均值就會偏高。
# 這個零值處理，和讀取有事件的列同樣重要。

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
from scipy.stats import poisson
from plotly.subplots import make_subplots
from gdms_toolkit.italy_models import get_forecast
cat = italy.experiment_catalog()
years = np.arange(1960, 2022)
large = cat.loc[cat.in_R & (cat.mb >= 5.0)]
annual = large.groupby("year").size().reindex(years, fill_value=0).to_numpy()
mean_count, var_count = annual.mean(), annual.var(ddof=1)
fig = go.Figure(go.Bar(x=years, y=annual, marker_color=ACCENT))
fig.add_hline(y=mean_count, line_dash="dash", annotation_text="全期年平均")
apply_layout(fig, title=f"HORUS 年計數：均值 {mean_count:.2f}，樣本變異數 {var_count:.2f}",
    xaxis_title="年份", yaxis_title="規模 ≥ 5.0 的事件數", height=400)
fig

# %% [markdown]
# 這條虛線描述整份資料的平均。
# 它不是預報模型的學習期估計。
# 圖中包含測試期，只用來認識計數波動。
# 製作正式預報時，不能把這條全期平均拿去用。
#
# ```{admonition} 定義：Poisson 分布（Poisson distribution）
# :class: definition
# Poisson 分布用一個非負期望數 $\Lambda$，指定整數計數 $N$ 的機率。
# 它的平均與變異數都等於 $\Lambda$。
# ```
#
# $$P(N=n)=e^{-\Lambda}\frac{\Lambda^n}{n!},
# \qquad n=0,1,2,\ldots.$$
#
# 式中的 $n!$ 是階乘，且 $0!=1$。
# 固定 $\Lambda$ 後，每個整數都有對應機率。
# 即使平均是 1.2，單年仍只能觀察整數。
# 下面把實際年計數的比例和 Poisson 比較。
# 為了對照形狀，Poisson 使用全期平均。
# 這是描述性比較，不是正式測試期檢驗。

# %% tags=["remove-input"]
k = np.arange(int(annual.max()) + 5)
fig = go.Figure(go.Bar(x=k, y=np.bincount(annual, minlength=len(k))/len(annual),
    name="HORUS 各計數所占比例", marker_color=ACCENT))
fig.add_trace(go.Scatter(x=k, y=poisson.pmf(k, mean_count), mode="lines+markers",
    name="同平均 Poisson", line_color=PALETTE[1]))
apply_layout(fig, title="相同平均數，計數分布仍可能不同", xaxis_title="一年內的事件數",
    yaxis_title="比例／機率", height=390)
fig

# %% [markdown]
# **過度離散（overdispersion）**指計數變異大於模型的基準變異。
# 對固定平均的 Poisson，基準就是平均本身。
# 實際樣本的均值與變異數不會完全相等。
# 因此，兩個估計值略有差異並不奇怪。
# 若差異很大，固定率假設就值得檢查。
#
# 年與年的真實率不同，也能增加變異。
# 事件群聚、觀測條件改變，都可能參與其中。
# 這張圖不能單獨證明哪一種原因。
# 第 16 章會正式處理計數變異的選擇。
#
# ## 4.2 分布回答多少，過程還回答何時
#
# **Poisson 過程（Poisson process）**描述事件隨時間出現的模型。
# 互不重疊時間段的計數彼此獨立。
# **齊次（homogeneous）**表示每單位時間的率固定。
# 若率為 $\lambda$ 顆／天，長度 $T$ 天的期望數是 $\Lambda=\lambda T$。
#
# **非齊次（nonhomogeneous）**表示率是預先給定的時間函數。
# 例如每天有不同的 $\lambda(t)$。
# 只要函數不由過程中新發生的事件改寫，仍可具有獨立增量。
# 第 11 章才處理事件發生後更新率的模型。
# 不要把「率會變」都當成同一種過程。
#
# $$\Lambda=\int_{t_1}^{t_2}\lambda(t)\,dt.$$
#
# 率曲線下面的面積，才是期間內的期望數。
# 峰值高但持續很短，總期望數仍可能小。
# 積分也能延伸到面積與規模。
# 若密度的單位是顆／天／平方公里／規模單位，
# 對指定窗、格、箱積分後，就得到顆數。
#
# $$\Lambda_{jk}=\int_{t_1}^{t_2}\int_{A_j}\int_{B_k}
# \lambda(t,x,y,m)\,dm\,dA\,dt.$$
#
# $A_j$ 是第 $j$ 個空間格，$B_k$ 是第 $k$ 個規模箱。
# 密度值本身不能直接拿來當格內機率。
# 改變格子大小時，積分範圍也必須改變。

# %% tags=["remove-input"]
show_diagram("d04_density_to_counts", caption="先對窗、格與規模箱積分，再為格內計數指定分布。")

# %% [markdown]
# ## 4.3 從期望數換成至少一次的機率
#
# Poisson 模型下，零次的機率是 $e^{-\Lambda}$。
# 所以至少一次的機率為
#
# $$P(N\ge1)=1-e^{-\Lambda}.$$
#
# 取 $\Lambda=0.2$，這個機率約為 0.181。
# 取 $\Lambda=2$，機率約為 0.865。
# 期望數雖然能超過 1，機率仍不會超過 1。
# 只有期望數很小時，兩個數才會相近。
#
# 同一個期望數，也可以搭配別的計數分布。
# 那些分布的零次機率未必是 $e^{-\Lambda}$。
# 因此，上式包含了 Poisson 假設。
# 地圖只交出期望數，還沒有交出全部機率規則。
#
# ## 4.4 讓模型重複產生資料
#
# **蒙地卡羅模擬（Monte Carlo simulation）**用重複隨機抽樣近似模型的結果分布。
# 先固定率，再模擬許多年的計數。
# 模型不會每次都交出平均值。
# 抽樣次數增加，模擬比例通常更接近理論機率。
#
# **亂數種子（random seed）**決定可重現的抽樣序列。
# 固定種子方便核對圖形，不會消除隨機性。
# 下面保留一千次計數，作為模型計數分布的模擬樣本。
# 它們是模型生成的資料，不是新增的歷史紀錄。

# %% tags=["remove-input"]
rng = np.random.default_rng(20260913)
simulated = rng.poisson(mean_count, size=1000)
fig = go.Figure(go.Histogram(x=simulated, histnorm="probability", xbins=dict(start=-.5, size=1),
    name="1,000 次模擬", marker_color=ACCENT))
fig.add_trace(go.Scatter(x=k, y=poisson.pmf(k, mean_count), name="理論機率", line_color=PALETTE[1]))
apply_layout(fig, title="固定同一平均數，重複產生一年計數", xaxis_title="模擬年計數",
    yaxis_title="比例／機率", height=380)
fig

# %% [markdown]
# **經驗分位數（empirical quantile）**是排序後對應指定累積比例的值。
# 例如 95% 分位數，表示約九成五樣本不超過該值。
# 離散資料會有很多同值，比例可能跨過門檻。
# 第 6 章會把這件事連回尾端機率。
#
# ## 4.5 做出空間均勻 Poisson 預報
#
# **空間均勻 Poisson（Spatially Uniform Poisson，SUP）**把相同面積分到相同期望數。
# 模型先用學習期的目標數估計總率。
# 再依期間長度、空間面積與規模比例分配。
# 177 格面積相同，所以各格的分配比例都是 $1/177$。
#
# 規模比例先採用已發表的 $b=1.084$。
# **Gutenberg–Richter 律（GR law）**描述大規模事件的數量隨門檻提高而下降。
# 第 7 章會解釋與估計 $b$，本章先把它當作給定設定。
# 以 $\pi_k$ 表示第 $k$ 箱的 GR 機率質量。
# 它以規模至少 5 的全部事件作分母。
# 本站保留上方尾端的質量，沒有在 7.5 截止處重新正規化。
# 因此，25 箱的權重合計略小於 1。
#
# $$\Lambda_{jk}=\widehat\lambda\,\Delta t\frac{1}{177}\pi_k,
# \qquad \sum_k\pi_k=1-10^{-b(7.5-5.0)}.$$
#
# $\widehat\lambda$ 必須與 $\Delta t$ 使用相同時間單位。
# 工具箱按天計算學習期長度，避免閏年換算誤差。
# 格與規模箱全部加總，得到規模落在輸出範圍的期望數。
# 它比所有規模至少 5 的期望數略小。
# 不同的尾端處理會改變總量，必須先交代再比較。
# 下一章會解釋為何用「數量除以時間」估計率。

# %% tags=["remove-input"]
sup = get_forecast("SUP", "testing")
windows = italy.forecast_windows()
targets = italy.target_events(cat)
first = targets.loc[(targets.t_days >= windows.iloc[0].t1) & (targets.t_days < windows.iloc[0].t2)]
fig = plot_forecast_map(sup[0].sum(axis=1), targets=first,
    title=f"SUP 第一窗：總期望 {sup[0].sum():.3f} 顆；只疊圖、不評分")
fig

# %% tags=["remove-input"]
fig = plot_forecast_map(sup.sum(axis=(0, 2)), targets=targets,
    title=f"SUP 測試期逐窗合計：總期望 {sup.sum():.2f} 顆；只疊圖、不評分")
fig

# %% [markdown]
# 兩圖都呈現均勻空間分配。
# 色階由各圖數值決定，不能只靠深淺比較總量。
# 請先讀圖名的期望數與預報期間。
# 紅點只協助確認位置與範圍。
# 數量、位置與規模的正式檢驗留到第 16–18 章。
#
# ## 4.6 均勻基準保留什麼、捨去什麼
#
# SUP 保留學習期的總率與指定規模比例。
# 它捨去歷史事件在空間上的不均勻分布。
# 新發生的事件也不會改變它的固定率。
# 因此，SUP 很容易解釋，也容易重做。
# 如果另一個模型勝過它，還要問優勢來自哪個部分。
# 可能是總量更接近，也可能是位置安排更好。
#
# 把某一個格拆成四個等面積小格，
# 在均勻模型下，每個小格得到原本四分之一的期望數。
# 四格加總仍回到原值，這叫期望數的可加性。
# 若改成至少一次的機率，則不能直接相加。
# 兩格都發生事件的情況，會被重複計算。
# Poisson 獨立格下，可先加期望數再換算聯集機率。
#
# 同樣地，十年合計圖也失去發生的時間排列。
# 某一窗特別高、其他窗低的模型，
# 可能與每窗固定率的模型有相同十年總量。
# 所以總量相同，並不代表預報相同。
# 後面的滾動模型會讓這個差別變得具體。
# 閱讀時應同時查看一窗圖與逐窗合計。
#
# 本書的四十個窗每窗長 91.31 天。
# 它們合計 3652.4 天，止於 2021 年 12 月 31 日 09:36。
# 到該年結束仍有 0.6 天未涵蓋。
# 這段尾端沒有符合條件的目標地震。
# 因此，本版目標總數沒有受影響。
# 期望數仍依實際窗長計算，不補成整十年。
# 這種小差異應由規格交代，不靠圖名猜測。
#
# 「Poisson」同時牽涉計數分布與獨立增量。
# 只有平均數的預報檔，無法證明兩項都成立。
# 看到預報檔中的 RATE 欄，先確認它是期望顆數。
# 看到論文說採 Poisson 評分，再確認其分布假設。
# 把輸出與假設分開讀，才能公平比較不同模型。
#
# ## 4.7 本章填入的規格欄位
#
# SUP 的預報規則至此齊全，評分細節仍待後續章節。
# 這張卡讓新模型只需說明改了哪些欄位。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("目錄／區域／深度", "HORUS；S=CPTI15、R=177 格；深度 ≤ 40 km"),
    ("三期", "1960–1989 前置；1990–2011 學習；2012–2021 測試"),
    ("窗與資料截止", "40 個 91.31 天窗；各窗只用發報前歷史，SUP 率固定於學習期"),
    ("門檻", "輸入 m0=2.45；目標取整規模 ≥ 5.0；完整度待第 8 章"),
    ("輸出", "每窗 × 177 格 × 25 規模箱；箱寬 0.1、範圍 [5.0,7.5)"),
    ("SUP 模型", "Poisson 計數；學習期率；等面積均分；GR b=1.084"),
    ("評分", "待填（第 6、16–18 章）"),
])))

# %% [markdown]
# 下一章 {doc}`05_likelihood_estimation` 追問：這個率怎麼估，還有多少不確定性？
#
# ## 參考資料與延伸閱讀
#
# - MIT OpenCourseWare，[Discrete Stochastic Processes，第 2 章](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/resources/mit6_262s11_chap02/)。免費講義，從計數與等待時間理解 Poisson 過程。
# - Reinhart（2018），[自激發時空點過程綜述](https://arxiv.org/abs/1708.02647)。免費作者稿；第 1–2 節比較外部率變化與事件觸發。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。對照 SUP 基準與本站共用的研究規格。
