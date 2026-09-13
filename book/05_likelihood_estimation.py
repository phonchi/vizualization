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
# # 5. 統計工具箱 II：概似與估計
#
# 學習期的目標事件數除以年數，就得到平均年率。
# 這個算式很熟悉，但背後仍有模型假設。
# 如果把率換成另一個數字，資料會有多不相容？
# 本章用同一份 HORUS 學習期回答。
# 最後區分「率估得不準」與「未來本來就會波動」。
#
# ## 5.1 把資料固定，讓參數移動
#
# 第 4 章固定期望數，計算各種計數的機率。
# 現在方向反過來：資料已觀察到，率還未知。
# **參數（parameter）**是控制模型分布的設定值。
# 本章的參數是年率 $\lambda$，單位為顆／年。
#
# ```{admonition} 定義：概似（likelihood）
# :class: definition
# 概似把已觀察的資料固定，將其模型機率視為參數的函數。
# 它用來比較不同參數對同一份資料的支持程度。
# ```
#
# 若學習期長度為 $T$ 年，總數為 $N$，
# Poisson 模型給出
#
# $$L(\lambda;N)=e^{-\lambda T}\frac{(\lambda T)^N}{N!}.$$
#
# 這裡的 $N$ 固定，改變的是橫軸上的 $\lambda$。
# 圖上的高點表示較能解釋這份計數。
# 概似並沒有給「參數為某值」一個機率。
# 把曲線面積縮成 1，也不會自動改變這件事。
# 如果要談參數機率，需要另加推論架構。
# 相關的貝氏觀點放在附錄 A。

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
cat = italy.experiment_catalog()
learn = italy.target_events(cat, "learning")
years = np.arange(1990, 2012)
annual = learn.groupby("year").size().reindex(years, fill_value=0).to_numpy()
T = float(len(years))
N = len(learn)
rate_hat = N / T
rates = np.linspace(.02, max(3.0, rate_hat*2.4), 500)
log_likelihood = poisson.logpmf(N, rates*T)
fig = go.Figure(go.Scatter(x=rates, y=log_likelihood, line_color=ACCENT))
fig.add_vline(x=rate_hat, line_dash="dash", annotation_text=f"最高點：{rate_hat:.3f} 顆／年")
apply_layout(fig, title=f"HORUS 學習期：{T:g} 年共 {N} 顆目標地震", xaxis_title="候選年率（顆／年）",
    yaxis_title="對數概似", height=420)
fig

# %% [markdown]
# **對數概似（log-likelihood）**是概似的自然對數，記為 $\ell$。
# 對數保持大小順序，所以最高點不變。
# 它把很多機率的乘積換成加總，較容易計算。
# 自然對數的底是 $e$；本書的概似都採這個慣例。
#
# $$\ell(\lambda)=N\ln\lambda-\lambda T+N\ln T-\ln N!.$$
#
# 最後兩項不隨 $\lambda$ 改變。
# 若只比較這份資料的不同率，可以暫時省略。
# 若要交付完整分數，就要先說明保留了哪些項。
# 不同資料或不同計分格式之間，不能任意沿用省略。
#
# ## 5.2 最高點就是最大概似估計
#
# **最大概似估計（maximum likelihood estimation，MLE）**選取使概似最大的參數。
# 對率微分後，斜率為 $N/\lambda-T$。
# 讓斜率等於零，就得到
#
# $$\widehat\lambda=\frac{N}{T}.$$
#
# 對這個簡單模型，平均率就是 MLE。
# 帽子表示這是由資料算出的估計值。
# 它不等於宣稱真實率已被完全知道。
# 若學習期沒有任何事件，最高點會落在零。
# 這是參數邊界；常見的大樣本近似會比較脆弱。
#
# 本章把 1990–2011 視為 22 個日曆年單位。
# 正式 SUP 快取則先使用實際天數算率。
# 兩者的換算需保留自己的時間單位。
# 不能拿顆／年的估計，直接乘 91.31 天。
#
# 一條尖的曲線表示偏離最高點後，支持程度很快下降。
# 一條平的曲線表示多個率仍可解釋資料。
# 但曲線外觀不是模型正確性的證據。
# 錯誤的獨立性假設，也能產生尖銳的最高點。
#
# ## 5.3 若知道總數，只想估事件的時間形狀
#
# 假設某段長度固定的期間已知有 $N$ 個事件。
# 現在關心它們偏早發生，還是偏晚發生。
# **機率密度（probability density）**用曲線面積表示連續值落在區間內的機率。
# 密度值可以大於 1，整段積分才必須等於 1。
#
# 令 $f(t;\theta)$ 是觀測窗內正規化的時間密度。
# $\theta$ 控制它的形狀。
# 若給定總數後，各事件時間可視為獨立抽樣，
# 固定總數的對數概似可寫成
#
# $$\ell(\theta\mid N)=\sum_{i=1}^{N}\ln f(t_i;\theta).$$
#
# 若事件依時間排序，密度另有一個 $N!$ 因子。
# 它不影響 $\theta$ 的最佳值。
# 這裡的關鍵是先條件於總數，不再估多少顆。
# 它不同於同時評估時間與總數的完整模型。
#
# 例如原始衰減形狀是 $h(t;\theta)$。
# 若只觀察 $a$ 到 $b$，就必須使用
#
# $$f(t;\theta)=\frac{h(t;\theta)}{\int_a^b h(u;\theta)\,du}.$$
#
# 分母可能也依賴參數，所以估計時不能丟掉。
# 省略它，等於在比較不同總質量的曲線。
# 第 10 章估計餘震衰減，就會用到這一步。
# 若事件仍互相觸發，獨立密度乘積只能作近似。
# 第 11 章會給能隨歷史更新的完整表達。
#
# ## 5.4 三種「不確定」先分清楚
#
# **抽樣分布（sampling distribution）**是重複取得資料時，估計值形成的分布。
# **標準誤（standard error）**是這個分布的標準差。
# 它描述估計值有多容易隨樣本變動。
# 地震計數本身的標準差，描述的是另一個量。
#
# 在固定率 Poisson 模型下，
# $\operatorname{Var}(N)=\lambda T$。
# 將計數除以 $T$，可得
#
# $$\operatorname{Var}(\widehat\lambda)=\frac{\lambda}{T},
# \qquad \widehat{\operatorname{SE}}(\widehat\lambda)=\frac{\sqrt N}{T}.$$
#
# 觀察時間變長，率的標準誤通常變小。
# 但未來一年的計數，仍有平均等於變異數的波動。
# **預測不確定性（predictive uncertainty）**描述未來觀測值的未知。
# 它包含事件隨機性；率未知時還需考慮參數不確定性。
#
# | 看到的範圍 | 回答的問題 | 單位 |
# |---|---|---|
# | 年計數分布 | 某年可能有幾顆？ | 顆 |
# | 率估計的抽樣分布 | 重做估計會得到哪些率？ | 顆／年 |
# | 率的標準誤 | 上述估計值有多分散？ | 顆／年 |
#
# ## 5.5 用 bootstrap 重做估計
#
# **拔靴法（bootstrap）**用重抽資料或擬合模型產生新資料，再重做同一估計。
# 這是一種近似抽樣分布的方法。
# 下面各做 2,000 次，只存小型年度計數陣列。
# 兩條路的差別是「如何產生新資料」。
#
# **非參數 bootstrap（nonparametric bootstrap）**在 22 個年度計數中，有放回地重抽 22 年。
# 某個年份可能重複出現，另一些年份可能沒被抽到。
# 每次都重新算平均率。
# 這保留年度計數的經驗形狀，但假設年度可以交換。
#
# **參數 bootstrap（parametric bootstrap）**則從已擬合的 Poisson 率模擬新總數。
# 每次模擬相同的 22 年長度，再除以 22。
# 它直接繼承 Poisson 的獨立與變異假設。
# 兩種結果不同，可能反映假設不同。

# %% tags=["remove-input"]
rng = np.random.default_rng(20260914)
n_boot = 2000
boot_year = rng.choice(annual, size=(n_boot, len(annual)), replace=True).mean(axis=1)
boot_poisson = rng.poisson(rate_hat*T, size=n_boot) / T
fig = go.Figure()
for values, label, color in [(boot_year, "整年重抽", ACCENT), (boot_poisson, "Poisson 模型重抽", PALETTE[1])]:
    fig.add_trace(go.Histogram(x=values, histnorm="probability", opacity=.55, name=label,
        xbins=dict(start=-.025, end=4, size=.10), marker_color=color))
fig.add_vline(x=rate_hat, line_dash="dash")
fig.update_layout(barmode="overlay")
apply_layout(fig, title=f"{n_boot:,} 次 bootstrap 的率估計", xaxis_title="年率估計（顆／年）",
    yaxis_title="比例", height=420)
fig

# %% tags=["remove-input"]
lo, hi = np.quantile(boot_year, [.025, .975])
plo, phi = np.quantile(boot_poisson, [.025, .975])
display(Markdown(f"整年重抽的中央 95% 百分位區間：**{lo:.2f}–{hi:.2f} 顆／年**。"
    f"Poisson 重抽為 **{plo:.2f}–{phi:.2f} 顆／年**；"
    f"Poisson 解析標準誤為 **{np.sqrt(N)/T:.3f} 顆／年**。"))

# %% [markdown]
# **百分位區間（percentile interval）**取重抽估計值的兩個分位數作端點。
# 這裡取 2.5% 與 97.5%，留下中間約 95%。
# 它是估計方法的近似區間，沒有保證這一次包住真值。
# 在頻率學派的解釋下，涵蓋率談的是重複實驗的表現。
#
# 整年重抽保留一年內的群聚，卻沒有保留跨年相依。
# 若一個序列跨過新年，兩部分可能被分開抽走。
# 只有 22 年，也未必看過所有重要變化。
# 因此，這張圖展示方法及假設，不能當成完整誤差審計。
# 正式研究還需檢查重抽單位與資料的時間結構。
#
# ## 5.6 讀一個估計結果的順序
#
# 假設論文報出一個率與它的區間。
# 先找觀察期間及事件門檻，再讀小數位。
# 若資料只保留規模至少 5 的事件，
# 那個率就不能直接解釋成所有地震的年率。
# 規模門檻改變，估計的對象也跟著改變。
#
# 再找區間的計算方式。
# Poisson 公式、年度重抽與完整序列模擬，
# 可能給出不同寬度，因為保留的變異來源不同。
# 較窄的區間未必代表較好的研究。
# 它也可能來自較強的假設或遺漏了相依性。
#
# 最後檢查資料能否支持模型的精細程度。
# 本章只估一個平均率，所以一個總數就足夠。
# 若要估各地、各時段都不同的率，
# 資料就要支持更多彼此不同的設定。
# 程式成功算出參數，不等於每個參數都估得穩。
# 後面的 ETAS 章會再遇到參數互相補償的例子。
#
# 本章的區間沒有處理目錄漏報與規模誤差。
# bootstrap 重抽已經看見的資料，
# 不會自動找回沒被記錄的地震。
# 第 8 章先檢查目錄完整度，正是要補上這一層。
#
# ## 5.7 本章填入的規格欄位
#
# 這一章沒有更換實驗範圍或預報模型。
# 新增的是讀懂估計值與區間的方法。
# 之後比較模型時，還要問參數由哪一期資料估計。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("估計資料", "學習期 1990–2011；R 內取整規模 ≥ 5.0"),
    ("本章方法", "Poisson 年率 MLE；年度／模型 bootstrap 示範"),
    ("完整度", "待填（第 8 章）"),
    ("正式評分細節", "待填（第 6、16–18 章）"),
])))

# %% [markdown]
# 有了估計與模擬，下一章 {doc}`06_simulation_tests_scores` 就能問：觀測結果是否出乎模型預期？
#
# ## 參考資料與延伸閱讀
#
# - Reinhart（2018），[自激發時空點過程綜述](https://arxiv.org/abs/1708.02647)。免費作者稿；第 3 節把概似、估計與診斷放在同一架構。
# - Daley 與 Vere-Jones（2003），[An Introduction to the Theory of Point Processes, Volume I](https://doi.org/10.1007/b97277)。第 7 章延伸到點過程概似；全文可經館藏取得。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。閱讀參數表時，對照學習資料、固定設定與估計值三種來源。
