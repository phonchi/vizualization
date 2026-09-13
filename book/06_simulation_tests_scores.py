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
# # 6. 統計工具箱 III：模擬、檢定與分數
#
# 如果模型預期十年約十多顆，最後卻觀察到更多顆，
# 這個差距只是隨機波動，還是模型需要重看？
# 上一章學會用資料估計率。
# 本章固定那個率，讓模型生成許多可能的未來。
# 再把真正觀測的位置，放進模擬分布中。
#
# ## 6.1 先說清楚要檢查哪一句話
#
# **虛無假設（null hypothesis）**是檢定暫時採用的資料生成規則。
# 本章的假設是：測試期間總數服從 Poisson 分布，
# 其期望數由學習期估計率與測試長度決定。
# 這包含分布形式與數值，不能只寫「預報正確」。
#
# **統計量（test statistic）**是用資料計算、供檢定比較的摘要。
# 本章選整個測試期間的事件總數 $N$。
# 它回答數量是否相容，沒有檢查位置與規模配置。
# 同樣總數的兩張地圖，可能給出完全不同的位置。
#
# **顯著水準（significance level）**是事前選定的拒絕門檻，記為 $\alpha$。
# 本章示範 $\alpha=0.05$。
# 它控制假設成立時，檢定錯誤拒絕的頻率上限。
# 離散計數下，實際頻率可能小於這個上限。
#
# 下面讀取 SUP 快取的測試期總期望數。
# 快取使用第 3 章同一批發報窗。
# 觀測也按相同窗、格與規模箱計數。
# 這讓分子分母都回答同一份規格。

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
from gdms_toolkit.italy_models import get_forecast
cat = italy.experiment_catalog()
rate = get_forecast("SUP", "testing")
targets = italy.target_events(cat, "testing")
omega = italy.bin_targets(targets)
expected, observed = float(rate.sum()), int(omega.sum())
rng = np.random.default_rng(20260915)
n_sim = 3000
simulated = rng.poisson(expected, size=n_sim)
fig = go.Figure(go.Histogram(x=simulated, histnorm="probability", xbins=dict(start=-.5, size=1),
    marker_color=ACCENT, name="模型模擬總數"))
fig.add_vline(x=observed, line_color=QUAKE_COLOR, annotation_text=f"實際總數 {observed}")
apply_layout(fig, title=f"SUP 期望 {expected:.2f} 顆：{n_sim:,} 次計數模擬",
    xaxis_title="測試窗合計事件數", yaxis_title="模擬比例", height=420)
fig

# %% [markdown]
# 圖上的每次抽樣代表整組測試窗的總數。
# 它不是從實際地震中重抽幾顆。
# 3000 次模擬只是近似計數分布。
# 模擬還把估計率當成固定已知。
# 上一章的參數不確定性，尚未被整合進來。
#
# ## 6.2 從排序位置走到尾端機率
#
# **分位數分數（quantile score）**描述觀測統計量在模擬分布中的位置。
# 本站採用「模擬值小於或等於觀測值」的比例：
#
# $$\widehat q=\frac{\#\{N^{(r)}\le N_{\rm obs}\}}{B}.$$
#
# $B$ 是模擬次數，$N^{(r)}$ 是第 $r$ 次結果。
# $q$ 接近 1，表示觀測值偏大。
# $q$ 接近 0，表示觀測值偏小。
# 它不是「模型正確的機率」。
#
# **p 值（p-value）**是在虛無假設下，出現至少同樣極端結果的機率。
# 「極端」必須先由檢定問題定義。
# 只問是否太多，是上尾問題。
# 只問是否太少，是下尾問題。
# 同時關心兩個方向，則是雙尾問題。

# %% tags=["remove-input"]
show_diagram("d06_quantile_pvalue", caption="先決定哪一端算極端，再把觀測值的排序位置換成尾端機率。")

# %% [markdown]
# Poisson 計數是整數，有許多同值。
# 本章用包含觀測值的兩個尾端：
#
# $$p_{\rm upper}=P(N\ge N_{\rm obs}),\qquad
# p_{\rm lower}=P(N\le N_{\rm obs}).$$
#
# 若 $q=P(N\le N_{\rm obs})$，
# 上尾是 $1-q+P(N=N_{\rm obs})$。
# 直接寫成 $1-q$，會漏掉觀測值那一格。
# 兩個含等號尾端的和，因此可能大於 1。
#
# 本章採用的雙尾規則是
#
# $$p_{\rm two}=\min\{1,\,2\min(p_{\rm lower},p_{\rm upper})\}.$$
#
# 離散分布還有其他雙尾定義，結果未必相同。
# 所以報 p 值時，必須連規則一起報。
# 這個等尾規則與第 16 章的兩端檢查銜接。
# 我們同時列模擬比例與解析機率供對照。

# %% tags=["remove-input"]
q_hat = np.mean(simulated <= observed)
upper_hat = np.mean(simulated >= observed)
lower = float(poisson.cdf(observed, expected))
upper = float(poisson.sf(observed-1, expected))
p_two = min(1.0, 2*min(lower, upper))
display(pd.DataFrame({
    "量": ["下尾／分位數分數", "上尾（含等號）", "本章雙尾規則"],
    "模擬近似": [q_hat, upper_hat, min(1., 2*min(q_hat, upper_hat))],
    "Poisson 解析值": [lower, upper, p_two],
}))

# %% [markdown]
# 很小的尾端機率，可能只由幾次模擬決定。
# 若一次也沒抽到，近似比例雖是零，真機率未必為零。
# 增加模擬可以降低抽樣誤差，但不會修正錯誤模型。
# 這裡有 Poisson 解析分布，因此能直接核對尾端。
# 後面空間檢驗較複雜，才更依賴模擬。
#
# ## 6.3 拒絕與未拒絕各能說到哪裡
#
# 若 p 值低於事前門檻，就拒絕這組假設。
# 含義是：這個計數在該模型下相當極端。
# 原因可能是率偏低，也可能是變異太小。
# 單一總數無法把這兩種原因完全分開。
# 地震序列也會讓多個事件彼此相關。
#
# 若沒有拒絕，只能說尚未找到足夠的不相容證據。
# 「未拒絕」不等於「證實」。
# 資料少、變異大時，很不同的模型也可能通過。
# **統計功效（statistical power）**是指定替代情況下，檢定能拒絕虛無假設的機率。
# 第 18 章會用模擬展示這個能力。
#
# 假設是在看結果前選定的。
# 若看到失敗才換分布，再把通過當成驗證，
# 就低估了選擇規則帶來的影響。
# 第 16 章會事先並列不同計數變異假設。
# 它們回答不同條件下的問題，不能互相替代。
#
# ## 6.4 分數回答：這份預報給結果多少支持
#
# **對數分數（logarithmic score）**是預報給實際結果的機率取自然對數。
# 本書採數值越大越好的方向。
# 若某種結果的機率很低，發生後就得到較低分。
# 它獎勵適當配置機率，沒有要求每次都猜中。
#
# 對 Poisson 總數而言，完整分數是
#
# $$S(\Lambda,N)=N\ln\Lambda-\Lambda-\ln N!.$$
#
# $N\ln\Lambda$ 讓事件出現時的較高率得到支持。
# $-\Lambda$ 則讓無限制提高總量付出代價。
# 即使只用總數，也不能只保留第一項。
# 否則把所有率提高，就會永遠得分更好。
#
# 下面比較兩個事先標為示意的期望數。
# A 是 SUP 的期望數；B 是它的兩倍。
# B 用來看懂分數，不是新訓練的預報模型。
# 觀測總數相同，兩條曲線的分數就能相減。

# %% tags=["remove-input"]
expected_a, expected_b = expected, 2*expected
possible = np.arange(max(observed+10, int(expected_b+5*np.sqrt(expected_b))))
fig = go.Figure()
for value, label, color in [(expected_a, "A：SUP", ACCENT), (expected_b, "B：兩倍期望（示意）", PALETTE[1])]:
    fig.add_trace(go.Scatter(x=possible, y=poisson.logpmf(possible, value), name=label, line_color=color))
fig.add_vline(x=observed, line_color=QUAKE_COLOR)
apply_layout(fig, title="相同結果，在兩份計數預報下的對數分數", xaxis_title="觀察到的總數",
    yaxis_title="對數分數（越大越好）", height=420)
fig

# %% [markdown]
# 圖上不同總數，可能由不同模型得較高分。
# 模型不會因為這次勝出，就被證明永遠較好。
# 若 B 是看過本次總數後才選的，
# 其分數差只能作回顧描述，不能當前瞻成效。
# 因此，正式比較要固定規則，再累積未見資料。
#
# **概似比（likelihood ratio）**將同一結果在兩個模型下的概似相除。
# 取對數後，就變成分數之差：
#
# $$\ln\frac{L_B}{L_A}=S_B-S_A
# =N\ln\frac{\Lambda_B}{\Lambda_A}-(\Lambda_B-\Lambda_A).$$
#
# 差值大於零，表示 B 對這份結果給較多支持。
# 差值等於零，表示這項計分無法區分兩者。
# 這是相對比較，並沒有檢查兩者是否都很差。
# 上半章的一致性檢定與這裡的比較，角色不同。
#
# **資訊增益（information gain）**以相對對數分數描述新預報增加的支持。
# 若換了基準，增益也會改變。
# 第 18 章會擴充到事件位置、規模與逐事件平均。
# 那時還要處理區間與事件相依，不能只看正負號。
#
# ## 6.5 把三個問題放回同一張桌上
#
# 概似估計問「哪些參數支持學習資料」。
# 一致性檢定問「測試資料是否出乎模型預期」。
# 模型比較問「相同結果在哪張預報下較受支持」。
# 三者可能使用同一條對數公式，用途卻不相同。
#
# 例如把測試期總數直接除以測試長度，
# 可以得到回顧性的最佳固定率。
# 這個數字可以作診斷，卻不能回填到已發出的預報。
# 資料的時間角色，比公式看起來是否熟悉更重要。
#
# ## 6.6 一個分數需要哪些陪伴資訊
#
# 若只看到「分數提高」，先問比較用了哪批事件。
# 兩個模型必須在同一組時間、空間與規模箱上計算。
# 少算一些沒有事件的區域，也可能提高分數。
# 但那等於改了問題，不能算成模型進步。
#
# 同樣地，先選出表現最好的一窗再報分數，
# 會隱藏其他窗的結果。
# 本書的正式比較使用事先固定的全部測試窗。
# 逐窗圖用來找差異出現在哪裡，
# 不能再據此刪掉不利於某模型的窗。
#
# 總數檢驗也沒有替地震逐顆認證。
# 如果一段序列貢獻很多目標事件，
# 總數雖然增加，獨立的資訊未必等比例增加。
# 第 16 章會先看逐窗累積，
# 第 18 章再討論比較區間的相依限制。
#
# 因此，一份可讀的評估至少交代模型、資料與規則。
# 數值保留到適當精度，尾端方向寫清楚。
# 模擬結果附上次數，解析結果註明分布。
# 讀者才有辦法重算，也才知道結論的範圍。
#
# ## 6.7 本章填入的規格欄位
#
# 評分欄先分成一致性與比較。
# 正式的空間、規模與數量檢驗名稱留待後面。
# 目前只用總數，避免太早混入所有模型。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("共同實驗規格", "沿用第 3–4 章的資料、窗、格與規模箱"),
    ("一致性", "統計量 → 模型分布 → 分位數與指定尾端；α=0.05"),
    ("比較", "相同觀測的對數分數差；基準先固定"),
    ("正式檢驗", "待填（第 16–18 章）"),
    ("完整度與規模比例", "待填（第 7–8 章）"),
])))

# %% [markdown]
# 三個統計工具章至此完成。
# 下一章 {doc}`07_gr_bvalue` 回到地震目錄：相同總數如何分到大小地震？
#
# ## 參考資料與延伸閱讀
#
# - Zechar（2010），[地震預報評估入門](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf)。CORSSA 免費教材，先讀檢驗目的，再看分數的定義。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。對照一致性檢驗與模型比較兩組結果，留意各自使用的分布假設。
