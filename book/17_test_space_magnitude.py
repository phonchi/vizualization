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
# # 17. 檢驗 II：位置與規模對不對？S、M、cL、L
#
# 25 顆目標地震，可以全擠在一個地區。
# 它們也可以分散到許多地區。
# 第 16 章只數總量，看不出這個差別。
# 本章保留同一批事件，再看落點與規模。
#
# 一個格箱同時寫著期望數與觀測數。
# 本章先用這兩欄算分數。
# 再讓模型產生模擬目錄，作為比較對象。
# 每項檢驗都要回答兩件事：
# 哪些量固定？哪些量重新抽樣？

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy, italy_models, csep_teaching as csep
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, apply_layout, show_diagram

models = ("SUP", "PPE", "EEPAS", "ETAS")
colors = dict(zip(models, [ACCENT, "#008300", "#4a3aa7", "#eda100"]))
cat = italy.experiment_catalog()
targets = italy.target_events(cat, "testing")
windows = italy.forecast_windows()
omega = italy.bin_targets(targets, windows)
forecasts = {m: italy_models.get_forecast(m, "testing") for m in models}
assert all(a.shape == omega.shape for a in forecasts.values())
assert int(omega.sum()) == len(targets)

# %% [markdown]
# ## 17.1 從一個格箱讀懂三種分數
#
# $\Lambda=0.2$、$\omega=1$ 是一個合成例子。
# 第一欄給期望數，第二欄給實際顆數。
# 逐箱 Poisson 對數概似稱為 POLL。
# 英文是 Poisson log-likelihood。
# 它衡量這個計數在指定率下的相容程度。
#
# $$\operatorname{POLL}(\Lambda,\omega)
# =-\Lambda+\omega\ln\Lambda-\ln(\omega!).$$
#
# $-\Lambda$ 會對配置過多期望數扣分。
# 有事件的箱還加上 $\omega\ln\Lambda$。
# 事件落在極低率的箱，會帶來較大扣分。
# $\ln(\omega!)$ 保留同箱多顆的計數因素。
# 空箱也有分數，因此計算必須涵蓋全部箱。
#
# 一整張預報的分數，是所有 POLL 相加。
# 聯合 Poisson 對數概似稱為 jPOLL。
# 英文是 joint Poisson log-likelihood。
# 本章用 $b$ 簡寫一個窗、格、規模箱。
#
# $$\operatorname{jPOLL}=\sum_b\operatorname{POLL}(\Lambda_b,\omega_b).$$
#
# $\omega=1$ 與 $\omega=2$ 也能另作處理。
# 若只記「有事件」，兩者都寫成一。
# 二元對數概似稱為 BILL。
# 英文是 binary log-likelihood。
# 它用箱內是否有事件來評分。
#
# $$z_b=\mathbf{1}(\omega_b>0),\qquad
# \operatorname{BILL}_b=
# \begin{cases}
# \ln(1-e^{-\Lambda_b}),&z_b=1,\\
# -\Lambda_b,&z_b=0.
# \end{cases}$$
#
# $1-e^{-\Lambda_b}$ 是箱內至少一顆的機率。
# BILL 將同箱的多顆事件合併成一次佔用。
# 這會改變叢集對分數的影響。
# 它也會捨去同箱顆數的部分資訊。
# 論文的二元檢驗與本章的計數檢驗須分開讀。

# %% tags=["remove-input"]
toy_rate = np.repeat(0.2, 3)
toy_omega = np.arange(3)
display(HTML(pd.DataFrame({"期望數": toy_rate, "觀測顆數": toy_omega,
    "POLL": csep.poll(toy_rate, toy_omega), "BILL": csep.bill(toy_rate, toy_omega)}).to_html(index=False, float_format=lambda x: f"{x:.3f}")))

# %% [markdown]
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d16_csep_tests_flow", caption="N、S、M、cL、L 各自回答不同問題。")

# %% [markdown]
# ## 17.2 固定總數，再抽事件落在哪裡
#
# 25 顆事件可以依預報比例重新分配。
# 條件化（conditioning）是固定已指定的量。
# 這裡固定總數，讓各箱的計數改變。
# 多項分布（multinomial distribution）描述這種分配。
# 各箱機率相加必須等於一。
#
# $$p_b=\frac{\Lambda_b}{\sum_a\Lambda_a},\qquad
# (\omega_1,\ldots,\omega_B)\mid N
# \sim\operatorname{Multinomial}(N;p_1,\ldots,p_B).$$
#
# $N p_b$ 是固定總數後的格箱期望數。
# 工具箱先將整張預報縮放到這個總數。
# 然後依 $p_b$ 抽出每顆事件所在的箱。
# 各箱共享固定總數，計數會彼此牽動。
# 某箱多分一顆，其他箱可分的顆數就減少。
#
# 一組模擬會產生一個 jPOLL 分數。
# 重複抽樣便得到參考分布。
# 分位數分數（quantile score）記作 $q$。
# 它是模擬分數小於或等於觀測值的比例。
# 本章採低尾判斷，因為低分代表較差的相容性。
#
# $N$ 固定後，總量的錯誤會被移除。
# 這項設定刻意隔開前章的數量問題。
# 因此，同時報告 N-test 很重要。
# 讀者才能看見「總量不足、分配尚可」的組合。
#
# ## 17.3 先把時間加總，再看位置
#
# 一個空間格可以包含多個規模箱。
# 邊際化（marginalization）是把不看的維度加總。
# 空間檢驗（spatial test，S-test）只看位置。
# 本章先合計所有時間窗，再合計規模箱。
# 每次模擬固定目標總數，再分配到空間格。

# %% tags=["remove-input"]
n_sim = 1000
n_bins = omega.size
estimated_dense_mb = n_sim * n_bins * 8 / 1e6
# Full dense simulations would exceed 1 GB; toolkit samples occupied bins only.
display(Markdown(f"模擬設定：每項 {n_sim:,} 次。每份模擬目錄使用相同的預報規格及評分方式。"))
space_mag_omega = omega.sum(axis=0)
space_mag_rates = {m: a.sum(axis=0) for m, a in forecasts.items()}
s_results = {m: csep.s_test(a, space_mag_omega, n_sim=n_sim, seed=1701) for m, a in space_mag_rates.items()}
from plotly.subplots import make_subplots

def distribution_figure(results, title):
    result_fig = make_subplots(rows=2, cols=2, subplot_titles=[f"{m} · q={results[m]['quantile']:.3f}" for m in models],
                               horizontal_spacing=0.15, vertical_spacing=0.23)
    for i, m in enumerate(models):
        r, k = divmod(i, 2)
        result_fig.add_trace(go.Histogram(x=results[m]["simulated"], nbinsx=35,
            marker_color=colors[m], name=m, showlegend=False), row=r+1, col=k+1)
        result_fig.add_vline(x=results[m]["statistic"], line_color=QUAKE_COLOR,
                            line_dash="dash", row=r+1, col=k+1)
    result_fig.update_xaxes(title_text="對數概似分數")
    result_fig.update_yaxes(title_text="模擬次數")
    return apply_layout(result_fig, title=title, height=660, bargap=0.08)

fig = distribution_figure(s_results, "S-test：只比較空間分配")
display(Markdown(f"圖：每模型模擬 {n_sim:,} 次。紅虛線為 {int(omega.sum())} 顆目標地震的分數。標題 q 使用含等號的低尾比例。"))
fig

# %% [markdown]
# SUP 的空間分數幾乎擠成一條線。
# 退化（degeneracy）指分布集中到少數取值。
# 等面積 SUP 給每格相同期望數。
# 若事件各占一格，換位置不會改變分數。
# 有多顆落在同格時，階乘項才會改分數。
#
# SUP 的低分因此主要反映同格多顆。
# 它無法指出哪一條構造帶更值得注意。
# 若把評分改成只保留逐事件對數率，
# 均勻模型在固定總數下甚至完全同分。
# 讀者要先查統計量，才解讀分位數。
#
# $q=0$ 表示有限次模擬中沒有更低分。
# 它仍受到模擬次數與同分處理影響。
# 浮點運算也可能改變幾乎相等的排序。
# 因此，接近退化的分布要連圖一起讀。
# 增加小數位數不會增加位置辨識力。
#
# ## 17.4 把位置加總，再看規模
#
# 同一個規模箱可以收進各地的事件。
# 規模檢驗（magnitude test，M-test）只看規模。
# 本章把時間與空間加總，保留規模箱。
# 模擬仍固定總數，只重新分配規模。
# 這樣能分開查位置分配與大小比例。

# %% tags=["remove-input"]
m_results = {m: csep.m_test(a, omega, n_sim=n_sim, seed=1702) for m, a in forecasts.items()}
fig = distribution_figure(m_results, "M-test：只比較規模分配")
display(Markdown(f"圖：每模型模擬 {n_sim:,} 次，保留 {omega.shape[-1]} 個規模箱。紅虛線為目標事件的分數。"))
fig

# %% [markdown]
# 幾張相似的規模分布，可能來自共同設定。
# SUP、PPE 與本站 ETAS 使用同一個 GR 斜率。
# 它們的 M-test 因而可能很接近。
# EEPAS 的規模核則帶入來源事件的尺度。
# 檢驗應核對實際輸出，避免從名字猜結果。
#
# 一張 M-test 圖只保留規模邊際分布。
# 它會把「大震放錯地區」的配對資訊加總掉。
# 空間圖也會把規模差異加總掉。
# 兩項各自合適，仍需檢查聯合分配。
# 下一項檢驗便保留原來的格箱。

# %% [markdown]
# ## 17.5 cL 與 L，相差一個總數條件
#
# 原始陣列保留窗、格與規模箱。
# 條件概似檢驗稱為 cL-test。
# 英文是 conditional likelihood test。
# 它固定總數，檢查完整格箱分配。
# 本章保留時間，所以也檢查事件集中在哪些窗。

# %% tags=["remove-input"]
cl_results = {m: csep.cl_test(a, omega, n_sim=n_sim, seed=1703) for m, a in forecasts.items()}
fig = distribution_figure(cl_results, "cL-test：固定總數，保留時間、位置與規模")
display(Markdown(f"圖：每模型模擬 {n_sim:,} 次，使用 {omega.size:,} 個時空規模箱。紅虛線為目標事件的分數。"))
fig

# %% [markdown]
# 另一種模擬會先抽總數，再抽格箱。
# 概似檢驗（likelihood test，L-test）採這個流程。
# 它使用原始率，不先縮放到觀測總數。
# 因此，總數與格箱分配都會影響它的分數。
# N-test 失敗時，L-test 也可能受到拖累。
#
# 下表把四種檢驗放在一起。
# S、M 與 cL 的參考目錄總數相同。
# L 的參考目錄則各有不同總數。
# 同一欄的低分位數具有相同方向。
# 不同欄保留的資料內容則各不相同。

# %% tags=["remove-input"]
l_results = {m: csep.l_test(a, omega, n_sim=n_sim, seed=1704) for m, a in forecasts.items()}
display(HTML(pd.DataFrame([
    {"模型": m, "S 的 q": s_results[m]["quantile"], "M 的 q": m_results[m]["quantile"],
     "cL 的 q": cl_results[m]["quantile"], "L 的 q": l_results[m]["quantile"]} for m in models
]).to_html(index=False, float_format=lambda x: f"{x:.3f}")))

# %% [markdown]
# 一列結果若只有 L 很低，先回查數量。
# 若 S 很低，先回查空間分配與同格叢集。
# 若 M 很低，先回查大小地震的比例。
# 若 cL 很低，再查時間與聯合分配。
# 這些順序用來定位問題，無法單獨證明成因。
#
# 一組小分位數也會受到多重判斷影響。
# 正式檢驗應事先列出主檢驗與總水準。
# 本章展示分數分布，不另選最有利的門檻。
# 原始分位數可以保留完整訊息。
# 第 16 章的 Bonferroni 方法可處理檢驗家族。
#
# 論文中的 BILL 採箱內有無事件評分。
# 本章的 S、M、cL、L 使用 POLL 系列。
# 因此，本站表格無法直接替代論文表格。
# 比較前還要核對目錄版本與時間窗。
# 只有統計量名稱相近，並不足以視為同一實驗。
#
# 25 顆觀測只提供有限資訊。
# 未被拒絕的模型，仍可能靠寬鬆分布留下。
# 本章回答的是各模型是否與資料相容。
# 下一章改用同一顆地震比較兩個模型。
# 那個問題需要明確寫出比較基準。
#
# 一份模擬目錄可以想成重新裝球。
# S-test 的盒子只標地區名稱。
# M-test 的盒子只標規模範圍。
# cL-test 的盒子還標發報窗。
# 每次裝球的機率由對應預報決定。
# 這個具體操作就是各圖的共同骨架。
#
# 一個盒子裝了多顆，仍須保留顆數。
# 若改成只記盒子有沒有球，
# 評分對象就轉成前面的 BILL。
# 模擬與觀測必須採用同一套記錄方式。
# 只改觀測端的分數會破壞比較。
# 換統計量時，也要一併重做參考分布。
#
# 一個空間熱區的判斷需要兩層證據。
# 先看觀測是否落在模型配置較高率的地區。
# 再看模型自己常產生哪些位置組合。
# S-test 把第二層放進模擬分布。
# 模型若配置非常尖銳的熱區，
# 落在熱區外的事件便可能受到較大扣分。
#
# 一張圖的直條高低是模擬頻數。
# 高條表示那段分數較常被模型產生。
# 紅線則只代表真實目標目錄的一個分數。
# 讀者應在每個小圖內比較兩者。
# 跨模型的原始分數尺度可能不同。
# 分位數才把位置轉成可讀的尾端比例。
#
# ## 17.6 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("共同目標", "HORUS；R 內、深度 ≤ 40 km、取整 Mw ≥ 5.0"),
    ("期間與輸出", "2012–2021；40 窗 × 177 格 × 25 規模箱"),
    ("數量", "N-test；Poisson 與外部負二項假設分開報告"),
    ("位置", "S-test；先合計時間，再合計規模；固定總數"),
    ("規模", "M-test；合計時間與空間；固定總數"),
    ("聯合一致性", "cL 固定總數；L 另抽總數；均使用 jPOLL"),
    ("模擬與分數", f"每項 {n_sim} 次；固定種子；q 含同分；低尾 α=0.05"),
    ("比較與基準", "待填（第 18 章）"),
])))

# %% [markdown]
# 一致性欄已填完，比較欄仍保留空位。
# 下一章：{doc}`18. 哪張預報比較有資訊 <18_test_comparison>`。
#
# ## 參考資料與延伸閱讀
#
# - Zechar（2010），[地震預報評估入門](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf)。
#   CORSSA 免費教材。先讀一致性，再讀模型比較。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。
#   [免費機構全文](https://www.earth-prints.org/handle/2122/17084)。先對照研究區與測試期，再讀檢驗。
# - Bayona 等（2022），[前瞻預報與二元評分](https://doi.org/10.1093/gji/ggac018)。
#   免費全文。對照計數與箱內有無事件的差異。
