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
# # 16. 檢驗 I：數量對不對？N-test
#
# 25 顆目標地震，已經落在測試期。
# 第 15 章的 EEPAS 圖給了期望數。
# 本章先把地圖收起來，只比總數。
# 數量檢驗（number test，N-test）問：
# 模型給的計數分布，容得下觀測嗎？
#
# 一張預報寫著「期望 14 顆」。
# 這個數字仍缺少上下波動的範圍。
# 一致性（consistency）指觀測符合指定模型的程度。
# 本章會固定預報，改看兩種計數假設。

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
# ## 16.1 先把同一批事件數清楚
#
# $\omega_{wjk}$ 是一個窗、格、規模箱的觀測數。
# 下標 $w$ 代表預報窗。
# 把三個方向加總，就得到目標總數。
# $\Lambda_{wjk}$ 是同一位置的期望數。
# 兩個陣列必須使用相同邊界。
#
# $$N_{\rm obs}=\sum_{w,j,k}\omega_{wjk},\qquad
# \mu=\sum_{w,j,k}\Lambda_{wjk}.$$
#
# 25 顆是 HORUS 2024 年版的測試結果。
# 論文使用的目錄版本有 27 顆。
# 本章的檢驗一律使用本站目錄。
# 表中的 ETAS 也採本站近似。
# 它只含背景與已知事件的第一代貢獻。
#
# 40 個預報窗各有自己的資料截止。
# 累積曲線把窗內期望數一路相加。
# 這叫逐窗累積（cumulative count）。
# 曲線的末端檢查十年總量。
# 曲線中途的斜率則顯示增加速度。

# %% tags=["remove-input"]
observed_by_window = omega.sum(axis=(1, 2))
fig = go.Figure()
for m, a in forecasts.items():
    fig.add_scatter(x=np.arange(1, len(windows)+1), y=a.sum(axis=(1, 2)).cumsum(),
                    name=m, line=dict(color=colors[m]))
fig.add_scatter(x=np.arange(1, len(windows)+1), y=observed_by_window.cumsum(),
                name="目標地震", line=dict(color=QUAKE_COLOR, shape="hv"))
apply_layout(fig, title="逐窗累積：哪段時間拉開差距？", xaxis_title="預報窗序號", yaxis_title="累積顆數")
display(Markdown(f"圖：共 {len(windows)} 窗。紅線累積到 {int(omega.sum())} 顆。模型線為期望數。"))
fig

# %% [markdown]
# 紅線的跳升對應目標事件集中出現。
# 線條在序列期間分開，值得追查。
# 序列（sequence）是一群時空接近的地震。
# 這裡先按曆年分解，避免事後挑窗。
# 2016 年另列，只作描述性比較。
# 曆年分組本身尚未判定地震的觸發關係。
#
# 下表將每顆事件保留在原目錄內。
# 「其他年份」也完整列入總數。
# 若只刪掉密集那年再評分，題目就變了。
# 序列分解應與完整結果一起報告。
# 未來實驗要在看分數前固定這張圖。

# %% tags=["remove-input"]
year_counts = targets.groupby(targets["time"].dt.year).size()
summary = pd.DataFrame({"組別": ["2016 年", "其他年份", "全部"],
    "目標顆數": [int(year_counts.get(2016, 0)), int(len(targets)-year_counts.get(2016, 0)), len(targets)]})
display(HTML(summary.to_html(index=False)))

# %% [markdown]
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d16_csep_tests_flow", caption="先分清數量與分配，再選擇要固定和模擬的量。")

# %% [markdown]
# ## 16.2 兩個尾端，各問一件事
#
# $N\sim\operatorname{Poisson}(\mu)$ 固定均值為預報總數。
# Poisson 分布的變異數也等於 $\mu$。
# 離散尾端機率（discrete tail probability）會累加整數機率。
# 本章的兩個尾端都包含觀測值。
#
# $$\delta_1=P(N\ge N_{\rm obs}),\qquad
# \delta_2=P(N\le N_{\rm obs}).$$
#
# $\delta_1$ 小，表示觀測數偏多。
# 模型很少生成這麼多事件，顯示低估。
# $\delta_2$ 小，則表示觀測數偏少。
# 雙尾檢驗（two-sided test）同時檢查兩端。
# 本章用總顯著水準 $\alpha=0.05$。
# 每端的門檻是 $\alpha/2$。
#
# $\delta_1+\delta_2$ 可能大於一。
# 兩個尾端都算到 $N=N_{\rm obs}$。
# 這是整數計數的正常結果。
# 讀表時先看尾端方向，再看門檻。
# 不要把一個尾端直接當成雙尾 p 值。
#
# 下圖的線段是中央預測區間。
# 預測區間（prediction interval）涵蓋模型下的計數波動。
# 它的單位是顆數。
# 它與參數估計的信賴區間回答不同問題。
# 離散邊界的判定，仍以尾端機率為準。

# %% tags=["remove-input"]
alpha = 0.05
poisson_results = {m: csep.n_test(a, omega, alpha=alpha) for m, a in forecasts.items()}
fig = go.Figure()
for i, (m, r) in enumerate(poisson_results.items()):
    lo, hi = r["interval"]
    fig.add_scatter(x=[lo, hi], y=[m, m], mode="lines+markers", name=m,
                    line=dict(color=colors[m], width=5))
fig.add_vline(x=int(omega.sum()), line_color=QUAKE_COLOR, line_dash="dash")
apply_layout(fig, title="Poisson 計數的中央預測區間", xaxis_title="十年目標顆數", showlegend=False)
display(Markdown(f"圖：線段涵蓋中央 {100*(1-alpha):.0f}% 的分位點範圍。紅虛線為觀測 {int(omega.sum())} 顆。"))
fig

# %% tags=["remove-input"]
display(HTML(pd.DataFrame([
    {"模型": m, "期望顆數": r["n_forecast"], "δ1": r["delta1"], "δ2": r["delta2"], "判定": "未拒絕" if r["passed"] else "拒絕"}
    for m, r in poisson_results.items()]).to_html(index=False, formatters={"期望顆數": lambda x: f"{x:.2f}", "δ1": lambda x: f"{x:.3g}", "δ2": lambda x: f"{x:.3g}"})))

# %% [markdown]
# 表中的四個期望數都低於觀測總數。
# Poisson 檢驗也都指出數量不一致。
# 此處的證據針對「預報加上計數假設」。
# 它尚未指出錯在背景率或觸發項。
# 下一步要分開追查均值與波動。
#
# 一條滾動更新的率曲線用到了新歷史。
# 本章把已產生的率陣列視為固定輸入。
# 再用 Poisson 計數作參考分布。
# 這是網格率的教學檢查。
# 它與完整 ETAS 序列模擬有不同假設。
# 完整模擬還會讓未來事件互相觸發。
#
# ## 16.3 同一個均值，容許更大的波動
#
# 67.76 是論文採用的十年計數變異數。
# 作者從 CPTI15 的歷史地震估計。
# 期間是 1882–2011，採不重疊十年窗。
# 該數值含主震與餘震。
# 本章將它當成外部指定的示範值。
# 它尚未針對本站 HORUS 版本重估。
#
# 變異數大於均值，稱為過度離散。
# 英文是 overdispersion。
# 負二項分布（negative binomial distribution）
# 可用均值與較大的變異數指定計數波動。
# 本章保留各模型的均值，只換波動假設。

# %% [markdown]
# $$v>\mu,\qquad r=\frac{\mu^2}{v-\mu},\qquad
# \theta=\frac{\mu}{v},\qquad N\sim\operatorname{NB}(r,\theta).$$
#
# $r$ 是形狀參數，$\theta$ 是此處的機率參數。
# 這套參數化給出均值 $\mu$、變異數 $v$。
# $r$ 可以是正實數。
# 負二項版沿用相同的兩個尾端。
# 因此兩圖的差異來自分布寬度。

# %% tags=["remove-input"]
variance = 67.76  # Biondini et al. (2023), p. 1688, CPTI15 1882–2011.
nb_results = {m: csep.n_test_nb(a, omega, variance=variance, alpha=alpha) for m, a in forecasts.items()}
fig = go.Figure()
for m, r in nb_results.items():
    fig.add_scatter(x=r["interval"], y=[m, m], mode="lines+markers", name=m,
                    line=dict(color=colors[m], width=5))
fig.add_vline(x=int(omega.sum()), line_color=QUAKE_COLOR, line_dash="dash")
apply_layout(fig, title="負二項計數的中央預測區間", xaxis_title="十年目標顆數", showlegend=False)
display(Markdown(f"圖：外部指定變異數 {variance:.2f}。中央 {100*(1-alpha):.0f}% 區間下，{sum(r['passed'] for r in nb_results.values())} 個模型未被拒絕。"))
fig

# %% [markdown]
# 較寬的區間涵蓋了觀測數。
# 「通過」在這裡表示未被檢驗拒絕。
# 各模型的期望總數完全沒有增加。
# 因此，通過本身無法證明率已抓準。
# 更寬的分布也會降低分辨錯誤的能力。
#
# 同一批結果先跑 Poisson，再看負二項。
# 這個順序用來看清假設的影響。
# 若因失敗才換分布，必須記下這項選擇。
# 正式實驗應事先決定變異數的來源。
# 測試期只負責提供觀測，不負責調寬區間。
#
# ## 16.4 多看幾次，也要多算一次門檻
#
# 四個模型各做一次檢驗，共有四次機會。
# 多重檢定（multiple testing）指同時作多個判斷。
# 即使各模型合適，也可能偶然出現拒絕。
# Bonferroni 校正會把總水準分給各檢驗。
# 英文名稱是 Bonferroni correction。

# %% [markdown]
# $$\alpha_{\rm each}=\frac{\alpha}{K},\qquad
# \alpha_{\rm tail}=\frac{\alpha}{2K}.$$
#
# $K$ 是事先列出的檢驗數。
# 這個界限不要求各次檢驗獨立。
# 若兩種計數假設都列為正式判斷，
# 檢驗家族也必須把兩者算入。
# 本章將負二項視為假設敏感度展示。
# 圖中仍顯示未校正區間，方便比較分布。
#
# 數量表應同時留下尾端機率與設定。
# 只留「通過」兩字，讀者無法重查。
# 本章填入數量的一致性規格。
# 位置與規模會在下一章分開檢查。
# 兩張總數相同的地圖，仍可能指向不同地區。
#

# %% [markdown]
# ## 16.5 從一列結果寫出可核對的判斷
#
# 表中的 SUP 列，可分三次讀完。
# 先讀期望顆數，知道模型中心在哪裡。
# 再讀上尾機率，衡量觀測偏多的程度。
# 最後才依既定門檻寫下判定。
# 這個順序讓數字與判斷各有位置。
#
# 一個很小的上尾機率，有明確的條件。
# 它假定 SUP 的均值就是計數分布中心。
# 它也假定十年總數服從 Poisson 分布。
# 若重做許多次這樣的實驗，
# 超過本站觀測數的結果很少出現。
# 這句話說的是模型下的重複實驗。
# 它沒有給出模型為真的機率。
#
# 負二項表保留完全相同的模型中心。
# 新的計數分布在兩側容許更大波動。
# 同一個觀測值因而變得較常見。
# 這項比較指出波動假設的重要性。
# 它並未找出低估總數的物理原因。
# 背景率與觸發貢獻仍須各自檢查。
#
# 累積曲線還保留另一種線索。
# 如果差距集中在少數跳升處，
# 後續研究可以檢查那些窗的更新方式。
# 如果差距每年都穩定擴大，
# 後續研究可以先查整體率的尺度。
# 兩者都是提出下一個問題的方法。
# 圖形本身尚不足以確定成因。
#
# 一個被拒絕的模型仍可能提供空間資訊。
# 例如模型把多數期望數放在正確地區，
# 卻把每格的率都乘上太小的常數。
# 數量檢驗會看到整體不足。
# 下一章的條件檢驗會另看分配比例。
# 把這兩項拆開，才能說清楚改哪裡。
#
# 一個未被拒絕的模型也可能位置很差。
# 只要它給出夠合適的總數分布，
# 數量檢驗就可能保留它。
# 因此，本章的合格範圍只有總量。
# 規格卡會留下位置與規模的空欄。
# 讀者可以沿空欄追到下一章的證據。
#
# 四條曲線、兩種區間各有用途。
# 曲線協助定位差距發生的時間。
# 區間協助衡量差距是否罕見。
# 尾端機率保留判斷的細節。
# 正式報告應讓三者使用同一份觀測。
# 如此才能避免圖、表與結論各說各話。
#

# %% [markdown]
# ## 16.6 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("共同目標", "HORUS；R 內、深度 ≤ 40 km、取整 Mw ≥ 5.0"),
    ("檢驗期間", "2012–2021；40 個滾動窗；保留全部目標事件"),
    ("數量一致性", "N-test；δ1、δ2；α=0.05；雙尾"),
    ("計數假設", "Poisson 主展示；負二項 v=67.76 為外部敏感度示範"),
    ("多重判斷", "先列 K，再用 Bonferroni 分配水準"),
    ("位置／規模一致性", "待填（第 17 章）"),
    ("比較與基準", "待填（第 18 章）"),
])))

# %% [markdown]
# 同一個總數，可以分配到不同格箱。
# 下一章：{doc}`17. 位置與規模的檢驗 <17_test_space_magnitude>`。
#
# ## 參考資料與延伸閱讀
#
# - Zechar（2010），[地震預報評估入門](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf)。
#   CORSSA 免費教材。先讀一致性，再讀模型比較。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。
#   [免費機構全文](https://www.earth-prints.org/handle/2122/17084)。先對照研究區與測試期，再讀檢驗。
# - Bayona 等（2022），[混合預報的前瞻檢驗](https://doi.org/10.1093/gji/ggac018)。
#   免費全文。對照負二項與多重判斷的用途。
