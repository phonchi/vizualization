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
# # 7. 大小地震的比例：GR 律與 b 值
#
# 規模 5 的地震，通常比規模 6 多。
# 第 6 章用總數比較預報。
# 本章把同一個總數拆到規模軸上。
# 大小地震的比例，決定每箱分到多少。
#
# HORUS 的 `mb` 欄保留一位小數。
# 本章選學習期、測試區內的事件。
# 規模下限先設為名目 2.5。
# 完整度的檢查留給下一章。
# 讀完後，你能解釋斜率與估計區間。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

# %% [markdown]
# ## 7.1 先數每個規模以上有幾顆
#
# 下圖的橫軸是規模門檻。
# 累積次數（cumulative count）是超過門檻的筆數。
# 規模 5 的點，也包含規模 6 的事件。
# 因此，相鄰點會共用許多資料。
#
# 下圖的縱軸採對數刻度。
# 規模—頻率分布（frequency–magnitude distribution）描述各規模的數量。
# 圖上的下降幅度，代表大小事件的比例。
# 最右端只有少數事件，跳動會較大。

# %% tags=["remove-input"]
cat = italy.experiment_catalog()
learn = cat.loc[cat.in_R & (cat.period == "learning")]
m = learn.loc[learn.mb >= 2.5, "mb"].to_numpy()
dm, lower = 0.1, 2.5
b_aki = np.log10(np.e) / (m.mean() - (lower - dm / 2))
b_exact = np.log1p(dm / (m.mean() - lower)) / (dm * np.log(10))
thresholds = np.round(np.arange(lower, 7.1, dm), 1)
counts = np.array([(m >= v - 1e-8).sum() for v in thresholds])
fig = go.Figure(go.Scatter(x=thresholds[counts > 0], y=counts[counts > 0],
    mode="markers", name="HORUS 累積次數", marker_color=ACCENT))
fig.add_trace(go.Scatter(x=thresholds,
    y=len(m) * 10 ** (-b_exact * (thresholds - lower)),
    mode="lines", name="離散 GR 擬合", line=dict(color="#4a3aa7")))
apply_layout(fig, title="學習期的規模—頻率分布",
    xaxis_title="取整規模門檻 Mw", yaxis_title="累積事件數", yaxis_type="log")
display(Markdown(f"圖 1｜HORUS 真實資料共 {len(m):,} 筆。"
    f"期間為 {italy.SPEC.learning[0]}–{italy.SPEC.learning[1]}。"
    f"線的離散估計 b＝{b_exact:.3f}。"))
fig

# %% [markdown]
# ## 7.2 一級規模，數量差多少？
#
# $$N(M\ge m)=10^{a-bm}.$$
#
# 這條式子是 Gutenberg–Richter 律，簡稱 GR 律。
# GR 律用直線近似規模與對數次數。
# $b$ 值（b-value）是直線下降的斜率大小。
# $a$ 值（a-value）是這條直線的截距。
#
# $b=1$ 時，門檻升一級，數量剩十分之一。
# $b$ 較大時，大地震所占比例較小。
# 這個比較固定了較低門檻的事件數。
# 若總活動量同時改變，大震數仍可能增加。
#
# $a$ 會隨觀測期間長短改變。
# 同一地區觀測兩倍時間，通常記到更多事件。
# 把次數除以年數，直線便上下平移。
# $b$ 描述比例，$a$ 則包含觀測範圍的資訊。
# 兩者各自回答不同的數量問題。
#
# $$\frac{N(M\ge m+1)}{N(M\ge m)}=10^{-b}.$$
#
# 這個比值只留下 $b$。
# 讀者可以先算比值，再看直線是否合理。
# GR 律適用的是一段規模範圍。
# 右端的少量大震，較難單獨確定斜率。
# 左端的漏記事件，也會改變線形。
#
# ## 7.3 從比例走到 Aki 估計式
#
# $$P(M\ge m\mid M\ge m_0)=e^{-\beta(m-m_0)},\qquad
# \beta=b\ln 10.$$
#
# 這裡先把規模視為連續量。
# 超額規模（excess magnitude）是 $M-m_0$。
# 指數分布（exponential distribution）用單一率描述超額量。
# $\beta$ 是規模方向的衰減係數。
# 它的倒數就是平均超額規模。
#
# $$s(m)=\beta e^{-\beta(m-m_0)},\qquad m\ge m_0.$$
#
# $s(m)$ 是規模機率密度（probability density）。
# 密度在規模區間的面積，就是該區間的機率。
# 第 9 章會用這個面積分配預報數。
# 此處的指數分布描述規模。
# 第 4 章的等待時間則是另一個隨機量。
#
# $$\ell(\beta)=N\ln\beta-\beta\sum_i(m_i-m_0).$$
#
# 這個對數概似沿用第 5 章的做法。
# 各筆規模在模型中視為獨立抽樣。
# 對 $\beta$ 微分，再令結果為零。
# 最佳值使模型均值等於樣本均值。
#
# $$\widehat\beta=\frac{1}{\overline m-m_0},\qquad
# \widehat b=\frac{\log_{10}e}{\overline m-m_0}.$$
#
# 這就是 Aki 的最大概似估計式。
# 平均規模越接近下限，估出的 $b$ 越大。
# 其意義是：事件更集中在小規模端。
# 這個結果來自每一筆規模。
# 累積圖主要負責呈現結果。

# %% [markdown]
# ## 7.4 一位小數也要進入公式
#
# `mb=2.5` 代表一個寬度為 0.1 的箱。
# 分箱修正（binning correction）把記錄精度納入估計。
# 下限箱的左邊界是 2.45。
# 連續近似便使用這個邊界。
#
# $$\widehat b_{\rm Aki}
# =\frac{\log_{10}e}{\overline m-(2.5-0.1/2)}.$$
#
# 式中的半箱修正保留了箱子的寬度。
# 若直接減去 2.5，分母會偏小。
# 同一批資料便會得到較大的 $b$。
# 計算前須確認使用的是箱心還是箱界。
#
# $$\widehat b_{\rm disc}=
# \frac{1}{\Delta m\ln 10}
# \ln\left(1+\frac{\Delta m}{\overline m-m_{\min}}\right).$$
#
# 這是離散指數模型的精確概似解。
# 此處 $m_{\min}=2.5$，$\Delta m=0.1$。
# 離散模型（discrete model）把機率放在規模格點。
# 每升一箱，機率乘上固定的比例。
# 「精確」指這個模型的求解方式。
# 有限樣本仍有抽樣誤差。

# %% tags=["remove-input"]
display(Markdown(f"同一批資料的平均規模為 {m.mean():.4f}。"
    f"半箱修正得到 b＝{b_aki:.4f}。"
    f"離散公式得到 b＝{b_exact:.4f}。"
    f"未修正則為 {np.log10(np.e)/(m.mean()-lower):.4f}。"))

# %% [markdown]
# ## 7.5 區間回答哪一種不確定性？
#
# 下圖把整段學習期分成數段。
# 每一段都使用相同的規模下限。
# 自助抽樣（bootstrap）是從現有樣本反覆抽回資料。
# 每次重新估計，便得到一批 $b$ 值。
# 本章沿用第 5 章的百分位數區間。
#
# 每個區間取抽樣結果的中間 95%。
# 這個信賴區間（confidence interval）描述估計的抽樣波動。
# 此處固定了區域、下限與規模尺度。
# 目錄漏記或尺度換算的誤差，另需證據。
# 事件逐筆重抽，也採用了獨立抽樣近似。

# %% tags=["remove-input"]
rng = np.random.default_rng(707)
n_sim = 600
# 上界估算：若一次保留所有重抽樣規模，需要以下記憶體。
bytes_full = n_sim * len(m) * 8
# 實際逐次重抽，只保留估計值。
def boot_b(values):
    estimates = np.empty(n_sim)
    for i in range(n_sim):
        sample = rng.choice(values, size=len(values), replace=True)
        estimates[i] = np.log1p(dm / (sample.mean() - lower)) / (dm * np.log(10))
    return np.quantile(estimates, [0.025, 0.975])
whole_ci = boot_b(m)
rows = []
for start, end in [(1990, 1995), (1995, 2000), (2000, 2005), (2005, 2012)]:
    values = learn.loc[(learn.year >= start) & (learn.year < end) & (learn.mb >= lower), "mb"].to_numpy()
    bh = np.log1p(dm / (values.mean() - lower)) / (dm * np.log(10))
    lo, hi = boot_b(values)
    rows.append((f"{start}–{end-1}", bh, lo, hi, len(values)))
r = pd.DataFrame(rows, columns=["period", "b", "lo", "hi", "n"])
fig = go.Figure(go.Scatter(x=r.period, y=r.b, mode="markers",
    error_y=dict(type="data", array=r.hi-r.b, arrayminus=r.b-r.lo),
    marker=dict(color=ACCENT, size=9), customdata=r.n,
    hovertemplate="%{x}<br>b=%{y:.3f}<br>N=%{customdata}<extra></extra>"))
fig.add_hline(y=italy.SPEC.b_value, line_dash="dot", line_color="#555555")
apply_layout(fig, title="固定規模下限後，各時段的 b 值",
    xaxis_title="學習期內的子期間", yaxis_title="b 值與 95% bootstrap 區間")
display(Markdown(f"圖 2｜每段重抽 {n_sim} 次。"
    f"整段估計為 {b_exact:.3f}，區間 [{whole_ci[0]:.3f}, {whole_ci[1]:.3f}]。"
    f"虛線是論文固定值 {italy.SPEC.b_value:.3f}。"
    "每次重新抽樣，都重新估計一次 b 值。"))
fig

# %% [markdown]
# 圖中的點高低不同，先回看各段資料量。
# 抽樣波動、目錄品質與真實活動都可能參與。
# 估計偏差（estimation bias）是估計平均值偏離真值。
# 漏掉小震會抬高觀測平均規模。
# 代入公式後，便可能把 $b$ 壓低。
#
# 圖中的短線只處理固定樣本規則下的波動。
# 若同一序列中的規模有相依性，區間可能過窄。
# 可改用整段序列或時間區塊重抽。
# 這需要先定義區塊與比較目的。
# 本章保留簡單版本，明列其抽樣假設。
#
# 1.084 是義大利實驗採用的固定值。
# Biondini 等人用學習期目標數校準它。
# 本頁則用小規模以上的樣本估計斜率。
# 兩個數字的資料範圍與估計目的不同。
# 讀到差距時，應先對照方法，再談解釋。
#
# 四張快取預報繼續使用固定值 1.084。
# 本章的自估值用來理解估計流程。
# 若改用自估值，規模箱的比例就會改變。
# 那會形成另一份需要重新檢驗的預報。
# 下一章先處理更前面的資料下限問題。
#
# ### 地震剛發生後：b-positive 在改變什麼？
#
# **b-positive** 使用時間相鄰事件的正規模差，估計GR斜率。先算 $d_i=m_i-m_{i-1}$，只保留超過預定差值門檻的正差。它的動機是：大震後小震容易漏記，比剛記錄的事件更大的下一顆，通常比較容易被看見。
#
# 例如連續規模為3.0、2.8、3.4，相鄰差是−0.2與0.6；正差法保留0.6，絕對差法則可能兩個都用。兩種選樣不同，不能把所有「規模差估計」都叫b-positive。這是方法示意，不是把本章HORUS估計偷偷換成另一個數字。
#
# 在理想獨立指數規模模型下，正差超過門檻後的超額仍為指數分布。實際漏測、離散分箱及相鄰配對會改變條件與不確定性。Tinti與Gasperini（2024）比較多種差值方法，沒有建立正差法在所有漏測情境都較好的結論。附錄B.5列出假設與差別；第8章回到短期完整度。
#
# ## 7.6 把一個估計值讀完整
#
# 一列報告可以同時列出規模下限與樣本數。
# 讀者先看資料選了哪些事件。
# 同樣的斜率若來自不同區域，回答的問題就不同。
# 本頁把地點固定在測試區內。
# 收集區外圍的事件留給預報模型使用。
# 這個選擇讓估計對準要描述的區域。
#
# 一條累積曲線上的點彼此重疊。
# 規模較大的事件會出現在許多門檻中。
# 若直接把這些點當成獨立觀測做迴歸，
# 估計區間就可能低估資料的共同波動。
# 逐筆規模概似避開了這種重複計數。
# 圖形仍能幫助讀者找出彎曲的位置。
#
# 規模下限升高後，剩下的事件會變少。
# 若目錄原先漏掉小震，斜率可能先改變。
# 若樣本原先已經完整，估計仍會隨機跳動。
# 因此，穩定的數字需要足夠事件支撐。
# 下一章會把這個觀察變成完整度判準。
#
# 同一個平均值也可能來自不同分布。
# 例如，小規模與大規模混合後，均值恰好相同。
# Aki 公式便會給出相同估計。
# 累積圖卻可能顯示不同的彎曲形狀。
# 所以，報告估計值時應保留資料圖。
# 均值與圖形提供互補的資訊。
#
# 規模量測的誤差也會沿公式傳遞。
# 假設所有規模都被往上移動，
# 固定下限以上的平均超額便會增加。
# 實際選入的樣本也可能跟著改變。
# 讀者應先確認尺度一致，再比較斜率。
# HORUS 的均一化規模正是為了這項需求。
#
# 一個較低的斜率，也可能出現在單一叢集。
# 叢集（cluster）是時間或位置相近的一群事件。
# 某段序列若包含較多大震，便會改變整段比例。
# 它與地殼狀態的關係，需要另外的分析設計。
# 本章的分段圖只整理已發生的目錄。
# 它提供描述，並未給出下一次大震的日期。
#
# 預報中的規模比例還要乘上總期望數。
# 即使兩個模型使用相同的斜率，
# 它們仍可能把活動放在不同位置與時間。
# 第 9 章先固定規模比例，改變空間分配。
# 第 10 章再觀察事件如何成群出現。
# 讀者便能逐項看出模型增加了什麼資訊。
#
#
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("規模尺度與樣本", "HORUS Mw；學習期 R 內，名目規模 ≥ 2.5"),
    ("教學自估 b", f"離散估計 {b_exact:.3f}；95% 區間 [{whole_ci[0]:.3f}, {whole_ci[1]:.3f}]"),
    ("四張預報固定 b", f"{italy.SPEC.b_value:.3f}；Biondini 等（2023）"),
    ("完整度 Mc", "待填（第 8 章）"),
    ("一致性檢驗細節", "待填（第 16 章）"),
])))

# %% [markdown]
# 名目 2.5 是這次估計的起點。
# 它是否適合整段目錄，需要下一步檢查。
# 接著讀 {doc}`從哪個規模開始相信目錄 <08_completeness>`。
#
# ## 參考資料與延伸閱讀
#
# - Tinti 與 Gasperini（2024），[分箱規模的 b 值估計](https://cris.unibo.it/handle/11585/980514)。免費出版全文。先讀分箱估計式，再讀區間比較。
# - SeismoStats，[b 值估計文件](https://seismostats.readthedocs.io/latest/user/estimate_b.html)。免費方法說明。可對照連續與離散公式。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://www.earth-prints.org/handle/2122/17084)。免費機構全文。第 2 節交代 HORUS 資料。表 3 列出固定參數。
