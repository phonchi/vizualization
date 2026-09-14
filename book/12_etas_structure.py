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
# # 12. ETAS：把叢集寫成分支過程
#
# 一顆地震之後，常接著多顆地震。
# 第 11 章已把新事件接進條件強度。
# 本章再問：新事件能再引發下一批嗎？
#
# 一棵家族樹可以表達這件事。
# 分支過程（branching process）讓個體產生後代。
# 後代也遵循相同的繁殖規則。
# ETAS 用這個結構描述地震叢集。
# 其英文全名為 Epidemic-Type Aftershock Sequence。
# 中文可稱流行型餘震序列模型。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy, italy_models
from gdms_toolkit.viz import ACCENT, PALETTE, SEQUENTIAL, QUAKE_COLOR, apply_layout, plot_forecast_map, show_diagram

# %% [markdown]
# ## 12.1 從一條親子連線開始
#
# 樹上的一條線連接兩顆模擬事件。
# 親代（parent）是直接引發該事件的來源。
# 子代（offspring）是這條線另一端的事件。
# 世代（generation）記錄事件離根部幾層。
# 背景事件列為第零代。
#
# 第一代事件也能產生第二代。
# 因此，餘震可以成為其他地震的親代。
# ETAS 先讓每顆事件遵守相同規則。
# 模型再由規模決定平均產能。
# 模型無須先挑出最大的一顆。
#
# 目錄中的一列只有時間、位置與規模。
# HORUS 沒有記錄真實的親子連線。
# 下圖的線是生成資料時留下的資訊。
# 真實事件的親代須以機率表達。
# 第 13 章會計算這些機率。

# %% tags=["remove-input"]
show_diagram("d12_etas_branching", caption=f"分支示意：背景為第 {0} 代。連線表示模擬中指定的親子關係。")

# %% [markdown]
# ### 從 Hawkes 到帶規模的時空 ETAS
#
# 第11章的Hawkes只記事件時間。加入震央位置，就能描述某地震如何改變附近的率；再加入每顆事件的規模，便得到**帶標記的時空點過程**（marked spatio-temporal point process）。**標記**（mark）是事件附帶的觀測量，本章用地震規模作標記。
#
# ETAS是這類帶標記Hawkes的一個特例。它仍把過去事件的非負貢獻線性相加，但讓較大親代增加更多產能，也可影響空間範圍。Omori時間核、GR新事件規模密度與規模依賴產能，是ETAS在一般Hawkes架構中作出的具體選擇。
#
# 這裡有兩個規模角色。親代規模 $m_i$ 決定其觸發能力；新事件規模 $m$ 則從 $s(m)$ 抽取。標準版本的這次規模抽樣不依過去歷史，但抽到的規模會影響更後面的事件。不能把這個假設簡化成「規模與整份目錄無關」。
#
# 下式正是這條連結的終點：括號內是對規模積分後的時空強度，外面乘上新規模密度。閱讀論文看到 marked Hawkes 時，可以逐項對回同一套構造。
#
# ## 12.2 五個組成，各回答一個問題
#
# 下面的式子先寫背景，再加歷史貢獻。
# 背景率（background rate）描述獨立移入的活動。
# $\mu_0(x,y)$ 是指定的參考率密度。
# $\nu$ 調整背景的大小。
#
# $$
# \lambda^*(t,x,y,m)=s(m)\left[\nu\mu_0(x,y)
# +\sum_{t_i<t}\kappa(m_i)g(t-t_i)h_i(x,y)\right].
# $$
#
# 式中的五個組成如下。
# 密度（density）表示每單位範圍的份量。
# 密度乘上小區間，才近似該區間的份量。
#
# | 組成 | 回答的問題 | 讀法 |
# |---|---|---|
# | $\nu\mu_0$ | 沒有新增事件時呢？ | 背景提供基本率 |
# | $\kappa(m_i)$ | 直接後代有多少？ | 規模越大，產能越高 |
# | $g(\tau)$ | 後代何時出現？ | 時間核分配等待時間 |
# | $h_i(x,y)$ | 後代出現在哪裡？ | 空間核分配位置 |
# | $s(m)$ | 後代有多大？ | GR 密度分配規模 |
#
# 時間差 $\tau=t-t_i$ 以天計。
# 核（kernel）是分配事件貢獻的函數。
# 本章的時間核在正時間積分為一。
# 空間核在整個平面積分為一。
# 規模核在輸入門檻以上積分為一。
#
# $$
# \begin{aligned}
# \kappa(m_i)&=K e^{\alpha(m_i-m_0)},\\
# g(\tau)&=\frac{p-1}{c}(1+\tau/c)^{-p},\quad \tau>0,\\
# s(m)&=\beta e^{-\beta(m-m_0)},\quad \beta=b\ln10,\\
# d_i^2&=D^2e^{\gamma(m_i-m_0)},\\
# h_i(x,y)&=\frac{q-1}{\pi d_i^2}
# \left(1+\frac{(x-x_i)^2+(y-y_i)^2}{d_i^2}\right)^{-q}.
# \end{aligned}
# $$
#
# $D$ 的單位是公里；$c$ 的單位是天。
# $p>1$ 與 $q>1$ 讓兩個核可積分。
# $\alpha$ 控制產能隨規模增加的速度。
# $\gamma$ 控制空間尺度增加的速度。
# $b$ 控制大小地震的比例。
#
# 一顆大事件會增加更多總貢獻。
# 它的空間分布也可能更寬。
# 這兩件事分別由 $\alpha$ 與 $\gamma$ 控制。
# 較寬的分布會把貢獻分到較遠處。
# 中心峰值因此未必同時升高。

# %% tags=["remove-input"]
P = italy_models.PARAMS["ETAS"]
tau = np.geomspace(0.001, 365, 240)
r = np.linspace(0, 50, 240)
fig = make_subplots(rows=1, cols=2, subplot_titles=("等待時間密度", "距震央的面密度"))
for p, color in [(1.042, ACCENT), (1.3, PALETTE[2])]:
    g = (p-1)/P["c"] * (1+tau/P["c"])**(-p)
    fig.add_trace(go.Scatter(x=tau, y=g, name=f"p={p:g}", line_color=color), row=1, col=1)
for m, color in [(3., ACCENT), (5., PALETTE[2])]:
    d2 = P["D"]**2*np.exp(P["gamma"]*(m-2.5))
    h = (P["q"]-1)/(np.pi*d2)*(1+r*r/d2)**(-P["q"])
    fig.add_trace(go.Scatter(x=r, y=h, name=f"輸入規模 {m:.1f}", line_color=color), row=1, col=2)
fig.update_xaxes(type="log", title_text="等待天數", row=1, col=1)
fig.update_yaxes(type="log", title_text="每天", row=1, col=1)
fig.update_xaxes(title_text="距離（km）", row=1, col=2)
fig.update_yaxes(type="log", title_text="每平方公里", row=1, col=2)
apply_layout(fig, height=440, legend=dict(orientation="h", y=-0.25))
display(Markdown(f"時間核固定 c={P['c']:g} 天。空間核固定 q={P['q']:g}。右圖是面密度，環帶機率還須乘上環帶面積。"))
fig

# %% [markdown]
# ## 12.3 K 的數字，要連同公式一起讀
#
# Biondini 表 3 列出 $K=0.029$。
# 原文把時間項寫成 $K_B(\tau+c)^{-p}$。
# 此處用 $K_B$ 標示那套係數。
# 本站預報快取也採這套寫法。
#
# $$
# K=K_B\frac{c^{1-p}}{p-1}.
# $$
#
# 換算後的 $K$ 才是門檻事件的平均產能。
# 正規化（normalization）把密度積分調為一。
# 把同一數字直接搬到正規化核旁邊，
# 會改變整個模型的觸發量。
#
# 名目門檻也須一起核對。
# 快取使用取整規模 $2.5$ 作產能參考。
# 實驗規格的有效輸入邊界為 $2.45$。
# 下方換算固定採快取的參考規模。
# 比較論文時，要同時核對核與門檻。

# %% tags=["remove-input"]
beta = italy.SPEC.b_value*np.log(10)
k_normal = P["K"]*P["c"]**(1-P["p"])/(P["p"]-1)
n_reference = k_normal*beta/(beta-P["alpha"])
display(Markdown(f"依原文時間核換算，K={k_normal:.3f}。在無上限 GR 假設下，n={n_reference:.3f}。這是參數推算值，並非從本章目錄重新估計。"))

# %% [markdown]
# ## 12.4 平均一顆，會接出幾顆？
#
# 假設每顆事件平均產生半顆直接後代。
# 十顆事件便平均接出五顆第一代。
# 那五顆又平均接出兩顆半第二代。
# 半顆描述重複實驗的平均數。
# 每次抽到的事件數仍是整數。
#
# 分支比（branching ratio）記為 $n$。
# 它是對親代規模平均後的直接後代數。
# 將產能乘上規模密度，再積分即可得到。
#
# $$
# n=\int_{m_0}^{\infty}\kappa(m)s(m)\,dm
# =\frac{K\beta}{\beta-\alpha},\qquad \alpha<\beta.
# $$
#
# $n<1$ 稱次臨界（subcritical）。
# 每代的期望規模逐漸縮小。
# $n=1$ 稱臨界（critical）。
# $n>1$ 稱超臨界（supercritical）。
# 後者具有家族持續增長的可能。
#
# 一顆背景事件的全家族平均大小為
#
# $$
# 1+n+n^2+\cdots=\frac{1}{1-n},\qquad n<1.
# $$
#
# 這個和包含背景事件本身。
# 它涵蓋無限時間與整個空間。
# 三個月內、義大利區內的數量較受限。
# 因此，這個和只能幫助理解世代累積。
# 上一節由參數推算的 $n=1.579$ 屬超臨界，不能代入 $1/(1-n)$ 當作平均家族大小。
# 為了看清次臨界家族如何逐代縮小，下圖另用 $n=0.55$ 的合成設定。

# %% tags=["remove-input"]
rng = np.random.default_rng(120913)
n_sim, n_bins = 1, 10
memory_bytes = n_sim*n_bins*8
n_demo, initial = 0.55, 80
counts = [initial]
for _ in range(n_bins-1):
    counts.append(int(rng.poisson(n_demo*counts[-1])))
fig = go.Figure(go.Bar(x=np.arange(n_bins), y=counts, name="一次模擬", marker_color=ACCENT))
fig.add_trace(go.Scatter(x=np.arange(n_bins), y=initial*n_demo**np.arange(n_bins), name="理論平均", line_color=PALETTE[2]))
apply_layout(fig, xaxis_title="世代", yaxis_title="事件數", height=390)
display(Markdown(f"合成分支示例：起始 {initial} 顆，n={n_demo:.2f}。每顆採相同產能，只展示世代關係。"))
fig

# %% [markdown]
# ## 12.5 規模上限與「七釘二估」
#
# 當 $\alpha=\beta$，大規模端仍持續貢獻。
# 無上限積分因而發散。
# 若模型指定最大規模 $M_{\max}$，
# 就要一起改用截斷後的規模密度。
# 令 $L=M_{\max}-m_0$，可得
#
# $$
# n=\frac{K\beta L}{1-e^{-\beta L}},\qquad \alpha=\beta.
# $$
#
# 這個有限值依賴指定的上限。
# 截斷（truncation）保留指定範圍，並重新分配機率。
# 只截掉圖的右端，並未完成模型截斷。
# 模型須在抽樣與積分時使用相同上限。
#
# 九個參數可以分成兩組。
# simplETAS原文把固定形狀視為可檢驗的平均地殼地震模型。它採 $b=1$、$\alpha=\beta$，以自相似與經驗證據為動機；時間核選代表值，空間核則減少難以區分的參數補償。固定值並非宣稱全球每個序列都相同，仍須用未參與估計的事件檢查預報。
#
# 本站採用的是精簡估計的思路，形狀值改取Biondini義大利實驗。因此本章的「七釘二估」不等於執行simplETAS原作。原作的固定值、理由、背景地圖與本站差別列於附錄C；第13章的估計示範還會額外提高輸入門檻。
#
# 七個形狀參數是 $b,c,p,D,\gamma,\alpha,q$。
# 兩個大小參數是 $\nu,K$。
# 本書把固定前組、估計後組稱為「七釘二估」。
# 這借用 simplETAS 的精簡估計思路。
# 本書的固定值取自義大利研究。
# 這套示範與 simplETAS 的完整設定有別。
#
# 固定的七個值仍含研究假設。
# 資料較少時，固定形狀能減少估計負擔。
# 形狀若不合，兩個大小參數便會代為吸收差異。
# 下一章會把這個代價寫在估計圖旁。

# %% [markdown]
# ### 從一個格箱，讀回整條公式
#
# 一個方格只占空間核的一部分。
# 預報先對方格內的位置積分。
# 接著對三個月與目標規模箱積分。
# 三個積分共同決定該格箱的期望數。
#
# 輸入地震的規模可能低於目標門檻。
# 模型仍可讓它對目標規模箱貢獻少量率。
# 規模核決定這份貢獻有多少。
# 所以，輸入事件與預報目標分屬兩個集合。
#
# 圖上的每個點是一筆抽樣結果。
# 圖上的線則是反覆抽樣的平均。
# 平均線下降時，單次抽樣仍可能上升。
# 因此，一代突然增加可以出現在次臨界模型。
# 判斷分支比需要整體模型與資料。
#
# 世代圖只追蹤有限層數。
# 最後一欄之後仍可能有少量後代。
# 本圖把後面的世代省略以方便閱讀。
# 理論平均線也只畫出相同的層數。
#
# 一個背景事件可以恰好沒有後代。
# 另一個背景事件可能接出很大的家族。
# 平均家族大小並未指定每棵樹的大小。
# 這也是計數變異可能偏大的來源。
#
# $\nu$ 調高時，更多背景事件進入系統。
# $K$ 調高時，每個既有事件增加更多後代。
# 兩種改動都能提高總率。
# 但它們產生的時間與位置分布不同。
# 資料中的間隔與距離可幫助區分兩者。
#
# 一份預報若只列出期望數，
# 便還沒完整指定未來目錄的抽樣方式。
# 完整 ETAS 會讓新抽到的事件繼續繁殖。
# 固定率的 Poisson 抽樣則維持既定率。
# 兩者即使平均接近，計數分布也可能不同。
#
# 地圖邊緣的事件會把後代分到區外。
# 收集區內、測試區外的事件也可向內貢獻。
# 這兩個方向一起構成邊界效應。
# 本書保留收集區歷史來描述向內的貢獻。
# 格內積分則只留下測試區接到的部分。

# %% [markdown]
# ## 12.6 先看模型交出的預報
#
# 下圖把各窗的預報加成一張地圖。
# 本站 ETAS 近似採背景加第一代觸發。
# 第一代的親代限於各次發報前已知事件。
# 發報後才出生的後代，要等下次更新才加入。
#
# 藍色方格表示期望數，紅點表示目標地震。
# 深色格與紅點重疊，只是位置對照。
# 本章只疊圖不評分。
# 正式評分在第 16–18 章。

# %% tags=["remove-input"]
assert (italy.FORECAST_DIR/"ETAS_testing.npy").exists(), "請先提供預報快取"
forecast = italy_models.get_forecast("ETAS", "testing")
targets = italy.target_events(italy.experiment_catalog(), "testing")
fig = plot_forecast_map(forecast.sum(axis=(0,2)), targets=targets, title="ETAS：背景＋第一代，測試期合計")
apply_layout(fig)
display(Markdown(f"{forecast.shape[0]} 窗合計期望 {forecast.sum():.2f} 顆；疊上 {len(targets)} 顆目標地震。各窗歷史分別更新。只疊圖不評分，評分在第 16–18 章。"))
fig

# %% [markdown]
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
("ETAS 結構", "背景、產能、時間核、空間核、規模核"),
("本站預報差異", "背景＋發報前歷史的第一代觸發"),
("形狀參數", "採 Biondini 表 3；K 須連同時間核讀取"),
("兩參數估計示範", "待填（第 13 章）"),
("數量檢驗", "待填（第 16 章）"),
("位置與規模檢驗", "待填（第 17 章）"),
])))

# %% [markdown]
# $\nu$ 與 $K$ 的估計接在下一章：
# {doc}`ETAS 在義大利 <13_etas_italy_forecast>`。
#
# ## 參考資料與延伸閱讀
#
# - Jalilian（2019），[時空 ETAS 模型與估計](https://www.jstatsoft.org/article/view/v088c01)。免費全文。先讀條件強度，再讀親代機率。
# - Ogata 與 Zhuang（2006），[時空 ETAS 的空間核](https://doi.org/10.1016/j.tecto.2005.10.016)。[免費作者全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。對照規模如何改變影響範圍。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。[免費機構全文](https://www.earth-prints.org/handle/2122/17084)。表 3 提供參數。附錄列出模型公式。
# - Mancini 與 Marzocchi（2023 線上發表），[simplETAS](https://doi.org/10.1785/0220230199)。[免費作者程式](https://github.com/smancini2/simplETAS)。對照精簡參數的動機。
