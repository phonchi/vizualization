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
# # 13. 預報檢驗 II：新增資訊是否有幫助
#
# {doc}`一致性檢驗 <17_testing_consistency>`讓我們分清總量、位置和規模的偏差，但仍可能留下好幾個與資料相容的模型。現在把問題縮小：與原本使用的背景預報相比，加入 ETAS 觸發或 EEPAS 時間變化，對同一批未來事件有沒有幫助？
#
# 這是一個成對比較。兩個模型必須面對相同的區域、時間窗、規模門檻與目錄版本，否則分數差也可能來自題目不同。接著分別看資訊增益、校準與決策。「分數較好」「機率可信」和「對某個決定有用」雖然相關，卻不是相同問題。
#
# ## 13.1 用同一份觀測比較兩張預報
#
# 令候選模型為 $Z$，基準為 $1$。把上一章的對數概似相減，再除以觀測事件數 $N>0$，得到每事件資訊增益：
#
# $$\mathrm{IGPE}=\frac{\ln L_Z-\ln L_1}{N}
# =\frac{\Lambda_1-\Lambda_Z}{N}
# +\frac1N\sum_{n=1}^{N}\ln\frac{\lambda_Z(n)}{\lambda_1(n)}.$$ (eq:igpe)
#
# $\lambda_Z(n)$ 表示模型在第 $n$ 個事件處的率密度，兩模型須用相同測度；$\Lambda_Z$ 是對應總量積分，在固定 Poisson 預報中就是期望總數。第一項處理總量差，第二項比較真正有事件處的率。只保留第二項，可能讓到處抬高率的模型得到不合理優勢。
#
# 正值表示候選在這份資料上給出較高對數分數；負值表示基準較高。它永遠是相對於指定基準的結果。比均勻模型好，與比已經很有資訊的平滑地震度模型好，難度不同。
#
# 事件不同，每事件貢獻也不同，所以增益需要搭配不確定性。常見 $t$ 區間以事件貢獻的樣本變異數估計標準誤，近似把總量差當常數。對高度叢集資料，事件貢獻未必獨立，這個區間可能過窄；序列或時間區塊的重抽樣、以及模型模擬是可考慮的替代方法。
#
# 下圖用合成目錄示範點估計和區間。請先看區間是否跨零，再看寬度；跨零表示證據不足以判定方向，不表示已證明兩模型等價。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.special import gammaln

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

GREEN = "#1baf7a"
GREY = "#8a8a8a"
rng = np.random.default_rng(20240403)

NY = NX = 40                       # 40 × 40 空間格，共 1600 箱
WATER = 1e-8                       # water level：避免 log(0)（見 18.8）


def norm_field(a):
    """把非負場正規化成機率場（總和為 1）。"""
    return a / a.sum()


# 真實率場：對白雜訊做高斯平滑再取指數 → 高度不均勻但空間連續
z = gaussian_filter(rng.standard_normal((NY, NX)), 2.0, mode="wrap")
truth = norm_field(np.exp(1.6 * z / z.std()))

# 四個模型：基準是模糊掉的真實場，三個候選各自朝不同方向偏離
uniform = np.full((NY, NX), 1.0 / (NY * NX))
p_bench = norm_field(gaussian_filter(truth, 2.6, mode="wrap"))   # 基準：偏模糊
p_a = norm_field(0.85 * truth + 0.15 * p_bench)                  # 甲：保留細節
p_b = norm_field(p_bench * np.exp(rng.normal(0, 0.45, (NY, NX))))  # 乙：基準加雜訊
p_c = norm_field(0.60 * uniform + 0.40 * p_bench)                # 丙：摻均勻成分

# 各模型的期望總數略有差異 → 修正項 (N̂₁ − N̂_Z)/N 才有作用
N_TRUE = 130
NHAT = {"基準": 128.0, "模型甲": 121.0, "模型乙": 134.0, "模型丙": 118.0}
PROB = {"基準": p_bench, "模型甲": p_a, "模型乙": p_b, "模型丙": p_c}

n_obs = rng.poisson(N_TRUE)
cells = rng.choice(NY * NX, size=n_obs, p=truth.ravel())


def igpe(name, bench="基準", conf=0.95):
    """回傳 (IGPE, 半寬)：完全照 18.1 的式子實作。"""
    lam_z = NHAT[name] * PROB[name].ravel()[cells] + WATER
    lam_1 = NHAT[bench] * PROB[bench].ravel()[cells] + WATER
    d = np.log(lam_z) - np.log(lam_1)                     # dₙ = X_Z − X₁
    n = d.size
    s2 = d @ d / (n - 1) - d.sum() ** 2 / (n**2 - n)      # 變異數估計式
    ig = (NHAT[bench] - NHAT[name]) / n + d.mean()
    half = stats.t.ppf(0.5 + conf / 2, n - 1) * np.sqrt(s2 / n)
    return ig, half


names = ["模型甲", "模型乙", "模型丙"]
tt = {m: igpe(m) for m in names}
colors = [GREEN, PALETTE[3], QUAKE_COLOR]

fig = go.Figure()
for m, c in zip(names, colors):
    g, e = tt[m]
    fig.add_trace(go.Scatter(
        x=[m], y=[g], error_y=dict(type="data", array=[e], width=10),
        mode="markers+text", marker=dict(size=13, color=c),
        text=[f"{g:+.2f} ± {e:.2f}"], textposition="middle right",
        showlegend=False))
fig.add_hline(y=0, line_dash="dash", line_color=GREY,
              annotation_text="基準模型")
apply_layout(fig,
             title=(f"T-test 判讀：{n_obs} 顆目標地震，"
                    "誤差棒跨零表示尚無法判定增益方向"),
             yaxis_title="每地震資訊增益 IGPE（相對基準）",
             xaxis_title="", hovermode="closest",
             yaxis_range=[-0.55, 0.55])
fig

# %% [markdown]
# 模型間的差異可能比目錄帶來的不確定性還小。把格子切得更多，不會增加獨立地震序列；因此誤差棒的解讀要回到資料如何產生，而非只看網格大小。
#
# ## 13.2 比較基準與擬合代價
#
# 文獻中的資訊分數有不同基準與單位。Kagan 類資訊分數常以均勻空間模型為參考，使用以二為底的對數；IGPE 常用自然對數作成對比較。閱讀數字前，先對齊基準、總量處理及對數底數。
#
# 若參數是在同一批被評分資料上估計，較複雜的模型通常有更多改善擬合的機會。AICc 型修正以參數數量和樣本數近似扣除這種樂觀偏差。可寫為
#
# $$\mathrm{IGPEc}=\mathrm{IGPE}-\frac{c(k_Z,N)-c(k_1,N)}{N},\qquad
# c(k,N)=k+\frac{k(k+1)}{N-k-1},\quad N>k+1.$$ (eq:igpec)
#
# 這是指定 AICc 慣例下的修正，並非所有點過程、相依資料或邊界參數都自動符合其近似條件。它也無法替看過資料後挑區域、挑目標或挑模型的所有選擇付清代價。時間上獨立的測試與事先固定規則仍然重要。
#
# ## 13.3 如果使用者需要一個警報區域
#
# 有些問題關心有限資源應優先放在哪裡。我們可把預報率由高到低排序，逐漸擴大警報範圍，觀察漏掉多少事件。
#
# Molchan 圖橫軸是警報佔用的時空比例 $\tau$，縱軸是漏報比例 $\nu$。完全不警報位於 $(0,1)$，全部警報位於 $(1,0)$；曲線越靠左下，表示用較少警報範圍捕捉較多事件。
#
# 參考模型很關鍵。如果只用面積作分母，總在原本就很活躍的地方警報，可能看似很有技巧。用長期活動率加權警報範圍，才是在問是否超越既有空間資訊。下面同時看曲線與所用基準；不要只從某一個事後選定門檻報告最佳命中率。
#
# 圖中藍色面積是曲線上方、完整單位方框內的面積，記為 AS。它將不同警報範圍下的命中比例加以平均，值越大表示這種排序越能優先涵蓋事件；在圖中所選測度下，理論隨機參考線的完整面積為 0.5。
#
# %% tags=["remove-input"]
counts = np.bincount(cells, minlength=NY * NX)          # 每格的事件數
n_cells = NY * NX


def molchan(score):
    """回傳 (tau, nu)：按分數由大到小放寬門檻，掃出整條軌跡。"""
    order = np.argsort(-score.ravel(), kind="stable")
    hits = np.cumsum(counts[order])
    tau = np.arange(1, n_cells + 1) / n_cells
    nu = 1.0 - hits / counts.sum()
    return np.r_[0.0, tau], np.r_[1.0, nu]               # 補上端點 (0,1)


def area_skill(tau, nu):
    """AS(1) = ∫₀¹ [1 − ν] dτ，梯形法（numpy 2 已移除 trapz）。"""
    h = 1.0 - nu
    return float(np.sum(0.5 * (h[1:] + h[:-1]) * np.diff(tau)))



tau_a, nu_a = molchan(p_a)                               # 有技能的模型甲
tau_r, nu_r = molchan(rng.random(n_cells))               # 純亂數分數
as_a, as_r = area_skill(tau_a, nu_a), area_skill(tau_r, nu_r)
i20 = np.searchsorted(tau_a, 0.20)                       # τ = 0.2 的命中率
hit20 = 1.0 - nu_a[i20]

fig = go.Figure()
fig.add_trace(go.Scatter(x=np.r_[tau_a, 1.0, 0.0], y=np.r_[nu_a, 1.0, 1.0],
                         fill="toself", fillcolor="rgba(42,120,214,0.16)",
                         line=dict(width=0), hoverinfo="skip",
                         name=f"面積技能分數 AS = {as_a:.2f}"))
fig.add_trace(go.Scatter(x=tau_a, y=nu_a, mode="lines", name="模型甲",
                         line=dict(color=ACCENT, width=2.8)))
fig.add_trace(go.Scatter(x=tau_r, y=nu_r, mode="lines", name="亂數分數",
                         line=dict(color=PALETTE[3], width=1.6)))
fig.add_trace(go.Scatter(x=[0, 1], y=[1, 0], mode="lines",
                         name="隨機參考線 ν = 1 − τ",
                         line=dict(color=GREY, dash="dash")))
fig.add_trace(go.Scatter(
    x=[0.20], y=[nu_a[i20]], mode="markers+text", showlegend=False,
    marker=dict(size=11, color=QUAKE_COLOR),
    text=[f"  τ=0.20 時抓到 {hit20:.0%}"], textposition="middle right"))
apply_layout(fig,
             title=(f"Molchan 圖：藍色面積即 AS = {as_a:.2f}"
                    f"（亂數分數 {as_r:.2f}，理論隨機值 0.50）"),
             xaxis_title="警報涵蓋的時空體積比例 τ",
             yaxis_title="漏報率 ν", hovermode="closest",
             xaxis_range=[0, 1], yaxis_range=[0, 1])
fig

# %% [markdown]
# ## 13.4 沒看出差異，可能是樣本不夠
#
# 統計功效描述：在指定替代模型真的較好時，檢驗有多大機會辨認出來。與上一章的顯著水準相對，功效關心漏掉真差異的機會。
#
# 地震預報常面對少量目標大地震，以及同一序列中高度相關的事件。再細的地圖也無法創造缺少的觀測。下面先用合成資料改變事件數，看看相同模型差異如何影響被檢出的比例。
#
# 網格的選擇也會影響功效，但不能概括為越粗或越細越好。過細的箱可能讓檢驗對某些空間差異不敏感，過粗又會抹掉有用結構。Khawaja 等人（2023）提供地震空間檢驗的具體研究，提醒我們將解析度當成實驗設計，而非繪圖偏好。
#
# %% tags=["remove-input"]
LO, HI = -3.0, 3.0
CDF_LO, CDF_HI = stats.norm.cdf(LO), stats.norm.cdf(HI)


def bin_probs(n_cell, kind):
    """回傳 (真實機率 p_true, 均勻模型的預報機率 q_fore)。"""
    if kind == "equal_width":
        edges = np.linspace(LO, HI, n_cell + 1)
    else:                                   # equal_rate：真實分布的分位數
        edges = stats.norm.ppf(
            np.linspace(CDF_LO, CDF_HI, n_cell + 1))
    p_true = np.diff(stats.norm.cdf(edges))
    q_fore = np.diff(edges)                 # 均勻模型：機率正比於箱寬
    return p_true / p_true.sum(), q_fore / q_fore.sum()


def joint_poll(omega, lam):
    """Σ_j [−λ_j + ω_j ln λ_j − ln ω_j!]，omega 形狀 (nsim, ncell)。"""
    return (omega * np.log(lam) - gammaln(omega + 1.0)).sum(1) - lam.sum()


def sim_ll(n_ev, prob, lam, nsim, gen, chunk=400):
    """分塊模擬，避免 (nsim × ncell) 的大陣列吃光記憶體。"""
    out = np.empty(nsim)
    for s in range(0, nsim, chunk):
        m = min(chunk, nsim - s)
        out[s:s + m] = joint_poll(
            gen.multinomial(n_ev, prob, size=m).astype(float), lam)
    return out


def power(n_ev, n_cell, kind, nsim=1500, alpha=0.05, seed=7):
    """S-test 對「空間均勻模型」的統計功效。"""
    gen = np.random.default_rng(seed)
    p_true, q_fore = bin_probs(n_cell, kind)
    lam = n_ev * q_fore + WATER             # 率正規化到 N_obs（S-test 慣例）
    crit = np.quantile(sim_ll(n_ev, q_fore, lam, nsim, gen), alpha)
    return float(np.mean(sim_ll(n_ev, p_true, lam, nsim, gen) < crit))


N_LIST = [4, 8, 16, 32, 64, 128, 256]
pw_w = [power(n, 20, "equal_width") for n in N_LIST]
pw_r = [power(n, 20, "equal_rate") for n in N_LIST]


def first_reach(xs, ys, target=0.9):
    for x, y in zip(xs, ys):
        if y >= target:
            return x
    return None


n90_w, n90_r = first_reach(N_LIST, pw_w), first_reach(N_LIST, pw_r)
tail = (f"等期望率分箱在 {n90_r} 顆達 0.9"
        + (f"，等寬分箱要 {n90_w} 顆" if n90_w else "，等寬分箱始終未達"))

fig = go.Figure()
fig.add_trace(go.Scatter(x=N_LIST, y=pw_r, mode="lines+markers",
                         name="等期望率分箱（資料驅動）",
                         line=dict(color=GREEN, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=N_LIST, y=pw_w, mode="lines+markers",
                         name="等寬分箱（標準網格）",
                         line=dict(color=QUAKE_COLOR, width=2.8),
                         marker=dict(size=8)))
fig.add_hline(y=0.9, line_dash="dot", line_color=GREY,
              annotation_text="功效 0.9")
apply_layout(fig,
             title=f"同一批資料、同一個檢驗，只換分箱方式：{tail}",
             xaxis_title="事件數 N（20 個箱）", yaxis_title="統計功效",
             xaxis_type="log", yaxis_range=[0, 1.02], hovermode="x")
fig

# %% [markdown]
# 下面再比較不同網格下的效果。兩張圖都使用指定模擬設定，沒有一個樣本數門檻可直接套到所有地區。網格應用訓練資料或事先規則決定，不能看過未來震央才調到最容易顯著的位置。
#
# %% tags=["remove-input"]
CELL_LIST = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
N_FIX = 64
pc_w = [power(N_FIX, c, "equal_width", nsim=1200) for c in CELL_LIST]
pc_r = [power(N_FIX, c, "equal_rate", nsim=1200) for c in CELL_LIST]
k_peak = int(np.argmax(pc_w))

fig = go.Figure()
fig.add_trace(go.Scatter(x=CELL_LIST, y=pc_r, mode="lines+markers",
                         name="等期望率分箱（資料驅動）",
                         line=dict(color=GREEN, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=CELL_LIST, y=pc_w, mode="lines+markers",
                         name="等寬分箱（標準網格）",
                         line=dict(color=QUAKE_COLOR, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(
    x=[CELL_LIST[k_peak]], y=[pc_w[k_peak]], mode="markers+text",
    showlegend=False, textposition="middle right",
    marker=dict(size=13, color=QUAKE_COLOR, symbol="circle-open",
                line=dict(width=3)),
    text=[f"  峰值 {pc_w[k_peak]:.2f} @ {CELL_LIST[k_peak]} 箱"]))
apply_layout(fig,
             title=(f"解析度愈高，檢驗愈弱：固定 {N_FIX} 顆事件，"
                    "等寬分箱的功效達峰後單調下滑"),
             xaxis_title="箱數 N_cell", yaxis_title="統計功效",
             xaxis_type="log", yaxis_range=[0, 1.02], hovermode="x")
fig

# %% [markdown]
# ## 13.5 說百分之五時，長期是否真的接近百分之五
#
# 資訊增益較高的模型，仍可能系統性高估機率。校準問的是另一件事：在許多次預報中，所有被報為相近機率的情況，實際發生比例是否相近？
#
# 可靠度圖把預報機率分組，橫軸畫各組預報平均，縱軸畫實際比例。接近對角線是相容的校準表現；樣本少的組則應有很大的不確定性，尤其地震高機率情境本來就稀少。
#
# 二元目標常用 Brier 分數，$S_B=K^{-1}\sum_{i=1}^{K}(P_i-y_i)^2$，越小越好。它對極端錯誤的懲罰比對數分數溫和；兩者衡量的敏感度不同，可以互相補充。
#
# 下面比較校準、過度自信及資訊不足的示意預報。沿對角線並不是唯一目標：永遠報長期平均的模型可能校準，卻分不出何時較危險。好的預報還要在資料支持下，區分不同情境。
#
# %% tags=["remove-input"]
K_BIN = 20                                   # 預報值離散成 20 個格點
N_CASE = 40000
centers = (np.arange(K_BIN) + 0.5) / K_BIN
gen = np.random.default_rng(1999)


def snap(x):
    """把機率吸附到 20 個格點之一（讓 Murphy 分解成為精確恆等式）。"""
    return centers[np.clip((x * K_BIN).astype(int), 0, K_BIN - 1)]


truth_p = snap(gen.beta(1.5, 5.0, N_CASE))   # 各場合的真實機率
y = (gen.random(N_CASE) < truth_p).astype(float)
s_bar = y.mean()

MODELS = {
    "校準良好": truth_p,
    "過度自信": snap(np.clip(s_bar + 2.2 * (truth_p - s_bar), 0.002, 0.998)),
    "過度保守": snap(np.clip(s_bar + 0.40 * (truth_p - s_bar), 0.002, 0.998)),
}
MCOLOR = {"校準良好": GREEN, "過度自信": QUAKE_COLOR, "過度保守": PALETTE[3]}


def brier_parts(P, y):
    """回傳 (S_B, REL, RES, UNC)；P 已吸附到格點，故分解為恆等式。"""
    sb = float(np.mean((P - y) ** 2))
    rel = res = 0.0
    for v in np.unique(P):
        sel = P == v
        w, o = sel.mean(), y[sel].mean()
        rel += w * (v - o) ** 2
        res += w * (o - y.mean()) ** 2
    return sb, rel, res, float(y.mean() * (1 - y.mean()))


fig = go.Figure()
fig.add_trace(go.Scatter(x=[0, 0.8], y=[0, 0.8], mode="lines",
                         name="完美校準", line=dict(color=GREY, dash="dash")))
lines = []
for name, P in MODELS.items():
    sb, rel, res, unc = brier_parts(P, y)
    assert abs(sb - (rel - res + unc)) < 1e-10          # Murphy 恆等式
    vs = [v for v in np.unique(P) if (P == v).sum() >= 50]   # 樣本太少不畫
    fig.add_trace(go.Scatter(
        x=vs, y=[y[P == v].mean() for v in vs], mode="lines+markers",
        name=name, line=dict(color=MCOLOR[name], width=2.6),
        marker=dict(size=7)))
    lines.append(f"{name} S_B={sb:.4f}")
fig.add_hline(y=s_bar, line_dash="dot", line_color=GREY,
              annotation_text=f"長期基準率 {s_bar:.2f}")
apply_layout(fig,
             title="可靠度圖與 Brier 分數：" + "；".join(lines),
             xaxis_title="預報機率", yaxis_title="實際發生比例",
             hovermode="closest", xaxis_range=[0, 0.85],
             yaxis_range=[0, 0.85])
fig

# %% [markdown]
# ## 13.6 預報如何改變一個決定
#
# 假設某種防護行動每次花費 $C$，若不行動而事件發生會損失 $L$；先採一個簡化模型，假設行動可以避免這項損失。當事件機率為 $P$，不行動的預期損失為 $PL$，因此 $P>C/L$ 時行動較合算。
#
# 這個門檻沒有普遍的固定數字。低成本的準備工作與高成本的撤離，面對同一個預報可能採取不同決定；差別來自後果與成本，並非機率互相矛盾。
#
# 價值分數 $V$ 比較使用預報與只用長期平均時的支出，再以完美預報能節省的幅度作尺度。$V=0$ 表示與長期基準持平，$V=1$ 表示達到完美預報的節省幅度，負值則表示使用預報反而支出更多。
#
# 下面比較不同 $C/L$ 的相對價值。它用來理解預報技巧如何轉成決策價值，不能直接當成現實撤離規則。真實情境還有防護效果不完全、行動延遲、不同人承擔不同損失等因素，需要使用者與專業機構共同評估。
#
# %% tags=["remove-input"]
ALPHAS = np.linspace(0.01, 0.99, 99)


def value_curve(P, y, alphas):
    """對每個 α = C/L，用最佳門檻 p* = α 行動，回傳價值分數 V。"""
    s = y.mean()
    out = np.empty_like(alphas)
    for i, a in enumerate(alphas):
        act = P > a                                   # 18.7.1 的最佳門檻
        e_f = a * act.mean() + np.mean((~act) & (y > 0))
        e_cl, e_pf = min(a, s), a * s
        out[i] = (e_cl - e_f) / (e_cl - e_pf)
    return out


fig = go.Figure()
peaks = []
for name, P in MODELS.items():
    v = value_curve(P, y, ALPHAS)
    j = int(np.argmax(v))
    peaks.append(f"{name} 峰值 {v[j]:.2f}／全段平均 {v.mean():+.2f}")
    fig.add_trace(go.Scatter(x=ALPHAS, y=v, mode="lines", name=name,
                             line=dict(color=MCOLOR[name], width=2.8)))
fig.add_hline(y=0, line_dash="dash", line_color=GREY,
              annotation_text="與長期基準持平")
fig.add_vline(x=s_bar, line_dash="dot", line_color=GREY,
              annotation_text=f"C/L = 基準率 {s_bar:.2f}")
apply_layout(fig,
             title="成本–損失價值曲線：" + "；".join(peaks),
             xaxis_title="成本損失比 C/L（＝最佳行動門檻 p*）",
             yaxis_title="價值分數 V", yaxis_range=[-0.4, 0.55],
             hovermode="x")
fig

# %% [markdown]
# ## 13.7 把評分放回預報問題
#
# 我們現在有三種互補的讀法：資訊增益比較同一資料上的相對分數，可靠度檢查機率是否與長期比例相容，成本損失模型則問資訊能否改變特定決定。沒有一個分數能代替其他所有問題。
#
# 評估前還須對齊目錄版本、規模尺度、深度和區域。如果事件落在零率格，對數分數會是負無窮；若模型要設背景底線，底線必須在測試前固定並計入總量，不能事後為了避免低分而補上。
#
# ```{admonition} 想看評分與決策的推導
# :class: dropdown
#
# {doc}`附錄 E <appendix_e_testing>`從概似差推到 IGPE，列出其區間近似、AICc 條件、Molchan 面積分數、Brier 分解與成本損失代數。特別注意零事件時間窗不能定義每事件增益，以及事件相依性對標準誤的影響。
# ```
#
# 經過比較，未必會出現一個在所有尺度上都最好的模型。某模型掌握短期叢集，另一個較能描繪長期空間結構，可能都有保留的理由。{doc}`模型組合 <19_ensembles>`將從這種互補性出發，但仍須注意：學出來的權重，也必須接受未來資料的檢驗。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Theory of CSEP Tests：Forecast comparison tests](https://docs.cseptesting.org/getting_started/theory.html#forecast-comparison-tests) — pyCSEP 開發團隊，官方文件（免費）。從資訊增益 IGPE 的定義與 T-test 範例開始，對照本章的基準模型、率修正與信賴區間，理解「分數較高」與「差異顯著」的區別。
# - [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。接著閱讀警報式預測與誤差圖的部分，比較不同評估方法回答的問題，銜接本章的 Molchan 圖與面積技能分數。
# - [Statistical power of spatial earthquake forecast tests](https://doi.org/10.1093/gji/ggad030) — Asim M. Khawaja 等（2023），*Geophysical Journal International*；[免費機構典藏全文](https://gfzpublic.gfz.de/pubman/item/item_5015770_1)。以空間檢驗說明樣本量與網格如何改變統計功效，正好延伸本章「透過檢驗不等於模型有辨識力」的討論；可接著比較等寬網格與 Quadtree 的設計。
# - [Enhancing the Statistical Evaluation of Earthquake Forecasts—An Application to Italy](https://doi.org/10.1785/0220240209) — Jonas R. Brehmer、Kristof Kraus、Tilmann Gneiting、Marcus Herrmann、Warner Marzocchi（2025），*Seismological Research Letters*（出版社全文可能需訂閱；[免費作者稿](https://arxiv.org/abs/2405.10712)）。從義大利預報案例延伸本章的模型比較、可靠度圖與校準，重點是評分函數如何對應預報目標，以及如何分開檢查校準與鑑別力。
