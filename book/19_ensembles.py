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
# # 14. 模型組合：利用不同尺度的資訊
#
# {doc}`上一章 <18_testing_comparison>`讓我們看到，模型的價值取決於它回答的問題與比較基準。若 ETAS 在地震剛發生後較敏感，EEPAS 提供較長尺度的變化，PPE 保留背景空間資訊，我們是否能把它們放在同一份預報中？
#
# 組合的動機是互補，而不是模型越多越好。兩個幾乎一樣的模型平均後，仍提供幾乎相同的資訊；兩個不同的模型也可能共享同一缺陷。這一章先從加權平均開始，再理解乘法修正，最後分清楚率的平均與完整機率分布的混合。
#
# ## 14.1 同一段歷史，不同時間尺度
#
# 下圖用合成曲線畫出短期尖峰、較慢的中期變化與兩者的平均。圖中的目標地震和高峰位置是示意設定，不能視為模型已提前預知主震。請看新事件出現後，不同方法各在什麼時間尺度調整率。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import minimize

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

GREEN = "#1baf7a"
GREY = "#8a8a8a"

t = np.linspace(0, 10, 1000)                     # 年；目標大地震發生於 t = 8
bg = 0.010
short = np.full_like(t, bg)                      # 短期模型：前震觸發的尖峰
for tf, amp in [(7.62, 0.7), (7.87, 2.5)]:       # 兩個前震
    m = t > tf
    short[m] += amp / ((t[m] - tf) * 365 + 3) ** 1.1 * 30
medium = bg + 0.12 * np.exp(-0.5 * ((t - 8.3) / 1.8) ** 2)   # 中期模型：緩坡
mix = 0.5 * short + 0.5 * medium

fig = go.Figure()
for y, name, color, dash in [(short, "短期模型（STEP／ETAS 型）", PALETTE[1], None),
                             (medium, "中期模型（EEPAS 型）", GREEN, None),
                             (mix, "五五凸組合", ACCENT, "dash")]:
    fig.add_trace(go.Scatter(x=t, y=y, mode="lines", name=name,
                             line=dict(color=color, dash=dash, width=2)))
fig.add_vline(x=8, line_dash="dot", line_color=QUAKE_COLOR,
              annotation_text="目標大地震")
apply_layout(fig, title="同一個地震、兩種機率軌跡（示意）：尖峰 vs 緩坡",
             xaxis_title="時間（年）", yaxis_title="發生率密度（相對值）",
             yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 加權平均讓預報保留兩種時間變化，但也會降低各自最突出的尖峰。這樣的折衷是否有用，要看實際事件如何分布。
#
# 下面改用很多時空格比較兩模型率比。有些格子兩者相近，少數格子卻差異很大。合成資料用來展示這種分歧；這些分歧讓我們有機會學習如何組合，但還不能證明組合已帶來有用的資訊。
#
# %% tags=["remove-input"]
gen = np.random.default_rng(1973)
NC = 20000                                   # 兩萬個時空格
z1, z2 = gen.standard_normal(NC), gen.standard_normal(NC)
burst = np.maximum(0.0, gen.standard_normal(NC) - 1.0)      # 少數格子在叢集中
lam_short = np.clip(10.0 ** (-5.0 + 1.30 * z1 + 2.40 * burst), 1e-9, None)
lam_long = np.clip(10.0 ** (-5.0 + 1.05 * (0.5 * z1 + 0.87 * z2)), 1e-9, None)
logratio = np.log10(lam_short / lam_long)

span = logratio.max() - logratio.min()
frac_mid = float(np.mean(np.abs(logratio) < 1.0))

fig = go.Figure(go.Histogram(x=logratio, nbinsx=90,
                             marker_color=ACCENT, opacity=0.85,
                             name="率比"))
for xv, lab, col in [(-1.0, "比值 0.1", GREY), (1.0, "比值 10", GREY)]:
    fig.add_vline(x=xv, line_dash="dash", line_color=col,
                  annotation_text=lab)
apply_layout(fig,
             title=(f"兩模型率比的分布：跨越 {span:.1f} 個數量級，"
                    f"{frac_mid:.0%} 落在 0.1–10 之間"),
             xaxis_title="log₁₀（短期模型率 ÷ 長期模型率）",
             yaxis_title="時空格數", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# ## 14.2 用非負權重平均
#
# 最簡單的組合為 $\lambda_\pi=\sum_i\pi_i\lambda_i$，其中 $\pi_i\ge0$、$\sum_i\pi_i=1$。逐格看，組合率介於成分最小與最大值之間；對全區積分，總量也滿足 $\Lambda_\pi=\sum_i\pi_i\Lambda_i$。
#
# 若每個成分原本就有相同總量，調權重便只改變率的分配，不改變預期總數。不過相同總數不等於總數正確，也不保證 N-test 透過。若權重會隨時間或位置變動，上述總量結論則需重新檢查，不能把變動權重直接移出積分。
#
# 可以先用等權重作一個透明的基準，再用訓練資料估計權重。候選模型越多，未知權重也越多；有限的大地震樣本容易讓權重對某一段序列過度反應。
#
# ## 14.3 較好的權重可能在中間，也可能在端點
#
# 以兩模型為例，$\lambda_r=(1-r)\lambda_1+r\lambda_2$，$0\le r\le1$。對固定資料的點過程對數概似而言，這個目標函數是凹的，便於尋找全域最優；但凹性不保證最大值在中間。若模型二在所有觀測事件處都不提供有效補充，最佳權重可以是零。
#
# 真正的互補是：模型二雖然整體分數較差，仍在模型一很弱的部分給出有用資訊。這時平均後可能改善對數分數，因為它減少了一些事件落在極低率區的情況。
#
# 下面用兩個等總量模型產生合成資料，畫出權重與概似的關係。這組例子的最大值在內部；它示範一種可能情況，並不是任何兩模型混合都會更好的定理。端點導數條件與退化情況見 {doc}`附錄 E <appendix_e_testing>`。
#
# %% tags=["remove-input"]
rng = np.random.default_rng(42)
tt = np.linspace(0, 100, 2000)                    # 100 年
LAM0 = 150.0                                      # 兩模型共同的期望總數

lam_A = 4 * (0.3 + 0.25 * np.sin(tt / 6) ** 2)    # 模型 A：緩慢起伏的背景
lam_B = 4 * (0.05 + 1.6 * np.exp(-0.5 * ((tt - 60) / 4) ** 2))   # 模型 B：叢集尖峰
lam_A *= LAM0 / np.trapezoid(lam_A, tt)
lam_B *= LAM0 / np.trapezoid(lam_B, tt)

lam_true = 0.7 * lam_A + 0.3 * lam_B              # 真實率：七三混合
lam_max = lam_true.max()                          # 由真實率抽目錄（thinning）
cand = rng.uniform(0, 100, rng.poisson(lam_max * 100))
ev = cand[rng.random(cand.size) < np.interp(cand, tt, lam_true) / lam_max]

rs = np.linspace(0, 1, 201)


def loglik(r):
    lam = (1 - r) * lam_A + r * lam_B
    return np.sum(np.log(np.interp(ev, tt, lam))) - np.trapezoid(lam, tt)


lnL = np.array([loglik(r) for r in rs])
k = int(np.argmax(lnL))
r_hat = float(rs[k])

# 19.3.1 的兩個端點條件，直接用資料驗算
rho = np.interp(ev, tt, lam_B) / np.interp(ev, tt, lam_A)
delta = float(np.mean(np.log(rho)))               # = 模型 B 相對 A 的 IGPE
assert rho.mean() > 1 and (1 / rho).mean() > 1    # ⇒ 最佳權重必在內部
assert 0 < r_hat < 1

fig = go.Figure(go.Scatter(x=rs, y=lnL - lnL[k], mode="lines",
                           line=dict(color=ACCENT, width=2.8),
                           name="ln L(r) − max"))
fig.add_vline(x=r_hat, line_dash="dash", line_color=GREEN,
              annotation_text=f"最佳權重 r = {r_hat:.2f}")
for xv, lab in [(0.0, f"純 A：{lnL[0] - lnL[k]:.1f}"),
                (1.0, f"純 B：{lnL[-1] - lnL[k]:.1f}")]:
    fig.add_trace(go.Scatter(x=[xv], y=[lnL[int(xv * 200)] - lnL[k]],
                             mode="markers+text", showlegend=False,
                             marker=dict(size=11, color=QUAKE_COLOR),
                             text=[lab], textposition="top center"))
apply_layout(fig,
             title=(f"凸組合的對數概似：{ev.size} 顆事件，"
                    f"B 單獨較差（IGPE = {delta:+.2f}）卻仍分到 "
                    f"{r_hat:.0%} 權重"),
             xaxis_title="模型 B 的權重 r（0 = 純 A，1 = 純 B）",
             yaxis_title="相對對數概似 ln L(r) − max", hovermode="x")
fig

# %% [markdown]
# 曲線的峰值告訴我們這份訓練目錄偏好的權重。它能比兩個端點高，雖然組合在任何一格都沒有超出成分率的範圍：逐格率的界線不等於整份非線性分數的界線。
#
# Rhoades 與 Gerstenberger（2009）在加州資料上混合 STEP 與 EEPAS，也發現訓練期的最佳權重保留兩者。這提供繼續測試的理由，還不能把該權重當成其他區域或未來期間的固定答案。
#
# ## 14.4 讓新資訊修正既有背景
#
# 如果我們已有一個相當好的背景地圖，也可以讓其他資料只回答「哪裡該提高、哪裡該降低」。乘法組合以正修正因子調整基準率，常見形式為
#
# $$\lambda_H(j,k)=\lambda_1(j,k)
# \exp\!\left[a_0+\sum_{i=2}^{I}f_i\bigl(\lambda_i(j,\cdot)\bigr)\right].$$ (eq:mulhyb)
#
# 這裡 $j,k$ 是空間及規模箱，$\lambda_i(j,\cdot)$ 對規模求和，$f_i$ 把其他模型的空間訊號轉成修正，$a_0$ 控制整體總量。不同論文對轉換與約束的定義不同，不能只對照參數名稱。
#
# 乘法可以把局部的率提高到所有成分之上，也能壓得更低，因而需要注意極低率與正規化。若新增訊號只是背景已經包含的資訊，再乘一次可能重複強化同一結構。
#
# 下圖比較兩種運算的幾何效果。加法平均保留折衷，乘法修正則更容易集中；集中是否有幫助，必須等真正的觀測出現才知道。
#
# %% tags=["remove-input"]
xs = np.linspace(0, 100, 801)                     # 一條空間剖線


def gauss(c, w, h):
    return h * np.exp(-0.5 * ((xs - c) / w) ** 2)


lam1 = 0.020 + gauss(25, 7, 0.55) + gauss(70, 9, 0.30)   # 基準（目錄平滑）
lam2 = 0.015 + gauss(45, 8, 0.65) + gauss(72, 6, 0.22)   # 共軛（大地測量）
A0, U, V = -2.9, 12.0, 0.70                       # 擬合期得到的乘法參數
lam_add = 0.5 * lam1 + 0.5 * lam2                 # 加法：五五凸組合
lam_mul = lam1 * np.exp(A0 + U * np.log1p(lam2) ** V)     # 乘法 hybrid

lo, hi = np.minimum(lam1, lam2), np.maximum(lam1, lam2)
assert np.all((lam_add >= lo - 1e-12) & (lam_add <= hi + 1e-12))   # 逐格內插
below = lam_mul < lo
j_worst = int(np.argmax(np.log10(lo / lam_mul)))
fac_below = float(lo[j_worst] / lam_mul[j_worst])
fac_above = float(np.max(lam_mul / hi))

fig = go.Figure()
fig.add_trace(go.Scatter(x=np.r_[xs, xs[::-1]], y=np.r_[hi, lo[::-1]],
                         fill="toself", fillcolor="rgba(138,138,138,0.18)",
                         line=dict(width=0), hoverinfo="skip",
                         name="兩成分之間（加法的值域）"))
for y, name, col, dash in [(lam1, "基準模型 λ₁", PALETTE[1], None),
                           (lam2, "共軛模型 λ₂", GREEN, None),
                           (lam_add, "加法：0.5 λ₁ + 0.5 λ₂", ACCENT, "dash"),
                           (lam_mul, "乘法 hybrid", QUAKE_COLOR, None)]:
    fig.add_trace(go.Scatter(x=xs, y=y, mode="lines", name=name,
                             line=dict(color=col, dash=dash, width=2.4)))
fig.add_trace(go.Scatter(
    x=[xs[j_worst]], y=[lam_mul[j_worst]], mode="markers+text",
    showlegend=False, marker=dict(size=13, color=QUAKE_COLOR,
                                  symbol="x", line=dict(width=2)),
    text=[f"  目標地震落在這裡：比兩者都低 {fac_below:.1f} 倍"],
    textposition="middle right"))
apply_layout(fig,
             title=(f"加法逐格內插、乘法可外推："
                    f"乘法在 {below.mean():.0%} 的格子低於兩個成分，"
                    f"最高處則超出 {fac_above:.1f} 倍"),
             xaxis_title="空間剖線位置（格）",
             yaxis_title="格點期望數 λ", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# ## 14.5 回溯改善能不能延續
#
# 如果組合模型的參數空間包含原基準，最大化同一份資料的概似，理論上就不會比基準差。因此訓練分數上升本身不是充分證據。
#
# Bayona 等人（2022）將乘法 hybrid 放進加州的獨立前瞻期，發現先前回溯優勢並未普遍延續。下面保留文獻比較的重點：先分清擬合與測試的時期，再看增益與不確定性，不能只挑最高的一根柱子。
#
# %% tags=["remove-input"]
labels = ["HKJ ⊗ Neokinema", "HKJ ⊗ PI", "HKJ ⊗ Shen"]
retro = [0.25, 0.25, 0.50]
prosp = [-0.42, -0.71, -0.68]

fig = go.Figure()
fig.add_trace(go.Bar(x=labels, y=retro, name="回溯（2006–2010，擬合期內，IGPEc）",
                     marker_color=GREEN,
                     text=[f"{v:+.2f}" for v in retro], textposition="outside"))
fig.add_trace(go.Bar(x=labels, y=prosp, name="前瞻（2011–2020，獨立測試，IGPE）",
                     marker_color=QUAKE_COLOR,
                     text=[f"{v:+.2f}" for v in prosp], textposition="outside"))
fig.add_hline(y=0, line_color=GREY,
              annotation_text="基準模型 HKJ")
apply_layout(fig,
             title=("乘法 hybrid 相對基準的每地震資訊增益："
                    "同一批模型、同一個測試區，符號全數翻轉"),
             yaxis_title="每地震資訊增益（相對 HKJ）",
             xaxis_title="", hovermode="x", barmode="group",
             yaxis_range=[-0.95, 0.75])
fig

# %% [markdown]
# 乘法模型若把率集中到某些歷史熱區，未來事件轉移到別處就可能承受較大損失。乘法不一定較差，但模型越靈活，越需要足夠資料與獨立檢驗。
#
# ## 14.6 權重究竟在學什麼
#
# 一種權重策略依每個模型自己的分數配比；另一種直接最佳化整體組合。後者會注意某模型能否補足其他模型，而不只看單獨排名。Herrmann 與 Marzocchi（2023）提供以 logistic 迴歸資訊建立組合的具體方法。
#
# 這個連結可以從勝算理解：若事件的 log-odds 對幾個對數率呈線性，取指數後就是各個率的乘冪。因此 logistic 係數可以描述訊號的條件貢獻。但把係數轉為非負權重是特定方法的設計，不是所有 logistic 模型都必須丟棄截距或負係數。
#
# 若每隔一段時間重新學權重，還要決定保留多久的歷史。短時間窗較能反映近期活動，也更容易受單一序列影響；長時間窗較穩定，卻可能混入已改變的背景。下面用模擬展示不同記憶長度的權重變化。更新規則本身也是模型的一部分，必須事先固定和測試。
#
# %% tags=["remove-input"]
gen = np.random.default_rng(20260403)
NCELL, NWIN = 240, 416                            # 240 格 × 416 週（8 年）
years = np.arange(NWIN) / 52.0
ug = np.linspace(0, 1, NCELL)
BASE = 1.3                                        # 每視窗期望事件數


def cell_bump(c, w, h=1.0):
    return h * np.exp(-0.5 * ((ug - c) / w) ** 2)


bg_true = 0.35 + cell_bump(0.30, 0.06) + cell_bump(0.72, 0.09, 0.7)
bg_wrong = 0.35 + cell_bump(0.38, 0.11) + cell_bump(0.66, 0.16, 0.7)
bg_true /= bg_true.sum()
bg_wrong /= bg_wrong.sum()
BURSTS = [(2.1, 0.30), (5.4, 0.71)]               # （年，位置）

phi_A = np.tile(BASE * bg_true, (NWIN, 1))        # 模型甲：長期平滑，時間不變
phi_B = np.empty((NWIN, NCELL))                   # 模型乙：叢集模型
lam_t = np.empty((NWIN, NCELL))                   # 真實率
for w in range(NWIN):
    clust = np.zeros(NCELL)
    for t_m, x0 in BURSTS:                        # Omori 型衰減的叢集
        if years[w] > t_m:
            clust += (8.0 / ((years[w] - t_m) * 365 + 4.0) ** 1.1
                      * cell_bump(x0, 0.025))
    tot = clust.sum()
    th = tot / (tot + 1.2)                        # 叢集成分佔比
    shape = clust / tot if tot > 0 else np.zeros(NCELL)
    phi_B[w] = BASE * ((1 - th) * bg_wrong + th * shape)
    lam_t[w] = BASE * ((1 - th) * bg_true + th * shape)

Y = (gen.poisson(lam_t) > 0).astype(float)        # 二元觀測
U = np.stack([np.log(phi_A + 1e-6), np.log(phi_B + 1e-6)], axis=-1)


def fit_beta(sl):
    """無截距的二元 logistic：最小化負對數概似（解析梯度）。"""
    X = U[sl].reshape(-1, 2)
    yv = Y[sl].ravel()

    def nll(b):
        return float(np.sum(np.logaddexp(0.0, X @ b) - yv * (X @ b)))

    def grad(b):
        return X.T @ (1.0 / (1.0 + np.exp(-(X @ b))) - yv)

    return minimize(nll, np.array([0.5, 0.5]), jac=grad,
                    method="L-BFGS-B").x


def to_weights(beta, tau=0.0):
    """19.7.2 的映射：w = exp(β) − exp(τ)，再正規化。"""
    w = np.where(beta > tau, np.exp(beta) - np.exp(tau), 0.0)
    return w / w.sum() if w.sum() > 0 else np.full(beta.size, 0.5)


out = list(range(52, NWIN, 3))                    # 第一年之後開始輸出
pi_all = np.array([to_weights(fit_beta(slice(0, w))) for w in out])
pi_win = np.array([to_weights(fit_beta(slice(w - 52, w))) for w in out])

fig = go.Figure()
fig.add_trace(go.Scatter(x=years[out], y=pi_all[:, 1], mode="lines",
                         name="方案 #1：全部歷史", line=dict(color=ACCENT, width=2.6)))
fig.add_trace(go.Scatter(x=years[out], y=pi_win[:, 1], mode="lines",
                         name="方案 #2：一年遺忘視窗",
                         line=dict(color=QUAKE_COLOR, width=2.6)))
for t_m, _ in BURSTS:
    fig.add_vline(x=t_m, line_dash="dot", line_color=GREY,
                  annotation_text="序列爆發")
fig.add_hline(y=0.5, line_dash="dash", line_color=GREY)
apply_layout(fig,
             title=(f"權重時序（叢集模型乙的份額）：遺忘視窗在 "
                    f"{pi_win[:, 1].min():.2f}–{pi_win[:, 1].max():.2f} "
                    f"之間擺盪，全歷史只從 {pi_all[0, 1]:.2f} 爬到 "
                    f"{pi_all[-1, 1]:.2f}"),
             xaxis_title="時間（年）", yaxis_title="模型乙的權重 π₂",
             yaxis_range=[-0.03, 1.05], hovermode="x")
fig

# %% [markdown]
# 權重跳動不一定表示物理機制突然切換，也可能只是資料不足。小事件可以提供更多樣本，但用它們學出的權重是否適用大目標，仍需要檢查規模尺度與完整度。
#
# ## 14.7 平均率會隱藏哪些不確定性
#
# 假設兩個模型對某格分別預報很低和很高的計數。只交出平均率，讀者看不出模型彼此分歧。若保留完整混合分布，則能同時表達各模型內的隨機波動與模型之間的差異。
#
# 對 Poisson 成分的分布混合，總變異數等於平均率加上成分均值的加權變異數。相反地，先把率平均再宣告是一個 Poisson，變異數只剩平均率。這兩種預報均值相同，尾部機率卻可能不同。
#
# 模型間差異不是全部認知不確定性：候選池可能漏掉重要機制，參數與權重也有估計誤差。組合能保留不同假設，但不能保證它們涵蓋真實地震行為。
#
# ```{admonition} 組合的數學性質
# :class: dropdown
#
# {doc}`附錄 E <appendix_e_testing>`推導凹性與內點條件、乘法修正的 KL 分解、混合 Poisson 的全變異數分解，以及 logistic 與乘冪的關係。尤其「未顯著分出勝負」不等於樣本分數恰好相同，不能用它保證混合有嚴格增益。
# ```
#
# 目前多數例子從區域地震目錄出發。如果改問一條特定斷層距離上次大破裂已經多久，我們就需要另一種時間記憶。{doc}`複發模型 <20_recurrence_models>`會把焦點轉到等待時間，讓長期地震發生模型與前面的點過程觀念接起來。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等（2024），*Reviews of Geophysics*（免費全文）。先讀模型發展與各國作業系統中組合模型的段落，瞭解為何需要結合不同預報，以及權重選擇如何與檢驗制度連在一起。
# - [Mixture Models for Improved Short-Term Earthquake Forecasting](https://doi.org/10.1785/0120080063) — David A. Rhoades、Matthew C. Gerstenberger（2009），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[SCEC 免費摘要](https://central.scec.org/node/3964)）。本章加法組合的原始案例，將 STEP 與 EEPAS 的不同時間尺度結合；重點是凸組合如何利用互補資訊，以及回溯估出的權重仍須接受前瞻檢驗。
# - [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。對照本章乘法 hybrid 的前瞻測試，觀察回溯資訊增益為何未必延續到未來；適合與前一篇的組合動機一起讀。
# - [Maximizing the forecasting skill of an ensemble model](https://doi.org/10.1093/gji/ggad020) — Marcus Herrmann、Warner Marzocchi（2023），*Geophysical Journal International*（免費全文）。本章 logistic 權重學習的主要來源；重點是直接最佳化整體組合的技巧，而非只依各成分模型的單獨成績配權重。
# - [Regional Earthquake Likelihood Models II: Information Gains of Multiplicative Hybrids](https://doi.org/10.1785/0120140035) — D. A. Rhoades、M. C. Gerstenberger、A. Christophersen、J. D. Zechar、D. Schorlemmer、M. J. Werner、T. H. Jordan（2014），*Bulletin of the Seismological Society of America*（出版社全文可能需訂閱；[Bristol 大學免費全文](https://research-information.bris.ac.uk/files/49915043/Rhoades_Hybrids_RELM_BSSA_2014.pdf)）。本章乘法 hybrid、保序轉換與正規化的直接來源；讀完建構方法後，對照上面的 Bayona 等人前瞻測試，區分擬合改善與未來預報表現。
#
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。
#
# - [Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759) — Mark Naylor、Francesco Serafini、Finn Lindgren、Ian G. Main（2023），Frontiers in Applied Mathematics and Statistics（免費全文）。對照模型間差異與單一模型內的參數不確定性；貝氏 ETAS 的後驗分布並不等同於模型組合權重。
#
