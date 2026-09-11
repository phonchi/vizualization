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
# # 10. Ψ：從震前活動看統計證據
#
# 在 {doc}`ETAS 估計 <14_etas_estimation>`，我們學會用觸發模型解釋事件之後的叢集，也知道一段活動是否異常，必須相對於某個參考模型來說。現在把觀察尺度拉長：在一場大地震之前數月到數年，周圍中小地震的活動是否包含額外資訊？
#
# 直覺上，人會想從震前紀錄找出徵兆。但先選定一場已發生的大地震，再回頭選時間窗和區域，很容易得到一段看似特殊的歷史。這一章用「前兆尺度增加」Ψ 作為具體案例，學習如何把觀察變成統計問題，並辨認哪些證據仍不足以支持預報。
#
# Ψ 指的是文獻在選定震前區域中辨識的規模水準與發生率增加。這個名稱描述一類目錄形態，不表示每場大地震都有可辨識徵兆，也不先排除 ETAS 型叢集的解釋。從量測圖開始，檢查搜尋和比較的方式，才能判斷下一章的 EEPAS 可以如何使用這些尺度關係。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 10.1 從規模時間圖讀出活動變化
#
# 先想像把同一區域的地震依時間排開。密集的點表示發生率較高，較高的點表示規模較大。兩種變化可以一起發生，也可以分開；只看事件數會失去規模資訊，只看最大規模又容易被一個事件主導。
#
# 累積規模異常 $C(t)$ 把兩種資訊放在同一張圖。給定分析窗 $[t_s,t_f]$，保留規模高於分析門檻 $m_c$ 的事件，定義
#
# $$C(t)=\sum_{t_s<t_i\le t}[M_i-(m_c-0.1)]-k(t-t_s),\qquad
# k=\frac{\sum_{t_s<t_i\le t_f}[M_i-(m_c-0.1)]}{t_f-t_s}.$$ (eq:cumag)
#
# 這裡 $m_c$ 是分析者選定的門檻，須高於目錄完整度 $M_c$。式中的 $0.1$ 沿用 Christophersen 等人（2024）的規模刻度；它不是對所有規模資料都適用的完整度修正。每個事件向上增加 $M_i-m_c+0.1$，兩事件之間則減去一段平均趨勢。時間若用年，$k$ 的單位就是每年的規模累積量。
#
# 所以，曲線持續下降表示累積速度低於整段平均，持續上升表示高於平均。頭尾回到零是定義的結果，不是物理系統恢復原狀。文獻把最小值位置定為一個候選起點；這個位置描述的是已選時窗內的相對變化，不能單靠它證明存在真正的變點。
#
# 圖中同時畫出事件跳升前後的值，候選起點取事件前的左極限；若起點恰好位於事件時刻，該事件歸入起點以後的區段。這樣分段計數才與曲線的量法一致。
#
# 以起點到主震的時間為 $T_P$，這段期間最大三個規模的平均為 $M_P$，包住選定事件的面積為 $A_P$，便把一張圖整理成三個尺度。先前期也取最大三個規模的平均，記為 $M^-$；率比 $r$ 則是前兆期的事件數除以期間長度，再除以先前期用相同方法算出的事件率。$M_P-M^-$ 比較規模水準，$r$ 比較發生頻率。
#
# 圖中 M.U. 是 magnitude units，即規模單位；$C(t)$ 累積的是規模超額，不是地震能量。以下是刻意加入活動變化的合成資料：先比較上方散佈圖的密度與高度，再看下方的 $C(t)$ 如何回應。
#
# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

rng = np.random.default_rng(15001)
MC = 4.0                                  # 分析用下限規模 m_c
T_S, T_M, T_F = 0.0, 20.0, 24.0           # 起點、主震、終點（年）
ONSET_TRUE = 13.0

t_prior = np.sort(rng.uniform(T_S, ONSET_TRUE, 22))          # 先前期
m_prior = MC + rng.exponential(0.30, t_prior.size)
t_pre = np.sort(rng.uniform(ONSET_TRUE, T_M, 27))            # 前兆期
m_pre = MC + 0.45 + rng.exponential(0.42, t_pre.size)
t_aft = T_M + np.sort(rng.exponential(0.55, 24))             # 餘震
m_aft = MC + rng.exponential(0.38, t_aft.size)
t_aft, m_aft = t_aft[t_aft < T_F], m_aft[t_aft < T_F]


def cumag(t, m, t0, t1, mc=MC):
    """回傳 (時間格點, C(t))；C 在每個事件前後各取一點以顯示跳躍。"""
    sel = (t > t0) & (t <= t1)
    te, me = t[sel], m[sel]
    exc = me - mc + 0.1
    k = exc.sum() / (t1 - t0)
    grid = np.concatenate([[t0], np.repeat(te, 2), [t1]])
    cum = np.concatenate([[0.0], np.repeat(np.cumsum(exc), 2)[:-1], [exc.sum()]])
    cum = np.concatenate([[0.0], np.repeat(np.concatenate([[0.0],
                          np.cumsum(exc)[:-1]]), 1)])
    # 重新以「事件前 / 事件後」兩點描出階梯
    grid = np.concatenate([[t0], np.repeat(te, 2), [t1]])
    before = np.concatenate([[0.0], np.cumsum(exc)[:-1]])
    cum = np.concatenate([[0.0], np.ravel(np.column_stack([before,
                          np.cumsum(exc)])), [exc.sum()]])
    return grid, cum - k * (grid - t0), k


t_all = np.concatenate([t_prior, t_pre])
m_all = np.concatenate([m_prior, m_pre])
g1, C1, k1 = cumag(t_all, m_all, T_S, T_M)                   # [t_s, t_M]
g2, C2, _ = cumag(t_aft, m_aft, T_M, T_F)                    # [t_M, t_f]
onset = g1[int(np.argmin(C1))]

pre_mask, pri_mask = t_all >= onset, t_all < onset
M_P = np.sort(m_all[pre_mask])[-3:].mean()
M_MINUS = np.sort(m_all[pri_mask])[-3:].mean()
r_ratio = (pre_mask.sum() / (T_M - onset)) / (pri_mask.sum() / (onset - T_S))
T_P_yr = T_M - onset

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.52, 0.48],
                    vertical_spacing=0.06)
fig.add_trace(go.Scatter(x=t_all[pri_mask], y=m_all[pri_mask], mode="markers",
                         name="先前期事件",
                         marker=dict(size=7, color=PALETTE[4], opacity=0.85)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=t_all[pre_mask], y=m_all[pre_mask], mode="markers",
                         name="前兆期事件",
                         marker=dict(size=8, color=ACCENT, opacity=0.9)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=t_aft, y=m_aft, mode="markers", name="餘震",
                         marker=dict(size=6, color=PALETTE[5], opacity=0.55)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=[T_M], y=[7.0], mode="markers", name="主震",
                         marker=dict(size=15, color=QUAKE_COLOR, symbol="star")),
              row=1, col=1)
for lvl, x0, x1, nm in ((M_MINUS, T_S, onset, "M⁻"), (M_P, onset, T_M, "M_P")):
    fig.add_trace(go.Scatter(x=[x0, x1], y=[lvl, lvl], mode="lines",
                             name=nm, line=dict(color="#666666", dash="dash",
                                                width=1.5)), row=1, col=1)
fig.add_trace(go.Scatter(x=g1, y=C1, mode="lines", name="C(t)：主震前",
                         line=dict(color=PALETTE[3], width=2)), row=2, col=1)
fig.add_trace(go.Scatter(x=g2, y=C2, mode="lines", name="C(t)：主震後",
                         line=dict(color=PALETTE[3], width=1.5, dash="dot")),
              row=2, col=1)
for r_ in (1, 2):
    fig.add_vline(x=onset, line_dash="dot", line_color="#1baf7a", row=r_, col=1)
fig.add_annotation(x=onset, y=float(C1.min()),
                   text=f"onset（T_P = {T_P_yr:.1f} 年）", showarrow=True,
                   ax=-55, ay=32, row=2, col=1)
fig.update_yaxes(title_text="規模", row=1, col=1)
fig.update_yaxes(title_text="C(t)（M.U.）", row=2, col=1)
fig.update_xaxes(title_text="時間（年）", row=2, col=1)
apply_layout(fig, height=560, hovermode="x",
             title=(f"Ψ 的量法（合成資料）：M_P = {M_P:.2f}、M⁻ = {M_MINUS:.2f}、"
                    f"M_P − M⁻ = {M_P - M_MINUS:.2f}、率比 r = {r_ratio:.1f}"))
fig

# %% [markdown]
# 下方曲線把前半段較低、後半段較高的活動合在一起；上方陰影則標出依最小值切出的兩段。圖中變化是模擬時設定的，因此它示範的是統計量的讀法，不是從真實資料證明震前訊號。
#
# 這也指出下一個問題：如果分析窗換了，平均趨勢就會改變，最小值可能跟著移動。必須連同選窗方法一起檢查，不能只保留最後最好看的曲線。
#
# %% [markdown]
# ## 10.2 分析區域為什麼會改變結論
#
# 假設先知道一場主震的位置，再向外擴大矩形；又把起始時間前後移動，直到圖上出現明顯的規模增加。即使每一次計算都正確，最後留下的圖也已經受過篩選。它回答「在哪一種選法下能看見變化」，而不是「未知主震到來之前能不能發現變化」。
#
# 這和前面估計章所說的選擇偏誤相同。多試幾種門檻、幾個時間窗，就多了幾次從雜訊中挑出形態的機會。若只把選中的結果與一次隨機序列比較，兩邊採用的分析程序不同，顯著性便沒有共同基準。
#
# Christophersen 等人（2024）將搜尋規則寫成演算法，讓真實目錄和對照目錄經過相同程序，方便重做搜尋並檢查選擇的影響。規則仍由人設定，而且仍以已知主震為條件，因此尚未變成前瞻預報。
#
# %% [markdown]
# ## 10.3 把搜尋規則固定下來
#
# 矩形搜尋的想法是依序改變區域大小和時間長度，對每個候選窗計算 $C(t)$，再檢查起點、規模增加與事件率是否符合事先設定的條件。符合者全部留下，而不只挑一個最清楚的例子。
#
# 這個做法讓「找得到多少組」本身也成為結果。搜尋區域較大可能納入不同構造帶；時間拉長可能跨越測網更替。前面讀過的目錄完整度與空間異質性，在這裡直接影響判讀。若演算法把觀測能力改善當成活動增加，後面的迴歸再精確也無法修補。
#
# 因此，閱讀演算法結果時先注意它的輸入、搜尋範圍及接受條件，再看找到哪些現象。矩形邊界與數值門檻的細節可查本章末的原論文。接著換一種量法，檢查結論是否依賴規模累積這個選擇。
#
# %% [markdown]
# ## 10.4 換一種量法：只看事件數
#
# $C(t)$ 同時使用時間與規模。我們也可以先問較單純的問題：後半段的事件率，是否高於前半段？圓形搜尋以震源附近的圓域代替矩形，搭配標準化的率差 $Z$ 比較活動強弱。
#
# 把兩段的事件率差除以估計的標準誤，便能讓事件數不同的時間窗有較接近的比較尺度。下面沿用前面的合成事件，標出兩段時間長度及事件數。率差相同時，樣本量也會影響標準化結果：事件越少，對「增加」的估計越不穩。
#
# $Z$ 的算式不直接使用規模，但輸入事件仍經過規模門檻篩選。改動規模會改變入選的事件，進而改變 $Z$。這是資料分析中常見的情形：統計量看似沒有使用某變數，前處理卻已經用到了它。
#
# %% tags=["remove-input"]
dT1, dT2 = onset - T_S, T_M - onset
N1, N2 = int(pri_mask.sum()), int(pre_mask.sum())
Z_demo = (N2 * dT1 - N1 * dT2) / np.sqrt(N2 * dT1 ** 2 + N1 * dT2 ** 2)

fig = go.Figure()
fig.add_vrect(x0=T_S, x1=onset, fillcolor=PALETTE[4], opacity=0.13,
              line_width=0, layer="below")
fig.add_vrect(x0=onset, x1=T_M, fillcolor=ACCENT, opacity=0.13,
              line_width=0, layer="below")
for tt in t_all[pri_mask]:
    fig.add_shape(type="line", x0=tt, x1=tt, y0=0, y1=1,
                  line=dict(color=PALETTE[4], width=1.6))
for tt in t_all[pre_mask]:
    fig.add_shape(type="line", x0=tt, x1=tt, y0=0, y1=1,
                  line=dict(color=ACCENT, width=1.8))
fig.add_trace(go.Scatter(x=[T_M], y=[1.15], mode="markers", name="主震",
                         marker=dict(size=15, color=QUAKE_COLOR, symbol="star")))
fig.add_annotation(x=(T_S + onset) / 2, y=1.55, showarrow=False,
                   text=f"先前期：ΔT₁ = {dT1:.1f} 年，N₁ = {N1}<br>"
                        f"速率 = {N1 / dT1:.2f} /年")
fig.add_annotation(x=(onset + T_M) / 2, y=1.55, showarrow=False,
                   text=f"前兆期：ΔT₂ = {dT2:.1f} 年，N₂ = {N2}<br>"
                        f"速率 = {N2 / dT2:.2f} /年")
fig.add_vline(x=onset, line_dash="dot", line_color="#1baf7a")
fig.add_annotation(x=onset, y=-0.28, text="onset", showarrow=False)
fig.update_yaxes(range=[-0.45, 1.95], showticklabels=False, title_text="")
apply_layout(fig, height=340, hovermode="closest",
             xaxis_title="時間（年）", showlegend=False,
             title=(f"Z 值的兩時窗定義：Z = {Z_demo:.2f}"
                    f"（≈ {Z_demo:.1f} 個標準差；公式完全不含規模）"))
fig

# %% [markdown]
# 這張圖比較的是標準化程度，不是大地震發生機率。就算 $Z$ 很大，也只表示相對於這個率差模型，兩段活動差異較明顯。若搜尋過很多圓心、半徑和時間窗，還須讓對照實驗重複同一輪搜尋，才能判斷所選最大值是否特別。
#
# %% [markdown]
# ## 10.5 對照實驗究竟排除了什麼
#
# 前面的量測找出一種形態，現在才開始追問來源。隨機化的作用是保留部分資料特徵、破壞另一部分，再看結果如何變化；每一種隨機化都對應不同問題。
#
# 若打亂規模與事件的配對，保留時間及位置，就能檢查規模排列是否重要。若打亂時間，則改變時間叢集。若先去除餘震再分析，則檢查結果是否高度依賴所用除叢規則。它們不是同一個虛無假設，也不能合併成「已排除所有隨機性」。
#
# 2024 年論文對物理模擬目錄安排這些比較。以下合成示意用相同想法呈現不同處理後的分數；它不重現論文的實測效應量。注意我們比較的是完整分析程序的結果，而不只是原始目錄外觀。
#
# %% tags=["remove-input"]
RNG3 = np.random.default_rng(15003)
MC_SCAN = 3.5                      # 掃描用下限規模（高於合成目錄完整度 3.0）
M_MIN, BETA = 3.0, np.log(10)      # 合成目錄的完整度與 GR 斜率（b = 1）
T_TOTAL, L_BOX, N_MS = 600.0, 400.0, 40      # 年、km、主震數

bg_t = RNG3.uniform(0, T_TOTAL, 6000)
bg_x, bg_y = RNG3.uniform(0, L_BOX, 6000), RNG3.uniform(0, L_BOX, 6000)
bg_m = M_MIN + RNG3.exponential(1 / BETA, 6000)

zx, zy = RNG3.uniform(70, L_BOX - 70, 5), RNG3.uniform(70, L_BOX - 70, 5)
zi = RNG3.integers(0, 5, N_MS)
ms_x, ms_y = zx[zi] + RNG3.normal(0, 25, N_MS), zy[zi] + RNG3.normal(0, 25, N_MS)
ms_t = np.sort(RNG3.uniform(150, T_TOTAL - 5, N_MS))

pt, px, py, pm, at, ax_, ay_, am = [], [], [], [], [], [], [], []
for k in range(N_MS):
    Tp = 10 ** RNG3.normal(1.05, 0.22)                 # 前兆時間（年）
    n_p = int(RNG3.integers(18, 34))
    pt.append(ms_t[k] - RNG3.uniform(0, Tp, n_p))
    px.append(ms_x[k] + RNG3.normal(0, 18, n_p))
    py.append(ms_y[k] + RNG3.normal(0, 18, n_p))
    pm.append(M_MIN + 0.55 + RNG3.exponential(1 / BETA, n_p))
    n_a = int(RNG3.integers(50, 90))
    at.append(ms_t[k] + 0.05 * ((1 - RNG3.random(n_a)) ** (-2.5) - 1))
    ax_.append(ms_x[k] + RNG3.normal(0, 12, n_a))
    ay_.append(ms_y[k] + RNG3.normal(0, 12, n_a))
    am.append(M_MIN + RNG3.exponential(1 / BETA, n_a))

cat_t = np.concatenate([bg_t] + pt + at)
cat_x = np.concatenate([bg_x] + px + ax_)
cat_y = np.concatenate([bg_y] + py + ay_)
cat_m = np.concatenate([bg_m] + pm + am)
is_aft = np.concatenate([np.zeros(bg_t.size, bool)]
                        + [np.zeros(a.size, bool) for a in pt]
                        + [np.ones(a.size, bool) for a in at])
keep_t = cat_t < T_TOTAL
cat_t, cat_x, cat_y, cat_m, is_aft = (a[keep_t] for a in
                                      (cat_t, cat_x, cat_y, cat_m, is_aft))

R_GRID = np.geomspace(6.0, 60.0, 10)
T_GRID = np.geomspace(2.0, 80.0, 10)


def scan_max_z(t, m, x, y):
    """對每個主震在 (R, T) 網格上掃描，回傳最大 Z 值（無合格窗則 nan）。"""
    out = np.full(N_MS, np.nan)
    keep = m >= MC_SCAN
    tk, mk, xk, yk = t[keep], m[keep], x[keep], y[keep]
    for k in range(N_MS):
        d = np.hypot(xk - ms_x[k], yk - ms_y[k])
        sel = (d <= R_GRID[-1]) & (tk < ms_t[k]) & (tk >= ms_t[k] - T_GRID[-1])
        if sel.sum() < 10:
            continue
        o = np.argsort(tk[sel])
        ts_, ms_m, ds_ = tk[sel][o], mk[sel][o], d[sel][o]
        best = -np.inf
        for R in R_GRID:
            in_r = ds_ <= R
            tr, mr = ts_[in_r], ms_m[in_r]
            if tr.size < 10:
                continue
            for T in T_GRID:
                t0 = ms_t[k] - T
                j = int(np.searchsorted(tr, t0))
                tw, mw = tr[j:], mr[j:]
                n = tw.size
                if n < 10:
                    continue
                exc = mw - MC_SCAN + 0.1
                before = np.concatenate(([0.0], np.cumsum(exc)[:-1]))
                C = before - (exc.sum() / T) * (tw - t0)
                i = int(np.argmin(C))
                n1, n2 = i, n - i
                if n1 < 3 or n2 < 3:
                    continue
                d1, d2 = tw[i] - t0, ms_t[k] - tw[i]
                if d1 <= 0 or d2 <= 0:
                    continue
                z = (n2 * d1 - n1 * d2) / np.sqrt(n2 * d1 ** 2 + n1 * d2 ** 2)
                best = max(best, z)
        if np.isfinite(best):
            out[k] = best
    return out


perm_m = RNG3.permutation(cat_m.size)
perm_t = RNG3.permutation(cat_t.size)
cases = {
    "原始": (cat_t, cat_m),
    "規模隨機化": (cat_t, cat_m[perm_m]),
    "時間隨機化": (cat_t[perm_t], cat_m),
}
colors = {"原始": QUAKE_COLOR, "規模隨機化": PALETTE[2], "時間隨機化": ACCENT}

fig = go.Figure()
for name, (tt, mm) in cases.items():
    for dec, dash, width in ((False, "solid", 2.4), (True, "dot", 1.4)):
        msk = ~is_aft if dec else np.ones(cat_t.size, bool)
        z = scan_max_z(tt[msk], mm[msk], cat_x[msk], cat_y[msk])
        z = np.sort(z[np.isfinite(z)])
        fig.add_trace(go.Scatter(
            x=z, y=np.arange(1, z.size + 1) / z.size, mode="lines",
            name=f"{name}{'（去餘震）' if dec else ''}",
            line=dict(color=colors[name], width=width, dash=dash)))
apply_layout(fig, height=440, hovermode="x",
             xaxis_title="每個主震掃描到的最大 Z 值",
             yaxis_title="累積比例",
             title=(f"三組隨機化對照（合成示意，非論文重製）："
                    f"{N_MS} 個主震、{cat_t.size} 個事件"))
fig

# %% [markdown]
# 若破壞時間排列後分數降低，合理的解釋是此分數依賴原來的時間結構。要進一步說這是某個物理前兆，還要比較能產生叢集的替代模型；要說它有預報價值，則必須在未參與選法的未來資料上評分。這幾個問題需要不同證據。
#
# %% [markdown]
# ## 10.6 一場主震可能對應多組尺度
#
# 保留全部合格時間窗後，同一主震可能有好幾組 $T_P$ 與 $A_P$。有的區域小、時間長；有的區域大、時間短。這可能是資料與判準還不足以決定唯一尺度，不一定表示演算法出錯。
#
# 想像在地圖上逐步放大搜尋圈：較遠的事件加入後，累積曲線的平均線與起點都可能改變。因此不能把某一次辨識得到的「前兆十年」當成一個不帶條件的地震屬性。較完整的結果應保留可接受範圍，並說明搜尋規則。
#
# 下面七個合成辨識點示範同一主震內的尺度差異，另放入 Rhoades、Rastin 與 Christophersen（2022）的[已發表對照](https://doi.org/10.3390/geosciences12090349)：2020 年愛琴海 M6.7 事件在原文 Figure 1 與 Figure 12 的兩組辨識，分別為 10,220 天、3,203 km²，以及 6,392 天、8,091 km²。
#
# 愛琴海兩點不是前面的合成目錄，也沒有與七個合成點合併估計斜率；圖中分別顯示合成擬合與真實兩點的連線。兩種資料都用來說明同一主震可以有不同辨識尺度，不據此估計跨主震的共同關係。
#
# %% tags=["remove-input"]
RNG4 = np.random.default_rng(15004)
logT_s = 3.5 + np.linspace(-0.45, 0.45, 7) + RNG4.normal(0, 0.03, 7)
logA_s = 3.9 - 1.0 * (logT_s - 3.5) + RNG4.normal(0, 0.07, 7)
b_syn = np.polyfit(logT_s, logA_s, 1)[0]

# 真實對照（愛琴海 M6.7，兩組合法識別）：T_P（天）、A_P（km²）
aeg_T = np.array([10220.0, 6392.0])
aeg_A = np.array([3203.0, 8091.0])
b_aeg = np.log10(aeg_A[1] / aeg_A[0]) / np.log10(aeg_T[1] / aeg_T[0])

xline = np.array([3.0, 4.1])
yline = (logA_s.mean() + logT_s.mean()) - xline          # 斜率 −1 的參考線

fig = go.Figure()
fig.add_trace(go.Scatter(x=10 ** xline, y=10 ** yline, mode="lines",
                         name="斜率 −1（等量抵換）",
                         line=dict(color="#666666", width=2, dash="dash")))
fig.add_trace(go.Scatter(x=10 ** logT_s, y=10 ** logA_s, mode="lines+markers",
                         name="合成：同一主震的 7 個識別",
                         line=dict(color=ACCENT, width=1.2, dash="dot"),
                         marker=dict(size=10, color=ACCENT)))
fig.add_trace(go.Scatter(x=aeg_T, y=aeg_A, mode="lines+markers",
                         name="真實：愛琴海 M6.7 的兩組識別",
                         line=dict(color=QUAKE_COLOR, width=1.6, dash="dot"),
                         marker=dict(size=13, color=QUAKE_COLOR, symbol="diamond")))
fig.update_xaxes(type="log", title_text="前兆時間 T_P（天）")
fig.update_yaxes(type="log", title_text="前兆面積 A_P（km²）")
apply_layout(fig, height=460, hovermode="closest",
             title=(f"多重識別沿抵換線分布：合成擬合斜率 {b_syn:.2f}、"
                    f"愛琴海兩點連線斜率 {b_aeg:.2f}"))
fig

# %% [markdown]
# 一個主震若提供很多辨識點，這些點共享同一份歷史和同一次主震，不能當成同樣多個獨立樣本。忽略這層相依性，會高估尺度關係的精確度，也會讓容易找到很多組的主震獲得過大權重。
#
# %% [markdown]
# ## 10.7 時間與空間如何互相補償
#
# 把同一主震的 $\log_{10}T_P$ 與 $\log_{10}A_P$ 畫在一起，常會看見向下的排列：擴大空間範圍，可以伴隨較短的時間範圍。這叫時空取捨，描述的是辨識方式與資料共同形成的關係。
#
# 在統計上，我們可以給每場主震自己的截距，再估計共享的組內斜率。這相當於先移除不同主震的平均水準，才問同一主震內的兩個尺度是否一起變動。它和把所有點直接混在一起回歸，回答的問題不同。
#
# 另一個細節是兩個尺度都有測量誤差。普通最小平方法把橫軸當成沒有誤差；當橫軸也有誤差時，斜率往往被拉向零。把橫、縱軸交換再倒回來，通常得不到原來那條線。以下用合成資料比較兩種回歸；{doc}`附錄 D <appendix_d_eepas>`說明差異從何而來。
#
# 這些現象提醒我們：兩條不同斜率不一定代表兩種物理機制。先弄清楚估計方法的假設，才有條件解讀尺度關係。
#
# %% tags=["remove-input"]
RNG5 = np.random.default_rng(15005)
B_TRUE = -1.0
SD_EPS, SD_DEL = 0.887, 1.425          # 兩個方向的量測誤差（見附錄 B）
n_pt = 140
Xs = RNG5.normal(0, 1.0, n_pt)                       # 真值 X*
Ys = B_TRUE * Xs
Xo = Xs + RNG5.normal(0, SD_EPS, n_pt)               # 觀測 log T_P
Yo = Ys + RNG5.normal(0, SD_DEL, n_pt)               # 觀測 log A_P
b_fwd = np.polyfit(Xo, Yo, 1)[0]
b_rev_inv = 1.0 / np.polyfit(Yo, Xo, 1)[0]

s_grid = np.linspace(0.0, 1.3, 40)
b_fwd_curve = B_TRUE / (1 + s_grid ** 2)
b_rev_curve = B_TRUE * (1 + s_grid ** 2)

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                    subplot_titles=("一次實現：兩條擬合線夾住真值",
                                    "衰減因子：誤差愈大，夾得愈開"))
xg = np.linspace(Xo.min(), Xo.max(), 2)
fig.add_trace(go.Scatter(x=Xo, y=Yo, mode="markers", name="觀測（雙向皆帶誤差）",
                         marker=dict(size=6, color=ACCENT, opacity=0.55)),
              row=1, col=1)
for slope, nm, col, dsh in ((B_TRUE, "真值 −1", "#666666", "dash"),
                            (b_fwd, f"正向 {b_fwd:.2f}", PALETTE[2], "solid"),
                            (b_rev_inv, f"反向換算 {b_rev_inv:.2f}",
                             QUAKE_COLOR, "solid")):
    fig.add_trace(go.Scatter(x=xg, y=slope * (xg - Xo.mean()) + Yo.mean(),
                             mode="lines", name=nm,
                             line=dict(color=col, width=2, dash=dsh)),
                  row=1, col=1)
fig.add_trace(go.Scatter(x=s_grid, y=b_fwd_curve, mode="lines",
                         name="正向（偏淺）", showlegend=False,
                         line=dict(color=PALETTE[2], width=2.5)), row=1, col=2)
fig.add_trace(go.Scatter(x=s_grid, y=b_rev_curve, mode="lines",
                         name="反向換算（偏陡）", showlegend=False,
                         line=dict(color=QUAKE_COLOR, width=2.5)), row=1, col=2)
fig.add_hline(y=B_TRUE, line_dash="dash", line_color="#666666", row=1, col=2)
for yv, txt in ((-0.56, "真實資料 −0.56"), (-3.03, "真實資料 −3.0")):
    fig.add_hline(y=yv, line_dash="dot", line_color="#999999", row=1, col=2)
    fig.add_annotation(x=0.05, y=yv, text=txt, showarrow=False, xanchor="left",
                       yshift=11, font=dict(size=11), row=1, col=2)
fig.update_xaxes(title_text="觀測 log₁₀ T_P（任意原點）", row=1, col=1)
fig.update_yaxes(title_text="觀測 log₁₀ A_P", row=1, col=1)
fig.update_xaxes(title_text="量測誤差標準差（真值標準差 = 1）", row=1, col=2)
fig.update_yaxes(title_text="估計斜率", range=[-4.2, 0.2], row=1, col=2)
apply_layout(fig, height=430, hovermode="closest",
             title="回歸稀釋：兩個方向各偏一邊，真值被夾在中間")
fig

# %% [markdown]
# 當對數斜率接近 $-1$，$\log A_P+\log T_P$ 近似固定，也就是 $A_PT_P$ 較穩定。因此可以檢查面積與時間的乘積是否較穩定，但不能單憑它斷言存在守恆定律。
#
# Rastin 等人（2021）還從 EEPAS 參數擬合觀察到類似取捨：把時間尺度改變後，空間尺度能補回部分預報表現。這提供另一條模型層次的線索；有限的目錄長度與搜尋邊界，也可能造成相似的補償。
#
# %% [markdown]
# ## 10.8 先分清楚組內與組間
#
# 同一主震可以用較大面積換較短時間，但不同主震之間，較大的主震又可能同時對應較大的面積與較長時間。於是組內是負相關，組間卻是正相關。
#
# 以下把每場主震用不同顏色表示。先沿同色點看，再比較各組中心，最後才看全部資料的擬合線。三條閱讀路徑回答三個不同問題；合在一起的線可能很平，卻不表示前兩種關係都不存在。
#
# 這是 Simpson 型聚合效應的教學例子。類似問題也會出現在不同構造區的 $b$ 值、不同年代的事件率，以及不同地震序列的預報分數中。跨區域或跨序列合併之前，先問哪些點共享同一背景。
#
# %% tags=["remove-input"]
RNG6 = np.random.default_rng(15006)
N_M6, offs = 8, np.array([-0.35, -0.175, 0.0, 0.175, 0.35])
logT_mean = np.linspace(3.0, 4.0, N_M6)
logA_mean = 0.30 + 0.65 * logT_mean + RNG6.normal(0, 0.05, N_M6)

fig = go.Figure()
allT, allA = [], []
for i in range(N_M6):
    tt = logT_mean[i] + offs + RNG6.normal(0, 0.02, offs.size)
    aa = logA_mean[i] - 1.0 * offs + RNG6.normal(0, 0.05, offs.size)
    allT.append(tt)
    allA.append(aa)
    fig.add_trace(go.Scatter(x=tt, y=aa, mode="lines+markers", showlegend=False,
                             line=dict(color=PALETTE[i % len(PALETTE)],
                                       width=1, dash="dot"),
                             marker=dict(size=7,
                                         color=PALETTE[i % len(PALETTE)])))
allT, allA = np.concatenate(allT), np.concatenate(allA)
b_pool, a_pool = np.polyfit(allT, allA, 1)
b_mean, a_mean = np.polyfit(logT_mean, logA_mean, 1)
xg = np.linspace(2.5, 4.5, 2)
fig.add_trace(go.Scatter(x=xg, y=a_pool + b_pool * xg, mode="lines",
                         name=f"全部識別擬合：斜率 {b_pool:+.2f}",
                         line=dict(color="#444444", width=3)))
fig.add_trace(go.Scatter(x=xg, y=a_mean + b_mean * xg, mode="lines",
                         name=f"每主震取平均：斜率 {b_mean:+.2f}",
                         line=dict(color=QUAKE_COLOR, width=3, dash="dash")))
fig.add_trace(go.Scatter(x=logT_mean, y=logA_mean, mode="markers",
                         name="各主震的平均",
                         marker=dict(size=13, color=QUAKE_COLOR, symbol="x")))
apply_layout(fig, height=470, hovermode="closest",
             xaxis_title="log₁₀ T_P（天）", yaxis_title="log₁₀ A_P（km²）",
             title=(f"Simpson 陷阱（合成示意）：主震內負相關、跨主震正相關，"
                    f"混在一起得到 {b_pool:+.2f}"))
fig

# %% [markdown]
# 圖中的組內線與組中心連線方向不同。若研究目標是跨主震尺度，應以主震為分析單位或使用階層模型；若研究目標是辨識不確定性，就保留主震內的變異。不能只選相關係數較高的那種整理方式。
#
# %% [markdown]
# ## 10.9 從尺度關係走向預報
#
# 經過前面的檢查，尺度關係較適合被看成帶散佈的經驗回歸，而不是一場地震對另一場地震的確定預告。EEPAS 沿用的核心形式是
#
# $$\begin{aligned}
# M_m &= a_M+b_M M_P+\epsilon_M,\\
# \log_{10}T_P &= a_T+b_T M_P+\epsilon_T,\\
# \log_{10}A_P &= a_A+b_A M_P+\epsilon_A.
# \end{aligned}$$ (eq:psi-scaling)
#
# $M_m$ 是主震規模，$M_P$ 是選定前兆群的規模摘要；$T_P$ 以天、$A_P$ 以平方公裡表達。改變單位會改變截距，不能只抄係數而省略單位。殘差 $\epsilon$ 保留每組資料偏離回歸線的程度。
#
# 一條回歸線提供中心趨勢，殘差提供不確定性。要成為預報，兩者缺一不可：只給中心值會把相當分散的關係畫成過於精確的時間和位置；只給很寬的範圍又可能不比長期背景模型更有資訊。
#
# 自動辨識研究支持繼續檢查這組尺度形式，同時顯示結果受選窗、分群及目錄限制影響。下一章據此建立可以估計、可以比較的機率核，不將某組回溯係數當成普適常數。
#
# %% [markdown]
# 其他震前研究也可以用同一套問題檢查：事件是如何選入的？分析窗是否在看過主震後調整？對照資料是否經過相同搜尋？哪些結果屬於單一主震內部，哪些才是跨主震關係？
#
# 這些問題不會直接給出「有」或「沒有」前兆的答案，卻能界定一張圖支持到哪裡。尤其不能把統計辨識出的起點，當成物理孕震過程的確切起始時間。
#
# %% [markdown]
# 物理模擬與不同構造環境的比較，可以進一步追問尺度從何而來。模擬中若出現類似形態，表示那套規則足以產生形態；還需改變輸入及對照模型，才能找出哪些機制不可或缺。模擬結果與真實地球的推論仍須分開。
#
# 對預報而言，更直接的考驗是能不能把現在已有的事件轉成未來機率，且不借用未來主震的位置與時間。{doc}`EEPAS 與 PPE <16_eepas_ppe>`將從這個限制出發：不先辨識哪一群事件是真正前兆，而讓每個已發生事件依尺度提供有限、可檢驗的貢獻。
#
# %% [markdown]
# ```{admonition} 延伸推導：量測與回歸的限制
# :class: dropdown
#
# {doc}`附錄 D <appendix_d_eepas>`推導 $C(t)$ 的端點性質、最小值與兩段速率的關係、測量誤差造成的回歸稀釋，以及組內與組間共變異數的分解。這些推導用來限定前面的解讀，不是額外的前兆判準。
# ```
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 2 節的 Ψ 現象與尺度關係，再讀限制與未解問題，掌握本章經驗觀察如何連到機率模型。
# - [Algorithmic Identification of the Precursory Scale Increase Phenomenon in Earthquake Catalogs](https://doi.org/10.1785/0220240233) — Annemarie Christophersen、David A. Rhoades、Sebastian Hainzl，2024，*Seismological Research Letters*；[免費機構典藏全文](https://gfzpublic.gfz.de/rest/items/item_5029405_4/component/file_5029659/content)。這是本章自動辨識與對照實驗的已發表來源，建議比較矩形、圓形搜尋與隨機化目錄的設計；辨識到統計現象仍須與前瞻預報能力分開判斷。
# - [Space–Time Trade-Off of Precursory Seismicity in New Zealand and California Revealed by a Medium-Term Earthquake Forecasting Model](https://doi.org/10.3390/app112110215) — Sepideh J. Rastin、David A. Rhoades、Annemarie Christophersen，2021，*Applied Sciences*；免費開放全文。研究時間與空間參數的取捨，適合延伸本章「同一事件的 Ψ 辨識不唯一」及參數解讀問題。
