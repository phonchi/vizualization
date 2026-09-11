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
# # 8. ETAS I：把每個事件接進預報
#
# 在 {doc}`12_clustering_laws`，單一 Omori 曲線遇到後續較大事件時，很難同時描述原來的衰減與新增活動。接著的估計章也提醒我們：一條曲線貼近資料，不代表其參數或假設已被資料確認。帶著這兩點，我們開始建立能描述重疊序列的模型。
#
# ETAS 的出發點是：每個已發生的事件，都可能增加後續一段時間的發生率。增加多少取決於事件規模，何時及何處出現後續活動則由衰減核描述。把所有事件的貢獻加起來，就不必事先指定哪一個是主震、哪些是餘震。
#
# 這一章沿著一份模擬目錄，先看率如何更新，再追蹤直接後代與整個家族的差別。模擬中的親子關係是我們生成資料時記錄的資訊；真實目錄沒有這項資訊。如何在不知道親代的情況下估計模型，留到下一章。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 8.1 一個事件發生後，模型更新了什麼？
#
# 先只考慮時間與規模。背景率描述模型尚未歸到先前事件觸發項的活動；每個事件再貢獻一條從自己的發震時刻開始的衰減曲線。因此，對目前為止的歷史 $H_t$，時間型 ETAS 寫成
#
# $$\lambda^*(t)=\mu+\sum_{i:t_i<t}\kappa(m_i)g(t-t_i).$$
#
# 這裡的 $\mu$ 是時間模型的背景率；加入空間後會改寫成 $\mu(x,y)$。$\kappa(m_i)$ 控制事件 $i$ 的平均直接後代數，$g$ 則使用前面的正規化 Omori 密度。求和只使用 $t$ 以前的事件，因此模型是在利用已知歷史更新未來的發生率。
#
# 想像第一個事件出現時，背景上方多了一條衰減曲線。過了一段時間，另一個事件出現，便再疊一條新的曲線；即使它在日常用語中被叫作餘震，模型也允許它產生後代。這就是「自激發」：事件改變後續事件的條件率。
#
# 這個機制和背景率隨環境改變並不相同。若活動由共同外部因素升高，未必能把所有叢集都解釋為事件彼此觸發。ETAS 提供一個可估計、可檢驗的統計模型，而不是僅憑事件接近便確認因果。
#
# 下圖使用固定種子的模擬目錄。上方是根據當時歷史計算的條件強度，下方是實際抽到的事件。先找出較大的事件，再看它之後強度如何跳升、衰減，以及其間新事件帶來的再次更新。
#
# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

M0, M_MAX = 3.0, 8.0                 # 模型門檻規模、規模上限
BETA = np.log(10.0)                  # GR 斜率（b = 1）
C_OM, P_OM = 0.01, 1.3               # Omori 核參數（天）
DM = M_MAX - M0                      # 截斷寬度


def A_from_n(n_ratio, alpha):
    """由分支比反推產能尺度 A（截斷 GR 版本，推導見 13.3）。"""
    d = BETA - alpha
    inner = (1 - np.exp(-d * DM)) / d if abs(d) > 1e-9 else DM
    return n_ratio * (1 - np.exp(-BETA * DM)) / (BETA * inner)


def draw_mag(rng, size):
    """截斷 GR 抽樣（反函數法，見 10.6）：m 落在 [M0, M_MAX]。"""
    u = rng.random(size)
    return M0 - np.log1p(-u * (1 - np.exp(-BETA * DM))) / BETA


def simulate_etas(mu, n_ratio, alpha, T, seed):
    """時間型 ETAS 的分支法模擬（虛擬碼見 13.6）。

    回傳 (目錄 DataFrame, 產能尺度 A)。DataFrame 欄位：
    t 時間、m 規模、parent 親代索引（背景事件為 −1）、
    gen 世代（背景為 0）、root 所屬家族的始祖索引。
    """
    rng = np.random.default_rng(seed)
    A = A_from_n(n_ratio, alpha)
    n_bg = rng.poisson(mu * T)                       # 第 0 代：背景事件
    t = list(rng.uniform(0, T, n_bg))
    m = list(draw_mag(rng, n_bg))
    parent, gen, root = [-1] * n_bg, [0] * n_bg, list(range(n_bg))
    todo = list(range(n_bg))
    while todo:                                      # 遞迴繁殖
        i = todo.pop()
        k = rng.poisson(A * np.exp(alpha * (m[i] - M0)))
        if k == 0:
            continue
        dt = C_OM * ((1 - rng.random(k)) ** (-1 / (P_OM - 1)) - 1)
        mk = draw_mag(rng, k)
        for d, mv in zip(dt, mk):
            if t[i] + d < T:                         # 落在窗外的直接丟棄
                t.append(t[i] + d)
                m.append(mv)
                parent.append(i)
                gen.append(gen[i] + 1)
                root.append(root[i])
                todo.append(len(t) - 1)
    cat = pd.DataFrame(dict(t=t, m=m, parent=parent, gen=gen, root=root))
    order = np.argsort(cat.t.to_numpy())          # 依時間排序並重編索引
    newpos = np.empty(len(cat), dtype=int)
    newpos[order] = np.arange(len(cat))
    cat = cat.iloc[order].reset_index(drop=True)
    par = cat.parent.to_numpy()
    cat["parent"] = np.where(par >= 0, newpos[np.maximum(par, 0)], -1)
    cat["root"] = newpos[cat.root.to_numpy()]
    return cat, A


def lam_star(tt, cat, mu, A, alpha):
    """時間型 ETAS 的條件強度曲線（逐事件疊加）。"""
    lam = np.full_like(tt, float(mu))
    for ti, mi in zip(cat.t.to_numpy(), cat.m.to_numpy()):
        msk = tt > ti
        lam[msk] += (A * np.exp(alpha * (mi - M0)) * (P_OM - 1) / C_OM
                     * (1 + (tt[msk] - ti) / C_OM) ** -P_OM)
    return lam


MU1, N1, ALPHA1, T1 = 1.0, 0.85, 1.2, 180.0
cat1, A1 = simulate_etas(MU1, N1, ALPHA1, T1, seed=10)
tt = np.linspace(0, T1, 3000)
lam1 = lam_star(tt, cat1, MU1, A1, ALPHA1)

fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                    row_heights=[0.45, 0.55], vertical_spacing=0.04)
fig.add_trace(go.Scatter(x=tt, y=lam1, mode="lines", name="條件強度 λ*(t)",
                         line=dict(color=ACCENT, width=1.4)), row=1, col=1)
fig.add_trace(go.Scattergl(x=cat1.t, y=cat1.m, mode="markers", name="模擬事件",
                           marker=dict(size=4 + (cat1.m - M0) * 3.5,
                                       color=QUAKE_COLOR, opacity=0.6)),
              row=2, col=1)
fig.update_yaxes(title_text="λ*（次/天）", type="log", row=1, col=1)
fig.update_yaxes(title_text="規模", row=2, col=1)
fig.update_xaxes(title_text="時間（天）", row=2, col=1)
apply_layout(fig, height=540, hovermode="x",
             title=f"模擬的時間型 ETAS 目錄（n = {N1}，α = {ALPHA1}，"
                   f"共 {len(cat1)} 個事件，其中背景 "
                   f"{int((cat1.gen == 0).sum())} 個）")
fig

# %% [markdown]
# 高強度區間較容易出現密集事件，但一次實現仍有隨機性：高率不保證馬上發生，低率也不表示不可能發生。背景事件同樣是隨機生成，不是固定間隔的輸入。
#
# 圖中的模型不需預先切出序列，卻能形成集中活動。接下來把同一個疊加想法延伸到地圖上，並把「多少、何時、何處」分清楚。
#
# %% [markdown]
# ## 8.2 把時間、位置與規模放在同一個模型
#
# 時間模型告訴我們某段時間是否活躍，空間模型還要回答活動集中在哪裡。本書使用的時空 ETAS 條件強度為
#
# $$\lambda^*(t,x,y,m)=s(m)\left[\mu(x,y)+\sum_{i:t_i<t}\kappa(m_i)g(t-t_i)f(x-x_i,y-y_i;m_i)\right].$$ (eq:etas-intensity)
#
# 閱讀這個式子時，先從方括號看起。背景場 $\mu(x,y)$ 可以在不同位置有不同的率；每個過去事件都放上一個以自身時空位置為中心的觸發核。將所有貢獻加起來後，再乘上規模密度 $s(m)$，把總率分配給不同規模。
#
# 這裡沿用一個明確假設：新事件規模從同一個 GR 分布抽出，不直接取決於親代規模。過去事件的規模 $m_i$ 會改變未來事件的**數量與位置分布**，但本式沒有讓它改變新事件自身的規模密度。這也是一種可以接受檢驗、可以被擴充的假設。
#
# 三個觸發部分是
#
# $$\kappa(m)=A e^{\alpha(m-m_0)},\qquad
# f(x,y;m)=\frac{q-1}{\pi D e^{\gamma(m-m_0)}}\left[1+\frac{x^2+y^2}{D e^{\gamma(m-m_0)}}\right]^{-q},$$
#
# 以及 {eq}`eq:omori-density` 的 $g(t)$。在本書的無限時間、全平面版本中，$p>1$、$q>1$，兩個密度各積分為 1。因此，把單一事件的觸發率對時間、位置和規模積分，就只剩 $\kappa(m)$：它才可以解讀成平均直接後代數。
#
# 參數可以順著這個分工記憶：$A,\alpha$ 管產量；$c,p$ 管時間；$D,\gamma,q$ 管空間。$D e^{\gamma(m-m_0)}$ 的單位是距離平方，不能當成距離本身。Ogata 與 Zhuang（2006）的擴充將 $\alpha$ 與 $\gamma$ 分開，讓產量和範圍不必以同樣速度增長。
#
# ```{admonition} 為什麼文獻中的 K 不能直接當成 A？
# :class: dropdown
#
# 若時間觸發項寫成 $K e^{\alpha(m-m_0)}(t+c)^{-p}$，且其他核已正規化，則 $K=A(p-1)c^{p-1}$。換時間單位或換 $p,c$，$K$ 都會改變；$A$ 才是這個正規化約定下的直接後代期望數。若還有未正規化空間核，則必須再包含其積分因子。
# ```
#
# 正規化是為了讓參數容易解讀，並不要求實際觀測區把所有後代都收進來。跑出研究區或發生在資料終點以後的事件仍可能存在，估計時須處理這些邊界。下一章會利用這個差別寫出概似。
#
# %% [markdown]
# ## 8.3 平均一個事件會接出多少後代？
#
# $\kappa(m)$ 回答的是「已知親代規模時」的平均直接後代數。若從模型的規模分布中任選一個事件，再對規模平均，就得到分支比 $n$：
#
# $$n=\mathbb E[\kappa(M)]=\frac{A\beta}{\beta-\alpha},\qquad \alpha<\beta.$$ (eq:branching-ratio)
#
# 這個結果把兩種相反效應放在一起。$\alpha$ 較大，大事件的產能增加得快；$\beta$ 較大，大事件變得稀少得快。只有後者下降得足夠快，無上界 GR 模型的平均產能才有限。積分的逐步計算見 {doc}`appendix_c_etas`。
#
# $\alpha<\beta$ 是這個**無上界規模模型的有限均值條件**，不等於說 $\alpha\ge\beta$ 時任何有限目錄上的條件強度都無法寫出。有限的已觀測事件仍各有有限觸發項。過程是否局部有限、能否長期平穩，以及某個期望是否存在，是不同問題，不能只用「模型發散」混在一起。
#
# 若規模有指定上限 $M_{\max}$，先將 GR 密度在有限區間重新正規化，便可得到有限的平均產能，即使 $\alpha=\beta$ 也一樣。這時分支比會依賴 $M_{\max}$ 和 $m_0$，截斷並不保證 $n<1$。報告分支比時，規模範圍是定義的一部分。
#
# 這個數字也不是一個事件一定生幾個後代。$n=0.7$ 可以由許多沒有後代的事件，和少數有大量後代的事件共同形成；下一節會看平均產能如何沿世代累積。
#
# %% [markdown]
# ## 8.4 從直接後代走到整個家族
#
# 先從一個背景事件開始。它的平均直接後代數是 $n$，下一代又各自有平均 $n$ 個後代，因此第 $k$ 代的平均數為 $n^k$。把整個家族都算入、沒有截去太晚或太遠的成員，$n<1$ 時得到
#
# $$\mathbb E[Y]=1+n+n^2+\cdots=\frac{1}{1-n}.$$
#
# $Y$ 包含最初的背景事件。分支比從 0.5 增加到 0.9，完整家族的期望大小便從 2 增加到 10。兩個參數都小於 1，活動量卻能有很大的差別；這就是接近臨界時模型對 $n$ 特別敏感的原因。
#
# 對具有常數背景率、有限均值的平穩時間模型，長期總率為 $\bar\lambda=\mu/(1-n)$，觸發事件在平均率中的佔比為 $n$。這不是說任意短窗中都恰有 $n$ 的比例是觸發事件，也不是說觀測到的「餘震比例」可直接當成同一數字。
#
# 在一般非退化分支模型下，$n<1$ 的家族以機率 1 終止且平均大小有限；$n=1$ 時家族也可幾乎必然終止，但平均大小無限；$n>1$ 則有正機率持續繁衍。這是跨越所有世代的性質。**超臨界不代表有限時間窗的概似積分必然無限**：若已觀測歷史有限、背景可積且核局部可積，有限窗內的概似仍可能完全有限。受影響的是長期有限均值平穩解與相關解讀。
#
# 下面固定其他設定，比較不同分支比的模擬。圖中只觀察有限期間，並非從無限過去開始的平穩樣本；所以模擬的事件比例與完整家族理論不需要精確相等。
#
# %% tags=["remove-input"]
fig = go.Figure()
MU3, ALPHA3, T3, N_REP = 0.4, 1.2, 365.0, 25
amp = []
for n_r, color in [(0.5, PALETTE[2]), (0.9, PALETTE[1])]:
    tot = bgs = 0
    for r in range(N_REP):                    # 平均倍率取多份實現
        c, _ = simulate_etas(MU3, n_r, ALPHA3, T3, seed=800 + r)
        tot += len(c)
        bgs += int((c.gen == 0).sum())
    amp.append(tot / bgs)
    cat_n, _ = simulate_etas(MU3, n_r, ALPHA3, T3, seed=800)   # 畫其中一份
    daily, _ = np.histogram(cat_n.t, bins=np.arange(0, int(T3) + 1))
    fig.add_trace(go.Scatter(
        x=np.arange(int(T3)), y=daily, mode="lines",
        name=f"n = {n_r}（此實現共 {len(cat_n)} 個，背景 "
             f"{int((cat_n.gen == 0).sum())} 個）",
        line=dict(color=color, width=1.1)))
apply_layout(fig, height=440, hovermode="x",
             xaxis_title="時間（天）", yaxis_title="每日事件數",
             yaxis_type="log",
             title=f"接近臨界的世界：{N_REP} 份模擬的平均放大倍率 "
                   f"{amp[0]:.2f}（n = 0.5，理論 2.0）對 "
                   f"{amp[1]:.2f}（n = 0.9，理論 10.0）")
fig

# %% [markdown]
# 分支比較高時，模擬中較容易出現較大的活動波動和延續較久的家族。但單次模擬不必處處比另一份目錄更活躍；圖中多次模擬平均也仍有抽樣誤差。
#
# 有限窗還會截掉後代。尤其 Omori 核接近 $p=1$ 時，等待時間尾端很長，觀測終點之外可能仍有不少事件。這解釋了為何用窗內總數除以背景數，往往低於完整家族的 $1/(1-n)$；不能只憑這項差異判定模擬失敗。
#
# 現在將最初那份模擬依世代拆開。背景、直接後代及更晚世代分別堆疊後，可以看出一個活動尖峰通常由哪些世代共同組成。
#
# %% tags=["remove-input"]
edges = np.arange(0, T1 + 1, 1.0)
gcap = np.minimum(cat1.gen.to_numpy(), 3)
fig = go.Figure()
frac = []
for g, name, color in [(0, "第 0 代（背景）", PALETTE[0]),
                       (1, "第 1 代", PALETTE[2]),
                       (2, "第 2 代", PALETTE[3]),
                       (3, "第 3 代以後", PALETTE[1])]:
    h, _ = np.histogram(cat1.t.to_numpy()[gcap == g], bins=edges)
    frac.append(int((gcap == g).sum()))
    fig.add_trace(go.Scatter(x=edges[:-1], y=h, mode="lines", name=name,
                             stackgroup="one", line=dict(width=0.6,
                                                         color=color),
                             fillcolor=color))
apply_layout(fig, height=440, hovermode="x",
             xaxis_title="時間（天）", yaxis_title="每日事件數（堆疊）",
             title=f"世代分解（n = {N1}）：背景 {frac[0]}、第 1 代 {frac[1]}、"
                   f"第 2 代 {frac[2]}、第 3 代以後 {frac[3]} 個")
fig

# %% [markdown]
# 底層是隨機背景事件，其期望率固定，但實際每天計數仍會起伏。上方各層是不同世代，同一段高活動可能同時包含第一代與更晚的事件。因此，用一條主震衰減曲線描述全部序列，會把多層傳遞合在同一組參數中。
#
# 這張圖能畫出世代，是因為模擬時已知生成關係。面對真實地震，我們通常只知道時間、位置與規模，並不知道家族樹。後面隨機除叢的任務，就是在模型條件下用機率表示這個未知結構。
#
# %% [markdown]
# ## 8.5 最大的事件不一定最先發生
#
# ETAS 的後代規模不受親代規模限制，因此小事件之後可以出現大事件。我們事後可能把較早的事件叫作前震；在事件發生當下，模型並沒有一個已知的「前震類別」。
#
# 要看出這點，先只問直接後代。給定親代規模 $m_1$，每個後代超過它的機率是 $e^{-\beta(m_1-m_0)}$。直接後代數為 Poisson，依規模獨立篩選後，超過親代的數目仍是 Poisson。因此在未截斷規模、涵蓋所有時間與位置的這個版本中，
#
# $$P(\text{a larger direct offspring}\mid m_1)=1-\exp\{-A e^{-(\beta-\alpha)(m_1-m_0)}\}.$$
#
# 它不是「未來一週出現更大事件」的完整預報：式子沒有指定有限時間窗，也沒有算後續世代或其他背景家族。若要問這些問題，就要把相應事件一起納入，不能直接套用單一 Poisson 公式。
#
# 這個模型示範的是：有些事後被稱為前震的模式，可以在一般叢集機制下出現。它沒有證明所有真實前震都出於同一機制，也沒有保證能事先認出哪些事件將被更大的事件跟隨。
#
# 下面將最初的模擬加上親子連線。圈出的始祖是家族最早的背景事件，未必是家族中規模最大的那個。先沿時間方向追蹤幾條連線，再找最大的成員位於哪一代。
#
# %% tags=["remove-input"]
tv, mv, pv = cat1.t.to_numpy(), cat1.m.to_numpy(), cat1.parent.to_numpy()
big = int(np.argmax(mv))
founder = int(cat1.root[big])                 # 該家族的始祖（背景事件）
fam_mask = cat1.root.to_numpy() == founder

edge_x, edge_y = [], []
for k in np.flatnonzero(pv >= 0):
    edge_x += [tv[pv[k]], tv[k], None]
    edge_y += [mv[pv[k]], mv[k], None]

fig = go.Figure()
fig.add_trace(go.Scattergl(x=edge_x, y=edge_y, mode="lines", name="親 → 子",
                           line=dict(color="#bbbbbb", width=0.7)))
for sel, name, color in [(~fam_mask, "其他事件", ACCENT),
                         (fam_mask, "最大事件所屬家族", QUAKE_COLOR)]:
    fig.add_trace(go.Scattergl(
        x=tv[sel], y=mv[sel], mode="markers", name=name,
        marker=dict(size=4 + (mv[sel] - M0) * 3.5, color=color, opacity=0.75)))
fig.add_trace(go.Scatter(
    x=[tv[founder]], y=[mv[founder]], mode="markers", name="該家族的始祖",
    marker=dict(size=16, color="rgba(0,0,0,0)", symbol="circle",
                line=dict(color="#333333", width=2))))
apply_layout(fig, height=470, hovermode="closest",
             xaxis_title="時間（天）", yaxis_title="規模",
             title=f"觸發樹：最大事件 M {mv[big]:.2f}（第 {int(cat1.gen[big])} "
                   f"代）與它的家族（共 {int(fam_mask.sum())} 個事件，"
                   f"始祖 M {mv[founder]:.2f}）")
fig

# %% [markdown]
# 同一個家族中可以有多層分支。最大的事件若出現在後面，前面的較小事件便具有「先小後大」的外觀。真實資料上，這種外觀本身不足以辨認觸發鏈；模擬只是讓我們看見一種可能的生成方式。
#
# 下一張圖重複生成許多家族，分別計算全部家族，以及至少有十個成員的家族中，「最大事件不是始祖」的比例。另畫一條忽略規模與產能關聯的平均場參考線。不同曲線的條件不同，不能互相當成同一機率的估計。
#
# %% tags=["remove-input"]
def theory_not_founder(n_ratio, n_theta=400):
    """平均場理論：解 w = (1−θ)exp[n(w−1)]，再對 θ ~ U(0,1) 平均。

    截斷 GR 之下 θ = P(M > m1) 幾乎就是 U(0,1)，故直接在 θ 上積分。
    """
    th = (np.arange(n_theta) + 0.5) / n_theta
    w = np.zeros_like(th)
    for _ in range(200):                      # n < 1 保證是壓縮映射
        w = (1 - th) * np.exp(n_ratio * (w - 1))
    return float(np.mean(1 - w / (1 - th)))


def simulate_family(alpha, A, rng, cap=200_000):
    """只追蹤族譜與規模（時間、空間與觀測窗都不影響這個統計量）。

    回傳 (始祖規模, 家族總成員數, 家族最大規模)。
    """
    m_root = float(draw_mag(rng, 1)[0])
    queue, size, m_top = [m_root], 1, m_root
    while queue:
        mp = queue.pop()
        k = rng.poisson(A * np.exp(alpha * (mp - M0)))
        if k == 0:
            continue
        kids = draw_mag(rng, k)
        size += k
        m_top = max(m_top, float(kids.max()))
        queue.extend(kids.tolist())
        if size > cap:
            break
    return m_root, size, m_top


n_grid = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9])
N_FAM, ALPHA_F = 12000, 1.2
emp_all, emp_big, theo = [], [], []
for j, n_r in enumerate(n_grid):
    rng = np.random.default_rng(3000 + j)
    A_f = A_from_n(n_r, ALPHA_F)
    fam = np.array([simulate_family(ALPHA_F, A_f, rng) for _ in range(N_FAM)])
    hit = fam[:, 2] > fam[:, 0] + 1e-12
    big_fam = fam[:, 1] >= 10
    emp_all.append(float(hit.mean()))
    emp_big.append(float(hit[big_fam].mean())
                   if big_fam.sum() >= 60 else np.nan)
    theo.append(theory_not_founder(n_r))

fig = go.Figure()
fig.add_trace(go.Scatter(x=n_grid, y=theo, mode="lines",
                         name="平均場理論（全部家族）",
                         line=dict(color="#888888", width=2, dash="dash")))
fig.add_trace(go.Scatter(x=n_grid, y=emp_all, mode="lines+markers",
                         name="模擬：全部家族",
                         line=dict(color=ACCENT, width=2),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=n_grid, y=emp_big, mode="lines+markers",
                         name="模擬：成員 ≥ 10 的家族",
                         line=dict(color=PALETTE[1], width=2),
                         marker=dict(size=8, symbol="square")))
apply_layout(fig, height=450, hovermode="x",
             xaxis_title="分支比 n", yaxis_title="家族最大事件不是始祖的比例",
             yaxis_range=[0, 1],
             title=f"「最大的不是第一個」：n = 0.5 時 {emp_all[2]:.0%}，"
                   f"n = 0.9 時 {emp_all[-1]:.0%}（大家族分別為 "
                   f"{emp_big[2]:.0%}、{emp_big[-1]:.0%}）")
fig

# %% [markdown]
# 全部家族包含只有始祖的孤立事件，而這些家族的最大事件必然是始祖。只選較大家族，相當於事後改變樣本；比例不同並不意外，也不能把選樣後的比例直接當成即時預報。
#
# 灰色平均場曲線把每個事件的直接後代數都改成 $\operatorname{Poisson}(n)$，不再由其規模決定。真正的 ETAS 則會讓大規模與較高產能相連，所以兩者不必重合。平均場公式只是一個數學對照；它在 $n\to1$ 時不一定趨近 1。
#
# 這個例子還有一個模擬限制：生成器對極大家族設有計算上限，觸及上限的家族只能得到部分資訊。涉及尾端機率的正式研究，必須追蹤這種截斷的比例及敏感性。這裡的圖用於辨認家族選樣與平均場假設的差別。
#
# %% [markdown]
# ## 8.6 模擬讓我們分開看不同參數
#
# 分支法模擬先生成背景事件，再為每個事件生成直接後代，後代又依同一規則繼續產生下一代。時間從 $g$ 抽樣，位置從 $f$ 抽樣，規模從 $s$ 抽樣。這個順序對應模型的統計生成方式，並不表示地震內部真的有可直接觀測的親子標籤。
#
# 如果只關心固定觀測窗，模擬可保留窗內事件，但不能忽略可能從窗外或較早歷史進入的觸發。計算前也要宣告規模上限及模擬起始條件。分支法的完整步驟與等待時間反函數放在附錄 C，原始 notebook 保留可重現的程式。
#
# 回到參數的直覺：$n$ 只告訴我們平均產能，不告訴我們產能分配在誰身上。下圖固定分支比，改變 $\alpha$ 並相應調整 $A$，比較較小與較大事件各自承擔多少觸發量。這樣才能隔離「總平均相同，但貢獻分配不同」的效果。
#
# %% tags=["remove-input"]
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.09,
                    subplot_titles=("α = 0.5：群震型（大小事件觸發力接近）",
                                    "α = 2.0：主震–餘震型（大事件主導）"))
share = []
for row, (al, color) in enumerate([(0.5, PALETTE[2]), (2.0, PALETTE[1])], 1):
    cat_a, _ = simulate_etas(mu=1.0, n_ratio=0.7, alpha=al, T=200.0, seed=23)
    par = cat_a.parent.to_numpy()
    trig = par >= 0                                  # 被觸發的事件
    cut = np.quantile(cat_a.m.to_numpy(), 0.95)      # 最大的 5% 事件
    top_idx = set(np.flatnonzero(cat_a.m.to_numpy() >= cut).tolist())
    by_top = sum(1 for p in par[trig] if int(p) in top_idx)
    share.append(by_top / max(1, trig.sum()))
    fig.add_trace(go.Scattergl(
        x=cat_a.t, y=cat_a.m, mode="markers",
        name=f"α = {al}（{len(cat_a)} 個事件）",
        marker=dict(size=4 + (cat_a.m - M0) * 3.5, color=color, opacity=0.7)),
        row=row, col=1)
fig.update_yaxes(title_text="規模", range=[M0 - 0.2, 7.4])
fig.update_xaxes(title_text="時間（天）", row=2, col=1)
apply_layout(fig, height=560, hovermode="closest",
             title=f"同樣的分支比 n = 0.7，只改 α：最大的 5% 事件"
                   f"直接觸發了 {share[0]:.0%}（α = 0.5）對 "
                   f"{share[1]:.0%}（α = 2.0）的被觸發事件")
fig

# %% [markdown]
# 兩個情境的理論平均產能相同，實現出的家族與時間分布卻不同。$\alpha$ 較大時，產能隨規模增加較快，少數較大事件可以佔較多直接觸發量；$\alpha$ 較小時，貢獻分布相對分散。
#
# 圖中最大的若干事件佔了多少實現觸發數，仍是有限樣本的統計量。它幫助理解參數，不能單靠這種視覺形狀唯一決定 $A$ 與 $\alpha$。下一章的概似面將顯示，這兩個參數往往能互相補償。
#
# %% [markdown]
# ## 8.7 從分支直覺回到模型界線
#
# ETAS 名稱中的 epidemic，借用了分支傳遞的直覺：一個事件可能引出下一代。但地震沒有對應感染人數上限的簡單「易感人口」，模型也沒有直接寫入應力守恆、斷層幾何或能量平衡。類比可以幫助理解世代累積，不能取代物理機制。
#
# 背景與觸發的分解同樣依賴模型。如果真實活動由慢滑移、流體或其他共同因素驅動，固定背景版本可能將部分變化歸給觸發項；若未觀測到較小或區域外的事件，部分觸發也可能歸入背景。這屬於統計識別問題，不能據此替每個事件指定永久不變的物理來源。
#
# Reinhart（2018）的回顧將這個架構放回更廣的時空點過程：其他領域也有自激發模型，並共享背景設定、邊界、估計與殘差診斷等問題。點過程教材則提供分支表示和條件強度背後的機率結構，讓我們知道哪些結論需要平穩性、有限均值或獨立標記。
#
# %% [markdown]
# ## 8.8 下一步：資料能辨認多少結構？
#
# 時間 ETAS 適合先理解序列更新；時空 ETAS 加入位置；背景可再容許空間異質性或隨時間變化，觸發核也可加入方向性。這些版本並不是越複雜越好：每增加一種結構，都要問資料是否支持、估計是否穩定，以及在未參與擬合的資料上是否改善預報。
#
# 目前我們能用指定參數生成目錄，知道每個符號如何改變率與家族。真實分析的方向恰好相反：只給一份目錄，要估出參數，還要處理不知道誰觸發誰、看不見窗外後代及早期漏測的問題。
#
# {doc}`14_etas_estimation` 就從這個反向問題開始。前面建立的「率是各個來源相加」會再次發揮作用：它既能寫出概似，也能把未知親代轉成機率權重，串起最大概似、EM 與貝氏不確定性。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01) — Abdollah Jalilian，2019，*Journal of Statistical Software*；[免費全文與程式](https://www.jstatsoft.org/article/view/v088c01)。先讀模型結構與參數說明，把本章的背景率、時間核、空間核對照到一套完整的時空 ETAS 表示式。
# - [ETAS: Epidemic-Type Aftershock Sequence](https://github.com/lmizrahi/etas) — Leila Mizrahi 等；免費作者程式庫與說明。可從目錄模擬範例追讀 ETAS 如何生成事件序列，對照本章的世代分解與 branching 模擬；程式所用參數慣例需先與本章核對。
# - [Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes](https://doi.org/10.1080/01621459.1988.10478560) — Yosihiko Ogata，1988，*Journal of the American Statistical Association*；出版社全文可能需訂閱，[免費研究機構全文](https://bemlar.ism.ac.jp/zhuang/Refs/Refs/ogata1988.pdf)。回到 ETAS 的經典來源，閱讀事件歷史如何進入條件強度，以及作者如何用資料比較不同叢集模型。
#
# - [Space–time ETAS models and an improved extension](https://doi.org/10.1016/j.tecto.2005.10.016) — Yosihiko Ogata、Jiancang Zhuang，2006，*Tectonophysics*；出版社全文可能需訂閱，[免費作者機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。這是本章時空 ETAS 結構的重要來源，說明如何把 Omori 衰減與餘震區尺度關係放進條件強度，並比較模型延伸。
#
# - [A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629) — Alex Reinhart，2018，*Statistical Science*。把 ETAS 放回一般自激發點過程，串起分支、邊界、估計與診斷；[免費作者預印本](https://arxiv.org/abs/1708.02647)。
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, Second Edition](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones，2003。第 6 章群集與標記模型、第 7 章條件強度與概似，提供本章分支直覺的嚴謹背景；出版社全文可能需訂閱。
# - [Seismicity Analysis through Point-process Modeling: A Review](https://doi.org/10.1007/s000240050275) — Yosihiko Ogata，1999，*Pure and Applied Geophysics*，155:471–507。從條件強度一路連到概似、ETAS 與殘差，適合把本章放回統計地震學的發展脈絡；出版社全文可能需訂閱。
