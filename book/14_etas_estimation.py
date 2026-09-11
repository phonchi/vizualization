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
# # 9. ETAS II：從目錄估計模型，再檢查它
#
# 前一章從參數出發生成目錄，看見事件率如何疊加、家族如何延伸。實際研究走的是反方向：我們只有時間、位置和規模，要從這些觀測估計背景、觸發與衰減，並判斷哪些結論有資料支持。
#
# 參數多只是其中一個困難。真實目錄沒有親子關係，研究區外的事件可能影響區內，大震後也可能漏掉小震。即使找到一組擬合很好的參數，還要問：其他參數是不是也一樣好？模型的不足會不會被某個參數吸收？
#
# 我們先把概似接到明確的觀測窗，接著利用模擬理解參數互相補償，再用機率權重處理未知親代。最後用診斷檢查模型及這些估計形成的預報分布。參數表本身並不足以回答這個問題。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 9.1 概似要對準我們真正觀察的範圍
#
# 前面點過程概似的兩部分，在 ETAS 中仍有相同意義：已發生事件的位置應有足夠的率，而整個觀測窗中模型所配置的活動量也必須合理。第一項對窗內事件求和，第二項對同一個窗積分，不能使用不同範圍。
#
# 設研究時間為 $[S,T]$、區域為 $\mathcal R$。較早或區域外已知的事件，可以當作條件歷史影響窗內的率，但不列入窗內事件的概似求和。這就像觀察一場已經進行中的序列：沒有從主震那刻開始收集資料，不代表主震對後面的影響消失了。
#
# 邊界的另一面是：窗內事件的部分後代可能發生在終點以後，或落在區域之外。正規化核的總面積是 1，但概似積分只能計入觀測窗的部分。對每個觸發源，時間佔比 $G_i$ 與空間佔比 $F_i$ 都介於 0 和 1，其貢獻為 $\kappa(m_i)G_iF_i$。把它們全部近似成 1，會改變估計問題。
#
# 本書的簡單可分離版本把新事件規模密度 $s(m)$ 乘在時空率外，因此在參數不共用、門檻固定的設定下，對數概似可分成規模與時空兩部分。這讓計算較方便，但不表示所有估計誤差在科學解讀上互不影響：分支比同時使用 $b$ 和產能參數，門檻與漏測也會影響兩部分。
#
# 完整的窗內積分及較早事件的時間修正見 {doc}`appendix_c_etas`。主文先記住：研究窗限制的是評分範圍，不應任意截斷可用歷史。
#
# %% [markdown]
# ## 9.2 最佳化解決計算，不能補足資料資訊
#
# ETAS 的時空參數通常沒有簡單的解析估計解，必須反覆調整參數，使對數概似增加。梯度與擬牛頓方法利用局部曲率加快搜尋；不同起始值有助於發現局部最優或平坦方向。程式回報成功，只表示達到其停止準則，不等於找到唯一或全域最佳答案。
#
# 計算量也來自事件歷史。若每個事件都與之前所有事件配對，工作量接近事件數的平方。提高輸入門檻可以減少事件，卻也刪去潛在觸發源；截斷很久以前或很遠的貢獻可加速，但應檢查這種近似是否改變預報。
#
# 誤差估計常從最優點附近的曲率開始。概似若近似二次、估計在參數空間內部且樣本足夠，觀測資訊矩陣的反矩陣可提供近似共變異數。平坦、彎曲或邊界附近的概似面不適合只報一組對稱標準誤。下一張圖讓我們直接看見這些條件。
#
# %% [markdown]
# ## 9.3 為什麼兩組參數都能描述同一份資料？
#
# 產能是 $A e^{\alpha(m-m_0)}$。若大部分事件的規模接近，提高 $A$、降低 $\alpha$，就可能在資料實際涵蓋的規模範圍產生近似的產能。資料能看見的是這個組合，未必能精確拆出每個部分。
#
# 同樣地，空間尺度是 $D e^{\gamma(m-m_0)}$。若親代規模幾乎相同，資料主要限制一個乘積，$D$ 和 $\gamma$ 就能互相補償。更大的規模跨度可能有幫助，更多樣本也有幫助；哪種增加更有效，仍取決於實際資訊結構。
#
# 先看一份時間 ETAS 模擬的 $A$–$\alpha$ 概似面。圖中固定其他參數，只改這兩項；因此它展示的是條件切面，不是已把所有未知量都考慮進去的聯合信賴區域。長而斜的等高線表示一組補償方向，短而集中才表示兩個參數都受約束。
#
# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.stats import kstest

from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, SEQUENTIAL, apply_layout

M0, BETA = 3.0, np.log(10)          # 門檻規模、GR 斜率（b = 1）
C_TRUE, P_TRUE = 0.01, 1.15         # Omori 核（天）
MU_TRUE, A_TRUE, AL_TRUE = 0.30, 0.20, 1.60
T_OBS = 1200.0                      # 觀測窗（天）


def simulate_etas(mu, A, alpha, T, seed, c=C_TRUE, p=P_TRUE):
    """時間型 ETAS 的分支法模擬（第 13 章的演算法）。"""
    rng = np.random.default_rng(seed)
    n_bg = rng.poisson(mu * T)
    t = list(rng.uniform(0, T, n_bg))
    m = list(M0 + rng.exponential(1 / BETA, n_bg))
    todo = list(range(len(t)))
    while todo:
        i = todo.pop()
        for _ in range(rng.poisson(A * np.exp(alpha * (m[i] - M0)))):
            dt = c * ((1 - rng.random()) ** (-1 / (p - 1)) - 1)
            if t[i] + dt < T:
                t.append(t[i] + dt)
                m.append(M0 + rng.exponential(1 / BETA))
                todo.append(len(t) - 1)
    order = np.argsort(t)
    return np.array(t)[order], np.array(m)[order]


def omori_matrices(t, T, c=C_TRUE, p=P_TRUE):
    """回傳 (G, G_int)：G[j,i] = g(t_j - t_i)，G_int[i] = ∫ g(t - t_i) dt。"""
    dt = t[:, None] - t[None, :]
    G = np.where(dt > 0, (p - 1) / c * (1 + np.abs(dt) / c) ** (-p), 0.0)
    return G, 1 - (1 + (T - t) / c) ** (1 - p)


t_syn, m_syn = simulate_etas(MU_TRUE, A_TRUE, AL_TRUE, T_OBS, seed=20240403)
G_syn, Gint_syn = omori_matrices(t_syn, T_OBS)
dm_syn = m_syn - M0

al_grid = np.linspace(1.00, 2.15, 46)
A_grid = np.linspace(0.04, 0.62, 46)
LL = np.empty((len(al_grid), len(A_grid)))
for ia, al in enumerate(al_grid):
    w = np.exp(al * dm_syn)
    S, W = G_syn @ w, (np.exp(al * dm_syn) * Gint_syn).sum()
    for ib, A in enumerate(A_grid):
        LL[ia, ib] = np.log(MU_TRUE + A * S).sum() - MU_TRUE * T_OBS - A * W

ia, ib = np.unravel_index(LL.argmax(), LL.shape)
# 山脊：固定「全目錄的總觸發量」A · Σ w_i G_i
tot = A_TRUE * (np.exp(AL_TRUE * dm_syn) * Gint_syn).sum()
A_ridge = tot / np.array([(np.exp(a * dm_syn) * Gint_syn).sum() for a in al_grid])

fig = go.Figure(go.Contour(
    x=A_grid, y=al_grid, z=LL - LL.max(), colorscale=SEQUENTIAL, reversescale=True,
    contours=dict(start=-30, end=0, size=1.5), colorbar=dict(title="Δ ln L")))
fig.add_trace(go.Scatter(x=A_ridge, y=al_grid, mode="lines",
                         name="總觸發量固定的方向",
                         line=dict(color="#444444", width=2, dash="dot")))
fig.add_trace(go.Scatter(x=[A_TRUE], y=[AL_TRUE], mode="markers", name="真值",
                         marker=dict(color=QUAKE_COLOR, size=13, symbol="x")))
fig.add_trace(go.Scatter(x=[A_grid[ib]], y=[al_grid[ia]], mode="markers",
                         name=f"網格 MLE（{A_grid[ib]:.3f}, {al_grid[ia]:.2f}）",
                         marker=dict(color=PALETTE[3], size=11,
                                     symbol="circle-open", line=dict(width=3))))
apply_layout(fig, title=f"A 與 α 的對數概似面：一條斜的山脊"
                        f"（合成目錄 N = {len(t_syn)}）",
             xaxis_title="A（基準產能）", yaxis_title="α（產能指數）",
             hovermode="closest", height=470)
fig

# %% [markdown]
# 沿著等高線較長的方向，改變兩個參數可以只造成很小的概似差異。若只拿最大點的 $\widehat\alpha$ 與另一研究比較，就會漏掉這個聯合不確定性。
#
# 再看空間核。下圖直接模擬已知親子距離，比較窄規模跨度與較寬跨度的情境；較寬跨度情境也有較多樣本，因此這個對照不能把改善完全歸因於跨度一項。左側看補償方向，右側看對 $D$ 最佳化後的 $\gamma$ 剖面概似。
#
# %% tags=["remove-input"]
Q_SP, D_SP, GAM_SP = 1.8, 2.0e-3, 1.0       # 空間核真值（D 單位 deg²）


def spatial_profile(n_pair, dm_max, seed=1988, n_grid=121):
    """回傳 (γ 網格, log10 D 網格, 對數概似面, 平均 Δm)。"""
    rng = np.random.default_rng(seed)
    dm = rng.exponential(1 / BETA, n_pair * 8)
    dm = dm[dm < dm_max][:n_pair]
    sig = D_SP * np.exp(GAM_SP * dm)
    r = np.sqrt(sig * ((1 - rng.random(len(dm))) ** (1 / (1 - Q_SP)) - 1))
    ga, lD = np.linspace(-2.0, 4.5, n_grid), np.linspace(-4.5, -1.2, n_grid)
    ll = np.empty((n_grid, n_grid))
    for i, g in enumerate(ga):
        s = 10 ** lD[None, :] * np.exp(g * dm)[:, None]
        ll[i] = np.sum(np.log((Q_SP - 1) / (np.pi * s))
                       - Q_SP * np.log1p((r ** 2)[:, None] / s), axis=0)
    return ga, lD, ll, dm.mean()


ga_n, lD_n, LL_n, dmbar = spatial_profile(250, 0.8)      # 窄規模跨度
ga_w, lD_w, LL_w, _ = spatial_profile(1500, 2.5)         # 寬規模跨度
iw, jw = np.unravel_index(LL_n.argmax(), LL_n.shape)
prof_n, prof_w = LL_n.max(axis=1), LL_w.max(axis=1)
ci_n, ci_w = ga_n[prof_n > prof_n.max() - 1.92], ga_w[prof_w > prof_w.max() - 1.92]
lD_ridge = lD_n[jw] + (ga_n[iw] - ga_n) * dmbar / np.log(10)   # D·e^{γΔm̄} 固定

fig = make_subplots(rows=1, cols=2, column_widths=[0.56, 0.44],
                    subplot_titles=("窄規模跨度：概似面是一條香蕉",
                                    "γ 的剖面概似：跨度決定寬度"))
fig.add_trace(go.Contour(x=lD_n, y=ga_n, z=LL_n - LL_n.max(), colorscale=SEQUENTIAL,
                         reversescale=True, showscale=False,
                         contours=dict(start=-25, end=0, size=1.5)), row=1, col=1)
fig.add_trace(go.Scatter(x=lD_ridge, y=ga_n, mode="lines", name="D·e^(γΔm̄) 固定",
                         line=dict(color="#444444", width=2, dash="dot")),
              row=1, col=1)
fig.add_trace(go.Scatter(x=[np.log10(D_SP)], y=[GAM_SP], mode="markers", name="真值",
                         marker=dict(color=QUAKE_COLOR, size=13, symbol="x")),
              row=1, col=1)
for ga_, pf_, lab, col in [
        (ga_n, prof_n, f"Δm ≤ 0.8（95% CI {ci_n.min():.2f}–{ci_n.max():.2f}）",
         PALETTE[1]),
        (ga_w, prof_w, f"Δm ≤ 2.5（95% CI {ci_w.min():.2f}–{ci_w.max():.2f}）",
         ACCENT)]:
    fig.add_trace(go.Scatter(x=ga_, y=pf_ - pf_.max(), mode="lines", name=lab,
                             line=dict(color=col, width=2)), row=1, col=2)
fig.add_hline(y=-1.92, line_dash="dash", line_color="#888888", row=1, col=2)
fig.update_xaxes(title_text="log₁₀ D（deg²）", row=1, col=1)
fig.update_yaxes(title_text="γ", row=1, col=1)
fig.update_xaxes(title_text="γ", row=1, col=2)
fig.update_yaxes(title_text="剖面 Δ ln L", range=[-12, 1], row=1, col=2)
apply_layout(fig, title="γ 與 D 只透過乘積 D·e^(γΔm) 被資料看見",
             hovermode="closest", height=440)
fig

# %% [markdown]
# 較平的剖面意味著多種 $\gamma$ 都能相容於資料，而不是說空間尺度真的隨時間大幅改變。圖中的 95% 區間採用大樣本概似比近似，碰到搜尋網格邊緣時也可能不完整。
#
# 這個例子還比真實目錄容易：模擬已知親子配對，實際上我們連這些配對都不知道。下一節把「不知道誰觸發誰」轉成模型中的潛在變數，接回先前學過的 EM 想法。
#
# %% [markdown]
# ## 9.4 用機率表達未知的親代
#
# 對某個事件 $j$，模型的總率由背景和各個先前事件的觸發貢獻相加。若背景在該時空點佔總率三成，那麼在指定模型及參數下，事件來自背景成分的機率就是三成；每個候選親代也按自己的率貢獻分配機率。
#
# 將對新事件規模積分後的時空強度記為 $\lambda^*_{ST}(t,x,y)$，便有
#
# $$\rho_{ij}=\frac{\kappa(m_i)g(t_j-t_i)f(x_j-x_i,y_j-y_i;m_i)}{\lambda^*_{ST}(t_j,x_j,y_j)},\qquad
# \phi_j=\frac{\mu(x_j,y_j)}{\lambda^*_{ST}(t_j,x_j,y_j)}.$$ (eq:rho-phi)
#
# $\rho_{ij}$ 是模型內事件 $j$ 的親代為 $i$ 的機率，$\phi_j$ 是背景機率，兩者滿足 $\phi_j+\sum_{i:t_i<t_j}\rho_{ij}=1$。共同的規模密度已消去，分母不能誤用尚未對規模積分的四維強度。
#
# Zhuang、Ogata 與 Vere-Jones（2002）的隨機除叢，用這些權重連結背景場估計與分支模型。背景率可由事件位置加權平滑，每個事件只以 $\phi_j$ 的權重貢獻；更新背景後再更新觸發參數與權重。平滑頻寬也影響哪些空間變化被分配給背景，須一併交代。
#
# Veen 與 Schoenberg（2008）將未知親代視為不完整資料，發展 EM 型估計。E 步先在目前參數下計算親代機率；M 步再依這些權重最大化期望完整資料概似。這相當於先用機率分配候選親代，再重新估計各成分，暫不指定每個事件的唯一來源。精確 EM 的概似不下降性依賴正確的 E/M 步；加入任意平滑、近似與截斷後，不能自動沿用這個保證。
#
# 這些機率是模型生成結構的推論，**不是物理因果被量測出來的機率**。漏測的小事件、區域外觸發或共同外部驅動，都可能改變分配。只在固定參數下反覆抽取親代，可呈現標籤不確定性，但不包含所有參數與模型不確定性。
#
# 下面用花蓮春季目錄展示背景權重的地圖。觸發參數採自公開 112 年委辦報告的既有示範設定，固定它們後只交替平滑背景及更新權重；這不是重新完成一套 ETAS 最大概似或完整 EM 估計。圖用短時間窗、經緯度平面距離，也未完整校正邊界，只適合學習權重的含義。
#
# %% tags=["remove-input"]
# CWA 112 年報告（詹忠翰等 2023）子計畫三的觸發參數
CWA = dict(A=0.6188, alpha=1.1733, c=0.0031, p=1.0616,
           D=5e-5, gamma=0.6786, q=1.5934)

cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
cat = (cat[(cat.depth <= 35) & (cat.ML >= 3.0)]
       .sort_values("time").reset_index(drop=True))
tt = (cat.time - cat.time.min()).dt.total_seconds().values / 86400.0
xx, yy, mm = cat.longitude.values, cat.latitude.values, cat.ML.values
N, T_WIN = len(tt), tt.max()

kap = CWA["A"] * np.exp(CWA["alpha"] * (mm - M0))
sig = CWA["D"] * np.exp(CWA["gamma"] * (mm - M0))
trig, CHUNK = np.zeros(N), 400
for s0 in range(0, N, CHUNK):                       # 觸發項只算一次
    s1 = min(s0 + CHUNK, N)
    d_t = tt[s0:s1, None] - tt[None, :]
    g = np.where(d_t > 0, (CWA["p"] - 1) / CWA["c"]
                 * (1 + np.abs(d_t) / CWA["c"]) ** (-CWA["p"]), 0.0)
    r2 = (xx[s0:s1, None] - xx[None, :]) ** 2 + (yy[s0:s1, None] - yy[None, :]) ** 2
    f = (CWA["q"] - 1) / (np.pi * sig) * (1 + r2 / sig) ** (-CWA["q"])
    trig[s0:s1] = (g * f * kap).sum(axis=1)

H_MIN, N_P = 0.03, 5                                 # ≈3 km，配 5 個最近鄰
dists, _ = cKDTree(np.c_[xx, yy]).query(np.c_[xx, yy], k=N_P + 1)
h = np.maximum(H_MIN, dists[:, N_P])
Kmat = np.empty((N, N), dtype=np.float32)            # Z_{h_k}(z_j - z_k)
for s0 in range(0, N, CHUNK):
    s1 = min(s0 + CHUNK, N)
    r2 = (xx[s0:s1, None] - xx[None, :]) ** 2 + (yy[s0:s1, None] - yy[None, :]) ** 2
    Kmat[s0:s1] = (np.exp(-r2 / (2 * h ** 2))
                   / (2 * np.pi * h ** 2)).astype(np.float32)

phi = np.full(N, 0.5)
for _ in range(15):                                  # 交替更新背景場與 φ
    u_hat = (Kmat @ phi.astype(np.float32)) / T_WIN
    phi = (u_hat / (u_hat + trig)).astype(float)
n_bg, main_i = phi.sum(), int(np.argmax(mm))

fig = go.Figure(go.Scattergl(
    x=xx, y=yy, mode="markers",
    marker=dict(size=3 + (mm - M0) * 2.6, color=phi, colorscale=SEQUENTIAL,
                cmin=0, cmax=1, showscale=True, line=dict(width=0),
                colorbar=dict(title="φ（背景機率）")),
    text=[f"ML {v:.1f}，φ = {p_:.3f}" for v, p_ in zip(mm, phi)],
    hoverinfo="text", name="事件"))
fig.add_trace(go.Scatter(x=[xx[main_i]], y=[yy[main_i]], mode="markers",
                         name=f"0403 主震 ML {mm[main_i]:.2f}（φ = {phi[main_i]:.2f}）",
                         marker=dict(color=QUAKE_COLOR, size=16, symbol="star")))
apply_layout(fig, title=f"隨機除叢：2024 春季目錄的背景機率 φ"
                        f"（N = {N}，Σφ = {n_bg:.0f}，叢集係數 "
                        f"{1 - n_bg / N:.2f}）",
             xaxis_title="經度（°E）", yaxis_title="緯度（°N）",
             hovermode="closest", height=560)
fig.update_yaxes(scaleanchor="x", scaleratio=1.0)
fig

# %% [markdown]
# 顏色表示這組設定下的背景機率，不能直接把淺色事件稱為已確認的某次主震後代。圖例中的 $\sum_j\phi_j$ 是固定模型下的期望背景標籤數；$1-\sum_j\phi_j/N$ 是這段目錄的平均觸發權重，不必等於長期分支比。
#
# 這裡的短窗可能讓背景和觸發互相吸收，距離近似及外部參數也會影響數值。這些限制沒有使圖失去教學用途：它清楚呈現「同一事件可能分到多種來源」的想法；只是正式預報仍需獨立估計與驗證。
#
# %% [markdown]
# ## 9.5 把模型當成時鐘，檢查它漏看了什麼
#
# 估計後不能只看概似值。前面介紹的時間變換，把原本快慢不一的條件率積分成一個新時鐘。如果模型正確，轉換後的事件間隔應與單位率 Poisson 過程相容。
#
# 一段物理時間內模型預期很多事件，新時鐘走得快；預期很少事件，則走得慢。將變換時間 $\tau_j$ 對事件序號 $j$ 作圖，正確模型下應大致沿斜率 1 的直線波動。偏離直線表示觀測活動與模型配置不同，但方向要配合座標讀：若事件序號增加得慢、變換時間卻增加很多，曲線會較陡，表示模型配置的活動比觀測多。
#
# 下圖用同一份合成目錄對照真參數 ETAS 與常數率 Poisson。後者將平均活動量攤在整段期間，可能匹配總數，卻無法同時描述背景空檔與叢集。這是「一階平均正確，但時間依賴結構不夠」的例子。
#
# 真正的診斷也要檢查變換間隔的分布、依賴性、空間殘差及規模分布。參數由同一資料估得時，檢定統計量須考慮估計效應，通常可透過模擬及重新擬合校準。
#
# %% tags=["remove-input"]
def transformed_time(t, m, mu, A, alpha, c=C_TRUE, p=P_TRUE):
    """τ_j = ∫₀^{t_j} λ*(u) du（時間型 ETAS）。"""
    d_t = t[:, None] - t[None, :]
    G_cum = np.where(d_t > 0, 1 - (1 + np.abs(d_t) / c) ** (1 - p), 0.0)
    return mu * t + A * (G_cum @ np.exp(alpha * (m - M0)))


tau_ok = transformed_time(t_syn, m_syn, MU_TRUE, A_TRUE, AL_TRUE)
tau_po = len(t_syn) / T_OBS * t_syn          # 錯設：忽略叢集的齊次 Poisson
jj = np.arange(1, len(t_syn) + 1)

fig = make_subplots(rows=1, cols=2, shared_yaxes=True,
                    subplot_titles=("正確模型：ETAS 真參數",
                                    "錯設模型：齊次 Poisson"))
for col, (tau, lab, color) in enumerate(
        [(tau_ok, "ETAS", ACCENT), (tau_po, "Poisson", PALETTE[1])], start=1):
    ks = kstest(1 - np.exp(-np.diff(np.r_[0, tau])), "uniform")
    fig.add_trace(go.Scatter(x=jj, y=tau, mode="lines",
                             name=f"{lab}（KS p = {ks.pvalue:.2g}）",
                             line=dict(color=color, width=2)), row=1, col=col)
    fig.add_trace(go.Scatter(x=jj, y=jj, mode="lines", showlegend=(col == 1),
                             name="y = x", line=dict(color="#888888", dash="dash")),
                  row=1, col=col)
    fig.update_xaxes(title_text="事件序號 j", row=1, col=col)
fig.update_yaxes(title_text="變換時間 τ", row=1, col=1)
apply_layout(fig, title="τ–j 診斷圖：對的模型貼直線，錯的模型走樓梯",
             hovermode="x", height=440)
fig

# %% [markdown]
# 圖中的 KS 數字只提供這個示範的分布對照，不等於模型已通過所有檢驗。左圖使用模擬真參數，右圖的平均率則由同一資料估計；兩者的參數條件並不相同，也沒有另檢查所有間隔的獨立性。
#
# 更重要的是，所謂「相對寧靜」必須相對某個模型定義。Poisson 認為不尋常的空檔，ETAS 可能視為叢集之間正常的間隔。若要把殘差解讀成新的前兆訊息，先要排除參考模型與資料條件不足的可能。
#
# %% [markdown]
# ## 9.6 看不見的事件會影響整段歷史
#
# 在規模統計章，漏測主要改變規模分布；在 ETAS 中，被漏掉的事件還是潛在觸發源。因此資料缺漏同時影響觀測數與歷史 $H_t$，不能只把少掉的事件數當成一般隨機缺值。
#
# 大事件後最初的小震若特別容易漏記，觀測資料可能使大事件的產能看起來較低，降低估計的 $\alpha$；早期率被壓低也可能改變 $c,p$。但多個參數、背景與邊界可以互相補償，偏差方向及大小仍依情境而定。
#
# Naylor et al.（2023）使用合成資料研究 ETAS 反演，明確呈現率相關的不完整性對估計的影響。貝氏方法可以表達參數的不確定性，卻不能在沒有合適觀測模型時自動消除漏測偏差。
#
# 下圖先生成已知參數的完整目錄，再刪掉大事件後短時間內的部分小事件，使用同一估計程式重新擬合。其他條件固定，才能把變化與這個特定漏測規則連起來。
#
# %% tags=["remove-input"]
def fit_mu_A_alpha(t, m, T, c=C_TRUE, p=P_TRUE):
    """固定 c、p，用 Nelder–Mead 估 (μ, A, α)。"""
    G, Gint = omori_matrices(t, T, c, p)
    dm = m - M0

    def nll(par):
        mu, A, al = np.exp(par[0]), np.exp(par[1]), par[2]
        if not 0.0 < al < 4.0:
            return 1e12
        w = np.exp(al * dm)
        lam = mu + A * (G @ w)
        if not np.all(lam > 0):
            return 1e12
        return -(np.log(lam).sum() - mu * T - A * (w * Gint).sum())

    res = minimize(nll, [np.log(0.3), np.log(0.2), 1.2], method="Nelder-Mead",
                   options=dict(maxiter=3000, xatol=1e-4, fatol=1e-3))
    return np.exp(res.x[0]), np.exp(res.x[1]), res.x[2]


def blind_after_big(t, m, cut, big=5.0, window=1.0):
    """模仿主震後的目錄失能：M≥big 事件後 window 天內刪掉 m < M0+cut 的事件。"""
    keep = np.ones(len(t), bool)
    for i in np.where(m >= big)[0]:
        keep &= ~((t > t[i]) & (t <= t[i] + window) & (m < M0 + cut))
    return t[keep], m[keep]


cuts, seeds = np.array([0.0, 0.5, 1.0, 1.5, 2.0]), (7, 20240403, 99)
alpha_hat = np.empty((len(seeds), len(cuts)))
A_hat, lost = np.empty_like(alpha_hat), np.empty_like(alpha_hat)
for si, sd in enumerate(seeds):
    t_c, m_c = simulate_etas(MU_TRUE, A_TRUE, AL_TRUE, T_OBS, seed=sd)
    for ci, cu in enumerate(cuts):
        t_b, m_b = blind_after_big(t_c, m_c, cu)
        _, A_hat[si, ci], alpha_hat[si, ci] = fit_mu_A_alpha(t_b, m_b, T_OBS)
        lost[si, ci] = 100 * (1 - len(t_b) / len(t_c))

fig = make_subplots(specs=[[{"secondary_y": True}]])
fig.add_trace(go.Scatter(x=cuts, y=alpha_hat.mean(axis=0), mode="lines+markers",
                         name="α̂（三個種子平均）",
                         error_y=dict(type="data", array=alpha_hat.std(axis=0)),
                         line=dict(color=ACCENT, width=2.5)), secondary_y=False)
fig.add_trace(go.Scatter(x=cuts, y=A_hat.mean(axis=0), mode="lines+markers",
                         name="Â（補償性上升）",
                         line=dict(color=PALETTE[1], width=2, dash="dot")),
              secondary_y=False)
fig.add_trace(go.Bar(x=cuts, y=lost.mean(axis=0), name="被刪掉的事件比例",
                     marker=dict(color="#cccccc"), opacity=0.55), secondary_y=True)
fig.add_hline(y=AL_TRUE, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"真值 α = {AL_TRUE}")
fig.update_yaxes(title_text="估計值", range=[0, 1.9], secondary_y=False)
fig.update_yaxes(title_text="刪掉的事件（%）", range=[0, 40], secondary_y=True)
apply_layout(fig, title=f"只刪掉 {lost[:, 2].mean():.0f}% 的早期小事件，"
                        f"α̂ 就從 {alpha_hat[:, 0].mean():.2f} 掉到 "
                        f"{alpha_hat[:, 2].mean():.2f}",
             xaxis_title="M≥5 事件後一天內被刪掉的規模範圍（m₀ 到 m₀+cut）",
             hovermode="x", height=460, barmode="overlay")
fig

# %% [markdown]
# 比較完整與刪減情境，可以看出少量、但集中在關鍵時段的缺漏如何改變參數。誤差棒是三個固定種子之間的散佈，僅呈現這組教學模擬的重複差異，不是完整的參數信賴區間。
#
# 實際資料若估得 $\alpha$ 偏小，應先檢查完整度與模型設定，而不是直接宣佈區域的觸發效率不同。同樣地，有限窗 Omori 擬合得到 $p<1$ 也不自動表示錯誤；它只不能直接套用本書要求 $p>1$ 的無限時間正規化核。背景錯設是可能原因之一，需要比較與診斷才能判斷。
#
# %% [markdown]
# ## 9.7 用較簡單的模型理解預報機率
#
# 完整 ETAS 需要在未來事件發生後持續更新率。若暫時只保留一個指定主震的影響，不讓新事件再觸發，問題會簡單許多。Reasenberg–Jones 型餘震率用
#
# $$R(t,M)=10^{a'+b(M_m-M)}(t+c)^{-p}$$ (eq:rj-rate)
#
# 描述規模至少 $M$ 的餘震率；$M_m$ 是指定主震規模，$t$ 是主震後時間，$a'$ 控制整體產能。這個式子直接把「門檻提高，事件更少」與「時間經過，率下降」連起來。
#
# 它與 ETAS 的單一觸發項有一個受限的代數關係：只留下指定主震、忽略背景與後續事件的再觸發、採未截斷 GR，並令 $\alpha=\beta$，便能得到同樣的率形式。這不代表 R–J 是完整自激發過程的等價預報，實際擬合整段序列的產能也不能直接當成 ETAS 的直接後代產能。
#
# 若給定參數後把 $R(t,M)$ 視為外生的非齊次 Poisson 率，未來 $[t_1,t_2]$ 至少一次事件的機率是
#
# $$P=1-\exp\left[-\int_{t_1}^{t_2}R(t,M)\,dt\right].$$
#
# 這裡的 Poisson 假設是計算的一部分。完整 ETAS 要計入窗內新事件帶來的再次觸發，通常使用條件模擬或相應機率方程；不能把一條固定歷史的率積分冒充完整級聯的機率。
#
# 下圖固定一組文獻參數，改變指定主震規模，計算同一目標門檻的未來一週機率。它示範率與機率如何轉換，不是當前事件的作業預報。
#
# %% tags=["remove-input"]
RJ_CA = dict(a=-1.67, b=0.91, p=1.08, c=0.05)      # 加州通用（62 個序列）
RJ_JP = dict(a=-1.83, b=0.85, p=1.3, c=0.3)        # Utsu (1969) 日本


def rj_prob(t0, Mm, M, win=7.0, a=-1.67, b=0.91, p=1.08, c=0.05):
    """R–J：從 t0 起 win 天內至少發生一次 M 以上事件的機率。"""
    integ = (10 ** (a + b * (Mm - M))
             * ((t0 + c) ** (1 - p) - (t0 + win + c) ** (1 - p)) / (p - 1))
    return 1 - np.exp(-integ)


t0s = np.logspace(-2, 2.5, 120)
fig = go.Figure()
for Mm, color, lab in [(7.2, QUAKE_COLOR, "M 7.2（0403 花蓮量級）"),
                       (6.4, ACCENT, "M 6.4（2025 大埔量級）"),
                       (5.5, PALETTE[2], "M 5.5")]:
    fig.add_trace(go.Scatter(x=t0s, y=[100 * rj_prob(s, Mm, 5.0, **RJ_CA)
                                       for s in t0s],
                             mode="lines", name=lab,
                             line=dict(color=color, width=2.2)))
p_ca = 100 * rj_prob(0.0, 6.0, 6.0, **RJ_CA)
p_jp = 100 * rj_prob(0.0, 6.0, 6.0, **RJ_JP)
fig.add_hline(y=p_ca, line_dash="dash", line_color="#666666",
              annotation_text=f"「被同等或更大事件跟隨」加州 {p_ca:.1f}%"
                              f"／日本 {p_jp:.1f}%")
apply_layout(fig, title="R–J 加州通用參數：「未來 7 天內 M≥5」的機率如何衰減",
             xaxis_title="主震後時間（天）", yaxis_title="7 天內機率（%）",
             xaxis_type="log", yaxis_type="log", hovermode="x", height=460)
fig

# %% [markdown]
# 規模增加先使率乘上一個因子，再經過 $1-e^{-\Lambda}$ 轉成機率。因此當機率較小時曲線近似成比例；機率接近 1 時便飽和，不能說所有曲線都只是上下平移。
#
# 模型給出的機率還要接受資料檢驗：許多被報為某一機率範圍的情境中，實際事件比例是否相符？這是後面預報評估的問題。單次發生或沒發生，都不足以證明機率正確或錯誤。
#
# 直接後代與全世代期望數之間，在次臨界、完整家族條件下有 $1/(1-n)$ 的放大關係；但不同世代會改變時間形狀，不能因此把實際 R–J 的所有參數逐項換成 ETAS 參數。附錄 C 分開推導單一觸發項的對照與全世代期望。
#
# %% [markdown]
# ## 9.8 把未知參數帶進預報
#
# 到目前為止，多數圖使用一組指定或估計參數。若概似面很平，只採最大點可能低估預報的不確定性。貝氏方法將先驗資訊與目錄概似結合，得到參數的聯合後驗分布；預報時將這些參數情境一起納入平均，避免只依賴一個點估計。
#
# Naylor、Serafini、Lindgren 與 Main（2023）介紹 ETAS.inlabru，利用 INLA 與 inlabru 的近似推論工具研究時間 ETAS。讀這篇文章時，先看合成案例如何分別改變背景、觸發和偵測條件，再看聯合後驗如何表現參數補償。演算法較快是實作優點，教學上更重要的是它把「哪些參數仍未知」帶到預報中。
#
# 另一種選擇是減少自由參數。simplETAS 固定部分參數，只估計資料較能支持的部分，以提高可重現性和作業可行性。這是在偏差、變異和計算之間取捨；固定值不會因為是文獻設定就變成所有地區的物理常數。若使用 $\alpha=\beta$，也須清楚交代規模截斷及相應分支比。
#
# 參數更少或使用貝氏方法，都不能免除外部檢驗。應區分三種不確定性：同一參數下事件本來就隨機、資料尚未決定的參數，以及未納入的模型或觀測機制。只有明確放入推論的部分，才會反映在輸出的區間內。
#
# %% [markdown]
# 前一章從模型生成叢集，這一章則處理資料不完整與觀測範圍有限時的估計問題。現在看到 ETAS 的背景率、親代機率或殘差圖時，我們知道它們是哪些資料與假設的共同結果，也知道可以如何進一步檢查。
#
# 接下來的 {doc}`15_psi_phenomenon` 把問題延伸到更長時間尺度：較大事件之前的活動變化，是否提供超出一般叢集模型的資訊？這仍是可檢驗的統計問題，不應先把事後排列解讀為已知的預示機制。ETAS 在這裡提供比較基準，幫助分辨新增方法究竟增加了什麼。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01) — Abdollah Jalilian，2019，*Journal of Statistical Software*；[免費全文與程式](https://www.jstatsoft.org/article/view/v088c01)。沿著資料整理、模型擬合與結果診斷讀一遍，對照本章的背景率估計及隨機除叢；數學方法與 R 範例可以分開閱讀。
# - [ETAS: Epidemic-Type Aftershock Sequence](https://github.com/lmizrahi/etas) — Leila Mizrahi 等；免費作者程式庫與說明。README 將參數估計、模擬及變動完整度的實作連回各自論文，適合延伸本章「不完整資料如何影響參數」的問題。
# - [The Effect of Declustering on the Size Distribution of Mainshocks](https://doi.org/10.1785/0220200231) — Leila Mizrahi、Shyam Nandan、Stefan Wiemer，2021，*Seismological Research Letters*；出版社全文可能需訂閱，[免費作者預印本](https://arxiv.org/abs/2012.09053)。用合成目錄檢查除叢後的估計偏差，適合思考本章的背景／觸發分類結果能否直接當成物理事實。
#
# - [SimplETAS: A Benchmark Earthquake Forecasting Model Suitable for Operational Purposes and Seismic Hazard Analysis](https://doi.org/10.1785/0220230199) — Simone Mancini、Warner Marzocchi，2023 年線上發表，*Seismological Research Letters*；全文可能需訂閱，[作者程式庫](https://github.com/smancini2/simplETAS)可免費閱讀。直接對照本章參數精簡的理由，思考固定哪些參數能減少估計困難，以及簡化模型應如何接受樣本外檢驗。
#
# - [Stochastic Declustering of Space-Time Earthquake Occurrences](https://doi.org/10.1198/016214502760046925) — Jiancang Zhuang、Yosihiko Ogata、David Vere-Jones，2002。親代機率及背景加權估計的原始來源；出版社全文可能需訂閱。
# - [Estimation of Space–Time Branching Process Models in Seismology Using an EM–Type Algorithm](https://doi.org/10.1198/016214508000000148) — Alejandro Veen、Frederic P. Schoenberg，2008。從未知分支結構理解 E 步與 M 步，並核對方法採用的近似；出版社全文可能需訂閱。
# - [Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759) — Mark Naylor、Francesco Serafini、Finn Lindgren、Ian G. Main，2023。免費開放全文；以合成案例連結聯合後驗、參數補償與率相關漏測。
# - [Earthquake Hazard After a Mainshock in California](https://doi.org/10.1126/science.243.4895.1173) — Paul A. Reasenberg、Lucile M. Jones，1989。R–J 預報的原始來源，對照固定主震率與完整自激發預報的差異；出版社全文可能需訂閱。
# - [臺灣地區112年中大型地震震源資訊之快速彙整與提供](https://scweb.cwa.gov.tw/webdata/PDF/Reports/MOTC-CWB-112-E-01.pdf) — 中華民國地球物理學會，2023，中央氣象署委辦報告，子計畫三。免費官方報告；本章除叢示範的外部觸發參數來源，示範本身沒有重現其完整估計流程。
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, Second Edition](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones，2003。第 7 章連結概似、條件強度與時間變換；出版社全文可能需訂閱。
# - [GDMS 地球物理資料管理系統](https://gdms.cwa.gov.tw/) — 中央氣象署。背景權重圖的目錄來源；本章使用既有春季快取進行教學示範。
# - [Seismicity Analysis through Point-process Modeling: A Review](https://doi.org/10.1007/s000240050275) — Yosihiko Ogata，1999，*Pure and Applied Geophysics*，155:471–507。從條件強度一路連到概似、ETAS 與殘差，適合把本章放回統計地震學的發展脈絡；出版社全文可能需訂閱。
