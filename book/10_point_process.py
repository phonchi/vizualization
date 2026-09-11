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
# # 4. 點過程：把歷史轉成可更新的發生率
#
# 在{doc}`foundation_randomness`，我們比較了平均率相同、事件排列卻不同的
# 序列，也看見各種模型對「剛發生一個事件」有不同反應。接下來需要把這些反應
# 寫成可以計算的量，才能從描述目錄走到預報。
#
# 這個量是條件強度。名稱雖然抽象，問題其實很具體：已經知道現在以前的事件，
# 接下來很短一段時間，預期會有多少新事件？我們先從這句話建立直覺，再看它
# 如何連到等待機率、概似與診斷。證明及抽樣公式放在{doc}`appendix_a_point_process`，
# 正文保留各步驟為什麼需要，以及公式使用時不能省略的條件。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 4.1 率會隨已知資訊改變
#
# 把時間 $t$ 以前已發生的事件記為 $H_t$，把截至時間 $t$ 的事件總數記為
# $N(t)$。條件強度 $\lambda^*(t)$ 是給定這段歷史後，接下來極短期間的平均
# 事件增量除以期間長度：
#
# $$\lambda^*(t)=\lim_{\Delta\downarrow0}
# \frac{E[N(t+\Delta)-N(t)\mid H_t]}{\Delta}.$$ (eq:cond-int)
#
# 星號提醒我們這個率依賴歷史。時間以天計時，率的單位就是次／天，所以率
# 可以大於一；它本身不是機率。前面的累積計數是一條階梯，這裡也不是直接
# 對那條階梯微分，而是在已知歷史下描述可能出現的增量。
#
# 在只由事件歷史和已知時間函式決定的模型裡，沒有新事件時，可以沿著目前
# 歷史算出一條未來率曲線；新事件一出現，就必須更新它。若模型還含未觀測的
# 外部隨機狀態，則需要一併處理那部分資訊，不能假定任意條件強度都只有事件
# 到來時才更新。這裡先考慮前一種情況。
#
# ## 4.2 下一個事件還沒來的機率
#
# 假設現在是 $t_0$，問下一個事件是否晚於 $t$。在「到 $t$ 都沒有新事件」
# 這條路徑上，把率記為 $\lambda^*_0(u)$，則
#
# $$P(T_{\rm next}>t\mid H_{t_0})
# =\exp\!\left[-\int_{t_0}^{t}\lambda^*_0(u)\,\mathrm du\right].$$
#
# 它把「每一小段都沒有事件」累積成整段的存活機率。若率始終不變，就回到
# 上一章的指數等待時間；若事件後的率逐漸衰減，就要把整條曲線積起來，
# 不能只拿此刻很高的率乘上整段時間。
#
# 下圖設定背景加上一個已知事件的衰減貢獻，並比較下一次等待時間的理論與
# 模擬。另一條曲線故意把起點的率當成常數，讓我們看見少掉「率會變」這一步
# 會產生什麼差別。這是下一個任意事件的等待示例，不是某個目標規模的大震預報。
#
# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

MU_D, KAPPA_D = 0.20, 3.0          # 背景率（次/天）、上一個事件的產能
C_D, P_D = 0.30, 1.40              # Omori 核參數（天）

def lam_demo(t):
    """單一觸發源的條件強度。"""
    return MU_D + KAPPA_D * (P_D - 1) / C_D * (1 + t / C_D) ** (-P_D)

def Lam_demo(t):
    """對應的 compensator（解析積分）。"""
    return MU_D * t + KAPPA_D * (1 - (1 + t / C_D) ** (-(P_D - 1)))

rng = np.random.default_rng(20)
n_sim = 40000
E = rng.exponential(1.0, n_sim)                 # 單位速率 Poisson 的第一個間隔
t_grid = np.linspace(0, 300, 300001)            # 用單調的 Λ 做數值反函數
t_sim = np.interp(E, Lam_demo(t_grid), t_grid)  # 解 Λ(t) = E

t_plot = np.logspace(-3, np.log10(60), 400)
S_theory = np.exp(-Lam_demo(t_plot))
S_emp = np.array([(t_sim > t).mean() for t in t_plot])
S_naive = np.exp(-lam_demo(0.0) * t_plot)       # 誤把 t=0 的率當常數
gap = np.abs(S_emp - S_theory).max()

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_plot, y=S_theory, mode="lines",
                         name="理論 exp(−∫λ*)",
                         line=dict(color=ACCENT, width=3)))
fig.add_trace(go.Scatter(x=t_plot, y=S_emp, mode="markers",
                         name=f"模擬經驗分布（{n_sim} 次）",
                         marker=dict(color=QUAKE_COLOR, size=5, opacity=0.7)))
fig.add_trace(go.Scatter(x=t_plot, y=S_naive, mode="lines",
                         name="若誤當成常數率 λ*(0)",
                         line=dict(color=PALETTE[3], width=2, dash="dot")))
apply_layout(fig, title="時間變換抽樣與固定率近似",
             xaxis_title="距上一事件的時間（天）",
             yaxis_title="P（下一事件還沒來）",
             xaxis_type="log", yaxis_type="log",
             yaxis_range=[-2.2, 0.05], height=430, hovermode="x")
fig

# %% [markdown]
# 這裡以同一條累積率做反函數抽樣，所以抽樣曲線接近理論曲線，是這種生成方法
# 預期的結果，不能當成對存活公式的獨立驗證。真正要比較的是固定起點率的
# 另一條曲線：它預期事件更快到來，顯示忽略等待期間的衰減會如何改變答案。
#
# 若我們改問「未來是否有規模超過門檻的事件」，事情會多一層：等待較大事件
# 的期間，較小事件仍可能先發生，並提高之後的率。對自激發模型，要整合這些
# 可能歷史，不能只沿「完全沒有新事件」的路徑積分再乘規模尾機率。給定率的
# Poisson模型有簡單換算，ETAS的目標事件預報則往往需要模擬未來目錄。
#
# ## 4.3 讓模型同時解釋事件與空白
#
# 有了率曲線，就可以比較不同設定對觀測的支持程度。一個模型若在事件發生
# 的時刻給出很低的率，觀測會使它受到質疑；但把所有時刻的率都拉高也不能
# 解決問題，因為還有許多時間沒有事件。時間點過程的對數概似把兩者合在一起：
#
# $$\ln L=\sum_{i=1}^{N}\ln\lambda^*(t_i)
# -\int_0^T\lambda^*(u)\,\mathrm du.$$ (eq:pp-loglik)
#
# 第一項看事件發生的位置，第二項累積整個觀察窗的率。積分必須到觀察截止
# 時間 $T$，不能只到最後一次事件，因為最後那段空白也是資料。自激發模型
# 的積分沿實際觀測歷史計算，通常是隨資料變動的量，不是事先固定的平均次數。
#
# 這個概似形式建立在具條件強度、局部可積且不爆發的簡單點過程等常見條件上。
# 估計參數的方法稍後再學；這裡先理解，模型必須同時對有事件與沒事件的地方
# 負責。換到空間與規模時也一樣，不能只在看見震央的區域積分。
#
# ## 4.4 同一份歷史，不同的更新方式
#
# 上一章比較了不同序列。現在反過來，固定一份合成序列，讓三種模型讀取完全
# 相同的歷史。固定率模型維持不變；自激發模型在事件後提高率，再逐漸衰減；
# 自我修正示例在事件後降低率，之後隨時間重新增加。圖中的曲線是三種假設的
# 反應方式，並不是分別從各自模型生成的三份目錄，也不是三個擬合結果。
#
# 對地震而言，自激發結構適合表達後續活動增加的現象；累積與釋放的想法則
# 能連到某些長期模型。適用哪一個，需要看問題尺度與資料，不能僅憑模型名稱
# 就認定它對應了唯一的物理機制。
#
# %% tags=["remove-input"]
rng = np.random.default_rng(31)
T_DEMO, M0, BETA = 60.0, 3.0, np.log(10.0)           # 門檻規模、GR 斜率（b = 1）

t_ev = np.sort(rng.uniform(0, T_DEMO, 16))
m_ev = M0 + rng.exponential(1 / BETA, t_ev.size)
grid = np.linspace(0, T_DEMO, 3000)

MU_B = 0.30                                          # 共同的基線率
lam_poi = np.full_like(grid, MU_B)

A_B, ALPHA_B, C_B, P_B = 0.45, 0.80, 0.30, 1.40      # self-exciting
lam_se = np.full_like(grid, MU_B)
for ti, mi in zip(t_ev, m_ev):
    msk = grid > ti
    lam_se[msk] += (A_B * np.exp(ALPHA_B * (mi - M0)) * (P_B - 1) / C_B
                    * (1 + (grid[msk] - ti) / C_B) ** (-P_B))

B_S = 0.30                                           # self-correcting
C_S = T_DEMO / t_ev.size                             # 加載與釋放長期平衡
n_past = np.searchsorted(t_ev, grid, side="right")
lam_sc = np.exp(np.log(MU_B) + B_S * (grid - C_S * n_past))

fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                    subplot_titles=("無記憶：Poisson λ* = μ",
                                    "正記憶：self-exciting（Hawkes／ETAS）",
                                    "負記憶：self-correcting（應力釋放）"))
for row, (lam, color) in enumerate(
        [(lam_poi, PALETTE[0]), (lam_se, PALETTE[1]), (lam_sc, PALETTE[2])], 1):
    fig.add_trace(go.Scatter(x=grid, y=lam, mode="lines", showlegend=False,
                             line=dict(color=color, width=1.8)), row=row, col=1)
    fig.add_trace(go.Scatter(x=t_ev, y=np.zeros_like(t_ev), mode="markers",
                             name="事件", showlegend=(row == 1),
                             marker=dict(color=QUAKE_COLOR, size=7,
                                         symbol="triangle-up")), row=row, col=1)
    fig.update_yaxes(title_text="λ*（次/天）", row=row, col=1)
fig.update_xaxes(title_text="時間（天）", row=3, col=1)
apply_layout(fig, title=f"同一份事件序列（{t_ev.size} 個），三種記憶結構",
             height=640, hovermode="x")
fig

# %% [markdown]
# 三張圖的事件時刻相同，率的跳動方向卻不同。這就是條件強度的用途：讓
# 「歷史怎麼影響下一步」從一句敘述變成一條可比較的曲線。ETAS會進一步
# 讓過去事件的規模影響跳升幅度，再把位置納入。
#
# ## 4.5 加上位置與規模
#
# 真實地震還有位置與規模。我們可以把 $\lambda^*(t,x,y,m)$ 理解成每單位
# 時間、面積、規模的條件事件率。空間座標若用公里，面積單位就是平方公里；
# 用經緯度建模時則要說清楚相應測度與換算，不能直接把平方度當成平方公里。
#
# 一個常用假設是把當次規模的密度分開：
#
# $$\lambda^*(t,x,y,m)=\lambda^*(t,x,y)s(m),\qquad
# \int_{m_0}^{\infty}s(m)\,\mathrm dm=1.$$
#
# 這表示給定事件發生的位置、時間與先前歷史，當次規模按同一個 $s(m)$ 抽取。
# 它不表示整份目錄的規模與時空完全獨立：一個過去的大事件仍可影響後面事件
# 發生的時間與位置。把兩種依賴分清楚，才能理解「活動變多」與「單次規模分布
# 改變」為什麼是兩個不同問題。
#
# 這裡的規模分布是模型假設，不是無須檢查的定律。下章會從目錄完整度與
# Gutenberg–Richter關係討論如何估計它。即使概似能分成規模與時空兩部分，
# 也只有在參數不共享、沒有連動限制時，才可以分別最佳化；
# {doc}`appendix_a_point_process`保留這個條件。
#
# ## 4.6 把模型預期的快慢放進時間軸
#
# 如果模型認為某段期間本來就應該很活躍，許多事件擠在那裡不一定算異常。
# 診斷需要先扣除模型所描述的快慢。定義累積條件率，也稱補償子：
#
# $$\Lambda(t)=\int_0^t\lambda^*(u)\,\mathrm du.$$
#
# 在自激發模型中，它沿觀測歷史累積，因此通常是隨機的。適當可積性條件下
# $E[N(t)]=E[\Lambda(t)]$，不能把每一條實現的 $\Lambda(t)$ 都寫成
# 確定性的 $E[N(t)]$ 曲線。
#
# 用 $\tau_i=\Lambda(t_i)$ 替事件重新標記時間。若強度模型正確、滿足時間
# 變換所需條件，變換後的間隔具有單位指數分布，並互相獨立；等價地，
#
# $$U_i=1-e^{-(\tau_i-\tau_{i-1})}\sim U(0,1).$$ (eq:time-rescale)
#
# 直覺是讓模型預期很快的時段展開、很慢的時段縮短。若還看見系統性的密集或
# 空白，代表模型可能留下未描述的結構。下圖用已知參數的模擬目錄，對照正確
# 結構與固定率錯設；兩者都看同一份資料。
#
# Q–Q 圖將理論與樣本在相同累積比例下的分位數配成一點，若各點接近對角線，
# 表示兩個分布在這些位置的數值相近。KS 統計量則是兩條累積分布之間最大的
# 垂直距離，數值越大，表示分布差異越明顯。
#
# %% tags=["remove-input"]
C_OM, P_OM = 0.30, 1.40                # 本章共用的 Omori 核參數（天）

def simulate_branching(mu, A, alpha, T, rng, c=C_OM, p=P_OM):
    """分支（cluster）法模擬時間型 ETAS，回傳 (時間, 規模)。"""
    t = list(rng.uniform(0, T, rng.poisson(mu * T)))       # 第 0 代：背景
    m = list(M0 + rng.exponential(1 / BETA, len(t)))
    todo = list(range(len(t)))
    while todo:
        i = todo.pop()
        for _ in range(rng.poisson(A * np.exp(alpha * (m[i] - M0)))):
            dt = c * ((1 - rng.random()) ** (-1 / (p - 1)) - 1)   # 反函數法
            if t[i] + dt < T:
                t.append(t[i] + dt)
                m.append(M0 + rng.exponential(1 / BETA))
                todo.append(len(t) - 1)
    order = np.argsort(t)
    return np.asarray(t)[order], np.asarray(m)[order]

def compensator(t, ts, ms, mu, A, alpha, c=C_OM, p=P_OM):
    """Λ(t) = ∫λ*：Omori 核的積分有解析形式。"""
    past = ts < t
    trig = (A * np.exp(alpha * (ms[past] - M0))
            * (1 - (1 + (t - ts[past]) / c) ** (-(p - 1))))
    return mu * t + trig.sum()

MU_R, ALPHA_R, T_R, N_RATIO = 0.30, 0.80, 400.0, 0.65
A_R = N_RATIO * (BETA - ALPHA_R) / BETA          # 由分支比反推產能尺度

t_cat, m_cat = simulate_branching(MU_R, A_R, ALPHA_R, T_R,
                                  np.random.default_rng(42))
N_cat = t_cat.size
j_idx = np.arange(1, N_cat + 1)
q_theory = (j_idx - 0.5) / N_cat

tau_ok = np.array([compensator(ti, t_cat, m_cat, MU_R, A_R, ALPHA_R)
                   for ti in t_cat])
tau_bad = N_cat * t_cat / T_R                    # 錯設：均勻 Poisson，總數校準

fig = make_subplots(rows=2, cols=2, vertical_spacing=0.13, horizontal_spacing=0.10,
                    subplot_titles=("正確模型：τ 對 j", "錯設模型：τ 對 j",
                                    "正確模型：U 的 Q–Q 圖",
                                    "錯設模型：U 的 Q–Q 圖"))
ks = {}
for col, (tau, name, color) in enumerate(
        [(tau_ok, "正確", PALETTE[2]), (tau_bad, "錯設", PALETTE[1])], 1):
    fig.add_trace(go.Scatter(x=j_idx, y=tau, mode="lines", showlegend=False,
                             line=dict(color=color, width=1.6)), row=1, col=col)
    fig.add_trace(go.Scatter(x=j_idx, y=j_idx, mode="lines", showlegend=False,
                             line=dict(color="#999999", width=1, dash="dash")),
                  row=1, col=col)
    u = np.sort(1 - np.exp(-np.diff(np.concatenate([[0.0], tau]))))
    ks[name] = max(np.max(j_idx / N_cat - u),
                   np.max(u - (j_idx - 1) / N_cat))
    fig.add_trace(go.Scattergl(x=q_theory, y=u, mode="markers", showlegend=False,
                               marker=dict(color=color, size=4, opacity=0.7)),
                  row=2, col=col)
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", showlegend=False,
                             line=dict(color="#999999", width=1, dash="dash")),
                  row=2, col=col)
    fig.update_xaxes(title_text="事件序號 j", row=1, col=col)
    fig.update_xaxes(title_text="理論分位數", row=2, col=col)
fig.update_yaxes(title_text="τ = Λ(t)", row=1, col=1)
fig.update_yaxes(title_text="U 的樣本分位數", row=2, col=1)
apply_layout(fig, height=680, hovermode="closest",
             title=f"隨機時間變換診斷（N = {N_cat}）："
                   f"KS 統計量 正確 {ks['正確']:.3f}／錯設 {ks['錯設']:.3f}")
fig

# %% [markdown]
# 上排橫軸是事件序號 $j$，縱軸是累積條件率 $\tau_j$。在一段活動比模型
# 預期密集的期間，很多新事件對應很少的模型時間，這一段曲線會偏平；相反的
# 安靜期則可能形成較陡的增量。不能只看曲線位於對角線哪一側，就把此前累積
# 偏差當成當下的活動變化。
#
# 下排的Q–Q圖延續{doc}`foundation_randomness`的間隔比較，只是先將事件搬到
# 模型所定義的時間軸，再與共同基準比較。接近對角線是符合模型的一個線索，卻不是完整
# 證明，因為還要檢查間隔間的相依。圖中的KS數字是描述性的分布差異；若參數
# 是用同一份資料估出來，正式檢定需要考慮估計與觀察窗的影響，常以重新擬合的
# 模擬來校準，不能直接照搬已知分布的臨界值。
#
# 殘差也可能反映目錄漏測、位置或規模處理改變，而非新的地震物理現象。
# 因此下一章會先回到觀測的完整度，再繼續組合模型。
#
# ## 4.7 從模型生成另一份可能的目錄
#
# 前面的圖用模擬比較模型。生成一份目錄時，先給定參數與起始歷史，
# 依模型產生下一次事件，更新歷史，再繼續往前。重複許多次，就能看見相同
# 條件下的不同未來，也能檢查某種診斷圖在模型正確時本來會波動多少。
#
# 一種方法稱為稀疏化：先用較高的率提出候選時刻，再按照真正的率與提案率
# 的比例保留候選點。像用較密的網先找可能位置，再做第二次篩選。這需要有效
# 的強度上界，而不要求所有積分都有解析解。本章的核在兩個事件之間單調遞減，
# 因此當下的強度就是接下來候選區間的有效上界；若換成先升後降的核，必須另找
# 覆蓋該區間的上界，不能直接沿用這個選擇。
#
# 下面的候選點中，有些保留，有些拒絕。保留事件會改變接下來的率；拒絕的
# 候選不是實際事件，不應作為觸發源加入歷史。
#
# %% tags=["remove-input"]
def lam_star_t(s, ts, ms, mu, A, alpha, c=C_OM, p=P_OM):
    """時間型 ETAS 的條件強度（ts 只含 s 之前的事件）。"""
    if ts.size == 0:
        return mu
    return mu + np.sum(A * np.exp(alpha * (ms - M0)) * (p - 1) / c
                       * (1 + (s - ts) / c) ** (-p))

def simulate_thinning(mu, A, alpha, T, rng, record=False):
    """Ogata thinning：回傳 (時間, 規模[, 提案紀錄])。"""
    t, ts, ms, rec = 0.0, [], [], []
    while True:
        at, am = np.asarray(ts), np.asarray(ms)
        lam_bar = lam_star_t(t, at, am, mu, A, alpha)      # 上界＝當下強度
        s = t + rng.exponential(1.0 / lam_bar)
        if s > T:
            break
        v = rng.uniform(0.0, lam_bar)
        accept = v <= lam_star_t(s, at, am, mu, A, alpha)
        if record:
            rec.append((t, s, v, lam_bar, accept))
        if accept:
            ts.append(s)
            ms.append(M0 + rng.exponential(1 / BETA))
        t = s
    out = (np.asarray(ts), np.asarray(ms))
    return out + (rec,) if record else out

t_th, m_th, rec = simulate_thinning(0.30, A_R, ALPHA_R, 25.0,
                                    np.random.default_rng(5), record=True)
gg = np.linspace(0, 25, 4000)
lam_g = np.array([lam_star_t(s, t_th[t_th < s], m_th[t_th < s],
                             0.30, A_R, ALPHA_R) for s in gg])

fig = go.Figure()
step_x, step_y = [], []
for t0, s, v, lb, acc in rec:
    step_x += [t0, s, None]
    step_y += [lb, lb, None]
fig.add_trace(go.Scatter(x=step_x, y=step_y, mode="lines", name="上界 λ̄（階梯）",
                         line=dict(color="#999999", width=1.4, dash="dot")))
fig.add_trace(go.Scatter(x=gg, y=lam_g, mode="lines", name="條件強度 λ*(t)",
                         line=dict(color=ACCENT, width=2)))
acc_pts = [(s, v) for _, s, v, _, a in rec if a]
rej_pts = [(s, v) for _, s, v, _, a in rec if not a]
for pts, name, color, sym in [(rej_pts, "提案被拒絕", PALETTE[3], "x"),
                              (acc_pts, "提案被接受", QUAKE_COLOR, "circle")]:
    px, py = zip(*pts)
    fig.add_trace(go.Scatter(x=list(px), y=list(py), mode="markers", name=name,
                             marker=dict(color=color, size=8, symbol=sym,
                                         opacity=0.85)))
apply_layout(fig, yaxis_type="log",
             title=f"Thinning：{len(rec)} 個提案點，接受 {len(acc_pts)} 個"
                   f"（接受率 {len(acc_pts) / len(rec):.0%}）",
             xaxis_title="時間（天）", yaxis_title="強度（次/天，對數軸）",
             height=460, hovermode="closest")
fig

# %% [markdown]
# 稀疏化讓我們從強度曲線往前生成事件。對可作分支表示的自激發模型，還有
# 另一條路：先產生背景事件，再讓每個事件依規則產生後代，後代又可以繼續
# 產生下一代。這條路直接顯示了ETAS名稱中「類流行病」的來源。
#
# 兩種方法只有在參數化、起始歷史、觀察邊界等設定一致時，才是在模擬同一個
# 模型。不能拿空歷史起步的一份目錄，與包含長期既往活動的一份目錄直接比較，
# 再把差異全歸因於演算法。
#
# 下圖在相同參數下比較兩種方法生成的事件數與間隔分布。兩者都從空歷史開始，
# 並只保留固定時間窗內的事件，所以比較的是這個共同設定。抽樣、核的反函式
# 與邊界處理細節集中在{doc}`appendix_a_point_process`。
#
# %% tags=["remove-input"]
MU_C, T_C, N_REP = 0.30, 60.0, 260
rng_b = np.random.default_rng(101)
rng_t = np.random.default_rng(202)

n_branch, n_thin, dt_branch, dt_thin = [], [], [], []
for _ in range(N_REP):
    tb, _mb = simulate_branching(MU_C, A_R, ALPHA_R, T_C, rng_b)
    tt, _mt = simulate_thinning(MU_C, A_R, ALPHA_R, T_C, rng_t)
    n_branch.append(tb.size)
    n_thin.append(tt.size)
    dt_branch.append(np.diff(tb))
    dt_thin.append(np.diff(tt))

n_branch, n_thin = np.array(n_branch), np.array(n_thin)
d_b = np.log10(np.concatenate(dt_branch) + 1e-6)
d_t = np.log10(np.concatenate(dt_thin) + 1e-6)
edges = np.arange(0, max(n_branch.max(), n_thin.max()) + 6, 5)
xs = np.linspace(-6, 2, 300)
ecdf_b = np.searchsorted(np.sort(d_b), xs) / d_b.size
ecdf_t = np.searchsorted(np.sort(d_t), xs) / d_t.size
d_max = np.abs(ecdf_b - ecdf_t).max()

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.11,
                    subplot_titles=("每份目錄的事件數 N", "事件間隔的經驗累積分布"))
for name, arr, ec, color in [("branching", n_branch, ecdf_b, PALETTE[0]),
                             ("thinning", n_thin, ecdf_t, PALETTE[1])]:
    h, _ = np.histogram(arr, bins=edges)
    fig.add_trace(go.Bar(x=(edges[:-1] + edges[1:]) / 2, y=h / arr.size,
                         name=name, legendgroup=name,
                         marker=dict(color=color, opacity=0.55)), row=1, col=1)
    fig.add_trace(go.Scatter(x=xs, y=ec, mode="lines", name=name,
                             legendgroup=name, showlegend=False,
                             line=dict(color=color, width=2)), row=1, col=2)
fig.update_xaxes(title_text="N", row=1, col=1)
fig.update_yaxes(title_text="相對次數", row=1, col=1)
fig.update_xaxes(title_text="log₁₀ 間隔（天）", row=1, col=2)
fig.update_yaxes(title_text="累積比例", row=1, col=2)
apply_layout(fig, barmode="overlay", height=430, hovermode="x",
             title=f"兩種模擬的一致性（各 {N_REP} 份目錄）："
                   f"平均 N＝{n_branch.mean():.1f} vs {n_thin.mean():.1f}，"
                   f"間隔分布最大差距 {d_max:.3f}")
fig

# %% [markdown]
# 兩種模擬得到的分布可用來核對實作是否相容；有限次模擬的接近仍不是等價性
# 的證明。分支表示的理論條件與演算法的數值近似，也需要另外確認。
#
# 生成目錄的思路還會出現在除叢。模擬時我們知道每個事件來自哪個親代；
# 真實目錄看不到這些關係，只能在模型下估計背景與不同觸發來源的權重。
# 這正是後面統計推論與ETAS估計之間的連結。
#
# ## 4.8 帶著共同語言回到資料
#
# 條件強度描述歷史如何改變率；沿著這個定義，
# 概似讓我們比較參數，時間變換幫助檢查遺漏的結構，模擬則提供可能的未來。
# 每一步都有自己的條件，並不是只要寫下一條率曲線，所有問題就自動解決。
#
# 接下來會一直用到三個不同的門檻：$M_c$ 是目錄完整度，$m_0$ 是模型納入
# 事件的下限，$m_T$ 是預報目標規模。它們可能不同，而且目的不同。
# 若漏掉了本來會影響後續活動的小事件，模型的歷史也會隨之改變。
#
# 因此在{doc}`11_catalog_completeness_b`，我們先暫停增加模型複雜度，
# 檢查規模分布與觀測門檻。之後的叢集律與估計，會把這些資料條件帶回來。
# 需要回查符號、概似推導或模擬公式時，可使用{doc}`appendix_a_point_process`，
# 第一次閱讀則可以直接往下一章走。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - Gallager, R. G.／MIT OpenCourseWare（2011），[Discrete Stochastic Processes，第2章](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/resources/mit6_262s11_chap02/)，免費講義。用等待時間與事件計數理解Poisson特例，適合回查條件強度之前的機率基礎。
# - Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I，第二版](https://doi.org/10.1007/b97277)。第7.2–7.6節是本章條件強度、標記、概似、時間變換與模擬的核心來源；全文需訂閱或館藏權限。
# - Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)；[免費作者預印本](https://arxiv.org/abs/1708.02647)。第2–3節將模型結構、估計與診斷串在一起，適合接到ETAS章。
# - Ogata, Y.（1988），[Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes](https://doi.org/10.1080/01621459.1988.10478560)；[研究機構提供的免費全文](https://bemlar.ism.ac.jp/zhuang/Refs/Refs/ogata1988.pdf)。對照地震模型及時間變換的實際用途，注意診斷相對於哪一個模型。
# - Jalilian, A.（2019），[ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01)，[免費全文](https://www.jstatsoft.org/article/view/v088c01)。延伸時空模型、補充歷史與邊界處理；理解文章的模型部分不必先學R。
