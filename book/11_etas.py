# %% [markdown]
# # 11. ETAS：把餘震邏輯寫成一條式子
#
# {doc}`第 5 章 <05_seismic>`與{doc}`第 10 章 <10_seismicity_stats>`
# 給了你兩條經驗律：GR 律（規模怎麼分布）與 Omori 律（餘震怎麼
# 衰減）。這一章要把它們縫成一個完整的隨機
# 過程模型——**ETAS**（Epidemic-Type Aftershock Sequence，
# 傳染型餘震序列模型，Ogata 1988）。
#
# 模型的名字已經把核心想法講完了：**地震像傳染病**。每個地震都是
# 一個感染源，會以一定的「繁殖力」觸發下一代地震；下一代又觸發
# 下下代。餘震的餘震、前震與主震，全部從這一條規則自動長出來。
# ETAS 是目前世界各國作業化地震預報（第 9、13 章）的主力引擎，
# 也是所有新模型都必須先贏過的基準。
#
# ## 11.1 從 Omori 到 ETAS：疊加這一步
#
# 修正 Omori 律描述**單一**主震後的餘震率 $n(t)=K/(t+c)^p$。
# 但真實序列裡，大餘震自己也會帶出一群餘震——1999 集集、2024
# 花蓮都看得到主震衰減途中疊上去的小尖峰。與其手動判斷「這是
# 主震的餘震還是餘震的餘震」，ETAS 的解法乾脆得多：**讓每一個
# 事件都掛上自己的 Omori 核，全部疊加起來。**
#
# 寫成式子，就是**條件強度函數（conditional intensity）**——
# 在時刻 $t$、已知全部歷史 $H_t$ 之下，單位時間的事件發生率：
#
# $$\lambda(t \mid H_t) = \mu + \sum_{i:\,t_i<t}
#   \underbrace{A\,e^{\alpha(m_i-m_0)}}_{\text{繁殖力}}\;
#   \underbrace{\frac{p-1}{c}\Bigl(1+\frac{t-t_i}{c}\Bigr)^{-p}}_{\text{Omori 時間核}}$$
#
# 逐項讀它：
#
# - $\mu$：**背景率**——不是被任何地震觸發、由板塊加載慢慢累積
#   出來的「自發」地震；
# - $A\,e^{\alpha(m_i-m_0)}$：規模 $m_i$ 的事件平均觸發的**直接
#   後代數**。$\alpha$ 控制「規模換算成觸發力」的效率：$\alpha$ 大
#   表示大地震壓倒性主導（典型主震–餘震型，日本經驗值 1.2–3.1）；
#   $\alpha$ 小表示大小地震觸發力差不多（群震型，0.35–0.85）；
# - Omori 核（已正規化成機率密度）：後代出生時間的分布；
# - 每個後代的**規模**則從 GR 律獨立抽出（$\beta = b\ln 10$）——
#   注意：**與親代規模無關**，所以後代完全可以比親代大。
#
# 最後一點值得停一秒。ETAS 的世界裡沒有「主震」這個概念，只有
# 「背景事件」與「被觸發事件」；我們慣稱的前震，只是恰好觸發了
# 更大後代的事件。第 9 章說「前震是事後標籤」，ETAS 把這句話
# 直接寫進了數學結構裡。
#
# ## 11.2 分支結構與臨界性
#
# 把「平均一個事件生幾個後代」對所有規模取期望，得到**分支比
# （branching ratio）**：
#
# $$n = \frac{A\,\beta}{\beta - \alpha} \qquad (\alpha < \beta)$$
#
# 它是整個模型最重要的一個數字，意義等同傳染病的基本再生數
# $R_0$：
#
# - $n < 1$（次臨界）：每代平均萎縮，任何序列終將熄滅，目錄中
#   平均有 $n$ 的比例是被觸發事件；
# - $n \to 1$（臨界）：序列可以拖得極長；
# - $n > 1$（超臨界）：級聯爆炸，模型發散——真實地球顯然不是
#   這樣，但不同地區、不同目錄擬出的 $n$ 多在 0.3–0.9 之間，**離臨界不遠**。
#
# 「地球在次臨界但接近臨界」這件事，正是為什麼一次大地震後的
# 序列又長又猛、卻終究會停。用模擬直接看 $n$ 的威力——下面
# 手刻一個最簡的**分支法模擬**（branching simulation）：先抽
# 背景事件，再讓每個事件按繁殖力生 Poisson 個後代、後代時間從
# Omori 核抽、規模從 GR 律抽，遞迴到沒有新事件為止：

# %% tags=["remove-input"]
import plotly.io as pio
pio.renderers.default = "notebook_connected"

# %% tags=["hide-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

M0, BETA = 3.0, np.log(10)          # 門檻規模、GR 斜率（b = 1）
C_OM, P_OM = 0.01, 1.1              # Omori 核參數（天）

def simulate_etas(mu, n_ratio, alpha, T, seed):
    """時間型 ETAS 的分支法模擬。回傳 (時間, 規模, 親代索引)。"""
    rng = np.random.default_rng(seed)
    A = n_ratio * (BETA - alpha) / BETA          # 由分支比反推繁殖力尺度
    t = list(rng.uniform(0, T, rng.poisson(mu * T)))      # 第 0 代：背景
    m = list(M0 + rng.exponential(1 / BETA, len(t)))
    parent = [-1] * len(t)
    todo = list(range(len(t)))
    while todo:
        i = todo.pop()
        for _ in range(rng.poisson(A * np.exp(alpha * (m[i] - M0)))):
            dt = C_OM * ((1 - rng.random()) ** (-1 / (P_OM - 1)) - 1)
            if t[i] + dt < T:
                t.append(t[i] + dt)
                m.append(M0 + rng.exponential(1 / BETA))
                parent.append(i)
                todo.append(len(t) - 1)
    return np.array(t), np.array(m), np.array(parent)

t_e, m_e, par = simulate_etas(mu=0.4, n_ratio=0.85, alpha=1.7, T=180, seed=7)

# λ(t) 曲線
tt = np.linspace(0, 180, 4000)
lam = np.full_like(tt, 0.4)
A_e = 0.85 * (BETA - 1.7) / BETA
for ti, mi in zip(t_e, m_e):
    mask = tt > ti
    lam[mask] += (A_e * np.exp(1.7 * (mi - M0)) * (P_OM - 1) / C_OM
                  * (1 + (tt[mask] - ti) / C_OM) ** -P_OM)

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.45, 0.55],
                    vertical_spacing=0.04)
fig.add_trace(go.Scatter(x=tt, y=lam, mode="lines", name="λ(t) 條件強度",
                         line=dict(color=QUAKE_COLOR, width=1.5)), row=1, col=1)
fig.add_trace(go.Scattergl(x=t_e, y=m_e, mode="markers", name="模擬事件",
                           marker=dict(size=4 + (m_e - M0) * 4, color=ACCENT,
                                       opacity=0.6)), row=2, col=1)
fig.update_yaxes(title_text="λ（次/天）", type="log", row=1, col=1)
fig.update_yaxes(title_text="規模", row=2, col=1)
fig.update_xaxes(title_text="時間（天）", row=2, col=1)
apply_layout(fig, title=f"模擬的 ETAS 目錄（n = 0.85，共 {len(t_e)} 個事件）",
             height=520, hovermode="x")
fig

# %% [markdown]
# 上圖是條件強度 $\lambda(t)$（對數軸），下圖是事件本身。每次較大
# 事件一出現，$\lambda$ 立刻跳升、再按 Omori 律衰減；衰減途中若又
# 冒出一個大事件，強度再度跳升。整份目錄裡沒有寫死任何「序列」，
# 但主震–餘震的樣子自己長了出來。
#
# 分支結構還可以直接畫出來——模擬時我們記錄了每個事件的親代，
# 把親子用線連起來：

# %% tags=["hide-input"]
big = np.argsort(m_e)[-1]           # 最大事件
fam = {big}                          # 找出它的整個家族（祖先＋後代）
j = big
while par[j] >= 0:                   # 祖先鏈
    j = par[j]
    fam.add(j)
changed = True                       # 後代（反覆掃描直到閉合）
while changed:
    changed = False
    for k in range(len(t_e)):
        if par[k] in fam and k not in fam:
            fam.add(k)
            changed = True

fig = go.Figure()
edge_x, edge_y = [], []
for k in range(len(t_e)):
    if par[k] >= 0:
        edge_x += [t_e[par[k]], t_e[k], None]
        edge_y += [m_e[par[k]], m_e[k], None]
fig.add_trace(go.Scattergl(x=edge_x, y=edge_y, mode="lines", name="親→子",
                           line=dict(color="#bbbbbb", width=0.7)))
is_fam = np.isin(np.arange(len(t_e)), list(fam))
for sel, name, color in [(~is_fam, "其他事件", ACCENT),
                         (is_fam, "最大事件的家族", QUAKE_COLOR)]:
    fig.add_trace(go.Scattergl(x=t_e[sel], y=m_e[sel], mode="markers",
                               name=name,
                               marker=dict(size=4 + (m_e[sel] - M0) * 4,
                                           color=color, opacity=0.75)))
apply_layout(fig, title="觸發樹：每個事件都有族譜",
             xaxis_title="時間（天）", yaxis_title="規模",
             hovermode="closest", height=480)
fig

# %% [markdown]
# 灰線是親子關係，紅色是模擬中最大事件所屬的整個家族。注意兩件事：
# 一、家族可以枝繁葉茂好幾「代」——二次、三次餘震不是特例而是
# 常態；二、最大的事件**不一定是家族的始祖**——它可能是某個較小
# 事件的後代。在傳統語彙裡，那個較小的祖先就叫「前震」。ETAS
# 不需要任何前震機制，前震自動出現。
#
# 再看分支比的效果。同樣的背景率，把 $n$ 從 0.5 調到 0.9：

# %% tags=["hide-input"]
fig = go.Figure()
for n_r, color in [(0.5, PALETTE[2]), (0.9, PALETTE[1])]:
    ts, ms, _ = simulate_etas(mu=0.4, n_ratio=n_r, alpha=1.7, T=365, seed=11)
    daily = pd.Series(1, index=pd.to_datetime("2024-01-01")
                      + pd.to_timedelta(ts, unit="D")).resample("D").sum()
    fig.add_trace(go.Scatter(x=np.arange(len(daily)), y=daily.values,
                             mode="lines", name=f"n = {n_r}（共 {len(ts)} 個）",
                             line=dict(color=color, width=1.2)))
apply_layout(fig, title="接近臨界的世界：分支比 n 對序列樣貌的影響",
             xaxis_title="時間（天）", yaxis_title="每日事件數",
             yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# $n=0.5$（綠）的世界裡，級聯短促，活動大致貼著背景率起伏；
# $n=0.9$（橘）的世界裡，同樣的背景事件三不五時引爆連鎖，單日
# 事件數可以衝高一兩個數量級。真實地殼比較像後者——這就是為什麼
# 「一半以上的地震是餘震」在多數地區都成立。
#
# ## 11.3 加上空間、擬合真實目錄
#
# 完整的**時空 ETAS** 再補兩塊：背景率變成空間場 $\mu(x,y)$
# （用平滑核從目錄估計），每個事件的觸發再乘上一個空間核——
# 餘震區面積隨親代規模指數放大（Utsu–Seki 律），距離上做冪次
# 衰減。參數增為八個：$(\mu, A, c, \alpha, p, D, q, \gamma)$。
#
# 估計靠點過程的最大概似法。對數概似只有兩項：
#
# $$\ln L = \sum_i \ln \lambda(t_i, x_i, y_i \mid H_{t_i})
#   - \int\!\!\!\int\!\!\!\int \lambda \, dt\, dx\, dy$$
#
# 第一項獎勵「事件真的發生的地方強度要高」，第二項懲罰「到處
# 都說高」。這個「$\sum\ln\lambda - \int\lambda$」骨架是一切
# 點過程模型的共同心臟，第 12 章的 EEPAS 也是同一副。
#
# 台灣已經有自己的時空 ETAS 參數，而且有兩份**互相獨立**的估計：
# 中央氣象署 112 年委辦計畫（1994–2021 目錄、$M_c=3.0$）得到
# $p = 1.06$、$\alpha = 1.17$；2025 年大埔地震快報（Hsieh et al.
# 2025）得到 $p = 1.04$、$\alpha = 1.04$。兩組數字高度一致，
# 互為交叉驗證——這就是「台灣的 ETAS 長什麼樣」的基準值。
# 大埔快報同時示範了作業化的工程面：參數**預先**用三十年目錄
# 訓練好（128 核跑一小時），地震來了之後每輪預報只需 8 核、
# 7 分鐘——先做重活、即時只做輕活，是所有作業化系統的共同
# 設計。它發布的「震後 10 天內 $M_L\ge5$ 機率 67.8%」對上實際
# 發生 8 次，這類數字怎麼算「準不準」，是第 15 章的主題。
#
# 擬合真實序列時，第 10 章的每個坑都會回來咬人。看 0403 花蓮
# 序列的餘震衰減與 Omori–Utsu 擬合：

# %% tags=["hide-input"]
from scipy.optimize import curve_fit

cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
main = cat.loc[cat.ML.idxmax()]
aft = cat[(cat.time > main.time) & (cat.ML >= 3.5)]
t_days = np.sort((aft.time - main.time).dt.total_seconds() / 86400)
t_days = t_days[t_days <= 60]

bins = np.logspace(-2, np.log10(60), 25)
counts, _ = np.histogram(t_days, bins=bins)
width = np.diff(bins)
centers = np.sqrt(bins[:-1] * bins[1:])
rate = counts / width

def omori(t, K, c, p):
    return K / (t + c) ** p

ok = rate > 0
popt, _ = curve_fit(lambda t, K, c, p: np.log(omori(t, K, c, p)),
                    centers[ok], np.log(rate[ok]),
                    p0=[100, 0.05, 1.1], maxfev=10000)
fig = go.Figure(go.Scatter(x=centers[ok], y=rate[ok], mode="markers",
                           name="觀測（ML≥3.5）",
                           marker=dict(color=ACCENT, size=8)))
fig.add_trace(go.Scatter(x=centers, y=omori(centers, *popt), mode="lines",
                         name=f"Omori–Utsu 擬合：p = {popt[2]:.2f}，"
                              f"c = {popt[1]:.3f} 天",
                         line=dict(color=QUAKE_COLOR, dash="dash")))
apply_layout(fig, title="0403 花蓮序列的餘震率衰減（主震後 60 天）",
             xaxis_title="主震後時間（天）", yaxis_title="餘震率（次/天）",
             xaxis_type="log", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 注意擬合出的 $p$ 略低於 1——先把這個數字記住，它馬上會變成
# 下面陷阱清單的活例子。同時回想第 10 章：池上序列的 $p$
# 從 0.92 到 1.39，只因時間窗不同；而擬到的 $c$ 值有多少是物理、
# 多少是主震後頭幾小時目錄不完整的偽裝，沒有人能單從這張圖回答。
# 參數估計還有幾個 ETAS 特有的陷阱，值得抄下來：
#
# - **參數之間高度相關**。文獻中有擬合出 $\gamma = 0.02$、
#   標準誤卻是 5.8 的真實案例（Jalilian 2019 的伊朗目錄）——
#   估計值比誤差小兩個數量級，等於這個參數完全沒被資料約束。
#   看 ETAS 參數表**必須**連標準誤一起看。
# - **早期不完整會壓低 $\alpha$**：大震剛過的高產出期被目錄
#   「咬掉」，繁殖力被系統性低估。修正不完整性之後 $\alpha$ 會
#   回升到接近 $\beta$——所以有一派模型（simplETAS，Mancini &
#   Marzocchi 2023）乾脆把描述叢集的七個參數（含 $\alpha=\beta$）全部
#   釘死在普世經驗值，只留背景率與繁殖力尺度兩個自由參數
#   （$\alpha=\beta$ 時分支比要搭配有限的最大規模才收斂——
#   這類自相似設定的標準配件），
#   在義大利從「日」到「四個世紀」的檢驗全部通過。參數多不如
#   參數穩，這是 bias–variance 取捨的地震學版本。
# - **$p<1$ 通常不是發現，是模型設定錯了**：把空間不均勻的
#   背景率硬塞進均勻假設，多出來的事件會被 Omori 核吸收，
#   把 $p$ 拉到 1 以下。上圖的擬合值正是一例：單一序列、混著
#   早期不完整與鄰近的背景活動，$p$ 就這樣被壓低了一點。
#
# 最後，ETAS 還免費送一個第 10 章遺留問題的解法。擬合好的模型
# 可以對每一對事件算出「$j$ 被 $i$ 觸發的機率」$\rho_{ij}$ 與
# 「$j$ 是背景事件的機率」$\phi_j$——這就是**隨機除叢**
# （stochastic declustering）：不再武斷二分主震／餘震，而是給
# 機率，還能用重抽樣把除叢的不確定性顯示出來。第 10 章除叢法
# 的選擇效應難題，在 ETAS 框架裡根本不會發生。
#
# ## 11.4 一條描述叢集的語言
#
# 這一章結束前，把一個容易誤會的定位講清楚：**ETAS 不是物理
# 定律，它是描述叢集的統計語言。**它的三個零件——GR、Omori、
# Utsu–Seki——全是先從資料歸納、再寫成公式的經驗律；空間核的
# 函數形式靠 AIC 挑選，不靠力學推導；不同構造區的參數各不相同，
# 沒有第一原理告訴你該是多少。它有物理上說得通的詮釋（例如
# $q=1.5$ 對應靜態應力隨距離三次方衰減），但那是事後的註解，
# 不是模型的來源。
#
# 正因為它只是語言，它的極限也很清楚：ETAS 的預報技巧幾乎全部
# 來自「大地震剛發生之後」的短期窗口。對「平靜期裡下一個大地震
# 什麼時候來」，ETAS 給的答案跟背景率差不多——約等於沒有答案。
# 想往前多看幾個月、幾年，需要完全不同的想法：不問「這個地震
# 會觸發什麼」，而問「這個地震**預示**著什麼」。
# {doc}`下一章 <12_eepas_ppe>`的 EEPAS 就是把時間之箭反過來的
# 模型——每一個小地震，都可能是未來更大地震的前兆。
