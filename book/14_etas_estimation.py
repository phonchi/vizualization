# %% [markdown]
# # 14. ETAS II：估計、診斷、變體與作業引擎
#
# {doc}`第 13 章 <13_etas_structure>`把 ETAS 的**結構**攤開了：條件強度
# 是背景率加上一群觸發核的疊加，觸發核由產能 $\kappa(m)=Ae^{\alpha(m-m_0)}$、
# 正規化 Omori 密度 $g(t)$ 與空間核 $f(x,y;m)$ 三塊相乘而成，規模則從
# GR 密度 $s(m)$ 獨立抽出；分支比 $n=A\beta/(\beta-\alpha)$ 決定整個
# 系統離臨界有多遠。
#
# 結構寫好了，模型還不能用——那八個參數
# $\theta=(\mu, A, \alpha, c, p, D, \gamma, q)$ 得從**資料**裡長出來，
# 而這件事比寫下模型難得多。本章回答四個問題：**怎麼估**（概似的
# 形狀、邊界效應、擬牛頓法、標準誤從哪來）、**估出來的數字能不能信**
# （有些參數根本沒被資料約束）、**怎麼知道模型錯了**（殘差與隨機
# 時間變換），以及**要不要少估幾個**（simplETAS 釘死七個參數的逐項
# 理由）。最後把兩個真實系統攤開：加州的 Reasenberg–Jones 公式其實
# 是 ETAS 的解析特例，而台灣已經有兩組互相獨立的八參數估計。
#
# 前面三條式子本章一律引用不重推：第 10 章的點過程對數概似
# {eq}`eq:pp-loglik`、隨機時間變換 {eq}`eq:time-rescale`，以及第 12 章
# 的正規化 Omori 密度 {eq}`eq:omori-density`。可分離性的拆解在
# {doc}`10.4 <10_point_process>`，分支比屬於第 13 章。
#
# 這章還埋著一條暗線。第 12 章留下的難題是：除叢演算法會在 b 值上
# 留下自己的指紋，而「主震」的定義本身不可驗證。ETAS 的答案在
# 14.4 節——**不做二分，只給機率**。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 14.1 概似與邊界效應
#
# ### 兩塊獨立的概似
#
# ETAS 的條件強度是可分離的（{doc}`10.4 <10_point_process>`）：
# $\lambda^*(t,x,y,m)=\lambda^*_\theta(t,x,y)\,s(m)$，
# $s(m)=\beta e^{-\beta(m-m_0)}$。代進時空標記版的
# {eq}`eq:pp-loglik`，規模的積分等於 1，整條對數概似裂成兩塊：
#
# $$l(\beta,\theta) = l_1(\beta) + l_2(\theta)$$
#
# $l_1$ 只含 $\beta$、有 Aki–Utsu 封閉解（第 11 章擁有），$l_2$ 只含
# 八個時空參數。這不只是省事：**$b$ 值的估計誤差不會直接污染
# $\alpha$**，兩者是分開的兩場戰役。但分離的是**估計**，不是**解讀**
# ——分支比同時吃 $\hat\beta$ 與 $\hat\alpha$，$b$ 值差 0.1，$n$ 就跟著跳。
#
# ### target 與 complementary：$\delta_i$ 是什麼
#
# $l_2$ 的第一項不是對「所有事件」求和，而是
#
# $$l_2(\theta) = \sum_{i} \delta_i \ln \lambda^*_\theta(t_i,x_i,y_i)
#   \;-\; \int_{t_s}^{t_s+T}\!\!\int_S \lambda^*_\theta
#   \,\mathrm{d}x\,\mathrm{d}y\,\mathrm{d}t$$
#
# $\delta_i=1$ 的是 **target event**（落在研究區 $S$ 與研究期
# $[t_s,t_s+T]$ 之內），$\delta_i=0$ 的是 **complementary event**
# ——研究期之前或研究區之外的事件（Jalilian 2019）。兩者待遇不同：
# complementary event **要進歷史 $H_t$**（它們確實會觸發區內事件），
# 但**不進第一項的求和**（我們沒有宣稱在預報它們）。
#
# 理由回到 {eq}`eq:pp-loglik` 的兩項——那是一場押注結算。你只在
# $S\times[t_s,t_s+T]$ 這個盒子裡下注，就只該為盒內事件領獎、為盒內
# 預期量付罰金。把盒外事件塞進第一項等於為沒下的注領獎；把它們從
# $H_t$ 拿掉更糟，研究期開頭那幾天的強度會被系統性低估，模型只好把
# 缺掉的活動硬塞進背景率——**這就是幾乎所有 ETAS 實作都要求一段
# 暖機期的原因。**
#
# ### 空間積分項：邊界會吃掉後代
#
# 把積分項逐塊展開，邊界效應就現形了：
#
# $$\begin{aligned}
# \int_{t_s}^{t_s+T}\!\!\int_S \lambda^*_\theta
#   &= T\int_S \mu(x,y)\,\mathrm{d}x\,\mathrm{d}y \\
#   &\quad + \sum_i \kappa(m_i)
#     \underbrace{\left[\int_{\max(t_i,t_s)}^{t_s+T} g(t-t_i)\,\mathrm{d}t\right]}
#       _{\displaystyle \equiv\, G_i \,\le\, 1}
#     \underbrace{\left[\int_S f(x-x_i,y-y_i;m_i)\,\mathrm{d}x\,\mathrm{d}y\right]}
#       _{\displaystyle \equiv\, F_i \,\le\, 1}
# \end{aligned}$$
#
# 兩個括號都小於 1，理由不同。$G_i<1$ 是**時間截斷**：事件 $i$ 的
# 後代有一部分生在研究期結束之後。由 {eq}`eq:omori-density` 直接積出
#
# $$G_i = 1 - \left(1 + \frac{t_s + T - t_i}{c}\right)^{1-p} ,$$
#
# 愈靠近研究期尾端的事件，$G_i$ 愈小——它來不及把後代生完。
# $F_i<1$ 則是**空間截斷**：後代有一部分掉到研究區外面。對距離邊界
# 很遠的事件 $F_i\approx1$，對貼著邊界的事件可以掉到 0.5 以下。
# **忽略 $F_i$ 會系統性高估懲罰項，最佳化只好把 $A$ 壓低來補償**
# ——這是 ETAS 實作中最常見也最隱蔽的 bug。
#
# $F_i$ 沒有初等封閉解（$S$ 是任意多邊形），但算得很準。標準做法叫
# **radial partitioning**：以事件 $i$ 為原點，把 $S$ 的邊界離散成
# $n_v$ 個節點，用徑向線段把 $S$ 切成一堆扇形，每個扇形轉極座標後
# **徑向積分有封閉解**（推導見附錄 A）：
#
# $$\int_0^{2\pi}\!\!\int_0^{R(\varphi)} f\,r\,\mathrm{d}r\,\mathrm{d}\varphi
#   = \frac{1}{2\pi}\int_0^{2\pi}
#     \left[1 - \left(1 + \frac{R(\varphi)^2}{D e^{\gamma(m_i-m_0)}}
#     \right)^{1-q}\right]\mathrm{d}\varphi$$
#
# 剩下的角度積分用數值求和。這比二維網格快兩三個數量級，而且在
# $r\to0$ 的奇異點附近不會出事——極座標的 Jacobian $r$ 剛好把 $f$ 在
# 原點的尖峰壓平。$R(\varphi)\to\infty$ 時右式趨近 1，正是應有的檢核。

# %% [markdown]
# ## 14.2 數值面：DFP、標準誤，與 $N^2$
#
# ### 為什麼是擬牛頓，不是牛頓
#
# 最大化 $l_2(\theta)$ 等價於最小化 $\xi(\theta)=-l_2(\theta)$。牛頓法
# 每步要算 $8\times8$ 的 Hessian，而 $l_2$ 每求值一次本身就是 $O(N^2)$
# 的雙重求和——解析二階導數要推 36 條式子、數值 Hessian 要多 36 次
# 求值。實務上一律改用**擬牛頓法**：只用梯度，把反 Hessian 逐步養
# 出來。ETAS 文獻的傳統選擇是 **Davidon–Fletcher–Powell（DFP）**
# （Ogata 的 Fortran 程式與 Jalilian 2019 的 R 套件都用它）：
#
# $$\theta_{k+1} = \theta_k - s_k\,H_k\,\nabla\xi(\theta_k)$$
#
# $H_k$ 是反 Hessian 的當前估計，$s_k$ 由線搜尋決定；每走一步就用
# 「走了多遠」與「梯度變了多少」更新 $H_k$（更新式見附錄 B）。
#
# 這裡有個常被忽略的副產品。收斂時
# $H_\infty \approx [\nabla^2\xi(\hat\theta)]^{-1}
# = [-\nabla^2 l_2(\hat\theta)]^{-1}$，而 $-\nabla^2 l_2$ 正是**觀測
# 資訊矩陣**。由 MLE 的漸近常態性，
#
# $$\mathrm{SE}(\hat\theta_j) = \sqrt{(H_\infty)_{jj}} .$$
#
# **所有 ETAS 論文表格裡的標準誤，都是這樣掉出來的。** 好處是不必
# 另外做 bootstrap；壞處是它建立在「概似在 $\hat\theta$ 附近近似
# 二次」上。參數貼邊界、或概似有一條長長的平坦山脊時（14.3 節），
# 這個數字只能當作「山脊有多寬」的粗略指標。而且 DFP 養出來的
# $H_\infty$ 只是**近似**——它從單位矩陣出發，收斂快時可能根本沒把
# 某些方向的曲率學到。
#
# ### 計算量與 $m_0$ 的取捨
#
# $l_2$ 每次求值都要對每個 target event 掃過它之前的所有事件：計算量
# $\propto N^2$。Jalilian（2019）的伊朗算例是 5970 個事件、4 次外層
# 迭代，18.37 分鐘；台灣大埔快報用 64,239 個事件、6 次迭代，128 核
# 約 1 小時。
#
# 這個平方律直接決定門檻規模 $m_0$ 的取法。GR 律給
# $N \propto 10^{-b\,m_0}$，所以
#
# $$\frac{\text{cost}(m_0 + \Delta)}{\text{cost}(m_0)} = 10^{-2b\Delta} .$$
#
# $b=1$ 時，門檻**抬高 0.5 個規模單位，計算量就掉到十分之一**；抬高
# 1.0 掉到百分之一。反過來說，把 $m_0$ 壓到 $M_c$ 的邊緣，代價是
# 計算量暴增**外加**不完整目錄污染參數（14.6 節）。這解釋了一個
# 乍看浪費的慣例：CWA 112 年報告明知台灣 1973 年後的完整度是
# $M_c=2.0$–$3.0$，卻**保守取上限 3.0**。
#
# ### 初始值敏感、不保證收斂
#
# Jalilian（2019）在文件裡直白寫著：`etas()` 對初始值敏感，預設值
# （$\mu=N/(4T|S|)$、$A=0.01$、$c=0.01$、$\alpha=1$、$p=1.3$、
# $D=0.01$、$q=2$、$\gamma=1$）只是粗估，**不保證收斂**。實務守則
# 四條：**多起點**（至少三組差很多的初始值，看是否收到同一點）；
# **檢查邊界**（$\hat c$ 貼下界、$\hat p$ 貼 1、$\hat q$ 貼 1 都是
# 模型或資料有問題的訊號）；**檢查 $\hat n$**（$\hat n\ge1$ 時概似的
# 積分項在穩態下發散，這組參數不能拿去模擬）；**檢查標準誤**
# ——下一節整節在講這件事。

# %% [markdown]
# ## 14.3 參數為什麼會說謊
#
# ### 一張把問題講完的表
#
# Jalilian（2019）用 ANSS 的伊朗目錄（1973–2016、$m_0=4.0$、5970 個
# 事件）擬合，得到：
#
# | 參數 | 估計值 | 標準誤 | 相對誤差 |
# |---|---|---|---|
# | $\beta$ | 5.6094 | 0.0453 | 0.8% |
# | $\mu$ | 0.5484 | 0.0133 | 2.4% |
# | $A$ | 0.1862 | 0.0519 | 28% |
# | $c$ | 0.0471 | 0.1093 | 232% |
# | $\alpha$ | 2.7071 | 0.0334 | 1.2% |
# | $p$ | 1.1548 | 0.0106 | 0.9% |
# | $D$ | 0.0160 | 0.1016 | 635% |
# | $q$ | 2.3234 | 0.0361 | 1.6% |
# | $\gamma$ | 0.0238 | **5.7553** | 24181% |
#
# 這張表是統計地震學最好的一堂課。$\beta$、$\alpha$、$p$、$q$ 被資料
# 釘得死死的（相對誤差 1% 上下）；$c$ 與 $D$ 的標準誤比估計值還大；
# 而 $\hat\gamma=0.0238$ 配上標準誤 5.7553——**估計值比誤差小兩個
# 數量級，這個參數完全沒被資料約束**。任何人拿 $\gamma\approx0.02$
# 去跟第 12 章那個 $\gamma=0.5\ln10\approx1.15$ 的理論值比對、宣稱
# 「伊朗的餘震區尺度律與 Utsu–Seki 不符」，都是把數值噪音當成物理
# 發現。
#
# ### 推導：$\gamma$ 與 $D$ 為什麼綁在一起
#
# 這不是壞運氣，是模型結構的必然。看空間核：
#
# $$f(x,y;m) = \frac{q-1}{\pi\,\sigma}
#   \left[1 + \frac{r^2}{\sigma}\right]^{-q},
#   \qquad \sigma \equiv D\,e^{\gamma\,\Delta m},
#   \quad \Delta m \equiv m - m_0 .$$
#
# **$D$ 與 $\gamma$ 只透過 $\sigma$ 進入模型**，沒有第二個管道。先做
# 一個極端假設：所有觸發源規模相同，$\Delta m \equiv \overline{\Delta m}$。
# 那麼 $(D,\gamma)\mapsto \sigma = D e^{\gamma\overline{\Delta m}}$ 是
# 一個從二維打到一維的映射：**沿著曲線
# $D(\gamma)=\sigma^\star e^{-\gamma\overline{\Delta m}}$，概似值一模
# 一樣。** 概似面不是山峰，是一條水平的山脊；Fisher 資訊矩陣奇異、
# 標準誤無窮大。$D$ 與 $\gamma$ **在數學上不可辨識**。
#
# 真實資料有規模散布，$\gamma$ 因此不是完全不可辨識——但可辨識的
# 程度**完全由規模的離散度決定**。令 $u_i=\ln\sigma_i=\ln D+\gamma\Delta m_i$，
# 記 $s_i = \partial \ell_i/\partial u_i$（$\ell_i$ 是第 $i$ 項的對數
# 概似貢獻）。連鎖律給出兩個分量的分數函數：
#
# $$\frac{\partial \ell}{\partial \ln D} = \sum_i s_i,
#   \qquad
#   \frac{\partial \ell}{\partial \gamma} = \sum_i \Delta m_i\, s_i ,$$
#
# 兩者只差一個權重 $\Delta m_i$。於是 $(\ln D,\gamma)$ 的資訊矩陣是
#
# $$I = \sum_i \iota_i
#   \begin{pmatrix} 1 & \Delta m_i \\
#                   \Delta m_i & \Delta m_i^2 \end{pmatrix},
#   \qquad \iota_i = -\,\mathbb{E}\!\left[\frac{\partial^2 \ell_i}
#     {\partial u_i^2}\right] ,$$
#
# 而 $\det I = \bigl(\sum_i \iota_i\bigr)^2
# \operatorname{Var}_{\iota}(\Delta m)$，其中 $\operatorname{Var}_\iota$
# 是以 $\iota_i$ 為權重的變異數（完整代數見附錄 C）。結論一行講完：
#
# $$\mathrm{SE}(\hat\gamma) \;\propto\;
#   \frac{1}{\sqrt{\operatorname{Var}_\iota(\Delta m)}} .$$
#
# **目錄的規模跨度愈窄，$\gamma$ 的標準誤愈大，而且是發散式地大。**
# 伊朗目錄 $m_0=4.0$，絕大多數事件落在 $M\,4.0$–$5.0$——5.7553 就是
# 這麼來的。
#
# ### 同一條引理也管 $\alpha$ 與 $A$
#
# 看產能：$\kappa(m)=Ae^{\alpha\Delta m}=\exp[\ln A + \alpha\Delta m]$。
# **這是同一個結構**：兩個參數以「截距 + 斜率 × $\Delta m$」的組合
# 進入模型，所以上面的引理原封不動適用。直觀說法：**資料只知道
# 「所有事件加起來一共觸發了多少後代」；要把總量拆成「基準產能 $A$」
# 與「規模加成 $\alpha$」，得靠不同規模的事件表現出不同產能——目錄裡
# 的事件規模都差不多時，這個拆分就沒有資訊來源。**
#
# $c$ 與 $p$ 的相關性來自別的機制。取對數：
# $\ln g(t)=\text{const}-p\ln(t+c)$。當 $t\gg c$ 時 $\ln(t+c)\approx\ln t$，
# $c$ 完全消失——**只有最早的那幾個事件攜帶 $c$ 的資訊**，而那正好是
# 目錄最不完整的時候（14.6 節）。這就是為什麼 $c$ 的標準誤在幾乎
# 每一份 ETAS 參數表裡都大得離譜。
#
# 先用合成資料把 $(A,\alpha)$ 的山脊畫出來。資料由已知參數
# （$\mu=0.30$、$A=0.20$、$\alpha=1.60$、$c=0.01$ 天、$p=1.15$）模擬，
# 固定 $\mu$、$c$、$p$，在 $(A,\alpha)$ 平面掃對數概似：

# %% tags=["hide-input"]
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
# 等高線不是同心圓，而是一條從左上斜到右下的長橢圓——這就是相關性的
# 幾何形狀。點線是「全目錄總觸發量固定」的方向，它幾乎貼著山脊長軸：
# 資料真正約束的是那個乘積，不是 $A$ 與 $\alpha$ 各自的值。網格 MLE
# 落在真值旁邊，但沿長軸走出去很遠，$\ln L$ 也只掉幾個單位。
# **只報 $\hat\alpha$ 不報 $\hat A$，讀者其實什麼都不知道。**
#
# 再看 $\gamma$–$D$。這次模擬「連親子配對都已知」的最理想情況——真實
# 估計連配對都不知道，只會更糟——並刻意把觸發源的規模跨度壓窄
# （$\Delta m \le 0.8$，模擬 $m_0$ 訂得高、目錄多為小事件的狀況）：

# %% tags=["hide-input"]
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
# 左圖的概似面沿著點線（$De^{\gamma\overline{\Delta m}}$ 固定的方向）
# 幾乎不下降——$D$ 走出去兩個數量級，只要 $\gamma$ 跟著補償，概似
# 幾乎一樣好。右圖把這件事量化：規模跨度從 2.5 縮到 0.8，$\gamma$ 的
# 95% 區間就從一條線變成一片荒原，連 0 都包進去了。$\hat D$ 反而
# 始終被估得不錯——因為資料約束的是**乘積**，而 $\overline{\Delta m}$
# 很小時乘積幾乎就是 $D$。
#
# 三條實務守則：**看 ETAS 參數表必須連標準誤一起看**（沒有標準誤的
# 參數表資訊量接近零）；**最好看聯合信賴區域**（上面兩張圖的橢圓
# 長軸都是斜的，投影到單一座標軸會丟掉最重要的訊息）；**跨研究比較
# 單一參數之前，先確認 $m_0$、規模尺度與核的正規化慣例都相同**
# （第 12 章證明過 $\gamma$ 的長度／面積慣例差一個因子 2，$A$ 與
# Ogata 的 $K$ 更是完全不能比）。

# %% [markdown]
# ## 14.4 背景率與隨機除叢
#
# 到目前為止我們一直假裝背景率 $\mu(x,y)$ 是已知的。它當然不是，而且
# 它的估計有一個雞生蛋的結構：**要估背景率，得先知道哪些事件是背景
# 事件；要知道哪些是背景事件，得先有一個模型。** ETAS 對這個循環的
# 解法，同時也是第 12 章那個除叢難題的答案。
#
# ### 推導：$\rho_{ij}$ 與 $\phi_j$
#
# 關鍵是把條件強度讀成**一堆獨立成分的疊加**。在時空點 $z=(t,x,y)$，
#
# $$\lambda^*(z) = \underbrace{\mu(x,y)}_{\lambda_0(z)}
#   + \sum_{i:\,t_i<t}
#   \underbrace{\kappa(m_i)\,g(t-t_i)\,f(x-x_i,y-y_i;m_i)}_{\lambda_i(z)} .$$
#
# ETAS 是**分支過程**：背景事件由速率 $\mu$ 的 Poisson 過程產生，
# 每個既有事件 $i$ 的直接後代由速率 $\lambda_i$ 的 Poisson 過程產生，
# 而且**給定歷史之後這些成分互相獨立**。疊加定理說：獨立 Poisson
# 疊起來仍是 Poisson，速率相加；觀測到的點來自第 $k$ 個成分的機率
# 就是速率的比例。在 $z_j$ 附近取體積 $\mathrm{d}z$ 的小盒，記
# $N_k(\mathrm{d}z)$ 為成分 $k$ 落在盒內的事件數、
# $N(\mathrm{d}z)=\sum_k N_k(\mathrm{d}z)$ 為盒內總事件數：
#
# $$\begin{aligned}
# P\bigl(N_k(\mathrm{d}z) = 1 \bigm| H_{t_j}\bigr)
#   &= \lambda_k(z_j)\,\mathrm{d}z + o(\mathrm{d}z), \\
# P\bigl(N(\mathrm{d}z) = 1 \bigm| H_{t_j}\bigr)
#   &= \lambda^*(z_j)\,\mathrm{d}z + o(\mathrm{d}z) .
# \end{aligned}$$
#
# 相除、令 $\mathrm{d}z\to0$（兩個成分同時落在同一小盒的機率是
# $O(\mathrm{d}z^2)$，消失）：
#
# $$
# \rho_{ij} = \frac{\kappa(m_i)\,g(t_j-t_i)\,
#   f(x_j-x_i,\,y_j-y_i;\,m_i)}{\lambda^*(t_j,x_j,y_j)},
# \qquad
# \phi_j = \frac{\mu(x_j,y_j)}{\lambda^*(t_j,x_j,y_j)}
# $$ (eq:rho-phi)
#
# $\rho_{ij}$ 讀作「事件 $j$ 被事件 $i$ 觸發的機率」，$\phi_j$ 讀作
# 「事件 $j$ 是背景事件的機率」。由 $\lambda^*$ 的定義，兩者自動滿足
# $\phi_j + \sum_{i:\,t_i<t_j}\rho_{ij}=1$。這條恆等式不是額外假設，
# 是分母的定義使然——**強度分解到哪裡，機率就分配到哪裡**。整個
# 隨機除叢的數學內容就這麼多：一個比例。
#
# ### 變頻寬核：背景率的無母數估計
#
# 有了 $\phi_j$，背景率就能用**加權**核密度估計反推——每個事件只以
# $\phi_j$ 的權重貢獻給背景場：
#
# $$\hat u(x,y) = \frac{1}{T}\sum_{j=1}^{N}
#   \phi_j \; Z_{h_j}(x-x_j,\, y-y_j),
#   \qquad
#   Z_h(\Delta x,\Delta y)=\frac{1}{2\pi h^2}
#     \exp\!\left[-\frac{\Delta x^2+\Delta y^2}{2h^2}\right] .$$
#
# 頻寬不是常數，而是**隨事件而變**：$h_j = \max\{h_{\min},\,r(j,n_p)\}$，
# 其中 $r(j,n_p)$ 是事件 $j$ 到第 $n_p$ 個最近鄰的距離。地震密的地方
# 頻寬小、解析度高；疏的地方頻寬自動放大，不會出現一堆孤立的針尖。
# $h_{\min}$ 防止重疊事件造成零頻寬，取值應是**定位誤差的量級**
# ——Jalilian 預設 $0.05^\circ\approx5.6$ km 配 $n_p=5$，CWA 112 年
# 報告取 3 km（對應台灣目錄 2–3 km 的定位誤差）配 $n_p=5$。這有物理
# 意義：**背景率場的解析度不該比定位精度更細**，否則畫出來的細節是
# 定位誤差的花紋。
#
# ### 迭代演算法
#
# 兩件事互相依賴，就交替更新（Zhuang et al. 2002；Jalilian 2019 的
# Algorithm 2）：**（一）初始化** $\theta$ 與 $\phi_j$（例如全設
# 0.5）；**（二）更新背景場**，用當前 $\phi_j$ 與變頻寬核算
# $\hat u(x,y)$；**（三）更新參數**，固定 $\hat u$ 用 DFP 最大化
# $l_2(\theta)$；**（四）更新機率**，用新的 $\hat\theta$ 與 $\hat u$
# 代進 {eq}`eq:rho-phi` 重算 $\phi_j$；**（五）回到第二步**，直到
# 兩者都不再動。收斂通常很快：伊朗 4 次、大埔 6 次、CWA 112 年報告
# 3–7 次。這是 EM 演算法的近親（$\phi_j$ 是隱藏標籤的後驗機率），
# 每步都不會讓概似變差，但**不保證收到全域最大**。
#
# 把總速率也做同樣分解，就得到**叢集係數**：
#
# $$\Lambda(x,y) \approx \hat u(x,y) + \frac{1}{T}\sum_{i}
#   \kappa(m_i)\,f(x-x_i,y-y_i;m_i),
#   \qquad
#   \omega(x,y) = 1 - \frac{\hat u(x,y)}{\Lambda(x,y)} .$$
#
# $\omega$ 接近 1 表示該處活動幾乎全是被觸發的，接近 0 表示幾乎全是
# 背景。它是 $\phi_j$ 的場版本，可以直接畫成地圖，用來分辨「這個熱區
# 是真的背景活躍，還是某次大震的餘波」。
#
# ### 這就是第 12 章那個難題的解法
#
# 回到第 12 章的破產點 $m_x$（{eq}`eq:mx`）。傳統除叢法的病根是
# **二分**：每個事件非主震即餘震，而「主震 = 叢集中最大事件」這條
# 規則會系統性壓低 b 值。隨機除叢在三個層面上不一樣：
#
# - **它給機率，不給標籤。** $\phi_j=0.3$ 就是 0.3，不會被四捨五入
#   成 0 或 1。
# - **它承認的二分是「背景 vs 被觸發」，不是「主震 vs 餘震」。**
#   ETAS 的後代可以比親代大。Mizrahi et al.（2021）的關鍵對照就在
#   這裡：以「未被觸發」定義主震時，b 值與全目錄**無顯著差異**；只有
#   套上「叢中最大」那條規則，30% 的 b 值偏差才冒出來。
# - **它可以把不確定性顯示出來。** 以 $\phi_j$ 為機率做 thinning 得到
#   一份背景子目錄，換個種子得到另一份；對一百份子目錄各算一次
#   b 值，散布就是**除叢不確定性的直接量測**。
#
# 代價要講清楚：$\phi_j$ **是模型的產物**，換一組參數或換一個空間核
# 就會變。隨機除叢沒有消滅主觀性，只是把它從「選哪個視窗」搬到
# 「選哪個模型」——差別是後者有概似值可以比較。
#
# 拿 2024 年春季的台灣目錄實作一次。用 CWA 112 年報告那七個觸發參數，
# 背景場則用上面的迭代法從資料本身估：

# %% tags=["hide-input"]
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
# 深色是高背景機率、淺色是高觸發機率。花蓮外海那一整片幾乎全是最淺
# 的顏色——它們是 0403 序列的後代，$\phi$ 接近 0；西部與南部零星的
# 深色點則是這段期間的背景活動。整份目錄 $\sum_j\phi_j$ 只有一百多個
# 背景事件，**其餘九成五都是被觸發的**。
#
# 三個但書。**其一**，這裡只用 121 天的窗，背景場是從這個窗自己估
# 出來的；正式作業流程用的是數十年的目錄（CWA 用 1994–2021），這裡
# 的數字只能當演示。**其二**，距離全部以「度」計，沒有做緯度校正，
# 也沒有考慮台灣餘震區明顯的**非等向性**——各向同性空間核是標準
# ETAS 的已知弱點（Ogata & Zhuang 2006 用一個 $2\times2$ 矩陣處理
# 橢圓形餘震區）。**其三**，也最有趣：$M_L\,7.19$ 的 0403 主震本身
# $\phi\approx0.98$，模型判它是**背景事件**。ETAS 不需要「主震」這個
# 詞——它只說這一顆沒有被誰觸發。

# %% [markdown]
# ## 14.5 診斷：模型錯在哪裡
#
# 參數估好了，下一個問題是模型對不對。點過程有一整套現成工具，全部
# 建立在第 10 章的隨機時間變換定理上。
#
# ### $\tau$–$j$ 圖：偵測相對寧靜的工具
#
# 定義**變換時間** $\tau_j=\int_{t_s}^{t_j}
# \lambda^{\rm temp}_{\hat\theta}(t\mid H_t)\,\mathrm{d}t$，其中
# $\lambda^{\rm temp}$ 是把空間積掉後的時間強度。若模型正確，
# $\{\tau_j\}$ 是單位速率 Poisson 過程，$\tau_j$ 對 $j$ 作圖應落在
# $y=x$ 附近。判讀規則兩條，而且都要加上「相對」兩個字：**低於直線
# 的區段是相對於模型的寧靜（quiescence）**，實際事件比模型預期的少；
# **高於直線是相對活化**。
#
# 「相對」不是修辭。它不是絕對的活動高低，而是**扣掉模型（尤其是
# 餘震正常衰減）之後剩下的部分**。這正是 Ogata & Zhuang（2006）造
# 時空 ETAS 的原始動機：Lomnitz & Nava（1983）批評說，大震前所謂的
# 「地震寧靜」多半只是前一場序列的餘震正常衰減造成的錯覺；要反駁
# 這個批評，就必須先有一個能描述**正常叢集**的參考模型。**ETAS
# 因此是一個前兆偵測工具——不是因為它能預測，而是因為它能把「正常」
# 定義清楚。**
#
# ### 三種殘差
#
# 更一般的工具是**加權殘差**。取時空區塊 $I\times B$ 與權重 $h$：
#
# $$R(I\times B; h) = \sum_i \delta_i\,
#   \mathbf{1}[\,t_i\in I,\ (x_i,y_i)\in B\,]\;h(t_i,x_i,y_i)
#   \;-\; \iiint_{I\times B} h\,\lambda_{\hat\theta}
#   \,\mathrm{d}x\,\mathrm{d}y\,\mathrm{d}t$$
#
# | 權重 $h$ | 名稱 | 特性 |
# |---|---|---|
# | $1$ | raw residual | 就是「觀測數 − 期望數」，被高強度區主導 |
# | $1/\lambda$ | reciprocal residual | 放大低強度區的偏差，對背景率誤設敏感 |
# | $1/\sqrt{\lambda}$ | Pearson residual | 變異數穩定化，兩端折衷，最常用 |
#
# 模型正確時三者都應在 0 附近沒有系統性偏離。實務上把 Pearson 殘差
# 畫成空間地圖最有用：一片系統性為正的區域代表模型在那裡**低估**了
# 活動——通常是背景率場沒把某個構造抓出來。
#
# 由 {eq}`eq:time-rescale`，$U_j = 1-e^{-(\tau_j-\tau_{j-1})}$ 應該是
# i.i.d. 的 $U(0,1)$，於是可以畫 Q–Q 圖、跑 Kolmogorov–Smirnov 檢定。
# 一個重要但書：**參數是從同一份資料估出來的，KS 的 p 值因此過於
# 樂觀**（標準 KS 假設分布完全指定）。嚴格做法是用參數 bootstrap
# 校正臨界值，或把資料切成訓練段與檢驗段。

# %% tags=["hide-input"]
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
# 左圖幾乎與 $y=x$ 重合，KS 檢定完全過關。右圖是同一份資料換成
# 「忽略叢集」的齊次 Poisson：曲線變成一道**樓梯**——長長的平段是
# 模型看不見的相對寧靜（實際上是正常的背景期，但 Poisson 把餘震的
# 量攤平到全部時間，所以在這裡預期太多），陡直的跳段是它完全沒預料
# 到的餘震爆發。KS 的 p 值差了三十個數量級。
#
# 這張圖同時示範了一件事：**$\tau$–$j$ 圖上的「寧靜」永遠是相對於
# 某個模型的**。同一段時間，用 Poisson 當參考模型看起來寧靜得不得了，
# 用正確的 ETAS 看就完全正常。任何宣稱偵測到「震前寧靜」的研究，
# 第一個該被問的問題是：**你的參考模型是什麼？**

# %% [markdown]
# ## 14.6 不完整性如何污染參數
#
# ETAS 的估計假設目錄**完整且時間穩態**。第 11 章已經證明兩件事都不
# 成立：$M_c$ 是隨時間、空間、尤其是**主震後時刻**變動的場。大震後
# 幾小時內小地震被尾波淹沒、測站飽和、人工檢視來不及，目錄會缺掉
# 一大塊——而這一塊剛好是模型最看重的一塊。
#
# **一、$\alpha$ 被壓低**，這是最嚴重也最不直觀的一個。缺漏最嚴重的
# 時刻是**大事件剛過**的時候，而那正是大事件展示產能的時刻。模型
# 看到的是：規模 6 的事件在頭一天只生了幾個看得見的後代，規模 3.5 的
# 事件在自己那個安靜的角落生了兩三個——**大事件的相對產能被系統性
# 低估**，$\hat\alpha$ 因此往下掉。Hainzl et al.（2013）指出，同時
# 處理不完整性與時變背景率之後，$\hat\alpha$ 會回升到接近 $\beta$。
# 這正是 simplETAS 敢把 $\alpha$ 直接釘在 $\beta$ 的依據之一。
#
# **二、$c$ 被高估。** 缺漏在觀測率上造出一段假的平緩期，而它的形狀
# 與 $c$ 的效果**無法區分**（12.2 節推導過）。測站愈密、用模板比對
# 補完目錄，$\hat c$ 就愈小——這是 $c$ 主要是資料人工品的最強證據。
#
# **三、$p$ 被扭曲，方向不定。** 只缺頭幾小時會讓早期率被壓低、使
# 衰減看起來比較平（$p$ 偏小）；若同時把擬合起點往後挪又會反過來。
# 這就是為什麼報告 $p$ 值**必須**同時報告擬合的起訖時間。
#
# ### $p<1$ 多半是背景率誤設的病徵
#
# 還有一個與不完整性無關的常見病因。Ogata & Zhuang（2006）在假設
# **均勻背景率**的那張表裡擬出 $\hat p<1$，作者直接判定：這代表均勻
# 背景率的假設不成立。機制很清楚——空間上不均勻的背景活動沒有地方
# 去，只好被 Omori 核吸收；核為了同時裝下「近期的真餘震」與「遠期
# 均勻散布的背景事件」，只能把尾巴拉平，也就是把 $p$ 壓到 1 以下。
# 而第 12 章證明過，$p\le1$ 時 $g(t)$ 根本不可正規化，
# {eq}`eq:omori-density` 失效，分支比公式也一起失效。
#
# > **$\hat p<1$ 幾乎從來不是「這個序列衰減特別慢」的發現，而是
# > 「模型設定錯了」的診斷訊號。** 第一個該檢查的是背景率——是不是
# > 用了常數 $\mu$？研究區是不是太大、把兩個構造區混在一起？第二個
# > 該檢查的是有沒有一個疊上來的次級序列被硬塞進同一個核。
#
# 用合成實驗把「$\alpha$ 被壓低」量出來。做法很乾淨：先用已知的
# $\alpha=1.60$ 模擬一份完整目錄，再模仿觀測網的失能——**每一個
# $M\ge5.0$ 事件之後的一天內，刪掉所有小於某個門檻的事件**（刪掉之後
# 分析者再也看不到它們，包括在歷史 $H_t$ 裡）——然後用同一套 MLE
# 重估 $(\mu, A, \alpha)$：

# %% tags=["hide-input"]
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
# 灰柱說明這個實驗有多溫和：刪掉一個規模單位時總事件數只少了一成
# 出頭，而且只發生在少數幾個大事件後的頭一天。但 $\hat\alpha$ 已經
# 從真值 1.60 附近掉到 1.0 上下，**低估了三分之一以上**。橘色點線是
# 同時發生的補償：$\hat A$ 反向上升——這正是 14.3 節那條山脊在作用，
# 資料堅持總觸發量不變，模型只能把「規模加成」換成「基準產能」。
#
# 所以，如果一份目錄的 $\hat\alpha$ 明顯低於 $\beta$，在下結論說
# 「這個區域的大地震觸發效率特別差」之前，先問有沒有做不完整性校正。
# 這也是 simplETAS 直接把 $\alpha$ 釘在 $\beta$ 的理由——與其估一個
# **已知會被系統性壓低**的參數，不如用理論值。

# %% [markdown]
# ## 14.7 Reasenberg–Jones 是 ETAS 的解析特例
#
# 1989 年，Reasenberg 與 Jones 在 *Science* 上發表了一條式子，成為
# 之後三十多年美國餘震預報的骨幹：主震規模 $M_m$ 之後 $t$ 時刻、
# 規模 $M$ 以上的餘震發生率
#
# $$
# R(t, M) = 10^{\,a' + b\,(M_m - M)}\,(t + c)^{-p}
# $$ (eq:rj-rate)
#
# 它看起來與 ETAS 沒什麼關係——沒有背景率、沒有疊加、沒有二次觸發、
# 參數只有四個。但它其實是 ETAS 在兩個假設下的**精確解析特例**。
#
# ### 推導
#
# **假設一：只看一個觸發源。** 忽略背景率與二次觸發，只留主震
# $M_m$ 這一項：
#
# $$\lambda^*(t, m) = \kappa(M_m)\,g(t)\,s(m)
#   = A\,e^{\alpha(M_m-m_0)}\;\frac{(p-1)c^{p-1}}{(t+c)^{p}}\;
#     \beta e^{-\beta(m-m_0)} .$$
#
# **對規模積分**，得到「規模 $\ge M$」的發生率：
#
# $$\begin{aligned}
# R(t,M) &= \int_M^{\infty} \lambda^*(t,m)\,\mathrm{d}m \\
#   &= A\,e^{\alpha(M_m-m_0)}\,\frac{(p-1)c^{p-1}}{(t+c)^{p}}
#      \int_M^{\infty}\beta e^{-\beta(m-m_0)}\,\mathrm{d}m \\
#   &= A\,(p-1)\,c^{\,p-1}\;
#      e^{\alpha(M_m-m_0)}\,e^{-\beta(M-m_0)}\;(t+c)^{-p} .
# \end{aligned}$$
#
# **假設二：$\alpha=\beta$**（自相似）。兩個指數立刻合併，$m_0$ 完全
# 消失：
#
# $$e^{\alpha(M_m-m_0)}e^{-\beta(M-m_0)}
#   \;\xrightarrow{\ \alpha=\beta\ }\;
#   e^{\beta(M_m - M)} = 10^{\,b\,(M_m - M)} ,$$
#
# 最後一步用的是本書的約定 $\beta=b\ln10$。於是
#
# $$R(t,M) = \underbrace{A\,(p-1)\,c^{\,p-1}}_{\displaystyle =\,10^{a'}}
#   \;10^{\,b(M_m-M)}\,(t+c)^{-p} ,$$
#
# 這就是 {eq}`eq:rj-rate`，並順帶給出兩套參數的精確對照
# $a' = \log_{10}[A(p-1)c^{\,p-1}]$。
#
# 三個結論值得記住。**其一，R–J 的 $b$ 同時扮演兩個角色**——既是
# 目標規模的 GR 斜率，也是主震規模的產能指數，因為 $\alpha=\beta$ 把
# 兩者黏在一起了。**其二，$M_m-M$ 這個差值形式不是巧合**，正是
# $\alpha=\beta$ 讓 $m_0$ 消掉的結果。**其三，R–J 沒有背景率**，
# 只適用於「主震剛過、餘震壓倒背景」的窗口；時間拉長必須另外加回
# 背景項（STEP 與紐西蘭的混成模型都這麼做）。
#
# 由第 10 章，$[t_1,t_2]$ 內至少發生一次的機率是
#
# $$P = 1 - \exp\!\left[-\int_{t_1}^{t_2} R(t,M)\,\mathrm{d}t\right],
#   \qquad
#   \int_{t_1}^{t_2}(t+c)^{-p}\,\mathrm{d}t
#   = \frac{(t_1+c)^{1-p}-(t_2+c)^{1-p}}{p-1} .$$
#
# 這裡偷偷用了非齊次 Poisson 的假設——把率當成外生的，忽略餘震彼此
# 的觸發。這是假設一的代價，也是 ETAS 相對於 R–J 的主要增益：ETAS
# 用蒙地卡羅模擬把級聯算進去，R–J 只能算期望率。
#
# Reasenberg & Jones 分析 **62 個加州序列**得到通用參數
# $a'=-1.67$、$b=0.91$、$p=1.08$、$c=0.05$ 天。把 $M=M_m$ 代進去
# （問「一個地震在一週內被同規模或更大的地震跟隨」的機率），得到
# **10.5%**；而 Jones（1985）統計南加州的實際前震比例是
# **$6.0\pm0.5\%$**。換上 Utsu（1969）的日本參數，同樣算法得到文獻
# 報的 **4.2%**（用這組已四捨五入的參數重算是 4.3%，差別來自捨入）。
#
# | 來源 | 參數 | 一週內被同等或更大事件跟隨 |
# |---|---|---|
# | R–J 加州通用（模型） | $a'=-1.67, b=0.91, p=1.08, c=0.05$ | 10.5% |
# | 南加州實測前震率 | Jones (1985) 的統計 | $6.0\pm0.5\%$ |
# | Utsu 日本參數（模型） | $a'=-1.83, b=0.85, p=1.3, c=0.3$ | 4.2% |
#
# 這張表值得盯著看一分鐘。**同一個模型換一組區域參數，答案差一倍
# 以上；模型與觀測之間又差一倍。** 參數在地化不是講究，是必要；而
# 「模型算出來的機率」與「歷史上的頻率」是兩個不同的東西，兩者的
# 差距本身就是研究題目（第 18 章的可靠度檢驗）。

# %% tags=["hide-input"]
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
# 三條線形狀一模一樣，只是上下平移——因為 $M_m$ 只透過
# $10^{b(M_m-M)}$ 這個**乘法因子**進入 {eq}`eq:rj-rate`，主震規模改變
# 的是機率的高度，不是衰減的形狀。M 7.2 主震剛過時，一週內再來一個
# M≥5 幾乎是確定的；一個月後掉到幾成；一年後回到個位數百分比。灰
# 虛線標出上面那兩個對照數字，**位置很低**——「下一個更大的地震」
# 始終是低機率事件，即使在序列最猛的時候。這是短期預報的宿命：機率
# 增益（相對背景可高上百倍）很大，絕對機率從來不高。
#
# 順帶澄清兩件常被弄混的事。**其一**，2005 年上線的 STEP
# （Gerstenberger et al. 2005）已經走入歷史，美國現行系統是**作業化
# 餘震預報（OAF）**，引擎是 R–J 或 ETAS，參數有 generic、
# sequence-specific 與貝氏結合三種來源；把 STEP 講成「USGS 現在在跑
# 的系統」是常見的錯誤。**其二**，R–J 的 $a'$ 是**dressed**（已含
# 級聯）的產能，ETAS 的 $A$ 是 **bare**（只算直接後代）的產能，兩者
# 差一個 $1/(1-n)$ 因子——附錄 D 用加州參數把這個換算做一遍，順便
# 回推出一個非常合理的分支比。

# %% [markdown]
# ## 14.8 參數精簡：simplETAS 的逐項理由
#
# ETAS 二十年來的主流方向是**加**：斷層幾何、深度分布、時空變動
# 參數、貝氏即時更新。Mancini & Marzocchi（2023）反其道而行：
# **如果把描述叢集的參數全部釘死，只留最少的自由度，還能用嗎？**
# 他們釘住七個參數，只估兩個明顯與區域有關的——總背景率與基準
# 產能 $A$。
#
# | 參數 | 釘死值 | 理由 |
# |---|---|---|
# | $\alpha$ | $\ln 10\,(=\beta)$ | 維持地震**自相似性**；可重現 Båth 定律；避免各向同性空間核造成的產能偏差（Helmstetter 2005；Hainzl 2008）；同時處理目錄不完整與時變背景率後，$\alpha$ 本來就接近 $\beta$（Hainzl 2013） |
# | $p$ | 1.15 | Utsu（1995）經驗範圍 $[0.9,1.4]$ 的中位數；加州全州目錄的期望範圍 1.0–1.4 |
# | $c$ | 0.005 天 | 經驗範圍 0.003–0.3 天，但 $c$ 會被早期不完整系統性高估、隨截切規模浮動；本文測試顯示它對時間表現的影響可忽略 |
# | $D$ | 1 km² | 全球估計橫跨數個數量級（加州 $<0.1$、隱沒帶 $>20$ km²）；取中間量級，也大致相當於現代目錄的平均水平定位誤差 |
# | $q$ | 1.5 | 與**靜態應力隨距離三次方衰減**一致 |
# | $\gamma$ | 1.5 | 資料豐富地區的估計典型落在 1.0–2.0 |
# | $\beta$ | $\ln 10$ | GR 律取 $b=1$ |
#
# 逐行讀會發現一個模式：**每一個被釘死的參數，都是前面幾節證明過
# 「很難估準」的那些**。$c$ 被早期不完整污染（14.6）、$D$ 與 $\gamma$
# 只有乘積被辨識（14.3）、$\alpha$ 被系統性壓低（14.6）。simplETAS
# 的邏輯不是「這些參數不重要」，而是**它們的估計噪音大於它們帶來的
# 擬合改善**。
#
# 這正是 **bias–variance 取捨**換上地震學的外衣。多留一個自由參數，
# 模型的**偏差**下降——它能貼合這個區域的特性；但**變異**上升——
# 那個參數的估計值帶著雜訊，而雜訊會被完整地帶進預報。當某個參數的
# 標準誤跟估計值同一個量級（回頭看伊朗表的 $c$、$D$、$\gamma$），
# 釘死一個「還算合理的常數」在均方誤差上幾乎必然更好。義大利的校正
# 結果只有兩個數字：$\nu = 18.27$ 次/年、$A=0.047$（$M_{\min}=3.95$）。
#
# 光是參數少不能證明模型好，得檢驗。他們的四尺度架構值得整套模仿：
# **（一）50 年、$M\ge3.95$**（HORUS 儀器目錄，深度 $\le30$ km）：
# 觀測累積曲線對上 10,000 次模擬的分布。**（二）392 年、$M\ge5.95$**
# （CPTI 歷史目錄，1630–2021）：**完全 out-of-sample**——參數只用近
# 50 年儀器資料估，卻要預測四個世紀的歷史地震數，這是全篇最強的
# 一項檢驗。**（三）序列日預報**（Amatrice–Visso–Norcia 2016）：
# pseudo-prospective，用增量 N-test 逐日評分。**（四）敏感度測試**：
# 擾動被釘住的七個參數，看模擬目錄的事件數與位置變化多大。跨越
# 「日」到「世紀」四個數量級全部通過——這是很強的主張，因為短期與
# 長期考驗模型的不同部位（短期考驗觸發核，長期考驗背景率與總量
# 守恆）。他們也用了 **Turing-style test**（Page & van der Elst
# 2018）：把合成目錄與真實目錄擺在一起，看專家能不能分辨。
#
# 作者列的三個用途都對後面章節有意義：**benchmark**（任何新模型都該
# 先贏過它，CSEP 跨區域比較才有共同基準）、**OEF**（目錄短的地區也
# 能做作業型預報）、**PSHA**（直接模擬任意長度的合成目錄，繞開除叢
# 與 Poisson 假設，正好接上第 12 章對除叢的批評）。
#
# 最後一個必須點破的數學細節。$\alpha=\beta$ 時，分支比的積分
# $n=A\beta\int_{m_0}^{\infty}e^{(\alpha-\beta)(m-m_0)}\,\mathrm{d}m$
# 的被積函數變成常數 1，積分**發散**——第 13 章那條
# $n=A\beta/(\beta-\alpha)$ 只在 $\alpha<\beta$ 時成立。所以任何
# $\alpha=\beta$ 的自相似參數化，**必須搭配一個有限的最大規模
# $M_{\max}$**，此時 $n=A\beta\,(M_{\max}-m_0)$ 才有限。這不是
# simplETAS 的瑕疵，是所有自相似設定的標準配件——但它提醒我們，
# $\alpha=\beta$ 是一條數學上的懸崖邊。

# %% [markdown]
# ## 14.9 台灣的 ETAS
#
# ### 兩組互相獨立的八參數估計
#
# 台灣已經有本土化的時空 ETAS 參數，而且有兩份**互相獨立**的估計：
# 中央氣象署 112 年（2023）委辦報告子計畫三（詹忠翰、謝銘哲、
# 呂奇祝），訓練窗 1994–2021；以及 2025 年大埔序列快報
# （Hsieh et al. 2025），訓練窗 1994–2024/12。兩者都取
# $M_c=M_L\,3.0$、都用 Cheng et al.（2015）的淺層地殼分區、都用
# 119.8–122.5°E × 21.7–25.6°N 的 81×201 格點：
#
# | 參數 | CWA 112 年（2023） | 2025 大埔快報 | 註記 |
# |---|---|---|---|
# | $\mu$ | 0.5098 | 0.5424 | 兩份的單位敘述不同，見下 |
# | $A$ | 0.6188 | 0.9166 | 基準產能 |
# | $c$ | 0.0031 天 | 0.0012 天 | 都遠小於 simplETAS 的 0.005 |
# | $\alpha$ | 1.1733 | 1.0408 | 都明顯低於 $\beta\approx2.30$ |
# | $p$ | 1.0616 | 1.0350 | 一致到小數第二位 |
# | $D$ | $5\times10^{-5}$ deg² | 原文寫 $D^2=0.0007$ deg² | **記號不同，不可直接比** |
# | $\gamma$ | 0.6786 | 0.3511 | 差近一倍 |
# | $q$ | 1.5934 | 2.4253 | 差很多 |
#
# 這張表本身就是 14.3 節的活教材。**先看一致的部分**：$p$ 差 0.03、
# $\alpha$ 差 0.13、$\mu$ 差 0.03。兩個團隊、兩段訓練窗、兩套實作，
# 卻收斂到幾乎同一組數字——這是很強的交叉驗證。**再看不一致的部分**：
# $\gamma$ 差近一倍、$q$ 差 0.83、$D$ 的記號根本不同（一份寫 $D$、
# 一份寫 $D^2$）。完全符合 14.3 節的預測——**被資料釘死的是 $p$ 與
# $\alpha$，鬆脫的是空間核那三個參數**。兩份報告都沒有附標準誤，
# 我們只能從「兩份獨立估計的離散度」反推不確定性；這是一種粗糙但
# 有效的 bootstrap。
#
# 兩個陷阱。**其一**，$\mu$ 在不同實作裡意義不同：Jalilian 的參數化
# 裡 $\tilde u(x,y)=\mu\,u(x,y)$，$\mu$ 是加速收斂用的鬆弛係數，不是
# 事件率密度。看到 $\mu=0.51$ 之前，先確認它是「次/天/deg²」還是一個
# 無量綱係數。**其二**，把這兩組 $(A,\alpha)$ 連同 $b=1$
# （$\beta=\ln10$）代進第 13 章的分支比公式，兩組都會得到 $n>1$。
# 這**不代表台灣的地殼超臨界**——更可能的解釋是慣例落差：$A$ 搭配的
# Omori 核有沒有正規化、$m_0$ 取的是 $M_c$ 還是 $M_c-\Delta M/2$、
# 台灣的 $b$ 值實際上是多少。這是很好的學生練習：算一遍，然後列出
# 所有可能讓 $n$ 回到 1 以下的慣例差異。**跨文獻搬參數之前先算一次
# 分支比，是最便宜的健康檢查。**
#
# ### 作業化的工程數字
#
# **CWA 112 年**：384 核約 52 分鐘完成完整參數迴歸；若沿用已訓練
# 參數、只做隨機除叢與速率計算，可降到 96 核。**2025 大埔**：128 核、
# 6 次 MLE 迭代約 **1 小時**完成訓練；沿用預訓練參數後，**每個即時
# 預報時窗只需 8 核、7 分鐘**。
#
# 這組對比是整章最實用的一段。**作業化的瓶頸在預訓練，不在即時
# 運算**——這解釋了為什麼所有作業系統的架構長得一樣：先花大錢把
# 八個參數訓練好、鎖起來，地震來了只更新目錄與背景率，跑輕量的
# 模擬。大埔的完整流程是：每小時更新 → 估 $\lambda$ 與 $\mu$ →
# 產生 **1,000 份合成目錄** → 用 GMM（Lin et al. 2011，含 $V_{S30}$
# 與孕震深度）算各站 PGA → 轉震度 → 統計超越機率。從目錄一路走到
# 場址震度機率，這是台灣目前最完整的端到端統計預報範例。
#
# ### 預報 vs 觀測
#
# 2025/01/20 大埔 $M_L\,6.4$（深度 15.8 km，大埔震度 6 弱）之後，
# 以震後 17:00 UTC 為起算的預報，對上實際發生次數：
#
# | 時窗 | $P(M_L\ge5)$ | $P(M_L\ge6)$ | 實際 $\ge5$ 次數 |
# |---|---|---|---|
# | 1 天 | 30.5% | 3.4% | 1 |
# | 3 天 | 42.1% | 6.2% | 1 |
# | 7 天 | 58.8% | 8.1% | 7 |
# | 10 天 | 67.8% | 11.2% | 8 |
#
# 怎麼讀這張表？**不能用單次結果判對錯**。10 天內 $P(M\ge6)=11.2\%$
# 而實際沒發生，這不算預報失敗——11.2% 的事件本來就有 88.8% 的機會
# 不發生。要判斷一組機率預報準不準，必須用**一整批**預報做一致性
# 檢定（N-test、S-test 等，第 18 章）。
#
# 最後三筆誠實面。大埔序列**震前沒有前兆**：源區沒有明顯高活動率，
# 即時目錄未顯示前驅訊號。序列呈四階段結構（主震 → 1/25
# $M_L\,5.7$ → 1/30 $M_L\,5.6$ → 之後），每隔約 5–6 天出現一次群集。
# 1/24 到 1/25 之間有一段「在 Omori 框架下不尋常的平靜」，但作者
# 謹慎地說**即時目錄的完整度與定位精度不足以判定那是不是真平靜**。
# 這句話點出一個教科書很少提的落差：**作業型預報用的是即時目錄，
# 研究用的是重定位目錄，兩者的統計性質不同。**

# %% [markdown]
# ## 14.10 常見誤解與陷阱
#
# **一、「參數估出來就是物理值。」** 看 14.3 那張伊朗表的
# $\hat\gamma = 0.0238 \pm 5.7553$。ETAS 參數彼此高度相關
# （$\gamma$–$D$ 只有乘積被辨識、$\alpha$–$A$ 同一個結構、$c$–$p$
# 只靠最早那幾個事件分辨），**一定要看標準誤，最好看聯合信賴區域**。
#
# **二、「$c$ 值是物理量。」** $c$ 主要反映早期餘震目錄不完整，隨
# 截切規模與測站密度浮動，補完目錄之後就縮小。simplETAS 乾脆把它釘
# 在 0.005 天。
#
# **三、「早期不完整只會影響 $c$。」** 不是。14.6 節的合成實驗顯示，
# 只刪掉一成出頭的早期小事件，$\hat\alpha$ 就從 1.60 掉到 1.0 附近。
# 缺漏最嚴重的時刻正是大事件展示產能的時刻，所以受害最深的是
# $\alpha$。
#
# **四、「$p<1$ 代表餘震衰減特別慢。」** 多半是模型設定錯了：空間
# 不均勻的背景率被硬塞進 Omori 核。而且 $p\le1$ 時
# {eq}`eq:omori-density` 根本不可正規化，分支比也一起失效。
#
# **五、「參數越多、模型越複雜就越好。」** simplETAS 是最好的反例：
# 釘死七個參數，在義大利從「日」到「四個世紀」四個尺度全部通過。
# 複雜模型必須先證明自己贏過這個基準。
#
# **六、「隨機除叢消滅了除叢的主觀性。」** 沒有。它把主觀性從「選
# 哪個視窗」搬到「選哪個 ETAS 模型與參數」——差別在於後者有概似值、
# AIC 與殘差圖可以比較。
#
# **七、「ETAS 的餘震一定比主震小。」** ETAS 對後代規模的唯一約束是
# 「從 GR 律獨立抽出」，時間上必須在親代之後，規模上**沒有上限**。
# 這正是它能自然重現前震現象的原因，也是它與傳統主震／餘震框架的
# 根本分歧。
#
# **八、「$\ln L$ 或 AIC 可以跨研究比較。」** 不行。只有在同一份
# 目錄、同一個 $m_0$、同一個研究區 $S$ 與同一段時間窗之下，$\ln L$
# 的差值才有意義。把 $S$ 縮到只包含事件密集區，$\ln L$ 會虛假地變
# 好看（10.2 節）。
#
# **九、「收斂了就代表估對了。」** DFP 收斂只代表梯度接近零。
# $\hat c$ 貼下界、$\hat p$ 貼 1、$\hat q$ 貼 1，最佳化一樣會回報
# 成功。**參數貼邊界是模型設定有問題的訊號**，不是估計成功的訊號。

# %% [markdown]
# ## 14.11 研究前沿：有人贏過 ETAS 了嗎
#
# ### 深度學習的第一回合：一場公開的挫敗
#
# DeVries et al.（2018）在 *Nature* 上發表用深度神經網路預測餘震
# **空間分布**的研究，一度被當成「機器學習終於攻進地震學」的標誌。
# 隔年 Mignan & Broccardo（2019）提出反駁：他們證明**一個遠更簡單的
# 參數化模型就能達到相當的資訊量**，神經網路的優勢主要來自檢驗設計
# 而非模型能力。這一回合讓社群卻步了幾年，也留下一條至今適用的
# 守則：**任何宣稱贏過既有模型的新方法，必須先與一個誠實的、調校過
# 的簡單基準比較**——不是與稻草人比較。
#
# ### 風向轉變：FERN 與神經點過程
#
# 近幾年情勢變了。Google 的 **FERN**（Zlydenko et al. 2023）在大規模
# 目錄上**略優於 ETAS**；**神經點過程**（Stockman et al. 2023；
# Dascher-Cousineau et al. 2023）在訓練資料充足時穩定優於 ETAS，
# 訓練還更快。
#
# 神經點過程的想法與第 10 章的框架完全相容：**保留
# {eq}`eq:pp-loglik` 這個概似骨架，只把 $\lambda^*$ 的函數形式換掉。**
# ETAS 把 $\lambda^*$ 寫成「背景 + 固定形狀的觸發核疊加」，神經點
# 過程改用 RNN 或 Transformer 把歷史 $H_t$ 編碼成隱藏狀態，再由它
# 輸出強度；訓練目標仍是最大化 $\sum\ln\lambda^*-\int\lambda^*$，
# 只是參數從八個變成幾十萬個。
#
# 它的優勢正好對應 ETAS 的三個結構性弱點：**不必假設觸發核的函數
# 形式**（Omori 的冪次、空間核的等向性、$\kappa$ 的指數形式全是可
# 協商的經驗選擇）；**對目錄的非平穩性極為適應**（ETAS 假設目錄完整
# 且穩態，而真實目錄的 $M_c$、測站密度、規模尺度都在變）；**能利用
# 機器學習產生的高解析度目錄**，即使那些目錄本身不完整。代價是可
# 解釋性沒了，而且神經模型要大量訓練資料，這對目錄短的地區是硬
# 約束——這正是 simplETAS 那條「參數少、可攜性高」路線的存在理由。
# 兩者不是誰取代誰，而是分別佔據資料光譜的兩端。
#
# ### 社群缺一個標準化的 ETAS 基準
#
# 最後一個前沿問題是工程問題，卻可能最重要：**沒有兩個 ETAS 實作是
# 相同的。** 邊界效應怎麼處理、complementary event 怎麼定、空間積分
# 怎麼近似、背景率用什麼頻寬、迭代收斂條件是什麼——每一項都有幾種
# 合理選擇，合起來足以讓兩個「同樣是 ETAS」的模型給出明顯不同的
# 預報。CSEP 第二階段因此把重心轉向開放軟體工具箱與可重現性套件，
# 並把「**建立一個真正標準化的 ETAS 基準版本**」點名為當前最值得
# 投資的社群工程。
#
# 這對「有人贏過 ETAS 了嗎」是致命的：**如果沒有標準的 ETAS，
# 「贏過 ETAS」這句話就沒有明確意義。** 現行所有作業化模型還有一個
# 共同限制值得反覆強調：它們都是叢集模型，**都無法以高機率預報未來
# 的大地震**。它們的強項是餘震，不是大地震。凡是宣稱能預報大地震的
# 系統，先問它有沒有通過前瞻檢驗與基準比較。

# %% [markdown]
# ## 14.12 附錄：本章推導細節
#
# ### A. 空間積分的極座標封閉式
#
# 空間核（記 $\sigma = De^{\gamma(m-m_0)}$）為
# $f(x,y;m)=\frac{q-1}{\pi\sigma}(1+r^2/\sigma)^{-q}$，$r^2=x^2+y^2$。
# 轉極座標，$\mathrm{d}x\,\mathrm{d}y=r\,\mathrm{d}r\,\mathrm{d}\varphi$，
# 令 $u=1+r^2/\sigma$（故 $r\,\mathrm{d}r=\tfrac{\sigma}{2}\mathrm{d}u$）：
#
# $$\begin{aligned}
# \int_0^{2\pi}\!\!\int_0^{R} f\,r\,\mathrm{d}r\,\mathrm{d}\varphi
#   &= \int_0^{2\pi}\frac{q-1}{\pi\sigma}\cdot\frac{\sigma}{2}
#      \int_1^{1+R^2/\sigma} u^{-q}\,\mathrm{d}u \,\mathrm{d}\varphi \\
#   &= \int_0^{2\pi}\frac{q-1}{2\pi}
#      \left[\frac{u^{1-q}}{1-q}\right]_1^{1+R^2/\sigma}\mathrm{d}\varphi \\
#   &= \frac{1}{2\pi}\int_0^{2\pi}
#      \left[1-\left(1+\frac{R^2}{\sigma}\right)^{1-q}\right]\mathrm{d}\varphi .
# \end{aligned}$$
#
# $R\to\infty$ 時（需 $q>1$）括號趨近 1，整式等於 1，正規化確認。在
# 有限區域上把 $R$ 換成隨方位角變動的邊界距離 $R(\varphi)$，就得到
# 14.1 節的 $F_i$。實作上把 $[0,2\pi)$ 切成 $n_v$ 段、每段取
# $R(\varphi)$ 為到多邊形邊界的距離，用梯形法求和——這就是
# **radial partitioning**。它比二維網格快兩三個數量級，原因是徑向
# （被積函數變化最劇烈的方向）已經被解析積掉了。
#
# ### B. DFP 更新式與標準誤
#
# 令 $\Delta\theta_k=\theta_{k+1}-\theta_k$、
# $y_k=\nabla\xi(\theta_{k+1})-\nabla\xi(\theta_k)$。DFP 的反 Hessian
# 更新為
#
# $$H_{k+1} = H_k
#   + \frac{\Delta\theta_k\,\Delta\theta_k^{\mathsf T}}
#          {\Delta\theta_k^{\mathsf T} y_k}
#   - \frac{H_k\,y_k\,y_k^{\mathsf T}\,H_k}{y_k^{\mathsf T} H_k y_k} .$$
#
# 兩個修正項各有職責：第二項注入新的曲率資訊（並保證滿足擬牛頓條件
# $H_{k+1}y_k=\Delta\theta_k$），第三項移除已經過時的舊曲率。只要
# $\Delta\theta_k^{\mathsf T}y_k>0$（線搜尋滿足 Wolfe 條件時成立），
# $H_{k+1}$ 保持正定。收斂時
# $H_\infty\approx[-\nabla^2 l_2(\hat\theta)]^{-1}$，於是
# $\mathrm{SE}(\hat\theta_j)=\sqrt{(H_\infty)_{jj}}$。
#
# **警告**：這是漸近結果，假設概似在 $\hat\theta$ 附近近似二次、
# $\hat\theta$ 在參數空間內部、樣本量夠大。14.3 節那條香蕉形山脊
# 三個條件全部踩線，所以「$\gamma$ 的標準誤 5.7553」不該讀成
# 「95% 區間是 $[-11,11]$」，而該讀成**這個方向上概似幾乎是平的**。
#
# ### C. 不可辨識性的完整條件
#
# 承 14.3 節。設 $\ell(\eta)=\sum_i \ell_i(u_i)$，
# $u_i=\eta_1+\eta_2\,\Delta m_i$（$\eta_1=\ln D,\ \eta_2=\gamma$；
# 或 $\eta_1=\ln A,\ \eta_2=\alpha$）。記
# $\iota_i=-\mathbb{E}[\partial^2\ell_i/\partial u_i^2]>0$、
# $W=\sum_i\iota_i$、$\bar{\Delta m}=\sum_i\iota_i\Delta m_i/W$。
# 資訊矩陣與其行列式為
#
# $$I = \begin{pmatrix}
#   \sum\iota_i & \sum\iota_i\Delta m_i \\
#   \sum\iota_i\Delta m_i & \sum\iota_i\Delta m_i^2
# \end{pmatrix},
# \qquad
# \det I = W^2 \operatorname{Var}_\iota(\Delta m) ,$$
#
# 其中 $\operatorname{Var}_\iota(\Delta m)
# =\sum_i\iota_i(\Delta m_i-\bar{\Delta m})^2/W$。反矩陣的對角元素給出
#
# $$\operatorname{Var}(\hat\eta_2)
#   = \frac{\sum\iota_i}{\det I}
#   = \frac{1}{W\operatorname{Var}_\iota(\Delta m)} .$$
#
# 三個推論。**其一**，所有事件規模相同時 $\operatorname{Var}_\iota=0$、
# $\det I=0$，參數**嚴格不可辨識**。**其二**，
# $\mathrm{SE}(\hat\eta_2)\propto[W\operatorname{Var}_\iota(\Delta m)]^{-1/2}$
# ——想把 $\gamma$ 或 $\alpha$ 估準，**加大規模跨度比加大事件數更
# 有效**（前者進變異數，後者只進 $W$，而變異數通常是瓶頸）。
# **其三**，兩個估計量的相關係數是
# $-\bar{\Delta m}/\sqrt{\overline{\Delta m^2}}$，**恆為負**——這解釋
# 了 14.6 節那張圖：$\hat\alpha$ 被壓低時 $\hat A$ 必然上升。
#
# ### D. R–J 的 $a'$ 與 ETAS 的 $A$：dressed 對 bare
#
# 14.7 節推出 $a'=\log_{10}[A(p-1)c^{p-1}]$。用加州通用參數反解：
#
# $$A_{\rm eff} = \frac{10^{a'}}{(p-1)c^{\,p-1}}
#   = \frac{10^{-1.67}}{0.08\times 0.05^{0.08}} \approx 0.34 .$$
#
# 這個數字偏大，原因是 R–J 是**用整段餘震序列擬合出來的**，而整段
# 序列包含餘震的餘震、餘震的餘震的餘震——全部世代。在分支過程裡，
# 一個事件的**全世代**期望後代數是
#
# $$\kappa_{\rm dressed} = \kappa\,(1 + n + n^2 + \cdots)
#   = \frac{\kappa}{1-n} ,$$
#
# 所以 $A_{\rm eff}=A_{\rm bare}/(1-n)$：R–J 的產能是 **dressed**
# （已含級聯），ETAS 的 $A$ 是 **bare**（只算直接後代）。
#
# 這可以做成自洽性檢查。$\alpha=\beta$ 時分支比要靠 $M_{\max}$ 截斷：
# $n=A_{\rm bare}\,\beta\,(M_{\max}-m_0)$。取加州的 $b=0.91$
# （$\beta=2.095$）、$m_0=3$、$M_{\max}=8$，代入
# $A_{\rm bare}=(1-n)A_{\rm eff}$：
#
# $$n = (1-n)\times 0.34 \times 2.095\times 5 = 3.56\,(1-n)
#   \quad\Longrightarrow\quad n = \frac{3.56}{4.56} \approx 0.78 .$$
#
# $n\approx0.78$、$A_{\rm bare}\approx0.075$——分支比落在文獻常見的
# 0.3–0.9 範圍內，而且相當接近臨界。**一條 1989 年的經驗公式，用
# ETAS 的語言重讀之後，內含的分支比竟然是合理的**，這是「R–J 是
# ETAS 特例」最漂亮的旁證。
#
# 但要標註限度：dressed／bare 的 $1/(1-n)$ 關係對**期望總數**是精確
# 的，對**時間形狀**不是——級聯會讓觀測到的衰減比 bare 的 $p$ 更平
# （Helmstetter & Sornette 的重整化結果），所以 R–J 的 $p=1.08$ 與
# ETAS 的 bare $p$ 也不該直接畫等號。這是量級檢查，不是參數轉換
# 公式。
#
# ---
#
# 回頭看這一章走過的路：從一條可分離的概似出發，經過邊界效應、
# 擬牛頓法、標準誤，走到「有些參數根本沒被資料約束」這個不太舒服
# 的事實；再用隨機除叢化解第 12 章的除叢難題，用隨機時間變換把
# 「模型錯了嗎」變成一張看得懂的圖；最後把 R–J、simplETAS 與台灣
# 的兩組參數擺在同一張桌子上。
#
# 壓成一句話：**ETAS 的價值不在於它估出來的八個數字，而在於它提供
# 了一個「正常」的定義。** 有了正常，異常才有意義——$\tau$–$j$ 圖上
# 的相對寧靜、殘差地圖上系統性為正的區塊、大埔序列那段「在 Omori
# 框架下不尋常的平靜」，全都是相對於這個參考模型才說得出口的話。
#
# 而這正好帶到下一個問題。ETAS 描述的是**觸發**——這個地震會引發
# 什麼。但地震序列裡還有另一類現象，時間之箭是反過來的：一群小地震
# 在數月到數年之後，被一場更大的地震跟上。那不是觸發，是**預示**。
# {doc}`第 15 章 <15_psi_phenomenon>`要處理的 Ψ 現象，就是把這條
# 反向的箭頭寫成可檢驗的迴歸關係——而它能被辨識出來的前提，正是
# 我們現在已經有了一個描述「正常叢集」的模型。
