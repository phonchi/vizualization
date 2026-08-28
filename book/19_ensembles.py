# %% [markdown]
# # 19. 模型組合：加法、乘法與權重學習
#
# {doc}`第 18 章 <18_testing_comparison>`留下一個很尷尬的處境。把
# 一群模型送進 CSEP，跑完 N／S／cL、跑完 T-test 與 Molchan，最常
# 見的結局不是「甲模型顯著勝出」，而是**一堆誤差棒都跨過零線**：
# 統計上誰也沒有顯著贏過誰，可是它們對「下一次地震落在哪裡」的
# 看法南轅北轍。決策者不能等三百年累積樣本，他明天早上就要一個
# 數字。
#
# 這一章的答案是：**也許不必選**。如果幾個模型各自捕捉到真實地震
# 活動的不同面向，把它們組合起來可能比任何單一模型都好——而且
# 19.3 節會證明，在一個乾淨的條件下，這句話是**定理**而不是希望。
#
# 但第 18 章同時也把帳單先開好了。組合模型多出來的每一個權重都是
# 一個參數，都要在 {eq}`eq:igpec` 的懲罰項裡付出代價；而且組合的
# 目標函數與評分標準是**同一條式子**——第 10 章的
# $\ln L=\sum_i\ln\lambda^*-\int\lambda^*$。這個巧合聽起來方便，
# 其實是本章最大的陷阱：**用來擬合的分數不能拿來當成績。** 這一章
# 因此分成三段：19.1–19.3 是加法家族與權重的最大概似，19.4–19.6
# 是乘法家族與它在真實前瞻測試中的翻車現場，19.7–19.8 回到加法，
# 但把問題從「用什麼運算子」換成「權重到底該學什麼」。
#
# 記號約定：組合中第 $i$ 個成分模型的權重寫 $\pi_i$（$\sum\pi_i=1$
# 時稱**凸組合**），乘法家族的基準模型寫 $\lambda_1$、保序轉換寫
# $f_i$、整體正規化參數寫 $a_0$。IGPE、IGPEc 與它們的變異數估計式
# 由第 18 章擁有，本章一律引用不重推。
#
# ## 19.1 互補性從何而來
#
# 「組合會更好」不是無條件成立的。把同一個模型複製十份再平均，
# 得到的還是它自己。組合的燃料是**互補性**，而互補性在地震預報裡
# 有兩個具體來源。
#
# **來源一：時間尺度。** 前面幾章的引擎剛好鋪滿了整條時間軸。
# ETAS 與 STEP（{doc}`第 13 章 <13_etas_structure>`、
# {doc}`第 14 章 <14_etas_estimation>`）吃的是 Omori–Utsu 衰減，
# 特徵時間是**小時到數週**；EEPAS（{doc}`第 16 章 <16_eepas_ppe>`）
# 吃的是前兆尺度增長，特徵時間是**數月到數十年**；PPE 與各種平滑
# 地震度模型吃的是長期空間分布，特徵時間是**數十年不變**。三者用
# 的是**同一份地震目錄**，抽取的卻是完全不同的訊號。
#
# **來源二：資料本身。** 目錄以外還有大地測量。GNSS 應變率、活動
# 斷層滑移率、地下水與地磁異常，這些量與地震目錄的相關性遠低於
# 「兩個都用目錄的統計模型」之間的相關性。19.10 節會說明，這才是
# 組合模型真正的成長空間。
#
# 先看第一種互補性長什麼樣子。對同一個目標地震，短期模型與中期
# 模型的機率軌跡是兩種完全不同的形狀——前者是事件之後的尖峰急衰，
# 後者是提前數年的緩坡：

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["hide-input"]
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
# 中期模型（綠）在事發前幾年就把機率抬起來，但抬得溫和；短期模型
# （橘）平時貼著背景，只在前震出現後的幾天內爆衝。凸組合（藍虛線）
# 兩邊都吃得到：平時繼承中期模型的緩坡，前震一來立刻繼承短期模型
# 的尖峰。注意縱軸是對數軸——兩條線在同一個時刻可以差好幾個數量級。
#
# ### 19.1.1 十二個數量級
#
# 「差好幾個數量級」不是誇飾。Rhoades 與 Gerstenberger（2009）把
# STEP 與 PPE 在同一套時空格上的預報率相除，比值**跨越十二個
# 數量級**，而其中約**一半落在 0.1 到 10 之間**。這兩個數字要一起
# 讀才有意義：分布的**主體**是溫和的（兩個模型大致同意），但**尾巴
# 極長**（在少數格子上一個模型認為機率是另一個的百萬倍）。下面用
# 合成資料重現這個形狀——短期模型的率是「背景 × 對數常態擾動 ×
# 偶發的叢集爆衝」，長期模型的率只有背景與一個相關性較低的擾動：

# %% tags=["hide-input"]
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
# 圖題裡的兩個數字由程式量出來，與文獻報告的「十二個數量級、一半
# 在 0.1–10」是同一個形狀。**差異大就是機會大**：如果兩個模型處處
# 一致，組合起來不會多出任何東西；正因為它們在尾巴上劇烈分歧，
# 才有「這一格該聽誰的」這個值得學習的問題。
#
# ### 19.1.2 投資組合類比，以及它的界線
#
# 這裡有一個很自然的類比：把每個模型當成一項資產，組合就是投資
# 組合。兩項資產的報酬率變異數分別是 $\sigma_1^2, \sigma_2^2$、
# 相關係數 $\rho$，權重 $\pi$ 與 $1-\pi$，組合的變異數是
#
# $$\sigma^2(\pi) = \pi^2\sigma_1^2 + (1-\pi)^2\sigma_2^2
# + 2\pi(1-\pi)\rho\,\sigma_1\sigma_2 .$$
#
# 對 $\pi$ 微分並令其為零，最小變異數的權重是
#
# $$\pi^{*} = \frac{\sigma_2^2 - \rho\,\sigma_1\sigma_2}
# {\sigma_1^2 + \sigma_2^2 - 2\rho\,\sigma_1\sigma_2} .$$
#
# 重點在 $\pi^*$ 什麼時候落在 $(0,1)$ 內。分子為正需要
# $\rho < \sigma_2/\sigma_1$；同理 $\pi^*<1$ 需要
# $\rho < \sigma_1/\sigma_2$。合起來就是
#
# $$\rho < \min\left(\frac{\sigma_1}{\sigma_2},\
# \frac{\sigma_2}{\sigma_1}\right)
# \quad\Longleftrightarrow\quad \pi^{*}\in(0,1).$$
#
# **只要兩項資產的相關性夠低，最佳配置一定是「兩個都持有」**，
# 即使其中一項單獨看又爛又波動。這正是本章要傳達的直覺。
#
# 但類比到此為止，三條界線必須講清楚。**其一，目標函數不同**。
# 投資組合最小化變異數，地震預報最大化期望對數分數；19.3 節會
# 證明後者有自己的內點最優條件，形式與上式不同。**其二，不能
# 放空**。$\pi_i<0$ 會讓組合的率在某些格子變成負數，機率式預報
# 沒有這種東西——這也是為什麼 19.7 節的權重映射一定要把負係數
# 截成零。**其三，樣本結構不同**。金融報酬率一年有兩百多個交易
# 日的觀測，地震預報一個測試期只有幾十顆目標地震，而且它們高度
# 叢集：**一個序列就可能主宰整個估計**（第 18 章的 18.8 節）。
# 分散化的邏輯成立，分散化的**估計精度**卻遠不如金融。
#
# ## 19.2 加法家族：三種混合與守恆定理
#
# 最早的實測從最簡單的形式開始。Rhoades 與 Gerstenberger（2009）
# 用加州 ANSS 目錄 1984–2004 的 152 個 $M \ge 5.0$ 事件，測試三種
# 加法混合。要先注意的是：**兩個母模型本身就已經是混合結構**。
#
# $$\lambda_{\rm STEP} = \max\bigl(\lambda_{\rm CLUST},\
# \lambda_{\rm STAT}\bigr), \qquad
# \lambda_{\rm EEPAS} = \lambda_{\Psi} + \mu_E\,\lambda_{\rm PPE}$$
#
# $\lambda_{\rm CLUST}$ 是 STEP 的時變群聚項、$\lambda_{\rm STAT}$
# 是它的靜態平滑地震度背景；$\lambda_{\Psi}$ 是 EEPAS 的前兆尺度
# 增長項、$\lambda_{\rm PPE}$ 是準靜態背景，混合權重
# $\mu_E \in [0,1]$ 代表「沒有前兆尺度增長的大地震比例」（記號
# 依本部規範由文獻的 $\mu$ 改寫，避免與 ETAS 背景率撞名）。
# **EEPAS 本身就是「時變項 + 背景項」的加法混合**——組合這件事
# 在模型內部早就發生了。
#
# 三個新混合模型是：
#
# $$\begin{aligned}
# \lambda_{\rm SE1} &= \lambda_{\rm CLUST} + q\,\lambda_{\rm EEPAS},
# & 0 &\le q \le 1, \\
# \lambda_{\rm SE2} &= (1-r)\,\lambda_{\rm STEP} + r\,\lambda_{\rm EEPAS},
# & 0 &\le r < 1, \\
# \lambda_{\rm SE3} &= \lambda_{\rm CLUST} + P(m)\,\lambda_{\rm EEPAS}. &&
# \end{aligned}$$
#
# 三者的差別值得一個一個看。**SE1 把 STEP 的靜態背景整個換掉**：
# 保留時變的 $\lambda_{\rm CLUST}$，用 EEPAS（自帶 PPE 背景）
# 取代 $\lambda_{\rm STAT}$。它是加法，但**不是凸組合**——權重
# $1$ 與 $q$ 加起來不等於一。**SE2 才是凸組合**，STEP 整個保留，
# 與 EEPAS 按 $(1-r):r$ 分配。**SE3 讓權重隨規模變**，下面單獨談。
#
# ### 19.2.1 凸組合為什麼不必再正規化
#
# 凸組合有一個常被順口帶過、其實需要一行證明的性質。設各成分
# 模型在測試區與測試期上的期望總數是
# $\Lambda_i = \int\lambda_i$，組合是
# $\lambda_\pi = \sum_i \pi_i\lambda_i$，$\pi_i\ge0$、
# $\sum_i\pi_i=1$。積分是線性算子，所以
#
# $$\begin{aligned}
# \Lambda_\pi &= \int \sum_i \pi_i \lambda_i
#  = \sum_i \pi_i \int \lambda_i = \sum_i \pi_i \Lambda_i .
# \end{aligned}$$
#
# 結論分兩層，不要混為一談。**一般情形**：組合的期望總數是各成分
# 期望總數的**加權平均**，因此必然落在
# $\bigl[\min_i\Lambda_i,\ \max_i\Lambda_i\bigr]$ 之內——**凸組合
# 在總數上不可能比最差的成分更離譜**。**特例**：若各成分模型都
# 已經校準到同一個區域平均地震率，$\Lambda_i \equiv \Lambda_0$，
# 則 $\Lambda_\pi = \Lambda_0\sum_i\pi_i = \Lambda_0$，**與權重
# 完全無關**。這時掃描權重不會動到總數，N-test 層次自動通過，
# 也不需要任何額外的正規化參數。
#
# 對照組是 SE1：$\Lambda_{\rm SE1} = \Lambda_{\rm CLUST}
# + q\,\Lambda_{\rm EEPAS}$ 隨 $q$ 單調遞增。它是加法，卻沒有
# 守恆保證——<strong>「加法」與「凸」是兩件事</strong>，而 19.4 節的乘法家族
# 之所以需要一個整體參數 $a_0$，正是因為它連加權平均這層保護
# 都沒有。
#
# ### 19.2.2 讓權重隨規模變：SE3 與奧坎剃刀
#
# SE3 的想法很漂亮。EEPAS 預報的是**主震**（獨立事件），STEP 的
# 群聚項預報的是**餘震**；規模愈大的事件愈可能是主震，所以
# EEPAS 的權重應該隨規模上升。用 Reasenberg 除叢給出的獨立性
# 機率做 logistic 迴歸：
#
# $$\ln\frac{P(m)}{1-P(m)} = a_l + b_l\,m,
# \qquad \hat a_l = -0.51,\quad \hat b_l = 0.19$$
#
# 代進去算一下就有感覺（$P$ 依本部規範一律寫成 $P(\cdot)$，
# 小寫 $p$ 保留給 Omori 指數）：
#
# | 規模 $m$ | $a_l+b_l m$ | $P(m)$ | 讀法 |
# |---|---|---|---|
# | 5.0 | $+0.44$ | 約 0.61 | 六成是主震 |
# | 6.0 | $+0.63$ | 約 0.65 | |
# | 7.0 | $+0.82$ | 約 0.69 | 七成是主震 |
#
# 結果呢？SE3 相對 SE1 的對數概似只增加 $0.1$ 到 $0.2$，卻多花了
# 一個參數。用第 18 章的懲罰算術：AIC 對每個參數收 $2$，
# $\Delta\mathrm{AIC} = -2\times0.15 + 2 = +1.7 > 0$——**扣完參數
# 之後 SE3 比 SE1 差**。作者最後選了較簡單的 SE1／SE2。這是一個
# 很乾淨的奧坎剃刀教材：**物理直覺對，不代表資料付得起這個直覺
# 的價錢**。31 個或 152 個目標地震，養不起太多自由度。
#
# ## 19.3 權重的最大概似：內點最優是定理
#
# 權重從哪裡來？跟這本書從第 10 章以來的一切一樣：**最大概似**。
# 目標函數就是 {doc}`第 10 章 <10_point_process>`的
# $\ln L = \sum_{i}\ln\lambda^*(t_i,x_i,y_i,m_i) - \int\lambda^*$，
# 把 $\lambda^*$ 換成 $\lambda_\pi$、對 $\pi$ 最佳化即可。網格版
# 用第 17 章的聯合 POLL，實務上以 Nelder–Mead 單體法求解。
#
# 這裡要先把一件事講死：**這個 $\ln L$ 同時是擬合的目標函數與
# 評分的標準**。所以擬合期的 $\ln L$ 上升**在數學上是必然的**，
# 它度量的是參數個數而不是技巧——回溯期一律要用第 18 章的
# {eq}`eq:igpec` 扣參數，理由 19.4.4 節會再具體化一次。
#
# ### 19.3.1 為什麼最佳點總是在內部
#
# 先看兩成分的凸組合 $\lambda_r = (1-r)\lambda_1 + r\lambda_2$。
# 對數概似對 $r$ 的形狀有一個結構性的答案：
#
# $$\ln L(r) = \sum_{n=1}^{N}\ln\bigl[(1-r)\lambda_1(n)
# + r\lambda_2(n)\bigr] - \bigl[(1-r)\Lambda_1 + r\Lambda_2\bigr]$$
#
# 第二項對 $r$ 是**線性**的；第一項是「對數 ∘ 仿射」的和，**對
# $r$ 嚴格凹**。凹函數加線性函數還是凹函數，所以 $\ln L(r)$ 在
# $[0,1]$ 上是凹的，最大值唯一，而且**沒有局部陷阱**——一維掃描
# 就能找到全域解，這是加法家族相對乘法家族的一個實務優勢。
#
# 接著問最大值落在哪裡。假設兩個模型都已校準到同一個期望總數
# $\Lambda_1=\Lambda_2$（19.2.1 節的特例），令第 $n$ 顆目標地震
# 所在時空格上的**率比**為
# $\varrho_n = \lambda_2(n)/\lambda_1(n)$。兩個端點的導數是
#
# $$\begin{aligned}
# \left.\frac{\mathrm{d}\ln L}{\mathrm{d}r}\right|_{r=0}
# &= \sum_{n=1}^{N}\bigl(\varrho_n - 1\bigr)
#  = N\left[\overline{\varrho} - 1\right], \\
# \left.\frac{\mathrm{d}\ln L}{\mathrm{d}r}\right|_{r=1}
# &= \sum_{n=1}^{N}\left(1 - \frac{1}{\varrho_n}\right)
#  = N\left[1 - \overline{\varrho^{-1}}\right],
# \end{aligned}$$
#
# 其中上劃線代表對 $N$ 顆目標地震取算術平均。所以**最佳權重嚴格
# 落在內部**的充要條件是
#
# $$\overline{\varrho} > 1 \quad\text{and}\quad
# \overline{\varrho^{-1}} > 1 .$$
#
# 這兩個條件看起來互相矛盾——一個說模型 2 的率平均比較高，另一個
# 說模型 1 的率平均比較高——但它們可以同時成立，因為算術平均對
# 大值敏感。事實上算術–幾何平均不等式直接給出：
#
# $$\overline{\varrho} \ \ge\ \mathrm{GM}(\varrho)
# = e^{\,\delta}, \qquad
# \overline{\varrho^{-1}} \ \ge\ \frac{1}{\mathrm{GM}(\varrho)}
# = e^{-\delta},$$
#
# 其中 $\delta = \frac{1}{N}\sum_n\ln\varrho_n$ **恰好就是模型 2
# 相對模型 1 的每地震資訊增益 IGPE**（總數已對齊，
# {eq}`eq:igpe` 的修正項為零）。於是有一個乾淨的結論：
#
# > **若兩個模型的期望總數相同、每地震資訊增益恰為零（$\delta=0$），
# > 而且它們的率比不是處處相等，則兩個不等式都嚴格成立，$\ln L(r)$
# > 的最大值必定落在 $(0,1)$ 的內部。**
#
# 換句話說：**當第 18 章的 T-test 判定兩個模型「統計上無法區分」
# 時，凸組合嚴格優於它們兩個。** 平手不是無話可說，平手是組合的
# 最佳時機。$\delta \ne 0$ 時結論不再自動成立，但由不等式可以讀出
# 條件——弱模型要進得了組合，它必須在**某些**地震上贏得夠多，
# 大到讓 $\overline{\varrho^{-1}}$ 突破 $1$。這就是 19.1 節那條
# 長尾巴的價值：**互補性的形式化定義，是率比分布的離散度**。
#
# 下面用合成資料把整條曲線跑出來。真實率是兩個模型的七三混合，
# 兩個模型都正規化到同樣的期望總數 $\Lambda_0=150$：

# %% tags=["hide-input"]
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
# 三件事值得停下來看。**第一，曲線是凹的、峰在內部**，兩端都明顯
# 低於最佳點，程式裡的兩個 `assert` 就是在驗算 19.3.1 的端點條件。
# **第二，模型 B 單獨評分是輸的**（圖題裡的 IGPE 為負），它在整條
# 時間軸上大部分時候都貼著近乎零的率，只有第 60 年附近爆一個尖峰
# ——但**它仍然拿到將近三成的權重**。原因就是那個尖峰：$\varrho_n$
# 在少數幾顆地震上大到把 $\overline{\varrho}$ 推過 $1$。
# **第三，因為兩個模型的期望總數都是 150，整條掃描過程中組合的
# 期望總數恆等於 150**——19.2.1 節的守恆定理在這裡是可以直接
# 驗算的。
#
# ### 19.3.2 2009 年加州實驗：一個參數換到兩倍
#
# 真實資料的結果比玩具還漂亮。加州 20 年、152 個 $M\ge5$ 事件，
# 相對於<strong>平穩均勻 Poisson 模型（SUP）</strong>的資訊率與平均機率增益
# $G=\exp(\text{資訊率})$：
#
# | 模型 | 資訊率 | 機率增益 $G$ |
# |---|---|---|
# | PPE | 1.67 | 5.3 |
# | EEPAS | 2.17 | 8.8 |
# | STEP | 約 2.1–2.4 | 8–12 |
# | SE1（CLUST + EEPAS） | 約 2.9–3.2 | 18–25 |
# | **SE2（STEP ⊕ EEPAS）** | **約 3.2–3.3** | **25–28** |
#
# 最佳解是 $r = 0.42$，也就是 $0.58\times$ STEP $+\ 0.42\times$
# EEPAS。相對 STEP 的平均機率增益是 **2.72 倍**，相對 EEPAS 也
# **超過兩倍**——而這只多擬合了**一個**參數。以第 18 章的懲罰
# 算術檢查：$N=152$、$n_p=1$ 時懲罰只有
# $\frac{1}{2N}\bigl(2+\frac{2}{150}\bigr)\approx0.0066$，相對
# 增益完全可以忽略。**這是本章唯一一個「便宜到不必猶豫」的組合。**
#
# 表裡還藏著本節第二個教訓。SE2（保留整個 STEP）勝過 SE1（丟掉
# STEP 的靜態背景 $\lambda_{\rm STAT}$）。而 $\lambda_{\rm STAT}$
# 單獨評分是全場最弱的成分——它就是一個比 PPE 還粗糙的平滑地震度
# 模型。**把最弱的成分拿掉，組合反而變差**：它在少數格子上帶著
# 別人沒有的特徵，那些格子剛好有目標地震。
#
# **不要只看總體排名就把模型踢出候選池。** 有效共軛的判準不是
# 「它自己很強」，而是「它與地震發生的相關方式，是其他成分沒有
# 捕捉到的」——這句話在 19.4 節還會以另一種形式回來。
#
# ## 19.4 乘法家族：把其他模型當修正因子
#
# 加法有個代價：**增益會被稀釋**。EEPAS 相對 PPE 有 0.5 的資訊
# 優勢，可是把 EEPAS 換進混合之後，SE1 相對對應的加法組合只剩
# 0.13。Aki（1981）、Utsu（1983）與 Imoto（2007）的**條件獨立
# 前兆理論**暗示了另一條路：如果前兆與地震的關係是**乘法**的，
# 那就該用乘法組合，才能保留完整的機率增益。
#
# Rhoades 等人（2014）在 RELM 五年實驗（加州，$M \ge 4.95$，
# 31 個目標地震）上實作了這個想法。RELM 的一階結論是
# Helmstetter 等人（2007）的平滑地震度模型（HKJ）最好，於是問題
# 變成：**能不能以 HKJ 為基準，把其他模型當修正因子乘上去？**
#
# $$
# \lambda_H(j,k) = \lambda_1(j,k)\,
# \exp\Bigl[a_0 + \sum_{i=2}^{n_i} f_i\bigl(\lambda_i(j,\cdot)\bigr)\Bigr],
# \qquad f_i(\lambda) = u_i\bigl[\log(1+\lambda)\bigr]^{v_i}
# $$ (eq:mulhyb)
#
# $j$ 是空間格、$k$ 是規模箱；$\lambda_1$ 是基準模型的格點期望數；
# $\lambda_i(j,\cdot)=\sum_k\lambda_i(j,k)$ 是第 $i$ 個**共軛模型**
# 對規模求和後的空間期望數；$a_0$ 是整體正規化參數，$u_i \ge 0$
# 與 $v_i > 0$ 是形狀參數（記號依本部規範由文獻的 $a, b_i, c_i$
# 改寫）。對照組的加法 hybrid 是
# $\lambda_H(j,k)=\sum_i a_i\lambda_i(j,k)$，$a_i\ge0$。
#
# 這個式子有四個設計決定，每一個都值得單獨拆開。
#
# ### 19.4.1 為什麼共軛模型要先對規模求和
#
# 隱含假設是：**共軛模型除了基準模型已有的規模資訊之外，不再提供
# 額外的規模資訊**。這不只是方便，在一個條件下它是**恆等式**。
# 假設所有模型都遵守空間不變、$b$ 值相同的 Gutenberg–Richter 律
# （第 11 章），則每個模型的格點期望數都可以分離成
#
# $$\lambda_i(j,k) = \lambda_i(j,\cdot)\, s_k,
# \qquad s_k = \int_{\text{bin }k} s(m)\,\mathrm{d}m,$$
#
# 規模密度 $s(m)=\beta e^{-\beta(m-m_0)}$ 與空間格無關。於是
#
# $$\frac{\lambda_i(j,k)}{\lambda_1(j,k)}
# = \frac{\lambda_i(j,\cdot)\,s_k}{\lambda_1(j,\cdot)\,s_k}
# = \frac{\lambda_i(j,\cdot)}{\lambda_1(j,\cdot)}$$
#
# **與 $k$ 無關**——乘數就算想依賴規模也沒有材料可依賴。所以
# {eq}`eq:mulhyb` 的乘數只寫成 $j$ 的函數，不是簡化，是這個假設
# 下的必然。實務意義：**規模分布由基準模型全權負責，共軛只負責
# 調整空間分布**。反過來說，如果某個共軛模型真的帶有規模資訊
# （例如它預測某些區域的 $b$ 值特別低），這個框架**接不住**。
#
# ### 19.4.2 為什麼只用排序、不用數值
#
# $f_i(\lambda)=u_i[\log(1+\lambda)]^{v_i}$ 在 $u_i\ge0$、
# $v_i>0$ 時對 $\lambda$ 單調非遞減，是**保序（order-preserving）
# 轉換**。求和後的空間期望數因此被當成 Zechar 與 Jordan（2008）
# 意義下的 **alarm function**（第 18 章的 Molchan 圖用的就是同一個
# 東西）：重要的是空間格之間的**排序**，不是絕對值。
#
# 這一步的報酬很大。共軛模型可以被任意單調重標度——換單位、換
# 校正基準、乘一個未知常數——排序不變，只要重新擬合 $u_i, v_i$
# 就能吸收掉尺度的差異。於是**共軛模型不必是一個完整的預報模型**：
# GNSS 應變率地圖、活動斷層滑移率、地震圖樣指標（PI）、甚至
# 一個「這格有沒有前兆異常」的二元變數，都能塞進 $f_i$。
# 這是本框架對「多源觀測資料整合」最直接的貢獻。
#
# ### 19.4.3 為什麼每個共軛只給兩個參數
#
# 作者自己下的比喻是**多元迴歸**：每加一個解釋變數只多付兩個
# 參數，就能解釋更多變異。兩個參數已經夠彈性——在
# $\lambda \ll 1$ 的範圍內 $\log(1+\lambda)\approx\lambda$，於是
# $f_i \approx u_i\lambda^{v_i}$，$v_i=1$ 近似線性、$v_i>1$ 上凸、
# $v_i<1$ 下凹、$u_i=0$ 退化成常數。
#
# 而節制的理由是算術。目標地震只有 **31 個**。第 18 章的懲罰項
# 在 $n_i=2$（兩個共軛，$n_p = 2\times2+1 = 5$ 個參數）時是
#
# $$\frac{1}{2N}\Bigl(2n_p + \frac{n_p+1}{N-n_p-1}\Bigr)
# = \frac{1}{62}\Bigl(10 + \frac{6}{25}\Bigr) \approx 0.165 .$$
#
# 而文獻報告的回溯 IGPEc 也才 0.25–0.35。**懲罰吃掉了一半以上的
# 名目增益**——這不是保守，這是這個樣本數下的現實。與多元迴歸
# 同理：**觀測要多、解釋變數要少**。
#
# ### 19.4.4 退化保證：為什麼一定要用 IGPEc
#
# 這一點最關鍵。在 {eq}`eq:mulhyb` 裡取 $a_0=0$ 且所有 $u_i=0$，
# 乘數恆等於
# $\exp[0+\sum_i 0]=1$，於是 $\lambda_H \equiv \lambda_1$。
# **參數空間包含了「hybrid 就是基準模型」這個點**。最大概似是在
# 整個參數空間上取上界，所以
#
# $$\ln L_H^{\max} \ \ge\ \ln L\bigl(a_0=0, u_i=0\bigr)
# = \ln L_1 .$$
#
# **不等式沒有例外**：擬合期的 hybrid 對數概似**在數學上不可能
# 低於基準**，就算共軛模型是純亂數也一樣。所以回溯期的
# $\Delta\ln L>0$ **不含任何資訊**——它度量的是參數個數。第 18 章
# 的 {eq}`eq:igpec` 存在的唯一理由就是這件事：回溯評估必須用
# IGPEc（扣 AICc 懲罰），前瞻評估才可以用未懲罰的 IGPE。
#
# ### 19.4.5 回溯成績單
#
# 扣完懲罰之後的回溯結果仍然漂亮：
#
# | 組合 | 區域 | $\Delta\ln L$ 或 IGPEc |
# |---|---|---|
# | HKJ ⊗ Neokinema（大地測量） | 全加州 | $\Delta\ln L=11.4$，IGPEc 約 0.25 |
# | HKJ ⊗ PI（圖樣指標） | 全加州 | $\Delta\ln L=11.2$，IGPEc 約 0.25 |
# | HKJ ⊗ Shen（大地測量） | 南加州 | IGPEc 超過 0.5 |
# | HKJ ⊗ Neokinema ⊗ PI | 全加州 | IGPEc $=0.35\pm0.17$ |
# | HKJ ⊗ Shen ⊗ PI | 南加州 | IGPEc $=0.79\pm0.27$ |
#
# 三個重點。**其一，乘法大勝加法**：同一批模型的最佳加法 hybrid
# 只有 $\Delta\ln L=4.8$，乘法是 11.4。**其二，邊際報酬遞減**：
# 第三個模型帶來的額外增益（$\Delta$IGPEc $=0.09$ 與 $0.22$）
# 都不顯著，而全加州三模型的 $0.35\pm0.17$ 也只是勉強達到 95%。
# **其三，異質性帶來增益**：增益最大的組合是**資料來源差異最大**
# 的組合。用大地測量當共軛特別有效，因為基準 HKJ 只用地震目錄；
# 反過來 ALM 與 Ebel 這類「粗糙版的平滑地震度」幾乎沒有增益——
# **它們和基準講的是同一件事**。19.3.2 節那句話在這裡以另一個
# 姿態出現：有效的共軛不是「單獨很強」的模型，是「講了基準沒講
# 的話」的模型。
#
# ## 19.5 加法只能內插，乘法可以外推
#
# 加法與乘法的本質差異可以用一行值域論證講完，而這一行論證同時
# 說明了它們各自的失敗模式。
#
# **加法：逐格內插。** 設 $\pi_i\ge0$、$\sum_i\pi_i=1$。對任意
# 時空格 $j$，
#
# $$\begin{aligned}
# \lambda_\pi(j) &= \sum_i \pi_i\lambda_i(j)
#  \ \le\ \sum_i \pi_i \max_{i'}\lambda_{i'}(j)
#  = \max_{i'}\lambda_{i'}(j), \\
# \lambda_\pi(j) &\ \ge\ \sum_i \pi_i \min_{i'}\lambda_{i'}(j)
#  = \min_{i'}\lambda_{i'}(j) .
# \end{aligned}$$
#
# **凸組合的率永遠夾在各成分之間**，等號只在所有成分同值或權重
# 完全集中時成立。
#
# **乘法：可以外推。** {eq}`eq:mulhyb` 的乘數
# $m(j)=\exp[a_0+\sum_i f_i]$ 的值域是整個 $(0,\infty)$。只要
# $m(j) > \max_i \lambda_i(j)/\lambda_1(j)$，hybrid 在該格就高於
# **所有**成分；只要 $m(j) < \min_i \lambda_i(j)/\lambda_1(j)$，
# 它就低於**所有**成分。前提是「資料整體支持這個乘數」，但
# **結構上沒有任何東西攔著它**。
#
# ### 19.5.1 一個一定要澄清的誤讀
#
# 「加法只能內插」是**逐格**的敘述，不是**分數**的敘述。很多人
# 會由此推論「所以凸組合的資訊增益不可能超過最好的成分」——
# **這是錯的**，19.3.1 節的內點最優定理就是反例：兩個平手的模型
# 組合起來嚴格更好。
#
# 為什麼逐格內插不會傳遞到分數？因為分數是**對數的和**，而對數
# 是凹函數。由 Jensen 不等式，對每一顆目標地震
#
# $$\ln\lambda_\pi(n) = \ln\Bigl(\sum_i\pi_i\lambda_i(n)\Bigr)
# \ \ge\ \sum_i \pi_i\ln\lambda_i(n),$$
#
# 兩邊對 $N$ 顆地震平均，並假設各成分期望總數相同，就得到
#
# $$\mathrm{IGPE}(\lambda_\pi)\ \ge\ \sum_i \pi_i\,\mathrm{IGPE}(\lambda_i).$$
#
# **凸組合的資訊增益不低於各成分資訊增益的加權平均。** 這是
# 19.10 節那句「加權平均不保證大贏，但從來不落後」的嚴格版本
# ——它是 Jensen 不等式，不是經驗觀察。而**上界不存在**：組合
# 可以嚴格超過每一個成分，因為不等式左邊的 $\max$ 是逐格取的，
# 不同的地震可以由不同的成分來救。
#
# 乘法家族沒有這條 Jensen 保護。下面把兩者畫在同一條空間剖線上：

# %% tags=["hide-input"]
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
# 灰色帶是「兩個成分之間」的區域。**藍虛線（加法）永遠在帶子裡**，
# 程式裡的 `assert` 就是在驗算這件事；**紅線（乘法）大方地走出去**
# ——在共軛模型認為安靜的區段被壓到兩者之下，在兩者都認為活躍的
# 區段又被推到兩者之上。
#
# **這正是乘法的全部價值，也是它的全部風險。** 走出去的方向對了
# 就是外推的紅利，方向錯了就是災難：紅色叉號標的那一格，乘法給
# 的率比兩個成分都低了好幾倍，**目標地震一旦落在那裡，hybrid 錯得
# 比任何一個成分都離譜**。加法輸的是天花板（最差不過被稀釋成
# 平庸），**乘法輸的是地板**。
#
# ## 19.6 開獎：十年前瞻測試
#
# 故事講到這裡都很美好——回溯測試裡的美好。2014 年那篇論文的作者
# 群把 16 個乘法 hybrid 全部送進 CSEP，加上 6 個原始 RELM 模型，
# 接受 2011-01-01 至 2020-12-31 整整十年的**前瞻**測試（Bayona
# 等人，2022）。目標事件 40 個 $M \ge 4.95$，含 2016 年 Hawthorne
# 群震與 2019 年 Ridgecrest 序列。開獎結果：

# %% tags=["hide-input"]
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
# 回溯時 $+0.25$ 到 $+0.5$ 的增益——**而且是已經扣過 AICc 懲罰的
# 數字**——前瞻時全數翻負，**沒有任何一個 hybrid 顯著贏過基準
# 模型**。順帶一提，這十年的一致性檢驗也不好看：所有模型都**高估**
# 了地震數而未通過 Poisson N-test（作者判讀為 2011–2020 是加州
# 相對平靜的十年，而非模型全錯），空間上只有 KAGAN 通過 Poisson
# S-test；改用第 17 章的二元概似之後才有幾個模型過關。以 HKJ 為例，
# Poisson 下的 7652 個空格加上少數幾個多震格，那些多震格貢獻了
# **39%** 的懲罰，改成二元概似後降到 **18%**。
#
# 這是統計地震學少見的、完整走完「提出 → 送測 → 十年後開獎」流程
# 的案例。作者給出四個原因，每一個都是一堂課。
#
# **原因一：權重是用 31 個事件擬合的。** 有效樣本數是**地震顆數**，
# 不是格子數——格子有數萬個。小樣本擬出的參數與信賴區間會隨時間
# 大幅漂移。
#
# **原因二：測試期沒有大型斷層上的地震。** 這對以斷層與大地測量
# 資料為共軛的模型格外不利——它們的空間排序押在斷層帶上，而這
# 十年的地震沒有按照這個劇本走。
#
# **原因三：共軛模型自身退化。** BIRD、PI、EBEL-C 相對 HKJ 的
# IGPE 從 $-0.70／-0.31／-1.64$ 掉到 $-1.62／-3.38／-2.95$。
# 共軛變差，hybrid 當然跟著變差——**組合模型不會比它的候選池
# 更耐用**。
#
# **原因四：乘法結構會放大錯誤。** 這一條需要推導，而且推完之後
# 它比前三條深得多。
#
# ### 19.6.1 乘法放大機制的推導
#
# 所有 hybrid 共用 HKJ 為基準、規模分布相同、總量已正規化，
# **唯一的差別就是空間分布**。把 {eq}`eq:mulhyb` 寫成
# $\lambda_H(j)=\lambda_1(j)\,m(j)$，代進第 18 章的
# {eq}`eq:igpe`：
#
# $$\mathrm{IGPE} = \frac{\hat N_1 - \hat N_H}{N}
# + \frac{1}{N}\sum_{n=1}^{N}\ln m(j_n) .$$
#
# 現在用上「總量已正規化」這個條件 $\hat N_H = \hat N_1$。定義
# 基準模型自己的**空間機率分布** $w_j = \lambda_1(j)/\hat N_1$
# （$\sum_j w_j = 1$），則正規化條件就是
#
# $$\hat N_H = \sum_j \lambda_1(j)\,m(j) = \hat N_1
# \quad\Longleftrightarrow\quad \sum_j w_j\,m(j) = 1 .$$
#
# **乘數在基準測度下的加權平均恰為 1。** 對凹函數 $\ln$ 用
# Jensen 不等式：
#
# $$\sum_j w_j \ln m(j) \ \le\ \ln\Bigl(\sum_j w_j\,m(j)\Bigr)
# = \ln 1 = 0 ,$$
#
# 等號只在 $m\equiv1$（hybrid 退化成基準）時成立。這一行的意思很
# 重：**如果目標地震是按照基準模型的空間分布散落的，hybrid 的
# 期望資訊增益必為負。** 乘法 hybrid 一出生就欠著一筆赤字，它
# 必須靠「地震確實偏好高乘數的格子」把赤字賺回來。
#
# 赤字有多大？令 $q_j = w_j m(j)$ 為 hybrid 正規化後的空間分布
# （正規化條件保證 $\sum_j q_j=1$），則
# $\ln m(j)=\ln(q_j/w_j)$，於是
#
# $$\sum_j w_j\ln m(j) = -\,D_{\mathrm{KL}}(w\,\|\,q) \ \le\ 0,$$
#
# **赤字精確等於「基準分布到 hybrid 分布」的 Kullback–Leibler
# 散度**。再把目標地震真正的空間分布記作 $p$，則 IGPE 的期望值是
#
# $$\begin{aligned}
# \mathbb{E}\,[\mathrm{IGPE}]
# &= \sum_j p_j \ln\frac{q_j}{w_j} \\
# &= \sum_j p_j\ln\frac{p_j}{w_j} - \sum_j p_j\ln\frac{p_j}{q_j}
#  = D_{\mathrm{KL}}(p\,\|\,w) - D_{\mathrm{KL}}(p\,\|\,q).
# \end{aligned}$$
#
# 這是全節的鑰匙：**hybrid 贏，若且唯若它在 KL 意義下比基準更
# 接近真實的空間分布。** 而放大機制就在右邊第二項——
# $D_{\mathrm{KL}}(p\,\|\,q)$ 對「$q_j \to 0$ 而 $p_j>0$」的格子
# **沒有上界**。乘法把某格的率壓到接近零時，只要真實分布在那格
# 不是零，這一項就會被推向無窮大；而左邊的 $D_{\mathrm{KL}}(p\|w)$
# 是基準模型的固有品質，是個固定的有限數。
#
# **加法沒有這個問題**：凸組合的 $q_j \ge \pi_1 w_j$，被基準本身
# 從下方托住，$D_{\mathrm{KL}}(p\|q) \le D_{\mathrm{KL}}(p\|w)
# + \ln(1/\pi_1)$ 有界。這就是「加法輸天花板、乘法輸地板」的
# 代數版本。
#
# 最後把兩段歷史對起來。2006–2010 的回溯期，幾乎所有模型都輕鬆
# 通過 S-test，$q$ 相當接近 $p$，乘法佔盡便宜；2011–2020 的前瞻
# 期，多數模型（含 HKJ、BIRD、PI）都明確**未通過** S-test，
# $D_{\mathrm{KL}}(p\|q)$ 暴增，乘法就反過來吃虧。**同一個結構
# 在資料合作時放大優點，在資料不合作時放大缺點。**
#
# ## 19.7 權重要學什麼
#
# 前瞻失敗並沒有終結組合模型——它把問題問得更精確了。運算子
# （加法或乘法）其實不是重點，**權重怎麼學**才是。Herrmann 與
# Marzocchi（2023）在義大利作業型預報系統 OEF-Italy 上點破了一個
# 微妙但根本的錯誤。
#
# 一般化的記法是 $f_E = g[\{f_i,\pi_i\}]$，$g[\cdot]$ 是組合
# 運算子。最常見的形式是加權平均
# $\bar f(x)=\sum_{i=1}^{m}f_i(x)\pi_i$；當 $\pi_i$ 被詮釋成
# 「$f_i$ 是真模型的機率」時，這就是**貝氏模型平均**。
#
# 過去 OEF-Italy 用的是 **SMA（Score Model Averaging）**：按
# **各模型自己的歷史表現**分配權重。作者主張這不對。權重應該
# **直接最大化組合後整體的表現**——這兩件事不是同一件事。差別
# 在哪裡？SMA 只問「這個模型有多好」，不問「它好在哪裡、跟別人
# 好的地方重不重疊」。一個表現中等但**與眾不同**的模型，SMA 給
# 它中等的權重；而以組合表現為目標的擬合會給它高權重，因為它
# 補上了別人的空缺。反過來，一個表現不錯但與最強模型高度重複的
# 模型，SMA 照樣給它份額，白白稀釋。**19.3.2 節那個「最弱的成分
# 不該踢掉」的故事，在 SMA 底下根本不可能發生。**
#
# ### 19.7.1 logistic 迴歸為什麼是乘法組合
#
# 實作的工具是**多元 logistic 迴歸**。把觀測二元化：時空格 $j$
# 內有沒有至少一個目標事件，$Y_j \in \{0,1\}$，
#
# $$P(Y_j = 1 \mid \mathbf{u}_j) = \frac{1}{1+e^{-g(\mathbf{u}_j)}},
# \qquad g(\mathbf{u}_j) = \beta_0 + \sum_{i=1}^{m}\beta_i u_{i,j}$$
#
# **關鍵細節在自變數**：$u_{i,j} = \ln\phi_{i,j}$，也就是各模型
# 預報率**取對數**。把 logit 的定義展開就看得出這是什麼結構：
#
# $$\begin{aligned}
# \ln\frac{P_j}{1-P_j} &= \beta_0 + \sum_i \beta_i \ln\phi_{i,j}, \\
# \frac{P_j}{1-P_j} &= e^{\beta_0}\prod_{i=1}^{m}\phi_{i,j}^{\,\beta_i} .
# \end{aligned}$$
#
# **在 log-odds 上做線性組合，等價於在率的尺度上做乘冪組合。**
# 目標事件極稀時 $P_j \ll 1$、$P_j/(1-P_j)\approx P_j$，於是
# logistic 迴歸骨子裡就是一個**幾何加權的乘法組合**——它跟
# 19.4 節的乘法 hybrid 是表兄弟（差別在連結函數：hybrid 用
# $\exp[u_i\lambda^{v_i}]$，logistic 用 $\prod\phi_i^{\beta_i}$）。
# 這件事本身就值得警覺：19.6.1 節的放大機制對它同樣適用。
#
# ### 19.7.2 為什麼一定要丟掉截距
#
# 這是本節最精采的轉折：**該用的是係數，不是這個模型本身**。
# 直接把擬好的 logistic 模型當組合預報，實測結果是**負的**
# （相對 SMA 的 CumIGPE 為 $-0.06$ 到 $-0.10$，比 SMA 還差）。
# 病灶在 $\beta_0$，三個理由：
#
# **理由一：自我強化偏誤。** $\beta_0$ 是 baseline log-odds，它
# 吸收的是各模型在目標格上**絕對**率水準的偏移。若過去這些模型
# 在真的發生地震的格子上給的率偏高（也就是表現好），擬合會給出
# $\beta_0>0$。拿這個 $\beta_0$ 回頭去組合**同一批**模型，等於
# 「因為你過去表現好，所以我再幫你把機率調高一次」——同一份
# 證據被用了兩次。
#
# **理由二：估計必然延遲。** $\beta$ 只能用「預報視窗已經結束、
# 觀測已經到齊」的過去資料估計。短期模型在序列爆發期間預報變化
# 極快，一個延遲數週的 $\beta_0$ 會嚴重失準——**而且錯的方向
# 剛好是最糟的**：序列剛開始時最需要調高，$\beta_0$ 卻還停留在
# 平靜期的水準。
#
# **理由三：稀有事件偏誤。** 本例的目標事件比例約 **0.006%**。
# 稀有事件的 logistic 迴歸本來就會讓截距有系統性偏誤（King 與
# Zeng，2001）；而斜率係數對這個偏誤穩健得多。
#
# 解法：**只取相對資訊，丟掉截距**。把係數映射成權重：
#
# $$w_i = \begin{cases}
# e^{\beta_i} - e^{\tau}, & \beta_i > \tau \\
# 0, & \text{otherwise}
# \end{cases},
# \qquad \pi_i = \frac{w_i}{\sum_{i'} w_{i'}} .$$
#
# 為什麼是指數映射？因為 $e^{\beta_i}$ 恰好是該模型的 **odds
# ratio**：$u_i$ 每增加一單位（即 $\phi_i$ 乘上 $e$），勝算就乘
# 上 $e^{\beta_i}$。所以 $e^{\beta_i}$ 直接讀作「這個模型的率
# 提高一個 $e$ 倍，能把事件的勝算放大多少」——**這正是「它對
# 組合有多少貢獻」的自然度量**。取 $\tau=0$ 時，$\beta_i\le0$
# 的模型（率愈高反而愈不發生地震，無法解釋觀測）直接歸零；減去
# $e^\tau$ 是為了讓權重在門檻處**平滑**歸零，而不是從 $1$ 突然
# 跳到 $0$。最後用這組 $\pi_i$ 回頭做**加權平均**組合——於是
# 完全擺脫對 $\beta_0$ 的依賴，也拿回 19.5.1 節的 Jensen 保護。
#
# ### 19.7.3 遺忘視窗：權重該記多久
#
# 第二個自由度是**用多久的過去資料擬合**。兩種方案：**#1 用全部
# 歷史**（權重隨時間收斂、愈來愈鈍），**#2 只用最近一段**（權重
# 會「遺忘」，反映近期表現）。作者實測不同視窗長度，相對 SMA 的
# 累積資訊增益是：
#
# | 遺忘視窗 | 相對 SMA 的 CumIGPE |
# |---|---|
# | 半年 | 0.049 |
# | **1.5 年** | **0.097**（最佳） |
# | 2 年 | 0.074 |
# | 全部歷史 | $0.064 \pm 0.017$ |
#
# 下面用合成序列把「全歷史 vs 一年遺忘視窗」的差別跑出來。設定：
# 240 個空間格、8 年、每週一個預報視窗；模型甲是時間不變的長期
# 平滑模型（背景畫對了），模型乙是叢集模型（會追序列，但它自己
# 的背景畫歪了）。真實率平時等於甲的背景、序列爆發時轉向叢集
# 形狀。權重用**無截距的二元 logistic** 擬合，再套 19.7.2 的映射：

# %% tags=["hide-input"]
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
# 兩條線的性格完全不同。<strong>紅線（一年遺忘視窗）</strong>在平靜期把權重
# 幾乎全給模型甲（模型乙的背景畫歪了，$\beta_2$ 掉到零以下直接
# 被截掉），序列一爆發就在幾週內把權重整個交給模型乙，序列衰減
# 完又交還回去。<strong>藍線（全部歷史）</strong>則是一條無法回頭的階梯：每
# 經過一次序列就往上跳一階，然後永遠記得——它到第八年還掛在
# 半途，既不適合平靜期也不適合爆發期。
#
# 但紅線也暴露了遺忘機制的代價，值得盯著看：**權重轉移是落後
# 的**。序列在第 2.1 年爆發，紅線要到大約第 2.3 年才追上去——因為
# 一年視窗裡此刻只有兩個月的序列資料，訊號還被前十個月的平靜期
# 稀釋著。這正是 19.7.2 節「理由二」的圖像版：**能用來學權重的
# 資料，永遠是已經過去的資料。**
#
# 而視窗長度本身是一個自由度，也就是一個過度配適的風險來源。
# 作者很誠實地報告：文中多項改良「單獨用都有效、一起用就沒效」。
# **調參的自由度愈多，前瞻期崩塌的風險愈大**——19.6 節剛付過
# 這個學費。
#
# ### 19.7.4 用門檻以下的小地震學權重
#
# 最後一招最實用。OEF-Italy 的預報目標是 $M_L \ge 3.95$，
# 8993 個 $0.1° \times 0.1°$ 空間格、2005–2020 共 6227 個重疊
# 預報視窗，每個模型約**五千六百萬個樣本**——聽起來很多，但其中
# **99.994% 是無事件格**，有效樣本數（第 18 章的老規矩）是地震
# 顆數，只有幾百顆。
#
# 解法：**預報門檻維持 $M_L\ge3.95$，但用 $M_L \ge 2.95$ 的事件
# 來擬合權重**。低一個規模單位，事件數依 Gutenberg–Richter 律
# 大約多十倍，權重估計立刻靈敏得多。成績（相對比較的基準各異，
# 要分開讀）：
#
# | 組合 vs 基準 | 用 $M_L\ge3.95$ 擬合 | 用 $M_L\ge2.95$ 擬合 |
# |---|---|---|
# | 加權平均 vs SMA | $0.064\pm0.017$（全歷史） | $0.112\pm0.019$ |
# | 加權平均 vs **最佳單一模型** | $0.002\pm0.019$（不顯著） | $\mathbf{0.050\pm0.019}$ |
#
# 第二列才是重點。用目標門檻的事件擬合時，組合相對最佳單一模型
# （ETAS-LM）的增益是 $0.002\pm0.019$——**教科書等級的平手**；
# 改用小地震擬合權重之後變成 $0.050\pm0.019$，**顯著勝過最佳
# 單一模型**。而這是在「候選模型彼此很像」的不利條件下拿到的
# （三個都是統計式群聚模型，背景率決定法也雷同）。
#
# **同樣的資料、同樣的運算子，只換了權重的學法，結論從平手變成
# 顯著獲勝。** 這句話是本章的收束。
#
# ## 19.8 集成的副產品：認知不確定性
#
# 加權平均還有一個乘法給不了、而且常被忽略的好處：**它保留了
# 各成分預報的離散度**。
#
# 把組合視為一個混合分布：以機率 $\pi_i$ 選中成分 $i$，再由該
# 成分的 Poisson 分布抽出格內事件數 $\omega$。混合分布的一二階
# 動差是
#
# $$\begin{aligned}
# \mathbb{E}[\omega] &= \sum_i \pi_i\lambda_i = \bar\lambda, \\
# \mathbb{E}[\omega^2] &= \sum_i \pi_i\bigl(\lambda_i + \lambda_i^2\bigr),
# \end{aligned}$$
#
# 第二行用了 Poisson 的 $\mathbb{E}[\omega^2]=\lambda+\lambda^2$。
# 相減得到
#
# $$\mathrm{Var}[\omega]
# = \underbrace{\bar\lambda}_{\text{aleatory}}
# + \underbrace{\sum_i \pi_i\bigl(\lambda_i - \bar\lambda\bigr)^2}
# _{\text{epistemic}} .$$
#
# **總變異乾淨地拆成兩塊**：第一項是 Poisson 本身的隨機性（就算
# 我們知道真實的率也躲不掉），第二項是**成分之間的加權變異數**
# ——模型彼此不同意的程度，也就是**認知不確定性**。
#
# 兩個推論。**其一**，這一項正是第 17 章負二項分布的來源：組合
# 預報相對 Poisson 的**過度離散**，數值上就等於模型間的加權
# 變異數。理論與實作在這裡接上了。**其二**，乘法家族與任何
# 「collapse 成單一分布」的組合方式，交出的是一個 $\lambda_H$，
# 變異數就只剩 $\lambda_H$，**認知不確定性項恆等於零**——不是
# 它比較小，是這個資訊在建模的當下就被丟掉了。
#
# 所以集成的價值不只是「平均比較保險」。**它同時是一個量測工具：
# 量測我們對自己的預報有多沒把握。** 對作業型預報而言，這往往
# 比點估計更值錢。
#
# ## 19.9 常見誤解與陷阱
#
# **誤解 1：「格子有幾萬個，樣本很充足。」** 有效樣本數是**目標
# 地震顆數**，不是格子數。2014 年那批 hybrid 是用 31 個事件擬合
# 的，2009 年的混合模型用 152 個，OEF-Italy 五千六百萬個樣本裡
# 99.994% 是空格。IGPE 的分母是 $N$、標準誤的分母是 $\sqrt N$、
# Student-$t$ 的自由度是 $N-1$，**通通都是地震顆數**。這是第 18
# 章的規矩，在組合模型裡只會更嚴重，因為每多一個成分就多一個
# 參數。
#
# **誤解 2：「加法只能內插，所以組合分數不可能超過最好的成分。」**
# 逐格內插是對的，推論是錯的。19.3.1 節證明了兩個平手的模型
# 組合起來嚴格更好，19.5.1 節解釋了原因：分數是**對數**的和，
# Jensen 不等式只給下界，不給上界。
#
# **誤解 3：「組合模型的 $\ln L$ 比較高，所以它比較好。」**
# 19.4.4 節的退化保證讓這件事**在數學上必然發生**——參數空間
# 包含「組合就是基準」這個點。回溯期一律用 IGPEc，而且要記得
# 第 18 章的警告：**AICc 懲罰不足以救回過擬合**，回溯評估永遠
# 只是 sanity check。
#
# **誤解 4：「最弱的成分該踢掉。」** STEP 的靜態背景項單獨看是
# 全場最弱，拿掉它組合反而變差。有效成分的判準是「它講了別人
# 沒講的話」，不是「它自己講得多好」。
#
# **誤解 5：「回溯擬出的最佳權重可以直接拿去作業。」** 權重會
# 漂移，而且漂移的幅度與樣本數成反比。19.7.3 節的圖裡，同一組
# 模型的權重在八年內從 0 擺盪到 1——**權重是狀態變數，不是常數**。
#
# **誤解 6：「等權重最保險，不必學。」** 等權重也是一個選擇，
# 而且是一個假設「所有候選模型同樣有用」的強假設。它會把明顯
# 較差的模型硬拉進來稀釋——這正是 SMA 與最大化整體表現之間
# 那條分界線的另一面。
#
# **誤解 7：「logistic 迴歸擬得很好，直接拿它當組合預報就行。」**
# 實測是負增益。19.7.2 節的三個理由（自我強化偏誤、估計延遲、
# 稀有事件偏誤）都指向同一個結論：**該用的是係數，不是模型**。
#
# **誤解 8：「乘法比加法強，這是已經驗證的規律。」** 回溯上是，
# 前瞻上不是。而且 19.6.1 節的推導說明這不是運氣問題：乘法
# 結構在基準空間表現不佳時會**系統性地**放大錯誤。
#
# ## 19.10 研究前沿與未解問題
#
# **相關的用加法、獨立的用乘法。** Bayona 等人（2022）在結尾
# 提出這個方向：模型之間高度相關時用加法（享受 Jensen 的下界
# 保護），概念與資料來源真正獨立時才用乘法（賺外推的紅利）。
# 要強調的是作者的措辭是「**值得進一步研究**」——這是一個
# 假說，不是已確立的經驗法則。要把它變成法則，需要一個把
# 「模型間相關性」量化的定義，而這個定義本身還沒有共識：用預報
# 率的空間相關？用 $\ln\lambda$ 在目標地震上的相關？兩者給出的
# 答案不一樣。
#
# **貝氏模型平均與 M-open 的困境。** 把 $\pi_i$ 詮釋成「$f_i$ 是
# 真模型的機率」在數學上很誘人，但它預設真模型**在候選池裡**。
# 地震預報顯然不是這種情形——所有候選模型都是錯的。這種設定下
# 純貝氏的權重會隨資料累積而集中到單一模型上，反而失去分散化的
# 好處；如何在「真模型不在池裡」的前提下給出有原則的權重，是
# 統計學本身的開放問題。
#
# **多樣性天花板。** 這是本章最重要的但書。OEF-Italy 的作者自評
# 增益「相對溫和（IGPE 約 0.05）」，並直指原因：**三個候選模型
# 太像**。這與 2014 年「概念差異愈大、增益愈大」是跨越九年的
# 同一個結論。所以組合模型的成長空間不在演算法，在候選池：
# **統計模型 + 物理模型、目錄資料 + 非目錄資料。**
#
# **對台灣的具體啟示。** 三件事。**其一**，最想要的新血是
# **GNSS 應變率與活動斷層滑移率**——19.4.2 節說明了為什麼它們
# 可以直接當共軛輸入：框架只要求排序，不要求它是完整的預報模型。
# 台灣的 GNSS 網密度在全球名列前茅，這是現成的資產。**其二**，
# 台灣測試區面積小、目標事件少，19.4.3 節的懲罰算術會比加州更
# 兇——**參數預算必須更省**，19.3 節那種「只擬合一個權重」的
# 凸組合是最合理的第一步。**其三**，台灣的地震序列群聚性極強
# （車籠埔、池上、大埔），第 17 章的二元概似與負二項 N-test
# 幾乎是必需品，否則少數幾個多震格就會主宰所有比較。
#
# **還沒有好答案的問題。** 權重的估計不確定性怎麼傳遞到組合
# 預報的不確定性上（19.8 節的分解假設 $\pi_i$ 是已知常數，
# 實際上它們是估計量）？遺忘視窗的長度該不該隨活動度自適應，
# 還是這只是又一個過度配適的入口？以及一個更根本的：**組合模型
# 的權重本身是否應該送進 CSEP 凍結測試**——目前的慣例是凍結
# 模型、不凍結權重更新規則，而 19.7.3 節的圖說明權重更新規則
# 對結果的影響，可能比模型本身還大。
#
# ## 19.11 附錄：本章推導細節
#
# ### A. 凸組合的凹性與內點最優
#
# 令 $\lambda_r(n) = (1-r)\lambda_1(n)+r\lambda_2(n)$。對 $r$ 求
# 一階與二階導數：
#
# $$\begin{aligned}
# \frac{\mathrm{d}\ln L}{\mathrm{d}r}
# &= \sum_{n=1}^{N}\frac{\lambda_2(n)-\lambda_1(n)}{\lambda_r(n)}
#  - (\Lambda_2-\Lambda_1), \\
# \frac{\mathrm{d}^2\ln L}{\mathrm{d}r^2}
# &= -\sum_{n=1}^{N}\frac{\bigl[\lambda_2(n)-\lambda_1(n)\bigr]^2}
#    {\lambda_r(n)^2} \ \le\ 0 .
# \end{aligned}$$
#
# 二階導數處處非正，且只在 $\lambda_1\equiv\lambda_2$ 時為零，
# 所以 $\ln L$ 嚴格凹、最大值唯一。取 $\Lambda_1=\Lambda_2$ 並
# 代入 $r=0$ 與 $r=1$ 即得 19.3.1 節的兩個端點條件。至於 AM–GM
# 的部分：對正數 $\varrho_n$，
#
# $$\frac{1}{N}\sum_n \varrho_n \ \ge\
# \Bigl(\prod_n\varrho_n\Bigr)^{1/N}
# = \exp\Bigl(\frac{1}{N}\sum_n\ln\varrho_n\Bigr) = e^{\delta},$$
#
# 等號僅在所有 $\varrho_n$ 相等時成立。把 $\varrho_n$ 換成
# $1/\varrho_n$ 得到 $\overline{\varrho^{-1}} \ge e^{-\delta}$。
# $\delta=0$ 時兩者同時 $\ge 1$，且只要率比不是處處相等就是嚴格
# 大於，於是兩個端點導數分別嚴格為正、嚴格為負，最大值必在內部。
#
# ### B. 乘法 hybrid 的 KL 分解
#
# 沿用 19.6.1 節的記號：$w_j=\lambda_1(j)/\hat N_1$、
# $m(j)=\lambda_H(j)/\lambda_1(j)$、$q_j=w_j m(j)$。正規化條件
# $\hat N_H=\hat N_1$ 給出 $\sum_j q_j = \sum_j w_j m(j)=1$，
# 所以 $q$ 確實是機率分布。$\ln m(j) = \ln q_j - \ln w_j$，於是
# 對任意分布 $p$，
#
# $$\begin{aligned}
# \sum_j p_j\ln m(j)
# &= \sum_j p_j \ln\frac{q_j}{w_j}
#  = \sum_j p_j\ln\frac{p_j}{w_j} - \sum_j p_j\ln\frac{p_j}{q_j} \\
# &= D_{\mathrm{KL}}(p\|w) - D_{\mathrm{KL}}(p\|q).
# \end{aligned}$$
#
# 取 $p=w$ 得 $\sum_j w_j\ln m(j) = -D_{\mathrm{KL}}(w\|q)\le0$，
# 即 19.6.1 節的赤字。加法的有界性：凸組合
# $q_j = \sum_i\pi_i w^{(i)}_j \ge \pi_1 w_j$，代入
#
# $$D_{\mathrm{KL}}(p\|q) = \sum_j p_j\ln\frac{p_j}{q_j}
# \ \le\ \sum_j p_j\ln\frac{p_j}{\pi_1 w_j}
# = D_{\mathrm{KL}}(p\|w) + \ln\frac{1}{\pi_1},$$
#
# 只要保留給基準模型的權重 $\pi_1$ 不為零，損失就有上界。
# **乘法沒有這個下托，所以沒有這個上界。**
#
# ### C. 混合 Poisson 的變異數分解
#
# 這是全變異數定律的特例。令 $I$ 是取值 $i$ 機率 $\pi_i$ 的隱藏
# 指標，$\omega \mid I=i \sim \mathrm{Poisson}(\lambda_i)$：
#
# $$\begin{aligned}
# \mathrm{Var}[\omega]
# &= \mathbb{E}\bigl[\mathrm{Var}(\omega\mid I)\bigr]
#  + \mathrm{Var}\bigl[\mathbb{E}(\omega\mid I)\bigr] \\
# &= \sum_i\pi_i\lambda_i
#  + \sum_i\pi_i(\lambda_i-\bar\lambda)^2
#  = \bar\lambda + \mathrm{Var}_\pi(\lambda).
# \end{aligned}$$
#
# 第一項用了 Poisson 的變異數等於期望值。過度離散指標
# $\mathrm{Var}/\mathbb{E} = 1 + \mathrm{Var}_\pi(\lambda)/\bar\lambda$
# 恆 $\ge 1$，等號僅在所有成分同率時成立——**組合預報必然過度
# 離散**，第 17 章的負二項不是權宜之計，是結構上的必然。
#
# ### D. logistic 係數與乘冪組合
#
# 由 $\ln\frac{P}{1-P}=\beta_0+\sum_i\beta_i\ln\phi_i$ 兩邊取
# 指數即得 $\frac{P}{1-P}=e^{\beta_0}\prod_i\phi_i^{\beta_i}$。
# 把某個 $\phi_i$ 乘上 $c>0$，勝算就乘上 $c^{\beta_i}$；取
# $c=e$ 便得到 odds ratio $e^{\beta_i}$，這就是 19.7.2 節權重
# 映射的來源。兩個邊界值得記住：$\beta_i=0$ 時該模型對勝算完全
# 沒有影響，映射後權重為 $e^0-e^0=0$；$\beta_i<0$ 時模型的率
# 愈高、事件反而愈不發生，映射直接截成零。**權重的正性不是
# 人為約束，是這個映射的自然結果。**
#
# ---
#
# 這一章的四篇關鍵文獻剛好排成一個完整的科學迴圈：先發現互補性
# 帶來增益（2009），再發明更強的組合文法（2014），然後被前瞻
# 測試打回原形（2022），最後把教訓消化成更好的方法（2023）。
# 注意這個迴圈能轉起來的前提：**有一個獨立的、事先約定規則的
# 測試機制，願意讓漂亮的想法難堪。**
#
# 也注意這一章反覆出現的同一句話：**組合的天花板不在演算法，
# 在候選池的多樣性。** 而到目前為止，第二部的候選池全部來自
# 同一口井——地震目錄。ETAS 讀它的短期叢集，EEPAS 讀它的中期
# 前兆尺度，PPE 讀它的長期空間分布，三個模型爭論的是同一份
# 資料該怎麼解讀。要真正擴大候選池，必須引進**目錄以外的物理**。
#
# 最直接的一條路，是回到地震的成因本身：斷層上的應力累積與釋放
# 是一個有記憶的過程，「上一次錯動到現在過了多久」本身就是資訊，
# 而這個資訊在整份目錄的統計裡是看不見的。
# {doc}`第 20 章 <20_recurrence_models>`就來處理這一類模型：
# 更新過程、危害函數、複發間隔的分布，以及一個貫穿全章的老問題
# ——**地震到底是不是週期性的？**
