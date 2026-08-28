# %% [markdown]
# # 21. PSHA：從發生率到地動危害
#
# {doc}`第 20 章 <20_recurrence_models>`的最後留下一個沒答完的問題。那一章算得出
# 「這條斷層未來 50 年破裂的機率」，但工程師從來沒問過這個問題。工程師問的是：
# **這塊地在未來 50 年裡，地表加速度超過 0.3 g 的機率是多少**？房子不在乎震央在
# 哪、規模多大、是哪一條孕震構造破裂的——房子只在乎腳下的地面晃多大。
#
# 這兩個問題之間隔著一整套計算：把所有可能的震源、所有可能的規模與距離，以及「同樣的
# 地震在同樣的距離上，地動仍然會差好幾倍」這件事，全部連同各自的機率積分起來。這套
# 計算叫**機率式地震危害度分析**（Probabilistic Seismic Hazard Analysis, PSHA），骨架
# 由 Cornell (1968) 奠定，至今仍是全世界建築耐震規範、核電廠選址與地震保險費率背後的
# 引擎。它也是這門學問通往社會最寬的一座橋，而橋的每一塊磚都是前面某一章的產品：規模
# 分布是第 11 章的 GR 律，除叢的代價是第 12 章的破產點，短期率的暴衝是第 13、14 章的
# ETAS，斷層的長期記憶是第 20 章的 BPT。這一章要做的，是把它們**接進同一個積分**。
#
# 本章以 Baker 的教學白皮書（Baker 2013, *Probabilistic Seismic Hazard Analysis*,
# White Paper v2.0.1；已被 Baker, Bradley & Stafford 2021 的教科書 *Seismic Hazard and
# Risk Analysis* 取代）為骨架，擁有**危害積分**、**截斷 GR 分布**、**反聚合**與**非齊次
# Poisson 過程**；hazard function {eq}`eq:hazard-def`、BPT 密度 {eq}`eq:bpt-pdf`、更新
# 過程條件機率 {eq}`eq:renewal-prob` 沿用第 20 章，GR 密度 {eq}`eq:gr-density` 沿用
# 第 11 章，破產點 {eq}`eq:mx` 沿用第 12 章，一律引用不重推。
#
# 先預告本章最反直覺的一句話：**在高地動強度端，危害的主要來源不是「很大的地震」，
# 而是「地動預測方程的殘差剛好很大」**。它會在 21.4 節被推導出來，在 21.7 節被量化
# 成一個具體的百分比。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 21.1 決定論的兩個致命問題
#
# 在 PSHA 之前（以及至今在某些場合仍然並用）的作法是**決定論式地震危害度分析**
# （DSHA）：挑一個「最壞情境」——通常是最近的那條斷層發生它能發生的最大地震——算出
# 這個情境下的地動，照著設計。聽起來比機率保守、穩健、也好溝通。它有兩個問題，而且
# 兩個都致命。
#
# ### 問題一：「最壞情境」沒有客觀定義
#
# 把「最壞」寫成可執行的規則，你立刻需要做一連串沒有標準答案的選擇：**哪一條斷層**
# （最近的那條規模上限可能很小，遠一點的可能大得多，哪個組合的地動大要算過才知道）；
# **規模取多少**（$m_{\max}$ 本身是由斷層長度與尺度律推出的估計值，不確定度常達 0.3
# 個規模單位）；**距離取多少**（破裂面是一個面不是一個點，取最近距離、破裂中心或投影
# 距離結果都不同）；以及最要命的**地動取哪一個**——給定規模與距離，地動預測方程給的
# 是一個**分布**，取中位數、中位數加一個標準差、還是加兩個？21.4 節會算出一個標準差
# 就是約 1.8 倍的地動，這個決定對設計值的影響比「選哪條斷層」還大。
#
# 每一個選擇都合理、都可辯護，而不同組合可以差好幾倍。「最壞情境」**不是從資料讀出來
# 的事實，是一串主觀選擇的乘積**。 DSHA 的問題不在於保守，而在於它把主觀性藏進一個
# 看起來很確定的數字裡。
#
# ### 問題二：完全忽略發生率
#
# 更根本的一項。DSHA 從頭到尾沒問過「這個情境多久發生一次」，於是一條一萬年動一次的
# 大斷層，與一條一百年動一次的小斷層平起平坐——只要前者算出來的地動大，它就決定了設計
# 值。但對一棟設計壽命 50 年的建築，兩者的意義天差地遠。**沒有發生率，就無法回答任何
# 與時間有關的問題**，而所有工程決策——設計壽命、保險費率、補強的投資報酬——本質上都
# 與時間有關。PSHA 的回答很直接：**不挑情境**，把所有可能的地震與地動連同各自的發生率
# 全部積分起來，輸出一條「地動強度 vs 年超越率」的曲線；挑不挑情境留給使用者在曲線上
# 自己決定，而且是在知道每個選擇對應多少年率的情況下決定。

# %% [markdown]
# ## 21.2 Baker 五步驟與危害積分
#
# Baker (2013) 把 PSHA 拆成**五個步驟**（原文用字是 "PSHA is composed of five
# steps"）：(1) 找出所有可能產生破壞性地動的**震源**；(2) 描述每個震源的**規模分布**；
# (3) 描述震源到場址的**距離分布**；(4) 給定規模與距離，預測**地動強度的分布**；
# (5) 用**全機率定理**把上述不確定性合併起來。
#
# 文獻上也常見「四步驟」的說法（震源幾何 → 規模–頻率分布 → GMPE → 危害積分），那是
# 把第 1、3 步併成「震源幾何與距離分布」的濃縮版。兩種講法都對，但**不要把四步驟說成
# 是 Baker 的**。以下把第五步完整推一遍。
#
# ### 從單一震源的一次地震開始
#
# 固定一個震源 $\text{src}$。它以年率 $\nu_{\rm src}$ 產生規模大於等於 $m_{\min}$ 的
# 地震（$m_{\min}$ 是工程上認定「可能造成結構破壞」的下限，實務上取 4.5 到 5.0，21.5
# 節會檢驗這個選擇的影響），每一次地震帶著兩個隨機**標記**：規模 $M$（密度 $f_M(m)$，
# 步驟 2）與源–站距離 $R$（密度 $f_R(r)$，步驟 3）。先問一個較簡單的問題：**這個震源的
# 下一次地震，讓場址的地動強度 $IM$ 超過 $x$ 的機率是多少**？步驟 4 給的是條件機率
# $P(IM>x\mid m,r)$，用全機率定理對 $(M,R)$ 積分：
#
# $$\begin{aligned}
# P(IM>x)
#   &= \int\!\!\int P(IM>x \mid M=m, R=r)\,f_{M,R}(m,r)\,\mathrm{d}m\,\mathrm{d}r \\
#   &= \int\!\!\int P(IM>x \mid m,r)\,f_M(m)\,f_R(r)\,\mathrm{d}m\,\mathrm{d}r .
# \end{aligned}$$
#
# 第二個等號用掉一個**假設**：$M$ 與 $R$ 獨立。這在點震源近似下成立，但對延伸破裂面
# 並不成立——大地震的破裂面比較長，離場址的最近距離統計上比較小；嚴格作法是寫成條件
# 密度 $f_R(r\mid m)$，把破裂尺度律放進去。本章的玩具模型用點震源，所以維持獨立寫法。
#
# ### 從機率到年率
#
# 震源的地震序列假設為速率 $\nu_{\rm src}$ 的 Poisson 過程，每一次地震以機率
# $P(IM>x)$ 造成超越。**Poisson 過程的獨立稀疏化仍是 Poisson 過程**（附錄 C），
# 速率為原速率乘上保留機率，所以「這個震源造成超越」本身也是 Poisson 過程，年率
# $\lambda_{\rm src}(x)=\nu_{\rm src}\,P(IM>x)$。最後把所有震源加起來——**獨立
# Poisson 過程疊加仍是 Poisson，速率相加**——於是
#
# $$\lambda_{IM}(x) \;\equiv\; \lambda(IM>x)
#   \;=\; \sum_{\rm src}\nu_{\rm src}
#   \int\!\!\int P(IM>x\mid m,r)\,f_M(m)\,f_R(r)\,\mathrm{d}m\,\mathrm{d}r$$ (eq:psha)
#
# 這就是 **PSHA 主方程**（危害積分），本章其餘所有內容都是它的註腳。右邊四個東西恰好
# 對應前四個步驟：$\sum_{\rm src}$ 是步驟 1、$f_M$ 是步驟 2、$f_R$ 是步驟 3、
# $P(IM>x\mid m,r)$ 是步驟 4，而「全部乘起來積分」就是步驟 5。
#
# 有四件事值得看清楚。**左邊是率不是機率**，單位是「次／年」，可以大於 1，轉成機率需要
# 額外假設（21.6 節）。**超越率是 $x$ 的遞減函數**，把 $x$ 從小掃到大畫出來就是**危害
# 曲線**。**這條式子沒有時間**：$\nu_{\rm src}$ 是常數，今年與明年一模一樣，這是傳統
# PSHA 的**時不變**性質，也是 21.8 節要動手術的地方。**兩種不確定性被混在一起了**：
# $f_M$、$f_R$ 與 GMPE 殘差是「即使模型正確結果仍然隨機」，而 $\nu_{\rm src}$ 該取多少、
# GMPE 該選哪一條則是「我們不知道真相」——21.11 節會把這兩種分開。

# %% [markdown]
# ## 21.3 步驟 2：截斷 GR 分布
#
# {eq}`eq:psha` 裡的 $f_M(m)$ 從哪來？第 11 章已經證明過，GR 律等價於「規模服從指數
# 分布」，密度是 {eq}`eq:gr-density`：$s(m)=\beta e^{-\beta(m-m_0)}$，$\beta=b\ln 10$。
#
# 但這個分布有一個 PSHA 不能接受的性質：**它沒有上界**。指數分布給 $M=12$ 一個很小但
# 非零的密度，而地球上不存在規模 12 的地震——斷層的長度、寬度與可用的應變能都有物理
# 上限。把無上界的分布餵進 {eq}`eq:psha`，高強度端會被一批不可能的地震撐起來。
#
# ### 推導截斷指數密度
#
# 處理方式是**條件化**：把分布限制在 $[m_{\min}, m_{\max}]$ 上，問「已知規模落在區間
# 內，它的密度長什麼樣」。直接套條件密度的定義：
#
# $$\begin{aligned}
# f_M(m)
#   &= \frac{s(m)}{P\bigl(m_{\min}\le M\le m_{\max}\bigr)}
#      \qquad (m_{\min}\le m\le m_{\max}) \\[2pt]
#   &= \frac{\beta e^{-\beta(m-m_{\min})}}
#           {\displaystyle\int_{m_{\min}}^{m_{\max}}
#            \beta e^{-\beta(u-m_{\min})}\,\mathrm{d}u}
#    = \frac{\beta e^{-\beta(m-m_{\min})}}
#           {\Bigl[-e^{-\beta(u-m_{\min})}\Bigr]_{m_{\min}}^{m_{\max}}} \\[2pt]
#   &= \frac{\beta\,e^{-\beta(m-m_{\min})}}{1-e^{-\beta(m_{\max}-m_{\min})}} .
# \end{aligned}$$ (eq:trunc-gr)
#
# 分母 $1-e^{-\beta(m_{\max}-m_{\min})}$ 就是「未截斷的指數分布落在區間內的機率質量」，
# 除掉它等於把被砍掉的尾巴**按比例重新分配**回區間內。這個重分配不是無害的——21.5 節
# 會看到它留下一個不到一成的痕跡。檢查兩個極限：$m_{\max}\to\infty$ 時分母趨於 1，退回
# {eq}`eq:gr-density`；$m_{\max}\to m_{\min}$ 時分母 $\approx\beta(m_{\max}-m_{\min})$，
# 密度趨於區間上的均勻分布——區間窄到看不出指數的彎曲。
#
# ### 換底：率與門檻的一致性
#
# 一個實務上常出錯的細節。震源的活動度通常以「$M\ge m_{\rm ref}$ 的年率」報告（例如
# $m_{\rm ref}=5.0$），但積分要用的是「$M\ge m_{\min}$ 的年率」。換算直接來自 GR 律的
# 累積形式：$\nu(m_{\min}) = \nu(m_{\rm ref})\cdot 10^{\,b\,(m_{\rm ref}-m_{\min})}$。
# 把 $m_{\min}$ 從 5.0 降到 4.5（$b=1$）會讓事件數變成 $10^{0.5}\approx3.16$ 倍。
# **這批多出來的地震幾乎全是小地震，對高強度端沒有貢獻**——這句話聽起來理所當然，
# 但它是可以量化檢驗的宣稱，21.5 節就檢驗它。
#
# ### 一個必須誠實面對的裂縫
#
# {eq}`eq:trunc-gr` 只適用於**面震源**與背景地震活動：一片區域內的地震規模分布是連續
# 的冪次。但對**斷層源**（孕震構造），第 20 章的更新過程框架建立在**特徵地震**假設上
# ——同一段斷層反覆產生規模幾乎相同的破裂，才談得上「複發間隔」。兩個模型不相容，
# 實務上的處理是混合式的：斷層源在特徵規模附近放一個窄峰、在其下接一段截斷 GR，權重
# 由滑移率平衡決定；接縫處會出現一個**凸起**，而它是真實現象還是模型拼接的產物至今
# 沒有定論（第 20 章 20.9 節誤解 3）。**這是整套計算裡最主觀的一步。**

# %% [markdown]
# ## 21.4 步驟 4：地動預測方程與 $\varepsilon$
#
# 現在處理 {eq}`eq:psha` 裡最後一項，也是最重要的一項：$P(IM>x\mid m,r)$。
#
# **地動預測方程**（ground motion prediction equation, GMPE；舊稱衰減律）是一條經驗
# 迴歸式：
#
# $$\ln IM = \overline{\ln IM}(m,r,\boldsymbol\theta)
#   \;+\; \sigma_{\ln IM}\,\varepsilon, \qquad \varepsilon\sim N(0,1),$$
#
# 其中 $\overline{\ln IM}$ 是規模、距離與其他預測變數 $\boldsymbol\theta$（震源機制、
# 場址 $V_{S30}$、上／下盤效應……）的函數，由數千到數萬筆強震紀錄迴歸而得。
#
# 為什麼取**對數常態**？地動強度是正的而常態分布會給出負值，取對數把值域搬到整條實線；
# 更重要的是地動的物理是**乘性**的——震源輻射、幾何擴散、非彈性衰減、場址放大，每一層
# 都是對前一層乘上一個因子，乘性效應的對數是加性的，中央極限定理因此作用在 $\ln IM$ 上
# 而不是 $IM$ 上。經驗上強震紀錄的殘差直方圖也確實接近對數常態。
#
# ### 超越機率
#
# 有了這個假設，$P(IM>x\mid m,r)$ 就是一行事：
#
# $$\begin{aligned}
# P(IM>x\mid m,r)
#   &= P\bigl(\ln IM > \ln x \mid m,r\bigr)
#    = P\!\left(\varepsilon >
#      \frac{\ln x - \overline{\ln IM}(m,r)}{\sigma_{\ln IM}}\right) \\
#   &= 1-\Phi\!\left(\frac{\ln x - \overline{\ln IM}(m,r)}{\sigma_{\ln IM}}\right).
# \end{aligned}$$ (eq:gmpe-exceed)
#
# $\Phi$ 是標準常態的累積分布函數。把 {eq}`eq:gmpe-exceed` 與 {eq}`eq:trunc-gr` 一起
# 代進 {eq}`eq:psha`，危害積分就完全可算了。
#
# ### $\varepsilon$ 是什麼
#
# 中間那個量本身有名字。定義
# $\varepsilon=[\ln IM-\overline{\ln IM}(m,r)]/\sigma_{\ln IM}$——**這一次地震在這個
# 場址產生的地動，比同規模同距離的平均值高出幾個標準差**。$\varepsilon=0$ 是「完全
# 平均」的一次地震，$\varepsilon=2$ 是「運氣很差」的一次。它在 21.7 節會成為反聚合的
# 第三個維度，在工程選波時則是挑選設計地震歷時的關鍵指標。
#
# 順帶一提，$\sigma_{\ln IM}$ 在現代 GMPE 裡還會再拆成**事件間**（同一次地震所有測站
# 共有的偏差，來自震源）與**事件內**（不同測站之間的差異，來自路徑與場址）兩部分；
# 這對多場址相關性分析很重要，但單場址的 {eq}`eq:psha` 只看得到兩者的總和。
#
# ### 一個標準差就是 1.8 倍
#
# 現在講本章最反直覺的數字。典型的 $\sigma_{\ln IM}$ 落在 **0.5 到 0.6** 個自然對數
# 單位。換回線性尺度：
#
# $$\frac{IM(\varepsilon=1)}{IM(\varepsilon=0)}
#   = \frac{e^{\overline{\ln IM}+\sigma_{\ln IM}}}{e^{\overline{\ln IM}}}
#   = e^{\sigma_{\ln IM}} \approx e^{0.57} \approx 1.77 .$$
#
# **同樣的地震、同樣的距離、同樣的場址，一個標準差就是約 1.8 倍的地動；兩個標準差
# 就是約 3.1 倍。** 而 $\pm1\sigma$ 只涵蓋 68% 的情形——「M6.5、距離 10 km」這個條件，
# 只能把 PGA 框在一個上下差三倍多的區間裡。
#
# 對照震源那一端：$f_M$ 與 $f_R$ 當然也很寬，但它們是「哪一種地震會發生」的不確定性，
# 而且我們對它們有相當的掌握（GR 律、斷層幾何）。**地動這一層的不確定性，是即使把
# 地震本身完全指定清楚也消不掉的**。
#
# 這件事有一個直接後果，也是整章最重要的一句話：**在高強度端，危害主要不是由「規模很
# 大的地震」貢獻的，而是由「規模中等但 $\varepsilon$ 很大的地震」貢獻的**——中等規模
# 地震的發生率高出好幾個數量級，即使需要 $\varepsilon=2$ 或 3 才能達到目標地動，它們的
# 「率 × 機率」乘積仍可以贏過稀有的大地震。21.7 節會把它量化。
#
# ### 把三層畫出來
#
# {eq}`eq:psha` 的被積函數是三個東西的乘積：震源年率 $\nu_{\rm src}$（常數）、規模密度
# $f_M(m)$、超越機率 $P(IM>x\mid m,r)$。三層各自的形狀完全不同，而它們的乘積決定了
# 危害。用玩具模型畫出來——**一個場址，附近兩條斷層**：近而小（10 km 外，
# $m_{\max}=6.0$，$M\ge5$ 年率 0.10）與遠而大（25 km 外，$m_{\max}=8.0$，年率 0.05）。
# GMPE 用 Baker 白皮書示範的 Cornell 型式（係數見附錄 D）。

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import invgauss, norm

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

# 玩具模型：名稱、距離（km）、M≥5 年率（1/yr）、最大規模
FAULTS = [("近而小的斷層", 10.0, 0.10, 6.0),
          ("遠而大的斷層", 25.0, 0.05, 8.0)]
BETA = np.log(10.0)          # β = b·ln10，取 b = 1
M_REF, M_MIN, SIGMA = 5.0, 5.0, 0.57


def gmpe_ln(m, r):
    """Cornell 型 GMPE 的 ln(PGA) 平均值（PGA 單位 g、r 單位 km）。"""
    return -0.152 + 0.859 * m - 1.803 * np.log(r + 25.0)


def mag_pdf(m, m_max, m_min=M_MIN):
    """截斷 GR 密度 f_M(m)，即式 (eq:trunc-gr)。"""
    return BETA * np.exp(-BETA * (m - m_min)) / (
        1.0 - np.exp(-BETA * (m_max - m_min)))


def eps0(x, m, r):
    """達到 IM = x 所需的 ε（標準化殘差門檻）。"""
    return (np.log(x) - gmpe_ln(m, r)) / SIGMA


def p_exceed(x, m, r):
    """P(IM > x | m, r)，即式 (eq:gmpe-exceed)。"""
    return 1.0 - norm.cdf(eps0(x, m, r))


def source_nu(nu_ref, m_min=M_MIN):
    """把 M≥M_REF 的年率換底成 M≥m_min 的年率（固定 GR 的 a 值）。"""
    return nu_ref * 10.0 ** (M_REF - m_min)


def source_rate(x, r, nu_ref, m_max, m_min=M_MIN, n=600):
    """單一震源對 λ_IM(x) 的貢獻（對規模數值積分）。"""
    m = np.linspace(m_min, m_max, n)
    integ = p_exceed(x, m, r) * mag_pdf(m, m_max, m_min)
    return source_nu(nu_ref, m_min) * np.trapezoid(integ, m)


def hazard(x, m_min=M_MIN):
    """總危害 λ_IM(x)：兩個震源相加。"""
    return sum(source_rate(x, r, nu, mx, m_min) for _, r, nu, mx in FAULTS)


X_DEMO = 0.5                                        # 示範用的地動門檻（g）

fig = make_subplots(
    rows=1, cols=3,
    subplot_titles=("① 規模密度 f_M(m)（截斷 GR）",
                    f"② 超越機率 P(IM>{X_DEMO} g | m, r)",
                    "③ 被積函數 ν·f_M·P（面積 = λ）"))
areas = []
for i, (name, r, nu, mx) in enumerate(FAULTS):
    m = np.linspace(M_MIN, mx, 400)
    f_m, p_x = mag_pdf(m, mx), p_exceed(X_DEMO, m, r)
    integ = source_nu(nu) * f_m * p_x
    areas.append(np.trapezoid(integ, m))
    for col, y in [(1, f_m), (2, p_x), (3, integ)]:
        fig.add_trace(go.Scatter(x=m, y=y, mode="lines", name=name,
                                 legendgroup=name, showlegend=(col == 1),
                                 fill="tozeroy" if col == 3 else None,
                                 line=dict(color=PALETTE[i], width=2.2)),
                      row=1, col=col)
for c in (1, 2, 3):
    fig.update_xaxes(title_text="規模 m", row=1, col=c)
fig.update_yaxes(title_text="f_M(m)（1/規模）", row=1, col=1)
fig.update_yaxes(title_text="超越機率", type="log", row=1, col=2)
fig.update_yaxes(title_text="率密度（次/年/規模）", type="log", row=1, col=3)
apply_layout(fig, height=400, hovermode="x", title=(
    f"圖 1　危害積分的三層分解（x = {X_DEMO} g）：兩塊面積 {areas[0]:.2e}"
    f" + {areas[1]:.2e} = {sum(areas):.2e} 次/年"))
fig

# %% [markdown]
# 三個面板正好對應 {eq}`eq:psha` 的三層，而形狀彼此拮抗。**面板 ①** 的規模密度隨 $m$
# 指數遞減（GR 律）；兩條斷層的起點不同，是因為截斷正規化的分母不同——近斷層只能到
# 6.0，被截掉的尾巴多，剩下的質量被抬高。**面板 ②** 的超越機率隨 $m$ 急遽遞增（縱軸是
# 對數）：規模每多 0.5 級，達到 0.5 g 所需的 $\varepsilon$ 就少一大截；遠斷層整條低於
# 近斷層，因為距離衰減把它的平均地動壓下去了。**面板 ③** 是兩者相乘——**遞減乘上遞增，
# 得到一個有峰的函數**，這就是「主控規模」的來源，曲線下面積即該震源對危害的貢獻。
#
# 面板 ③ 已經預告了 21.7 節的全部內容：**危害不是由最大的地震貢獻的，也不是由最常見
# 的地震貢獻的，而是由那個乘積最大的區間貢獻的**，而那個區間會隨 $x$ 移動。

# %% [markdown]
# ## 21.5 危害曲線的讀法
#
# 現在把 $x$ 從小掃到大，每個 $x$ 都做一次圖 1 的積分，畫出危害曲線。這是 PSHA 的
# 最終輸出，也是全世界耐震規範背後的那張圖。

# %% tags=["hide-input"]
x_grid = np.logspace(-2, 0.3, 90)
lam_each = {name: np.array([source_rate(x, r, nu, mx) for x in x_grid])
            for name, r, nu, mx in FAULTS}
lam_tot = sum(lam_each.values())

LAM_475 = 1.0 / 475.0
pga_475 = np.exp(np.interp(np.log(LAM_475), np.log(lam_tot[::-1]),
                           np.log(x_grid[::-1])))
frac_far = lam_each[FAULTS[1][0]] / lam_tot
i_cross = int(np.argmax(frac_far > 0.5))

fig = go.Figure()
for i, (name, _, _, _) in enumerate(FAULTS):
    fig.add_trace(go.Scatter(x=x_grid, y=lam_each[name], mode="lines",
                             name=name,
                             line=dict(color=PALETTE[i], width=2, dash="dot")))
fig.add_trace(go.Scatter(x=x_grid, y=lam_tot, mode="lines", name="總危害",
                         line=dict(color=ACCENT, width=3)))
fig.add_hline(y=LAM_475, line=dict(color=QUAKE_COLOR, width=1.4, dash="dash"),
              annotation_text="λ = 1/475 /yr（50 年 10%）",
              annotation_position="bottom left")
fig.add_vline(x=pga_475, line=dict(color=QUAKE_COLOR, width=1, dash="dot"))
fig.add_annotation(x=np.log10(x_grid[i_cross]),
                   y=np.log10(lam_tot[i_cross]), text="遠斷層反超",
                   showarrow=True, arrowhead=2, ax=-45, ay=-30,
                   font=dict(size=11))
fig.update_xaxes(title_text="地動強度 x：PGA（g）", type="log")
fig.update_yaxes(title_text="年超越率 λ_IM(x)（次/年）", type="log")
apply_layout(fig, height=440, hovermode="x", title=(
    f"圖 2　玩具 PSHA 危害曲線：1/475 對應 PGA = {pga_475:.3f} g，"
    f"遠斷層在 {x_grid[i_cross]:.2f} g 之上反超"))
fig

# %% [markdown]
# **橫著讀**：選一個年率（規範常用 1/475 或 1/2475），水平線與總危害曲線的交點就是設計
# 地動，玩具模型的數字印在圖標題裡。**豎著讀**：選一個地動值讀出它的年率——這在既有
# 結構評估裡更常用（「這棟樓撐得住 0.3 g，那 0.3 g 多久超越一次？」）。
#
# **看兩條虛線的相對位置**，這是危害曲線最有教育意義的部分。低強度端幾乎完全由近而小
# 的斷層主導：它常動，那個強度對它是家常便飯。往右走，近斷層的曲線掉得愈來愈快——
# $m_{\max}=6.0$ 是一道硬牆，超過某個強度之後它只能靠愈來愈極端的 $\varepsilon$ 硬撐；
# 遠而大的斷層掉得慢得多，最後在圖上標出的那個 PGA 之上**反超**。這個交叉是普遍現象
# 而非巧合：**每一條震源在危害曲線上都有自己的「勢力範圍」，邊界由 $m_{\max}$ 與距離
# 共同決定**。這正是 21.7 節要處理的問題。
#
# ### $m_{\min}$ 敏感度：小地震到底重不重要
#
# 21.3 節說「$m_{\min}$ 從 5.0 降到 4.5 會多出兩倍多的地震，但它們對高強度端沒有
# 貢獻」。直接檢驗：固定 GR 的 $a$ 值（也就是固定 $M\ge5$ 的年率），只改 $m_{\min}$，
# 看三條危害曲線差多少。

# %% tags=["hide-input"]
x_g2 = np.logspace(-2, 0.2, 80)
fig = go.Figure()
curves = {}
for i, mmin in enumerate([4.5, 5.0, 5.5]):
    y = np.array([hazard(x, mmin) for x in x_g2])
    curves[mmin] = y
    fig.add_trace(go.Scatter(x=x_g2, y=y, mode="lines",
                             name=f"m_min = {mmin:.1f}",
                             line=dict(color=PALETTE[i], width=2.4)))
lo = np.array([curves[m][0] for m in (4.5, 5.0, 5.5)])
hi = np.array([curves[m][-1] for m in (4.5, 5.0, 5.5)])
fig.add_vline(x=x_g2[-1], line=dict(color="#bbb", width=1))
fig.update_xaxes(title_text="PGA（g）", type="log")
fig.update_yaxes(title_text="年超越率 λ_IM(x)（次/年）", type="log")
apply_layout(fig, height=430, hovermode="x", title=(
    f"圖 3　m_min 的影響：低強度端（{x_g2[0]:.2f} g）三者差 "
    f"{lo.max()/lo.min():.1f} 倍，高強度端（{x_g2[-1]:.2f} g）"
    f"只差 {hi.max()/hi.min():.2f} 倍"))
fig

# %% [markdown]
# 結果非常乾淨。低強度端三條曲線分得很開——把 $m_{\min}$ 從 5.5 降到 4.5，年超越率
# 差了將近一個數量級，因為那裡的危害幾乎全由小地震包辦。往右走三條線迅速收攏，到圖的
# 右端只剩下圖標題印出的那個倍數。
#
# 而那個殘餘差距**不是**小地震的貢獻，是 {eq}`eq:trunc-gr` 的正規化分母造成的：
# $m_{\min}$ 改變會改變 $1-e^{-\beta(m_{\max}-m_{\min})}$，等於把「被 $m_{\max}$ 截掉
# 的機率質量」以不同比例重分配回區間內，連帶輕微改動大規模端的密度。**真正由 $M<5$
# 地震直接貢獻的高強度危害，比這個數字還小。** 所以實務上 $m_{\min}$ 取 4.5 到 5.0
# 是有根據的：更小的地震「不被認為有破壞結構的能力」，略掉它們不是為了省計算量，而是
# 因為它們真的不影響答案的那一端。但反過來也要記住，**如果你關心的是低強度端**（設備
# 功能性、非結構構件，或把危害曲線接到損失曲線去算年化期望損失），$m_{\min}$ 的選擇
# 就會直接改變你的答案。

# %% [markdown]
# ## 21.6 回歸期不是週期
#
# 危害曲線的縱軸是年率，但所有規範、新聞與溝通場合講的都是「475 年回歸期」或
# 「50 年 10% 超越機率」。這兩個說法之間的橋只有一條式子，而它是全工程界被誤解最深的
# 一段。
#
# **回歸期**（return period）定義為發生率的倒數：$T_R=1/\lambda_{IM}(x)$。就這樣，
# 沒有別的。年率 0.01 對應回歸期 100 年。Baker 特別強調它的正確全名是 **mean return
# period**——**平均**回歸期，因為它就是超越事件之間平均間隔的長度。
#
# 記號提醒：本章的 $T_R$ 是**地動超越**的平均間隔，第 20 章的 $T_r$ 是**單一斷層破裂**
# 的平均複發時間。兩者形式相似、意義完全不同；一個場址的地動可以被好幾條斷層超越，
# 所以 $T_R$ 通常遠短於任何一條斷層的 $T_r$。
#
# ### 從率到時窗機率
#
# 要把年率換成「未來 $t$ 年內至少超越一次的機率」，必須額外假設超越事件在時間上的
# 分布。21.2 節已經論證過超越事件構成速率 $\lambda_{IM}$ 的 Poisson 過程，所以直接
# 寫得出來。令 $N(t)$ 為 $[0,t]$ 內的超越次數：
#
# $$\begin{aligned}
# P\bigl(N(t)=k\bigr) &= \frac{(\lambda_{IM}t)^{k}}{k!}\,e^{-\lambda_{IM}t}, \\
# P\bigl(\text{at least one}\bigr) &= 1-P\bigl(N(t)=0\bigr)
#   = 1-e^{-\lambda_{IM}t}.
# \end{aligned}$$
#
# 反解：給定時窗 $t$ 與目標機率 $P$，得 $\lambda_{IM}=-\ln(1-P)/t$、
# $T_R=t/[-\ln(1-P)]$。代入工程上最常見的一組數字（$P=0.10$、$t=50$ 年）：
#
# $$T_R = \frac{50}{-\ln 0.9} = \frac{50}{0.10536} = 474.6 \approx 475 .$$
#
# 「50 年 10% 超越機率 ≈ 475 年回歸期」這句話的全部由來，**就是這一行**。同樣地，
# $P=0.02$、$t=50$ 年給出 2475 年（另一個規範常用值）。這兩個數字不是從地質資料算出
# 來的，是從 Poisson 假設換算出來的。
#
# Baker 給的三條「為什麼用 Poisson」的理由是：數學形式簡單、多數情況與觀測相符、更
# 複雜的模型通常不顯著改變最終結果。第三條值得咀嚼——它不是說 Poisson 是對的，而是說
# **在目前的不確定度下，換模型帶來的改變被其他不確定性淹沒了**（第 20 章 20.4 節量化
# 過：當「距上次事件的時間」約在半個複發期附近，BPT 與 Poisson 的差異小到不影響決策）。
# 但這句辯護有有效範圍——對剛破裂完或明顯逾期的斷層，差異就是一階效應，21.8 節會畫。
#
# ```{admonition} 「475 年才來一次，我這輩子遇不到」
# :class: warning
# 這句話錯三次。**第一，回歸期是平均間隔，不是週期**——Poisson 過程沒有記憶，
# 「明年就超越」的機率每一年都是 $1-e^{-1/475}\approx0.21\%$，不會因為剛超越過就變低。
# **第二，機率並不小**：50 年裡至少遇上一次是 10%，而房子的壽命通常不只 50 年。
# **第三，475 年講的不是「大地震發生」的間隔，是「這個場址的地動超越設計值」的間隔**
# ——一個場址的地動可以被許多不同斷層、許多不同規模的地震超越，這正是 {eq}`eq:psha`
# 那個求和號的意思。
# ```
#
# Baker 本人的建議更乾脆：**乾脆只報年率，不要報回歸期**，免得「期」這個字誘導出週期
# 性的聯想。這是一個難得的、由原作者親口說出的溝通建議。

# %% [markdown]
# ## 21.7 反聚合：危害是誰貢獻的
#
# {eq}`eq:psha` 把所有情境加總，這是 PSHA 相對於 DSHA 的核心優點，但同時是代價：
# **算完之後，「這個危害水準主要來自哪一種地震」反而看不到了**——而工程師需要知道，
# 因為選設計地震歷時、做情境演練、判斷長短週期需求，都得有一個具體的地震當對象。
# **反聚合**（deaggregation；disaggregation 亦通用）就是把積分拆回去。
#
# ### 由貝氏定理推導
#
# 問題是：**已知地動超越了 $x$，造成這次超越的地震規模是 $m$ 的機率是多少**？這是一個
# 條件機率，分母是已知的事件「$IM>x$」，直接套定義：
#
# $$P(M=m \mid IM>x) \;=\;
#   \frac{\lambda\bigl(IM>x,\; M=m\bigr)}{\lambda\bigl(IM>x\bigr)} .$$
#
# 分母就是主方程 {eq}`eq:psha`；**分子只是「不對 $m$ 積分」的同一個式子**。寫成連續
# 密度的形式：
#
# $$f_{M\mid IM>x}(m) = \frac{\displaystyle\sum_{\rm src}\nu_{\rm src}\,f_M(m)
#   \int P(IM>x\mid m,r)\,f_R(r)\,\mathrm{d}r}{\lambda_{IM}(x)} .$$
#
# 分母的存在保證它積分為 1。這裡有一個容易滑過去的要點：上式的左邊是
# **機率**，右邊是兩個**率**相除——量綱在相除時消掉了，這正是為什麼反聚合可以在
# 「率」的世界裡做貝氏更新而不必先轉成機率。
#
# ### 推廣到三個維度
#
# 同樣的道理可以不對 $r$ 積分，得到聯合反聚合 $f_{M,R\mid IM>x}(m,r)$；還可以再往下拆
# 一層。注意 {eq}`eq:gmpe-exceed` 裡的 $P(IM>x\mid m,r)$ 其實是 $\varepsilon$ 的一個尾
# 機率 $\int_{\varepsilon_0}^{\infty}\varphi(e)\,\mathrm{d}e$，其中
# $\varepsilon_0(x,m,r)=[\ln x-\overline{\ln IM}(m,r)]/\sigma_{\ln IM}$。把這個積分也
# 攤開，就得到**三維反聚合**：
#
# $$f_{M,R,\varepsilon\mid IM>x}(m,r,e)
#   = \frac{\displaystyle\sum_{\rm src}\nu_{\rm src}\,
#           f_M(m)\,f_R(r)\,\varphi(e)\,
#           \mathbf{1}\{e>\varepsilon_0(x,m,r)\}}{\lambda_{IM}(x)} ,$$
#
# $\varphi$ 是標準常態密度、$\mathbf 1\{\cdot\}$ 是指示函數。三個維度各回答一個問題：
# $M$ 說「多大的地震」、$R$ 說「多遠」、$\varepsilon$ 說「地動要多不走運」。USGS 的
# 標準反聚合輸出就是把這三個量的邊際分布與平均值一起報出來。用玩具模型算兩個危害水準：

# %% tags=["hide-input"]
X_LO, X_HI = 0.2, 1.0
m_edges = np.arange(5.0, 8.01, 0.5)
e_edges = np.array([-np.inf, 0.0, 1.0, 2.0, 3.0, np.inf])
e_labels = ["ε < 0", "0–1", "1–2", "2–3", "ε > 3"]


def deagg_mag(x, lo, hi, n=300):
    """M 落在 [lo, hi) 的超越率貢獻。"""
    tot = 0.0
    for _, r, nu, mx in FAULTS:
        if lo >= mx:
            continue
        m = np.linspace(lo, min(hi, mx), n)
        tot += source_nu(nu) * np.trapezoid(
            p_exceed(x, m, r) * mag_pdf(m, mx), m)
    return tot


def deagg_eps(x, lo, hi, n=1500):
    """ε 落在 [lo, hi) 的超越率貢獻（首尾開放，各箱總和等於 λ）。"""
    tot = 0.0
    for _, r, nu, mx in FAULTS:
        m = np.linspace(M_MIN, mx, n)
        w = np.clip(norm.cdf(hi) - norm.cdf(np.maximum(eps0(x, m, r), lo)),
                    0.0, None)
        tot += source_nu(nu) * np.trapezoid(w * mag_pdf(m, mx), m)
    return tot


def deagg_means(x, n=2000):
    """反聚合的 (平均規模, 平均距離, 平均 ε)。"""
    num_m = num_r = num_e = den = 0.0
    for _, r, nu, mx in FAULTS:
        m = np.linspace(M_MIN, mx, n)
        p = p_exceed(x, m, r)
        w = source_nu(nu) * mag_pdf(m, mx) * p
        den += np.trapezoid(w, m)
        num_m += np.trapezoid(w * m, m)
        num_r += np.trapezoid(w * r, m)
        num_e += np.trapezoid(w * norm.pdf(eps0(x, m, r))
                              / np.maximum(p, 1e-300), m)   # E[ε|ε>ε₀]
    return num_m / den, num_r / den, num_e / den


fig = make_subplots(rows=1, cols=2, column_widths=[0.55, 0.45],
                    subplot_titles=("規模反聚合", "ε 反聚合"))
stats_txt, eps_share = [], {}
for i, x in enumerate([X_LO, X_HI]):
    cm = np.array([deagg_mag(x, a, b)
                   for a, b in zip(m_edges[:-1], m_edges[1:])])
    ce = np.array([deagg_eps(x, a, b)
                   for a, b in zip(e_edges[:-1], e_edges[1:])])
    eps_share[x] = 100 * ce[3:].sum() / ce.sum()
    fig.add_trace(go.Bar(x=[f"{a:.1f}–{b:.1f}" for a, b in
                            zip(m_edges[:-1], m_edges[1:])],
                         y=100 * cm / cm.sum(), name=f"PGA > {x} g",
                         legendgroup=str(x), marker_color=PALETTE[i]),
                  row=1, col=1)
    fig.add_trace(go.Bar(x=e_labels, y=100 * ce / ce.sum(),
                         name=f"PGA > {x} g", legendgroup=str(x),
                         showlegend=False, marker_color=PALETTE[i]),
                  row=1, col=2)
    mm, rr, ee = deagg_means(x)
    stats_txt.append(f"{x} g 的平均 M={mm:.2f}、R={rr:.1f} km、ε={ee:.2f}")
fig.update_xaxes(title_text="規模區間", row=1, col=1)
fig.update_xaxes(title_text="ε 區間（地動殘差的標準差數）", row=1, col=2)
fig.update_yaxes(title_text="貢獻比例（%）", row=1, col=1)
fig.update_yaxes(title_text="貢獻比例（%）", row=1, col=2)
apply_layout(fig, height=430, barmode="group", hovermode="x",
             title=f"圖 4　反聚合：{stats_txt[0]}｜{stats_txt[1]}")
fig

# %% [markdown]
# 兩張圖說了兩件不同的事。
#
# **左圖（規模）**：危害水準從 0.2 g 提高到 1.0 g，貢獻明顯往大規模端搬，平均規模與
# 平均距離都上升（圖標題的數字由程式帶入）——這就是那條著名的定性規律，**危害水準
# 愈高，主控的震源愈大、愈遠**。Baker 白皮書的算例是同一個現象：他的兩斷層例子在
# PGA = 0.376 g 時 $P(M=6.5)=0.77$、$P(M=7.5)=0.23$，到 1 g 時變成 0.58 / 0.42。
#
# 但要看清楚**搬的幅度**：在 1.0 g，貢獻最大的規模區間**仍然是 5.5–6.0**，只是
# $M\ge6.5$ 的份額從個位數升到三分之一以上。重心確實移動了，但沒有整個換手——中等規模
# 地震因為發生率高，即使需要極端的地動殘差，仍守住大半的貢獻。
#
# **右圖（$\varepsilon$）**：這才是本章那句反直覺結論的證據。在 0.2 g，貢獻集中在
# $\varepsilon$ 介於 0 到 2 之間——大致是「地動略高於平均」的普通地震。到了 1.0 g，
# $\varepsilon>2$ 的份額暴增：程式算出在 0.2 g 只佔約 11%，在 1.0 g 佔約 80%。
#
# 把兩張圖合起來，結論就完整了：**高強度端的危害，主要來自「中等規模、距離不遠、但
# 地動殘差高達兩三個標準差」的地震。** 稀有的不只是地震本身，更是地動——這是 21.4 節
# 那個 $\sigma_{\ln IM}\approx0.57$ 的直接後果。
#
# 兩個實務意涵。第一，**選設計地震波時不能只看 $M$ 與 $R$**，還必須挑 $\varepsilon$
# 匹配的紀錄，否則選出來的波在目標週期上的強度會系統性偏低（這正是條件平均譜 CMS 這類
# 方法要處理的問題）。第二，**高強度端的危害對 GMPE 的 $\sigma$ 極度敏感**：$\sigma$
# 估大一點尾巴就厚一點，2475 年回歸期的設計地動就抬高一截，而它恰好是 GMPE 裡最難估、
# 各家差異最大的參數之一。

# %% [markdown]
# ## 21.8 三種記憶都要裝
#
# {eq}`eq:psha` 裡的 $\nu_{\rm src}$ 是常數，所以傳統 PSHA 的危害曲線今年與明年一模
# 一樣。但這本書前面十幾章都在說：**地震率不是常數**。到這裡三種記憶結構已經全部到齊，
# 而且各自有明確的時間尺度與適用對象：
#
# | 記憶 | 模型 | 時間尺度 | 章 |
# |---|---|---|---|
# | 無記憶 | 齊次 Poisson | 長期平均 | 10 |
# | 正記憶（自我激發） | ETAS | 分鐘–月 | 13、14 |
# | 負記憶（自我抑制） | BPT 更新過程 | 數十–數百年 | 20 |
#
# **一個完整的危害模型，三種都要裝。** 而它們接進 {eq}`eq:psha` 的方式是同一個：
# **把常數 $\nu_{\rm src}$ 換成時間的函數 $\nu_{\rm src}(t)$**。這個替換的數學基礎，
# 就是把齊次 Poisson 過程換成非齊次的。
#
# ### 齊次到非齊次 Poisson 的替換
#
# 令 $\nu(t)$ 為隨時間變化的瞬時率。把時窗 $[0,T]$ 切成 $n$ 個等長小段 $\Delta=T/n$，
# 每段裡假設率近似為常數 $\nu(t_i)$ 且各段獨立。「整段時窗都沒有事件」的機率是
#
# $$\begin{aligned}
# P\bigl(N(T)=0\bigr)
#   &= \prod_{i=1}^{n}\bigl[1-\nu(t_i)\,\Delta + o(\Delta)\bigr]
#    = \exp\left\{\sum_{i=1}^{n}
#      \ln\bigl[1-\nu(t_i)\Delta+o(\Delta)\bigr]\right\} \\
#   &= \exp\left\{-\sum_{i=1}^{n}\nu(t_i)\,\Delta + o(1)\right\}
#   \;\xrightarrow[n\to\infty]{}\;
#   \exp\left[-\int_{0}^{T}\nu(t)\,\mathrm{d}t\right].
# \end{aligned}$$
#
# 第二行用了 $\ln(1-u)=-u+O(u^2)$，而 $\sum O(\Delta^2)=O(1/n)\to0$。於是
#
# $$P\bigl(\text{exceedance in }[0,T]\bigr)
#   \;=\; 1-\exp\left[-\int_{0}^{T}\lambda_{IM}(x;t)\,\mathrm{d}t\right],$$ (eq:nhpp)
#
# 其中時變的危害積分是 $\lambda_{IM}(x;t)=\sum_{\rm src}\nu_{\rm src}(t)\int\!\!\int
# P(IM>x\mid m,r)f_M(m)f_R(r)\,\mathrm{d}m\,\mathrm{d}r$。{eq}`eq:nhpp` 就是
# 21.6 節那條時窗機率的推廣：$\lambda t$ 換成 $\int\lambda\,\mathrm{d}t$。整個
# 時變危害的框架，數學上就只有這一步。
#
# **要看清楚哪一層變了、哪一層沒變。** 變的只有最外層的 $\nu_{\rm src}(t)$，規模分布、
# 距離分布、GMPE 三層原封不動——這等於假設「餘震序列期間的地震，其規模–距離–地動關係
# 與平時一樣」。這不見得對（餘震的深度、機制、破裂尺度都可能與背景地震有系統性差異），
# 但它讓整套計算可以複用，代價是把時變性全部壓進一個純量。
#
# 兩種記憶各自怎麼提供 $\nu_{\rm src}(t)$。**短期用 ETAS**：$\nu(t)=\lambda^*(t\mid
# H_t)$，第 13、14 章的條件強度——微妙之處在於 $\lambda^*$ 條件於歷史，而未來還沒發生
# 的餘震會再觸發自己的餘震，嚴格算 {eq}`eq:nhpp` 需要對未來所有可能歷史取期望，實務上
# 就是跑數萬條模擬目錄再平均，這也是時變危害系統遠比傳統 PSHA 昂貴的原因。**長期用
# BPT**：$\nu(t)=h(t-t_{\rm last})$，$h$ 是第 20 章的危害函數 {eq}`eq:hazard-def`。
#
# 先看短期那一端：大地震之後第一週，整條超越機率曲線怎麼動。

# %% tags=["hide-input"]
WEEKS = 365.25 / 7.0                          # 一年幾週
BOOST = 200.0                                 # 餘震期率放大倍數（示意）


def weekly_prob(x, boost=1.0):
    """一週內超越機率：把 ν 乘上 boost，套式 (eq:nhpp)。"""
    return 1.0 - np.exp(-hazard(x) * boost / WEEKS)


x_g3 = np.logspace(-2, 0.0, 70)
p_bg = np.array([weekly_prob(x) for x in x_g3])
p_af = np.array([weekly_prob(x, BOOST) for x in x_g3])
i_ref = int(np.argmin(np.abs(x_g3 - 0.2)))

fig = go.Figure()
for y, name, color in [(p_bg, "平時（長期背景率）", ACCENT),
                       (p_af, f"大地震後第一週（率 ×{BOOST:.0f}，示意）",
                        QUAKE_COLOR)]:
    fig.add_trace(go.Scatter(x=x_g3, y=y, mode="lines", name=name,
                             line=dict(color=color, width=2.6)))
fig.add_vline(x=x_g3[i_ref], line=dict(color="#888", width=1, dash="dot"))
fig.update_xaxes(title_text="PGA（g）", type="log")
fig.update_yaxes(title_text="一週內超越機率", type="log")
apply_layout(fig, height=430, hovermode="x", title=(
    f"圖 5　時變危害：在 {x_g3[i_ref]:.2f} g，一週內超越機率由 "
    f"{p_bg[i_ref]:.1e} 升到 {p_af[i_ref]:.1e}"))
fig

# %% [markdown]
# 整條曲線往上平移。以 0.2 g 為例（圖上的灰色虛線），平時一週內超越的機率是萬分之六
# 左右，序列期間變成一成一——**升了兩個數量級**。放大倍數只是示意（真實的 ETAS 率增益
# 隨主震規模、距離與時間衰減，要由第 13、14 章的模型算），但量級是真實的。而大地震之後
# 最需要做決定的正是那幾天——搜救人員能不能進入受損建築、災民能不能回家、哪些結構要先
# 貼紅單。**這些決定面對的危害不是那條長期平均曲線，而是上面那條紅線。** 這正是第 22
# 章作業化系統存在的理由，也是 21.10 節要接上的地方。
#
# ### 長期那一端：BPT 讓危害重新分配
#
# 短期的效應是「整條曲線上移」，方向單一而直觀。長期的效應完全不同：它**有正有負**。
# 把第 20 章的 {eq}`eq:renewal-prob` 套在 BPT 上，看一條斷層的 50 年條件破裂機率隨
# 「距上次破裂的時間」怎麼變，並與 Poisson 的定值對照。

# %% tags=["hide-input"]
T_R_FAULT, CV_FAULT, WINDOW = 300.0, 0.50, 50.0     # 年、aperiodicity、時窗
bpt = invgauss(mu=CV_FAULT**2, scale=T_R_FAULT / CV_FAULT**2)

t_since = np.linspace(0.0, 900.0, 900)
p_bpt = (bpt.sf(t_since) - bpt.sf(t_since + WINDOW)) / np.maximum(
    bpt.sf(t_since), 1e-300)
p_poi = 1.0 - np.exp(-WINDOW / T_R_FAULT)
i_x = int(np.argmax(p_bpt > p_poi))

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_since, y=p_bpt, mode="lines",
                         name=f"BPT（T_r={T_R_FAULT:.0f} 年, c_v={CV_FAULT}）",
                         line=dict(color=PALETTE[0], width=2.6)))
fig.add_trace(go.Scatter(x=t_since, y=np.full_like(t_since, p_poi),
                         mode="lines", name=f"Poisson（定值 {p_poi:.3f}）",
                         line=dict(color="#8a8a8a", width=1.6, dash="dot")))
fig.add_vline(x=t_since[i_x], line=dict(color=QUAKE_COLOR, width=1.2,
                                        dash="dash"))
fig.add_annotation(x=t_since[i_x], y=p_poi, text="交叉", showarrow=True,
                   arrowhead=2, ax=48, ay=-34, font=dict(size=11))
fig.add_vrect(x0=0, x1=t_since[i_x], fillcolor=PALETTE[2], opacity=0.07,
              line_width=0)
fig.update_xaxes(title_text="距上次破裂的時間 T（年）")
fig.update_yaxes(title_text=f"未來 {WINDOW:.0f} 年破裂機率")
apply_layout(fig, height=430, hovermode="x", title=(
    f"圖 6　BPT vs Poisson 的 {WINDOW:.0f} 年條件機率：T < {t_since[i_x]:.0f} 年"
    f"（{t_since[i_x]/T_R_FAULT:.2f} T_r）時 BPT 較低，長期趨於 {p_bpt[-1]:.3f}"))
fig

# %% [markdown]
# Poisson 是一條水平線——無記憶，距上次破裂多久完全不影響。BPT 則從零開始爬：
# **剛破裂完的斷層，時變模型給的機率比時間獨立模型低得多**（圖上淺綠色那一段）。交叉點
# 大約落在 0.44 個平均複發期；之後 BPT 超過 Poisson，最後停在一個高原，高度約為
# Poisson 的兩倍——這正是第 20 章推導的漸近危害率 $1/(2T_rc_v^2)$ 換算成 50 年機率的
# 結果（條件機率在趨於高原之前會先略微過衝，圖上看得到）。
#
# 這張圖說明了一件政策上很重要的事：**把時變模型放進國家級危害度圖，不是「一律加碼」
# 的安全係數，而是一次重新分配。** 剛破裂過的斷層附近，設計地動會**下修**；沉寂已久的
# 斷層附近會**上修**。使用者必須同時接受兩端——只挑上修那一半來用是不誠實的，這也是
# 各國官方模型導入 BPT 時都非常謹慎、通常只用在資料最好的少數幾條斷層上的原因。
#
# 三種記憶合起來的圖像是：**Poisson 給骨架，BPT 在長期尺度上把危害在斷層之間重新分配，
# ETAS 在事件之後的短期尺度上把整條曲線暫時抬起來。** 三者的時間尺度差了四到六個數量
# 級，所以實務上是分層疊加而不是寫成同一條式子——而那個難題第 20 章 20.10 節已經說過
# 它還沒被解決。

# %% [markdown]
# ## 21.9 台灣：TEM 的兩代模型
#
# 台灣有自己的國家級 PSHA。**台灣地震模型**（Taiwan Earthquake Model, TEM）由學界團隊
# 建立，至今有兩代正式發布。
#
# **TEM PSHA2015**（Wang et al. 2016, *TAO*）的結構是標準的「斷層源 + 面震源」兩層：
# 斷層源採用 TEM 地質團隊辨識的 **38 條孕震構造**（Shyu et al. 2016 的構造資料庫），
# 面震源則有 **28 個淺層面震源**、**4 個隱沒帶板塊間震源**與 **12 個隱沒帶板塊內震源**，
# 共三類背景地震活動，各自配上對應構造環境的 GMPE。場址條件上，它只評估到**工程岩盤**
# （$V_{S30}=760$ m/s），因此**忽略場址放大效應**——這是讀 2015 版危害圖時最重要的一條
# 但書，真實地表在軟弱沖積層上的地動可以比岩盤值高出可觀的倍數。
#
# 定性的結論圖像：**危害最高的區域在西南部與東部縱谷**；西部人口稠密的都會之中，以
# <strong>台南（短週期）與台中（長週期）</strong>最受關注——原文據此指出台南的低樓層建物與台中的
# 高樓層建物的耐震設計特別重要（低樓層建物自然週期短、高樓層長）。本節不引用任何具體
# 危害數值：那些數字有明確的適用條件（哪一版模型、哪個場址條件、哪個回歸期、哪個
# 週期），抽出來單獨引用幾乎必然被誤讀。
#
# **TEM PSHA2020**（Chan et al. 2020, *Earthquake Spectra*）的四項更新，每一項都對應
# 本書前面的一個主題：
#
# 1. **三維幾何**。構造資料庫擴充並納入新辨識的三維構造，不再把斷層當成簡單平面。這
#    直接影響 {eq}`eq:psha` 裡的 $f_R$——距離分布是由幾何算出來的。
# 2. **多構造同時破裂**。允許相鄰構造連動，因而能產生比單一構造上限更大的地震。1999
#    集集與 2024 花蓮都提醒過這件事的重要性。
# 3. **場址放大**。新一組 GMPE 並加入場址放大係數，補上 2015 版留下的空白。
# 4. **對斷層源導入 BPT**。這是本章與第 20 章的交會點：每一條孕震構造掛上一個
#    {eq}`eq:bpt-pdf` 的複發分布，用 {eq}`eq:renewal-prob` 算條件機率，再把結果當成
#    {eq}`eq:psha` 裡那條斷層的 $\nu_{\rm src}$。**這是台灣官方危害度圖裡唯一時間相依
#    的成分**，也就是圖 6 那條藍線在真實模型裡的角色。
#
# 另外，無法歸屬到特定構造的地殼地震，同時採用**面震源**與**平滑核**兩種模型——這是
# 「不知道該用哪一種就兩種都用、給權重」的 logic tree 思維（21.11 節）。下一代模型
# （PSHA2025）正在發展中；本書不引用它的任何數字或結論，因為它尚未正式出版。

# %% [markdown]
# ## 21.10 把 OEF 接上 PSHA
#
# 21.8 節推的 {eq}`eq:nhpp` 是數學；這一節講它在真實世界裡被接上了什麼。
#
# ### PBEE：危害只是第一格
#
# 風險評估的主流框架是 **PBEE**（performance-based earthquake engineering, Cornell &
# Krawinkler 2000），把三個機率描述串起來：**場址地震危害**（本章的 $\lambda_{IM}$）、
# **易損性**（給定地動，結構損壞到什麼程度的機率）、**後果與曝險**（給定損壞，損失多少
# 的機率）。三者相乘積分，得到「損失超過某門檻的年發生率」。傳統 PBEE 用的是古典
# PSHA：齊次 Poisson、只算主震、並假設兩次主震之間有足夠時間修復，所以風險度量是
# **時不變**的——這正是 21.8 節動手術的對象。把 OEF 接進 PBEE 的經典嘗試是 **Yeo &
# Cornell (2009)** 的餘震危害分析：地震發生改用**非齊次 Poisson 過程**，均值函數由
# （Reasenberg–Jones 改良的）Omori 律給出，於是可以估出「單位時間內超過某地動強度門檻
# 的**時變率**」，這就是 {eq}`eq:nhpp` 的具體實作。但接上去之後冒出兩個不在原本框架裡
# 的新問題：**序列期間來不及修復，損失會在多次事件中累積**（需要「已受損結構」的易損性
# 模型）；以及**曝險本身也可能時變**（疏散、遷移）。兩者至今都沒有成熟的解法。
#
# ### 兩個真實的系統
#
# **義大利：OELF 與 Mantis-K。** 把 OEF 的地震率直接接上全國住宅建物清冊與易損性模型，
# 每日（或每個 $M\ge3.5$ 事件之後）產出**預期不可使用建物數、倒塌建物數、傷亡人數**，
# 時間範圍與 OEF 預報一致。這套引擎叫 **Mantis-K**（Iervolino et al. 2015），是
# **OELF**（operational earthquake loss forecasting）的原型。
#
# **紐西蘭：從時變危害到強制補強。** Canterbury 序列（2010 M7.1 Darfield → 2011 M6.2
# Christchurch）大幅改變了未來數十年的預期危害，於是為 Christchurch 重建開發了一套
# **時變地震危害模型**（Gerstenberger et al. 2014, 2016）：年步長、最長 **50 年**、
# 時變與時不變元件取極大值，第二階段用**結構化專家判斷**決定各模型權重。2016 年
# Kaikōura 地震之後，這套模型被用來計算「地動超越建築規範要求的**機率增益**」，直接
# 支援了中紐西蘭**無筋磚造建築強制補強**的政策決定。
#
# 請留意這條因果鏈的完整性：ETAS 類模型算出短期地震率 → {eq}`eq:nhpp` 把它變成地動
# 超越機率 → 與規範設計值比較得到機率增益 → 政府據以要求業主補強。**一條統計模型的
# 輸出，最後變成了法律義務。** 這是本書所有內容裡，機率走得離社會最近的一次。
#
# 這一切的定調文件是 ICEF 報告（Jordan et al. 2011），它明確要求：**OEF 必須與 PSHA
# 的長期預報一致地提供完整的危害描述——地動超越機率，而不只是地震發生機率**。這句話
# 的技術內容就是 {eq}`eq:nhpp`：短期預報若只給「M5 以上發生機率」，使用者沒辦法拿它做
# 任何工程決策；只有換算成地動超越機率，它才能與規範、與既有的危害曲線、與 PBEE 的
# 其他兩格對話。

# %% [markdown]
# ## 21.11 兩種不確定性與 logic tree
#
# 21.2 節結尾埋了一個伏筆：{eq}`eq:psha` 裡混著兩種性質完全不同的不確定性。
#
# **偶然不確定性**（aleatory variability）：即使模型完全正確，結果仍然隨機。「下一次
# 地震規模多大」（$f_M$）、「震央在哪」（$f_R$）、「地動殘差多大」（$\varepsilon$）
# 都屬於這一類，而它們**已經被積分進 {eq}`eq:psha` 裡了**——危害積分做的事，本質上
# 就是把偶然不確定性積掉。
#
# **認知不確定性**（epistemic uncertainty）：我們不知道真相是什麼。$\nu_{\rm src}$ 該
# 取多少？$m_{\max}$ 是 7.0 還是 7.3？該用哪一條 GMPE？$\sigma_{\ln IM}$ 該多大？斷層
# 該怎麼分段？$c_v$ 取 0.3 還是 0.7（第 20 章）？這些**不能**被積分掉，因為它們不是
# 隨機的——真相只有一個，只是我們不知道。
#
# 這個區分有可操作的後果：**偶然不確定性變大會讓危害曲線的尾巴變厚（一條曲線改變
# 形狀）；認知不確定性變大會讓你得到很多條不同的曲線。** 表達後者的標準工具是 **logic tree**：對每一個無法由資料決定的選擇，列出所有
# 合理選項、各給一個權重，然後把整棵樹跑完，輸出一**族**曲線（通常報告成平均危害曲線
# 加上分位數帶）。這棵樹可以長得非常大——**UCERF3 有 5760 個 logic-tree 分支**（第 20
# 章）。這個數字本身就是認知不確定性規模的量尺：一個需要 5760 個分支才能表達的模型，
# 等於在說「這件事我們真的不太確定」。
#
# 兩個常見陷阱。**權重不是機率**：它是專家對「這個選項是對的」的信念度，不是任何頻率
# 意義下的機率。**平均危害曲線不對應任何一個分支**：每一條分支曲線都代表一個「內部
# 一致」的世界觀，平均之後那條曲線可能不是任何一個專家會簽名的模型。

# %% [markdown]
# ## 21.12 常見誤解與陷阱
#
# **誤解 1：「475 年回歸期表示我這輩子遇不到。」** 21.6 節的 admonition 已經拆過三次。
# 補一個算式：一棟樓若使用 100 年，期間至少超越一次設計地動的機率是
# $1-e^{-100/475}\approx19\%$——**大約每五棟就有一棟，在它的壽命內會遇上一次
# 「475 年一遇」的地動**。這個數字通常比大眾的直覺高一個檔次。
#
# **誤解 2：「高強度端的危害來自最大的地震。」** 21.7 節的 $\varepsilon$ 反聚合直接
# 反駁：在 1.0 g，約八成的貢獻來自 $\varepsilon>2$ 的情形，而主控規模區間仍然是
# 5.5–6.0。**稀有的是地動殘差，不只是地震。** 這也意味著 $m_{\max}$ 估得準不準，
# 對高強度端的影響常常不如 $\sigma_{\ln IM}$ 大。
#
# **誤解 3：「除叢只是資料清理，不影響結果。」** 第 12 章 12.8 節列出除叢對 PSHA 的
# **三重後果**：**(a)** 主震的定義本身不可驗證——不同演算法留下的主震數差 6.1 倍；
# **(b)** 丟掉餘震會低估危害，而餘震一樣致災；**(c)** 被壓低的 $b$ 值會在 $m>m_x$ 之上
# **高估**危害，這是破產點 {eq}`eq:mx` 的直接後果。「一個低估、一個高估剛好抵銷」的
# 辯護不成立：兩個誤差只在 $m=m_x$ 這**單一**規模上抵銷，其他每個規模都是錯的，方向還
# 隨規模翻轉。除叢仍有它的理由（{eq}`eq:psha` 的 Poisson 假設要求事件獨立），但要記住
# **除叢是為了滿足模型的假設，不是為了描述地球**。
#
# **誤解 4：「危害度就是風險。」** 不是。**危害**（hazard）是「地面會怎麼晃」，
# **風險**（risk）是「晃了之後會損失什麼」。一片沒有人的荒野可以有極高的危害與接近零
# 的風險；一座蓋在中等危害區的老舊高密度市區可以有極高的風險。PSHA 只算 PBEE 三格裡的
# **第一格**，把它當成風險評估，等於把易損性與曝險兩層當成常數。
#
# **誤解 5：「時變模型一定比時間獨立模型保守。」** 圖 6 已經反駁：剛破裂完的斷層，
# BPT 給的機率**低於** Poisson。時變模型是重新分配，不是加碼。另外，PSHA 的機率也很難
# 前瞻檢驗——2475 年回歸期的曲線在人類時間尺度上收集不到足夠的超越事件，而「拿許多個
# 場址當樣本」的替代方案又要求場址獨立，同一次地震卻會同時影響一大片場址。
#
# **誤解 6：「危害曲線是觀測到的。」** 它是**模型的輸出**，不是量測值，底下疊著目錄
# 完整度（第 11 章 $M_c$）、除叢（第 12 章）、GR 律外插到沒觀測過的規模、GMPE 從有限
# 紀錄迴歸再外插、斷層幾何與滑移率的地質判讀、複發模型的選擇——每一層都有不確定性，
# 而它們最後全被壓進一個看起來很權威的數字裡。

# %% [markdown]
# ## 21.13 研究前沿與未解問題
#
# **用 ETAS 模擬合成目錄，繞開除叢與 Poisson 假設。** 這是目前最有希望的一條路，而且
# 同時解決誤解 3 的兩端：用全目錄擬合 ETAS（不必區分主震與餘震），模擬數十萬年的合成
# 目錄，再把每一條丟進 {eq}`eq:psha` 的前四步統計超越次數。這樣一來，**除叢不必做**
# （ETAS 只依賴全目錄的 GR 律）、**Poisson 假設不必要**（叢集自然出現在模擬裡）、
# **時空相關性自動涵蓋**；代價是計算量與 ETAS 參數本身的不確定性被搬到台前。Mizrahi
# et al. (2021) 對除叢的批評與這條路線正好接得上：與其在一個任意的分類上做統計，不如用
# 一個能同時描述背景與觸發的點過程模型。
#
# **UCERF3-ETAS：把短期與長期裝進同一份模擬。** 加州的做法是把時變的長期複發率（第 20
# 章的更新過程）餵進 ETAS 當背景率，讓短期觸發與長期記憶共存在同一份模擬目錄裡，概念上
# 是本章那張「三種記憶」表格的完整實現。但第 20 章 20.10 節的難題仍在：兩層之間沒有真正
# 的耦合——大地震既觸發餘震又消耗應力，還是沒有被同一條式子描述。
#
# **時變危害的作業化。** 目前全世界只有少數幾個系統真的每天在跑時變危害，多半停留在
# 研究或有限授權的階段。障礙不在數學（{eq}`eq:nhpp` 很簡單），而在於模擬的計算量、
# 即時目錄的品質（大震後 $M_c$ 會急遽劣化，第 11 章），以及最麻煩的一項——**每天都在
# 變的設計地動，工程與法規體系該怎麼使用**。規範的邏輯建立在「一個固定的數字」上。
# 這不是統計問題，是制度問題，第 22 章會再回到它。
#
# **$\varepsilon$ 該不該截斷。** 21.7 節顯示高強度端的危害幾乎全由大 $\varepsilon$ 撐起。
# 物理上 $\varepsilon$ 應該有上限，但截在 3 還是 4 個標準差，對核設施那種 $10^{-5}$ 到
# $10^{-6}$ 年率的評估影響巨大，而經驗資料在那個尾巴上幾乎是空的。
#
# **台灣的具體缺口有兩個。** 其一，孕震構造的複發參數（$T_r$、$c_v$）依賴古地震探槽，
# 而台灣的探槽數量與定年精度離加州、日本都還有距離——第 20 章 20.7 節那個「四個帶
# 誤差的數字」的困境在台灣更嚴重。其二，台灣的危害度評估慣用視窗法除叢，而 Mizrahi
# et al. (2021) 的批評完全適用；用 ETAS 全目錄模擬重做一次台灣的背景率，是一件既可行、
# 又能直接檢驗的事。

# %% [markdown]
# ## 21.14 附錄：本章推導細節
#
# ### A. 截斷 GR 的累積分布、抽樣與平均規模
#
# 由 {eq}`eq:trunc-gr` 積分得累積分布
# $F_M(m)=\bigl[1-e^{-\beta(m-m_{\min})}\bigr]\big/\bigl[1-e^{-\beta(m_{\max}
# -m_{\min})}\bigr]$。這個式子在模擬時直接可用：由第 10 章的反函數法，取
# $U\sim\mathrm{Unif}(0,1)$ 則 $m=m_{\min}-\beta^{-1}\ln[1-U(1-e^{-\beta(m_{\max}
# -m_{\min})})]$ 就是一個截斷 GR 的抽樣。平均規模由分部積分得 $E[M]=m_{\min}+1/\beta
# -(m_{\max}-m_{\min})e^{-\beta(m_{\max}-m_{\min})}/[1-e^{-\beta(m_{\max}-m_{\min})}]$，
# $m_{\max}\to\infty$ 時退化成指數分布的 $m_{\min}+1/\beta$（$b=1$ 時即
# $m_{\min}+0.434$）。
#
# ### B. 危害積分的 $\varepsilon$ 形式
#
# 把 21.7 節的 $P(IM>x\mid m,r)$ 展開成 $\varepsilon$ 的積分，{eq}`eq:psha` 有一個等價
# 寫法 $\lambda_{IM}(x)=\sum_{\rm src}\nu_{\rm src}\iiint_{\varepsilon_0(x,m,r)}^{\infty}
# \varphi(e)f_M(m)f_R(r)\,\mathrm{d}e\,\mathrm{d}m\,\mathrm{d}r$。被積函數不再含 $x$
# ——$x$ 只出現在積分下限裡。這個寫法讓反聚合的三個維度完全對稱，並清楚顯示 **PSHA
# 本質上是在 $(m,r,\varepsilon)$ 這個三維情境空間上，對一個由 $x$ 決定的區域做積分**；
# 改變 $x$ 就是移動那個區域的邊界。
#
# ### C. 為什麼超越事件是 Poisson 過程
#
# 21.2 節用了兩條性質，這裡補上理由。
#
# **稀疏化。** 設 $\{t_i\}$ 是速率 $\nu$ 的 Poisson 過程，每個點獨立地以機率 $p$ 保留。
# 在 $[0,T]$ 上 $N\sim\mathrm{Poisson}(\nu T)$、$K\mid N\sim\mathrm{Bin}(N,p)$，於是
#
# $$P(K=k) = \sum_{n\ge k}\frac{(\nu T)^{n}e^{-\nu T}}{n!}
#   \binom{n}{k}p^{k}(1-p)^{n-k}
#   = \frac{(\nu Tp)^{k}e^{-\nu T}}{k!}\sum_{j\ge 0}\frac{[\nu T(1-p)]^{j}}{j!}
#   = \frac{(\nu pT)^{k}}{k!}e^{-\nu pT},$$
#
# 其中令 $j=n-k$。所以 $K\sim\mathrm{Poisson}(\nu pT)$。本章的 $p$ 是 $P(IM>x)$，由
# $(M,R,\varepsilon)$ 這些**獨立的標記**決定、與事件時間無關，稀疏化的獨立性因此成立。
#
# **疊加。** 獨立 Poisson 過程之和仍是 Poisson、速率相加（Poisson 分布對可加性封閉的
# 直接結果），這就是 {eq}`eq:psha` 那個 $\sum_{\rm src}$ 的正當性。
#
# **哪裡會壞掉。** 兩條性質都要求**震源之間、事件之間獨立**：餘震破壞後者（這是除叢的
# 動機），應力轉移破壞前者。所以 Poisson 骨架不是自然律，是為了讓計算可行而付出的
# 代價——21.13 節那條「用 ETAS 模擬合成目錄」的路線，正是要把這個代價省掉。
#
# ### D. 玩具模型的設定
#
# 本章所有圖使用同一組參數，列在這裡以便重現。
#
# | 量 | 近而小的斷層 | 遠而大的斷層 |
# |---|---|---|
# | 距離 $r$（km） | 10 | 25 |
# | $M\ge5$ 年率（1/yr） | 0.10 | 0.05 |
# | $m_{\max}$ | 6.0 | 8.0 |
#
# | 共用參數 | 值 |
# |---|---|
# | $b$（$\beta=b\ln10$） | 1.0 |
# | $m_{\min}$（基準） | 5.0 |
# | $\sigma_{\ln IM}$ | 0.57 |
# | GMPE | $\ln PGA=-0.152+0.859m-1.803\ln(r+25)$ |
#
# 這條 GMPE 是 Baker 白皮書示範用的 Cornell 型式。可以驗證：代入 $m=6.5$ 與
# $r=3,10,30$ km 得到 $\overline{\ln PGA}=-0.5765,\,-0.9788,\,-1.7937$，與白皮書的算例
# 逐位相符。積分一律用 `np.trapezoid` 在規模上做等距數值積分；距離分布退化成點震源，
# 所以 $f_R$ 的積分消失。本章沒有任何亂數，所有圖都是確定性計算，重跑必然一致。

# %% [markdown]
# 回頭看，這一章其實只做了一件事：把「地震會發生」翻譯成「地面會晃」。翻譯的工具是
# 一個積分 {eq}`eq:psha`，而積分的每一層都是前面某一章的產品。這是全書結構上最收斂的
# 一刻——第 11 章的 GR 律、第 12 章的除叢代價、第 13 章的觸發、第 20 章的斷層記憶，
# 在這裡全部匯進同一條式子。
#
# 但也正因為如此，這一章是「一個數字算得出來，不等於它站得住腳」這句話的終極考場。
# 危害曲線是這本書裡最有影響力的輸出——它決定了鋼筋要多粗、樓要蓋多高、保費要收多少
# ——同時也是假設疊得最厚的一個。這不是要你不信任規範：規範是目前為止把已知的東西整理
# 得最好的結果，每一個假設都寫在文件裡供人檢查。真正該記住的是**規範是模型的輸出，而
# 模型會更新**——TEM 從 2015 到 2020 換掉了場址放大、加進多構造破裂與斷層記憶，就是
# 最好的證據。懂得掀開引擎蓋的人，才有資格參與下一次更新。
#
# 最後還有一件事沒做完。本章的時變危害停在數學上：{eq}`eq:nhpp` 告訴我們怎麼把時變的
# 地震率換成地動超越機率。但一個系統要 24 小時運轉、在大地震後很短的時間內產出結果、
# 還要把那個結果講給不懂機率的人聽，需要的遠不只是一條式子。{doc}`第 22 章
# <22_operational_systems>`就來處理這一段：從 STEP 到現行的作業化系統，模型上線之後
# 會冒出哪些不再是統計的問題。
