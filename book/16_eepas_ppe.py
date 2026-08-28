# %% [markdown]
# # 16. EEPAS 與 PPE：把 Ψ 變成可運轉的預報
#
# {doc}`第 15 章 <15_psi_phenomenon>`把 Ψ（前兆尺度增加）攤開來看：
# 大地震之前，震源區的中小地震會在**規模水準與發生率上同時抬升**；
# 從累積規模異常 $C(t)$ 量出前兆規模 $M_P$、前兆時間 $T_P$、前兆面積
# $A_P$ 之後，三條迴歸把它們全部繫在主震規模上。那一章也把壞消息
# 講完了：辨識要先知道主震在哪、多大，本質是回溯的；同一個主震有
# 好幾組同樣合格卻差到十倍的辨識；47 個主震裡有 13 個一組都找不到。
#
# 於是問題變成：**一個無法事前辨識的前兆現象，還能拿來做預報嗎？**
#
# **EEPAS**（Every Earthquake a Precursor According to Scale，
# Rhoades & Evison 2004）的答案是一步漂亮的側身——**那就不要辨識**。
# 假設每一個地震都是某個未來更大地震的前兆，只是它「預告」的尺度
# 由它自己的規模決定。這一步把 Ψ 從一種回溯的敘事變成一台可以自動
# 運轉、可以被概似函數評分、可以送進 CSEP 擂台的機器。
#
# 這一章從那一步側身開始，把整台機器拆開：三個機率核為什麼是常態、
# 對數常態與二維常態（16.2）；把它們加起來為什麼還需要一個正規化
# 函數 $\eta(m)$ 與一個補償函數 $\Delta(m)$（16.3、16.4）；餘震權重
# $w_i$ 其實就是第 14 章的背景機率（16.5）；基準模型 PPE 為什麼是
# 整個領域繞不開的參照（16.6）。然後是家族、補償、二十年成績單，
# 以及一份誠實到有點刺眼的失敗清單。
#
# 前面幾章的式子本章一律引用不重推：第 10 章的點過程對數概似
# {eq}`eq:pp-loglik`、第 12 章的正規化 Omori 密度
# {eq}`eq:omori-density`、第 14 章的背景機率 {eq}`eq:rho-phi`，
# 以及第 15 章的三條 Ψ 迴歸與時空取捨。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 16.1 放棄辨識這一步
#
# ### 從「哪些地震是前兆」到「每個地震都是前兆」
#
# 第 15 章的三條迴歸寫成以前兆規模為自變數的形式：
#
# $$M_m = a_M + b_M M_P,\qquad
#   \log_{10} T_P = a_T + b_T M_P,\qquad
#   \log_{10} A_P = a_A + b_A M_P$$
#
# 要用它們做預報，得先有一個 $M_P$；要有 $M_P$，得先辨識出一段前兆
# 序列；要辨識前兆序列，得先知道主震——這個循環在事前是打不開的。
#
# EEPAS 的做法是把 $M_P$ 這個「群體的統計量」換成 $m_i$ 這個「單一
# 地震的規模」。目錄裡每一個 $m_i\ge m_0$ 的地震，都被當成一次
# $M_P=m_i$ 的觀測；它對未來的貢獻，就是把三條迴歸的**條件分布**
# 攤在時間、規模、空間上。不需要框選、不需要判準、不需要知道主震
# 在哪。**辨識這一步被整個刪掉了。**
#
# 這是全章最該記住的一句話，值得再說一次：EEPAS 從來不需要做出
# Ψ 辨識，所以第 15 章列出的三個困難——辨識不唯一、28% 找不到、
# 判準含有事後知識——**它一個都不用面對**。
#
# ### 代價：一個誠實的「弱預報模型」
#
# 天下沒有白吃的午餐。原作者自己在 2004 年的論文裡用了
# 「weak forecasting model」這個詞，並且給了量化的比較：EEPAS 對
# **單一地震**能達到的最大概似比只有 29.7；而如果允許人工辨識出
# 明顯的前兆群震，同樣的比較可以到 114。**放棄辨識，把預報強度
# 削掉了大約四分之三。**
#
# 削掉的東西是資訊。人工辨識用到了「這一群地震看起來像一個前兆
# 序列」這個判斷，而 EEPAS 只用到「這個地震的規模是 $m_i$」。前者
# 的資訊量顯然大得多——問題在於它不能事前取得。
#
# 換來的是三件事，而且每一件都是可運轉的預報模型的必要條件：
# **可自動化**（輸入是目錄，輸出是率密度網格，中間沒有人的判斷）、
# **可檢驗**（有明確的概似函數，可以做 AIC、資訊增益、CSEP 一致性
# 檢定）、**可否證**（下面 16.7 節那個 ERDEEP 的故事就是一個假設
# 被自己的檢驗推翻的乾淨案例）。
#
# ### 一個必然的系統性偏移
#
# 這一步替換有一個可預期的副作用，值得先講清楚，否則讀參數表時
# 會困惑。Ψ 的 $M_P$ 定義是**前兆期最大三個地震的規模平均**，而
# EEPAS 把**每一個**地震的規模都當成一次 $M_P$。前兆群本身近似
# 一個 GR 集合，小地震的數量遠多於大地震，所以 EEPAS 看到的
# 「$M_P$ 樣本」平均而言比 Ψ 辨識用的小了不少。
#
# 迴歸式因此必須自我補償：同樣要預告一個 $M_m$ 的主震，$m_i$ 小了
# 就得把 $a_M$ 調大；同樣要指向同一個時刻，$m_i$ 小了就得把 $a_T$
# 調大。2004 年紐西蘭的擬合正是如此——$a_M$ 比 Ψ 迴歸大 0.5、
# $a_T$ 大 0.14。**這不是矛盾，是可預期的統計後果。** 也因此，
# 第 15 章那組 $a_M=3.16$、$a_T=1.36$ 的數字**不能直接搬進 EEPAS
# 當參數用**；Ψ 迴歸提供的是函數形式，數值一律逐區以最大概似重估。

# %% [markdown]
# ## 16.2 三核為何是常態、對數常態、二維常態
#
# 學生最常問的問題是：時間核為什麼不是 Omori 那種冪次律？答案
# 樸素得幾乎令人失望——**因為它們是三條線性迴歸的殘差分布**。
#
# ### 起點：一條線性迴歸的條件分布
#
# 取最一般的簡單線性迴歸
#
# $$Y = a + bX + \epsilon,\qquad \epsilon \sim N(0,\sigma^2),$$
#
# 其中 $\epsilon$ 與 $X$ 獨立。給定 $X=x$，$Y$ 只剩下 $\epsilon$ 這
# 一個隨機來源，於是
#
# $$Y \mid X=x \;\sim\; N\!\left(a+bx,\ \sigma^2\right),
#   \qquad
#   p(y\mid x)=\frac{1}{\sigma\sqrt{2\pi}}
#     \exp\!\left[-\frac12\left(\frac{y-a-bx}{\sigma}\right)^2\right].$$
#
# 這條式子看起來平凡，但它是本章所有核函數的唯一來源。接下來的
# 工作只是把三條 Ψ 迴歸逐一代進去，注意 $Y$ 是誰、$X$ 是誰，以及
# **$Y$ 有沒有取過對數**。
#
# ### 規模核：直接就是常態
#
# 取 $Y=M_m$、$X=M_P\to m_i$、殘差標準差記為 $\sigma_M$：
#
# $$g(m\mid m_i)=\frac{1}{\sigma_M\sqrt{2\pi}}
#   \exp\!\left[-\frac12\left(\frac{m-a_M-b_M m_i}{\sigma_M}\right)^2\right]$$
#
# 讀法：一個規模 $m_i$ 的地震，「預告」的主震規模是一個以
# $a_M+b_M m_i$ 為中心的常態分布。$b_M=1$ 且 $a_M\approx1$ 時，
# 中心就落在「比自己大一個規模單位」——與紐西蘭辨識前兆群震的
# 老經驗一致。**這是 EEPAS 與 ETAS 最刺眼的分歧**：ETAS 的後代
# 規模服從單調下降的 GR 密度（最可能比親代小），EEPAS 的目標規模
# 分布**有一個高於 $m_i$ 的峰值**。
#
# ### 時間核：變數變換與那個 $\ln 10$
#
# 時間的迴歸左邊是**對數**：$Y=\log_{10}T_P$。所以直接套用上面的
# 結論只能得到「$\log_{10}$ 前兆時間服從常態」，而我們要的是**前兆
# 時間本身**的密度。這中間差一個雅可比，而它正是各篇論文公式不
# 一致的那個地方，值得逐步做完。
#
# 記 $\tau=t-t_i$ 為等待時間（天），$u=\log_{10}\tau$，
# $\mu_{T,i}=a_T+b_T m_i$。由迴歸，
#
# $$U \sim N(\mu_{T,i},\ \sigma_T^2),\qquad
#   \varphi_U(u)=\frac{1}{\sigma_T\sqrt{2\pi}}
#     \exp\!\left[-\frac12\left(\frac{u-\mu_{T,i}}{\sigma_T}\right)^2\right].$$
#
# 變換是 $\tau=10^{u}$，即 $u=\ln\tau/\ln 10$，這是單調可微的一對一
# 映射，於是
#
# $$\begin{aligned}
# f(\tau)
#   &= \varphi_U\bigl(u(\tau)\bigr)\left|\frac{\mathrm{d}u}{\mathrm{d}\tau}\right|,
#   \qquad
#   \frac{\mathrm{d}u}{\mathrm{d}\tau}
#     = \frac{\mathrm{d}}{\mathrm{d}\tau}\frac{\ln\tau}{\ln 10}
#     = \frac{1}{\tau\ln 10}, \\[2pt]
#   &= \frac{1}{\tau\,\sigma_T\ln 10\,\sqrt{2\pi}}
#      \exp\!\left[-\frac12\left(
#        \frac{\log_{10}\tau-a_T-b_T m_i}{\sigma_T}\right)^2\right].
# \end{aligned}$$
#
# 補上因果性（前兆必須在目標之前）就得到 EEPAS 的時間核：
#
# $$f(t\mid t_i,m_i)=\frac{H(t-t_i)}{(t-t_i)\,\sigma_T\ln 10\,\sqrt{2\pi}}
#   \exp\!\left[-\frac12\left(
#     \frac{\log_{10}(t-t_i)-a_T-b_T m_i}{\sigma_T}\right)^2\right]$$
#
# $H$ 是單位階梯函數。**分母那個 $(t-t_i)\,\sigma_T\ln 10$ 就是
# 雅可比**：$(t-t_i)$ 來自 $\mathrm{d}u/\mathrm{d}\tau$ 的 $1/\tau$，
# $\ln 10$ 來自「迴歸用的是以 10 為底的對數」。檢核很容易：
# $\int_0^\infty f\,\mathrm{d}\tau=\int_{-\infty}^{\infty}\varphi_U\,
# \mathrm{d}u=1$，正規化成立。若把 $\ln 10$ 漏掉，整條密度會被放大
# $\ln 10\approx2.3026$ 倍，積分等於 2.3026 而不是 1。
#
# 這不是假想的風險。**2004 年的原始論文式 (6) 就漏了這個 $\ln 10$**，
# 2007 年以後的論文才補上；Biondini et al.（2023）在附錄 A 重新完整
# 推導整組公式時明白寫道，EEPAS 系列論文之間「有些含有 typo，使得
# 各篇公式不完全一致」。實作之前務必核對版本——這是本章第一條
# 實務守則。
#
# 對數常態有三個特徵時間，讀參數表時常常搞混（推導見附錄 A）：
# 中位數 $10^{\mu_{T,i}}$、眾數 $10^{\mu_{T,i}-\sigma_T^2\ln 10}$、
# 平均 $10^{\mu_{T,i}}\exp(\sigma_T^2\ln^2\!10/2)$。$\sigma_T$ 一大，
# 三者可以差好幾倍——**論文說的「前兆時間」幾乎都是中位數**。
#
# ### 空間核：從面積迴歸到二維常態
#
# 第三條迴歸的左邊是**面積**的對數。面積不是位置，所以這裡多一層
# 轉譯：把「前兆活動散布在一塊面積 $A_P$ 的區域內」翻譯成一個圓
# 對稱二維常態，讓它的尺度平方正比於 $A_P$。取每軸標準差
# $\sigma_i$，令 $\sigma_i^2=\sigma_A^2\,10^{b_A m_i}$：
#
# $$h(x,y\mid x_i,y_i,m_i)=\frac{1}{2\pi\sigma_A^2 10^{b_A m_i}}
#   \exp\!\left[-\frac{(x-x_i)^2+(y-y_i)^2}
#                    {2\sigma_A^2 10^{b_A m_i}}\right]$$
#
# 式子裡看不到 $a_A$，因為**它的角色被 $\sigma_A$ 吸收了**。轉換
# 常數可以直接算出來：若「最小涵蓋矩形」大致是每軸 $\pm k\sigma_i$
# 的方框，其面積為 $(2k\sigma_i)^2$，令它等於
# $A_P=10^{a_A+b_A m_i}$ 就得到
#
# $$\sigma_A=\frac{10^{a_A/2}}{2k}.$$
#
# 理論上 $k\approx1.5$–$2$（即比例常數 $1/(2k)\approx0.25$–$0.33$），
# 而 2004 年實際擬合出的比例是 0.37——**EEPAS 攤開的機率雲比
# 最小涵蓋矩形略寬一點**，合理，因為它還要吸收位置本身的不確定性。
#
# 把三個核畫在一起，整個模型的尺度感就出來了：

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import norm

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

# 紐西蘭 AMC 目錄的 EEPAS_1F 擬合值（Rhoades et al. 2022），全章示範用
P = dict(a_M=1.00, b_M=1.00, sig_M=0.20,
         a_T=1.97, b_T=0.35, sig_T=0.20,
         b_A=0.59, sig_A=0.51)
LN10 = np.log(10.0)
DAY_YR = 365.25


def g_mag(m, mi, p=P):
    """規模核：常態，中心 a_M + b_M·m_i。"""
    return norm.pdf(m, p["a_M"] + p["b_M"] * mi, p["sig_M"])


def f_time(tau_d, mi, p=P):
    """時間核：對數常態（tau_d 以天為單位），含 ln 10 的雅可比。"""
    mu = p["a_T"] + p["b_T"] * mi
    z = (np.log10(tau_d) - mu) / p["sig_T"]
    return np.exp(-0.5 * z ** 2) / (tau_d * p["sig_T"] * LN10 * np.sqrt(2 * np.pi))


def sigma_space(mi, p=P):
    """空間核每軸標準差（km）。"""
    return p["sig_A"] * 10 ** (0.5 * p["b_A"] * mi)


mi_show = [4.0, 5.0, 6.0]
mgrid = np.linspace(3.0, 8.5, 400)
tgrid = np.logspace(1.0, 5.0, 400)              # 10 天到 10 萬天
rgrid = np.linspace(0.0, 160.0, 400)            # km

fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.075,
                    subplot_titles=("規模核 g(m | mᵢ)",
                                    "時間核 f(t | tᵢ, mᵢ)",
                                    "空間核（徑向）h(r | mᵢ)"))
for k, mi in enumerate(mi_show):
    col = PALETTE[k]
    sig = sigma_space(mi)
    fig.add_trace(go.Scatter(x=mgrid, y=g_mag(mgrid, mi), mode="lines",
                             name=f"mᵢ = {mi:.0f}", legendgroup=f"m{k}",
                             line=dict(color=col, width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=tgrid / DAY_YR, y=f_time(tgrid, mi) * DAY_YR,
                             mode="lines", showlegend=False, legendgroup=f"m{k}",
                             line=dict(color=col, width=2)), row=1, col=2)
    fig.add_trace(go.Scatter(x=rgrid, y=2 * np.pi * rgrid
                             * np.exp(-rgrid ** 2 / (2 * sig ** 2))
                             / (2 * np.pi * sig ** 2),
                             mode="lines", showlegend=False, legendgroup=f"m{k}",
                             line=dict(color=col, width=2)), row=1, col=3)
med_yr = [10 ** (P["a_T"] + P["b_T"] * mi) / DAY_YR for mi in mi_show]
fig.update_xaxes(title_text="目標規模 m", row=1, col=1)
fig.update_xaxes(title_text="等待時間（年，對數軸）", type="log", row=1, col=2)
fig.update_xaxes(title_text="距離 r（km）", row=1, col=3)
fig.update_yaxes(title_text="機率密度", row=1, col=1)
apply_layout(fig, title=f"三個核隨前兆規模的變化（中位等待時間 "
                        f"{med_yr[0]:.1f} / {med_yr[1]:.1f} / {med_yr[2]:.0f} 年；"
                        f"每軸空間尺度 {sigma_space(4.0):.1f} / "
                        f"{sigma_space(5.0):.1f} / {sigma_space(6.0):.0f} km）",
             hovermode="closest", height=400)
fig

# %% [markdown]
# 三張圖說的是同一件事：**規模愈大的前兆，預告的主震愈大、要等得
# 愈久、可能發生的範圍愈廣**。空間那一欄畫的是徑向機率密度
# $2\pi r\,h(r)$（對 $r$ 積分為 1），所以峰值不在原點——最可能的
# 距離是 $r=\sigma_i$，這是二維常態的標準結果。
#
# ### 對照組：EEPAS 的時間核 vs ETAS 的時間核
#
# 這一章最有教學力的一張圖，是把 EEPAS 的對數常態時間核與第 12 章
# 的正規化 Omori 密度 {eq}`eq:omori-density` 放在同一組座標軸上。
# 兩者都是**正規化的等待時間密度**，可以直接比：

# %% tags=["hide-input"]
C_OM, P_OM = 0.01, 1.15                        # Omori 核（天）：ETAS 的典型值
t_cmp = np.logspace(-2.0, 5.0, 600)            # 0.01 天到 10 萬天
g_omori = (P_OM - 1) / C_OM * (1 + t_cmp / C_OM) ** (-P_OM)
f_eepas = f_time(t_cmp, 5.0)
med_ee = 10 ** (P["a_T"] + P["b_T"] * 5.0)
med_om = C_OM * (2 ** (1 / (P_OM - 1)) - 1)
cross = t_cmp[np.argmin(np.abs(np.log(f_eepas / g_omori))[t_cmp > 1.0])
              + int((t_cmp <= 1.0).sum())]

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_cmp, y=g_omori, mode="lines",
                         name=f"ETAS：Omori 冪次（c = {C_OM}, p = {P_OM}）",
                         line=dict(color=PALETTE[1], width=2.5)))
fig.add_trace(go.Scatter(x=t_cmp, y=f_eepas, mode="lines",
                         name="EEPAS：對數常態（mᵢ = 5）",
                         line=dict(color=ACCENT, width=2.5)))
fig.add_vline(x=cross, line_dash="dot", line_color="#666666")
fig.add_annotation(x=np.log10(cross), y=-3.2, text=f"交叉 ≈ {cross / DAY_YR:.1f} 年",
                   showarrow=True, ax=55, ay=-20)
fig.update_xaxes(type="log", title_text="等待時間（天，對數軸）")
fig.update_yaxes(type="log", title_text="機率密度（每天）",
                 range=[-9.5, 1.5])
apply_layout(fig, title=f"兩種時間核：冪次重尾 vs 對數常態峰"
                        f"（中位等待 {med_om:.3f} 天 vs "
                        f"{med_ee / DAY_YR:.1f} 年）",
             hovermode="x", height=460)
fig

# %% [markdown]
# 這張圖把「短期看觸發、中長期看前兆」這句話變成幾何事實。
#
# **Omori 核（橘）在 $t\to0$ 發散**，一半的機率質量落在頭
# 幾分鐘之內；它從第一天起就單調下降，成為一條直線（在對數—對數
# 座標上冪次律就是直線），永遠沒有峰值。**對數常態核（藍）在
# $t\to0$ 趨近於零**——它的指數項壓過了 $1/\tau$ 的發散，前幾天
# 幾乎不貢獻任何機率；它先爬升、在十幾年處達到峰值、然後兩邊
# 對稱地衰減（在對數時間軸上）。
#
# 兩條線在數年的尺度上交叉。交叉點左邊，任何叢集現象都應該優先用
# 觸發解釋；交叉點右邊，Omori 核已經衰減到可以忽略，而 EEPAS 核
# 才剛剛進入它的黃金時段。**這條交叉線就是 16.9 節那個「三個月
# ETAS 勝、六個月以上 EEPAS 勝」實證結果的幾何解釋。**
#
# 三個延伸觀察。**其一**，兩者的形狀差異不是模型偏好，而是來源
# 不同：Omori 是直接對餘震衰減曲線擬合出來的經驗律，對數常態是
# 一條**線性迴歸的殘差**經過變數變換的結果。**其二**，冪次律沒有
# 特徵尺度（自相似），對數常態有——峰值位置隨 $m_i$ 移動，這正是
# 「尺度相稱」四個字的數學內容。**其三**，兩條密度都對時間積分
# 為 1，但它們前面乘的**總量**完全不同：ETAS 乘的是產能
# $\kappa(m_i)$，EEPAS 乘的是 $\eta(m_i)w_i$——下一節就處理這件事。

# %% [markdown]
# ## 16.3 總率密度與 $\eta(m)$：為什麼需要正規化
#
# ### 疊加
#
# 把每個過去地震的貢獻寫成三個核的乘積，
# $\lambda_i=w_i\,f(t\mid t_i,m_i)\,g(m\mid m_i)\,h(x,y\mid x_i,y_i,m_i)$，
# 再加上一個背景項，就得到 EEPAS 的**總率密度**（對時間、規模、
# 面積三重微分後的期望個數）：
#
# $$\lambda(t,m,x,y)=\mu_E\,\lambda_{\rm PPE}(t,m,x,y)
#   +\sum_{t_i\ge t_0,\ m_i\ge m_0}
#     \eta(m_i)\,w_i\,f\,g\,h$$ (eq:eepas-rate)
#
# 這裡 $\mu_E$ 是**混合權重**（第 15 章與文獻寫作 $\mu$；本部把
# $\mu$ 保留給 ETAS 的背景率，見 10.8 節的記號表），可以讀成
# 「沒有可辨識前兆的目標地震所佔的比例」。它同時是模型的誠實
# 儀表板：時變項找不到東西時，最大概似會自動把 $\mu_E$ 推高、
# 轉而依賴背景。
#
# 剩下 $\eta(m_i)$ 這個因子。它不是自由參數，而是被一個要求逼出來
# 的——**長期平均的規模分布必須回到 GR 律**。
#
# ### 為什麼不能省略
#
# 先看省略會怎樣。假設 $\eta\equiv1$，且輸入地震本身服從 GR：
# 規模密度 $\propto\beta e^{-\beta(\nu-m_0)}$，$\beta=b\ln 10$。
# 那麼長期平均的規模率密度是
#
# $$\bar\lambda(m)\ \propto\ \int_{m_0}^{\infty}
#   g(m\mid \nu)\,\beta e^{-\beta(\nu-m_0)}\,\mathrm{d}\nu .$$
#
# 換元 $u=a_M+b_M\nu$（$\mathrm{d}\nu=\mathrm{d}u/b_M$），指數項變成
# $e^{-(\beta/b_M)u}\times$ 常數。也就是說，**輸出的 GR 斜率是
# $b/b_M$，不是 $b$**。$b_M=0.65$ 時輸出的 b 值變成 1.54——模型會
# 系統性地少預報大地震。$\eta$ 的職責就是把這個扭曲扳回來。
#
# ### 推導
#
# 要求：$\bar\lambda(m)\propto e^{-\beta(m-m_0)}$。試探解取
# $\eta(\nu)=c_0\,e^{-\beta(b_M-1)\nu}$，代入：
#
# $$\begin{aligned}
# \bar\lambda(m)
#   &\propto \int_{m_0}^{\infty} c_0 e^{-\beta(b_M-1)\nu}\,
#      g(m\mid\nu)\,\beta e^{-\beta(\nu-m_0)}\,\mathrm{d}\nu \\
#   &= c_0\,\beta\,e^{\beta m_0}\int_{m_0}^{\infty}
#      e^{-\beta b_M \nu}\,g(m\mid\nu)\,\mathrm{d}\nu .
# \end{aligned}$$
#
# 換元 $u=a_M+b_M\nu$，於是 $e^{-\beta b_M\nu}=e^{-\beta(u-a_M)}$，
# 而 $g$ 變成以 $u$ 為中心、寬 $\sigma_M$ 的常態在 $m$ 處取值：
#
# $$\bar\lambda(m)\ \propto\ \frac{c_0\beta e^{\beta(m_0+a_M)}}{b_M}
#   \int_{a_M+b_M m_0}^{\infty} e^{-\beta u}\,
#   \phi_{\sigma_M}(m-u)\,\mathrm{d}u .$$
#
# **先把積分下限放寬到 $-\infty$**（放寬造成的缺口就是下一節的
# $\Delta(m)$）。令 $v=m-u$，$V\sim N(0,\sigma_M^2)$：
#
# $$\int_{-\infty}^{\infty} e^{-\beta u}\phi_{\sigma_M}(m-u)\,\mathrm{d}u
#   = e^{-\beta m}\,\mathbb{E}\!\left[e^{\beta V}\right]
#   = e^{-\beta m}\,e^{\sigma_M^2\beta^2/2},$$
#
# 最後一步用的是常態的**動差生成函數**
# $\mathbb{E}[e^{sV}]=e^{s^2\sigma^2/2}$。整理：
#
# $$\bar\lambda(m)\ \propto\
#   \frac{c_0\beta}{b_M}\,e^{\beta(m_0+a_M)}\,
#   e^{\sigma_M^2\beta^2/2}\,e^{-\beta m}.$$
#
# 指數部分已經是純粹的 $e^{-\beta m}$——**GR 律恢復了**。把前面那
# 兩個常數吸收進 $c_0$，並選定使總率守恆的比例常數，就得到
#
# $$\eta(m)=\frac{b_M\,(1-\mu_E)}{\mathbb{E}(w)}\,
#   \exp\!\left[-\beta\left(a_M+(b_M-1)m
#     +\frac{\sigma_M^2\beta}{2}\right)\right]$$
#
# 三個因子各有職責，值得逐一點名。$e^{-\beta(b_M-1)m}$ 修正斜率
# （$b_M<1$ 時它隨 $m$ **上升**，因為小前兆被過度代表了，要把大
# 前兆的權重補回去）；$e^{-\beta\sigma_M^2\beta/2}$ **正好抵掉動差
# 生成函數吐出來的 $e^{+\sigma_M^2\beta^2/2}$**——規模核愈寬，它
# 把機率質量往高規模端搬得愈多，就要扣得愈多；
# $b_M(1-\mu_E)/\mathbb{E}(w)$ 讓總率守恆：時變項承擔 $1-\mu_E$
# 的份額、背景項承擔 $\mu_E$，而除以平均權重 $\mathbb{E}(w)$ 是
# 因為每個貢獻已經先被 $w_i$ 打了折。
#
# ### 退化情形
#
# 取 $b_M=1$、$\mu_E=0$、$w_i\equiv1$（於是 $\mathbb{E}(w)=1$）：
#
# $$\eta=\exp\!\left[-\beta\left(a_M+\frac{\sigma_M^2\beta}{2}\right)\right]
#   \ =\ \text{常數}.$$
#
# 這是實務上最常見的設定（16.10 節那張表裡四組參數有三組把 $b_M$
# 釘在 1），此時 $\eta$ 退化成一個純粹的尺度因子，只影響總量、
# 不影響形狀。**留意最後那項必須是 $\sigma_M^2\beta/2$ 而不是
# $\sigma_M^2/2$**——文獻裡兩種寫法都出現過，而只有前者能與上面的
# 動差生成函數對消。這是本章第二條實務守則的實例。
#
# 用數值積分把「有沒有 $\eta$」的差別畫出來，取 $b_M=0.65$ 讓效應
# 看得見：

# %% tags=["hide-input"]
B_GR = 1.0
BETA = B_GR * LN10
M0_IN, MU_IN = 2.0, 8.0                        # 前兆規模積分範圍
BM_DEMO, AM_DEMO, SM_DEMO = 0.65, 3.16, 0.35   # 取 b_M ≠ 1 才看得出效應

nu = np.linspace(M0_IN, MU_IN, 2400)
m_out = np.linspace(4.5, 7.5, 160)
w_gr = BETA * np.exp(-BETA * (nu - M0_IN))                    # 輸入 GR 密度
eta_nu = np.exp(-BETA * (AM_DEMO + (BM_DEMO - 1) * nu
                         + SM_DEMO ** 2 * BETA / 2))          # 正規化函數
gmat = norm.pdf(m_out[:, None], AM_DEMO + BM_DEMO * nu[None, :], SM_DEMO)
lam_raw = np.trapezoid(gmat * w_gr[None, :], nu, axis=1)
lam_eta = np.trapezoid(gmat * eta_nu[None, :] * w_gr[None, :], nu, axis=1)
slope_raw = -np.polyfit(m_out, np.log10(lam_raw), 1)[0]
slope_eta = -np.polyfit(m_out, np.log10(lam_eta), 1)[0]

fig = go.Figure()
for y, lab, col, dash in [
        (lam_raw / lam_raw[0], f"η ≡ 常數（等效 b = {slope_raw:.2f}）",
         PALETTE[1], "solid"),
        (lam_eta / lam_eta[0], f"乘上 η(m)（等效 b = {slope_eta:.2f}）",
         ACCENT, "solid"),
        (10 ** (-B_GR * (m_out - m_out[0])), f"輸入的 GR 律（b = {B_GR:.2f}）",
         "#666666", "dash")]:
    fig.add_trace(go.Scatter(x=m_out, y=y, mode="lines", name=lab,
                             line=dict(color=col, width=2.5, dash=dash)))
fig.update_yaxes(type="log", title_text="長期平均規模率密度（相對值）")
apply_layout(fig, title=f"η(m) 的職責：把長期 FMD 的斜率從 b/b_M 扳回 b"
                        f"（示範用 b_M = {BM_DEMO}）",
             xaxis_title="目標規模 m", hovermode="x", height=430)
fig

# %% [markdown]
# 橘線明顯比虛線陡：$b_M=0.65$ 時，未正規化的長期分布等效 b 值
# 接近 $1/0.65\approx1.54$，模型會在 M7 那一端少預報一個數量級。
# 藍線與灰虛線幾乎重合——**$\eta$ 讓 EEPAS 的長期平均行為與輸入
# 目錄的 GR 律一致**。
#
# 這件事在哲學上比在數值上更重要。EEPAS 是一個**時變**模型，但它
# 的時變只能是「把機率在時空中重新分配」，不能憑空製造或消滅
# 地震。$\eta$ 就是那條守恆律的執行者：**模型可以說某年某地機率
# 高三倍，但它不能說整個地區的長期地震率不是 GR 律說的那樣。**

# %% [markdown]
# ## 16.4 $\Delta(m)$：規模不完整的補償
#
# 上一節推導 $\eta$ 時偷渡了一步：把積分下限從 $a_M+b_M m_0$ 放寬到
# $-\infty$。現實裡放寬不了——**只有 $m_i\ge m_0$ 的地震進得了目錄**，
# 比門檻小的前兆存在，但我們看不見。缺掉的那一塊要補回來。
#
# ### 缺口的大小
#
# 定義補償函數為「實際拿得到的貢獻」除以「理論上應有的貢獻」：
#
# $$\Delta(m)=\frac{\displaystyle\int_{a_M+b_M m_0}^{\infty}
#     e^{-\beta u}\,\phi_{\sigma_M}(m-u)\,\mathrm{d}u}
#   {\displaystyle\int_{-\infty}^{\infty}
#     e^{-\beta u}\,\phi_{\sigma_M}(m-u)\,\mathrm{d}u}.$$
#
# 分子分母的被積函數相同，只差積分下限，所以 $\Delta(m)$ 就是一個
# **尾機率**——只要認出被積函數正比於哪個機率密度即可。
#
# ### 推導：指數傾斜
#
# 把被積函數的指數部分配方（完整代數見附錄 C）：
#
# $$\begin{aligned}
# -\beta u-\frac{(u-m)^2}{2\sigma_M^2}
#   &= -\frac{1}{2\sigma_M^2}
#      \left[u^2-2\bigl(m-\sigma_M^2\beta\bigr)u+m^2\right] \\
#   &= -\frac{\bigl(u-(m-\sigma_M^2\beta)\bigr)^2}{2\sigma_M^2}
#      + \text{（與 }u\text{ 無關的項）}.
# \end{aligned}$$
#
# 也就是說，$e^{-\beta u}\phi_{\sigma_M}(m-u)$ 正比於一個
# $N(m-\sigma_M^2\beta,\ \sigma_M^2)$ 的密度。這叫**指數傾斜**
# （exponential tilting）：乘上 $e^{-\beta u}$ 不改變常態的形狀，
# 只把中心往左搬 $\sigma_M^2\beta$。與 $u$ 無關的那一項在分子分母
# 中相同，直接消掉。於是
#
# $$\Delta(m)=P\bigl(U\ge a_M+b_M m_0\bigr),
#   \qquad U\sim N\!\left(m-\sigma_M^2\beta,\ \sigma_M^2\right),$$
#
# $$\Delta(m)=1-\Phi\!\left(\frac{a_M+b_M m_0-m+\sigma_M^2\beta}{\sigma_M}\right)
#   =\Phi\!\left(\frac{m-a_M-b_M m_0-\sigma_M^2\beta}{\sigma_M}\right).$$
#
# 補償的做法是把每個前兆的貢獻除以 $\Delta(m)$：
#
# $$\lambda(t,m,x,y)=\mu_E\lambda_{\rm PPE}
#   +\sum_{t_i\ge t_0,\ m_i\ge m_0}
#     \eta(m_i)\,\frac{w_i\,f\,g\,h}{\Delta(m)} .$$
#
# ### 怎麼讀這個數字
#
# $\Delta$ 的自變數是**目標規模** $m$，不是前兆規模。它的意義是：
# 「要預報一個規模 $m$ 的地震，理論上該有的前兆裡，有多大比例的
# 規模在 $m_0$ 之上、因而被目錄看得見？」
#
# 拿一組真實參數試算（$a_M=1.0$、$b_M=1$、$\sigma_M=0.24$、
# $b=1.08$ 即 $\beta=2.49$、$\sigma_M^2\beta=0.14$）：目標比門檻高
# 2.0 個規模單位時 $\Delta\approx0.9998$，幾乎沒有缺口；高 1.5 時
# $\Delta\approx0.93$；高 1.2 時掉到 $\approx0.59$——**近半的前兆
# 貢獻在門檻底下，模型看不見。**
#
# 這解釋了一條實作規範：**$m_T-m_0$ 不應遠小於 2.0**。理由現在是
# 定量的而不是口訣——$b_M\approx1$、$a_M\approx1$ 意味著典型前兆
# 比目標小約一個單位，而前兆規模本身還有 $\pm2\sigma_M$ 的散布，
# 再加上 $\sigma_M^2\beta$ 這個傾斜位移，兩個單位剛好把整個分布
# 包進來。門檻抬得太高，$\Delta$ 一小，$1/\Delta$ 就會把少數看得見
# 的前兆放大到不合理的權重——**補償函數不是萬靈丹，它只在缺口
# 不大時可信。**

# %% [markdown]
# ## 16.5 $w_i$：餘震降權，其實是第 14 章的老朋友
#
# 三個核裡沒有 $w_i$，因為它不來自 Ψ，而來自一個實務問題：**EEPAS
# 的理論基礎是「獨立事件的前兆活動」，但目錄裡大部分地震是餘震。**
# 一場 M7 之後的幾千個餘震，如果每一個都放下一個「預告更大地震」
# 的機率包裹，模型會在剛剛破裂完的區域堆出一座虛假的高峰。
#
# ### 定義
#
# 解法是先用一個 ETAS 型的餘震模型 $\lambda'$ 把目錄過一遍，
# 再取
#
# $$w_i=\frac{\nu\,\lambda_0(t_i,m_i,x_i,y_i)}
#            {\lambda'(t_i,m_i,x_i,y_i)} .$$
#
# 分母 $\lambda'$ 是「背景 + 被觸發」的總率密度，分子
# $\nu\lambda_0$ 是其中的背景成分（$\nu$ 在這裡是餘震模型裡背景項
# 的比例常數，與 18.4 節 Molchan 圖的漏報率 $\nu$ 無關）。所以
#
# > $w_i$ **就是「事件 $i$ 是獨立事件的機率」**；而一群事件的平均
# > 權重 $\overline{w}$ **就是這群事件裡獨立事件的比例** $p_I$。
#
# ### 這正是隨機除叢
#
# 把它與第 14 章的 {eq}`eq:rho-phi` 並排看：$\phi_j=\mu(x_j,y_j)/
# \lambda^*(t_j,x_j,y_j)$，同樣是「背景成分 ÷ 總強度」。**$w_i$ 與
# $\phi_j$ 是同一個量**，都是疊加定理給出的成分歸屬機率，都不做
# 二分、只給機率。第 14 章那句「強度分解到哪裡，機率就分配到哪裡」
# 原封不動適用。
#
# 三個差異值得標出來。**其一**，背景項不同：ETAS 用變頻寬核從資料
# 估出 $\mu(x,y)$，EEPAS 直接用 PPE 當 $\lambda_0$——省一次迭代，
# 代價是背景場的形狀被 PPE 的平滑距離綁死。**其二**，餘震模型多了
# 一條非標準的限制：**餘震至少比親代小 $\delta$ 個規模單位**
# （2004 年擬合出 $\delta=0.7$）。第 14 章強調過標準 ETAS 對後代
# 規模**沒有上限**，這正是它能自然重現前震的原因；EEPAS 刻意把
# 這條路堵住，因為「比親代大的後續事件」正是它要留給前兆機制去
# 解釋的東西。**兩個模型在同一個現象上劃出了不同的管轄權。**
# **其三**，$w_i$ 只在**擬合與預報時**當權重用，不回頭改寫目錄。
#
# ### 一個誠實的麻煩
#
# $w_i$ 在理論上無懈可擊，在實證上卻反覆惹麻煩。2007 年南加州那
# 篇把目標規模下修到 M4.95，結果發現**等權重（$w_i\equiv1$）明顯
# 優於餘震降權**，$\ln L$ 差距約 20。2011 年日本的分解更刺眼：
# 獨立事件（$w_i=1$）的平均機率增益 2.3，而**純餘震（$w_i=0$）也
# 有 1.7**。
#
# 作者自己的診斷值得整段抄下來：目標規模一低，目標集裡就混進大量
# 餘震；最大概似為了預報這些餘震，會把參數推去遷就它們，於是
# 「目前的實作已偏離了模型的理論基礎」。**最佳化目標與模型理念
# 不一致時，參數會說謊**——這是第 14 章那句「模型評分不等於模型
# 正確」在另一個模型上的重演。

# %% [markdown]
# ## 16.6 PPE 與 SUP：基準模型的三個角色
#
# {eq}`eq:eepas-rate` 裡的 $\lambda_{\rm PPE}$ 到現在還沒定義。它是
# **PPE**（Proximity to Past Earthquakes，「與過去地震的鄰近性」，
# 改寫自 Jackson & Kagan 1999），一個樸素到近乎無聊、卻誰都繞不開
# 的模型。
#
# ### 三個因子
#
# PPE 是完全可分離的：$\lambda_{\rm PPE}=f_0(t)\,g_0(m)\,h_0(x,y)$。
#
# **空間**：把過去每個震央攤上一個平滑核再疊加，
#
# $$h_{0i}(r_i)=a\,(m_i-m_c)
#   \left[\frac{1}{\pi\,(d^2+r_i^2)}+s\right] .$$
#
# $d$ 是**平滑距離**（義大利擬合出約 30 km、紐西蘭 AMC 約 5 km），
# $s$ 是一個小常數，代表「遠離所有歷史震央的地方仍然可能發震」的
# 底線。權重 $(m_i-m_c)$ 讓大地震的貢獻更重。核的尾巴是
# $r^{-2}$，**在無窮遠處不可積**（$\int_0^R 2\pi r\,\mathrm{d}r/
# [\pi(d^2+r^2)]=\ln(1+R^2/d^2)$ 隨 $R$ 對數發散），所以 $a$ 是由
# 「在有限的預報區域上正規化」定出來的，不是一個獨立的物理量。
# 這個重尾是刻意的：**未來地震落在歷史震央附近，但「附近」的
# 尺度不該被高斯尾巴切斷。**
#
# **規模**：$g_0(m)=\beta e^{-\beta(m-m_c)}$，就是 GR 律。
#
# **時間**：$f_0(t)=1/(t-t_0)$，$t_0$ 是目錄起點。這一項最容易被
# 誤讀。它不是「率隨時間衰減」，而是**把累積的貢獻除以經過的時間
# ，得到一個滾動平均率**。所以 PPE 的率密度在兩次地震之間隨時間
# 緩慢下降（分母變大），一有新地震就往上跳一階（分子加一項）。
# **PPE 因此不是一個時間獨立模型，而是一個會自我更新的長期率
# 估計器**——這也是為什麼把它當基準時，它並不是稻草人。
#
# ### SUP：刻度的原點
#
# 還有一個更樸素的模型：**SUP**（stationary uniform Poisson），
# 時間齊次、空間均勻、規模服從 GR，率由目錄總數定。它幾乎不可能
# 贏過任何東西，但它扮演**共同刻度原點**的角色。EEPAS 文獻的評分
# 一律寫成
#
# $$I_M=\frac{\mathrm{AIC}_{\rm SUP}-\mathrm{AIC}_M}{2N},$$
#
# 也就是「模型 $M$ 相對 SUP 的每地震資訊增益」。任意兩個模型的
# 比較就是 $I_A-I_B$，SUP 消掉。**用同一個爛模型當原點，不同論文
# 的數字才能相減。**
#
# ### 為什麼「贏過平滑地震度多少」是第一個該問的問題
#
# PPE 幾乎沒有物理，卻是地震預報界最重要的一類模型，因為它扮演
# 三個角色：{eq}`eq:eepas-rate` 的**背景項**、$w_i$ 分子裡的
# **獨立事件參照**，以及最重要的——**所有花俏模型都必須贏過的
# 基準線**。
#
# 理由很硬：地震的空間分布極度不均勻，光是「把歷史震央抹平」就能
# 得到巨大的資訊增益。一個模型若只報告「我比均勻分布好很多」，
# 它報告的可能只是台灣的地震集中在東部這件常識。真正的問題是
# **扣掉這個常識之後還剩多少**——這就是第 18 章 IGPE 的精神。
# 用台灣長期目錄的 M≥5 事件把 PPE 的空間項畫出來：

# %% tags=["hide-input"]
cat_long = load_taiwan_catalog(min_ml=5.0)
D_KM, S_BG, MC_PPE = 15.0, 1e-4, 5.0           # 平滑距離、遠域常數、參照門檻
lons = np.arange(119.0, 123.5, 0.1)
lats = np.arange(21.0, 26.0, 0.1)
LON, LAT = np.meshgrid(lons, lats)
dens = np.zeros_like(LON)
ev = cat_long[["longitude", "latitude", "ML"]].to_numpy()
for lo, la, mi in ev:
    r2 = ((LON - lo) * 111 * np.cos(np.radians(la))) ** 2 + ((LAT - la) * 111) ** 2
    dens += (mi - MC_PPE + 0.1) * (1 / (np.pi * (D_KM ** 2 + r2)) + S_BG)

fig = go.Figure(go.Heatmap(x=lons, y=lats, z=np.log10(dens), colorscale="Blues",
                           colorbar=dict(title="log₁₀ 相對率")))
apply_layout(fig, title=f"PPE 的空間項 h₀（台灣 1973–2025，M ≥ {MC_PPE:.0f} 共 "
                        f"{len(ev)} 筆，平滑距離 d = {D_KM:.0f} km，"
                        f"動態範圍 {np.ptp(np.log10(dens)):.1f} 個數量級）",
             xaxis_title="經度（°E）", yaxis_title="緯度（°N）",
             yaxis_scaleanchor="x", hovermode="closest", height=560)
fig

# %% [markdown]
# 東部外海與西南部亮起來，中央山脈與西部平原暗下去——這就是
# 「地震會發生在以前發生過地震的地方附近」的定量版本。整張圖的
# 動態範圍有好幾個數量級，而這幾個數量級**全部是 EEPAS 必須先
# 免費繼承、然後才能開始談自己貢獻了什麼的東西**。
#
# 兩個但書。$d=15$ km 是示範值，正式應用要與 $\beta$、$s$ 一起用
# 最大概似估；而這裡只畫了空間項，完整的 PPE 還要乘上
# $f_0(t)g_0(m)$。

# %% [markdown]
# ## 16.7 EEPAS 家族與兩個漂亮的簡化
#
# 二十年下來，EEPAS 長出了一整族變體。每一個變體都對應一個被指出
# 的缺陷，這張表因此也是一份缺陷清單：
#
# | 變體 | 加了什麼 | 主要結論 |
# |---|---|---|
# | 標準 EEPAS | 基準版本 | 紐西蘭擬合、加州獨立測試皆勝 PPE |
# | STEP-EEPAS | 與短期叢集模型凸組合 | 0.42/0.58 配比，24 小時預報增益 > 2 |
# | EAS | 把預期主震的餘震一併積分 | 比 EEPAS 再進步 $\Delta\ln L/N\approx0.1$ |
# | ERDEEP | 令 $a_T,\sigma_A$ 隨長期地震率變 | 正式測試輸給標準版，假設被放棄 |
# | Janus | 與 ETAS 凸組合，權重 $\pi_{\rm ETAS}$ | horizon = 0 仍勝純 ETAS |
# | LEEPAS | 補償 time-lag 造成的前兆缺失 | lag 15 年仍有 $IG\approx0.5$ |
# | FLEEPAS | 固定 lead time 擬合 | lead time 越短，$a_T$ 越小 |
# | FLCEEPAS | 固定 lead time 且補償 | 兩種修正併用 |
#
# 其中 **EAS**（Rhoades 2009）補的是一個邏輯漏洞：EEPAS 預報的是
# **主震**，但它宣稱的率密度要對**所有**超過門檻的地震負責——包括
# 那些尚未發生的主震自己的餘震。EAS 的補法是把餘震模型對「所有
# 可能的未來主震」積分：
#
# $$\lambda_A(t_a,m_a,x_a,y_a)=\int\!\!\int\!\!\iint
#   f_A\,g_A\,h_A\;\lambda_M(t,m,x,y)\,\mathrm{d}y\,\mathrm{d}x
#   \,\mathrm{d}m\,\mathrm{d}t$$
#
# 這個四重積分看起來很嚇人，但其中兩重可以被兩個漂亮的近似消掉。
#
# ### 簡化一：長期預報不必區分主震與餘震的時間分布
#
# 時間那一重是一個卷積。令 $s=t_a-t$：
#
# $$\int_{t_i}^{t_a} f_A(t_a-t)\,f_{1i}(t)\,\mathrm{d}t
#   = \int_0^{t_a-t_i} f_A(s)\,f_{1i}(t_a-s)\,\mathrm{d}s .$$
#
# 關鍵在**兩個時間尺度差了三個數量級**：$f_A$ 是 Omori 核，質量
# 集中在數天到數月；$f_{1i}$ 是 EEPAS 的對數常態核，特徵尺度是
# 數年到數十年（16.2 節那張對照圖就是這件事的圖像）。把 $f_{1i}$
# 在 $t_a$ 附近展開：
#
# $$f_{1i}(t_a-s)=f_{1i}(t_a)-s\,f_{1i}'(t_a)+O(s^2),$$
#
# 代回並用 $\int_0^\infty f_A(s)\,\mathrm{d}s=1$：
#
# $$\int f_A(s)f_{1i}(t_a-s)\,\mathrm{d}s
#   \approx f_{1i}(t_a)-\bar{s}\,f_{1i}'(t_a)
#   = f_{1i}(t_a)\left[1-\bar{s}\,
#     \frac{\mathrm{d}\ln f_{1i}}{\mathrm{d}t}\Big|_{t_a}\right],$$
#
# 其中 $\bar{s}=\int s\,f_A(s)\,\mathrm{d}s$ 是平均餘震延遲。相對
# 誤差是 $\bar{s}\,|\mathrm{d}\ln f_{1i}/\mathrm{d}t|$；對數常態在
# 峰值附近的對數導數量級是 $1/\tau$（$\tau$ 為特徵等待時間），
# 所以誤差約為 $\bar{s}/\tau$——月比上年，**百分之幾**。於是
#
# $$\int f_A(t_a\mid t)\,f_{1i}(t)\,\mathrm{d}t \approx f_{1i}(t_a).$$
#
# 一句話：**在中長期預報的尺度上，「主震的時間」與「它的餘震的
# 時間」是同一件事**，時間那一重積分可以直接拿掉。
#
# ### 簡化二：兩個二維常態的卷積仍是二維常態
#
# 空間那一重更漂亮。EAS 的餘震空間核**不是冪次律，而是二維常態**
# （變異數 $\sigma_V^2 10^{m}$，符合 Utsu 的餘震面積關係），而
# EEPAS 的 $h_{1i}$ 也是二維常態。兩個圓對稱二維常態的卷積可以
# 逐軸做，每軸都是一維常態卷積（配方法見附錄 D）：
#
# $$\int \phi_{\sigma_2}(x_a-x)\,\phi_{\sigma_1}(x-x_i)\,\mathrm{d}x
#   = \phi_{\sqrt{\sigma_1^2+\sigma_2^2}}(x_a-x_i).$$
#
# 兩軸相乘就是
#
# $$\iint h_A\,h_{1i}\,\mathrm{d}y\,\mathrm{d}x
#   = \frac{1}{2\pi\bigl(\sigma_V^2 10^{m}+\sigma_A^2 10^{b_A m_i}\bigr)}
#     \exp\!\left[-\frac{(x_a-x_i)^2+(y_a-y_i)^2}
#       {2\bigl(\sigma_V^2 10^{m}+\sigma_A^2 10^{b_A m_i}\bigr)}\right].$$
#
# **變異數相加。** 這條式子有一個立即可讀的物理結論：
#
# > 預期餘震的空間散布**一定**大於預期主震的空間散布。
#
# 因為它必須同時吸收兩份不確定性——**未來主震會落在哪裡**
# （$\sigma_A^2 10^{b_A m_i}$），以及**餘震會散在主震周圍多遠**
# （$\sigma_V^2 10^{m}$）。這不是模型選擇，是機率論的結果。
#
# 順帶一提，EAS 刻意用常態而非冪次律當餘震空間核，作者給的理由
# 值得記：遠距觸發的事件，主震可能是扳機但不是主因；**它如果早就
# 在孕育中，長期模型本來就該預期得到它**，不需要靠餘震項再算一次。
#
# 把 {eq}`eq:eepas-rate` 的機制畫成圖，三個核相乘的結果就是每個
# 地震在未來時空中放下的一個「機率包裹」：

# %% tags=["hide-input"]
# 示意用 Evison–Rhoades 的 Ψ 迴歸值（第 15 章），凸顯尺度隨規模變化
PSI = dict(a_M=3.16, b_M=0.65, a_T=1.36, b_T=0.40, sig_T=0.35,
           b_A=0.35, sig_A=1.0)
events = [(1.0, 20.0, 4.0), (3.0, 55.0, 5.0), (5.5, 80.0, 5.8)]
tg = np.linspace(0.05, 25.0, 320)              # 年
xg = np.linspace(0.0, 110.0, 240)              # km（一維剖面）
TT, XX = np.meshgrid(tg, xg)
dens = np.zeros_like(TT)
for t0, x0, mi in events:
    dt_d = np.clip((TT - t0) * DAY_YR, 1e-3, None)
    z = (np.log10(dt_d) - PSI["a_T"] - PSI["b_T"] * mi) / PSI["sig_T"]
    f = np.exp(-0.5 * z ** 2) / (dt_d * PSI["sig_T"] * LN10 * np.sqrt(2 * np.pi))
    f[TT <= t0] = 0.0
    sx = PSI["sig_A"] * 10 ** (0.5 * PSI["b_A"] * mi)
    h = np.exp(-0.5 * ((XX - x0) / sx) ** 2) / (sx * np.sqrt(2 * np.pi))
    dens += f * h * 10 ** (-B_GR * (PSI["b_M"] - 1) * mi)      # η(m_i) 的斜率項

fig = go.Figure(go.Heatmap(x=tg, y=xg, z=np.sqrt(dens), colorscale="Blues",
                           showscale=False))
fig.add_trace(go.Scatter(
    x=[e[0] for e in events], y=[e[1] for e in events], mode="markers+text",
    name="過去的地震",
    text=[f"M{e[2]:.1f} → 預告 M{PSI['a_M'] + PSI['b_M'] * e[2]:.1f}，"
          f"中位等待 {10 ** (PSI['a_T'] + PSI['b_T'] * e[2]) / DAY_YR:.1f} 年"
          for e in events],
    textposition="middle left",
    marker=dict(size=[8, 12, 16], color=QUAKE_COLOR)))
apply_layout(fig, title="每個地震在未來時空放下一個機率包裹（顏色為率密度的平方根）",
             xaxis_title="時間（年）", yaxis_title="位置（km，一維剖面）",
             hovermode="closest", height=470)
fig

# %% [markdown]
# 三個地震（紅點）各自貢獻一團往未來延伸的密度。規模愈大，包裹
# 愈晚打開、攤得愈開、權重也愈大。M4 的包裹兩三年內就過期，M5.8
# 的要十幾年後才完全展開。把一個地區幾十年目錄裡**所有**地震的
# 包裹疊起來，就是 EEPAS 的預報地圖——真實的模型只是把這張圖從
# 「時間 × 一維位置」換成「時間 × 經度 × 緯度 × 規模」。

# %% [markdown]
# ## 16.8 時間完整度：$p(T,L,m)$ 與凸組合補償
#
# $\Delta(m)$ 補的是**規模**方向的缺口。還有一個更麻煩的缺口在
# **時間**方向：M7 的前兆時間中位數是十幾年，而多數地區的可用
# 目錄只有二、三十年——**目錄開始之前的前兆，永遠不會被看到**。
#
# ### 兩個時間尺度
#
# 定義兩個量：**time-lag** $T$ 是目錄結束到預報時刻的落後，
# **catalogue span** $L$ 是目錄涵蓋的長度。一個前兆若要被用上，
# 它到目標地震的等待時間必須落在 $[T,\ T+L]$ 之內——比 $T$ 短的
# 前兆還沒進目錄，比 $T+L$ 長的前兆發生在目錄開始之前。
#
# 於是一個規模 $\nu$ 的前兆能被觀測到的機率，就是它的對數常態
# 時間核落在那個窗內的質量：
#
# $$\Phi\!\left(\frac{\log_{10}(T+L)-a_T-b_T\nu}{\sigma_T}\right)
#  -\Phi\!\left(\frac{\log_{10}T-a_T-b_T\nu}{\sigma_T}\right).$$
#
# 把它對所有可能的前兆規模加權平均（權重是「有多少這種前兆」
# $10^{-b\nu}$ 乘上「它對目標規模 $m$ 貢獻多少」$\eta(\nu)g(m\mid\nu)$）：
#
# $$c(T,L,m)=\int_{m_0}^{m_u}
#   \left[\Phi\!\left(\frac{\log_{10}(T+L)-a_T-b_T\nu}{\sigma_T}\right)
#        -\Phi\!\left(\frac{\log_{10}T-a_T-b_T\nu}{\sigma_T}\right)\right]
#   \eta(\nu)\,g(m\mid\nu)\,10^{-b\nu}\,\mathrm{d}\nu$$
#
# $$p(T,L,m)=\frac{c(T,L,m)}{c(0,\infty,m)}$$
#
# 分母取 $T=0$、$L=\infty$，方括號等於 1，代表「目錄無限長且即時
# 更新」的理想情形。所以 **$p$ 就是完整度：實際拿得到的前兆貢獻
# 佔理論總量的比例**，與 $\Delta(m)$ 完全平行，只是換了一個維度。
#
# ### 補償：兩個端點的凸組合
#
# 知道缺了多少，還要決定怎麼補。有兩種極端立場：
#
# **模型 A（放大背景項）**：缺掉的前兆資訊已經無從得知，把缺掉的
# 那部分率密度攤回 PPE 的長期平滑場。這是**保守**的立場——不敢
# 宣稱看不見的前兆長在哪裡。
#
# **模型 B（放大時變項）**：看得見的前兆是全體的隨機樣本，把時變
# 項整個除以 $p$ 放大回去。這是**進取**的立場——相信空間格局是
# 對的，只是強度被削弱了。
#
# 真相在兩者之間，於是取凸組合：
#
# $$\lambda_C=\phi\,\lambda_A+(1-\phi)\,\lambda_B,\qquad 0\le\phi\le1$$
#
# （這裡的 $\phi$ 是兩個端點模型的混合權重，與第 14 章隨機除叢的
# 背景機率 $\phi_j$ 無關；凸組合的一般性質屬於第 19 章。）$\phi$
# 由最大概似定，而**它本身就是一個可讀的診斷量**：$\phi$ 接近 0
# 表示「看得見的前兆確實有代表性」，接近 1 表示「缺掉的部分不能
# 靠外推補」。
#
# 把 $p$ 對 lag 畫出來：

# %% tags=["hide-input"]
def completeness(T_yr, L_yr, m_target, p=P, b=1.0, m0=2.45, mu_max=8.05,
                 n_nu=1600):
    """p(T, L, m)：前兆時間完整度（T、L 以年為單位）。"""
    nu_ = np.linspace(m0, mu_max, n_nu)
    mu_T = p["a_T"] + p["b_T"] * nu_
    wgt = (np.exp(-b * LN10 * (p["b_M"] - 1) * nu_)
           * norm.pdf(m_target, p["a_M"] + p["b_M"] * nu_, p["sig_M"])
           * 10 ** (-b * nu_))
    T_d = np.atleast_1d(T_yr) * DAY_YR
    hi = norm.cdf((np.log10(T_d[:, None] + L_yr * DAY_YR) - mu_T) / p["sig_T"])
    lo = norm.cdf((np.log10(T_d[:, None]) - mu_T) / p["sig_T"])
    return (np.trapezoid((hi - lo) * wgt, nu_, axis=1)
            / np.trapezoid(wgt, nu_))


L_CAT = 25.0                                    # 目錄長度（年）
T_yr = np.logspace(-1.0, 1.6, 140)              # lag 0.1–40 年
fig = go.Figure()
for k, m_t in enumerate([5.0, 6.0, 7.0]):
    pv = completeness(T_yr, L_CAT, m_t)
    fig.add_trace(go.Scatter(x=T_yr, y=pv, mode="lines",
                             name=f"目標 M{m_t:.0f}（lag 15 年時 "
                                  f"{completeness(15.0, L_CAT, m_t)[0]:.1%}）",
                             line=dict(color=PALETTE[k], width=2.5)))
fig.add_vline(x=15.0, line_dash="dot", line_color="#666666")
fig.update_xaxes(type="log", title_text="time-lag T（年，對數軸）")
fig.update_yaxes(title_text="前兆完整度 p(T, L, m)", range=[0, 1.05])
apply_layout(fig, title=f"目錄落後愈久，看得見的前兆愈少"
                        f"（目錄長度 L = {L_CAT:.0f} 年）",
             hovermode="x", height=430)
fig

# %% [markdown]
# 三條曲線都從接近 1 掉到接近 0，但**掉的位置隨目標規模右移**：
# 目標愈大、前兆時間愈長，能容忍的落後也愈久。M5 在 lag 十幾年時
# 只剩百分之幾——這正是 LEEPAS 論文報告的量級。
#
# 然後是整個 EEPAS 文獻裡最反直覺的一個結果：**lag 拉到 15 年、
# M5 的前兆完整度只剩約 5% 時，LEEPAS 相對 PPE 的資訊增益仍達
# 0.5，甚至高於 lag = 1 年時的值。**
#
# 怎麼可能？兩個機制。**其一**，lag 大時，能用上的前兆只剩「等待
# 時間特別長」的那一群，而那一群恰恰對應**最大的目標地震**——
# 樣本被篩選了，但篩得對。**其二**，補償把權重放回正確的量級，
# 於是剩下的少數前兆不會因為「數量少」而被模型當成「訊號弱」。
#
# 這件事的實務意義很大：**做對補償，預報視窗可以從數個月拉到
# 數十年而不太損失資訊**。紐西蘭把 EEPAS 用在基督城重建的 50 年
# 危害度模型、Kaikōura 之後的 100 年模型，靠的就是這一族方法。

# %% [markdown]
# ## 16.9 二十年成績單
#
# ### 增益隨目標規模上升
#
# 日本本土的正式 CSEP 測試把目標切成三個規模級距分別擬合，得到
# EEPAS 相對 PPE 的每地震資訊增益 $I_{\rm EEPAS}-I_{\rm PPE}$ 為
# 0.24、0.42、1.02。取指數就是機率增益：

# %% tags=["hide-input"]
BRACKET = [("3.95 < M < 4.45", 1040, 0.24),
           ("4.45 < M < 4.95", 396, 0.42),
           ("4.95 < M < 9.05", 148, 1.02)]
gains = [np.exp(d) for _, _, d in BRACKET]

fig = go.Figure(go.Bar(
    x=[b[0] for b in BRACKET], y=gains, marker_color=PALETTE[:3],
    text=[f"{g:.2f}" for g in gains], textposition="outside",
    hovertext=[f"N = {b[1]}，ΔI = {b[2]:.2f}" for b in BRACKET]))
fig.add_hline(y=1.0, line_dash="dash", line_color="#888888",
              annotation_text="與基準持平")
for lab, val, col in [("跨級距移植（低→高）", 0.74, PALETTE[4]),
                      ("跨級距移植（高→低）", 0.61, PALETTE[5])]:
    fig.add_hline(y=val, line_dash="dot", line_color=col,
                  annotation_text=f"{lab} {val:.2f}",
                  annotation_position="bottom right")
fig.update_yaxes(title_text="機率增益 exp(I_EEPAS − I_PPE)", range=[0, 3.2])
apply_layout(fig, title=f"日本本土 CSEP：增益隨目標規模單調上升"
                        f"（{gains[0]:.2f} → {gains[1]:.2f} → {gains[2]:.2f}）",
             xaxis_title="目標規模級距", hovermode="closest", height=430)
fig

# %% [markdown]
# 增益從 1.27 爬到 2.77，**完全符合模型理念**——目標愈大、Ψ 訊號
# 愈清楚。2.77 已屬 EEPAS 歷來應用的高端。
#
# 但同一份研究也給了兩條煞車。**第一條**是圖上那兩條點線：把低
# 規模級距的參數套到高規模級距，增益掉到 0.74；反過來只有 0.61
# ——**比基準模型還糟**。參數不可跨規模級距移植，這是很硬的實證
# 結論。**第二條**是離散度：個別地震的增益從 **0.05 到 20** 都有。
# 平均增益 2.77 是一個集體統計量，**它不承諾任何單一地震**。
#
# ### 短期看觸發、中長期看前兆
#
# 義大利的擬前瞻實驗（HORUS 目錄，學習期 1990–2011、測試期
# 2012–2021、$m_T=5.0$）把 ETAS 與 EEPAS 放上同一張表，按預報
# 視窗長度排開：**三個月內 ETAS 勝；六個月到十年 EEPAS 勝**，
# 而 EEPAS 的優勢**只在五年與十年期達到統計顯著**（主震＋餘震
# 目標集；只取主震集時任何區間都不顯著）。這條時間尺度的分界線，
# 就是 16.2 節那張核對照圖上的交叉點。
#
# 更深刻的是 Janus 混合實驗
# $\lambda=(1-\pi_{\rm ETAS})\lambda_{\rm EEPAS}
# +\pi_{\rm ETAS}\lambda_{\rm ETAS}$。短時間尺度上 ETAS 遠勝，
# 這不意外；意外的是**即使把預報視窗縮到零，混合模型仍勝過純
# ETAS，IGPE 約 0.1**。零視窗意味著沒有任何「等待時間」可言，
# 純粹比較同一瞬間的率密度場——而 EEPAS 仍然貢獻了資訊。結論很
# 強：**觸發串級與前兆尺度增加是兩種大致獨立的可預報性來源**，
# 正確的做法不是二選一，而是按時間尺度加權混合（第 19 章）。
# 紐西蘭約六個月、加州約兩年之後，EEPAS 成為混合模型的主成分。
#
# ### 參數不可跨區域移植，而且有物理原因
#
# 2009 年 EAS 那篇同時擬合了加州與日本關東，得到一組經典數字：
# 加州 $a_T\approx2.05$–$2.11$、$\sigma_A\approx0.72$–$0.85$；
# 關東 $a_T\approx1.18$–$1.31$、$\sigma_A\approx1.76$–$2.25$。
# 換算下來，**加州的平均前兆時間約為關東的 6 倍，而前兆面積
# （正比於 $\sigma_A^2$）約為關東的 1/6**。而關東 M>4 的地震率
# 也恰好約為加州的 6 倍。
#
# 這個 6 倍不是巧合。用物理式地震模擬器系統性地降低斷層滑移
# 速率，前兆時間會成反比拉長——**前兆時間與構造加載速率成反比**。
# 穩定大陸內部（例如澳洲）因此前兆時間長到超出任何目錄的長度，
# 資料點明顯偏離迴歸線。這同時解釋了為什麼「全球 EEPAS」到今天
# 仍是一個未解的挑戰。
#
# 順帶說明：那個 6 倍與 1/6 同時出現，也是第 15 章時空取捨的另一
# 種面貌——$\sigma_A^2$ 與 $10^{a_T}$ 的乘積近似守恆。
#
# ### 最強的一次獨立驗證
#
# 最後是 2004 年那次移植。模型在紐西蘭 1965–2000 的資料上擬合，
# 然後**參數一個不動**（只重估 b 值與 PPE 的三個基線參數）搬到
# 加州 1975–2001 做獨立測試：EEPAS 預期 37.9 個目標地震（實際
# 40 個），PPE 只預期 29.5 個，$\Delta\ln L=34.8$，**累積概似比
# 到 2001 年達到 $10^{15}$**。
#
# 同一篇還做了一組乾淨的隨機化檢驗：把目標地震的時間、位置、
# 規模**分別**隨機打亂各 1000 次，EEPAS 的 $\ln L$ 中位數都低於
# 真實目錄。**三個維度各自都貢獻了預報能力**，不是只有空間平滑
# 在起作用。
#
# 這兩件事合起來，是「EEPAS 不是過擬合」最強的證據——也是任何
# 宣稱有預報能力的模型都應該提供的兩樣東西：**一次真正獨立的
# 移植測試，以及一組把訊號打散的對照實驗**。

# %% [markdown]
# ## 16.10 參數與典型值
#
# 四組已發表的擬合並列。先看規模核與時間核：
#
# | 參數 | 紐西蘭 AMC | 南加州 | 關東 | 義大利 |
# |---|---|---|---|---|
# | $a_M$ | 1.00 | 1.00 | — | 1.23 |
# | $b_M$ | 1.0（固定） | 1.0（固定） | — | 1.0（固定） |
# | $\sigma_M$ | 0.20 | 0.58 | — | 0.24 |
# | $a_T$ | 1.97 | 1.49 | 1.18–1.31 | 2.71 |
# | $b_T$ | 0.35 | 0.48 | — | 0.32 |
# | $\sigma_T$ | 0.20 | 0.81 | — | 0.15 |
#
# 再看空間核、混合權重與設定：
#
# | 參數 | 紐西蘭 AMC | 南加州 | 關東 | 義大利 |
# |---|---|---|---|---|
# | $b_A$ | 0.59 | 0.61 | — | 0.51 |
# | $\sigma_A$ | 0.51 | 0.66 | 1.76–2.25 | 1.00 |
# | $\mu_E$ | 0.36 | 0（擬合值） | — | 0.16 |
# | $b$ | — | 0.96 | — | 1.084 |
# | $m_0$／$m_T$ | —／— | 2.45／4.95 | 2.45／4.75 | 2.45／4.95 |
#
# 出處：紐西蘭 AMC 目錄 2006–2020 的 EEPAS_1F（Rhoades et al.
# 2022，$N=147$；PPE 的 $a=0.55$、$d=5.26$ km）；南加州 1981–2004
# （Rhoades 2007）；關東（Rhoades 2009，該文只報告尺度參數）；
# 義大利 1990–2011 學習期、主震＋餘震目標集（Biondini et al. 2023；
# PPE 的 $a=0.62$、$d\approx30$ km）。破折號表示原文未列。
#
# 讀這兩張表要帶三個提醒。
#
# **一、$\sigma_A$ 不能跨研究直接比。** 它總是與 $10^{b_A m_i}$
# 綁在一起出現，$b_A$ 不同、目錄的規模分布不同，$\sigma_A$ 的
# 數值就不可比——這與第 14 章 $D$ 和 $\gamma$ 只透過乘積被辨識是
# 同一類問題。
#
# **二、南加州那一欄的異常值是有故事的。** $\sigma_M=0.58$、
# $\sigma_T=0.81$、$b_T=0.48$ 都異常大，$a_M$、$\sigma_A$ 異常小。
# 原因在 16.5 節講過：$m_T=4.95$ 遠低於該目錄的最大地震，目標集
# 裡混進大量餘震（1992 年 Landers 之後半年內 50 km 內就有 19 個
# M≥5），最佳化只好把核撐寬去裝它們。**參數表裡的離群值往往不是
# 區域差異，而是設定問題的病徵。**
#
# **三、$\mu_E=0$ 不是壞消息。** 它的意思是「幾乎每個目標地震都有
# 可辨識尺度的前兆活動」。紐西蘭 1965–2000 與南加州都擬出 0，
# 紐西蘭 AMC 目錄卻是 0.36、義大利是 0.16——**這是區域性質，不是
# 普世常數**。
#
# 最後補三條操作規範，都是踩過坑才寫下來的：**delay 通常設 50 天**
# （不設的話參數會被短期叢集主導）；**warm-up 期須對 $M>m_T$
# 完整、擬合期須對 $M>m_0$ 完整**；**lead time 越短，擬合出的
# $a_T$ 越小**——這是系統性偏差，不是區域差異，FLEEPAS 就是為了
# 處理它而生的。

# %% [markdown]
# ## 16.11 常見誤解與陷阱
#
# **一、「EEPAS 能預測下一個大地震。」** 不能。它輸出的是率密度，
# 機率增益 2–3 倍意味著把年機率從 1% 提高到 2–3%，仍然是低機率
# 事件；而且單一地震的增益分布極寬（日本：0.05 到 20）。原作者
# 自己用的詞是「弱預報模型」。
#
# **二、「機率增益大就是好模型。」** 增益永遠是相對於某個基準。
# 換基準、換規模級距、換時間視窗，數字就變。同一個 EEPAS，跨
# 規模級距套用參數後增益掉到 0.61——比基準還糟。**報告增益而不
# 報告基準與目標集定義，等於沒有報告。**
#
# **三、「EEPAS 是 ETAS 的長期版。」** 不是。因果敘事相反（觸發
# vs 醞釀的證據）、預報方向相反（更小的後續事件 vs 更大的後續
# 事件）、核函數族不同（冪次重尾 vs 對數常態峰）、規模核的形狀
# 不同（單調下降 vs 峰值高於 $m_i$）。EEPAS 用到 ETAS 的地方只有
# 算 $w_i$。
#
# **四、「$\mu_E=0$ 表示模型有問題。」** 相反，那表示時變項幾乎
# 解釋了全部——「幾乎每個目標地震都有可辨識尺度的前兆活動」。
# 該擔心的是 $\mu_E$ 接近 1：那才是時變項什麼都沒找到。
#
# **五、「通過 CSEP N-test 表示模型抓對了地震率。」** 不一定。
# 義大利實驗中測試期的目標地震率（約 2.7/年）是學習期（約 1.2/年）
# 的兩倍多，**所有模型都低估了實際個數**，是靠信賴區間較寬的
# 負二項 N-test 才通過。**通過檢定與校準良好是兩件事**（第 17、
# 18 章）。
#
# **六、「公式照抄論文就好。」** 危險。2004 年原文的時間核漏了
# $\ln 10$；$\eta$ 的退化式在不同篇裡差一個 $\beta$；Biondini
# et al.（2023）明白指出各篇公式互有出入而在附錄重新推導。
# **實作前務必核對版本，並用「密度積分是否為 1」自我檢查。**
#
# **七、「模型擬合得好就代表理念正確。」** 2007 年南加州與 2011
# 年日本都示範了反例：目標規模一低，最佳化就會把參數推去遷就
# 餘震，讓模型「偏離了它所依據的理論基礎」。
#
# **八、「$\Delta(m)$ 與 $p(T,L,m)$ 是技術細節。」** 它們是本章
# 最重要的方法學進展。沒有補償，短目錄會系統性低估前兆貢獻；
# 有了補償，同一份目錄能支撐的預報視窗從數月變成數十年。

# %% [markdown]
# ## 16.12 研究前沿與未解問題
#
# 2022 年的二十年回顧列出七項未解挑戰，每一項都對應本章的某個
# 缺口。
#
# **一、Ψ 的物理機制仍然未知。** 這是最根本的一項。Ψ 可以在只放
# 入摩擦定律與斷層交互作用的物理式模擬器產生的合成目錄裡被辨識
# 出來，EEPAS 在那些合成目錄上也一樣有效——**它是某種普遍的力學
# 或幾何後果，但沒有人能說出是哪一種**。這與 Omori 律被確立時的
# 處境一模一樣：先有可用的經驗律，機制留給後人。
#
# **二、把長期地震率納入模型。** ERDEEP 試過令 $a_T$、$\sigma_A$
# 隨長期地震率變化，在正式 CSEP 測試中輸給標準版，假設被放棄。
# 後續釐清了：地震率**確實**控制 $a_T$（前兆時間與滑移速率成
# 反比），時空取捨也**確實**存在，**但兩者無關**。怎麼把前者納入
# 而不引進後者，仍未解決。
#
# **三、三維版本。** 目前的空間核是平面的，沒有深度。台灣、日本、
# 紐西蘭這種隱沒帶環境，深度顯然攜帶資訊。
#
# **四、target-oriented 的缺失前兆補償。** $\Delta$ 與 $p$ 都是
# 對「平均目標」補償的；理想上補償應該隨目標的規模、位置逐格
# 調整。
#
# **五、隨時空變化的目錄完整度。** 本章的 $m_0$ 是一個常數，而
# 第 11 章證明過 $M_c$ 是隨時間、空間、尤其是主震後時刻變動的
# 場。把 $M_c(t,x,y)$ 接進 $\Delta$ 與 $\eta$ 是一項未完成的工程。
#
# **六、時空取捨的形式化。** EEPAS 的 $f$ 與 $h$ 是**相乘**的，
# 也就是假設時間與位置獨立；第 15 章從兩條獨立路徑都證明它們
# 負相關。目前的權宜補救是沿取捨線取三組參數混合（零新增自由
# 參數，獨立測試期勝出），作者自承「如何最佳地把取捨納入 EEPAS
# 仍未解決」。這是最適合當作課堂 open problem 的一個缺口。
#
# **七、全球 EEPAS。** 前兆時間與加載速率成反比、參數不可跨區域
# 移植，兩件事合起來意味著單一套全球參數不可能work——除非把
# 加載速率當成一個顯式的協變量寫進模型。
#
# 對台灣，本章留下一個明確的座標：短期（ETAS）已有本土參數與
# 作業化實測（第 14 章），長期（PSHA）也有國家級模型（第 21 章），
# **中期這一段目前是空白**——而 EEPAS 的台灣在地化工作正在進行
# 中。要填上這一格，本章列出的門檻都得先逐一檢查：目錄需完整到
# 目標規模以下約兩個單位（16.4 節的 $\Delta$）、b 值與餘震比例
# 可能隨規模變化（可能需要分級距擬合）、複雜構造可能需要分區。

# %% [markdown]
# ## 16.13 附錄：本章推導細節
#
# ### A. 對數常態的三個特徵時間
#
# 承 16.2 節，$U=\log_{10}\tau\sim N(\mu_T,\sigma_T^2)$。
#
# **中位數**：$\log_{10}$ 是單調函數，中位數可以直接搬過去，
# $\mathrm{med}(\tau)=10^{\mu_T}$。
#
# **平均**：寫 $\tau=e^{U\ln 10}$，用常態的動差生成函數
# $\mathbb{E}[e^{sU}]=e^{s\mu_T+s^2\sigma_T^2/2}$ 取 $s=\ln 10$：
#
# $$\mathbb{E}[\tau]=10^{\mu_T}\exp\!\left(\frac{\sigma_T^2\ln^2 10}{2}\right).$$
#
# **眾數**：對 $\ln f(\tau)=-\ln\tau-\frac{(\log_{10}\tau-\mu_T)^2}
# {2\sigma_T^2}+\text{const}$ 微分並令其為零。以
# $u=\log_{10}\tau$ 表示，$\mathrm{d}u/\mathrm{d}\tau=1/(\tau\ln10)$：
#
# $$-\frac{1}{\tau}-\frac{u-\mu_T}{\sigma_T^2}\cdot\frac{1}{\tau\ln 10}=0
#   \quad\Longrightarrow\quad
#   u=\mu_T-\sigma_T^2\ln 10 ,$$
#
# 故 $\mathrm{mode}(\tau)=10^{\mu_T-\sigma_T^2\ln 10}$。
#
# 三者的順序恆為 眾數 < 中位數 < 平均。以紐西蘭 AMC 的
# $\sigma_T=0.20$ 為例，平均比中位數大 $\exp(0.04\times2.651/2)
# \approx1.055$ 倍，差別很小；但 2007 年南加州的 $\sigma_T=0.81$
# 讓這個倍數變成 $\exp(0.656\times2.651/2)\approx2.4$ 倍。
# **論文報告的「前兆時間」幾乎都是中位數 $10^{a_T+b_T m}$**，
# 引用時要確認。
#
# ### B. $\eta$ 推導裡用到的高斯積分
#
# 承 16.3 節。要證
# $\int_{-\infty}^{\infty}e^{-\beta u}\phi_\sigma(m-u)\,\mathrm{d}u
# =e^{-\beta m}e^{\sigma^2\beta^2/2}$。令 $v=m-u$，
# $\mathrm{d}u=-\mathrm{d}v$，積分上下限對調：
#
# $$\begin{aligned}
# \int_{-\infty}^{\infty}e^{-\beta u}\phi_\sigma(m-u)\,\mathrm{d}u
#   &= \int_{-\infty}^{\infty}e^{-\beta(m-v)}\phi_\sigma(v)\,\mathrm{d}v \\
#   &= e^{-\beta m}\int_{-\infty}^{\infty}e^{\beta v}\phi_\sigma(v)\,\mathrm{d}v
#    = e^{-\beta m}\,\mathbb{E}\!\left[e^{\beta V}\right],
# \end{aligned}$$
#
# $V\sim N(0,\sigma^2)$。動差生成函數 $\mathbb{E}[e^{sV}]
# =e^{s^2\sigma^2/2}$ 給出 $e^{\sigma^2\beta^2/2}$。
#
# **這個因子的物理讀法值得記**：規模核的寬度 $\sigma$ 讓一部分
# 機率質量被搬到高規模端，而高規模端在 GR 權重下「稀有」，所以
# 加權平均反而被放大。$\eta$ 裡的 $-\beta\sigma_M^2\beta/2$ 就是
# 把這份虛增扣回去。$\sigma_M$ 愈大，扣得愈多。
#
# ### C. $\Delta(m)$ 的配方法
#
# 承 16.4 節。展開指數：
#
# $$\begin{aligned}
# -\beta u-\frac{(u-m)^2}{2\sigma_M^2}
#   &= -\frac{1}{2\sigma_M^2}\left[2\sigma_M^2\beta u+u^2-2mu+m^2\right] \\
#   &= -\frac{1}{2\sigma_M^2}
#      \left[u^2-2\bigl(m-\sigma_M^2\beta\bigr)u+m^2\right] \\
#   &= -\frac{1}{2\sigma_M^2}
#      \left[\bigl(u-(m-\sigma_M^2\beta)\bigr)^2
#        -\bigl(m-\sigma_M^2\beta\bigr)^2+m^2\right].
# \end{aligned}$$
#
# 最後兩項與 $u$ 無關，在 $\Delta$ 的分子分母中相同而消去。剩下
# 的部分是 $N(m-\sigma_M^2\beta,\ \sigma_M^2)$ 的核，於是
#
# $$\Delta(m)=P\bigl(U\ge a_M+b_M m_0\bigr)
#   =\Phi\!\left(\frac{m-a_M-b_M m_0-\sigma_M^2\beta}{\sigma_M}\right),$$
#
# 最後一步用了常態的對稱性 $1-\Phi(z)=\Phi(-z)$。
#
# 兩個極限檢核。$m\to\infty$ 時 $\Delta\to1$（大目標的前兆都在
# 門檻之上，沒有缺口）；$m_0\to-\infty$ 時 $\Delta\to1$（目錄
# 無限完整）。兩者都是應有的行為。
#
# ### D. 兩個高斯的卷積
#
# 承 16.7 節。逐軸計算：
#
# $$\begin{aligned}
# \int\phi_{\sigma_2}(x_a-x)\,\phi_{\sigma_1}(x-x_i)\,\mathrm{d}x
# \end{aligned}$$
#
# 指數部分為
# $-\frac{(x_a-x)^2}{2\sigma_2^2}-\frac{(x-x_i)^2}{2\sigma_1^2}$。
# 令 $\sigma_*^{-2}=\sigma_1^{-2}+\sigma_2^{-2}$，即
# $\sigma_*^2=\sigma_1^2\sigma_2^2/(\sigma_1^2+\sigma_2^2)$，
# 對 $x$ 配方：
#
# $$-\frac{1}{2\sigma_*^2}
#   \left(x-\sigma_*^2\left[\frac{x_a}{\sigma_2^2}
#     +\frac{x_i}{\sigma_1^2}\right]\right)^2
#   -\frac{(x_a-x_i)^2}{2(\sigma_1^2+\sigma_2^2)} .$$
#
# 第一項對 $x$ 積分給出 $\sigma_*\sqrt{2\pi}$，與前面的
# $1/(2\pi\sigma_1\sigma_2)$ 相乘後恰為
# $1/\sqrt{2\pi(\sigma_1^2+\sigma_2^2)}$（因為
# $\sigma_*/(\sigma_1\sigma_2)=1/\sqrt{\sigma_1^2+\sigma_2^2}$）；
# 第二項與 $x$ 無關，直接留下。合起來
#
# $$\int\phi_{\sigma_2}(x_a-x)\,\phi_{\sigma_1}(x-x_i)\,\mathrm{d}x
#   =\phi_{\sqrt{\sigma_1^2+\sigma_2^2}}(x_a-x_i).$$
#
# 兩軸獨立相乘就是 16.7 節的結果。**變異數相加**這件事有一個更
# 快的看法：卷積對應獨立隨機變數相加，$Z=X+E$，而獨立變數的
# 變異數本來就相加。配方法只是把它算了一遍。
#
# ### E. 從 $a_A$ 到 $\sigma_A$ 的幾何常數
#
# 承 16.2 節。Ψ 的 $A_P$ 是「涵蓋前兆、主震與餘震的最小矩形
# 面積」，而 EEPAS 的 $h$ 是一個圓對稱二維常態。要把兩者接起來，
# 得先選一個約定：矩形大致對應每軸 $\pm k\sigma_i$ 的方框。於是
#
# $$A_P\approx(2k\sigma_i)^2=4k^2\sigma_A^2\,10^{b_A m_i}
#   \quad\text{而}\quad A_P=10^{a_A+b_A m_i},$$
#
# 兩式的 $10^{b_A m_i}$ 對消，得 $\sigma_A=10^{a_A/2}/(2k)$。
# $k=1.5$ 給比例常數 0.33，$k=2$ 給 0.25——這就是「理論上約
# 0.25–0.33」的來源。2004 年實際擬合出 0.37，對應
# $k\approx1.35$：**EEPAS 攤開的機率雲比最小涵蓋矩形稍寬**。
#
# 這個轉換也解釋了 $\sigma_A$ 為什麼不可跨研究直接比：它吸收了
# $a_A$、吸收了矩形與高斯之間的幾何約定，還與 $b_A$ 綁在一起。
# **要比較兩個研究的空間尺度，該比的是某個共同 $m_i$ 下的
# $\sigma_A\,10^{b_A m_i/2}$，不是 $\sigma_A$ 本身。**
#
# ---
#
# 回頭看這一章走過的路：從「放棄辨識」這一步側身出發，用一條
# 平凡的線性迴歸條件分布長出三個核，用一個守恆要求長出 $\eta$，
# 用兩個「看不見的東西」長出 $\Delta$ 與 $p$，再用一個舊朋友
# （第 14 章的背景機率）長出 $w_i$。整台機器沒有一個零件是為了
# 好看而加的——**每一個都是被一個具體的缺陷逼出來的**。
#
# 壓成一句話：**EEPAS 的價值不在於它預報得多準，而在於它示範了
# 如何把一個爭議中的前兆現象，轉化成可自動化、可檢驗、可否證的
# 預報模型。** 它的作者把 28% 找不到 Ψ、單一辨識不唯一、機制
# 未知、參數會為了遷就餘震而說謊，全部寫進了論文——這在前兆
# 研究的歷史上並不常見。
#
# 但這也把一個問題推到了眼前。這一章一路在用「資訊增益 0.42」
# 「機率增益 2.77」「通過負二項 N-test」這些說法，卻沒有定義過
# 其中任何一個。$\ln L$ 差 34.8 到底算大還是小？「所有模型都低估
# 了地震數，但靠較寬的信賴區間通過檢定」是什麼意思？兩個模型的
# 增益差 0.1，可以說誰比較好嗎？{doc}`第 17 章 <17_testing_consistency>`
# 就從最基本的問題開始：**給定一份預報與一份觀測，怎麼判斷它們
# 相不相容？**
