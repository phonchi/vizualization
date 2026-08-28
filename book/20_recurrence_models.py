# %% [markdown]
# # 20. 時間相依複發模型：Weibull、BPT 與應力釋放
#
# {doc}`第 19 章 <19_ensembles>`把一堆模型組合起來，但那些模型有一個
# 共同的時間尺度：天到年。ETAS 看的是「上一個地震剛剛發生」，EEPAS
# 看的是「這幾年的活動變密了」，PPE 看的是「過去五十年哪裡常震」。
# 沒有一個模型在問這個問題：**這條斷層上一次破裂是 1906 年，現在
# 是不是快到了**？
#
# 這個問題屬於另一個世界。它的時間單位是**世紀**，它的資料不是三十五萬
# 筆的儀器目錄而是三到五次的古地震記錄，它的使用者不是防災中心而是
# 建築規範的制定者。它也是台灣官方危害度模型 TEM PSHA2020 真正在算的
# 東西——那份模型對每一條孕震構造都掛了一個 **Brownian passage time**
# 分布，用來描述「斷層的記憶」。
#
# {doc}`第 10 章 <10_point_process>`的 10.3 節已經把地圖畫好了。那一節把
# 條件強度依「記憶的方向」分成三類——無記憶的 Poisson、正記憶的
# self-exciting、負記憶的 self-correcting——並且在最後補了第四類：
# **更新過程**（renewal process），$\lambda^*$ 只依賴「距上一個事件多久」，
# 記得上一次、忘記上上次。這一章就是那一節的**量化展開**：把三種記憶
# 結構從定性的曲線形狀，推成可以拿去算機率、可以餵進 PSHA 積分的公式。
#
# 分工上，本章只重推**新東西**：hazard function、BPT 的首達時間、
# 應力釋放模型的參數意義、b 值與 aperiodicity 的橋樑。條件強度
# {eq}`eq:cond-int`、概似 {eq}`eq:pp-loglik`、時間變換 {eq}`eq:time-rescale`
# 全部沿用第 10 章，不重推。
#
# 一句話預告本章的結論：**選複發分布時該看的不是機率密度 $f(t)$，
# 而是危害函數 $h(t)$**。三個候選分布的 $f$ 看起來幾乎一樣，$h$ 卻在
# 長時間端分道揚鑣到物理意義完全相反——一個趨於無窮、一個趨於零、
# 一個趨於有限值。而斷層上的問題，恰好都發生在長時間端。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 20.1 為何長期模型需要「記憶」
#
# 先把虛無假設的性質講清楚，才知道要反對什麼。
#
# Poisson 過程的間隔服從指數分布 $f(t)=\lambda e^{-\lambda t}$，
# 存活函數 $S(t)=P(T>t)=e^{-\lambda t}$。問一個具體問題：**已經等了 $s$
# 年沒有地震，再等 $t$ 年還是沒有地震的機率是多少**？
#
# $$\begin{aligned}
# P(T>s+t \mid T>s)
#   &= \frac{P(T>s+t,\;T>s)}{P(T>s)} \\
#   &= \frac{P(T>s+t)}{P(T>s)}
#      \qquad (\{T>s+t\}\subset\{T>s\}) \\
#   &= \frac{e^{-\lambda(s+t)}}{e^{-\lambda s}}
#    = e^{-\lambda t} \;=\; P(T>t).
# \end{aligned}$$
#
# 右邊完全不含 $s$。**已經等了多久，對未來沒有任何影響**——這就是
# 無記憶性（memorylessness），而指數分布是唯一具有這個性質的連續分布。
#
# 這條性質在短期預報裡無害甚至好用（第 21 章的 PSHA 就整章建立在
# 它上面）。但把它套到單一斷層上，它說的是：**1906 年梅山地震之後，
# 這條斷層在 1907 年再破裂一次的機率，跟它在 2026 年破裂的機率一樣。**
#
# 這與彈性回跳理論（Reid 1910）正面衝突。斷層的物理圖像是：板塊運動
# 以大致固定的速率把彈性應變能存進岩體，能量累積到岩石強度的上限就
# 破裂、把應變放掉、然後從低點重新開始累積。**一次大破裂之後，斷層
# 需要時間重新裝填**。剛破裂完的斷層應該是安全的；等得愈久，愈危險。
# 指數分布不但沒有這個性質，它還說「危險程度恆定」——連方向都不對。
#
# ### 兩種記憶，兩個時間尺度
#
# 這裡要小心一個容易搞混的地方。{doc}`第 13 章 <13_etas_structure>`的
# ETAS 也是有記憶的模型，但它的記憶方向**相反**：一個地震剛發生完，
# ETAS 說接下來的風險**上升**（餘震），而彈性回跳說接下來的風險
# **下降**（應力放掉了）。兩者不是矛盾，是尺度不同：
#
# - **短期（分鐘～月）**：看到的是應力重分布的正效應。主震把週邊的
#   庫倫應力推高，觸發一大群餘震。這是 self-exciting，第 13、14 章。
# - **長期（數十～數百年）**：看到的是應力耗竭。同一段斷層面上的
#   應變能真的被放掉了，要等板塊重新加載。這是 self-correcting，本章。
#
# 10.3 節那句話值得重讀：**地震觸發地震**與**地震消耗應力**，在點過程
# 的語言裡是同一條句型的正反兩讀。本章要做的就是把「反讀」那一半
# 寫成可以計算的東西。

# %% [markdown]
# ## 20.2 更新過程的核心式
#
# 更新過程的設定極簡：斷層上連續兩次特徵地震的**間隔** $T_1,T_2,\dots$
# 是獨立同分布的隨機變數，共同的機率密度是 $f(t)$。整個 fault-based
# 家族的差異，全在於**選哪一個 $f$**。
#
# ### 條件複發機率
#
# 實務上要回答的問題長這樣：距上次破裂已經過了 $T$ 年，**未來 $\Delta T$
# 年內破裂的機率是多少**？（$T$ 是已知的、確定的等待時間；$\Delta T$
# 通常取 30 年或 50 年，因為那是建築物的設計壽命。）
#
# 令 $X$ 為這一輪的複發間隔（隨機）。已知的資訊是 $X>T$，要求的是
# $P(T< X \le T+\Delta T)$ 在這個條件下的值。直接套條件機率定義：
#
# $$\begin{aligned}
# P\bigl(T \le t \le T+\Delta T \,\big|\, t > T\bigr)
#   &= \frac{P\bigl(T < X \le T+\Delta T,\; X > T\bigr)}{P(X>T)} \\
#   &= \frac{P\bigl(T < X \le T+\Delta T\bigr)}{P(X>T)} \\
#   &= \frac{\displaystyle\int_T^{T+\Delta T} f(u)\,\mathrm{d}u}
#           {\displaystyle\int_T^{\infty} f(u)\,\mathrm{d}u} .
# \end{aligned}$$ (eq:renewal-prob)
#
# 第二個等號成立是因為事件 $\{T<X\le T+\Delta T\}$ 本身就蘊含 $X>T$。
# 分母是**存活函數** $S(T)=1-F(T)$。
#
# {eq}`eq:renewal-prob` 就是所有「未來 30 年破裂機率 xx%」這類新聞
# 標題背後的式子。加州聖安地列斯、日本南海海槽、台灣的孕震構造，
# 官方發布的數字全部由它算出，差別只在 $f$ 的選擇與 $f$ 的參數。
#
# 順帶檢查一致性：把指數分布代進去，分子 $=e^{-\lambda T}-
# e^{-\lambda(T+\Delta T)}$、分母 $=e^{-\lambda T}$，相除得
# $1-e^{-\lambda\Delta T}$——不含 $T$，正是 20.1 節的無記憶性。
#
# ### 該看的是 hazard，不是 pdf
#
# {eq}`eq:renewal-prob` 有一個更好用的等價寫法。把 $\Delta T$ 縮到無窮小、
# 除以 $\Delta T$ 取極限，得到**危害函數**（hazard function，又稱瞬時
# 破裂率）：
#
# $$h(t) \;=\; \lim_{\Delta\to 0}\frac{1}{\Delta}\,
#   P\bigl(t < X \le t+\Delta \,\big|\, X>t\bigr)
#   \;=\; \frac{f(t)}{1-F(t)} \;=\; \frac{f(t)}{S(t)}$$ (eq:hazard-def)
#
# 對照 {eq}`eq:cond-int` 會發現這其實是同一個東西：**更新過程的條件強度
# 就是 hazard function**，只是它的自變數不是絕對時間 $t$ 而是「距上次
# 事件的時間」。也就是說，$\lambda^*(t)=h(t-t_{\text{last}})$。這是把
# 本章接回第 10 章框架的那一根釘子——寫得出 $h$，概似
# {eq}`eq:pp-loglik`、殘差 {eq}`eq:time-rescale` 全部自動可用。
#
# $h$ 與 $S$ 之間有一個互推關係。由 $S'(t)=-f(t)$ 得
# $h(t)=-S'(t)/S(t)=-\frac{\mathrm{d}}{\mathrm{d}t}\ln S(t)$，積分後
#
# $$S(t) = \exp\left[-\int_0^t h(u)\,\mathrm{d}u\right],
#   \qquad
#   P\bigl(t\le T+\Delta T \mid t>T\bigr)
#   = 1-\exp\left[-\int_T^{T+\Delta T} h(u)\,\mathrm{d}u\right].$$
#
# 於是 {eq}`eq:renewal-prob` 完全由 $h$ 決定：**知道 $h$ 就知道一切，
# 而且知道的方式比 $f$ 直觀得多**。
#
# 為什麼堅持看 $h$ 不看 $f$？因為 $f(t)$ 回答的是「從破裂那一刻算起，
# 下一次破裂落在 $t$ 附近的機率密度」——但我們**已經知道**它還沒發生。
# 條件在「還沒發生」之後，$f$ 的絕對高度就沒有意義了，有意義的只有
# $f$ 相對於剩餘機率 $S$ 的比例。長時間端 $f$ 與 $S$ 通常都趨於零，
# **兩個趨於零的量相除，結果可以是零、可以是常數、可以是無窮**——
# 下一節的三個分布正好把這三種結局各佔一個。

# %% [markdown]
# ## 20.3 三個分布
#
# 文獻上最常用的三個複發時間分布是 Weibull、對數常態與 BPT。三者都是
# 正值、單峰、右偏，都由兩個參數決定（一個定尺度、一個定形狀）。
#
# ### Weibull
#
# $$f(t\mid k,\theta) = \frac{k}{\theta}\left(\frac{t}{\theta}\right)^{k-1}
#   \exp\left[-\left(\frac{t}{\theta}\right)^{k}\right],\qquad t>0$$
#
# 記號說明：文獻常把形狀參數寫成 $\beta$ 或 $b$，兩個都與本書已占用的
# 符號撞名（$\beta=b\ln 10$、$b$ 是 GR 斜率）。本書一律用 $k$ 表形狀、
# $\theta$ 表尺度。這不是潔癖——**把 Weibull 的形狀參數叫做 $b$，
# 在同一章裡又出現 GR 的 $b$ 值，是教科書等級的地雷**。
#
# Weibull 的存活函數可以直接積出來（這是它受歡迎的主因）：
# $S(t)=\exp[-(t/\theta)^k]$。於是 hazard 是
#
# $$h(t) = \frac{f(t)}{S(t)}
#   = \frac{\frac{k}{\theta}(t/\theta)^{k-1}
#           \exp[-(t/\theta)^{k}]}{\exp[-(t/\theta)^{k}]}
#   = \frac{k}{\theta}\left(\frac{t}{\theta}\right)^{k-1}.$$
#
# 指數項整個消掉，剩下一條純冪次。形狀參數 $k$ 的物理意義因此非常乾淨：
#
# - $k=1$：$h(t)=1/\theta$ 為常數，Weibull 退化成指數分布，**無記憶**，
#   回到 Poisson。
# - $k<1$：$h$ 單調**遞減**——剛發生完最危險，之後愈來愈安全。這是
#   **短期叢集**的簽名（餘震序列就長這樣）。
# - $k>1$：$h$ 單調**遞增**——等愈久愈危險。這是**準週期**行為，
#   fault-based 模型的典型情形。
#
# ### 對數常態
#
# $$f(t\mid m,\sigma) = \frac{1}{\sqrt{2\pi}\,\sigma t}
#   \exp\left[-\frac{(\ln t-m)^{2}}{2\sigma^{2}}\right]$$
#
# 也就是 $\ln T\sim N(m,\sigma^2)$。它是重尾分布，對「超過平均複發期
# 很多」的長間隔給出偏高的機率——這聽起來是保守（安全）的選擇，但
# 20.4 節會看到它其實給出物理上最不合理的 hazard。
#
# 它的存活函數寫得出來但積不成初等函數：$S(t)=\Phi(-z)$，其中
# $z=(\ln t-m)/\sigma$、$\Phi$ 是標準常態的累積分布。於是
#
# $$h(t) = \frac{\varphi(z)}{\sigma t\,\Phi(-z)},
#   \qquad \varphi(z)=\tfrac{1}{\sqrt{2\pi}}e^{-z^{2}/2}.$$
#
# 大 $t$ 時 $z\to\infty$，用 Mills 比 $\Phi(-z)\sim\varphi(z)/z$ 得
# $h(t)\sim z/(\sigma t)=(\ln t-m)/(\sigma^{2}t)$。分母的 $t$ 贏過分子的
# $\ln t$，所以 $h(t)\to 0$。**等愈久，愈不會發生**——記住這個結論，
# 20.4 節要拿它做文章。
#
# ### BPT：從布朗鬆弛振盪器推起
#
# 前兩個分布是「借來的」——Weibull 來自材料疲勞、對數常態來自乘性
# 誤差累積，套到斷層上是類比而非推導。BPT 不一樣：**它是從一個明確的
# 物理模型推出來的**，這是它在 UCERF2／UCERF3／TEM PSHA2020 中被採用
# 的主因。
#
# 物理模型叫 **Brownian Relaxation Oscillator**。斷層上的應力寫成
#
# $$X(t) \;=\; \rho\,t \;+\; \sigma\,W(t)$$
#
# 兩項各有所指：$\rho$ 是**定值加載率**（板塊運動送進來的應力，
# 單位時間增加 $\rho$），$W(t)$ 是**標準布朗運動**、$\sigma$ 是它的
# 強度——它代表所有沒被模型寫出來的干擾：鄰近地震的應力轉移、
# 慢滑移、孔隙壓變化、斷層面幾何的不均勻。事件發生的條件是
# **$X$ 首次觸及破裂門檻 $a$**；一旦觸及，應力歸零、時鐘重設。
#
# 所以複發間隔就是布朗運動加漂移的**首達時間**（first passage time）：
#
# $$\tau = \inf\{t>0:\;X(t)\ge a\}.$$
#
# 這個分布有閉式解，叫**逆高斯分布**（inverse Gaussian）。完整推導
# 放在 20.11 節附錄 A（反射原理 ＋ Girsanov 測度變換），結果是
#
# $$f(\tau) = \frac{a}{\sigma\sqrt{2\pi\tau^{3}}}
#   \exp\left[-\frac{(a-\rho\tau)^{2}}{2\sigma^{2}\tau}\right].$$
#
# 現在換成有物理意義的參數。平均值是 $T_r=a/\rho$（門檻除以加載率，
# 完全符合直覺），變異係數是 $c_v=\sigma/\sqrt{a\rho}$。把 $a=\rho T_r$、
# $\sigma^2=a\rho\,c_v^2=\rho^2T_r c_v^2$ 代回去、約掉 $\rho$：
#
# $$f(t \mid T_r, c_v) = \sqrt{\frac{T_r}{2\pi c_v^{2}t^{3}}}\;
#   \exp\left[-\frac{(t-T_r)^{2}}{2\,T_r\,c_v^{2}\,t}\right]$$ (eq:bpt-pdf)
#
# 這就是 **BPT 分布**（Brownian passage time，Matthews et al. 2002）。
# 兩個參數的意義都很實在：
#
# - $T_r$ 是**平均複發時間**。文獻常寫 $\mu$，與本書的背景率撞名，
#   本書一律寫 $T_r$。
# - $c_v$ 是 **aperiodicity**（非週期度），恰好等於複發間隔的變異係數
#   $\mathrm{SD}/\mathrm{mean}$。文獻常寫 $\alpha$，與 ETAS 的產能指數
#   撞名，本書一律寫 $c_v$。
#
# $c_v$ 讀起來就是<strong>「噪音相對於加載的比例」</strong>：$c_v=\sigma/\sqrt{a\rho}$，
# 干擾愈大、門檻愈低、加載愈慢，斷層就愈不規律。極限行為很乾淨——
# $c_v\to 0$ 時 $\sigma\to0$，應力是一條直線，每 $T_r$ 年準時破裂一次
# （完全週期）；$c_v\to\infty$ 時漂移相對可忽略，破裂時間趨近隨機。
# 實務研究多落在 $0.3\le c_v\le 0.7$。
#
# 先把三個分布放在**同一個平均複發期**下畫出來，看看它們的 $f$ 有多像：

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import optimize, special, stats

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

T_R, CV = 100.0, 0.50          # 平均複發期（年）與 aperiodicity


def make_three(T_r=T_R, cv=CV):
    """在同一 (平均, 變異係數) 下建構 BPT／Weibull／對數常態。"""
    bpt = stats.invgauss(mu=cv**2, scale=T_r / cv**2)          # 見 20.11 附錄 B
    k = optimize.brentq(                                        # 解 Weibull 形狀
        lambda x: np.sqrt(special.gamma(1 + 2 / x)
                          / special.gamma(1 + 1 / x) ** 2 - 1) - cv, 0.2, 50)
    wbl = stats.weibull_min(k, scale=T_r / special.gamma(1 + 1 / k))
    s = np.sqrt(np.log(1 + cv**2))                              # 對數常態
    lgn = stats.lognorm(s, scale=np.exp(np.log(T_r) - s**2 / 2))
    return bpt, wbl, lgn, k, s


bpt, wbl, lgn, K_W, S_LN = make_three()
tg = np.linspace(1.0, 400.0, 1200)

fig = go.Figure()
for d, name, color in [(bpt, f"BPT（T_r={T_R:.0f}, c_v={CV:.2f}）", PALETTE[0]),
                       (wbl, f"Weibull（k={K_W:.2f}）", PALETTE[1]),
                       (lgn, f"對數常態（σ={S_LN:.3f}）", PALETTE[2])]:
    fig.add_trace(go.Scatter(x=tg, y=d.pdf(tg), mode="lines", name=name,
                             line=dict(color=color, width=2.2)))
fig.add_vline(x=T_R, line=dict(color="#888", width=1, dash="dot"))
fig.add_annotation(x=T_R, y=max(bpt.pdf(tg)) * 1.02, text=f"T_r = {T_R:.0f} 年",
                   showarrow=False, yshift=8, font=dict(size=11))
fig.update_xaxes(title_text="距上次破裂的時間 t（年）")
fig.update_yaxes(title_text="機率密度 f(t)（1/年）")
apply_layout(fig, height=430,
             title=(f"圖 1　同一平均複發期與同一 c_v={CV:.2f} 下的三個分布："
                    "pdf 幾乎重合"))
fig

# %% [markdown]
# 三條線幾乎疊在一起。峰值位置差幾年、尾巴粗細略有不同，但如果你手上
# 只有三到五個歷史複發間隔（20.7 節會說明這是常態），**這三條曲線在
# 統計上完全無法區分**。任何以「哪個分布擬合得比較好」為名的模型選擇，
# 在這種樣本數下都是自欺。
#
# 這正是為什麼要看 hazard。

# %% [markdown]
# ## 20.4 hazard 形狀決定一切
#
# 把同樣三個分布換成 $h(t)=f(t)/S(t)$ 再畫一次。這是本章最重要的一張圖。

# %% tags=["hide-input"]
th = np.linspace(2.0, 400.0, 1500)
H_ASYM = 1.0 / (2.0 * T_R * CV**2)                # BPT 的 hazard 漸近值

fig = go.Figure()
for d, name, color in [(bpt, "BPT", PALETTE[0]),
                       (wbl, f"Weibull（k={K_W:.2f}>1）", PALETTE[1]),
                       (lgn, "對數常態", PALETTE[2])]:
    fig.add_trace(go.Scatter(x=th, y=d.pdf(th) / d.sf(th), mode="lines",
                             name=name, line=dict(color=color, width=2.2)))
fig.add_trace(go.Scatter(
    x=th, y=np.full_like(th, H_ASYM), mode="lines",
    name=f"BPT 漸近線 1/(2·T_r·c_v²)={H_ASYM:.4f}",
    line=dict(color=PALETTE[0], width=1.4, dash="dash")))
fig.add_trace(go.Scatter(
    x=th, y=np.full_like(th, 1.0 / T_R), mode="lines",
    name=f"Poisson（1/T_r={1/T_R:.3f}，水平）",
    line=dict(color="#8a8a8a", width=1.4, dash="dot")))
fig.add_vline(x=T_R / 2, line=dict(color="#bbb", width=1))
fig.add_annotation(x=T_R / 2, y=0.055, text="半個複發期", showarrow=False,
                   xshift=42, font=dict(size=11))
fig.update_xaxes(title_text="距上次破裂的時間 t（年）")
fig.update_yaxes(title_text="危害函數 h(t)（1/年）", range=[0, 0.075])
apply_layout(fig, height=460,
             title="圖 2　同樣三個分布的 hazard：長時間端分道揚鑣")
fig

# %% [markdown]
# 同樣三個分布，pdf 幾乎重合，hazard 卻走向三個完全不同的命運。
#
# **Weibull（$k>1$）單調遞增、永不飽和**。$h(t)=\frac{k}{\theta}
# (t/\theta)^{k-1}$ 是純冪次，$t\to\infty$ 時 $h\to\infty$。物理上的
# 讀法是：只要還沒破，危險度就無上限地攀升。對一條沉寂了十倍複發期
# 的斷層，Weibull 會給出荒謬的高機率。
#
# **對數常態的 hazard 先升後降，最後趨於零**。20.3 節推過
# $h(t)\sim(\ln t-m)/(\sigma^2 t)\to 0$。白話：**等愈久愈不會發生**。
# 對一條沉寂太久的斷層，對數常態的回答是「它大概不會再動了」——
# 除非有獨立證據說斷層已停止活動，否則這在物理上說不通。它讓
# 「最危險的斷層」看起來最安全，方向完全錯了。
#
# **BPT 的 hazard 趨於有限值** $1/(2T_rc_v^{2})$。圖上的虛線就是它。
# 完整推導在 20.11 節附錄 C，結論的物理讀法很漂亮：一條沉寂極久的
# 斷層，並不會愈來愈危險到無窮（因為布朗噪音已經把「應力到底累積
# 到哪」的資訊消磨光了），也不會變安全，而是**退化成一個 Poisson
# 過程**——只是它的等效速率 $1/(2T_rc_v^2)$ 通常高於長期平均速率
# $1/T_r$。三種行為裡，只有這一種在物理上站得住腳。
#
# 圖上還有一件事值得注意：$c_v=0.5$ 時 BPT 的漸近速率
# 是 $1/(2\times100\times0.25)=0.02$ 次／年，剛好是 Poisson 速率
# $1/T_r=0.01$ 的**兩倍**。「沉寂很久的斷層退化成 Poisson」不等於
# 「回到虛無假設」——它退化成一個**兩倍危險的** Poisson。
#
# ### 什麼時候選錯分布會出事
#
# 看圖 2 的左半段：$t$ 在半個複發期附近（$t\approx 50$ 年），三條線
# 擠在一起，數值差異小到不影響任何決策。Polidoro et al. (2013) 的
# 觀察就是這件事：**當「距上次事件的時間」約為複發期的一半時，
# 各模型結果差不多；時間拖得愈久，選錯分布的代價愈大**。
#
# 這對實務有直接的意義。台灣多數孕震構造的複發期是數百到數千年，
# 而可靠的古地震記錄往往只往回追幾百年，$T/T_r$ 常常落在 0.1–0.5
# 的區間——正是三個模型都差不多的區間。**在這個區間裡爭論該用哪個
# 分布，是把力氣花錯地方；真正該擔心的是 $T_r$ 本身估得準不準。**
# 反過來，對於已經明顯「逾期」的斷層（$T>T_r$），分布選擇就變成
# 一階效應，而那正是最需要小心的情境。

# %% [markdown]
# ### $c_v$ 對形狀的影響
#
# 固定 $T_r$、只動 aperiodicity，BPT 的 hazard 會從「幾乎週期性的
# 尖峰」變成「幾乎平坦的 Poisson」。這是實務上唯一需要調的旋鈕：

# %% tags=["hide-input"]
tc = np.linspace(2.0, 300.0, 1200)
fig = go.Figure()
for i, cv in enumerate([0.2, 0.5, 0.8]):
    d = stats.invgauss(mu=cv**2, scale=T_R / cv**2)
    asym = 1.0 / (2 * T_R * cv**2)
    fig.add_trace(go.Scatter(x=tc, y=d.pdf(tc) / d.sf(tc), mode="lines",
                             name=f"c_v = {cv:.1f}（漸近 {asym:.3f}）",
                             line=dict(color=PALETTE[i], width=2.2)))
    fig.add_trace(go.Scatter(x=tc, y=np.full_like(tc, asym), mode="lines",
                             showlegend=False,
                             line=dict(color=PALETTE[i], width=1, dash="dash")))
fig.add_trace(go.Scatter(x=tc, y=np.full_like(tc, 1 / T_R), mode="lines",
                         name=f"Poisson 1/T_r = {1/T_R:.2f}",
                         line=dict(color="#8a8a8a", width=1.4, dash="dot")))
fig.add_vline(x=T_R, line=dict(color="#bbb", width=1))
fig.update_xaxes(title_text="距上次破裂的時間 t（年）")
fig.update_yaxes(title_text="h(t)（1/年）", range=[0, 0.16])
apply_layout(fig, height=440,
             title=f"圖 3　aperiodicity 決定 BPT 的記憶強度（T_r={T_R:.0f} 年）")
fig

# %% [markdown]
# $c_v=0.2$ 的曲線在前 50 年幾乎貼著零（**斷層剛破裂完，非常安全**），
# 然後在 $T_r$ 附近急遽衝高，最後停在 $1/(2T_rc_v^2)=0.125$——是 Poisson
# 速率的 12.5 倍。這是強記憶：時鐘走得準，快到期就非常危險。
#
# $c_v=0.8$ 幾乎是一條平線，早期就已經接近漸近值。這是弱記憶：噪音
# 大到「上次什麼時候破」幾乎沒有資訊量。**$c_v$ 愈大，模型愈接近
# Poisson——但這句話只在早期成立**，長期端 $c_v$ 大反而讓漸近速率
# 更低（$0.0078$ 對比 $c_v=0.2$ 的 $0.125$）。
#
# 這也解釋了為什麼危害度模型對 $c_v$ 這麼敏感、為什麼它總是被放進
# logic tree 的分支裡：從 0.2 換到 0.8，同一條斷層的 30 年破裂機率
# 可以差一個數量級，而**沒有任何資料能把 $c_v$ 釘死**。

# %% [markdown]
# ### 用模擬把 BPT 的推導驗證一次
#
# {eq}`eq:bpt-pdf` 是從布朗運動首達推出來的，那就直接模擬布朗運動、
# 量它的首達時間，看直方圖對不對得上理論曲線。這一步也順便把
# $\rho$、$\sigma$、$a$ 三個物理量與 $(T_r,c_v)$ 的換算驗證掉。

# %% tags=["hide-input"]
rng = np.random.default_rng(2020)
A_TH = 1.0                                  # 破裂門檻（無量綱化）
RHO = A_TH / T_R                            # 加載率：門檻 / 平均複發期
SIG = CV * np.sqrt(A_TH * RHO)              # 噪音強度：c_v = σ/√(aρ)
DT, NSTEP, NPATH = 0.05, 6000, 4000         # 時間步、步數（至 300 年）、路徑數

dW = rng.normal(0.0, np.sqrt(DT), size=(NPATH, NSTEP))
X = np.cumsum(RHO * DT + SIG * dW, axis=1)
hit = X >= A_TH
has_hit = hit.any(axis=1)
tau = (hit.argmax(axis=1) + 1) * DT
tau_ok = tau[has_hit]

fig = make_subplots(rows=1, cols=2, column_widths=[0.52, 0.48],
                    subplot_titles=("應力軌跡 X(t)=ρt+σW(t) 與首達",
                                    "首達時間分布 vs 理論 BPT"))
tt = np.arange(1, NSTEP + 1) * DT
for j in range(12):                          # 只畫 12 條，其餘用於統計
    n_show = (hit.argmax(axis=1)[j] + 1) if has_hit[j] else NSTEP
    fig.add_trace(go.Scatter(x=tt[:n_show], y=X[j, :n_show], mode="lines",
                             showlegend=False, opacity=0.65,
                             line=dict(color=PALETTE[j % len(PALETTE)],
                                       width=1)), row=1, col=1)
    if has_hit[j]:
        fig.add_trace(go.Scatter(x=[tau[j]], y=[A_TH], mode="markers",
                                 showlegend=False,
                                 marker=dict(color=QUAKE_COLOR, size=7,
                                             symbol="triangle-up")),
                      row=1, col=1)
fig.add_trace(go.Scatter(x=tt, y=RHO * tt, mode="lines", name="平均加載 ρt",
                         line=dict(color="#555", width=1.6, dash="dash")),
              row=1, col=1)
fig.add_hline(y=A_TH, line=dict(color=QUAKE_COLOR, width=1.2, dash="dot"),
              row=1, col=1)

fig.add_trace(go.Histogram(x=tau_ok, histnorm="probability density",
                           nbinsx=60, name="模擬首達時間",
                           marker=dict(color=ACCENT, opacity=0.55)),
              row=1, col=2)
fig.add_trace(go.Scatter(x=tg, y=bpt.pdf(tg), mode="lines",
                         name="理論 BPT pdf",
                         line=dict(color=QUAKE_COLOR, width=2.4)), row=1, col=2)
fig.update_xaxes(title_text="時間（年）", range=[0, 300], row=1, col=1)
fig.update_yaxes(title_text="應力 X(t)（門檻 = 1）", row=1, col=1)
fig.update_xaxes(title_text="首達時間 τ（年）", range=[0, 300], row=1, col=2)
fig.update_yaxes(title_text="密度", row=1, col=2)
apply_layout(fig, height=430, title=(
    f"圖 4　{NPATH} 條布朗軌跡的首達時間："
    f"平均 {tau_ok.mean():.1f} 年（理論 {T_R:.0f}）、"
    f"c_v {tau_ok.std()/tau_ok.mean():.3f}（理論 {CV:.2f}）"))
fig

# %% [markdown]
# 左圖是十二條應力軌跡。灰色虛線是「沒有噪音時」的加載直線，
# $\rho t$ 在 $t=T_r=100$ 年準時抵達門檻。實際軌跡在直線附近亂走，
# 有的提早撞線、有的拖到兩百多年——**同一組參數、同一條斷層，
# 複發間隔的離散完全來自那個 $\sigma W(t)$ 項**。
#
# 右圖把四千條軌跡的首達時間做成直方圖，疊上 {eq}`eq:bpt-pdf` 的理論
# 曲線。兩者吻合，平均值與 $c_v$ 都對得上（圖標題的數字由程式帶入）。
# 模擬平均略低於理論值，是因為離散時間步只在格點上檢查是否越界，
# 會漏掉步與步之間的短暫越界，這個誤差隨 $\Delta t$ 縮小而消失；
# 另外有少數軌跡在 300 年內未觸及門檻而被截掉。
#
# 這張圖的教學價值在於它把「分布」還原成「過程」。BPT 不是有人挑了
# 一個好看的偏態分布，而是**一個明確物理機制的必然後果**——你只要
# 接受「應力線性累積 ＋ 隨機擾動 ＋ 門檻破裂」，逆高斯分布就自己
# 跑出來了。

# %% [markdown]
# ### 參數與典型值
#
# 以下是**文獻報告值**與換算關係，供對照與健全性檢查。
#
# | 量 | 符號 | 典型值／範圍 | 備註 |
# |---|---|---|---|
# | 平均複發時間 | $T_r$ | $10^2$–$10^4$ 年 | 由滑移率與特徵位移推得 |
# | aperiodicity | $c_v$ | 0.3–0.7 | UCERF 系列的慣用範圍 |
# | Weibull 形狀 | $k$ | $>1$（準週期） | $k=1$ 即 Poisson |
# | 加載率 | $\rho$ | $a/T_r$ | 門檻除以平均複發期 |
# | 噪音強度 | $\sigma$ | $c_v\sqrt{a\rho}$ | 未建模擾動的總和 |
#
# | 對照關係 | 式子 | $b=1$ 時 |
# |---|---|---|
# | Zöller 橋（20.5 節） | $c_v=\sqrt{b/(3-b)}$ | $c_v\approx0.71$ |
# | BPT hazard 漸近值 | $1/(2T_rc_v^2)$ | $1/T_r$ 的 $1/(2c_v^2)$ 倍 |
# | Weibull 變異係數 | $\sqrt{\Gamma(1+2/k)/\Gamma(1+1/k)^2-1}$ | $c_v=0.5\Rightarrow k\approx2.10$ |
#
# 台灣端的落點：TEM PSHA2020（Chan et al. 2020）對孕震構造源導入 BPT
# 以描述斷層記憶，是時變危害度在台灣的**官方實作**；第 21 章會把它
# 接進危害積分裡。

# %% [markdown]
# ## 20.5 從 b 值到 aperiodicity：Zöller 橋
#
# 20.4 節留了一個尷尬的問題：$c_v$ 這麼關鍵，卻沒有資料能定它。
# 一條斷層的三到五個歷史間隔，估變異係數的標準誤大到毫無意義。
# 文獻慣用的 0.3–0.7 說穿了是**專家判斷**，不是量測。
#
# Zöller et al. (2008) 提出了一條出乎意料的捷徑：**用儀器目錄的 b 值
# 去推 $c_v$**。
#
# $$c_v \;=\; \sqrt{\frac{b}{3-b}},\qquad 0<b<3$$
#
# 論證的骨架是這樣。斷層上的應力不只被大地震清空，也被斷層面上的
# 中小地震一點一點卸掉，同時被斷層外的地震一點一點加上去。假設這
# **兩種效應大致相抵**，剩下的淨擾動主要來自斷層面上中小地震的
# 規模分布——而那個分布由 GR 律描述，斜率就是 $b$。地震矩與規模的
# 關係讓「一次事件釋放多少應力」的分布有一條 $-(2b/3+1)$ 的冪次尾巴
# （$3/2$ 來自 $\log M_0 = 1.5M+\text{const}$），把它的擾動累積起來，
# 就得到布朗噪音的等效強度，再換算成 $c_v$。
#
# 這條式子的性質值得逐一讀：
#
# - $b$ 愈大 ⇒ 小地震相對愈多 ⇒ 應力以細碎的方式釋放 ⇒ 擾動相對
#   加載愈大 ⇒ $c_v$ 愈大 ⇒ 斷層愈不規律。方向合理。
# - $b=1$（全球與台灣的長期平均，第 11 章）給出
#   $c_v=\sqrt{1/2}\approx0.71$——**恰好落在文獻慣用範圍 0.3–0.7 的
#   上緣**。兩條完全獨立的路線給出相容的數字，這件事本身就值得注意。
# - $b\to3$ 時 $c_v\to\infty$，模型退化成純隨機；$b\to0$ 時
#   $c_v\to0$，退化成完美週期。兩個極限都合理。
#
# 為什麼這座橋重要？因為它是**第 11 章（b 值，儀器目錄，數萬筆
# 小地震，時間尺度數十年）與長期危害（古地震，三五筆記錄，時間
# 尺度數千年）之間罕見的量化連結**。這兩端平常各說各話：搞統計
# 地震學的人算 b 值，搞古地震的人數探槽層位，兩邊的資料在數量級
# 上差了四個級數。Zöller 橋說：**你手上那份儀器目錄的斜率，可以
# 回頭餵給複發時間模型。**
#
# 當然要提醒代價。「斷層內卸載與斷層外加載相抵」是一個很強的假設，
# 而且無法在單一斷層上驗證；$b$ 值本身還有第 11 章那一整章的估計
# 陷阱（$M_c$ 選錯、除叢造成的偏誤、樣本數）。所以合理的用法不是
# 「用它取代專家判斷」，而是**把它當成 logic tree 上的一個分支**：
# 讓「由 b 值推導的 $c_v$」與「文獻慣用值」互相制衡。

# %% [markdown]
# ## 20.6 自我修正過程
#
# 更新過程有一個結構性的限制：它只記得**上一次**。真實斷層系統顯然
# 不是這樣——一次大破裂在鄰段留下的應力變化，會影響好幾個世代的
# 行為。要描述這件事，就得回到第 10 章的條件強度框架。
#
# **應力釋放模型**（stress release model, Vere-Jones 1978）的想法是
# 把彈性回跳直接翻譯成 $\lambda^*$。設應力水準
# $X(t)=X(0)+\rho t-S(t)$，其中 $\rho$ 是加載率、$S(t)$ 是到 $t$ 為止
# 的**累積釋放量**（實務上取累積地震矩或 Benioff 應變）。假設破裂
# 率隨應力水準指數成長，就得到
#
# $$\lambda^*(t) \;=\; \exp\bigl\{a_s + b_s\,[\,t - c_s\,S(t)\,]\bigr\}$$
#
# （文獻原本用 $a,b,c$，與 GR 的 $b$ 值、Omori 的 $c$ 撞名，本書一律
# 改寫成 $a_s,b_s,c_s$。）三個參數各有所指：$a_s$ 定基線水準、
# $b_s$ 定「應力對率的敏感度」、$c_s$ 定「單位釋放量折算成多少
# 加載時間」。長期平衡的條件是 $c_s\times(\text{單位時間平均釋放量})=1$。
#
# **Linked SRM**（Liu et al. 1998）把區域切成若干子區，加入應力轉移：
#
# $$\lambda_i^*(t) = \exp\left\{a_{s,i} + b_{s,i}
#   \left[\,t - \sum_j c_{s,ij}\,S_j(t)\,\right]\right\}$$
#
# $c_{s,ij}$ 是「子區 $j$ 的釋放對子區 $i$ 的應力影響」。對角項為正
# （自己放掉自己的應力），**離對角項可正可負**——這是模型能容納
# 庫倫應力觸發與應力陰影的地方。
#
# ### 與 ETAS 同框架、反號
#
# 10.3 節已經指出，self-exciting 與 self-correcting 可以寫成同一條
# 句型 $\lambda^*(t)=\Phi[\eta_0\pm\sum_i w(t-t_i,m_i)]$，ETAS 取
# 加號、應力釋放取減號。這裡把「反號」量化。
#
# 兩個模型都在事件發生的瞬間讓 $\lambda^*$ 乘上一個因子：
#
# - **ETAS**：事件 $i$ 讓強度**加上** $\kappa(m_i)g(0^+)$。因為
#   Omori 核在 $t\to0^+$ 極陡，這是一個很大的正跳，跳完隨即衰減。
#   跳幅隨規模指數成長（$e^{\alpha(m_i-m_0)}$）。
# - **應力釋放**：事件讓 $S$ 增加 $\Delta S_i$，強度**乘上**
#   $e^{-b_sc_s\Delta S_i}<1$。這是一個乘性的下跳，跳完之後以
#   $e^{b_st}$ 的速率緩慢爬回來。
#
# 差別不只在符號，還在**時間尺度感**與**恢復方式**：ETAS 是「瞬間
# 暴衝、冪次衰減」，應力釋放是「瞬間掉落、指數回升」。
#
# 最乾淨的量化對照是**間隔的變異係數**——也就是本章的主角 $c_v$。
# 用同一個長期平均速率模擬三種過程，量它們的 $c_v$：

# %% tags=["hide-input"]
def sim_hawkes(rng, T, mu0, K, tau):
    """指數核 Hawkes，thinning（10.6 節）。事件間 λ* 遞減，上界取事件後值。"""
    ev, t, s, t_last = [], 0.0, 0.0, 0.0
    while True:
        ub = mu0 + s
        t += rng.exponential(1.0 / ub)
        s *= np.exp(-(t - t_last) / tau)
        t_last = t
        if t >= T:
            return np.array(ev)
        if rng.random() <= (mu0 + s) / ub:
            ev.append(t)
            s += K


def sim_stress_release(rng, T, a_s, b_s, c_s):
    """應力釋放：事件間 Λ 可解析反轉，直接用反函數法（10.6 節），不需 thinning。"""
    ev, t, lam = [], 0.0, np.exp(a_s)
    while True:
        x = np.log1p(b_s * rng.exponential() / lam) / b_s
        t += x
        if t >= T:
            return np.array(ev)
        lam *= np.exp(b_s * x - b_s * c_s)     # 加載 e^{b_s x}、釋放 e^{-b_s c_s}
        ev.append(t)


rng = np.random.default_rng(717)
T_SIM, RATE = 3000.0, 1.0                       # 年、目標長期速率（次/年）
MU0, N_BR, TAU_K = 0.40, 0.60, 0.30             # Hawkes：μ0/(1-n) = 1.0
A_S, B_S, C_S = 0.0, 1.5, 1.0                   # 應力釋放：c_s × 速率 = 1

t_poi = np.cumsum(rng.exponential(1 / RATE, size=int(T_SIM * RATE * 1.4)))
t_poi = t_poi[t_poi < T_SIM]
t_haw = sim_hawkes(rng, T_SIM, MU0, N_BR / TAU_K, TAU_K)
t_srm = sim_stress_release(rng, T_SIM, A_S, B_S, C_S)

names = ["無記憶 Poisson", "正記憶 self-exciting", "負記憶 self-correcting"]
cvs = [np.diff(e).std() / np.diff(e).mean() for e in (t_poi, t_haw, t_srm)]

W0, W1 = 900.0, 1000.0                          # 放大窗（年）
gz = np.linspace(W0, W1, 2000)
lam_h = MU0 + sum((N_BR / TAU_K) * np.exp(-(gz - ti) / TAU_K) * (gz > ti)
                  for ti in t_haw[(t_haw > W0 - 8) & (t_haw < W1)])
n_past = np.searchsorted(t_srm, gz, side="right")
lam_s = np.exp(A_S + B_S * (gz - C_S * n_past))

fig = make_subplots(rows=1, cols=2, column_widths=[0.62, 0.38],
                    subplot_titles=("同一速率下的 λ*（對數軸，100 年窗）",
                                    "間隔的變異係數 c_v"))
for lam, nm, col, ev in [(lam_h, names[1], PALETTE[1], t_haw),
                         (lam_s, names[2], PALETTE[2], t_srm)]:
    fig.add_trace(go.Scatter(x=gz, y=lam, mode="lines", name=nm,
                             line=dict(color=col, width=1.9)), row=1, col=1)
    sel = ev[(ev > W0) & (ev < W1)]
    fig.add_trace(go.Scatter(x=sel, y=np.full_like(sel, 0.05), mode="markers",
                             showlegend=False,
                             marker=dict(color=col, size=6,
                                         symbol="triangle-up")), row=1, col=1)
fig.add_hline(y=RATE, line=dict(color="#8a8a8a", width=1.4, dash="dot"),
              row=1, col=1)
fig.add_trace(go.Bar(x=names, y=cvs, showlegend=False,
                     marker=dict(color=[PALETTE[0], PALETTE[1], PALETTE[2]]),
                     text=[f"{c:.2f}" for c in cvs], textposition="outside"),
              row=1, col=2)
fig.add_hline(y=1.0, line=dict(color="#8a8a8a", width=1.4, dash="dot"),
              row=1, col=2)
fig.update_xaxes(title_text="時間（年）", row=1, col=1)
fig.update_yaxes(title_text="λ*（次/年，對數軸）", type="log", row=1, col=1)
fig.update_yaxes(title_text="c_v = SD/mean", range=[0, 2.2], row=1, col=2)
apply_layout(fig, height=430, title=(
    f"圖 5　同一長期速率 {RATE:.0f} 次/年、{T_SIM:.0f} 年模擬："
    f"三種記憶結構的 c_v = {cvs[0]:.2f} / {cvs[1]:.2f} / {cvs[2]:.2f}"))
fig

# %% [markdown]
# 左圖用**對數縱軸**，因為兩條曲線的動態範圍差太多才擺得下。橘線
# （self-exciting）是尖峰急衰：每個事件把 $\lambda^*$ 推高一到兩個
# 數量級，然後指數衰回背景。綠線（self-correcting）是規律的鋸齒：
# 事件讓它掉一格，之後穩定爬升，幅度小得多。灰色虛線是兩者共同的
# 長期平均速率。
#
# 右圖是本章的量化重點。三個過程的長期速率完全一樣，間隔的**變異
# 係數**卻分成三檔：
#
# - Poisson 的 $c_v\approx1$（指數分布的 $c_v$ 恆為 1，這是理論值）。
# - self-exciting 的 $c_v$ 遠**大於** 1——叢集造成大量極短間隔與
#   少數極長間隔。這正是 20.3 節 Weibull $k<1$ 描述的情形。
# - self-correcting 的 $c_v$ 明顯**小於** 1，落在 0.6 附近——**恰好
#   就是 BPT 文獻慣用的 0.3–0.7 範圍**。
#
# 這個對應不是巧合，它是本章與第 10 章之間的縫合線：**aperiodicity
# $c_v$ 不只是 BPT 的一個參數，它是「記憶方向」的量尺**。$c_v>1$
# 是正記憶、$c_v=1$ 是無記憶、$c_v<1$ 是負記憶。你拿到一串複發間隔，
# 第一件該算的事就是它的 $c_v$；這個數字落在 1 的哪一側，就決定了
# 你該去第 13 章還是留在第 20 章。

# %% [markdown]
# ## 20.7 資料稀少的困境
#
# 現在講這一章最不舒服的部分。
#
# 前面所有的數學都預設了一件事：我們有辦法從資料裡估出 $f$ 的參數。
# 對儀器目錄這沒問題——第 14 章估 ETAS 用的是幾萬筆事件。但
# fault-based 模型的資料量是這樣的：**一條斷層的歷史破裂記錄，
# 可能只有 3 到 5 次。**
#
# 這些記錄從哪來？古地震探槽（開挖斷層帶、辨識被錯開的地層、
# 碳十四定年）、歷史文獻、地形位移。每一次的年代都帶著數十到
# 數百年的誤差。三到五個帶大誤差的間隔，要估兩個參數。
#
# 後果很直接。$T_r$ 的相對標準誤大約是 $c_v/\sqrt{n}$，$n=4$、
# $c_v=0.5$ 就是 25%；而 $c_v$ 本身的標準誤更差，粗估是
# $c_v/\sqrt{2n}\approx 0.18$——也就是說「$c_v=0.5$」這個估計的
# 95% 區間大概涵蓋 0.15 到 0.85，**幾乎橫跨整個文獻慣用範圍**。
# 至於「Weibull 還是對數常態還是 BPT」，圖 1 已經回答了：
# **在這種樣本數下，三個分布在統計上無法區分。分布的選擇由假設
# 決定，不由資料決定。**
#
# 這與{doc}`第 18 章 <18_testing_comparison>`那條原則是同一件事的
# 另一個版本：**有效樣本數是地震顆數，不是資料點數**。第 18 章的
# 情境是「二十萬個時空箱、12 顆目標地震，有效樣本數是 12」；這裡
# 是「幾千年的探槽剖面、4 個複發間隔，有效樣本數是 4」。兩者都
# 在提醒同一件事：**看起來資訊很多的資料集，可以只攜帶極少的
# 統計自由度。**
#
# ### logic tree：把無知寫進模型
#
# 既然資料無法決定，就必須誠實地把「不知道」表達出來。這就是
# **logic tree**（邏輯樹）的用途：對每一個無法由資料決定的選擇
# （用哪個分布、$c_v$ 取多少、$T_r$ 取多少、斷層分段怎麼切），
# 列出所有合理選項、各給一個權重，然後把整棵樹跑完，得到的不是
# 一條危害曲線而是一族曲線。
#
# 這棵樹可以長得非常大。**UCERF3 有 5760 個 logic-tree 分支**——
# 這個數字本身就是一把量尺：它量的是**認知不確定性**（epistemic
# uncertainty，「我們不知道真相是什麼」）的規模，而不是隨機不確定性
# （aleatory，「即使知道真相，結果仍然隨機」）。一個需要 5760 個
# 分支才能表達的模型，等於在說「這件事我們真的不太確定」。
#
# 把這件事說白：**時間相依複發模型的最大不確定性，不在數學裡，
# 在資料裡。** 本章推的每一條式子都是精確的；把精確的式子套在
# 四個帶誤差的數字上，得到的結論依然模糊。這不是式子的錯，也不是
# 可以靠更好的統計方法解決的問題——它需要更多的古地震探槽。

# %% [markdown]
# ## 20.8 其他家族（點到為止）
#
# 為了讓地圖完整，另外兩條線值得知道名字。
#
# **隱藏馬可夫與半馬可夫模型**（HMM／semi-Markov）。想法是：斷層
# 系統有若干個**觀測不到的狀態**，最自然的解釋就是「應力水準」
# （例如低／中／高三態）。每個狀態有自己的地震發生率，狀態之間
# 依轉移機率跳動；我們只看得到地震序列，看不到狀態，所以要用
# EM 之類的方法同時估狀態序列與參數。它比更新過程多一層彈性
# （記憶不再只有「上一次」），代價是參數更多、可辨識性更差——
# 在 20.7 節那種資料量下通常無法使用，多用於區域尺度的目錄而非
# 單一斷層。
#
# **加速矩釋放**（Accelerating Moment Release, AMR）。觀察是：某些
# 大地震之前數年到數十年，區域內的累積 Benioff 應變
# $\varepsilon(t)=\sum_i\sqrt{E_i}$（$E_i$ 為地震能量）會偏離線性、
# 呈現加速。經驗式寫成
#
# $$\varepsilon(t) = A - B\,(t_f - t)^{m},\qquad 0.1\le m\le 0.5$$
#
# $t_f$ 是預測的大地震時間、$m$ 是冪次。它的吸引力在於能同時給出
# 時間與規模；問題在於**回溯擬合幾乎總是成功**（$t_f$、$A$、$B$、
# $m$ 四個自由參數，加上「區域半徑」與「時間窗」兩個通常被自由
# 挑選的量），前瞻表現則遠不如宣稱。它是{doc}`第 9 章
# <09_forecasting_intro>`那套「回溯漂亮、前瞻不行」故事的標準案例，
# 放在這裡供辨識，不建議當作預報工具。

# %% [markdown]
# ## 20.9 常見誤解與陷阱
#
# **誤解 1：「BPT 給的 30 年機率比 Poisson 高，所以它比較保守。」**
# 不一定。圖 2 顯示，$t$ 小於約半個複發期時 BPT 的 hazard **低於**
# Poisson——剛破裂完的斷層，時變模型給的機率比時間獨立模型**低**。
# 1999 集集之後的車籠埔就是這種情形。時變模型不是「一律加碼」的
# 安全係數，它是**重新分配**：把危害從剛破裂的斷層搬到沉寂已久的
# 斷層。用它之前必須接受這個重新分配的兩端，只取一端是不誠實的。
#
# **誤解 2：「複發間隔的樣本數就是探槽找到的事件數。」** 估一個
# 分布的兩個參數，用的是**間隔**不是**事件**：$n$ 次古地震只給出
# $n-1$ 個完整間隔。而且最後一個間隔通常是**設限資料**（censored，
# 只知道「已經超過 $T$ 年」還沒結束），不能當成一次觀測直接丟進
# 樣本平均——否則會系統性低估 $T_r$。四次古地震其實只有三個間隔，
# 這時候多算一個或少算一個都是 33% 的差別。
#
# **誤解 3：「特徵地震模型與 GR 律是相容的，只是尺度不同。」**
# 兩者有真實的張力。**特徵地震模型**假設一段斷層反覆產生規模幾乎
# 相同的破裂（這是整個更新過程框架的前提——沒有「特徵規模」就沒有
# 「複發間隔」可言）；GR 律則說規模分布是連續的冪次，沒有任何
# 特徵尺度。把兩者放在一起做危害度計算，會在特徵規模附近出現
# **凸起**（bump），而那個凸起是模型拼接的產物還是真實現象，
# 至今沒有定論。第 21 章的 PSHA 必須同時使用兩種震源（構造源用
# 特徵模型、面震源用 GR），拼接處的處理是那套計算裡最主觀的一步。
#
# **誤解 4：「$c_v$ 是斷層的物理性質，可以量出來。」** $c_v$ 是
# **模型參數**，不是可觀測量。20.7 節算過它的標準誤幾乎橫跨整個
# 合理範圍；20.5 節的 Zöller 橋提供另一條估計路線，但依賴一個
# 無法驗證的假設。實務上 $c_v$ 的值來自專家判斷與 logic tree，
# 把它寫成「這條斷層的 aperiodicity 是 0.5」而不加不確定度，
# 是把假設偽裝成量測。
#
# **誤解 5：「時間相依模型取代了 Poisson 模型。」** 沒有。目前
# 世界各國的官方危害度圖，主體仍然是時間獨立的（第 21 章）；
# 時變成分只加在**少數幾條資料最好的斷層**上。原因就是 20.7 節：
# 對絕大多數構造，我們連 $T_r$ 都估不準，遑論 $c_v$。BPT 是
# TEM PSHA2020 的一個組件，不是它的骨架。

# %% [markdown]
# ## 20.10 研究前沿與未解問題
#
# **用力學取代參數分布。** 本章所有模型的共同弱點是：複發分布的
# 形狀是**假設**的。物理式地震模擬器想繞過這一步——RSQSim 之類的
# 系統直接在三維斷層網路上解 rate-and-state 摩擦律，跑出數十萬年
# 的合成地震目錄，複發時間分布是**輸出**而不是輸入。UCERF3-ETAS
# 則走另一條路，把時變的長期複發率餵進 ETAS 當背景率，讓短期觸發
# 與長期記憶共存在同一份模擬裡。兩者都還在驗證階段，共同的難題
# 也一樣：模擬器裡的摩擦參數與斷層幾何，本身也是猜的，只是把
# 不確定性從「哪個分布」搬到了「哪組摩擦參數」。
#
# **把兩種記憶裝進同一個模型。** 這是本章開頭那個張力的正面版本。
# 短期是 self-exciting、長期是 self-correcting，兩者物理上都對，
# 但寫成單一 $\lambda^*$ 卻異常困難：$\lambda^*$ 得在事件後先暴衝
# （餘震）、再掉到低於背景（應力耗竭）、再緩慢回升，而且三段的
# 時間尺度差了四到六個數量級。目前的作法多半是**把兩層疊起來**
# （長期 renewal 決定背景率、短期 ETAS 疊在上面），但這樣寫的
# 兩層之間沒有真正的耦合——大地震既觸發餘震又消耗應力，這件事
# 沒有被同一條式子描述。
#
# **中間那段仍然是空白。** 10.3 節說過：短期看觸發、長期看耗竭，
# **月到年的灰色地帶兩者都不管用**。EEPAS（第 16 章）與 Ψ（第 15
# 章）試圖填補，但至今沒有一個模型在這個尺度上有壓倒性的前瞻證據。
#
# **驗證的根本困難。** 這是最誠實也最令人洩氣的一點：一個複發期
# 一千年的模型，**在人類的時間尺度上無法被前瞻檢驗**。第 17、18
# 章那套 CSEP 方法論需要在測試期內看到足夠多的目標事件；對長期
# 複發模型，測試期得以千年計。目前的替代方案是「拿整個斷層系統
# 的許多條斷層當作一次實驗的許多個樣本」，但那要求各斷層彼此獨立
# ——而應力轉移正好破壞這個假設。**這一章的模型，可能是全書最
# 難被證偽的一類。**

# %% [markdown]
# ## 20.11 附錄：本章推導細節
#
# ### A. 布朗運動加漂移的首達時間
#
# 目標：證明 $X(t)=\rho t+\sigma W(t)$ 首次觸及 $a>0$ 的時間
# $\tau$ 的密度是 $f(\tau)=\frac{a}{\sigma\sqrt{2\pi\tau^3}}
# \exp[-(a-\rho\tau)^2/(2\sigma^2\tau)]$。
#
# **第一步：無漂移的情形（反射原理）。** 設 $\rho=0$，
# $M_t=\max_{s\le t}\sigma W(s)$。反射原理說：每一條在 $[0,t]$ 內
# 觸及過 $a$、終點在 $a$ 以下的路徑，都可以沿首達時刻之後的部分
# 對水平線 $a$ 做鏡射，一對一地對應到一條終點在 $a$ 以上的路徑。
# 因此
#
# $$P(M_t\ge a) = 2\,P\bigl(\sigma W(t)\ge a\bigr)
#   = 2\left[1-\Phi\!\left(\frac{a}{\sigma\sqrt t}\right)\right].$$
#
# 對 $t$ 微分（注意 $\{M_t\ge a\}=\{\tau\le t\}$）：
#
# $$\begin{aligned}
# f_0(t) &= \frac{\mathrm{d}}{\mathrm{d}t}
#   \left[2-2\Phi\!\left(\frac{a}{\sigma\sqrt t}\right)\right]
#  = -2\,\varphi\!\left(\frac{a}{\sigma\sqrt t}\right)\cdot
#    \left(-\frac{a}{2\sigma t^{3/2}}\right) \\
#  &= \frac{a}{\sigma\sqrt{2\pi t^{3}}}\,
#    \exp\left[-\frac{a^{2}}{2\sigma^{2}t}\right].
# \end{aligned}$$
#
# **第二步：加上漂移（測度變換）。** 令 $Q$ 為使 $X$ 無漂移的測度，
# $P$ 為真實測度。Girsanov 定理給出在 $\mathcal{F}_t$ 上的
# Radon–Nikodym 導數
#
# $$\left.\frac{\mathrm{d}P}{\mathrm{d}Q}\right|_{\mathcal F_t}
#   = \exp\left[\frac{\rho}{\sigma^{2}}X(t)
#     - \frac{\rho^{2}t}{2\sigma^{2}}\right].$$
#
# 關鍵在於：在事件 $\{\tau=t\}$ 上，路徑的終點被釘死了——
# $X(\tau)=a$。所以這個導數在該事件上是一個**只依賴 $t$ 的常數**
# $\exp[\rho a/\sigma^{2}-\rho^{2}t/(2\sigma^{2})]$。於是
#
# $$\begin{aligned}
# f(t) &= \exp\left[\frac{\rho a}{\sigma^{2}}
#         - \frac{\rho^{2}t}{2\sigma^{2}}\right] f_0(t) \\
#      &= \frac{a}{\sigma\sqrt{2\pi t^{3}}}\exp\left[
#         -\frac{a^{2}}{2\sigma^{2}t} + \frac{\rho a}{\sigma^{2}}
#         - \frac{\rho^{2}t}{2\sigma^{2}}\right] \\
#      &= \frac{a}{\sigma\sqrt{2\pi t^{3}}}\exp\left[
#         -\frac{a^{2}-2\rho a t+\rho^{2}t^{2}}{2\sigma^{2}t}\right]
#       = \frac{a}{\sigma\sqrt{2\pi t^{3}}}\exp\left[
#         -\frac{(a-\rho t)^{2}}{2\sigma^{2}t}\right].
# \end{aligned}$$
#
# 第三個等號只是把三項通分成 $2\sigma^2 t$ 的完全平方。代入
# $T_r=a/\rho$ 與 $c_v^2=\sigma^2/(a\rho)$ 即得 {eq}`eq:bpt-pdf`。
#
# 順帶算兩個矩：$E[\tau]=a/\rho=T_r$（直觀：平均而言噪音互相抵消，
# 需要 $a/\rho$ 的時間走到門檻），$\mathrm{Var}[\tau]=a\sigma^2/\rho^3$，
# 故 $\mathrm{SD}/\mathrm{mean}=\sqrt{a\sigma^2/\rho^3}\big/(a/\rho)
# =\sigma/\sqrt{a\rho}=c_v$。**aperiodicity 正是變異係數**，
# 這不是定義而是計算結果。
#
# ### B. 與 SciPy 參數化的對照
#
# 逆高斯的標準參數化是 $\mathrm{IG}(\mu_{\rm IG},\lambda_{\rm IG})$，
# 密度 $\sqrt{\lambda/(2\pi t^3)}\exp[-\lambda(t-\mu)^2/(2\mu^2t)]$，
# 平均 $\mu_{\rm IG}$、變異數 $\mu_{\rm IG}^3/\lambda_{\rm IG}$。
# 對照 {eq}`eq:bpt-pdf` 得 $\mu_{\rm IG}=T_r$、
# $\lambda_{\rm IG}=T_r/c_v^{2}$。
#
# SciPy 的 `invgauss(mu, scale)` 定義為 $X=\text{scale}\cdot Y$、
# $Y\sim\mathrm{IG}(\texttt{mu},1)$，對應
# $\mu_{\rm IG}=\texttt{mu}\times\texttt{scale}$、
# $\lambda_{\rm IG}=\texttt{scale}$。所以本章的呼叫方式是
# `invgauss(mu=c_v**2, scale=T_r/c_v**2)`——圖 1 的建構函式就是
# 這樣寫的，而圖 4 的模擬（直接跑布朗運動）是對這個對照的獨立驗證。
#
# ### C. BPT 的 hazard 漸近值
#
# 目標：證明 $t\to\infty$ 時 $h(t)\to\lambda_{\rm IG}/(2\mu_{\rm IG}^2)
# =1/(2T_rc_v^{2})$。以下用 $\mu,\lambda$ 簡寫。
#
# 逆高斯的累積分布是
#
# $$F(t)=\Phi(u)+e^{2\lambda/\mu}\,\Phi(-v),\qquad
#   u=\sqrt{\tfrac{\lambda}{t}}\Bigl(\tfrac{t}{\mu}-1\Bigr),\quad
#   v=\sqrt{\tfrac{\lambda}{t}}\Bigl(\tfrac{t}{\mu}+1\Bigr),$$
#
# 故 $S(t)=\Phi(-u)-e^{2\lambda/\mu}\Phi(-v)$。
#
# **關鍵恆等式**。展開平方：
# $u^{2}/2=\frac{\lambda t}{2\mu^{2}}-\frac{\lambda}{\mu}
# +\frac{\lambda}{2t}$、
# $v^{2}/2=\frac{\lambda t}{2\mu^{2}}+\frac{\lambda}{\mu}
# +\frac{\lambda}{2t}$，兩者相差恰好 $2\lambda/\mu$。因此
#
# $$e^{2\lambda/\mu}\,\varphi(v)
#  = e^{2\lambda/\mu}\cdot\tfrac{1}{\sqrt{2\pi}}e^{-u^{2}/2-2\lambda/\mu}
#  = \varphi(u).$$
#
# 指數項乾淨地互相抵消——這是整個推導能收斂的原因。
#
# **取極限**。$t\to\infty$ 時 $u,v\to\infty$，用 Mills 比
# $\Phi(-x)=\varphi(x)/x\cdot[1+O(x^{-2})]$：
#
# $$\begin{aligned}
# S(t) &\approx \frac{\varphi(u)}{u}
#   - e^{2\lambda/\mu}\frac{\varphi(v)}{v}
#  = \varphi(u)\left(\frac{1}{u}-\frac{1}{v}\right)
#  = \varphi(u)\,\frac{v-u}{uv} \\
#  &= \varphi(u)\,\frac{2\sqrt{\lambda/t}}{uv}
#  \;\approx\; \varphi(u)\,\frac{2\sqrt{\lambda/t}\;\mu^{2}}{\lambda t}
#  = \frac{2\mu^{2}\varphi(u)}{\sqrt{\lambda}\;t^{3/2}},
# \end{aligned}$$
#
# 其中用了 $v-u=2\sqrt{\lambda/t}$ 與 $uv\to\lambda t/\mu^{2}$。
#
# 分子同樣用 $u$ 改寫。注意
# $u^{2}=\frac{\lambda}{t}\bigl(\frac{t}{\mu}-1\bigr)^{2}
# =\frac{\lambda(t-\mu)^{2}}{\mu^{2}t}$，正好是 pdf 指數項的兩倍，
# 所以
#
# $$f(t)=\sqrt{\frac{\lambda}{2\pi t^{3}}}\,e^{-u^{2}/2}
#      =\sqrt{\frac{\lambda}{t^{3}}}\;\varphi(u).$$
#
# 相除，$\varphi(u)$ 與 $t^{3/2}$ 全部消掉：
#
# $$h(t)=\frac{f(t)}{S(t)}
#  \;\longrightarrow\;
#  \frac{\sqrt{\lambda}\,t^{-3/2}\varphi(u)}
#       {2\mu^{2}\varphi(u)\,\lambda^{-1/2}t^{-3/2}}
#  = \frac{\lambda}{2\mu^{2}}
#  = \frac{T_r/c_v^{2}}{2T_r^{2}} = \frac{1}{2T_rc_v^{2}} .$$
#
# 對照另外兩個分布：Weibull 的 $h$ 是純冪次，$k>1$ 時發散；
# 對數常態的 $h\sim(\ln t-m)/(\sigma^{2}t)\to0$。三種結局的來源
# 都是「$f$ 與 $S$ 誰趨於零得比較快」——BPT 的兩者以完全相同的
# $t^{-3/2}\varphi(u)$ 速率趨零，所以商數收斂到常數。
#
# ### D. 更新過程的概似
#
# 為完整性補一句：給定 $n$ 個完整間隔 $t_1,\dots,t_n$ 與一個
# 設限的殘餘時間 $T_c$（距最後一次破裂至今），對數概似是
#
# $$\ln L = \sum_{i=1}^{n}\ln f(t_i) + \ln S(T_c).$$
#
# 這正是 {eq}`eq:pp-loglik` 在「$\lambda^*(t)=h(t-t_{\rm last})$」
# 這個特例下的樣子——把 $\ln\lambda^*$ 的求和寫成 $\sum\ln h(t_i)$、
# 把 $-\int\lambda^*$ 寫成 $\sum\ln S(t_i)+\ln S(T_c)$，兩者合併
# 即得上式。**設限項 $\ln S(T_c)$ 不可省略**（20.9 節誤解 2）：
# 省掉它等於假裝「最後一次破裂到現在都不算數」，會系統性低估
# $T_r$。實作上就是把第 10 章那套概似最大化直接搬過來用，本章
# 不重複。

# %% [markdown]
# 回頭看，這一章其實只換了一個視角：不看 $f$，看 $h$。換完之後，
# 三個長得幾乎一樣的分布立刻分出高下——一個在物理上荒謬（對數常態
# 說沉寂愈久愈安全）、一個在數學上失控（Weibull 的 hazard 沒有
# 上界）、一個剛好合理（BPT 退化成一個兩倍危險的 Poisson）。
# 而 BPT 之所以合理，是因為它不是被挑出來的，是被**推導**出來的。
#
# 但這一章也留下了一個不太體面的結局。20.7 節說得很白：所有精確的
# 數學，最後都套在三到五個帶著世紀級誤差的數字上。這使得時間相依
# 複發模型成為全書最奇特的一類——**理論最乾淨，證據最薄弱**。
# 它之所以還是被寫進官方危害度模型，不是因為證據夠強，而是因為
# 「假裝斷層沒有記憶」這個替代方案更糟。
#
# 這也正好把我們帶到最後一段路。本章算出來的東西——「這條斷層
# 未來 50 年破裂的機率」——還不是工程師要的答案。工程師要的是
# 「這個場址未來 50 年地表加速度超過 0.3 g 的機率」。從前者到後者，
# 中間隔著震源幾何、規模分布、距離分布、地動預測方程，以及一個
# 把全部可能情境加總起來的積分。
#
# {doc}`第 21 章 <21_psha>`就來走這座橋。你會看到本章的 BPT 以
# 一個具體的身分出現：**TEM PSHA2020 對台灣每一條孕震構造掛上的
# 「斷層記憶」模型**，也是台灣官方危害度圖裡唯一時間相依的成分。
