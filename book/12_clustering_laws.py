# %% [markdown]
# # 12. 目錄統計 II：叢集的經驗律與除叢
#
# {doc}`第 11 章 <11_catalog_completeness_b>`把目錄的「靜態」性質
# 整理乾淨了：完整規模 $M_c$ 是隨時間、空間與主震後時刻變動的場，
# GR 斜率 $b$ 有離散、不完整與規模尺度三重偏差。那一章談的是
# 「一顆一顆地震的規模怎麼分布」；這一章要問另一半：**這些地震
# 在時間與空間上怎麼擠在一起。**
#
# 叢集有三條老得可以進博物館的經驗律：Omori（1894）的餘震衰減、
# Utsu 與 Seki（1955）的餘震區尺度、Båth（1965）的最大餘震規模差。
# 三者年紀加起來超過三百歲，卻仍是現代預報模型的骨架——
# {doc}`第 13 章 <13_etas_structure>`的 ETAS 條件強度函數逐項拆開，
# 就是這三條律再加上 GR。所以這章不是考古，是在鋪第 13 章的地基。
#
# 這章還要處理一件更尷尬的事：**除叢（declustering）**。長期危害度
# 分析傳統上要把目錄拆成「主震」與「餘震」，只留前者。這個看似
# 技術性的前處理會系統性扭曲 b 值——而且是**演算法造的**。第 11 章
# 結尾那句警告（報告 b 值必須同時報告除叢法）在 12.7 節會被完整
# 推導出來，包括一個可以用紙筆算出的**破產點** $m_x$。會用到的舊
# 式子只有兩條：第 11 章的 $\beta=b\ln10$，以及第 10 章的
# $\ln L=\sum_i\ln\lambda^*-\int\lambda^*$。
#
# ## 12.1 從 Omori 到 Omori–Utsu：為什麼非要 $p>1$ 不可
#
# 大森房吉在 1894 年看濃尾地震的餘震，寫下地震學第一條定量經驗律：
# 餘震率隨時間**雙曲線**衰減，$n(t)=K/(t+c)$。六十七年後宇津德治
# （Utsu 1961）發現指數不必是 1，推廣成今天通用的**修正 Omori 律**：
#
# $$n(t) = \frac{K}{(t + c)^{p}}$$
#
# 三個參數分工明確：$K$ 是產能尺度（帶量綱，單位「次 × 天$^{p-1}$」，
# 所以**不同 $p$ 的 $K$ 不能互相比較**），$c$ 是讓 $t=0$ 不發散的
# 平緩期，$p$ 是衰減陡度。
#
# 進 ETAS 之前必須先正規化：ETAS 把觸發拆成「生幾個後代」乘上
# 「後代什麼時候出生」，後者必須積分為 1。所以先問 $n(t)$ 積分多少：
#
# $$\begin{aligned}
# \int_0^{\infty} \frac{K}{(t+c)^{p}}\,\mathrm{d}t
# &= K\left[\frac{(t+c)^{1-p}}{1-p}\right]_0^{\infty} \\
# &= \frac{K}{1-p}\Bigl[\lim_{t\to\infty}(t+c)^{1-p} - c^{1-p}\Bigr]
#  = \frac{K\,c^{1-p}}{p-1}\,.
# \end{aligned}$$
#
# 收斂條件全壓在那個極限上：$1-p<0$ 時 $(t+c)^{1-p}\to0$，積分有限；
# $p=1$ 時原函數是 $\ln(t+c)$，發散；$p<1$ 時發散得更快。所以
# **$p>1$ 不是經驗觀察，而是「餘震總數有限」的數學等價條件**。
# 把 $n(t)$ 除以自己的積分，得到本書統一使用的正規化 Omori 密度：
#
# $$
# g(t) = \frac{p-1}{c}\left(1 + \frac{t}{c}\right)^{-p}
# $$ (eq:omori-density)
#
# 驗算：$\int_0^\infty g=\frac{p-1}{c}\cdot\frac{c}{p-1}
# \bigl[-(1+t/c)^{1-p}\bigr]_0^\infty=1$。這條式子的擁有權在本章，
# 第 13、14 章一律引用 {eq}`eq:omori-density` 不重推。
#
# 那 $p\le1$ 怎麼辦？現實中擬合出 $p<1$ 太常見了。答案是**必須綁
# 一個時間上界**：
#
# $$\begin{aligned}
# N(0,T) &= \frac{K}{1-p}\bigl[(T+c)^{1-p} - c^{1-p}\bigr],
#   \qquad p \neq 1, \\
# N(0,T) &= K\ln\!\left(1 + \frac{T}{c}\right), \qquad p = 1 .
# \end{aligned}$$
#
# $p<1$ 時 $N(0,T)$ 隨 $T^{1-p}$ 發散，$p=1$ 時隨 $\ln T$ 發散——
# 都是無限多個餘震。在工程上這很要命：**任何「這場序列總共會有
# 幾個餘震」的估計，$p\le1$ 時都必須明講積分到哪一天。** 這也是
# 作業型系統（第 22 章）一律給「未來 7 天」「未來 30 天」這種有界
# 視窗數字的原因。先把 0403 花蓮序列擬合一次，當後面的共同素材：

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["hide-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit, minimize

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, apply_layout

GREEN = "#1baf7a"
cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
main = cat.loc[cat.ML.idxmax()]
MC_SEQ = 3.5                      # 本章分析 0403 序列統一用的門檻
aft = cat[cat.time > main.time]
t_aft = np.sort((aft[aft.ML >= MC_SEQ].time - main.time)
                .dt.total_seconds().to_numpy() / 86400)
# 起算時刻 S 由 Mc(t) 反解（見 12.2）：Mc(t) = Mm − 4.5 − 0.76 log10 t
S_START = 10 ** ((float(main.ML) - 4.5 - MC_SEQ) / 0.76)


def omori_rate(t, K, c, p):
    """未正規化的 Omori–Utsu 率 n(t) = K /(t+c)^p。"""
    return K / (t + c) ** p


def omori_integral(S, T, c, p):
    """∫_S^T (t+c)^{-p} dt，p = 1 時取對數形式。"""
    if abs(p - 1.0) < 1e-9:
        return np.log((T + c) / (S + c))
    return ((T + c) ** (1 - p) - (S + c) ** (1 - p)) / (1 - p)


def omori_mle(ts, S, T):
    """點過程 MLE：K 以剖面法消去，只對 (log c, p) 數值最佳化。

    回傳 (K, c, p, sigma_p, N)；sigma_p 取數值 Hessian 反矩陣。
    """
    ts = ts[(ts >= S) & (ts <= T)]
    n = len(ts)

    def nll(theta):
        c, p = np.exp(theta[0]), theta[1]
        integ = omori_integral(S, T, c, p)
        if integ <= 0 or not np.isfinite(integ):
            return 1e12
        return -(n * np.log(n / integ) - p * np.sum(np.log(ts + c)) - n)

    bounds = ((np.log(1e-4), np.log(3.0)), (0.3, 2.5))
    best = None
    for c0 in (0.002, 0.02, 0.1, 0.5):
        for p0 in (0.7, 1.0, 1.4):
            res = minimize(nll, [np.log(c0), p0], method="L-BFGS-B",
                           bounds=bounds)
            if best is None or res.fun < best.fun:
                best = res
    x, h, e = best.x, 1e-4, np.eye(2) * 1e-4
    hess = np.array([[(nll(x + e[i] + e[j]) - nll(x + e[i] - e[j])
                       - nll(x - e[i] + e[j]) + nll(x - e[i] - e[j]))
                      / (4 * h * h) for j in range(2)] for i in range(2)])
    try:
        sig_p = float(np.sqrt(abs(np.linalg.inv(hess)[1, 1])))
    except np.linalg.LinAlgError:
        sig_p = np.nan
    c_hat, p_hat = float(np.exp(x[0])), float(x[1])
    return n / omori_integral(S, T, c_hat, p_hat), c_hat, p_hat, sig_p, n


# 方法一：分箱後對數線性最小平方（教科書作法）
bins = np.logspace(-2, np.log10(60), 25)
counts, _ = np.histogram(t_aft[t_aft <= 60], bins=bins)
centers = np.sqrt(bins[:-1] * bins[1:])
rate_obs = counts / np.diff(bins)
ok = rate_obs > 0
p_lsq, _ = curve_fit(lambda t, K, c, p: np.log(omori_rate(t, K, c, p)),
                     centers[ok], np.log(rate_obs[ok]),
                     p0=[100, 0.05, 1.1], maxfev=20000)
# 方法二：點過程最大概似（不分箱）
K_ml, c_ml, p_ml, sp_ml, n_ml = omori_mle(t_aft, S_START, 60.0)

grid = np.logspace(np.log10(S_START), np.log10(60), 200)
fig = go.Figure(go.Scatter(x=centers[ok], y=rate_obs[ok], mode="markers",
                           name=f"觀測（ML ≥ {MC_SEQ}，分箱）",
                           marker=dict(color=ACCENT, size=8)))
fig.add_trace(go.Scatter(
    x=grid, y=omori_rate(grid, *p_lsq), mode="lines",
    name=f"分箱最小平方：p = {p_lsq[2]:.2f}，c = {p_lsq[1]:.3f} 天",
    line=dict(color=QUAKE_COLOR, dash="dash")))
fig.add_trace(go.Scatter(
    x=grid, y=omori_rate(grid, K_ml, c_ml, p_ml), mode="lines",
    name=f"點過程 MLE：p = {p_ml:.2f} ± {sp_ml:.2f}，c = {c_ml:.3f} 天",
    line=dict(color=GREEN)))
fig.add_vline(x=S_START, line_dash="dot", line_color="#888888",
              annotation_text=f"MLE 起算 S = {S_START:.3f} 天")
apply_layout(fig, title=f"0403 花蓮序列的餘震率衰減（主震 ML {main.ML:.2f}，"
                        f"主震後 60 天，N = {n_ml}）",
             xaxis_title="主震後時間（天）", yaxis_title="餘震率（次/天）",
             xaxis_type="log", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 兩條擬合線都給出 $p<1$。先別急著宣布「花蓮的餘震衰減特別慢」
# ——第 13 章會說明 $p<1$ 通常是模型設定錯了；而 12.2 與 12.3 會
# 說明，光是「從哪一刻開始算」與「算到哪一天」這兩個純粹的分析
# 選擇，就足以讓同一份資料的 $\hat p$ 移動一倍以上。
#
# ## 12.2 c 值：物理，還是缺漏事件的化名？
#
# $c$ 是三個參數裡最沒共識的一個。「物理派」說：主震破裂結束後，
# 應力重分布、孔隙壓擴散、餘滑都需要時間，餘震率理應有平緩期。
# 「假象派」只有一句話，但很致命：**主震後幾小時內，小餘震根本
# 記不到。**
#
# 把這句話寫成數學。設真實餘震率是純冪次的 $n(t)=Kt^{-p}$，而目錄
# 只收錄規模在時變完整度 $M_c(t)$ 以上的事件。依 GR 律，門檻 $m_0$
# 以上被記到的比例是 $P(m\ge M_c(t))=10^{-b[M_c(t)-m_0]}$。主震後的
# 完整度普遍可寫成 $M_c(t)=M_m-a_1-a_2\log_{10}t$（Helmstetter et al.
# 2005 的南加州版取 $a_1=4.5$、$a_2=0.76$），代入得
#
# $$\begin{aligned}
# 10^{-b\,M_c(t)}
# &= 10^{-b\,(M_m - a_1)}\cdot 10^{\,b\,a_2\log_{10}t}
#  = 10^{-b\,(M_m - a_1)}\; t^{\,b\,a_2}, \\
# n_{\rm obs}(t) &\;\propto\; t^{-p}\cdot t^{\,b\,a_2}
#   = t^{-(p - b\,a_2)} .
# \end{aligned}$$
#
# 這一步是全節重點。取 $b=1$、$a_2=0.76$，觀測到的早期斜率是
# $p-0.76$：一條真實斜率 $p=1.1$ 的衰減，在目錄還沒補齊的期間看
# 起來只有 0.34——**幾乎是平的**。而「先平一段、再轉成 $t^{-p}$」
# 正好就是 Omori–Utsu 函數的形狀，轉折點在 $t\approx c$。所以
# **早期缺漏的事件會被 $c$ 完整吸收，$\hat c$ 大約就是 $M_c(t)$
# 落到分析門檻的那一刻**。
#
# 這可以直接驗算。0403 花蓮主震 $M_L\,7.19$、分析門檻 3.5，解
# $M_m-4.5-0.76\log_{10}t=3.5$ 得 $t\approx0.086$ 天（約兩小時）；
# 12.1 的 MLE 給的 $\hat c$ 是同一量級的幾十分鐘。兩個數字沒有理由
# 這麼接近，除非它們在描述同一件事。看一次主震後前幾天的規模–時間
# 圖，「浮出水面」的過程一目了然：

# %% tags=["hide-input"]
t_days_all = (aft.time - main.time).dt.total_seconds() / 86400
fig = go.Figure(go.Scattergl(
    x=t_days_all, y=aft.ML, mode="markers",
    marker=dict(size=4, color=ACCENT, opacity=0.5), name="餘震"))
tt = np.logspace(-3, 1, 100)
fig.add_trace(go.Scatter(
    x=tt, y=float(main.ML) - 4.5 - 0.76 * np.log10(tt), mode="lines",
    line=dict(color=QUAKE_COLOR, dash="dash"),
    name="Mc(t) = Mm − 4.5 − 0.76 log10 t"))
fig.add_hline(y=MC_SEQ, line_dash="dot", line_color=GREEN,
              annotation_text=f"分析門檻 {MC_SEQ}")
fig.add_vline(x=S_START, line_dash="dot", line_color=GREEN,
              annotation_text=f"兩線相交：t = {S_START:.3f} 天")
apply_layout(fig, title="主震後的早期不完整：小餘震要過一段時間才「浮出來」",
             xaxis_title="主震後時間（天，對數軸）", yaxis_title="規模 ML",
             xaxis_type="log", yaxis_range=[2.8, 7.4], hovermode="closest")
fig

# %% [markdown]
# 主震後最初一小時內，目錄裡幾乎只剩 $M\,4$ 以上的餘震——不是小的
# 沒發生，是記不到。紅色虛線與綠色門檻線的交點就是上面算出的 $S$；
# 在它左邊做任何衰減擬合，統計的是儀器不是地震。台灣自己的短期
# 完整度研究（Tsai et al. 2012）用堆疊序列得到同樣的定性結論：
# $M\,7.4$–$7.9$ 的序列起初 $M_c\approx4.0$，隨 $\log t$ 回落到 2.8
# 左右，再降到約 2.4 的平台。係數與南加州版不同，形狀一樣。
#
# $c$ 到底是不是物理量？誠實的答案是：**目前無法從目錄本身分辨**。
# 文獻有一個標準診斷——把分析門檻一路抬高，若 $\hat c$ 隨之縮小，
# 表示原來的 $c$ 主要來自缺漏事件（抬高門檻等於只留「本來就記得到」
# 的那一群），而物理上的成核時間不該隨你選的門檻改變。多數目錄做
# 這個測試都會看到 $\hat c$ 縮小一到兩個數量級；典型報告值因此橫跨
# $10^{-4}$ 到 $10^{-1}$ 天，而且與觀測網密度高度相關。結論可以講得
# 很硬：**把 $c$ 當成模型的正規化參數，不要當物理常數解讀**，除非
# 你能先證明目錄在 $t<c$ 這段期間是完整的。近年的作業型模型（如
# simplETAS）乾脆把 $c$ 釘死在普世值。
#
# ## 12.3 p 值：擬合方法、時間窗，與一條台灣的迴歸線
#
# 12.1 的圖上兩條線來自兩種完全不同的作法。
#
# **分箱最小平方**：切箱、算率、取對數、跑迴歸。直觀而且畫得出來，
# 但有三個系統性問題：(1) 取對數後的最小平方隱含殘差同方差，但每箱
# 計數是 Poisson，尾端的箱只有兩三個事件、對數殘差方差大得多，卻被
# 同等加權；(2) **空箱會被丟掉**（$\log0$ 沒定義），丟掉的一定是率低
# 的箱，等於系統性抬高尾巴、壓低 $p$；(3) 箱寬與箱界全是分析者的選擇。
#
# **點過程最大概似**：直接用第 10 章的骨架，把 $\lambda^*(t)=K(t+c)^{-p}$
# 代入，觀測窗取 $[S,T]$：
#
# $$\ln L(K,c,p) = N\ln K - p\sum_{i=1}^{N}\ln(t_i + c)
#   - K\!\int_S^T\!(t+c)^{-p}\,\mathrm{d}t .$$
#
# 對 $K$ 微分等於零可直接解出 $\hat K=N\big/\!\int_S^T(t+c)^{-p}
# \mathrm{d}t$，代回去把三維最佳化降成二維——這叫**剖面概似**，
# 附錄 12.12 有代數。不分箱、不丟資料、誤差由 Hessian 給。代價是它
# 逼你明確宣告 $[S,T]$ 並假設這段期間完整。$S$ 從哪來？12.2 已經
# 給了答案：由 $M_c(t)$ 反解。這是本章兩節間最實用的一條連結。
#
# 至於 $T$——2022 池上 $M_L\,6.8$ 序列是最乾淨的在地案例
# （Chen et al. 2024）：同一份目錄、同一個程式，取全期 30 天得
# $p=0.92$，取前 12 天得 1.30，只取前 6 天得 **1.39**（完整式子是
# $n(t)=27.7/(t+0.02)^{0.92}$）。為什麼？因為真實序列不是單一個
# Omori 核。窗一旦拉長到涵蓋一個大餘震，那個大餘震自己的次級序列
# 就疊上來，把尾巴撐高，單一核唯一能做的事就是把 $p$ 調小。把 0403
# 花蓮的窗長從 3 天掃到 60 天，這個機制看得清清楚楚：

# %% tags=["hide-input"]
windows = np.unique(np.round(np.logspace(np.log10(3), np.log10(60), 22), 2))
scan = pd.DataFrame(
    [(T, *omori_mle(t_aft, S_START, float(T))[1:4]) for T in windows],
    columns=["T", "c", "p", "sp"])
t_big = np.sort((aft[aft.ML >= 6.0].time - main.time)
                .dt.total_seconds().to_numpy() / 86400)
t_big = t_big[(t_big > 1) & (t_big <= 60)]

fig = go.Figure()
fig.add_trace(go.Scatter(
    x=np.r_[scan["T"], scan["T"][::-1]],
    y=np.r_[scan.p + scan.sp, (scan.p - scan.sp)[::-1]],
    fill="toself", fillcolor="rgba(31,119,180,0.18)", line=dict(width=0),
    hoverinfo="skip", name="± 1σ（Hessian）"))
fig.add_trace(go.Scatter(x=scan["T"], y=scan.p, mode="lines+markers",
                         name="點過程 MLE 的 p̂", line=dict(color=ACCENT),
                         marker=dict(size=6)))
fig.add_hline(y=1.0, line_dash="dash", line_color="#888888",
              annotation_text="教科書的 p = 1")
for tb in t_big:
    fig.add_vline(x=tb, line_dash="dot", line_color=QUAKE_COLOR)
apply_layout(fig, title=f"同一序列、同一程式：p̂ 從 {scan.p.min():.2f} 走到 "
                        f"{scan.p.max():.2f}（起算 S = {S_START:.3f} 天，"
                        f"紅線＝ML ≥ 6.0 的大餘震）",
             xaxis_title="擬合窗長 T（天）", yaxis_title="估計的 p 值",
             hovermode="x")
fig

# %% [markdown]
# 圖形被紅色虛線切成兩段。窗長從 3 天拉到 18 天，$\hat p$ 從 0.9
# 一路爬到 1.58，因為窗愈長、愈看得到真正的長期尾巴。然後在 4 月
# 22–23 日那一串 $M_L\,6.1$–$6.2$ 大餘震進窗的瞬間，$\hat p$ **崩到
# 0.66**——單一 Omori 核描述不了疊了次級序列的活動史，於是它用
# 「衰減得很慢」來解釋多出來的事件。之後隨著次級序列自己衰減下去，
# $\hat p$ 才慢慢爬回 0.8 附近。同一份資料、同一個估計法，$\hat p$
# 的變動範圍是它自己標準誤（0.02–0.14）的**十倍以上**。順帶一提，
# 那幾個崩塌點的 $\hat c$ 同時被推到最佳化的下界（$10^{-4}$ 天）：
# 模型在試圖把自己變成純冪次。**參數貼著邊界，是模型設定不對的
# 警訊，不是發現。**
#
# 有沒有可能 $p$ 其實有系統性結構，只是被這些雜訊蓋住了？
# Tsai et al. (2012) 用 1973–2011 的 CWBSN 目錄（深度 < 35 km，集集
# 改用 $M_w\,7.6$）做了很小心的檢驗：先處理長期與主震後兩層完整度，
# 再用**三種不同的除叢法**堆疊餘震序列，得到
#
# $$p(M_m) \simeq (0.11 \pm 0.01)\,M_m + (0.38 \pm 0.02),
#   \qquad R^2 \approx 0.70\text{–}0.87 .$$
#
# 三種除叢法結果幾乎一致——這是穩健性論證的教科書模板。分區後，
# 東部（隱沒帶、海洋地殼）是 $p\approx0.11M_m+0.36$，西部（造山帶、
# 陸殼）是 $p\approx0.07M_m+0.50$；西部斜率略小，但兩區的 $t$ 檢定
# 都在 0.995 以上拒絕「斜率為零」。代進去看量級：$M_m=5.0$ 給
# $p=0.93$，$M_m=7.0$ 給 $p=1.15$。也就是「$p\approx1$」在台灣只對
# $M_m\approx5.6$ 的主震成立。作者同時很誠實地指出：係數在不同目錄
# 之間**並未**與任何構造特徵明顯相關，所以這條線目前只是現象。
#
# ```{admonition} 報告 p 值的最小清單
# :class: tip
# 一個 $p$ 值單獨出現時沒有意義。至少要同時給出：
# **（1）** 擬合區間 $[S,T]$，並說明 $S$ 怎麼決定的；
# **（2）** 分析門檻 $M_c$，以及它是否高於序列的短期完整度；
# **（3）** 估計法（分箱迴歸／點過程 MLE）與誤差怎麼算的；
# **（4）** 序列怎麼圈出來的（除叢法與參數）；
# **（5）** 主震規模 $M_m$——因為 $p$ 依賴它。
# 少任何一項，別人都沒辦法判斷你的 $p$ 和他的 $p$ 是不是同一件事。
# ```
#
# ## 12.4 產能與餘震區尺度：$\gamma$ 是怎麼來的
#
# Omori 律說「什麼時候」，還缺「幾個」與「在哪裡」。
#
# **幾個**：宇津的產能律說，規模 $m$ 的事件平均觸發的直接後代數隨
# 規模指數成長，本書統一寫成 $\kappa(m)=A\,e^{\alpha(m-m_0)}$。
# $\alpha$ 是「規模換算成觸發力」的效率。注意 $A$ 的意義綁定在
# 「搭配正規化核 {eq}`eq:omori-density`」這個約定上，**不能**跟文獻
# 中 Ogata 式寫法的 $K$ 直接比較（$K$ 吸收了 $c^{p-1}/(p-1)$ 這一團）。
# 這是跨論文比參數最常見的翻車點。
#
# **在哪裡**：Utsu 與 Seki（1955）量了一批日本地震的餘震分布面積，
# 得到 $\log_{10}A_{\rm aft}=m+4.0$；文獻另一個常見的等價寫法是餘震區
# 長度 $\log_{10}L=0.5\,m-1.8$。兩式的**常數項依單位與資料集而異**
# （面積與長度各用什麼單位、餘震區怎麼圈，各家不同），互相換算對不
# 太起來是正常的。真正穩健也真正重要的是**斜率**：面積 $\propto10^{m}$、
# 長度 $\propto10^{0.5m}$。這與規模定義自洽——規模每加一級地震矩約增
# $10^{1.5}$ 倍，若應力降大致固定，破裂長度就該增 $10^{0.5}$ 倍。
#
# 現在做本節要的推導。ETAS 空間核必須寫成「與規模無關的形狀」乘上
# 「隨規模伸縮的尺度」，否則模型不自相似。設尺度為 $L(m)$，由
# Utsu–Seki：
#
# $$\begin{aligned}
# L(m) &= L_0 \cdot 10^{\,0.5(m - m_0)}
#       = L_0 \exp\!\bigl[0.5\ln 10\,(m-m_0)\bigr] \\
#      &\equiv L_0\, e^{\gamma (m - m_0)}
#   \quad\Longrightarrow\quad \gamma = 0.5\ln 10 \approx 1.15 .
# \end{aligned}$$
#
# **$\gamma$ 不是自由發明的參數，它是 Utsu–Seki 律換一個底數寫出來
# 的結果。** 把尺度塞進反冪次核，就是本書統一使用的形式
#
# $$f(x,y;m) = \frac{q-1}{\pi D e^{\gamma(m-m_0)}}
#   \left[1 + \frac{r^2}{D e^{\gamma(m-m_0)}}\right]^{-q},$$
#
# 其中 $r$ 是到親代的震央距，$q>1$ 是遠場衰減冪次（正規化的收斂
# 條件，推法與 12.1 的 $p>1$ 一模一樣）。
#
# 這個值與文獻對得上嗎？Zhuang et al. (2004) 用診斷圖擬合日本目錄，
# 經驗值正是 $\tilde\gamma=0.5\ln10$；Ogata & Zhuang (2006) 對 JMA
# 三個資料集的正式估計落在 0.80–1.33；simplETAS 釘死在 1.5，並註明
# 「資料豐富地區的估計典型落在 1.0–2.0」。同一數量級，但**散得比你
# 以為的寬**——附錄 12.12 說明為什麼（提示：$De^{\gamma(m-m_0)}$
# 到底是長度還是面積，各家寫法不同，差一個因子 2）。
#
# 還有一段歷史值得記：Ogata (1998) 原本用**同一個** $\alpha$ 同時
# 決定「生幾個」與「散多遠」。Zhuang et al. (2004) 用隨機除叢診斷，
# 發現「親代規模對子代距離眾數」的斜率明顯比 $\hat\alpha$ 平緩；
# Ogata & Zhuang (2006) 因此把兩者拆開，參數由七個變八個。**產能與
# 空間尺度是兩件事，不要用一個參數兼差。**
#
# 到這裡，叢集的三個問題都有答案了：$\kappa(m)$ 管幾個、
# {eq}`eq:omori-density` 的 $g(t)$ 管什麼時候、$f(x,y;m)$ 管在哪裡。
# 三者相乘，就是第 13 章條件強度函數裡的那一項觸發核。
#
# ## 12.5 Båth 定律：「最壞會多壞」的三種算法
#
# 防災現場最常被問的問題只有一個：**接下來還會不會有更大的？**
# Båth (1965) 給了一個好記到危險的答案：最大餘震比主震約小
# $\Delta_1 \equiv M_m - M_{\rm max,aft} \approx 1.2$。這一節把三種
# 估計法完整推一次，再用極值論證說明為什麼**它們的平均值對單一
# 序列幾乎都沒有預測力**。
#
# ### 估計一：常數 $\Delta \approx 1.2$
#
# 最簡單，也最常被寫進新聞稿。台灣的系統性檢驗（Chan & Wu 2013）
# 用 706 個序列（$M_L\ge4.0$、深度 < 30 km、1993–2011，以 Wu & Chiao
# 2006 的 5 km／3 天雙連結除叢）算出 $\bar{\Delta}_1=1.20$，與 Båth
# 原值完全一致——但**標準差高達 0.73**。這代表 68% 區間橫跨 1.5 個
# 規模單位，對應能量差約 180 倍。這種區間對疏散決策沒有幫助。
#
# ### 估計二：Utsu 的 $1/\beta$ 下界
#
# 這條界線的推導很漂亮。假設序列裡 $N+1$ 個事件（含被我們叫做主震
# 的那一個）的規模都是同一個 GR 分布的獨立抽樣，即 $M_i-m_0\sim
# \mathrm{Exp}(\beta)$，$\beta=b\ln10$。由大到小排成順序統計量
# $M_{(1)}>M_{(2)}>\dots$，指數分布的**無記憶性**給了一個特別的結構
# （Rényi 表示）：順序統計量之間的間距互相獨立，而且
#
# $$D_k \equiv M_{(k)} - M_{(k+1)} \sim \mathrm{Exp}(k\beta) .$$
#
# 取 $k=1$：最大與次大之差 $D_1\sim\mathrm{Exp}(\beta)$，因此
#
# $$\mathbb{E}[D_1] = \frac{1}{\beta} = \frac{1}{b\ln 10}
#   \approx \frac{0.434}{b},$$
#
# 而且**與序列裡有幾個事件無關**（附錄 12.12 有完整證明）。這就是
# Utsu 的 $1/\beta$。$b=1$ 時它只有 0.43，遠小於 Båth 的 1.2——換句
# 話說，真實主震**比「一群同分布事件裡最大的那個」還要大得多**，
# $1/\beta$ 因此被當成 $\Delta_1$ 的下界。
#
# 這條界線可以直接檢驗。純 i.i.d. 假設下 $P(D_1>1/\beta)=e^{-1}
# \approx36.8\%$；Chan & Wu 在 119 個餘震數 > 50 的台灣序列上得到
# 56% 滿足 $\Delta_1>1/\beta$，再限定主震 $M\ge6.0$ 則跳到
# **86%（12/14）**。也就是：**主震愈大，序列愈不像「一堆同分布
# 事件」，這條下界愈有效。**
#
# ### 估計三：GR 外推的 $m^* = a/b$
#
# 第三種作法不看實測的最大餘震，改用序列自己的 GR 律外推。定義
# 「推論最大餘震」$m^*$ 為期望累積數降到 1 的規模：
#
# $$\log_{10} N(\ge m^*) = a - b\,m^* = 0
#   \quad\Longrightarrow\quad m^* = \frac{a}{b},
#   \qquad \Delta M^* = M_m - \frac{a}{b} .$$
#
# 優點是用了**整個**規模分布而非單一極值樣本，方差小得多；缺點是
# $a$ 值依賴 $M_c$、時間窗與空間範圍——$a$ 不是序列的「性質」，
# 是一個**計數**。
#
# 1999 集集是最著名的案例（Lee et al. 2013）。主震 $M_w\,7.65$，
# 序列 GR 擬合 $\log_{10}N(\ge m)=6.4-0.84m$，給出 $m^*=6.4/0.84
# =7.62$，於是 $\Delta M^*=0.03$；而實測最大餘震是 $M_w\,6.70$，
# $\Delta_1=0.95$。兩個數字差了三十倍。這不是計算錯誤，而是在說
# 一件很具體的事：**集集的小餘震多到 GR 律預期最大餘震應該幾乎跟
# 主震一樣大**（主震後 1000 天內區內記到 42,952 個 $m\ge2.0$ 的
# 餘震）。若要讓加州平均的 $\Delta M^*=1.11$ 成立，主震規模得是
# 8.73 才對得上。三種方法在 Chan & Wu 的四個台灣個案上各有勝負：
# Båth 常數在大埔、南澳、瑞里最準，集集則是 $1/\beta$ 最準；GR 法
# 在 119 個序列上的平均殘差最小（0.13）。沒有一種全勝。
#
# ### 極值論證：為什麼平均值救不了個案
#
# 設 $N$ 個餘震的規模是 i.i.d. 的 $\mathrm{Exp}(\beta)$（相對於
# $M_c$），最大值 $X_{(1)}$ 的分布是
#
# $$P\bigl(X_{(1)} \le x\bigr) = \bigl(1 - e^{-\beta x}\bigr)^{N} .$$
#
# **這不是指數分布**，$N$ 大時它趨近 Gumbel 分布。用間距表示可以
# 直接算出動差：$X_{(1)}=\sum_{k=1}^{N}D_k$，各項獨立且
# $D_k\sim\mathrm{Exp}(k\beta)$，所以
#
# $$\begin{aligned}
# \mathbb{E}\bigl[X_{(1)}\bigr]
#   &= \frac{1}{\beta}\sum_{k=1}^{N}\frac{1}{k} = \frac{H_N}{\beta}
#    \approx \frac{\ln N + 0.5772}{\beta}, \\
# \mathrm{Var}\bigl[X_{(1)}\bigr]
#   &= \frac{1}{\beta^2}\sum_{k=1}^{N}\frac{1}{k^2}
#   \xrightarrow[N\to\infty]{} \frac{\pi^2}{6\beta^2} .
# \end{aligned}$$
#
# 兩條結論直接掉出來。
#
# **第一，$\Delta_1$ 必然隨餘震數下降。** 主震規模固定時，
# $\Delta_1=(M_m-M_c)-H_N/\beta$，隨 $\ln N$ 遞減。從 $N=5$ 到
# $N=50$，$H_{50}-H_5=4.499-2.283=2.216$，除以 $\beta=\ln10$ 得
# **0.96 個規模單位**的預期下降。Chan & Wu 觀測到的是 1.20 → 0.74，
# 下降 0.46——方向一致、量級同級。觀測值較小是因為餘震多的序列
# 主震通常也較大，兩個效應部分抵銷（他們限定主震 $M\ge6.0$ 後
# $\bar\Delta_1$ 回升到 1.26，正是這個混淆的直接證據）。
#
# **第二，$\sigma=0.73$ 幾乎不需要任何物理來解釋。** 極限方差給的
# 標準差是 $\pi/(\sqrt6\,\beta)\approx1.283/\beta\approx0.56$
# （$b=1$）。也就是說，**即使地球完全均一、每個序列的 $b$ 值都一樣**，
# 單一序列的最大餘震規模本身就有 0.56 的標準差；再加上主震規模、
# $M_c$、$b$ 值的變異，0.73 就滿了。「平均 1.2」是一個統計量的位置
# 參數，不是一條可以套在你手上這場序列的預測。把兩件事畫在一起：

# %% tags=["hide-input"]
# 左：台灣長期目錄的實測 ΔM（Gardner–Knopoff 視窗法）
cat_long = load_taiwan_catalog()
sub = cat_long[cat_long.ML >= 3.5].sort_values("time").reset_index(drop=True)
t_num = sub.time.astype("int64").to_numpy() / 86400e9
lat, lon, ml = (sub[c].to_numpy() for c in ("latitude", "longitude", "ML"))


def gk_window(m):
    """Gardner & Knopoff (1974) 的除叢時空窗（公里、天）。"""
    L = 10 ** (0.1238 * m + 0.983)
    T = 10 ** (0.032 * m + 2.7389) if m >= 6.5 else 10 ** (0.5409 * m - 0.547)
    return L, T


is_aft = np.zeros(len(sub), dtype=bool)
bath = []
for i in np.argsort(-ml):
    if ml[i] < 6.0:
        break
    if is_aft[i]:
        continue                       # 自己是更大事件的餘震，跳過
    L, T = gk_window(ml[i])
    sl = slice(i + 1, np.searchsorted(t_num, t_num[i] + T, side="right"))
    dist = np.hypot((lat[sl] - lat[i]) * 111,
                    (lon[sl] - lon[i]) * 111 * np.cos(np.radians(lat[i])))
    idx = np.arange(sl.start, sl.stop)[(dist <= L) & (ml[sl] < ml[i])]
    is_aft[idx] = True
    if len(idx) >= 5:
        bath.append(ml[i] - ml[idx].max())
bath = np.array(bath)

# 右：i.i.d. GR 零假設下 ΔM 隨餘震數的變化
rng = np.random.default_rng(20240403)
BETA, DM_MAX = np.log(10.0), 2.5       # b = 1；主震高出 Mc 2.5 個規模單位
n_grid = np.unique(np.round(np.logspace(0, 3, 22)).astype(int))
sim = np.array([np.percentile(DM_MAX - rng.exponential(
    1 / BETA, size=(3000, int(n))).max(axis=1), [5, 50, 95]) for n in n_grid])
mean_u = np.array([DM_MAX - sum(1 / np.arange(1, n + 1)) / BETA
                   for n in n_grid])
n_zero = float(np.exp(BETA * DM_MAX - 0.5772))

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.11,
                    subplot_titles=(f"台灣 1973–2025 實測（ML ≥ 6 主震，"
                                    f"{len(bath)} 個序列）",
                                    "i.i.d. GR 零假設：ΔM 隨餘震數 N 變化"))
fig.add_trace(go.Histogram(x=bath, xbins=dict(size=0.2), marker_color=ACCENT,
                           opacity=0.85, name="實測 ΔM"), row=1, col=1)
fig.add_vline(x=1.2, line_dash="dash", line_color=QUAKE_COLOR, row=1, col=1)
fig.add_vline(x=float(bath.mean()), line_color=GREEN, row=1, col=1)
fig.add_trace(go.Scatter(
    x=np.r_[n_grid, n_grid[::-1]], y=np.r_[sim[:, 2], sim[::-1, 0]],
    fill="toself", fillcolor="rgba(31,119,180,0.18)", line=dict(width=0),
    hoverinfo="skip", name="5–95% 區間"), row=1, col=2)
fig.add_trace(go.Scatter(x=n_grid, y=mean_u, mode="lines",
                         line=dict(color=ACCENT),
                         name="E[ΔM]（主震規模固定）"), row=1, col=2)
fig.add_trace(go.Scatter(x=n_grid, y=np.full(len(n_grid), 1 / BETA),
                         mode="lines", line=dict(color=GREEN, dash="dot"),
                         name=f"Utsu 下界 1/β = {1/BETA:.2f}"), row=1, col=2)
fig.add_hline(y=1.2, line_dash="dash", line_color=QUAKE_COLOR, row=1, col=2)
fig.add_hline(y=0.0, line_color="#888888", row=1, col=2)
fig.update_xaxes(title_text="ΔM（主震 − 最大餘震）", row=1, col=1)
fig.update_yaxes(title_text="序列數", row=1, col=1)
fig.update_xaxes(title_text="餘震數 N（對數軸）", type="log", row=1, col=2)
fig.update_yaxes(title_text="ΔM", row=1, col=2)
apply_layout(fig, title=f"Båth 定律的兩張臉：實測平均 {bath.mean():.2f}"
                        f"（紅虛線＝Båth 的 1.2），"
                        f"理論在 N ≈ {n_zero:.0f} 時歸零",
             height=460, hovermode="closest")
fig

# %% [markdown]
# 左圖先給一個難堪的結果：視窗法算出的平均 $\Delta M$ 只有 0.3 上下，
# 離 Båth 的 1.2 很遠。這**不是**在推翻 Båth，而是 12.6 節的預告——
# Gardner–Knopoff 對 $M\,6.5$ 以上主震的時間窗長達九百多天、空間窗
# 七十多公里，在台灣這種密集的構造區裡，窗內幾乎一定會有另一場獨立
# 的大地震（雙主震、重疊序列都是常態），於是它被當成「餘震」，
# $\Delta M$ 被系統性壓扁。序列怎麼定義，答案就跟著變；改用 Chan & Wu
# 的 5 km／3 天連結法，同一個地區得到的是 $1.20\pm0.73$。
#
# 右圖是零假設下的理論曲線（$b=1$、主震固定高出 $M_c$ 2.5 個規模
# 單位）。三件事值得盯著看：一、藍線隨 $\ln N$ 下滑，正是
# 「$\Delta_1$ 隨餘震數變小」的來源；二、5–95% 的帶寬**始終超過一個
# 規模單位**，這是不可壓縮的極值變異，不是量測誤差；三、藍線在
# $N\approx178$ 處穿過零——一個主震若有一百多個 $M_c$ 以上的餘震，
# 純 i.i.d. 的 GR 模型就預期「最大餘震該和主震一樣大」。集集有四萬
# 多個餘震，這正是 $\Delta M^*=0.03$ 的來源。而綠色虛線（Utsu 的
# $1/\beta$）是**平的**，因為上面證明了 $D_1$ 的期望與 $N$ 無關。
#
# 第三點還有一層更深的意思，12.7 節會把它兌現：**當 $N$ 大到讓藍線
# 降到零以下，「主震是序列裡最大的」這件事本身就變成一個選擇條件**，
# 而不是資料的性質。
#
# 防災上怎麼講才誠實？Wang et al. (2016) 回顧 30 個台灣 $M>5$ 地震，
# 把主震與最大餘震的四個差值全部對主震規模做迴歸，得到四個**否定
# 結果**：規模差的迴歸是 $\Delta M=-1.01+0.31M_m$（方向甚至與 Båth
# 相反），但作者明講相關性太低、不可用於預測；時間差 $\Delta T$ 與
# 深度差 $\Delta H$ 也與 $M_m$ 無明確相關。唯一有實用價值的是
# **條件式、有界**的那一條：$M_m>6.5$ 時最大餘震的震央距 $\Delta D$
# 幾乎都在 **35 公里以內**（$M_m<6.5$ 時可散到 85 公里）。範圍有界
# 就能劃警戒區——這才是能寫進防災 SOP 的結論。
#
# ## 12.6 除叢的四個家族
#
# 上面每一段計算都偷偷用了一個操作：把目錄拆成「序列」。這件事的
# 正式名字是**除叢**，方法分成四個家族。
#
# **一、視窗法（Gardner–Knopoff 1974）。** 最老、最常用：對每個候選
# 主震依規模開一個時空窗，窗內比它小的事件全部標記為餘震。本章用的
# 版本是
#
# $$\begin{aligned}
# \log_{10} L &= 0.1238\,m + 0.983 \quad (L \text{ in km}), \\
# \log_{10} T &= \begin{cases}
#   0.5409\,m - 0.547, & m < 6.5,\\
#   0.032\,m + 2.7389, & m \ge 6.5,
# \end{cases}\quad (T \text{ in days}).
# \end{aligned}$$
#
# 優點是完全可重現、快、只要規模與位置。缺點在 12.5 的左圖已經看到：
# 窗是規模的**確定性函數**，與這場序列實際的活躍程度無關，而且對大
# 地震開得極寬（$m=7.3$ 時 $L\approx77$ km、$T\approx933$ 天）。台灣
# 的機率式危害度評估傳統上就是用這一家。
#
# **二、連結法（Reasenberg 1985）。** 不開固定窗，改判斷「這個事件
# 和已有的叢集有沒有互動」：互動半徑由主震的應力場尺度給，向前看的
# 時間長度由「下一個事件的等待時間」在某個信心水準下決定，隨叢集內
# 最大事件動態更新。台灣慣用 Wu & Chiao (2006) 的時空雙連結版本，
# 連結參數 5 km／3 天。通病是**串接（chaining）**：$A$ 連 $B$、
# $B$ 連 $C$，於是 $A$ 與 $C$ 被歸為同一叢，即使相隔很遠。在活躍區
# 的長目錄上，單純的單連結可以一路串下去，把幾十年的活動縫成一個
# 超級叢集；實作時必須加限制，否則會退化。
#
# **三、最近鄰距離（Zaliapin & Ben-Zion）。** 把時間、空間、規模三者
# 合成一個標量。對事件對 $(i,j)$ 定義
#
# $$\eta_{ij} = t_{ij}\; r_{ij}^{\,d}\; 10^{-b\,m_i},
#   \qquad t_{ij} = t_j - t_i > 0 ,$$
#
# 而 $t_{ij}\le0$ 時定義 $\eta_{ij}=\infty$。$r_{ij}$ 是震央距，$d$ 是
# 震央分布的碎形維度（常取 1.5–1.6），$b$ 是 GR 斜率。每個事件 $j$
# 只保留 $\eta_j=\min_{i<j}\eta_{ij}$ 與對應的「親代」$i$。邏輯是：
# $10^{-b m_i}$ 依 GR 律把「等一個 $m_i$ 級事件平均要多久」折算進去，
# 於是大事件後面的短間隔就不算稀奇；剩下的 $t\cdot r^d$ 是時空體積。
# **背景事件的 $\eta$ 大、被觸發事件的 $\eta$ 小**，$\log\eta$ 的直方
# 圖因此呈雙峰，谷底就是切點 $\eta_0$。優點：只有兩個參數且都能從
# 資料估、不需要「主震」概念、叢集是一棵樹（每個事件只有一個親代），
# 天生不串接。缺點：$\eta_0$ 隨 $M_c$、$d$、$b$ 移動，完整度差的目錄
# 上兩個峰會糊在一起，切點變成主觀選擇。
#
# **四、隨機除叢（ETAS）。** 前三家都在回答「這個事件屬於哪一叢」。
# 第四家換一個問法：**「這個事件是背景事件的機率有多大？」** 擬合好
# 的 ETAS 可以算出「$j$ 被 $i$ 觸發」的機率 $\rho_{ij}$ 與「$j$ 是
# 背景」的機率 $\phi_j=1-\sum_i\rho_{ij}$。不做二分，只給機率，而且
# 能用重抽樣把除叢本身的不確定性顯示出來。這是第 14 章的主題，這裡
# 先記住它存在——因為 12.7 節要用它當對照組。
#
# ## 12.7 除叢的選擇效應：一個可以算出來的破產點
#
# 除叢後的目錄有個大家都喜歡的性質：背景事件近似 Poisson 過程，年率
# 穩定，可以直接餵給危害度積分。既然這麼好用，問題出在哪裡？
#
# Mizrahi et al. (2021) 拿 1980 年起的加州目錄，套上**五類**常用除叢法
# （Reasenberg、三種視窗法、Zaliapin 最近鄰、兩種 ETAS 隨機除叢），
# 比較除叢前後的規模分布。結果是：**除叢後的「主震」b 值最多下降
# 30%**，各法散布在 0.73–1.00 之間（全目錄是 1.01），中間沒有明顯
# 間隙；不同方法留下的主震數相差 **6.1 倍**。
#
# 先看這個下降的邏輯後果。若全目錄與主震目錄都被硬套上 GR 直線：
#
# $$\begin{aligned}
# \log_{10} N(m) &= a - b\,m, \\
# \log_{10} N_{\rm main}(m) &= a_{\rm main} - b_{\rm main}\,m ,
# \end{aligned}$$
#
# 當 $b_{\rm main}\neq b$ 兩條線斜率不同，必定在某個規模相交。令兩式
# 相等：$a-b\,m=a_{\rm main}-b_{\rm main}\,m$，即 $a-a_{\rm main}
# =(b-b_{\rm main})\,m$，於是
#
# $$
# m_x = \frac{a - a_{\rm main}}{b - b_{\rm main}}
# $$ (eq:mx)
#
# 再看主震占比。由兩式相減，
#
# $$r(m) = \frac{N_{\rm main}(m)}{N(m)}
#   = 10^{\,(a_{\rm main}-a) + (b - b_{\rm main})\,m} .$$
#
# 當 $b_{\rm main}<b$（實際觀測到的方向），指數隨 $m$ **單調遞增**：
# $r(m)$ 在 $m=m_x$ 剛好等於 1，再往上就**大於 1**——模型預測
# 「主震數多於總地震數」，在觀測上不可能。$m_x$ 因此是一個可以用
# 紙筆算出的**破產點**：低於它，危害度被低估（餘震被丟掉了）；高於
# 它，危害度被高估（$b$ 被壓低，大地震的相對頻率被抬高）。加州資料
# 算出的 $m_x$ 落在 **6.9 到 8.8** 之間——正好是工程上最關心的規模
# 區間。畫出來是這樣：

# %% tags=["hide-input"]
a_full, b_full = 5.60, 1.01              # 示意值：全目錄
b_main, mc_ref, frac_at_mc = 0.78, 3.6, 0.13   # 除叢後主震；mc 處主震占比
a_main = np.log10(frac_at_mc) + a_full - (b_full - b_main) * mc_ref
m_x = (a_full - a_main) / (b_full - b_main)

mgrid = np.linspace(3.0, 9.5, 200)
fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                    subplot_titles=("兩條 GR 直線必定相交",
                                    "主震占比 r(m) 在 m_x 之上超過 1"))
fig.add_trace(go.Scatter(x=mgrid, y=a_full - b_full * mgrid, mode="lines",
                         line=dict(color=ACCENT),
                         name=f"全目錄：b = {b_full:.2f}"), row=1, col=1)
fig.add_trace(go.Scatter(x=mgrid, y=a_main - b_main * mgrid, mode="lines",
                         line=dict(color=QUAKE_COLOR),
                         name=f"除叢後主震：b = {b_main:.2f}"), row=1, col=1)
fig.add_trace(go.Scatter(
    x=mgrid, y=10 ** ((a_main - a_full) + (b_full - b_main) * mgrid),
    mode="lines", line=dict(color=GREEN), name="主震占比 r(m)"), row=1, col=2)
fig.add_hline(y=1.0, line_dash="dash", line_color=QUAKE_COLOR, row=1, col=2)
for col in (1, 2):
    fig.add_vline(x=m_x, line_dash="dash", line_color=GREEN, row=1, col=col)
fig.update_xaxes(title_text="規模 m", row=1, col=1)
fig.update_xaxes(title_text="規模 m", row=1, col=2)
fig.update_yaxes(title_text="log10 N(≥ m)", row=1, col=1)
fig.update_yaxes(title_text="r(m) = N_main / N", range=[0, 2.2], row=1, col=2)
apply_layout(fig, title=f"除叢的破產點（合成示意）：m_x = ({a_full:.2f} − "
                        f"{a_main:.2f}) / ({b_full:.2f} − {b_main:.2f}) "
                        f"= {m_x:.2f}",
             height=430, hovermode="x")
fig

# %% [markdown]
# 這張圖是**合成示意**，$a$、$b$ 取自加州量級的代表值；重點不在數字
# 而在幾何：只要 $b_{\rm main}<b$，兩條線一定會交，而交點右邊的區域，
# 模型在說一件不可能的事。
#
# 機制其實很單純，而且與 12.5 的極值論證是同一件事的兩面。**多數
# 除叢演算法把「叢集中最大的事件」定義為主震。** 一個小地震要成為
# 某一叢中最大的那個，本來就比較難；因此小地震被移除的比例高於大
# 地震，剩下的相對頻率往大規模端傾斜，斜率自然變小。這是**選擇
# 效應**，不是物理。數學上的根本理由 12.5 已經證過：一組 i.i.d. 指數
# 變數的**最大值**的分布不是指數分布（它趨近 Gumbel）。所以「主震
# 規模服從 GR 律」本來就不自洽——Lombardi (2003) 與 Zhuang & Ogata
# (2006) 早指出，全目錄的 GR 律只在 $m\to\infty$ 的漸近意義下對主震
# 成立。
#
# 這個說法要怎麼證明，而不只是講得通？Mizrahi et al. 設計了一個非常
# 乾淨的對照。第一步：用 ETAS 模擬 2000 份**合成目錄**，模擬時所有
# 事件（不分主震餘震）的規模都抽自**同一個** GR 分布，$b$ 值已知。
# 若「主震 b 值較低」是物理性質，合成資料上不該出現任何下降。結果
# **照樣下降**，幅度甚至比真實資料還大。第二步——這才是關鍵——他們
# 刻意做了兩種 ETAS 除叢：
#
# - **ETAS-main**：把叢集中最大的事件定義為主震（沿用傳統定義，好與
#   其他四類方法比較）→ b 值明顯下降；
# - **ETAS-background**：把「未被觸發的事件」定義為主震（這才是 ETAS
#   原生的語意，背景事件可以是任何規模）→ **b 值與全目錄無顯著
#   差異**，合成與真實資料皆然。
#
# 同一個模型、同一份目錄、同一套機率，只換「主震」的定義，一個下降
# 30%、一個不動。這組對照把因果講得再清楚不過：**問題不在「除叢」，
# 而在「主震 ＝ 叢集中最大事件」這個定義。** ETAS 世界裡的餘震可以
# 比它的親代更大（唯一要求是「發生在後」），這與傳統主震定義根本
# 不相容。至於為什麼合成資料的效應**比真實資料更強**？因為所有除叢
# 法都假設餘震空間分布等向——這對合成目錄成立、對真實地震不成立，
# 所以合成目錄的叢集比較好抓，小事件被剔得更乾淨。
#
# ## 12.8 三重後果，與 PSHA 的接口
#
# 把上面兩節合起來，Mizrahi et al. 對機率式地震危害度評估（PSHA）
# 提出三條指控：
#
# 1. **主震的定義本身不可驗證。** 沒有任何觀測可以判定「這個事件是不
#    是主震」——它是演算法的輸出，不是地震的屬性。不同演算法留下的
#    主震數差 6.1 倍，就是這句話的量化版本。
# 2. **丟掉餘震會低估危害度。** 餘震一樣致災：集集的嘉義餘震、2024
#    花蓮 4 月 22–23 日那一串 $M\,6$ 級事件，對已經受損的結構都是實質
#    威脅。把它們從目錄裡刪掉，等於假設它們不存在。
# 3. **被壓低的 b 值會在 $m>m_x$ 之上高估危害度**，這是 {eq}`eq:mx`
#    的直接後果。
#
# 最容易被拿來辯護的是「反正一個低估、一個高估，剛好抵銷」。原作者
# 不接受，理由很硬：**兩個錯誤只在 $m=m_x$ 這一個規模上剛好抵銷**，
# 其他每個規模都是錯的，而且錯的方向隨規模翻轉。除叢仍有它成立的
# 理由——PSHA 的 Poisson 假設要求事件獨立，而叢集顯然不獨立。真正
# 的問題是：**除叢是為了滿足模型的假設，不是為了描述地球。** 記住
# 這句話，就會知道除叢後的目錄能拿來做什麼（估背景率）、不能拿來
# 做什麼（估 b 值、估短期危害）。看一次除叢前後的累積事件數：

# %% tags=["hide-input"]
kept = sub[~is_aft]
fig = go.Figure()
for data, name, color in [
        (sub, f"原始目錄（ML ≥ 3.5，{len(sub)} 筆）", ACCENT),
        (kept, f"除叢後（{len(kept)} 筆，約 "
               f"{100 * len(kept) / len(sub):.0f}%）", GREEN)]:
    fig.add_trace(go.Scattergl(x=data.time, y=np.arange(1, len(data) + 1),
                               mode="lines", name=name,
                               line=dict(color=color)))
for day, label in [("1999-09-21", "集集"), ("2024-04-03", "0403 花蓮")]:
    fig.add_vline(x=day, line_dash="dot", line_color=QUAKE_COLOR,
                  annotation_text=label)
apply_layout(fig, title="除叢前後的累積事件數（Gardner–Knopoff 視窗法，"
                        "僅以 ML ≥ 6 為主震）",
             xaxis_title="時間", yaxis_title="累積事件數", hovermode="x")
fig

# %% [markdown]
# 原始目錄（藍）在集集與 0403 花蓮處有明顯的階梯跳升——那是餘震
# 序列；除叢後（綠）曲線平滑許多，接近等速率的直線，正是 Poisson
# 假設想要的樣子。看起來很成功。但現在你知道要問下一個問題了：
# **那條綠線背後的規模分布，斜率被改成什麼樣子了？**
#
# 台灣的接口很直接：本地的 PSHA 同樣慣用 Gardner–Knopoff 視窗法除叢，
# 所以上面每一條指控都適用。第 21 章討論危害度積分時會回到這裡，
# 並看看「直接用 ETAS 模擬長目錄、繞開除叢與 Poisson 假設」這條替代
# 路線能走多遠。
#
# ## 12.9 參數與典型值
#
# 先是本章三條經驗律的參數本身：
#
# | 參數 | 意義 | 典型值 | 註記 |
# |---|---|---|---|
# | $p$ | Omori 衰減指數 | 0.9–1.4 | 正規化需 $p>1$；隨窗長劇烈變動 |
# | $c$ | Omori 平緩期 | $10^{-4}$–$10^{-1}$ 天 | 多半吸收了早期不完整 |
# | $\alpha$ | 產能指數 | 0.35–3.1 | 群震型小、主震–餘震型大 |
# | $\gamma$ | 空間尺度指數 | 0.8–1.5 | $=0.5\ln10$ 由 Utsu–Seki 導出 |
# | $q$ | 空間遠場冪次 | 1.5–1.8 | 正規化需 $q>1$ |
#
# 台灣的 $p$ 值有好幾個版本，差異全部來自分析選擇：
#
# | 序列／來源 | $p$ | 條件 | 出處 |
# |---|---|---|---|
# | 1999 集集 | 1.05 | $c_0=6\times10^{-3}$ 天 | Lee et al. 2013 |
# | 2022 池上（6 天） | 1.39 | 同一份資料 | Chen et al. 2024 |
# | 2022 池上（12 天） | 1.30 | 同一份資料 | Chen et al. 2024 |
# | 2022 池上（30 天） | 0.92 | $c=0.02$ 天 | Chen et al. 2024 |
# | $p(M_m)$ 全台 | $0.11M_m+0.38$ | 三種除叢法一致 | Tsai et al. 2012 |
# | $p(M_m)$ 東部 | $0.11M_m+0.36$ | 隱沒／海洋地殼 | Tsai et al. 2012 |
# | $p(M_m)$ 西部 | $0.07M_m+0.50$ | 造山／陸殼 | Tsai et al. 2012 |
# | ETAS 全台 | 1.04–1.06 | 兩份獨立估計 | 見第 13 章 |
#
# Båth 定律的各種分組，混亂本身就是重點：
#
# | 分組 | $\bar{\Delta}_1$ | $\sigma$ | 序列數 | 出處 |
# |---|---|---|---|---|
# | 全部序列 | 1.20 | 0.73 | 706 | Chan & Wu 2013 |
# | 餘震數 > 50 | 0.74 | 0.52 | 119 | 同上 |
# | 主震 $M\ge6.0$ | 1.26 | 0.43 | 14 | 同上 |
# | 集集（實測） | 0.95 | — | 1 | Lee et al. 2013 |
# | 集集（GR 外推） | 0.03 | — | 1 | Lee et al. 2013 |
#
# ## 12.10 常見誤解與陷阱
#
# **一、「$p\approx1$ 是普世常數」。** Tsai et al. (2012) 在台灣量到
# $p$ 隨主震規模線性上升，三種除叢法都成立；Chen et al. (2024) 在同
# 一個序列上量到 $p$ 隨時間窗從 1.39 掉到 0.92。$p\approx1$ 是很好的
# **起始猜測**，不是常數。
#
# **二、「$c$ 是物理常數」。** 12.2 證明了主震後的時變完整度會在觀測
# 率上造出一段偽平緩期，形狀與 $c$ 的效果無法區分。
#
# **三、報告 $p$ 值卻不報時間窗。** 這等於沒有報告。不報 $M_c$、
# 不報 $S$、不報除叢法也一樣。
#
# **四、「Båth 說 1.2，所以最大餘震會是 $M_m-1.2$」。** $\sigma=0.73$
# 代表 68% 區間橫跨 1.5 個規模單位；12.5 證明其中至少 0.56 是**不可
# 壓縮的極值變異**，就算地球完全均一也消不掉。
#
# **五、「主震的 b 值天生比較低」。** ETAS-main／ETAS-background 的
# 對照顯示，下降只在套用「主震 ＝ 叢中最大」這條規則時出現。
#
# **六、「除叢後的目錄是乾淨的，可以直接拿去做各種統計」。** 除叢的
# 設計目的只有一個：讓事件近似獨立以滿足 Poisson 假設。它**不保證**
# 保留規模分布，實際上還會系統性破壞它。
#
# **七、「$\gamma$（或 $c$、$D$、$K$）可以跨論文直接比較」。** 不行。
# $\gamma$ 要看它乘的是長度還是面積（差一個因子 2）；$c$ 與 $D$ 要看
# $M_c$ 取多少；$K$ 要看核有沒有正規化。更麻煩的是 ETAS 參數彼此高度
# 相關（$\gamma$ 與 $D$ 在等向核裡只有乘積被辨識），**一定要連標準誤
# 一起看**。
#
# **八、「擬合出 $p<1$，代表這個序列衰減特別慢」。** 通常不是。更常
# 見的原因是模型設定錯了：把空間不均勻的背景活動、或一個疊上來的次級
# 序列，硬塞進單一 Omori 核。12.3 的掃描圖就是活例——順帶提醒，那裡
# 的 $\hat c$ 貼在下界上時最佳化一樣「成功」，**參數貼邊界是模型設定
# 有問題的訊號**。
#
# ## 12.11 研究前沿與未解問題
#
# **繞開除叢。** 最乾脆的一條路：既然「主震」的定義不可驗證，就不要
# 用它。直接用 ETAS 對**全目錄**做率估計，只需要全目錄的 GR 律，不必
# 對「被任意挑出來的大事件」假設任何規模分布；再用數十萬次情境模擬
# 涵蓋時空叢集。Nandan et al. (2019) 在加州的前瞻式實驗中，ETAS 全面
# 勝過平滑地震活動度與應變率模型。台灣的 PSHA 目前仍是視窗法除叢，
# 這是一個明顯可做的在地化題目。
#
# **最近鄰方法的現代化。** Zaliapin 的 $\eta$ 把三個維度壓成一個標量，
# 代價是切點 $\eta_0$ 的主觀性。近年方向包括：把雙峰分布正式擬合成
# 混合模型而不是目視找谷底、讓 $d$ 與 $b$ 隨區域變動、把 $\eta$ 的
# 兩個分量（重整化的時間與空間）分開，用二維而非一維的切割。
#
# **不需要除叢的率估計。** 第 11 章的 b-positive 開了一條有趣的路：
# 用相鄰事件的**規模差**取代規模本身，繞開時變完整度。同樣的想法能
# 不能用到率估計上（「a-positive」）？也就是找一個對目錄缺漏不敏感
# 的統計量來估地震活動度水準。這個方向還很新。
#
# **$c$ 值論戰。** 需要的是同時具備高密度觀測網與可靠早期目錄的資料。
# 日本 Hi-net、加州的模板比對（template matching）補完目錄，都讓
# $\hat c$ 一路縮小；但縮到什麼程度會停，或者根本不會停，目前沒有
# 定論。台灣有全球密度最高的觀測網之一，很有條件回答這個問題。
#
# **$p(M_m)$ 的機制。** 多碎形應力活化模型預測了這條線，四套獨立目錄
# 也都測到了，但係數在不同目錄之間與構造特徵沒有明顯相關。現象穩固、
# 機制不明——這是統計地震學相當典型的狀態。
#
# ## 12.12 附錄：本章推導細節
#
# ### A. $p\le1$ 時的截斷正規化
#
# 若要在 $p\le1$ 時仍使用機率密度形式，必須先選定上界 $T$：
#
# $$g_T(t) = \frac{(t+c)^{-p}}
#   {\displaystyle\int_0^{T}(s+c)^{-p}\,\mathrm{d}s}
#   = \frac{(1-p)\,(t+c)^{-p}}{(T+c)^{1-p} - c^{1-p}},
#   \qquad 0 \le t \le T ,$$
#
# $p=1$ 時分母換成 $\ln[(T+c)/c]$。**$g_T$ 依賴 $T$，所以由它導出的
# 分支比（第 13 章）也依賴 $T$**——這是為什麼 $p<1$ 的擬合結果不能
# 直接餵進標準 ETAS 的分支比公式。
#
# ### B. $\gamma$ 的長度／面積慣例
#
# 本書的空間核把 $De^{\gamma(m-m_0)}$ 放在 $r^2$ 的位置上，嚴格說它的
# 量綱是**面積**。若堅持這個讀法，由 $A_{\rm aft}\propto10^{m}$ 直接
# 比對指數會得到 $\gamma=\ln10\approx2.30$；而 12.4 的推導把
# $De^{\gamma(m-m_0)}$ 讀成**長度**尺度，得到 $0.5\ln10\approx1.15$。
# 兩者剛好差一個因子 2。文獻的實測值（JMA 三個資料集 0.80–1.33、
# Zhuang et al. 2004 的 $0.5\ln10$、simplETAS 釘死的 1.5、「典型
# 1.0–2.0」）落在這兩個理論值之間，正好覆蓋這個因子 2 的區間。這不是
# 巧合，而是各家寫法不一致的直接後果。**實務守則：看到別人報的
# $\gamma$，先確認它乘的是 $r$ 還是 $r^2$，再決定能不能跟你的比。**
# 同樣的檢查也適用於 $D$（度$^2$ 還是 km$^2$）。
#
# ### C. Rényi 表示與最大值的動差
#
# 設 $X_1,\dots,X_N$ 獨立同分布於 $\mathrm{Exp}(\beta)$，順序統計量為
# $X_{(1)}>X_{(2)}>\dots>X_{(N)}$。由無記憶性可證：間距
# $D_k=X_{(k)}-X_{(k+1)}$（$k=1,\dots,N-1$）以及 $D_N=X_{(N)}$ 彼此
# **獨立**，且 $D_k\sim\mathrm{Exp}(k\beta)$。直觀理由：把 $N$ 個變數
# 想成 $N$ 個獨立的「壽命」，最小值是 $N$ 個競爭者的第一個到期時間，
# 服從 $\mathrm{Exp}(N\beta)$；由無記憶性，扣掉這一段之後剩下的 $N-1$
# 個仍是獨立的 $\mathrm{Exp}(\beta)$，遞迴即得。由此立刻得到本章用到
# 的三個結果：
#
# $$\begin{aligned}
# \mathbb{E}\bigl[X_{(1)} - X_{(2)}\bigr] &= \frac{1}{\beta}, \\
# \mathbb{E}\bigl[X_{(1)}\bigr]
#   &= \sum_{k=1}^{N}\frac{1}{k\beta} = \frac{H_N}{\beta}, \\
# \mathrm{Var}\bigl[X_{(1)}\bigr]
#   &= \sum_{k=1}^{N}\frac{1}{(k\beta)^2}
#   \xrightarrow[N\to\infty]{} \frac{\pi^2}{6\beta^2} .
# \end{aligned}$$
#
# 第一式是 Utsu 的 $1/\beta$，且**與 $N$ 無關**；第二式給出 $\Delta_1$
# 隨 $\ln N$ 遞減；第三式給出 $\sigma\to\pi/(\sqrt6\beta)\approx
# 1.283/\beta$。三者在 $b=1$ 時分別是 0.43、$(\ln N+0.58)/2.303$
# 與 0.56。
#
# ### D. Omori MLE 的剖面概似
#
# 由 12.3 的對數概似對 $K$ 微分：
#
# $$\begin{aligned}
# \frac{\partial \ln L}{\partial K}
#   &= \frac{N}{K} - \int_S^T (t+c)^{-p}\,\mathrm{d}t = 0
# \quad\Longrightarrow\quad
# \hat{K}(c,p) = \frac{N}{I(c,p)} , \\
# \ln L_{\rm prof}(c,p) &= N\ln\frac{N}{I(c,p)}
#   - p\sum_{i=1}^{N}\ln(t_i + c) - N ,
# \end{aligned}$$
#
# 其中 $I(c,p)=\bigl[(T+c)^{1-p}-(S+c)^{1-p}\bigr]/(1-p)$（$p=1$ 時為
# $\ln[(T+c)/(S+c)]$）。三維最佳化因此降成二維，收斂穩定得多。$p$ 的
# 標準誤取觀測資訊矩陣（$-\partial^2\ln L_{\rm prof}$ 的數值 Hessian）
# 反矩陣對角元素開根號——本章程式就是這樣算的。注意這是**漸近**標準
# 誤，在 $\hat c$ 貼邊界時完全不可信。
#
# ---
#
# 這一章把叢集的三條老經驗律推到了它們該有的精確度，也把「除叢」這個
# 看似無害的前處理拆開來看了一次。回頭看，每一節的結論都指向同一個
# 方向：**只要你先在資料上劃一條「主震／餘震」的線，後面所有統計量都
# 會帶著那條線的指紋。** $p$ 帶著時間窗的指紋、$c$ 帶著完整度的指紋、
# $\Delta M$ 帶著序列定義的指紋、$b$ 帶著除叢法的指紋。
#
# 出路只有一條：**不要劃那條線。** 與其先做一個不可驗證的分類再統計，
# 不如用一個同時描述背景與觸發的模型，讓「這是不是餘震」變成模型輸出
# 的一個機率，而不是輸入的一個假設。{doc}`第 13 章 <13_etas_structure>`
# 要做的正是這件事：把本章的 $\kappa(m)$、$g(t)$、$f(x,y;m)$ 與第 11 章
# 的 GR 律縫成一條條件強度函數，讓「地震觸發地震」變成一個可以擬合、
# 可以模擬、可以檢驗的數學物件。屆時你會發現，本章辛苦推導的每一條
# 式子都在那條函數裡佔一個位置——包括那個一直在製造麻煩的 $c$。
