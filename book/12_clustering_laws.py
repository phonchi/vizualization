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
# # 6. 目錄統計 II：叢集的經驗律與除叢
#
# 前一章把注意力放在規模軸：小地震有多少、哪些事件可能漏記，以及如何估計大小事件的比例。現在將同一批地震放回時間與空間。即使兩份目錄的規模分布相同，事件也可能一份分散、一份成群；只看 $b$ 值，無法分辨這個差異。
#
# 最容易觀察的叢集是大震後的餘震序列。活動通常先密集，接著逐漸減少，卻未必平順下降，中間仍可能再出現較大的事件。要描述這段過程，我們依次問：發生率怎麼變、較大的事件會帶來多少活動、活動分布在哪裡，以及最大餘震能知道多少。
#
# 這些問題產生了幾條廣泛使用的經驗律。它們用少數量化關係整理現象，也構成下一個預報模型的基礎；適用範圍與觀測條件仍須逐一確認。
#
# ## 6.1 餘震率為什麼先高、後低？
#
# Omori–Utsu 律用一條緩慢衰減的曲線描述主震後時間 $t$ 的發生率：
#
# $$n(t)=K(t+c)^{-p},\qquad t\ge0,\quad K>0,\ c>0.$$
#
# 這裡 $n(t)$ 表示餘震率，與後面 ETAS 的純量分支比 $n$ 不同。$K$ 調整整體活動量，$c$ 控制最初的轉折，$p$ 控制較晚時間的衰減。若時間用天，$c$ 也用天，而 $K$ 的單位會隨 $p$ 改變。因此比較兩份文獻的 $K$ 前，還要確認時間單位與參數寫法。
#
# 曲線下方某段面積，代表那段時間內的期望餘震數。只要 $c>0$、觀測時間有限，這個面積對有限的 $p$ 都是有限的；**有限窗的 Omori 擬合不要求 $p>1$**。$p\le1$ 可以描述有限期間的緩慢衰減，並不自動表示擬合無效。
#
# 下一步建立 ETAS 時，我們想把「預期產生幾個直接後代」和「這些後代何時出現」分開。若時間密度要定義在整個 $[0,\infty)$，便需要 $p>1$，此時
#
# $$g(t)=\frac{p-1}{c}\left(1+\frac{t}{c}\right)^{-p},\qquad t\ge0.$$ (eq:omori-density)
#
# 它的總面積是 1。這個條件限制的是無限時間正規化，不是說真實序列一定有 $p>1$。若使用有限觸發期限，則須改用相應的截斷核與產能定義；附錄 B 列出兩種積分。
#
# 下面用花蓮序列展示同一份資料的兩種擬合。散點是分箱後的率，曲線分別由分箱迴歸與點過程概似得到；估計起點採用外地短期完整度經驗式提供的教學設定，並非對此序列完整性的獨立認證。先看整體衰減，再留意曲線未能描述的局部活動。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["remove-input"]
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
        at_boundary = any(abs(x[i] - edge) < 1e-5
                          for i, pair in enumerate(bounds) for edge in pair)
        valid = (best.success and not at_boundary
                 and np.all(np.linalg.eigvalsh(hess) > 0))
        sig_p = float(np.sqrt(np.linalg.inv(hess)[1, 1])) if valid else np.nan
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
# 兩條曲線使用的資料表示方式和擬合區間不同，得到的參數也可能不同。分箱圖很適合看趨勢，但不能把曲線差異完全歸因於演算法優劣。單一衰減曲線若遇到後續較大事件，也可能需要拉平尾端來描述新增的活動。
#
# ## 6.2 早期平緩可能來自觀測
#
# 還記得前一章的短期不完整嗎？主震後波形重疊時，小事件較難辨認；觀測到的餘震率是實際發生率經過偵測過程篩選的結果。因此，一段先平後降的曲線，可能包含真實時間機制，也可能部分由漏測造成。
#
# Helmstetter et al.（2005）在南加州研究中使用隨時間下降的完整度關係。下圖將其形式 $M_c(t)=M_m-4.5-0.76\log_{10}t$ 疊在花蓮目錄上，$t$ 以天計。它在這裡是跨區域套用的參考曲線，不能因為外觀相近就當成臺灣的偵測模型。
#
# 這個例子要建立的直覺是：若最初只能辨識較大事件，觀測率會被壓低；隨偵測能力恢復，小事件逐漸加入，會減弱圖上看見的衰減。$c$ 和 $p$ 都可能吸收這個觀測效應。若要把它們解讀為物理量，需要更多觀測證據。
#
# %% tags=["remove-input"]
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
# 散點與參考曲線的相對位置，有助於提出「哪一段可能不完整」的問題，卻不能證明曲線下方每個沒記到的事件都確實發生過。實際分析可比較不同門檻、不同起算時間，也可藉助更完整的波形重偵測目錄。
#
# ```{admonition} 漏測如何改變看見的斜率？
# :class: dropdown
#
# 在 $M_c(t)>m_0$ 的期間，若真實規模符合 GR，理想的硬門檻偵測比例為 $e^{-\beta[M_c(t)-m_0]}$。將對數時間形式代入，比例與 $t^{0.76b}$ 成正比；乘在 $t^{-p}$ 上，便得到表觀斜率 $p-0.76b$。當 $M_c(t)\le m_0$，偵測比例應封頂為 1。這是理想化示範，不代表漏測能被單一 $c$ 精確吸收。
# ```
#
# ## 6.3 擬合時間窗也是問題的一部分
#
# 分箱後取對數很直觀，卻會遇到空箱無法取對數、低計數箱的變異較大，以及箱界選擇的影響。點過程概似則直接利用事件時間，並把沒有事件的時段透過積分項一起計算。這樣就保留了「沒有事件」的資訊。
#
# 不過，不分箱並不等於沒有假設。仍須選定 $[S,T]$、規模門檻及研究區域，並確認這段資料能以指定的率模型描述。給定 $c,p$，最佳的 $K$ 會讓模型在窗內的積分符合事件數；真正要辨認的是發生時間如何分布。這個剖面概似的計算放在附錄 B。
#
# 下圖固定起點、逐步延長終點，再重新估計參數。閱讀時把參數變化和新進入時間窗的事件對照。這不是同時獨立觀測到許多個 $p$；每個估計都使用大量重疊資料，目的是檢查對分析範圍的敏感性。
#
# %% tags=["remove-input"]
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
# 若新加入的時間段包含另一串活動，單一衰減曲線的估計會跟著改變。參數貼近數值邊界也值得進一步檢查，但可能的原因包括資料資訊不足、參數限制與模型錯設，不能只用一種解釋。
#
# 模型還需要考慮後續的觸發：除了最早那個事件，後續事件也可能成為新的觸發源。ETAS 會利用這個結構處理重疊序列；在開始疊加之前，先補上事件數量與空間範圍。
#
# ## 6.4 大事件影響多少事件、分布多遠？
#
# 一個較大的事件，通常與較多餘震及較廣的活動區域相關。這是兩個問題，不能只用「影響比較大」混在一起。產能律描述平均直接後代數：
#
# $$\kappa(m)=A e^{\alpha(m-m_0)}.$$
#
# $A$ 是門檻規模事件的平均直接後代數，$\alpha$ 則控制親代規模增加後的產能加成。這個解釋須搭配面積為 1 的時間與空間核；若核未正規化，前面的係數就不是同一個 $A$。
#
# 空間分布可用一個隨規模伸縮、遠處逐漸衰減的密度表示。本書採用
#
# $$f(x,y;m)=\frac{q-1}{\pi\sigma(m)}\left[1+\frac{r^2}{\sigma(m)}\right]^{-q},\qquad \sigma(m)=D e^{\gamma(m-m_0)},\quad q>1.$$
#
# 其中 $r^2=x^2+y^2$，所以 $\sigma$ 與 $D$ 是**距離平方的尺度**。$\gamma$ 控制活動區域如何隨規模擴張，$q$ 控制遠端尾巴。Utsu–Seki 型關係提供「大事件活動區域較大」的經驗動機，但不直接決定所有空間核的參數。
#
# 這裡有一個容易忽略的單位問題：若餘震區面積確實與 $10^m$ 成正比，對應上述面積尺度會得到 $\gamma=\ln10$；若改用長度尺度 $L\propto10^{m/2}$，長度指數才是 $\ln10/2$。實際擬合可以偏離這些簡化關係，不能把兩個指數直接互換。Ogata 與 Zhuang（2006）將產能和空間尺度分開估計，正是為了讓資料分辨「產生幾個」與「分布多遠」。
#
# ## 6.5 最大餘震為什麼仍很不確定？
#
# 如果一次序列已有許多事件，人們自然會問是否還可能出現更大的。Båth 定律以許多序列的平均規模差概括這個問題：主震與最大餘震通常相差約 1.2。但這是受序列定義與選樣影響的經驗摘要，不是對任何一次地震都成立的上限。
#
# 先考慮一個更簡單的統計情境：固定有 $N$ 個獨立 GR 規模樣本。抽樣次數越多，最大值通常越大；即使 $b$ 完全不變，最大餘震也會隨可用事件數改變。其分布為
#
# $$P(M_{\max}\le m\mid N)=\left[1-e^{-\beta(m-m_0)}\right]^N,\qquad m\ge m_0.$$
#
# 這個公式描述的是固定樣本數、未依規模篩選的獨立樣本。不能無條件搬到 ETAS 家族：在 ETAS 中，家族大小會受到成員規模影響。它的教學用途是展示極值的不穩定，而不是直接給出完整餘震預報。
#
# 另一個常見參考值是 $1/\beta$。它是理想獨立指數樣本中，最大與次大規模差的**期望值**，不是每個序列必須遵守的下界。若改用 GR 累積曲線下降到期望數 1 的位置 $m^*=a/b$，得到的也不是確定最大值；在 Poisson 計數假設下，仍有 $1-e^{-1}$ 的機率至少出現一個超過該門檻的事件。
#
# 下面左圖由目錄以指定視窗規則選出序列，右圖則是上述理想獨立抽樣的計算。兩者回答不同問題：左圖讓我們看見序列定義的影響，右圖隔離了抽樣極值本身的變化。
#
# %% tags=["remove-input"]
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
# 左圖和文獻中的平均規模差不一致時，不能立刻解讀成地區物理差異。應先核對主震選擇、空間窗、時間窗、完整度與序列是否重疊。本圖只對指定的較大候選主震使用視窗分類，是教學示範，並未重現 Chan 與 Wu（2013）的完整選樣程序。
#
# 右圖把主震規模固定後，增加餘震數，預期最大餘震會增大，規模差便下降。陰影反映隨機極值的變異，說明知道平均關係仍不足以準確預測個別序列。附錄 B 推導這個分布與其動差。
#
# ## 6.6 要不要先把餘震刪掉？
#
# 要計算前面的規模差，得先決定哪些事件屬於同一個序列。這正是除叢要處理的問題。除叢試圖把短期相互依賴的活動與背景活動分開，但背景如何定義，取決於研究目的和方法。
#
# 視窗法依候選主震的規模設定時間與距離範圍，將符合條件的事件標成同群。它容易重現，但固定視窗不會自動配合每次序列的形狀。連結法則依事件之間的關係擴充叢集，較有彈性，也可能因一段段串接而合併原本相距甚遠的活動。
#
# 最近鄰方法將時間、距離與親代規模合成距離分數，為每個事件尋找最接近的先前事件。分數若呈現可分辨的群體，可作為分類依據；雙峰不是所有目錄都保證出現，樹狀表示也不會自動消除長鏈。
#
# 隨機除叢換一種方式處理：先建立背景加觸發的模型，再為每個事件分配背景與候選親代的機率。它避免一定要選出唯一親代，卻把不確定性轉移到模型與參數上。Zhuang、Ogata 與 Vere-Jones（2002）把這個想法用在時空 ETAS；後面的估計章會從機率權重說明它如何運作。
#
# 選擇方法前先問分析要保留什麼。研究短期餘震預報，叢集正是重要訊息；研究背景活動率，可能需要分離短期活動；研究規模分布，則必須檢查選樣是否改變大、小事件的相對比例。除叢不是所有分析都需要的清理步驟。
#
# ## 6.7 分類會如何改變規模分布？
#
# 若每群只留下最大事件，小事件被排除的機率通常較高。即使原先所有事件遵循同一規模分布，留下來的子目錄也可能有較低的 $b$。Mizrahi et al.（2021）使用觀測與合成目錄比較不同方法，指出這種選擇效應必須與物理解釋分開。
#
# 一個簡單的檢查就能揭示外推問題。把全目錄與子目錄都擬合成 GR 直線，交點為
#
# $$m_x=\frac{a-a_{\rm main}}{b-b_{\rm main}},\qquad b\ne b_{\rm main}.$$ (eq:mx)
#
# 若子目錄斜率較小，往足夠大的規模外推，就可能預測「子目錄事件數比全目錄還多」。這當然不符合子集合關係，表示兩條直線不能同時無限外推。它不是某個物理臨界規模，也不能單靠交點判定完整危害度積分的偏差方向。
#
# 下圖以合成參數示範這個交叉。先追蹤兩條累積曲線，再讀它們的比值；重點是模型外推是否保有原本的集合關係。
#
# %% tags=["remove-input"]
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
# 在交點以前，較少事件和較平斜率可能部分抵銷；交點以後，比值超過 1 暴露了外推的不一致。實際應用還必須限制可信規模範圍，並檢查規模尾端、研究區域與其他模型假設。
#
# 長期危害度研究常以 Poisson 模型作為近似，因此會考慮除叢。但一條除叢後較平滑的累積曲線，不能證明增量獨立或等待時間服從指數分布。是否適用這個近似，仍需要檢驗；如果問題本來就關心短期群聚，直接使用叢集模型通常更符合目的。
#
# 下面回到原始目錄，比較這個教學視窗法保留前後的累積次數。它能呈現被移除的活動集中在哪裡，不能單獨證明留下來的事件是物理上彼此無關的「真正主震」。
#
# %% tags=["remove-input"]
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
# 曲線的階梯減少，表示短期集中活動被部分移除了。接下來應檢查規模分布和殘差是否也改變，以及這些改變是否符合原本研究目的。只評估曲線是否變平，容易忽略選樣造成的其他後果。
#
# 我們現在已知道餘震率、產能與空間分布可以怎麼描述，也看見硬性切分序列的侷限。把每個事件都視為可能的觸發源，是下一個模型的自然方向。不過，加入更多參數後，如何知道資料足以辨認它們？先到 {doc}`foundation_inference` 建立估計與不確定性的直覺，再在 {doc}`13_etas_structure` 把這些關係組合成 ETAS。
#
# 需要重算本章的正規化、極值分布或概似時，可查閱 {doc}`appendix_b_catalog`。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS；免費官方說明。先看「How the Aftershock Forecasts Work」，把餘震產能、發生率衰減與規模分布三件事分開理解，再回頭閱讀本章的經驗律。
# - [The Centenary of the Omori Formula for a Decay Law of Aftershock Activity](https://doi.org/10.4294/jpe1952.43.1) — Tokuji Utsu、Yosihiko Ogata、Ritsuko S. Matsu'ura，1995，*Journal of Physics of the Earth*；[免費全文](https://www.jstage.jst.go.jp/article/jpe1952/43/1/43_1_1/_article)。這篇回顧整理 Omori 公式、參數擬合與早期漏測問題，尤其適合追讀本章對 c 值與 p 值的解釋。
# - [The Effect of Declustering on the Size Distribution of Mainshocks](https://doi.org/10.1785/0220200231) — Leila Mizrahi、Shyam Nandan、Stefan Wiemer，2021，*Seismological Research Letters*；出版社全文可能需訂閱，[免費作者預印本](https://arxiv.org/abs/2012.09053)。以觀測與合成目錄研究除叢如何改變規模分布，對應本章「除叢也是一種資料選擇」的討論。
#
# - [New Empirical Tests of the Multifractal Omori Law for Taiwan](https://doi.org/10.1785/0120110237) — Ching-Yi Tsai、Guy Ouillon、Didier Sornette，2012，*Bulletin of the Seismological Society of America*；全文可能需訂閱。延伸閱讀臺灣餘震的長短期目錄完整度處理，以及如何用不同除叢方法比較結果；本章未重新估計該文的 p–規模迴歸。
#
# - [Importance of small earthquakes for stress transfers and earthquake triggering](https://doi.org/10.1029/2004JB003286) — Agnès Helmstetter、Yan Y. Kagan、David D. Jackson，2005。對照短期不完整與觸發產能；[免費機構版本](https://escholarship.org/uc/item/7q2277sq)。本章將其南加州完整度關係用作跨區教學參考，未當成臺灣校準結果。
# - [Maximum magnitudes in aftershock sequences in Taiwan](https://doi.org/10.1016/j.jseaes.2013.05.006) — Chung-Han Chan、Yih-Min Wu，2013，*Journal of Asian Earth Sciences*。對照 Båth、GR 與最大餘震的統計摘要，留意序列選取規則；出版社全文可能需訂閱。
# - [Space–time ETAS models and an improved extension](https://doi.org/10.1016/j.tecto.2005.10.016) — Yosihiko Ogata、Jiancang Zhuang，2006。分清產能和空間尺度；[免費機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。
# - [Stochastic Declustering of Space-Time Earthquake Occurrences](https://doi.org/10.1198/016214502760046925) — Jiancang Zhuang、Yosihiko Ogata、David Vere-Jones，2002。以機率分配背景與候選親代的原始論文；出版社全文可能需訂閱。
# - [GDMS 地球物理資料管理系統](https://gdms.cwa.gov.tw/) — 中央氣象署。圖中臺灣目錄的來源；本章視窗除叢是局部教學實作，未重現所有正式研究流程。
