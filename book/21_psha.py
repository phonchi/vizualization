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
# # 16. PSHA：從地震發生走到場址搖晃
#
# {doc}`複發模型 <20_recurrence_models>`回答了一段斷層未來發生事件的機率，但同一個地震對不同地方的影響不一樣。距離、破裂方向、地質與場址條件都會改變地動。現在把讀者的位置固定在一個場址，問：未來搖晃超過某個強度的機會有多大？
#
# 機率式地震危害分析（PSHA）把可能震源、規模、距離與地動變異合起來，形成一條危害曲線。我們先以時間不變的 Poisson 震源模型理解計算，再分清楚何時能延伸到時間相依的情況。
#
# 本章的地動例子都是教學用模型，不提供任何真實場址的設計值。重點是知道每一層輸入回答什麼問題，以及它如何影響最後的曲線。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 16.1 為什麼不能只選一場代表地震
#
# 一場近處中型地震可能帶來明顯搖晃，一場遠處大型地震也可能如此。只選「最大地震」會忽略距離，只選「最常見地震」又可能低估少見但強烈的尾部。
#
# 情境地震仍很有用：它讓我們檢查某次指定破裂的可能後果。PSHA 增加的是發生頻率的資訊，把許多情境按其發生率加權。兩種分析可以互補，沒有一種圖能單獨回答所有工程問題。
#
# %% [markdown]
# ## 16.2 把事件率與搖晃機率接起來
#
# Baker 的入門教材將工作分為五步：辨認震源、描述規模、描述震源到場址距離、預估地動，以及整合前述不確定性。先從一個震源看：每年發生 $\nu_s$ 次目標事件，每次只有一部分會讓場址的地動指標 $IM$ 超過門檻 $x$。
#
# 把所有震源相加，就得到年超越率：
#
# $$\lambda_{IM}(x)=\sum_s\nu_s\int\!\!\int
# P(IM>x\mid m,r,s)\,f_{M,R\mid s}(m,r)\,dm\,dr.$$ (eq:psha)
#
# $\nu_s$ 的單位是每年，內部積分是一次事件超越門檻的機率，所以結果也是每年。若規模與距離在給定震源後獨立，聯合密度才可寫成 $f_Mf_R$；有限破裂尺度等因素可能使它們相依。
#
# 這條式子運用全期望，將各種情境加權平均。單純計算長期平均超越次數，不要求所有事件都是 Poisson；但稍後把年率轉成「至少一次」機率，便需要額外的發生過程假設。
#
# %% [markdown]
# ## 16.3 給大規模端一個明確邊界
#
# 前面 GR 律描述規模越大越少。工程危害分析還要指定使用範圍 $m_{\min}\le m\le m_{\max}$，並重新正規化：
#
# $$f_M(m)=\frac{\beta e^{-\beta(m-m_{\min})}}
# {1-e^{-\beta(m_{\max}-m_{\min})}},\qquad \beta=b\ln10.$$ (eq:trunc-gr)
#
# 分母補回截去上尾後失去的質量。$\nu_s$ 也必須對應同一規模範圍；不能一邊改門檻，一邊保留原來的率而不檢查。
#
# 最大規模通常牽涉構造與資料限制，並不是從少量觀測中精確知道的硬邊界。不同合理 $m_{\max}$ 可以形成危害模型的分支，尤其在較低超越率的尾部比較結果。推導、抽樣與平均規模見 {doc}`附錄 F <appendix_f_hazard>`。
#
# %% [markdown]
# ## 16.4 同樣規模與距離，搖晃仍不同
#
# 地動預估式描述給定規模、距離等條件後的地動分布。常見近似是假設 $\ln IM$ 為常態，中心為 $\mu_{\ln IM}(m,r)$，標準差為 $\sigma_{\ln IM}$，所以
#
# $$P(IM>x\mid m,r)=1-\Phi\!\left(\frac{\ln x-\mu_{\ln IM}(m,r)}{\sigma_{\ln IM}}\right).$$ (eq:gmpe-exceed)
#
# $\Phi$ 是標準常態累積分布。門檻若在中位數以上，超越機率仍不為零；同樣地，不應把中位數低於門檻的所有事件直接刪掉。許多可能事件的尾部累積後，也可能貢獻顯著危害。
#
# 常見的 $\varepsilon$ 是標準化殘差，表示地動偏離對數中心多少個標準差。它不是一個新震源，也不是可預知的「加強係數」。
#
# 下面依序畫出規模密度、指定門檻的條件超越機率，以及兩者乘上事件率後的貢獻。最後一格曲線下的面積，正是各震源的年超越率；兩個震源各固定一個示意距離。
#
# %% tags=["remove-input"]
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


def source_nu(nu_ref, m_max, m_min=M_MIN):
    """在相同有限上界內換規模門檻，保持未正規化規模率不變。"""
    return nu_ref * (np.exp(-BETA * m_min) - np.exp(-BETA * m_max)) / (
        np.exp(-BETA * M_REF) - np.exp(-BETA * m_max))


def source_rate(x, r, nu_ref, m_max, m_min=M_MIN, n=600):
    """單一震源對 λ_IM(x) 的貢獻（對規模數值積分）。"""
    m = np.linspace(m_min, m_max, n)
    integ = p_exceed(x, m, r) * mag_pdf(m, m_max, m_min)
    return source_nu(nu_ref, m_max, m_min) * np.trapezoid(integ, m)


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
    integ = source_nu(nu, mx) * f_m * p_x
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
# 讀這張圖時，先固定一個規模和距離，再想像重複很多次事件所得的搖晃分布。PSHA 積分的就是每種情境跨越門檻的比例，按其發生率與情境機率加權。
#
# %% [markdown]
# ## 16.5 危害曲線是一整組門檻的回答
#
# 逐漸提高地動門檻 $x$，能超過它的事件只會變少，年超越率應隨之不增加。下面沿橫軸選一個搖晃門檻，再讀縱軸，便得到這個教學場址的年超越率。
#
# %% tags=["remove-input"]
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
# 同一條曲線的不同位置，可能由不同地震群貢獻。較低門檻可能收到許多中小事件，較高門檻則更依賴近處、較大規模或地動偏高的事件。
#
# 下面改變最低規模門檻，看看哪些部分受到影響。這是在檢查截斷是否漏掉重要貢獻，而不是預先斷言小地震對工程一定不重要。
#
# %% tags=["remove-input"]
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
# 門檻敏感度應與震源率的定義一起檢查。若刪去小地震後仍用全部事件率重新正規化，便在無意間把它們的率搬到較大規模，得到的差異不再只是「忽略小地震」。
#
# %% [markdown]
# ## 16.6 回歸期不表示每隔固定時間發生
#
# 在平穩模型下，年超越率的倒數 $T_R=1/\lambda_{IM}(x)$ 常稱為回歸期。它是長期頻率的表達，不是週期，也不是距離下一次事件的倒數計時。
#
# 若超越事件服從固定率 Poisson 過程，$T$ 年內至少一次超越的機率才是
#
# $$P(N_T\ge1)=1-e^{-\lambda_{IM}(x)T}.$$
#
# 這個換算來自零事件機率。它不是把「每年機率」直接乘年數；只有乘積很小時線性近似才合理。更不能把任何自激發模型的平均數都代進去，因為相同平均數可以伴隨不同的零事件機率。
#
# ```{admonition} 回歸期、等待時間與機率
# :class: dropdown
#
# Poisson 假設下，上式可以反解年率 $\lambda=-\ln(1-P)/T$。更新過程則應使用上一章的條件機率；一般叢集過程需由完整計數模型或模擬計算。{doc}`附錄 F <appendix_f_hazard>`說明各自的適用條件。
# ```
#
# %% [markdown]
# ## 16.7 哪些地震在貢獻這段危害
#
# 危害曲線加總了很多可能事件。有時我們想倒回來問：在超過這個門檻的事件中，哪些規模與距離佔比較大？這就是反聚合。
#
# 每個規模距離格的「事件率 × 超越機率」就是它的貢獻；再除以總超越率，可畫成貢獻比例。這類似條件分布，但條件是指定門檻的超越事件，而非下一次地震一定會有某個規模。
#
# 下圖把總率拆回規模及地動殘差 $\varepsilon$，並在圖題列出相應的平均距離。提高門檻後，主要貢獻區可能移動。反聚合因此能協助選擇代表情境，卻不能把最亮的一格解讀為唯一會造成搖晃的地震。
#
# %% tags=["remove-input"]
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
        tot += source_nu(nu, mx) * np.trapezoid(
            p_exceed(x, m, r) * mag_pdf(m, mx), m)
    return tot


def deagg_eps(x, lo, hi, n=1500):
    """ε 落在 [lo, hi) 的超越率貢獻（首尾開放，各箱總和等於 λ）。"""
    tot = 0.0
    for _, r, nu, mx in FAULTS:
        m = np.linspace(M_MIN, mx, n)
        w = np.clip(norm.cdf(hi) - norm.cdf(np.maximum(eps0(x, m, r), lo)),
                    0.0, None)
        tot += source_nu(nu, mx) * np.trapezoid(w * mag_pdf(m, mx), m)
    return tot


def deagg_means(x, n=2000):
    """反聚合的 (平均規模, 平均距離, 平均 ε)。"""
    num_m = num_r = num_e = den = 0.0
    for _, r, nu, mx in FAULTS:
        m = np.linspace(M_MIN, mx, n)
        p = p_exceed(x, m, r)
        w = source_nu(nu, mx) * mag_pdf(m, mx) * p
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
# 圖上各組百分比加總為一百，表示相對貢獻，不是各格的絕對發生機率。左側比較規模，右側比較地動殘差；要知道這種搖晃多常出現，仍須回到危害曲線。
#
# %% [markdown]
# ## 16.8 當發生率隨時間改變
#
# 大地震之後，未來一週的危害可能與長期平均不同。先採一個最簡單的延伸：假設未來超越率 $\lambda_{IM}(x;t)$ 是已知、確定的時間函數，且事件為非齊次 Poisson，則
#
# $$P(\text{at least one exceedance in }[0,T])
# =1-\exp\!\left[-\int_0^T\lambda_{IM}(x;t)\,dt\right].$$ (eq:nhpp)
#
# 與固定率相比，只是把率乘時間改成時間積分。下面用率放大的 Poisson 示意，呈現短時間窗曲線如何上升；它沒有模擬 ETAS 的後續世代，因此不是完整 ETAS 危害計算。
#
# 對 ETAS，未來地震會再觸發事件，必須傳遞整份未來歷史。較直接的途徑是模擬目錄及各事件地動，再統計時間窗內至少一次超越的比例。一般不能把平均 ETAS 率代入上式，也不能把已發生目錄上的強度路徑當成預先已知的非齊次 Poisson 率。
#
# %% tags=["remove-input"]
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
# 率提高時，固定短時間窗的超越機率也會增加，但這張示意圖的放大倍數不是任何真實序列的估計。要轉成場址可用的產品，需要經過完整模型與地動驗證。
#
# 長期端也有時間資訊。下面回到上一章的 BPT，畫出距上次事件不同年齡下，未來固定時間窗的下一次破裂機率。這與將整個區域率乘上一個短期放大係數，是不同的建模操作。
#
# %% tags=["remove-input"]
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
# BPT 圖回答的是下一次指定破裂，不是直接給出至少一次地動超越。若時間窗可能包含多次破裂，或每次破裂未必超過指定地動門檻，還要整合相應的更新歷史與地動分布。
#
# %% [markdown]
# ## 16.9 把方法帶回臺灣
#
# 臺灣的危害模型整合地震目錄、孕震構造、地質與測地資訊，以及合適的地動預估式。TEM PSHA2020 提供已發表的具體例子；本章末列出正式論文和計畫說明，供讀者對照每一類輸入位於哪一步。
#
# 這裡先保留方法上的連結。到了第二部，我們會逐一認識臺灣的地動、GNSS 等觀測，再問這些資料能限制模型哪一部分。這樣便不必在還沒理解資料來源前，先記住一套地方參數或結果。
#
# %% [markdown]
# ## 16.10 危害還不是損失
#
# 相同搖晃對不同建築、裝置與使用者會有不同後果。PSHA 描述地動，後面還需結構反應、損傷與損失模型，才走到工程風險。不能把「某地超越率較高」直接當成「每個人的損失都較大」。
#
# 短期預報與危害模型的連結，則使使用者能夠根據當前事件歷史重新理解可能搖晃。這也帶來更新、驗證與溝通的需求：每次資料修訂可能改變事件率，地動模型和場址條件又帶來另一層不確定性。
#
# %% [markdown]
# 模型通常區分自然變異與知識不足。給定模型下每場事件地動不同，是一層變異；不同地動式、最大規模或斷層幾何，是另一層模型差異。兩層應分別傳遞，不能把所有差異都塞進同一個標準差。
#
# 邏輯樹可以保留多個合理分支，但其曲線範圍取決於候選假設及權重，也不是全部未知的保證範圍。讀危害圖時，應同時知道它對哪些條件敏感。
#
# %% [markdown]
# 從目錄到模型、從模型到評分，再從事件走到場址地動，統計預報的主要環節已經連起來。真正持續發布時，每一環還要在資料未齊、會修訂且使用者需要及時資訊的情況下運作。
#
# {doc}`作業化預報 <22_operational_systems>`將把這些環節接成一個持續更新的流程，也讓我們準備好在第二部閱讀臺灣的實際觀測。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Introduction to Probabilistic Seismic Hazard Analysis](https://scits.stanford.edu/sites/g/files/sbiybj22081/files/media/file/baker_2013_intro_psha_v2_0.pdf) — Jack W. Baker（2013），White Paper Version 2.0（免費教材）。本章 PSHA 推導的入門來源，依五個步驟走過震源、規模、距離、地動與積分；接著讀反聚合及回歸期，理解危害曲線背後有哪些事件在貢獻。
# - [TEM 計畫概要與研究成果](https://tem.tw/TEM2020/portfolio-overview.html) — Taiwan Earthquake Model 計畫團隊，TEM2020 專頁（免費）。先看孕震構造、測地資料與危害度評估如何分工，再循頁面列出的論文了解臺灣模型的輸入來源。
# - [Probabilistic seismic hazard assessment for Taiwan: TEM PSHA2020](https://doi.org/10.1177/8755293020951587) — Chung-Han Chan 等（2020），*Earthquake Spectra*（出版社全文可能需訂閱；[中央大學免費全文](https://www.gep.ncu.edu.tw/storage/thesis/2020/2020%20Chung-Han%20Chan_ES2.pdf)）。本章臺灣案例的正式文獻，說明孕震構造、地震目錄、地動預估式與場址效應如何納入 PSHA；可對照教材檢查每一類不確定性出現在哪一步。
#
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。
#
