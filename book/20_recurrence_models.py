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
# # 15. 複發模型：已經等了多久，也是資訊
#
# {doc}`模型組合 <19_ensembles>`連結了短期觸發與較長期的活動分布。現在把空間範圍縮到一段斷層，把時間拉長到數十年至數百年：若知道它上次大破裂的時間，這項資訊應如何影響未來機率？
#
# 這個問題不能只用區域目錄的平均率回答。它需要定義哪些破裂算同一類事件、如何估計古地震年代，以及兩次破裂間隔有多不規則。資料通常比儀器目錄少得多，所以直覺與不確定性都很重要。
#
# 我們從等待時間開始，先理解「已經等到現在」這個條件，再比較幾種複發分布。最後把它們接回點過程，讓下一章能將事件發生模型轉成場址的地動危害。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 15.1 等待會不會改變下一步
#
# 想像兩種簡化機制。第一種每年都有相同的小機率發生事件，不管之前等多久；第二種事件之後需要重新累積，短時間內重複較不容易。前者適合用指數間隔描述，後者則需要能表達間隔形狀的分布。
#
# 指數分布有無記憶性：已經等了 $s$ 年，再等 $t$ 年仍無事件的條件機率，與從現在重新開始等 $t$ 年相同。這是一個明確可檢查的假設，不是對所有斷層物理的結論。
#
# 彈性應變累積與釋放提供另一種動機：同一破裂面的再次大破裂可能需要較長恢復時間。但不能由此說剛破裂後整個區域安全，因為周圍仍可能出現大量餘震，甚至其他段落的大事件。區域觸發與同一斷層段複發，分析物件不同。
#
# 在點過程語言中，兩者都可使用條件強度；差別在於歷史以何種方式進入模型。更新過程只記距上次指定事件的時間，而 ETAS 累積許多過去事件的觸發作用。
#
# %% [markdown]
# ## 15.2 從間隔分布到條件機率
#
# 令 $X$ 是一次完整複發間隔，其累積分布為 $F$、密度為 $f$。我們已經知道這次間隔超過 $T$ 年，想問接下來 $\Delta T$ 年內發生的機率：
#
# $$P(T<X\le T+\Delta T\mid X>T)
# =\frac{F(T+\Delta T)-F(T)}{1-F(T)}.$$ (eq:renewal-prob)
#
# 分子是原分布在未來時間窗中的面積；分母去掉那些本來會更早發生、但已被目前觀測排除的情況。因此不能只看密度曲線未來那一小段，還要條件在「已等待到今天」。
#
# 把時間窗縮短，就得到危害函數，也就是瞬時發生率：
#
# $$h(t)=\frac{f(t)}{S(t)},\qquad S(t)=1-F(t).$$ (eq:hazard-def)
#
# $S(t)$ 是存活函數，此處指間隔仍未結束的機率。$h(t)$ 的單位為每年，並不是一個介於零與一之間的機率。在更新過程中，絕對時間的條件強度是 $\lambda^*(t)=h(t-t_{\rm last})$；事件發生後，年齡重新從零計算。
#
# 因此可以沿用前面的點過程概似與殘差方法。不過，斷層資料很少時，能否辨認 $h$ 的形狀，往往比會不會代公式更困難。
#
# %% [markdown]
# ## 15.3 三種間隔形狀，三種建模想法
#
# Weibull 用形狀參數 $k$ 控制等待效應。$k=1$ 回到無記憶的指數分布；$k>1$ 時危害函數隨等待增長。它的形式簡單，但持續增長到很長時間是否合理，需要配合使用範圍判斷。
#
# 對數常態適合描述倍數式間隔變異。它的密度有一個峰，右側保留長尾；等得非常久仍沒有事件時，條件分布會逐漸偏向長間隔那一群，所以危害函數不必一直上升。
#
# BPT（Brownian passage time）則從帶隨機擾動的載入量出發。平均向上累積，到門檻就記一次破裂；不同路徑因擾動而在不同時間首次碰到門檻。用平均間隔 $T_r$ 與間隔變異係數 $c_v$ 表示，其密度為
#
# $$f(t)=\sqrt{\frac{T_r}{2\pi c_v^2t^3}}
# \exp\!\left[-\frac{(t-T_r)^2}{2c_v^2T_rt}\right],\qquad t>0.$$ (eq:bpt-pdf)
#
# $c_v$ 小表示間隔較集中，較大則表示較不規則。這是一個可產生複發分布的理想化機制，不是直接量到斷層全部應力。下面將三個分布調到相同平均與變異係數，先比較密度：光看密度，可能不容易察覺它們對很長等待時間的不同看法。
#
# %% tags=["remove-input"]
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
# 三條曲線的主要機率質量可能相近；當我們已等過這個主要區段，真正影響未來條件機率的卻是尾部。接著用危害函數看同一組模型。
#
# %% [markdown]
# ## 15.4 尾部如何改變長等待的解讀
#
# 下面比較 $h(t)$。固定其他條件後，哪條曲線在目前等待時間附近較高，哪個模型就給更高的短時間窗條件機率。不要把整張密度圖的峰值直接當成「現在最危險」的判準。
#
# %% tags=["remove-input"]
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
# Weibull 在 $k>1$ 的設定下持續上升；對數常態最終下降；BPT 的長時間危害率趨向有限正值 $1/(2T_rc_v^2)$。這些是各自尾部假設的數學後果，不代表資料已經看過足夠長時間來驗證。
#
# 因此，「已超過平均間隔」並不足以決定風險升降。平均只是一個摘要，需要配合分布形狀、上次事件時間及參數不確定性。若目錄只記錄三四次破裂，尾部差異通常難由該段資料單獨決定。
#
# %% [markdown]
# ### 間隔越不規則，單一平均值越不夠
#
# 下面固定 BPT 平均間隔，只改變 $c_v$。這樣就能將「平均多久一次」與「每次相差多大」分開看；兩者都會影響當前等待時間下的條件機率。
#
# %% tags=["remove-input"]
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
# 較小的 $c_v$ 把間隔集中在平均附近，較大的 $c_v$ 留下更寬的可能範圍。相同平均複發時間可以給出不同的未來機率，因此引用參數時不能只報 $T_r$。
#
# %% [markdown]
# ### 看見分布背後的路徑
#
# 接下來以帶漂移的隨機載入路徑示意 BPT 的來源。每條線在碰到門檻之前都有不同起伏；首次碰到門檻的時間，才是間隔樣本。
#
# %% tags=["remove-input"]
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
# 路徑說明的是為什麼同樣平均載入也會有不同破裂時刻。有限步長的模擬可能漏掉兩取樣時刻之間的碰撞，因此只是近似驗證；解析首達時間與參數對照放在 {doc}`附錄 F <appendix_f_hazard>`。
#
# 這也說明更新假設的強度：每次破裂後都重新抽一個同分布間隔，便忽略了更早的載入歷史與其他斷層互動作用。模型是否需要放寬，取決於資料與問題尺度。
#
# %% [markdown]
# ```{admonition} 閱讀 BPT 的參數
# :class: dropdown
#
# $T_r$ 與時間軸使用相同單位；$c_v=\mathrm{sd}(X)/\mathbb E[X]$ 無量綱。不同軟體的 inverse Gaussian 參數不一定直接使用這兩個量。{doc}`附錄 F <appendix_f_hazard>`提供對照，避免相同數字代表不同分布。
# ```
#
# %% [markdown]
# ## 15.5 規模統計能否補足間隔資料
#
# 前面 GR 律描述事件大小，這裡的 $c_v$ 描述事件間隔；它們不是同一個量。但若加上應力累積、規模與釋放量的關係等物理假設，就可能推導兩者之間的連結。
#
# 這類模型的教學價值，是讓我們看見額外資訊如何透過假設進入推論。它不能支持「知道區域 $b$ 值就知道每條斷層 $c_v$」這種直接換算。區域內事件的來源很多，與同一斷層段的重複破裂並非同一抽樣單位。
#
# 實務上應把此類關係當成可比較的建模假設，與直接間隔資料、地質限制及其他分布一起檢查，而不是用一個看似精確的公式消除資料不足。
#
# %% [markdown]
# ## 15.6 讓過去釋放量留在歷史中
#
# 更新過程在事件後重設年齡。另一種模型則保留累積與釋放的狀態，例如用
#
# $$\lambda^*(t)=\exp\!\left[a_s+b_st-c_s\sum_{t_i<t}R_i\right],\qquad b_s,c_s>0,$$
#
# 表示自我修正強度。$R_i$ 是第 $i$ 次事件的釋放量代理，$b_s$ 描述隨時間增加的趨勢，$c_s$ 控制一次釋放如何降低後續率。指數確保強度非負。
#
# 相較 ETAS 的事件後上升，這裡是事件後下降；相較更新過程，它記住的不只是最後一次時間，還包括累積釋放。下面的示意曲線讓這兩種歷史反應可以直接比較。它們是不同尺度的理想化模型，不能只憑正負號決定哪個較符合所有地震。
#
# %% tags=["remove-input"]
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
# 圖上的下跳由模型指定的釋放量產生，上升段則來自持續載入。若要套到真實資料，還需說明如何由規模換算 $R_i$、研究區是否漏掉重要事件，以及是否需要鄰近區域的互動作用。
#
# %% [markdown]
# ## 15.7 少量古地震資料，如何保留不確定性
#
# 儀器目錄通常有明確時間戳記；古地震年代常是一個區間，且可能漏掉破裂。平均間隔和 $c_v$ 的不確定性因此不只是樣本小，還包含年代、對應事件與斷層分段的選擇。
#
# 若觀測只到今天，而下一次還沒有發生，最後這段等待是右設限資料。它仍提供「間隔至少這麼長」的資訊，不能因為不是完整間隔就丟掉。若觀測開始時距前次事件的年齡未知，起點也需要另行處理；Daley 與 Vere-Jones 的更新過程章特別區分這些抽樣條件。
#
# 邏輯樹可將不同合理假設分支保留，再計算其結果差異。但分支權重不是從少量資料中自動得到的真理，候選樹也可能漏掉重要情境。相比只提供一個最可能值，呈現機率對這些假設的敏感度更能反映知識範圍。
#
# %% [markdown]
# 若不同斷層之間有相互作用，或間隔分布隨時間改變，單一更新模型可能不足。可以擴充套件成耦合、隱藏狀態或時變模型，但每次擴充套件都增加資料需求。先問現有觀測能否區分這些假設，再決定模型要複雜到哪裡。
#
# %% [markdown]
# 從區域觸發到斷層複發，我們已經能用不同歷史資訊描述事件何時、何地及多大。防災與工程接著關心的是：「這些可能事件，會讓某個場址搖到什麼程度？」
#
# {doc}`下一章的 PSHA <21_psha>`將把事件模型與地動分布接起來。這個轉換需要整合許多可能地震，而不是只挑一個最像平均情況的事件。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Related Distributions：Hazard Function 與 Survival Function](https://itl.nist.gov/div898/handbook/eda/section3/eda362.htm) — NIST／SEMATECH，*e-Handbook of Statistical Methods*（免費）。先看機率密度、存活函數與危害函數的關係，釐清本章「已經等到現在，再發生的瞬時率」與一般發生機率的差別。
# - [A Brownian model for recurrent earthquakes](https://doi.org/10.1785/0120010267) — Mark V. Matthews、William L. Ellsworth、Paul A. Reasenberg（2002），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[USGS 免費摘要](https://pubs.usgs.gov/publication/70024443)）。BPT 複發模型的核心論文，從帶有布朗擾動的載入過程推到首達時間，銜接本章的物理直覺與長時間危害率。
# - [Modeling the earthquake occurrence with time-dependent processes: a brief review](https://doi.org/10.1007/s11600-019-00284-4) — Ourania Mangira、Christos Kourouklas、Dimitris Chorozoglou、Aggelos Iliopoulos、Eleftheria Papadimitriou（2019），*Acta Geophysica*（全文可能需訂閱）。將本章的更新過程、複發分布與應力釋放放回時間相依模型的全貌，閱讀時比較各模型適用的時間尺度、資料與假設。
#
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。
#
