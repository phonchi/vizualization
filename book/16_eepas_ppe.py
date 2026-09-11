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
# # 11. EEPAS 與 PPE：讓尺度關係成為機率預報
#
# {doc}`上一章 <15_psi_phenomenon>`留下了一個實際限制：Ψ 的回溯辨識需要知道主震，而真正發預報時不能使用這個資訊。EEPAS 的建模選擇，是讓每一個已發生的地震都對未來提供一小部分率，而不先判定它是否屬於某個前兆群。
#
# 想像在每顆已知地震周圍放上一片會隨時間改變的機率密度；規模不同，這片分布的時間、空間及目標規模尺度也不同。把所有貢獻加起來，再搭配長期背景，便得到可以在任何發布時刻計算的預報。
#
# 這個模型假設不表示每顆地震真的會被更大地震跟隨。它的價值要由整份預報在未來資料上的表現判斷。我們先讀懂單一事件留下的分布，再理解如何相加、如何處理看不到的資料，最後回到比較基準。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 11.1 從一群事件到一個事件的貢獻
#
# Ψ 的 $M_P$ 是前兆群中最大三個規模的平均。EEPAS 用個別事件規模 $m_i$ 作為輸入，將回歸的中心與散佈轉成未來目標事件的分布。這是新的模型假設，不能把群體統計量與個別觀測視為完全相同。
#
# 較小的輸入事件很多，每顆分到的貢獻有限；較大的事件較少，卻指向不同的未來尺度。模型不必對某顆地震宣佈「這就是前震」，仍可讓一片區域的預報率隨新資料逐漸改變。
#
# 替換輸入定義後，回歸係數也須重新估計。將人工辨識得到的參數直接套入 EEPAS，等於忽略兩種樣本的差異。前一章的經驗關係提供建模方向，參數與預報能力則留給資料檢查。
#
# %% [markdown]
# ## 11.2 規模、時間和位置的三片分布
#
# 先回想基礎迴歸：$Y=a+bX+\epsilon$ 若假設 $\epsilon$ 為常態，給定 $X$ 後，$Y$ 就有一個以回歸線為中心的常態分布。EEPAS 將這個想法用在目標規模、等待時間的對數與空間尺度。
#
# 規模核 $g(m\mid m_i)$ 描述目標規模的分布；時間核 $f(\tau\mid m_i)$ 描述正等待時間 $\tau=t-t_i$；空間核 $h(x,y\mid x_i,y_i,m_i)$ 則把貢獻攤在事件周圍。原始模型使用常態規模核、對數常態時間核及圓對稱二維常態空間核。
#
# 先看時間：如果 $\log_{10}\tau$ 呈常態，表示不確定性是倍數式的。等候一年到兩年，和十年到二十年，在對數軸上相隔相同距離；直接對天數假設常態，則是另一種誤差結構。
#
# $$f(\tau\mid m_i)=\frac{1}{\tau\sigma_T\ln 10\sqrt{2\pi}}
# \exp\!\left[-\frac{(\log_{10}\tau-a_T-b_Tm_i)^2}{2\sigma_T^2}\right],\qquad \tau>0.$$
#
# $\tau$ 若用天，密度的單位是每天。分母的 $\tau\ln10$ 來自座標轉換，讓曲線對天數積分後仍等於一。中心等待時間的中位數是 $10^{a_T+b_Tm_i}$ 天，不等於平均數或眾數。
#
# 下面三張核圖應一起讀：沿不同顏色追蹤輸入規模變大後，預報往哪個目標規模、哪個時間與多大的空間範圍移動。這些形狀來自模型設定；相乘還假設給定輸入事件後三者可分離，並不是資料已證明它們獨立。
#
# %% tags=["remove-input"]
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
# 三個核分別控制「哪種規模、多久以後、在哪裡」，核外的權重才控制「總共多少」。為了看清時間尺度的差異，下面只比較 EEPAS 的時間核與 ETAS 的 Omori 核，先不比較它們預報的總事件數。
#
# %% tags=["remove-input"]
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
# 在 $c>0$ 的 Omori 核中，剛發生後的密度有限且最高，之後逐漸衰減；EEPAS 的對數常態時間核則從零附近上升，到正等待時間才達到峰值。兩者把資訊放在不同時段，這是後面模型組合的動機之一。
#
# 交叉點只屬於這組示意參數，不能當成區分觸發與其他機制的物理界線。Omori 核的長尾也可能延伸很久，對數常態密度在對數時間軸上並非左右對稱；真正對稱的是 $\log_{10}\tau$ 的常態密度。
#
# %% [markdown]
# ## 11.3 把貢獻相加之前，先校準總量
#
# 有了單一事件的三個核，可以寫出 EEPAS 的基本率密度：
#
# $$\lambda(t,m,x,y)=\mu_E\lambda_{\rm PPE}(t,m,x,y)
# +\sum_{t_0\le t_i<t,\ m_i\ge m_0}\eta(m_i)w_i
# f(t-t_i\mid m_i)g(m\mid m_i)h(x,y\mid x_i,y_i,m_i).$$ (eq:eepas-rate)
#
# $\mu_E$ 控制背景項，$w_i$ 控制個別事件的權重，$\eta(m_i)$ 調整不同輸入規模的總貢獻。時間、面積與規模核分別正規化，不代表整體長期率已經正確；把更多輸入事件加進來，總量仍會變動。
#
# 我們希望模型累積很多年之後，能維持與目錄相容的規模頻率分布。假如每顆小地震與大地震都分配相同總量，小地震龐大的數量便會主導結果。$\eta$ 的作用是配合 GR 輸入分布與規模核，平衡這個效應。
#
# 下面比較不做這個調整與做了調整後的長期規模分布。先看曲線斜率，再看接近輸入門檻的彎曲：斜率問題與漏掉小地震貢獻的問題，需要分開處理。
#
# %% tags=["remove-input"]
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
# 正規化使輸出的長期尺度與指定 GR 律相容；有限輸入門檻仍會留下低目標規模端的缺口。這不是靠把整張預報乘上一個常數就能完全修正，因為不同目標規模缺少的比例不同。
#
# %% [markdown]
# ## 11.4 目錄看不到的事件去了哪裡
#
# 前面完整度章談的是觀測：小事件可能未被目錄記錄。在 EEPAS 裡，這也變成預報問題，因為未記錄事件原本可能對未來有所貢獻。目標越接近輸入門檻，缺少這部分來源的影響通常越明顯。
#
# 模型以 $\Delta(m)$ 描述：在指定的 GR 與常態核假設下，目標規模 $m$ 的長期貢獻有多少比例由 $m_i\ge m_0$ 的事件提供。$\Delta$ 接近一時缺口較小，接近零時則代表補償高度依賴模型外推。
#
# 補償不會重新找回漏測的個別事件，也不能還原它們實際的位置與時間。它只修正假設模型下的平均比例。因此輸入門檻須配合目錄品質，不能把很大的補償係數當成「低品質資料也能照常使用」。完整推導與適用條件放在 {doc}`附錄 D <appendix_d_eepas>`。
#
# %% [markdown]
# ## 11.5 避免一串餘震被重複解讀
#
# 若每個事件都給完整權重，一場主震後的大量餘震可能被解讀成很多份新的中期資訊。這時可以沿用 {doc}`ETAS 估計 <14_etas_estimation>`的隨機除叢想法，讓較可能屬於背景的事件權重較高。
#
# 在採用 ETAS 背景機率的版本中，可取 $w_i=\phi_i$；另一些版本使用不同除叢規則或等權重。權重是模型推論，並不是觀測已揭露某事件的真實身分。若 ETAS 參數、完整度或背景率改變，權重也可能改變。
#
# 這樣做連起短期與中期模型：先用觸發模型辨認已能解釋的叢集，再問剩餘資訊是否改善較長期的預報。權重是否有幫助，仍要透過固定其他條件的比較，而非因為用了較複雜方法就預設較好。
#
# %% [markdown]
# ## 11.6 先和「附近以前發生過地震」比較
#
# 預報若只和空間均勻模型比較，可能只是重新發現地震集中在活動構造附近。更有意義的基準，是保留過去地震的空間資訊，卻不使用 Ψ 的時間變化。
#
# PPE（Proximity to Past Earthquakes）依過去事件的位置與尺度建構平滑地震度，再搭配規模分布及長期率。它代表一個具體問題：「在已知哪裡常發生地震之後，新的時間資訊還增加了多少？」SUP 則是空間均勻的簡單參考，用來顯示最基本的地理集中性。
#
# 下面用臺灣公開目錄畫出 PPE 型空間項。這是一個保留過去震央資訊的教學實作，尚未擬合成正式臺灣預報；它讓讀者先看見平滑背景所使用的資訊。平滑距離和尾部會改變未知地區分到多少率，稍後會直接影響對數分數。
#
# %% tags=["remove-input"]
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
# 核的峰值高，不代表整份預報較好。它可能增加近處事件的分數，也可能讓稍遠的事件落在極低率區。評估時必須連同沒有事件的位置以及整體總量一起計算，不能只看圖上幾個命中的高峰。
#
# %% [markdown]
# ## 11.7 當目標也包含餘震
#
# 基本模型把注意力放在中長期較大事件。如果預報目標包含它們的餘震，就要把「未來可能的主震」與「主震之後的後代」接起來。這和 ETAS 的世代累積有相通之處，但不是把兩個模型的率任意相加。
#
# 長時間窗下，若主震與餘震的時間差相對於預報尺度很短，可以考慮近似合併時間分布。空間上則需對未知的主震位置積分：目標餘震的位置等於主震位置加上相對位移，對應兩個空間分布的卷積。
#
# 在兩者都是獨立二維常態的簡化情況下，卷積仍是常態，變異數相加；推導見附錄 D。為了先建立疊加的直覺，下圖回到基本模型，畫出三個不同規模事件在未來時間與位置的貢獻。它只畫時間空間剖面，並未實作主震與餘震的卷積。
#
# %% tags=["remove-input"]
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
# 每個標記後方逐漸出現一片分布，較大輸入事件的時間與空間尺度也不同。亮區是多個貢獻疊加後的相對率，不是已知的未來震央；若要加入目標餘震，才進一步使用上面說明的卷積。
#
# %% [markdown]
# ## 11.8 新開始的目錄，如何做長期預報
#
# 即使規模門檻足夠低，目錄開始之前的事件仍然看不到。對數常態時間核可能讓很久以前的事件持續貢獻；因此短目錄會漏掉一部分長時間來源。這是時間完整度，和規模完整度是兩個不同方向的缺口。
#
# 想像今天要發布未來一年的預報：往前有三年資料或三十年資料，在 EEPAS 裡不一定等價。大目標常對應更長的時間尺度，因此缺口也可能隨目標規模改變。
#
# Rhoades 等人（2020）用固定目錄前置時間與補償模型研究這件事。下面固定可用目錄長度 $L$，改變目錄末端到預報時刻的落後時間 $T$。圖顯示資料越陳舊時，可見貢獻比例如何依目標規模改變。重新估參數可能讓模型部分適應短目錄，但也可能透過空間時間尺度補償，掩蓋原來缺少的資訊。
#
# %% tags=["remove-input"]
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
apply_layout(fig, title=f"前兆可見比例取決於時間窗與目標規模"
                        f"（目錄長度 L = {L_CAT:.0f} 年）",
             hovermode="x", height=430)
fig

# %% [markdown]
# 閱讀這些曲線時，請同時看歷史長度和目標規模，不要只找一個全區通用的「至少幾年」。固定長度的可用時間窗沿時間核移動時，可見比例可能先增加再下降；它不保證隨落後時間單調減少。當資料不足時，較合理的預報可能需要提高背景的作用，並如實表達額外不確定性。前置時間與補償規則都應在測試前固定。
#
# %% [markdown]
# ## 11.9 成績必須連同題目一起讀
#
# 下面保留 Rhoades（2011）日本本土研究的回溯比較，三根柱代表不同目標規模的 EEPAS 相對 PPE 增益，虛線則為跨規模級距移植參數的結果。這些數字是為 CSEP 前瞻測試準備模型時的回溯分析，不是後來前瞻實驗已取得的成績。
#
# %% tags=["remove-input"]
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
# 一些研究中，EEPAS 在較長時間尺度提供了超出平滑背景的資訊；另一些版本或設定的改善較有限。這和前面的核形狀一致，但核的形狀本身不能保證預報技巧。
#
# 尤其要區分擬合資料、保留測試資料和真正前瞻資料。先看到整份目錄再調時間尺度，與在每次發布當下只能使用過去資料，是不同的實驗。較好的回溯分數可以支持進一步測試，不能直接當成未來保證。
#
# 區域間的測網、規模尺度、構造背景與目錄長度不同，參數也未必可直接移植。比起尋找一組通用的「標準 EEPAS 參數」，更重要的是說清楚輸入、目標、可用歷史及比較基準。
#
# %% [markdown]
# ```{admonition} 閱讀參數時，先辨認它控制哪一部分
# :class: dropdown
#
# $a_M,b_M,\sigma_M$ 控制目標規模；$a_T,b_T,\sigma_T$ 控制對數等待時間；$b_A,\sigma_A$ 控制空間範圍。$\mu_E$ 是背景係數，$w_i$ 是事件權重。不同版本的符號及正規化可能不同，因此比較數值前須先對齊核函數、單位與門檻。完整核、$\eta$、$\Delta$ 及卷積推導見 {doc}`附錄 D <appendix_d_eepas>`。
# ```
#
# %% [markdown]
# 回到我們一路追問的問題：新的地震究竟帶來什麼資訊？ETAS 把它轉成後續觸發，EEPAS 把它轉成較長尺度的分布，PPE 則保留它對空間活動度的資訊。三種解讀可以從同一目錄開始，卻必須對同一未來目標交出可比較的預報。
#
# %% [markdown]
# 這些模型已經能產生預報，接著需要一套共同的評估方式：觀測到的事件數、位置和規模，是否符合模型原先的預期？這是 {doc}`預報一致性檢驗 <17_testing_consistency>`的起點。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 3–5 節，建立 EEPAS、模型組合與缺失前兆補償的全貌，再回到本章的正規化推導。
# - [Long-range Earthquake Forecasting with Every Earthquake a Precursor According to Scale](https://doi.org/10.1007/s00024-003-2434-9) — David A. Rhoades、Frank F. Evison，2004，*Pure and Applied Geophysics*；全文可能需訂閱。這是 EEPAS 原始論文，重點是如何把尺度關係轉成每個事件對未來地震率的貢獻，而不是先判定哪個事件必然是前兆。
# - [Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123) — Emanuele Biondini、David A. Rhoades、Paolo Gasperini，2023，*Geophysical Journal International*；出版社全文可能需訂閱，[免費機構典藏全文](https://www.earth-prints.org/handle/2122/17084)。將本章公式連到義大利目錄的實際應用，閱讀 PPE、ETAS 與 EEPAS 的比較時，特別留意學習期、測試期和預報時間窗。
#
# - [The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data](https://doi.org/10.3390/e22111264) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2020，*Entropy*；免費開放全文。說明目錄開始前遺漏的事件如何影響 EEPAS，適合延伸本章的 lead time 與時間完整度補償推導。
# - [Long-range earthquake forecasting allowing for aftershocks](https://doi.org/10.1111/j.1365-246X.2008.04083.x) — D. A. Rhoades，2009，*Geophysical Journal International*；[出版社網頁全文](https://academic.oup.com/gji/article/178/1/244/644120)可免費閱讀。閱讀 EEPAS 如何加入預報事件的餘震貢獻，對照本章 EAS 延伸與「降低輸入餘震權重」的不同角色。
#
# - [Application of a long-range forecasting model to earthquakes in the Japan mainland testing region](https://doi.org/10.5047/eps.2010.08.002) — David A. Rhoades（2011），Earth, Planets and Space（免費全文）。本章日本規模級距增益圖的直接來源；摘要及結果明確區分回溯擬合與預定的前瞻測試。
