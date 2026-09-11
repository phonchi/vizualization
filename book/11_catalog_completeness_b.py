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
# # 5. 目錄統計 I：完整度與規模分布
#
# 前一章用條件強度描述「此刻有多容易發生地震」。要把這個想法用在目錄上，還有一個很實際的問題：模型使用的是記錄下來的地震，未必涵蓋地下實際發生的所有事件。
#
# 想像兩段長度相同的目錄，後一段的小地震多出許多。這可能是活動變頻繁，也可能只是測站增加了。若直接把事件數增加解讀成地殼改變，後面的預報會從一開始就追錯問題。規模分布可以幫我們分辨這兩種可能，因為漏測通常先影響小事件。
#
# 我們沿著同一份目錄，先看不同規模的事件各有多少，再問哪些規模可信，最後估計大小地震的相對比例。這個比例用 $b$ 值表示；它也是後面 ETAS 與危害度模型分配規模機率時需要的資訊。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 5.1 先看目錄記錄了什麼
#
# 把地震依規模排序，逐一計算「規模至少這麼大」的事件數，就得到累積規模頻度分布，常簡稱 FMD。同一事件會同時出現在多個累積門檻中，因此圖上的相鄰點不是互相獨立的觀測。
#
# 下面沿用臺灣公開目錄作為讀圖例子，將不同年代的事件數除以各自的年數，畫成年發生率。先看曲線在小規模端是否彎曲，再看較大規模端是否仍有差異。這兩種差異需要不同解釋：前者可能與漏測有關，後者也可能包含序列活動、樣本波動或規模尺度改變。不能只看曲線高低就判定哪個年代比較危險。
#
# 數十年間，目錄的觀測條件也在變。測站、儀器動態範圍、事件定位與規模計算都可能改變。Wang et al.（2015）的回顧整理了臺灣規模統計的觀測背景；本章先利用這個案例學會辨認資料條件，第二部再回到觀測網本身。
#
# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

DM = 0.1            # 台灣 CWA 目錄的規模格距 ΔM
DELTA = DM / 2      # 半格 δ
L10 = np.log(10.0)

cat_long = load_taiwan_catalog()
cat_long = cat_long.assign(M=np.round(cat_long.ML.to_numpy() / DM) * DM)

ERAS = {"1973–1987.5（類比 M_D）": ("1973", "1987-06-01"),
        "1987.6–1991.2（數位 M_D）": ("1987-06-01", "1991-03-01"),
        "1994–2011（連續記錄 M_L）": ("1994", "2012"),
        "2012–2025（觀測網成熟）": ("2012", "2026")}
BINS = np.arange(1.95, 7.05, DM)
CENTERS = np.round(BINS[:-1] + DELTA, 2)


def maxc(m, corr=0.2):
    """最大曲率法（MAXC）＋慣用的 0.2 保守修正。"""
    counts, _ = np.histogram(np.asarray(m, float), bins=BINS)
    return round(float(CENTERS[np.argmax(counts)]) + corr, 2)


mags = np.arange(2.0, 7.01, DM)
fig = go.Figure()
era_mc = {}
for (label, (t0, t1)), color in zip(ERAS.items(), PALETTE):
    sub = cat_long[(cat_long.time >= t0) & (cat_long.time < t1)]
    yrs = (sub.time.max() - sub.time.min()).days / 365.25
    rate = np.array([(sub.M >= m - 1e-9).sum() / yrs for m in mags])
    era_mc[label] = maxc(sub.M.to_numpy())
    fig.add_trace(go.Scatter(x=mags, y=rate, mode="lines", name=label,
                             line=dict(color=color, width=2)))
    j = int(np.argmin(np.abs(mags - era_mc[label])))
    fig.add_trace(go.Scatter(x=[mags[j]], y=[rate[j]], mode="markers",
                             showlegend=False,
                             marker=dict(color=color, size=11,
                                         symbol="diamond",
                                         line=dict(color="white", width=1))))
apply_layout(fig,
             title="分年代累積 FMD 與各年代 Mc（菱形＝MAXC＋0.2："
                   + "、".join(f"{v:.1f}" for v in era_mc.values()) + "）",
             xaxis_title="規模 ML", yaxis_title="年發生率 N（M ≥ ML）",
             yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 讀這張圖時，要把「曲線往下彎」與「確定不完整的門檻」分開。前者是資料呈現的形狀，後者還需要一套估計方法。圖中的菱形採用 MAXC 加修正量，只是初步摘要；稍後會看到其他方法可能給出不同結果。
#
# 這份公開目錄本身有收錄下限。即使後來的觀測網能記錄更小的事件，我們也無法只靠這份擷取後的資料估計那一段能力。較大規模事件通常較容易偵測，但不同年代曲線相近，也不等於規模尺度與資料品質已完全一致。
#
# ### 同樣叫「規模」，數字不一定能直接混用
#
# $M_L$ 主要由局部地震波振幅估計，$M_w$ 則由地震矩換算；它們的物理定義與校準方式不同。將兩種尺度混在一份規模分布中，可能改變曲線的斜率。跨目錄比較時，應先確認使用的尺度、適用規模區間與轉換來源，不能只把欄位名稱統一。
#
# 尺度轉換也會改變斜率。若在研究範圍內有 $M_2=a_sM_1+d_s$、$a_s>0$，同一條規模頻度曲線換成 $M_2$ 表示時，其斜率變為 $b_2=b_1/a_s$。轉換式本身的誤差仍須另外傳遞；兩個方向分別做的迴歸也不必互為倒數。換算過程見附錄。報告 $b$ 值時，應同時寫出規模尺度。
#
# %% [markdown]
# ## 5.2 用一個斜率描述大小地震的比例
#
# 許多足夠完整的目錄，在累積次數的對數圖上有一段近似直線。Gutenberg–Richter 律把這個規律寫成
#
# $$N(\ge m)=10^{a-bm}.$$
#
# $b$ 決定線有多陡。若 $b=1$，規模門檻提高一級，事件數約剩十分之一；若 $b$ 較大，大事件相對更少。這描述的是整個規模分布，不是下一個事件一定有多大。$a$ 則與整段資料的活動量、期間長度及研究區域有關。
#
# 為了把次數接回前一章的標記點過程，先選定輸入門檻 $m_0$，將累積次數除以 $N(\ge m_0)$。得到的比例就是一個事件超過規模 $m$ 的機率：
#
# $$P(M\ge m)=10^{-b(m-m_0)}=e^{-\beta(m-m_0)},\qquad \beta=b\ln10.$$
#
# 因此，門檻以上的「超額規模」$M-m_0$ 服從指數分布，其密度為
#
# $$s(m)=\beta e^{-\beta(m-m_0)},\qquad m\ge m_0.$$ (eq:gr-density)
#
# 這裡的指數分布是規模分布，和前面 Poisson 過程的指數等待時間是不同隨機變數。模型可以同時使用兩者，也可以只用其中一種。ETAS 就保留這個簡單的規模分布，讓事件時間具有叢集結構。
#
# ```{admonition} 次數與年率怎麼換算？
# :class: dropdown
#
# 若期間長度為 $T$ 年，則年率 $\nu(\ge m)=N(\ge m)/T$。在對數圖上，這會把截距減去 $\log_{10}T$，不改變 $b$。若改以門檻 $m_0$ 為原點，寫成 $N(\ge m)=10^{a_0-b(m-m_0)}$，則 $a_0=a-bm_0$。比較截距以前，時間、面積、規模尺度與門檻都要一致。
# ```
#
# GR 律是一段規模區間內的統計近似。小規模端可能受漏測影響，極大規模端也可能需要截斷或其他尾端模型。下一步是先決定哪些事件能用來估計這段直線。
#
# %% [markdown]
# ## 5.3 從哪個規模開始相信這條分布？
#
# 完整度規模 $M_c$ 表示目錄從哪個規模以上可近似視為收錄齊全。它是需要估計的資料性質；$m_0$ 是分析者設定的模型輸入門檻；$m_T$ 是預報要回答的目標門檻。三者不能互換。例如資料可能從某個小規模起較完整，模型卻選擇較高的輸入門檻，而預報只關心更大的事件。
#
# 完整度會隨空間、年代及主震後時間改變。偏遠地區測站較疏，可能漏掉更多小震；大震後波形重疊，也會讓原本能辨認的小事件暫時消失。若研究期間跨越這些變化，單一 $M_c$ 至多是方便的近似。
#
# 估計 $M_c$ 的方法，其實是在用不同證據回答同一個觀測問題。MAXC 尋找非累積直方圖的峰值，再使用經驗修正；它容易理解，但峰值會受分箱與偵測曲線形狀影響。$b$ 值穩定法逐步抬高門檻，找估計斜率不再明顯改變的區間。分布適合度方法則比較門檻以上的經驗分布與擬合的指數分布，檢查兩者是否相容。
#
# 下面把三種方法放在同一張圖上。它們使用同一批事件，因此差異不來自地震活動，而來自判準。圖中的 KS 型方法用模擬校準比較，不應把「沒有拒絕」解讀為證明沒有漏測。大樣本能發現很小的偏離；小樣本則可能根本沒有足夠能力區分模型。
#
# 圖上的三個標示都是各方法選出的完整度規模門檻；KS 方法旁的數值也是門檻，不是檢定的 p 值。
#
# %% tags=["remove-input"]
def b_exact(m, mc, dm=DM):
    """離散精確式（Tinti & Mulargia 1987）＋ Shi & Bolt (1982) 標準差。"""
    m = np.asarray(m, float)
    m = m[m >= mc - 1e-9]
    b = np.log1p(dm / (m.mean() - mc)) / (L10 * dm)
    return b, L10 * b ** 2 * np.sqrt(m.var(ddof=1) / len(m)), len(m)


def mc_bstability(m, grid, k_win=5):
    """b 值穩定度判準（Cao & Gao 2002；Woessner & Wiemer 2005）。"""
    bs, ss = np.array([b_exact(m, c)[:2] for c in grid]).T
    for i, c in enumerate(grid):
        if i + k_win < len(grid) and \
                abs(bs[i + 1:i + 1 + k_win].mean() - bs[i]) < ss[i]:
            return round(float(c), 2), bs, ss
    return float(grid[-1]), bs, ss


def mc_ks(m, grid, rng, nsim=400, nmax=2000, p_th=0.1):
    """Clauset/Mizrahi 式的模擬 p 值程序；nmax 為每個門檻的樣本上限。"""
    for c in grid:
        s = m[m >= c - 1e-9]
        if len(s) < 100:
            break
        if len(s) > nmax:
            s = rng.choice(s, nmax, replace=False)
        n, k = len(s), np.round((s - c) / DM).astype(int)
        r = 10 ** (-b_exact(s, c)[0] * DM)
        grid_k = np.arange(int(k.max()) + 1)
        f_mod = 1 - r ** (grid_k + 1)
        d_obs = np.abs(np.bincount(k, minlength=len(grid_k)).cumsum() / n
                       - f_mod).max()
        sim = np.sort(rng.geometric(1 - r, size=(nsim, n)) - 1, axis=1)
        f_sim = np.stack([np.searchsorted(row, grid_k, side="right")
                          for row in sim]) / n
        if float((np.abs(f_sim - f_mod).max(axis=1) >= d_obs).mean()) >= p_th:
            return round(float(c), 2)
    return float(grid[-1])


m94 = cat_long[(cat_long.time >= "1994") & (cat_long.time < "2012")].M.to_numpy()
grid_mc = np.round(np.arange(2.0, 4.51, DM), 2)
mc_a = maxc(m94)
mc_b, b_curve, s_curve = mc_bstability(m94, grid_mc)
mc_c = mc_ks(m94, grid_mc, np.random.default_rng(11))

counts, _ = np.histogram(m94, bins=BINS)
fig = go.Figure(go.Bar(x=CENTERS, y=counts, marker_color=ACCENT,
                       opacity=0.75, name="非累積 FMD", showlegend=False))
for x, name, color in [(mc_a, f"MAXC＋0.2 = {mc_a:.1f}", PALETTE[1]),
                       (mc_b, f"b 值穩定度 = {mc_b:.1f}", PALETTE[2]),
                       (mc_c, f"KS 模擬法：Mc = {mc_c:.1f}", PALETTE[3])]:
    fig.add_vline(x=x, line_color=color, line_dash="dash")
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name=name,
                            line=dict(color=color, dash="dash")))
apply_layout(fig,
             title=f"三種 Mc 估計法（台灣 1994–2011，N = {len(m94):,}）："
                   f"彼此差 {max(mc_a, mc_b, mc_c) - min(mc_a, mc_b, mc_c):.1f}"
                   f" 個規模單位",
             xaxis_title="規模 ML", yaxis_title="事件數",
             xaxis_range=[1.9, 5.5], yaxis_type="log", hovermode="x",
             legend=dict(orientation="h", x=0, y=1.02, yanchor="bottom"),
             margin=dict(l=60, r=20, t=100, b=40))
fig

# %% [markdown]
# 三條門檻線的差異，是完整度不確定性的一部分。方法對直方圖峰值、斜率穩定與分布形狀各有不同敏感度，不能只挑出最符合預期的一條。
#
# 接著把門檻逐格抬高，重新估計 $b$。這張曲線能把「門檻選在哪裡」與「估計有多穩定」放在一起讀，也能看出提高門檻所付出的樣本數代價。
#
# %% tags=["remove-input"]
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=np.r_[grid_mc, grid_mc[::-1]],
    y=np.r_[b_curve + s_curve, (b_curve - s_curve)[::-1]],
    fill="toself", fillcolor="rgba(42,120,214,0.20)", line=dict(width=0),
    hoverinfo="skip", name="±σ（Shi & Bolt）"))
fig.add_trace(go.Scatter(x=grid_mc, y=b_curve, mode="lines+markers",
                         name="b(Mc)", line=dict(color=ACCENT),
                         marker=dict(size=5)))
fig.add_vrect(x0=mc_b, x1=mc_b + 0.5, fillcolor=PALETTE[2], opacity=0.13,
              line_width=0, annotation_text="穩定度判準的 L = 0.5 視窗",
              annotation_position="top left")
for x, color in [(mc_a, PALETTE[1]), (mc_b, PALETTE[2]), (mc_c, PALETTE[3])]:
    fig.add_vline(x=x, line_color=color, line_dash="dash")
i_b = int(np.argmin(np.abs(grid_mc - mc_b)))
apply_layout(fig,
             title=f"b(Mc) 曲線與穩定度判準帶（1994–2011）：低門檻端由 "
                   f"{b_curve[0]:.2f} 爬升到平台 {b_curve[i_b]:.2f}",
             xaxis_title="截取門檻 Mc", yaxis_title="估計的 b 值",
             hovermode="x")
fig

# %% [markdown]
# 若低門檻漏掉許多小事件，留下來的事件平均規模會偏大，估計的 $b$ 往往偏低。門檻提高後可能出現平臺，但平臺本身仍不足以證明分布完全正確。右側只剩少數事件時，誤差也會增加。
#
# 因此，門檻選擇不是越低或越高越好。應選擇有資料支持的完整區間，再報告合理門檻變動下的結果。即使這一步做好了，目錄中規模通常只有一位或兩位小數，估計公式還需要配合它的記錄方式。
#
# %% [markdown]
# ## 5.4 為什麼平均規模可以估計 $b$？
#
# 在連續、完整的指數模型下，超額規模的平均值是 $1/\beta$。因此，樣本平均超出門檻越多，就表示分布尾端較重，$b$ 越小。最大概似估計把這個直覺寫成
#
# $$\widehat b_{\rm cont}=\frac{1}{\ln10\,(\bar m-m_0)}.$$
#
# 這裡的 $m_0$ 必須是實際納入樣本的下限。若選在估計的完整度上，就是 $m_0=\widehat M_c$。公式本身很短，真正重要的是分母的含義，以及資料是否符合完整、連續規模的假設。概似微分的過程移到 {doc}`appendix_b_catalog`，不影響這裡的讀圖。
#
# 假如規模只記到 0.1，門檻那一格便有非零機率，連續密度不再精確描述資料。直接把格點值當連續觀測，會把平均超額規模估得過小，因而使 $b$ 偏高。這與漏測小事件造成的方向相反，所以兩種偏差甚至可能互相抵銷，產生一個看似合理的答案。
#
# 下圖特意用已知 $b$ 的模擬目錄，把離散化和不完整放在同一個例子中。先觀察低門檻端兩種估計是否都偏離真值，再看高門檻端的差異。它能區分「資料少記了事件」和「公式沒有配合資料格式」這兩種問題。
#
# %% tags=["remove-input"]
from scipy.stats import norm

rng = np.random.default_rng(42)
B_TRUE = 1.0
m_raw = rng.exponential(np.log10(np.e) / B_TRUE, 400_000) + 2.0
m_raw = np.round(m_raw / DM) * DM                     # 分箱到 0.1
obs = m_raw[rng.random(m_raw.size)
            < norm.cdf(m_raw, loc=3.0, scale=0.25)]   # 偵測率平滑上升

cutoffs = np.round(np.arange(2.6, 4.41, DM), 2)
b_aki = np.array([np.log10(np.e) / (obs[obs >= c - 1e-9].mean() - c)
                  for c in cutoffs])
b_exa = np.array([b_exact(obs, c)[0] for c in cutoffs])

fig = go.Figure()
fig.add_trace(go.Scatter(x=cutoffs, y=b_aki, mode="lines+markers",
                         name="Aki 連續式（未修 binning）",
                         line=dict(color=PALETTE[1]), marker=dict(size=5)))
fig.add_trace(go.Scatter(x=cutoffs, y=b_exa, mode="lines+markers",
                         name="離散精確式", line=dict(color=ACCENT),
                         marker=dict(size=5)))
fig.add_hline(y=B_TRUE, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"真值 b = {B_TRUE:.1f}")
apply_layout(fig,
             title=f"合成不完整＋分箱目錄（N = {len(obs):,}）：Aki 高原 "
                   f"{b_aki.max():.3f}，精確式高原 {b_exa.max():.3f}",
             xaxis_title="截取門檻（假設的 Mc）", yaxis_title="估計的 b 值",
             hovermode="x")
fig

# %% [markdown]
# 低門檻端的偏差來自偵測機率，即使用離散資料的正確公式，也無法自動補回沒有記到的事件。門檻足夠高後，兩條曲線的差異才主要反映離散化。
#
# 圖中真值是模擬時設定的參數，因此可以衡量偏差；真實目錄沒有這條參考線。這也是先用模擬理解方法，再回頭分析觀測資料的理由。
#
# %% [markdown]
# ## 5.5 讓估計公式配合規模格距
#
# 令記錄格距為 $\Delta M$，納入的最小格點為 $m_0$。在離散 GR 模型下，規模只可能是 $m_0,m_0+\Delta M,m_0+2\Delta M,\ldots$。每提高一格，機率乘上同一個比例，所以這是一個幾何分布。
#
# 令半格寬 $\delta=\Delta M/2$，對應的最大概似估計為
#
# $$\widehat b=\frac{1}{2\delta\ln10}\ln\left(1+\frac{2\delta}{\bar m-m_0}\right).$$ (eq:b-exact)
#
# 這個公式假設格距固定、門檻以上完整，且規模分布符合離散 GR。若所有事件恰好落在最小格點，分母為零，表示無法從資料估得有限的斜率，不能把它當成一般的數值誤差略過。格距趨近零時，它回到前節的連續式。
#
# 常見的 Utsu 半格修正，把連續式的分母換成 $\bar m-m_0+\Delta M/2$。它是小格距時的近似，而離散式直接使用格距資訊。Tinti 與 Gasperini（2024）比較了這些估計及區間，提醒我們：公式的差異必須連同樣本大小與記錄格距一起判斷。
#
# ```{admonition} 為什麼離散式裡的門檻是最小格點？
# :class: dropdown
#
# 若四捨五入後保留 $m\ge m_0$，被保留的原始連續值可能從 $m_0-\Delta M/2$ 開始。對每一格積分並條件化到保留區間後，格點索引服從幾何分布。離散公式已包含這個結構，不應再任意把 $m_0$ 減半格，否則會重複修正。不同擷取規則須重新確認其支撐集。
# ```
#
# %% [markdown]
# ## 5.6 一個估計值還不夠
#
# 從同一個分布反覆抽樣，得到的 $b$ 不會完全相同。在理想的連續指數模型下，大樣本標準誤約為 $b/\sqrt N$：樣本數增加四倍，隨機誤差才約減半。Shi–Bolt 型近似則利用樣本規模的離散程度估計誤差。兩者都不會自動包含門檻選擇、漏測或模型設定造成的系統差異。
#
# $b$ 是平均超額規模的非線性函數，因此平均規模上的對稱區間，換成 $b$ 後通常不對稱。這在小樣本時尤其需要注意。下圖從花蓮序列抽取子樣本，再重抽樣估計 $b$，比較所得分布與不同區間近似。圖中的重抽樣把觀測規模視為可交換樣本；若存在規模相關、時變完整度或選樣效應，這種 bootstrap 不會包含所有不確定性。
#
# 讀圖時先找分布中心，再看左右尾端，而不要只比較哪一條誤差棒最短。短區間可能來自額外假設，不一定代表資料更有資訊。離散區間的代數與適用條件整理於附錄。
#
# %% tags=["remove-input"]
cat24 = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
main = cat24.loc[cat24.ML.idxmax()]
win = cat24[(cat24.time > main.time + pd.Timedelta(days=1))
            & (cat24.time <= main.time + pd.Timedelta(days=30))]
pool = np.round(win.ML.to_numpy() / DM) * DM
pool = pool[pool >= 3.5 - 1e-9]

rng_bs = np.random.default_rng(2024)
samp = rng_bs.choice(pool, 60, replace=False)
b_hat, sig_sb, n_sub = b_exact(samp, 3.5)
boot = np.array([b_exact(rng_bs.choice(samp, n_sub, replace=True), 3.5)[0]
                 for _ in range(2000)])

c_val = 10 ** (DM * b_hat)
sN = np.sqrt(c_val / n_sub)
b1 = np.log((c_val + sN) / (1 + sN)) / (L10 * DM)
b2 = np.log((c_val - sN) / (1 - sN)) / (L10 * DM)

fig = go.Figure(go.Histogram(x=boot, xbins=dict(size=0.02),
                             marker_color=ACCENT, opacity=0.8,
                             name="bootstrap 分布"))
fig.add_vline(x=b_hat, line_color=QUAKE_COLOR,
              annotation_text=f"點估計 {b_hat:.3f}")
for x, dash, color, name in [(b_hat - sig_sb, "dot", PALETTE[1],
                              "Shi & Bolt ±σ"),
                             (b_hat + sig_sb, "dot", PALETTE[1], ""),
                             (b1, "dash", PALETTE[2], f"b1 = {b1:.3f}"),
                             (b2, "dash", PALETTE[2], f"b2 = {b2:.3f}")]:
    fig.add_vline(x=x, line_dash=dash, line_color=color, annotation_text=name)
apply_layout(fig,
             title=f"bootstrap 的 b 值分布（0403 花蓮餘震子樣本 N = {n_sub}）："
                   f"σ1 = {b_hat - b1:.3f} 對 σ2 = {b2 - b_hat:.3f}",
             xaxis_title="估計的 b 值", yaxis_title="次數", hovermode="x")
fig

# %% [markdown]
# 這個例子呈現非線性換算後的區間形狀。比較兩段目錄的 $b$ 時，應同時考慮各自的區間與共同資料處理方式，不能只憑兩個點估計一高一低就認定活動機制改變。
#
# 目前仍有一個問題未處理：主震後完整度變化最快，卻正好是最想即時估計規模分布的時候。下一種方法就是為了減輕這種影響而提出。
#
# %% [markdown]
# ## 5.7 從規模差減輕短期漏測的影響
#
# 主震剛發生時，小餘震可能被較強波形掩蓋，完整度門檻隨時間下降。此時用整段序列的一個 $M_c$，容易把偵測能力變化誤當成 $b$ 變化。b-positive 等方法改看相鄰或適當配對事件的規模差，希望選出較不受當時偵測門檻影響的資訊。
#
# 其數學起點很簡單：若兩個規模各自獨立服從同一個指數分布，它們的差是左右對稱的 Laplace 分布；保留大於指定正門檻的差，再減去該門檻，其尾端仍是指數形狀。這讓規模差也能用來估計同一個 $b$。
#
# 但目錄被漏測篩選後，觀測配對不一定仍滿足理想獨立假設。正差法、絕對差法與尋找下一個足夠大事件的配對法，保留的資訊各不相同。修剪門檻提高，可以減少某些不完整影響，卻也會減少可用配對。沒有一種取差方式能保證完全繞開資料品質問題。
#
# 下面用時變偵測能力的合成序列比較三種取差方式。橫軸是保留差值的門檻，並不是完整度規模 $M_c$。先看曲線在門檻增加時是否接近真值，再留意高門檻端的波動。Tinti 與 Gasperini（2024）的比較也以這種方式區分估計公式、配對規則與不完整情境。
#
# 圖中正差和負差都取自時間順序上相鄰的事件，負差以其大小表示。不重疊配對則將第 1、2 個事件配成一對，第 3、4 個配成下一對，再取各對的絕對差，每個事件只用一次。本圖沒有實作前面提到的「尋找下一個足夠大事件」配對法。
#
# ```{admonition} 相鄰差值為什麼不算獨立樣本？
# :class: dropdown
#
# $m_{i+1}-m_i$ 與 $m_{i+2}-m_{i+1}$ 共用 $m_{i+1}$，即使原規模獨立，兩個差仍相關。取絕對值或只留正差後，相關結構還會改變。因此不能把所有保留差值的數量直接當成獨立樣本數，照搬原本的標準誤。附錄從理想指數模型說明這個差別。
# ```
#
# %% tags=["remove-input"]
def b_trim(d, dc):
    """trimmed 差分估計式。"""
    d = np.asarray(d, float)
    d = d[d >= dc - 1e-9]
    return np.log1p(DM / (d.mean() - dc)) / (L10 * DM)


rng_af = np.random.default_rng(7)
N_SYN, C_OM, T_OM, M_MAIN, M0 = 200_000, 0.01, 30.0, 5.6, 2.0
u_om = rng_af.random(N_SYN)
t_af = np.sort(C_OM * ((1 + T_OM / C_OM) ** u_om - 1))     # Omori（p = 1）
k_af = rng_af.geometric(1 - 10 ** (-B_TRUE * DM), size=N_SYN) - 1
m_af = np.round(M0 + k_af * DM, 2)
mc_t = np.maximum(M_MAIN - 4.5 - 0.76 * np.log10(np.maximum(t_af, 1e-4)), M0)
m_obs = m_af[rng_af.random(N_SYN) < norm.cdf(m_af, loc=mc_t, scale=0.5)]

b_plain = b_exact(m_obs, M0)[0]
dd = np.diff(m_obs)                                    # 取法 (A)：相鄰
half = len(m_obs) // 2 * 2
d_pair = np.abs(m_obs[1:half:2] - m_obs[0:half:2])     # 取法 (B)：配對
trims = np.round(np.arange(0.1, 1.01, 0.1), 2)

fig = go.Figure()
for name, arr, color in [("相鄰正差", dd[dd > 0], PALETTE[0]),
                         ("相鄰負差的大小", -dd[dd < 0], PALETTE[1]),
                         ("不重疊配對絕對差", d_pair, PALETTE[2])]:
    fig.add_trace(go.Scatter(x=trims, y=[b_trim(arr, dc) for dc in trims],
                             mode="lines+markers", name=name,
                             line=dict(color=color), marker=dict(size=6)))
fig.add_hline(y=B_TRUE, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"真值 b = {B_TRUE:.1f}")
fig.add_hline(y=b_plain, line_dash="dot", line_color="#888888",
              annotation_text=f"直接用規模的精確式 = {b_plain:.3f}")
apply_layout(fig,
             title=f"trimming 門檻掃描（合成餘震序列，觀測到 {len(m_obs):,} "
                   f"個事件）：三種取差方式一起收斂",
             xaxis_title="trimming 門檻 ΔM'c", yaxis_title="估計的 b 值",
             hovermode="x")
fig

# %% [markdown]
# 這個特定漏測情境中，修剪規模差能減輕偏差；這是合成資料上的方法比較，不能直接推論所有地區的目錄都適用同一門檻。高門檻下若只剩很少配對，估計也可能重新變得不穩定。
#
# 至此我們已經看過完整度、格距、樣本量與配對規則。還有一個更根本的檢查：資料真的支持用一條 GR 分布描述嗎？
#
# %% [markdown]
# ## 5.8 先檢查分布，再討論變化
#
# 估計公式即使精確，所用的模型仍可能不適合資料。若不同時段、不同構造區或不同規模尺度被混在一起，整體規模分布可能彎曲；此時單一 $b$ 只能是摘要，不宜把它解讀成共同的物理參數。
#
# 分布適合度檢查會把觀測累積分布與擬合模型比較。因為參數通常從同一批資料估得，直接使用「參數已知」的 KS 臨界值並不正確；可採用相應的 Lilliefors 型校準，或模擬、重估參數後再比較統計量。若連 $M_c$ 都由資料挑選，校準程序也應反映這個選擇。
#
# 這和「兩段目錄的 $b$ 是否不同」是不同問題。前者問分布形式是否合理，後者問差異是否超出隨機波動。若反覆掃描許多時間窗，只展示差異最大的那一段，還會有多重比較與事後挑選問題。前瞻性預報需要在看到結果前固定選窗和判定規則。
#
# 因此，比較 $b$ 值時應一起交代目錄期間、研究區域、規模尺度、格距、門檻、樣本數和是否除叢。這些條件界定了這個數字所描述的資料與分析方式。
#
# %% [markdown]
# 現在可以回到一開始的問題：後一段目錄小震變多，究竟是觀測變好了，還是規模分布改變？我們有了分開檢查的工具，也知道僅憑一個 $b$ 值無法回答全部問題。
#
# 接下來保留這些對資料的認識，將事件放回時間與空間。同樣多、同樣規模分布的地震，可能均勻散開，也可能集中在少數序列。{doc}`12_clustering_laws` 從這個差異出發，說明餘震如何衰減、最大事件為何不容易預測，以及切分序列會如何反過來改變規模統計。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Estimate magnitude of completeness](https://seismostats.readthedocs.io/latest/user/estimate_mc.html) — SeismoStats 開發團隊；免費官方文件。用公式與短例子比較 MAXC、KS 與 b 值穩定法，閱讀時可對照本章對門檻選擇及方法假設的討論。
# - [Estimate b-value](https://seismostats.readthedocs.io/latest/user/estimate_b.html) — SeismoStats 開發團隊；免費官方文件。整理傳統與差分式估計法，適合查閱規模離散化、完整度門檻與各方法適用條件。
# - [The estimation of b-value of the frequency–magnitude distribution and of its 1σ intervals from binned magnitude data](https://doi.org/10.1093/gji/ggae159) — S. Tinti、P. Gasperini，2024，*Geophysical Journal International*；[免費出版版全文與補充資料](https://cris.unibo.it/handle/11585/980514)。深入比較分箱規模資料的估計式與不確定度，適合讀完本章連續、離散估計推導後，進一步理解哪些修正有統計依據。
#
# - [b-Values Observations in Taiwan: A Review](https://doi.org/10.3319/TAO.2015.04.28.01%28T%29) — Jeen-Hwa Wang、Kou-Cheng Chen、Pei-Ling Leu、Jeng-Hsin Chang，2015，*Terrestrial, Atmospheric and Oceanic Sciences*；[免費全文](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1361-2015042801t)。整理臺灣 b 值的時空變化、構造背景及估計研究，適合對照本章的臺灣目錄例子；閱讀前兆相關解釋時，仍需分辨回溯觀察與前瞻檢驗。
#
# - [B-Positive: A Robust Estimator of Aftershock Magnitude Distribution in Transiently Incomplete Catalogs](https://doi.org/10.1029/2020JB021027) — Nicholas J. van der Elst，2021，*Journal of Geophysical Research: Solid Earth*。規模差方法的原始來源，重點閱讀短期完整度的假設與限制；出版頁提供全文取得資訊。
# - [GDMS 地球物理資料管理系統](https://gdms.cwa.gov.tw/) — 中央氣象署。本章臺灣目錄圖的資料來源；下載期間、尺度與篩選條件保留於原始 notebook。
