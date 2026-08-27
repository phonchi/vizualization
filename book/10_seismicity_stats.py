# %% [markdown]
# # 10. 地震目錄統計的進階課
#
# {doc}`第 5 章 <05_seismic>`你已經算過 GR 律的 b 值、畫過 Omori
# 衰減，也在章末聽過
# 一句警告：「b 值分析的困難不在算式，而在它對取樣選擇非常敏感。」
# 這一章就是把那句警告展開。所有預報模型——ETAS、EEPAS、PSHA——
# 的輸入都是地震目錄與它的統計參數；地基歪一度，上層建築就斜一片。
#
# 這章要回答四個問題：目錄從哪一筆開始可以信（Mc）？b 值的誤差
# 到底多大（以及 b-positive 是什麼）？Omori 律有哪些教科書沒講的
# 皺褶？最大餘震能不能預測（Båth 定律）？最後談「除叢」——一個
# 看似技術性、實則影響深遠的選擇。
#
# 本章的示意圖大多使用一份新的資料：**1973–2025 年的台灣長期目錄**
# （中央氣象署公開目錄，約 35 萬筆、ML≥2，由本站 `gdms_toolkit`
# 提供載入函式）。五十年的目錄才裝得下這一章要講的坑。
#
# ## 10.1 完整規模 Mc：一切統計的地基
#
# 完整規模 $M_c$（magnitude of completeness）是「目錄從這個規模
# 以上大致收錄齊全」的門檻。第 5 章提過它一次；這裡要講清楚的是：
# **Mc 不是一個數字，而是一個隨時間、空間、甚至主震後時刻變動的場。**
#
# 先看時間。台灣目錄的偵測能力在五十年間有幾次跳躍式的提升——
# 1973 年儀器目錄起點、1991 年改用模擬 Wood–Anderson 的 $M_L$、
# 1994 年即時資料由觸發式改為連續記錄（年偵測數 4 千 → 2 萬）、
# 2012 年觀測網升級成熟（再翻倍到 4 萬，多為 M2 以下微震）。把長期目錄按年代切開，
# 各畫一條規模–頻率分布（FMD），差異一目了然：

# %% tags=["remove-input"]
import plotly.io as pio
pio.renderers.default = "notebook_connected"

# %% tags=["hide-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

cat_long = load_taiwan_catalog()

eras = {
    "1973–1987（類比延時規模）": ("1973", "1988"),
    "1994–2011（連續記錄）": ("1994", "2012"),
    "2012–2025（觀測網成熟）": ("2012", "2026"),
}
mags = np.arange(2.0, 7.01, 0.1)
fig = go.Figure()
for (label, (t0, t1)), color in zip(eras.items(), PALETTE):
    sub = cat_long[(cat_long.time >= t0) & (cat_long.time < t1)]
    yrs = (sub.time.max() - sub.time.min()).days / 365.25
    n_per_yr = [(sub.ML >= m).sum() / yrs for m in mags]
    fig.add_trace(go.Scatter(x=mags, y=n_per_yr, mode="lines+markers",
                             name=label, line=dict(color=color),
                             marker=dict(size=4)))
apply_layout(fig, title="不同年代的規模–頻率分布：Mc 隨觀測網進步下降",
             xaxis_title="規模 ML", yaxis_title="年發生率 N（ML ≥ M）",
             yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 三條曲線的高規模端（M≥4.5）大致重合——大地震誰都漏不掉，年發生
# 率是穩定的。分歧全在低規模端：曲線偏離 GR 直線、往下彎的位置，
# 就是各年代的 Mc。1970 年代約在 2.5 以上，近年已壓到 2 附近
# （本目錄的收錄下限就是 ML 2.0，再低看不到）。
# **這也表示：任何跨年代的地震活動度比較，若不先固定一個共同的
# Mc（取各年代的最大值），看到的「地震變多」多半只是儀器變好。**
#
# 空間上同樣如此。台灣本島測站密集，Mc 可低到 1.5–2.0；外海
# 測站涵蓋差，Mc 高達 2.5–3.2（Chan & Wu 2013）——最低與最高
# 差達 1.7 個規模單位，對應事件數差約 50 倍。對全目錄套用單一 Mc，等於讓
# 外海的不完整資料污染整體統計。
#
# 估計 Mc 的常用方法有三類（SeismoStats 等現代軟體都有實作）：
#
# 1. **最大曲率法（MAXC）**：取非累積 FMD 的峰值規模，再加上
#    保守修正 0.2——模擬研究證實這個「加 0.2」的老慣例確有道理；
# 2. **b 值穩定法**：利用「Mc 取太低會低估 b」的性質，逐步抬高
#    門檻直到 b 值估計趨於穩定；
# 3. **擬合優度法（KS 距離）**：找最低的門檻使門檻以上的規模
#    分布通過指數分布檢定。
#
# 沒有哪個方法全勝；認真的研究會同時報告多種方法的結果與差異。
#
# ## 10.2 b 值的不確定性：從 Aki 公式到 b-positive
#
# 第 5 章用的 Aki (1965) 最大概似式 $b = \log_{10}e\,/\,(\bar{M}-M_c)$
# 是在「規模連續、目錄完整」兩個假設下推導的。實際目錄兩個假設
# 都不成立，各埋了一個方向相反的偏差：
#
# - **規模是離散的**（記到 0.1）：$\bar{M}$ 被系統性壓低 →
#   b 被**高估**。修正有 Utsu 半格修正（把 $M_c$ 下修半個 bin），
#   以及理論上正確的離散精確式（Tinti & Mulargia 1987）——Utsu
#   其實只是精確式的二階近似，bin 到 0.5（如震度換算的規模）時
#   會失效。
# - **目錄是不完整的**：$M_c$ 取低了，小事件缺漏、分布尾巴被削 →
#   b 被**低估**。
#
# 兩個偏差有時剛好抵銷，讓未修正的公式「看起來很準」——這是
# 假象，不是驗證（Tinti & Gasperini 2024）。
#
# 不確定度也要升級。Aki 的 $\sigma_b = b/\sqrt{N}$ 過於樂觀，
# 實務標準是 Shi & Bolt (1982) 用資料實際離散度算的標準差；而且
# b 值估計的分布**天生不對稱**（右尾較長），嚴謹的區間估計不該寫成
# 對稱的 $b\pm\sigma$。拿第 5 章的花蓮前後對照重算一次，這次把
# 離散修正與誤差一起給：

# %% tags=["hide-input"]
def b_exact(m, mc, dm=0.1):
    """離散精確式（Tinti & Mulargia 1987）＋ Shi & Bolt (1982) 標準差。"""
    m = np.asarray(m[m >= mc], dtype=float)
    p = 1 + dm / (m.mean() - mc)
    b = np.log(p) / (np.log(10) * dm)
    sigma = np.log(10) * b ** 2 * np.sqrt(m.var(ddof=1) / len(m))
    return b, sigma, len(m)

cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
main = cat.loc[cat.ML.idxmax()]
for label, sel in [("主震前（3/1–4/2）", cat.time < main.time),
                   ("餘震期（30 天）", (cat.time >= main.time) &
                    (cat.time <= main.time + pd.Timedelta(days=30)))]:
    b, s, n = b_exact(cat[sel].ML.to_numpy(), mc=3.5)
    print(f"{label}：N = {n:4d}，b = {b:.2f} ± {s:.2f}")

# %% [markdown]
# 兩個 b 值的差距與各自的誤差相比如何？現在你可以自己判斷第 5 章
# 那個「主震前後 b 值有差」的觀察，證據力到底有多強。文獻給這類
# 宣稱設了兩道正式關卡：先用 Lilliefors 檢定確認規模分布真的是
# 指數（GR 律這個前提本身可以檢定！），再檢定 b 值變化是否顯著
# 超過估計不確定度。有一篇論文標題就叫〈How to be fooled searching
# for significant variations of the b-value〉（Marzocchi et al.
# 2020）——光是標題就值得抄在實驗室牆上。
#
# 還有一個更陰險的問題：**大地震剛過後的短期不完整（STAI）**。
# 主震後幾小時到幾天，小餘震被大事件的波形淹沒，Mc 短暫飆高，
# 這段期間混進統計就會拉低 b 值。van der Elst (2021) 的
# **b-positive** 方法給了一個漂亮的迴避：不用規模本身，改用
# **相鄰事件的規模差**——若規模服從指數分布，規模差服從 Laplace
# 分布，取正差再套最大概似即可估 b；而規模差對「完整度隨時間
# 變化」遠不敏感。用一個合成實驗看它的表現：

# %% tags=["hide-input"]
rng = np.random.default_rng(42)
b_true, mc_true = 1.0, 3.0
m_all = rng.exponential(np.log10(np.e) / b_true, 300_000) + 2.0
# 偵測機率隨規模平滑上升（模擬不完整目錄）
from scipy.stats import norm
obs = m_all[rng.random(m_all.size) < norm.cdf(m_all, loc=mc_true, scale=0.25)]

cutoffs = np.arange(2.6, 4.41, 0.1)
b_aki = [np.log10(np.e) / (obs[obs >= c].mean() - c) for c in cutoffs]

d = np.diff(obs)          # 觀測順序下的相鄰規模差
dpos = d[d > 0]
b_pos = {dc: np.log10(np.e) / (dpos[dpos >= dc].mean() - dc)
         for dc in (0.2, 0.6)}   # 兩種修剪門檻

fig = go.Figure()
fig.add_trace(go.Scatter(x=cutoffs, y=b_aki, mode="lines+markers",
                         name="傳統 MLE（隨截取門檻變動）",
                         line=dict(color=ACCENT)))
fig.add_hline(y=b_pos[0.2], line_color="#1baf7a", line_dash="dot",
              annotation_text=f"b-positive（修剪 0.2）= {b_pos[0.2]:.2f}")
fig.add_hline(y=b_pos[0.6], line_color="#1baf7a",
              annotation_text=f"b-positive（修剪 0.6）= {b_pos[0.6]:.2f}")
fig.add_hline(y=b_true, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text="真值 b = 1.0")
apply_layout(fig, title="合成不完整目錄：傳統估計 vs b-positive",
             xaxis_title="截取門檻（假設的 Mc）", yaxis_title="估計的 b 值",
             hovermode="x")
fig

# %% [markdown]
# 合成目錄的真實 b 值是 1.0，但完整度在規模 3.0 附近平滑劣化。
# 傳統估計（藍線）在門檻取低時嚴重低估，要把門檻抬到真實 Mc 之上
# 約 0.2–0.4 才回到真值——代價是丟掉大量資料。b-positive（綠線）
# 完全不需要知道 Mc 在哪：淺修剪（0.2）就把四成的低估縮到約一成，
# 修剪門檻加大到 0.6 後殘餘偏差幾乎消失。注意它並非一步到位——
# 這個合成情境是「觀測網長期偵測能力不足」，不是 b-positive 原本
# 針對的主震後短期不完整，殘餘偏差是已知的限制（Lippiello &
# Petrillo 2024）；而「加大修剪就收斂」正是近年系統性檢驗的結論。
#
# 順帶一提，同一份系統性檢驗（Tinti & Gasperini 2024）還發現
# 「非取**正**差不可」其實證據不足——正差、負差、絕對差表現相當，
# 關鍵在修剪。方法學還在演進中，這正是統計地震學活著的證據。
#
# ## 10.3 Omori 律再訪：c 值、p 值與時間窗
#
# 第 5 章的 Omori 圖用 $1/t$ 參考線帶過。完整的修正 Omori 律
# （Omori–Utsu）是
#
# $$n(t) = \frac{K}{(t + c)^{p}}$$
#
# 多出來的兩個參數各藏一個故事。
#
# **c 值**：讓 $t=0$ 不發散的「平緩期」，典型值幾分鐘到幾小時。
# 它是物理（主震後應力重分布需要時間）還是假象（早期餘震被波形
# 淹沒、記不到）？這場論戰至今未決——而 10.1 節的短期不完整
# 正好提供了「假象派」的彈藥：主震剛過時 Mc 可以飆到 4 以上，
# 隨 $\log t$ 才慢慢回落。看看 0403 花蓮主震後的頭幾天：

# %% tags=["hide-input"]
aft = cat[cat.time > main.time].copy()
t_days = (aft.time - main.time).dt.total_seconds() / 86400
fig = go.Figure(go.Scattergl(
    x=t_days, y=aft.ML, mode="markers",
    marker=dict(size=4, color=ACCENT, opacity=0.5), name="餘震"))
tt = np.logspace(-3, 1, 100)
fig.add_trace(go.Scatter(
    x=tt, y=7.2 - 4.5 - 0.76 * np.log10(tt), mode="lines",
    line=dict(color=QUAKE_COLOR, dash="dash"),
    name="Mc(t) 南加州型經驗式"))
apply_layout(fig, title="主震後的早期不完整：小餘震要過一陣子才「浮出來」",
             xaxis_title="主震後時間（天，對數軸）", yaxis_title="規模 ML",
             xaxis_type="log", yaxis_range=[2.8, 7.4], hovermode="closest")
fig

# %% [markdown]
# 主震後最初的一小時內，目錄裡幾乎只有 M4 以上的餘震——不是小的
# 沒發生，是記不到。紅色虛線是 Helmstetter 等人的南加州經驗式
# $M_c(t) = M_m - 4.5 - 0.76\log_{10}t$；台灣也有自己的版本
# （Tsai et al. 2012）。
# 在這條線以下做任何統計，都是在統計儀器而不是地震。
#
# **p 值**：教科書說 $p\approx 1$，但它既不是常數、也對分析選擇
# 敏感。台灣的證據有兩條：
#
# - Tsai et al. (2012) 用台灣目錄發現 $p$ 隨主震規模**線性增加**：
#   $p \simeq 0.11\,M_m + 0.38$——主震愈大、餘震衰減愈快，且用
#   三種除叢法做都成立；
# - 2022 池上序列（Chen et al. 2024）：同一份資料，取前 6 天擬合
#   $p = 1.39$，取 30 天 $p = 0.92$。**報告 p 值而不報告時間窗，
#   等於沒有報告。**
#
# ## 10.4 Båth 定律：最大餘震有多大？
#
# Båth (1965) 的經驗律說：最大餘震比主震約小 $\Delta M \approx 1.2$。
# 這條律在防災上太重要了——它回答「接下來最壞會多壞」。但它是
# 「定律」嗎？用長期目錄自己算一次。做法：以簡化的 Gardner–Knopoff
# 時空窗連結餘震，對每個 ML≥6.0 的主震找出窗內最大餘震
# ——這個「簡化」馬上會自己變成教材：

# %% tags=["hide-input"]
sub = cat_long[cat_long.ML >= 3.5].reset_index(drop=True)
t_num = sub.time.astype("int64").to_numpy() / 86400e9   # 天
lat, lon, ml = (sub[c].to_numpy() for c in ("latitude", "longitude", "ML"))

def gk_window(m):
    """Gardner & Knopoff (1974) 的除叢時空窗（公里、天）。"""
    L = 10 ** (0.1238 * m + 0.983)
    T = 10 ** (0.032 * m + 2.7389) if m >= 6.5 else 10 ** (0.5409 * m - 0.547)
    return L, T

is_aft = np.zeros(len(sub), dtype=bool)
bath = []
order = np.argsort(-ml)  # 由大到小處理
for i in order:
    if ml[i] < 6.0:
        break
    if is_aft[i]:
        continue                       # 自己是更大事件的餘震，跳過
    L, T = gk_window(ml[i])
    dist = np.hypot((lat - lat[i]) * 111,
                    (lon - lon[i]) * 111 * np.cos(np.radians(lat[i])))
    inwin = (dist <= L) & (t_num > t_num[i]) & (t_num <= t_num[i] + T)
    inwin &= ml < ml[i]
    is_aft |= inwin
    if inwin.sum() >= 5:               # 至少 5 個餘震才計入
        bath.append(ml[i] - ml[inwin].max())

fig = go.Figure(go.Histogram(x=bath, xbins=dict(size=0.2),
                             marker_color=ACCENT, opacity=0.85))
fig.add_vline(x=1.2, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text="Båth：ΔM = 1.2")
fig.add_vline(x=float(np.mean(bath)), line_color="#1baf7a",
              annotation_text=f"本目錄平均 {np.mean(bath):.2f}")
apply_layout(fig, title=f"主震與最大餘震的規模差（台灣 1973–2025，"
                        f"ML≥6 主震共 {len(bath)} 個）",
             xaxis_title="ΔM（主震 − 最大餘震）", yaxis_title="序列數",
             hovermode="x")
fig

# %% [markdown]
# 圖上先冒出一個驚訝：平均只有 0.3 左右，離教科書的 1.2 很遠！
# 先別急著宣布 Båth 錯了——換幾組視窗參數重算，平均值就在 0.3
# 到 0.5 之間跳動。問題出在「餘震」的定義：視窗法把窗內鄰近的
# 獨立大地震（台灣常見的雙主震與重疊序列）也算成了餘震，ΔM 被
# 系統性壓扁——序列怎麼定義，答案就跟著變。改用逐事件連結的
# 除叢法，台灣的系統性研究（Chan & Wu 2013，706 個序列）給出 $\bar{\Delta}_1 = 1.20 \pm 0.73$：**標準差 0.73 意味著
# 「平均 1.2」對單一序列幾乎沒有預測力**。更麻煩的是 $\Delta_1$
# 還隨餘震數量變：只取餘震數大於 50 的序列，平均掉到 0.74。
# 極端案例是集集地震——若用序列自身的 GR 律外推「推論最大餘震」
# $m^* = a/b$，得到的 $\Delta M^*$ 只有 0.03（Lee et al. 2013）：
# 集集的餘震多到 GR 律預期「最大餘震應該幾乎跟主震一樣大」
# （實測差 0.95）。
#
# 教訓與第 5 章的 b 值如出一轍：**平均值成立，不代表可以拿來預測
# 個案**。防災上真正有用的反而是條件式的陳述，例如台灣資料顯示
# 主震 M>6.5 時，最大餘震幾乎都落在 35 公里內（Wang et al. 2016）
# ——範圍有界，就能劃設警戒區。
#
# ## 10.5 除叢：把背景與叢集拆開
#
# 上面的 Båth 計算偷偷用了一個重要操作：**除叢（declustering）**
# ——把目錄拆成「背景地震」與「餘震（被觸發者）」。長期危害度
# 分析（第 16 章）傳統上只想要背景率，於是各種除叢法應運而生：
# Gardner–Knopoff 時空窗、Reasenberg 連結法、最近鄰距離法……
# 除叢後的目錄有個好性質：背景地震近似 Poisson 過程，累積數
# 隨時間近似直線：

# %% tags=["hide-input"]
kept = sub[~is_aft]
fig = go.Figure()
for data, name, color in [(sub, "原始目錄（ML≥3.5）", ACCENT),
                          (kept, "除叢後（近似背景活動）", "#1baf7a")]:
    fig.add_trace(go.Scattergl(x=data.time, y=np.arange(1, len(data) + 1),
                               mode="lines", name=name,
                               line=dict(color=color)))
for yr, label in [("1999-09-21", "集集"), ("2024-04-03", "花蓮")]:
    fig.add_vline(x=yr, line_dash="dot", line_color=QUAKE_COLOR)
apply_layout(fig, title="除叢前後的累積事件數（簡化視窗法，僅以 ML≥6 為主震；紅線＝集集、0403 花蓮）",
             xaxis_title="時間", yaxis_title="累積事件數", hovermode="x")
fig

# %% [markdown]
# 原始目錄（藍）在集集與 0403 花蓮處有明顯的階梯跳升——那是餘震
# 序列；除叢後（綠）曲線平滑許多，接近等速率的直線。看起來很成功？
#
# 但除叢藏著一個近年才被看清的陷阱。Mizrahi et al. (2021) 用五類
# 常用除叢法處理加州目錄，發現**除叢後的「主震」b 值比全目錄低了
# 最多 30%**——而且拿「b 值設計上完全均一」的合成目錄重跑，一樣
# 出現下降。原因是純粹的**選擇效應**：多數演算法把「叢集中最大的
# 事件」定義為主震，而小地震本來就不容易成為一群地震中最大的那個，
# 於是被除掉的比例偏高，規模分布的斜率就被壓平了。「主震的 b 值
# 天生比較低」這句流傳甚廣的話，至少有一大半是演算法造的。
#
# 這件事的後果直通危害度評估：被壓低的 b 值外插到高規模端，會在
# 工程上最關心的規模區間（加州的例子是 M6.9–8.8）**高估**大地震
# 的相對頻率；同時「只算主震」又**低估**了餘震的危害。兩個錯不會
# 抵成一個對。也因此，近年的趨勢是乾脆**繞開除叢**：與其在任意的
# 主震定義上做統計，不如用能同時描述背景與觸發的模型——這正是
# 下一章 ETAS 的出場理由。
#
# ## 10.6 目錄是一段行政史
#
# 把這一章倒過來讀，你會發現每一節都在說同一件事：**目錄的統計
# 性質，一半來自地球，另一半來自觀測系統與分析者的選擇。** Mc 隨
# 測站建設變、b 值隨規模尺度換代跳（台灣 1973–1987 年約 0.83、
# 1994 年後約 0.99——不是地體構造變了，是規模定義換了）、p 值隨
# 時間窗變、Båth 隨餘震數變、除叢法本身會製造統計假象。
#
# 所以在台灣用目錄做任何長期統計之前，先記住幾個該切一刀的年份：
# **1973**（儀器目錄起點）、**1991**（改用 $M_L$）、**1994**
# （連續記錄，多數現代研究的起點）、**2012**（微震偵測能力再升級）。
# 以及一條實用的換算式：台灣區域 $M_w = 0.87\,M_L + 0.23$
# （AutoBATS）——$M_L\,7.0$ 其實只是 $M_w\,6.3$，把 $M_L$ 目錄
# 直接餵給以 $M_w$ 校準的模型，b 值與觸發率都會系統性偏差。
#
# 好消息是：這些坑都有人踩過、也都有對策——這正是這一章存在的
# 意義。現在地基打穩了，可以開始蓋房子。{doc}`下一章 <11_etas>`
# 登場的是現代地震預報的主力引擎：ETAS。它把你在第 5 章與這一章
# 看到的兩條經驗律
# ——GR 與 Omori——縫成一條條件強度函數，讓「地震觸發地震」
# 變成一個可以擬合、可以模擬、可以檢驗的數學物件。
