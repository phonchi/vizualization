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
# # 17. 作業化預報：讓資料、模型與使用者接得起來
#
# {doc}`上一章 <21_psha>`把事件發生模型轉成場址地動，但一次離線計算還不是持續可用的預報。真實地震發生後，位置和規模可能先有初值，幾分鐘或更久之後再修訂；小事件可能尚未補齊，使用者卻已經需要判斷接下來的情境。
#
# 作業化預報把前面各章接起來：資料進來，按照事先定義的規則更新模型，產生指定時間窗的產品，儲存發布版本，最後等觀測累積後檢驗。模型的統計表現、資料品質與數字是否能被正確理解，都會影響這條流程。
#
# ## 17.1 從發布時刻回頭看
#
# 想像在今天上午發布未來一天的預報。先確定資料截止時刻、可用目錄與品質標記；接著估計或更新模型，對未來歷史積分或模擬，再轉成事件數、機率或地動產品。
#
# 發布後保留版本十分重要。後來完整目錄的地震數，可能不是當時系統看得到的地震數。若回溯重算時只用最新目錄，得到的是另一個問題的答案，不能假裝完全重現當時預報能力。
#
# ## 17.2 用 STEP 理解短期更新
#
# Gerstenberger 等人（2005）的 STEP 是一個歷史案例：它結合背景地震度與餘震叢集模型，估計未來短時間窗的事件及搖晃。它讓我們看見前面 GR、Omori 衰減及地動計算如何在同一產品中連起來。
#
# 資料少時，需要借用區域經驗；序列逐漸累積後，才有條件估計更細的序列或空間特徵。以下示意不同複雜度的反應，重點是資料增加如何支持更具體的模型，而不是最複雜版本在所有時刻都最好。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.integrate import quad
from scipy.stats import gamma as gamma_dist

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

rng = np.random.default_rng(22)

N_MIN_FIT = 100                              # 序列自身參數的資料門檻
N_SEQ = np.logspace(0, 3.4, 400)             # 累積的完整規模以上餘震數


def _switch(n, n_half, width=0.22):
    """在 log10(n) 上的平滑切換函數（示意用，非 AIC 權重公式）。"""
    return 1.0 / (1.0 + np.exp(-(np.log10(n) - np.log10(n_half)) / width))


score = np.vstack([
    np.ones_like(N_SEQ),                            # 通用參數：永遠可用
    3.2 * _switch(N_SEQ, 1.6 * N_MIN_FIT),          # 序列自身參數
    4.0 * _switch(N_SEQ, 7.0 * N_MIN_FIT),          # 參數容許空間變化
])
w = score / score.sum(axis=0)
n_cross = N_SEQ[np.argmax(w[1] > w[0])]             # 序列層超越通用層之處

fig = go.Figure()
for k, lab in enumerate(["通用參數（區域先驗）", "序列自身參數",
                         "參數容許空間變化"]):
    fig.add_trace(go.Scatter(x=N_SEQ, y=w[k], mode="lines", name=lab,
                             line=dict(width=0.6, color=PALETTE[k]),
                             stackgroup="one"))
fig.add_vline(x=N_MIN_FIT, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"{N_MIN_FIT} 個事件門檻")
apply_layout(fig,
             title=f"STEP 三元素的權重隨序列累積而移轉（示意；"
                   f"序列自身參數約在 {n_cross:.0f} 個事件處超越通用層）",
             xaxis_title="序列累積的完整規模以上餘震數（log）",
             xaxis_type="log", yaxis_title="權重（三者和為 1）",
             hovermode="x", height=420)
fig

# %% [markdown]
# 多個叢集可能在時空中重疊。若每個分量都包含相同背景，直接相加會重複計入。STEP 的特定組合規則處理這個問題；它與 ETAS 把各事件直接後代強度相加的生成結構不同，不能只把運運算元號互換。
#
# ```{admonition} 為什麼有些模型取最大值
# :class: dropdown
#
# 對兩個非負分量 $a,b$，有 $\max(a,b)\le a+b$，差值為 $\min(a,b)$。這說明重疊貢獻如何改變總量，但沒有證明取最大值就是正確機率模型。是否重複計算，取決於兩分量各自代表什麼；模型設計與觀測檢驗需一起考慮。
# ```
#
# ## 17.3 各國經驗提供不同角度
#
# Mizrahi 等人（2024）的回顧整理義大利、紐西蘭與美國的作業經驗。這些案例使用不同資料與制度，適合比較它們如何處理共同問題，而不是排成哪一國最好。
#
# 義大利案例凸顯持續測試與模型組合；紐西蘭案例強調公開產品和長序列中的溝通；美國的公開餘震預報示範機率、事件數範圍及模型資訊如何一起呈現。STEP 是歷史模型名稱，不能拿來概括後來所有 USGS OAF 產品。
#
# 下面用合成序列示範其中一個共同環節：從較寬的區域先驗開始，隨事件累積更新生產力倍數的後驗分布。曲線逐漸集中只表示這個假設模型下參數資訊增加，不保證模型誤差也同時消失。各國最新產品請從章末官方入口檢視，分清文獻案例與現行設定。
#
# %% tags=["remove-input"]
A_GEN, C_OM, P_OM = 30.0, 0.05, 1.08     # 通用參數下的餘震生產力（合成）
THETA_TRUE = 0.60                        # 本序列真實的生產力倍數
ALPHA0, BETA0 = 2.0, 2.0                 # 通用先驗 Gamma：均值 1、很寬

EDGES = np.array([0.0, 0.25, 0.5, 1.0, 3.0, 7.0, 30.0])
E_j = np.array([quad(lambda t: A_GEN * (t + C_OM) ** -P_OM, a, b)[0]
                for a, b in zip(EDGES[:-1], EDGES[1:])])
n_j = rng.poisson(THETA_TRUE * E_j)                    # 逐段觀測到的餘震數
alpha_t = ALPHA0 + np.concatenate([[0], np.cumsum(n_j)])
beta_t = BETA0 + np.concatenate([[0], np.cumsum(E_j)])

theta = np.linspace(0.02, 3.0, 500)
snap = [(0, "t = 0（通用先驗）"), (1, "0.25 天後"), (3, "1 天後"),
        (5, "7 天後"), (6, "30 天後")]
fig = go.Figure()
widths = []
for (k, lab), color in zip(snap, PALETTE):
    dist = gamma_dist(alpha_t[k], scale=1.0 / beta_t[k])
    widths.append(dist.ppf(0.95) - dist.ppf(0.05))
    fig.add_trace(go.Scatter(x=theta, y=dist.pdf(theta), mode="lines",
                             name=lab, line=dict(color=color, width=2.2)))
fig.add_vline(x=THETA_TRUE, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"本序列真值 {THETA_TRUE:.2f}")
apply_layout(fig,
             title=f"貝氏後驗逐日變窄：90% 區間寬度從 {widths[0]:.2f} "
                   f"收斂到 {widths[-1]:.2f}（合成序列）",
             xaxis_title="生產力倍數 θ（1 = 區域通用值）",
             yaxis_title="後驗密度", hovermode="x", height=430)
fig

# %% [markdown]
# ## 17.4 資料不完整時，預報會怎麼改變
#
# 主震之後波形重疊，小事件可能被漏掉；位置和規模修訂又會改變觸發貢獻。若把不完整目錄當成完整資料，系統可能低估活動，或將修訂造成的變化誤讀成地震過程改變。
#
# 因此作業系統除了模型參數，也需保留資料品質、處理延遲及更新時間。這和前面的完整度與參數不確定性直接相連：目錄不是一份永遠固定、完全正確的事件清單。
#
# ## 17.5 把相同數字說清楚
#
# 發布一個機率時，至少要讓讀者知道事件定義、時間窗、地理範圍和比較背景。只寫「風險增加十倍」會漏掉原來是千分之一還是十分之一；只寫「很低」也無法讓使用者理解某項行動是否值得。
#
# 下圖把小機率換成許多相同機率情況下的事件數量，以點陣表達。這是頻率框架的示意，不是在一張圖內測量不同說法的心理效應。機率溝通研究的意義在於，表達方式本身可以被測試，不能只靠發布者覺得文字已經清楚。
#
# %% tags=["remove-input"]
# Kaikōura M7.8 後兩週（2016-11-28 起算）GeoNet 的未來 7 天預報
KAIKOURA = [("M5.0–5.9", 98, 5.6), ("M6.0–6.9", 41, 0.53),
            ("M≥7.0", 5, 0.05)]
N_DOTS = 100                             # 點陣的點數
PER_DOT = 100_000 // N_DOTS              # 每個點代表幾個「地方」
idx = np.arange(N_DOTS)
gx, gy = idx % 10, -(idx // 10)

fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.06,
                    subplot_titles=[f"{lab}：{p}%（平均 {mu} 個）"
                                    for lab, p, mu in KAIKOURA])
for col, (lab, p, mu) in enumerate(KAIKOURA, start=1):
    hit = idx < p * N_DOTS // 100
    for mask, color, name in [(~hit, "#d5dde6", "不會發生"),
                              (hit, ACCENT, "會發生")]:
        fig.add_trace(go.Scatter(x=gx[mask], y=gy[mask], mode="markers",
                                 name=name, showlegend=(col == 1),
                                 marker=dict(size=11, color=color)),
                      row=1, col=col)
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False, scaleanchor="x")
apply_layout(fig,
             title=f"頻率框架：{N_DOTS * PER_DOT:,} 個「機率與此地相同」的"
                   f"地方裡，一週內會發生的有幾個"
                   f"（每點 = {PER_DOT:,} 個地方）",
             hovermode="closest", height=380)
fig

# %% [markdown]
# 當產品同時提供事件數範圍與至少一次機率，它們回答的是不同問題。範圍可以呈現序列活動有多大波動，機率則對應指定規模與時間內是否發生。兩者都應標示適用條件，不能讓讀者把小事件數量誤當成大地震機率。
#
# 有時情境敘述可以協助理解：較常見情況是活動逐漸下降，較少見情況可能包含相似規模事件，更少見但需準備的情況可能是更大事件。下面依 2024 年回顧所整理的 Kaikōura 歷史情境作教學重繪，以當時列出的三組機率說明產品形式。它不是現在對該區域發布的預報；衝擊軸只是敘事順序，並非量測損失。
#
# %% tags=["remove-input"]
SCEN = [("情境一：餘震照預報衰減", 70, 1,
         "無新增重大災害；序列逐步平息"),
        ("情境二：M7.0–7.8", 25, 2,
         "區域性破壞；隱沒帶可能參與、局部海嘯"),
        ("情境三：大於主震", 5, 3,
         "廣域破壞；含 M8 以上板塊介面破裂的可能")]

fig = go.Figure()
for (lab, p, imp, desc), color in zip(SCEN, PALETTE):
    fig.add_trace(go.Scatter(x=[p], y=[imp], mode="markers+text",
                             name=f"{lab}（{p}%）", text=[f" {p}%"],
                             textposition="middle right",
                             marker=dict(size=18 + 34 * imp / 3, color=color,
                                         opacity=0.85),
                             hovertext=[desc], hoverinfo="text"))
fig.update_yaxes(tickvals=[1, 2, 3],
                 ticktext=["序列平息", "區域性破壞", "廣域破壞"],
                 title_text="潛在衝擊（序數，示意）", range=[0.4, 3.6])
apply_layout(fig,
             title=f"Kaikōura 歷史情境（教學重繪）：機率合計 "
                   f"{sum(s[1] for s in SCEN)}%",
             xaxis_title="未來一年的機率（%）", xaxis_range=[-4, 88],
             hovermode="closest", height=430, showlegend=True)
fig

# %% [markdown]
# ## 17.6 降低中的機率，仍可能高於背景
#
# 一場序列的活動會隨時間衰減，但「比昨天低」與「已回到長期背景」不同。下面使用文獻整理的 Canterbury 歷史機率點，配合合成衰減曲線作示意；曲線經過這些點是設定結果，不是獨立擬合或當前預報。它也顯示新事件出現後，原先下降的機率可能再次升高。
#
# 發布端若只宣佈下降，讀者可能誤以為警戒理由已經消失；若只說仍高於背景，又可能忽略確實下降的趨勢。連續產品需要保持相同的時間窗、門檻與語言，讓不同發布版本可以比較。
#
# %% tags=["remove-input"]
T_DARFIELD, T_CHCH, T_M57 = 0.0, 0.47, 5.45      # 事件時間（年）
T_BEFORE = 5.075                                  # 2015-10-01 起算的一年
MU_BG, C_YR, P_YR = 0.05, 0.02, 1.10              # 背景率與 Omori（年）
P_BEFORE, P_AFTER = 0.50, 0.63                    # 公布的兩個機率


def _omori_int(t1, t2, t0):
    """單一母事件對 [t1, t2] 期望數的形狀因子。"""
    lo = max(t1, t0)
    return 0.0 if t2 <= lo else quad(
        lambda s: (s - t0 + C_YR) ** -P_YR, lo, t2)[0]


# 由 P(一年內至少一次) 反解兩個生產力：先解主序列，再解 M5.7
shape_b = (_omori_int(T_BEFORE, T_BEFORE + 1, T_DARFIELD)
           + 0.6 * _omori_int(T_BEFORE, T_BEFORE + 1, T_CHCH))
K_B = (-np.log(1 - P_BEFORE) - MU_BG) / shape_b
lam_base2 = (MU_BG + K_B * _omori_int(T_M57, T_M57 + 1, T_DARFIELD)
             + 0.6 * K_B * _omori_int(T_M57, T_M57 + 1, T_CHCH))
K_7 = (-np.log(1 - P_AFTER) - lam_base2) / _omori_int(T_M57, T_M57 + 1, T_M57)


def p_year(t):
    """在時刻 t 發布的「未來一年至少一次 M5.0–5.9」機率。"""
    lam = (MU_BG + K_B * _omori_int(t, t + 1, T_DARFIELD)
           + 0.6 * K_B * _omori_int(t, t + 1, T_CHCH))
    if t >= T_M57:
        lam += K_7 * _omori_int(t, t + 1, T_M57)
    return 1 - np.exp(-lam)


t_a = np.linspace(0.02, T_M57 - 1e-6, 160)
t_b = np.linspace(T_M57, 6.4, 100)
fig = go.Figure()
for tt, name, color in [(t_a, "M5.7 之前", ACCENT),
                        (t_b, "M5.7 之後", QUAKE_COLOR)]:
    fig.add_trace(go.Scatter(x=tt, y=[100 * p_year(s) for s in tt],
                             mode="lines", name=name,
                             line=dict(color=color, width=2.5)))
for t_ev, lab in [(T_DARFIELD, "M7.1"), (T_CHCH, "M6.2"), (T_M57, "M5.7")]:
    fig.add_vline(x=t_ev, line_dash="dot", line_color="#888888",
                  annotation_text=lab)
fig.add_trace(go.Scatter(x=[T_BEFORE, T_M57],
                         y=[100 * P_BEFORE, 100 * P_AFTER],
                         mode="markers", name="公布值（校準點）",
                         marker=dict(size=11, symbol="diamond",
                                     color="#4a3aa7")))
apply_layout(fig,
             title=f"序列第六年：一個 M5.7 把未來一年機率從 "
                   f"{100 * p_year(T_BEFORE):.0f}% 推到 "
                   f"{100 * p_year(T_M57):.0f}%（合成，校準至公布值）",
             xaxis_title="主震後時間（年）",
             yaxis_title="未來一年至少一次 M5.0–5.9 的機率（%）",
             hovermode="x", height=430)
fig

# %% [markdown]
# ## 17.7 預報要與使用問題一起設計
#
# 一份預報可以降低意外感，卻未必讓所有人採取相同行動。結構檢查、搜救部署、交通管理與家庭準備的成本和後果不同，正如前面成本損失模型所示。使用者需要的時間窗、地理範圍與呈現方式也不同。
#
# 與使用者共同設計產品，不表示可以改變統計結論來配合期待；它表示先弄清楚決策需要什麼資訊，再選擇適當摘要。機率與不確定性應可以追溯到同一套模型和資料，不在不同受眾之間變成互相矛盾的訊息。
#
# 若某種新現象尚無可靠統計模型，例如資料稀少的特殊序列，可以進行有結構的專家徵詢。但專家判斷仍需標示依據與分歧，不等同於經前瞻驗證的機率模型。新的觀測進來後，也應有可檢查的更新方式。
#
# ## 17.8 從預報回到觀測
#
# 到這裡，第一部的各章已經連成一個迴圈：先認識目錄能記錄什麼，再用統計描述活動，建立模型，檢查預報，最後在持續更新中面對資料和溝通限制。任何一環出現新問題，都可能需要回到觀測來源。
#
# 地震目錄之外，地下水、地磁、地動與 GNSS 各自量到不同物理量。它們可能幫助理解過程、約束模型，也可能受到天氣、儀器或環境變化影響。不能因為某個訊號出現在地震前，就跳過前面學過的替代解釋和獨立檢驗。
#
# {doc}`第二部：臺灣地球物理觀測 <01_overview>`將依觀測機制認識這些資料。接著仍要問：看到的是什麼量、在什麼尺度變化、有哪些其他原因、又能為地震研究增加什麼資訊？
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS，官方餘震預報說明（免費）。先看預報表、時間窗、機率與事件數範圍的讀法，再閱讀模型參數頁籤說明，對照本章從統計模型到發布產品的流程。
# - [Earthquake forecasts](https://www.geonet.org.nz/earthquake/forecast/) — GeoNet，官方預報說明（免費）。用較少公式說明機率預報的用途與限制，適合比較本章紐西蘭案例如何把數字連到不同使用者的需求。
# - [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等（2024），*Reviews of Geophysics*（免費全文）。本章各國系統與機率溝通的主要回顧來源，優先讀義大利、紐西蘭與美國的制度比較及使用者共同設計；論文記錄的是出版時的狀態，後續產品可配合前兩項官網閱讀。
# - [Real-time forecasts of tomorrow's earthquakes in California](https://doi.org/10.1038/nature03622) — Matthew C. Gerstenberger、Stefan Wiemer、Lucile M. Jones、Paul A. Reasenberg（2005），*Nature*（全文可能需訂閱）。本章 STEP 系統的原始論文，說明如何結合背景地震率與叢集模型，產生未來 24 小時的強震動機率圖；可與上面的 USGS OAF 說明比較不同時期的預報產品。
