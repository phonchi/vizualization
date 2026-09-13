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
# # 8. 從哪個規模開始相信目錄：完整度 Mc
#
# 名目規模 2.5 是前章的估計下限。
# 若小地震漏記，前章的平均規模會偏高。
# 本章先看直方圖，再決定哪些規模可用。
# 完整度規模（magnitude of completeness）記為 $M_c$。
# 它是指定時空內可近似收錄齊全的下限。
#
# HORUS 跨越數十年的觀測環境。
# 同一個下限，需要接受年代與地點的檢查。
# 本章詳講一個容易重現的方法。
# 另外兩法用來比較證據與取捨。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

# %% [markdown]
# ## 8.1 先看最高的規模箱
#
# 下圖每個柱子只收一個規模箱。
# 非累積次數（non-cumulative count）就是箱內筆數。
# 完整樣本常在小規模端有較多事件。
# 若偵測開始漏掉小震，柱高可能轉而下降。
#
# 最高的柱子提供一個候選位置。
# 最大曲率法（maximum curvature，MaxC）用這個峰估下限。
# 名稱來自累積曲線的轉折。
# 實作時，先找非累積直方圖的峰。
# 本頁再加上 0.2 的經驗修正。
#
# $$\widehat M_c^{\rm nom}=m_{\rm peak}+0.2.$$
#
# 上標表示名目規模，也就是箱心。
# 修正量（correction term）是額外提高下限的幅度。
# 0.2 用來降低峰值偏低的影響。
# 這是經驗選擇，效果仍依目錄而變。
# 它也需要與其他完整度證據對照。

# %% tags=["remove-input"]
cat = italy.experiment_catalog()
learning_all = cat.loc[cat.in_R & (cat.period == "learning")]
# 零值峰單獨報告；MaxC 的操作樣本限正規模。
learning = learning_all.loc[learning_all.mb > 0]
dm = 0.1
def maxc(values, min_n=100):
    values = np.asarray(values)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < min_n:
        return np.nan
    bins, counts = np.unique(np.round(values, 1), return_counts=True)
    return round(float(bins[np.argmax(counts)]) + 0.2, 1)
mc = maxc(learning.mb)
values, counts = np.unique(learning.mb.to_numpy(), return_counts=True)
fig = go.Figure(go.Bar(x=values, y=counts, marker_color=ACCENT, name="正規模樣本"))
fig.add_vline(x=mc, line_dash="dot", line_color="#555555")
apply_layout(fig, title="學習期規模直方圖：從峰值找候選下限",
    xaxis_title="取整規模 Mw", yaxis_title="各箱事件數", xaxis_range=[0, 5])
display(Markdown(f"圖 1｜正規模樣本 {len(learning):,} 筆。"
    f"另有 {(learning_all.mb == 0).sum():,} 筆零規模，未進入本圖估計。"
    f"MaxC 加修正的名目下限為 {mc:.1f}。"))
fig

# %% [markdown]
# 直方圖外的零值，需要另外追查其意義。
# 本頁先把零值獨立列出，限定分析正規模。
# 這項操作使估計結果有清楚的適用樣本。
# 它也表示，本頁無法判斷零規模附近的完整性。
#
# 一個柱子的峰，也可能來自資料匯整方式。
# 選擇效應（selection effect）是收錄規則改變樣本組成。
# 若原始資料已有規模下限，峰可能貼著下限。
# 把資料先切到 2.5，再估完整度，便有此風險。
# 所以，本頁先讀取完整的實驗目錄。
#
# ## 8.2 一個總體數字會遮住年代差異
#
# 下圖每個點代表一年的估計。
# 移動窗（moving window）是沿時間平移的資料區間。
# 圖中的另一條線使用三年窗。
# 每個窗都截止於橫軸所標的年份。
# 它能減少單一年份的抽樣起伏。

# %% tags=["remove-input"]
years = np.arange(1960, 2012)
r_all = cat.loc[cat.in_R & (cat.year <= 2011) & (cat.mb > 0)]
yearly = np.array([maxc(r_all.loc[r_all.year == y, "mb"]) for y in years])
rolling = np.array([maxc(r_all.loc[(r_all.year >= y-2) & (r_all.year <= y), "mb"])
    if y >= 1962 else np.nan for y in years])
fig = go.Figure(go.Scatter(x=years, y=yearly, mode="markers", name="逐年 MaxC＋修正",
    marker=dict(color=ACCENT, size=5)))
fig.add_trace(go.Scatter(x=years, y=rolling, mode="lines", name="向前三年窗",
    line=dict(color="#4a3aa7")))
ref_years = [1960, 1981, 1990, 2003, 2005, 2011]
ref_mc = [4.0, 3.0, 2.5, 2.1, 1.8, 1.8]
fig.add_trace(go.Scatter(x=ref_years, y=ref_mc, mode="lines", name="文獻名目下限",
    line=dict(color="#555555", dash="dash", shape="hv")))
apply_layout(fig, title="不同年代，需要不同的規模下限",
    xaxis_title="窗的最後一年", yaxis_title="名目完整度 Mw")
display(Markdown(f"圖 2｜每窗至少 {100} 筆正規模才估計。"
    f"缺點表示資料不足。實線使用 {3} 年窗。"
    f"虛線依 Lolli（2020），由 Biondini（2023）引述。"))
fig

# %% [markdown]
# 虛線給出文獻對義大利本土的判斷。
# 1960 年起使用名目 4.0。
# 1981 年起降到 3.0，1990 年起降到 2.5。
# 2003 年起為 2.1，2005 年起為 1.8。
# 這些數字來自目錄研究，並非本頁重新擬合。
#
# 實線與虛線可能落在不同高度。
# 總體 MaxC 會偏向資料量大的年代。
# 後期小震多，整段直方圖的峰就偏向後期。
# 因此，整段估出的低門檻，不能代表早期表現。
# 讀者應先選期間，再解讀完整度。
#
# 三年窗的相鄰點共用大部分事件。
# 線形變平順，表示統計摘要較平滑。
# 點與點的差距仍含相依性。
# 不能把每個點當成獨立的年代證據。
# 發報前的分析，也應只使用截止前的窗。
#
# ## 8.3 同一年代，位置也有差別
#
# 下圖沿用實驗的方格。
# 每格只使用學習期內的正規模事件。
# 顏色表示各格的 MaxC 候選下限。
# 灰色格表示樣本不足，暫不給估計。
# 門檻高低與資料量需要一起閱讀。

# %% tags=["remove-input"]
cells = italy.testing_cells()
cell_mc = np.array([maxc(learning.loc[learning.cell == j, "mb"]) for j in cells.cell])
cell_n = np.array([(learning.cell == j).sum() for j in cells.cell])
features = []
for row in cells.itertuples():
    ring = [[float(x), float(y)] for x, y in zip(row.corner_lon, row.corner_lat)]
    features.append(dict(type="Feature", id=int(row.cell),
        geometry=dict(type="Polygon", coordinates=[ring])))
geojson = dict(type="FeatureCollection", features=features)
fig = go.Figure(go.Choropleth(geojson=geojson, locations=cells.cell,
    z=np.zeros(len(cells)), colorscale=[[0,"#dddddd"],[1,"#dddddd"]],
    showscale=False, marker_line_width=0.5, hoverinfo="skip"))
ok = np.isfinite(cell_mc)
fig.add_trace(go.Choropleth(geojson=geojson, locations=cells.cell[ok], z=cell_mc[ok],
    colorscale=SEQUENTIAL, colorbar_title="名目 Mc", marker_line_width=0.5,
    customdata=cell_n[ok], hovertemplate="格 %{location}<br>Mc=%{z:.1f}<br>N=%{customdata}<extra></extra>"))
fig.update_geos(fitbounds="locations", visible=False)
apply_layout(fig, title="學習期各格的完整度候選值", height=520)
display(Markdown(f"圖 3｜共 {len(cells)} 格，{ok.sum()} 格可估計。"
    f"其中 {(cell_mc[ok]>2.5).sum()} 格高於名目 {2.5:.1f}。"
    f"這是資料診斷圖，顏色未表示預報率。"))
fig

# %% [markdown]
# 灰色格仍屬於測試區。
# 資料不足（insufficient sample）表示方法缺少估計依據。
# 把它填成零，會誤導讀者以為觀測更完整。
# 本頁因此保留缺值，並列出可估計格數。
#
# 深色格提示較高的候選下限。
# 位置差異可能涉及測站分布與事件組成。
# MaxC 圖本身無法拆開這些原因。
# 若某格高於輸入門檻，應列為資料風險。
# 整體門檻通過比較，也不能替每格背書。
#
# ## 8.4 大震後，完整度可能暫時變高
#
# 下圖從 2009 年 L'Aquila 主震起算。
# 主震（mainshock）是序列中作為參照的大事件。
# 地震波重疊時，小事件較難辨認。
# 短期不完整（short-term incompleteness）指這段暫時漏記。
# 圖中每點是一筆實際記錄的規模。

# %% tags=["remove-input"]
candidate = cat.loc[(cat.year == 2009) & (cat.mb >= 5)]
main = candidate.loc[candidate.mb.idxmax()]
radius = 50.0
distance = np.hypot(cat.x_km-main.x_km, cat.y_km-main.y_km)
tau = (cat.time-main.time).dt.total_seconds()/86400
seq = cat.loc[(distance <= radius) & (tau > 0) & (tau <= 30) & (cat.mb > 0)].copy()
seq["tau"] = tau.loc[seq.index]
fig = go.Figure(go.Scatter(x=seq.tau, y=seq.mb, mode="markers",
    marker=dict(color=ACCENT, size=4, opacity=0.45), name="局部目錄事件"))
fig.add_hline(y=2.5, line_dash="dot", line_color="#555555")
apply_layout(fig, title="L’Aquila 主震後：小規模端隨時間怎樣變？",
    xaxis_title="主震後天數（對數刻度）", xaxis_type="log", yaxis_title="取整規模 Mw")
display(Markdown(f"圖 4｜主震規模 {main.mb:.1f}。"
    f"取震央周圍 {radius:.0f} km、其後 {30} 天。"
    f"圖中有 {len(seq):,} 筆正規模記錄。虛線為名目 {2.5:.1f}。"))
fig

# %% [markdown]
# 圖的左端如果缺少小事件，便值得追查。
# 但稀疏的點也可能只是較少事件的結果。
# 要估早期下限，可比較更短窗的規模分布。
# 每個短窗都需要足夠樣本。
# 波形重新偵測的資料能再提供外部證據。
#
# 固定的規模線只是分析選擇。
# 第 10 章會提高門檻，並改變擬合起點。
# 這種敏感度分析（sensitivity analysis）觀察選擇改變後的結果。
# 若參數大幅變動，資料範圍就需要重新檢視。
# 單一曲線的吻合程度不足以決定原因。
#
# ## 8.5 三種方法，各看哪項證據？
#
# 下表把峰值、斜率與分布分開。
# $b$ 穩定法（b-value stability）找提高下限後斜率穩定的位置。
# KS 距離（Kolmogorov–Smirnov distance）量測兩條累積分布的最大差距。
# 兩法都以門檻以上近似 GR 為出發點。
#
# | 方法 | 觀察的量 | 選下限的方式 | 主要取捨 |
# |---|---|---|---|
# | MaxC 加修正 | 非累積峰值 | 峰值再提高 | 快速；需檢查峰的來源 |
# | b 穩定法 | 多個下限的 b | 變動落入容許誤差 | 需較多事件支撐 |
# | KS 比較 | 經驗與 GR 分布 | 模擬後比較距離 | 需校準估參與分箱效應 |
#
# KS 模擬要沿用觀測的規模精度。
# 每次模擬也應重新估計分布參數。
# 這樣才把估參造成的貼合納入比較。
# 第 6 章的「未拒絕」仍只代表相容。
# 它提供統計證據，完整性仍要結合目錄研究。
#
# 一個小樣本可能讓三法都缺少辨別力。
# 一個大樣本則能顯出很細微的分布偏離。
# 因此，方法比較應連同樣本數報告。
# 本頁只把 MaxC 數字填入診斷欄。
# 文獻門檻則保留在實驗設定欄。
#
# ## 8.6 把三種下限放在同一把尺上
#
# 名目 2.5 對應箱左界 2.45。
# 輸入門檻（input threshold）$m_0$ 使用這個左界。
# 目標門檻（target threshold）$m_T$ 決定預報的地震集合。
# 完整度則描述觀測資料的收錄能力。
# 比較時，三者都須先換成相同表示法。

# %% tags=["remove-input"]
show_diagram("d03_thresholds", caption=f"三門檻示意｜輸入名目規模 {2.5:.1f}；目標名目規模 {italy.SPEC.mT:.1f}。")

# %% [markdown]
# 1990 年起的文獻下限是名目 2.5。
# 換成箱界後，$m_0=2.45\ge M_c=2.45$。
# 這項確認適用於採用的整體文獻設定。
# 局部高下限與早期序列仍保留診斷標記。
# 本頁的 MaxC 低值，也不會自動降低輸入門檻。
#
# 1960–1989 年的暖身期，使用較高完整度下限。
# 暖身期（warm-up period）提供正式學習前的事件歷史。
# EEPAS 的長時間核會用到較早的小震。
# 時間核（temporal kernel）描述事件影響隨時間的分布。
# 漏掉早期小震，可能使後續貢獻不足。
#
# 完整度的補償公式需要假設遺漏的分布。
# 第 15 章會說明補償在模型中的位置。
# 歷史開始前的缺失，也有另外的時間問題。
# 讀者應先記住哪些年代較難收錄小震。
# 資料範圍的紀錄，是解讀補償的起點。
#
# 一份完整度紀錄應同時保存四項資訊。
# 第一項是區域，第二項是資料期間。
# 第三項是估計方法，第四項是規模精度。
# 這四項決定別人能否重做同一個判斷。
# 單寫一個下限，會失去它的適用範圍。
#
# 本頁的格圖也提示後續分析的順序。
# 先回查高下限格的直方圖與事件年份。
# 再看峰值是否來自少量集中事件。
# 若需要合併鄰格，須重新註明空間範圍。
# 合併能增加樣本，也會平均掉地點差異。
# 這項取捨應隨分析目的決定。
#
# 規格卡保留文獻值與診斷值兩欄。
# 後續作者便能知道預報使用哪個設定。
# 讀者也能沿資料圖追查這個設定的弱點。
# 完整度估計因此成為可重查的資料紀錄。
#
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("學習期 Mc：文獻設定", "名目 2.5；箱左界 2.45，1990 年起"),
    ("MaxC：整段診斷", f"正規模樣本；名目 {mc:.1f}，非逐格完整性保證"),
    ("輸入門檻核對", "m0＝2.45 ≥ 文獻 Mc 箱界 2.45；局部風險另列"),
    ("暖身期", "1960–1989；早期小震不完整"),
    ("PPE 空間分配", "待填（第 9 章）"),
])))

# %% [markdown]
# 名目下限已經有資料依據。
# 接著把歷史事件的位置轉成預報。
# 下一章是 {doc}`鄰近過去地震：PPE 預報 <09_ppe_forecast>`。
#
# ## 參考資料與延伸閱讀
#
# - SeismoStats，[完整度三法](https://seismostats.readthedocs.io/latest/user/estimate_mc.html)。免費文件。先比較各法使用的統計量。
# - Lolli 等（2020），The homogenized instrumental seismic catalog (HORUS) of Italy from 1960 to present。刊於 SRL，91，3208–3222。本文是各年代完整度的資料來源。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://www.earth-prints.org/handle/2122/17084)。免費機構全文。第 2 節交代 HORUS 資料。表 3 列出固定參數。
