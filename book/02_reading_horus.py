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
# # 2. 讀一份地震目錄：HORUS 與義大利
#
# HORUS 目錄的一列長這樣：1960 年 1 月 4 日 09:20，北緯 43.13、東經 13.17，
# 深度 0 km，Mw 3.94。這一列告訴我們四件事：什麼時候、在哪裡、多深、多大。
# 整份目錄有 49 萬列。上一章看過一張預報圖；那張圖的每一個數字，
# 都是從這種列算出來的。
#
# 這一章只做一件事：把 HORUS 逐欄讀懂。讀完你能說出三個決定的理由。
# 為什麼要統一規模？為什麼只留 40 km 以內的地震？為什麼要畫兩個區域？
# 章末會立起一張「規格卡」，先填好資料側的欄位，其餘留給後面的章。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import HTML, display

from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

horus = italy.load_horus()                 # 完整目錄，不做任何篩選
cat = italy.experiment_catalog()           # 套用規格卡資料側篩選後的工作目錄
cells = italy.testing_cells()
poly = italy.collection_polygon()

n_horus = len(horus)
last_day = horus.time.max().strftime("%Y-%m-%d")
row = horus.iloc[1]                        # 開場那一列（第 0 列是一顆 290 km 的深震）

# %% [markdown]
# ## 2.1 一列資料有哪些欄位
#
# HORUS 的全名是 Homogenized instrumental seismic catalogue of Italy，
# 由 Lolli 等人（2020）編製，收錄 1960 年至今的義大利儀器地震。
# **儀器目錄**（instrumental catalogue）指每一筆都來自地震儀的記錄，
# 而不是歷史文獻。下表列出原始檔的欄位，範例值取自開場那一列。

# %% tags=["remove-input"]
columns = pd.DataFrame([
    ("Year … Se", "發震時刻（UTC），年月日時分秒六欄", "—", row.time.strftime("%Y-%m-%d %H:%M:%S")),
    ("Lat、Lon", "震央的緯度、經度", "度", f"{row.lat:.4f}、{row.lon:.4f}"),
    ("Depth", "震源深度", "km", f"{row.depth:.1f}"),
    ("Mw", "均一化後的矩規模", "—", f"{row.mw:.2f}"),
    ("sigMw", "Mw 的估計不確定度（一個標準差）", "—", f"{row.sig_mw:.2f}"),
    ("Geo-Ita", "震央是否在義大利國界內（* 為是）", "旗標", "*" if row.in_italy else ""),
    ("Geo-CPTI15", "震央是否在 CPTI15 多邊形內（* 為是）", "旗標", "*" if row.in_cpti15 else ""),
    ("Ev. type", "事件種類；x 表示爆炸等非地震事件", "代碼", "（空白＝地震）"),
    ("Iside n.", "對應 INGV 即時目錄 ISIDe 的事件編號", "編號", "（早期為空）"),
], columns=["欄位", "意義", "單位", "範例值"])
columns

# %% [markdown]
# 四個核心欄位是時間、經緯度、深度、規模。另外兩個旗標欄位替我們先判斷過
# 「震央在不在某個區域內」。Geo-CPTI15 判斷用的多邊形，就是這個實驗的
# **收集區**（collection region, S）：提供輸入事件的範圍，2.6 節細講。
# 這一章後面會自己重算一次旗標，看看差在哪裡。
#
# 目錄裡也有不是地震的事件。本章的每一個決定都是一次篩選；
# 篩完留下的那份資料叫**工作目錄**（working catalogue）。
# 第一個篩選就是用 Ev. type 欄把爆炸等事件排除。

# %% tags=["remove-input"]
n_non_eq = int((~horus.is_earthquake).sum())
n_in_S = int(horus.in_cpti15.sum())
print(f"HORUS 2024 年 6 月版：{n_horus:,} 筆，最後一筆 {last_day}")
print(f"其中非地震事件 {n_non_eq:,} 筆；Geo-CPTI15 旗標為 * 的 {n_in_S:,} 筆")

# %% [markdown]
# ## 2.2 規模：為什麼要統一，為什麼要取整
#
# 開場那一列的規模寫 Mw 3.94，不確定度 0.23。**矩規模**
# （moment magnitude, Mw）由地震釋放的地震矩換算而來，是目前最通用的規模尺度。
# 但 1960 年代的儀器測不出地震矩；當年報的是別種規模，例如近震規模 $M_L$。
# 不同規模在同一顆地震上會差零點幾個單位，直接混用會把資料弄髒。
#
# **均一化**（homogenization）就是把各種來源的規模，用經驗關係換算成同一種 Mw。
# HORUS 的作者做了這件事，並附上每一筆的不確定度 sigMw。
# 全目錄 sigMw 的中位數是 0.2，換句話說，規模的第二位小數沒有意義。
#
# 因此 Biondini 等人（2023）把所有規模取整到 0.1，公式是
#
# $$m_{b} = \frac{\lfloor 10\,m_{\rm raw} + 0.5 \rfloor}{10}.$$
#
# **取整規模**（binned magnitude, $m_b$）是本站往後所有計算用的規模。
# 取整有一個常被忽略的後果：門檻的位置會移動。
# 原始 Mw 4.95 取整後變成 5.0，所以「取整後 $\ge 5.0$」等於「原始 $\ge 4.95$」。
# 同理，名目門檻 2.5 對應原始 2.45。第 3 章談門檻時會直接沿用這兩個數字。

# %% tags=["remove-input"]
edge_cases = cat[(cat.mw >= 4.95) & (cat.mw < 5.0)]
print(f"工作目錄中原始 Mw 落在 4.95–4.99 的地震：{len(edge_cases)} 筆，"
      f"取整後都算 Mw 5.0；其中 {int(edge_cases.in_R.sum())} 筆在 177 格內（測試區 R，2.6 節）")

# %% [markdown]
# ## 2.3 從年計數尋找完整度的線索
#
# 下圖上半是 1960–2021 年收集區內 $M_w \ge 4.0$ 的規模–時間圖。
# 下半是三個最小規模下的每年事件數。看下半圖：$M_w \ge 4.0$ 的年計數
# 六十年來大致持平；$M_w \ge 2.5$ 的年計數卻從 1960 年代起漲了好幾倍。

# %% tags=["remove-input"]
big = cat[cat.mb >= 4.0]
years = np.arange(1960, 2022)
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.45, 0.55],
                    vertical_spacing=0.06,
                    subplot_titles=("收集區 S 內 Mw ≥ 4.0 的地震", "每年事件數（對數尺度）"))
fig.add_trace(go.Scattergl(x=big.time, y=big.mb, mode="markers", name="Mw ≥ 4.0",
                           marker=dict(size=4, color=ACCENT, opacity=0.45)), row=1, col=1)
decade_note = {}
for i, m in enumerate([2.5, 3.0, 4.0]):
    yc = cat[cat.mb >= m].groupby("year").size().reindex(years, fill_value=0)
    decade_note[m] = (yc.loc[1960:1969].mean(), yc.loc[2010:2019].mean())
    fig.add_trace(go.Scatter(x=pd.to_datetime(years, format="%Y"), y=yc.values, mode="lines",
                             name=f"Mw ≥ {m:.1f}", line=dict(color=PALETTE[i], width=2)),
                  row=2, col=1)
y_top = cat[cat.mb >= 2.5].groupby("year").size().max()
for yr, m, k in [(1981, 3.0, 1.5), (1990, 2.5, 3.0)]:
    fig.add_vline(x=pd.Timestamp(f"{yr}-01-01"), line_dash="dot",
                  line_color="rgba(90,100,120,0.6)", row=2, col=1)
    fig.add_annotation(x=pd.Timestamp(f"{yr}-01-01"), y=np.log10(y_top * k), xanchor="left",
                       text=f"{yr} 起 Mw ≥ {m} 完整", showarrow=False, font=dict(size=11),
                       row=2, col=1)
fig.update_yaxes(title_text="取整規模 Mw", range=[3.9, 6.9], row=1, col=1)
fig.update_yaxes(title_text="每年事件數", type="log", row=2, col=1)
r25, r40 = decade_note[2.5], decade_note[4.0]
apply_layout(fig, height=640, hovermode="x",
             title=(f"Mw ≥ 4.0 年均 {r40[0]:.0f} → {r40[1]:.0f} 筆；"
                    f"Mw ≥ 2.5 年均 {r25[0]:.0f} → {r25[1]:.0f} 筆（1960 年代 → 2010 年代）"))
fig

# %% [markdown]
# 小震計數的增加與測網改善的歷史相符，但這張圖無法單獨分開
# 觀測能力與真實活動的變化。目錄沒有記到的地震叫**漏報**（missed events）。
# 從某個規模起，目錄幾乎記到了全部的地震，這個規模叫**完整度規模**
# （magnitude of completeness, $M_c$）。
#
# Lolli 等人（2020）給 HORUS 的完整度是：$M_w \ge 4.0$ 自 1960 年起、
# $\ge 3.0$ 自 1981 年、$\ge 2.5$ 自 1990 年、$\ge 2.1$ 自 2003 年、
# $\ge 1.8$ 自 2005 年。
# 這些年代來自文獻中的完整度研究，不是由本圖的曲線平坦程度推得。
# 完整目錄仍可有很強的真實率變動。1990 年開始學習，是為了對齊
# 名目 2.5 的輸入門檻；更早的暖身資料仍需注意小震不完整。
# 第 8 章會自己估 $M_c$，並確認模型用的輸入門檻不低於它。

# %% [markdown]
# ## 2.4 深度：為什麼只留 40 km 以內
#
# HORUS 的第一列其實是一顆深震：1960 年 1 月 3 日，深度 290 km，Mw 6.34，
# 位置在第勒尼安海下方。義大利有兩種地震。一種是陸上淺層斷層的地震，
# 另一種在南方海域下方數百公里深。下圖是收集區內 0–100 km 的深度直方圖。

# %% tags=["remove-input"]
in_S = horus[horus.is_earthquake & horus.in_cpti15 & (horus.time < "2022-01-01")]
n_deep = int((in_S.depth > 40).sum())
n_zero = int((in_S.depth == 0).sum())
n_ten = int((in_S.depth == 10).sum())
n_neg = int((in_S.depth < 0).sum())
big_all = in_S[in_S.mb >= 5.0]
n_big_deep = int((big_all.depth > 40).sum())
edges = np.arange(0, 101, 1.0)
counts, _ = np.histogram(in_S.depth.clip(lower=0), bins=edges)
fig = go.Figure()
fig.add_trace(go.Bar(x=edges[:-1] + 0.5, y=counts, width=0.9, marker_color=ACCENT,
                     name="每 1 km 的事件數"))
fig.add_vrect(x0=40, x1=100, fillcolor="gray", opacity=0.12, line_width=0,
              annotation_text=f"深度 > 40 km：{n_deep:,} 筆不進工作目錄",
              annotation_position="top right")
fig.add_annotation(x=10.5, y=np.log10(n_ten), text=f"10.0 km：{n_ten:,} 筆", showarrow=True,
                   ax=60, ay=-10)
fig.add_annotation(x=0.5, y=np.log10(n_zero), text=f"0.0 km：{n_zero:,} 筆", showarrow=True,
                   ax=60, ay=-30)
fig.update_yaxes(type="log", title_text="事件數（對數尺度）")
fig.update_xaxes(title_text="深度（km）")
apply_layout(fig, height=420, hovermode="x",
             title=f"收集區 S 內 1960–2021 年地震的深度分布（{len(in_S):,} 筆）")
fig

# %% [markdown]
# 兩個尖峰值得停下來看。深度剛好 10.0 km 的有幾萬筆，剛好 0.0 km 的也有一萬多筆。
# 這通常是定位程序給的**固定深度**（fixed depth）。資料不夠決定深度時，
# 程式就填一個預設值。這種欄位值是假的，只能當「淺層」讀，不能拿來算深度分布。
# 另有一筆深度是 −0.2 km，代表震源在海平面以上，那是火山區的地震。
#
# Biondini 等人（2023）只保留深度 $\le 40$ km 的地震。這是一次**截斷**
# （truncation）：超出範圍的資料整筆不用。被截掉的地震多在國界外的海域下方，
# 比例見下方輸出。模型要學的是陸上淺層地震。

# %% tags=["remove-input"]
print(f"S 內 Mw ≥ 5.0 的地震共 {len(big_all)} 顆，其中 {n_big_deep} 顆深度超過 40 km，被截斷")
deep = in_S[in_S.depth > 40]
print(f"深度 > 40 km 的 {len(deep):,} 筆中，Geo-Ita 旗標為國界外的占 {(~deep.in_italy).mean():.0%}")
print(f"深度為負值的地震：{n_neg} 筆")

# %% [markdown]
# ## 2.5 位置：經緯度與公里
#
# 目錄用經緯度記位置，模型卻用公里算距離。緯度一度大約 111 km，
# 但經度一度的長度隨緯度變短：在北緯 42 度只有約 83 km。
# 一個在經緯度上看起來是正方形的格子，在地面上其實是扁的長方形。
#
# 所以實驗把經緯度換成**投影座標**（projected coordinates）。
# 這是一套把地球表面攤平、以公里為單位的平面座標。本站用義大利官方的
# EPSG:7794 投影，工具箱的 `lonlat_to_km()` 負責換算。
# 測試區 R（2.6 節）的 177 個方格是在這套公里座標上切的。每格邊長 $30\sqrt{2}\approx 42.4$ km，
# 面積 1,800 km²。換回經緯度畫在地圖上時，格子會略呈梯形。

# %% [markdown]
# ## 2.6 兩個區域：收集區 S 與測試區 R
#
# 這個實驗畫了兩個區域，一大一小，小的完全在大的裡面。
# **測試區**（testing region, R）是被評分的範圍。預報只對這 177 格出數字，
# 目標地震也只算這 177 格內的。**收集區**（collection region, S）是提供
# 輸入事件的範圍。模型計算時，S 內所有符合條件的地震都算歷史。

# %% tags=["remove-input"]
show_diagram("d02_regions_s_r",
             caption="收集區 S 包住測試區 R。R 邊上的格會受到 R 外地震的影響，所以輸入事件要從更大的 S 收。")

# %% [markdown]
# 為什麼要兩個區域？因為地震會互相影響，而影響不認格線。
# R 邊上那一格的活動，有一部分來自 R 外幾十公里的地震。
# 如果輸入事件也只收 R 內的，邊上的格就會被低估，這叫**邊界效應**
# （edge effect）。解法是把輸入範圍放大到 S。
#
# S 的邊界取自 CPTI15。**CPTI15**（Catalogo Parametrico dei Terremoti Italiani）
# 是義大利的參數化歷史地震目錄，由 Rovida 等人（2020）編製，
# 以歷史文獻與儀器資料重建每顆地震的參數。Biondini 等人用它看 1600–1959 年，
# 用 HORUS 看 1960 年以後。
# CPTI15 的涵蓋多邊形有 14 個頂點，HORUS 的 Geo-CPTI15 旗標就是對它判斷的。

# %% tags=["remove-input"]
ring_lon, ring_lat = [], []
for _, c in cells.iterrows():
    ring_lon += list(c.corner_lon) + [None]
    ring_lat += list(c.corner_lat) + [None]
fig = go.Figure()
fig.add_trace(go.Scattergeo(lon=big.lon, lat=big.lat, mode="markers", name="Mw ≥ 4.0（1960–2021）",
                            marker=dict(size=3, color=ACCENT, opacity=0.35), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(lon=ring_lon, lat=ring_lat, mode="lines", name="測試區 R：177 格",
                            line=dict(color="rgba(36,48,64,0.55)", width=0.8), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(lon=poly.lon, lat=poly.lat, mode="lines", name="收集區 S：CPTI15 多邊形",
                            line=dict(color=PALETTE[1], width=2.2), hoverinfo="skip"))
fig.update_geos(fitbounds="locations", visible=False, showcountries=True,
                countrycolor="rgba(120,130,150,0.5)", showcoastlines=True,
                coastlinecolor="rgba(120,130,150,0.5)", showland=True,
                landcolor="rgba(238,242,247,0.6)", projection_type="mercator")
n_S_cat, n_R_cat = len(cat), int(cat.in_R.sum())
apply_layout(fig, height=560, hovermode="closest", margin=dict(l=10, r=10, t=50, b=10),
             legend=dict(orientation="h", y=-0.02),
             title=f"工作目錄：S 內 {n_S_cat:,} 筆，其中 R 內 {n_R_cat:,} 筆")
fig

# %% [markdown]
# 圖上藍點是 1960–2021 年 $M_w \ge 4.0$ 的震央。S 內、R 外的藍點大多在海上；
# 它們不會被評分，卻會當作輸入事件進入模型。
#
# 177 格是怎麼選的？Biondini 等人（2023）沿用 Gasperini 等人（2021）的規則。
# 一格要留下，格內必須有至少一顆陸上的 $M \ge 4.0$ 地震。
# 依據是 CPTI15 的 1600–1959 年紀錄，加上 HORUS 的 1960–2021 年紀錄。
# 離島等不相連的格再剔除。請注意「1960–2021 年」這段：
# 選格用到了 2021 年前的資訊，包含測試期本身。
# 所以本站的實驗是**既有研究區的重演**。實驗設計並沒有全部在測試期前凍結。
# 這件事不影響資料側的讀法，但第 3 章談擬前瞻時要記得。
#
# 最後一個技術細節：多邊形頂點的順序。工具箱讀進來的 CPTI15 頂點是順時針排列，
# 最後一點重複第一點。順時針可以從有號面積為負看出來（下方輸出）。
# EEPAS 的原始程式要求順時針；Jalilian（2019）的 ETAS R 套件則要求反時針，
# 且首點不可重複。同一個多邊形餵給兩套軟體，要先各自整理一次頂點。

# %% tags=["remove-input"]
x, y = poly.x_km.to_numpy(), poly.y_km.to_numpy()
signed_area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
corner_lon = np.concatenate(cells.corner_lon.to_list())
corner_lat = np.concatenate(cells.corner_lat.to_list())
n_corner_in = int(italy.in_collection_region(corner_lon, corner_lat).sum())
print(f"CPTI15 多邊形：{len(poly) - 1} 個頂點，有號面積 {signed_area / 1e3:,.0f} 千 km²"
      f"（負值＝順時針）；面積約 {abs(signed_area) / 1e6:.2f} 百萬 km²")
print(f"177 格的 {len(corner_lon)} 個角點全部落在 S 內：{n_corner_in == len(corner_lon)}")

# %% [markdown]
# ## 2.7 旗標可以信嗎？自己重算一次
#
# Geo-CPTI15 是 HORUS 的作者算好的旗標。我們有多邊形的頂點，可以自己判斷
# 每個震央在不在裡面，再和旗標對照。下表是對照結果。

# %% tags=["remove-input"]
own = italy.in_collection_region(horus.lon, horus.lat)
table = pd.crosstab(horus.in_cpti15.map({True: "旗標：在 S 內", False: "旗標：在 S 外"}),
                    pd.Series(own, name="").map({True: "自算：在 S 內", False: "自算：在 S 外"}))
table.index.name = ""
table

# %% tags=["remove-input"]
diff = horus[horus.in_cpti15 & ~own].copy()
dx, dy = italy.lonlat_to_km(diff.lon, diff.lat)
px, py = poly.x_km.to_numpy(), poly.y_km.to_numpy()
dist = []
for xi, yi in zip(dx, dy):
    best = np.inf
    for i in range(len(px) - 1):
        ax, ay, bx, by = px[i], py[i], px[i + 1], py[i + 1]
        t = np.clip(((xi - ax) * (bx - ax) + (yi - ay) * (by - ay)) / ((bx - ax) ** 2 + (by - ay) ** 2), 0, 1)
        best = min(best, np.hypot(xi - (ax + t * (bx - ax)), yi - (ay + t * (by - ay))))
    dist.append(best)
diff["dist_km"] = dist
diff["cell"] = italy.cell_of(dx, dy, cells)
n_diff, max_dist = len(diff), max(dist)
n_diff_in_R = int((diff.cell >= 0).sum())
fig = go.Figure()
fig.add_trace(go.Scattergeo(lon=poly.lon, lat=poly.lat, mode="lines", name="收集區 S",
                            line=dict(color=PALETTE[1], width=2.2), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(lon=ring_lon, lat=ring_lat, mode="lines", name="測試區 R",
                            line=dict(color="rgba(36,48,64,0.45)", width=0.6), hoverinfo="skip"))
fig.add_trace(go.Scattergeo(
    lon=diff.lon, lat=diff.lat, mode="markers", name="旗標在內、自算在外",
    marker=dict(size=8, color=PALETTE[3], line=dict(color="white", width=0.8)),
    text=[f"{t:%Y-%m-%d} Mw {m:.1f}，距邊界 {d:.1f} km" for t, m, d in zip(diff.time, diff.mw, diff.dist_km)],
    hovertemplate="%{text}<extra></extra>"))
fig.update_geos(fitbounds="locations", visible=False, showcountries=True,
                countrycolor="rgba(120,130,150,0.5)", showcoastlines=True,
                coastlinecolor="rgba(120,130,150,0.5)", showland=True,
                landcolor="rgba(238,242,247,0.6)", projection_type="mercator")
apply_layout(fig, height=520, hovermode="closest", margin=dict(l=10, r=10, t=50, b=10),
             legend=dict(orientation="h", y=-0.02),
             title=f"兩種判斷不一致的 {n_diff} 筆：全部距 S 邊界 {max_dist:.1f} km 以內，"
                   f"{n_diff_in_R} 筆在 R 內")
fig

# %% [markdown]
# 49 萬筆裡只有 23 筆不一致，而且全貼著多邊形的邊線（最遠距離見圖說）。
# 差異來自兩邊的多邊形版本或判斷方式略有不同。沿著邊線，這種事很正常。
# 這 23 筆沒有一筆在 R 內，對目標地震毫無影響。對輸入事件的影響也可以忽略。
#
# 這個對照給了一個讀目錄的習慣。別人算好的旗標可以用，但要知道它怎麼算。
# 最好再用自己的方法重算一次。本站的工作目錄沿用 HORUS 的旗標。

# %% [markdown]
# ## 2.8 本章填入的規格欄位
#
# 一份預報實驗要先說清楚十一件事，本站把它們寫成一張規格卡。
# 這一章填好了資料側的五個欄位；其餘欄位標出會在哪一章填。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("目錄與規模尺度", f"HORUS 2024 年 6 月版（Lolli et al. 2020），Mw 均一化，規模取整到 0.1；全目錄 {n_horus:,} 筆"),
    ("時間範圍", "1960–2021 年（工作目錄只取到 2021 年底）"),
    ("深度", f"≤ 40 km（S 內截掉 {n_deep:,} 筆深震）；排除非地震事件"),
    ("收集區 S", f"CPTI15 多邊形，{len(poly) - 1} 個頂點；工作目錄 {n_S_cat:,} 筆"),
    ("測試區 R", f"177 個 30√2 km 方格（EPSG:7794 公里座標，共 {cells.area_km2.sum():,.0f} km²）；"
                f"R 內 {n_R_cat:,} 筆。選格規則用到 2021 年前的資訊，屬既有研究區的重演"),
    ("三段期間", "待填（第 3 章）"),
    ("發報時間軸", "待填（第 3 章）"),
    ("門檻", "待填（第 3 章給角色、第 8 章給 M_c）"),
    ("規模箱", "待填（第 3 章）"),
    ("輸出格式", "待填（第 3 章）"),
    ("評分方式", "待填（第 16–18 章）"),
])))

# %% [markdown]
# 資料側說清楚了，題目本身還沒有。學習哪幾年？測試哪幾年？每幾個月出一張預報？
# 目標是多大的地震？下一章把這些一次說清楚：{doc}`標準預報實驗規格 <03_experiment_spec>`。
#
# ## 參考資料與延伸閱讀
#
# - Lolli, B., Randazzo, D., Vannucci, G. 與 Gasperini, P.（2020），[The Homogenized Instrumental Seismic Catalog (HORUS) of Italy from 1960 to Present](https://doi.org/10.1785/0220200148)，*Seismological Research Letters* 91(6), 3208–3222。HORUS 的原始論文。先讀規模均一化的做法與完整度表，再看[目錄網站](http://horus.bo.ingv.it/)的欄位說明。目錄可免費下載，論文全文需訂閱或館藏權限。
# - Rovida, A., Locati, M., Camassi, R., Lolli, B. 與 Gasperini, P.（2020），[The Italian earthquake catalogue CPTI15](https://doi.org/10.1007/s10518-020-00818-y)，*Bulletin of Earthquake Engineering* 18, 2953–2984。出版社全文可能需訂閱；[波隆那大學典藏](https://cris.unibo.it/retrieve/0eca64cc-4a09-474e-8e16-5d56cbbe643c/BEEE-D-19-00711_R1.pdf)有作者稿。歷史目錄如何從文獻重建參數。讀第 2 節就能理解本章的收集區多邊形從何而來。
# - Biondini, E., Rhoades, D. A. 與 Gasperini, P.（2023），[Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123)，*Geophysical Journal International* 234, 1681–1700。本站主線的來源論文；本章對應「Application to Italy」一節，可對照圖 1、圖 2 與規模取整公式。
# - Naylor, M. 等（2023），[Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759)，*Frontiers in Applied Mathematics and Statistics*，開放取用。第 3 節用合成資料示範漏報如何扭曲估計，適合在讀第 8 章之前先建立直覺。
