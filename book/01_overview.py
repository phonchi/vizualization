# %% [markdown]
# # 1. 台灣地球物理觀測總覽
#
# 中央氣象署在全台佈設了多個觀測網，GDMS 彙整了其中七個測網的資料。
# 本章先載入**真實的測站清單**（由 GDMS 取得，內建於本專案
# `data/stations/`），看看它們的規模與分布。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %%
import folium
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from gdms_toolkit import load_stations, NETWORKS
from gdms_toolkit.viz import NETWORK_COLORS, apply_layout

st = load_stations()          # 全部測網
st.groupby("network").size().rename("測站數").to_frame().join(
    pd.Series(NETWORKS, name="測網名稱")
)

# %% [markdown]
# 各測網觀測的物理量與取樣方式整理如下：
#
# | 測網 | 觀測量 | 儀器 | 取樣 |
# |---|---|---|---|
# | CWASN | 地震波（速度） | 短週期／寬頻地震儀 | 100 Hz 連續 |
# | TSMIP | 強地動（加速度） | 強震儀 | 觸發／連線站連續 |
# | GW | 地下水位、氣壓、水溫 | 觀測井（水位計） | 1 秒 |
# | MAGNET | 地磁三分量 X, Y, Z | 磁力儀 | 1 秒 |
# | GNSS 系列 | 衛星定位原始觀測 | GNSS 接收儀 | 30 秒（RINEX） |
#
# ## 測站分布互動地圖
#
# 每個測網一個圖層（右上角可切換），顏色全書一致：
# 藍＝地下水、紫＝地磁、橘＝地震、紅＝強震、青／綠／洋紅＝GNSS。

# %%
m = folium.Map(location=[23.7, 121.0], zoom_start=7, tiles="cartodbpositron")

for net, name in NETWORKS.items():
    fg = folium.FeatureGroup(name=f"{net}（{name}）",
                             show=net not in ("TSMIP",))  # TSMIP 533 站預設先關
    sub = st[st.network == net].dropna(subset=["lat", "lon"])
    for _, r in sub.iterrows():
        folium.CircleMarker(
            location=[r.lat, r.lon], radius=4,
            color=NETWORK_COLORS[net], fill=True, fill_opacity=0.8, weight=1,
            tooltip=f"[{net}] {r.station_code} {r.chinese_station_name}",
        ).add_to(fg)
    fg.add_to(m)

folium.LayerControl(collapsed=False).add_to(m)
m

# %% [markdown]
# 幾個值得注意的分布特徵：
#
# - 強震站（TSMIP，紅）幾乎鋪滿全島，都會區特別密。它的任務是記錄強烈
#   搖晃，供工程設計與地震速報使用，所以哪裡有人、哪裡就要有站。
# - 地下水觀測井（GW，藍）只有少數幾口有連續即時資料釋出，集中在平原
#   含水層。井的數量少，但每口的取樣密（1 秒），資訊量其實很大。
# - 地磁站（MAGNET，紫）刻意遠離都市，因為車輛與電力設施的電磁干擾會
#   淹沒地磁訊號。東部縱谷沿線那串測站，底下正是歐亞與菲律賓海板塊的
#   縫合帶。
# - GNSS（青／綠）沿著活動斷層與花東縱谷加密，那裡是台灣變形最快的地方。
#
# ## 各測網測站數量

# %%
counts = (st.groupby("network").size()
            .reindex(NETWORKS.keys()).rename("count").reset_index())
fig = px.bar(counts, x="count", y="network", orientation="h",
             color="network", color_discrete_map=NETWORK_COLORS,
             labels={"count": "測站數", "network": "測網"})
apply_layout(fig, title="GDMS 各測網測站數", showlegend=False, hovermode="y")
fig.update_traces(marker_line_width=0)
fig

# %% [markdown]
# ## 測網是怎麼長出來的？
#
# 測站清單裡的 `start_time`（設站日期）記錄了每個測網的擴建史。
# 把各測網的累積測站數畫成時間曲線：

# %%
st["start_year"] = pd.to_datetime(st.start_time, errors="coerce").dt.year
fig = go.Figure()
for net in NETWORKS:
    years = st.loc[st.network == net, "start_year"].dropna().sort_values()
    if years.empty:
        continue
    fig.add_trace(go.Scatter(x=years, y=list(range(1, len(years) + 1)),
                             mode="lines", name=net,
                             line=dict(color=NETWORK_COLORS[net], width=2)))
apply_layout(fig, title="各測網累積測站數（依 GDMS 登錄的設站日期）",
             xaxis_title="年", yaxis_title="累積測站數")
fig

# %% [markdown]
# 台灣觀測網的擴建，幾乎每一波都跟在一場災害地震後面：
#
# - 1999 集集地震（M7.6）之後，強震觀測與地震速報系統大幅擴充，也催生
#   了跨部會的地震前兆觀測計畫。地下水井與地磁站的現代化佈設多從 2000
#   年代中期開始。
# - 2016 美濃、2018 花蓮地震之後，GNSS 與強震網在斷層帶周邊繼續加密。
# - 曲線上的「階梯」對應一次次計畫性佈建。要注意 GDMS 登錄的是資料開始
#   釋出的日期，實際建站通常更早。
#
# 一條看似平凡的連續時間序列，背後是十幾年的儀器維運與經費投入。這也
# 提醒我們一件做資料分析時很容易忽略的事：序列裡的變化，不見得是地球
# 的變化。儀器更換、測站遷移、韌體升級，都會在資料裡留下看起來像訊號
# 的假階變。分析前先弄清楚測站的沿革，往往比急著找異常更重要。
