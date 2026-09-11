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
# # 18. 臺灣地球物理觀測總覽
#
# 前一部從地震目錄建立了預報模型：過去的事件可以告訴我們，未來哪些地方的
# 發生率較高、餘震如何衰減，以及一份預報應該怎樣接受檢驗。現在把視線移向
# 產生資料的儀器。模型中的一個點，要經過測站記錄、訊號辨識、定位與規模
# 估計，才成為目錄的一列。
#
# 臺灣的觀測資料讓這個過程更具體：地震波、地下水、地磁與地表位移都可以
# 連續觀測，同一場地震在各類資料中卻有不同的形狀。地震儀可能記下幾十秒的劇烈振盪，
# 觀測井可能出現水位偏移，GNSS 則可能在震前後的座標之間留下階變。要把它們
# 放在一起，先要知道每種儀器對什麼物理量有反應。
#
# 這一部沿著觀測機制讀圖，再將讀到的現象帶回統計問題。下列測站清單是本書
# 儲存的 GDMS 資料快照；它描述教材資料來源，不是今日全部運作測站的即時盤點。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
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
# ## 18.1 測網各自在量什麼
#
# 一個測站的價值，不能只用每秒產生幾個數字衡量。取樣快有助於觀察短暫變化，
# 但可解析的訊號還受儀器頻帶、校正、所在環境和資料處理影響。先把本書用到的
# 資料產品分開：
#
# | 測網／產品 | 直接記錄的量 | 可以研究的現象 | 讀圖前要確認 |
# |---|---|---|---|
# | CWASN 波形 | 地震儀數位訊號 | 震波到時、頻率及傳播 | 儀器響應、頻道方向 |
# | TSMIP 強震紀錄 | 加速度感測器訊號 | 場址搖晃、工程需求 | 校正、動態範圍與場址 |
# | GW 觀測井 | 水位、氣壓、水溫 | 水文與應變相關反應 | 井構造、含水層及抽水 |
# | MAGNET 地磁 | 磁場分量 | 日變化、磁擾與局部變化 | 方向、參考站及外源場 |
# | GNSS 衛星觀測 | 偽距、載波相位等 | 解算後的速度與位移 | 解算方式及參考框架 |
#
# 這些量之間存在物理連結，但沒有一個換算比例能適用全臺所有測站。例如同樣
# 的地殼應變，兩口井可能因含水層性質不同而呈現不同水位反應；同樣的斷層滑移，
# 不同位置的 GNSS 站也會往不同方向移動。
#
# ## 18.2 把測站放回地圖
#
# 地圖右上角可以切換測網圖層。先看每一種測網覆蓋了哪些區域，再比較不同
# 測網能否在相近位置互相對照。顏色只是區分資料類型，點的密集程度不是地震
# 風險的直接量測。
#
# %% tags=["remove-input"]
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
# 地圖上的空白有兩種可能：那裡較少布站，或本書儲存的清單未包含可用資料。
# 兩者都不能解讀成當地沒有地震。陸上與外海的幾何差異尤其重要：震源若被測站
# 包圍，到時資料通常比較能限制位置；測站若都在同一側，某些方向的不確定性
# 就較難縮小。這正是前一部「目錄完整度與定位品質」在觀測端的來源。
#
# 強震觀測關心有人居住的場址會如何搖晃，地下水觀測受可用井及含水層限制，
# 地磁觀測還要避開局部電磁干擾。布站目的不同，因此不能要求不同測網在地圖上
# 完全重合。跨觀測分析必須保留這個差異，不能把相距很遠的站當成同一個位置。
#
# ## 18.3 站數不等於有效樣本數
#
# 下圖統計清單內各測網的站數。它適合盤點可查詢的站碼，卻沒有計入每站的
# 缺測天數、不同年代的運作狀況，也沒有扣掉鄰站共享的訊號。十個同步受到
# 同一場磁暴影響的站，不能當作十次獨立的地震證據。
#
# %% tags=["remove-input"]
counts = (st.groupby("network").size()
            .reindex(NETWORKS.keys()).rename("count").reset_index())
fig = px.bar(counts, x="count", y="network", orientation="h",
             color="network", color_discrete_map=NETWORK_COLORS,
             labels={"count": "測站數", "network": "測網"})
apply_layout(fig, title="GDMS 各測網測站數", showlegend=False, hovermode="y")
fig.update_traces(marker_line_width=0)
fig

# %% [markdown]
# ## 18.4 測網與資料一起改變
#
# 把清單中的日期累積起來，可以看到測網資料涵蓋範圍如何隨時間擴大。這裡使用
# `start_time` 欄位，應讀作清單所記錄的起始日期；若要判定儀器實際建置、
# 資料上線或換代的日期，還需要測站沿革檔案。
#
# 先留意曲線突然上升的年份，再想它會如何影響統計：更多站可能讓小地震
# 比較容易被記錄，也可能讓定位更穩定。年事件數上升，因此有地球活動與觀測
# 能力兩種解釋。只靠累積站數圖，還不能把某個轉折歸因於單一政策或災害。
#
# %% tags=["remove-input"]
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
# 一條長時間序列同時保留了地球與觀測系統的歷史。儀器更換、測站遷移、
# 資料格式調整可能造成階變；真正的同震位移也會造成階變。形狀相似，原因
# 卻不同，測站紀錄也是科學資料的一部分。
#
# 對統計模型而言，這表示觀測機制可能隨時間改變。模型若把新增測站帶來的
# 偵測改善，全解釋成地下活動增加，就會把資料取得的歷史誤寫成地震的歷史。
# 反過來，知道測網的變化，也讓我們能選擇較一致的時段、調整門檻，或明確
# 建模偵測機率。
#
# 接下來的{doc}`資料產品與品質 <02_download>`會從這張地圖往下一層走：
# 同一個測站可能提供原始訊號、整理後的序列或事件目錄，選哪一種產品，取決於
# 我們要回答的問題。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **入門・免費文章**：蕭乃祺（2019），〈[臺灣地震測報的發展](https://www.ntsec.edu.tw/liveSupply/detail.aspx?a=6829&cat=6841&lid=16154&p=1&print=1)〉，《科學研習》58 卷 6 期。
#   優先讀「基礎」中的觀測網介紹，對照本章地圖上的地震、GNSS 與地下水測站；留意歷史測站數不等於目前仍在運作的站數。
# - **資料來源・免費網站**：中央氣象署／FDSN，〈[Central Weather Administration Seismographic Network（T5）](https://www.fdsn.org/networks/detail/T5/)〉。
#   查閱地震測網的代碼、營運機構與資料引用方式，理解測網代碼與單一測站代碼的差別。這是地震測網的登錄資料，不代表本站全部地球物理測網。
# - **系統背景・免費官方文章**：中央氣象署地震測報中心，〈[資訊服務－地球物理資料管理系統](https://scweb.cwa.gov.tw/zh-tw/page/twenty/129)〉，《20 週年專刊》。
#   說明不同觀測系統的資料如何整合、儲存與供人查詢，可接著思考本章的測網沿革為何會影響資料涵蓋範圍。
