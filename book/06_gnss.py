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
# # 23. GNSS：從衛星訊號到地表位移
#
# 地震儀擅長記錄震波造成的振盪；若要知道一場地震之後測站停在什麼新位置，
# 還需要地表位移資料。GNSS 利用衛星訊號估計接收站座標，經過精密解算與
# 適當的參考框架，能讓我們比較長期速度、地震前後的位移及較緩慢的變形。
#
# 這些變化的形狀對應不同問題。震間速度告訴我們地表如何持續變形；同震
# 階變約束地震造成的位移；震後曲線可能含有餘滑、黏彈性鬆弛及水文負載等
# 成分。單一曲線通常不能唯一識別地下機制，需要多站空間型態與物理模型。
#
# 臺灣位於板塊聚合區，但板塊相對運動速率、跨某條基線的縮短率，以及某一
# 測站的速度，並不是同一個數字。沒有指定參考框架與地理範圍，就不能把
# 一個「每年幾公分」的速率套用到全島。
#
# ## 23.1 衛星原始觀測還不是位移
#
# 本書使用的 GDMS GNSS 產品是 RINEX 觀測檔，儲存偽距與載波相位等資料。
# 它們受到衛星軌道、時鐘、大氣、天線與多路徑影響，需經解算才能得到座標。
# 本章不將原始觀測偽裝成位移成果；先看產品與品質，再說明一份解算序列
# 應怎樣讀。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input", "remove-output"]
import gzip
import tarfile
from collections import Counter

import pandas as pd
import plotly.graph_objects as go

from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, apply_layout

tgz = CACHE_DIR / "edu-gnss-hualien2024.tgz"
with tarfile.open(tgz) as tar:
    names = sorted(m.name for m in tar.getmembers() if m.name.endswith(".gz"))
print("\n".join(names))

# %% [markdown]
# 一份 RINEX 檔的標頭記錄版本、測站、觀測量種類及取樣間隔。它相當於
# 儀器資料的身分說明。標頭中也可能有近似座標，供解算初始化或位置查核，
# 但那不是每個時刻重新估計的精密位置。
#
# 檔名與欄位解析、ECEF 座標換算移到{doc}`附錄 G <appendix_g_observations>`。
# 主文留下與科學解讀直接相關的一件事：檔案裡有很多數字，不表示地表位移
# 已經被測成很多個時間點。
#
# %% tags=["remove-input", "remove-output"]
with tarfile.open(tgz) as tar:
    member = next(m for m in tar.getmembers() if "hual0930" in m.name)
    raw = gzip.decompress(tar.extractfile(member).read()).decode("ascii")
lines = raw.splitlines()
print("\n".join(lines[:20]))

# %% [markdown]
# ## 23.2 品質從可用觀測開始
#
# 下面顯示 HUAL 站某日各曆元記錄的衛星數。本書資料為三十秒間隔；衛星數
# 變化可能反映星座移動、遮蔽或觀測可用性。突然減少值得檢查，但它不是
# 地表突然移動的證據。
#
# 定位能力也受衛星幾何影響。即使衛星數一樣，分散在天空各方向與集中在
# 同一方向，對座標的約束也不同。這與上一章地震定位的測站幾何很相似：
# 觀測數量重要，觀測從哪個方向來同樣重要。
#
# %% tags=["remove-input"]
def parse_epoch_sats(lines):
    """從 RINEX 2 觀測檔抓出（時間, 衛星數）序列。"""
    out = []
    for ln in lines:
        # 曆元行：' 24  4  2  0  0  0.0000000  0 17G05G13...'
        if len(ln) > 32 and ln[28] == "0" and ln[:3].strip().isdigit():
            try:
                yy, mo, dd, hh, mi = (int(ln[1:3]), int(ln[4:6]), int(ln[7:9]),
                                      int(ln[10:12]), int(ln[13:15]))
                nsat = int(ln[29:32])
                out.append((pd.Timestamp(2000 + yy, mo, dd, hh, mi), nsat))
            except ValueError:
                continue
    return pd.DataFrame(out, columns=["time", "nsat"]).set_index("time")

sats = parse_epoch_sats(lines)
fig = go.Figure(go.Scattergl(x=sats.index, y=sats.nsat, mode="lines",
                             line=dict(color=ACCENT, width=1.2)))
apply_layout(fig, title=f"HUAL 站 2024/04/02 可見衛星數（30 秒取樣，共 {len(sats)} 曆元）",
             yaxis_title="衛星數", showlegend=False)
fig

# %% [markdown]
# 這張圖回答的是觀測供應的變化。要評估座標品質，還需要檢視解算殘差、
# 相位周跳、幾何和正式誤差等資訊。不能把衛星較多的時段直接當作所有
# 分量都同樣準確，也不能由衛星數反推出毫米級精度。
#
# ## 23.3 位移需要一個參考
#
# 地心直角座標可以換成經緯度與橢球高，也可以將相對位移轉為東、北、上
# 分量。幾何換算改變的是表示方式，不會提升原始座標的精度。因此把檔頭
# 近似座標轉到小數點後很多位，仍然不會得到同震位移。
#
# %% tags=["remove-input", "remove-output"]
import math
import re

def ecef_to_lla(x, y, z):
    """地心直角座標（公尺）轉為經緯度（度）與橢球高（公尺）。"""
    a, e2 = 6378137.0, 6.69437999014e-3     # WGS84
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - e2))
    h = 0.0
    for _ in range(6):
        N = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return math.degrees(lat), math.degrees(lon), h

m = re.search(r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+APPROX POSITION", raw)
xyz = [float(v) for v in m.groups()]
lat, lon, h = ecef_to_lla(*xyz)
print(f"ECEF: {xyz}")
print(f"→ 緯度 {lat:.5f}, 經度 {lon:.5f}, 橢球高 {h:.1f} m")

from gdms_toolkit import load_stations
hual = load_stations("GNSS").query("station_code == 'HUAL'").iloc[0]
print(f"測站清單 HUAL：緯度 {hual.lat}, 經度 {hual.lon}")

# %% [markdown]
# ## 23.4 讀一份已解算座標序列
#
# Nevada Geodetic Laboratory 等研究機構提供解算好的位置與速度產品。
# 閱讀前先確認參考框架、單位、取樣間隔與產品版本；同一站換到不同框架，
# 長期趨勢可能不同。將多站放在同一框架下，才能比較它們的相對運動。
#
# 閱讀東、北、上三分量圖時，可以先看長期斜率，再看週期起伏與突變。
# 斜率是測站在該框架中的速度；季節成分可能包含水文負載和其他環境效應；
# 階變則需對照地震與裝置更換紀錄。許多站同時改變，是值得追查的空間證據，
# 但參考框架或解算流程的共同誤差也可能跨站出現。
#
# | 序列形狀 | 可能關心的過程 | 還需要的對照 |
# |---|---|---|
# | 長期近線性趨勢 | 區域地表速度 | 參考框架、跨站速度場 |
# | 突然偏移 | 同震位移或裝置階變 | 事件時間、維運紀錄、鄰站 |
# | 震後逐漸變緩 | 餘滑或其他鬆弛過程 | 更長序列、空間分布與模型 |
# | 季節性往返 | 地表負載等週期效應 | 降雨、水文及背景資料 |
#
# 每日解能比較震前後座標是否不同，卻不能從一天一點的圖判定幾秒內如何
# 移動。研究快速同震過程需要高頻 GNSS 解及相應處理；原始觀測每三十秒
# 一筆，也不等於已經有三十秒精密位移。
#
# 慢滑移可以在沒有相應強烈地震波的情況下逐漸累積位移。GNSS 是重要觀測
# 方式之一，但還有其他方式；應變、傾斜或其他地球物理資料也能提供約束。
# 是否看得到，取決於滑移位置、幅度、時間尺度及觀測網的敏感度。
#
# 本章沒有建立花蓮事件的精密 GNSS 解，因此下一章的共同時間軸只比較
# 實際具備的波形、目錄、地下水及地磁，不以近似座標補出一條不存在的
# 位移曲線。這個區分讓{doc}`花蓮案例 <07_case_hualien2024>`可以清楚
# 說明哪些證據已在手上，哪些需要另外取得。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **資料實作・免費網站**：Nevada Geodetic Laboratory，〈[Plug and Play GPS Data Products](https://geodesy.unr.edu/PlugNPlayPortal.php)〉。
#   從測站清單進入位置時間序列與資料格式說明，延伸本章「原始衛星觀測」和「已解算座標」的區別；比較位移前先確認參考框架及單位。
# - **格式查詢・免費檔案**：IGS／RTCM RINEX Working Group，〈[RINEX](https://igs.org/wg/rinex/)〉。
#   依手上檔案的版本選擇規格，查閱標頭與觀測量代碼；詳細解析移到附錄 G，檔頭的近似座標不等於逐時刻解算的位置序列。
# - **資料方法導讀・免費文章**：Blewitt, G., Hammond, W. C., & Kreemer, C.（2018），〈[Harnessing the GPS Data Explosion for Interdisciplinary Science](https://doi.org/10.1029/2018EO104623)〉，*Eos*；[免費全文](https://eos.org/science-updates/harnessing-the-gps-data-explosion-for-interdisciplinary-science)。
#   瞭解大量 GNSS 資料如何整理成可研究的速度與位移產品，特別看季節變化、地震階變及裝置更動如何影響時間序列的解讀。
