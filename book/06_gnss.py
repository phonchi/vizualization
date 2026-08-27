# %% [markdown]
# # 6. GNSS 地殼變形觀測
#
# ## 6.1 為什麼 GNSS 與地震有關？
#
# GNSS（GPS 等衛星定位系統）連續站以毫米級精度追蹤測站位置：
#
# - **震間（inter-seismic）**：板塊持續推擠，台灣以每年約 8 公分的速度
#   縮短——應變正在累積。
# - **同震（co-seismic）**：地震瞬間測站「跳」到新位置，位移場直接
#   反映斷層滑移分布。
# - **震後（post-seismic）**：餘滑與黏彈性鬆弛造成的緩慢衰減運動。
# - **慢地震（SSE）**：不放地震波的緩慢滑移事件，只有 GNSS 看得到。
#
# ## 6.2 GDMS 提供的是「原始觀測檔」（RINEX）
#
# 這裡要先建立一個重要觀念：GDMS 釋出的 GNSS 資料是 **RINEX 觀測檔**
# ——接收儀記錄的衛星訊號原始觀測量（虛擬距離、載波相位），
# **不是**現成的位置時間序列。從 RINEX 到毫米級座標需要專業軟體
# （GAMIT/GLOBK、Bernese、PRIDE PPP-AR 等）解算。
#
# 本章目標：看懂 RINEX 檔案的結構與內容。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %%
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
# 檔名規則（RINEX 2）：`hual0930.24o.gz` → 測站 `hual`（花蓮）、
# 年積日 `093`（2024 年第 93 天 = 4/2）、`24o` = 2024 年觀測檔。
#
# ## 6.3 打開一個 RINEX 檔看看

# %%
with tarfile.open(tgz) as tar:
    member = next(m for m in tar.getmembers() if "hual0930" in m.name)
    raw = gzip.decompress(tar.extractfile(member).read()).decode("ascii")
lines = raw.splitlines()
print("\n".join(lines[:20]))

# %% [markdown]
# 標頭裡的重要欄位：`APPROX POSITION XYZ`（測站概略座標，地心直角座標）、
# `TYPES OF OBSERV`（記錄的觀測量：C1 虛擬距離、L1/L2 載波相位…）、
# `INTERVAL`（取樣間隔 30 秒）。
#
# 標頭之後是一個個「曆元（epoch）」：每 30 秒一筆，記錄當下收到哪些
# 衛星（G=GPS、R=GLONASS）、各觀測量的值。
#
# ## 6.4 簡單解析：每個曆元收到幾顆衛星？

# %%
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
# 衛星數隨星座幾何在 8–20 顆之間起伏。衛星數與幾何分布（DOP）
# 直接影響定位精度。
#
# ## 6.5 從標頭座標反推測站位置
#
# 標頭裡的 `APPROX POSITION XYZ` 是測站的地心直角座標（ECEF），單位公尺。
# 把它轉成我們熟悉的經緯度，就能和第 1 章的測站清單對照，確認自己讀對了
# 檔案。轉換公式是大地測量的標準流程（WGS84 橢球，迭代求緯度）：

# %%
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
# 從 RINEX 標頭算出的位置，和測站清單裡登錄的經緯度吻合到小數點後好幾位。
# 這種交叉檢查看似瑣碎，卻是處理陌生資料格式時最實在的一步：先確認自己
# 讀進來的東西「位置對得上」，再談後面的分析。
#
# ## 6.6 想要位置時間序列？用解算好的成果
#
# 要看同震位移，需要的是每日一點的座標序列，不是原始觀測檔。學術研究
# 通常直接使用專業機構解算好的成果：
#
# - **中研院地球所 GPS Lab**（<https://gps.earth.sinica.edu.tw/>）：台灣
#   全網每日座標解，可申請下載。
# - **Nevada Geodetic Laboratory**（<https://geodesy.unr.edu/>）：全球
#   （含台灣多站）每日解，網頁直接可抓。
#
# 拿到解算後的東西向、南北向、垂直向三分量序列（單位毫米，一天一點），
# 就能看到 2024/4/3 那天花蓮周邊測站在幾秒內跳了數十公分。從原始觀測檔
# 自己解算到毫米級座標，需要 GAMIT/GLOBK、Bernese 或 PRIDE PPP-AR 這類
# 專業軟體，超出本課範圍，但知道路在哪裡，需要時才找得到。
#
# ```{admonition} 延伸工具
# :class: tip
# 想用 Python 讀 RINEX，可以看 `georinex`；想自己嘗試 PPP 精密單點定位，
# 開源的 PRIDE PPP-AR 是不錯的起點。
# ```
