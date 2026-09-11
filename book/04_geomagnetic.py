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
# # 21. 地磁觀測：共同背景與局部變化
#
# 地下水讓我們看到，一種觀測可以同時回應多個原因。地磁的情況更明顯：
# 磁力儀量到的是地表總磁場，其中包括地核、地殼，以及電離層與磁層電流的
# 貢獻。測站就在花蓮，不表示曲線的主要變化也來自花蓮地下。
#
# ## 21.1 從來源理解時間尺度
#
# 主磁場的長期變化、岩石磁化造成的空間差異，以及日變化和磁擾，會疊加在
# 同一份紀錄。地震相關研究另外提出應力改變岩石磁化、地下流體造成電流等
# 機制。機制上可能產生磁場，不代表每個測站都量得到，更不代表訊號必然
# 在地震之前出現。Johnston（1997）的綜述因此同時討論物理來源和觀測限制。
#
# 新城（XCG）的紀錄呈現日常背景，多站比較與磁場分量比值則讓我們進一步
# 辨認資料中的變化。這個順序讓候選解釋先有參照，不必
# 一開始就把起伏分成「正常」和「前兆」。
#
# ## 21.2 測站與分量
#
# 本書使用的 IAGA-2002 資料含 X（向北）、Y（向東）、Z（垂直向下）分量。
# 下表是教材儲存的測站清單。分析前先確認方向、單位與有效欄位，站數不能
# 當作本時段全部可用秒資料的數量。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import gdms_toolkit as gt
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import PALETTE, QUAKE_COLOR, apply_layout

mags = gt.load_stations("MAGNET")
mags[["station_code", "chinese_station_name", "lat", "lon", "location_county"]]

# %% tags=["remove-input"]
tgz = CACHE_DIR / "edu-mag-hualien2024.tgz"
mag = gt.read_geomagnetic(tgz, "XCG", resample="1min")
mag.describe().round(1)

# %% [markdown]
# ## 21.3 三分量不一定一起變化
#
# 此處使用一分鐘平均的新城紀錄。全磁力可由三分量平方和開根號計算，但
# 推算值與獨立感測器量到的全磁力不是同一項觀測；本書以 `F_calc` 明確區分。
# 先看 X、Y、Z 各自的日常變化，紅線只標示主震時刻。
#
# %% tags=["remove-input"]
comps = [("X", "X（北向，nT）"), ("Y", "Y（東向，nT）"), ("Z", "Z（垂直，nT）")]
fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                    subplot_titles=[t for _, t in comps])
for i, (c, _) in enumerate(comps, 1):
    fig.add_trace(go.Scattergl(x=mag.index, y=mag[c], name=c,
                               line=dict(color=PALETTE[i - 1], width=1.2)), i, 1)
eq_t = "2024-04-02 23:58:11"
for row in (1, 2, 3):
    fig.add_vline(x=eq_t, line_dash="dash", line_color=QUAKE_COLOR, row=row, col=1)
apply_layout(fig, height=640, showlegend=False,
             title="新城（XCG）地磁三分量，1 分鐘平均（UTC）")
fig

# %% [markdown]
# 每天相似的起伏，與日變化背景相容；不同分量的形狀未必相同，因為磁場
# 是向量。較不規則的日期值得再與 GFZ 的 Kp 指數對照，但要判定某個起伏
# 是否為磁暴，不能只憑曲線振幅。Kp 描述大範圍地磁活動，也不是本站每一段
# 局部干擾的完整紀錄。
#
# ## 21.4 把不同日期對齊
#
# 下面以 UTC 時刻對齊每日 X 分量，並扣掉各日平均值。這會凸顯日內形狀，
# 同時移除每日平均的差異。如果研究的是持續數天的偏移，就不能只用這張圖，
# 因為處理本身已經拿掉了部分目標訊號。
#
# %% tags=["remove-input"]
z = mag["X"].copy()
daily = z.groupby(z.index.date)
fig = go.Figure()
for i, (day, s) in enumerate(daily):
    hours = s.index.hour + s.index.minute / 60
    fig.add_trace(go.Scattergl(
        x=hours, y=s - s.mean(), mode="lines", name=str(day),
        line=dict(width=1, color="rgba(42,120,214,0.35)"), showlegend=False))
apply_layout(fig, title="每日 X 分量疊圖（各日去平均）：日變化的形狀",
             xaxis_title="UTC 時（台灣時間 = UTC+8）",
             yaxis_title="ΔX（nT）", hovermode=False)
fig

# %% [markdown]
# 多日曲線的相似部分勾勒出背景；偏離背景的日期則提示還有其他
# 來源。它們可以引導後續查證，卻尚不是對地震原因的判定。
#
# ## 21.5 參考站讓共同成分有機會被辨認
#
# 若兩站共享相近的外源場，相減能減弱共同變化。本書先用池上（CSG）當
# 新城的參考站。效果取決於外源場的空間差異、地方地下導電構造、儀器與
# 取樣一致性，因此相減後留下的不一定全是局部地殼訊號。
#
# %% tags=["remove-input"]
mag_csg = gt.read_geomagnetic(CACHE_DIR / "edu-mag-csg.tgz", "CSG",
                              resample="1min")
diff = (mag["F_calc"] - mag_csg["F_calc"]).dropna()
fig = go.Figure(go.Scattergl(x=diff.index, y=diff - diff.mean(), mode="lines",
                             line=dict(color=PALETTE[6], width=1.2)))
fig.add_vline(x=eq_t, line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="全磁力兩站差：新城（XCG）－ 池上（CSG），去平均",
             yaxis_title="ΔF（nT）", showlegend=False)
fig

# %% [markdown]
# 比較差分圖與原始圖時，先看哪些變化被減弱、哪些仍在。某段起伏若在兩站
# 都有，差分後較小是預期結果；若只在一站有，差分仍會保留它。但「只在一站」
# 可能來自地下局部來源，也可能是附近的人為干擾，還不足以分辨兩者。
#
# ## 21.6 參考站選擇本身也是模型選擇
#
# 換用卑南（TTN）後，殘差散佈可能改變。這反映了背景估計的不確定性：
# 不同參考站未必得到相同結果。比較時也要注意兩站共同有資料的時段，
# 否則標準差差異可能只是使用了不同日期。
#
# %% tags=["remove-input"]
mag_ttn = gt.read_geomagnetic(CACHE_DIR / "edu-mag-ttn.tgz", "TTN",
                              resample="1min")
diff_ttn = (mag["F_calc"] - mag_ttn["F_calc"]).dropna()
print(f"減 CSG 後殘差標準差：{(diff - diff.mean()).std():.2f} nT")
print(f"減 TTN 後殘差標準差：{(diff_ttn - diff_ttn.mean()).std():.2f} nT")

# %% [markdown]
# 較小的殘差表示這個差分組合在該資料窗中更平穩，卻不能證明背景已完全
# 去除；目標訊號也可能被共同減掉。若在看過地震時間之後反覆挑選參考站，
# 直到得到漂亮的候選異常，便增加了選擇偏差。較好的做法是在訓練時段決定
# 比較規則，再用其他時段檢查效果。
#
# ## 21.7 同一個比值，可以回答不同問題
#
# 最後看看垂直與水平磁場的比值 Z/H。下面直接從原始場分量算出每日中位數，
# 主要反映磁場方向及其緩慢變化。它容易計算，卻不是文獻裡先選定頻帶、
# 再比較擾動功率的 ULF 極化指標。先分清這個差異，才能知道這張圖支持
# 哪些問題。
#
# %% tags=["remove-input"]
zh = (mag.Z / mag.H)
zh_daily = zh.groupby(zh.index.date).median()
fig = go.Figure(go.Scatter(x=list(map(str, zh_daily.index)), y=zh_daily.values,
                           mode="lines+markers",
                           line=dict(color=PALETTE[6], width=2)))
fig.add_vline(x="2024-04-03", line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="新城（XCG）每日 Z/H 中位數",
             xaxis_title="日期", yaxis_title="Z / H", showlegend=False)
fig

# %% [markdown]
# 若 Z/H 大致穩定，表示本圖呈現的場方向沒有很大的每日變化；它不能直接
# 排除特定頻帶內的微弱擾動。頻帶分析要先定義觀測量、去趨勢與濾波
# 方式，再檢查每個頻帶的背景和不確定性。
#
# 取樣率也設定了界線。一秒取樣的 Nyquist 頻率是 0.5 Hz；一分鐘取樣則是
# 約 0.0083 Hz。後者仍可能保留更慢的變化，但無法完整分析較高頻率的 ULF
# 訊號。平均與降取樣怎樣改變資訊，見{doc}`附錄 G <appendix_g_observations>`。
#
# USGS 對若干知名電磁前兆宣稱的重新檢查，指出全球磁場變化與資料問題可以
# 造成誤判。這些案例提醒我們，需要把來源識別與預報能力分別驗證。即使確認
# 了同震磁場變化，也只能先說明地震伴隨的效應，不能倒推成震前預報證據。
#
# {doc}`下一章 <05_seismic>`回到與地震最直接相關的震波，看看連續訊號
# 如何整理成前一部模型使用的事件目錄。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **背景查詢・免費網站**：GFZ，〈[Kp index](https://kp.gfz.de/en/)〉。
#   查閱全球地磁活動指數與說明，對照本章的擾動日；先檢查太空天氣背景，再討論單站或兩站差值中的局部異常。
# - **入門・免費官方文章**：USGS Geomagnetism Program，〈[Overview](https://www.usgs.gov/programs/geomagnetism/science/overview)〉。
#   閱讀 USGS 對若干磁場前兆宣稱的重新檢驗，瞭解全球背景與資料問題如何影響來源判斷；這不是對所有候選效應的單一總結。
# - **核心綜述・全文可能需訂閱**：Johnston, M. J. S.（1997），〈[Review of electric and magnetic fields accompanying seismic and volcanic activity](https://doi.org/10.1023/A:1006500408086)〉，*Surveys in Geophysics*；[USGS 免費摘要](https://www.usgs.gov/publications/review-electric-and-magnetic-fields-accompanying-seismic-and-volcanic-activity)。
#   對照壓磁效應、流體相關電磁效應與觀測頻帶，特別區分同震訊號和震前訊號；看到地震伴隨的磁場變化，不等於能據此提前預測地震。
