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
# # 20. 地下水位：含水層如何回應外界擾動
#
# 上一章提醒我們，資料欄位不等於它背後的物理原因。地下水位正是一個例子：
# 井中的水面升高，可能因降雨補注、氣壓降低、抽水停止，也可能與地殼應變或
# 地震引起的滲透性改變有關。要研究地震，必須先理解平常的水位為什麼會動。
#
# ## 20.1 一口井如何感受到地震
#
# 在合適的承壓含水層與水力連通條件下，井水位變化能反映地下壓力水頭的變化。
# 含水層受壓縮時，孔隙與流體的反應可能使水位上升；但開放含水層、井內儲水、
# 排水速度及局部地質會改變這個對應。因此不能把所有地下水位都當作孔隙壓力
# 的直接量測，更不能把水位公分數直接換成地表位移公分數。
#
# 震波經過時可能造成短暫振盪；震後也可能出現較持久的水位偏移。兩者涉及的
# 時間尺度不同，未必同時出現。USGS 的地下水觀測整理顯示，不同井對同一場
# 地震可以有不同反應，甚至在離震源很遠的地方也可能受震波擾動。
#
# 在討論地震之前，我們先找兩種平常就存在的激勵：潮汐與氣壓。知道輸入何時
# 改變，觀察水位如何跟著變，有助於理解這口井對外界擾動的敏感度。
#
# ## 20.2 從測站清單認識觀測井
#
# 下表列本書資料中的觀測井。站碼能把時間序列接回位置與沿革，但含水層性質
# 仍需井的地質與施工資料，不能只用經緯度推斷。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import gdms_toolkit as gt
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import PALETTE, ACCENT, QUAKE_COLOR, apply_layout

wells = gt.load_stations("GW")
wells[["station_code", "chinese_station_name", "lat", "lon",
       "elevation", "location_county", "start_time"]]

# %% [markdown]
# ## 20.3 先看一段完整的日常變化
#
# 以下使用本書儲存的 2024 年 3–4 月壯圍站（TUN）資料。原始取樣間隔為
# 一秒，長時間圖先用一分鐘平均，讓日夜尺度的起伏容易閱讀。平均會減弱
# 短暫振盪，因此這張圖適合看背景，不能替代同震短窗分析。
#
# 同一份資料同時有水位、氣壓與水溫。把三者並列，是為了讓可能的解釋有機會
# 互相競爭：水位變了，當時氣壓是否也變？只有一欄突然跳動，還是所有欄位
# 一起中斷？
#
# %% tags=["remove-input", "remove-output"]
tgz = CACHE_DIR / "edu-gw-hualien2024.tgz"
gw = gt.read_groundwater(tgz, "TUN", resample="1min")
gw.info()
gw.head(3)

# %% [markdown]
# ## 20.4 水位與氣壓放在一起看
#
# 先追蹤重複出現的起伏，再看較慢的趨勢。紅線標示花蓮主震的發震時刻；
# 先用整段背景判斷紅線附近的變動是否與平常不同，再討論它是否為地震
# 反應。水位、氣壓與水溫各自使用不同縱軸，不能按曲線高度比較強弱。
#
# %% tags=["remove-input"]
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                    subplot_titles=("水位（cm）", "氣壓（hPa）"))
fig.add_trace(go.Scattergl(x=gw.index, y=gw.water_level_cm, name="水位",
                           line=dict(color=PALETTE[0], width=1.5)), 1, 1)
fig.add_trace(go.Scattergl(x=gw.index, y=gw.pressure_hPa, name="氣壓",
                           line=dict(color=PALETTE[1], width=1.5)), 2, 1)
# 2024/4/3 花蓮地震（UTC 23:58 於 4/2）
eq_t = "2024-04-02 23:58:11"
for row in (1, 2):
    fig.add_vline(x=eq_t, line_dash="dash", line_color=QUAKE_COLOR, row=row, col=1)
fig.add_annotation(x=eq_t, y=1, yref="paper", text="M7.2 花蓮地震",
                   showarrow=False, font=dict(color=QUAKE_COLOR), xanchor="left")
apply_layout(fig, height=520, title="壯圍（TUN）觀測井，2024/03–04（時間為 UTC）",
             showlegend=False)
fig

# %% [markdown]
# 若水位與氣壓在某些時段呈反向變化，這與承壓含水層的氣壓反應相容，
# 但光憑視覺同步還不能估出一個通用的物理係數。降雨、抽水、潮汐和延遲反應
# 可能同時存在。後面的簡單迴歸是一個描述起點，完整分析需要把這些因素
# 及其不確定性一起考慮。
#
# ## 20.5 規律起伏提供了什麼資訊
#
# 放大到一週後，可以辨認接近半日尺度的起伏。重複訊號讓我們能
# 多次比較：同樣的外界激勵，水位的振幅與延遲是否穩定？若要把變化歸於
# 地震，先要建立這種平常反應的範圍。
#
# %% tags=["remove-input"]
week = gw.loc["2024-03-10":"2024-03-17"]
fig = px.line(week, y="water_level_cm",
              labels={"Time": "時間（UTC）", "water_level_cm": "水位（cm）"},
              color_discrete_sequence=[ACCENT])
fig.update_traces(line_width=1.5)
apply_layout(fig, title="一週的水位：接近半日尺度的起伏", showlegend=False)
fig

# %% [markdown]
# 接近一天兩次的起伏提示半日週期成分，潮汐是其中的重要候選來源。
# 不過 M2 月球半日潮與 S2 太陽半日潮的頻率很近，氣壓也可能有日週期；
# 僅憑一週曲線數出兩個波峰，還不能把各來源分離。
#
# ## 20.6 頻譜把「多久重複一次」放到橫軸
#
# 下面將同一條水位序列表示成不同頻率的功率。橫軸改成週期後，越靠近
# 十二小時的峰，表示序列含有越明顯的半日尺度變化。這是時間圖的另一種
# 摘要，不是增加了新的觀測；峰值的寬度與高度還會受資料長度、缺測處理
# 和估計設定影響。
#
# %% tags=["remove-input"]
import numpy as np
from scipy import signal

x = gw.water_level_cm.interpolate().dropna()
freq, psd = signal.welch((x - x.mean()).to_numpy(), fs=1 / 60, nperseg=2 ** 14)
period_hr = 1 / freq[1:] / 3600
fig = go.Figure(go.Scattergl(x=period_hr, y=psd[1:], mode="lines",
                             line=dict(color=ACCENT, width=1.5)))
for p, name in [(12.42, "M2 太陰半日潮"), (12.0, "S2 太陽半日潮"),
                (25.82, "O1"), (24.0, "K1/S1")]:
    fig.add_vline(x=p, line_dash="dot", line_color="#999")
    fig.add_annotation(x=np.log10(p), y=1.02, yref="paper", text=name,
                       showarrow=False, font=dict(size=11), textangle=-30)
apply_layout(fig, title="水位功率譜（Welch 法）",
             xaxis=dict(type="log", title="週期（小時）",
                        range=[np.log10(3), np.log10(200)]),
             yaxis=dict(type="log", title="PSD"), hovermode="x")
fig

# %% [markdown]
# 接近十二與二十四小時的峰，值得與理論潮汐及氣壓頻譜對照。頻率相近
# 不等於原因已確定，還要比較相位、跨時段穩定性及外界激勵。若要把水位
# 換成含水層應變，需要另外校準水文地質響應，不能把幾公分水位變化直接
# 解讀成幾公分地殼形變。
#
# ## 20.7 不同井的差異也有科學意義
#
# 以下用相同的頻譜設定比較壯圍與赤山（CHI）。顯示的是半日頻帶內的
# 峰值指標，便於相對比較；它不是經過完整潮汐分析估得的 M2 振幅。頻譜
# 指標的單位與限制見{doc}`附錄 G <appendix_g_observations>`。
#
# %% tags=["remove-input"]
def m2_amplitude(df):
    """半日頻帶 PSD 峰值的平方根（cm/√Hz），不是水位振幅。"""
    x = df.water_level_cm.interpolate().dropna().to_numpy()
    f, p = signal.welch(x - x.mean(), fs=1 / 60, nperseg=2 ** 14)
    per = 1 / f[1:] / 3600
    band = (per > 11.5) & (per < 13)
    return float(np.sqrt(p[1:][band].max()))

chi = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "CHI",
                          resample="1min")
for name, df in [("壯圍 TUN", gw), ("赤山 CHI", chi)]:
    print(f"{name}：半日頻帶峰值指標（cm/√Hz） = {m2_amplitude(df):.1f}")

# %% [markdown]
# 兩口井的指標不同，可能與含水層的承壓程度、透水性、井內儲水、背景干擾
# 及儀器條件有關。若想分辨這些原因，下一步應比較已知潮汐輸入與水位輸出，
# 而非只用較高的峰值判定哪口井比較好。跨站差異能幫助設計研究，但不能
# 單靠一個數字反推出含水層是否封閉。
#
# ## 20.8 用迴歸描述氣壓反應
#
# 這裡先將相鄰有效紀錄取差分，再迴歸水位變化對氣壓變化。這會弱化部分
# 慢趨勢，但不保證排除潮汐，也沒有建模延遲效應。輸出的斜率單位是 cm/hPa，
# 應稱為氣壓響應斜率；水文地質學常用的無量綱氣壓效率，還需壓力水頭換算
# 及符號約定。
#
# %% tags=["remove-input"]
def barometric_efficiency(df):
    d = df[["water_level_cm", "pressure_hPa"]].dropna().diff().dropna()
    slope = np.polyfit(d.pressure_hPa, d.water_level_cm, 1)[0]
    return float(slope)

for name, df in [("壯圍 TUN", gw), ("赤山 CHI", chi)]:
    be = barometric_efficiency(df)
    print(f"{name}：氣壓響應斜率 = {be:+.3f} cm/hPa")

# %% [markdown]
# 氣壓反應可以幫助解釋水位，也能提醒我們模型有多簡化。如果背景關係隨
# 季節改變，拿一段時間的斜率去扣除全年資料，殘差就可能出現系統性偏差。
# 因此「扣除背景」應理解為建立可檢查的背景模型，而非保證剩下的全是地震。
#
# ## 20.9 缺測欄位與仍在工作的感測器
#
# 花蓮井（HWA）在這段教材資料中的水位欄為缺值，但水溫仍有紀錄。下圖保留
# 水溫，是為了區分欄位的可用性。只看缺值無法確定故障原因，必須再查維運
# 或資料釋出紀錄；它更不能證明那段時間水位沒有反應。
#
# %% tags=["remove-input"]
hwa = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "HWA",
                          resample="10min")
print("水位有效點數：", hwa.water_level_cm.notna().sum())
print("上層水溫有效點數：", hwa.temp_upper_C.notna().sum())
fig = px.line(hwa.dropna(subset=["temp_upper_C"]), y="temp_upper_C",
              labels={"Time": "時間（UTC）", "temp_upper_C": "上層水溫（°C）"},
              color_discrete_sequence=[PALETTE[1]])
fig.update_traces(line_width=1.5)
apply_layout(fig, title="花蓮（HWA）井：水位欄缺測，水溫仍有紀錄",
             showlegend=False)
fig

# %% [markdown]
# 地下水研究也能回答前兆以外的問題。辨認一口井的潮汐與氣壓反應、瞭解
# 地震後含水層是否改變、量化觀測不到訊號時的上限，都能增加對地下系統的
# 認識。這些目標需要不同資料窗，也支持不同程度的結論。
#
# 接下來的{doc}`地磁觀測 <04_geomagnetic>`會遇到相似的分離問題：一條
# 看似反映局部變化的曲線，可能主要由遠方的共同來源驅動。那時多站比較會比單站放大
# 更有幫助。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **入門・免費網站**：USGS，〈[How does an earthquake affect groundwater levels and water quality in wells?](https://www.usgs.gov/faqs/how-does-earthquake-affect-groundwater-levels-and-water-quality-wells)〉。
#   先分清震波經過時的水位振盪、震後水位偏移與水質變化，再回頭判斷本章各口井究竟觀測到了哪一種反應。
# - **觀測機制・免費官方教材**：USGS（2003），〈[Earthquakes—Rattling the Earth's Plumbing System](https://pubs.usgs.gov/fs/fs-096-03/)〉，Fact Sheet 096-03。
#   用不同地點的觀測說明地震如何擾動地下水，適合延伸本章「同一場地震，不同井反應不同」的討論。
# - **核心綜述・全文需訂閱**：Roeloffs, E. A.（1988），〈[Hydrologic precursors to earthquakes: A review](https://doi.org/10.1007/BF00878996)〉，*Pure and Applied Geophysics*。
#   重點讀潮汐響應、承壓含水層與應變的關係，以及如何排除氣壓、降雨與抽水影響；文獻整理的候選前兆不能直接當成已驗證的預測方法。
#
# - **頻譜單位・免費檔案**：SciPy 開發團隊，〈[welch](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html)〉。對照本章半日頻帶指標，特別檢視功率譜密度與功率譜的單位差異；PSD 峰值開根號不能直接當作水位振幅。
