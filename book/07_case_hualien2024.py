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
# # 24. 花蓮案例：同一場地震，不同觀測回答什麼
#
# 前面分別讀過各類觀測，現在把它們放回 2024 年 4 月 3 日花蓮 $M_L\,7.2$
# 地震。當地發震時間約為上午 7 時 58 分，對應 UTC 4 月 2 日 23 時 58 分。
# 這個時間轉換是共同時間軸的起點，但不同站的震波到時仍會晚於發震時刻。
#
# 本章不把四類資料當作四張獨立的贊成票。先問每一類資料理應對什麼過程
# 敏感，再看它是否在相應時間尺度出現可辨認的變化。地震目錄與波形並非
# 彼此獨立，因為目錄本來就由波形建構；地下水和地磁也有各自共享的背景。
#
# | 資料 | 本章測站／範圍 | 主要比較對象 |
# |---|---|---|
# | 地震目錄 | 本書春季目錄 | 主震後事件的時空分布 |
# | 波形 | HWA 花蓮站 | 起振、持續時間及分量 |
# | 地下水 | TUN 壯圍等井 | 振盪與較持久的水位偏移 |
# | 地磁 | XCG 新城 | 主震時窗與日常背景 |
#
# 站點位置不同，所以水位和磁場不能被視為震源處的直接量測。GNSS 需要
# 另有已解算位移產品才可加入，本章不包含這項成果。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import gdms_toolkit as gt
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import PALETTE, ACCENT, QUAKE_COLOR, SEQUENTIAL, apply_layout

EQ_UTC = pd.Timestamp("2024-04-02 23:58:11", tz="UTC")

# %% [markdown]
# ## 24.1 從事件分布建立背景
#
# 先看主震後三十天的目錄點，顏色表示距主震的天數。這是一個時間切片，
# 不是以親代機率判定出的 ETAS 家族；圖中也可能包含區域背景活動。
# 把它稱為主震後事件分布，比先假設所有點都由主震觸發更精確。
#
# %% tags=["remove-input"]
cat_csv = CACHE_DIR / "catalog_2024spring.csv"
if cat_csv.exists():
    cat = pd.read_csv(cat_csv, parse_dates=["time"])
else:
    cat = gt.gdms_earthquake_catalog(gt.GDMSSession(), "2024-03-01", "2024-06-30",
                                     min_ml=3)
    cat.to_csv(cat_csv, index=False)
main = cat.loc[cat.ML.idxmax()]
print(f"主震：{main.time} UTC（當地時間 4/3 07:58） ML {main.ML} "
      f"深度 {main.depth} km ({main.latitude:.3f}N, {main.longitude:.3f}E)")

# 餘震時空分布：主震後 30 天
aft = cat[(cat.time >= main.time) &
          (cat.time <= main.time + pd.Timedelta(days=30))].copy()
aft["days"] = (aft.time - main.time).dt.total_seconds() / 86400
fig = px.scatter_map(aft, lat="latitude", lon="longitude",
                     size=2 ** aft.ML / 8, color="days",
                     color_continuous_scale=SEQUENTIAL,
                     zoom=7, center=dict(lat=24.0, lon=121.6),
                     map_style="carto-positron", height=560,
                     labels={"days": "主震後天數"})
apply_layout(fig, title="主震後 30 天事件分布（顏色＝發生時間）", hovermode="closest")
fig

# %% [markdown]
# 密集點帶提示活動區域，也顯示後續事件未必集中在單一位置。要將它解釋
# 為某一斷層上的破裂，還需深度、震源機制與重定位結果。Zheng et al.（2024）
# 的研究進一步結合波形與大地測量資料建立震源模型，示範的是更完整的
# 推論，不是單憑震央圖描線。
#
# ## 24.2 波形讓事件時間具體可見
#
# 三分量波形在同一時間軸上顯示強烈震動。主震紅線與到站訊號之間的差異
# 包含傳播時間；繪圖為了瀏覽而減少顯示點數，不適合從這張圖讀取精確峰值
# 或逐樣本到時。物理振幅仍需儀器響應校正。
#
# %% tags=["remove-input"]
st = gt.read_waveform(CACHE_DIR / "edu-wave-hualien2024.mseed")
fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                    subplot_titles=[tr.id for tr in st])
for i, tr in enumerate(st, 1):
    x = pd.to_datetime(tr.times("timestamp"), unit="s", utc=True)
    step = 10  # 100 Hz → 10 Hz 顯示用降頻
    fig.add_trace(go.Scattergl(x=x[::step], y=tr.data[::step], mode="lines",
                               line=dict(color=PALETTE[i - 1], width=0.8)), i, 1)
for row in (1, 2, 3):
    fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR,
                  row=row, col=1)
apply_layout(fig, height=600, showlegend=False, hovermode=False,
             title="花蓮站三分量波形（UTC；紅線＝發震時刻）")
fig

# %% [markdown]
# ## 24.3 地下水：先分清振盪與偏移
#
# 地下水可能在震波經過時振盪，也可能震後停在不同水位，或逐漸恢復。
# 這些形狀涉及不同機制與頻率，應分別檢查。只比較震前震後平均值，可能
# 漏掉短暫振盪；只看短窗差分，則可能漏掉持續偏移。
#
# 先用一分鐘平均的壯圍資料看多日背景。若日常潮汐與氣壓造成的變化已經
# 很大，就必須把候選同震反應與這些背景一起比較，不能只因它接近紅線
# 就認定有關。
#
# %% tags=["remove-input"]
gw = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "TUN",
                         start="2024-03-30", end="2024-04-06", resample="1min")
fig = go.Figure(go.Scattergl(x=gw.index, y=gw.water_level_cm, mode="lines",
                             line=dict(color=PALETTE[0], width=1.5)))
fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="壯圍（TUN）水位，主震前後各三天",
             yaxis_title="水位（cm）", showlegend=False)
fig

# %% [markdown]
# 在這個多日縱軸下，主震附近未呈現明顯的大幅水位階變。但一分鐘平均
# 會減弱短暫振盪，這個判讀只適用於目前顯示尺度。接著回到一秒資料，
# 檢查主震前五分鐘到後十五分鐘的短窗。
#
# %% tags=["remove-input"]
gw_s = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "TUN",
                           start="2024-04-02", end="2024-04-03", resample=None)
win = gw_s.loc[EQ_UTC - pd.Timedelta("5min"): EQ_UTC + pd.Timedelta("15min")]
fig = go.Figure(go.Scattergl(x=win.index, y=win.water_level_cm, mode="lines",
                             line=dict(color=PALETTE[0], width=1.5)))
fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="壯圍（TUN）水位，主震前後 20 分鐘（1 秒取樣）",
             yaxis_title="水位（cm）", showlegend=False)
fig

# %% [markdown]
# 這段短窗沒有呈現像地震波形那樣明顯的瞬間擾動。這是對本井、本資料窗
# 及本圖解析度的描述，還不是「水位完全沒有反應」的統計證明。
#
# 下面另外比較四口井的一秒差分標準差：震前使用十分鐘、震後使用五分鐘。
# 比值可描述短時間波動是否增加，但兩窗長度不同，且觀測有時間相關；
# 不能把比值接近一直接當成接受「無反應」假設的檢定。
#
# %% tags=["remove-input"]
def wobble(station):
    s = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", station,
                            start="2024-04-02", end="2024-04-03",
                            resample=None).water_level_cm.diff()
    before = s.loc[EQ_UTC - pd.Timedelta("10min"): EQ_UTC].std()
    after = s.loc[EQ_UTC: EQ_UTC + pd.Timedelta("5min")].std()
    return before, after

for st_code in ["TUN", "DON", "NAB", "CHI"]:
    b, a = wobble(st_code)
    print(f"{st_code}：震前 std={b:.4f}  震後 std={a:.4f}  比值={a / b:.1f}×")

# %% [markdown]
# 若比值沒有明顯放大，只表示這項摘要未顯示突出的高頻波動增加。
# 持續水位偏移主要影響階變附近的一次差分，較慢的反應也可能被這個指標
# 忽略。因此應將它與原始曲線一起看；完整的同震反應分析仍不可省略。
#
# 要給出可檢測上限，可以在具有相似背景的非事件時窗注入不同幅度與持續
# 時間的模擬訊號，再用固定流程檢查檢出率。這一步目前尚未執行，所以本章
# 不聲稱已量出最小可偵測水位或排除了某種物理機制。井的敏感度和實際
# 地殼應變也未在這裡校準，不能只用震央距離替它們指定數值。
#
# ## 24.4 地磁：把紅線與日常背景放在一起
#
# 下圖仍使用一分鐘平均。總磁力與垂直分量的日常起伏提供背景輪廓；它們
# 來自同一組分量紀錄，不能算成兩次獨立的異常證據。
#
# %% tags=["remove-input"]
mag = gt.read_geomagnetic(CACHE_DIR / "edu-mag-hualien2024.tgz", "XCG",
                          start="2024-03-30", end="2024-04-06", resample="1min")
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                    subplot_titles=("全磁力 F（由三分量計算）", "垂直分量 Z"))
fig.add_trace(go.Scattergl(x=mag.index, y=mag.F_calc, mode="lines",
                           line=dict(color=PALETTE[6], width=1.2)), 1, 1)
fig.add_trace(go.Scattergl(x=mag.index, y=mag.Z, mode="lines",
                           line=dict(color=PALETTE[2], width=1.2)), 2, 1)
for row in (1, 2):
    fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR,
                  row=row, col=1)
apply_layout(fig, height=520, showlegend=False,
             title="新城（XCG）地磁，主震前後各三天（nT）")
fig

# %% [markdown]
# 目前的圖未顯示一個可以直接歸因於地震的明顯變化。要再往前判斷，應沿用
# 地磁章的方法：查太空天氣、比較參考站、指定頻帶，並考慮儀器與局部干擾的
# 可能。本圖沒有完成這些識別步驟，因此也不宣稱排除所有震磁效應。
#
# ## 24.5 共同時間軸能比較什麼
#
# 四列圖把發震前後的紀錄對齊，但波形涵蓋的視窗比其他序列短。空白區域
# 表示這裡沒有展示該產品，不能視作零訊號。每列的單位和縱軸也不同，
# 只能比較各列相對自身背景的變化，不能把曲線高度當成共同效果量。
#
# %% tags=["remove-input"]
tr = st.select(channel="HHZ")[0].copy()
tr.detrend("demean")
wave_t = pd.to_datetime(tr.times("timestamp"), unit="s", utc=True)

t0, t1 = EQ_UTC - pd.Timedelta("12h"), EQ_UTC + pd.Timedelta("12h")
gw_w = gw.loc[t0:t1]
mag_w = mag.loc[t0:t1]
cat["time_utc"] = cat.time.dt.tz_localize("UTC")   # 目錄時間本就是 UTC
cat_w = cat[(cat.time_utc >= t0) & (cat.time_utc <= t1)]

fig = make_subplots(rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                    subplot_titles=("地震（ML≥3）", "波形 HWA.HHZ（±30 分鐘）",
                                    "地下水位 TUN（cm）", "地磁 F XCG（nT）"))
fig.add_trace(go.Scatter(x=cat_w.time_utc, y=cat_w.ML, mode="markers",
                         marker=dict(color=QUAKE_COLOR, size=6, opacity=0.7)), 1, 1)
fig.add_trace(go.Scattergl(x=wave_t[::10], y=tr.data[::10], mode="lines",
                           line=dict(color=PALETTE[1], width=0.8)), 2, 1)
fig.add_trace(go.Scattergl(x=gw_w.index, y=gw_w.water_level_cm, mode="lines",
                           line=dict(color=PALETTE[0], width=1.5)), 3, 1)
fig.add_trace(go.Scattergl(x=mag_w.index, y=mag_w.F_calc, mode="lines",
                           line=dict(color=PALETTE[6], width=1.5)), 4, 1)
for row in range(1, 5):
    fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR,
                  row=row, col=1)
fig.update_xaxes(range=[t0, t1], row=4, col=1)
apply_layout(fig, height=800, showlegend=False, hovermode=False,
             title="2024/4/3 花蓮地震：四類觀測同一時間軸（UTC）")
fig

# %% [markdown]
# 目錄中的群集與波形中的強震動，在這組展示資料中最容易辨認；它們分別
# 描述事件的發生與到站震動。地下水和地磁則需更細緻的背景模型與敏感度
# 分析，才能判定有多小的效應仍可能存在。辨認各類觀測對哪些問題有資訊，
# 並不需要替它們排價值高低。
#
# 這個案例以已知地震時刻回看資料，適合研究伴隨效應與設計後續分析，
# 不能單靠它建立震前預報技巧。若想把某個水位或磁場指標加入預報模型，
# 還需要知道非地震時段出現相同指標的頻率，以及加入指標後是否比既有
# 基準提供更多預測資訊。
#
# {doc}`下一章 <08_explore_ideas>`把這兩個問題接起來：從可解釋的物理
# 反應，走到有對照、有不確定性、也能被新資料檢查的統計證據。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **案例背景・免費報告**：臺灣地震科學中心（TEC）（2024），〈[2024 M7.2 花蓮地震](https://tec.earth.sinica.edu.tw/specialEQ/pdf/20240403hwalien.pdf)〉。
#   先看主震位置、區域構造與觀測摘要，再與本章的事件時間軸對照；這是事件報告，閱讀時要留意版本與資料來源。
# - **對照資料・免費網站**：USGS（2024），〈[Hualien earthquake：事件頁 us7000m9g4](https://earthquake.usgs.gov/earthquakes/eventpage/us7000m9g4/executive)〉。
#   比較官方事件參數、震源機制與地動產品；USGS 與本章使用的規模尺度、定位結果及時間表示可能不同，先核對定義再比較數字。
# - **核心案例・開放取用論文**：Zheng et al.（2024），〈[Thrust-dominated unilateral rupture of a blind listric fault associated with the 2024 Hualien earthquake](https://doi.org/10.1038/s41598-024-82971-x)〉，*Scientific Reports*。
#   看作者如何結合地震波與大地測量資料推估破裂過程，特別注意不同觀測對斷層幾何的約束；這是超越本章時間序列對照、進一步建立震源模型的例子。
# - **臺灣比較案例・論文**：Wang et al.（2016），〈[Studies on Aftershocks in Taiwan: A Review](https://doi.org/10.3319/TAO.2016.09.12.01)〉，*Terrestrial, Atmospheric and Oceanic Sciences*；[期刊文章與全文入口](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1498-2016091201t)。
#   把花蓮序列放進臺灣其他地震的背景中，對照餘震分布與觸發機制；這篇早於 2024 年，提供的是比較框架，並非花蓮事件的分析結果。
