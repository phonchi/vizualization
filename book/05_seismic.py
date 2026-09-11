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
# # 22. 地動觀測：波形如何成為地震目錄
#
# 地下水和地磁的變化需要仔細辨認來源；地震儀則直接記錄震波造成的地面運動。
# 不過，波形上有一段震動，仍不等於我們已知道震源在哪裡、規模多大。前一部
# 統計模型使用的每個事件，都經過了從訊號到估計量的過程。
#
# 這章以花蓮站的連續紀錄和本書的 2024 春季目錄互相對照。波形保留訊號
# 隨時間變動的細節；目錄將事件整理成時間、位置、深度與規模。把兩者一起讀，
# 可以看清楚目錄提供了什麼便利，又在哪裡受到觀測能力限制。
#
# ## 22.1 一場地震在三個方向留下的紀錄
#
# 以下是花蓮站 HWA 的三分量波形。原始數位振幅以 counts 表示，必須配合
# 儀器響應才能換成速度或加速度等物理量。圖中的高峰能幫我們定位強烈震動
# 發生的時段，不能直接拿它與別種感測器的峰值比較。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

import gdms_toolkit as gt
from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, SEQUENTIAL, apply_layout

st = gt.read_waveform(CACHE_DIR / "edu-wave-hualien2024.mseed")
print(st)

# %% [markdown]
# 垂直分量標為 HHZ，兩個水平分量為 HH1 與 HH2，後者的實際方位應查測站
# metadata。三條曲線不同，可能來自震波的極化、入射方向、場址以及儀器響應。
# 先比較起振時間與持續長度，再比較各分量的形狀；不要把最大的一條直接
# 當成震源朝向的答案。
#
# %% tags=["remove-input"]
_ = st.plot(size=(1000, 500))

# %% [markdown]
# 強震段把背景振幅壓縮得幾乎看不見。這是共用縱軸的視覺效果，不代表
# 主震前沒有微弱訊號。放大時間窗、選擇合適頻帶能幫助辨認到時，但每一種
# 處理也會改變圖形。
#
# ## 22.2 到時差把時間連回空間
#
# 下面的垂直分量已去均值與線性趨勢，再做 0.5–10 Hz 帶通濾波。紅線是
# 目錄中的發震時刻，不是本站的 P 波到時。兩者之間的延遲，主要來自震波
# 從震源傳到測站所需的時間。濾波與取樣細節保留在原始檔及
# {doc}`附錄 G <appendix_g_observations>`。
#
# %% tags=["remove-input"]
tr = st.select(channel="HHZ")[0].copy()
tr.detrend("demean").detrend("linear")
tr.filter("bandpass", freqmin=0.5, freqmax=10)

t = pd.to_datetime(tr.times("timestamp"), unit="s", utc=True)
sub = slice(int(27.5 * 60 * 100), int(31 * 60 * 100))  # 主震前後幾分鐘
fig = go.Figure(go.Scattergl(x=t[sub], y=tr.data[sub], mode="lines",
                             line=dict(color=ACCENT, width=1)))
fig.add_vline(x="2024-04-02 23:58:11", line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="HWA.HHZ 帶通 0.5–10 Hz（紅線＝發震時刻）",
             yaxis_title="未校正振幅（counts）", showlegend=False, hovermode=False)
fig

# %% [markdown]
# P 波和 S 波的傳播速度不同，因此同站的兩個到時差能約束震源距離。
# 但單站不能給出唯一震央方向，地殼速度也不是處處相同。地震定位需要
# 多站到時、速度模型與誤差評估，得到的震央和深度都是估計結果。
#
# 未移除儀器響應的圖仍以 counts 表示。若要與工程上的地動速度、加速度
# 或震度對照，還要做物理校正，不能直接把這個縱軸當成公尺每秒。
#
# ## 22.3 時頻圖：同一訊號的另一個切面
#
# 時頻圖把每個時間窗裡的頻率成分分開顯示。它能幫助辨認短暫寬頻訊號
# 與較持續的窄頻干擾，但時間窗越長，頻率分辨通常較細、時間定位則較粗。
# 讀圖時也要考慮這個取捨。
#
# %% tags=["remove-input"]
_ = tr.spectrogram(log=True, wlen=10, dbscale=True,
                   title="HWA.HHZ spectrogram")

# %% [markdown]
# 強烈震動通常在多個頻帶留下能量，後續短暫訊號則可能對應其他事件。
# 僅靠時頻圖還不能逐一命名地震，需結合多站到時或已建立的目錄。
#
# ## 22.4 從波形走到事件資料
#
# 下面接回本書使用的春季目錄。每一列的時間、位置與規模已由觀測流程估計，
# 因此可以拿來畫規模–時間圖、做頻率統計。取得方式放在原始檔；此處更重要的
# 是看清本次時間範圍、規模門檻及資料版本，避免把這份子目錄當成完整的全臺
# 歷史紀錄。
#
# %% tags=["remove-input"]
cat_csv = CACHE_DIR / "catalog_2024spring.csv"
if cat_csv.exists():
    cat = pd.read_csv(cat_csv, parse_dates=["time"])
else:
    gdms = gt.GDMSSession()
    cat = gt.gdms_earthquake_catalog(gdms, "2024-03-01", "2024-06-30", min_ml=3)
    cat.to_csv(cat_csv, index=False)
print(f"2024/03–06 ML≥3 共 {len(cat)} 筆")
cat.tail(3)

# %% [markdown]
# ## 22.5 規模–時間圖把群集顯示出來
#
# 橫軸為 UTC，縱軸為 $M_L$，每一個點是一筆目錄事件，顏色表示深度。
# 主震附近的密集點列讓我們看到短時間內事件率的變化；高低分布則提供規模
# 分布的線索。兩者是不同問題，不能由「點很多」直接推成「下一震一定更大」。
#
# %% tags=["remove-input"]
fig = px.scatter(cat, x="time", y="ML", color="depth",
                 color_continuous_scale=SEQUENTIAL,
                 labels={"time": "時間", "ML": "規模 ML", "depth": "深度（km）"},
                 hover_data={"latitude": ":.2f", "longitude": ":.2f"})
fig.update_traces(marker=dict(size=5, opacity=0.6))
fig.add_vline(x="2024-04-02 23:58", line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="規模–時間圖（紅線＝花蓮主震；目錄時間為 UTC，"
                       "當地時間為 4/3 07:58）", hovermode="closest")
fig

# %% [markdown]
# 主震後密集的事件，將前一部的叢集模型變成具體觀測。之後活動可能衰減，
# 也可能因較大餘震出現而再次升高。ETAS 讓各事件都能加入觸發貢獻，正是
# 為了描述這種重疊。
#
# 「前震」則是一種在後續較大事件出現後才確定的關係標籤。可以事先計算
# 某事件後出現更大地震的機率，卻不能在看到一個事件時就確定它將成為前震。
#
# ## 22.6 空間分布提供幾何線索
#
# 把目錄點放回地圖，先比較事件密集區、深度和海陸位置。符號大小隨規模
# 增加而放大，並非地震破裂範圍的比例尺。
#
# %% tags=["remove-input"]
fig = px.scatter_map(cat, lat="latitude", lon="longitude",
                     size=2 ** cat.ML / 10, color="depth",
                     color_continuous_scale=SEQUENTIAL,
                     hover_data={"time": True, "ML": True},
                     zoom=6.3, center=dict(lat=23.7, lon=121.2),
                     map_style="carto-positron", height=600)
apply_layout(fig, title="2024/03–06 ML≥3 震央分布（符號隨規模放大、顏色＝深度）",
             hovermode="closest")
fig

# %% [markdown]
# 密集事件帶可提示活動構造與震源區的空間範圍，但震央是震源位置在地表的
# 投影，餘震也可能出現在主要破裂面周圍。要確定斷層幾何及滑移分布，還需
# 重定位、震源機制、波形和地表位移等資料，不能直接把點雲邊界畫成破裂面。
#
# ## 22.7 規模–頻率圖應怎樣讀
#
# 下面重畫前一部的 Gutenberg–Richter 關係。門檻越高，累積事件數越少；
# 如果某段範圍近似直線，斜率描述那個範圍內大小事件的相對比例。
# 圖上的最小平方法直線只是描述性參考，不是本書推薦的正式 $b$ 值估計，
# 因為累積計數彼此相依。正式估計與完整度處理請回看
# {doc}`規模與完整度 <11_catalog_completeness_b>`。
#
# %% tags=["remove-input"]
mags = np.arange(3, cat.ML.max() + 0.1, 0.1)
N = [(cat.ML >= m).sum() for m in mags]
fig = go.Figure(go.Scatter(x=mags, y=N, mode="markers",
                           marker=dict(color=ACCENT, size=7)))
# 以 3.5–5.5 區間做最小二乘擬合
sel = (mags >= 3.5) & (mags <= 5.5)
b, a = np.polyfit(mags[sel], np.log10(np.array(N)[sel]), 1)
fig.add_trace(go.Scatter(x=mags, y=10 ** (a + b * mags), mode="lines",
                         line=dict(color=QUAKE_COLOR, dash="dash"),
                         name=f"擬合 b = {-b:.2f}"))
apply_layout(fig, title="Gutenberg–Richter：累積次數 vs 規模",
             xaxis_title="規模 ML", yaxis_title="N（ML ≥ M）",
             yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 低規模端變平，可能反映收錄門檻與小事件漏測；不能把曲線轉折直接等同
# 已經確認的 $M_c$。大規模端則事件少，單一事件就能明顯改變尾部形狀。
# 本圖把整個春季合併，還混合了主震前後不同活動與偵測條件。
#
# ## 22.8 衰減是平均趨勢，不是每日的時間表
#
# 下面以固定日窗數事件，並畫一條 $1/t$ 參考線。參考線沒有估計最佳 $p$，
# 也沒有包含各次餘震的二次觸發，計數亦未另做親代分類。它的用途是讓我們看見量級如何隨時間下降，
# 以及哪些日子明顯偏離簡單的平均輪廓。
#
# %% tags=["remove-input"]
main = cat.loc[cat.ML.idxmax()]
aft = cat[cat.time >= main.time].copy()
aft["day"] = ((aft.time - main.time).dt.total_seconds() // 86400).astype(int) + 1
daily_n = aft[aft.day.between(1, 30)].groupby("day").size()

fig = go.Figure(go.Scatter(x=daily_n.index, y=daily_n.values, mode="markers",
                           marker=dict(color=ACCENT, size=8)))
fig.add_trace(go.Scatter(x=daily_n.index, y=daily_n.iloc[0] / daily_n.index,
                         mode="lines", line=dict(color=QUAKE_COLOR, dash="dash"),
                         name="~ 1/t 參考線"))
apply_layout(fig, title="主震後每日事件數（ML≥3）：衰減參考",
             xaxis_title="主震後天數", yaxis_title="當日事件數",
             xaxis_type="log", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 事件數不會每天單調下降。較大餘震可以開啟新的群集，早期漏測也可能
# 壓低最初的計數。圖中的第一個日窗從主震時刻開始，計數包含主震本身；
# 若要正式擬合餘震衰減，必須先定義是否排除它，以及空間區域與完整門檻。
#
# ## 22.9 同一份目錄的兩個 $b$ 值
#
# 最後比較主震前與主震後三十天的規模資料。下列估計使用 $M_L\ge3.5$，
# 並對 0.1 規模刻度作修正。這是一個敏感度示範：顯示的樣本數是各時段
# 全部收錄事件數，真正進入估計的，是其中超過門檻的事件。兩個時段的樣本
# 大小與偵測條件並不相同。
#
# %% tags=["remove-input"]
def b_value(magnitudes, mc):
    """Aki (1965) 最大概似法估 b 值。"""
    m = magnitudes[magnitudes >= mc]
    return np.log10(np.e) / (m.mean() - (mc - 0.05))

pre = cat[cat.time < main.time].ML
post = aft[aft.day.between(1, 30)].ML
print(f"主震前（3/1–4/2）：{len(pre)} 筆，b = {b_value(pre, 3.5):.2f}")
print(f"餘震期（30 天）：{len(post)} 筆，b = {b_value(post, 3.5):.2f}")

# %% [markdown]
# 即使兩個 $b$ 值不同，仍有抽樣、門檻、尺度與空間組成等替代解釋。
# 不確定性估計也必須對應資料的相依結構：若要評估整個目錄流程，不能
# 無條件把群集事件當作獨立樣本重抽。可比較完整度較一致的範圍，或用
# 合適的分塊／模型模擬來檢查差異是否穩定。
#
# 這些圖呈現了波形如何整理成點過程資料。波形與目錄讓我們認識地震的
# 震動與發生分布，卻不完整描述地表最後移到了哪裡。
# {doc}`下一章的 GNSS <06_gnss>`補上這個幾何觀點。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **操作入門・免費教材**：ObsPy 開發團隊，〈[Tutorial](https://docs.obspy.org/tutorial/)〉。
#   依序看 Reading Seismograms、Filtering Seismograms 與 Plotting Spectrograms，把本章的讀檔、濾波及時頻圖流程接起來。
# - **工具論文**：Beyreuther et al.（2010），〈[ObsPy: A Python Toolbox for Seismology](https://doi.org/10.1785/gssrl.81.3.530)〉，*Seismological Research Letters*；[出版學會免費全文](https://www.seismosoc.org/Publications/SRL/SRL_81/srl_81-3_es/)。
#   瞭解 ObsPy 如何把不同格式的地震資料接到共同的處理流程；論文用來認識設計背景，程式語法以目前官方檔案為準。
# - **臺灣規模統計・論文**：Wang et al.（2015），〈[b-Values Observations in Taiwan: A Review](https://doi.org/10.3319/TAO.2015.04.28.01%28T%29)〉，*Terrestrial, Atmospheric and Oceanic Sciences*。
#   延伸本章的規模－次數圖與震前震後比較，留意目錄品質、取樣區間和構造差異如何影響 $b$ 值；不同的估計值本身還不足以證明前兆。
# - **臺灣餘震統計・論文**：Wang et al.（2016），〈[Studies on Aftershocks in Taiwan: A Review](https://doi.org/10.3319/TAO.2016.09.12.01)〉，*Terrestrial, Atmospheric and Oceanic Sciences*；[期刊文章與全文入口](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1498-2016091201t)。
#   對照臺灣不同序列的餘震分布、Omori 衰減與觸發機制，思考本章花蓮個案的結果能推廣到什麼範圍。
