# %% [markdown]
# # 7. 案例分析：2024/4/3 花蓮地震 M7.2
#
# 2024 年 4 月 3 日上午 7 時 58 分（台灣時間；UTC 4/2 23:58），
# 花蓮外海發生規模 7.2 地震——台灣自 921 之後最大的地震。
# 本章把前面學過的**四類資料放到同一個事件上**，完整走一遍
# 「多參數對照分析」的流程。
#
# | 資料 | 測站 | 距震央 |
# |---|---|---|
# | 地震目錄 | 全台 | — |
# | 波形 | HWA 花蓮氣象站 | ~12 km |
# | 地下水位 | TUN 壯圍（宜蘭） | ~95 km |
# | 地磁 | XCG 新城（花蓮） | ~30 km |

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %%
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
# ## 7.1 主震與餘震序列

# %%
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
apply_layout(fig, title="主震後 30 天餘震分布（顏色＝發生時間）", hovermode="closest")
fig

# %% [markdown]
# 餘震帶從震央向北北東延伸約 60 公里——大致就是斷層破裂面在地表的投影。
#
# ## 7.2 強烈的地動：主震波形

# %%
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
# ## 7.3 地下水位：預期會有同震反應，實際呢？
#
# 教科書上，大地震會在地下水位留下兩種痕跡：震波經過時水位像液面地震儀
# 一樣振盪（稱為水震盪，hydroseismogram），以及震後水位停在跟震前不同
# 的高度（同震階變，co-seismic step）。壯圍井距震央約 95 公里，理論上有
# 機會看到。先看主震前後各三天（1 分鐘平均）：

# %%
gw = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "TUN",
                         start="2024-03-30", end="2024-04-06", resample="1min")
fig = go.Figure(go.Scattergl(x=gw.index, y=gw.water_level_cm, mode="lines",
                             line=dict(color=PALETTE[0], width=1.5)))
fig.add_vline(x=str(EQ_UTC), line_dash="dash", line_color=QUAKE_COLOR)
apply_layout(fig, title="壯圍（TUN）水位，主震前後各三天",
             yaxis_title="水位（cm）", showlegend=False)
fig

# %% [markdown]
# 水位就是照著原本的潮汐節奏起伏，紅線那一刻看不出明顯的跳動。放大到
# 1 秒解析度、只看地震前後這 20 分鐘，再確認一次：

# %%
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
# 誠實的結果：**看不到水震盪，也沒有可辨識的階變**。水位在地震前後就是
# 平順地延續原本每分鐘約 0.01 cm 的潮汐趨勢，起伏幅度不到 0.3 cm。
#
# 用數字確認這個判斷，別只靠肉眼。比較「地震前」與「地震後」5 分鐘內
# 1 秒差分的標準差，如果有水震盪，震後這個值應該明顯放大：

# %%
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
# 四口井的比值全都接近 1.0，也就是說地震前後的擾動程度沒差別，沒有水震盪。
# 這不代表理論錯了。距離約 100 公里、規模 7.2，靜態體應變本來就小，換算
# 成水位可能只有毫米量級，剛好埋在這幾口井約 0.04 cm 的量測噪訊底下。要
# 在近場、井的封閉性又好，再加上乾淨的儀器，才比較有機會抓到教科書上那種漂亮
# 的階變。
#
# 這是本課程最想讓你帶走的一課：**「應該會有訊號」和「資料裡真的看得到
# 訊號」是兩回事**。真正的分析要能算出一個門檻，誠實回答「以這口井的噪訊
# 水準，多大的訊號才看得出來」，而不是憑印象宣稱看到了什麼。
#
# ## 7.4 地磁場有反應嗎？

# %%
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
# 一樣看不出與地震明顯相關的變化。每天規律的日變化，遠遠大過任何可能的
# 震磁訊號。這和地下水那節得到的是同一個教訓：震磁效應至今仍有爭論，
# 部分原因就是它（如果存在）實在太小，很難從背景裡分出來。負面結果也是
# 結果，把它老實呈現出來，是科學態度的一部分。
#
# ## 7.5 四合一總覽

# %%
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
# 這張圖就是本課程的縮影，四類資料在同一時間軸上，清楚程度天差地遠：
#
# - 地震目錄：主震後餘震密集出現，訊號明確。
# - 波形：主震的巨大振幅，一眼就看得到。
# - 地下水：在這幾口井、這個距離下，沒有可解析的同震訊號。
# - 地磁：日變化淹沒一切，看不出地震相關的變化。
#
# 換句話說，這場地震留給我們的、真正乾淨明確的觀測，是波形和餘震。地下水
# 和地磁的「前兆」或「同震反應」，在這個案例裡並沒有出現在資料上。這不是
# 失敗的分析，而是一次誠實的分析：它告訴我們，用哪些資料、在什麼條件下，
# 才有機會看到什麼樣的訊號。
#
# 帶著這個經驗回頭想「地震前兆」這個問題，會比一開始務實很多。想更進一步，
# 第 8 章整理了一套從「單一事件的巧合」走向「統計上站得住腳的證據」的方法。
