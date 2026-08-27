# %% [markdown]
# # 3. 地下水位觀測
#
# ## 3.1 為什麼地下水與地震有關？
#
# 地下水位是**孔隙壓力**的直接量測。地殼應力改變時，含水層被壓縮或
# 舒張，水位隨之升降。因此理論上：
#
# - **同震階變（co-seismic step）**：地震瞬間應力重新分布，水位在幾秒到
#   幾分鐘內階梯式升降——這是最確定、最容易觀察到的訊號。
# - **震前異常（?）**：若應力在震前緩慢累積，水位「可能」出現異常趨勢。
#   這是有爭議的開放問題，也是氣象署佈設地震前兆觀測井的動機。
# - 此外水位還忠實記錄**固體潮**（月球與太陽的引潮力壓縮含水層）與
#   **氣壓負載**——這些「已知訊號」是檢驗儀器靈敏度的天然標準源，
#   也是尋找異常前必須先扣除的背景。
#
# ## 3.2 認識觀測井
#
# GDMS 釋出六口即時觀測井，1 秒取樣，同時記錄水位、氣壓與水溫：

# %% tags=["remove-input"]
import plotly.io as pio
pio.renderers.default = "notebook_connected"

# %%
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
# ## 3.3 讀取資料
#
# 使用第 2 章下載的 2024 年 3–4 月資料（`edu-gw-hualien2024.tgz`）。
# 原始資料是 1 秒取樣（每天 86,400 筆），教學上先重取樣成 1 分鐘平均。
# 我們選**壯圍（TUN，宜蘭）**——六口井中資料完整、又離 2024 花蓮地震
# 震央最近的一口。

# %%
tgz = CACHE_DIR / "edu-gw-hualien2024.tgz"
gw = gt.read_groundwater(tgz, "TUN", resample="1min")
gw.info()
gw.head(3)

# %% [markdown]
# ## 3.4 兩個月的水位長什麼樣？

# %%
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
# 讀圖重點：
#
# - 水位有明顯的週期性振盪，這是固體潮與氣壓一起造成的，不是雜訊。
# - 紅色虛線標的是 2024/4/3 花蓮 M7.2 地震。地震會不會在水位上留下痕跡？
#   第 7 章會把這一刻放大到 1 秒解析度來看，答案可能跟你以為的不一樣。
# - 水位與氣壓大致反相：氣壓升、水位降。這個反應的強弱叫「氣壓效率」，
#   是含水層的一個重要性質，3.8 節會實際算出來。
#
# ## 3.5 放大看固體潮

# %%
week = gw.loc["2024-03-10":"2024-03-17"]
fig = px.line(week, y="water_level_cm",
              labels={"Time": "時間（UTC）", "water_level_cm": "水位（cm）"},
              color_discrete_sequence=[ACCENT])
fig.update_traces(line_width=1.5)
apply_layout(fig, title="一週的水位：一天兩次的固體潮清晰可見", showlegend=False)
fig

# %% [markdown]
# 一天兩個波峰，主要是半日潮的 M2、S2 分量。潮汐訊號的振幅與相位反映
# 含水層的彈性性質，有研究主張這些性質在大地震前後會改變，這是地震前兆
# 研究裡一條重要的線索（第 8 章會談方法）。
#
# ## 3.6 頻譜：把週期成分拆開看

# %%
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
# 半日潮（約 12 小時）與全日潮（約 24 小時）的尖峰清楚跳出來。一口好的
# 觀測井，本身就是一台天然的應變儀，能量測到月球與太陽引潮力造成的、
# 只有幾公分的地殼形變。
#
# ## 3.7 不同的井，反應差很多
#
# 同一段時間、同樣的引潮力，不同的井卻有不同的潮汐振幅。把壯圍（TUN）
# 和赤山（CHI，屏東）放在一起比：

# %%
def m2_amplitude(df):
    """用 Welch 功率譜在 M2 頻帶（約 12.4 小時）估潮汐振幅。"""
    x = df.water_level_cm.interpolate().dropna().to_numpy()
    f, p = signal.welch(x - x.mean(), fs=1 / 60, nperseg=2 ** 14)
    per = 1 / f[1:] / 3600
    band = (per > 11.5) & (per < 13)
    return float(np.sqrt(p[1:][band].max()))

chi = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "CHI",
                          resample="1min")
for name, df in [("壯圍 TUN", gw), ("赤山 CHI", chi)]:
    print(f"{name}：M2 潮汐振幅指標 = {m2_amplitude(df):.1f}")

# %% [markdown]
# 兩口井的潮汐振幅差了將近一倍。這不是儀器好壞的問題，而是含水層本身
# 的差異。承壓（受壓）含水層被上覆不透水層封住，地殼一被引潮力壓縮，
# 壓力幾乎立刻反映到水位，潮汐訊號就強；非承壓（自由液面）含水層則像
# 開口的杯子，形變容易從液面洩掉，潮汐反應相對弱。看一口井的潮汐有多
# 清楚，某種程度就能反推它「封」得多好。
#
# ## 3.8 把氣壓的影響算出來：氣壓效率
#
# 前面看到水位和氣壓反相。到底多強？直接對「水位」和「氣壓」做迴歸會被
# 共同的潮汐與長期趨勢污染，所以標準做法是先各取一階差分（看每分鐘的
# 變化量），再迴歸。斜率就是氣壓效率。

# %%
def barometric_efficiency(df):
    d = df[["water_level_cm", "pressure_hPa"]].dropna().diff().dropna()
    slope = np.polyfit(d.pressure_hPa, d.water_level_cm, 1)[0]
    return float(slope)

for name, df in [("壯圍 TUN", gw), ("赤山 CHI", chi)]:
    be = barometric_efficiency(df)
    print(f"{name}：氣壓效率 = {be:+.3f} cm/hPa")

# %% [markdown]
# 氣壓效率是水文地質學家描述含水層的基本參數之一，也是尋找地震相關
# 異常前一定要先扣掉的背景。氣壓每天都在變，如果不先把它對水位的貢獻算
# 清楚、扣乾淨，很容易把一場鋒面過境誤認成「地殼異常」。
#
# ## 3.9 缺測不代表沒資訊
#
# 花蓮（HWA）井的水位在這段期間全是缺測值 9999，讀進來後整欄是 NaN。
# 但同一個檔案裡的水溫欄位還在動：

# %%
hwa = gt.read_groundwater(CACHE_DIR / "edu-gw-hualien2024.tgz", "HWA",
                          resample="10min")
print("水位有效點數：", hwa.water_level_cm.notna().sum())
print("上層水溫有效點數：", hwa.temp_upper_C.notna().sum())
fig = px.line(hwa.dropna(subset=["temp_upper_C"]), y="temp_upper_C",
              labels={"Time": "時間（UTC）", "temp_upper_C": "上層水溫（°C）"},
              color_discrete_sequence=[PALETTE[1]])
fig.update_traces(line_width=1.5)
apply_layout(fig, title="花蓮（HWA）井：水位計故障，但水溫仍持續記錄",
             showlegend=False)
fig

# %% [markdown]
# 這是實務上很常見的情況：一口井的某個感測器壞了，其他感測器還好好的。
# 讀資料時把缺測值轉成 NaN、而不是直接丟掉整口井，就能保住還能用的部分。
# 判斷資料能不能用，永遠要看清楚是哪個欄位、哪段時間缺，而不是看到
# 「有一口井」就當它整口都可信。
