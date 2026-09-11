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
# # 26. 臺灣展望：把觀測、預報與驗證接起來
#
# 上一章把物理觀測轉成可檢查的問題，也重新連回前一部的機率預報。現在
# 回到臺灣：我們已有什麼樣的資料與已發表案例，它們支持哪些判斷，接下來
# 又需要累積什麼證據？
#
# 這個問題不能只用測站數回答。臺灣地震頻繁、陸上觀測密集，但大規模目標
# 事件仍然少；外海、早期年代與大震剛過的時段，資料品質也不同。許多小震
# 不會自動變成同樣多次獨立的大震預報實驗。前面介紹的完整度、相依性與
# 統計功效，到了這裡都是研究設計的實際限制。
#
# 本章依序看長期目錄、空間完整度與大埔預報案例，再用花蓮事件討論空間
# 基準的作用。最後將觀測與模型放進同一個驗證流程。文獻描述的成果以其
# 研究年代與資料範圍為準，不把單一研究當作全臺服務現況的盤點。
#
# ## 26.1 長期目錄同時記錄了兩種歷史
#
# 下圖使用本書儲存的臺灣長期目錄，上方畫較大事件，下方畫全部收錄事件的
# 年數量。圖中規模欄位沿用資料產品的標示；跨年代的原始量測定義仍需
# 核對，不能因欄名相同就假定所有事件都已在同一尺度。
#
# 這份資料有自己的起始年代與收錄門檻，因此其中的事件總數不應與官方
# 其他版本目錄直接相減。判斷是否漏了資料之前，先對齊年代、區域、深度
# 及規模條件。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

cat = load_taiwan_catalog()
yearly = cat.set_index("time").resample("YE").size()

# 兩個標竿事件直接由目錄取出（目錄時間為 UTC，0403 主震落在 04-02 UTC）
def biggest(t0, t1):
    win = cat[(cat.time >= t0) & (cat.time < t1)]
    return cat.loc[win.ML.idxmax()]

chichi = biggest("1999-09-20", "1999-09-22")
hualien = biggest("2024-04-02", "2024-04-04")

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.58, 0.42],
                    vertical_spacing=0.05)
big = cat[cat.ML >= 5.0]
fig.add_trace(go.Scattergl(x=big.time, y=big.ML, mode="markers", name="ML ≥ 5",
                           marker=dict(size=4, color=ACCENT, opacity=0.5)),
              row=1, col=1)
fig.add_trace(go.Bar(x=yearly.index, y=yearly.values, name="年事件數（全目錄）",
                     marker_color=PALETTE[2]), row=2, col=1)
for ev, label in [(chichi, "集集"), (hualien, "0403 花蓮")]:
    for r in (1, 2):
        fig.add_vline(x=ev.time, line_dash="dot", line_color=QUAKE_COLOR,
                      row=r, col=1)
    fig.add_annotation(x=ev.time, y=7.45, text=f"{label} ML {ev.ML:.1f}",
                       showarrow=False, font=dict(color=QUAKE_COLOR),
                       row=1, col=1)
fig.add_annotation(x="1994-06-01", y=yearly.max() * 0.92, xanchor="right",
                   text="1990 年代測網與記錄方式變動", showarrow=False,
                   font=dict(size=11), row=2, col=1)
fig.update_yaxes(title_text="規模 ML", range=[4.8, 7.7], row=1, col=1)
fig.update_yaxes(title_text="年事件數", row=2, col=1)
apply_layout(fig,
             title=f"台灣長期目錄總覽（1973–2025，ML ≥ 2 共 {len(cat):,} 筆；"
                   f"其中 ML ≥ 5 共 {len(big):,} 筆）",
             height=580, hovermode="x")
fig

# %% [markdown]
# 1999 與 2024 年的群集是地震活動的歷史；跨年代的記錄密度變化還包含
# 儀器、測網與目錄編製的歷史。某一年度事件增加，可能同時有多個原因，
# 不能由總數圖把增幅全部分配給其中一項。
#
# Chang et al.（2016）整理臺灣歷史目錄時，面對的正是這個問題。不同年代
# 使用不同規模定義，部分事件還受尺度飽和影響。作者結合較可靠的寬頻
# 規模與經驗轉換，並明確訂定資料來源優先順序。這個案例說明，光記住一條
# 轉換式還不夠：可比較性要靠資料處理與證據建立，不能由欄位名稱保證。
#
# 線性轉換本身有散佈與適用範圍，已失去解析能力的飽和區也不能靠代數
# 反解恢復。若把一份歷史研究的轉換直接套到新的即時目錄，應先檢查資料
# 來源和估計範圍是否一致。完整的規模討論回看
# {doc}`規模與完整度 <11_catalog_completeness_b>`。
#
# 觀測篇的各類產品也在這裡各有位置：波形與矩張量幫助限制震源，強震紀錄
# 連到場址地動，GNSS 和其他連續資料提供變形及環境背景。它們互相補充，
# 但不是每個事件、每個位置都有同樣完整的一套資料。
#
# %% [markdown]
# ## 26.2 同一張臺灣地圖，觀測能力並不均勻
#
# Chan & Wu（2013）分析當時的臺灣目錄，指出陸上與外海的完整度有明顯
# 差異。原因之一是測站幾何：小事件在較遠或涵蓋不足的區域較難被穩定
# 記錄。這種觀測差異也會影響統計估計。
#
# 下面以本書目錄計算格網內的最大曲率指標，作為教學比較。不同顏色表示
# 規模直方圖轉折的位置，不是地下應力的直接地圖。每格只在事件數達到
# 程式設定的門檻時估計，因此空白格也不能解讀成零風險。
#
# %% tags=["remove-input"]
DEG, DM = 0.2, 0.1
MIN_N = 50                       # 每格至少要有這麼多事件才估 Mc

lon_e, lat_e = np.arange(119.4, 123.21, DEG), np.arange(21.4, 26.01, DEG)
mag_e = np.arange(1.95, 7.05, DM)
shallow = cat[(cat.time >= "1994") & (cat.depth <= 30)]

# 三維直方圖：一次掃完，不對 35 萬筆做逐格比對
H, _ = np.histogramdd((shallow.longitude.to_numpy(), shallow.latitude.to_numpy(),
                       shallow.ML.to_numpy()), bins=(lon_e, lat_e, mag_e))
n_cell = H.sum(axis=2)
mag_c = np.round(mag_e[:-1] + DM / 2, 2)
mc_grid = np.where(n_cell >= MIN_N, mag_c[np.argmax(H, axis=2)] + 0.2, np.nan)

lon_c, lat_c = lon_e[:-1] + DEG / 2, lat_e[:-1] + DEG / 2
LO, LA = np.meshgrid(lon_c, lat_c, indexing="ij")
on_land = (LO > 120.0) & (LO < 121.9) & (LA > 22.0) & (LA < 25.2)  # 粗略本島方框
mc_land, mc_sea = np.nanmean(mc_grid[on_land]), np.nanmean(mc_grid[~on_land])

fig = go.Figure(go.Heatmap(x=lon_c, y=lat_c, z=mc_grid.T, colorscale="Blues",
                           zmin=2.2, zmax=3.5, colorbar=dict(title="Mc"),
                           hovertemplate="%{x:.1f}°E %{y:.1f}°N<br>"
                                         "Mc=%{z:.1f}<extra></extra>"))
for ev, name in [(chichi, "1999 集集"), (hualien, "2024 0403 花蓮")]:
    fig.add_trace(go.Scatter(x=[ev.longitude], y=[ev.latitude], mode="markers",
                             name=name,
                             marker=dict(symbol="x", size=12, color=QUAKE_COLOR,
                                         line=dict(width=2))))
apply_layout(fig,
             title=f"Mc 空間場（0.2° 網格，1994 年後，深度 ≤ 30 km，"
                   f"每格 ≥ {MIN_N} 筆；本島平均 {mc_land:.2f}、"
                   f"其餘 {mc_sea:.2f}）",
             xaxis_title="經度", yaxis_title="緯度",
             yaxis_scaleanchor="x", hovermode="closest", height=560)
fig

# %% [markdown]
# 這份公開資料在 $M_L\ge2$ 收錄，圖上的估計又加了固定的 0.2 修正。
# 所以即使真實觀測能力可以看到更小事件，本圖也無法由未收錄的資料推回
# 那段能力。資料已經截斷，演算法多算幾次也無法補回未收錄的資訊。
#
# 圖可用來找值得進一步檢查的區域差異，卻不能直接當作已驗證的逐格
# 完整門檻。事件數、時間窗、群集及分箱方式都會影響估計。Chan & Wu
# 使用的目錄、鄰域和年代與本圖不同，兩者可以比較問題與方向，不能把
# 色階數字一對一當成重現結果。
#
# 對預報的影響很具體：若某區小事件較容易漏掉，模型可能低估當地觸發
# 活動，再把無法解釋的部分分配到背景率。選擇門檻時，需要在保留事件數
# 和資料一致性之間取捨；這個取捨應在模型比較之前固定。
#
# %% [markdown]
# ## 26.3 大埔案例：從事件率走到場址地動
#
# Hsieh et al.（2025）研究大埔地震序列，將時空 ETAS 與地動模型連結。
# 流程先以歷史資料估計模型，預報時接入即時目錄，再模擬可能的未來事件。
# 每一份模擬目錄都可以轉成指定場址的地動，最後統計超越某個震度的比例。
#
# 這個案例把前一部的幾個概念接成了實際流程：條件強度描述事件的發生，
# 模擬處理未來群集的隨機性，地動模型將震源條件連到場址；即時資料的
# 完整度則影響整條流程能否正確更新。已有這樣的案例是重要進展，但它
# 仍需跨事件、跨時段的檢驗，不能單靠一場序列證明普遍校準。
#
# 下面保留論文中同一發報時刻對不同未來窗的機率及觀測摘要，起點為
# 2025 年 1 月 20 日 17:00 UTC。長條是「至少一次」，右軸點線是實際
# 次數；兩種量放在一起是為了閱讀資料，不能用雙軸高度比較準不準。
#
# | 未來窗 | $P(M_L\ge5)$ | $P(M_L\ge6)$ | $P(M_L\ge7)$ | 觀測 $M_L\ge5$ 次數 |
# |---|---|---|---|---|
# | 1 天 | 30.5% | 3.4% | 0.1% | 1 |
# | 3 天 | 42.1% | 6.2% | 0.2% | 1 |
# | 7 天 | 58.8% | 8.1% | 0.3% | 7 |
# | 10 天 | 67.8% | 11.2% | 0.3% | 8 |
#
# 第三列依原文 §3.2 與結論採七天；原 Table 2 的該列印作四天，存在
# 文字不一致，來源處理記於{doc}`附錄 G <appendix_g_observations>`。
#
# %% tags=["remove-input"]
win_label = ["1 天", "3 天", "7 天", "10 天"]
p_ge5 = [30.5, 42.1, 58.8, 67.8]      # P(ML ≥ 5)，%
p_ge6 = [3.4, 6.2, 8.1, 11.2]         # P(ML ≥ 6)，%
p_ge7 = [0.1, 0.2, 0.3, 0.3]          # P(ML ≥ 7)，%
obs_ge5 = [1, 1, 7, 8]                # 實際發生的 ML ≥ 5 次數
obs_ge6 = 0

fig = make_subplots(specs=[[{"secondary_y": True}]])
for name, vals, color in [("P(ML ≥ 5)", p_ge5, PALETTE[0]),
                          ("P(ML ≥ 6)", p_ge6, PALETTE[1]),
                          ("P(ML ≥ 7)", p_ge7, PALETTE[3])]:
    fig.add_trace(go.Bar(x=win_label, y=vals, name=name, marker_color=color),
                  secondary_y=False)
fig.add_trace(go.Scatter(x=win_label, y=obs_ge5, mode="lines+markers",
                         name="實際 ML ≥ 5 次數",
                         line=dict(color=QUAKE_COLOR, width=2.5, dash="dot"),
                         marker=dict(size=11, symbol="diamond")),
              secondary_y=True)
fig.update_yaxes(title_text="至少發生一次的機率（%）", range=[0, 80],
                 secondary_y=False)
fig.update_yaxes(title_text="實際發生次數", range=[0, 10], showgrid=False,
                 secondary_y=True)
apply_layout(fig,
             title=f"2025 大埔序列：ETAS 預報機率 vs 實際觀測"
                   f"（10 天內實際 ML ≥ 5 共 {obs_ge5[-1]} 次、"
                   f"ML ≥ 6 共 {obs_ge6} 次）",
             xaxis_title="自震後 17:00 UTC 起算的時窗",
             barmode="group", hovermode="x", height=470)
fig

# %% [markdown]
# 同一發報時刻下，較長視窗的「至少一次」機率不小於較短視窗，是事件集合
# 包含關係的結果。它不是額外四次獨立的預報成功：一、三、七、十天的窗
# 彼此重疊，也共享同一個序列。
#
# 實際出現八次 $M_L\ge5$ 事件，不能直接與 67.8% 比大小，因為後者沒有
# 告訴我們事件數為八的機率。若要檢查計數是否偏低，需要模型的完整計數
# 分布；若檢查「至少一次」，觀測就應轉成是否發生的二元結果，並在更多
# 預報中評估。相反地，某門檻一次也沒發生，不會自動使低機率預報失敗。
#
# 論文也指出即時目錄的完整度與定位品質，限制了部分相對平靜現象的解讀。
# 這與本部的觀測課題完全相連：少了一筆小震，可能同時改變事件率、觸發
# 歸屬和下一次預報。因此應儲存每次發報當時真正可用的目錄版本，而非
# 只儲存事後整理最完整的一份。
#
# %% [markdown]
# ## 26.4 用花蓮看空間基準的價值
#
# 回到本部的花蓮案例。只用主震之前的地震資料，也能建立一張地震活動
# 較集中在哪裡的空間圖。下面沿用{doc}`EEPAS 與 PPE <16_eepas_ppe>`
# 介紹的空間項，將主震之前的目錄轉成率密度，再把主震位置標上去。
#
# 這裡只展示空間摘要，沒有發布完整的時空規模預報，也沒有執行前瞻檢驗。
# 它可以幫助理解基準模型如何利用歷史空間分布，不能視為當時已發出的
# 官方警報。
#
# %% tags=["remove-input"]
D_KM, S_BG, MC_PPE = 15.0, 1e-4, 5.0        # 沿用 16.6 節的 h0 參數
pre = cat[(cat.time < hualien.time) & (cat.ML >= MC_PPE)]

lons, lats = np.arange(119.0, 123.51, 0.1), np.arange(21.0, 26.01, 0.1)
LON, LAT = np.meshgrid(lons, lats)
dens = np.zeros_like(LON)
for lo, la, mi in pre[["longitude", "latitude", "ML"]].to_numpy():
    r2 = ((LON - lo) * 111 * np.cos(np.radians(LAT))) ** 2 \
         + ((LAT - la) * 111) ** 2
    dens += (mi - MC_PPE + 0.1) * (1 / (np.pi * (D_KM ** 2 + r2)) + S_BG)

# 0403 震央落在全島率密度的哪個百分位？
j = int(np.argmin(np.abs(lats - hualien.latitude)))
i = int(np.argmin(np.abs(lons - hualien.longitude)))
pct = (dens < dens[j, i]).mean() * 100

fig = go.Figure(go.Heatmap(x=lons, y=lats, z=np.log10(dens), colorscale="Blues",
                           colorbar=dict(title="log₁₀ 相對率"),
                           hovertemplate="%{x:.1f}°E %{y:.1f}°N<extra></extra>"))
fig.add_trace(go.Scatter(x=[hualien.longitude], y=[hualien.latitude],
                         mode="markers", name=f"0403 花蓮 ML {hualien.ML:.1f}",
                         marker=dict(symbol="x", size=15, color=QUAKE_COLOR,
                                     line=dict(width=3))))
apply_layout(fig,
             title=f"只用 0403 之前的目錄（{len(pre):,} 個 ML ≥ "
                   f"{MC_PPE:.0f} 事件）：震央落在率密度第 {pct:.0f} 百分位",
             xaxis_title="經度", yaxis_title="緯度",
             yaxis_scaleanchor="x", hovermode="closest", height=560)
fig

# %% [markdown]
# 震央若落在較高率密度區，表示這個事件與空間基準的排序相容。
# 但單一落點還不足以評估整張地圖；也要統計模型把多少率放在沒有事件
# 的地方，以及對其他事件是否同樣有效。
#
# 空間模型即使不隨時間變化，仍可相對均勻空間基準提供資訊。它沒有精確
# 指出某一天，並不表示資訊量為零。新模型可能改善空間、時間、規模或
# 總量中的一項，評估時應說清楚改進在哪裡，並在相同目標上比較。
#
# ETAS 則在既有背景上，隨新事件加入而更新觸發貢獻。主震發生後，模型
# 通常會提高周邊短期事件率；有前震時也可能在較大地震前提高機率。
# 因此不能概括成「ETAS 只知道大震之後」。它的限制是：不會替某個
# 尚未發生的大震指定必然的時間與位置，預報效果仍取決於當時歷史與
# 模型設定。
#
# 地下水、地磁或 GNSS 若提供額外指標，應在這些基準之外接受檢驗。
# 同震反應的物理解釋和震前預報改善是兩個研究成果，兩者都值得做，但
# 不能用前者代替後者。
#
# %% [markdown]
# ## 26.5 時間尺度不同，需要的證據也不同
#
# 下圖是一張工具用途的示意圖。預警關心已發生地震的傳播與通訊；短期預報
# 關心未來事件；中長期模型與危害評估則服務不同期間的風險問題。圖中的
# 時間範圍只是教學定位，不能讀成各系統固定適用的邊界，也不是全臺所有
# 服務的現況認證。
#
# %% tags=["remove-input"]
YR = 365.25
tools = [
    ("地震預警 EEW", np.log10(3 / 86400 / YR), np.log10(60 / 86400 / YR),
     "#1baf7a", "已發生地震的傳播與通訊"),
    ("短期預報 ETAS/OAF", np.log10(1 / YR), np.log10(90 / YR),
     "#2a78d6", "事件率與餘震機率：已發表在地案例"),
    ("中期預報", np.log10(0.25), np.log10(20),
     "#e34948", "較長窗口的模型：需要在地驗證"),
    ("長期危害 TEM PSHA", np.log10(10), np.log10(500),
     "#4a3aa7", "地震來源與場址地動的長期評估"),
]
fig = go.Figure()
for name, lo, hi, color, status in tools:
    fig.add_trace(go.Bar(y=[name], x=[hi - lo], base=[lo], orientation="h",
                         marker_color=color, opacity=0.85, name=status,
                         text=status, textposition="inside",
                         insidetextanchor="middle",
                         hovertemplate=f"{name}<br>{status}<extra></extra>"))
fig.add_annotation(x=np.log10(1.5), y=3.62, showarrow=False,
                   text="資料品質與驗證需求貫穿各種工具",
                   font=dict(size=11, color="#555"))
fig.update_xaxes(title_text="時間尺度（年，log₁₀）", tickvals=[-6, -4, -2, 0, 2],
                 ticktext=["秒–分", "小時", "天", "年", "百年"])
apply_layout(fig, title="地震資訊工具的時間尺度示意（非固定適用界線）",
             showlegend=False, hovermode="closest", height=420)
fig

# %% [markdown]
# 短期案例能讓我們檢查目錄延遲、模型更新與場址地動的連結；更長期間的
# 研究則要處理歷史資料一致性、長期率與較稀少的大事件。這些工作不只是
# 把同一條曲線往外延伸，因為影響預報的背景、資料與使用情境都可能改變。
#
# EEPAS 的臺灣在地化工作正在進行中。此處不引用未發表成果，也不由既有
# 案例推定它已經建立在地的前瞻技巧。對任何新的中期方法，需要累積的
# 是有版本紀錄、可與基準比較的預報，及足以限制不確定性的測試資料。
#
# 同樣地，一份文獻不足以判定全臺缺乏某種研究文化或所有技術都已就緒。
# 較能推進工作的問法是：目前這項研究是否儲存了發報紀錄？哪些目標和
# 基準已對齊？結果最受哪個序列影響？這些問題可以由資料與研究設計逐項
# 回答。
#
# %% [markdown]
# ## 26.6 讓臺灣案例累積成可比較的證據
#
# 比較模型之前，先讓它們預報同一件事。區域、深度、規模尺度、目標門檻、
# 期間與目錄版本必須一致；模型需要哪種輸入目錄，可以不同，但用來評分的
# 目標不能為各模型臨時更換。時間獨立模型也不必一律對除叢目錄評分，
# 要看它原本預報的是全部事件還是背景事件。
#
# 網格大小則影響問題的解析度與檢驗功效。格子太細、事件又少時，很難區分
# 相近模型；自適應或多解析度網格可以作為候選設計，但應在測試前決定，
# 不能承諾換一種網格就一定抓得到錯誤模型。幾何聚合時期望數可以相加，
# 事件相依性卻不會因此消失。
#
# 同樣地，將計數改成「有或沒有」會減少多事件格的權重，但不自動解除
# 跨格或跨時窗的叢集依賴。若用 Poisson、負二項或二元計分，應一起說明
# 分布假設及哪些資料特徵可能使它失準。這些差異可回看
# {doc}`一致性檢驗 <17_testing_consistency>`與
# {doc}`模型比較 <18_testing_comparison>`。
#
# 評估不能只看一個總分。逐段檢查不同序列、不同區域與資料品質，
# 可以看出優勢是否只由少數事件帶動。未拒絕模型可能表示它相容，也可能
# 是資料不足以分辨；觀測密集不會免除功效問題。若把測試規則與不確定性
# 一起公開，新的研究才能接續同一組問題，而不是每次重新挑一個有利案例。
#
# %% [markdown]
# ## 26.7 科學產品與使用者之間
#
# 本部從儀器出發，前一部從預報目標出發，兩條線最後都要面對使用者。
# 地震目錄、短期事件機率、場址地動與地震預警是不同產品。即使共用觀測網，
# 它們的更新時機、提前時間與適用決策仍不同，發布時不能只用「地震資訊」
# 一個名稱帶過。
#
# 有用的機率資訊需要連同時間窗、區域、目標門檻與更新時間呈現。相對背景
# 升高很多的機率，絕對值仍可能很小；對是否值得採取行動的判斷，也取決於
# 行動成本與可能損失。這正是{doc}`作業化系統 <22_operational_systems>`
# 和{doc}`比較與決策 <18_testing_comparison>`已介紹的連結。
#
# 觀測研究也需要對應的溝通精度。「在本井、本時段未辨認出明顯反應」與
# 「地下水沒有地震訊號」是不同範圍的主張；「單一序列的預報有合理表現」
# 也不同於「已證實全臺適用」。說清資料範圍和不確定性，讓研究成果可以
# 被後續工作使用，而不是讓每次新的觀測都被迫在成功與失敗之間二選一。
#
# %% [markdown]
# ## 26.8 把這兩條線一起帶走
#
# 本書先從機率目標建立統計工具，再回到臺灣觀測，看資料如何形成、有哪些
# 背景，以及同一事件在不同儀器中留下什麼。模型需要觀測告訴它哪些量可以
# 被估計；觀測需要模型幫助區分正常變動、資料限制與新的現象。
#
# 接下來可以：為一口井建立可跨時段檢查的背景模型，為地磁指標
# 選擇有物理意義的頻帶與對照站，核對一份目錄的完整度，或儲存一次機率
# 預報並與後續資料比較。每一步都先定義問題、保留來源，再讓未參與選擇的
# 資料檢查結果。
#
# 這些工作不要求從第一天就解決下一場大地震的時間與位置。它們要求的是
# 每個結論都比原來更清楚：它解釋了哪種現象，比哪個基準增加了資訊，
# 還有哪些不確定性值得繼續觀測。主文在這裡告一段落；需要回查推導、
# 數值方法或資料處理時，可以進入第三部附錄。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [最近地震](https://www.cwa.gov.tw/V8/C/E/index.html) — 交通部中央氣象署，官方資料頁（持續更新；免費）。從實際地震報告檢視發震時間、震央、深度、規模與震度，對照本章對資料尺度與目錄用途的區分；單筆事件報告與經整理的研究目錄各有不同用途。
# - [An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized Mw magnitudes](https://doi.org/10.1186/s40623-016-0414-4) — Wen-Yen Chang、Kuei-Pao Chen、Yi-Ben Tsai（2016），*Earth, Planets and Space*（免費全文）。本章規模均一化與歷史目錄的主要來源，重點是跨年代、跨測網的規模如何對齊；使用轉換結果時，一併檢視文章頁連結的勘誤。
# - [Fast report: performance of the ETAS model in forecasting aftershock occurrence and site-specific ground-shaking intensity for the 2025 Dapu, Taiwan, earthquake sequence](https://doi.org/10.1007/s44195-025-00097-7) — Ming-Che Hsieh、Min-Hsuan Chang、Yu-Chen Tai、Chun-Te Chen、Ting-Ying Lu（2025），*Terrestrial, Atmospheric and Oceanic Sciences*（免費全文）。沿著本章大埔案例閱讀即時目錄、ETAS 模擬與場址地動的完整流程，並留意單一序列的結果能支持哪些判斷、仍有哪些跨序列驗證待做。
#
# - [Maximum magnitudes in aftershock sequences in Taiwan](https://doi.org/10.1016/j.jseaes.2013.05.006) — Chung-Han Chan、Yih-Min Wu（2013），*Journal of Asian Earth Sciences*；全文可能需訂閱。讀 §2 的完整度與測網分布，理解陸上、外海的觀測差異；其目錄與估計範圍不同於本書格網示範。
# - [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen（2022），*Geosciences*，免費全文。認識 PPE 空間基準與 EEPAS 的關係，延伸本章「空間資訊也有預報價值」的討論；不把其他地區的結果直接套到臺灣。
# - [b-Values Observations in Taiwan: A Review](https://doi.org/10.3319/TAO.2015.04.28.01%28T%29) — Jeen-Hwa Wang、Kou-Cheng Chen、Pei-Ling Leu、Jeng-Hsin Chang（2015），*Terrestrial, Atmospheric and Oceanic Sciences*。回查臺灣測網與規模定義的歷史，對照長期目錄圖中的記錄密度變化；歷史差異不應全部解讀成地震活動改變。
# - [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Mizrahi et al.（2024），*Reviews of Geophysics*，免費全文。閱讀基準比較、前瞻測試與使用者共同設計三部分，理解研究成果如何成為可累積、可溝通的預報證據。
