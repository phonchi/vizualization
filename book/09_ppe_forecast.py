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
# # 9. 鄰近過去地震：PPE 預報
#
# 一張均勻預報把相同面積分到相同數量。
# 第 4 章的 SUP 採用這項假設。
# 本章改看過去震央的位置。
# 鄰近過去地震模型（Proximity to Past Earthquakes）簡稱 PPE。
# PPE 把較多期望數放在歷史大震附近。
#
# 下圖先展示來源事件，接著拆開計算步驟。
# 第 7 章的 GR 律負責規模比例。
# 本章新增位置權重與時間累積。
# 最後把預報寫成能交換的十欄檔案。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

# %% tags=["remove-input"]
from pathlib import Path
from gdms_toolkit import italy_models
from gdms_toolkit.viz import plot_forecast_map
cat = italy.experiment_catalog()
cells = italy.testing_cells()
windows = italy.forecast_windows()
cache_path = italy.FORECAST_DIR / "PPE_testing.npy"
assert cache_path.exists(), "請先由工具箱維護者準備 PPE 快取"
forecast = italy_models.get_forecast("PPE", "testing")
assert forecast.shape == (len(windows), len(cells), 25)
sources = cat.loc[cat.mb >= italy.SPEC.mT]
first_sources = sources.loc[sources.t_days < windows.iloc[0].t1 - 50]
P = italy_models.PARAMS["PPE"]

# %% [markdown]
# ## 9.1 哪些過去地震提供位置？
#
# 下圖的藍點是收集區內的歷史大震。
# 來源事件（source event）是提供預報貢獻的舊事件。
# 紅點則保留給要檢驗的目標地震。
# 本圖只畫來源，所以使用藍色。
# 收集區外圍的來源，也能影響測試區邊緣。

# %% tags=["remove-input"]
fig = go.Figure(go.Scattergeo(lon=sources.lon, lat=sources.lat, mode="markers",
    marker=dict(color=ACCENT, size=5), customdata=np.c_[sources.mb, sources.year],
    hovertemplate="Mw %{customdata[0]:.1f}<br>%{customdata[1]:.0f} 年<extra></extra>",
    name="全期候選來源"))
poly = italy.collection_polygon()
fig.add_trace(go.Scattergeo(lon=poly.lon, lat=poly.lat, mode="lines",
    line=dict(color="#555555"), name="收集區 S"))
fig.update_geos(fitbounds="locations", visible=False)
apply_layout(fig, title="PPE 的來源位置：全期回顧圖", height=500)
display(Markdown(f"圖 1｜全期共 {len(sources)} 顆名目 Mw ≥ {italy.SPEC.mT:.1f} 候選來源。"
    f"首窗截止前僅 {len(first_sources)} 顆。"
    f"全期圖含後來事件；各窗會另依截止時間篩選。"))
fig

# %% [markdown]
# 全期地圖用來認識資料覆蓋範圍。
# 第一張預報只使用當時已可取得的來源。
# 後來發生的事件，等符合截止條件才加入。
# 把全期地圖一次餵給首窗，會用到未來資訊。
# 每窗的來源清單因此是預報的一部分。
#
# ## 9.2 一個震央如何變成一片權重？
#
# $$k_d(r)=\frac{1}{\pi(d^2+r^2)}.$$
#
# $r$ 是離來源震央的公里距離。
# 核（kernel）是把一個位置展開成鄰近權重的函數。
# 帶寬（bandwidth）$d$ 控制權重鋪開的尺度。
# $d$ 較小時，震央附近的峰較尖。
# $d$ 較大時，遠近差別較平緩。
#
# 下圖是穿過震央的一條剖面。
# 它呈現權重隨距離的變化。
# 剖面（profile）只沿一條線讀取二維函數。
# 因此，線下的面積與地圖總權重不同。
# 地圖積分還要計入周圍各圈的面積。

# %% tags=["remove-input"]
radius = np.linspace(0, 160, 321)
fig = go.Figure()
for d, color in [(15., "#4a3aa7"), (P["d"], ACCENT), (60., "#1baf7a")]:
    fig.add_trace(go.Scatter(x=radius, y=1/(np.pi*(d*d+radius*radius)),
        name=f"d＝{d:.0f} km", line=dict(color=color)))
apply_layout(fig, title="單一來源的空間權重剖面", xaxis_title="離震央距離（km）",
    yaxis_title="核權重（km⁻²）")
display(Markdown(f"圖 2｜公式示範。藍線使用論文 d＝{P['d']:.0f} km。"
    f"其他 {2} 條線只展示帶寬效果。"))
fig

# %% tags=["remove-input"]
show_diagram("d09_kernel_smoothing", caption=f"核疊加示意｜PPE 的固定帶寬為 {P['d']:.0f} km。")

# %% [markdown]
# ## 9.3 先分清權重與機率密度
#
# $$\int_0^L k_d(r)\,2\pi r\,dr
# =\ln\left(1+\frac{L^2}{d^2}\right).$$
#
# 半徑越大，這個積分仍會繼續增加。
# 因此，此核在整個平面上沒有單位總面積。
# 它在 PPE 中提供的是空間權重。
# 正規化（normalization）是除以總量，使積分成為一。
# 要得到機率密度，必須明說正規化的區域。
#
# $$h_0(x,y)=s+\sum_{i\in I(u)}a(m_i-m_T)k_d(r_i).$$
#
# $I(u)$ 是發報時刻 $u$ 可用的來源集合。
# $a$ 在本章是 PPE 的尺度參數。
# 它與第 7 章的 GR 截距用途不同。
# $s$ 是各處都有的微小基底權重。
# 它讓遠離已知震央的位置仍分到期望數。
#
# $m_i-m_T$ 是來源規模超過目標門檻的量。
# 較大來源在這個公式中有較大的權重。
# 恰好等於門檻的來源，這一項為零。
# 所以，來源顆數與有效權重總和是兩件事。
# 讀者應把來源清單與加權公式一起看。
#
# $$H_j=\int_{R_j}h_0(x,y)\,dx\,dy,\qquad
# q_j=\frac{H_j}{\sum_l H_l}.$$
#
# $H_j$ 是第 $j$ 格的空間總權重。
# $q_j$ 是以測試區為條件的空間機率。
# 它們的總和分別是權重總量與一。
# 工具箱保留 $H_j$，再乘時間與規模項。
# 本頁用 $q_j$ 解釋比例，不改動快取數值。
#
# ## 9.4 把面積、時間與規模乘起來
#
# 一個方格內的權重可能隨位置改變。
# 格內積分近似（cell integration approximation）用小塊面積加總。
# 工具箱在每格放置等距中點。
# 每個中點的權重乘上代表面積。
# 全部相加，就近似 $H_j$。
#
# $$H_j\approx\sum_{v=1}^{V}h_0(x_v,y_v)\Delta A_v.$$
#
# 中點足夠密時，平滑函數較容易估準。
# 若核很尖，格心一個值可能錯過高峰。
# 這就是格內取樣與空間解析度的關係。
# 本章讀取既有快取，保留共同積分設定。
# 公式用來理解每個陣列元素的來源。
#
# $$G_k=e^{-\beta(m_k-m_T)}-e^{-\beta(m_{k+1}-m_T)}.$$
#
# $G_k$ 是第 7 章密度在一個規模箱的面積。
# 本站各模型採用相同的 $b=1.084$。
# 規模箱上端有限，因此所有箱的總和略小於一。
# 它表示保留範圍內的期望數比例。
# 省略的較大規模尾端仍屬於原始 GR 模型。
#
# $$\Lambda_{jk}=\ln\left(\frac{t_2}{t_1}\right)H_jG_k.$$
#
# $t_1,t_2$ 是從 1960 年起算的天數。
# 對應時間率為 $1/t$，積分得到對數比。
# 這個時間尺度從目錄起點算起。
# 它與某顆地震剛發生後的餘震時間不同。
# 第 10 章會另外引入主震後時間。
#
# ## 9.5 參數固定，歷史仍會增加
#
# 首窗的 $H_j$ 使用發報前的來源清單。
# 準時間相依（quasi time-dependent）指歷史新增時更新位置權重。
# PPE 的參數固定，來源卻能逐窗增加。
# 窗內沿用發報時刻的來源集合。
# 新增事件在後續發報時才進入模型。
#
# 本站 PPE 快取也保留 50 天延遲。
# 這與論文附錄的來源截止式一致。
# PPE 每窗取發報前 50 天以前的來源。
# 這個截止條件與 EEPAS 相同。
# 發報時刻與來源截止時刻都應列入紀錄。
#
# 首窗圖的顏色表示各格期望數。
# 所有規模箱先在格內加總。
# 同一顏色代表同一圖例刻度。
# 若兩張圖的色階範圍不同，應讀圖例再比較。

# %% tags=["remove-input"]
first = forecast[0].sum(axis=1)
fig = plot_forecast_map(first, title="PPE：測試期首窗預報")
apply_layout(fig, height=520)
display(Markdown(f"圖 4｜發報日 {windows.iloc[0].start:%Y-%m-%d}。"
    f"窗長 {windows.iloc[0].t2-windows.iloc[0].t1:.2f} 天。"
    f"全區期望數 {first.sum():.3f}。資料截止另提前 {50} 天。"))
fig

# %% [markdown]
# ## 9.6 十欄格式與重新分格
#
# 一列十欄檔案描述一個空間格與規模箱。
# CSEP 是地震預報可比較性合作研究。
# 英文為 Collaboratory for the Study of Earthquake Predictability。
# 十欄格式（ten-column format）列出邊界、期望數與啟用旗標。
# 下表將同一列拆成兩段，方便閱讀。
#
# 177 個原始方格使用投影公里座標。
# 把它們改成 0.1 度格，需要重新分格。
# 重新分格（regridding）把原格的數量分到新格。
# 此處按重疊面積比例配置。
# 它採用了原格內均勻分布的近似。
#
# $$\Lambda'_{gk}=\sum_j W_{gj}\Lambda_{jk}.$$
#
# $W_{gj}$ 是原格 $j$ 分給新格 $g$ 的比例。
# 若每一欄的和為一，總期望數就守恆。
# 單純改寫經緯度外接框，不會完成這項配置。
# 因此，本頁呼叫正式重新分格的匯出函式。
# 匯出後仍須核對總數，不能只看副檔名。

# %% tags=["remove-input"]
total = forecast.sum(axis=0)
# 僅匯出快取預報；不重新計算 PPE。
root = Path(italy.__file__).resolve().parents[1]
path = root / "book/_static/downloads/09_ppe_testing_01deg.tsv"
italy.write_csep_10col_regridded(total, path)
exported = pd.read_csv(path, sep="\t")
assert exported.shape[1] == 10
assert np.isfinite(exported.RATE).all() and (exported.RATE >= 0).all()
relative_difference = (exported.RATE.sum()/total.sum()-1)*100
labels = ["經度下界", "經度上界", "緯度下界", "緯度上界", "深度下界 km",
          "深度上界 km", "規模下界", "規模上界", "期望數", "啟用旗標"]
row = exported.iloc[0]
for start in [0, 5]:
    display(HTML(pd.DataFrame({"欄位": exported.columns[start:start+5],
        "意義": labels[start:start+5], "首列值": row.iloc[start:start+5].to_numpy()})
        .to_html(index=False)))
display(Markdown(f"匯出 {len(exported):,} 列。原始總數 {total.sum():.6f}，"
    f"匯出總數 {exported.RATE.sum():.6f}，相差 {relative_difference:+.3f}%。"
    "面積權重逐來源格正規化；總期望數守恆，邊界仍是取樣近似。"))

# %% [markdown]
# 上面的差額是實際匯出檢查結果。
# 面積取樣會帶來近似誤差。
# 每個來源格的面積權重會再正規化。
# 欄和為一，讓總期望數完整分配。
# 總量守恆不表示邊界位置已經精確。
# 後續評分仍使用原始 177 格快取。
#
# [下載十欄示範檔](_static/downloads/09_ppe_testing_01deg.tsv)。
# 檔案包含欄名，數值列依十欄順序排列。
# 檢驗平台的讀取器還可能要求額外資訊。
# 時間窗與模型版本應隨檔案一起交付。
# 重新分格的詳細設定可接著查附錄 E。
#
# ## 9.7 把整段預報與目標事件疊在一起
#
# 下圖加總測試期的全部預報窗。
# 藍色格表示期望數，紅點表示實際目標地震。
# 紅點是事後加入的觀測圖層。
# 預報本身仍遵守每窗資料截止。
# 位置看來接近，只是提出檢驗問題的起點。

# %% tags=["remove-input"]
targets = italy.target_events(cat, "testing")
fig = plot_forecast_map(total.sum(axis=1), targets=targets, title="PPE：測試期合計與目標地震")
apply_layout(fig, height=520)
display(Markdown(f"圖 5｜{len(windows)} 窗合計期望 {total.sum():.2f} 顆。"
    f"事後疊上 {len(targets)} 顆目標地震。"
    "只疊圖不評分，評分在第 16–18 章。"))
fig

# %% [markdown]
# 紅點與深色格的相對位置，還受格子大小影響。
# 較大的格能包含更多不同位置。
# 第 16–18 章會用固定規格比較數量與位置。
# 本章先保留預報，讓檢驗面對明確的數值。
# 讀者現在已能把一格的期望數拆成三個因素。
#
# 一列預報資料也可以當成計算清單。
# 先固定時間窗，再選一個空間格。
# 接著指出規模箱的左右界。
# 這三項都確定後，期望數才有清楚的對象。
# 若把時間窗加長，通常會累積更多期望數。
# 若只改規模箱，改變的則是規模分配比例。
#
# 首窗與十年合計也有不同的用途。
# 首窗讓讀者檢查發報當時能使用的歷史。
# 合計圖則整理多次滾動更新的結果。
# 兩者的單位都是期望數，涵蓋期間卻不同。
# 比較格色以前，先讀圖說中的期間與總數。
# 這個閱讀順序可避免把累積量誤當成瞬時率。
#
# 固定參數使這份實驗可以逐窗重演。
# 每次只需要核對新的資料截止與輸入事件。
# 若發報規則改成更短的延遲，來源清單會變動。
# 若帶寬改變，空間分配也會跟著改變。
# 這些改動都應留下新的模型版本。
# 本章保留共同快取，讓後續比較有一致的起點。
#
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("PPE 來源", "S 內名目 Mw ≥ 5.0；每窗資料截止提前 50 天"),
    ("參數", f"a＝{P['a']}；d＝{P['d']:.0f} km；s＝{P['s']:.1e}；b＝{italy.SPEC.b_value}"),
    ("預報陣列", f"{forecast.shape}；每元素為窗 × 格 × 規模箱期望數"),
    ("交換格式", f"0.1° 十欄示範；總數差 {relative_difference:+.3f}%，待工具箱修正"),
    ("檢驗", "待填（第 16 章）"),
])))

# %% [markdown]
# PPE 已使用歷史位置。
# 接著看大震後幾天內的活動如何改變。
# 下一章是 {doc}`地震會互相引發 <10_clustering_laws>`。
#
# ## 參考資料與延伸閱讀
#
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://www.earth-prints.org/handle/2122/17084)。免費機構全文。第 2 節交代 HORUS 資料。表 3 列出固定參數。
# - Rhoades、Rastin 與 Christophersen（2022），[EEPAS 二十年回顧](https://doi.org/10.3390/geosciences12090349)。免費全文。先讀 PPE 的空間角色，再比較其他核。
# - [公開義大利模型程式](https://github.com/phonchi/EEPAS)。可追查參數與資料檔。程式實作需與論文公式一起核對。
