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
# # 13. ETAS 在義大利：估兩個參數、出一張三個月預報
#
# 一列 HORUS 事件提供時間、位置與規模。
# 第 12 章把這些欄位接進 ETAS。
# 本章先用學習期資料估兩個參數。
# 接著讀取快取中的三個月預報。
#
# 兩張成果使用不同的設定。
# 估計圖使用提高門檻的教學子目錄。
# 預報圖使用 Biondini 的既定參數。
# 本站近似為**背景＋第一代觸發**。
# 本章會逐項說明差別。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy, italy_models
from gdms_toolkit.viz import ACCENT, PALETTE, SEQUENTIAL, QUAKE_COLOR, apply_layout, plot_forecast_map, show_diagram

# %% [markdown]
# ## 13.1 先列出誰參與估計
#
# $M\ge4.0$ 是本次估計示範的門檻。
# HORUS 在早期也較能完整記錄這些事件。
# 學習期取 1990–2011 年。
# 測試期事件全數留到預報圖才使用。
#
# 測試區內的學習事件提供概似的加總項。
# 收集區內更早的事件提供歷史。
# 補充事件（complementary events）是只提供歷史的事件。
# 它們可以位於測試區外或學習期前。
#
# 本次概似中的目標事件（target events），
# 指接受條件強度評估的學習事件。
# 它們採示範門檻 $4.0$。
# 全書預報檢驗的目標門檻仍為 $5.0$。
# 兩種用途必須連同門檻一起讀。
#
# 一個區外事件也可能在區內引發後代。
# 因此，估計保留收集區內的歷史事件。
# 每個加總項只使用該事件之前的歷史。
# 學習期後段事件無法解釋前段事件。

# %% tags=["remove-input"]
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from time import perf_counter
cat = italy.experiment_catalog()
cells = italy.testing_cells()
P = italy_models.PARAMS["ETAS"]
source = cat.loc[(cat.mb >= 4.0) & (cat.period != "testing")].sort_values("t_days").reset_index(drop=True)
obs = source.loc[(source.period == "learning") & source.in_R].copy()
t_start, t_end = italy.year_start_days(1990), italy.year_start_days(2012)
area, duration = cells.area_km2.sum(), t_end-t_start
radius, m_ref = 100.0, 4.0
background = len(obs)/(area*duration)
display(Markdown(f"估計使用 {len(obs)} 顆區內學習事件；歷史來源共 {len(source)} 顆。來源門檻為 M≥{m_ref:.1f}，最遠納入 {radius:.0f} km。"))

# %% [markdown]
# ## 13.2 固定形狀，只改兩個高度
#
# $\nu$ 控制背景大小，$K$ 控制觸發大小。
# 固定形狀估計（fixed-shape estimation）只調整振幅。
# 本例釘住 $b,c,p,D,\gamma,\alpha,q$。
# 七個值皆採 Biondini 表 3。
#
# $B_j$ 是事件位置上的參考背景。
# $A_j$ 是取 $K=1$ 時的歷史觸發貢獻。
# $B_R$ 與 $A_R$ 是學習時空範圍的積分。
# 第 11 章的概似因此縮成
#
# $$
# \ell(\nu,K)=\sum_j\log(\nu B_j+K A_j)
# -\nu B_R-K A_R+\mathrm{const}.
# $$
#
# 規模密度的參數已固定。
# 它在本次最佳化中只增加常數。
# 最大概似估計（maximum likelihood estimation），
# 簡寫 MLE，選出使上式最大的參數。
#
# 本例的 $K$ 乘在 $(\tau+c)^{-p}$ 前面。
# 這與義大利原文的時間核慣例相同。
# 但本例把產能參考規模提高到 $4.0$。
# 空間尺度也相對這個門檻計算。
# 因此，估計值只描述本次簡化模型。
# 它不直接替換標準預報的參數。
#
# 一百公里以外的觸發項設為零。
# 鄰近搜尋樹（k-d tree）先找出距離內的事件。
# 程式再排除未來的來源事件。
# 補償子的空間積分也採相同限距。
# 這讓概似兩項使用同一個模型。

# %% tags=["remove-input"]
fit_started = perf_counter()
sxy = source[["x_km", "y_km"]].to_numpy()
oxy = obs[["x_km", "y_km"]].to_numpy()
st = source.t_days.to_numpy()
sm = source.mb.to_numpy()
d2 = P["D"]**2*np.exp(P["gamma"]*(sm-m_ref))
prod = np.exp(P["alpha"]*(sm-m_ref))
neighbors = cKDTree(sxy).query_ball_point(oxy, radius)
A = np.zeros(len(obs))
for j, (tj, xy, near) in enumerate(zip(obs.t_days, oxy, neighbors)):
    ids = np.asarray(near, dtype=int)
    ids = ids[st[ids] < tj]
    r2 = np.sum((sxy[ids]-xy)**2, axis=1)
    spatial = (P["q"]-1)/(np.pi*d2[ids])*(1+r2/d2[ids])**(-P["q"])
    A[j] = np.sum(prod[ids]*(tj-st[ids]+P["c"])**(-P["p"])*spatial)
# 積分用固定極座標樣點；只保留落在 R 的份量。
n_quad = 1024
u = (np.arange(n_quad)+0.5)/n_quad
angle = np.arange(n_quad)*np.pi*(3-np.sqrt(5))
radial_mass = 1-(1+radius**2/d2)**(1-P["q"])
rad = np.sqrt(d2[:,None]*((1-radial_mass[:,None]*u)**(1/(1-P["q"]))-1))
points = sxy[:,None,:]+rad[:,:,None]*np.stack([np.cos(angle),np.sin(angle)],axis=1)
centers = np.column_stack(((cells.x0+cells.x1)/2,(cells.y0+cells.y1)/2))
nearest = cKDTree(centers).query(points.reshape(-1,2))[1]
offsets = np.abs(points.reshape(-1,2)-centers[nearest])
half = np.column_stack(((cells.x1-cells.x0)/2,(cells.y1-cells.y0)/2))
in_region = np.all(offsets <= half[nearest]+1e-8,axis=1).reshape(len(source),n_quad)
space_mass = radial_mass*in_region.mean(axis=1)
lower = np.maximum(t_start-st,0)
upper = t_end-st
time_mass = ((lower+P["c"])**(1-P["p"])-(upper+P["c"])**(1-P["p"] ))/(P["p"]-1)
A_R = np.sum(prod*time_mass*space_mass)
B_R = background*area*duration
B = np.full(len(obs),background)
def objective(v):
    rate = v[0]*B+v[1]*A
    value = -np.log(rate).sum()+v[0]*B_R+v[1]*A_R
    grad = np.array([B_R-np.sum(B/rate), A_R-np.sum(A/rate)])
    return value, grad
fit = minimize(objective,[0.5,0.02],jac=True,bounds=[(1e-8,1.),(1e-8,None)],method="L-BFGS-B",options={"ftol":1e-12,"gtol":1e-7})
assert fit.success, fit.message
nu_hat, K_hat = fit.x
fit_seconds = perf_counter()-fit_started
assert fit_seconds < 30
nu_grid = np.linspace(max(0.001,nu_hat*0.4),min(1.,nu_hat*1.6),55)
k_grid = np.linspace(max(1e-6,K_hat*0.35),K_hat*1.65,55)
ll = np.array([[-objective([v,k])[0] for k in k_grid] for v in nu_grid])
fig = go.Figure(go.Contour(x=k_grid,y=nu_grid,z=ll+fit.fun,colorscale=SEQUENTIAL,colorbar_title="相對最大值"))
fig.add_trace(go.Scatter(x=[K_hat],y=[nu_hat],mode="markers",marker=dict(color=PALETTE[2],symbol="x",size=12),name="MLE"))
apply_layout(fig,xaxis_title="K（本例門檻與時間核）",yaxis_title="ν",height=420,hovermode="closest")
display(Markdown(f"本次估計 ν={nu_hat:.3f}、K={K_hat:.4f}。每個來源用 {n_quad} 個固定空間樣點積分；這是限距近似的概似面。"))
fig

# %% [markdown]
# 藍色等高線表示不同參數的概似差異。
# 沿著較長方向移動，概似可能變化較小。
# 不可辨識性（non-identifiability）表示不同參數難以區分。
# 本圖可以檢視參數互相補償的方向。
# 是否足以區分，還要看資料量與曲面寬度。
#
# 圖上的最大點只解決已指定的估計問題。
# 固定參數、限距與積分近似都影響位置。
# 空間樣點若太少，邊界份量也會改變。
# 因此，本圖用來理解估計步驟。
# 正式應用還須檢查積分收斂與參數穩定性。
#
# 背景高度與觸發高度都能提高事件位置的率。
# 補償子則扣除整段時空範圍的總率。
# 兩項一起使用，才會限制任意抬高預報。
# 只看事件位置的高率，會漏掉這份代價。

# %% [markdown]
# ## 13.3 一顆事件，可以分給多個來源
#
# 事件 $j$ 的總率可拆成背景與各親代貢獻。
# 隨機除叢（stochastic declustering）用機率分配來源。
# 令 $a_{ij}$ 為來源 $i$ 的單位觸發貢獻。
# 便可寫成
#
# $$
# \phi_j=\frac{\nu B_j}{\nu B_j+K A_j},\qquad
# \rho_{ij}=\frac{K a_{ij}}{\nu B_j+K A_j}.
# $$
#
# $\phi_j$ 是背景機率。
# $\rho_{ij}$ 是由事件 $i$ 觸發的機率。
# 同一事件滿足 $\phi_j+\sum_i\rho_{ij}=1$。
# 這些機率由指定模型計算。
# 它們表達統計歸屬，物理因果仍須另找證據。
#
# <details><summary>把單一主震固定，會留下什麼？</summary>
#
# 一顆指定主震可以獨自提供歷史觸發項。
# Reasenberg–Jones 模型可簡寫為 R–J。
# 它用主震後的衰減率描述餘震活動。
# 固定主震，並略去後續事件的再觸發，
# 就得到與 R–J 相近的單序列形式。
#
# ETAS 允許後續事件再產生觸發項。
# 因此，每次新事件都可能再次改變率。
# 兩種形式的產能係數使用慣例可能不同。
# 比較數值前，須先把時間與規模項對齊。
#
# </details>

# %% [markdown]
# ## 13.4 同一批起點，讓後代繼續繁殖
#
# 兩百條合成目錄可以顯示世代差異。
# 本例設定一顆已知的起始事件。
# 先抽背景事件與它的直接後代。
# 第一代近似在此停止。
# 完整分支版本再讓這批新事件繁殖。
#
# 兩個版本共用同一批初始抽樣。
# 配對模擬（paired simulation）讓共同部分保持一致。
# 計數差異因而來自新增的繁殖步驟。
# 本例只比較時間與規模。
# 參數為合成設定，未擬合義大利資料。
#
# 計算先估計所需的計數陣列。
# 程式每次只保留一條事件清單。
# 模擬持續追蹤後代，直到本次預報窗結束。
# 模擬未把測試期真實事件當作未來輸入。

# %% tags=["remove-input"]
sim_started = perf_counter()
rng = np.random.default_rng(130913)
n_sim, n_bins = 200, 20
memory_bytes = n_sim*n_bins*8
T, b_syn, alpha_syn, k_syn, c_syn, p_syn = 91.31, 1.0, 1.0, 0.20, 0.05, 1.3
beta_syn, m_min, m_max = b_syn*np.log(10), 2.5, 7.5
edges = np.linspace(0,T,n_bins+1)
first_counts = np.zeros((n_sim,n_bins),dtype=np.int64)
full_counts = np.zeros_like(first_counts)
max_events = 0
def mags(size):
    return m_min-np.log1p(-rng.random(size)*(1-np.exp(-beta_syn*(m_max-m_min))))/beta_syn
def children(t,m):
    mass = 1-(1+(T-t)/c_syn)**(1-p_syn)
    number = rng.poisson(k_syn*np.exp(alpha_syn*(m-m_min))*mass)
    waits = c_syn*((1-rng.random(number)*mass)**(1/(1-p_syn))-1)
    return list(zip(t+waits,mags(number)))
for run in range(n_sim):
    nb = rng.poisson(0.08*T)
    events = list(zip(rng.uniform(0,T,nb),mags(nb)))+children(0.,5.)
    first_counts[run] = np.histogram([e[0] for e in events],edges)[0]
    cursor = 0
    while cursor < len(events):
        events.extend(children(*events[cursor]))
        cursor += 1
        assert len(events) < 50000, "事件清單超過記憶體保護上限"
    max_events = max(max_events,len(events))
    full_counts[run] = np.histogram([e[0] for e in events],edges)[0]
sim_seconds = perf_counter()-sim_started
assert sim_seconds < 20
fig = go.Figure()
for values,label,color in [(first_counts,"背景＋已知事件第一代",ACCENT),(full_counts,"允許未來事件再繁殖",PALETTE[2])]:
    fig.add_trace(go.Histogram(x=values.sum(axis=1),name=label,opacity=0.65,marker_color=color,xbins=dict(start=0,end=100,size=5)))
apply_layout(fig,barmode="overlay",xaxis_title="一窗模擬事件數（M≥2.5）",yaxis_title="目錄份數",height=400,legend=dict(orientation="h",y=-0.25))
display(Markdown(f"合成目錄各 {n_sim} 條。平均數：第一代 {first_counts.sum(axis=1).mean():.2f}，再繁殖 {full_counts.sum(axis=1).mean():.2f}。"))
fig

# %% [markdown]
# 綠色分布多出的事件來自未來的後代。
# 這張圖說明第一代近似省略了什麼。
# 它的計數門檻與義大利目標門檻不同。
# 兩張預報的差異大小須用各自設定計算。
#
# 每次正式發報只使用當時已知的歷史。
# 下一個窗可以納入剛發生的真實事件。
# 這種逐窗更新仍遵守資料截止時間。
# 十年合計圖因此是一串預報的總和。
# 它也保留了每次更新前的不同資訊量。

# %% [markdown]
# ## 13.5 讀取義大利預報快取
#
# 下圖選取接近 2016 年第三季的標準窗。
# 標準窗長依規格固定，並非逐季日曆切割。
# 圖說會列出快取採用的實際起迄。
# 紅點是該窗結束後才能知道的目標事件。
#
# 每次發報時，模型先固定已知歷史。
# 快取使用名目輸入門檻 $2.5$。
# 參數使用 Biondini 表 3 的值。
# 本章剛估出的兩個值只用於估計示範。
#
# 藍色表示方格內各規模箱的期望數總和。
# 深色區反映已知事件與背景的共同貢獻。
# 本章只疊圖不評分。
# 正式評分在第 16–18 章。

# %% tags=["remove-input"]
assert (italy.FORECAST_DIR/"ETAS_testing.npy").exists(), "請先提供預報快取"
forecast = italy_models.get_forecast("ETAS","testing")
windows = italy.forecast_windows()
idx = int(np.argmin(np.abs(windows.t1-italy.decimal_days(pd.Timestamp("2016-07-01")))))
w = windows.iloc[idx]
targets = italy.target_events(cat,"testing")
in_window = targets.loc[(targets.t_days>=w.t1)&(targets.t_days<w.t2)]
fig = plot_forecast_map(forecast[idx].sum(axis=1),targets=in_window,title="ETAS：背景＋第一代，接近 2016 Q3 的窗")
apply_layout(fig)
display(Markdown(f"第 {idx+1} 窗：{w.start:%Y-%m-%d %H:%M} 至 {w.end:%Y-%m-%d %H:%M} UTC。期望 {forecast[idx].sum():.3f} 顆，疊上 {len(in_window)} 顆目標事件。只疊圖不評分，評分在第 16–18 章。"))
fig

# %% tags=["remove-input"]
fig = plot_forecast_map(forecast.sum(axis=(0,2)),targets=targets,title="ETAS：背景＋第一代，測試期合計")
apply_layout(fig)
display(Markdown(f"{forecast.shape[0]} 窗合計期望 {forecast.sum():.2f} 顆，目標事件 {len(targets)} 顆。這是逐窗更新後的合計。只疊圖不評分，評分在第 16–18 章。"))
fig

# %% [markdown]
# ### 從圖例讀到可比較的預報
#
# 一個紅點代表一次事後觀測。
# 發報時，模型只能使用它出現之前的事件。
# 讀圖時要把這兩個時間角色分開。
# 預報顏色才保有當時資訊下的意義。
#
# 三個月圖只累積該窗的期望數。
# 十年圖則累積所有窗的期望數。
# 同一種深藍色在兩張圖上可能代表不同數量。
# 因此，比較顏色前要先讀右側色條。
# 滑鼠移到方格時，可查看原始期望數。
#
# 紅點集中在深色格時，位置配置可能有幫助。
# 模型同時還要交代總共預期幾顆。
# 預報的數量與位置是兩個檢驗面向。
# 後續章節會分開處理這兩個問題。
#
# 兩參數估計已完成從目錄到概似的連結。
# 快取地圖接著完成從期望數到圖像的連結。
# 兩段連結使用的門檻與近似已列在規格卡。
# 讀者可以據此核對每個數字的來源。

# %% [markdown]
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
("目錄與範圍", "HORUS；S 提供歷史，R 接受預報；深度 ≤40 km"),
("估計示範", f"1990–2011；M≥4.0；ν={nu_hat:.3f}、K={K_hat:.4f}；限距近似"),
("正式預報參數", "Biondini 表 3；名目輸入 2.5；非本章重估值"),
("本站 ETAS 差異", "背景＋已知歷史第一代；未含窗內新事件再繁殖"),
("輸出", "每窗×177 格×25 規模箱；2012–2021 滾動更新"),
("數量檢驗", "待填（第 16 章）"),
("位置與規模檢驗", "待填（第 17 章）"),
])))

# %% [markdown]
# 三個月之外，地震活動還有其他尺度。
# 下一章觀察大地震之前的活動變化：
# {doc}`Ψ 現象 <14_psi_precursory_scale>`。
#
# ## 參考資料與延伸閱讀
#
# - Jalilian（2019），[時空 ETAS 模型與估計](https://www.jstatsoft.org/article/view/v088c01)。免費全文。先讀條件強度，再讀親代機率。
# - Ogata 與 Zhuang（2006），[時空 ETAS 的空間核](https://doi.org/10.1016/j.tecto.2005.10.016)。[免費作者全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。對照規模如何改變影響範圍。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。[免費機構全文](https://www.earth-prints.org/handle/2122/17084)。表 3 提供參數。附錄列出模型公式。
# - Mancini 與 Marzocchi（2023 線上發表），[simplETAS](https://doi.org/10.1785/0220230199)。[免費作者程式](https://github.com/smancini2/simplETAS)。讀固定參數如何降低估計負擔。
