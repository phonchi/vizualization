# %% [markdown]
# # 12. EEPAS 與 PPE：前兆尺度增加與中長期預報
#
# ETAS 的預報技巧集中在大地震**之後**的短期窗口；對「下一個大地震
# 什麼時候來」，它幾乎無話可說。這一章的模型把目光轉向大地震
# **之前**——它建立在一個累積了數十年觀測的經驗現象上：大地震
# 發生前，震源區的中小地震常常會先「變多、變大」。
#
# 這個現象叫 **Ψ（前兆尺度增加，precursory scale increase）**；
# 把它變成可運轉的預報模型，就是 **EEPAS**（Every Earthquake a
# Precursor According to Scale——「每個地震都是與其尺度相稱的
# 前兆」，Rhoades & Evison 2004）。EEPAS 是目前世界上少數真正
# 上線運轉的**中期**（數月到數十年）預報模型，也是紐西蘭公開
# 地震預報的中期成分。這一章同時要教它的忠實配角 **PPE**——
# 一個樸素到近乎無聊、卻誰都繞不開的基準模型。
#
# ## 12.1 Ψ 現象：大震之前，中小地震先變多變大
#
# Ψ 的辨識程序是這樣的：對一個已知的大地震，框選一個時空範圍，
# 畫出範圍內的規模–時間圖，再畫**累積規模異常**曲線 $C(t)$——
# 把每個地震的規模超額累加起來、扣掉平均趨勢。若震前存在一段
# 「規模與發生率同時抬升」的時期，$C(t)$ 會先下探再回升，
# 最低點就是 Ψ 的**起始點（onset）**。用合成資料畫出典型樣貌：

# %% tags=["remove-input"]
import plotly.io as pio
pio.renderers.default = "notebook_connected"

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit import load_taiwan_catalog
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

rng = np.random.default_rng(42)
MC = 4.0
# 平時期（0–14 年）：低發生率、規模上限低；前兆期（14–20 年）：率與規模同時抬升
t_bg = np.sort(rng.uniform(0, 14, 26))
m_bg = MC + rng.exponential(0.35, len(t_bg))
t_pre = np.sort(rng.uniform(14, 20, 30))
m_pre = MC + 0.4 + rng.exponential(0.45, len(t_pre))
t_all = np.concatenate([t_bg, t_pre, [20.0]])
m_all = np.concatenate([m_bg, m_pre, [7.0]])          # 20 年時主震 M7

# 累積規模異常 C(t)：規模超額的累積，扣掉線性平均趨勢
excess = m_all[:-1] - MC - 0.1
k = excess.sum() / 20.0
tc = np.linspace(0, 20, 400)
C = np.array([excess[t_all[:-1] <= t].sum() for t in tc]) - k * tc
onset = tc[np.argmin(C)]

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.55, 0.45],
                    vertical_spacing=0.05)
fig.add_trace(go.Scatter(x=t_all[:-1], y=m_all[:-1], mode="markers",
                         name="中小地震", marker=dict(size=7, color=ACCENT,
                                                  opacity=0.75)), row=1, col=1)
fig.add_trace(go.Scatter(x=[20], y=[7.0], mode="markers", name="主震",
                         marker=dict(size=14, color=QUAKE_COLOR,
                                     symbol="star")), row=1, col=1)
fig.add_trace(go.Scatter(x=tc, y=C, mode="lines", name="累積規模異常 C(t)",
                         line=dict(color=PALETTE[3], width=2)), row=2, col=1)
for r in (1, 2):
    fig.add_vline(x=onset, line_dash="dot", line_color="#1baf7a", row=r, col=1)
fig.add_annotation(x=onset, y=float(C.min()), text="Ψ onset", showarrow=True,
                   ax=-40, ay=30, row=2, col=1)
fig.update_yaxes(title_text="規模", row=1, col=1)
fig.update_yaxes(title_text="C(t)", row=2, col=1)
fig.update_xaxes(title_text="時間（年）", row=2, col=1)
apply_layout(fig, title="Ψ 現象示意（合成資料）：主震前規模與發生率同時抬升",
             height=520, hovermode="x")
fig

# %% [markdown]
# 從一次 Ψ 辨識可以量出三個量：前兆規模 $M_P$（前兆期最大三個
# 地震的平均規模）、前兆時間 $T_P$（onset 到主震）、前兆面積
# $A_P$。Evison & Rhoades 從四個地區的 47 個大地震歸納出三條
# 迴歸關係——**前兆的尺度處處跟著主震規模走**：
#
# $$M_m = 3.16 + 0.65\,M_P,\qquad
#   \log_{10} T_P = 1.36 + 0.40\,M_P,\qquad
#   \log_{10} A_P = 1.88 + 0.35\,M_P$$
#
# 讀出兩個直覺：前兆群中最大的地震通常比主震小約一個規模單位
# （$b_M=0.65$）；主震愈大、前兆愈早開始——換算下來 M6 的
# 前兆時間約四、五年，**M7 約十五年**。這就是「中長期」的由來。
#
# 但誠實的敘述必須馬上加上三個保留。第一，Ψ 的辨識傳統上是
# **手工**的：研究者已知主震在哪、多大，一邊調整框選範圍一邊
# 「最大化尺度增加」——第 8 章警告過的事後選擇偏誤全套到齊，
# 連原作者群近年都自承這種程序「可能被視為在滿足預設好的結論」
# （Christophersen et al. 2024）。第二，改用公開判準的自動演算法
# 重做後發現：Ψ 的辨識**不唯一**——同一主震平均能找出約 3 組
# 同樣合格、但 $T_P$ 與 $A_P$ 可差到十倍的辨識（模擬目錄中
# 甚至差兩個數量級）；而且 47 個
# 主震中有 13 個（28%）連一組都找不到。第三，物理機制至今未明。
#
# 那 Ψ 還算數嗎？算。同一份研究在物理模擬器（RSQSim——只放入
# 摩擦定律與斷層交互作用）生成的合成目錄上做了乾淨的對照實驗：
# 把發生時間隨機打散，訊號消失；把餘震拿掉，訊號仍在。而 Ψ
# 從來不是被寫進模擬器的規則——**它自己長了出來**，尺度關係
# 還與真實資料一致。結論的正確措辭是：
# 主震前確實存在一種不同於餘震衰減的時空叢集，但它的參數化
# 不唯一、機制未知——與 Omori 律被確立時的處境一模一樣。
#
# ## 12.2 EEPAS：每個地震都是與其尺度相稱的前兆
#
# Ψ 不能事前辨識（你不知道主震在哪），這似乎宣判了它對預報
# 無用。EEPAS 的核心創意是一步漂亮的側身：**那就不要辨識。**
# 假設**每一個**地震都是某個未來更大地震的前兆，它「預告」的
# 主震規模、等待時間、發生範圍，全部由它自己的規模 $m_i$ 按
# 三條 Ψ 迴歸給定——迴歸的殘差是常態的，所以三條迴歸直接
# 變成三個機率核：
#
# | Ψ 迴歸（回溯統計） | EEPAS 機率核（前瞻預報） | 直覺 |
# |---|---|---|
# | $M_m = a_M + b_M M_P$ | 規模核：常態，中心 $a_M+b_M m_i$ | 預告的主震多大 |
# | $\log T_P = a_T + b_T M_P$ | 時間核：對數常態，中位數 $\propto 10^{a_T+b_T m_i}$ | 要等多久 |
# | $\log A_P = a_A + b_A M_P$ | 空間核：二維常態，變異數 $\sigma_A^2\,10^{b_A m_i}$ | 落在多大範圍 |
#
# 總發生率密度是背景項加上所有過去地震的貢獻：
#
# $$\lambda(t,m,x,y) = \mu\,\lambda_0(t,m,x,y)
#   + \sum_{t_i \ge t_0} \eta(m_i)\, w_i\,
#     f(t \mid t_i, m_i)\, g(m \mid m_i)\, h(x,y \mid x_i, y_i, m_i)$$
#
# $\eta$ 是讓長期平均仍符合 GR 律的正規化函數；$w_i$ 是餘震
# 降權（用一個 ETAS 型模型算出「這個地震是獨立事件的機率」）；
# $\mu$ 是混合權重——「沒有可辨識前兆的地震」佔的比例，同時
# 也是模型的誠實儀表板：時變項找不到東西時，擬合就會把 $\mu$
# 推高、轉而依賴背景。直觀地說，**每個小地震都在未來的時空中
# 放下一個「機率包裹」**，包裹的大小、遠近、何時打開，由它的
# 規模決定：

# %% tags=["hide-input"]
a_T, b_T, sigT = 1.36, 0.40, 0.35            # 時間核（天，log10）
b_A, sigA = 0.35, 1.0                         # 空間核（km）
events = [(1.0, 20.0, 4.0), (3.0, 55.0, 5.0), (5.5, 80.0, 5.8)]

tg = np.linspace(0, 25, 300)                  # 年
xg = np.linspace(0, 110, 220)                 # km
TT, XX = np.meshgrid(tg, xg)
dens = np.zeros_like(TT)
for t0, x0, m in events:
    dt_days = np.clip((TT - t0) * 365.25, 1e-3, None)
    mu_T = a_T + b_T * m
    f = (np.exp(-0.5 * ((np.log10(dt_days) - mu_T) / sigT) ** 2)
         / (dt_days * sigT * np.log(10) * np.sqrt(2 * np.pi)))
    f[TT <= t0] = 0
    sig_x = sigA * 10 ** (b_A * m / 2)
    h = np.exp(-0.5 * ((XX - x0) / sig_x) ** 2) / (sig_x * np.sqrt(2 * np.pi))
    dens += f * h * 10 ** (0.5 * m)           # 權重隨規模放大（示意）

fig = go.Figure(go.Heatmap(x=tg, y=xg, z=np.sqrt(dens), colorscale="Blues",
                           showscale=False))
fig.add_trace(go.Scatter(x=[e[0] for e in events], y=[e[1] for e in events],
                         mode="markers+text", name="過去的地震",
                         text=[f"M{e[2]:.1f}" for e in events],
                         textposition="middle left",
                         marker=dict(size=[8, 12, 16], color=QUAKE_COLOR)))
apply_layout(fig, title="EEPAS 的機制：每個地震在未來放下一個機率包裹",
             xaxis_title="時間（年）", yaxis_title="位置（km，一維剖面）",
             hovermode="closest", height=460)
fig

# %% [markdown]
# 三個地震（紅點）各自貢獻一團往未來延伸的機率密度：規模愈大，
# 包裹愈大、愈晚達到高峰、空間上攤得愈開。M4 的包裹一兩年內
# 就過期；M5.8 的包裹要十幾年後才完全打開。把一個地區幾十年
# 目錄裡**所有**地震的包裹疊起來，就是 EEPAS 的預報地圖。
#
# 這一步側身的代價與收穫都要講明。收穫是可操作性：不需要人工
# 辨識，全自動、可檢驗。代價是預報强度被稀釋——原作者自己稱
# EEPAS 是「弱預報模型」：它給的不是警報，而是把目標地震發生
# 機率相對基準抬高數倍的率密度。但這個「數倍」經得起最嚴格的
# 考驗：模型在紐西蘭擬合後，**參數一個不動**（僅重估 b 值與 PPE 基線）移植到加州做獨立
# 測試，26 年下來相對基準模型的累積概似比達到 $10^{15}$——
# 這種數字不可能靠過擬合湊出來。
#
# ## 12.3 PPE：樸素的基準
#
# 上式裡的背景項 $\lambda_0$ 是 **PPE**（Proximity to Past
# Earthquakes，「與過去地震的鄰近性」）：未來地震的發生率，
# 正比於把過去每個震央攤上一個平滑核之後的疊加，規模服從 GR 律，
# 不對時間叢集做任何建模。翻成白話：**「地震會發生在以前發生過地震的地方
# 附近」**——第 9 章那張平滑地圖的正式版。用台灣長期目錄的
# M≥5 事件畫出來：

# %% tags=["hide-input"]
cat_long = load_taiwan_catalog(min_ml=5.0)
d_km, s_bg = 15.0, 1e-4                       # 平滑距離、遠域常數（泛用示意值）
step = 0.1
lons = np.arange(119.0, 123.5, step)
lats = np.arange(21.0, 26.0, step)
LON, LAT = np.meshgrid(lons, lats)
dens = np.full_like(LON, 0.0)
ev = cat_long[["longitude", "latitude", "ML"]].to_numpy()
for lo, la, m in ev:
    r2 = ((LON - lo) * 111 * np.cos(np.radians(la))) ** 2 + ((LAT - la) * 111) ** 2
    dens += (m - 5.0 + 0.1) * (1 / (np.pi * (d_km ** 2 + r2)) + s_bg)  # 權重∝規模超額

fig = go.Figure(go.Heatmap(x=lons, y=lats, z=np.log10(dens), colorscale="Blues",
                           colorbar=dict(title="log₁₀ 相對率")))
apply_layout(fig, title=f"PPE 式平滑地震度（台灣 1973–2025，M≥5 共 {len(ev)} 筆，"
                        f"平滑距離 {d_km:.0f} km）",
             xaxis_title="經度", yaxis_title="緯度",
             yaxis_scaleanchor="x", hovermode="closest", height=560)
fig

# %% [markdown]
# 東部外海與西南部亮起來——過去五十年地震密集的地方。PPE 幾乎沒有
# 時間結構、沒有任何物理，卻是地震預報界最重要的一類模型，
# 因為它扮演三個角色：EEPAS 的背景項、餘震權重的參照，以及
# 最重要的——**所有花俏模型都必須贏過的基準線**。一個模型宣稱
# 有預報能力，第一個問題永遠是：「你贏過平滑地震度多少？」
#
# ## 12.4 二十年旅程：EEPAS 學到的事
#
# EEPAS 自 2004 年起被移植到加州、日本、義大利，並成為
# 紐西蘭官方混成預報模型的中期成分（基督城震後重建的 50 年
# 危害度模型、Kaikōura 之後的 100 年模型都有它）。二十年的
# 成績單與教訓，濃縮成五條：
#
# 1. **增益隨目標規模上升。** 日本本土的正式 CSEP 測試中，
#    EEPAS 相對 PPE 的機率增益在 M4.0–4.5 為 1.27、M4.5–5.0 為
#    1.52、M5 以上為 2.77——目標愈大、Ψ 訊號愈清楚，完全符合
#    模型理念。但注意增益是集體統計量：單一地震的增益從 0.05 到
#    20 都有，**平均增益不承諾任何個案**。
# 2. **短期看觸發、中長期看前兆。** 義大利的擬前瞻實驗把 ETAS
#    與 EEPAS 放上同一張表：三個月內 ETAS 勝，**六個月到十年
#    EEPAS 勝**（顯著性集中在五年與十年期）。更深刻的是 Janus 混合實驗：即使預報視窗縮到
#    零，EEPAS+ETAS 的混合仍勝過純 ETAS——表示**觸發串級與
#    前兆尺度增加是兩種大致獨立的可預報性來源**，正確做法不是
#    二選一，而是按時間尺度加權混合（第 14 章）。
# 3. **參數不可跨脈絡移植。** 跨規模級距套用參數，增益會掉到
#    比基準還糟（日本：0.61）；跨區域更不行——加州的平均前兆
#    時間約為日本關東的 6 倍，前兆面積卻只有 1/6。背後有物理：
#    模擬實驗顯示**前兆時間與構造加載速率成反比**，穩定大陸
#    內部（如澳洲）的前兆時間長到超出目錄長度。
# 4. **資料限制要建模，不能忽略。** 目錄太短，早年的長前兆沒
#    被記到，模型就會系統性偏差。近年的補償方法（LEEPAS／
#    FLEEPAS 家族）把「目錄能看到多少前兆」寫成公式修正回去，
#    讓預報視窗能從數月拉到數十年而不損失太多資訊。
# 5. **假設要能被自己的檢驗推翻。** 曾有一個變體假設「地震率
#    同時控制前兆時間與面積」，在正式前瞻測試中輸給標準版，
#    假設被放棄——統計地震學少見的乾淨否證案例，值得敬佩。
#
# ## 12.5 時空取捨：同一個前兆的兩種讀法
#
# Ψ 還藏著一個微妙的結構。前兆時間與前兆面積各自對主震規模的
# 迴歸都不錯，但**彼此**的相關性很差；而兩者的乘積 $A_P T_P$
# 對主震規模的解釋力反而最高。原因是 Ψ 的辨識存在**時空取捨
# （space-time trade-off）**：同一個主震，可以用「大面積、
# 短前兆時間」辨識，也可以用「小面積、長前兆時間」辨識，兩者
# 一樣合法——多組辨識沿著 $A_P \times T_P \approx$ 常數的
# 抵換線分布。EEPAS 的擬合也獨立看到同一件事：$\sigma_A^2$ 與
# $10^{a_T}$ 的乘積近似守恆，把參數沿抵換線挪動，預報表現
# 幾乎不變。
#
# 這件事有三層教學價值。第一，它解釋了為什麼不同團隊、不同
# 實作擬合出不同的 $a_T$、$\sigma_A$ 卻得到幾乎相同的概似——
# **參數不可辨識性**不是 bug，是模型面在這個方向上天生平坦。
# 第二，它是一堂活的統計課：研究者對取捨斜率做正反兩個方向的
# 迴歸，得到 $-0.56$ 與 $-0.33$（後者換算回同一座標系約為 $-3$）
# ——兩個估計都因量測誤差被**回歸稀釋**而偏淺，真值被夾在中間，
# 與理論的 $-1$ 相容；而 $T_P$ 與 $A_P$ 的相關性在「同一主震內」
# 是負的、「跨主震」是正的——不分群畫散布圖會得到完全相反的
# 結論（Simpson 悖論的地震版）。第三，取捨可以被利用：沿抵換線
# 取三組參數混合，**不增加任何自由參數**，在獨立測試期顯著
# 勝出——比任何「加參數後變好」都更有說服力。
#
# ## 12.6 反思：把第 8 章的功課寫成模型
#
# 回頭看，EEPAS 像是把第一部{doc}`第 8 章 <08_explore_ideas>`的
# 方法論功課逐條寫成了數學：
# 事後選擇偏誤？——那就放棄事後辨識，讓每個地震自動貢獻。
# 判準要事先定好？——判準就是三條公開的迴歸式與最大概似。
# 要用獨立資料檢驗？——紐西蘭擬合、加州獨立測試，再送進
# CSEP 前瞻擂台。要誠實面對失敗？——28% 找不到 Ψ、單一辨識
# 不唯一、機制未知，全部寫進論文。它不是因為「相信前兆」而
# 值得教，而是因為它示範了**如何把一個爭議中的前兆現象，
# 轉化成可檢驗、可否證、可上線的預報模型**。
#
# 對台灣，這一章留下一個明確的座標：短期（ETAS）台灣已有本土
# 參數與作業化實測（第 11 章），長期（PSHA）也有國家級模型
# （第 16 章），**中期這一段目前是空白**——而 EEPAS 的在地化
# 工作，正由台灣的研究團隊進行中。填上這一格之後，
# 台灣才有材料做下一章的事：把不同時間尺度的模型組裝成一個
# 隨時間演化的作業化系統。{doc}`下一章 <13_step_oef>`我們就去看
# 世界各國是怎麼把模型接上社會的——從 STEP 到作業化地震預報。
