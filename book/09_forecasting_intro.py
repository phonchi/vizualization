# %% [markdown]
# # 9. 地震預報的問題設定
#
# 第一部的結尾（{doc}`第 8 章 <08_explore_ideas>`）留下一個誠實但不
# 舒服的結論：想在單一觀測量裡找到可靠的地震前兆，非常困難；前兆
# 研究的歷史是一座墓園。那麼科學就對地震束手無策了嗎？
#
# 不是。過去三十年，地震科學走通了另一條路——放棄「說出下一場大
# 地震的時間地點」，改為回答一個謙虛得多、但答得出來的問題：
#
# > **在這個區域、這段時間內，發生某規模以上地震的機率是多少？**
#
# 這條路叫**地震預報（earthquake forecasting）**。它的材料不是神祕
# 的前兆訊號，而是你在{doc}`第 5 章 <05_seismic>`已經摸過的東西：
# 地震目錄本身。
#
# 這一章是第二部的問題設定。它不推導任何模型，只做四件事：把
# 「預測／預報／預警」三個詞分乾淨；說清楚我們憑什麼能預報；
# 展示一份預報實際上長什麼樣、有哪四種輸出形式；最後給你一張
# 十五章的地圖，以及讀它們的方式。
#
# ## 9.1 預測、預報、預警：三個詞的嚴格定義
#
# 中文日常把這三個詞混用，但在地震科學裡它們指三件不同的事，混淆
# 是最常見的溝通災難來源。
#
# | 詞 | 問的問題 | 時間尺度 | 現況 |
# |---|---|---|---|
# | 預測 prediction | 某時某地將發生規模 X 的地震 | — | **做不到** |
# | 預報 forecast | 此區、此期間、M≥X 的機率是多少 | 天～數十年 | **多國作業化運轉** |
# | 預警 early warning | 地震已發生，破壞性震波幾秒後到 | 秒 | **台灣世界前段班** |
#
# **預測**要求確定性，而且範圍要窄到可以據此疏散（Main 1999）。它的
# 失敗不是沒試過：美國在 Parkfield 斷層段守了十幾年，預測的地震遲到
# 十年；各種前兆方案（第 8 章那座墓園）在事後檢驗中一一倒下。台灣
# 中央氣象署自己的業務回顧給出的數字同樣誠實——短期前兆觀測中
# 「可視為成功發現前兆的比例」在**兩成以下**（蕭乃祺，約 2019）。
#
# **預警**是台灣民眾最熟悉的那個：手機在搖晃前幾秒響起。它的物理很
# 單純——地震**已經發生**，只是破壞性的 S 波還在路上，搶的是震波速
# 與電磁波速之間的幾秒到幾十秒。2018 年花蓮地震後 17 秒發布、20 秒
# 觸達手機，這是台灣真正領先世界的部分。但**預警不預測任何還沒發生
# 的事**，它與本書第二部要講的東西完全無關。
#
# **預報**問的則是「地震還沒發生時，接下來發生的機率」。它給的不是
# 答案而是機率，時間尺度從幾天到幾十年。這三者的物理、方法、難度
# 完全不同：台灣在預警上領先全球，在預報上則和所有國家一樣——
# 只能給機率。
#
# 那機率預報憑什麼給得出來？憑地震目錄裡最強、最穩定的一個訊號。
#
# ## 9.2 可預報性住在哪裡
#
# 地震的行為落在兩個極端之間。一端是**完全隨機**：像公正銅板一樣
# 互不相關，此時最好的預報就是長期平均率，任何模型都不可能贏過
# 它。另一端是**完全決定論**：只要量夠準就能算出下一次破裂，此時
# 預測是可能的，只是我們量不夠準。
#
# 真實地震落在中間，而且我們相當清楚**可預報性住在哪一格**：時空
# 叢集行為（Omori 1894；Kagan & Jackson 1991）與叢集內的規模頻率
# 分布（Gutenberg & Richter 1944）。這句話值得反覆強調——我們能
# 預報的其實是**叢集**，不是「地震」。
#
# 叢集不是印象，是可以量的。如果地震像 Poisson 過程那樣互不相關，
# 事件間隔時間應該服從指數分布。拿 2024 年春天的目錄（第 5 章用過
# 的那份，含 0403 花蓮主震與餘震序列）檢查：

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["hide-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
t_sorted = cat.time.sort_values()
dt_hours = t_sorted.diff().dt.total_seconds().dropna() / 3600
mean_dt = dt_hours.mean()

# 過度離散：每日事件數的變異數 / 平均
daily = t_sorted.dt.floor("D").value_counts().sort_index()
fano = daily.var() / daily.mean()

edges = np.arange(0, 12.25, 0.25)
hist, _ = np.histogram(dt_hours, bins=edges, density=True)
centers = edges[:-1] + 0.125

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.11,
                    subplot_titles=("事件間隔分布", "每日事件數：離散程度"))
fig.add_trace(go.Bar(x=centers, y=hist, name="觀測", marker_color=ACCENT,
                     opacity=0.75), row=1, col=1)
fig.add_trace(go.Scatter(x=centers, y=np.exp(-centers / mean_dt) / mean_dt,
                         mode="lines", name="同平均率的 Poisson",
                         line=dict(color=PALETTE[3], width=2.5)), row=1, col=1)
n_max = int(daily.max())
obs_hist, obs_edges = np.histogram(daily.values, bins=np.arange(0, n_max + 5, 4))
from scipy import stats as _st
pois = _st.poisson.pmf(np.arange(0, n_max + 4), daily.mean())
pois_binned = [pois[a:b].sum() for a, b in zip(obs_edges[:-1], obs_edges[1:])]
ctr = (obs_edges[:-1] + obs_edges[1:]) / 2
fig.add_trace(go.Bar(x=ctr, y=obs_hist / obs_hist.sum(), name="觀測（每日）",
                     marker_color=ACCENT, opacity=0.75, showlegend=False),
              row=1, col=2)
fig.add_trace(go.Scatter(x=ctr, y=pois_binned, mode="lines+markers",
                         name="Poisson 預期", line=dict(color=PALETTE[3]),
                         showlegend=False), row=1, col=2)
fig.update_xaxes(title_text="與前一事件的間隔（小時）", row=1, col=1)
fig.update_yaxes(title_text="機率密度", type="log", row=1, col=1)
fig.update_xaxes(title_text="每日事件數", row=1, col=2)
fig.update_yaxes(title_text="比例", type="log", row=1, col=2)
apply_layout(fig, title=f"地震不是隨機撒點（每日事件數的變異數／平均 = "
                        f"{fano:.0f}，Poisson 應為 1）",
             height=420, hovermode="x")
fig

# %% [markdown]
# 左圖：觀測（藍）在短間隔端遠高於同平均率的 Poisson 預期（黃線）
# ——大量地震擠在前一個地震後的幾分鐘到一小時內。右圖用另一個角度
# 說同一件事：Poisson 過程的事件數變異數應該等於平均（Fano 因子
# 為 1），實測卻是它的幾十倍。**過度離散就是叢集的指紋**，而這個
# 指紋正是所有預報模型的技巧來源。
#
# 叢集有兩條百年經驗律撐腰，你在第 5 章都算過：
#
# - **Omori–Utsu 律**：餘震率隨主震後時間衰減，$n(t)=K/(t+c)^p$，
#   $p$ 典型接近 1；
# - **Gutenberg–Richter 律**：$\log_{10}N=a-bM$，$b$ 典型接近 1，
#   規模每小一級、數量約多十倍。
#
# 把兩條合起來，就能回答「接下來這一週，這附近發生 M≥5 的機率」
# ——Omori 給「還會有多少餘震」，GR 給「其中多大比例是大的」。
# **現代地震預報的核心引擎，本質上就是這兩條定律的各種精緻組合。**
# 這也立刻說明了一件事：現行模型最擅長的是叢集（餘震），而不是
# 憑空冒出的大地震。第二部大半的篇幅在把這兩條律推到極致；
# {doc}`第 15 <15_psi_phenomenon>`、{doc}`16 章 <16_eepas_ppe>`
# 則會介紹一條不同的路——不問「這個地震會觸發什麼」，改問
# 「這個地震預示著什麼」。
#
# ## 9.3 前震、主震、餘震只是事後標籤
#
# 這裡藏著一個反直覺、但對防災至關重要的觀念。
#
# 目前**沒有任何已知的物理量**能在當下區分「這是餘震」還是「這是
# 更大地震的前震」（Felzer et al. 2004）；這個分類只能在序列結束後
# 回頭貼上。序列進行中，永遠無法排除更大的還在後面。
#
# 統計上這不是杞人憂天：1980–2019 年間，全球 M≥6 地震有超過一成
# 在 60 天、100 公里的時空窗內被**更大**的地震跟上（Taroni 2023）。
# 2016 年熊本地震是最痛的一課——M6.5 發生後，日本氣象廳依慣例
# 發布「餘震」預報，28 小時後真正的主震 M7.3 才來。此後 JMA 認定
# 原有程序失效，修改了整套發布方式。
#
# 所以機率預報不只是「算得出來」，它還是**唯一誠實的說法**：序列
# 進行中，科學能說的是「接下來一週再來一個 M≥6 的機率是百分之
# 幾」，而不是「餘震會慢慢變小，請放心」。這個觀念會在
# {doc}`第 13 章 <13_etas_structure>`得到數學形式——ETAS 的世界裡
# 根本沒有「主震」這個概念，只有背景事件與被觸發事件，而前震
# 是模型自己長出來的湧現行為。
#
# ## 9.4 一份預報長什麼樣：四種輸出
#
# 「地震預報」這個詞底下其實有四種不同的輸出，它們可以互相換算，
# 但溝通時混用會出事。從最原始到最貼近使用者：
#
# 1. **率密度** $\lambda(t,x,y,m)$：單位時間、單位面積、單位規模的
#    期望事件數。這是模型真正輸出的東西，也是
#    {doc}`第 10 章 <10_point_process>`整章要建立的語言。
# 2. **期望數** $\Lambda$：把率密度對一個時空規模盒子積分，
#    $\Lambda=\iiint\lambda\,dt\,dA\,dm$。
# 3. **機率** $P$：至少發生一次的機率。若盒內事件近似 Poisson，
#
#    $$P(N\ge1)=1-P(N=0)=1-e^{-\Lambda}$$
#
#    這條式子是預報溝通的主力，注意兩件事：$\Lambda$ 小時
#    $P\approx\Lambda$（機率與期望數幾乎相等），$\Lambda$ 大時 $P$
#    飽和趨近 1（再多的期望數也只能說「幾乎確定」）。
# 4. **地動超越機率**：使用者真正在乎的量——「我家會不會晃到某個
#    程度」。要從第 3 項再走一步，得配上地動預測方程，這是
#    {doc}`第 21 章 <21_psha>`的主題。
#
# 用一個合成情境把四者一次算給你看：

# %% tags=["hide-input"]
rng = np.random.default_rng(9)
days = np.linspace(0.01, 30, 400)
K, c_om, p_om = 12.0, 0.05, 1.10           # Omori 參數（示意）
lam_t = K / (days + c_om) ** p_om          # 率密度（每天，M≥4）
cum = np.concatenate([[0], np.cumsum(np.diff(days) * lam_t[:-1])])
win = 7.0
Lam_win = np.array([K / (p_om - 1) * ((t + c_om) ** (1 - p_om)
                                      - (t + win + c_om) ** (1 - p_om))
                    for t in days])        # 未來 7 天的期望數
P_win = 1 - np.exp(-Lam_win)               # 未來 7 天至少一次
P_shake = 1 - np.exp(-Lam_win * 0.12)      # 乘上「該事件造成場址超越門檻」的比例

fig = make_subplots(rows=2, cols=2, vertical_spacing=0.14, horizontal_spacing=0.11,
                    subplot_titles=("① 率密度 λ(t)", "② 未來 7 天期望數 Λ",
                                    "③ 未來 7 天機率 P = 1 − e^(−Λ)",
                                    "④ 場址震度超越機率"))
for (y, r, c, col) in [(lam_t, 1, 1, ACCENT), (Lam_win, 1, 2, PALETTE[2]),
                       (P_win, 2, 1, PALETTE[1]), (P_shake, 2, 2, QUAKE_COLOR)]:
    fig.add_trace(go.Scatter(x=days, y=y, mode="lines",
                             line=dict(color=col, width=2.5), showlegend=False),
                  row=r, col=c)
fig.update_yaxes(type="log", row=1, col=1); fig.update_yaxes(type="log", row=1, col=2)
fig.update_xaxes(title_text="主震後天數", row=2, col=1)
fig.update_xaxes(title_text="主震後天數", row=2, col=2)
apply_layout(fig, title="同一份預報的四種輸出（合成情境）",
             height=560, hovermode="x")
fig

# %% [markdown]
# 四條曲線描述的是**同一件事**，只是換了單位與問法。注意 ③ 與 ②
# 在早期幾乎重合、在後期一起下降——因為 $\Lambda$ 小的時候
# $1-e^{-\Lambda}\approx\Lambda$。而 ④ 永遠低於 ③：不是每一個
# M≥4 地震都會讓你家超過某個震度門檻。
#
# 現實中的預報也是這樣呈現的。紐西蘭 GeoNet 的公開預報表每一格
# 同時給**平均期望數、可能範圍、以及至少發生一次的機率**三個數字
# ——正是上面的 ② 與 ③。為什麼要三個一起給、而不是只給機率？
# 這是機率溝通的實證研究結果，{doc}`第 22 章 <22_operational_systems>`
# 會展開。
#
# ## 9.5 誰在算、誰能發：作業化與權威性
#
# 當這樣的機率被**自動化、定期產出，並由權責機關發布**，就叫做
# **作業化地震預報（Operational Earthquake Forecasting, OEF）**。
# 這個詞是 2009 年義大利 L'Aquila 地震（$M_w$ 6.3，約 300 人罹難）
# 之後，國際地震預報委員會（ICEF）給出的方向（Jordan et al. 2011）：
# 與其在「能不能預測」上空轉，不如把**當下已知的機率**用權威、
# 透明、常態化的方式交到社會手上。
#
# 定義裡有兩個關鍵字，值得現在就記住：
#
# - **作業化（operational）**：借自氣象學——自動、及時、持續運轉，
#   不是研究者事後發表論文。
# - **權威性（authoritative）**：資訊來自依法負有職責的機關。
#   **權威性不能自我宣稱，必須由制度授予**——這句話把 OEF 從純
#   技術問題變成制度問題。
#
# 世界上目前有三種代表性做法（義大利、紐西蘭、美國），它們算機率
# 的能力相近，分歧全在「要不要、怎麼給社會」。這是
# {doc}`第 22 章 <22_operational_systems>`的主場；這裡先記住一個
# 結論：**從模型走到系統，中間有七步，而第三步之後的瓶頸都不再
# 是統計。**
#
# ## 9.6 淘汰機制：CSEP 的三條紀律
#
# 讀到這裡你應該會想問：模型百百種，怎麼知道誰的機率值得信？
#
# 第 8 章的 8.2 節其實已經給過答案的雛形——**規則要在看到答案之前
# 定好**（第 3 關）、**要用沒有地震的時段檢驗誤報率**（第 4 關）。
# 把這兩條紀律放大成國際建制，就是 **CSEP**（Collaboratory for the
# Study of Earthquake Predictability）。它的三條紀律是：
#
# 1. **事前註冊**：模型參數、預報格式、目標資料來源，全部在觀測
#    之前釘死，達成零自由度的可否證陳述；
# 2. **前瞻檢驗**：模型對**未來**的地震給出正式預報，不是回頭
#    解釋過去；
# 3. **第三方執行**：由獨立測試中心用事先議定的統計檢驗打分數；
#    **「開發者自己相信自己的模型」不構成任何證據**，而且模型
#    換一個地區使用，必須重新檢驗。
#
# 這一刀切掉了第 8 章講的事後選擇偏誤：你不可能對還沒發生的地震
# p-hacking。CSEP 二十年來累積的最重要教訓也很直白——不少模型在
# 回溯測試中大放異彩，一到前瞻測試就原形畢露。你會在
# {doc}`第 19 章 <19_ensembles>`看到一個完整的十年案例：同一批
# 模型、同一個統計量，回溯 $+0.5$、前瞻 $-0.7$，符號整個翻轉。
#
# 檢驗的完整工具箱在{doc}`第 17 <17_testing_consistency>`與
# {doc}`第 18 章 <18_testing_comparison>`。第二部每次講到任何模型，
# 我們都會問同一個問題：<strong>它前瞻檢驗過了嗎？成績如何？</strong>這是這個
# 領域與前兆墓園最大的差別——它建立了淘汰機制。
#
# ## 9.7 第二部的地圖
#
# 接下來十四章分成五個群組，涵蓋從秒到數十年的時間尺度：

# %% tags=["hide-input"]
groups = [
    ("地基：語言與資料", "10 點過程｜11 完整度與 b 值｜12 叢集律與除叢",
     -2.6, 0.5, PALETTE[0]),
    ("短期：觸發", "13 ETAS I｜14 ETAS II｜22 作業化系統",
     -2.4, -0.6, PALETTE[1]),
    ("中期：前兆", "15 Ψ 現象｜16 EEPAS 與 PPE",
     -0.6, 1.0, PALETTE[2]),
    ("長期：複發與危害", "20 複發模型｜21 PSHA",
     0.7, 2.7, PALETTE[6]),
    ("橫貫：怎麼判斷好壞", "17–18 檢驗｜19 組合｜23 台灣",
     -2.6, 2.7, "#888888"),
]
fig = go.Figure()
for i, (name, chapters, lo, hi, color) in enumerate(groups):
    fig.add_trace(go.Bar(y=[name], x=[hi - lo], base=[lo], orientation="h",
                         marker_color=color, opacity=0.85 if i < 4 else 0.35,
                         text=chapters, textposition="inside",
                         insidetextanchor="middle", showlegend=False,
                         textfont=dict(size=11)))
fig.update_xaxes(title_text="時間尺度（年，log₁₀）",
                 tickvals=[-2.6, -1.4, 0, 1, 2.7],
                 ticktext=["天", "月", "年", "十年", "數百年"])
apply_layout(fig, title="第二部的五個群組：各自負責哪一段時間尺度",
             height=380, hovermode="closest")
fig

# %% [markdown]
# 逐章路線圖：
#
# | 章 | 主題 | 你會學到的核心工具 |
# |---|---|---|
# | 10 | 點過程 | 條件強度 $\lambda^*$、概似 $\sum\ln\lambda^*-\int\lambda^*$、模擬 |
# | 11 | 完整度與 b 值 | $M_c$ 三法、Aki MLE 與離散修正、差分法 |
# | 12 | 叢集律與除叢 | Omori–Utsu、Utsu–Seki、Båth、除叢的選擇效應 |
# | 13 | ETAS I | 分支比 $n$、世代分解、branching 模擬 |
# | 14 | ETAS II | 參數估計與診斷、隨機除叢、R–J 特例 |
# | 15 | Ψ 現象 | 累積規模異常 $C(t)$、自動辨識、時空取捨 |
# | 16 | EEPAS 與 PPE | 三個機率核、正規化 $\eta$、基準模型 |
# | 17 | 檢驗 I | Poisson 似然、N/M/S/cL、負二項與二元 |
# | 18 | 檢驗 II | 資訊增益、Molchan、統計功效、成本–損失 |
# | 19 | 模型組合 | 凸組合、乘法 hybrid、權重學習 |
# | 20 | 複發模型 | hazard function、BPT、應力釋放 |
# | 21 | PSHA | 危害積分、反聚合、時變危害 |
# | 22 | 作業化系統 | STEP 三層、三國系統、機率溝通 |
# | 23 | 台灣 | 在地基準值、四個缺口、如何做台灣的 CSEP |
#
# **怎麼讀這一部**。三個提醒：
#
# 第一，**數學是主體，程式碼不是**。每章都有完整推導——不是列出
# 結果，而是交代每一步怎麼來。所有圖都由摺疊的程式產生，你不需要
# 讀它們；但如果好奇某個數字怎麼算出來的，展開就看得到。
#
# 第二，**記號在第 10 章統一定案**。地震統計文獻的符號互相打架
# （同一個 $\mu$ 在三個模型裡指三件不同的事），本書為此定了一套
# 統一寫法並在 10.8 節列表。讀論文時記得回頭對照。
#
# 第三，**每章末尾有三個固定區塊**：常見誤解與陷阱、研究前沿與
# 未解問題、以及推導附錄。前兩者是散文，可以當休息；附錄放長代數，
# 第一次讀可以跳過。
#
# ## 9.8 常見誤解
#
# **「地震預報就是預測，只是講得含糊。」** 不是。預報給的是可以被
# 檢驗、可以被淘汰的機率陳述；預測給的是無法用統計檢驗的單一斷言。
# 兩者的科學地位完全不同——一個有前瞻測試制度，一個沒有。
#
# **「機率增益 10 倍，那應該很準了吧。」** 機率增益是**相對**量。
# 背景日機率萬分之一，抬高十倍還是千分之一。第二部會反覆看到這個
# 落差：增益很大、絕對機率很小，兩者同時為真。
#
# **「模型通過了檢驗，所以它是對的。」** 通過檢驗只代表資料不足以
# 拒絕它。{doc}`第 18 章 <18_testing_comparison>`會給你一個震撼的
# 案例：一個宣稱「地球上每個地方機率都一樣」的模型，通過了標準的
# 空間檢驗。
#
# **「台灣有地震預警，所以離預報不遠了。」** 兩者的技術路線毫無
# 交集（9.1 節）。預警強不代表預報強。
#
# ## 9.9 研究前沿：這個領域現在往哪裡走
#
# 三個方向會在後續章節逐一展開，這裡先給座標。
#
# **機器學習能不能贏過 ETAS？** 深度學習預測餘震位置的早期宣稱被
# 證明與遠更簡單的參數化資訊量相當（Mignan & Broccardo 2019）；
# 近年神經點過程（把條件強度換成 RNN 或 Transformer）開始能匹敵
# 甚至略勝，且對目錄非平穩性適應力更好——見
# {doc}`14.11 <14_etas_estimation>`。
#
# **b 值能不能當前震判別？** 前震紅綠燈（FTLS）主張追蹤大震後
# $b$ 值的時變可以一階判別「這是正常餘震序列還是前兆序列」。原理
# 有實驗室與野外支持，但 $b$ 值估計本身的偏差來源太多——見
# {doc}`11.11 <11_catalog_completeness_b>`。
#
# **當模型不存在時怎麼辦？** 2016 年 Kaikōura 之後，紐西蘭觀測到
# 隱沒帶廣泛的慢滑移，而當時沒有任何現成模型可用。他們的做法是
# 結構化專家徵詢——見 {doc}`22.8 <22_operational_systems>`。
#
# ## 從懷疑到制度
#
# 第一部教你對「看起來像前兆的異常」保持懷疑。第二部要教你的，是
# 對「看起來很準的模型」保持同樣的懷疑——只是這一次，懷疑有了
# 正式的工具：基準模型、前瞻檢驗、資訊增益、統計功效。
#
# 地震預報這門學問最了不起的地方，不在於它能算出機率，而在於它
# 建立了一套讓錯誤的機率無所遁形的制度。
#
# {doc}`下一章 <10_point_process>`我們先學這套制度所使用的語言：
# 點過程。它是後面十三章共用的數學骨架——條件強度、概似、模擬、
# 殘差診斷，全部從那裡長出來。
