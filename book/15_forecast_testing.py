# %% [markdown]
# # 15. 預報好不好，怎麼算？CSEP 檢驗
#
# 前面六章反覆出現同幾個詞：前瞻測試、資訊增益、N-test。這一章
# 正式打開這個工具箱。它可能是第二部最重要的一章——因為模型
# 會過時，但「怎麼判斷一個預報值不值得信」這套方法，你會用一輩子。
#
# 先復習動機。地震預測研究曾長期背負罵名，根本病因不是想法
# 不夠聰明，而是**結果無法重現、技巧無法獨立驗證**：模型在
# 自己人的回溯測試裡百戰百勝，換一份資料就原形畢露。CSEP 的
# 解方（第 9 章）是把預報變成**零自由度的可否證陳述**：模型
# 參數、預報格式、目標資料來源，全部在觀測**之前**釘死，再由
# 第三方測試中心打分數。這一章講的就是打分數的方法。
#
# ## 15.1 似然：共同貨幣
#
# 標準的網格化預報把測試區切成空間格 × 規模箱，模型對每箱交出
# 期望地震數 $\lambda$。假設箱內事件數服從 Poisson 分布，單一
# 箱觀測到 $\omega$ 個事件的對數似然是
#
# $$\ln P(\omega \mid \lambda) = -\lambda + \omega \ln\lambda - \ln(\omega!)$$
#
# 把全部箱加總，就得到「這個模型賦予**實際發生的那些地震**多少
# 機率」的總分。這就是整個檢驗科學的原子：後面每一種檢驗——
# 一致性檢驗、比較檢驗、資訊增益——都是這條式子的變奏。似然
# 提供了跨模型的**共同貨幣**：不管你的模型是 ETAS、EEPAS 還是
# 神經網路，最後都被換算成同一種分數。
#
# 而有了共同貨幣，就能定義這幾章一直在用的**每地震資訊增益**：
# 兩個模型的對數似然差除以地震數，$I = (\ln L_X - \ln L_Y)/N$；
# 取指數就是機率增益。看到任何資訊增益數字，永遠先問一句：
# **「相對於哪個基準？」**——相對「均勻隨機」的增益與相對
# 「當今最佳模型」的增益，是完全不同量級的成就。
#
# ## 15.2 一致性檢驗：N、M、S、L
#
# 一致性檢驗問的是「模型與觀測資料合不合」。CSEP 的設計哲學
# 是**把預報拆開來、一次檢查一個面向**——每個檢驗都刻意丟掉
# 其他資訊：
#
# | 檢驗 | 檢什麼 | 刻意抽掉什麼 |
# |---|---|---|
# | **N-test** | 總地震數對不對 | 空間、規模 |
# | **M-test** | 規模分布對不對 | 空間、總數 |
# | **S-test** | 空間圖形對不對 | 規模、總數 |
# | **cL-test** | 空間×規模聯合結構 | （不正規化總數） |
#
# （原始的 L-test 把全部混在一起檢，結果對總數過度敏感——率
# 抓錯就整體被拒絕，看不出是哪裡錯——因此實務上被 cL-test
# 取代。）觀測分數沒有解析分布，做法是模擬：按預報抽一萬份
# 合成目錄、算出分數的分布，看觀測值落在第幾分位。以最直觀的
# N-test 為例：

# %% tags=["remove-input"]
import plotly.io as pio
pio.renderers.default = "notebook_connected"

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from scipy import stats

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

n_exp = 20.0                                     # 模型的期望地震數
k = np.arange(0, 45)
pmf = stats.poisson.pmf(k, n_exp)
lo, hi = stats.poisson.ppf([0.025, 0.975], n_exp)

fig = go.Figure(go.Bar(x=k, y=pmf, name="模型的預報數分布（Poisson）",
                       marker_color=ACCENT, opacity=0.75))
fig.add_vrect(x0=lo - 0.5, x1=hi + 0.5, fillcolor="#1baf7a", opacity=0.12,
              line_width=0, annotation_text="95% 區間")
fig.add_vline(x=24, line_color="#1baf7a", line_width=3,
              annotation_text="觀測 A：24（通過）")
fig.add_vline(x=38, line_color=QUAKE_COLOR, line_width=3,
              annotation_text="觀測 B：38（不一致）")
apply_layout(fig, title="N-test：觀測總數落在模型預報數分布的哪裡？",
             xaxis_title="目標地震數", yaxis_title="機率", hovermode="x")
fig

# %% [markdown]
# 模型說期望 20 個。觀測 24 個（綠線）落在 95% 區間內——與
# 模型一致；觀測 38 個（紅線）落在遠端尾巴——模型顯著低估了
# 地震率。其他檢驗同理，只是分數換成空間或規模維度的似然。
#
# 兩個必須跟檢驗本身一起教的細節：
#
# - **Poisson 假設會被叢集打爆。**真實地震強烈叢集（第 9 章
#   第一張圖），總數的變異數遠大於 Poisson 所允許；一個群震
#   序列塞爆少數幾格，就能讓 S-test 對所有模型亮紅燈。因此
#   實務上有兩個放寬版：總數改用**負二項分布**（變異數由歷史
#   資料估計）、空間似然改用**二元版**（只問這格「有沒有」，
#   不問「有幾顆」）。
# - **檢驗是診斷，不是判決。**CSEP 自己的立場是：不因為某個
#   檢驗失敗就正式淘汰模型，而是把分位數分數當診斷指標——
#   指出模型與資料在哪個面向不合，值得回去追究。
#
# ## 15.3 比較檢驗：T-test 與資訊增益
#
# 一致性檢驗是「模型 vs 資料」；真正決定排名的是「模型 vs
# 模型」。**T-test** 對每一顆觀測地震，比較兩個模型在該地震
# 所在箱的對數似然，平均起來（並修正兩模型期望總數的差異）
# 得到每地震資訊增益 IGPE，再用 Student-t 給信賴區間。判讀
# 規則簡單到可以畫成一張圖：

# %% tags=["hide-input"]
models = ["模型甲", "模型乙", "模型丙"]
ig = [0.45, 0.12, -0.38]
err = [0.20, 0.25, 0.22]
colors = ["#1baf7a", PALETTE[3], QUAKE_COLOR]
fig = go.Figure()
for m, g, e, c in zip(models, ig, err, colors):
    fig.add_trace(go.Scatter(x=[m], y=[g], error_y=dict(type="data", array=[e]),
                             mode="markers", name=m,
                             marker=dict(size=12, color=c)))
fig.add_hline(y=0, line_dash="dash", line_color="#888",
              annotation_text="基準模型")
apply_layout(fig, title="T-test 判讀：誤差棒跨過零線就是平手",
             yaxis_title="每地震資訊增益（相對基準）",
             hovermode="closest", showlegend=False)
fig

# %% [markdown]
# 模型甲顯著優於基準（區間整個在零之上）；模型乙的誤差棒跨過
# 零線——**統計上無法區分，就是平手**，不管點估計多好看；
# 模型丙顯著較差。搭配的 W-test（無母數版本，Rhoades et al.
# 2011）用來確認 T-test 的結論不依賴常態假設。
#
# 還有一族完全不同語言的檢驗：**警報式**。把預報率用一個門檻
# 二值化——高於門檻的格子「發警報」——就得到命中／漏報／
# 誤報的列聯表；把門檻從嚴掃到鬆，畫出**Molchan 圖**：漏報率
# 對警報覆蓋比例的軌跡。

# %% tags=["hide-input"]
rng = np.random.default_rng(42)
n_cells = 2000
lam = rng.lognormal(0, 1.6, n_cells)             # 真實率場（高度不均勻）
lam /= lam.sum()
events = rng.choice(n_cells, size=200, p=lam)    # 事件按真實率落格

def molchan(alarm_score):
    order = np.argsort(-alarm_score)
    tau = np.arange(1, n_cells + 1) / n_cells
    hits = np.cumsum(np.isin(order, events, assume_unique=False)
                     * np.array([np.sum(events == i) for i in order]))
    return tau, 1 - hits / len(events)

tau1, nu1 = molchan(lam + rng.normal(0, lam.std() * 0.5, n_cells))  # 有技能
tau2, nu2 = molchan(rng.random(n_cells))                            # 亂猜

fig = go.Figure()
fig.add_trace(go.Scatter(x=tau1, y=nu1, mode="lines", name="有技能的模型",
                         line=dict(color=ACCENT, width=2.5)))
fig.add_trace(go.Scatter(x=tau2, y=nu2, mode="lines", name="無資訊模型",
                         line=dict(color=PALETTE[3], width=1.5)))
fig.add_trace(go.Scatter(x=[0, 1], y=[1, 0], mode="lines", name="隨機猜測參考線",
                         line=dict(color="#888", dash="dash")))
apply_layout(fig, title="Molchan 圖：用多小的警報範圍，抓到多高比例的地震？",
             xaxis_title="警報覆蓋的時空比例 τ", yaxis_title="漏報率 ν",
             hovermode="closest")
fig

# %% [markdown]
# 好模型的軌跡壓向左下角：用兩成的警報面積就抓到約三分之二的
# 地震（整條軌跡的面積技能分數約 0.76，遠高於隨機的 0.5）；無資訊模型貼著對角線——警報覆蓋多少比例，就只抓到
# 多少比例，與亂猜無異。把軌跡上方的面積積分成單一分數（area
# skill score，隨機參考為 0.5），就能一次涵蓋所有門檻。警報式
# 語言對決策情境特別有用——「要不要提高警戒」本質上就是一個
# 二值化問題——而且 Molchan 圖可以用長期活動度加權警報體積，
# 讓「只在本來就常地震的地方發警報」不再佔便宜。
#
# ## 15.4 檢驗的檢驗：統計功效
#
# 到這裡你可能以為：模型通過全部檢驗＝好模型。這一節要打破
# 這個想法，用的是近年最發人深省的一個實驗（Khawaja et al.
# 2023）：一個宣稱「地球上每個地方地震機率都一樣」的**空間
# 均勻模型**，在標準的 0.1° 全球網格上，**通過了 S-test**。
#
# 怎麼會？因為全球網格有六百多萬格，測試期只有幾百顆目標
# 地震——平均一萬格才攤到一顆。格子細到每顆地震各佔一格時，
# 似然分數只由「有幾格中了一顆」決定，「中在哪裡」的資訊被
# 解析度稀釋殆盡，檢驗就失去了拒絕爛模型的力氣。這個「力氣」
# 有正式名字：**統計功效（statistical power）**——正確拒絕
# 錯誤模型的機率。作者估計，要在慣用的全球高解析度網格上
# 得到有力的 S-test，需要約三萬顆地震，相當於**三百年**的
# 觀測。解方不是等三百年，而是改變設計：用資料驅動的多解析度
# 網格（地震多的地方細、少的地方粗），**八顆**地震就能達到
# 最大功效。網格是研究者自己選的自由變數——檢驗設計本身
# 就是實驗設計。
#
# 由此得到本章最重要的三句警語：
#
# 1. **「沒被拒絕」不等於「模型好」**，往往只代表檢驗沒有力氣
#    拒絕。嚴謹的說法是「在此網格與樣本數下，資料不足以拒絕
#    該模型」，並附上功效評估。回想{doc}`第 8 章 <08_explore_ideas>`
#    的誤報率精神：檢驗只有在「拒絕」時才攜帶強資訊。
# 2. **少數大震主宰似然。**一個序列擠進幾格，就能翻轉整份
#    排名——義大利的模型排名曾僅因測試期是否涵蓋一個序列
#    而整個改寫。樣本數是地震顆數，不是格子數（第 14 章）。
# 3. **回溯成績永遠只是 sanity check。**上一章的十年開獎已經
#    示範過：同一批模型、同一統計量，回溯 +0.5、前瞻 −0.7。
#    沒有前瞻測試，就沒有結論。
#
# 再補三個實作陷阱：任何含對數的分數，只要一顆地震落在預報率
# 為零的箱，分數就是 $-\infty$——爛掉的是整份預報，實務上要
# 給零率箱設底線率；跨模型比較前要對齊規模尺度（$M_w$ vs
# $M_L$）、深度範圍與除叢方式，不對齊就沒有可比性；資料驅動
# 的網格會繼承歷史盲點——紐西蘭 Canterbury 序列的起點，恰好
# 落在按歷史資料建格的低解析度區。
#
# ## 15.5 pyCSEP：把檢驗變成人人可用的工具
#
# 以上所有工具——四種一致性檢驗、T/W-test、負二項與二元
# 變體、Molchan 與 area skill score、多解析度網格——都已被
# 實作在開源套件 **pyCSEP** 裡。這件事的意義超過方便：CSEP
# 第一階段的檢驗程式封在測試中心裡，外人無法檢視；把它重寫
# 成開源套件，等於把「裁判的尺」也攤在陽光下。現代的標竿
# 做法更進一步：論文附上**可重現性套件**——程式碼、資料、
# 凍結的執行環境，一道指令重跑全部圖表。可重現性不是附錄，
# 是方法論的一部分。
#
# 對台灣，這一章的接點很具體：台灣測試區小、目標地震少，
# 正是統計功效警告的低功效情境——多解析度網格幾乎是必需品；
# 台灣序列叢集極強（車籠埔、池上、大埔），負二項與二元似然
# 也是；而任何台灣模型宣稱有價值之前，都有一個現成的全球
# 基準模型（GEAR1）可以投影過來當對照組——順帶一提，在
# 加州、紐西蘭、義大利這三個全世界研究最透徹的地區，區域
# 模型都沒能穩定打敗這個全球模型。「我們的模型是為本地量身
# 打造的，所以一定比較準」——這句話，現在你知道該怎麼
# 檢驗了。
#
# {doc}`下一章 <16_psha>`我們離開「會不會有地震」，走向工程師
# 真正的問題：「地表會晃多大？」——從預報到危害度，PSHA。
