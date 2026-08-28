# %% [markdown]
# # 15. Ψ：前兆尺度增加的解剖
#
# {doc}`第 14 章 <14_etas_estimation>`結束時留下一句話：ETAS 描述的是
# **觸發**——這個地震會引發什麼；但地震序列裡還有另一類現象，時間之箭
# 是反過來的——一群中小地震在數月到數十年之後，被一場更大的地震跟上。
# 那不是觸發，是**預示**。
#
# 這一章就專門解剖那個反向的箭頭。它有名字：**Ψ（前兆尺度增加，
# precursory scale increase）**——大地震發生前，其震源區內中小地震的
# 「規模水準」與「發生率」會**同時**出現一次階梯式的抬升。這個現象由
# Evison 與 Rhoades 在 1990–2000 年代從四個地區的 47 個大地震歸納出來，
# 後來成為 EEPAS 預報模型的全部地基（{doc}`第 16 章 <16_eepas_ppe>`）。
#
# 但本章不是「介紹一個前兆」。Ψ 是這本教材裡最好的**方法論標本**：
# 它的原始辨識程序是手工的、事後的、以「最大化尺度增加」為目標函數的
# ——第一部{doc}`第 8 章 <08_explore_ideas>`警告過的事後選擇偏誤全套
# 到齊；而二十年後，同一群作者把程序寫成演算法、加上隨機化對照、
# 並公開承認前輩（也是自己）的做法「可能被視為在滿足一個預設好的
# 結論」（Christophersen et al. 2024）。**偏誤與矯正發生在同一條研究
# 線上，而且都被寫進論文。** 這種標本不多見。
#
# 本章的順序是：先把現象量出來（$C(t)$ 曲線與四個 Ψ 變數），再檢查
# 量法本身有幾個自由度，然後看兩套自動化演算法怎麼把自由度攤在陽光
# 下，接著用隨機化對照問「訊號是不是真的」，最後處理三個統計陷阱
# ——不唯一性、時空取捨與回歸稀釋、Simpson 反轉——並把七條尺度關係
# 的重估結果攤開。
#
# 本章要用到的前置：第 11 章的目錄完整度 $M_c$ 與 GR 律、第 12 章的
# 餘震衰減與除叢、第 14 章的「什麼叫正常叢集」。Ψ 之所以能被辨識成
# 一種**不同於餘震衰減**的叢集，前提正是我們已經有了一個描述正常
# 叢集的參考模型。
#
# 本章擁有三組式子，其他章一律引用不重推：累積規模異常 $C(t)$
# {eq}`eq:cumag`、Ψ 的三條尺度迴歸 {eq}`eq:psi-scaling`，以及時空
# 取捨與回歸稀釋（15.7 節）。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% [markdown]
# ## 15.1 現象與量法：累積規模異常
#
# ### 為什麼需要一條曲線
#
# 「大地震前中小地震變多變大」這句話，光看規模–時間散布圖是量不出來
# 的。人眼很容易在任何一段隨機序列裡看出「這裡比較密」。要把它變成
# 可計算的量，需要一個**把規模與發生率同時吃進去、並且自帶零基準**
# 的統計量。Evison 與 Rhoades 用的就是**累積規模異常**（cumulative
# magnitude anomaly，文獻慣稱 cumag）。
#
# 設在時窗 $[t_s, t_f]$ 內，取一個高於目錄完整度的下限規模 $m_c$
# （注意：$m_c$ 是**分析時自選的截切門檻**，與第 11 章那個目錄性質
# $M_c$ 不是同一個東西，要求只是 $m_c \ge M_c$）。定義
#
# $$C(t) = \sum_{t_s < t_i \le t}\bigl(M_i - m_c - 0.1\bigr)
#   \;-\; k\,(t - t_s),
#   \qquad
#   k = \frac{\displaystyle\sum_{t_s < t_i \le t_f}\bigl(M_i - m_c - 0.1\bigr)}
#            {t_f - t_s}$$ (eq:cumag)
#
# 三個元件各有職責。**$M_i - m_c - 0.1$** 是第 $i$ 個地震的「規模超額」
# ——目錄規模一般以 0.1 為級距，$m_c$ 這一格的事件中心值是 $m_c$，
# 減去 $m_c + 0.1$ 讓最小的合格事件貢獻約為零，於是曲線的上升幾乎
# 完全來自「比門檻大的部分」。**求和項**是規模超額的累積，單調不減。
# **$k(t-t_s)$** 是整段時窗的**平均規模累積速率**，單位是每年多少個
# 規模單位（magnitude units per year，M.U. yr$^{-1}$）；$m_c = 5.0$ 時
# 一年發生一個 $M\,5.9$ 就相當於 $1$ M.U. yr$^{-1}$。
#
# ### 三個性質的證明
#
# {eq}`eq:cumag` 有三個一眼看不出、但一算就清楚的性質，Ψ 的整套量法
# 都建立在它們之上。
#
# **性質一：$C(t_s) = C(t_f) = 0$。** 起點端，求和的範圍是空集合，
# 第一項為 0，第二項 $k(t_s - t_s) = 0$，故 $C(t_s) = 0$。終點端把
# $k$ 的定義代回去：
#
# $$\begin{aligned}
# C(t_f) &= \sum_{t_s < t_i \le t_f}(M_i - m_c - 0.1)
#   - \frac{\sum_{t_s < t_i \le t_f}(M_i - m_c - 0.1)}{t_f - t_s}\,(t_f - t_s)\\
#   &= \sum_{t_s < t_i \le t_f}(M_i - m_c - 0.1)
#      - \sum_{t_s < t_i \le t_f}(M_i - m_c - 0.1) \;=\; 0 .
# \end{aligned}$$
#
# 換句話說，**$k$ 這個扣除項是被「讓曲線頭尾接回零」這個要求唯一
# 決定的**。$C(t)$ 因此不是一條自由的累積曲線，而是一條被兩端釘死的
# 橋——它只能量「相對於整段平均趨勢的偏離」，量不到絕對水準。
#
# **性質二：每個地震處往上跳一格。** 在沒有事件的區間，
# $\mathrm{d}C/\mathrm{d}t = -k < 0$（只要時窗內至少有一個超過門檻的
# 事件，$k > 0$），曲線以固定斜率下滑；在 $t = t_i$ 處，求和項瞬間
# 增加 $M_i - m_c - 0.1 \ge 0$，曲線垂直往上跳這麼多。**跳的高度由
# 規模決定，兩跳之間下滑的長度由等待時間決定**——一條曲線同時編碼了
# 規模與發生率，這正是我們要的。
#
# **性質三：最小值即 onset。** 把 $C$ 寫成
# $C(t) = \bigl[\bar\rho(t_s,t) - k\bigr](t - t_s)$，其中
# $\bar\rho(t_s,t)$ 是 $[t_s,t]$ 內的平均規模累積速率。$C$ 在 $t^*$
# 取到最小值，代表 $[t_s, t^*]$ 這一段的平均速率相對於全段最低、而
# $[t^*, t_f]$ 那一段最高——$t^*$ 就是「速率換檔」的位置。定義
# **$t^*$ 為 Ψ 的起始點（onset）**。
#
# ### 從一次辨識量出四個變數
#
# 有了 onset，四個 Ψ 變數就都能讀出來（記主震時間為 $t_M$、規模為
# $M_m$）：
#
# - **前兆時間** $T_P = t_M - t^*$，慣以「天」為單位。
# - **前兆規模** $M_P$：前兆期 $(t^*, t_M)$ 內**最大三個地震規模的
#   平均**。取三個而不是一個，是為了降低單一事件的隨機性。
# - **先前規模水準** $M^-$：先前期 $(t_s, t^*)$ 內最大三個規模的平均。
# - **前兆面積** $A_P$：包住前兆（原始程序還要求包住主震與餘震）的
#   最小矩形面積，慣以 km$^2$ 為單位。
#
# 尺度增加的**強度**有兩種量法，兩者互補：規模式的 $M_P - M^-$，以及
# 純計數式的**率比** $r$——前兆期平均發生率對先前期平均發生率的比值。
# 先用合成資料把整張圖畫出來（真實範例見 Christophersen et al. 2024
# 的 Figure 1，1999 年加州 Hector Mine $M\,7.1$）：

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

rng = np.random.default_rng(15001)
MC = 4.0                                  # 分析用下限規模 m_c
T_S, T_M, T_F = 0.0, 20.0, 24.0           # 起點、主震、終點（年）
ONSET_TRUE = 13.0

t_prior = np.sort(rng.uniform(T_S, ONSET_TRUE, 22))          # 先前期
m_prior = MC + rng.exponential(0.30, t_prior.size)
t_pre = np.sort(rng.uniform(ONSET_TRUE, T_M, 27))            # 前兆期
m_pre = MC + 0.45 + rng.exponential(0.42, t_pre.size)
t_aft = T_M + np.sort(rng.exponential(0.55, 24))             # 餘震
m_aft = MC + rng.exponential(0.38, t_aft.size)
t_aft, m_aft = t_aft[t_aft < T_F], m_aft[t_aft < T_F]


def cumag(t, m, t0, t1, mc=MC):
    """回傳 (時間格點, C(t))；C 在每個事件前後各取一點以顯示跳躍。"""
    sel = (t > t0) & (t <= t1)
    te, me = t[sel], m[sel]
    exc = me - mc - 0.1
    k = exc.sum() / (t1 - t0)
    grid = np.concatenate([[t0], np.repeat(te, 2), [t1]])
    cum = np.concatenate([[0.0], np.repeat(np.cumsum(exc), 2)[:-1], [exc.sum()]])
    cum = np.concatenate([[0.0], np.repeat(np.concatenate([[0.0],
                          np.cumsum(exc)[:-1]]), 1)])
    # 重新以「事件前 / 事件後」兩點描出階梯
    grid = np.concatenate([[t0], np.repeat(te, 2), [t1]])
    before = np.concatenate([[0.0], np.cumsum(exc)[:-1]])
    cum = np.concatenate([[0.0], np.ravel(np.column_stack([before,
                          np.cumsum(exc)])), [exc.sum()]])
    return grid, cum - k * (grid - t0), k


t_all = np.concatenate([t_prior, t_pre])
m_all = np.concatenate([m_prior, m_pre])
g1, C1, k1 = cumag(t_all, m_all, T_S, T_M)                   # [t_s, t_M]
g2, C2, _ = cumag(t_aft, m_aft, T_M, T_F)                    # [t_M, t_f]
onset = g1[int(np.argmin(C1))]

pre_mask, pri_mask = t_all > onset, t_all <= onset
M_P = np.sort(m_all[pre_mask])[-3:].mean()
M_MINUS = np.sort(m_all[pri_mask])[-3:].mean()
r_ratio = (pre_mask.sum() / (T_M - onset)) / (pri_mask.sum() / (onset - T_S))
T_P_yr = T_M - onset

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.52, 0.48],
                    vertical_spacing=0.06)
fig.add_trace(go.Scatter(x=t_all[pri_mask], y=m_all[pri_mask], mode="markers",
                         name="先前期事件",
                         marker=dict(size=7, color=PALETTE[4], opacity=0.85)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=t_all[pre_mask], y=m_all[pre_mask], mode="markers",
                         name="前兆期事件",
                         marker=dict(size=8, color=ACCENT, opacity=0.9)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=t_aft, y=m_aft, mode="markers", name="餘震",
                         marker=dict(size=6, color=PALETTE[5], opacity=0.55)),
              row=1, col=1)
fig.add_trace(go.Scatter(x=[T_M], y=[7.0], mode="markers", name="主震",
                         marker=dict(size=15, color=QUAKE_COLOR, symbol="star")),
              row=1, col=1)
for lvl, x0, x1, nm in ((M_MINUS, T_S, onset, "M⁻"), (M_P, onset, T_M, "M_P")):
    fig.add_trace(go.Scatter(x=[x0, x1], y=[lvl, lvl], mode="lines",
                             name=nm, line=dict(color="#666666", dash="dash",
                                                width=1.5)), row=1, col=1)
fig.add_trace(go.Scatter(x=g1, y=C1, mode="lines", name="C(t)：主震前",
                         line=dict(color=PALETTE[3], width=2)), row=2, col=1)
fig.add_trace(go.Scatter(x=g2, y=C2, mode="lines", name="C(t)：主震後",
                         line=dict(color=PALETTE[3], width=1.5, dash="dot")),
              row=2, col=1)
for r_ in (1, 2):
    fig.add_vline(x=onset, line_dash="dot", line_color="#1baf7a", row=r_, col=1)
fig.add_annotation(x=onset, y=float(C1.min()),
                   text=f"onset（T_P = {T_P_yr:.1f} 年）", showarrow=True,
                   ax=-55, ay=32, row=2, col=1)
fig.update_yaxes(title_text="規模", row=1, col=1)
fig.update_yaxes(title_text="C(t)（M.U.）", row=2, col=1)
fig.update_xaxes(title_text="時間（年）", row=2, col=1)
apply_layout(fig, height=560, hovermode="x",
             title=(f"Ψ 的量法（合成資料）：M_P = {M_P:.2f}、M⁻ = {M_MINUS:.2f}、"
                    f"M_P − M⁻ = {M_P - M_MINUS:.2f}、率比 r = {r_ratio:.1f}"))
fig

# %% [markdown]
# 上圖三件事值得停下來看。第一，$C(t)$ **確實從零出發、回到零**——
# 這不是巧合，是 $k$ 的定義逼出來的（性質一）。第二，最小值出現在
# 「下滑段結束、上升段開始」的轉折，而不是在任何一個大事件上；
# **onset 是速率換檔點，不是最大事件的時刻**。第三，主震之後另有
# 一條 $C(t)$ 曲線（虛線）——原始程序在 $[t_s, t_M]$ 與 $[t_M, t_f]$
# 各畫一條，後者用來檢查餘震期的行為，不參與 onset 的判定。
#
# 順帶一提斜率的讀法：$C$ 曲線任一段的斜率加上 $k$，就是那一段的
# 平均規模累積速率。原論文的圖上畫了一個「量角器」，把斜率直接翻譯
# 成 M.U. yr$^{-1}$——這是一個很好的視覺化設計，值得學起來。

# %% [markdown]
# ## 15.2 手工辨識與 hindsight 偏誤
#
# ### demarcation：一個有三個自由度的程序
#
# Evison 與 Rhoades 把「為某個主震辨識 Ψ」的程序稱作 **demarcation
# （劃界）**。實際操作是這樣的：先設一個初始下限規模 $m_c$（通常比
# 主震規模低 2–3 個單位、但仍高於目錄完整度），再框一個南北–東西向
# 的矩形，涵蓋主震與餘震震央；然後**逐步調整矩形的位置、形狀與大小，
# 使尺度增加最大化**。作者還記下一條經驗：「先前期的長度與前兆時間
# 相當時，尺度增加顯示得最清楚。」
#
# 把這段話翻成統計語言，問題就無所遁形了：
#
# - **三個自由度**（位置、形狀、大小）**被連續調整**；
# - **目標函數是「最大化尺度增加」**——也就是最大化我們想看到的那個
#   訊號本身；
# - **調整時已知答案**（主震的位置、時間、規模都在手上）；
# - 還有第四層：**內建的尺度直覺**。作者自承，由於他們在找的是一個
#   「長期」孕震過程，**短前兆時間根本不被考慮**——$M\,5$–$6$ 的主震
#   或許可以接受一年的前兆時間，$M\,7$–$8$ 則「絕不可能」。
#
# 這正是第 8 章講的**事後選擇偏誤**的教科書形態：搜尋空間很大、選擇
# 準則在看過資料之後才定、而且沒有人記錄「試了幾個框才成功」。
#
# ### 原作者的自我檢舉
#
# 最誠實、也最值得整段引用的評語，來自二十年後的同一群人。
# Christophersen, Rhoades & Hainzl（2024, SRL）在討論演算法版與人工版
# 的差異時寫道：由於演算法完全不預設任何尺度觀念、允許前兆時間遠比
# 人工版短、前兆面積遠比人工版大或小，因此——**在 Evison 與 Rhoades
# 的研究中，那些尺度關係「可能會被視為在滿足一個預設好的結論」
# （might have been regarded as fulfilling a preconceived outcome），
# 而對本文的演算法識別則不能這樣說**。
#
# 這句話的份量在於：它不是外部批評者寫的，而是原始尺度關係的共同作者
# 自己寫的，而且寫在同一條研究線的延續論文裡。
#
# ### 演算法版的價值不在「更準」
#
# 一個常見的誤讀是「自動化比較客觀，所以結果比較可信」。不對。
# 15.3 節會看到，演算法本身也充滿人訂的門檻（$r > 3$、$0.4$、
# 中央 20%…），換一組門檻就換一組結果，作者自己也說「這些演算法不是
# 唯一的，也不是不可爭議的」。**自動化真正買到的是三樣別的東西**：
#
# 1. **可重現**——判準寫成程式碼，別人能重跑、能質疑、能改門檻重跑。
# 2. **允許失敗**——13 個主震找不到任何合格的 Ψ，這個數字被寫進論文，
#    而不是被靜靜移出樣本（15.6 節）。
# 3. **允許醜陋的答案**——不預設「$M\,7$ 不可能只有一年前兆時間」，
#    於是出現了前兆時間極短、甚至前兆面積小於震源面積的識別。
#
# 第三點特別重要。人工程序的「尺度直覺」等於在資料上先疊了一層結論；
# 演算法拿掉這層，尺度關係**仍然出現**——這才是真正的檢驗。

# %% [markdown]
# ## 15.3 自動化 I：rectangular 演算法
#
# ### 輸入與流程
#
# rectangular 演算法沿用人工程序的矩形框，但把「調整」換成**網格
# 搜尋**，並把所有通過門檻的解**全部保留**（而不是停在第一個好看的
# 解）。輸入是：一份主震清單（時間 $t_M$、位置、規模 $M_m$）、對應的
# 目錄、深度範圍、最大搜尋距離 $R$、回溯時間 $T$、下限規模 $m_c$。
# 流程可以壓成八步：
#
# 1. 若 $m_c < M_m - 2.5$，把它拉高到 $M_m - 2.5$（避免門檻離主震
#    規模太遠，前兆群被小事件淹沒）。
# 2. 選出主震前 $T$ 年內、南北向與東西向距離**皆**小於 $R$ 的事件，
#    構成初選集 $S$（這是一個以主震為中心的正方形，不是圓）。
# 3. 令 $R_{\max}$ 等於 $S$ 中任一事件到主震的最大南北或東西距離加上
#    $\Delta = 0.5$ km。
# 4. 若 $S$ 至少有 10 個事件、且 $S$ 內最大規模 $\le M_m - 0.5$，
#    就算 $C(t)$ 的最小值時間 $T_{\min}$（即 onset）。
# 5. 若 onset 之後至少有 3 個前兆，**把矩形縮到剛好包住這些前兆**
#    （南北、東西方向各自縮，下限為 $R_{\max}/2$），得到子集。
# 6. 逐步提高 $m_c$，直到再提高會使 onset 之前的事件少於 3 個；若總
#    事件數仍 $\ge 10$，計算 $T_{\min}$、$M_P$、$M^-$ 與率比 $r$，
#    並施加接受條件。
# 7. $R_{\max} \leftarrow 0.95\,R_{\max}$，重複 4–6，最多 80 次。
# 8. $T \leftarrow 0.95\,T$，重複 2–7，最多 240 次。
#
# ### 接受條件與去重規則
#
# 一個候選解只有全部通過下表才會被存下來。這張表就是本章的參數表
# ——**Ψ 的「有沒有」完全取決於這幾個數字**：
#
# | 判準 | 值 | 為什麼 |
# |---|---|---|
# | 事件數 | $\ge 10$ | 太少的話 $C(t)$ 只是雜訊 |
# | 窗內最大規模 | $\le M_m - 0.5$ | 排除窗內另有大事件的情形 |
# | onset 位置 | 落在時間範圍中央 20% | 先前期與前兆期長度相當 |
# | 尺度增加（規模） | $M_m - M_P > 0.4$ 且 $M_P - M^- > 0.4$ | 前兆要比主震小、比背景大 |
# | 尺度增加（率） | $r > 3$ | 發生率至少抬升三倍 |
#
# 另外兩條幾何與去重規則：矩形的**長寬比限制在 0.5–2.0**（避免退化
# 成細長條；實際上多數解是正方形）；去重時，先在「起始時間、$T_P$、
# $M_P$ 皆相同」的解裡留下 $r$ 最大者，再在剩下「$T_P$ 與 $r$ 皆
# 相同」的解裡留下 $A_P$ 最小者。
#
# ### 兩點與人工程序的差異
#
# 這裡最容易被跳過、卻最該講清楚的是：**演算法不是人工程序的忠實
# 複製，它在空間上一邊更嚴、一邊更鬆**。
#
# **更嚴的一邊**：演算法**強制把主震震央放在前兆區的正中心**（步驟 2
# 的正方形以主震為心）。人工程序允許主震落在前兆區內的任何位置。
#
# **更鬆的一邊**：演算法**完全不管餘震落在哪裡**。人工程序要求前兆區
# 必須包住餘震震央，用意是讓前兆面積至少不小於主震震源面積。拿掉這
# 條之後，**演算法原則上允許前兆面積比震源還小**——這種在物理直覺上
# 「不好看」的解，正是 15.2 節說的「允許醜陋答案」。
#
# 最後一句作者的自白值得抄下來：**這個演算法不做全域最佳化**。它是在
# 一個給定的大時空體積裡，找出一個「看起來有 Ψ」的子集；改變體積的
# 大小，就得到不同的識別。**多重識別不是 bug，是設計的直接後果**
# ——而它恰好讓 15.6、15.7 兩節的研究成為可能。

# %% [markdown]
# ## 15.4 自動化 II：circular 演算法與 Z 值
#
# ### 一個完全不含規模資訊的強度量測
#
# 第二套演算法換了三件事：區域改用**以主震為心的圓**、在（半徑 $R$、
# 回溯時間 $T$）上做對數網格搜尋、而且**每個主震只保留一個最顯著的
# 解**。「最顯著」由 Habermann（1981）的 **$Z$ 值**定義。
#
# $Z$ 要回答的問題是：兩段時間內觀測到的地震數，能不能用「同一個
# 泊松率」解釋？設第 1 段（先前期）長度 $\Delta T_1$、事件數 $N_1$，
# 第 2 段（前兆期）長度 $\Delta T_2$、事件數 $N_2$。速率的自然估計是
# $\hat\rho_i = N_i / \Delta T_i$；在泊松假設下 $\mathrm{Var}(N_i)
# = \mathbb{E}[N_i] \approx N_i$，故
#
# $$\mathrm{Var}(\hat\rho_i) \approx \frac{N_i}{\Delta T_i^2} .$$
#
# 兩段獨立，把差除以標準差：
#
# $$\begin{aligned}
# \frac{\hat\rho_2 - \hat\rho_1}
#      {\sqrt{\mathrm{Var}(\hat\rho_1) + \mathrm{Var}(\hat\rho_2)}}
# &= \frac{\dfrac{N_2}{\Delta T_2} - \dfrac{N_1}{\Delta T_1}}
#         {\sqrt{\dfrac{N_1}{\Delta T_1^2} + \dfrac{N_2}{\Delta T_2^2}}}
#  = \frac{\dfrac{N_2\Delta T_1 - N_1 \Delta T_2}{\Delta T_1 \Delta T_2}}
#         {\dfrac{1}{\Delta T_1 \Delta T_2}
#          \sqrt{N_1 \Delta T_2^2 + N_2 \Delta T_1^2}} \\[4pt]
# &= \frac{N_2 \Delta T_1 - N_1 \Delta T_2}
#         {\sqrt{N_2 \Delta T_1^2 + N_1 \Delta T_2^2}} \;\equiv\; Z .
# \end{aligned}$$
#
# （分子分母同乘 $\Delta T_1 \Delta T_2$ 即得。）在 $N_1, N_2$ 都大的
# 極限下，若兩段真的來自同一個隨機（泊松）過程，$Z$ 服從標準差為 1
# 的高斯分布——所以 $Z$ 可以直接當「幾個標準差」讀。$Z > 0$ 代表
# 第 2 段的速率高於第 1 段。
#
# ### 為什麼 $Z$ 是關鍵的對照
#
# 請注意 $Z$ 的公式裡**只有事件數與時間長度，一個規模都沒有**。
# 這一點極其重要，值得單獨成段：
#
# **在 $Z$ 之前，所有 Ψ 強度的量測都吃規模資訊**——$C(t)$ 的每一跳
# 是規模超額，$M_P - M^-$ 直接是規模差，連率比 $r$ 都依賴 $m_c$ 這個
# 規模門檻決定誰被算進去。因此 rectangular 演算法（規模式強度）與
# circular 演算法（純計數強度）是**兩條在設計上互相獨立的量測路徑**。
# 兩者得到彼此一致的尺度關係，比同一條路徑上做兩次實驗有力得多。
#
# 用 15.1 節那份合成資料把 $Z$ 的兩時窗定義畫出來：

# %% tags=["hide-input"]
dT1, dT2 = onset - T_S, T_M - onset
N1, N2 = int(pri_mask.sum()), int(pre_mask.sum())
Z_demo = (N2 * dT1 - N1 * dT2) / np.sqrt(N2 * dT1 ** 2 + N1 * dT2 ** 2)

fig = go.Figure()
fig.add_vrect(x0=T_S, x1=onset, fillcolor=PALETTE[4], opacity=0.13,
              line_width=0, layer="below")
fig.add_vrect(x0=onset, x1=T_M, fillcolor=ACCENT, opacity=0.13,
              line_width=0, layer="below")
for tt in t_all[pri_mask]:
    fig.add_shape(type="line", x0=tt, x1=tt, y0=0, y1=1,
                  line=dict(color=PALETTE[4], width=1.6))
for tt in t_all[pre_mask]:
    fig.add_shape(type="line", x0=tt, x1=tt, y0=0, y1=1,
                  line=dict(color=ACCENT, width=1.8))
fig.add_trace(go.Scatter(x=[T_M], y=[1.15], mode="markers", name="主震",
                         marker=dict(size=15, color=QUAKE_COLOR, symbol="star")))
fig.add_annotation(x=(T_S + onset) / 2, y=1.55, showarrow=False,
                   text=f"先前期：ΔT₁ = {dT1:.1f} 年，N₁ = {N1}<br>"
                        f"速率 = {N1 / dT1:.2f} /年")
fig.add_annotation(x=(onset + T_M) / 2, y=1.55, showarrow=False,
                   text=f"前兆期：ΔT₂ = {dT2:.1f} 年，N₂ = {N2}<br>"
                        f"速率 = {N2 / dT2:.2f} /年")
fig.add_vline(x=onset, line_dash="dot", line_color="#1baf7a")
fig.add_annotation(x=onset, y=-0.28, text="onset", showarrow=False)
fig.update_yaxes(range=[-0.45, 1.95], showticklabels=False, title_text="")
apply_layout(fig, height=340, hovermode="closest",
             xaxis_title="時間（年）", showlegend=False,
             title=(f"Z 值的兩時窗定義：Z = {Z_demo:.2f}"
                    f"（≈ {Z_demo:.1f} 個標準差；公式完全不含規模）"))
fig

# %% [markdown]
# ### circular 演算法的搜尋設定
#
# 套用在 RSQSim 模擬目錄上時，設定是：跳過最前面 5000 年（暖機期），
# 取所有 $M \ge 7.0$ 的事件為候選主震；對每個主震，取距離 $r \le R$、
# 規模 $m \ge 5.0$、時間落在 $[t_M - T, t_M)$ 的事件；掃描範圍
# $1 \le R \le 60$ km（60 個對數空間格）與 $1 \le T \le 1000$ 年
# （50 個對數時間格）。**含有「與主震規模相差半個單位以內」的事件的
# 時空窗一律排除**——因為 Ψ 的既有識別都有「前兆規模遠小於主震規模」
# 這個特徵，不排除就會讓另一個大事件冒充前兆。最後保留 $Z$ 最大的
# 那一個解。

# %% [markdown]
# ## 15.5 訊號是真的嗎：三組隨機化對照
#
# ### 先把措辭釘死
#
# 這一節是全章最有力、也最容易被引用錯的一段。**三組對照都是在
# RSQSim 物理模擬目錄上做的，不是在真實目錄上做的。** RSQSim
# （rate-and-state earthquake simulator）是一個只放進速率–狀態摩擦
# 定律、斷層應力交互作用與破裂準則的物理模擬器；作者用中紐西蘭地殼
# 斷層網跑了 **20,000 年**，得到 **495,626 個地震**。真實目錄（四個
# 地區、47 個主震）太短、事件太少，無法支撐這種等級的隨機化實驗——
# 這正是要用模擬器的理由，也是引用時必須帶上的限定條件。
#
# ### 三組對照與各自要排除的解釋
#
# 對這份模擬目錄跑 circular 演算法，比較四種版本的最大 $Z$ 值累積
# 分布：
#
# 1. **原始目錄** —— 基準，$Z$ 最大。
# 2. **時間隨機化**（把所有事件的發生時間重新洗牌，保留位置與規模）
#    —— $Z$ 最小、**訊號消失**。要排除的解釋是：「$Z$ 只是有限樣本的
#    隨機起伏」。既然打散時間後訊號就沒了，$Z$ 的高值必定來自時間
#    結構。
# 3. **規模隨機化**（把規模重新洗牌，保留時間與位置）—— $Z$ **減弱
#    但仍接近原始**，且原始目錄平均而言仍顯著較高。要排除的解釋是：
#    「訊號只是規模分布的假象」。結果顯示訊號**主要**來自時空叢集
#    （洗掉規模後叢集仍在，所以 $Z$ 仍高），規模結構只提供額外的一小
#    部分。
# 4. **去餘震**（用 Gardner & Knopoff 1974 的視窗法除叢，事件數從
#    481,948 降到 153,409，砍掉超過三分之二）—— 各曲線的 $Z$ 都變小
#    （窗內事件數少了，統計顯著性自然下降），但**四條曲線的相對位置
#    不變**，原始與規模隨機化版本仍顯著。要排除的解釋是：「Ψ 只是餘震
#    序列的偽影」——例如前一次大地震的餘震尾巴，剛好落在後一次大地震
#    的前兆窗裡。
#
# ### 結論的正確講法
#
# 把三組合起來，能說的**恰好只有這一句**：
#
# > 在物理模擬目錄中，大地震之前存在一種**不同於餘震衰減的時空叢集**；
# > 它依賴事件的時間結構（時間隨機化後消失），不完全依賴規模結構
# > （規模隨機化後仍在），也不是餘震的偽影（去餘震後仍在）。
#
# 不能說的有三句：不能說「Ψ 被證明是真的前兆」（模擬目錄不是真實
# 地球）；不能說「Ψ 可以用來預報」（辨識需要事先知道主震在哪，
# 15.6 節會看到連辨識本身都不唯一）；也不能說「真實目錄通過了同樣的
# 對照」（沒做，也做不了）。
#
# 下面用一份**合成**目錄把這個實驗的骨架跑一次——不是重製論文結果，
# 而是示範這套對照怎麼在課堂上自己動手做。合成目錄裡植入了 40 個
# 主震，每個主震前有一群「規模與發生率同時抬升」的前兆，之後有一串
# 帶長尾的餘震；我們對每個主震在 $(R, T)$ 網格上掃描、每個窗用
# $C(t)$ 定 onset、算 $Z$，取最大值：

# %% tags=["hide-input"]
RNG3 = np.random.default_rng(15003)
MC_SCAN = 3.5                      # 掃描用下限規模（高於合成目錄完整度 3.0）
M_MIN, BETA = 3.0, np.log(10)      # 合成目錄的完整度與 GR 斜率（b = 1）
T_TOTAL, L_BOX, N_MS = 600.0, 400.0, 40      # 年、km、主震數

bg_t = RNG3.uniform(0, T_TOTAL, 6000)
bg_x, bg_y = RNG3.uniform(0, L_BOX, 6000), RNG3.uniform(0, L_BOX, 6000)
bg_m = M_MIN + RNG3.exponential(1 / BETA, 6000)

zx, zy = RNG3.uniform(70, L_BOX - 70, 5), RNG3.uniform(70, L_BOX - 70, 5)
zi = RNG3.integers(0, 5, N_MS)
ms_x, ms_y = zx[zi] + RNG3.normal(0, 25, N_MS), zy[zi] + RNG3.normal(0, 25, N_MS)
ms_t = np.sort(RNG3.uniform(150, T_TOTAL - 5, N_MS))

pt, px, py, pm, at, ax_, ay_, am = [], [], [], [], [], [], [], []
for k in range(N_MS):
    Tp = 10 ** RNG3.normal(1.05, 0.22)                 # 前兆時間（年）
    n_p = int(RNG3.integers(18, 34))
    pt.append(ms_t[k] - RNG3.uniform(0, Tp, n_p))
    px.append(ms_x[k] + RNG3.normal(0, 18, n_p))
    py.append(ms_y[k] + RNG3.normal(0, 18, n_p))
    pm.append(M_MIN + 0.55 + RNG3.exponential(1 / BETA, n_p))
    n_a = int(RNG3.integers(50, 90))
    at.append(ms_t[k] + 0.05 * ((1 - RNG3.random(n_a)) ** (-2.5) - 1))
    ax_.append(ms_x[k] + RNG3.normal(0, 12, n_a))
    ay_.append(ms_y[k] + RNG3.normal(0, 12, n_a))
    am.append(M_MIN + RNG3.exponential(1 / BETA, n_a))

cat_t = np.concatenate([bg_t] + pt + at)
cat_x = np.concatenate([bg_x] + px + ax_)
cat_y = np.concatenate([bg_y] + py + ay_)
cat_m = np.concatenate([bg_m] + pm + am)
is_aft = np.concatenate([np.zeros(bg_t.size, bool)]
                        + [np.zeros(a.size, bool) for a in pt]
                        + [np.ones(a.size, bool) for a in at])
keep_t = cat_t < T_TOTAL
cat_t, cat_x, cat_y, cat_m, is_aft = (a[keep_t] for a in
                                      (cat_t, cat_x, cat_y, cat_m, is_aft))

R_GRID = np.geomspace(6.0, 60.0, 10)
T_GRID = np.geomspace(2.0, 80.0, 10)


def scan_max_z(t, m, x, y):
    """對每個主震在 (R, T) 網格上掃描，回傳最大 Z 值（無合格窗則 nan）。"""
    out = np.full(N_MS, np.nan)
    keep = m >= MC_SCAN
    tk, mk, xk, yk = t[keep], m[keep], x[keep], y[keep]
    for k in range(N_MS):
        d = np.hypot(xk - ms_x[k], yk - ms_y[k])
        sel = (d <= R_GRID[-1]) & (tk < ms_t[k]) & (tk >= ms_t[k] - T_GRID[-1])
        if sel.sum() < 10:
            continue
        o = np.argsort(tk[sel])
        ts_, ms_m, ds_ = tk[sel][o], mk[sel][o], d[sel][o]
        best = -np.inf
        for R in R_GRID:
            in_r = ds_ <= R
            tr, mr = ts_[in_r], ms_m[in_r]
            if tr.size < 10:
                continue
            for T in T_GRID:
                t0 = ms_t[k] - T
                j = int(np.searchsorted(tr, t0))
                tw, mw = tr[j:], mr[j:]
                n = tw.size
                if n < 10:
                    continue
                exc = mw - MC_SCAN - 0.1
                before = np.concatenate(([0.0], np.cumsum(exc)[:-1]))
                C = before - (exc.sum() / T) * (tw - t0)
                i = int(np.argmin(C))
                n1, n2 = i, n - i
                if n1 < 3 or n2 < 3:
                    continue
                d1, d2 = tw[i] - t0, ms_t[k] - tw[i]
                if d1 <= 0 or d2 <= 0:
                    continue
                z = (n2 * d1 - n1 * d2) / np.sqrt(n2 * d1 ** 2 + n1 * d2 ** 2)
                best = max(best, z)
        if np.isfinite(best):
            out[k] = best
    return out


perm_m = RNG3.permutation(cat_m.size)
perm_t = RNG3.permutation(cat_t.size)
cases = {
    "原始": (cat_t, cat_m),
    "規模隨機化": (cat_t, cat_m[perm_m]),
    "時間隨機化": (cat_t[perm_t], cat_m),
}
colors = {"原始": QUAKE_COLOR, "規模隨機化": PALETTE[2], "時間隨機化": ACCENT}

fig = go.Figure()
for name, (tt, mm) in cases.items():
    for dec, dash, width in ((False, "solid", 2.4), (True, "dot", 1.4)):
        msk = ~is_aft if dec else np.ones(cat_t.size, bool)
        z = scan_max_z(tt[msk], mm[msk], cat_x[msk], cat_y[msk])
        z = np.sort(z[np.isfinite(z)])
        fig.add_trace(go.Scatter(
            x=z, y=np.arange(1, z.size + 1) / z.size, mode="lines",
            name=f"{name}{'（去餘震）' if dec else ''}",
            line=dict(color=colors[name], width=width, dash=dash)))
apply_layout(fig, height=440, hovermode="x",
             xaxis_title="每個主震掃描到的最大 Z 值",
             yaxis_title="累積比例",
             title=(f"三組隨機化對照（合成示意，非論文重製）："
                    f"{N_MS} 個主震、{cat_t.size} 個事件"))
fig

# %% [markdown]
# 三條實線的相對位置重現了論文的定性結果：**原始最右（$Z$ 最大）、
# 規模隨機化居中且靠近原始、時間隨機化最左**。虛線是去餘震版本——
# 整體左移（窗內事件變少、顯著性下降），但**三者的先後次序不變**。
#
# 一個容易忽略的技術細節值得點出來：規模隨機化之所以會讓 $Z$ 下降，
# 機制是**掃描前先做了 $m_c$ 截切**。前兆群的規模被人為抬高，截切後
# 通過門檻的比例特別高；洗掉規模之後，通過比例回到全域平均，於是
# 前兆窗的計數 $N_2$ 縮水。**如果不做截切，規模隨機化對純計數的 $Z$
# 會完全沒有影響**——這是一個很好的提醒：所謂「與規模無關的統計量」，
# 只有在資料前處理也與規模無關時才真的與規模無關。

# %% [markdown]
# ## 15.6 不唯一性：一個主震有幾個 Ψ
#
# ### 誠實的數字
#
# 把 rectangular 演算法（真實目錄用 $R = 200$ km 上限；模擬目錄用
# $T = 1500$ 年、$m_c = 5.0$）跑完，得到的數字必須整組一起講：
#
# - **真實目錄**：47 個主震產生 **109 個識別，但只來自其中 34 個**
#   ——**13 個主震（28%）連一個合格的 Ψ 都找不到**。
# - **模擬目錄**：376 個主震產生 3113 個識別，來自其中 369 個
#   （7 個找不到）。
# - 找得到時，**真實主震平均約 3 個識別、模擬主震平均約 9 個**；
#   單一主震的最多識別數分別是 **10 個**與 **24 個**。
# - 同一主震內 $M_P$ 的變動幅度（最大減最小）中位數約 **0.2**（真實）
#   ／**0.5**（模擬）。
# - 同一主震內 $T_P$ 與 $A_P$ 的**最大最小比值**：真實資料**多半小於
#   10**；模擬資料**中位數約 100**——也就是差兩個數量級。
#
# 真實與模擬的差距有可解釋的來源：模擬目錄長得多（20,000 年 vs 數十
# 年），容得下更長的 $T_P$；模擬區域的地理範圍較小，且真實資料有定位
# 誤差（會把 $A_P$ 撐大）。
#
# ### 這句話該怎麼講
#
# 於是「這場地震的前兆期是 $N$ 年」這種說法，**在方法上就站不住**。
# 正確的講法是：
#
# > 訊號存在，但它的**參數化不唯一**。同一個主震可以有多組同樣合格的
# > $(T_P, A_P, M_P)$，而且它們沿著一條 $A_P \times T_P \approx$ 常數的
# > 線分布。
#
# **「不唯一」比「不存在」更難教，也更誠實。** 大眾論述習慣把前兆
# 問題化約成有／沒有的二分，但這裡的真實情況是第三種：有，可是量不
# 出唯一的值。
#
# 上限的截斷效應也要一起講。作者明白承認：$R \le 200$ km 與
# $T \le 1500$ 年這兩個**人訂的上限，截斷了 $A_P$ 與 $T_P$ 的取值**；
# 若把上限放寬，就會多出一些 $A_P$ 與 $T_P$ 更大的識別。任何從有限
# 時空窗口估出來的尺度關係，都先天低估大事件那一端——這與第 8 章的
# 「你看到的樣本是被誰選出來的」是同一個問題。
#
# 下面畫出同一個主震的多重識別在 $(\log_{10} T_P, \log_{10} A_P)$
# 平面上的樣子，並疊上一個真實對照：2020/10/30 愛琴海 $M\,6.7$ 的
# **兩組合法識別**（Rhoades et al. 2022 的 20 年回顧）：

# %% tags=["hide-input"]
RNG4 = np.random.default_rng(15004)
logT_s = 3.5 + np.linspace(-0.45, 0.45, 7) + RNG4.normal(0, 0.03, 7)
logA_s = 3.9 - 1.0 * (logT_s - 3.5) + RNG4.normal(0, 0.07, 7)
b_syn = np.polyfit(logT_s, logA_s, 1)[0]

# 真實對照（愛琴海 M6.7，兩組合法識別）：T_P（天）、A_P（km²）
aeg_T = np.array([10220.0, 6392.0])
aeg_A = np.array([3203.0, 8091.0])
b_aeg = np.log10(aeg_A[1] / aeg_A[0]) / np.log10(aeg_T[1] / aeg_T[0])

xline = np.array([3.0, 4.1])
yline = (logA_s.mean() + logT_s.mean()) - xline          # 斜率 −1 的參考線

fig = go.Figure()
fig.add_trace(go.Scatter(x=10 ** xline, y=10 ** yline, mode="lines",
                         name="斜率 −1（等量抵換）",
                         line=dict(color="#666666", width=2, dash="dash")))
fig.add_trace(go.Scatter(x=10 ** logT_s, y=10 ** logA_s, mode="lines+markers",
                         name="合成：同一主震的 7 個識別",
                         line=dict(color=ACCENT, width=1.2, dash="dot"),
                         marker=dict(size=10, color=ACCENT)))
fig.add_trace(go.Scatter(x=aeg_T, y=aeg_A, mode="lines+markers",
                         name="真實：愛琴海 M6.7 的兩組識別",
                         line=dict(color=QUAKE_COLOR, width=1.6, dash="dot"),
                         marker=dict(size=13, color=QUAKE_COLOR, symbol="diamond")))
fig.update_xaxes(type="log", title_text="前兆時間 T_P（天）")
fig.update_yaxes(type="log", title_text="前兆面積 A_P（km²）")
apply_layout(fig, height=460, hovermode="closest",
             title=(f"多重識別沿抵換線分布：合成擬合斜率 {b_syn:.2f}、"
                    f"愛琴海兩點連線斜率 {b_aeg:.2f}"))
fig

# %% [markdown]
# 愛琴海那一對是很好的教材，因為兩組都是**發表過的合法識別**：
# 一組 $M_P = 5.6$、$T_P = 10{,}220$ 天、$A_P = 3203$ km$^2$；另一組
# $M_P = 5.7$、$T_P = 6392$ 天、$A_P = 8091$ km$^2$。面積差 2.5 倍、
# 時間差 1.6 倍，方向相反。要誠實地說：這一對連線的斜率約
# $-2$，比等量抵換線陡——**單一一對點量不出抵換斜率**，這正是下一節
# 要做統計的理由。

# %% [markdown]
# ## 15.7 時空取捨與回歸稀釋
#
# ### 帶個別截距的迴歸
#
# 抵換是**主震內部**的現象：對同一個主震，前兆面積大的識別前兆時間
# 短、反之亦然。要量它，就不能把所有識別混在一起做一條迴歸——那樣會
# 把「主震之間的差異」與「主震內部的抵換」攪在一起（15.8 節會看到
# 攪在一起的後果有多嚴重）。正確的做法是**給每個主震一個自己的截距**：
#
# $$\log_{10} A_P = \sum_{i=1}^{n_m} a_i I_i + b\,\log_{10} T_P ,
#   \qquad\text{以及反向}\qquad
#   \log_{10} T_P = \sum_{i=1}^{n_m} a_i I_i + b\,\log_{10} A_P$$
#
# $n_m$ 是主震數、$I_i$ 是第 $i$ 個主震的指示函數（該識別屬於它為 1、
# 否則 0）、$a_i$ 是第 $i$ 個主震自己的截距、$b$ 是共用斜率。幾何上
# 這就是**一組平行線，每個主震一條**；$b$ 只由「同一主震內各識別之間
# 的相對位置」決定，主震之間的水準差異被 $a_i$ 全部吸收掉。
#
# （記號警告：這裡的 $b$ 是迴歸斜率，與第 11 章 GR 律的 $b$ 值毫無
# 關係。本節之內 $b$ 一律指斜率。）
#
# 擬合結果：
#
# | 迴歸方向 | 真實資料 $b$ | 模擬資料 $b$ |
# |---|---|---|
# | $\log_{10} A_P$ 對 $\log_{10} T_P$ | $-0.56 \pm 0.14$（$R^2 = 0.66$） | $-0.77 \pm 0.01$（$R^2 = 0.63$） |
# | $\log_{10} T_P$ 對 $\log_{10} A_P$ | $-0.33 \pm 0.08$（$R^2 = 0.62$） | $-0.67 \pm 0.01$（$R^2 = 0.62$） |
#
# 兩個方向都是負的——抵換確實存在。但兩個數字**不能直接比較**，因為
# 它們住在不同的座標系。
#
# ### 把兩個方向換到同一個座標系
#
# 反向迴歸給的是「$\log T_P$ 每變動一單位，$\log A_P$ 變動多少」的
# **倒數**。要換回 $\log A_P$ 對 $\log T_P$ 的座標系，取倒數：
#
# $$b_{\text{反向}\to\text{同座標}} = \frac{1}{-0.33} \approx -3.0 .$$
#
# 於是兩個估計是 $-0.56$ 與 $-3.0$，**理論上的等量抵換值 $-1$ 恰好
# 落在兩者之間**。這不是巧合，而是量測誤差的必然。
#
# ### 完整推導：回歸稀釋
#
# 設真實的線性關係是 $Y = \alpha + b\,X^\ast$，其中 $X^\ast$ 是**無法
# 直接觀測的真值**（真正的 $\log T_P$），我們觀測到的是帶誤差的
#
# $$X = X^\ast + \epsilon, \qquad
#   \mathbb{E}[\epsilon] = 0,\quad
#   \epsilon \perp X^\ast,\quad \mathrm{Var}(\epsilon) = \sigma_\epsilon^2 .$$
#
# 最小二乘的斜率估計是 $\hat b = \widehat{\mathrm{Cov}}(X,Y)
# / \widehat{\mathrm{Var}}(X)$。取機率極限，分子分母各算一次：
#
# $$\begin{aligned}
# \mathrm{Cov}(X, Y)
#   &= \mathrm{Cov}\bigl(X^\ast + \epsilon,\; \alpha + b X^\ast\bigr)
#    = b\,\mathrm{Var}(X^\ast) + b\,\mathrm{Cov}(\epsilon, X^\ast)
#    = b\,\mathrm{Var}(X^\ast), \\
# \mathrm{Var}(X)
#   &= \mathrm{Var}(X^\ast) + \sigma_\epsilon^2 .
# \end{aligned}$$
#
# 因此
#
# $$\hat b \;\xrightarrow{\;p\;}\;
#   b \cdot \frac{\mathrm{Var}(X^\ast)}{\mathrm{Var}(X^\ast) + \sigma_\epsilon^2}
#   \;\equiv\; b\,\lambda , \qquad 0 < \lambda \le 1 .$$
#
# $\lambda$ 稱為**衰減因子**（attenuation factor）。結論是：
# **只要自變數帶量測誤差，最小二乘的斜率就系統性地偏向零**（絕對值
# 偏小，即「偏淺」）。誤差愈大、$\lambda$ 愈小、偏得愈多；而且這是
# **偏誤**，不是變異——樣本再大也不會消失。
#
# 反向那條呢？把 $Y$ 當自變數（也帶誤差 $\delta$，
# $\mathrm{Var}(\delta) = \sigma_\delta^2$）：
#
# $$\hat b' = \frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}
#   = \frac{b\,\mathrm{Var}(X^\ast)}
#          {b^2 \mathrm{Var}(X^\ast) + \sigma_\delta^2}
#   \;\;\Longrightarrow\;\;
#   \frac{1}{\hat b'} = b\left(1
#     + \frac{\sigma_\delta^2}{b^2\,\mathrm{Var}(X^\ast)}\right).$$
#
# 括號裡大於 1，所以 $|1/\hat b'| > |b|$——**換回同一座標系之後，反向
# 迴歸的斜率系統性地偏陡**。兩個方向一個偏淺、一個偏陡，於是
#
# $$\left|\hat b\right| \;\le\; |b| \;\le\; \left|1/\hat b'\right| ,$$
#
# **真值必被夾在兩者之間**。代入真實資料：$0.56 \le |b| \le 3.0$，
# 而理論預期的等量抵換 $|b| = 1$ 舒服地落在區間內。作者的原話是「等量
# 抵換線大約落在兩條擬合線的中間」。
#
# 用模擬把這個夾擠現象跑出來——真值固定為 $-1$，看兩個方向怎麼各偏
# 一邊：

# %% tags=["hide-input"]
RNG5 = np.random.default_rng(15005)
B_TRUE = -1.0
SD_EPS, SD_DEL = 0.887, 1.425          # 兩個方向的量測誤差（見附錄 B）
n_pt = 140
Xs = RNG5.normal(0, 1.0, n_pt)                       # 真值 X*
Ys = B_TRUE * Xs
Xo = Xs + RNG5.normal(0, SD_EPS, n_pt)               # 觀測 log T_P
Yo = Ys + RNG5.normal(0, SD_DEL, n_pt)               # 觀測 log A_P
b_fwd = np.polyfit(Xo, Yo, 1)[0]
b_rev_inv = 1.0 / np.polyfit(Yo, Xo, 1)[0]

s_grid = np.linspace(0.0, 1.3, 40)
b_fwd_curve = B_TRUE / (1 + s_grid ** 2)
b_rev_curve = B_TRUE * (1 + s_grid ** 2)

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                    subplot_titles=("一次實現：兩條擬合線夾住真值",
                                    "衰減因子：誤差愈大，夾得愈開"))
xg = np.linspace(Xo.min(), Xo.max(), 2)
fig.add_trace(go.Scatter(x=Xo, y=Yo, mode="markers", name="觀測（雙向皆帶誤差）",
                         marker=dict(size=6, color=ACCENT, opacity=0.55)),
              row=1, col=1)
for slope, nm, col, dsh in ((B_TRUE, "真值 −1", "#666666", "dash"),
                            (b_fwd, f"正向 {b_fwd:.2f}", PALETTE[2], "solid"),
                            (b_rev_inv, f"反向換算 {b_rev_inv:.2f}",
                             QUAKE_COLOR, "solid")):
    fig.add_trace(go.Scatter(x=xg, y=slope * (xg - Xo.mean()) + Yo.mean(),
                             mode="lines", name=nm,
                             line=dict(color=col, width=2, dash=dsh)),
                  row=1, col=1)
fig.add_trace(go.Scatter(x=s_grid, y=b_fwd_curve, mode="lines",
                         name="正向（偏淺）", showlegend=False,
                         line=dict(color=PALETTE[2], width=2.5)), row=1, col=2)
fig.add_trace(go.Scatter(x=s_grid, y=b_rev_curve, mode="lines",
                         name="反向換算（偏陡）", showlegend=False,
                         line=dict(color=QUAKE_COLOR, width=2.5)), row=1, col=2)
fig.add_hline(y=B_TRUE, line_dash="dash", line_color="#666666", row=1, col=2)
for yv, txt in ((-0.56, "真實資料 −0.56"), (-3.03, "真實資料 −3.0")):
    fig.add_hline(y=yv, line_dash="dot", line_color="#999999", row=1, col=2)
    fig.add_annotation(x=0.05, y=yv, text=txt, showarrow=False, xanchor="left",
                       yshift=11, font=dict(size=11), row=1, col=2)
fig.update_xaxes(title_text="觀測 log₁₀ T_P（任意原點）", row=1, col=1)
fig.update_yaxes(title_text="觀測 log₁₀ A_P", row=1, col=1)
fig.update_xaxes(title_text="量測誤差標準差（真值標準差 = 1）", row=1, col=2)
fig.update_yaxes(title_text="估計斜率", range=[-4.2, 0.2], row=1, col=2)
apply_layout(fig, height=430, hovermode="closest",
             title="回歸稀釋：兩個方向各偏一邊，真值被夾在中間")
fig

# %% [markdown]
# ### 抵換的物理圖像與 $A_P T_P$ 近似守恆
#
# 斜率 $-1$ 在對數座標上的意思是 $\log A_P + \log T_P =$ 常數，也就是
#
# $$A_P \times T_P \approx \text{常數} \quad (\text{對給定的主震}) .$$
#
# 直觀圖像：**最早的前兆傾向緊貼震源，晚期的前兆散布到震源周圍的
# 廣大範圍。** 想辨識得早，就得把框畫大（否則早期那些近源事件太少，
# 撐不起一次尺度增加）；想把框畫小，就只能接受晚一點才辨識得出來。
# 要嘛「早而遠」、要嘛「晚而近」——資料不允許同時要短前兆時間又要小
# 前兆面積。對預報而言這很殘酷：**小面積、短時間的預報最有價值，
# 而它恰好是抵換線上最不可得的那一端。**
#
# 這條抵換還有第二條完全獨立的證據鏈。Rastin et al.（2021）根本不做
# Ψ 辨識，改用 EEPAS 模型當量測工具：把時間尺度參數與空間尺度參數
# 輪流固定在偏離最佳值的位置、重新擬合另一個，在紐西蘭與加州都得到
# 與斜率 $-1$ 一致的關係（$\sigma_A^2 \times 10^{a_T}$ 近似守恆）。
# **兩條方法完全不同的路徑指向同一件事**，這是抵換為真的最強論據。
# 他們還把抵換線上的三組參數混成 hybrid 模型，**不增加任何自由參數**
# 就在獨立測試期勝出——細節屬於第 16 章。
#
# 最後一個保留意見要講：作者自己指出，抵換也**可能**是統計副作用。
# 若 $\log T_P$ 與 $\log A_P$ 各自與 $M_m$ 相關，而識別程序在給定
# $M_m$ 下有取捨，兩者就會呈負相關——這與「物理上存在時空守恆量」
# 觀測上難以區分。目前沒有辦法判定。

# %% [markdown]
# ## 15.8 Simpson 陷阱：同一組資料，兩個相反的結論
#
# 把 $T_P$ 與 $A_P$ 的關係換一種算法，符號會整個翻面：
#
# - 用**全部識別**擬合（不分主震）：真實資料的斜率是
#   $-0.05 \pm 0.07$、模擬資料是 $-0.59 \pm 0.01$——**負的或近乎零**。
# - 用**每個主震取平均**再擬合：真實資料變成 $+0.65 \pm 0.08$
#   （$R^2 = 0.65$）、模擬資料 $+0.82 \pm 0.02$——**正的**，而且比
#   ER2004 的 $+0.57 \pm 0.12$（$R^2 = 0.34$）還陡、還準。
#
# 這是教科書級的 **Simpson 悖論**：**同一組資料，因為「在主震內比」
# 還是「跨主震比」，相關性可以整個翻面。** 兩個結論都對，只是在回答
# 不同的問題：
#
# - **主震內**問的是「對這一場地震，我可以怎麼在時間與空間之間換」
#   ——答案是負相關（抵換）。
# - **跨主震**問的是「大地震的前兆是不是又久又廣」——答案是正相關
#   （尺度關係）。
#
# 混在一起畫散布圖，兩個效應互相抵銷，於是真實資料得到 $-0.05$：
# **看起來像「$T_P$ 與 $A_P$ 沒有關係」，而事實是它們有兩種方向相反、
# 大小接近的關係。** 這也順帶解釋了一個二十年的舊謎：ER2004 的
# $T_P$–$A_P$ 迴歸 $R^2$ 只有 $0.34$，遠低於其他關係——因為人工識別
# 也混入了抵換。
#
# 順帶一提作者去除抵換的技巧（附錄 C 有完整說明）：他們不是單純取
# 算術平均，而是**先用帶個別截距的迴歸把主震內的抵換扣掉**，再把第
# $i$ 個主震的代表值定義為「它那條平行線在全體 $\overline{\log T_P}$
# 處的高度」。這是一個很漂亮的操作——用一個模型把不想要的變異挪走，
# 再對殘餘做統計。
#
# 畫出來就一目瞭然：

# %% tags=["hide-input"]
RNG6 = np.random.default_rng(15006)
N_M6, offs = 8, np.array([-0.35, -0.175, 0.0, 0.175, 0.35])
logT_mean = np.linspace(3.0, 4.0, N_M6)
logA_mean = 0.30 + 0.65 * logT_mean + RNG6.normal(0, 0.05, N_M6)

fig = go.Figure()
allT, allA = [], []
for i in range(N_M6):
    tt = logT_mean[i] + offs + RNG6.normal(0, 0.02, offs.size)
    aa = logA_mean[i] - 1.0 * offs + RNG6.normal(0, 0.05, offs.size)
    allT.append(tt)
    allA.append(aa)
    fig.add_trace(go.Scatter(x=tt, y=aa, mode="lines+markers", showlegend=False,
                             line=dict(color=PALETTE[i % len(PALETTE)],
                                       width=1, dash="dot"),
                             marker=dict(size=7,
                                         color=PALETTE[i % len(PALETTE)])))
allT, allA = np.concatenate(allT), np.concatenate(allA)
b_pool, a_pool = np.polyfit(allT, allA, 1)
b_mean, a_mean = np.polyfit(logT_mean, logA_mean, 1)
xg = np.linspace(2.5, 4.5, 2)
fig.add_trace(go.Scatter(x=xg, y=a_pool + b_pool * xg, mode="lines",
                         name=f"全部識別擬合：斜率 {b_pool:+.2f}",
                         line=dict(color="#444444", width=3)))
fig.add_trace(go.Scatter(x=xg, y=a_mean + b_mean * xg, mode="lines",
                         name=f"每主震取平均：斜率 {b_mean:+.2f}",
                         line=dict(color=QUAKE_COLOR, width=3, dash="dash")))
fig.add_trace(go.Scatter(x=logT_mean, y=logA_mean, mode="markers",
                         name="各主震的平均",
                         marker=dict(size=13, color=QUAKE_COLOR, symbol="x")))
apply_layout(fig, height=470, hovermode="closest",
             xaxis_title="log₁₀ T_P（天）", yaxis_title="log₁₀ A_P（km²）",
             title=(f"Simpson 陷阱（合成示意）：主震內負相關、跨主震正相關，"
                    f"混在一起得到 {b_pool:+.2f}"))
fig

# %% [markdown]
# 每一條虛線是同一個主震的多重識別（斜率 $-1$ 的抵換），紅色叉是
# 各主震的平均（沿正斜率排列）。灰色實線是把全部點混在一起的擬合
# ——它幾乎是水平的，正好落在真實資料 $-0.05$ 的量級上。**不分群畫
# 散布圖，就會得到「沒有關係」這個雙重錯誤的結論。**

# %% [markdown]
# ## 15.9 七條尺度關係的重估
#
# ### ER2004 的原始版本
#
# 先把被重估的對象列清楚。Evison & Rhoades（2004）從 47 個主震的
# 人工識別擬合出七條關係。**單位約定**：$T_P$ 以**天**為單位；
# $A_P$ 以 **$10^3$ km$^2$** 為單位（這個換算不是原文明寫的，是由
# 截距推得的——同一條 $A_P$–$M_P$ 關係在 Rhoades & Evison 2004 的
# EEPAS 論文中寫成 $\log_{10} A_P = 1.88 + 0.35 M_P$（km$^2$），
# 與下表的 $-1.12$ 恰好相差 $3.00$，正是 $10^3$；而 $T_P$ 的截距
# $1.36$ 在兩處完全相同，故單位不變）。
#
# | 關係 | $a \pm \sigma_a$ | $b \pm \sigma_b$ | RSE | $R^2$ |
# |---|---|---|---|---|
# | $\log_{10} A_P = a + b M_m$ | $-2.44 \pm 0.43$ | $0.48 \pm 0.06$ | 0.27 | 0.54 |
# | $\log_{10} T_P = a + b M_m$ | $0.08 \pm 0.38$ | $0.52 \pm 0.06$ | 0.23 | 0.65 |
# | $\log_{10}(A_P T_P) = a + b M_m$ | $-2.36 \pm 0.56$ | $1.00 \pm 0.08$ | 0.35 | **0.75** |
# | $\log_{10} T_P = a + b \log_{10} A_P$ | $3.11 \pm 0.11$ | $0.57 \pm 0.12$ | 0.32 | 0.34 |
# | $\log_{10} A_P = a + b M_P$ | $-1.12 \pm 0.31$ | $0.35 \pm 0.07$ | 0.29 | 0.48 |
# | $\log_{10} T_P = a + b M_P$ | $1.36 \pm 0.25$ | $0.40 \pm 0.04$ | 0.23 | 0.65 |
# | $M_m = a + b M_P$ | $3.16 \pm 0.35$ | $0.65 \pm 0.06$ | 0.32 | 0.71 |
#
# 三個讀法。**其一**，最後三列就是 EEPAS 三個機率核所繼承的
# $(a_A, b_A)$、$(a_T, b_T)$、$(a_M, b_M)$——第 16 章的模型完全建在
# 這三行上。把它們單獨標出來：
#
# $$M_m = a_M + b_M M_P, \qquad
#   \log_{10} T_P = a_T + b_T M_P, \qquad
#   \log_{10} A_P = a_A + b_A M_P$$ (eq:psi-scaling)
#
# **其二**，$b_M = 0.65 < 1$ 的意思是：前兆群裡最大的那幾個地震通常
# 比主震小約一個規模單位，而且主震愈大、差距愈大。**其三**，把
# $M_m = 7.0$ 代入 $\log_{10} T_P = 0.08 + 0.52 M_m$ 得
# $T_P \approx 10^{3.72} \approx 5300$ 天 $\approx 14$–15 年——這就是
# 「$M\,7$ 的前兆時間約十五年」這句廣為流傳的話的出處，也是 EEPAS
# 被歸類為**中長期**模型的原因。$A_P T_P$ 那一列的斜率恰好是
# $1.00 \pm 0.08$，$R^2$ 也是七條裡最高的 $0.75$——這是「乘積比單項
# 穩」的第一個線索（但見 15.10 節第五條的警告）。
#
# ### 演算法版：每主震取平均（真實目錄、rectangular）
#
# 重估版有四種資料集（真實／模擬 × 全部識別／每主震取平均）。這裡列
# **與 ER2004 最可比的一組**——真實目錄、每主震取平均（$n_m = 34$，
# 殘差自由度 32）：
#
# | 關係 | $a \pm \sigma_a$ | $b \pm \sigma_b$ | RSE | $R^2$ |
# |---|---|---|---|---|
# | $\log_{10} A_P = a + b M_m$ | $-3.92 \pm 0.82$ | $0.62 \pm 0.12$ | 0.40 | 0.45 |
# | $\log_{10} T_P = a + b M_m$ | $-0.86 \pm 0.56$ | $0.58 \pm 0.08$ | 0.27 | 0.60 |
# | $\log_{10}(A_P T_P) = a + b M_m$ | $-4.78 \pm 1.24$ | $1.21 \pm 0.18$ | 0.61 | 0.57 |
# | $\log_{10} T_P = a + b \log_{10} A_P$ | $2.89 \pm 0.05$ | $0.65 \pm 0.08$ | 0.26 | 0.65 |
# | $\log_{10} A_P = a + b M_P$ | $-2.89 \pm 0.56$ | $0.59 \pm 0.10$ | 0.39 | 0.50 |
# | $\log_{10} T_P = a + b M_P$ | $0.13 \pm 0.38$ | $0.54 \pm 0.07$ | 0.26 | 0.65 |
# | $M_m = a + b M_P$ | $2.19 \pm 0.30$ | $0.84 \pm 0.06$ | 0.20 | **0.88** |
#
# ### 三個結論
#
# **一、所有以規模為自變數的關係，斜率都是正的，且在參數不確定度內
# 與 ER2004 沒有顯著差異。** 這是重估最重要的結果：一個**不預設任何
# 尺度觀念**的演算法，仍然復現了人工程序歸納出的斜率。
#
# **二、擬合品質常常持平甚至更好。** 對照兩張表的 RSE 與 $R^2$：
# $T_P$–$M_m$、$T_P$–$M_P$ 兩列的 RSE 與 $R^2$ 與 ER2004 相當；
# $M_m$–$M_P$ 的 $R^2$ 從 $0.71$ 升到 $0.88$、RSE 從 $0.32$ 降到
# $0.20$；$T_P$–$A_P$ 的 $R^2$ 從 $0.34$ 升到 $0.65$（15.8 節解釋了
# 為什麼）。
#
# **三、截距有系統性位移，而且方向可解釋。** 真實資料的 $T_P$ 平均
# 比 ER2004 **小約 3 倍**、$A_P$ 小到 **4–5 倍**；模擬資料反過來，
# $T_P$ **大 2–3 倍**、$A_P$ **小超過一個數量級**。作者歸因於目錄的
# 長度與空間範圍：模擬目錄長得多，容得下更長的 $T_P$；模擬區域地理
# 範圍較小、且真實資料有定位誤差（會把 $A_P$ 撐大）。
#
# 還有一個必須配著讀的警告：**用「全部識別」擬合時，多數關係的
# $R^2$ 極低**——真實資料的 $A_P$–$M_m$ 只有 $0.13$，模擬資料只有
# $0.006$。抵換與跨主震尺度關係方向相反、互相抵銷（15.8 節），
# 必須先按主震平均才看得到訊號。**同一份資料，$R^2$ 可以從 0.006
# 變到 0.45，差別只在你有沒有分群。**

# %% [markdown]
# ## 15.10 常見誤解與陷阱
#
# **一、「Ψ 是一個可以事前偵測的前兆。」** 不是。整套辨識程序都需要
# **先知道主震的位置、時間與規模**——$m_c$ 由 $M_m$ 決定、搜尋區以主震
# 為心、接受條件用到 $M_m - M_P > 0.4$。Ψ 是一個**回溯**現象。
# 它之所以仍有預報價值，是因為 EEPAS 走了完全不同的一步：放棄辨識，
# 假設每個地震都是某個尺度上的前兆（第 16 章）。
#
# **二、「這場地震的前兆期是 $N$ 年。」** 方法上站不住。同一主震平均
# 有約 3 組（真實）到 9 組（模擬）同樣合格的識別，$T_P$ 的最大最小
# 比在模擬目錄中位數達 100 倍。要講就講區間，並且說明用的是哪個演算法
# 與哪組門檻。
#
# **三、「自動化 = 客觀 = 更可信。」** 演算法只是把主觀性從「調框框」
# 搬到「訂門檻」——$r > 3$、兩個 $0.4$、中央 20%、長寬比 0.5–2.0、
# $R \le 200$ km、$T \le 1500$ 年，每一個都是人訂的，作者自己說這些
# 演算法「不是唯一的，也不是不可爭議的」。真正的進步是**判準被公開、
# 可重跑、可質疑**，而不是主觀性消失了。
#
# **四、「$Z$ 值與規模無關，所以它是純淨的計數統計量。」** 只對一半。
# $Z$ 的**公式**不含規模，但送進 $Z$ 的**事件集合**是用 $m_c$ 篩過的。
# 15.5 節的合成實驗就示範了：規模隨機化之所以能讓 $Z$ 下降，正是透過
# 這道截切。任何「無關某變數」的宣稱，都要連同前處理一起檢查。
#
# **五、「$A_P T_P$ 對 $M_m$ 的 $R^2 = 0.75$，最高，所以它抓到了物理。」**
# 危險。兩個**各自**與 $M_m$ 正相關、彼此**負相關**的量，取乘積（對數
# 相加）之後散布必然縮小、$R^2$ 必然上升——這是代數，不是發現。
# **高 $R^2$ 不等於發現物理，有時只等於座標選得好。** 這一條要和第 8
# 章的多重比較與選擇性報告串起來講。
#
# **六、「時間隨機化實驗證明了真實地震有前兆。」** 措辭錯了兩層。
# 第一，那組實驗是在 **RSQSim 模擬目錄**上做的；第二，它證明的是
# 「$Z$ 訊號依賴時間結構」，不是「前兆存在於真實地球」。正確的講法
# 見 15.5 節末的引文框。
#
# **七、「28% 找不到 Ψ 表示那些地震沒有前兆。」** 也不能這樣說。
# 找不到的原因至少有四種混在一起：目錄在那個時期不完整、目錄長度不夠
# （前兆時間超出目錄）、真的沒有尺度增加、或者門檻設得太嚴。**「演算
# 法在這組門檻下找不到」與「不存在」是兩件事。**
#
# **八、「不分群的散布圖可以看出 $T_P$ 與 $A_P$ 有沒有關係。」**
# 15.8 節整節在反駁這一條。真實資料混在一起是 $-0.05$（看起來沒關係），
# 分群後主震內是 $-0.56$、跨主震是 $+0.65$。**畫散布圖之前先問：
# 這些點是不是來自不同的群？**

# %% [markdown]
# ## 15.11 研究前沿與未解問題
#
# ### 物理機制缺席，而作者自比 Omori
#
# 2024 那篇論文最好的一段自我定位是這個類比：Omori 在 1894 年發現
# 1891 年濃尾地震後有感餘震的發生率隨時間約以 $1/T$ 衰減；一個世紀
# 之後，這個現象在無數大地震上被確認、參數被當成構造環境的函數研究、
# 統計模型（第 12–14 章的 Omori–Utsu 與 ETAS）被建立起來——**但物理
# 機制始終沒被寫下**，直到 Dieterich（1994）用速率–狀態摩擦提出一個
# 本構律，才把應力變化與時間相依的成核連起來。
#
# **Ψ 目前處在同一個階段**：現象被確立、參數被量化、統計模型（EEPAS）
# 上線運轉，物理機制仍未寫下。這不是 Ψ 的獨有困境，而是統計地震學的
# 常態——**多數有用的定律（Omori、GR、Ψ）都是先有經驗規律，機制在
# 幾十年後才補上，或至今沒補上。**
#
# ### RSQSim 自己長出了 Ψ
#
# 前沿裡最令人振奮的一條：RSQSim 這個只放進速率–狀態摩擦、斷層應力
# 交互作用與破裂準則的模擬器，**沒有人把 Ψ 寫進它的規則，Ψ 卻自己
# 出現了**，而且尺度關係與真實資料一致（15.9 節）。作者的推論是：
# 尺度關係很可能是**已知物理的產物**，反映大地震前的應力累積。
#
# 但要小心不要多說一步。**湧現不等於機制被理解**——我們只知道「這些
# 規則足以產生 Ψ」，不知道「Ψ 對應規則裡的哪一段」。不過這打開了一
# 條可操作的研究路徑：**系統性地改變模擬器的輸入參數，看 Ψ 怎麼變。**
# 這是真實資料永遠做不到的實驗。
#
# ### 加載速率控制前兆時間
#
# 沿著這條路徑已經有一個結果。Christophersen, Rhoades & Colella
# （2017）發現：把模擬器裡斷層的滑移速率降低時，**前兆時間與滑移速率
# 的降幅成反比**。與此呼應的觀測是澳洲——板塊內部、應變速率比板塊
# 邊界低約三個數量級——的資料點明顯偏離尺度關係，$T_P$ 比預期大約
# **10 倍**。
#
# 這代表 **Ψ 的尺度關係會隨構造環境改變**，不是普適常數。對台灣這種
# 高應變速率的板塊邊界環境，含意是前兆時間應該偏短端。但這只是尺度
# 上的推論——台灣的在地化工作仍在進行中，本書不對台灣的 Ψ 參數作任何
# 具體宣稱。
#
# ### 如何把取捨形式化進 EEPAS
#
# 最後一個明確的開放問題。EEPAS 的時間核 $f$ 與空間核 $h$ 是**相乘**
# 的，也就是假設時間與位置獨立；但 15.7 節的兩條獨立證據鏈都說它們
# 負相關。目前唯一的補救是 Rastin et al.（2021）的權宜做法：沿抵換線
# 取三組參數平均，零新增自由度、在獨立測試期勝出。作者自承
# **「如何最佳地把抵換納入 EEPAS 仍未解決」**——這是一個適合當課堂
# open problem 的缺口，而且它的形狀很清楚：需要一個帶負相關的二維
# 時空聯合核，取代現在的乘積形式。
#
# 其餘尚未解決的還有：把長期地震率納入模型、三維（含深度）版本、
# 隨時空變化的目錄完整度，以及全球尺度的 Ψ 統計。

# %% [markdown]
# ## 15.12 附錄：本章推導細節
#
# ### A. $C(t)$ 的最小值為什麼是「速率換檔點」
#
# 承 15.1 節性質三。定義 $[t_s, t]$ 內的平均規模累積速率
#
# $$\bar\rho(t) = \frac{1}{t - t_s}\sum_{t_s < t_i \le t}(M_i - m_c - 0.1) ,$$
#
# 則 $C(t) = [\bar\rho(t) - k](t - t_s)$。設 $t^\ast$ 是最小值點，
# 於是對任意 $t$，$C(t) \ge C(t^\ast)$。取 $t = t_f$（此時
# $C(t_f) = 0$，性質一）得 $C(t^\ast) \le 0$，即
#
# $$\bar\rho(t^\ast) \;\le\; k .$$
#
# 另一方面，把 $[t_s, t_f]$ 拆成兩段，總量守恆給出
#
# $$\bar\rho(t^\ast)(t^\ast - t_s)
#   + \rho_{\text{後}}\,(t_f - t^\ast) = k\,(t_f - t_s) ,$$
#
# 其中 $\rho_{\text{後}}$ 是 $[t^\ast, t_f]$ 的平均速率。既然
# $\bar\rho(t^\ast) \le k$，必有 $\rho_{\text{後}} \ge k$。
# **$t^\ast$ 把時窗切成「低於全段平均」與「高於全段平均」兩塊，而且
# 這個切法在所有切點裡讓落差最大**——這就是「速率換檔點」的精確意思。
#
# 兩個實務推論。**其一**，$C(t)$ 的最小值位置**依賴 $t_s$ 與 $t_f$
# 的選擇**：改變時窗就改變 $k$，也就改變 onset。這正是 rectangular
# 演算法要掃描 $T$、並要求「onset 落在中央 20%」的原因——後者等於
# 要求先前期與前兆期長度相當，使 $k$ 不被單邊主導。**其二**，
# $C(t)$ 量的是**相對**於全段平均的偏離，所以它偵測不到「整段都在
# 緩慢抬升」這種沒有換檔的變化。
#
# ### B. 回歸稀釋的兩個分支與夾擠區間
#
# 承 15.7 節。記 $v = \mathrm{Var}(X^\ast)$、$\sigma_\epsilon^2$ 為
# 自變數誤差、$\sigma_\delta^2$ 為應變數誤差，真斜率 $b$。兩個方向的
# 極限分別是
#
# $$\hat b \to \frac{b}{1 + \sigma_\epsilon^2 / v},
#   \qquad
#   \frac{1}{\hat b'} \to b\left(1
#     + \frac{\sigma_\delta^2}{b^2 v}\right).$$
#
# 兩個誤差相等且 $|b| = 1$ 時（記 $s^2 = \sigma^2/v$），兩支恰好是
# $-1/(1+s^2)$ 與 $-(1+s^2)$——在對數尺度上**對稱地夾住 $-1$**，這正是
# 圖中右панель的兩條曲線。真實資料的兩個估計不對稱（$-0.56$ 與
# $-3.0$），反解出
#
# $$\frac{\sigma_\epsilon^2}{v} = \frac{1}{0.56} - 1 \approx 0.79,
#   \qquad
#   \frac{\sigma_\delta^2}{b^2 v} = 3.0 - 1 = 2.0 ,$$
#
# 也就是 $A_P$ 的相對量測誤差約為 $T_P$ 的 2.5 倍。這在物理上說得通
# ——$A_P$ 由震央的**空間分布**決定，直接受定位誤差影響（15.9 節作者
# 用同一個理由解釋真實目錄的 $A_P$ 偏大）。
#
# 若要給一個點估計而非區間，標準做法是取兩支的幾何平均（等價於
# **reduced major axis** 迴歸）：
#
# $$b_{\rm RMA} = -\sqrt{0.56 \times 3.0} \approx -1.30 .$$
#
# 與理論的 $-1$ 同量級，但不完全相等——**要誠實說「相容」，不要說
# 「證實」**。
#
# ### C. 從兩個方向的斜率反推主震內的相關係數
#
# 對任何最小二乘配對，正向與反向斜率的乘積等於相關係數的平方：
#
# $$\hat b_{Y|X}\,\hat b_{X|Y}
#   = \frac{S_{XY}}{S_{XX}} \cdot \frac{S_{XY}}{S_{YY}}
#   = \frac{S_{XY}^2}{S_{XX}S_{YY}} = \hat r^2 .$$
#
# 帶個別截距的迴歸（15.7 節）等價於先對每個主震做組內去中心化、再做
# 普通最小二乘，所以同一恆等式在**組內**平方和上成立。代入真實資料：
#
# $$\hat r_{\text{組內}} = -\sqrt{0.56 \times 0.33} \approx -0.43 .$$
#
# 這個數字值得記住：**主震內 $\log T_P$ 與 $\log A_P$ 的相關係數只有
# $-0.43$，抵換是一個明確但相當鬆散的趨勢。** 表中那兩個 $R^2$
# （$0.66$、$0.62$）看起來高得多，但它們是**含 34 個主震截距**的整體
# 模型解釋力，絕大部分來自主震之間的水準差異，不是來自抵換本身。
# **看到高 $R^2$ 要先問：那是誰貢獻的？**
#
# ### D. 去除抵換的平均：作者的操作
#
# 15.8 節提到的「每主震取平均」不是算術平均。令 $a_i, b$ 是
# $\log_{10}A_P$ 對 $\log_{10}T_P$ 那條帶個別截距迴歸的擬合值，
# $\overline{\log_{10}T_P}$ 是**全體識別**的 $\log_{10}T_P$ 平均。
# 第 $i$ 個主震的代表值定為
#
# $$\bigl[\overline{\log_{10}A_P}\bigr]_i
#   = a_i + b\,\overline{\log_{10}T_P} ,$$
#
# 亦即「第 $i$ 條平行線在共同橫座標處的高度」；$T_P$ 的代表值以反向
# 迴歸同法求得。這樣做的用意很明確：**若直接取算術平均，各主震識別
# 在抵換線上的分布位置不同（有的偏早而遠、有的偏晚而近），平均值會
# 帶著這個位置差異；改用共同橫座標求截距，等於把所有主震拉到同一個
# 參考點再比較。** 這是一個通用技巧——用模型把已知的干擾變異挪走，
# 再對殘餘做統計——值得從 Ψ 這個例子學起來。
#
# ---
#
# 回頭看這一章走過的路：一條被兩端釘死的 $C(t)$ 曲線、四個從曲線上
# 讀出的變數、一套把手工調框換成網格搜尋的演算法、一個完全不含規模
# 的 $Z$ 值、三組在物理模擬目錄上做的隨機化對照，最後是三個統計陷阱
# ——不唯一、回歸稀釋、Simpson 反轉。
#
# 壓成一句話：**Ψ 的正確講法不是「有沒有前兆」，而是「有一種不同於
# 餘震衰減的震前時空叢集，但它的參數化不唯一、機制未知」。** 這句話
# 比「大地震前會有前兆」弱得多，卻是資料真正支持的全部。而 Ψ 這條
# 研究線最值得學的其實不是結論，是它示範了一套矯正程序：**把判準寫成
# 程式碼、允許失敗被記錄、允許醜陋的答案出現、設隨機化對照、再拿
# 獨立資料前瞻測試。**
#
# 但這一整章也留下一個看似致命的問題。Ψ 的辨識需要事先知道主震在哪、
# 多大——那它對預報還有什麼用？28% 的主震連一個識別都沒有，剩下的
# 又各有好幾組互相矛盾的參數，這樣的現象怎麼可能變成模型？
#
# {doc}`第 16 章 <16_eepas_ppe>`要講的，就是那一步漂亮的側身：
# **既然辨識不出來，那就不要辨識。** 假設每一個地震都是某個未來更大
# 地震的前兆，它預告的規模、時間與範圍全部由 {eq}`eq:psi-scaling`
# 那三條迴歸按它自己的規模給定——迴歸的殘差是常態的，於是三條迴歸
# 直接變成三個機率核。這一步繞開了本章列出的**全部**困難，代價是
# 預報強度被稀釋。那筆交易划不划算，下一章見分曉。
