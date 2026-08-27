# %% [markdown]
# # 16. 從預報到危害度：PSHA
#
# 第二部到目前為止都在回答「**會不會有地震**」。但對工程師、
# 保險公司、建築規範的制定者來說，這不是最終的問題。房子不在乎
# 震央在哪、規模多大——房子只在乎**腳下的地表會晃多大**。
#
# 把「地震發生的機率」轉換成「地動強度的機率」，就是**機率式
# 地震危害度分析**（Probabilistic Seismic Hazard Analysis,
# PSHA，Cornell 1968）。你家的建築耐震設計、核電廠選址、
# 地震保險費率，背後全是這套計算。它也是地震預報這門學問
# 通往社會的最寬的一座橋——而橋的每一塊磚，都是前面幾章
# 教過的東西。
#
# 本章以 Baker 的教學白皮書（*Probabilistic Seismic Hazard
# Analysis*, White Paper v2.0.1, 2013；已擴充為 Baker, Bradley
# & Stafford 2021 教科書）為骨架，只講觀念，搭配一個玩具計算。
#
# ## 16.1 為什麼要「機率式」
#
# 傳統的替代方案是決定論式的：「假設最壞情境——最近的斷層
# 發生最大地震——算出地動，照著設計。」聽起來穩健，實則有
# 兩個致命問題：「最壞情境」沒有客觀定義（最大規模？最近
# 距離？地動要取平均還是平均加幾個標準差？），而且它完全
# 忽略**發生率**——一條一萬年動一次的大斷層與一條一百年
# 動一次的小斷層，對設計的意義天差地遠。
#
# PSHA 的回答是：不挑情境，**把所有可能的地震、所有可能的
# 地動，連同各自的機率，全部積分起來**。Baker 把計算拆成
# 五個步驟：
#
# 1. 找出所有可能產生破壞性地動的**震源**（斷層或面震源）；
# 2. 描述每個震源的**規模分布**——就是 GR 律（第 5、10 章）；
# 3. 描述震源到場址的**距離分布**；
# 4. 給定規模與距離，用**地動預測方程（GMPE**）描述地動強度
#    的分布——經驗迴歸式，殘差取對數常態，典型的 σ 約 0.5–0.6
#    個自然對數單位（意思是：同樣的地震、同樣的距離，一個
#    標準差就是約 1.8 倍的 PGA——**地動的不確定性比震源的
#    不確定性還大**，這是 PSHA 最反直覺的事實之一）；
# 5. 用全機率定理把以上全部積分，得到「場址地動超越某強度」
#    的年發生率。
#
# ## 16.2 玩具 PSHA：一個場址、兩條斷層
#
# 用最小的例子把五步驟走一遍：一個場址，附近有兩條斷層——
# 近而小（10 km 外，最大規模 6.0，活動較頻繁）與遠而大
# （25 km 外，最大規模 8.0，活動較稀）。規模服從截斷 GR
# （b=1），GMPE 用 Baker 白皮書示範的 Cornell 型式。對每個
# 地動門檻 $x$，把兩條斷層、所有規模的貢獻積分起來：

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from scipy.stats import norm

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

# 兩條斷層：名稱、距離 (km)、M≥5 年發生率、最大規模
faults = [("近而小的斷層", 10.0, 0.10, 6.0),
          ("遠而大的斷層", 25.0, 0.05, 8.0)]
B = np.log(10)                                    # b = 1（自然對數版）
M_MIN = 5.0

def gmpe_lnpga(m, r):
    """Cornell 型 GMPE（Baker 白皮書示範式；PGA 單位 g）。"""
    return -0.152 + 0.859 * m - 1.803 * np.log(r + 25)

SIGMA = 0.57

def exceed_rate(x, dist, rate5, m_max):
    """單一震源對「PGA > x」年率的貢獻（對規模數值積分）。"""
    m = np.linspace(M_MIN, m_max, 300)
    # 截斷指數分布的規模密度
    fm = B * np.exp(-B * (m - M_MIN)) / (1 - np.exp(-B * (m_max - M_MIN)))
    p_ex = 1 - norm.cdf((np.log(x) - gmpe_lnpga(m, dist)) / SIGMA)
    return rate5 * np.trapezoid(p_ex * fm, m)

x_grid = np.logspace(-2, 0.3, 80)
total = np.zeros_like(x_grid)
fig = go.Figure()
for (name, dist, rate5, m_max), color in zip(faults, PALETTE):
    lam = np.array([exceed_rate(x, dist, rate5, m_max) for x in x_grid])
    total += lam
    fig.add_trace(go.Scatter(x=x_grid, y=lam, mode="lines", name=name,
                             line=dict(color=color, dash="dot")))
fig.add_trace(go.Scatter(x=x_grid, y=total, mode="lines", name="總危害",
                         line=dict(color=ACCENT, width=3)))
fig.add_hline(y=1 / 475, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text="1/475 年（50 年 10%）")
apply_layout(fig, title="玩具 PSHA 危害曲線：年超越率 vs 地動強度",
             xaxis_title="PGA（g）", yaxis_title="年超越率 λ(PGA > x)",
             xaxis_type="log", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 這就是 **危害曲線（hazard curve）**——PSHA 的最終輸出。
# 讀法：垂直找一個年率（例如工程規範常用的 1/475），水平對到
# 的 PGA 就是設計地動。注意兩條斷層的角色變化：低強度端由
# 近而小的斷層主導（它常動），愈往高強度端，遠而大的斷層佔比
# 愈高、並在最高強度端反超（只有夠大的地震配上運氣不好的地動
# 殘差，才晃得出超過 1 g）。
#
# 這裡順便拆掉全工程界最常被誤解的一個詞。「475 年回歸期」
# **不是**「每 475 年來一次」，它只是年超越率 1/475 的倒數——
# 平均間隔。把年率換算成「50 年內至少超越一次的機率」要再
# 借一個假設（通常是 Poisson）：$P = 1 - e^{-\lambda t}$，
# 代入 $\lambda = 1/475$、$t = 50$ 年，得到約 10%——這就是
# 「50 年 10% 超越機率 ≈ 475 年回歸期」的全部由來。Baker
# 甚至建議乾脆只報年率、不報回歸期，免得誤會。
#
# ```{admonition} 「475 年才來一次，我這輩子遇不到」
# :class: warning
# 這句話錯兩次。第一，回歸期是平均值，不是週期——明年就超越
# 的機率每年都在。第二，475 年講的是「這個場址的地動超越
# 設計值」的平均間隔，不是「大地震發生」的間隔；一個場址的
# 地動可以被許多不同斷層的地震超越。
# ```
#
# ## 16.3 反聚合：危害是誰貢獻的？
#
# PSHA 把所有情境加總，是它的優點也是代價：算完之後，「危害
# 主要來自哪種地震？」反而看不到了。**反聚合（deaggregation）**
# 就是把積分拆回去——問「已知 PGA 超越了 $x$，它是規模多大、
# 距離多遠的地震造成的機率各是多少」（一個貝氏定理的直接
# 應用）。用玩具模型算兩個危害水準的規模貢獻：

# %% tags=["hide-input"]
m_edges = np.arange(5.0, 8.01, 0.5)
fig = go.Figure()
for x_val, color in [(0.2, PALETTE[0]), (1.0, PALETTE[1])]:
    contrib = []
    for m_lo, m_hi in zip(m_edges[:-1], m_edges[1:]):
        c = 0.0
        for name, dist, rate5, m_max in faults:
            if m_lo >= m_max:
                continue
            m = np.linspace(m_lo, min(m_hi, m_max), 60)
            fm = (B * np.exp(-B * (m - M_MIN))
                  / (1 - np.exp(-B * (m_max - M_MIN))))
            p_ex = 1 - norm.cdf((np.log(x_val) - gmpe_lnpga(m, dist)) / SIGMA)
            c += rate5 * np.trapezoid(p_ex * fm, m)
        contrib.append(c)
    contrib = np.array(contrib) / np.sum(contrib) * 100
    fig.add_trace(go.Bar(x=[f"{a:.1f}–{b:.1f}" for a, b in
                            zip(m_edges[:-1], m_edges[1:])],
                         y=contrib, name=f"PGA > {x_val} g",
                         marker_color=color))
apply_layout(fig, title="反聚合：不同危害水準下，各規模區間的貢獻",
             xaxis_title="規模區間", yaxis_title="貢獻比例（%）",
             barmode="group", hovermode="x")
fig

# %% [markdown]
# 規律一目瞭然：**危害水準愈高，主控的震源愈大、愈遠**。
# PGA > 0.2 g 的危害幾乎全由近斷層的中等地震包辦；PGA > 1 g 時，
# 遠斷層的大地震（M≥6.5）已搶下超過三分之一的貢獻，且愈往
# 高強度端佔比愈高。反聚合是 PSHA 與工程設計之間的翻譯器——
# 選設計地震波、做情境演練，都從這張圖出發。
#
# ## 16.4 台灣地震模型（TEM）
#
# 台灣有自己的國家級 PSHA。**台灣地震模型（Taiwan Earthquake
# Model, TEM**）由學界團隊建立，兩代正式發布：
#
# - **TEM PSHA2015**（Wang et al. 2016, *TAO*）：採用地質調查
#   辨識的 38 條孕震構造，加上淺層、隱沒帶板塊內、隱沒帶板塊
#   間三類背景地震活動，各配對應的 GMPE。結論的定性圖像：
#   **危害最高的區域在西南部與東部縱谷**；西部人口稠密都會中
#   以台南（短週期）與台中（長週期）最受關注——原文特別點名
#   台南的低樓層建物與台中的高樓層建物的耐震設計。
# - **TEM PSHA2020**（Chan et al. 2020, *Earthquake Spectra*）：
#   構造資料庫擴充並納入三維幾何、允許**多構造同時破裂**
#   （0403 花蓮與 1999 集集都提醒過我們這件事）、加入場址
#   放大效應，並對斷層源導入 **BPT（Brownian passage time）**
#   更新模型——斷層有「記憶」：剛破裂過的斷層短期內再破裂
#   的機率低，隨應力重新累積而升高。下一代模型（PSHA2025）
#   正在發展中。
#
# BPT 值得多停一秒，因為它是這一章與前面所有章節思想上的
# 交會點：GR 律與 Poisson 假設說地震「無記憶」，ETAS 說地震
# 短期內「彼此招引」（叢集），BPT 說單一斷層長期上「自我
# 抑制」（週期性）——三種記憶結構，各自主宰不同的時間尺度
# 與空間對象。一個完整的危害模型，三種都要裝。
#
# ## 16.5 時變危害：把 OEF 接上 PSHA
#
# 傳統 PSHA 是**時不變**的：Poisson 假設下，今年與明年的危害
# 曲線一模一樣。但第 11–13 章整整三章都在說：短期內地震率
# 可以比背景高出百倍千倍。大地震剛過的城市，重建決策面對的
# 危害顯然不是那條長期平均曲線。把 OEF 的時變地震率餵進
# PSHA 的積分（把齊次 Poisson 換成非齊次 Poisson），就得到
# **時變危害**：

# %% tags=["hide-input"]
def weekly_exceed(x, boost=1.0):
    """一週內超越機率：背景率 × boost（主震後餘震率放大）。"""
    lam_wk = sum(exceed_rate(x, d, r5 * boost, mm) / 52.18
                 for _, d, r5, mm in faults)
    return 1 - np.exp(-lam_wk)

x_grid2 = np.logspace(-2, 0.2, 60)
fig = go.Figure()
for boost, name, color in [(1.0, "平時（長期背景）", ACCENT),
                           (200.0, "大地震後第一週（餘震率放大，示意）",
                            QUAKE_COLOR)]:
    fig.add_trace(go.Scatter(x=x_grid2,
                             y=[weekly_exceed(x, boost) for x in x_grid2],
                             mode="lines", name=name,
                             line=dict(color=color, width=2.5)))
apply_layout(fig, title="時變危害示意：主震後一週的超越機率曲線整條上移",
             xaxis_title="PGA（g）", yaxis_title="一週內超越機率",
             xaxis_type="log", yaxis_type="log", hovermode="x")
fig

# %% [markdown]
# 主震後第一週，各強度的超越機率整條上移約兩個數量級——
# 平時每週萬分之一的搖晃機率，序列期間變成百分之幾。這不是
# 學術演習：ICEF 報告（Jordan et al. 2011）明確要求 OEF 必須
# 「與 PSHA 的長期預報一致地」輸出完整的危害描述——**地動
# 超越機率，而不只是地震發生機率**。實作上的代表作有兩個：
# 義大利把 OEF 地震率接上建物清冊與易損性模型，每日產出
# 預期損失（OELF）；紐西蘭為 Canterbury 重建建立了 50 年
# 時變危害模型，Kaikōura 之後更直接用它計算「地動超越規範
# 設計值的機率增益」，支撐了中紐西蘭無筋磚造建築的**強制
# 補強**政策。第 13 章的機率溝通難題，在這裡有了最實在的
# 出口：機率變成了法規與工程決策。
#
# ## 16.6 反思：規範背後全是模型
#
# 走完這一章，回頭看「一個數字算得出來，不等於它站得住腳」
# 這句話——PSHA 是它的終極考場。一條危害曲線背後疊著多少
# 層假設：目錄的完整性與除叢（第 10 章——記得除叢會偏移
# b 值、而 b 值直接進入危害積分）、GR 律外插到沒觀測過的
# 大規模、GMPE 從有限的強震紀錄迴歸並外插、Poisson 或 BPT
# 的發生模型、斷層幾何與滑移率的地質判讀。每一層都有不確定
# 性（除叢與 b 值的細節見{doc}`第 10 章 <10_seismicity_stats>`），
# 而最終那個「475 年回歸期的設計地動」以工程規範的姿態
# 出現時，所有的不確定性都藏進了一個看起來很權威的數字裡。
# 這不是要你不信任規範——而是要你知道：規範是**模型的輸出**,
# 模型會更新（TEM 兩代已發布、第三代進行中，就是證據），而懂得掀開引擎蓋的
# 人，才有資格參與更新它。
#
# 第二部的工具至此全部到齊：統計地基、短中長期模型、組合、
# 檢驗、危害。{doc}`最後一章 <17_taiwan_outlook>`，我們把所有
# 東西帶回家——用台灣的
# 目錄、台灣的地震、台灣的預報現況，做一次總整理，並誠實
# 地問：台灣的下一步是什麼？
