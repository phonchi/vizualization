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
# # 1. 地震預報從什麼問題開始
#
# 一次明顯的地震之後，接下來幾天還會不會有較大的地震，是很自然的疑問。
# 我們也可能在沒有地震的日子，想知道某個地方未來幾十年的地動危害。
# 這兩個問題的時間尺度不同，但都需要把有限的觀測轉成對未來的判斷。
# 這個網站從這裡出發：先理解地震活動中有哪些能利用的統計規律，再學習如何
# 把規律寫成模型、如何檢查它，最後回到臺灣的物理觀測，看看資料還能增加什麼資訊。
#
# 我們不需要一開始就知道點過程或ETAS。先把問題問清楚，比先選一個複雜模型
# 更重要。「未來可能有地震」幾乎沒有提供資訊；「指定區域接下來七天內，至少
# 發生一次指定規模以上地震的機率」則界定了要預報的事件，也讓之後的評估有依據。
#
# ## 1.1 預報的是尚未發生的活動
#
# 日常說的地震預測，常指在很窄的時間、地點與規模範圍內，事先指出一場地震。
# 目前沒有能穩定達成這種要求的方法。統計預報改以機率描述未來活動，仍然需要
# 具體的區域、期間與目標，並把不確定性明確寫進答案；檢查標準並不因此降低。
#
# 地震預警的起點不同：地震已經發生，系統利用較早收到的資訊，估計震波將對
# 其他地方造成的影響。預警與預報都使用地震觀測，也都涉及不確定性，但前者
# 處理已開始的事件，後者處理未來可能發生的活動。這裡主要學習後者。
#
# 機率預報不是單看一場地震來決定對錯。若許多可比較的情境都報出相近機率，
# 我們希望實際發生比例與它相符，也希望它比只用長期平均的基準提供更多資訊。
# 這兩種要求會在預報檢驗篇分開說明；現在先記住，模型必須接受沒有參與建模的
# 未來資料檢查，才知道它是否有可延續的預報能力。
#
# ## 1.2 地震目錄為什麼提供線索
#
# 地震活動並非每一天、每個地方都一樣。有些差異來自長期的區域背景，有些
# 跟剛發生的事件有關。餘震序列尤其容易觀察：活動一開始很密集，之後通常
# 逐漸減少，但中間也可能再出現活躍期。這讓「過去剛發生什麼」成為可能有用的資訊。
#
# 先看一份真實目錄的概要。下圖使用2024年春季的臺灣目錄，包含花蓮主震與
# 後續活動。左圖比較相鄰事件間隔，右圖比較每日事件數；參考曲線所表示的
# 固定率Poisson模型，暫時可以理解成「始終以相同平均率、互不影響地出現」。
# 下一章會先介紹目錄的構成，再逐步說明這個模型的假設。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gdms_toolkit.download import CACHE_DIR
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

cat = pd.read_csv(CACHE_DIR / "catalog_2024spring.csv", parse_dates=["time"])
window_start = pd.Timestamp("2024-03-01", tz=cat.time.dt.tz)
window_end = pd.Timestamp("2024-07-01", tz=cat.time.dt.tz)
cat = cat.loc[(cat.time >= window_start) & (cat.time < window_end)]
t_sorted = cat.time.sort_values()
dt_hours = t_sorted.diff().dt.total_seconds().dropna() / 3600
mean_dt = dt_hours.mean()

# 過度離散：每日事件數的變異數 / 平均
calendar = pd.date_range("2024-03-01", "2024-06-30", freq="D", tz=t_sorted.dt.tz)
daily = t_sorted.dt.floor("D").value_counts().reindex(calendar, fill_value=0)
fano = daily.var() / daily.mean()

edges = np.arange(0, 12.25, 0.25)
hist_counts, _ = np.histogram(dt_hours, bins=edges)
hist = hist_counts / (len(dt_hours) * np.diff(edges))
centers = edges[:-1] + 0.125

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.11,
                    subplot_titles=("事件間隔分布", "每日事件數：離散程度"))
fig.add_trace(go.Bar(x=centers, y=hist, name="觀測", marker_color=ACCENT,
                     opacity=0.75), row=1, col=1)
fig.add_trace(go.Scatter(x=centers, y=np.exp(-centers / mean_dt) / mean_dt,
                         mode="lines", name="同平均間隔的指數參考",
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
apply_layout(fig, title=f"固定率能描述這份目錄嗎？（每日事件數變異數／平均 = "
                        f"{fano:.0f}，Poisson 應為 1）",
             height=420, hovermode="x")
fig

# %% [markdown]
# 左圖只顯示前十二小時，但密度以全部間隔正規化；右圖也保留觀察期間中
# 沒有事件的日期。很短的間隔較多、每日計數起伏較大，表示固定平均率沒有充分描述這份目錄。
# 但這不是僅靠一張圖就能證明事件互相觸發：研究區內的背景差異、活動率變動與
# 偵測能力變化，也可能影響圖形。圖的作用是指出需要解釋的現象，讓我們知道
# 接下來要比較哪些模型。
#
# 地震學常用兩組經驗關係整理這些現象。一組描述事件之後活動如何隨時間衰減，
# 另一組描述小地震與大地震的數量比例。它們分別會帶到Omori–Utsu律與
# Gutenberg–Richter律。先不用背公式，可以把它們理解成「有多活躍」和
# 「各種大小佔多少」兩個互補問題。模型還需要說明位置、歷史更新及不確定性，
# 並不是把兩條曲線相乘就完成全部預報。
#
# ## 1.3 一份預報需要說清楚什麼
#
# 假設我們關心某個固定區域未來七天的活動。首先可以報出期望事件數，它是
# 在模型下重複許多次這種情境的平均，不保證單次觀察會等於這個數字。
# 也可以報出至少一次的機率，或事件數可能落在哪一個範圍。
# 這些量彼此相關，但不相同，尤其不能只知道平均數就唯一決定所有機率。
#
# 在Poisson計數假設下，期望數為 $\Lambda$ 時，至少一次的機率是
# $P(N\ge1)=1-e^{-\Lambda}$。這是有條件的換算，不是所有地震模型通用的
# 捷徑。後面自激發模型允許新事件繼續產生後續活動，未來機率常需要模擬或
# 其他計算，不能只把目前的率固定住。
#
# 再往使用者的問題靠近一步，就要問某個場址是否會晃到指定程度。相同規模
# 事件發生在不同距離、深度與傳播環境，地動可能不同，因此還需要地動模型。
# 這是後面PSHA從事件活動連到場址危害的動機。
#
# 下圖用一個給定衰減率的非齊次Poisson示例展示這個關係。它不是花蓮資料的
# 擬合，也不是正式預報。這裡假設每次事件獨立地有12%的機率使場址地動超越
# 門檻，則超越事件的期望數是 $0.12\Lambda$，再換算為至少一次的機率
# $1-e^{-0.12\Lambda}$。要縮小的是事件率，不是把原本「至少一次」的機率乘以12%。
#
# %% tags=["remove-input"]
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
                                    "④ 場址至少一次超越機率"))
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
# 上排先描述活動率與七天內的平均次數，下排則轉成事件與場址地動的機率。
# 期望數可以大於一，機率則不會；只有期望數很小時，Poisson的至少一次機率
# 才近似等於期望數。兩種量的縱軸也不相同，不能把曲線的高度直接當成同一個數字。
#
# 這個差別對理解預報很重要。相對於很小的背景機率增加許多倍，絕對機率仍然
# 可能很小；反過來，活動逐漸下降也不表示已經沒有風險。預報應交代比較基準、
# 目標期間與不確定性，讓讀者知道數字變動究竟代表什麼。
#
# ## 1.4 為什麼還需要檢驗與制度
#
# 模型能重現已經看過的曲線，是有用的第一步，卻可能部分來自針對這段資料
# 調整了許多設定。若要判斷它能否用於未來，就要事先訂好預報及評分方式，
# 保留未參與選模的資料，並與合理基準比較。
#
# CSEP的前瞻測試工作提供了這類共同評估的架構。預先固定規則能減少看到結果
# 後調整說法的空間，但不是偏誤自動消失的保證：反覆挑選測試結果、忽略失敗的
# 預報或只展示有利區域，仍然會讓評估失真。後面會把一致性、相對表現與實際
# 決策價值分開，逐步看清它們各自回答的問題。
#
# 當預報需要持續更新、及時發布並供社會使用，就進入作業化地震預報的範圍。
# 除了模型，還需要穩定資料、可追溯的更新程式及清楚的溝通。Jordan等人的
# 報告與Mizrahi等人的綜述將這些環節放在一起討論；我們會在模型與檢驗學完後
# 回到這裡，而不在第一章要求你同時記住整套系統。
#
# ## 1.5 先理解資料，再理解模型
#
# 後面的學習會沿著同一個問題持續往前。先從目錄的計數、等待時間與位置看見
# 現象，接著用隨機模型建立比較基準，再引入條件強度，把歷史資訊轉成可更新的
# 發生率。完整度、規模分布與叢集律提供建模材料，估計方法則說明資料能讓我們
# 知道多少。到了ETAS與EEPAS，我們就能比較兩種不同的預報想法，而不只記住公式。
#
# 學會模型之後，再問它對未來有沒有幫助、如何組合、如何連到地動危害與公開
# 資訊。第二部回到臺灣的地下水、地磁、地動與GNSS，把相同的統計判斷用在
# 不同物理觀測上。長推導放在第三部，正文會保留理解所需的結果與限制。
#
# 現在先把模型名稱放在一旁，從{doc}`foundation_catalog`的一列地震資料開始。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Can you predict earthquakes?](https://www.usgs.gov/faqs/can-you-predict-earthquakes) — USGS；免費官方說明。先讀這份短問答，釐清「預測特定地震」與「估計發生機率」的差別，再回看本章的三種用語。
# - [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS；免費官方說明。用實際預報介面解釋期望次數、規模門檻與時間窗，適合對照本章「一份預報長什麼樣」。
# - [Operational Earthquake Forecasting: State of Knowledge and Guidelines for Utilization](https://doi.org/10.4401/ag-5350) — Thomas H. Jordan 等，2011，*Annals of Geophysics*；[免費出版版全文](https://www.research.ed.ac.uk/files/10773224/AnnalsofGeoPhys.pdf)。這是 OEF 的重要報告，建議先讀摘要與建議，瞭解為什麼機率模型、公開溝通與權責制度必須一起考慮。
#
# - [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等，2024，*Reviews of Geophysics*；免費開放全文。從第 2 節的領域導覽讀起，再比較義大利、紐西蘭與美國系統，瞭解本章所說的基準比較、前瞻檢驗與透明溝通如何進入實務。
