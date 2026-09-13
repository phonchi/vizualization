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
# # 14. 大地震之前，小地震先變多變大？Ψ 現象
#
# 一張規模對時間的圖，可以同時看疏密與高度。
# 第 13 章關注地震之後的觸發。
# 本章把目光移到大地震之前。
# 部分回溯研究曾觀察到活動尺度增加。
#
# Ψ 是前兆尺度增加現象的簡稱。
# 英文為 precursory scale increase phenomenon。
# 它描述選定區域內，地震率與規模一起提高。
# 這項觀察後來啟發 EEPAS 模型。
# 本章先說清楚量法與證據的用途。

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
# ## 14.1 先看一份刻意改變活動的目錄
#
# 下圖的前半段事件較疏，規模也較小。
# 後半段由程式刻意增加事件率與規模。
# 合成目錄（synthetic catalogue）是依設定生成的資料。
# 它讓讀者知道變化確實由何處加入。
#
# 一條累積曲線能把兩種資訊放在一起。
# 累積規模異常（cumulative magnitude anomaly）記為 $C(t)$。
# 每顆事件先增加一份規模超額。
# 曲線再扣除整段期間的平均累積趨勢。
#
# $$
# \begin{aligned}
# C(t)&=\sum_{t_s<t_i\le t}[m_i-(m_a-0.1)]-k(t-t_s),\\
# k&=\frac{\sum_{t_s<t_i<t_f}[m_i-(m_a-0.1)]}{t_f-t_s}.
# \end{aligned}
# $$
#
# $m_a$ 是本章的分析門檻（analysis threshold）。
# 程式只保留達到這個門檻的事件。
# 此符號另設，以免混同目錄完整度 $M_c$。
# $0.1$ 沿用原研究的規模刻度。
# 它讓恰在門檻上的事件也有正貢獻。
#
# 上式用 $[t_s,t_f]$ 表示選定的分析期間。
# 本章目錄排除兩個端點上的事件。
# 因此，兩個求和在終點有相同總量。
# $k$ 的單位是每天的規模累積量。
# $C(t)$ 的單位是規模單位。
#
# 一顆事件到來時，曲線向上跳。
# 兩顆事件之間，曲線以固定斜率下降。
# 持續上升表示累積速度高於全段平均。
# 較多事件或較大事件都能拉高曲線。
# 這條曲線累積規模超額，並非能量。

# %% tags=["remove-input"]
rng = np.random.default_rng(140913)
t0, change, tf, ma = 0., 600., 1000., 2.5
n_before, n_after = 60, 100
t = np.r_[rng.uniform(0,change,n_before),rng.uniform(change,tf,n_after)]
m = np.r_[ma+rng.exponential(0.22,n_before),ma+0.3+rng.exponential(0.3,n_after)]
order = np.argsort(t); t,m = t[order],m[order]
def cumulative_anomaly(times,magnitudes,start,end,threshold):
    excess = np.asarray(magnitudes)-threshold+0.1
    cumulative = np.r_[0,np.cumsum(excess)]
    grid = np.r_[start,np.repeat(times,2),end]
    steps = np.r_[0,np.column_stack((cumulative[:-1],cumulative[1:])).ravel(),cumulative[-1]]
    slope = cumulative[-1]/(end-start)
    values = steps-slope*(grid-start)
    assert np.allclose(values[[0,-1]],0,atol=1e-10)
    return grid,values,slope
grid,cv,slope = cumulative_anomaly(t,m,t0,tf,ma)
onset = grid[np.argmin(cv)]
fig = make_subplots(rows=2,cols=1,shared_xaxes=True,subplot_titles=("事件時間與規模","扣除全段平均後的 C(t)"))
fig.add_trace(go.Scatter(x=t,y=m,mode="markers",name="合成事件",marker_color=ACCENT),row=1,col=1)
fig.add_trace(go.Scatter(x=grid,y=cv,name="C(t)",line_color=ACCENT),row=2,col=1)
fig.add_vline(x=change,line_dash="dash",line_color=PALETTE[2])
fig.update_yaxes(title_text="規模",row=1,col=1)
fig.update_yaxes(title_text="規模累積量",row=2,col=1)
fig.update_xaxes(title_text="天",row=2,col=1)
apply_layout(fig,height=510,showlegend=False)
display(Markdown(f"合成資料共 {len(t)} 顆。程式在第 {change:.0f} 天改變活動。C(t) 最低點在第 {onset:.1f} 天；兩個日期來自不同的定義。"))
fig

# %% [markdown]
# ## 14.2 曲線最低點能告訴什麼
#
# 圖上的最低點是一個候選起點。
# 候選起點（candidate onset）是量法提出的分段位置。
# 它由整段曲線決定，仍須搭配接受條件。
# 原研究還檢查活動率與規模差異。
# 單憑最低點不足以完成 Ψ 辨識。
#
# $C(t_s)=C(t_f)=0$ 是公式的結果。
# 因為扣除的直線剛好連接總累積量。
# 因此，尾端回到零是作圖方式的安排。
# 讀者應把注意力放在中間的相對變化。
#
# 終點 $t_f$ 也參與平均斜率 $k$ 的計算。
# 若今天尚未到達終點，$k$ 仍可能改變。
# 新事件加入後，最低點的位置也可能移動。
# 所以，事後畫好的曲線帶有完整時窗資訊。
# 它適合描述已發生的活動。

# %% [markdown]
# ## 14.3 把一段活動整理成三個尺度
#
# 一組候選事件可以整理成三個數字。
# 前兆規模（precursor magnitude）記為 $M_P$。
# 它取選定前兆群中最大三顆的平均規模。
# 主震規模（mainshock magnitude）記為 $M_m$。
# 主震在此是回溯分析選定的較大事件。
#
# 前兆時間（precursor time）記為 $T_P$。
# 它是候選起點到主震的時間差，以天計。
# 前兆面積（precursor area）記為 $A_P$。
# 它是依搜尋規則圈出的活動面積，以平方公里計。
#
# 三個最大規模只摘要較大的前兆事件。
# $M_P$ 並未使用全部規模的平均。
# 事件數增加時，最大幾顆也可能隨之增大。
# 因此，規模摘要與取樣數量要一起理解。
#
# 一個矩形量法會包住指定的事件群。
# 原始定義也可納入主震與餘震的分布。
# 圓形量法則用搜尋半徑決定面積。
# 兩種量法得到的 $A_P$ 未必相同。
# 閱讀尺度關係前，須先核對圈選規則。
#
# 跨多個地震序列，可以分別畫三張散布圖。
# 尺度關係（scaling relation）連結不同大小的量。
# 常見的迴歸形式如下。
#
# $$
# \begin{aligned}
# M_m&=a_m+b_mM_P+\epsilon_m,\\
# \log_{10}(T_P)&=a_t+b_tM_P+\epsilon_t,\\
# \log_{10}(A_P)&=a_a+b_aM_P+\epsilon_a.
# \end{aligned}
# $$
#
# $\epsilon$ 是殘差（residual），即點與迴歸線的差。
# 同一個 $M_P$ 可以對應不同的主震規模。
# 時間與面積也保留這種散布。
# 下一章會把散布轉成機率密度。
#
# 三條式子的係數用小寫下標區分。
# 它們描述已圈選的事件群。
# EEPAS 的係數改用單顆事件建立預報。
# 兩者概念相連，估計對象卻不同。
# 因此，本章不直接搬入 EEPAS 的參數。

# %% tags=["remove-input"]
mp = np.linspace(3.0,5.2,36)
responses = [1.1+mp+rng.normal(0,.23,len(mp)),
             1.6+.45*mp+rng.normal(0,.22,len(mp)),
             1.2+.65*mp+rng.normal(0,.25,len(mp))]
labels = ["主震規模 M_m","log₁₀（T_P／天）","log₁₀（A_P／km²）"]
fig = make_subplots(rows=1,cols=3,subplot_titles=("規模","時間","面積"))
for col,(yy,label) in enumerate(zip(responses,labels),1):
    coef = np.polyfit(mp,yy,1)
    fig.add_trace(go.Scatter(x=mp,y=yy,mode="markers",marker_color=ACCENT,showlegend=False),row=1,col=col)
    fig.add_trace(go.Scatter(x=mp,y=np.polyval(coef,mp),line_color=PALETTE[2],showlegend=False),row=1,col=col)
    fig.update_xaxes(title_text="M_P",row=1,col=col)
    fig.update_yaxes(title_text=label,row=1,col=col)
apply_layout(fig,height=390,margin=dict(l=60,r=20,t=55,b=60))
display(Markdown(f"三組合成迴歸各有 {len(mp)} 個點。點與係數皆為示意，未使用義大利研究的迴歸結果。"))
fig

# %% [markdown]
# ## 14.4 橫軸也有誤差，斜率會怎樣
#
# $M_P$ 只取幾顆事件，因此帶有取樣變動。
# 測量誤差也會改變橫軸位置。
# 回歸稀釋（regression dilution）是誤差使斜率趨平的現象。
# 它在古典的獨立加性誤差條件下成立。
#
# 設觀測橫軸為 $X=X^*+u$。
# $X^*$ 表示原來的尺度，$u$ 表示獨立誤差。
# 若真實關係的斜率為 $b_*$，則
#
# $$
# b_{\mathrm{obs}}\approx b_*
# \frac{\operatorname{Var}(X^*)}
# {\operatorname{Var}(X^*)+\operatorname{Var}(u)}.
# $$
#
# 分母多了一份誤差變異。
# 因此，觀測斜率在這組假設下較接近零。
# 若圈選造成誤差相依，結果可能不同。
# 迴歸線的斜率須連同量測與選樣方式判讀。
#
# 時間軸上的十倍差異，在對數軸上相差一格。
# 這讓很長與很短的前兆時間可放在同一張圖。
# 把天改成年會平移截距。
# 因此，公式旁的單位也是模型設定的一部分。

# %% [markdown]
# ## 14.5 L’Aquila：一張明標回溯的圖
#
# 下圖以 2009 年 L’Aquila 事件為中心。
# 中心位置由那顆已知的大事件決定。
# 分析期間從 2003 年起，到該事件之前。
# 半徑與門檻固定後，再畫出目錄與曲線。
#
# 回溯分析（retrospective analysis）使用已知結局整理過去。
# 這張圖正是回溯示意。
# 本章未以此圖宣稱完成前瞻辨識。
# 曲線上的起伏只描述這個選定範圍。
#
# 一個圓圈會排除圈外事件。
# 若把半徑放大，加入的事件會改變平均斜率。
# 若把起始日期前移，曲線也會重新計算。
# 這些設定都要與圖一起保留。

# %% tags=["remove-input"]
cat = italy.experiment_catalog()
main_candidates = italy.target_events(cat,"learning")
main = main_candidates.loc[(main_candidates.time>=pd.Timestamp("2009-04-06"))&(main_candidates.time<pd.Timestamp("2009-04-07"))].sort_values("mb").iloc[-1]
radius = 50.
start = italy.decimal_days(pd.Timestamp("2003-01-01"))
end = float(main.t_days)
distance = np.hypot(cat.x_km-main.x_km,cat.y_km-main.y_km)
local = cat.loc[(distance<=radius)&(cat.t_days>start)&(cat.t_days<end)&(cat.mb>=2.5)].sort_values("t_days")
grid,cv,_ = cumulative_anomaly(local.t_days.to_numpy(),local.mb.to_numpy(),start,end,2.5)
fig = make_subplots(rows=2,cols=1,shared_xaxes=True,subplot_titles=("L’Aquila 周圍事件：回溯示意","C(t)：回溯示意"))
for mask,color,label in [(~(local.in_R & (local.mb>=5)),ACCENT,"輸入事件"),(local.in_R & (local.mb>=5),QUAKE_COLOR,"目標事件")]:
    fig.add_trace(go.Scatter(x=local.loc[mask,"time"],y=local.loc[mask,"mb"],mode="markers",marker=dict(color=color,size=4),name=label),row=1,col=1)
fig.add_trace(go.Scatter(x=[main.time],y=[main.mb],mode="markers",marker=dict(color=QUAKE_COLOR,size=11,symbol="star"),name="選定的目標主震"),row=1,col=1)
dates = pd.Timestamp("1960-01-01")+pd.to_timedelta(grid,unit="D")
fig.add_trace(go.Scatter(x=dates,y=cv,line_color=ACCENT,name="C(t)"),row=2,col=1)
fig.update_yaxes(title_text="取整 Mw",row=1,col=1)
fig.update_yaxes(title_text="規模累積量",row=2,col=1)
fig.update_xaxes(title_text="日期",row=2,col=1)
apply_layout(fig,height=540,legend=dict(orientation="h",y=-.18))
display(Markdown(f"回溯示意：中心為 {main.time:%Y-%m-%d}、Mw {main.mb:.1f} 的事件。半徑 {radius:.0f} km，門檻 {2.5:.1f}，主震前共 {len(local)} 顆。紅星未納入 C(t)。"))
fig

# %% [markdown]
# ## 14.6 從漂亮案例走到可檢驗規則
#
# 圖中的一個最低點可以提出候選時段。
# 前瞻檢驗（prospective test）先凍結規則，再等待新資料。
# 凍結內容包括搜尋區域、時間長度與門檻。
# 每次符合規則的候選都要保留。
#
# 一份候選清單也須保留沒有大地震的案例。
# 只展示後來發生大地震的案例，
# 會看不到候選出現後落空的頻率。
# 預報實驗正是把這些結果放回同一份帳本。
#
# 兩個搜尋窗可能圈到同一批密集事件。
# 較大的面積可能搭配較短的時間。
# 時空抵換（space–time trade-off）描述這種補償關係。
# 它提醒讀者，辨識出的尺度可能依賴搜尋方式。
#
# 一份完整紀錄應列出所有搜尋設定。
# 原始目錄版本、規模取整也應保存。
# 若更換設定後才挑出清楚的圖，
# 新設定就需要另一批資料來檢驗。
#
# 合成圖已知改變活動的真正起點。
# 義大利圖只知道觀測到的事件。
# 兩張圖的證據角色因而不同。
# 前者檢查量法如何回應已知改變。
# 後者展示量法在真實資料上的樣子。
#
# 三條尺度迴歸提供一個可延伸的想法。
# 小事件可以對未來多種尺度給出少量貢獻。
# 下一章把每顆事件的貢獻攤成機率密度。
# 模型再把所有貢獻加成一份網格預報。

# %% [markdown]
# ## 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
("本章規格更新", "無；本章建立 Ψ 的量法與證據角色"),
("真實資料圖", "L’Aquila 回溯示意；未作 Ψ 前瞻辨識"),
("EEPAS 三核與參數", "待填（第 15 章）"),
("數量檢驗", "待填（第 16 章）"),
("位置與規模檢驗", "待填（第 17 章）"),
])))

# %% [markdown]
# 每顆事件如何分配未來貢獻，接在下一章：
# {doc}`EEPAS 在義大利 <15_eepas_italy_forecast>`。
#
# 累積規模曲線的水平段表示沒有新的入選事件。
# 向上的階梯來自新增事件，其高度取決於採用的規模權重。
# 先逐段讀出這兩種變化，再與時間、位置及搜尋條件對照，才能看懂回溯圖。
#
# ## 參考資料與延伸閱讀
#
# - Rhoades 等（2022），[EEPAS 的二十年回顧](https://doi.org/10.3390/geosciences12090349)。免費全文。第 2 節介紹 Ψ 與尺度關係。
# - Christophersen、Rhoades 與 Hainzl（2024），[Ψ 的演算法辨識](https://doi.org/10.1785/0220240233)。[免費機構全文](https://gfzpublic.gfz.de/rest/items/item_5029405_4/component/file_5029659/content)。式 1–2 定義累積曲線。再讀搜尋條件與對照設計。
# - Rastin 等（2021），[前兆活動的時空抵換](https://doi.org/10.3390/app112110215)。免費全文。閱讀面積與時間如何互相補償。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。[免費機構全文](https://www.earth-prints.org/handle/2122/17084)。表 3 提供參數。附錄列出模型公式。
