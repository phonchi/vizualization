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
# # 11. 讓率隨歷史更新：條件強度與點過程概似
#
# 一顆地震在今天下午發生，明天的率就可能改變。第 10 章只替一顆主震畫衰減曲線，本章把同樣的想法掛到每顆事件上。模型在每次收到新事件後，都能重新計算下一刻的率。
#
# **歷史**（history） $H_t$ 是時刻 $t$ 以前已發生的事件集合。它包含時間、位置與規模；不包含當下才要評估的事件。若把事件本身提前放進歷史，模型就會用答案解釋答案。
#
# ## 11.1 下一小段時間，率是多少？
#
# **條件強度**（conditional intensity）是已知歷史後，下一小段時間的即時事件率。暫時不分位置與規模，可以寫成
#
# $$P\{N(t,t+dt]=1\mid H_t\}=\lambda^*(t)dt+o(dt).$$
#
# $dt$ 是很短的時間，$o(dt)$ 代表比 $dt$ 更快縮小的餘項。星號提醒讀者：這個率依賴歷史。它的單位是每日事件數，數值可以大於一；只有乘上一小段時間後，才近似該段的單事件機率。
#
# 同樣的定義可以延伸至位置與規模。$\lambda^*(t,x,y,m)$ 的單位是每一天、每平方公里、每規模單位的事件數。給定固定的密度曲面，須積分到指定格箱，才得到第 3 章規格卡上的期望數。
#
# ## 11.2 固定時間表與事件回饋
#
# 下圖用三個人為指定的事件，示範兩種率。灰線是一開始就排好的時間表。藍線則在事件發生後加入新的衰減項。這三個時間是合成示意，沒有從模型抽樣，也沒有擬合義大利資料。
#
# **自激發**（self-exciting）指事件會增加後續事件的條件率。對照之下，**非齊次 Poisson 過程**（inhomogeneous Poisson process）使用事先指定、可隨時間改變的率。率曲線都會起伏，但只有前者會因新事件而更新。
#
# ### 用一個模型名字連起來：Hawkes 過程
#
# **線性 Hawkes 過程**（linear Hawkes process）把背景率和每顆過去事件的觸發貢獻相加：
#
# $$\lambda^*(t)=\mu+\sum_{t_i<t}\phi(t-t_i),\qquad \mu>0.$$
#
# $\mu$ 是背景每日率；$\phi(\tau)$ 是單一事件在經過 $\tau$ 天後增加的率。本章使用非負、因果的核：$\phi(\tau)\ge0$，且在 $\tau\le0$ 時取零。因果在這裡指只往未來提供貢獻，並不等於已證明地震間的物理因果。
#
# 「線性」指各事件的貢獻直接相加。一顆事件加一份核，兩顆就加兩份；核本身可以是彎曲的指數或冪律。不是所有依賴歷史的點過程都是這個模型：若以非線性函數轉換總貢獻，或讓事件降低率，就需要不同形式與條件。
#
# 下圖的藍線就是把指定事件時間代入這個公式。核取 $\phi(\tau)=1.5e^{-\tau/1.3}$（$\tau>0$）；背景為每天 $0.3$。因此每顆事件剛出現後，新增率約為每天 $1.5$，再逐漸消退。事件時間是手動指定的，這張曲線沒有模擬 Hawkes 目錄。
#
# 核的面積還有另一個意思：$n=\int_0^\infty\phi(\tau)\,d\tau$ 是每顆事件的平均直接後代數。此圖的面積是 $1.5\times1.3=1.95$，所以這組視覺示範不是平穩、次臨界的參數設定。有限歷史的曲線仍可畫出；不能據此宣稱它代表穩定長期地震率。
#
# 若背景固定，非負核可積且 $n<1$，標準線性Hawkes具有有限均值的平穩版本。第12章將用分支過程解釋這個條件，並把單一核擴成與規模有關的ETAS。這裡先記住：核高度描述即時增加，核面積描述整個未來的直接後代平均數。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.optimize import minimize
from IPython.display import HTML, display
from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, apply_layout, show_diagram
te=np.array([2.,5.,5.8])
t=np.linspace(0,12,1201)
base=.3
rate=np.full_like(t,base)
for ti in te:
    rate += np.where(t>ti,1.5*np.exp(-np.maximum(t-ti,0)/1.3),0.)
fig=go.Figure()
fig.add_scatter(x=t,y=.5+.25*np.sin(t/2),name='事先指定的率',line_color='#64748b')
fig.add_scatter(x=t,y=rate,name='依已發生事件更新',line_color=ACCENT)
for ti in te:
    fig.add_vline(x=ti,line_dash='dot',line_color='#64748b')
apply_layout(fig,title='合成示意：同一條時間軸，兩種更新規則',xaxis_title='時間（天）',yaxis_title='事件數／天')
fig

# %% [markdown]
# 事件出現以前，藍線不能先上升。若將 $t_i=t$ 的事件納入求和，曲線會在評估事件時使用自身的貢獻。正確做法是以左側已知歷史計算 $\lambda^*(t_i)$，收到事件後才更新下一刻的率。
#
# 圖中的時間核用指數衰減，方便看清楚跳升與回落。地震預報常改用第 10 章的 Omori 核。選哪個核，是關於記憶能持續多久的假設；條件強度這個概念本身沒有指定核的形狀。
#
# %% tags=["remove-input"]
show_diagram("d11_hawkes_process", caption="線性 Hawkes：非負觸發核沿已發生事件相加；每次新事件只改變其後的率。")

# %% [markdown]
# ## 11.3 一顆顆事件的分數，還要扣掉什麼？
#
# 假設模型在每個已發生事件處都給很高的率。這樣的模型是否一定好？若把全天的率都提高，它也能做到。因此對數概似必須同時獎勵命中位置，並扣掉整段預期活動量。
#
# 在觀測窗 $[S,T]$，時間點過程的對數概似為
#
# $$\ell=\sum_{i:S<t_i\le T}\ln\lambda^*(t_i)-\int_S^T\lambda^*(u)\,du.$$ (eq:point-process-loglik)
#
# 第一項在事件時間取值。第二項在整個觀測窗累積。沒有事件的時段也會增加第二項；模型若在安靜時段給很高的率，就會付出代價。第 5 章 Poisson 概似的「事件項減去期望數」，在這裡有了時間解析度。
#
# 率的量綱看似讓取對數不容易理解。實際推導先把時間切成小格，對每格發生或未發生寫機率，再令格寬縮小。與參數無關的格寬常數在比較時消去。比較兩個模型仍須使用相同時間單位與觀測範圍。
#
# 位置與規模一起評估時，事件項代入完整的 $(t_i,x_i,y_i,m_i)$，積分項也要涵蓋相同的時間、空間與規模範圍。收集區外或觀測窗前事件若已知，可以貢獻歷史；只有評估集合內事件才放進事件項。
#
# ## 11.4 用模型的時鐘看活動
#
# **補償子**（compensator）是沿已知歷史累積的條件率：
#
# $$A(t)=\int_S^t\lambda^*(u)\,du.$$
#
# 當模型認為地震頻繁，$A(t)$ 走得快；模型認為安靜時，它走得慢。這個新時鐘能把快慢不同的時段放在同一尺度。若指定模型正確且符合時間變換條件，事件在新時鐘中的間隔應像平均一的指數分布。
#
# 下圖回到 L'Aquila 2009 的同一個 50 km 範圍。仍使用規模至少 3.0、主震後 1–60 天的事件。先用單主震 Omori 模型估參數，再將實際累積計數與 $A(t)$ 疊在一起。這是該簡化模型的診斷，還沒有建立完整自激發模型。

# %% tags=["remove-input"]
cat=italy.experiment_catalog()
x,y=italy.lonlat_to_km(13.38,42.35)
near=np.hypot(cat.x_km-x,cat.y_km-y)<=50
cand=cat.loc[near & (cat.time>=pd.Timestamp('2009-04-06')) & (cat.time<pd.Timestamp('2009-04-07'))]
main=cand.loc[cand.mb.idxmax()]
tau=(cat.time-main.time).dt.total_seconds().to_numpy()/86400
keep=(np.hypot(cat.x_km-main.x_km,cat.y_km-main.y_km)<=50)&(cat.mb>=3)&(tau>=1)&(tau<60)
ts=np.sort(tau[keep])
def integral(a,b,c,p):
    if abs(p-1)<1e-7:
        return np.log((b+c)/(a+c))
    return ((b+c)**(1-p)-(a+c)**(1-p))/(1-p)
def nll(z):
    c,p=np.exp(z[0]),z[1]
    return p*np.log(ts+c).sum()+len(ts)*np.log(integral(1,60,c,p))
fits=[minimize(nll,[np.log(c),p],method='L-BFGS-B',bounds=[(np.log(1e-4),np.log(10)),(.3,2.5)]) for c,p in [(.01,.8),(.1,1.2),(1.,1.5)]]
valid=[r for r in fits if r.success and np.isfinite(r.fun)]
if not valid:
    raise RuntimeError('診斷模型未收斂')
fit=min(valid,key=lambda r:r.fun)
c,p=np.exp(fit.x[0]),fit.x[1]
tgrid=np.linspace(1,60,400)
a=len(ts)/integral(1,60,c,p)
fig=go.Figure()
fig.add_scatter(x=np.r_[1,ts,60],y=np.r_[0,np.arange(1,len(ts)+1),len(ts)],line_shape='hv',name='HORUS 累積計數',line_color='#64748b')
fig.add_scatter(x=tgrid,y=a*integral(1,tgrid,c,p),name='擬合模型補償子',line_color=ACCENT)
apply_layout(fig,title=f"L'Aquila：{len(ts)} 顆事件的累積診斷",xaxis_title='距主震（天）',yaxis_title='累積事件數／累積率')
fig

# %% [markdown]
# 兩條曲線在終點相遇，是因為振幅由總數估得。中途若出現系統性差距，才顯示模型沒有描述的時間結構。即使線條接近，也只檢查了累積時間形狀，沒有檢查位置、規模或間隔獨立性。
#
# 同一份資料既估參數又診斷，參考分布也受到估計影響。不能直接用完全已知參數的理論分布宣稱顯著。若需要正式檢定，可從模型模擬，再對每份目錄重估參數，重做相同診斷。詳細條件見附錄 A。
#
# ## 11.5 發報當下與事後評分
#
# 2016 年某個發報時刻記為 $t_0$。當時只知道 $H_{t_0}$，還不知道窗內會出現哪些事件。若模型允許新事件繼續觸發，未來率本身也是隨機的。計算完整預報時，必須平均這些可能歷史。
#
# 等到預報窗結束，實際事件都已知。可以沿著實際歷史，每次用事件前的率計算概似。這時積分得到的補償子依賴窗內事件，通常不是發報當下就知道的固定數字。
#
# 兩個量因此要分清楚：發報期望數是給定起點歷史、對未來可能路徑取平均；事後補償子是沿實際路徑累積率。它們有相關的理論關係，但不能逐次視為同一個數值。
#
# 本書第 13 章採用背景加已知歷史第一代的 ETAS 教學近似。固定發報時刻後，不在窗內補進未來後代。因此輸出比較容易計算，同時少了完整級聯。這是特定實作的近似，並不是條件強度的定義限制。
#
# ### 率如何改變，決定是哪一種點過程
#
# Hawkes描述事件之後增加活動，**更新過程**（renewal process）則常描述同一斷層反覆破裂。後者只記距離上次指定事件多久；每次事件把年齡歸零。指數等待時間的特例是齊次Poisson，其他等待分布可有不同的年齡效應。
#
# **Cox 過程**（Cox process）以外生隨機率生成條件Poisson事件。它也會形成叢集，例如共同環境使某段期間的率偏高。觀察到成群事件，因而不足以單憑圖形認定Hawkes觸發。事件歷史可以更新對未知外生率的推測，這與事件本身回饋生成機制要分開。
#
# **自我修正過程**（self-correcting process）提供相反的回饋示例：無事件時率逐漸上升，事件發生後率下降。應力釋放模型把這個想法連到持續載入與地震釋放，但需另外指定物理對應。附錄A比較各模型使用哪些資訊，附錄F給出具體公式；這些模型不是互相可替換的命名方式。
#
# ## 11.6 把一段方法文字翻成操作
#
# 論文若寫「以條件強度最大化概似」，可以依序找三件事。首先找事件求和的集合：評估哪些日期、位置與規模？再找歷史集合：窗前與區外事件是否提供貢獻？最後找積分範圍：是否與觀測事件的集合一致？
#
# 例如窗前一顆地震不在事件求和內，仍可增加窗內強度及積分。反過來，區內一段沒有事件的時間也要保留積分。這兩個細節都會改變估計，不能因為資料表沒有對應的一列，就把它們省略。
#
# 這種讀法讓公式直接對回資料處理。即使還不能自行推導全部點過程理論，也能判斷一個概似是否使用未來資訊、是否漏掉邊界來源，以及作者究竟在估哪一個範圍的活動。
#
# ## 本章填入的規格欄位
#
# 本章沒有新增預報參數。需要記住的是更新順序：先用過去評估當下，再收進新事件。下一章會把背景、產能、時間、空間與規模組成完整公式；每個部分都能對回本章的率與積分。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([('歷史條件','Hₜ 只含 t 以前的事件'),('預報輸出','固定發報歷史後的格箱期望數'),('診斷輸出','沿實際歷史累積的補償子'),('模型參數','待填（第 12–13 章）'),('正式檢驗','待填（第 16–18 章）')])))

# %% [markdown]
# 下一章：{doc}`12_etas_structure`。推導與診斷補充見 {doc}`appendix_a_point_process`。
#
# 例如某段十天沒有事件，概似仍包含這十天的積分。
# 模型若在這段時間持續給很高的率，就會付出較大的積分代價。
# 事件項與積分項一起計算，才能讓有地震與沒有地震的時段都提供資訊。
#
# ## 參考資料與延伸閱讀
#
# - Reinhart（2018），[自激發時空點過程回顧](https://arxiv.org/abs/1708.02647)。免費作者稿；先讀條件強度與分支表示，對照兩種率的差別。
# - Jalilian（2019），[ETAS 套件論文](https://www.jstatsoft.org/article/view/v088c01)。免費全文；對照歷史、觀測區域及診斷流程。
# - Ogata（1988），[地震統計模型與點過程殘差](https://doi.org/10.1080/01621459.1988.10478560)。理解補償子如何成為模型檢查工具。
