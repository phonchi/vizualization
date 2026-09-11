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
# # 2. 從地震目錄提出統計問題
#
# 上一章把預報寫成一個有區域、期間與規模門檻的問題。接下來先不選模型，
# 而是想像你收到一份地震清單：每一列有時間、經緯度、深度和規模。
# 這些欄位看起來很熟悉，卻還沒有告訴我們什麼是一次「增加」、一群「餘震」，
# 或一段值得注意的「安靜」。這些判斷都需要比較的對象與觀察的尺度。
#
# 我們會沿用一份合成目錄，練習把同一份資料問成不同的問題。這份資料刻意安排了
# 若干事件群，位置與規模也是模擬的；它的用途是讓每一步可以看清楚，並不代表
# 臺灣某場地震。後面先用它比較隨機性與叢集，再把問題帶到模型和真實觀測。
#
# ## 2.1 一列資料如何成為一次事件
#
# 連續地震波形是儀器隨時間記錄的振動；地震目錄則是經過偵測、定位、規模估計
# 與整理後的事件清單。兩者的關係很像錄音與逐字稿：逐字稿方便分析，卻是經過
# 判讀的產品。微弱事件可能漏掉，兩筆紀錄可能需要合併，定位與規模也可能更新。
# 因此「目錄沒有事件」至少有兩種解釋：那裡沒有發生符合條件的地震，或它沒有
# 被這套觀測與處理流程記錄下來。
#
# 做統計前先界定一次分析的範圍。時間窗要同時包含有事件與沒有事件的時段；
# 空間區域不能看到震央後才沿著它們描邊；規模要使用已宣告的尺度。
# $M_L$ 與 $M_w$ 都描述地震大小，但定義與估計方式不同，不能因為數字相近
# 就不加處理地混在一起。這些選擇決定我們究竟在研究哪一群事件。
#
# 下表是共用合成目錄的前幾筆。時間從觀察開始算起，位置使用平面公里座標，
# $m$ 是教學用規模標記，沒有套用到任何真實地區的規模轉換。
#
# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from gdms_toolkit.teaching import learning_catalog
from gdms_toolkit.viz import setup_plotly, apply_layout, ACCENT, PALETTE, QUAKE_COLOR
setup_plotly()
catalog = learning_catalog()
catalog.head(6).round(3).rename(columns={"day": "時間（天）", "x_km": "東西（km）", "y_km": "南北（km）", "m": "規模標記"})

# %% [markdown]
# ## 2.2 同一份資料，可以有三種讀法
#
# 如果我們關心「最近是不是比較活躍」，可以沿時間累計事件數；若關心「事件是否
# 接連發生」，就看相鄰事件的間隔；若要問「較大的事件集中在哪裡」，則要同時
# 看位置與規模。這幾種圖各自保留不同資訊，可以一起閱讀。
#
# 在下圖左上，累積曲線每遇到一個事件就上升一階。陡的區段表示事件密集，
# 平的區段表示暫時沒有新事件。右上把事件位置保留下來，圓點大小表示規模；
# 這份合成資料的位置與規模獨立於時間抽取，刻意沒有加入和時間群相連的空間群。
# 圖中央較密來自設定的位置密度，並不是事件彼此觸發的證據；真實序列則可能
# 同時包含時空相依。因此不能只憑時間聚集，就推論空間上也有對應的事件群。
# 下方再把時間壓成每段期間的事件數，觀察我們因此得到和失去了什麼。
#
# %% tags=["remove-input"]
fig = make_subplots(rows=2, cols=2, subplot_titles=("事件逐次累積", "位置與規模標記", "每 2 天計數", "每 15 天計數"))
fig.add_trace(go.Scatter(x=np.r_[0,catalog.day,120], y=np.r_[0,np.arange(1,len(catalog)+1),len(catalog)], line_shape="hv", line_color=ACCENT, showlegend=False),row=1,col=1)
fig.add_trace(go.Scatter(x=catalog.x_km,y=catalog.y_km,mode="markers",marker=dict(size=5+5*(catalog.m-3),color=catalog.day,colorscale="Blues"),showlegend=False),row=1,col=2)
for col,width in [(1,2),(2,15)]:
    edges=np.arange(0,120+width,width)
    counts,_=np.histogram(catalog.day,bins=edges)
    fig.add_trace(go.Bar(x=edges[:-1]+width/2,y=counts,width=width*.9,marker_color=ACCENT,showlegend=False),row=2,col=col)
    fig.update_xaxes(title_text="時間（天）",row=2,col=col)
    fig.update_yaxes(title_text="事件數",row=2,col=col)
fig.update_xaxes(title_text="時間（天）",row=1,col=1)
fig.update_yaxes(title_text="累積事件數",row=1,col=1)
fig.update_xaxes(title_text="東西（km）",row=1,col=2)
fig.update_yaxes(title_text="南北（km）",row=1,col=2)
apply_layout(fig,title=f"同一份合成目錄：{len(catalog)} 個事件，觀察 120 天",height=670)
fig

# %% [markdown]
# 改成較寬的時間箱後，短暫的活躍期被平均進較長的區間。這不代表寬箱錯了：
# 若問題是每月需要多少處理量，月計數很有用；若問題是一次事件後幾小時的變化，
# 月計數就太粗。統計方法必須配合問題的尺度，而不是選一個最能讓圖出現峰值的
# 箱寬。把箱子的起點平移，也可能改變最高峰落在哪一箱。
#
# 空間分格同樣如此。地震發生在行政區邊界兩側時，按行政區計數會把鄰近事件
# 分開；一個大格內也可能混有不同斷層與觀測條件。Reinhart 的綜述從這個
# 聚合問題說明點過程的動機：保留事件的位置與時間，可以避免一開始就把所有
# 訊息壓進任意格網。但保留原始位置也不會自動消除定位誤差。
#
# ## 2.3 次數需要搭配觀察的機會
#
# 一個區域一年有很多事件，另一區只有少數事件，還不能直接說前者比較活躍。
# 兩區面積是否一樣？觀測時間是否一樣長？能記到的最小規模是否相同？
# 次數的比較需要搭配「有多少機會被觀察」，也就是觀測時間、面積及偵測範圍。
#
# 最簡單的平均率是事件數除以觀測時間，單位可以是次／天。它把不同長度的
# 觀測期放到相同尺度，卻還沒有處理活動隨時間變化的情況。如果所有事件都集中
# 在第一週，四個月的平均率仍會把這段活躍攤平。下一章會比較平均率完全相同、
# 排列卻很不同的三份目錄，讓這個限制更具體。
#
# 另一個容易忽略的機會差異來自偵測。假設地下的活動沒有變，但一段期間較小的
# 事件比較難被記錄，曲線也會變平。下面直接從同一份合成資料移除部分較小事件，
# 並把完整與缺漏的版本放在一起。這是人工設定的漏測示範，不是某套儀器的偵測模型。
#
# %% tags=["remove-input"]
rng = np.random.default_rng(20260912)
eligible = (catalog.day.between(40,80)) & (catalog.m<3.6)
keep = ~eligible | (rng.random(len(catalog))<.3)
observed=catalog.loc[keep]
fig=go.Figure()
for frame,name,color in [(catalog,"原始合成事件",PALETTE[0]),(observed,"人工漏測後",PALETTE[1])]:
    fig.add_trace(go.Scatter(x=np.r_[0,frame.day,120],y=np.r_[0,np.arange(1,len(frame)+1),len(frame)],mode="lines",line_shape="hv",name=name,line_color=color))
fig.add_vrect(x0=40,x1=80,fillcolor="gray",opacity=.12,line_width=0)
apply_layout(fig,title=f"觀測過程也會改變曲線：人工移除 {len(catalog)-len(observed)} 筆小事件",xaxis_title="時間（天）",yaxis_title="累積事件數",height=380)
fig

# %% [markdown]
# 陰影區內兩條曲線開始分開；離開陰影區後，累積的差距仍然留下。
# 只看到下方曲線的人，可能把斜率下降解讀成地下變安靜了。知道偵測流程的人，
# 則會先檢查資料是否完整。完整度會直接影響我們如何解釋活動變化，
# 需要和發生率一起考慮。{doc}`11_catalog_completeness_b`會進一步
# 說明如何從規模分布判斷可用的門檻。
#
# ## 2.4 從看見差異到提出可檢驗的問題
#
# 現在可以把「最近地震好像變多」拆成比較清楚的問題。先固定區域、規模門檻
# 與時間窗，再選一個比較基準，例如同一區過去的平均活動，或能描述餘震衰減的
# 模型。觀察到的差異有多大，是資料提供的結果；在基準下這種差異多常出現，
# 才是統計推論要回答的部分。
#
# 這個順序也讓我們分清三件事。描述是把已經發生的樣貌整理出來；解釋是提出
# 為什麼出現這個樣貌的假設；預報則必須只用當時已知的資料，對未來提出可評估
# 的機率。三者會互相幫助，但一張漂亮的回顧圖不能直接充當預報成績。
# 例如先找到主震、再往前挑一段最活躍的期間，很適合產生研究假說，卻需要
# 在沒有參與挑選的資料上檢查，才能知道它是否提供了額外的預報資訊。
#
# 閱讀一張地震統計圖時，不妨先沿著這個順序理解：它選了哪些事件，保留了
# 什麼資訊，和什麼基準比較，最後想回答哪一個問題。這些線索比先記住模型
# 名稱更有用，因為之後每種模型都必須面對同樣的資料選擇。
#
# 合成目錄有些時段擠、有些時段疏。但「擠在一起」究竟超過隨機
# 波動多少？只要是隨機資料，不就也會偶爾連續發生嗎？
# {doc}`foundation_randomness`接著建立這個比較基準。
#
# ## 參考資料與延伸閱讀
#
# - Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I, 第二版](https://doi.org/10.1007/b97277)。第2、3章從事件計數與平穩性建立描述語言，第5章介紹矩與事件對；出版社全文需訂閱或館藏權限。
# - Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)；[免費作者預印本](https://arxiv.org/abs/1708.02647)。先讀第1節對分箱及空間聚合的討論，再看為何保留事件時間與位置有助建模。
# - Naylor, M. 等（2023），[Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759)，開放取用。第3節以合成資料說明歷史缺漏與短期不完整如何改變估計；本章漏測圖只示範這個問題，未重製其實驗。
