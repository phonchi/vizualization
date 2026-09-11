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
# # 3. 隨機性與叢集：平均之外還有什麼
#
# {doc}`foundation_catalog`讓我們看見，累積曲線的斜率會隨時間改變。
# 但即使發生機制始終相同，隨機事件也不會等間隔出現。要判斷某一段密集事件
# 是否需要額外解釋，必須先知道：在一個簡單、明確的隨機模型下，資料本來就會
# 有多不整齊？
#
# 先沿用前章的觀察範圍，把同樣數目的事件排成三種序列。一份接近等間隔，
# 一份均勻隨機撒點，另一份是前章的群集目錄。它們的平均率相同，所以只報
# 「每天平均幾次」無法分辨。這次我們把注意力放在事件之間的距離。
#
# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from gdms_toolkit.teaching import learning_catalog
from gdms_toolkit.viz import setup_plotly,apply_layout,ACCENT,PALETTE,QUAKE_COLOR
setup_plotly()
catalogs={k:learning_catalog(k) for k in ["regular","poisson","cluster"]}
names={"regular":"近規律", "poisson":"均勻隨機（固定總數）", "cluster":"分群（固定總數）"}
fig=make_subplots(rows=3,cols=1,shared_xaxes=True,subplot_titles=list(names.values()))
for row,(k,frame) in enumerate(catalogs.items(),1):
    fig.add_trace(go.Scatter(x=frame.day,y=np.ones(len(frame)),mode="markers",marker=dict(symbol="line-ns",size=15,color=PALETTE[row-1],line=dict(width=2,color=PALETTE[row-1])),showlegend=False),row=row,col=1)
    fig.update_yaxes(visible=False,row=row,col=1)
fig.update_xaxes(title_text="時間（天）",row=3,col=1)
apply_layout(fig,title=f"每份都是 {len(catalogs['cluster'])} 個事件／120 天，排列方式卻不同",height=410)
fig

# %% [markdown]
# 中間的隨機撒點也有短暫的密集和空白。因此，看到幾個點靠近，
# 還不足以排除隨機性。下方分群序列則在更多地方同時出現短間隔與較長空白，
# 暗示只用平均率可能漏掉結構。這裡固定了總數，方便比較排列；中間序列對應
# 齊次 Poisson 過程「已知這段時間共有這些事件」的條件分布，而不是一次完整
# 的無條件 Poisson 模擬。
#
# ## 3.1 Poisson 模型為什麼是起點
#
# 假設每個很短的時間段都有相同的小機會發生事件，而且不重疊期間的事件數
# 互相獨立，就得到齊次 Poisson 過程。它只有一個率參數 $\mu$；觀察長度
# 為 $T$ 的時間窗時，事件數的平均與變異數都是 $\mu T$：
#
# $$E[N(T)]=\mu T,\qquad \operatorname{Var}[N(T)]=\mu T.$$
#
# 這個關係使模型可以檢查。假設許多可比較的時間窗有相近的觀測條件，如果
# 事件數的變動遠比平均大，固定率、獨立事件這一組假設就值得重新檢視。
# 但不能立刻跳到「證明餘震互相觸發」：各窗的率不同、漏測程度不同，甚至
# 混入不同地區，也可能讓變異增加。
#
# Poisson 是有用的起點，並不是預設所有地震都互不影響。對某些長時間平均
# 或經過特定篩選的問題，它可能是合理近似；對正在演化的餘震序列，固定率通常
# 太粗略。模型是否合用，取決於它要回答的問題及要保留的尺度。
#
# ## 3.2 從等待時間看見另一面
#
# 把相鄰事件的時間相減，就得到事件間隔。對齊次 Poisson 過程，下一個間隔
# $W$ 服從指數分布，$P(W>w)=e^{-\mu w}$。這表示短間隔常見，但長空白也
# 有可能發生；「平均每兩天一次」並不表示第二天一定輪到下一次。
#
# 指數分布還有一個特別性質：如果你已經等了一段時間，剩下的等待時間分布
# 仍然相同，這叫無記憶性。它是這個模型的假設，不是「隨機」兩字的必然意思。
# 若一次大破裂後需要重新累積條件，距上次事件多久可能影響下一次的率；
# 若事件後容易接連出現新事件，剛發生過反而可能縮短下一次等待。
#
# 下面只比較三份目錄中完整觀測到的相鄰間隔；觀察開始前與結束後的兩段沒有
# 完整事件對，所以不把它們混成普通間隔。虛線是同平均率、未固定總數的Poisson
# 指數參考，延續{doc}`09_forecasting_intro`用過的比較基準，這次改畫成累積比例。
# 固定總數會使有限窗的精確間隔分布略有不同；此處用虛線提供近似參考，
# 精確關係見{doc}`appendix_a_point_process`。
#
# %% tags=["remove-input"]
fig=make_subplots(rows=1,cols=2,subplot_titles=("相鄰間隔的累積比例", "箱寬改變，離散程度也改變"))
for idx,(k,frame) in enumerate(catalogs.items()):
    gaps=np.sort(np.diff(frame.day))
    fig.add_trace(go.Scatter(x=gaps,y=np.arange(1,len(gaps)+1)/len(gaps),name=names[k],line_color=PALETTE[idx]),row=1,col=1)
    widths=np.array([1,2,5,10,20])
    fanos=[]
    for width in widths:
        counts=np.histogram(frame.day,bins=np.arange(0,120+width,width))[0]
        fanos.append(counts.var(ddof=1)/counts.mean())
    fig.add_trace(go.Scatter(x=widths,y=fanos,name=names[k],showlegend=False,line_color=PALETTE[idx],mode="lines+markers"),row=1,col=2)
gap_grid=np.linspace(0,max(np.diff(frame.day).max() for frame in catalogs.values()),300)
fig.add_trace(go.Scatter(x=gap_grid,y=1-np.exp(-(80/120)*gap_grid),name="Poisson 指數參考（未固定總數）",line=dict(color="#666666",dash="dash")),row=1,col=1)
fig.add_hline(y=1,line_dash="dash",line_color="#666666",row=1,col=2)
fig.update_xaxes(title_text="間隔（天）",row=1,col=1)
fig.update_yaxes(title_text="累積比例",row=1,col=1)
fig.update_xaxes(title_text="箱寬（天）",row=1,col=2)
fig.update_yaxes(title_text="計數變異數／平均",row=1,col=2)
apply_layout(fig,title="以等待時間與分箱計數描述同一組序列",height=430)
fig

# %% [markdown]
# 左圖在很短間隔就迅速上升的曲線，代表有較多緊接著出現的事件；接近規律的
# 序列則把多數間隔集中在狹窄範圍。右圖的變異數／平均稱為 Fano 因子。
# 同一份資料在不同箱寬得到不同結果，說明「聚集有多強」必須連同尺度一起說。
# 對這裡的均勻隨機序列，使用樣本變異數計算時，各箱寬的Fano期望值都是1，
# 即使事件總數固定也一樣；右圖的水平虛線標出這個基準。單一實現仍會上下波動，
# 尤其寬箱只剩少量時間箱時更明顯，不能把偏離1直接當成正式檢定或物理觸發的證明。
#
# 我們也可以問任一事件附近，是否比基準更常找到其他事件。這類問題屬於
# **二階統計**：平均率描述一個位置附近的活動量，事件對則描述兩個位置或
# 時刻共同出現的方式。若某個區域原本就有較高背景率，附近事件自然較多，
# 所以事件對的比較要先處理這種不均勻性。{doc}`appendix_a_point_process`提供矩與成對計數的形式，
# 這裡先保留「平均」與「相依結構」是不同資訊這個觀念。
#
# ## 3.3 同樣成群，可能來自不同機制
#
# 假設某一段時間事件一起變多。一個解釋是外部條件改變，使所有事件都比較
# 容易發生；另一個解釋是第一批事件提高了後續事件的發生率。兩者都能畫出
# 成群的目錄，但對下一步的預報，以及參數的物理解讀，意義並不相同。
#
# 在 **Cox 過程**中，發生率本身是隨機的；給定那條率之後，事件依 Poisson
# 機制產生。它可以表達未完全觀測到的環境變化。非齊次 Poisson 則把變動的率
# 視為給定的確定函式，這兩者不能只因曲線都會起伏就混為一談。
#
# **群集過程**先產生一些中心，再在中心附近產生事件。開頭的分群示例就是
# 為了展示這個構造，但中心數和每群大小被固定，不能直接當成擬合好的地震
# 模型。**Hawkes 自激發過程**進一步讓每一個事件都可能提高後續活動，
# 新出現的事件也能繼續影響下一代。這個遞迴的結構將在 ETAS 章成為主角。
#
# 與它們不同，**更新過程**在每次事件後重新抽下一個間隔，記住的是距離
# 上次事件的時間；**自我修正模型**則讓事件發生後的率下降，表達某種累積與
# 釋放的想法。這些模型對歷史資訊如何影響未來，各自作了不同假設。同一個地區在短期餘震與長期斷層複發的問題中，也可能需要不同描述。
#
# ```{admonition} 多等一會兒，是否就更接近下一次？
# :class: dropdown
#
# 指數間隔的無記憶性可直接由存活機率看出：
# $P(W>s+t\mid W>s)=e^{-\mu(s+t)}/e^{-\mu s}=e^{-\mu t}$。
# 一般間隔分布未必滿足這個等式。另外，在隨機時刻開始等待，較容易落在
# 較長的間隔內；「從一次事件之後開始計時」與「任選今天開始等」因此可能
# 得到不同的等待分布。這是觀察方式造成的差別，並不表示模型自相矛盾。
# ```
#
# ## 3.4 從目錄的樣貌走向可更新的預報
#
# 比較模型時，要看它們如何利用已經發生的事，不能只看是否畫得出密集的一團點。若新來一個事件，固定率模型不改變預報；自激發模型
# 會提高近期活動；更新模型則重新計算距上次事件的時間。
#
# 這也把「統計觸發」與「物理因果」分開了。自激發模型成功描述資料，表示
# 過去事件對後續活動有預報資訊；但未觀測的共同原因可能仍然存在。要進一步
# 主張應力傳遞或流體變化，還需要物理觀測與其他證據。第二部會回到這個問題。
#
# 下一步需要一個共同的量，讓這些不同假設都能回答「已知現在以前的資料，
# 接下來一小段時間有多活躍」。{doc}`10_point_process`會把這個量寫成
# 條件強度，並說明如何從它連到機率與模型檢查。
#
# ## 參考資料與延伸閱讀
#
# - Gallager, R. G.／MIT OpenCourseWare（2011），[Discrete Stochastic Processes，第2章](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/resources/mit6_262s11_chap02/)，免費講義。從計數、指數等待時間與無記憶性三個角度理解Poisson過程。
# - Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I, 第二版](https://doi.org/10.1007/b97277)。第4章討論更新過程，第5、8章處理矩與二階結構，第6章比較Cox與群集過程；全文需訂閱或館藏權限。
# - Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)；[免費作者預印本](https://arxiv.org/abs/1708.02647)。第1–2節幫助區分外部變動、群集與自激發，並說明相似數學結構如何出現在不同應用中。
