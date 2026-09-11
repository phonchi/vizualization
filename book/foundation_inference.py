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
# # 7. 從描述到估計：資料能讓我們知道多少
#
# 前面已經介紹目錄的完整度、規模分布與餘震衰減。它們讓我們知道模型應該
# 留意哪些現象，卻還沒有解決一個共同問題：同樣是一條逐漸下降的曲線，
# 不同參數可能看起來很像，我們怎麼選？選出的參數又能相信到什麼程度？
#
# 在組成 ETAS 之前，先用最簡單的事件率估計走一遍。這樣後面遇到許多參數、
# 看不到的觸發關係，以及貝氏推論時，就能辨認它們分別在處理哪一種困難。
#
# ## 7.1 先把模型與觀測分開
#
# 假設一段固定觀察窗內事件依齊次 Poisson 過程發生，率為 $\mu$。
# 我們觀察到的事件數 $N$ 是資料，$\mu$ 是待估的模型參數。即使 $\mu$
# 完全相同，重新觀察另一段同樣長的期間，也可能得到不同的 $N$。
# 這種差異來自模型本身的隨機性，不是估計程式失敗。
#
# 反過來，拿到一個 $N$，也不能唯一倒推出真正的 $\mu$。較高的率可能碰巧
# 生出較少事件，較低的率也可能碰巧生出較多事件。估計是在這些可能性之間
# 衡量資料支持程度，並不是從資料中讀出一個沒有誤差的答案。
#
# 下面固定觀察時間與真實率，重複模擬許多份目錄。左圖每一個可能的事件數都
# 是相同模型能產生的結果；右圖把事件數除以時間，得到每次觀察的估計率。
#
# %% tags=["remove-input"]
import numpy as np
from scipy.stats import gamma,nbinom
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from gdms_toolkit.viz import setup_plotly,apply_layout,ACCENT,PALETTE
setup_plotly()
rng=np.random.default_rng(20260913)
T_INFER=120.0
MU_TRUE=80/120
counts=rng.poisson(MU_TRUE*T_INFER,3000)
fig=make_subplots(rows=1,cols=2,subplot_titles=("同一模型的重複觀察", "每次用 N/T 估計發生率"))
fig.add_trace(go.Histogram(x=counts,marker_color=ACCENT,showlegend=False),row=1,col=1)
fig.add_trace(go.Histogram(x=counts/T_INFER,marker_color=PALETTE[1],showlegend=False),row=1,col=2)
fig.add_vline(x=MU_TRUE,line_dash="dash",row=1,col=2)
fig.update_xaxes(title_text="觀測事件數",row=1,col=1)
fig.update_xaxes(title_text="估計率（次／天）",row=1,col=2)
fig.update_yaxes(title_text="模擬次數")
apply_layout(fig,title=f"真實率固定為 {MU_TRUE:.3f} 次／天，重複 {len(counts)} 次觀察",height=390)
fig

# %% [markdown]
# 這個示例沿用前面相同的時間長度與平均事件數，但刻意換成已知的 Poisson
# 模型，讓估計的來源能單獨看清楚。前面的群集目錄並不滿足這個假設，不能
# 直接套用這裡的誤差範圍，把事件都當成互不相關。
#
# ## 7.2 概似是在比較哪些參數較能解釋資料
#
# 點過程章已經看過概似的兩部分：模型要在事件發生處給足夠的率，也要對沒有
# 事件的整段觀察負責。在固定率的例子中，對數概似簡化成
# $\ell(\mu)=N\log\mu-\mu T$，略去了不影響參數比較的項。
# 把這條曲線推到最高點，得到 $\hat\mu=N/T$。
#
# 因此，熟悉的「次數除以時間」正是這個模型的最大概似估計。
# 但概似不是參數的機率密度：橫軸不同的 $\mu$ 是我們拿來比較的設定，
# 不是已經賦予機率的隨機變數。曲線很尖，表示偏離最高點很快會與資料不合；
# 曲線很平，表示許多設定仍然難以區分。
#
# 在複雜模型裡，平坦的方向可能來自參數互相補償。例如增加背景活動、減少
# 觸發活動，有時仍能描出類似的總率。這叫辨識困難：最佳化程式可以回傳
# 很多小數位，卻不代表資料真的有那麼多資訊。後面 ETAS 估計會把這個問題
# 畫成具體的參數關係。
#
# ## 7.3 看不到觸發關係，還能估計嗎
#
# 假設一個事件附近有兩個更早的事件，我們不知道它是背景活動，還是由其中
# 一個先前事件所觸發。如果先硬分一個家族，再估計模型，就把不確定的分類
# 當成了確定資料。另一條路是把家族關係視為**潛在變數**：它參與模型，
# 但沒有被直接觀測。
#
# EM 的想法是交替處理這兩層。先用目前模型，為每種可能來源分配權重；
# 再用這些權重更新參數；新參數又會改變下一輪的權重。直覺上，某事件與
# 一個較早的大事件很接近時，該來源可能得到較高權重，但背景仍可能佔一部分。
#
# 這個過程沒有替我們發現唯一、可直接當成物理事實的家譜。它是在模型假設下
# 處理缺少的資訊。EM 也不保證找到全域最佳解；初始化、停止條件與可辨識性
# 仍然重要。Veen 與 Schoenberg 將這個想法用在地震分支模型，後面的
# {doc}`14_etas_estimation`會把權重寫成背景與各觸發源的相對貢獻。
#
# ## 7.4 貝氏推論保留一整組可能的參數
#
# 最大概似先找資料最支持的設定。貝氏推論則先明確指定參數的先驗分布，
# 再以概似更新，得到後驗分布：
#
# $$\text{posterior}\ \propto\ \text{likelihood}\times\text{prior}.$$
#
# 先驗可以表達已有知識，也可以限制不合理的參數區域；它不是免費增加資料。
# 資料很少或參數難辨識時，先驗的影響可能很大，因此需要說明選擇並檢查敏感度。
# 後驗區間則是在這套模型與先驗條件下，描述參數仍有哪些可能值。
#
# 用前面的固定率例子，選擇 Gamma 先驗後可以直接算出後驗。下圖用同樣的
# 模擬真實率，比較短期與較長觀察所得到的曲線。虛線是生成資料時的率，
# 不是估計時可以偷看的答案；增加資料通常會縮小不確定性，但單次估計不保證
# 一定更接近虛線。這裡先驗的影響較弱，後驗外觀會接近標準化的概似曲線；
# 差別不在於曲線一定長得不同，而在於後驗已把先驗納入，具有參數機率的意義。
# 單純把概似縮放成面積為一，不會在未指定先驗的情況下自動賦予它這種意義。
#
# %% tags=["remove-input"]
alpha_prior,beta_prior=1.0,1.0
mu_grid=np.linspace(.01,1.6,500)
fig=go.Figure()
observations=[]
for i,horizon in enumerate([15.,120.]):
    n=int(rng.poisson(MU_TRUE*horizon))
    observations.append((horizon,n))
    density=gamma.pdf(mu_grid,a=alpha_prior+n,scale=1/(beta_prior+horizon))
    fig.add_trace(go.Scatter(x=mu_grid,y=density,name=f"{horizon:g} 天，觀察 {n} 次",line_color=PALETTE[i]))
fig.add_vline(x=MU_TRUE,line_dash="dash")
apply_layout(fig,title="資料更新後的率參數分布（Gamma 先驗的教學示例）",xaxis_title="率 μ（次／天）",yaxis_title="後驗密度",height=390)
fig

# %% [markdown]
# ```{admonition} 這個例子的後驗怎麼算？
# :class: dropdown
#
# 以shape–rate記法，若 $\mu\sim\operatorname{Gamma}(a,b)$，其密度含
# $\mu^{a-1}e^{-b\mu}$；概似含 $\mu^N e^{-T\mu}$。兩者相乘後，
# $\mu\mid N\sim\operatorname{Gamma}(a+N,b+T)$。這是特定模型的共軛結果，
# 並不表示ETAS的多參數後驗也能直接寫成Gamma分布。
# ```
#
# 地震模型裡，後驗常需要數值近似。Naylor 等人的研究用 INLA 近似時間 ETAS
# 的後驗，讓參數間的相關與目錄缺漏的影響更容易被檢查。這裡要學的重點是
# 為什麼保留參數分布有用；計算方法與具體適用範圍放在 ETAS 章和附錄。
#
# ## 7.5 參數知道得更清楚，未來仍會波動
#
# 參數區間回答「率可能是多少」，預測分布回答「下一段期間可能發生幾次」。
# 就算率完全已知，未來事件數仍會波動；率還不確定時，則要把不同率下的
# 未來結果一起考慮。這是後驗預測的直覺，而不是把後驗平均代進模型就結束。
#
# 在Gamma–Poisson例子中，令 $h$ 為未來期間長度，總變異數公式給出：
#
# $$\operatorname{Var}(N_{\rm future}\mid\text{data})
# =hE(\mu\mid\text{data})+h^2\operatorname{Var}(\mu\mid\text{data}).$$
#
# 第一項是即使知道率仍有的事件波動，第二項來自率未知。這個分解只針對
# 目前的Poisson模型；自激發目錄還有事件間相依，不能不加修改地沿用。
# 模型結構是否合理、測站是否漏測，也沒有因為畫出一條後驗曲線就被解決。
#
# ## 7.6 回到資料，檢查模型留下什麼
#
# 有了參數與區間，下一步仍是回到資料。模型是否把密集期與安靜期都描述好？
# 殘差是否還有成群結構？換一段沒有參與估計的期間，它是否比簡單基準更好？
# 這三個問題分別指向擬合、診斷與預報評估，不應以同一個分數概括。
#
# 重抽樣與模擬能幫助檢查估計的變動，但抽樣方式必須配合資料。把互相觸發的
# 每一筆地震當成獨立資料做普通bootstrap，可能低估不確定性；以模型模擬
# 整份目錄，又會依賴模型是否合適。選擇方法時，需要交代它保留了哪些相依
# 結構，而不是把「有算誤差」當作充分保證。
#
# 接下來閱讀模型時，要分清事件怎麼產生、資料如何限制參數，以及
# 參數未知如何影響未來預報。接下來在{doc}`13_etas_structure`，將背景、
# 規模與衰減組合成一個會持續更新的事件模型；再用這一章的思路理解它的估計。
#
# ## 參考資料與延伸閱讀
#
# - Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I, 第二版](https://doi.org/10.1007/b97277)。第7章連結條件強度、概似、模擬與診斷；全文需訂閱或館藏權限。
# - Veen, A. 與 Schoenberg, F. P.（2008），[Estimation of Space–Time Branching Process Models in Seismology Using an EM–Type Algorithm](https://doi.org/10.1198/016214508000000148)。第4節說明如何將未知家族關係當成不完整資料；出版社全文可能需訂閱。
# - Zhuang, J., Ogata, Y. 與 Vere-Jones, D.（2002），[Stochastic Declustering of Space-Time Earthquake Occurrences](https://doi.org/10.1198/016214502760046925)。第3–5節連結背景率、觸發機率與隨機除叢，幫助理解為何分類也有不確定性；出版社全文可能需訂閱。
# - Naylor, M. 等（2023），[Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759)，開放取用。第3節比較隨機實現、資料長度、歷史缺漏與短期不完整對後驗的影響。
# - Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)；[免費作者預印本](https://arxiv.org/abs/1708.02647)。第3節將估計、不確定性與診斷放在同一個方法架構中。
