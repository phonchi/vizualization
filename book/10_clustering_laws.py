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
# # 10. 地震會互相引發：Omori、Båth 與產能
#
# 2009 年 L'Aquila 地震之後，附近的地震密集出現。把時間拉長，活動通常逐漸減少。第 9 章的 PPE 只看過去地震的鄰近程度，本章開始追問：距離一顆大地震多久，會怎樣改變後續活動？
#
# **主震**（mainshock）是序列中作為參考的主要大事件。**餘震**（aftershock）是其後鄰近地區的相關活動。**群震**（swarm）則是一段密集活動，沒有明顯支配全局的主震。這些名稱整理觀察；單靠先後與距離，仍不能確定哪顆地震觸發哪一顆。
#
# L'Aquila 2009 與 Umbria–Marche 1997 都在學習期內。下圖從 HORUS 取出兩段序列。參考事件定為指定日期附近、指定地區最大的事件。篩選半徑為 50 km、規模至少 3.0，時間為主震後第 1 至第 60 天。這是可重做的教學篩選，沒有替每顆事件認定親代。
#
# ## 10.1 先把計數換成率
#
# 時間箱的寬度不同，計數不能直接比較。每箱事件數除以箱寬，得到每日事件率。圖上越晚的箱越寬，能收集較稀疏的事件。沒有事件的箱仍保留於資料，只因對數軸不能畫零而不顯示。
#
# 兩條序列即使都逐漸下降，局部起伏也會不同。後續較大事件可能帶來新一波活動。把整段都歸給第一顆主震，會將再觸發混進同一條衰減曲線。

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
from gdms_toolkit.viz import ACCENT, apply_layout
cat = italy.experiment_catalog()

def sequence(date, lon, lat):
    x, y = italy.lonlat_to_km(lon, lat)
    near = np.hypot(cat.x_km-x, cat.y_km-y) <= 50
    day = pd.Timestamp(date)
    candidates = cat.loc[near & (cat.time >= day) & (cat.time < day+pd.Timedelta(days=1))]
    main = candidates.loc[candidates.mb.idxmax()]
    tau = (cat.time-main.time).dt.total_seconds().to_numpy()/86400
    keep = (np.hypot(cat.x_km-main.x_km, cat.y_km-main.y_km) <= 50) & (cat.mb >= 3.0) & (tau >= 1) & (tau < 60)
    return main, np.sort(tau[keep])

sequences = {name: sequence(*args) for name,args in {
    "L'Aquila 2009": ('2009-04-06',13.38,42.35),
    'Umbria–Marche 1997': ('1997-09-26',12.89,43.03)}.items()}
bins = np.geomspace(1,60,15)
fig = go.Figure()
for (name,(main,ts)),color in zip(sequences.items(),[ACCENT,'#64748b']):
    counts,_ = np.histogram(ts,bins)
    yy = np.where(counts>0,counts/np.diff(bins),np.nan)
    fig.add_scatter(x=np.sqrt(bins[:-1]*bins[1:]),y=yy,mode='lines+markers',name=f'{name}，{len(ts)} 顆',line_color=color)
apply_layout(fig,title='HORUS：同一篩選規則下的兩段活動',xaxis_title='距參考主震（天）',yaxis_title='事件數／天',xaxis_type='log',yaxis_type='log')
fig

# %% [markdown]
# ## 10.2 用 Omori–Utsu 律描述下降
#
# **Omori–Utsu 律**是一個餘震率的經驗模型。以 $\tau$ 表示主震後天數，寫成
#
# $$r(\tau)=A(\tau+c)^{-p},\qquad A>0,\ c>0.$$
#
# $A$ 調整活動量，$c$ 是早期轉折的時間尺度，$p$ 決定較晚的衰減速度。這裡用 $A$ 避免與後面的正規化產能 $K$ 混淆。時間改用小時，$A$ 與 $c$ 也須一起轉換。
#
# **冪律**（power law）指數量依某個變數的次方改變。當 $\tau$ 遠大於 $c$，時間加倍，率約乘 $2^{-p}$。若 $p=1$，第 20 天的率約為第 10 天的一半。這是平均曲線的關係，單次序列的計數仍會波動。
#
# 若只研究第 $S$ 至第 $T$ 天，將 $(\tau+c)^{-p}$ 除以該窗積分，就得到總面積為一的時間密度。給定窗內事件數，使用第 5 章的密度概似可估 $c,p$：
#
# $$\ell(c,p)=-p\sum_i\ln(\tau_i+c)-N\ln I(c,p),\quad
# I(c,p)=\int_S^T(\tau+c)^{-p}\,d\tau.$$
#
# 圖中曲線直接使用事件時間擬合。估計完成後取 $\widehat A=N/I$，使窗內積分等於觀測數。這個配合是估計方法帶來的，不是模型通過驗證的證據。分箱只負責呈現，不參與這次估計。

# %% tags=["remove-input"]
def integ(a,b,c,p):
    if abs(p-1) < 1e-7:
        return np.log((b+c)/(a+c))
    return ((b+c)**(1-p)-(a+c)**(1-p))/(1-p)

main, ts = sequences["L'Aquila 2009"]
def objective(theta):
    c,p = np.exp(theta[0]),theta[1]
    return p*np.log(ts+c).sum()+len(ts)*np.log(integ(1,60,c,p))
fits = [minimize(objective,[np.log(c),p],method='L-BFGS-B',bounds=[(np.log(1e-4),np.log(10)),(.3,2.5)]) for c,p in [(.01,.8),(.1,1.2),(1.,1.5)]]
valid = [r for r in fits if r.success and np.isfinite(r.fun)]
if not valid:
    raise RuntimeError('Omori 最佳化未收斂')
fit = min(valid,key=lambda r:r.fun)
c_hat,p_hat = np.exp(fit.x[0]),fit.x[1]
a_hat = len(ts)/integ(1,60,c_hat,p_hat)
tgrid=np.geomspace(1,60,300)
counts,_=np.histogram(ts,bins)
fig=go.Figure()
fig.add_scatter(x=np.sqrt(bins[:-1]*bins[1:]),y=np.where(counts>0,counts/np.diff(bins),np.nan),mode='markers',name='HORUS 分箱率',marker_color='#64748b')
fig.add_scatter(x=tgrid,y=a_hat/(tgrid+c_hat)**p_hat,name='有限窗 MLE',line_color=ACCENT)
apply_layout(fig,title=f"L'Aquila：c={c_hat:.3g} 天，p={p_hat:.3f}",xaxis_title='距主震（天）',yaxis_title='事件數／天',xaxis_type='log',yaxis_type='log')
fig

# %% [markdown]
# $c$ 的搜尋範圍為 $10^{-4}$ 至 10 天，$p$ 為 0.3 至 2.5。若最佳解靠近邊界，資料可能沒有充分約束該參數。尤其只保留一天以後的事件，很難看見最早期轉折，不能把估得的 $c$ 當成精確物理時間。
#
# 主震剛發生時，重疊波形會使小地震漏記。這會讓觀測率看似比較平坦，也會改變 $c,p$。提高規模門檻或延後起算，可減少部分影響，同時也會失去資料。第 8 章的完整度問題，因此直接進入參數解讀。
#
# ## 10.3 從衰減率拆出時間密度
#
# 模型還需要區分「總共幾顆」與「何時出現」。若時間範圍延伸至無限，正規化的 Omori 密度為
#
# $$g(\tau)=\frac{p-1}{c}(1+\tau/c)^{-p},\qquad \tau>0,\ p>1.$$ (eq:omori-density)
#
# 密度總面積為一，每一小段面積代表直接後代落在該時段的比例。乘上產能，才得到事件率。有限窗擬合允許 $p\le1$；無限時間密度卻需要 $p>1$。兩者的積分範圍不同，因此沒有矛盾。完整積分見附錄 B。
#
# **產能**（productivity）是某顆事件平均產生的直接後代數。以輸入門檻 $m_0$ 為基準，常寫成
#
# $$\kappa(m)=K e^{\alpha(m-m_0)}.$$
#
# $K$ 是門檻規模事件的產能；$\alpha$ 表示產能增加速度。規模增加一級，產能乘上 $e^{\alpha}$。這是模型中的直接後代，與半徑內所有後續事件的計數有別。後者還包含背景及更後面的世代。
#
# ## 10.4 最大餘震與序列整理
#
# **Båth 律**（Båth's law）描述許多序列中，主震與最大餘震的平均規模差約為 1.2。這是跨序列的經驗摘要，不能解讀為每次都差 1.2，或下一顆地震必小於主震。若後來出現更大的地震，原先的主震名稱也可能改變。
#
# 最大值會隨抽樣數改變。同樣的規模分布，抽到更多事件時，最大值通常更大。因此 Båth 律、GR 律與產能應一起閱讀。選定主震門檻、餘震範圍與截止時間，也會改變這個差值；短目錄可能還沒看見最大餘震。
#
# **除叢**（declustering）是將目錄中的群集活動分離或降低權重。Gardner–Knopoff 方法依主震規模給出時空窗，再標記窗內事件。它容易操作，但窗界線是一套規則。留下的事件不會自動成為已證實獨立的背景。
#
# 本書四模型的主線使用包含餘震的目錄。此處不先刪除一批事件。第 13 章會介紹用模型給親代機率的方法，讓同一事件能以不同程度歸入背景與觸發活動。兩種方法回答的問題不同，結果也不必一致。
#
# ## 10.5 讀曲線時，檢查三個尺度
#
# 第 10 天有十顆事件，單看這個數字很難判斷活動強弱。必須同時知道規模門檻、空間範圍與時間箱寬。把門檻從 3.0 降低，通常會收進更多小地震；把半徑擴大，則可能混入別的序列。
#
# 因此比較兩條 Omori 曲線前，先讓這三項設定一致。即使設定相同，目錄完整度也可能不同。1997 年與 2009 年的儀器觀測條件不能由圖形相似就視為相同；這張圖適合比較形狀，不足以單獨歸因地區的物理差異。
#
# 一條平滑曲線還隱藏了事件的順序。若第 20 天再出現大事件，第 21 天的活動可能相應升高。單主震公式只會持續下降，便可能把新一波活動解讀成較小的衰減指數。下一章讓每顆事件都提供貢獻，正是為了表達這種情況。
#
# $\widehat A=N/I$ 也提示另一個限制：總數被估計步驟配合掉了。若將同一份序列再拿來問總數是否正確，這個檢查提供的資訊很少。更有意義的問題是另一段未參與估計的時間，或同窗內的詳細事件時間能否被模型描述。
#
# 最後，產能的平均值允許非常不同的實現。平均直接後代數是一，不表示每顆來源剛好生出一顆。可能多數沒有後代，少數形成較大的家族。這種家族結構讓計數比獨立事件更容易成群，也為第 16 章的過度離散檢驗提供背景。
#
# 同一個家族也可能同時提供多顆目標地震。將這些事件當成許多次互不相關的成功，會誇大證據。往後的模型比較會保留全部事件，同時明確檢查序列集中造成的相依性。
#
# ## 本章填入的規格欄位
#
# 本章補上時間衰減與產能的語言，沒有更改預報實驗門檻。序列示範使用規模 3.0、50 km、1–60 天；這組診斷範圍不替代全國預報卡。Omori 曲線目前只描述指定序列；下一章才讓每一顆事件都能更新率。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([('本章新增','Omori 衰減與產能的模型假設'),('主線輸入／目標','m₀ = 2.45；mT = 5.0'),('序列示範','學習期；Mw ≥ 3.0；50 km；主震後 1–60 天'),('正式檢驗','待填（第 16–18 章）')])))

# %% [markdown]
# 例如把相同的事件窗切成較寬的時間箱，長條高度會改變。
# 概似估計仍使用每顆事件的時間，因此不隨畫圖的箱寬改變。
# 讀擬合圖時，先分清顯示用分箱與估計用樣本，就能理解兩者各自的角色。
#
# 下一章：{doc}`11_conditional_intensity`。
#
# ## 參考資料與延伸閱讀
#
# - Utsu、Ogata 與 Matsu'ura（1995），[Omori 百年回顧](https://doi.org/10.4294/jpe1952.43.1)。免費全文；先讀開頭的率曲線，再看早期漏測如何影響參數。
# - Zhuang 等（2012），[Seismicity declustering](https://www.corssa.org/en/articles/theme_iii/)。CORSSA 免費教學；比較時空窗與模型權重的假設。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。對照包含餘震與主震資料集的差別。
