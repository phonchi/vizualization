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
# # 15. EEPAS：每顆地震都是它尺度上的前兆
#
# 一顆規模 3.5 的地震，可以提供哪些未來資訊？第 14 章的 Ψ 尺度關係，連結前兆規模、等待時間與活動面積。本章將這三種尺度轉成密度，再把所有輸入事件的貢獻加起來。
#
# **EEPAS**（Every Earthquake a Precursor According to Scale）意為每顆地震都是其尺度上的前兆。模型把每顆事件當成一份微弱線索，按規模分配未來活動。它不先判定哪一顆是真前兆，也不承諾每顆小震後都出現指定大震。
#
# ## 15.1 先讀一顆事件的三個核
#
# **條件密度**（conditional density）描述給定資訊後，另一個量的分布。這裡給定事件 $i$ 的位置、時間與規模，再問未來事件的規模、時間與位置。第 9 章的核是平滑貢獻；本章將三個密度各當成一個核。
#
# 規模核 $g(m\mid m_i)$ 是常態密度，中心在 $a_M+b_Mm_i$，標準差為 $\sigma_M$。輸入規模改變，整條曲線的中心跟著移動。曲線下某個規模箱的面積，才是落入該箱的比例。
#
# 時間核使用 $\log_{10}\tau$，其中 $\tau=t-t_i$ 以天計。**對數常態分布**（lognormal distribution）指一個正值變數取對數後服從常態。對數時間的中心是 $a_T+b_Tm_i$，標準差是 $\sigma_T$。
#
# 取對數能將不同數量級放在同一尺度。一百天與一千天相差一個十倍，對數軸上相差一格。這個轉換不表示一天和十天有相同機率；從對數密度換回每日密度時，還需要除以 $\tau\ln10$。
#
# 空間核 $h_i(x,y)$ 是以輸入震央為中心的圓對稱二維常態。每一座標方向的標準差為 $\sigma_i=\sigma_A10^{b_Am_i/2}$，單位是公里。較大輸入事件會把貢獻分配到較寬範圍。
#
# 下圖指定一顆規模 3.5 的合成事件。參數來自 Biondini 等（2023）表 3 的等權重版本。三張小圖分別看規模、時間及空間，不代表三份不同預報。空間面板畫穿過震央的一條截線，面積積分仍須在二維平面計算。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import norm
from IPython.display import HTML, display
from gdms_toolkit import italy, italy_models
from gdms_toolkit.viz import ACCENT, apply_layout, plot_forecast_map
P=italy_models.PARAMS['EEPAS']
mi=3.5
m=np.linspace(2.5,6.5,350)
tau=np.geomspace(100,50000,350)
x=np.linspace(-60,60,350)
mt=P['aT']+P['bT']*mi
sig=P['sA']*10**(P['bA']*mi/2)
fig=make_subplots(rows=1,cols=3,subplot_titles=('規模密度','時間密度','空間密度截線'))
fig.add_scatter(x=m,y=norm.pdf(m,loc=P['aM']+P['bM']*mi,scale=P['sM']),row=1,col=1,line_color=ACCENT,showlegend=False)
fig.add_scatter(x=tau,y=norm.pdf(np.log10(tau),loc=mt,scale=P['sT'])/(tau*np.log(10)),row=1,col=2,line_color=ACCENT,showlegend=False)
fig.add_scatter(x=x,y=np.exp(-x*x/(2*sig*sig))/(2*np.pi*sig*sig),row=1,col=3,line_color=ACCENT,showlegend=False)
fig.update_xaxes(title_text='未來規模',row=1,col=1)
fig.update_xaxes(title_text='距輸入事件（天）',type='log',row=1,col=2)
fig.update_xaxes(title_text='東西距離（km）',row=1,col=3)
fig.update_yaxes(title_text='每規模單位',row=1,col=1)
fig.update_yaxes(title_text='每天',row=1,col=2)
fig.update_yaxes(title_text='每 km²',row=1,col=3)
apply_layout(fig,title=f'合成輸入 m={mi:.1f}：三個尺度的密度')
fig

# %% [markdown]
# ## 15.2 相乘之前，先說假設
#
# **聯合密度**（joint density）描述幾個量同時落入指定範圍的可能性。EEPAS 採用可分離的乘積核：給定輸入事件後，時間、規模及位置的貢獻寫成 $f_i g_i h_i$。這是模型的條件獨立假設，不是資料已證明三個量互不相關。
#
# 三核的公式集中寫成
#
# $$\begin{aligned}
# g_i(m)&=\frac{1}{\sigma_M\sqrt{2\pi}}
#  e^{-(m-a_M-b_Mm_i)^2/(2\sigma_M^2)},\\
# f_i(t)&=\frac{1}{\tau\sigma_T\ln10\sqrt{2\pi}}
#  e^{-(\log_{10}\tau-a_T-b_Tm_i)^2/(2\sigma_T^2)},\\
# h_i(x,y)&=\frac{1}{2\pi\sigma_i^2}
#  e^{-[(x-x_i)^2+(y-y_i)^2]/(2\sigma_i^2)}.
# \end{aligned}$$
#
# 時間核在 $\tau\le0$ 時為零。每個核在其完整範圍內積分為一；輸出只取三個月、177 格與目標規模箱，通常只能收到其中一部分。落在區外的質量不能搬回區內，除非另行宣告重新正規化。
#
# 單顆事件的原始貢獻還須乘 $\eta(m_i)$。**規模權重** $\eta$ 調整不同大小輸入事件的總貢獻，使長期規模率與指定 GR 關係相容。它包含背景比例的設定，不等同於每顆事件實際產生的後代數。
#
# **門檻補償** $\Delta(m)$ 表示在理想規模模型下，只收錄 $m_i\ge m_0$ 時保留下來的比例。將規模核除以 $\Delta(m)$，補回模型預期的缺失貢獻。它依賴指定尺度關係，不能讓實際漏記事件重新出現。推導與小分母風險見附錄 D。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import show_diagram
show_diagram("d15_eepas_kernels", caption="一顆輸入事件在規模、時間、空間三個方向提供貢獻。")

# %% [markdown]
# ## 15.3 所有事件加起來，再加背景
#
# 第 9 章的 PPE 描述鄰近過去地震的背景形狀。以 $\lambda_0$ 表示 PPE 率密度，EEPAS 在發報時刻 $t_0$ 的形式是
#
# $$\lambda(t,x,y,m\mid H_{t_0})=
# \mu_E\lambda_0(t,x,y,m)+
# \sum_{t_i\le t_0-d}\eta(m_i)
# \frac{f_i(t)g_i(m)h_i(x,y)}{\Delta(m)},\qquad d=50\ \mathrm{days}.$$
#
# $\mu_E$ 控制背景份量。它與 $\eta$ 的正規化一起決定背景和時變部分，不能直接把某次有限窗的 $\mu_E$ 解讀為觀測事件的背景比例。本書採等權重 EEPAS-NW，每顆輸入事件的額外除叢權重都取一。
#
# 50 天的 delay 是事件進入發報歷史前的等待規則。發報時還不夠舊的事件，整個窗都不加入。下一次發報再重新判斷；不是把每顆事件的時間核向右平移 50 天。本站 PPE 快取亦採相同資料延遲，方便保持背景一致。
#
# **前置時間**（lead time）則描述預報起點與可用歷史截止之間的距離。文獻有不同實驗定義，閱讀時要對照時間圖。目錄起點以前的缺失歷史是另一件事：即使目前資料即時可用，1960 年前的長尾貢獻仍然未知。
#
# 本書使用 1960 年起的暖機資料，但早期小震不完整。$\Delta(m)$ 處理的是模型中的規模截斷，沒有自動處理不同年代的漏測。較長的時間核尤其需要早期歷史，這是解讀義大利預報時必須保留的限制。
#
# ## 15.4 參數從哪裡來？
#
# 本書採已發表的 EEPAS-NW 參數，沒有用測試期重估。規模組為 $(1.22,1,0.25)$，時間組為 $(2.55,0.35,0.150)$，空間組為 $(0.52,1.00)$，背景權重為 $0.18$。順序分別對應上文的截距、斜率與散布。
#
# 原研究採分階段擬合：先固定部分尺度形狀，估其餘參數；再交換固定與估計的部分；最後聯合調整。這能提供較好的起始值，仍不保證找到全域最佳。不同參數可能互相補償，因此參數差異不必等於預報差異。
#
# 下圖從快取讀取同一個三個月窗的 EEPAS 與 PPE。顏色是各格加總所有目標規模箱的期望數。先讀活動放在哪裡，以及時間窗有多長；兩張圖的色階若自動縮放，不能只比藍色深淺判斷總數。

# %% tags=["remove-input"]
for model in ('EEPAS', 'PPE'):
    if not (italy.FORECAST_DIR / f'{model}_testing.npy').exists():
        raise FileNotFoundError(f'{model} 預報快取缺失；請先完成全書預算')
eepas=italy_models.get_forecast('EEPAS','testing')
ppe=italy_models.get_forecast('PPE','testing')
windows=italy.forecast_windows()
cells=italy.testing_cells()
cat=italy.experiment_catalog()
targets=italy.target_events(cat,'testing')
w=windows.iloc[0]
fig=plot_forecast_map(eepas[0].sum(axis=1),cells=cells,title=f'EEPAS：{w.start:%Y-%m-%d} 起，期望數 {eepas[0].sum():.3f}')
fig

# %% tags=["remove-input"]
fig=plot_forecast_map(ppe[0].sum(axis=1),cells=cells,title=f'PPE：同一窗，期望數 {ppe[0].sum():.3f}')
fig

# %% [markdown]
# ## 15.5 回到同一份標準預報
#
# 十年合計圖將各次實際發報的格箱期望數相加。每次使用當時已知的歷史，因此這不是 2012 年一次發出的十年預報。歷史在更新，已發表參數維持固定；兩種「固定」必須分開讀。
#
# 紅點是測試期真實目標事件。此處只疊圖，不評分；評分在第 16–18 章。事件靠近高率格值得注意，卻還不足以判定模型勝出。正式比較還會考慮總數、規模分配及沒有事件的格子。

# %% tags=["remove-input"]
fig=plot_forecast_map(eepas.sum(axis=(0,2)),cells=cells,targets=targets,title=f'EEPAS 滾動窗合計：期望 {eepas.sum():.2f}；目標 {len(targets)} 顆（只疊圖）')
fig

# %% [markdown]
# ## 15.6 模型名稱背後，哪些事尚未知道？
#
# 單顆輸入事件的規模核有一個中心，不表示它已與某顆未來大震配對。許多輸入事件的核會重疊，同一個未來位置可能收到很多微小貢獻。模型可以因此產生可比較的預報，仍不會給出確定的前兆名單。
#
# 三尺度參數也不能只靠地圖反推。時間核變寬、空間核變窄，有時能形成相似的格箱總數。閱讀參數表時，要同時看採用的時間單位、目標規模、資料期間及固定項目；不同地區的數值不能直接互換。
#
# 第 16 章開始才用同一份目標目錄檢查四張預報。把規格固定下來，能讓比較聚焦於模型分配的期望數。若看過結果再改半徑、目標或參數，就應另留新的測試資料，而不能把修改後的同份比較當成原來的擬前瞻實驗。
#
# ## 本章填入的規格欄位
#
# EEPAS 與其他模型輸出相同形狀的期望數陣列。空間常態對方格使用解析積分，時間核使用累積分布差，規模箱則採細分中點積分。三種積分方法都服務同一件事：把密度轉成可以跟觀測計數對照的量。

# %% tags=["remove-input"]
display(HTML(italy.spec_card([('模型','EEPAS-NW；所有額外事件權重為 1'),('目錄／區域','HORUS；深度 ≤ 40 km；S 為 CPTI15；R 為 177 格'),('三期','暖機 1960–1989；學習 1990–2011；測試 2012–2021'),('輸入／目標','m₀ = 2.45；mT = 5.0；25 個規模箱'),('發報','每窗 91.31 天；歷史截止於窗首前 50 天'),('參數','Biondini 2023 表 3；測試期不重估'),('輸出','窗 × 177 格 × 25 箱的期望數'),('檢驗','待填（第 16–18 章）')])))

# %% [markdown]
# 下一章：{doc}`16_test_number`。公式推導見 {doc}`appendix_d_eepas`。
#
# ### 再讀延伸論文時，要分開哪三件事？
#
# **時間完整度**（temporal completeness）描述已記錄歷史保留多少預期前兆貢獻。
# 它與規模門檻補償不同；需要把各種可能輸入規模的時間核一起平均。
# 前置時間較短，大目標所需的早期歷史可能還不夠。
# 附錄 D 用原文的完整度比例，說明補到背景或放大時變訊號的差別。
#
# **目錄前置時間**（catalogue lead time）從目錄起點量到目標地震。
# **預報期限**（forecast horizon）從發報量到預報窗末。
# **輸入延遲**（input delay）由發報往回推資料截止，本站是 50 天。
# **歷史時間差**（time lag）從資料截止量到窗內目標。
# 因此 time lag 等於目標距發報的時間再加 50 天。
# 它不是對每顆目標都固定為 50 天，四個時間量須分開。
#
# **EAS 餘震延伸**把預報主震可能帶來的餘震也積分進來。
# 它不同於 EEPAS-W 對已知輸入餘震降權，也不同於只計已知事件第一代的 ETAS 近似。
# 本站仍使用等權 EEPAS-NW，沒有把這些延伸當成已執行的結果。
# 三種結構的公式與適用限制見 {doc}`appendix_d_eepas`。
#
# ## 參考資料與延伸閱讀
#
# - Rhoades、Rastin 與 Christophersen（2022），[EEPAS 二十年回顧](https://doi.org/10.3390/geosciences12090349)。免費全文；先讀三尺度如何變成預報，再看資料限制。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 應用](https://doi.org/10.1093/gji/ggad123)。[免費機構典藏](https://www.earth-prints.org/handle/2122/17084)。表 3 核對 NW 參數，附錄 A 核對正規化。
# - Rhoades 等（2020），[目錄前置時間的影響](https://doi.org/10.3390/e22111264)。免費全文；區分預報起點、歷史截止及缺失前兆。
# - [PyEEPAS 公開程式](https://github.com/phonchi/EEPAS)。對照公開義大利資料與實作；本章數值參數來自上述已發表研究。
