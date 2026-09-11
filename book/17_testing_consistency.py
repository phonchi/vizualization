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
# # 12. 預報檢驗 I：觀測是否符合模型
#
# {doc}`EEPAS 與 PPE <16_eepas_ppe>`讓我們看到同一份目錄可以產生不同預報。現在假設發布時間已到：模型都已交出未來一年各地、各規模的預期事件數。一年後真正觀測到的地震，應該如何與預報比較？
#
# 只看有沒有一場地震落在高率區並不夠。模型可能把全區的率都設得很高，於是任何事件都看似命中；也可能猜對總數，卻放錯位置。一致性檢驗將這些問題拆開，分別檢查總數與分布的差異。下一章才進一步問兩個模型誰提供較多資訊。
#
# ## 12.1 評分之前，先固定題目
#
# 前面 Ψ 的例子顯示，分析者在看過結果之後仍可改時間窗和區域，就很難判斷成功來自方法還是選擇。預報測試因此要先約定目標規模、時間與空間範圍、目錄來源和評分方式。
#
# CSEP 將這個想法做成可重現的測試程序。模型可以隨新資料更新，但更新規則也必須事先固定；「前瞻」允許參數更新，但每次更新都不能使用當時尚未取得的資訊。這個原則讓 ETAS、EEPAS 和較簡單的背景模型可以接受同樣的檢查。
#
# ## 12.2 一張預報圖說了多少事
#
# 網格化預報對空間格 $j$、規模箱 $k$ 提供期望事件數 $\Lambda_{jk}$，觀測則是計數 $\omega_{jk}$。如果某格預報 0.2 顆，不是說會觀測到五分之一顆地震，而是說在同樣條件下重複很多次，平均數為 0.2。
#
# 平均數沒有決定全部機率。每次都接近平均、或大多數時候沒有但偶爾成群，都可能有同樣期望數。若模型只交出一張均值圖，評分時還要指定計數分布和格子間的相依性。
#
# 另一種作法是交出多份模擬目錄。對 ETAS 而言，每次模擬都讓未來事件繼續觸發後代，便能保留事件數波動、時間先後和空間叢集。把模擬目錄再平均為一張圖，會失去其中部分資訊。
#
# 對歷史相依模型，期望數的概念可寫成
#
# $$\Lambda_{jk}=\mathbb E\!\left[\int_{T_0}^{T_1}\int_{S_j}\int_{m_k}^{m_{k+1}}
# \lambda^*(t,x,y,m)\,dm\,dx\,dy\,dt\ \middle|\ H_{T_0}\right].$$
#
# 外面的期望是對發布當下尚未知的歷史平均。它提醒我們：在觀測路徑上積分出來的補償子，不能直接當成發布時已知的固定期望數。
#
# ## 12.3 為什麼未發生的事件也會影響分數
#
# 先採一個便於理解的假設：各箱獨立，箱內計數服從 Poisson 分布。此時整份觀測的對數概似為
#
# $$\mathrm{jPOLL}=\sum_{j,k}\left[-\Lambda_{jk}+\omega_{jk}\ln\Lambda_{jk}
# -\ln(\omega_{jk}!)\right].$$ (eq:jpoll)
#
# 有事件的箱會透過 $\omega\ln\Lambda$ 區分高率與低率；每一個箱又都要付出 $-\Lambda$。因此把全區率一律抬高，不會無條件得到好分數。空白區域也是資料的一部分。
#
# 階乘項只跟觀測有關，在同一份目錄的兩模型比較中會抵消；但當我們拿觀測與很多不同模擬目錄比較時，每份目錄的計數不同，它就不能隨便省略。分數的絕對值也依賴網格和單位，應在相同設定內解讀；除以事件數仍不能自動消除網格尺度。
#
# ## 12.4 先問總數是否合理
#
# 如果模型預期一年約二十顆，卻看到兩百顆，可以先檢查總量。N-test 將所有格子相加，只看總事件數，暫時忽略位置和規模。在獨立 Poisson 設定下，總數仍為 Poisson，其均值是 $\sum_{j,k}\Lambda_{jk}$。
#
# 下面把模型的總數分布與觀測位置放在一起。落在中央表示這個總數不特別意外，落在尾部則表示模型高估或低估活動量。這只檢查一個面向，還沒有告訴我們地圖畫得對不對。
#
# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from scipy.special import gammaln

from gdms_toolkit.viz import (ACCENT, PALETTE, QUAKE_COLOR, SEQUENTIAL,
                              apply_layout)

N_fore = 20.0                                   # 模型的預報總數
kk = np.arange(0, 46)
pmf_pois = stats.poisson.pmf(kk, N_fore)
lo, hi = stats.poisson.ppf([0.025, 0.975], N_fore)

obs_a, obs_b = 24, 38                           # 兩個假想的觀測總數
d1_a = 1.0 - stats.poisson.cdf(obs_a - 1, N_fore)
d1_b = 1.0 - stats.poisson.cdf(obs_b - 1, N_fore)

fig = go.Figure(go.Bar(x=kk, y=pmf_pois, marker_color=ACCENT, opacity=0.75,
                       name=f"Poisson({N_fore:.0f}) 預報數分布"))
fig.add_vrect(x0=lo - 0.5, x1=hi + 0.5, fillcolor=PALETTE[2], opacity=0.10,
              line_width=0,
              annotation_text=f"95% 區間 [{lo:.0f}, {hi:.0f}]")
fig.add_vline(x=obs_a, line_color=PALETTE[2], line_width=3,
              annotation_text=f"觀測 A={obs_a}（δ₁={d1_a:.2f}，一致）")
fig.add_vline(x=obs_b, line_color=QUAKE_COLOR, line_width=3,
              annotation_text=f"觀測 B={obs_b}（δ₁={d1_b:.1e}，不一致）")
apply_layout(fig, title="N-test：觀測總數落在模型預報數分布的哪裡",
             xaxis_title="預報期間的目標地震數",
             yaxis_title="機率", hovermode="x")
fig

# %% [markdown]
# ## 12.5 總數相同，地圖仍可能不同
#
# 為了單看空間形狀，可以把每格率除以全區總率，得到事件落在各格的比例，再固定模擬事件總數為實際觀測數。這樣做刻意拿掉模型總量的影響，讓高分或低分主要來自位置安排。
#
# S-test 先對規模求和，檢查空間；M-test 先對空間求和，檢查規模分布；條件概似 cL-test 保留空間與規模的聯合分布，同樣條件在觀測總數上。這四種檢驗用來辨認偏差來源，通過它們仍不保證模型有用。
#
# | 檢驗 | 保留的問題 | 暫時不回答的問題 |
# |---|---|---|
# | N | 發生多少事件 | 事件在哪裡、多大 |
# | S | 事件落在哪裡 | 總量與規模分布 |
# | M | 規模如何分布 | 總量與位置 |
# | cL | 位置和規模的聯合安排 | 總量是否正確 |
#
# 分數高低需要一個比較尺度。我們在模型下產生很多份目錄，對每份目錄用同一規則評分，再把實際觀測分數放進這個分布。下面示範這個程序：若觀測分數落在低端，表示模型自己很少產生這麼不符合預報形狀的資料。
#
# 分位數不是「模型為真的機率」。它描述在指定模型和模擬規則下，觀測統計量的位置；換了統計量、網格或模擬模型，也可能改變它。
#
# %% tags=["remove-input"]
rng = np.random.default_rng(17)

nx, ny = 30, 20                                  # 600 個空間格
gx = np.linspace(120.0, 122.0, nx)
gy = np.linspace(21.9, 25.3, ny)
GX, GY = np.meshgrid(gx, gy, indexing="ij")


def blob(x0, y0, s, amp):
    return amp * np.exp(-((GX - x0) ** 2 + (GY - y0) ** 2) / (2 * s ** 2))


bg = 0.05
lam_model = bg + blob(121.5, 24.0, 0.30, 1.0) + blob(120.7, 22.8, 0.25, 0.8)
lam_true = lam_model + blob(121.8, 22.6, 0.17, 2.4)   # 模型沒預期到的一團
lam_model, lam_true = lam_model.ravel(), lam_true.ravel()

N_obs = 60
p_true = lam_true / lam_true.sum()
omega = rng.multinomial(N_obs, p_true)           # 一份「觀測」

lam_s = lam_model * N_obs / lam_model.sum()      # S-test 的正規化率場
log_lam_s = np.log(lam_s)
const = lam_s.sum()


def jpoll(counts):
    return counts @ log_lam_s - const - gammaln(counts + 1).sum(axis=-1)


n_sim = 10_000
sims = rng.multinomial(N_obs, lam_s / lam_s.sum(), size=n_sim)
scores_sim = jpoll(sims)
score_obs = jpoll(omega)
gamma = float(np.mean(scores_sim <= score_obs))

fig = go.Figure(go.Histogram(x=scores_sim, nbinsx=60, marker_color=ACCENT,
                             opacity=0.8, name="模擬分數（一萬份合成目錄）"))
fig.add_vline(x=score_obs, line_color=QUAKE_COLOR, line_width=3,
              annotation_text=f"觀測分數 {score_obs:.1f}（γ={gamma:.4f}，"
                              f"n={n_sim:,}）")
fig.add_vline(x=float(np.quantile(scores_sim, 0.05)), line_dash="dash",
              line_color="#888", annotation_text="5% 分位")
apply_layout(fig,
             title=f"S-test：觀測分數落在模擬分布的第 {gamma * 100:.2f} 百分位",
             xaxis_title="空間聯合對數概似 jPOLL",
             yaxis_title="模擬次數", hovermode="x", showlegend=False)
fig

# %% [markdown]
# 整體分數也可以拆回各格。接下來的圖把每格的貢獻分開，讓我們看到是否只有少數叢集區主導結果。這對地震尤其重要：一個強餘震序列可能在局部區域留下很多事件，卻未必代表所有區域的空間模式都同樣失準。
#
# %% tags=["remove-input"]
cell_obs = omega * log_lam_s - lam_s - gammaln(omega + 1)
cell_exp = sims.mean(axis=0) * log_lam_s - lam_s - gammaln(sims + 1).mean(axis=0)
deficit = cell_exp - cell_obs                    # 每格「拖累」總分多少
order = np.argsort(-deficit)
share = np.cumsum(deficit[order]) / deficit[deficit > 0].sum()
top5, top20 = share[4] * 100, share[19] * 100

fig = make_subplots(rows=1, cols=2, column_widths=[0.58, 0.42],
                    subplot_titles=("逐格分數缺口的空間分布",
                                    "累積貢獻：少數格主宰總分"))
fig.add_trace(go.Heatmap(z=deficit.reshape(nx, ny).T, x=gx, y=gy,
                         colorscale=SEQUENTIAL, colorbar=dict(x=0.46,
                                                              title="缺口"),
                         hovertemplate="經度 %{x:.2f}<br>緯度 %{y:.2f}"
                                       "<br>缺口 %{z:.2f}<extra></extra>"),
              row=1, col=1)
hit = omega > 0
fig.add_trace(go.Scatter(x=GX.ravel()[hit], y=GY.ravel()[hit], mode="markers",
                         marker=dict(color=QUAKE_COLOR, size=5,
                                     line=dict(width=0)),
                         name="有事件的格", showlegend=False),
              row=1, col=1)
fig.add_trace(go.Scatter(x=np.arange(1, 41), y=share[:40] * 100, mode="lines",
                         line=dict(color=ACCENT, width=2.5), showlegend=False),
              row=1, col=2)
fig.add_hline(y=top20, line_dash="dash", line_color=QUAKE_COLOR,
              annotation_text=f"前 5 格佔 {top5:.0f}%，"
                              f"前 20 格佔 {top20:.0f}%", row=1, col=2)
fig.update_xaxes(title_text="經度", row=1, col=1)
fig.update_yaxes(title_text="緯度", row=1, col=1)
fig.update_xaxes(title_text="按缺口排序的格子名次", row=1, col=2)
fig.update_yaxes(title_text="累積佔正缺口的百分比", row=1, col=2)
apply_layout(fig, title="總分的來源極度不均：幾格決定了整個檢驗結果",
             hovermode="closest", showlegend=False)
fig

# %% [markdown]
# ## 12.6 當事件數比 Poisson 更不穩定
#
# 前面叢集章已經提醒我們：平均率相同的序列，波動程度仍可能不同。Poisson 把變異數固定為均值，對成群活動可能給出太窄的計數區間。因此極端事件數可能反映分布假設不足，而不只是平均率錯誤。
#
# 負二項分布允許均值 $\Lambda$ 與額外離散程度分開調整。常用形式為
#
# $$\mathbb E[N]=\Lambda,\qquad \mathrm{Var}(N)=\Lambda+\frac{\Lambda^2}{r},\quad r>0.$$
#
# $r$ 越大越接近 Poisson，$r$ 較小則允許更大波動。下面在同一平均數下比較兩種分布：先看中央，再看尾部。同一個觀測數，在較寬分布裡不一定那麼意外。
#
# 額外離散參數也需要資料或模型支持，不能看過某次極端結果後才放寬到讓它透過。負二項改善的是邊際計數分布，並沒有自動恢復箱間相關或完整的觸發歷史；這些資訊仍可能需要目錄式模擬。
#
# %% tags=["remove-input"]
mu_nb = N_fore                                   # 與 N-test 圖同一個平均
sigma2_nb = 80.0                                 # 由歷史不重疊時段估得
nu_nb = mu_nb / sigma2_nb
tau_nb = mu_nb ** 2 / (sigma2_nb - mu_nb)
kk2 = np.arange(0, 61)
pmf_nb = stats.nbinom.pmf(kk2, tau_nb, nu_nb)
lo_nb, hi_nb = stats.nbinom.ppf([0.025, 0.975], tau_nb, nu_nb)

fig = go.Figure()
fig.add_trace(go.Bar(x=kk2, y=stats.poisson.pmf(kk2, mu_nb), opacity=0.65,
                     marker_color=ACCENT,
                     name=f"Poisson：μ={mu_nb:.0f}, σ²={mu_nb:.0f}"))
fig.add_trace(go.Bar(x=kk2, y=pmf_nb, opacity=0.65, marker_color=PALETTE[3],
                     name=f"負二項：μ={mu_nb:.0f}, σ²={sigma2_nb:.0f}"
                          f"（τ={tau_nb:.2f}, ν={nu_nb:.2f}）"))
fig.add_vline(x=obs_b, line_color=QUAKE_COLOR, line_width=3,
              annotation_text=f"觀測 B={obs_b}")
apply_layout(fig, title=f"同一個平均、不同的變異數："
                        f"Poisson 95% 區間 [{lo:.0f}, {hi:.0f}]，"
                        f"負二項 [{lo_nb:.0f}, {hi_nb:.0f}]",
             xaxis_title="預報期間的目標地震數", yaxis_title="機率",
             barmode="overlay", hovermode="x")
fig

# %% [markdown]
# ## 12.7 如果問題只是「這格有沒有事件」
#
# 另一種評估目標是不計重複次數，只記一格是否至少有一次事件。觀測變成 $y_{jk}=1$ 或 $0$，預報也必須轉為相同事件的機率 $P_{jk}$。
#
# 若該格計數確實服從 Poisson，則 $P_{jk}=1-e^{-\Lambda_{jk}}$，可用二元對數分數
#
# $$y_{jk}\ln P_{jk}+(1-y_{jk})\ln(1-P_{jk})$$
#
# 評估。這會降低一格內很多顆餘震對分數的重複影響，但代價是放棄計數資訊：一顆和一百顆都記為一。
#
# 下面比較相同預報率下，不同觀測數對兩種分數的影響。二元化改變了預報目標，不能把較容易通過檢驗當成原模型改善的證據。對一般自激發計數，不能只把均值代入 $1-e^{-\Lambda}$；至少一次機率應由該模型的計數分布或模擬估計。
#
# %% tags=["remove-input"]
lam_grid = np.logspace(-5, 0, 200)
fig = go.Figure()
for idx, w in enumerate([0, 1, 2, 3, 5]):
    poll = -lam_grid + w * np.log(lam_grid) - gammaln(w + 1)
    bill = np.where(w >= 1, np.log1p(-np.exp(-lam_grid)), -lam_grid)
    fig.add_trace(go.Scatter(x=lam_grid, y=poll - bill, mode="lines",
                             name=f"ω={w}",
                             line=dict(color=PALETTE[idx], width=2.5)))
diff3 = float((-1e-3 + 3 * np.log(1e-3) - gammaln(4))
              - np.log1p(-np.exp(-1e-3)))
fig.add_annotation(x=np.log10(1e-3), y=diff3, text=f"Λ=10⁻³, ω=3：{diff3:.1f}",
                   showarrow=True, arrowhead=2, ax=60, ay=-30)
fig.update_xaxes(type="log")
apply_layout(fig, title="POLL 減 BILL：只有多震格才有實質差異",
             xaxis_title="該格的預報期望數 Λ",
             yaxis_title="POLL − BILL（對數單位）", hovermode="x")
fig

# %% [markdown]
# 二元評分也不會自動使格子彼此獨立。若要把逐格分數當統計量，可以用保留叢集的模擬來校準；若把乘積直接當聯合機率，則另有獨立假設需要說明。
#
# ## 12.8 把多個診斷放在一起
#
# 同時做許多檢驗，至少遇到一次小分位數的機會會增加。若要控制整個檢驗家族的誤拒機率，可以使用 Bonferroni 等預先指定的方法。Bonferroni 用實際檢驗數，不要求各檢驗獨立；S 與 cL 相關，並不能直接把四個檢驗改算成兩個。若要利用相關性減少保守程度，必須用聯合模擬或有依據的校正程序。
#
# 另一方面，多個檢驗都沒有拒絕，也不保證模型具備有用的預報能力。樣本太少時，差異很大的模型可能都透過。比起只列透過與否，更有教學意義的讀法是：數量、空間、規模各在哪些地方和資料不一致？這些差異是否集中在一段序列或一個區域？
#
# ```{admonition} 分數與模擬的計算細節
# :class: dropdown
#
# {doc}`附錄 E <appendix_e_testing>`列出條件多項分布、網格與連續概似的關係、Poisson 與二元分數、負二項參數及蒙地卡羅誤差。主文只需掌握：同一規則同時用在觀測和模擬，才有可比較的分位數。
# ```
#
# 現在假設兩個模型都與資料大致相容，我們仍想知道增加時間變化或新資料是否值得。{doc}`下一章 <18_testing_comparison>`將固定同一份觀測，比較兩個預報各自提供多少資訊，並進一步檢查校準與實際使用價值。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - [Theory of CSEP Tests](https://docs.cseptesting.org/getting_started/theory.html) — pyCSEP 開發團隊，官方文件（免費）。先對照各檢驗的目標、模擬方式與分位數分數，再看 N、S、M 與 conditional L-test 的程式範例，可把本章公式接到實際檢驗流程。
# - [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。從「觀測是否與預報一致」出發整理檢驗方法，適合先建立觀念，再回頭理解本章為何要拆開事件數、空間與規模。
# - [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。本章 Poisson／二元似然、負二項計數與多重檢定的主要實例來源；閱讀時留意餘震叢集如何影響各種評分。
#
# - [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。
#
