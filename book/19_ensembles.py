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
# # 19. 把模型加起來：凸組合
#
# 一格由 ETAS 給出 0.2，EEPAS 給出 0.6。
# 若兩者各占一半，混合期望數就是 0.4。
# 這是合成數字，先讓相加的規則看得見。
# 本章再把相同規則用到義大利快取預報。
#
# 第 18 章的事件差值各有起伏。
# 兩個模型可能在不同事件上得分較好。
# 集合預報（ensemble forecast）會結合多個成分模型。
# 本章只採期望數的加法混合。
# 重點是加什麼、如何選權重，以及怎麼驗證。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, Markdown, display
from gdms_toolkit import italy, italy_models, csep_teaching as csep
from gdms_toolkit.viz import ACCENT, QUAKE_COLOR, apply_layout, show_diagram

models = ("SUP", "PPE", "EEPAS", "ETAS")
colors = dict(zip(models, [ACCENT, "#008300", "#4a3aa7", "#eda100"]))
cat = italy.experiment_catalog()
targets = italy.target_events(cat, "testing")
windows = italy.forecast_windows()
omega = italy.bin_targets(targets, windows)
forecasts = {m: italy_models.get_forecast(m, "testing") for m in models}
assert all(a.shape == omega.shape for a in forecasts.values())
assert int(omega.sum()) == len(targets)

# %% [markdown]
# ## 19.1 每個格箱，都用同一個加法
#
# $\pi=0.5$ 表示 ETAS 占一半權重。
# 權重（weight）是成分參與混合的比例。
# 凸組合（convex combination）要求權重非負。
# 所有成分的權重相加還必須等於一。
# 兩模型的式子因此很簡單。
#
# $$\Lambda_{\pi,b}=\pi\Lambda_{{\rm ETAS},b}
# +(1-\pi)\Lambda_{{\rm EEPAS},b},\qquad 0\le\pi\le1.$$
#
# $b$ 代表同一個窗、格與規模箱。
# 相加前，兩張預報必須使用相同單位。
# 窗長、區域邊界與規模範圍也必須相同。
# 一張每年的率，不能直接加到每窗期望數。
# 本章的快取已共用相同格箱。
#
# $\pi=0$ 完整保留 EEPAS。
# $\pi=1$ 完整保留 ETAS。
# 中間值在兩張期望數陣列之間移動。
# 每個格箱的混合值都落在兩個成分之間。
# 凸組合不會憑空造出更高的單箱率。

# %% tags=["remove-input"]
show_diagram("d19_convex_mix", caption=f"示意：混合的是 {len(['ETAS', 'EEPAS'])} 張期望數預報；權重沿同一組格箱套用。")

# %% [markdown]
# ## 19.2 「守恆」到底保留什麼？
#
# 把每格的混合率相加，就得到混合總量。
# 線性相加讓總量遵守同一個加權關係。
# 這是期望數的守恆關係。
# 它沒有要求兩個成分原本的總量相同。
#
# $$N_\pi=\pi N_{\rm ETAS}+(1-\pi)N_{\rm EEPAS}.$$
#
# $N_{\rm ETAS}=N_{\rm EEPAS}=N_0$ 時，
# 每個權重都保留共同總量 $N_0$。
# 如果兩個總量不同，混合總量也會隨權重變。
# 因此，「權重加起來等於一」只完成一半條件。
# 要固定總量，成分還要先有相同總量。
#
# $\Lambda_b/\sum_a\Lambda_a$ 可分開表示分配形狀。
# 正規化（normalization）把總和縮放到指定值。
# 若先各自縮放到共同 $N_0$，
# 混合只會改變格箱分配，總量保持固定。
# 那會形成另一個預報方案，應另作記錄。
#
# 下圖保留快取的原始期望數。
# 本章沒有把測試期觀測數拿來重新縮放。
# 紅線只提供事後數量對照。
# 各成分若都低於它，凸組合也會低於它。
# 混合能改變分配，但無法突破總量的上下界。

# %% tags=["remove-input"]
weights = np.linspace(0, 1, 51)
test_totals = np.array([p*forecasts["ETAS"].sum()+(1-p)*forecasts["EEPAS"].sum() for p in weights])
fig = go.Figure(go.Scatter(x=weights, y=test_totals, mode="lines", line=dict(color=ACCENT), name="混合期望數"))
fig.add_hline(y=int(omega.sum()), line_color=QUAKE_COLOR, line_dash="dash")
apply_layout(fig, title="原始率混合：總量在兩個端點之間", xaxis_title="ETAS 權重 π", yaxis_title="測試期目標顆數")
display(Markdown(f"圖：兩端期望數為 {test_totals[0]:.2f} 與 {test_totals[-1]:.2f}；紅線是 {int(omega.sum())} 顆觀測。"))
fig

# %% [markdown]
# 一個格箱的事件機率又有不同運算。
# Poisson 下至少一顆的機率為 $1-e^{-\Lambda}$。
# 應先混合期望數，再依這個公式換算。
# 直接平均兩個成分的機率，會得到另一個數字。
# 原因是指數轉換並非線性。
#
# 一份完整預報還可以混合兩個機率分布。
# 那相當於先按權重選一個成分模型。
# 再由被選中的模型生成整份目錄。
# 這種分布混合通常也不等於單一 Poisson。
# 本章只混合率，沿用既定的計數檢驗。
#
# ## 19.3 相加為何可能改善對數分數？
#
# 一顆事件落在某模型的極低率箱，
# 會給那個模型很大的對數扣分。
# 如果另一模型在此配置較高的率，
# 混合便能減輕這次扣分。
# 這是成分互補（complementarity）的具體例子。
#
# 對數曲線在兩個端點間向上彎。
# 凹性（concavity）指曲線高於端點連線。
# 對正率而言，下式給出混合的分數下界。
#
# $$\ln[\pi a+(1-\pi)b]\ge\pi\ln a+(1-\pi)\ln b.$$
#
# jPOLL 的負總率項是線性的。
# 固定觀測的階乘項也不隨權重改變。
# 所以同一份觀測上的混合分數滿足：
#
# $$\operatorname{jPOLL}(\Lambda_\pi,\omega)
# \ge\pi\operatorname{jPOLL}(\Lambda_{\rm ETAS},\omega)
# +(1-\pi)\operatorname{jPOLL}(\Lambda_{\rm EEPAS},\omega).$$
#
# 式中的兩個成分，都使用同一個觀測陣列。
# 右邊是兩個分數的加權平均。
# 它通常低於其中較高的那個分數。
# 因此，下界不保證混合勝過最佳成分。
# 讀公式時要看清楚比較的是哪個對象。

# %% [markdown]
# 兩張完全相同的預報沒有互補空間。
# 它們混合後仍是同一張圖。
# 兩張都漏掉相同地區，也很難靠相加補救。
# 互補要在不同事件或不同區域上有證據。
# 模型名稱不同，本身不構成互補證明。
#
# ETAS 會利用最近已知事件的觸發貢獻。
# EEPAS 則利用來源事件的較長時間尺度。
# 兩種訊息的作用時間可能不同。
# 本站 ETAS 仍是背景加第一代的近似。
# 本章的混合也一併繼承這項近似。
#
# ## 19.4 用學習期選好，才打開測試期
#
# 1990–2011 是這次選權重使用的資料。
# 候選權重事先列在固定網格上。
# 每個權重都計算學習期 jPOLL。
# 最高分的候選值記作固定權重。
# 若同分，本章採最小的 ETAS 權重。
#
# 2012–2021 的觀測不參與這個選擇。
# 測試曲線只用來看事後表現。
# 事後最佳（oracle optimum）指看過答案才選出的高點。
# 它可以展示選擇的樂觀程度。
# 它無法冒充當時能發布的預報。
#
# 下圖同時畫學習與測試兩條曲線。
# 兩者都換成相對 EEPAS 的每事件分數差。
# 這樣可減少不同期間顆數造成的尺度差。
# 但曲線高度仍反映不同的事件內容。
# 真正要追的是同一個固定權重在兩期的位置。

# %% tags=["remove-input"]
learning_windows = italy_models.learning_windows()
learning_targets = italy.target_events(cat, "learning")
learning_omega = italy.bin_targets(learning_targets, learning_windows)
learning_rates = {m: italy_models.get_forecast(m, "learning") for m in ("ETAS", "EEPAS")}
assert all(a.shape == learning_omega.shape for a in learning_rates.values())
assert int(learning_omega.sum()) == len(learning_targets)
learning_scores = np.array([csep.jpoll(p*learning_rates["ETAS"]+(1-p)*learning_rates["EEPAS"], learning_omega) for p in weights])
best_index = int(np.argmax(learning_scores))
chosen_weight = float(weights[best_index])
# The selection above uses learning data only; testing is diagnostic from here on.
testing_scores = np.array([csep.jpoll(p*forecasts["ETAS"]+(1-p)*forecasts["EEPAS"], omega) for p in weights])
learning_gain = (learning_scores-learning_scores[0])/learning_omega.sum()
testing_gain = (testing_scores-testing_scores[0])/omega.sum()
assert np.all(learning_scores >= (1-weights)*learning_scores[0]+weights*learning_scores[-1]-1e-8)
fig = go.Figure()
fig.add_scatter(x=weights, y=learning_gain, name="學習期", line=dict(color=ACCENT))
fig.add_scatter(x=weights, y=testing_gain, name="測試期（事後展示）", line=dict(color="#4a3aa7", dash="dash"))
fig.add_vline(x=chosen_weight, line_dash="dot", line_color="#666666")
apply_layout(fig, title="固定候選網格：在哪一期選權重？", xaxis_title="ETAS 權重 π", yaxis_title="相對 EEPAS 的每事件分數差")
display(Markdown(f"圖：共 {len(weights)} 個候選權重。學習期 {int(learning_omega.sum())} 顆選出 π={chosen_weight:.2f}；測試期 {int(omega.sum())} 顆只用於評估。"))
fig

# %% [markdown]
# 學習曲線的最高點回答擬合問題。
# 它已使用成分參數的學習資料。
# 所以即使權重只用學習期選出，
# 學習成績仍是偏樂觀的回溯結果。
# 評估選權重程序，需要新的事件資料。
#
# 測試曲線若在另一處最高，
# 只表示兩期偏好的組合可能不同。
# 它也可能受到有限事件的波動影響。
# 本章保留原先選出的權重，不追著高點改。
# 若要更新選擇規則，下一期再作獨立評估。
#
# 下表只列兩個端點與學習期選擇。
# 端點讓讀者看見是否超過單一成分。
# 表中的期望數則保留總量是否低估的線索。
# 較高資訊分數可以與數量不足同時存在。
# 因此，混合模型仍需前兩章的一致性檢驗。

# %% tags=["remove-input"]
mixed = chosen_weight*forecasts["ETAS"]+(1-chosen_weight)*forecasts["EEPAS"]
assert np.isclose(mixed.sum(), chosen_weight*forecasts["ETAS"].sum()+(1-chosen_weight)*forecasts["EEPAS"].sum())
display(HTML(pd.DataFrame([
    {"方案": label, "ETAS 權重": p, "測試期望數": float(a.sum()),
     "對 EEPAS 的 IGPE": csep.information_gain(a, forecasts["EEPAS"], omega)["igpe"]}
    for label, p, a in [("EEPAS", 0, forecasts["EEPAS"]), ("學習期選擇", chosen_weight, mixed), ("ETAS", 1, forecasts["ETAS"])]
]).to_html(index=False, float_format=lambda x: f"{x:.3f}")))

# %% tags=["remove-input"]
from gdms_toolkit.viz import plot_forecast_map
fig = plot_forecast_map(mixed.sum(axis=(0, 2)), targets=targets, title="固定學習期權重的混合預報")
display(Markdown(f"圖：π={chosen_weight:.2f}，十年總期望 {mixed.sum():.2f} 顆。藍色為每格期望數；{len(targets)} 個紅點為事後目標疊圖。此圖本身不評分。"))
fig

# %% [markdown]
# ## 19.5 前瞻退步，會推翻哪一個承諾？
#
# Bayona 等人的前瞻研究提供一個提醒。
# 部分乘法混合曾在回溯期改善分數。
# 到了新一期，優勢未能維持。
# 乘法混合（multiplicative hybrid）以相乘結合率場。
# 該研究的組合法與本章的凸組合不同。
# 但兩者都必須面對時間外的資料。
#
# [Bayona 等（2022）](https://doi.org/10.1093/gji/ggac018)
# 檢查的是預先提出模型的前瞻表現。
# 這項結果提醒讀者，擬合改善可能依賴特定期間。
# 它沒有證明所有混合都必然失敗。
# 也不能直接當作本站加法混合的檢驗結果。
#
# 一個正式方案應一起固定成分版本與候選集合。
# 還要固定選擇分數、學習期間與更新時點。
# 只宣布最後權重，無法重建選擇過程。
# 若參數、網格與門檻也一起挑選，
# 整套挑選程序都需要樣本外檢查。
#
# 一張混合圖仍有清楚的用途範圍。
# 它給出事件在格箱內的期望數。
# 使用者若要知道建築會搖多大，
# 還需要地震到地動的轉換。
# 下一章先讀懂原論文，再補上這一段連結。
#
# 固定權重也讓讀者能追查每格的來源。
# 某格變深，可以回看哪個成分原本較高。
# 某格變淡，則表示另一成分把加權平均拉低。
# 這種解釋只涉及預報率的組成。
# 它沒有替地震指定唯一的物理成因。
#
# 一個學到端點的權重也應照實保留。
# 端點代表這批學習資料偏好單一成分。
# 混合方法允許這個答案。
# 為了讓圖看起來像混合而移動權重，
# 反而會破壞事先公布的選擇規則。
#
# ## 19.6 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("成分", "本站 ETAS 第一代近似＋等權 EEPAS"),
    ("混合對象", "相同窗、格、規模箱內的期望數；原始總量"),
    ("權重條件", "0≤π≤1；EEPAS 權重為 1−π"),
    ("選擇資料", "1990–2011 學習期；成分參數也來自學習期"),
    ("選擇規則", f"{len(weights)} 個等距候選；最大 jPOLL；同分取最小 π"),
    ("固定權重", f"ETAS π={chosen_weight:.2f}"),
    ("評估資料", "2012–2021；測試曲線只作事後展示"),
    ("比較", "對 EEPAS 的 IGPE；與兩端點並列"),
])))

# %% [markdown]
# 混合完成後，仍要清楚說出它預報的量。
# 下一章：{doc}`20. 預報之後 <20_beyond_forecast>`。
#
# ## 參考資料與延伸閱讀
#
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。
#   [免費機構全文](https://www.earth-prints.org/handle/2122/17084)。先對照研究區與測試期，再讀檢驗。
# - Rhoades 與 Gerstenberger（2009），[短期預報的混合模型](https://doi.org/10.1785/0120080063)。
#   [免費摘要](https://central.scec.org/node/3964)。從不同時間尺度理解加法組合。
# - Bayona 等（2022），[乘法混合的前瞻評估](https://doi.org/10.1093/gji/ggac018)。
#   免費全文。重點讀回溯改善與前瞻表現的落差。
