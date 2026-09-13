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
# # 18. 檢驗 III：哪張預報比較有資訊？
#
# 同一顆地震，在兩張預報上有兩個率。
# 它若落在 A 的高率區、B 的低率區，
# A 就在這顆事件上得到較好的分數。
# 本章把每顆事件配成一對，再合計差異。
#
# 第 17 章的分位數檢查各模型自身。
# 本章改問兩個模型之間的差距。
# 基準模型（baseline model）是比較的參照。
# 先用 SUP，再換成 PPE。
# 每個結論都會帶著基準名稱。

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
# ## 18.1 一顆事件，配一個對數率比
#
# $\Lambda_A=0.4$、$\Lambda_B=0.2$ 是合成率。
# 同一個箱若出現事件，A 的率是 B 的兩倍。
# 對數率比（log rate ratio）把這個倍數取對數。
# 它在兩率相同時等於零。
# 以下一律使用自然對數。
#
# $$d_i=\ln\frac{\Lambda_{A,b_i}}{\Lambda_{B,b_i}}.$$
#
# $b_i$ 是第 $i$ 顆目標地震所在的格箱。
# 配對（pairing）指兩個分數來自同一顆事件。
# 因此，每個差值都具有共同的觀測條件。
# 事件的時間、位置與規模完全相同。
# 改變的只有提供率的模型。
#
# $\ln 2$ 表示這顆事件偏向 A。
# $\ln(1/2)$ 則表示它偏向 B。
# 大小相反的兩個差值會互相抵銷。
# 這個刻度方便把事件貢獻相加。
# 它沒有把分數換成「命中百分比」。
#
# 下圖依窗、格與規模箱排序事件。
# 同箱若有多顆，每顆都保留相同差值。
# 這個排序與工具箱輸出一致。
# 它保留逐事件貢獻，尚未加上總率修正。

# %% tags=["remove-input"]
against_sup = {m: csep.information_gain(forecasts[m], forecasts["SUP"], omega) for m in models if m != "SUP"}
fig = go.Figure()
for m, r in against_sup.items():
    fig.add_scatter(x=np.arange(1, r["n"]+1), y=r["per_event"], mode="lines+markers",
                    name=f"{m} / SUP", line=dict(color=colors[m]))
fig.add_hline(y=0, line_dash="dash", line_color="#666666")
apply_layout(fig, title="逐事件對數率比：哪些事件拉開差距？", xaxis_title="格箱排序後的事件序號", yaxis_title="ln(模型率 / SUP 率)")
display(Markdown(f"圖：各線都有 {int(omega.sum())} 個目標事件貢獻。相同格箱中的事件會重複計入。"))
fig

# %% [markdown]
# ## 18.2 把空箱的代價也算回來
#
# 把 A 的所有率乘大，事件上的分數會增加。
# 但模型也在沒有事件的箱放進更多期望數。
# 第 17 章的 jPOLL 會為這部分扣分。
# 所以兩模型的比較必須保留總率項。
#
# 每事件資訊增益稱為 IGPE。
# 英文是 information gain per earthquake。
# 它是兩張預報的對數概似差，除以事件數。
# 正值表示 A 在指定觀測上比 B 得分高。
#
# $$\operatorname{IGPE}_{A:B}
# =\frac{\operatorname{jPOLL}_A-\operatorname{jPOLL}_B}{N}
# =\overline d+\frac{N_B-N_A}{N}.$$
#
# $N_A$ 與 $N_B$ 是兩模型的期望總數。
# $N$ 是觀測總數，$\overline d$ 是配對差的平均。
# 同一份觀測的階乘項相消。
# 剩下的總率修正可能增加或減少 IGPE。
# 它的方向由兩模型的總量決定。
#
# 下表把平均差與總率項分開。
# 讀者可以直接檢查兩欄相加是否等於 IGPE。
# 這能防止把「全部率抬高」誤讀為技巧。
# 率的絕對尺度與分配比例都參與評分。

# %% tags=["remove-input"]
score_rows = []
for m, r in against_sup.items():
    direct = (csep.jpoll(forecasts[m], omega)-csep.jpoll(forecasts["SUP"], omega))/omega.sum()
    assert np.isclose(direct, r["igpe"])
    score_rows.append({"模型 / SUP": m, "平均對數率比": r["per_event"].mean(), "總率修正": r["rate_term"], "IGPE": r["igpe"]})
display(HTML(pd.DataFrame(score_rows).to_html(index=False, float_format=lambda x: f"{x:.3f}")))

# %% [markdown]
# ## 18.3 同樣的平均，也可能有不同區間
#
# 兩組差值可以具有相同平均。
# 一組每顆都接近平均，另一組起伏很大。
# 後者對抽到哪些事件會更敏感。
# 配對 t 區間（paired t interval）用差值的離散程度估計精度。
# 本章把預報率視為已固定。
#
# $$s_d^2=\frac{1}{N-1}\sum_i(d_i-\overline d)^2,\qquad
# \operatorname{IGPE}\ \pm\ t_{0.975,N-1}\frac{s_d}{\sqrt N}.$$
#
# $s_d$ 是事件差值的樣本標準差。
# $t_{0.975,N-1}$ 是 Student t 分布的分位點。
# 固定的總率修正移動區間中心。
# 它不改變這個公式使用的區間寬度。
# 樣本少或差值偏斜時，近似會更不穩。
#
# 零若落在區間內，差異仍難確定。
# 整段位於零以上，則支持 A 得分較高。
# 這是固定率與近似獨立事件下的判讀。
# 同一序列的事件可能一起偏向某個模型。
# 下節會把這層相依性另列出來。

# %% tags=["remove-input"]
def gain_figure(results, baseline):
    gain_fig = go.Figure()
    for m, r in results.items():
        gain_fig.add_scatter(x=[r["igpe"]], y=[m], mode="markers", name=m,
            marker=dict(color=colors[m], size=10),
            error_x=dict(type="data", array=[r["ci"][1]-r["igpe"]],
                         arrayminus=[r["igpe"]-r["ci"][0]], symmetric=False, color=colors[m]))
    gain_fig.add_vline(x=0, line_dash="dash", line_color="#666666")
    return apply_layout(gain_fig, title=f"相對 {baseline} 的 IGPE 與配對 t 區間",
                        xaxis_title="每事件資訊增益（自然對數）", showlegend=False)

fig = gain_figure(against_sup, "SUP")
display(Markdown(f"圖：{int(omega.sum())} 顆事件的近似 95% 配對 t 區間。區間未校正序列相依與多重比較。"))
fig

# %% [markdown]
# ## 18.4 配對有幫助，事件相依仍要處理
#
# 一顆事件若很容易被兩模型預報，
# 兩個分數可能一起升高。
# 共變異（covariance）描述兩個量共同起伏的方向。
# 配對差會把這個共同部分算進去。
# 把兩組分數當獨立樣本，便會漏掉它。
#
# $$\operatorname{Var}(X-Y)=\operatorname{Var}(X)
# +\operatorname{Var}(Y)-2\operatorname{Cov}(X,Y).$$
#
# 兩顆餘震也可能一起支持 ETAS。
# 相依性（dependence）指事件資訊彼此相關。
# 配對處理的是同一事件的兩個模型分數。
# 它沒有自動處理不同事件之間的關係。
# 因此，「已配對」仍不足以保證區間可靠。
#
# 一個序列若提供許多相近的差值，
# 有效資訊量可能小於事件顆數。
# 正式分析可按事先定義的序列作整塊重抽樣。
# 整塊重抽樣（block bootstrap）會一起抽取相關觀測。
# 區塊太少時，它自身也很難估準區間。
#
# 固定參數的模擬，只讓地震目錄改變。
# 地震隨機性（aleatory variability）是這種事件波動。
# 參數不確定性（parameter uncertainty）來自有限資料估計。
# 若要納入後者，還需在學習資料上重新估參數。
# 本章區間只反映固定預報下的事件差值波動。
#
# ## 18.5 把基準換成 PPE
#
# SUP 在等面積格上分配相同空間率。
# 勝過它，表示模型利用了均勻基準沒有的資訊。
# PPE 已利用過去地震的位置。
# 把基準換成 PPE，便是在更具體的對照下比較。
# 同一張預報的物理內容完全沒有改變。

# %% tags=["remove-input"]
against_ppe = {m: csep.information_gain(forecasts[m], forecasts["PPE"], omega) for m in models if m != "PPE"}
for m, r in against_ppe.items():
    sup_gain = 0.0 if m == "SUP" else against_sup[m]["igpe"]
    assert np.isclose(r["igpe"], sup_gain-against_sup["PPE"]["igpe"])
fig = gain_figure(against_ppe, "PPE")
display(Markdown(f"圖：仍使用同一批 {int(omega.sum())} 顆事件。中心值可由兩個對 SUP 的 IGPE 相減得到；區間則重算配對差。"))
fig

# %% [markdown]
# 兩個對 SUP 的區間不能直接相減。
# 它們共用事件，也共用 SUP 分數。
# 重新計算 A 與 PPE 的逐事件差，
# 才能保留正確的共變異結構。
# 基準更換後的點估計容易算，區間仍要重做。
#
# 一個正 IGPE 也沒有指定實際用途。
# 實用性（utility）要連到行動成本與損失。
# 例如，強震後的巡檢需要特定時間窗。
# 十年總評分可能掩蓋那個窗的不足。
# 第 20 章會把預報接到地動與決策問題。
#
# ## 18.6 少量事件，能抓出多大的錯誤？
#
# 四個合成箱已足以示範辨識能力。
# 真實模型把較多事件放進第一箱。
# 受檢模型的配置則較平均。
# 統計功效（statistical power）是指定真相下的拒絕比例。
# 它取決於差異大小、樣本數與檢驗方法。
#
# 下圖逐次增加合成期望總數。
# 每次先生成一份目錄，再對受檢模型做 cL-test。
# 重複這個流程，就能估計拒絕比例。
# 兩模型的總期望數相同，差別只有分配。
# 這些合成結果與義大利模型排名無關。

# %% tags=["remove-input"]
n_outer, n_inner = 200, 150
synthetic_bins = 4
memory_bytes = max(n_outer, n_inner) * synthetic_bins * 8
assert memory_bytes < 100_000_000
sample_sizes = [10, 30, 90]
power = []
for i, n in enumerate(sample_sizes):
    true_rate = n*np.array([0.7, 0.1, 0.1, 0.1])
    tested_rate = n*np.array([0.4, 0.3, 0.2, 0.1])
    power.append(csep.power_by_simulation(true_rate, tested_rate,
        n_sim=n_outer, n_inner=n_inner, seed=1805+i, alpha=0.05))
fig = go.Figure(go.Scatter(x=sample_sizes, y=power, mode="lines+markers", line=dict(color=ACCENT)))
apply_layout(fig, title="合成示範：事件變多後的檢驗功效", xaxis_title="合成目錄的期望顆數", yaxis_title="拒絕受檢模型的比例", yaxis_range=[0, 1])
display(Markdown(f"圖：每點 {n_outer} 次外層模擬，每次 {n_inner} 份參考目錄。外層比例的標準誤上界約 {0.5/np.sqrt(n_outer):.3f}。"))
fig

# %% [markdown]
# 較高的拒絕比例表示較有機會抓到指定錯誤。
# 它不代表任何真實模型為假的機率。
# 若把兩種合成配置改得更接近，
# 同樣事件數通常更難分辨差異。
# 因此，功效報告必須一起列出假定真相。
#
# 每個功效點本身也來自有限模擬。
# 小幅起伏可能只是抽樣波動。
# 增加外層次數會減少比例估計的波動。
# 增加內層次數則細化檢驗的分位數。
# 兩層次數解決的問題不同。
#
# 一份沒有拒絕的檢驗結果，應連功效一起想。
# 少量事件可能讓很多模型都留下。
# 這與模型之間完全沒有差異是兩件事。
# 可讀的結論應同時留下點估計、區間與假設。
# 讀者才能判斷證據到底支持多強的比較。
#
# 一顆事件的極端差值，值得回到目錄核對。
# 先查兩張率是否使用同一個發報窗。
# 再查震央是否落在同一個空間格。
# 最後核對規模取整與箱邊界。
# 這些是比較成立前的資料條件。
# 率若對錯箱，精緻的區間也無法補救。
#
# 兩個模型都可能通過一致性檢驗。
# 它們仍可以具有不同的平均對數分數。
# 一致性問的是模型能否容納觀測。
# 比較問的是誰把相對較多資訊放對位置。
# 兩類結論並列，才不會把排名當成全面合格。
#
# 一張排名表也要留下模型版本。
# 本站 ETAS 只含背景與第一代貢獻。
# EEPAS 則採等權重的教學設定。
# 這些版本決定了表中的率與總數。
# 換成完整後代模擬後，應重新檢驗。
# 舊排名無法直接移植到新的實作。
#
# ## 18.7 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("共同觀測", "HORUS 測試期；沿用同一個 ω 與格箱邊界"),
    ("比較分數", "IGPE=(jPOLL_A−jPOLL_B)/N；保留總率項"),
    ("基準", "主比較 SUP；另列 PPE 敏感度對照"),
    ("區間", "近似 95% 配對 t；固定預報；未校正序列相依"),
    ("多模型比較", "原始區間；正式家族門檻須預先指定"),
    ("模型差異", "ETAS 為背景加第一代；EEPAS 採等權重"),
    ("功效", f"四箱合成示範；外層 {n_outer}、內層 {n_inner} 次"),
])))

# %% [markdown]
# 兩張圖的優勢可能出現在不同事件。
# 下一章：{doc}`19. 把模型加起來 <19_ensembles>`。
#
# ## 參考資料與延伸閱讀
#
# - Zechar（2010），[地震預報評估入門](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf)。
#   CORSSA 免費教材。先讀一致性，再讀模型比較。
# - Biondini、Rhoades 與 Gasperini（2023），[義大利 EEPAS 實驗](https://doi.org/10.1093/gji/ggad123)。
#   [免費機構全文](https://www.earth-prints.org/handle/2122/17084)。先對照研究區與測試期，再讀檢驗。
# - Khawaja 等（2023），[空間預報檢驗的統計功效](https://doi.org/10.1093/gji/ggad030)。
#   [免費機構全文](https://gfzpublic.gfz.de/pubman/item/item_5015770_1)。對照樣本量與網格的影響。
