# %% [markdown]
# # 18. 預報檢驗 II：比較、警報、功效與價值
#
# {doc}`第 17 章 <17_testing_consistency>`把「模型 vs 資料」這一路
# 走完了：Poisson 對數似然 POLL、二元似然 BILL、負二項的變異數
# 放寬，以及 N／M／S／cL 四個一致性檢驗與背後的模擬程序。那一章
# 的結論是一句克制的話——一致性檢驗是**診斷矩陣**，它告訴你模型
# 在哪個面向與資料不合，卻不告訴你該相信誰。
#
# 因為「與資料不矛盾」是一道太低的門檻。同一批目錄底下，可以有
# 十個模型全部通過 N-test、全部通過 S-test，而它們對下一次地震
# 該落在哪裡的看法南轅北轍。決策者要的不是「這個模型沒被拒絕」，
# 而是「這十個裡面該用哪一個」，甚至更進一步：「用了它，我到底
# 省下多少錢」。這一章因此分四段：18.1–18.3 把散落的資訊增益寫法
# 收攏成同一條式子，補上變異數估計式與參數懲罰；18.4 換一種語言，
# 把機率式預報二值化成警報；18.5–18.6 是自省——檢驗本身有多少
# 力氣、說 5% 的時候真的是 5% 嗎；18.7 走到最後一哩，問統計學不
# 回答的問題：這份預報**值多少錢**。會用到的舊式子只有第 10 章的
# $\ln L=\sum_i\ln\lambda^*-\int\lambda^*$ 與第 17 章的
# $\mathrm{POLL}=-\lambda+\omega\ln\lambda-\ln(\omega!)$；POLL、
# BILL、負二項與模擬程序一律引用不重推。
#
# ## 18.1 資訊增益的統一結構
#
# 前面幾章至少出現過三種「資訊增益」的寫法：第 10 章談連續時空
# 點過程時寫成兩個對數似然的差除以事件數；第 16 章比較 EEPAS 與
# PPE 用的是同一個量；而 CSEP 的比較檢驗寫出來卻多了一個看起來
# 很突兀的「期望總數差」修正項。它們是同一件事，三行代數就看得
# 清楚。
#
# ### 18.1.1 從對數似然差到 IGPE
#
# 連續版。第 10 章的時空規模點過程對數似然是
#
# $$\ln L = \sum_{i=1}^{N}\ln\lambda(t_i,x_i,y_i,m_i) - \Lambda,
# \qquad \Lambda = \int\!\!\int\!\!\int\lambda\,\mathrm{d}t\,
# \mathrm{d}A\,\mathrm{d}m$$
#
# 模型 $X$ 相對模型 $Y$ 的**每地震資訊增益**是兩者相減再除以觀測
# 事件數：
#
# $$\begin{aligned}
# I(X,Y) &= \frac{\ln L_X - \ln L_Y}{N} \\
# &= \frac{1}{N}\sum_{i=1}^{N}
#    \bigl[\ln\lambda_X(i) - \ln\lambda_Y(i)\bigr]
#    - \frac{\Lambda_X - \Lambda_Y}{N} .
# \end{aligned}$$
#
# 後項是「兩個模型對總量的看法差多少」，前項是「每顆地震平均賺到
# 多少對數機率」。單位是 **per earthquake**——分母是地震顆數，不是
# 格子數，也不是天數。這件事在 18.8 節會再犯一次案。
#
# 現在做網格化。CSEP 把測試區切成空間格 $j$ × 規模箱 $k$，模型交出
# 每箱的期望數 $\lambda(j,k)$、觀測數為 $\omega(j,k)$，第 17 章的
# 聯合對數似然是
# $\mathrm{jPOLL}=\sum_{j,k}[-\lambda+\omega\ln\lambda-\ln(\omega!)]$。
# 把候選模型 $Z$ 與基準模型 $1$ 的 jPOLL 相減，逐項看：
#
# $$\begin{aligned}
# \sum_{j,k}\lambda_Z(j,k) &= \hat N_Z, \qquad
# \sum_{j,k}\lambda_1(j,k) = \hat N_1, \\
# \sum_{j,k}\omega(j,k)\ln\lambda_Z(j,k)
# &= \sum_{n=1}^{N}\ln\lambda_Z(j_n,k_n) \equiv \sum_{n=1}^{N}X_Z(n), \\
# \sum_{j,k}\ln\bigl(\omega(j,k)!\bigr) &= \text{兩模型相同，相減後消掉}.
# \end{aligned}$$
#
# 第二行是樞紐：**把「對所有箱加總、以觀測數加權」改寫成「對所有
# 地震加總」**。一個裝了 $\omega$ 顆地震的箱在總和裡出現 $\omega$
# 次，正好等於那 $\omega$ 顆地震各報一次到。$(j_n,k_n)$ 是第 $n$ 顆
# 地震所在的箱，$X_Z(n)$ 是模型 $Z$ 在那個箱的對數率。於是
#
# $$
# \mathrm{IGPE} = \frac{\mathrm{jPOLL}_Z - \mathrm{jPOLL}_1}{N}
# = \frac{\hat N_1 - \hat N_Z}{N}
# + \frac{1}{N}\sum_{n=1}^{N}\bigl[X_Z(n) - X_1(n)\bigr]
# $$ (eq:igpe)
#
# 那個「突兀的修正項」原來就是連續版裡的
# $-(\Lambda_X-\Lambda_Y)/N$，只是把積分寫成求和，$\hat N$ 就是
# $\Lambda$。**三種寫法是同一條式子在不同離散化下的樣子。** 方向
# 要一次講定，因為文獻兩種都有：本書一律採 Rhoades 等人（2011）
# 與 Bayona 等人（2022）的方向——**下標 1 是基準模型、$Z$ 是候選
# 模型，IGPE 為正代表候選贏**；pyCSEP 文件的基準方向相反。
#
# 搭配二元似然 BILL 的版本叫 **IGPA**（information gain per active
# bin），形式一模一樣，只把分母的 $N$ 換成 active bin 數 $M$、
# $X$ 換成 BILL 分數。$M\le N$，且只有「多顆地震擠進同一箱」時
# 兩者才有實質差別——**IGPE 與 IGPA 的落差本身就是叢集強度的
# 診斷量**。
#
# ### 18.1.2 變異數估計式與 Student-$t$ 區間
#
# IGPE 是個平均值，平均值就有標準誤，有標準誤就能做 $t$ 檢定——
# 這才是「T-test」這個名字的由來。令 $d_n = X_Z(n) - X_1(n)$，
# 文獻裡的變異數估計式長這樣：
#
# $$s^{2} = \frac{1}{N-1}\sum_{n=1}^{N}d_n^{2}
# - \frac{1}{N^{2}-N}\Bigl(\sum_{n=1}^{N}d_n\Bigr)^{2}$$
#
# 看起來很怪，其實就是樣本變異數的展開式。設
# $\bar d=\frac1N\sum d_n$：
#
# $$\begin{aligned}
# \frac{1}{N-1}\sum_{n}(d_n-\bar d)^2
# &= \frac{1}{N-1}\Bigl[\sum_n d_n^2 - N\bar d^{\,2}\Bigr] \\
# &= \frac{1}{N-1}\sum_n d_n^2
#    - \frac{N}{N-1}\cdot\frac{1}{N^2}\Bigl(\sum_n d_n\Bigr)^2
#  = \frac{1}{N-1}\sum_n d_n^2
#    - \frac{1}{N^2-N}\Bigl(\sum_n d_n\Bigr)^2 .
# \end{aligned}$$
#
# 兩式相同，信賴區間因此是標準的單樣本 $t$ 區間
# $\mathrm{IGPE}\pm t_{N-1,\,1-\alpha/2}\,s/\sqrt{N}$。這裡藏著一個
# **必須誠實交代的近似**：修正項 $(\hat N_1-\hat N_Z)/N$ 被當成常數
# 平移，沒有納入變異數；區間只反映「每顆地震對數率差」的抽樣變異，
# 兩模型 $\hat N$ 差很多時會偏窄。判讀規則簡單到可以印在便利貼上：
#
# | IGPE 與區間 | 判讀 |
# |---|---|
# | 區間整個在零之上 | 候選模型顯著較具資訊量 |
# | 區間跨過零 | **統計上無法區分，就是平手** |
# | 區間整個在零之下 | 基準模型顯著較好 |
#
# 中間那一列最常被誤讀：點估計 $+0.3$ 但區間 $[-0.1,0.7]$，不能寫
# 成「本模型優於基準」。另有無母數的孿生兄弟 **W-test**（Wilcoxon
# signed-rank 版本，Rhoades et al. 2011），用來確認結論不依賴常態
# 假設；$d_n$ 通常右尾很長，T 與 W 不一致時該相信 W。下面用完全
# 合成的例子把整套算式跑一次：真實率場、基準模型與三個候選模型都
# 由亂數生成，事件按真實率場抽樣，再照 {eq}`eq:igpe` 與上面的
# $s^2$ 算出誤差棒。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly

setup_plotly()

# %% tags=["hide-input"]
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.special import gammaln

from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, apply_layout

GREEN = "#1baf7a"
GREY = "#8a8a8a"
rng = np.random.default_rng(20240403)

NY = NX = 40                       # 40 × 40 空間格，共 1600 箱
WATER = 1e-8                       # water level：避免 log(0)（見 18.8）


def norm_field(a):
    """把非負場正規化成機率場（總和為 1）。"""
    return a / a.sum()


# 真實率場：對白雜訊做高斯平滑再取指數 → 高度不均勻但空間連續
z = gaussian_filter(rng.standard_normal((NY, NX)), 2.0, mode="wrap")
truth = norm_field(np.exp(1.6 * z / z.std()))

# 四個模型：基準是模糊掉的真實場，三個候選各自朝不同方向偏離
uniform = np.full((NY, NX), 1.0 / (NY * NX))
p_bench = norm_field(gaussian_filter(truth, 2.6, mode="wrap"))   # 基準：偏模糊
p_a = norm_field(0.85 * truth + 0.15 * p_bench)                  # 甲：保留細節
p_b = norm_field(p_bench * np.exp(rng.normal(0, 0.45, (NY, NX))))  # 乙：基準加雜訊
p_c = norm_field(0.60 * uniform + 0.40 * p_bench)                # 丙：摻均勻成分

# 各模型的期望總數略有差異 → 修正項 (N̂₁ − N̂_Z)/N 才有作用
N_TRUE = 130
NHAT = {"基準": 128.0, "模型甲": 121.0, "模型乙": 134.0, "模型丙": 118.0}
PROB = {"基準": p_bench, "模型甲": p_a, "模型乙": p_b, "模型丙": p_c}

n_obs = rng.poisson(N_TRUE)
cells = rng.choice(NY * NX, size=n_obs, p=truth.ravel())


def igpe(name, bench="基準", conf=0.95):
    """回傳 (IGPE, 半寬)：完全照 18.1 的式子實作。"""
    lam_z = NHAT[name] * PROB[name].ravel()[cells] + WATER
    lam_1 = NHAT[bench] * PROB[bench].ravel()[cells] + WATER
    d = np.log(lam_z) - np.log(lam_1)                     # dₙ = X_Z − X₁
    n = d.size
    s2 = d @ d / (n - 1) - d.sum() ** 2 / (n**2 - n)      # 變異數估計式
    ig = (NHAT[bench] - NHAT[name]) / n + d.mean()
    half = stats.t.ppf(0.5 + conf / 2, n - 1) * np.sqrt(s2 / n)
    return ig, half


names = ["模型甲", "模型乙", "模型丙"]
tt = {m: igpe(m) for m in names}
colors = [GREEN, PALETTE[3], QUAKE_COLOR]

fig = go.Figure()
for m, c in zip(names, colors):
    g, e = tt[m]
    fig.add_trace(go.Scatter(
        x=[m], y=[g], error_y=dict(type="data", array=[e], width=10),
        mode="markers+text", marker=dict(size=13, color=c),
        text=[f"{g:+.2f} ± {e:.2f}"], textposition="middle right",
        showlegend=False))
fig.add_hline(y=0, line_dash="dash", line_color=GREY,
              annotation_text="基準模型")
apply_layout(fig,
             title=(f"T-test 判讀：{n_obs} 顆目標地震，"
                    "誤差棒跨過零線就是平手"),
             yaxis_title="每地震資訊增益 IGPE（相對基準）",
             xaxis_title="", hovermode="closest",
             yaxis_range=[-0.55, 0.55])
fig

# %% [markdown]
# 三種結局。**模型甲**保留了真實率場的大部分細節，區間整個落在零線
# 上方——顯著較具資訊量。**模型丙**摻了六成均勻成分，把資訊攤平掉
# 了，顯著較差。真正值得停下來看的是**模型乙**：它是「基準乘上一個
# 獨立的隨機擾動」，與基準明顯不同，但那些不同全是雜訊。它的點估計
# 不是零，誤差棒卻大剌剌地跨過零線——這時唯一能寫的結論是**「統計
# 上無法區分」**，點估計的正負號在跨零的區間裡不攜帶資訊。
#
# 另外注意誤差棒的**寬度**與模型甲的整個增益是同一量級。這不是巧合
# ——$N$ 只有一百多顆，$\sqrt{N}$ 才十出頭。**樣本數是地震顆數，
# 不是一千六百個格子。** 這句話在 18.5 節會變成主旋律。
#
# ## 18.2 Kagan 資訊分數：對照一個固定的無資訊參考
#
# IGPE 是**成對**比較，每個數字都綁著一個基準模型；十個模型兩兩比
# 就是四十五對。Kagan（2009）的**資訊分數 $I_1$** 走另一條路：所有
# 模型都對照同一個固定的無資訊參考，於是可以直接排序。
#
# $$I_1 = \frac{1}{n}\sum_{i=1}^{n}\log_2\frac{\lambda(x_i)}{E}$$
#
# $\lambda(x_i)$ 是第 $i$ 顆觀測地震所在箱的預報率，$E>0$ 是參考的
# **均勻 Poisson 過程**在同一箱的率。單位是 **bit**：$I_1=1$ 表示
# 模型平均替每顆地震把機率提高一倍；用 2 為底是為了讓「資訊」這個
# 詞回到 Shannon 的原意。
#
# 與 IGPE 有兩個實質差別。其一，**參考模型固定 vs 成對比較**：
# $I_1$ 可以把十個模型放在同一張長條圖上排名，IGPE 的數字則只在
# 「相對誰」講清楚時才有意義。其二，**沒有期望總數修正項**——
# $I_1$ 定義裡的 $\lambda(x_i)$ 不正規化到觀測事件數，系統性高估
# 總率的模型會白白賺到分數，而 IGPE 用 $(\hat N_1-\hat N_Z)/N$ 把
# 這條後門堵起來。
#
# $I_1$ 的兩個限制都很硬。**它完全不管沒有地震的格子**：求和只跑過
# $n$ 顆觀測地震，模型在其餘幾千格報了多少率一概不問。**零率箱則給
# $-\infty$**：一顆地震落在 $\lambda=0$ 的箱，整個分數崩潰——這個
# 陷阱對**所有含對數的分數**都成立，標準補救是設一個 water level
# 底線率（見 18.8）。所以 $I_1$ 適合粗排名與快速篩選，真正的勝負
# 判定回到 IGPE 加信賴區間。
#
# ## 18.3 懲罰：IGPEc 與 AICc
#
# 這一節解決回溯評估躲不掉的作弊漏洞。把基準模型 $\lambda_1$ 乘上一
# 個可調的乘數造出 hybrid：只要參數空間包含「乘數恆等於 1」這個點，
# 最佳化後的 hybrid 對數似然**在數學上不可能低於**基準。也就是說，
# 在擬合資料上多加參數的模型永遠贏，贏多少只是運氣。回溯期的 IGPE
# 因此不是技巧的度量，是**參數個數的度量**。
#
# 修正方式是把 AIC 家族的懲罰搬進來。CSEP 文獻用小樣本修正版 AICc
# （Hurvich & Tsai, 1989），因為地震檢驗的 $N$ 動輒只有幾十：
#
# $$\mathrm{AICc} = -2\ln L + 2n_p + \frac{n_p+1}{N-n_p-1}$$
#
# （文獻習慣把參數個數記作 $p$；本書的 $p$ 保留給 Omori 指數，改寫
# 成 $n_p$，指候選模型**比基準多出來**的自由參數個數。）
#
# ### 18.3.1 三項分解的完整推導
#
# 定義 $\Delta=\mathrm{AICc}_H-\mathrm{AICc}_1$，令
# $\mathrm{IGPEc}=-\Delta/(2N)$——除以 2 抵銷 AICc 定義裡的
# $-2\ln L$，除以 $N$ 讓單位回到「每顆地震」，與 IGPE 直接可比。
# 展開：
#
# $$\begin{aligned}
# \mathrm{IGPEc}
# &= \frac{-1}{2N}\Bigl[\bigl(-2\ln L_H + 2n_p
#    + \tfrac{n_p+1}{N-n_p-1}\bigr) - \bigl(-2\ln L_1\bigr)\Bigr] \\
# &= \frac{\ln L_H - \ln L_1}{N}
#    - \frac{1}{2N}\Bigl(2n_p + \frac{n_p+1}{N-n_p-1}\Bigr).
# \end{aligned}$$
#
# 第一項就是未懲罰的 IGPE，套用 18.1 的拆解：
#
# $$
# \mathrm{IGPEc} = \underbrace{\frac{\hat N_1 - \hat N_H}{N}}_{\text{總量校正}}
# - \underbrace{\frac{1}{2N}\Bigl(2n_p + \frac{n_p+1}{N-n_p-1}\Bigr)}_{\text{複雜度懲罰}}
# + \underbrace{\frac{1}{N}\sum_{n=1}^{N}\bigl[X_H(n)-X_1(n)\bigr]}_{\text{空間命中}}
# $$ (eq:igpec)
#
# **資訊增益 = 總量校正 − 複雜度懲罰 + 空間命中。** 總量校正罰
# 「總數預報不準」，與空間分布無關；複雜度懲罰罰「參數太多、資料
# 太少」，注意分母有 $N$，且 $\frac{n_p+1}{N-n_p-1}$ 在
# $N\to n_p+1$ 時發散——這是 AICc 內建的煞車，明白告訴你「參數個數
# 逼近樣本數時，別再談模型比較了」；第三項才是真正的「模型知道
# 地震在哪裡」。具體感覺：$N=40$、$n_p=3$ 時懲罰是
# $\frac{1}{80}(6+\frac{4}{36})\approx0.076$，而文獻中回溯期 hybrid
# 對基準的 IGPEc 也才 0.25 左右——懲罰吃掉三分之一；$N=10$ 時懲罰
# 跳到 0.14，超過一半。**小樣本的比較檢驗，懲罰項不是修飾，是主角。**
#
# ### 18.3.2 回溯必用 IGPEc，前瞻才可用 IGPE
#
# 這是本書反覆出現的那條「回溯／前瞻斷裂線」的形式化：
#
# | 評估情境 | 目標事件與擬合資料 | 該用 |
# |---|---|---|
# | 回溯（retrospective） | 重疊 | **IGPEc** |
# | 前瞻（prospective） | 獨立 | IGPE |
#
# 懲罰項存在的唯一理由，是模型在**同一批資料**上調過參數；前瞻期的
# 目標地震在模型凍結之後才發生，這時 IGPE 本身就是無偏的技巧度量。
# 實證後果有多嚴重，加州的十年開獎示範過（見
# {doc}`第 19 章 <19_ensembles>`會再談的乘法 hybrid 家族）：十六個
# hybrid 在回溯期對基準有 $+0.25$ 到 $+0.5$ 的 IGPEc——**而且是已經
# 扣過 AICc 懲罰的數字**——到了前瞻期沒有任何一個顯著優於基準，
# 好幾個掉到 $-0.4$ 至 $-0.7$。
#
# 教訓比「要罰參數」更深一層：**AICc 懲罰不足以救回過擬合**。懲罰項
# 假設「多出來的參數各自貢獻一份獨立的自由度」，但真實的過擬合是
# 模型結構整體對某一段時期的地震分布產生了依賴，那不是 $n_p$ 抓得住
# 的。所以**回溯評估永遠只是 sanity check**。
#
# ## 18.4 警報式語言：從列聯表到面積技能分數
#
# 前面三節都在似然的世界裡。但決策者問的不是「這格的率是多少」，
# 而是「今天要不要提高警戒」。這是個二值問題，而翻譯器只有一個：
# **門檻**。給定 $\lambda_{\rm th}$，把預報率高於門檻的格子標成
# 「發布警報」，於是每一格都落進 2 × 2 列聯表：
#
# | | 觀測有地震 | 觀測無地震 |
# |---|---|---|
# | **發布警報** | $a$（hit 命中） | $b$（false alarm 誤報） |
# | **未發警報** | $c$（miss 漏報） | $d$（correct negative 正確不報） |
#
# 關鍵在下一句：**逐步改變 $\lambda_{\rm th}$，就得到一整族列聯表**。
# 門檻從高掃到低，$(a,b,c,d)$ 畫出一條軌跡，所有警報式檢驗都只是這
# 條軌跡的不同座標系。挑單一門檻報一組命中率是不夠的，因為門檻是
# 使用者的自由變數（18.7 節會證明門檻該由誰決定）。
#
# ### 18.4.1 ROC 與它的隱含參考模型
#
# 最老牌的座標是 **ROC**：以命中率 $H=a/(a+c)$ 對誤報率
# $F=b/(b+d)$ 作圖，曲線愈往左上角 $(0,1)$ 靠愈好，對角線 $H=F$
# 代表隨機。它的弱點很要命：對角線這條參考線**隱含假設地震在空間上
# 均勻分布**——$F$ 的分母是「沒有地震的格子數」，地震稀疏時
# $b+d\approx$ 全部格子，於是 $F\approx$ 警報面積比例，ROC 等於在跟
# 「按面積隨機發警報」比。但地震本來就集中在斷層帶上，**只要把警報
# 都發在花東縱谷，就能輕鬆打敗這條對角線**，而這不需要任何預報技巧。
#
# ### 18.4.2 Molchan 圖：兩個端點與一條對角線
#
# Molchan 圖換了兩個座標：
#
# $$\tau = \frac{a+b}{a+b+c+d}, \qquad \nu = \frac{c}{a+c}$$
#
# $\tau$ 是**警報所佔的時空體積比例**，$\nu$ 是**漏報率**，畫 $\nu$
# 對 $\tau$。這個座標比 ROC 誠實，因為兩軸都直接對應決策者關心的
# 東西：我要付出多少警戒成本，會漏掉多少地震。兩個端點值得手算
# 一次。**門檻極高**時沒有格子發警報，$a=b=0$，於是
# $\tau=0/(c+d)=0$、$\nu=c/(0+c)=1$，得到左上角 $(0,1)$：不發警報、
# 漏掉全部。**門檻極低**時全區發警報，$c=d=0$，於是
# $\tau=(a+b)/(a+b)=1$、$\nu=0/(a+0)=0$，得到右下角 $(1,0)$：全面
# 警戒、一顆不漏。任何模型的軌跡都必須從 $(0,1)$ 走到 $(1,0)$，
# **差別只在中間怎麼走。**
#
# 對角線 $\nu=1-\tau$ 的意義：一個「按時空體積比例隨機發警報」的
# 策略，每顆地震落進警報區的機率就是 $\tau$，漏報率期望值因此是
# $1-\tau$。軌跡壓在對角線**左下方**代表比隨機好；越到右上方就是
# **比亂猜還差**——過度自信的局部預報真的會這樣。
#
# ### 18.4.3 面積技能分數
#
# 挑一個門檻報一個點仍是主觀選擇。Zechar 與 Jordan（2008, 2010）把
# 整條軌跡積分成一個數：
#
# $$\mathrm{AS}(\tau) = \frac{1}{\tau}\int_{0}^{\tau}
# \bigl[1 - \nu(t)\bigr]\,\mathrm{d}t$$
#
# $1-\nu(t)$ 是命中率，所以 $\mathrm{AS}(\tau)$ 就是「警報體積從 0
# 掃到 $\tau$ 這段區間內的**平均命中率**」；取 $\tau=1$ 時，
# $\mathrm{AS}(1)=\int_0^1[1-\nu]\,\mathrm{d}t$ 正好是 Molchan 軌跡
# **上方**的面積，值域 $[0,1]$。隨機參考模型的分數要算清楚，代入
# $\nu(t)=1-t$：
#
# $$\mathrm{AS}(\tau) = \frac{1}{\tau}\int_0^{\tau} t\,\mathrm{d}t
# = \frac{\tau}{2}$$
#
# 所以隨機模型在 $\tau=1$ 時 $\mathrm{AS}=0.5$——這就是「隨機參考為
# 0.5」的來源，但它**只在整條軌跡（$\tau=1$）成立**。若只積到
# $\tau=0.2$，隨機基準是 $0.1$ 而不是 $0.5$；報告部分面積分數時務必
# 附上對應的隨機基準。
#
# ### 18.4.4 Molchan 勝過 ROC 的兩件事
#
# 第一，**參考模型可以自選**。ROC 的對角線是釘死的（面積均勻）；
# Molchan 圖的 $\tau$ 可以用任意參考測度定義，最有用的選擇是拿長期
# 地震活動度當權重——把「警報體積」定義成警報區的長期期望地震數
# 佔全區的比例。這麼一改，「只在花東縱谷發警報」就不再佔便宜，因為
# 縱谷本來就貢獻了大部分長期期望數。**要打敗這條新對角線，模型必須
# 提供長期活動度以外的資訊**——那才是真正的預報技巧。
#
# 第二，**座標直接對應決策成本**。$\tau$ 是你付的錢、$\nu$ 是你漏掉
# 的災害；18.7 節會證明成本–損失模型的期望支出在 $(\tau,\nu)$ 平面
# 上是一條直線。下面用 18.1 那個合成率場實際跑一次，警報分數用模型
# 甲的預報率，對照組是純亂數分數。

# %% tags=["hide-input"]
counts = np.bincount(cells, minlength=NY * NX)          # 每格的事件數
n_cells = NY * NX


def molchan(score):
    """回傳 (tau, nu)：按分數由大到小放寬門檻，掃出整條軌跡。"""
    order = np.argsort(-score.ravel(), kind="stable")
    hits = np.cumsum(counts[order])
    tau = np.arange(1, n_cells + 1) / n_cells
    nu = 1.0 - hits / counts.sum()
    return np.r_[0.0, tau], np.r_[1.0, nu]               # 補上端點 (0,1)


def area_skill(tau, nu):
    """AS(1) = ∫₀¹ [1 − ν] dτ，梯形法（numpy 2 已移除 trapz）。"""
    h = 1.0 - nu
    return float(np.sum(0.5 * (h[1:] + h[:-1]) * np.diff(tau)))


tau_a, nu_a = molchan(p_a)                               # 有技能的模型甲
tau_r, nu_r = molchan(rng.random(n_cells))               # 純亂數分數
as_a, as_r = area_skill(tau_a, nu_a), area_skill(tau_r, nu_r)
i20 = np.searchsorted(tau_a, 0.20)                       # τ = 0.2 的命中率
hit20 = 1.0 - nu_a[i20]

fig = go.Figure()
fig.add_trace(go.Scatter(x=np.r_[tau_a, 1.0, 0.0], y=np.r_[nu_a, 1.0, 1.0],
                         fill="toself", fillcolor="rgba(42,120,214,0.16)",
                         line=dict(width=0), hoverinfo="skip",
                         name=f"面積技能分數 AS = {as_a:.2f}"))
fig.add_trace(go.Scatter(x=tau_a, y=nu_a, mode="lines", name="模型甲",
                         line=dict(color=ACCENT, width=2.8)))
fig.add_trace(go.Scatter(x=tau_r, y=nu_r, mode="lines", name="亂數分數",
                         line=dict(color=PALETTE[3], width=1.6)))
fig.add_trace(go.Scatter(x=[0, 1], y=[1, 0], mode="lines",
                         name="隨機參考線 ν = 1 − τ",
                         line=dict(color=GREY, dash="dash")))
fig.add_trace(go.Scatter(x=[0.20], y=[nu_a[i20]], mode="markers+text",
                         marker=dict(size=11, color=QUAKE_COLOR),
                         text=[f"  τ=0.20 時抓到 {hit20:.0%}"],
                         textposition="middle right", showlegend=False))
apply_layout(fig,
             title=(f"Molchan 圖：藍色面積即 AS = {as_a:.2f}"
                    f"（亂數分數 {as_r:.2f}，理論隨機值 0.50）"),
             xaxis_title="警報涵蓋的時空體積比例 τ",
             yaxis_title="漏報率 ν", hovermode="closest",
             xaxis_range=[0, 1], yaxis_range=[0, 1])
fig

# %% [markdown]
# 模型甲的軌跡壓在對角線左下方：只用兩成的警報面積就抓到圖上標註的
# 那個比例，藍色填色區就是軌跡上方的面積即 $\mathrm{AS}(1)$。亂數
# 分數則貼著對角線起伏，AS 在 0.5 附近抖動——抖動幅度本身就是
# 「單一次實驗的 AS 有多不穩」的提醒。兩個實作細節：軌跡必須手動
# 補上端點 $(0,1)$，否則積分會少掉最左邊那一段；計算命中數別寫成
# 雙層迴圈，用 `np.bincount` 算出每格事件數再對排序後的格子做
# `np.cumsum`，就是 $O(n\log n)$ 的一行——格數上萬時這是能不能跑
# 的差別。
#
# ## 18.5 檢驗的檢驗：統計功效
#
# 到目前為止的印象可能是：跑完 N/S/cL、再跑 T-test 和 Molchan，模型
# 就算被嚴格審查過了。這一節要把這個印象打碎。Khawaja 等人（2023）
# 做了一個實驗，結論可以印成海報：**一個宣稱「地球上每個地方地震
# 機率都一樣」的空間均勻模型，在標準的 0.1° × 0.1° 全球網格上通過
# 了 S-test**。不是勉強通過，是在他們測試的**所有樣本數下**都通過。
# 這個模型的空間資訊量精確地等於零，而 CSEP 最主力的空間一致性檢驗
# 抓不到它。
#
# ### 18.5.1 功效的定義與地震學的困境
#
# 抓不到，不代表檢驗寫錯了。統計檢定有兩種錯誤：型一誤（$H_0$ 為真
# 卻拒絕）與型二誤（$H_0$ 為假卻不拒絕）。顯著水準 $\alpha$ 控制
# 前者，**功效**控制後者：
#
# $$\text{Power} = P(\text{正確拒絕 } H_0) = 1 - P(\text{型二誤})$$
#
# 一致性檢驗只把 $\alpha$ 釘在 0.05，從來沒有人保證功效；而功效低到
# 0 的檢驗，「不拒絕」這個結果攜帶的資訊量精確地等於零。
#
# 地震學算功效有個結構性困難：**沒有已知的真實模型**。醫學試驗可以
# 說「假設藥效是 5 mmHg，需要多少受試者」，地震學連「真實的空間率場
# 長什麼樣」都不知道；CSEP 比較的是 equipollent hypotheses——一群
# 地位平等、都不知道對不對的假說。迂迴作法是**拿一個模型當資料產生
# 器**：選一個公認不錯的模型 $\Lambda_1$ 模擬出許多份「觀測」目錄，
# 再拿去檢驗另一個模型 $\Lambda_2$，
#
# $$\text{Power of S-test}
# = \frac{\text{S-test 失敗次數}}{\text{總模擬次數}}$$
#
# 這個定義有循環性（結論依賴你選了誰當「真實」），但它是目前唯一
# 可操作的路。
#
# ### 18.5.2 為什麼高解析度網格會殺死功效
#
# 0.1° 的全球網格有六百四十八萬格；六年的 $M$ 5.95 以上觀測大約
# 六百五十顆——**平均一萬格才攤到一顆地震**。格子細到「幾乎每顆
# 地震各佔一格」時，觀測目錄的聯合對數似然幾乎完全由「有幾格中了
# 1 顆」決定，而不是由「中在哪一格」決定：因為 S-test 把率正規化到
# $N_{\rm obs}$，每格的 $\lambda$ 都極小，$\omega=1$ 的箱貢獻
# $\ln\lambda$、$\omega=0$ 的箱貢獻 $-\lambda$，後者加總起來對所有
# 模型都一樣是 $-N_{\rm obs}$。空間資訊被解析度稀釋掉了。反過來說：
# **在高解析度網格上，均勻模型只有在觀測目錄含叢集（某些格子拿到
# 2 顆以上）時才可能被拒絕**——S-test 於是退化成一個叢集偵測器。
# 這推翻了一個很符合直覺的想法：「解析度愈高，檢驗愈嚴格」，
# **恰好相反**。
#
# ### 18.5.3 一維合成實驗
#
# 這個現象在一維上就能重現。真實事件從標準常態抽樣、限制在
# $[-3,3]$ 內；被檢驗的模型是**空間均勻模型**；把 $[-3,3]$ 切成
# $N_{\rm cell}$ 個箱跑一個 S-test 型的檢定。關鍵在**怎麼切箱**：
# **等寬分箱**把 $[-3,3]$ 均分，是標準 CSEP 網格的一維類比；
# **等期望率分箱**讓每箱涵蓋**真實分布的相同機率質量**（取真實分布
# 的分位數當箱界），是資料驅動多解析度網格的一維類比——事件多的
# 地方切細，少的地方切粗。檢定統計量用第 17 章的聯合 POLL，虛無
# 分布由「從均勻模型抽 $N$ 顆事件」模擬得到，觀測目錄則從常態抽；
# 功效 = 觀測分數低於虛無分布 5% 分位數的比例。

# %% tags=["hide-input"]
LO, HI = -3.0, 3.0
CDF_LO, CDF_HI = stats.norm.cdf(LO), stats.norm.cdf(HI)


def bin_probs(n_cell, kind):
    """回傳 (真實機率 p_true, 均勻模型的預報機率 q_fore)。"""
    if kind == "equal_width":
        edges = np.linspace(LO, HI, n_cell + 1)
    else:                                   # equal_rate：真實分布的分位數
        edges = stats.norm.ppf(
            np.linspace(CDF_LO, CDF_HI, n_cell + 1))
    p_true = np.diff(stats.norm.cdf(edges))
    q_fore = np.diff(edges)                 # 均勻模型：機率正比於箱寬
    return p_true / p_true.sum(), q_fore / q_fore.sum()


def joint_poll(omega, lam):
    """Σ_j [−λ_j + ω_j ln λ_j − ln ω_j!]，omega 形狀 (nsim, ncell)。"""
    return (omega * np.log(lam) - gammaln(omega + 1.0)).sum(1) - lam.sum()


def sim_ll(n_ev, prob, lam, nsim, gen, chunk=400):
    """分塊模擬，避免 (nsim × ncell) 的大陣列吃光記憶體。"""
    out = np.empty(nsim)
    for s in range(0, nsim, chunk):
        m = min(chunk, nsim - s)
        out[s:s + m] = joint_poll(
            gen.multinomial(n_ev, prob, size=m).astype(float), lam)
    return out


def power(n_ev, n_cell, kind, nsim=1500, alpha=0.05, seed=7):
    """S-test 對「空間均勻模型」的統計功效。"""
    gen = np.random.default_rng(seed)
    p_true, q_fore = bin_probs(n_cell, kind)
    lam = n_ev * q_fore + WATER             # 率正規化到 N_obs（S-test 慣例）
    crit = np.quantile(sim_ll(n_ev, q_fore, lam, nsim, gen), alpha)
    return float(np.mean(sim_ll(n_ev, p_true, lam, nsim, gen) < crit))


N_LIST = [4, 8, 16, 32, 64, 128, 256]
pw_w = [power(n, 20, "equal_width") for n in N_LIST]
pw_r = [power(n, 20, "equal_rate") for n in N_LIST]


def first_reach(xs, ys, target=0.9):
    for x, y in zip(xs, ys):
        if y >= target:
            return x
    return None


n90_w, n90_r = first_reach(N_LIST, pw_w), first_reach(N_LIST, pw_r)
tail = (f"等期望率分箱在 {n90_r} 顆達 0.9"
        + (f"，等寬分箱要 {n90_w} 顆" if n90_w else "，等寬分箱始終未達"))

fig = go.Figure()
fig.add_trace(go.Scatter(x=N_LIST, y=pw_r, mode="lines+markers",
                         name="等期望率分箱（資料驅動）",
                         line=dict(color=GREEN, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=N_LIST, y=pw_w, mode="lines+markers",
                         name="等寬分箱（標準網格）",
                         line=dict(color=QUAKE_COLOR, width=2.8),
                         marker=dict(size=8)))
fig.add_hline(y=0.9, line_dash="dot", line_color=GREY,
              annotation_text="功效 0.9")
apply_layout(fig,
             title=f"同一批資料、同一個檢驗，只換分箱方式：{tail}",
             xaxis_title="事件數 N（20 個箱）", yaxis_title="統計功效",
             xaxis_type="log", yaxis_range=[0, 1.02], hovermode="x")
fig

# %% [markdown]
# 兩條線的差距就是重點。**同一批事件、同一個檢驗、同一個被檢驗的爛
# 模型**——只因為箱界畫在不同地方，一條很快達到高功效，另一條慢吞吞
# 地爬。分箱不是無關緊要的技術細節，它決定了你的檢驗有沒有牙齒。
# 直覺是：等寬分箱時 $[-3,-2]$ 這種尾巴箱在真實分布下幾乎沒有事件，
# 在均勻模型下卻分到跟中央一樣的率，這些箱幾乎不貢獻資訊卻稀釋了
# 統計量；等期望率分箱把箱界排到事件真正會出現的地方，每個箱都在
# 承載資訊。
#
# 原始的全球實驗給出的門檻更誇張：在慣用的 0.1° 全球網格上，要讓
# S-test 達到有力的功效需要 **32,000 顆以上的目標地震，約當 300 年
# 的 $M \ge 5.95$ 記錄**；改用資料驅動的多解析度網格，**8 顆**地震
# 就能達到最大功效。（這兩個數字出自 Khawaja et al. 2023 的全球
# 實驗，不是上面這個一維玩具的產物；玩具只負責重現機制。）把樣本數
# 固定、改掃解析度，就能直接看到功效隨格數下降的曲線。

# %% tags=["hide-input"]
CELL_LIST = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
N_FIX = 64
pc_w = [power(N_FIX, c, "equal_width", nsim=1200) for c in CELL_LIST]
pc_r = [power(N_FIX, c, "equal_rate", nsim=1200) for c in CELL_LIST]
k_peak = int(np.argmax(pc_w))

fig = go.Figure()
fig.add_trace(go.Scatter(x=CELL_LIST, y=pc_r, mode="lines+markers",
                         name="等期望率分箱（資料驅動）",
                         line=dict(color=GREEN, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=CELL_LIST, y=pc_w, mode="lines+markers",
                         name="等寬分箱（標準網格）",
                         line=dict(color=QUAKE_COLOR, width=2.8),
                         marker=dict(size=8)))
fig.add_trace(go.Scatter(x=[CELL_LIST[k_peak]], y=[pc_w[k_peak]],
                         mode="markers+text", showlegend=False,
                         marker=dict(size=13, color=QUAKE_COLOR,
                                     symbol="circle-open",
                                     line=dict(width=3)),
                         text=[f"  峰值 {pc_w[k_peak]:.2f} "
                               f"@ {CELL_LIST[k_peak]} 箱"],
                         textposition="middle right"))
apply_layout(fig,
             title=(f"解析度愈高，檢驗愈弱：固定 {N_FIX} 顆事件，"
                    "等寬分箱的功效達峰後單調下滑"),
             xaxis_title="箱數 N_cell", yaxis_title="統計功效",
             xaxis_type="log", yaxis_range=[0, 1.02], hovermode="x")
fig

# %% [markdown]
# 紅線的形狀就是 18.5.2 的機制：格數很少時空間資訊被壓縮掉（極端
# 情況只剩 1 格，S-test 退化成 N-test，而且因為模擬事件數被固定成
# $N_{\rm obs}$，必然通過）；格數一多，資訊又被稀釋。中間有個甜蜜
# 點，位置取決於樣本數——**這意味著網格解析度應該是實驗設計的一
# 部分，不是慣例**。
#
# 綠線則幾乎不隨格數衰減，這就是 **Quadtree 多解析度網格**的一維版
# 原理。Quadtree 是樹狀鋪磚法：整個地球（Mercator 投影）當根節點，
# 每個 cell 要嘛不分、要嘛切成四個子 cell，由每格最多容納的資料點數
# $N_{\max}$ 與最大縮放層級 $L_{\max}$ 控制，結果是**地震多的地方
# 細、地震少的地方粗**。網格之間的轉換有理論保證：**聚合**（把小格
# 的率相加）之後獨立 Poisson 的和仍是 Poisson，率為
# $\Lambda=\sum_i\lambda_i$，所以真模型聚合後仍會通過檢驗；
# **反聚合**（把大格的率均勻攤到小格）則會注入模型沒有的假資訊。
#
# 三句必須背下來的話：（1）**「沒被拒絕」不等於「模型好」**，嚴謹的
# 寫法是「在此網格與此樣本數下，資料不足以拒絕該模型」並附上功效；
# （2）**低功效的檢驗只有拒絕時才有意義**，這與醫學、心理學界的
# underpowered study 危機是同一個病；（3）**設計比資料更能改善
# 檢驗**——觀測資料無法加速累積（三百年！），但網格是我們自己選的
# 自由變數。
#
# ## 18.6 校準：說 5% 的時候真的是 5% 嗎
#
# 資訊增益衡量的是**鑑別力**（discrimination）：模型能不能把高風險
# 的地方跟低風險的地方分開。但一份預報還有另一種好，叫**校準**
# （calibration／reliability）：模型說 5% 的那些場合，長期下來真的
# 有 5% 發生嗎？兩件事互相獨立——一個模型可以鑑別力很好卻嚴重失準
# （永遠把機率報成實際值的三倍，但排序完全正確），也可以完美校準卻
# 毫無鑑別力（永遠報氣候平均值 20%）。鑑別力決定「該把資源投到
# 哪裡」，校準決定「該投多少」。
#
# ### 18.6.1 可靠度圖怎麼建
#
# 四步：（1）收集一大批「預報機率 $P_i$、二元觀測 $y_i\in\{0,1\}$」
# 的配對，網格預報裡一個配對就是一個時空箱；（2）把 $P_i$ 分箱；
# （3）每箱算出預報機率平均 $\bar P_k$ 與實際發生比例 $\bar o_k$；
# （4）畫 $\bar o_k$ 對 $\bar P_k$。**完美校準的模型落在 45° 對角線
# 上**；曲線在對角線**下方**代表高估、**上方**代表低估；比對角線
# **平緩**（斜率小於 1）代表**過度自信**——敢報極端值卻兌現不了，
# 比對角線**陡**代表**過度保守**。一定要同時畫**樣本數**：機率 0.9
# 那一箱只有三個樣本時，那個點的位置沒有任何意義。
#
# ### 18.6.2 Brier 分數與它的分解
#
# 配套的單一分數是 **Brier 分數**：預報機率與二元觀測的均方差
#
# $$S_B = \frac{1}{n}\sum_{i=1}^{n}\bigl(P_i - y_i\bigr)^2$$
#
# 愈小愈好，完美預報為 0。相對於對數分數，Brier 有一個很實際的
# 優勢：**永遠有限**。地震預報在低活動區常出現機率為零的箱，對數類
# 分數會給出 $-\infty$ 的無限懲罰，Brier 最壞也只罰 1；所以建議兩者
# **併看**。Brier 分數可以精確拆成三項（Murphy, 1973）：把預報值
# 分成 $K$ 組，第 $k$ 組有 $n_k$ 個樣本、預報值 $P_k$、實際發生比例
# $\bar o_k$，全體發生比例為 $\bar o$，則
#
# $$S_B = \underbrace{\sum_{k}\frac{n_k}{n}\bigl(P_k-\bar o_k\bigr)^2}_{\mathrm{REL}}
# - \underbrace{\sum_{k}\frac{n_k}{n}\bigl(\bar o_k-\bar o\bigr)^2}_{\mathrm{RES}}
# + \underbrace{\bar o\,(1-\bar o)}_{\mathrm{UNC}}$$
#
# **REL（校準項）**是可靠度圖上各點偏離對角線的加權平方距離，愈小
# 愈好，只跟校準有關。**RES（鑑別項）**是各組發生比例偏離氣候平均值
# 的加權平方距離，愈大愈好（前面是減號）；永遠報氣候平均值的模型
# RES 恰為 0。**UNC（不確定項）**只取決於事件本身的基準率，**與模型
# 無關**——所以跨區域比較 Brier 分數之前要先確認基準率相近，否則
# 你比的是地質不是模型。推導在 18.11 節。
#
# ### 18.6.3 proper scoring rule
#
# **評分函數（scoring rule）** $S(P,y)$ 是把「預報 $P$」與「觀測
# $y$」映到一個數的函數。它是 **proper** 的，若且唯若當預報者相信
# 真實機率是 $Q$ 時，報出 $Q$ 本身能讓期望分數最佳；以負向（愈小
# 愈好）的分數寫成
#
# $$\mathbb{E}_{y\sim Q}\bigl[S(Q,y)\bigr]
# \le \mathbb{E}_{y\sim Q}\bigl[S(P,y)\bigr]
# \quad\text{對所有 } P$$
#
# 嚴格不等式（$P\neq Q$ 時必然更差）稱為 strictly proper。實際意義
# 是：**proper 的分數讓「說實話」成為最佳策略**。對數分數與 Brier
# 分數都是 strictly proper（證明見 18.11）。反例：線性分數
# $S=-\frac1n\sum P_i y_i$（「事件發生時報愈高分愈高」）是
# **improper** 的——不管你真正相信的機率是多少，把 $P$ 一律報成 1
# 都能拿到最好的期望分數。用 improper 的分數排名機率式地震預報會
# 產生**系統性偏誤**，通常是獎勵過度自信的模型。下面用三個合成模型
# 跑一次可靠度圖，並在程式裡**數值驗證** Murphy 分解是恆等式。

# %% tags=["hide-input"]
K_BIN = 20                                   # 預報值離散成 20 個格點
N_CASE = 40000
centers = (np.arange(K_BIN) + 0.5) / K_BIN
gen = np.random.default_rng(1999)


def snap(x):
    """把機率吸附到 20 個格點之一（讓 Murphy 分解成為精確恆等式）。"""
    return centers[np.clip((x * K_BIN).astype(int), 0, K_BIN - 1)]


truth_p = snap(gen.beta(1.5, 5.0, N_CASE))   # 各場合的真實機率
y = (gen.random(N_CASE) < truth_p).astype(float)
s_bar = y.mean()

MODELS = {
    "校準良好": truth_p,
    "過度自信": snap(np.clip(s_bar + 2.2 * (truth_p - s_bar), 0.002, 0.998)),
    "過度保守": snap(np.clip(s_bar + 0.40 * (truth_p - s_bar), 0.002, 0.998)),
}
MCOLOR = {"校準良好": GREEN, "過度自信": QUAKE_COLOR, "過度保守": PALETTE[3]}


def brier_parts(P, y):
    """回傳 (S_B, REL, RES, UNC)；P 已吸附到格點，故分解為恆等式。"""
    sb = float(np.mean((P - y) ** 2))
    rel = res = 0.0
    for v in np.unique(P):
        sel = P == v
        w, o = sel.mean(), y[sel].mean()
        rel += w * (v - o) ** 2
        res += w * (o - y.mean()) ** 2
    return sb, rel, res, float(y.mean() * (1 - y.mean()))


fig = go.Figure()
fig.add_trace(go.Scatter(x=[0, 0.8], y=[0, 0.8], mode="lines",
                         name="完美校準", line=dict(color=GREY, dash="dash")))
lines = []
for name, P in MODELS.items():
    sb, rel, res, unc = brier_parts(P, y)
    assert abs(sb - (rel - res + unc)) < 1e-10          # Murphy 恆等式
    xs, ys = [], []
    for v in np.unique(P):
        sel = P == v
        if sel.sum() >= 50:                              # 樣本太少的點不畫
            xs.append(v)
            ys.append(y[sel].mean())
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers", name=name,
                             line=dict(color=MCOLOR[name], width=2.6),
                             marker=dict(size=7)))
    lines.append(f"{name} S_B={sb:.4f}（REL={rel:.4f}, RES={res:.4f}）")
fig.add_hline(y=s_bar, line_dash="dot", line_color=GREY,
              annotation_text=f"氣候基準率 {s_bar:.2f}")
apply_layout(fig,
             title="可靠度圖：UNC=" + f"{s_bar * (1 - s_bar):.4f}；"
                   + "；".join(lines),
             xaxis_title="預報機率", yaxis_title="實際發生比例",
             hovermode="closest", xaxis_range=[0, 0.85],
             yaxis_range=[0, 0.85])
fig

# %% [markdown]
# 三條線各自示範一種病。**校準良好**的模型貼著對角線，REL 幾乎為零。
# **過度自信**的模型敢把機率往兩端推，曲線比對角線平緩。**過度保守**
# 的模型縮在氣候基準率附近，曲線比對角線陡：它幾乎不出錯（REL 小），
# 但也幾乎不提供資訊（RES 小）。三個分數與分解由程式帶入圖題，並在
# cell 裡用 `assert` 檢查了
# $S_B=\mathrm{REL}-\mathrm{RES}+\mathrm{UNC}$；UNC 三者相同，它只跟
# 資料的基準率有關。
#
# 注意過度保守的模型有可能拿到比過度自信更好的 $S_B$。這是 Brier 的
# 性格：它對「大膽而錯」的懲罰重於「膽小而平庸」。這就是為什麼
# **不要追求一個總分**——Brier、log score、資訊增益、Molchan 各自
# 懲罰不同的性質，併看才能定位模型在哪個面向失敗。
#
# ## 18.7 預報的價值：成本–損失模型
#
# 前面六節的分數全是統計量。但主管機關真正要問的是：**用了這份預報，
# 我到底省下多少？** 一個 IGPE 顯著為正的模型，如果在使用者的成本
# 結構下毫無用處，那個顯著性就是學術裝飾。
#
# 成本–損失模型是決策理論裡最簡單、也最常用來銜接預報與決策的框架。
# 設定：有 $n$ 個決策機會（例如 $n$ 個「格子 × 時段」），每個機會可
# 選擇**行動**（防護、加固、提高警戒）或**不行動**；行動成本是 $C$，
# **不管事件有沒有發生都要付**；不行動而事件發生，損失是 $L$（假設
# $L>C$，否則永遠不值得行動）；行動之後事件發生只付 $C$。沿用 18.4
# 的列聯表，行動 = 發布警報，期望總支出是
#
# $$E = C\,(a+b) + L\,c$$
#
# $a+b$ 是行動次數（不論事件是否發生都付 $C$），$c$ 是「沒行動卻出
# 事」的次數；$d$ 完全不花錢，$a$ 只花 $C$。
#
# ### 18.7.1 最佳行動門檻
#
# 單一機會的問題：已知這次事件發生的機率是 $P$，該不該行動？
#
# $$\begin{aligned}
# \mathbb{E}[\text{支出} \mid \text{行動}] &= C, \\
# \mathbb{E}[\text{支出} \mid \text{不行動}] &= L\,P .
# \end{aligned}$$
#
# 行動划算的條件是 $C < L\,P$，也就是
#
# $$P > p^{*} \equiv \frac{C}{L}$$
#
# **最佳行動門檻就是成本–損失比。** （$p^*$ 依文獻慣例寫成小寫；
# 本書其餘地方的機率一律寫 $P(\cdot)$，$p$ 保留給 Omori 指數。）
#
# 這個結論簡單得有點不像話，但有三層很重的意思。**其一，門檻不是
# 模型的性質，是使用者的性質**：防護成本低廉的使用者（例如發個簡訊
# 提醒）該在很低的機率就行動，成本高昂的（例如疏散一座城市）要等到
# 機率很高——**所以預報應該發布機率，不要替使用者決定門檻**。
# 其二，這正是 18.4 那個「掃過所有門檻」的動機：不同使用者站在
# Molchan 軌跡上的不同點。其三，**機率必須校準門檻才有意義**——
# 模型說 0.2 但實際是 0.05 時，$p^*=0.1$ 的使用者會被騙去行動。
#
# ### 18.7.2 價值分數
#
# 把三個策略的期望支出並列（都除以 $n$）：**氣候基準**不看預報，
# 只知道長期基準率 $s=(a+c)/n$，要嘛永遠行動（每次付 $C$）要嘛永遠
# 不行動（每次期望損失 $L\,s$），理性的使用者取小的那個，
# $E_{\rm climate}/n=\min(C,\ L\,s)$；**完美預報**只在事件真的會發生
# 時行動，$E_{\rm perfect}/n=C\,s$；**實際預報**則是
# $E_{\rm forecast}/n=C\,\tau+L\,s\,\nu$（推導見下）。**價值分數**把
# 實際預報放在兩個端點之間做線性內插：
#
# $$V = \frac{E_{\rm climate} - E_{\rm forecast}}
# {E_{\rm climate} - E_{\rm perfect}}$$
#
# 讀法：$V=1$ 表示與完美預報同等價值；$V=0$ 表示與「不看預報只看
# 氣候值」一樣；$V<0$ 表示**用了這份預報比不用還糟**。$V$ 是 $C/L$
# 的函數——**同一份預報對不同使用者有不同價值**，所以標準作法是畫
# 一條 $V$ 對 $C/L$ 的曲線，而不是報單一數字。
#
# ### 18.7.3 $V$ 與 Molchan 圖是同一件事
#
# 把支出式子除以 $nL$，並令 $\alpha=C/L=p^*$：
#
# $$\begin{aligned}
# e_{\rm forecast} &\equiv \frac{E_{\rm forecast}}{nL}
#  = \frac{C(a+b)}{nL} + \frac{Lc}{nL}
#  = \alpha\,\frac{a+b}{n} + \frac{c}{n} \\
# &= \alpha\,\tau + \underbrace{\frac{c}{a+c}}_{\nu}
#    \cdot\underbrace{\frac{a+c}{n}}_{s}
#  = \alpha\,\tau + s\,\nu .
# \end{aligned}$$
#
# 這一行就是全節重點：**期望支出在 Molchan 座標 $(\tau,\nu)$ 上是
# 一個仿射函數**。因此等支出線是斜率 $-\alpha/s$ 的**直線**；最小化
# 支出 = 在 Molchan 軌跡上找一條斜率 $-\alpha/s$ 的支撐直線的切點；
# 而 $e_{\rm climate}=\min(\alpha,s)$ 恰是右下角 $(1,0)$ 與左上角
# $(0,1)$ 兩個端點的支出，$e_{\rm perfect}=\alpha s$ 對應理想點
# $(\tau,\nu)=(s,0)$。所以 $V$ 與 Molchan 圖攜帶的資訊完全相同，
# 只是換了座標：**Molchan 圖畫的是整條抵換曲線，$V$–$C/L$ 曲線畫的
# 是「對每個成本結構的使用者，這條抵換曲線上最好的那一點值多少
# 錢」。** 前者是模型開發者的語言，後者是使用者的語言。完整代數放在
# 18.11 節。

# %% tags=["hide-input"]
ALPHAS = np.linspace(0.01, 0.99, 99)


def value_curve(P, y, alphas):
    """對每個 α = C/L，用最佳門檻 p* = α 行動，回傳價值分數 V。"""
    s = y.mean()
    out = np.empty_like(alphas)
    for i, a in enumerate(alphas):
        act = P > a                                   # 18.7.1 的最佳門檻
        e_f = a * act.mean() + np.mean((~act) & (y > 0))
        e_cl, e_pf = min(a, s), a * s
        out[i] = (e_cl - e_f) / (e_cl - e_pf)
    return out


fig = go.Figure()
peaks = []
for name, P in MODELS.items():
    v = value_curve(P, y, ALPHAS)
    j = int(np.argmax(v))
    peaks.append(f"{name} 峰值 {v[j]:.2f}／全段平均 {v.mean():+.2f}")
    fig.add_trace(go.Scatter(x=ALPHAS, y=v, mode="lines", name=name,
                             line=dict(color=MCOLOR[name], width=2.8)))
fig.add_hline(y=0, line_dash="dash", line_color=GREY,
              annotation_text="與氣候基準持平")
fig.add_vline(x=s_bar, line_dash="dot", line_color=GREY,
              annotation_text=f"C/L = 基準率 {s_bar:.2f}")
apply_layout(fig,
             title="成本–損失價值曲線：" + "；".join(peaks),
             xaxis_title="成本損失比 C/L（＝最佳行動門檻 p*）",
             yaxis_title="價值分數 V", yaxis_range=[-0.4, 0.55],
             hovermode="x")
fig

# %% [markdown]
# 三條曲線都在 $C/L$ 等於基準率附近達到峰值，這不是巧合：$C/L=s$ 時
# $e_{\rm climate}$ 的兩個選項（永遠行動、永遠不行動）打平，氣候基準
# 最無能為力，預報的相對價值也就最大。往兩端走 $V$ 都掉下來——
# $C/L$ 極小時「永遠行動」本來就夠好，$C/L$ 極大時「永遠不行動」
# 夠好。**價值曲線兩端趨近於零是這個框架的內建結構，不是模型的
# 缺點**；對「檢查逃生路線」這種 $C/L$ 極小的行動，任何預報都幫不上
# 忙（本來就該永遠做），對「疏散城市」也一樣。**預報真正能改變決策
# 的，是中間那一段。**
#
# 三個模型的**峰值幾乎一樣、全段平均卻天差地遠**，這才是這張圖真正
# 的教學點。峰值一樣是因為三者的機率互為單調變換，**排序完全相同**，
# $C/L$ 恰好等於基準率時錯誤的校準剛好把門檻推到對的地方。分手發生
# 在離開甜蜜點之後：**過度自信**的模型在低 $C/L$ 時報得太低而該行動
# 時不行動、在高 $C/L$ 時報得太高而不該行動時行動，曲線兩端都掉到
# 零線以下——**用了它比不用還糟**；**過度保守**的模型機率全擠在基準
# 率附近，門檻一離開那個區間它要嘛永遠行動、要嘛永遠不行動，於是
# **乾淨地退化成氣候基準**（$V=0$）。兩種失準，兩種後果：**保守只是
# 浪費資訊，過度自信會主動造成損害**——這也解釋了為什麼作業型系統
# 寧可保守。
#
# ## 18.8 實作陷阱
#
# **$\log$ 的 $-\infty$ 與 water level。** 任何含對數的分數，只要有
# **一顆**地震落在預報率為零的箱，整份分數就是 $-\infty$。標準補救是
# 給所有箱設底線率 $\lambda_{\min}$，但這個補救本身是要交代的選擇：
# 訂 $10^{-8}$ 與訂 $10^{-4}$，對同一顆落在零率箱的地震分數差了
# $\ln 10^{4}\approx9.2$，足以翻轉整場比較；加了 water level 之後
# **必須重新正規化**，否則 $\hat N$ 被墊高，IGPE 的總量校正項會出現
# 假訊號。最誠實的作法是**把 $\lambda_{\min}$ 當成模型宣告的一部分**，
# 在預報凍結時一起釘死。（Brier 分數在這裡完全免疫。）
#
# **對齊：規模尺度、深度、除叢。** 台灣目錄常見 $M_L$、全球模型用
# $M_w$，兩者不是同一個量（第 11 章），正確作法是先出示轉換式再把率
# 換算過去——一個全球模型從 $M_w$ 5.95 外推到區域的 $M_L$ 4.95，
# 中間至少要經過 GR 外推與尺度換算兩步。深度也一樣：全球模型可能算到
# 70 km、區域模型只到 30–40 km，標準作法是**實證檢查**深部事件貢獻
# 可忽略才宣稱可比。除叢與否更常被忽略：用完整目錄去檢驗一個以除叢
# 目錄擬合的模型，它必然低估總數而被 N-test 拒絕——那不是模型錯，
# 是檢驗設定錯。**檢驗的前提是 apples-to-apples。**
#
# **資料驅動的網格會繼承歷史盲點。** 網格由歷史目錄建構，於是**歷史
# 上安靜的地方會被切成粗格**，而地震最愛在歷史上安靜的地方發生。
# 2010 年紐西蘭 Canterbury 序列的起始位置，在依 2006 年以前資料建立
# 的 Quadtree 網格上恰好落在**低解析度區**。這不是否定多解析度網格，
# 而是提醒：**任何資料驅動的設計都會把歷史的偏誤帶進未來**；緩解
# 方式是把大地測量應變率、活動斷層等非目錄資訊也納入建格依據。
#
# **有效樣本數是地震顆數，不是格子數。** 一份台灣的網格預報可能有
# 一萬個空間格 × 二十個規模箱 = 二十萬箱，測試期五年只有 12 顆目標
# 地震。有人會覺得「二十萬個資料點」，於是對統計顯著性很有信心——
# **有效樣本數是 12。** IGPE 的分母是 $N$、標準誤的分母是 $\sqrt N$、
# Student-$t$ 的自由度是 $N-1$，通通都是地震顆數。後果在文獻裡看得
# 到：義大利測試區 8 年只有 11 顆目標地震，該區模型排名與更早研究
# 不同的主因僅僅是測試期有沒有涵蓋一個特定序列——**單一序列可以
# 翻轉整份排名。** 附帶推論：改用二元似然時有效樣本數變成 active
# bin 數 $M\le N$，放寬 Poisson 假設的代價就是樣本數。
#
# ## 18.9 常見誤解與陷阱
#
# **誤解 1：「模型通過了所有一致性檢驗，所以它是好模型。」** 最強的
# 反例是均勻模型在 0.1° 全球網格上通過 S-test。檢驗沒有拒絕，往往只
# 代表**檢驗沒有力氣拒絕**；正確的說法是「在此網格與此樣本數下，
# 資料不足以拒絕該模型」，且應同時報告統計功效。
#
# **誤解 2：「likelihood 分數低就是模型爛。」** 少數叢集事件會主宰
# 整個分數。加州的評估裡，兩個群震序列擠在極少數格子裡，Poisson
# S-test 對這幾格的重罰幾乎決定了所有模型的成敗——這正是二元似然被
# 提出的動機。診斷方法：畫出對數似然的**空間分布圖**，看是不是幾格
# 紅到發黑、其餘一片平淡。
#
# **誤解 3：「解析度愈高，檢驗愈嚴格。」** 恰好相反，18.5 節的紅線
# 就是證據。格子細到每顆地震各佔一格時空間資訊被稀釋，S-test 失去
# 鑑別力；慣用的全球高解析度網格上要達到有力檢驗需要約 32,000 顆
# 地震（三百年），改用資料驅動的多解析度網格 8 顆就夠。
#
# **誤解 4：「資訊增益 = 0.5，所以模型好。」** 永遠要問「相對於什麼
# 基準」。相對於 SUP（最小資訊模型）的增益，與相對於當時最佳模型的
# 增益，完全不是同一個量級的成就。同樣地，Kagan $I_1$ 的參考是均勻
# Poisson 過程，與成對比較的 IGPE **不可混用**——兩個模型的 $I_1$
# 相減也不等於它們的 IGPE，因為 $I_1$ 沒有期望總數修正項。
#
# **誤解 5：「同時跑很多檢驗比較保險。」** 跑愈多檢驗，至少一個假
# 陽性的機率愈高。要做多重檢定校正（如 Bonferroni $\alpha/T$），而
# $T$ 應該是**有效獨立檢驗數**：模擬證據顯示 S-test 與 cL-test 高度
# 相關（合理，cL 含 S），兩者都與 N-test 幾乎獨立（相關係數約
# 0.01–0.03），所以有效上 $T=2$ 而不是 4。
#
# **誤解 6：「模型愈在地、資料解析度愈高，一定愈準。」** 在加州、
# 紐西蘭、義大利這三個全世界儀器最密的地震區，一個以全球資料訓練的
# 模型（GEAR1）在前瞻期分別排到第 2、第 1、第 3 名。**地方資料的
# 優勢是需要被證明的假設，不是公理。**「我們的模型是為台灣量身打造
# 的，所以一定比較準」——這句話現在你知道該怎麼檢驗了：把全球基準
# 投影到台灣測試區，做 T-test。
#
# **誤解 7：「檢驗不通過就要淘汰模型。」** CSEP 自己的立場是**不因為
# 某個檢驗失敗就正式 reject 模型**，而是把分位數分數當診斷指標。檢驗
# 的目的是**改進模型**，不是頒發及格證書。
#
# **誤解 8：「有信賴區間就代表統計嚴謹。」** IGPE 的 $t$ 區間把總量
# 校正項當成常數，沒有納入 $\hat N$ 的不確定性；而且 $d_n$ 的分布
# 右尾很長，$N$ 只有幾十顆時常態近似本來就勉強——這是為什麼要搭配
# 無母數的 W-test。**區間是估計，不是保證。**
#
# ## 18.10 研究前沿與未解問題
#
# **CSEP 的第二階段與可重現性套件。** 第一代 CSEP 把檢驗程式封在
# 測試中心裡，外人無法檢視——這在「可重現性」這個核心訴求上是個
# 尷尬的矛盾。第二階段的方向是把裁判的尺也攤在陽光下：開源的評估
# 工具、公開的實驗規格、社群共同維護的基準模型。與之配套的是
# **可重現性套件**——程式碼、資料、凍結的執行環境（容器映像與版本
# 鎖定），一道指令重跑全部圖表。這不是附錄，是方法論的一部分。
#
# **pyCSEP 生態系。** 本章與上一章的工具——四種一致性檢驗、
# T/W-test、負二項與二元變體、Molchan 與面積技能分數、Quadtree
# 多解析度網格、把全球基準模型投影到任意區域——都已被實作在開源
# 套件裡。要理解的不是 API，而是那條資料流：**區域（region）→
# 預報（forecast）→ 目錄（catalog）→ 評估（evaluation）**。（本書
# 環境未安裝該套件，所有示範一律用合成資料自己實作。）
#
# **一致性評分函數的理論化。** proper scoring rule 的框架近年才被
# 系統性地引進地震預報評估。開放的問題包括：目錄式預報該用什麼
# proper score？如何設計一個既 proper、又對零率箱穩健、還能反映決策
# 價值的分數？18.7 節的價值分數不是 proper 的（它取決於使用者的
# $C/L$），而 proper 的分數不直接反映價值——兩者能不能調和仍是開放
# 問題。
#
# **低功效區域的檢驗設計。** 台灣測試區面積小、$M\ge4.95$ 的目標
# 事件少，正落在 18.5 節警告的低功效區間。出路有三條：多解析度
# 網格、降低目標規模門檻換取樣本數（但要處理 $M_c$ 與叢集）、跨區域
# 聯合檢驗（用空間換時間）。哪一條最有效，目前沒有定論。
#
# ## 18.11 附錄：本章推導細節
#
# ### A. Brier 分數的 Murphy 分解
#
# 把樣本按預報值分成 $K$ 組，第 $k$ 組的預報值為 $P_k$、樣本數
# $n_k$、發生比例 $\bar o_k$。先算組內：
#
# $$\begin{aligned}
# \frac{1}{n_k}\sum_{i\in k}(P_k - y_i)^2
# &= \frac{1}{n_k}\sum_{i\in k}\bigl[(P_k-\bar o_k)
#    + (\bar o_k - y_i)\bigr]^2 \\
# &= (P_k-\bar o_k)^2
#    + \frac{2(P_k-\bar o_k)}{n_k}\sum_{i\in k}(\bar o_k - y_i)
#    + \frac{1}{n_k}\sum_{i\in k}(\bar o_k - y_i)^2 .
# \end{aligned}$$
#
# 中間項為零（因為 $\sum_{i\in k}y_i=n_k\bar o_k$），末項因為
# $y_i\in\{0,1\}$ 而等於 $\bar o_k(1-\bar o_k)$。所以
# $S_B=\sum_k\frac{n_k}{n}(P_k-\bar o_k)^2
# +\sum_k\frac{n_k}{n}\bar o_k(1-\bar o_k)$，第一項已經是 REL。處理
# 第二項，令 $\bar o=\sum_k\frac{n_k}{n}\bar o_k$：
#
# $$\begin{aligned}
# \sum_k \frac{n_k}{n}\bar o_k(1-\bar o_k)
# &= \bar o - \sum_k \frac{n_k}{n}\bar o_k^2 \\
# &= \bar o - \bar o^2
#    + \Bigl(\bar o^2 - \sum_k \frac{n_k}{n}\bar o_k^2\Bigr)
#  = \bar o(1-\bar o)
#    - \sum_k \frac{n_k}{n}\bigl(\bar o_k - \bar o\bigr)^2 ,
# \end{aligned}$$
#
# 最後一步用了 $\sum_k\frac{n_k}{n}(\bar o_k-\bar o)^2
# =\sum_k\frac{n_k}{n}\bar o_k^2-\bar o^2$。合併即得
# $S_B=\mathrm{REL}-\mathrm{RES}+\mathrm{UNC}$。**這個分解只有在
# 「預報值真的取離散值」時才精確**；若 $P_i$ 連續，分組是人為的，
# 分解會隨分箱方式漂移——這也是 18.6 節的程式刻意把機率吸附到 20 個
# 格點的原因。
#
# ### B. Brier 分數是 strictly proper
#
# 設真實機率為 $Q$、預報者報 $P$，單次事件的期望 Brier 分數是
#
# $$g(P) = Q(P-1)^2 + (1-Q)P^2 = P^2 - 2QP + Q .$$
#
# 對 $P$ 微分：$g'(P)=2P-2Q$ 在 $P=Q$ 唯一為零，$g''(P)=2>0$ 故為
# 嚴格最小，所以報出真實機率是唯一的最佳策略。（對數分數
# $g(P)=-Q\ln P-(1-Q)\ln(1-P)$ 同理，$g'(P)=\frac{P-Q}{P(1-P)}$。）
# 反例：線性分數 $g(P)=-QP$ 對 $P$ 是遞減的線性函數，最佳策略恆為
# $P=1$，與 $Q$ 無關——improper。
#
# ### C. 成本–損失的完整代數
#
# 沿用 18.7 的記號，全部除以 $nL$ 化成無量綱形式，令 $\alpha=C/L$、
# $s=(a+c)/n$、$\tau=(a+b)/n$、$\nu=c/(a+c)$。實際預報的支出是
# $e_{\rm f}=\frac{C(a+b)+Lc}{nL}=\alpha\tau+\frac{c}{n}
# =\alpha\tau+s\nu$。氣候基準只能全報或全不報：
#
# $$\begin{aligned}
# \text{全報：}\ \tau=1,\ \nu=0 &\Rightarrow e = \alpha, \\
# \text{全不報：}\ \tau=0,\ \nu=1 &\Rightarrow e = s, \\
# e_{\rm cl} &= \min(\alpha,\ s).
# \end{aligned}$$
#
# 完美預報只在事件發生處報（$\tau=s$、$\nu=0$），
# $e_{\rm pf}=\alpha s$。於是
#
# $$V(\alpha) = \frac{\min(\alpha,s) - (\alpha\tau + s\nu)}
# {\min(\alpha,s) - \alpha s}$$
#
# 分母展開：$\alpha<s$ 時為 $\alpha(1-s)$，$\alpha>s$ 時為
# $s(1-\alpha)$；兩者在 $\alpha=s$ 相接，值為 $s(1-s)$，正是分母的
# 最大值——這解釋了為什麼價值曲線的峰值出現在 $\alpha\approx s$。
# 分子中的 $\alpha\tau+s\nu$ 對固定 $\alpha$ 是 $(\tau,\nu)$ 的
# **仿射函數**，等值線是斜率 $-\alpha/s$ 的直線；最大化 $V$ 等於
# 最小化 $e_{\rm f}$，也就是在 Molchan 軌跡上找斜率 $-\alpha/s$ 的
# 支撐線切點。當 $\alpha=s$ 時斜率為 $-1$，切點正是 Molchan 軌跡上
# 「離對角線最遠」的那一點——而那個最遠距離也正是 Molchan 圖最傳統
# 的技能度量。**兩套語言在這裡完全對上。**
#
# ---
#
# 這兩章合起來是一整套「怎麼判斷一個預報值不值得信」的方法，值得把
# 層次再排一次：一致性檢驗問「模型與資料矛盾嗎」，比較檢驗問「哪個
# 模型更好」，功效分析問「這個檢驗有沒有力氣回答前兩個問題」，校準問
# 「機率數字本身可不可以直接拿去用」，成本–損失則問「拿去用之後省下
# 多少」。**五個層次，五種不同的失敗方式**；只報一個總分的評估，一定
# 漏掉其中四種。
#
# 但這一章也留下一個尷尬的處境。18.1 的判讀圖裡，模型乙的誤差棒跨過
# 零線——統計上與基準平手。這在真實的 CSEP 實驗裡是**最常見的結局**：
# 一群模型互有勝負卻誰也沒有顯著贏過誰，而它們對「下次地震在哪裡」
# 的看法各不相同。決策者不能等三百年累積樣本，他明天就要一個數字。
# 那該用哪一個？{doc}`第 19 章 <19_ensembles>`給出的答案是：**也許
# 不必選。** 如果幾個模型各自捕捉到真實地震活動的不同面向，把它們
# 組合起來可能比任何單一模型都好。但那一章開頭就會回到本章的結論：
# 組合模型多出來的每一個權重都是一個參數，都要在 {eq}`eq:igpec` 的
# 第二項裡付出代價；而組合的天花板不在演算法，在候選池的**多樣性**。
# 檢驗的尺已經在手上了，接下來要量的是組合。
