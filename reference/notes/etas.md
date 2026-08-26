# ETAS 系列論文筆記
> 主題：ETAS 短期叢集模型 | 來源：reference/ 下 4 篇 | 供第二部 11 章（ETAS）等使用

**閱讀前的兩個檔名提醒（重要）**

1. `reference/[2017] ETAS_R.pdf` 的實際出版年是 **2019**（JSS Vol. 88, Code Snippet 1）。
2. `reference/[2023] Calibrated ETAS.pdf` 的內文**不是**任何「calibrated ETAS」論文，而是 Mizrahi, Nandan & Wiemer (2021) 談 declustering 對 b 值影響的 SRL 論文。全文中沒有出現 "calibrat" 這個字根。本筆記依實際內容撰寫，並以實際出版年排序。

---

## [2006] Ogata & Zhuang — Space–time ETAS models and an improved extension

- **書目**：Yosihiko Ogata, Jiancang Zhuang（The Institute of Statistical Mathematics, Tokyo）。*Tectonophysics* **413** (2006), 13–23。doi:10.1016/j.tecto.2005.10.016。原始檔：`[2006] Space–time ETAS models and an improved extension.txt`

- **這篇在做什麼**：先把 ETAS 從「修正大森公式」一路推到「時空 ETAS」的發展史整理成一條線，再指出 Ogata (1998) 最佳時空模型的一個系統性偏差——空間叢集尺度與 productivity 被綁在同一個參數 $\alpha$ 上。作者鬆開這個綁定，多加一個參數 $\gamma$，用日本氣象廳（JMA）1926–1995 三個構造背景不同的資料集做 AIC 比較，證明新模型普遍更好。動機是實務的：要偵測「地震寧靜（quiescence）」這類前兆異常，得先有一個能描述「正常叢集」的統計參考模型，否則寧靜只是餘震衰減的錯覺（Lomnitz & Nava 1983 的批評）。

- **關鍵觀念與公式**

  **(1) 修正大森公式（modified Omori formula）**——一切的起點：

  $$\nu(t)=K\,(t+c)^{-p}$$

  $t$ 為主震後經過時間（天），$K$ 為產出率尺度，$c$ 為避免 $t\to0$ 發散的時間常數（也吸收了早期紀錄不完整），$p$ 控制衰減快慢。Ogata (1983) 用**點過程最大概似法**估計，而不是先分 bin 再做迴歸：

  $$\ln L(\theta)=\sum_{i=1}^{N}\ln \nu(t_i)-\int_{S}^{T}\nu(t)\,dt,\qquad \theta=(K,c,p)$$

  第一項是「事件真的發生在 $t_i$」的獎勵，第二項是「這段期間預期發生幾個」的懲罰。這個「$\sum\log\lambda-\int\lambda$」骨架在整個 ETAS 家族反覆出現，值得單獨教一次。

  **(2) 時間 ETAS**——把每個事件都當成新的觸發源，做修正大森核的加權疊加：

  $$\lambda_\theta(t)=\mu+\sum_{\{j:\,t_j<t\}} e^{\alpha (M_j-M_c)}\,\nu(t-t_j)$$

  $\mu$（次/天）是背景率，$M_c$ 為資料完整度門檻。$e^{\alpha(M_j-M_c)}$ 就是 **productivity law**：一個規模 $M_j$ 的事件比門檻事件多產出 $e^{\alpha(M_j-M_c)}$ 倍的後代。$K$ 因此是「已用 $e^{\alpha(M-M_c)}$ 標準化後」的產出率。$\alpha$（單位：magnitude$^{-1}$）量測「規模轉換成觸發能力」的效率。

  **(3) 時空 ETAS** 的一般形式（本文式 4）：

  $$\lambda_\theta(t,x,y)=\mu(x,y)+\sum_{\{j:\,t_j<t\}}\nu(t-t_j)\,g\!\left(x-x_j,\;y-y_j;\;M_j-M_c\right)$$

  Ogata (1998) 比較了三種空間核 $g$：高斯型（式 5，指數收斂、邊界銳利）、以及兩種反冪次型（式 6、7）。AIC 一律選中式 (7)：

  $$g(\Delta x,\Delta y; M_j-M_c)\;\propto\;\left[\frac{(\Delta x,\Delta y)\,S_j\,(\Delta x,\Delta y)^{\mathsf T}}{e^{\alpha (M_j-M_c)}}+d\right]^{-q}$$

  $S_j$ 是 $2\times2$ 正定對稱矩陣，描述餘震區的橢圓形狀（斷層走向、傾角、定位誤差都揉在裡面）；$d$ 是空間平滑常數；$q$ 控制遠場冪次衰減。分母的 $e^{\alpha(M_j-M_c)}$ 就是 **Utsu–Seki 定律**的模型化：餘震區面積隨規模指數增長（$\log_{10}A = M+4.0$，或等價的 $\log_{10}L=0.5M-1.8$）。

  三個結論（原文 R1–R3）值得抄進講義：叢集區在空間上**超出傳統餘震區**、邊界模糊、冪次衰減（R1）；可能存在近場（破裂面周邊）與遠場（動態觸發、應力場遷移）兩個成分（R2）；叢集尺度與 Utsu–Seki 公式一致（R3）。

  **(4) 本文的貢獻：把空間尺度從 productivity 解耦。** 一般化寫法（式 9）把觸發核拆成「大小函數 × 正規化時間密度 × 正規化空間密度」：

  $$g(x,y,M)=\kappa(M)\cdot\frac{(p-1)c^{p-1}}{(t+c)^p}\cdot\frac{1}{\pi\sigma(M)}\,h\!\left(\frac{(x,y)S(x,y)^{\mathsf T}}{\sigma(M)}\right)$$

  Ogata (1998) 假設 $\kappa(M)\propto\sigma(M)\propto e^{\alpha M}$——**同一個 $\alpha$ 同時決定「產幾個」和「散多遠」**。Zhuang et al. (2004) 用 stochastic declustering 做診斷，發現「親代規模 vs 子代距離眾數」的斜率明顯比 $\hat\alpha$ 平緩，證實這個綁定造成偏差。本文因此提出（式 10）：

  $$g(\Delta x,\Delta y;M_j-M_c)=e^{(\alpha-\gamma)(M_j-M_c)}\left[\frac{(\Delta x,\Delta y)S_j(\Delta x,\Delta y)^{\mathsf T}}{e^{\gamma (M_j-M_c)}}+d\right]^{-q}$$

  參數由七個變八個：$\theta=(\mu,K,c,\alpha,\gamma,p,d,q)$。$\alpha$ 專管產量、$\gamma$ 專管空間尺度。用最小平方擬合診斷圖，Zhuang et al. (2004) 得到經驗值（式 11）：

  $$\tilde\gamma = 0.5\log_e 10 \approx 1.15$$

  也就是「餘震區長度 $\propto 10^{0.5M}$」的直接翻譯。診斷圖（Fig. 2c）在改用式 (10) 後系統性偏差幾乎消失。

  **(5) 模型選擇用 AIC**：$\mathrm{AIC}=-2\ln L(\hat\theta)+2\dim(\theta)$，越小越好。作者給的實用門檻：多 $d$ 個參數時，**AIC 差 2 約等於 5% 顯著水準**。

  **(6) Stochastic declustering（隨機去叢集）**——不做二分法，只給機率：

  $$\rho_{i,j}=\frac{\nu(t_j-t_i)\,g(x_j-x_i,y_j-y_i;M_i-M_c)}{\lambda(t_j,x_j,y_j)},\qquad \phi_j=\frac{\mu(t_j,x_j,y_j)}{\lambda(t_j,x_j,y_j)},\qquad \rho_j=1-\phi_j=\sum_{i<j}\rho_{i,j}$$

  $\rho_{i,j}$ 是「事件 $j$ 由事件 $i$ 觸發」的機率，$\phi_j$ 是「事件 $j$ 是背景事件」的機率。以 $\phi_j$ 為機率做 **thinning**（稀疏化），就得到一份背景子過程；換一個亂數種子就得到另一份。作者強調這本質上是 **bootstrap 重抽樣**，優點正是它把不確定性顯示出來，而傳統去叢集法把不確定性藏起來了。

- **教學上可用的洞見**
  - 「傳染病比喻」：ETAS 的 E 就是 epidemic。每個地震都是感染源，$\kappa(M)=Ke^{\alpha(M-M_c)}$ 是基本再生數，$\nu(t)$ 是潛伏期分布，$g(x,y)$ 是傳播距離分布。學生對 COVID 的 $R_0$ 有直覺，這個橋很好走。
  - **$\alpha$ 的經典數值範圍（Ogata 1992，日本）**：**群震（swarm）型活動 $\alpha \in [0.35, 0.85]$；非群震活動 $\alpha \in [1.2, 3.1]$**。$\alpha$ 小 = 大小地震觸發能力差不多 = 群震；$\alpha$ 大 = 大地震主導 = 典型主震–餘震序列。這是 ETAS 少數能直接對應「活動型態分類」的參數。
  - 本文擬合 JMA 三個資料集（式 10，含空間變動背景率）得到的量級：$\hat p \approx 1.03$–$1.05$、$\hat q \approx 1.58$–$1.74$、$\hat\alpha \approx 1.1$–$1.65$、$\hat\gamma \approx 0.80$–$1.33$。
  - **警語一**：在 Table 1（假設 $\mu(x,y)=$ 常數）中出現 $\hat p<1$，作者直接說這代表「均勻背景率的假設不成立」。$p<1$ 常常不是物理，是模型設定錯了——因為背景事件被硬塞進餘震核裡。
  - **警語二**（結論段）：除了 $\alpha, \gamma, p, q$ 之外，**其餘參數的估計值都會隨 $M_c$ 的選擇而變**（尺度差異）。跨研究比較參數前，先確認 $M_c$ 是否相同。
  - 大地震的震央位於餘震區「邊緣」是常態（震央 = 破裂起始點），所以模型中應把大事件的位置換成**餘震重心（centroid）**，否則空間核會被系統性拉偏。

- **與台灣的關聯**：文中未見台灣相關內容（資料全部來自日本 JMA）。教學橋接（本筆記標註，非原文）：台灣的隱沒帶（東部外海）與陸內（西部麓山帶）在構造上正對應本文的 Region A / Region B 對照，可拿同樣的分區邏輯講「為何 ETAS 參數要分區估計」。

---

## [2019] Jalilian — ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data

- **書目**：Abdollah Jalilian（Razi University）。*Journal of Statistical Software*, January 2019, **Vol. 88, Code Snippet 1**。doi:10.18637/jss.v088.c01。原始檔：`[2017] ETAS_R.txt`（檔名年份與實際出版年不符）

- **這篇在做什麼**：把 Zhuang–Ogata 的 Fortran 程式重寫成 R 套件 `ETAS`，用 Zhuang et al. (2002) 的 stochastic declustering 迭代法同時估背景率 $u(x,y)$ 與八個參數 $\theta$。全文其實是一份寫得很完整的「時空 ETAS 教科書 + 實作手冊」，是四篇裡**公式與演算法最完整**的一篇，最適合當講義主幹。

- **關鍵觀念與公式**

  **(1) 條件強度函數的完整分解**（式 1、2）：

  $$\lambda_{\beta,\theta}(t,x,y,m\mid H_t)=\nu_\beta(m)\,\lambda_\theta(t,x,y\mid H_t)$$

  $$\nu_\beta(m)=\beta\exp\!\left[-\beta(m-m_0)\right],\qquad \beta = b\ln 10$$

  $$\lambda_\theta(t,x,y\mid H_t)=\tilde u(x,y)+\sum_{i:\,t_i<t}\kappa_{A,\alpha}(m_i)\,g_{c,p}(t-t_i)\,f_{D,\gamma,q}(x-x_i,y-y_i;m_i)$$

  其中 $H_t=\{(t_i,x_i,y_i,m_i)\in X: t_i<t\}$ 是 $t$ 之前的全部歷史。**符號逐一說明**：

  - $\tilde u(x,y)=\mu\,u(x,y)$：背景率，時間上恆定、空間上非均勻。$u(x,y)$ 是無母數的平滑函數，$\mu$ 是加速收斂用的鬆弛係數。
  - $\kappa_{A,\alpha}(m)=A\exp[\alpha(m-m_0)]$：規模 $m$ 事件的**期望直接後代數**。
  - $g_{c,p}(t-t_i)=\frac{p-1}{c}\left(1+\frac{t-t_i}{c}\right)^{-p}$：後代發生時間的**機率密度**（已正規化，積分為 1，需 $p>1$）。注意這裡 Omori 核被寫成密度，所有「量」都被吸進 $\kappa$——這是比 Ogata 寫法更乾淨的參數化，$A$ 與 $K$ 不能直接互相比較。
  - $f_{D,\gamma,q}(\cdot)=\frac{q-1}{\pi D e^{\gamma(m_i-m_0)}}\left(1+\frac{(x-x_i)^2+(y-y_i)^2}{D e^{\gamma(m_i-m_0)}}\right)^{-q}$：後代位置的機率密度（各向同性，需 $q>1$）。$D$ 是空間尺度常數（度$^2$），$\gamma$ 就是 Ogata & Zhuang (2006) 解耦出來的餘震區尺度指數。
  - $m_0$：完整度門檻規模；$\beta$：GR 律的自然對數版斜率。

  **(2) Branching ratio 與臨界性**——整篇最該背下來的一條式子。對「隨機抽一個規模」的事件取期望：

  $$n=\int_{m_0}^{\infty}\kappa_{A,\alpha}(m)\,\nu_\beta(m)\,dm=\frac{A\beta}{\beta-\alpha}$$

  $n<1$（等價於 $A\beta/(\beta-\alpha)<1$，且需 $\alpha<\beta$）時模型是**次臨界（subcritical）**，序列必然熄滅、存在穩態版本，$T$ 期間的期望事件數正比於 $T$。$n\to1$ 是臨界點，$n>1$ 則分支過程爆炸、概似函數的積分項發散。$n$ 的物理意義：**目錄中有多大比例的事件是被觸發的**（$n=0.5$ 意味著約一半的地震是餘震）。

  Iran 目錄的實例（見下）：$n=\dfrac{0.1862\times5.6094}{5.6094-2.7071}\approx 0.36$。

  **(3) 叢集係數（clustering coefficient）**——把「哪裡容易叢集」畫成地圖：

  $$\Lambda(x,y)\approx \mu u(x,y)+\frac1T\sum_{i:t_i<T}\kappa_{A,\alpha}(m_i)f_{D,\gamma,q}(x-x_i,y-y_i;m_i),\qquad \omega(x,y)=1-\frac{u(x,y)}{\Lambda(x,y)}$$

  $\omega$ 接近 1 = 該處活動幾乎都是觸發事件；接近 0 = 幾乎都是背景事件。

  **(4) 估計：概似可分離 + 迭代去叢集。** 對數概似

  $$l(\beta,\theta\mid H_T)=\underbrace{l_1(\beta\mid H_T)}_{\text{規模}}+\underbrace{l_2(\theta\mid H_T)}_{\text{時空}},\qquad l_2(\theta)=\sum_i \delta_i\log\lambda_\theta(t_i,x_i,y_i\mid H_{t_i})-\int_{t_{\rm start}}^{t_{\rm start}+T}\!\!\iint_S \lambda_\theta\,dx\,dy\,dt$$

  $\delta_i=1$ 表示 $i$ 是 target event，$\delta_i=0$ 表示是 **complementary event**（研究期之前、或研究區之外的事件）——後者不進概似的第一項，但要進歷史 $H_t$，這是處理**邊界效應（edge effect）**的標準做法。因為 $\int\nu_\beta\,dm=1$，$\beta$ 與 $\theta$ 完全分離，$\beta$ 有封閉解：

  $$\hat\beta=\frac{N'}{\sum_i \delta_i (m_i-m_0)}$$

  即 Aki–Utsu 估計式。$\theta$ 則靠 **Davidon–Fletcher–Powell（DFP）**擬牛頓法最小化 $\xi(\theta)=-l_2(\theta)$，收斂時的 $H_k$ 就是反 Hessian 的估計，也就是漸近共變異矩陣（標準誤由此而來）。空間積分用 **radial partitioning**（把 $S$ 的邊界切成 $n_v$ 個節點、以徑向線段分割、轉極座標）近似。

  外層迴圈（Algorithm 2）交替更新：估 $\theta$ → 由式 (3)(4) 更新去叢集機率 $\hat p_j$ → 用**變頻寬高斯核**重估背景率

  $$\hat u(x,y)=\frac1T\sum_{j=1}^{N}(1-\hat p_j)\,\varphi(x-x_j,y-y_j;h_j),\qquad h_j=\max\{h_{\min},\,r(j,n_p)\}$$

  $r(j,n_p)$ 是第 $j$ 事件到第 $n_p$ 個最近鄰的距離（密的地方頻寬小、疏的地方頻寬大），$h_{\min}$ 防止重疊事件造成零頻寬，建議取定位誤差量級（預設 $h_{\min}=0.05^\circ\approx5.56$ km，$n_p=5$）。通常 10 次以內收斂。

  **(5) 殘差分析與 transformed time**——模型檢驗的核心工具。一階時空殘差

  $$R(I\times B;h)=\sum_i \delta_i \mathbf 1[t_i\in I,(x_i,y_i)\in B]\,h(t_i,x_i,y_i)\lambda_{\hat\theta}(t_i,x_i,y_i\mid H_{t_i})-\iiint_{I\times B} h\,\lambda_{\hat\theta}\,dx\,dy\,dt$$

  取 $h=1$ 得 raw residual、$h=1/\lambda$ 得 reciprocal residual、$h=1/\sqrt\lambda$ 得 Pearson residual。模型正確時殘差應在 0 附近無系統性偏離。更直觀的是**時間變換（transformed time）**：

  $$\tau_j=\int_{t_{\rm start}}^{t_j}\lambda^{\rm temp}_{\hat\theta}(t\mid H_t)\,dt,\qquad U_j=1-\exp\!\left[-(\tau_j-\tau_{j-1})\right]$$

  若模型正確，$\{\tau_j\}$ 是**單位速率 Poisson 過程**，$\{U_j\}$ 是 i.i.d. $U(0,1)$。診斷做法：畫 $\tau_j$ vs $j$（應落在 $y=x$ 直線上）、畫 $U_j$ 的 Q–Q 圖、跑 Kolmogorov–Smirnov 檢定。**$\tau$–$j$ 圖低於直線的區段就是「相對於模型的寧靜（quiescence）」，高於直線就是「相對活化」**——這正是 Ogata & Zhuang (2006) 一開始要解決的問題，兩篇在此接上。

- **教學上可用的洞見**
  - **可直接拿來上課的實例（Iran 目錄，1973–2016，$m_0=4.0$，5970 事件，ANSS）**，4 次迭代收斂、耗時 18.37 分鐘：

    | 參數 | $\beta$ | $\mu$ | $A$ | $c$ | $\alpha$ | $p$ | $D$ | $q$ | $\gamma$ |
    |---|---|---|---|---|---|---|---|---|---|
    | 估計 | 5.6094 | 0.5484 | 0.1862 | 0.0471 | 2.7071 | 1.1548 | 0.0160 | 2.3234 | 0.0238 |
    | 標準誤 | 0.0453 | 0.0133 | 0.0519 | 0.1093 | 0.0334 | 0.0106 | 0.1016 | 0.0361 | **5.7553** |

  - **參數相關性的最佳教材就在這張表**：$\hat\gamma=0.0238$ 而標準誤 **5.7553**，比估計值大兩個數量級——這個參數完全沒被資料約束。原因是 $\gamma$ 與 $D$ 在各向同性核裡高度反相關（$De^{\gamma(m-m_0)}$ 只有乘積被辨識）。**看到 ETAS 參數表，一定要一起看標準誤，否則會把數值噪音當成物理發現。**
  - **警語（原文明講）**：`etas()` 對初始值敏感，預設值（$\mu=N/(4T|S|)$、$A=0.01$、$c=0.01$、$\alpha=1$、$p=1.3$、$D=0.01$、$q=2$、$\gamma=1$）只是粗估，**不保證收斂**；且模型假設目錄**完整且時間穩態**，若有不完整或趨勢/季節性，結果不可靠。
  - **計算成本**：概似計算約正比於 $N^2$。目錄從 3000 筆變成 30000 筆，時間變 100 倍——這解釋了為什麼實務上要設 $m_0$，而不是「有多少資料用多少」。
  - 這一版套件**不含模擬功能**（作者在結論列為未來工作），模擬要另外實作或用別的套件。

- **與台灣的關聯**：文中未見台灣相關內容（範例資料為伊朗、義大利、日本）。教學橋接（本筆記標註）：`ETAS` 套件的 `catalog()` + `etas()` 流程可直接餵入 CWA/AutoBATS 目錄，是研究生做台灣 ETAS 練習最低門檻的入口。

---

## [2021] Mizrahi, Nandan & Wiemer — The Effect of Declustering on the Size Distribution of Mainshocks

- **書目**：Leila Mizrahi, Shyam Nandan, Stefan Wiemer（Swiss Seismological Service, ETH Zürich）。*Seismological Research Letters* **92**(4), 2333–2342, July 2021。doi:10.1785/0220200231。原始檔：`[2023] Calibrated ETAS.txt`（**檔名與內容不符**，見文首提醒）

- **這篇在做什麼**：問一個很尖銳的問題——大家在 PSHA 裡算的「主震 b 值」到底有沒有意義？作者拿加州 1980 年起的目錄，用五類去叢集法（Reasenberg、三種 window 法、Zaliapin 最近鄰、以及兩種 ETAS-based）掃過各種參數設定，再用 **ETAS 模擬 2000 份「b 值已知且主震餘震共用同一分布」的合成目錄**做對照。結論是：去叢集後 b 值下降的現象，在合成資料上同樣重現，因此**主要是演算法產物，不是主震的物理特性**。

- **關鍵觀念與公式**

  GR 律與「主震 GR 律」：

  $$\log_{10}N(m)=a-b\,m,\qquad \log_{10}N_{\rm main}(m)=a_{\rm main}-b_{\rm main}\,m$$

  當 $b_{\rm main}\neq b$，兩條線必在某個規模交叉：

  $$m_x=\frac{a-a_{\rm main}}{b-b_{\rm main}}$$

  以及主震占比 $r(m)=N_{\rm main}(m)/N(m)$。**邏輯陷阱**：$b_{\rm main}<b$ 意味著在 $m>m_x$ 之上，「預期主震數」會大於「預期總地震數」——這在觀測上是不可能的。$m_x$ 因此是一個可計算的**破產點**：低於 $m_x$ 危害度被低估（餘震被丟掉），高於 $m_x$ 危害度被高估（b 值被壓低、大地震機率被抬高）。

  ETAS 在本文有雙重角色：(a) 生成合成目錄的**生成模型**（參數由 EM 演算法從加州目錄反演，Veen & Schoenberg 2008；Nandan et al. 2017）；(b) 提供兩種去叢集定義：
  - **ETAS-main**：把叢集中最大事件定義為主震（與其他方法可比）。
  - **ETAS-background**：把「未被觸發的事件」定義為主震（符合 ETAS 本身的語意，任何規模的背景事件都可以引發串級）。

- **教學上可用的洞見**
  - **量級數字**：加州目錄去叢集後 b 值最多下降約 **30%**；不同方法得到的 $b$ 落在 **0.73–1.00** 之間，中間沒有明顯間隙（表示這是連續的方法學光譜，不是兩群）。不同方法留下的主震數相差 **6.1 倍**。
  - **關鍵對照實驗**：合成目錄的所有規模都抽自同一個分布（主震餘震無差別），去叢集後 b 值**照樣**下降。若 b 值下降是主震的物理性質，這在合成資料上不該發生。
  - **唯一的例外**：ETAS-background 去叢集後的 b 值與全目錄**無顯著差異**（合成與真實資料皆然）。差異只在套用「最大規模者為主震」這條規則時才出現。這個對照非常乾淨，可以當成一堂課的高潮：**問題不在「去叢集」，而在「主震 = 叢集中最大事件」這個定義**。ETAS 的餘震可以比它的親代更大（只被要求「發生在後」），這與傳統主震定義根本不相容。
  - **為什麼合成資料的效應比真實資料更強**：所有去叢集法都假設餘震空間分布各向同性——這對合成資料成立（生成時就是各向同性），對真實資料不成立。所以合成目錄的叢集比較好抓、小事件被移除得更徹底。這是一個漂亮的「模型誤設如何影響方法表現」的例子。
  - **給 PSHA 的三條指控**（結論）：主震定義不可驗證；忽略餘震會低估危害度（餘震一樣會致災）；被壓低的 b 值在 $m>m_x$ 造成高估。作者不接受「一低一高剛好抵銷」的辯護——那只在 $m=m_x$ 這一點成立。
  - **正面出路**：直接用 ETAS 做全目錄預報，不做去叢集。ETAS 只依賴**全目錄的 GR 律**，不需要對「被任意挑出來的大事件」假設規模分布；且用數十萬次情境模擬自然涵蓋時空叢集。Nandan et al. (2019) 的加州 pseudo-prospective 實驗中，ETAS 全面勝過 smoothed seismicity 與應變率模型。
  - 順帶一個統計常識：一組 i.i.d. 隨機變數的**最大值**的分布，本來就不會是原分布（Pareto/GR 的最大值不是 Pareto）。所以「主震規模服從 GR 律」在數學上本來就是不自洽的假設（Lombardi 2003；Zhuang & Ogata 2006 指出只有 $m\to\infty$ 的漸近情形才成立）。

- **與台灣的關聯**：文中未見台灣相關內容（資料為加州 ANSS ComCat）。教學橋接（本筆記標註）：台灣的機率式地震危害度評估同樣慣用 Gardner–Knopoff 視窗法去叢集，這篇的批評可直接搬過來當「為什麼台灣也該試 ETAS-based 全目錄評估」的論據。

---

## [2023] Mancini & Marzocchi — SimplETAS: A Benchmark Earthquake Forecasting Model Suitable for Operational Purposes and Seismic Hazard Analysis

- **書目**：Simone Mancini（Scuola Superiore Meridionale, Naples）, Warner Marzocchi（University of Naples Federico II）。*Seismological Research Letters* **95**(1), 38–49（cite-as 標 2023，刊於 January 2024 期）。doi:10.1785/0220230199。程式碼：https://github.com/smancini2/simplETAS 。原始檔：`[2023] SimpleETAS.txt`

- **這篇在做什麼**：反其道而行。當多數研究把 ETAS 越做越複雜（斷層幾何、深度分布、時空變動參數、貝氏即時更新），本文問：**能不能把大部分參數釘死，只留最少的自由度，還能用？** 答案是把描述叢集的七個參數 $\{\alpha,p,c,D,\gamma,q,\beta\}$ 依物理與經驗釘住，只估**背景率 $\nu$** 與**產出率 $A$** 兩個明顯與區域相關的參數。用義大利資料檢驗，從「日」到「世紀」四個時間尺度都通過。

- **關鍵觀念與公式**

  條件強度（式 1–5，符號略異於前兩篇）：

  $$\lambda(t,x,y,m)=\left[\nu\,\mu(x,y)+\sum_{i:\,t_i<t}\kappa(m_i)\,g(t-t_i)\,f(x-x_i,y-y_i;m_i)\right]s(m)$$

  $$\kappa(m_i)=A\,e^{\alpha(M_i-M_{\min})},\qquad g(t-t_i)=\frac{c^{p-1}(p-1)}{(t-t_i+c)^p}$$

  $$f(\cdot)=\frac{q-1}{\pi D e^{\gamma(M_i-M_{\min})}}\left[1+\frac{(x-x_i)^2+(y-y_i)^2}{De^{\gamma(M_i-M_{\min})}}\right]^{-q},\qquad s(m)=\beta e^{-\beta(m-M_{\min})}$$

  $\nu$ 是總背景率（時空恆定的純量），$\mu(x,y)$ 是背景事件的空間機率密度（本文直接取自義大利國家危害度模型 MPS19），$M_{\min}=M_c-\Delta M/2$（$\Delta M$ 為規模 binning）。標準參數化共九個相關參數。

  **釘住每個參數的理由**（這段是全篇教學價值最高的部分）：

  | 參數 | 值 | 理由 |
  |---|---|---|
  | $\alpha$ | $\ln 10$（$=\beta$，取 $b=1$） | 維持地震**自相似性**；可重現 **Båth 定律**；避免各向同性空間核造成的 productivity 偏差（Helmstetter 2005; Hainzl 2008）；且在同時考慮目錄不完整與時變背景率後，$\alpha$ 本來就很接近 $\beta$（Hainzl 2013） |
  | $p$ | 1.15 | 經驗範圍 $[0.9,1.4]$（Utsu 1995）的中位數；加州全州目錄的期望範圍 1.0–1.4 |
  | $c$ | 0.005 天 | 經驗範圍 0.003–0.3 天，但 $c$ 會被**早期餘震不完整**系統性高估、隨截切規模浮動；本文測試顯示 $c$ 對時間表現影響可忽略 |
  | $D$ | 1 km$^2$ | 全球估計橫跨數個數量級（加州 $<0.1$ 到隱沒帶 $>20$ km$^2$）；1 km$^2$ 取中間量級，也大致相當於現代目錄的平均水平定位誤差 |
  | $q$ | 1.5 | 與**靜態應力隨距離三次方衰減**一致 |
  | $\gamma$ | 1.5 | 資料豐富地區的估計典型落在 1.0–2.0 |
  | $\beta$ | $\ln 10$ | GR 律 $b=1$ |

  義大利校正結果：$\nu = 18.27$ 次/年、$A=0.047$（$M_{\min}=3.95$）。

  **檢驗設計（值得模仿的實驗架構）**：
  1. **50 年、$M\ge3.95$**（HORUS 儀器目錄，深度 $\le30$ km）：觀測累積曲線 vs 10,000 次模擬的分布。
  2. **392 年、$M\ge5.95$**（CPTI 歷史目錄，1630–2021）：**完全 out-of-sample**——參數只用近 50 年儀器資料估，卻要預測四個世紀的歷史地震數。這是很強的外推檢驗。
  3. **AVN（Amatrice–Visso–Norcia, 2016）序列的日預報**：pseudo-prospective，用 incremental N-test（Zechar et al. 2010 的 $\delta_1,\delta_2$ 分位數分數）逐日評分。
  4. 敏感度測試：擾動被釘住的參數，模擬目錄的事件數與位置變化有限。

  用的是 **Turing-style test**（Page & van der Elst 2018）：把合成目錄與真實目錄擺在一起，看能不能分辨。

- **教學上可用的洞見**
  - **一句話總結**：ETAS 的九個參數裡，真正「屬於這個地區」的可能只有兩個（背景率與產出率）；其餘七個描述的是**地震叢集的普世行為**。支持這個假設的獨立證據：Stallone & Marzocchi (2019) 發現日本、南加州、義大利的叢集特徵高度相似；Chu et al. (2011) 發現標準 ETAS 在全球不同構造區都擬合得不錯。
  - **三個用途**（作者自己列的，也是教學上最好的動機）：
    1. **Benchmark**：任何新預報模型都該先贏過 simplETAS；有了共同基準，跨區域的 CSEP 實驗才能公平比較「真正的技巧增益」。
    2. **OEF**：目錄短、資料少的地區（大多數國家）也能做作業型地震預報。
    3. **PSHA**：直接模擬任意長度的合成目錄，**繞開去叢集與 Poisson 假設**——與上一篇 Mizrahi et al. 的批評正好接上。
  - **教學價值最高的一點**：這篇示範了「不是所有參數都要估」。參數少 → 概似函數簡化 → 估計快且穩 → 可攜性高。對研究生而言，這是一個關於 **bias–variance trade-off** 的地震學版本：多兩個自由度換來的擬合改善，可能還不如它帶來的估計噪音。
  - **一個要在課堂上點破的數學細節（本筆記標註，原文正文未討論）**：當 $\alpha=\beta$ 時，branching ratio 的積分 $A\beta\int_{M_{\min}}^{\infty}e^{(\alpha-\beta)(m-M_{\min})}dm$ **發散**——$n=A\beta/(\beta-\alpha)$ 的公式只在 $\alpha<\beta$ 時成立。$\alpha=\beta$ 的參數化必須搭配**有限的最大規模 $M_{\max}$**，$n$ 才有限。這是 $\alpha=\beta$ 這類「自相似」設定的標準注意事項，教學時值得與上面的 branching ratio 一起講。
  - 排除了 Etna 火山區資料——火山構造的觸發機制不是標準 ETAS 能描述的。這是「知道模型的適用邊界」的好示範。

- **與台灣的關聯**：文中未見台灣相關內容（資料為義大利 HORUS / CPTI）。教學橋接（本筆記標註）：simplETAS 的邏輯對台灣特別有吸引力——若七個參數真的普世，台灣只需估 $\nu$ 與 $A$；同時台灣有龜山島等火山/地熱區，正好對應本文排除 Etna 的處理。

---

## 跨篇綜合：這個主題教什麼

### 一、符號對照表（先講，學生一定會在這裡卡住）

四篇的符號互不相同，第一堂課就要把對照表發下去，並宣告本課的**正規寫法**。

| 概念 | Ogata & Zhuang (2006) | Jalilian (2019) / 本課正規 | Mancini & Marzocchi (2023) |
|---|---|---|---|
| 門檻規模 | $M_c$ | $m_0$ | $M_{\min}=M_c-\Delta M/2$ |
| 產出率尺度 | $K$（Omori 核未正規化） | $A$（Omori 核已正規化為密度） | $A$ |
| productivity 指數 | $\alpha$ | $\alpha$ | $\alpha$ |
| Omori 參數 | $c, p$ | $c, p$ | $c, p$ |
| 空間尺度常數 | $d$ | $D$ | $D$ |
| 空間規模指數 | $\gamma$（2006 起才與 $\alpha$ 解耦） | $\gamma$ | $\gamma$ |
| 空間衰減指數 | $q$ | $q$ | $q$ |
| GR 斜率 | $b$ | $\beta=b\ln10$ | $\beta=b\ln10$ |
| 背景率 | $\mu(x,y)$ | $\tilde u(x,y)=\mu u(x,y)$ | $\nu\,\mu(x,y)$ |

**最容易出錯的一點**：$K$ 與 $A$ 不能直接比較。Ogata 的 $\nu(t)=K(t+c)^{-p}$ 不是機率密度；Jalilian 的 $g_{c,p}$ 是密度（積分為 1），所有「量」都收在 $A$ 裡。跨論文比較 productivity 前，先確認時間核有沒有正規化、$m_0$ 是否相同。

### 二、教學順序建議：從 Omori 到 ETAS 的七步推進

每一步只加一個新想法，學生跟得上：

1. **大森–宇津律**：$\nu(t)=K(t+c)^{-p}$。單一主震、單一序列。畫 1999 集集或 2024 花蓮的餘震率 log–log 圖，直線就在那裡。
2. **點過程的語言**：從「分 bin 數個數再迴歸」升級到 $\lambda(t)$ 與 $\sum\log\lambda-\int\lambda$。強調這一步不是為了炫技，而是因為分 bin 會浪費資訊、且時間解析度低時偏差嚴重。
3. **Gutenberg–Richter 與 productivity**：$s(m)=\beta e^{-\beta(m-m_0)}$ 與 $\kappa(m)=Ae^{\alpha(m-m_0)}$。問學生：$\alpha$ 與 $\beta$ 誰大？答案決定了「大地震主導」還是「小地震群主導」，也決定模型穩不穩定。
4. **疊加（superposition）= ETAS 的核心一步**：餘震也會生餘震。把單一 Omori 核換成 $\mu+\sum_j \kappa(m_j)g(t-t_j)$。此時「二次餘震」自動出現，不需要額外機制——這是 ETAS 最漂亮的地方。
5. **Branching ratio 與臨界性**：$n=A\beta/(\beta-\alpha)$。用「傳染病 $R_0$」比喻，$n<1$ 疫情熄滅、$n>1$ 爆炸。順帶點出 $\alpha=\beta$ 需要 $M_{\max}$ 截斷。
6. **加上空間**：Utsu–Seki 定律 → $De^{\gamma(m-m_0)}$ 的尺度 → 反冪次核 $[\,\cdot\,]^{-q}$。此處引入 Ogata & Zhuang (2006) 的教訓：**$\alpha$ 與 $\gamma$ 必須分開**，否則診斷圖會出現系統性斜率偏差。
7. **背景率與 declustering**：$\mu(x,y)$ 從哪來？→ stochastic declustering 的 $\rho_{ij}$、$\phi_j$ 與變頻寬核 → 迭代估計（Jalilian 的 Algorithm 2）。到這裡整個模型閉合。

補充兩堂：**檢驗**（transformed time、Q–Q、KS、N-test）與**應用**（OEF、PSHA、benchmark）。

### 三、「ETAS 是描述叢集的語言，不是物理定律」——支撐這個觀點的四條證據

這是整章的中心論點，四篇論文各提供一塊拼圖：

1. **它的組件全都是經驗律，不是從物理推導出來的。** Omori (1894)、Utsu–Seki (1955)、Gutenberg–Richter (1944) 都是先從資料看出來、再寫成公式的。Ogata & Zhuang (2006) 說得很直白：模型形式「based on empirical laws」。ETAS 是把三條經驗律組裝成一個自洽的隨機過程，不是應力轉移或摩擦定律的推論。
2. **同一份資料可以有多種函數形式，靠 AIC 挑，不靠物理挑。** Ogata (1998) 比了三個空間核（式 5/6/7），選中的理由是 AIC 最小，不是哪個比較「對」。2006 年再多加一個參數也是同一套邏輯。函數形式是可協商的。
3. **不同構造區的參數不同，卻沒有第一原理告訴你該是多少。** 群震 $\alpha\in[0.35,0.85]$ vs 非群震 $\alpha\in[1.2,3.1]$；$D$ 從加州的 $<0.1$ km$^2$ 到隱沒帶的 $>20$ km$^2$。反過來說，simplETAS 又發現把它們釘在中間值也堪用——這正是「描述性語言」的特徵：足夠有彈性去描述，但沒有唯一正確答案。
4. **它對「主震」這種人為概念完全無感。** Mizrahi et al. 的 ETAS 模擬裡，餘震可以比親代大，「主震」只是事後挑出來的標籤。ETAS 只承認「背景 vs 被觸發」，不承認「主震 vs 餘震」。當你把外部定義硬套上去，就會製造出 30% 的 b 值偏差這種假訊號。

補一句平衡：這不代表 ETAS 沒有物理內涵。$q=1.5$ 對應靜態應力隨距離三次方衰減、遠場成分對應動態觸發，都是有物理根據的解讀。但這些是**事後的物理詮釋**，不是模型的來源。

### 四、模擬演算法怎麼教最直觀

四篇裡沒有一篇寫出完整的模擬虛擬碼（Jalilian 的套件當時甚至還沒有模擬功能），但這是學生最需要動手做的部分。**建議用「分支（branching）法」教，不要用 thinning 法教**——分支法直接對應 ETAS 的物理故事，thinning 法只對應數學。

**分支法（cluster / branching simulation）——推薦教法**

1. **第 0 代（背景）**：在時間 $[0,T]$、空間 $S$ 上，依速率 $\nu\mu(x,y)$ 抽一個非齊次 Poisson 過程，得到背景事件的時空位置；每個事件的規模從 GR 律 $s(m)$ 獨立抽出。
2. **繁殖**：對「上一代」的每一個事件 $i$：
   - 後代數 $N_i \sim \text{Poisson}\!\left(A e^{\alpha(m_i-m_0)}\right)$；
   - 每個後代的**時間延遲**從 $g_{c,p}$ 抽（反函數法：$t-t_i=c[(1-U)^{-1/(p-1)}-1]$，$U\sim U(0,1)$）；
   - **空間位移**從 $f_{D,\gamma,q}$ 抽（極座標：角度均勻、半徑用反函數法）；
   - **規模**從 $s(m)$ 抽，**與親代規模無關**。
3. **遞迴**：把新產生的這一代當成上一代，重複步驟 2，直到某一代沒有新事件（$n<1$ 保證必然發生）或超出時間窗。
4. **收尾**：丟掉落在 $[0,T]\times S$ 之外的事件；為了避免邊界效應，起始時間要往前多跑一段 burn-in。

教學上的好處：學生寫完這 30 行程式碼，會**親眼看到二次餘震自己長出來**——不需要在模型裡寫任何「二次餘震」的規則。這比任何解釋都有說服力。也可以順手做兩個實驗：把 $n$ 從 0.3 調到 0.95，看序列長度如何爆炸；把 $\alpha$ 從 0.5 調到 2.0，看活動型態如何從「群震」變成「主震–餘震」。

**Thinning 法**則留給進階：先找到 $\lambda$ 的上界 $\lambda^*$，用速率 $\lambda^*$ 的齊次 Poisson 過程提案，以機率 $\lambda(t)/\lambda^*$ 接受。它的價值在於與 stochastic declustering 是同一個運算的正反兩面（模擬時「稀疏化」提案點；去叢集時以 $\phi_j$ 稀疏化真實點），這個對偶關係本身就是一個很好的收尾。

### 五、常見誤解清單（可做成課堂 quiz）

1. **「$p<1$ 代表餘震衰減特別慢。」** 多半是模型設定錯了。Ogata & Zhuang 明講：假設均勻背景率卻擬合出 $\hat p<1$，代表背景率的空間非均勻性被錯誤地吸收進 Omori 核。而且 $p\le1$ 時 $g_{c,p}$ 根本不可正規化。
2. **「$c$ 值是物理量。」** $c$ 主要反映**早期餘震目錄不完整**——大震後幾小時內小地震被大波形淹沒、測不到。$c$ 會隨截切規模與測站密度浮動。simplETAS 乾脆把它釘在 0.005 天，因為它對預報表現影響很小。
3. **「早期不完整只會影響 $c$。」** 不是。早期的高活動期被「咬掉」一塊，會系統性**壓低 $\alpha$ 的估計**（大地震剛發生後的產出被低估最嚴重）。Hainzl et al. (2013) 指出，同時處理不完整性與時變背景率後，$\alpha$ 會回升到接近 $\beta$——這正是 simplETAS 敢設 $\alpha=\beta$ 的依據之一。
4. **「參數估出來就是物理值。」** 看 Iran 的 $\hat\gamma=0.0238\pm5.7553$。ETAS 參數之間高度相關（$\gamma$–$D$ 反相關、$\alpha$–$A$ 相關、$c$–$p$ 相關），**一定要看標準誤，最好看聯合信賴區域**。單獨報一個參數的點估計幾乎沒有意義。
5. **「主震 b 值比餘震低，是主震的物理特性。」** Mizrahi et al. 用「b 值已知的合成目錄」證明這主要是去叢集演算法的產物。而且「一組 i.i.d. 變數的最大值」的分布本來就不是原分布，所以「主震服從 GR 律」在數學上不自洽。
6. **「ETAS 是預測地震的模型。」** ETAS 給的是**時間相依的機率**（未來 $\Delta t$ 內某格子發生 $M\ge m$ 的期望數），不是「哪天會地震」。它的技巧絕大部分來自「剛發生大地震之後」的短期窗口。
7. **「參數越多、模型越複雜就越好。」** simplETAS 是最好的反例：釘死七個參數，在義大利從日到世紀四個尺度都通過檢驗。複雜模型必須證明自己贏過這個 benchmark 才值得。
8. **「ETAS 的餘震一定比主震小。」** ETAS 對後代規模的唯一約束是「從 GR 律獨立抽出」，時間上必須在親代之後，規模上**沒有上限**。所以 ETAS 天然能產生「餘震大於主震」的序列——這正是它能自然重現前震現象的原因，也是它與傳統主震/餘震框架的根本分歧。

### 六、與台灣的連結（本筆記標註，非上述四篇內容）

四篇均未涉及台灣（資料分別來自日本、伊朗/義大利/日本、加州、義大利）。同一份 reference/ 收藏中，`Taiwan/` 與 `Taiwan__[2025] Fast_Report_Dapu_ETAS` 等檔案才是把這些方法落地到台灣的材料；教學上建議的接法是：先用本章四篇建立模型與方法的骨架，再用台灣案例（1999 集集、2024 花蓮 M7.2、2025 大埔）做參數估計與預報檢驗的實作。
