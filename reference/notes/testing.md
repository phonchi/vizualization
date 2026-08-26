# 預報檢驗／CSEP 論文筆記
> 主題：地震預報的評估與檢驗 | 來源：reference/ 下 6 篇 | 供第二部 15 章（CSEP 檢驗）等使用

---

## [2022] Bayona Tests — Prospective evaluation of multiplicative hybrid earthquake forecasting models in California

- **書目**：J. A. Bayona, W. H. Savran, D. A. Rhoades, M. J. Werner (2022), *Geophysical Journal International* 229(3), 1736–1753. DOI: 10.1093/gji/ggac018（Open Access）

- **這篇在做什麼**：
  RELM 實驗（2006–2010）曾判定 Helmstetter 等人的 adaptive smoothed seismicity 模型（代號 HKJ）是加州最具資訊量的時間獨立模型。Rhoades 等人（2014）以 HKJ 為 baseline、乘上其他模型當 conjugate，造出 16 個 multiplicative hybrid，回溯測試中不少 hybrid 顯著贏過 HKJ。本文把這 16 個 hybrid 加上 6 個原始 RELM 模型丟進**真正的前瞻期**（2011-01-01 至 2020-12-31，40 個 M ≥ 4.95 目標地震）重測，並同時引入以 **binary（Bernoulli）likelihood** 為基礎的新版 CSEP 檢驗，以降低檢驗結果對地震叢集的敏感度。結論很反直覺：回溯期看到的 information gain 在前瞻期**全部消失**，沒有任何 hybrid 顯著優於 HKJ。

- **關鍵觀念與公式**：

  **1. 網格化預報與 Poisson 假設**
  CSEP grid-based forecast 把測試區切成 0.1° × 0.1° 的空間格 $j$、再切 0.1 magnitude unit 的規模箱 $k$，模型交出每個 $(j,k)$ 箱在預報期間的**期望地震數** $\lambda(j,k)$。傳統檢驗假設各箱互相獨立、且箱內事件數服從 Poisson：

  $$P(\omega \mid \lambda) = \frac{\lambda^{\omega}}{\omega!}\,e^{-\lambda}$$

  其中 $\omega$ 為觀測事件數、$\lambda$ 為預報期望數。取對數得到單一箱的 log-likelihood（文中稱 POLL, Poisson log-likelihood）：

  $$\mathrm{POLL} = \ln P(\omega\mid\lambda) = -\lambda + \omega\ln\lambda - \ln(\omega!)$$

  **2. 一致性檢驗（consistency tests）：各自檢什麼**
  - **N-test（number）**：只看**總數**。$N_{\text{fore}} = \sum_{j,k}\lambda(j,k)$，問觀測總數 $N_{\text{obs}} = \sum_{j,k}\omega(j,k)$ 有沒有落在預報數分布的 95%（或 Bonferroni 校正後的 97.5%）區間內。檢的是「地震率抓得對不對」。
  - **S-test（spatial）**：只看**空間分布**。先把 $\lambda(j,k)$ 對規模箱 $k$ 加總得到 $\lambda(j)$，再把總和**正規化成 $N_{\text{obs}}$**（刻意抽掉率的資訊，只留形狀），算 joint log-likelihood 後與模擬分布比較。檢的是「震央的空間圖形對不對」。
  - **M-test（magnitude）**：對稱地，把 $\lambda$ 對空間加總只留規模維度、正規化到 $N_{\text{obs}}$。本文**沒有報 M-test**，理由是所有模型都用 Gutenberg–Richter 分配規模，結果必然雷同、不具鑑別力。
  - **cL-test（conditional likelihood）**：空間 × 規模的聯合分布，**不**做率的正規化調整（條件在觀測事件數上），是 S-test 與 M-test 的合體：

  $$\mathrm{jPOLL} = \sum_{j=1}^{l}\sum_{k=1}^{m}\Big[-\lambda(j,k) + \omega(j,k)\ln\lambda(j,k) - \ln\big(\omega(j,k)!\big)\Big]$$

  - **L-test**：原始的完整 likelihood 檢驗，本文**刻意不用**，因為它對「期望地震數」過度敏感（模型只要率抓錯就整個被拒絕，混淆了 N 與空間兩個面向的訊息）。這正是 cL-test 存在的理由——把「數量對不對」和「形狀對不對」拆開來看。

  **3. 模擬與分位數**
  觀測 log-likelihood 沒有解析分布，作法是**模擬**：把各箱的率正規化後做成累積分布，固定要模擬 $N_{\text{obs}}$ 個事件，抽 $U[0,1)$ 隨機數把事件丟進對應的箱，重複 10,000 次，得到一組模擬 catalog 及其 joint log-likelihood 分布。觀測值落在下尾即判定 inconsistent。分位數分數等價於單尾檢定的 p-value。

  **4. 比較檢驗 T-test 與 information gain per earthquake（IGPE）**
  一致性檢驗是「模型 vs. 資料」，比較檢驗是「模型 vs. 模型」。以 benchmark 模型 1（此處為 HKJ）對照模型 $Z$：

  $$\mathrm{IGPE} = \frac{\hat N_1 - \hat N_Z}{N} + \frac{1}{N}\sum_{n=1}^{N}\big[X_Z(n) - X_1(n)\big]$$

  其中 $N$ 是觀測目標地震數；$\hat N_1, \hat N_Z$ 是兩模型的期望總數；$X_Z(n) = \ln\lambda_Z(j_n,k_n)$ 是模型 $Z$ 在「第 $n$ 個地震所在那個箱」的 log-likelihood 分數。第一項修正兩模型率的差異，第二項是**每顆地震平均賺到多少 log-likelihood**。
  變異數估計：

  $$s^{2} = \frac{1}{N-1}\sum_{n=1}^{N}\big(X_Z(n)-X_1(n)\big)^{2} - \frac{1}{N^{2}-N}\Big(\sum_{n=1}^{N}\big(X_Z(n)-X_1(n)\big)\Big)^{2}$$

  信賴區間為 $\mathrm{IGPE} \pm t\,s/\sqrt{N}$（$t$ 為自由度 $N-1$ 的 Student t 分位數）。IGPE 為正且區間不含 0 → $Z$ 顯著較具資訊量；為負且不含 0 → benchmark 勝；含 0 → 兩者統計上無法區分。
  另有 **IGPEc**（corrected），用 AICc 精神加上參數懲罰項：

  $$\mathrm{IGPEc} = \frac{\hat N_1 - \hat N_H}{N} - \frac{1}{2N}\Big(2p + \frac{p+1}{N-p-1}\Big) + \frac{1}{N}\sum_{n=1}^{N}\big[X_H(n) - X_1(n)\big]$$

  $p$ 為擬合參數個數（hybrid 為 3）。回溯評估時**必須**用 IGPEc，因為模型在同一批資料上調過參數；前瞻評估時目標事件與擬合資料獨立，才能用未懲罰的 IGPE。

  **5. 放寬 Poisson：negative binomial N-test**
  Poisson 的變異數等於平均，對真實（未去餘震的）地震活動而言**過度低估變異**。負二項分布（NBD）多一個參數、變異數大於平均，較能容納時空叢集：

  $$p(\omega \mid \tau, \nu) = \frac{\Gamma(\tau+\omega)}{\Gamma(\tau)\,\omega!}\,\nu^{\tau}(1-\nu)^{\omega},\qquad \mu = \tau\frac{1-\nu}{\nu},\quad \sigma^{2} = \tau\frac{1-\nu}{\nu^{2}}$$

  實作訣竅：**平均取自模型的期望數，變異數由歷史觀測估計**——本文用 ANSS 目錄 1932–2010 的不重疊 10 年期估出全加州 $\sigma^2_C \approx 314.21$（南加州 SHEN、WARD 測試區分別為 185.70、164.83）。

  **6. 放寬 Poisson：binary（Bernoulli）likelihood**
  核心想法：不要問「這格出現幾顆」，只問「這格**有沒有**地震」。給定 $\lambda$，$p_0 = e^{-\lambda}$、$p_1 = 1 - e^{-\lambda}$，令 $X_i \in \{0,1\}$ 表示箱 $i$ 是否 active：

  $$\mathrm{BILL} = X_i \ln\!\big(1-e^{-\lambda}\big) + (1-X_i)\ln\!\big(e^{-\lambda}\big)$$

  $$\mathrm{jBILL} = \sum_{j=1}^{l}\sum_{k=1}^{m}\Big[X(j,k)\ln\big(1-e^{-\lambda(j,k)}\big) + \big(1-X(j,k)\big)\ln\big(e^{-\lambda(j,k)}\big)\Big]$$

  性質：$\omega = 0$ 時 POLL = BILL；$\omega = 1$ 且 $\lambda \to 0$ 時兩者近似；**只有在 $\omega \ge 2$ 的箱才會出現實質差異**——也就是叢集發生的地方。binary S-test 把率正規化到 **active cell 數**（而非 $N_{\text{obs}}$），模擬時也固定 active cell 數。
  對應的比較檢驗改成 **IGPA（information gain per active bin）**，把 $N$ 換成 active bin 數 $M$：

  $$\mathrm{IGPA} = \frac{\hat N_1 - \hat N_Z}{M} + \frac{1}{M}\sum_{n=1}^{M}\big[X_Z(n) - X_1(n)\big]$$

  **7. 多重檢定問題**
  同一個模型跑 N/S/cL 多個檢驗，型一誤（false-inconsistency）率會膨脹。本文用 Bonferroni 校正 $\alpha_B = 0.05/2 = 0.025$，並用 HKJ 當資料產生器模擬 1000 次、算各檢驗分位數分數之間的相關係數來決定「有效獨立檢驗數」：結果 **S-test 與 cL-test 高度相關，兩者都與 N-test 幾乎獨立**（$R_{NS} = 0.01$、$R_{NcL} = 0.03$），所以有效上只算 2 個獨立檢驗。

- **教學上可用的洞見**：
  - **回溯 ≠ 前瞻**：這篇最有教學價值的一句話。回溯期 hybrid 對 HKJ 有 0.2–0.5 的 IGPE，前瞻期全部歸零。作者的解釋是「形成 hybrid 的擬合本身在時間上不穩定」——換句話說，回溯期的 gain 大半是過擬合的幻覺。
  - **少數叢集主宰 likelihood**：2016 Hawthorne 群震（3 顆）與 2019 Ridgecrest 序列（8 顆）擠在極少數格子裡，Poisson S-test 對這幾格的懲罰極重，幾乎主導了整個檢驗結果。教學時可以畫 log-likelihood 的空間分布圖（本文 Fig. 8），讓學生直接看到「幾格紅到發黑，其餘一片平淡」。
  - **binary likelihood 的適用邊界**：本文 40 顆地震幾乎每顆各佔一格，因此 binary 與 Poisson 的 cL-test 結果**幾乎沒差**；差異只在 S-test（空間維度、有叢集）顯現。這是很好的反例教材——新方法不是萬靈丹，要在「多顆事件落同一箱」的情境（更大格子、更低規模門檻、更長預報期）才會發揮作用。
  - **CSEP 的檢驗哲學**：作者明講「我們不因為模型沒過某一檢驗就正式 reject 它」，而是把分位數分數當作**診斷工具**，指出模型與資料哪裡不合、值得科學上追究。這句話應該原封不動搬進教材。
  - 一個有趣的次要結論：把小地震（M ≥ 2）做高解析度平滑，仍是加州 5–15 年尺度預報中大地震最穩健的方法。複雜的組合模型沒有打敗它。

- **與台灣的關聯**：文中未提及台灣。

---

## [2023] Bayona Tests — Are Regionally Calibrated Seismicity Models More Informative than Global Models? Insights from California, New Zealand, and Italy

- **書目**：J. A. Bayona, W. H. Savran, P. Iturrieta, M. C. Gerstenberger, K. M. Graham, W. Marzocchi, D. Schorlemmer, M. J. Werner (2023), *The Seismic Record* 3(2), 86–95. DOI: 10.1785/0320230006（Open Access, CC-BY）

- **這篇在做什麼**：
  檢驗一個大家默認的假設：「區域模型用了解析度更高的區域資料，所以一定比全球模型更有資訊量」。作者把全球模型 **GEAR1**（smoothed seismicity × interseismic strain rate 的乘法組合，以 1977–2013 全球 CMT 的 $M_w$ 5.767+ 淺震校準）投影到加州、紐西蘭、義大利三個 CSEP 測試區，與 19 個區域時間獨立模型做前瞻對決（2014-01-01 至 2022-01-01，M 4.95+）。結果：GEAR1 在紐西蘭排第 1、加州第 2、義大利第 3——全球模型表現「意外地好」。

- **關鍵觀念與公式**：
  - **跨規模／跨深度的可比性處理**：全球模型要下放到區域測試，必須先對齊定義。作者用全球 b = 1 把 GEAR1 的 $M_w$ 5.95+ 率外推到 4.95+；義大利因為校準用 $M_w$、測試用 $M_L$，率要除以 1.602 做規模尺度換算；深度範圍 GEAR1 到 70 km、區域模型只到 30–40 km，作者實證檢查「深部事件貢獻可忽略」（紐西蘭 40–70 km 只佔 6%、義大利 30 km 以下掛零）才宣稱可比。**這一段是實作教材的黃金素材**：檢驗的前提是 apples-to-apples，規模尺度、深度範圍、去餘震與否若不對齊，後面的統計全是空談。
  - **主力工具是 T-test / IGPE**（同 Rhoades et al. 2011，公式見上一篇），加上 S-test 與 binary S-test。
  - **「每顆地震」與「每個 active cell」的正規化**：三區的格子數（加州 7682、紐西蘭 6343、義大利 8993）與目標事件數（38、47、11）差距懸殊，直接比 joint log-likelihood 沒有意義。作者改報**每顆地震的 log-likelihood 分數**（Poisson）與**每個 active cell 的分數**（binary），才能跨區比較。
  - **Gini coefficient $G_c$ 量化「平滑 vs. 集中」**：$G_c$ 定義為 ROC 曲線與對角線之間面積的兩倍，取值 $[0,1)$。ROC 曲線在此是「把格子按預報率由大到小排序後的累積率曲線」。$G_c \to 0$ 表示預報幾乎均勻攤平，$G_c \to 1$ 表示極度集中在少數格子。
  - **關鍵發現：平滑程度存在最佳區間**。$G_c \approx 0.6$–$0.75$ 的模型（HKJ、紐西蘭的 GEAR1、NZHM）拿到最高的 likelihood 分數；過度平滑（SUP、TRIPLE_S，低 $G_c$）系統性損失資訊，過度集中（PI、ALM，高 $G_c$）則是「過度自信的局部預報」，一旦地震沒在你指定的格子裡發生就被重罰。建議用 adaptive smoothing 這類「中間程度」的平滑。

- **教學上可用的洞見**：
  - **「模型愈在地愈準」是需要被檢驗的假設，不是公理。** 加州、紐西蘭、義大利是全世界儀器覆蓋最密、研究最透徹的地震區，區域模型卻沒能穩定打敗一個 1977 年起的全球資料集所訓練出的模型。這對「我們的模型是為台灣量身打造的所以一定比較好」這類主張是很好的警鐘。
  - **樣本量小到什麼程度會讓結論不穩**：義大利 8 年只有 11 顆目標地震，而本文與 Taroni et al.（2018）的排名之所以不同，主要就因為 2012 Emilia 序列（單一序列貢獻 11 顆 M 4.95+）落在或不落在測試期內。**單一序列可以翻轉整份排名**——這是「統計功效不足」最生動的例子。
  - 空間預報技巧可以拆解：作者用 residual analysis 逐格比較 GEAR1 與 NZHM 在 Kaikōura 地震周邊的表現，發現 NZHM（有斷層資訊）在主震破裂區贏，GEAR1（有大地測量應變率 + 較新的地震資料）在 170 km 外的餘震區贏。**逐格殘差分析比一個總分更能教會學生模型為什麼贏或輸。**
  - 建議語：GEAR1 可以當作區域模型的 global reference / benchmark，任何新的區域模型都應該先證明自己贏得過它。

- **與台灣的關聯**：文中未提及台灣。但方法論上這篇是「台灣模型該怎麼被公允評估」的直接範本——若要主張某個台灣區域模型有價值，正確作法是把 GEAR1 投影到台灣測試區當 benchmark 做 T-test，而不是只報自家模型的絕對分數。

---

## [2023] Khawaja Stat_power_test — Statistical power of spatial earthquake forecast tests

- **書目**：A. M. Khawaja, S. Hainzl, D. Schorlemmer, P. Iturrieta, J. A. Bayona, W. H. Savran, M. Werner, W. Marzocchi (2023), *Geophysical Journal International* 233(3), 2053–2066. DOI: 10.1093/gji/ggad030

- **這篇在做什麼**：
  問一個尖銳的問題：CSEP 的 S-test 到底有沒有能力分辨好壞模型？作者做了一個全球預報實驗，測試 S-test 拒絕「空間完全無資訊的均勻預報模型」的統計功效。結果令人不安：在慣用的 0.1° × 0.1° 全球網格上，**均勻模型竟然通過 S-test**。作者釐清了兩個決定功效的因子（樣本數、網格解析度），並提出用 Quadtree 多解析度網格作為解方。

- **關鍵觀念與公式**：

  **1. 統計功效（statistical power）的定義與地震學的困境**

  $$\text{Power} = \Pr(\text{正確拒絕 } H_0)$$

  問題在於：地震學裡**沒有已知的真實模型**，CSEP 的 likelihood 檢驗比較的是 equipollent hypotheses。作者的迂迴作法是——**拿一個模型當資料產生器**：以 GEAR1 為 $\Lambda_1$ 模擬「觀測」目錄，再用這些目錄去檢驗另一個模型 $\Lambda_2$，重複多次，

  $$\text{Power of S-test} = \frac{\text{S-test 失敗次數}}{\text{總模擬次數}}$$

  **2. 為什麼高解析度網格會殺死功效（直覺）**
  0.1° × 0.1° 的全球網格有 648 萬格，而 6 年的 M 5.95+ 觀測只有 651 顆——平均**一萬格才一顆地震**。格子小到「幾乎每顆地震各佔一格」時，觀測 catalog 與模擬 catalog 的 log-likelihood 都由「有幾格中了 1 顆」決定，而不是由「中在哪裡」決定。空間資訊被解析度稀釋掉了，S-test 遂失去鑑別力。反過來說，均勻預報在高解析度網格上**只有在觀測目錄含叢集（某些格子拿到 2 顆以上）時才可能被拒絕**——這等於把空間檢驗變成了叢集偵測器。

  **3. 樣本數的量化門檻**
  以 GEAR1 為資料產生器，S-test 拒絕各競爭模型所需的事件數差異很大（表 1）：TEAM、SHIFT_GSRM2f 在 64 顆就被拒絕（power 0.75–0.98），KJSS 與 WHEEL 要 1024–2048 顆才達 power ≈ 1（因為 WHEEL 與 GEAR1 共用 KJSS 元件，兩者太像）；均勻模型在**所有**樣本數下 power 都是 0。作者估計：在慣用的全球高解析度網格上，要達到有力的檢驗需要 **32,000 顆以上地震，約當 300 年的 M ≥ 5.95 記錄**。

  **4. 一維合成實驗（最好的教學示意）**
  在 $[-3, 3]$ 上以標準常態抽事件、用均勻模型去預報，比較兩種分箱：
  - **uniform bins（等寬）**：功效在 $N_{\text{cell}} = 4$ 達到最高（僅約 0.3），之後隨格數單調下降。
  - **density-based bins（等期望率）**：功效隨格數收斂到約 0.85。
  固定 $N_{\text{cell}} = 20$ 時，density-based 只要 **15 顆**事件就達到 power 0.9，等寬分箱要 **38 顆**。同樣一批資料，只是換了分箱方式，檢驗力就差了一倍以上。

  **5. Quadtree 多解析度網格**
  Quadtree 是樹狀鋪磚法：整個地球（Mercator 投影）當根節點，每個 cell 要嘛不分、要嘛切成四個子 cell。以兩個參數控制——每格最多容納的資料點數 $N_{\max}$、最大縮放層級 $L_{\max}$——就能讓**地震多的地方細、地震少的地方粗**。在多解析度網格上，S-test 只要 **8 顆地震**就能達到最大功效。
  預報在網格間的轉換：聚合（aggregation）就是把小格的率相加；反聚合（de-aggregation）則把大格的率均勻攤到小格。理論保證：Poisson 相加仍是 Poisson（率為 $\Lambda = \sum\lambda_i$），所以**真模型聚合後仍會通過檢驗**（作者實測通過率 95.6–100%，符合理論的 97.5%），但**假模型會在新網格上被抓出來**，因為功效提高了。

  **6. 重新評估的結果**
  把全球模型逐一聚合到不同網格重測：均勻模型在 L8（65,536 格）開始失敗、WHEEL 在 L7、GEAR1 與 KJSS 在 L5（1024 格）失敗；在**所有多解析度網格上，沒有一個模型通過 S-test**。改用 binary S-test（較不受叢集影響）後有幾個模型多通過幾個網格，但多數多解析度網格仍全數失敗——說明**短期叢集不是模型失敗的唯一原因**，空間結構本身就不對。
  另一個有趣的極端：網格只剩 1 格時，S-test 退化成 N-test，而且因為 S-test 的正規化條件（模擬事件數固定為 $N_{\text{obs}}$）必然通過。

- **教學上可用的洞見**：
  - **「沒被拒絕 ≠ 模型好」的教科書級案例**：一個宣稱「地球上每個地方地震機率都一樣」的模型通過了 CSEP 空間檢驗。任何「我的模型通過所有一致性檢驗」的宣稱，都必須附上「這個檢驗在此設定下的功效是多少」。
  - **通過檢驗要問功效，被拒絕才是資訊**：低功效的檢驗只有拒絕時才有意義；不拒絕基本上什麼也沒說。這與醫學／心理學界的 underpowered study 危機是同一個問題（作者確實引了 Button et al. 2013）。
  - **設計比資料更能改善檢驗**：觀測資料量無法加速累積（300 年！），但網格是我們自己選的。同一批 651 顆地震，換成資料驅動的多解析度網格就能得到有力的檢驗。**教學訊息：檢驗設計是研究者的自由變數，應該當作實驗設計的一部分認真對待。**
  - S-test 有系統性偏好：它**偏袒平滑（uniform-ish）的模型**，懲罰那些敢於指定精確位置的模型。這與上一篇 Bayona 2023 的 Gini coefficient 發現互相呼應，但方向相反——值得在課堂上放在一起討論。
  - 作者的推薦：未來 CSEP 實驗應採用 **Quadtree 資料驅動多解析度網格，讓可用資料決定解析度**。

- **與台灣的關聯**：文中未提及台灣。惟台灣測試區面積小、目標地震數少，正是本文所警告的「低功效」情境；台灣的 CSEP 型實驗特別應該考慮多解析度網格。

---

## [2022] Savran pyCSEP — pyCSEP: A Python Toolkit for Earthquake Forecast Developers

- **書目**：W. H. Savran, J. A. Bayona, P. Iturrieta, K. M. Asim, H. Bao, K. Bayliss, M. Herrmann, D. Schorlemmer, P. J. Maechling, M. J. Werner (2022), *Seismological Research Letters* 93(5), 2858–2870. DOI: 10.1785/0220220033

- **這篇在做什麼**：
  介紹 pyCSEP——把過去封閉在 CSEP testing center 裡的評估程式碼，重寫成開源、模組化、物件導向的 Python 套件。四大模組：（1）地震目錄存取與處理、（2）機率式預報的表示、（3）統計檢驗、（4）視覺化與其他工具。文章附一個 reproducibility package，讀者用一道指令就能重現文中所有圖。

- **關鍵觀念與公式**：
  - **CSEP 的存在理由**：地震預測研究長期背負負面觀感（Geller 1997），根本問題是缺乏**可重現性**（無法重新產生當初發出的預報／檢驗結果）與**可複製性**（不同資料無法得到相同的模型技巧結論）。同儕審查不足以保證這兩件事。CSEP 的解方是**前瞻檢驗**：所有模型參數、預報規格、目標資料來源都必須在觀測之前明確定義，達成 **zero-degree-of-freedom 的獨立檢驗**。
  - **兩種預報表示法**：
    - **grid-based forecast**（`GriddedForecast`）：時空規模箱裡的期望率，用 Poisson likelihood 家族檢驗。
    - **catalog-based forecast**（`CatalogForecast`）：一大堆模擬目錄（例如 UCERF3-ETAS），用經驗分布捕捉預報自身的不確定性，**不需要 Poisson 假設**。
    這個二分法是理解整個檢驗生態系的關鍵：**預報的表示形式決定了可用的檢驗方法。**
  - **grid-based 一致性檢驗**（N / S / M / cL）與上文相同，此處補充 pyCSEP 的說明重點：S、M、cL 三個檢驗都把模擬事件數固定為 $N_{\text{obs}}$，以移除對預報率的依賴；cL-test 實質上是 S 與 M 的組合。
  - **比較檢驗 T-test 與 W-test**：pyCSEP 提供兩者，都基於 IGPE，

  $$\mathrm{IGPE} = \frac{1}{N}\sum_{i=1}^{N}\big(X_i - Y_i\big) - \frac{N_A - N_B}{N}$$

  （$X_i = \ln\lambda_A(k_i)$、$Y_i = \ln\lambda_B(k_i)$ 為兩模型在第 $i$ 顆地震所在箱 $k_i$ 的 log-likelihood；$N_A, N_B$ 為期望總數。符號寫法與 Bayona 2022 的 eq. 6 等價，只是 benchmark 的方向相反。）T-test 用 Student t 分布檢定 IGPE 是否顯著異於 0。**W-test 本組論文只點名「Rhoades et al. (2011) 的 T- 與 W-test」而未細述**；一般理解是 Wilcoxon signed-rank 的無母數版本，不假設每顆地震的資訊增益差呈常態，用以佐證 T-test 的結論穩健性——這一點在本組六篇中找不到文字支持，教材若要寫應標註出處為 Rhoades et al. (2011) 原文。
  - **catalog-based 檢驗**（Savran et al. 2020）：N-test 比對非 Poisson 的數量分布；M-test 比對逐箱的增量規模分布差；S-test 比對目標事件率的**幾何平均**；另有仿連續點過程 likelihood 的 pseudolikelihood test。這些是前述一致性檢驗的類比物，但**放寬了 Poisson 假設**。輸出的兩個統計量 $\delta_1$（預報目錄中事件數 ≥ 觀測的比例）與 $\delta_2$（≤ 觀測的比例）構成雙尾判準。
  - **reproducibility package** 的定義：包含重建論文所有圖表所需的程式碼、資料與其他實驗產物；本文的作法是（1）從 Zenodo 取回並驗證資料、（2）建立含指定 pyCSEP 版本的 Docker image、（3）一道指令重跑全部圖。

- **教學上可用的洞見**：
  - **pyCSEP 能做什麼（觀念層級）**：載入／過濾地震目錄 → 定義測試區與規模箱 → 載入 grid-based 或 catalog-based 預報 → 跑一致性檢驗與比較檢驗 → 畫出診斷圖（N-test 長條、S-test 直方圖、T-test 資訊增益帶信賴區間、預報地圖）。學生不必記 API，但必須理解「區域（region）→ 預報（forecast）→ 目錄（catalog）→ 評估（evaluation）」這條資料流。
  - **可重現性是方法論的一部分，不是附錄**。CSEP 這個社群之所以誕生，直接原因就是先前的預測研究無法被重現。教學時可以把「附上 reproducibility package」當作作業要求。
  - 圖 4/5 是很好的讀圖練習：N-test 圖中綠方塊 = 落在 95% 區間內、紅圓 = 不一致；T-test 圖中水平虛線 = benchmark、誤差棒跨過虛線 = 統計上無法區分。**「誤差棒跨線就是平手」比任何文字說明都有效。**
  - 顯著水準的選擇（$\alpha = 0.01$ 或 0.05）在歷史實驗中並不統一，pyCSEP 的立場是把一致性檢驗當**診斷工具**而非生死判決。

- **與台灣的關聯**：文中未提及台灣（歷史上的 CSEP testing center 設在加州、紐西蘭、義大利、日本、中國）。

---

## [2024] Graham New pyCSEP — New Features in the pyCSEP Toolkit for Earthquake Forecast Development and Evaluation

- **書目**：K. M. Graham, J. A. Bayona, A. M. Khawaja, P. Iturrieta, F. Serafini, E. Biondini, D. A. Rhoades, W. H. Savran, P. J. Maechling, M. C. Gerstenberger, F. Silva, M. J. Werner (2024), *Seismological Research Letters* 95(6), 3449–3463. DOI: 10.1785/0220240197

- **這篇在做什麼**：
  Savran et al. (2022) 的續篇，把社群這兩年新增的功能整理成一篇：新的權威目錄接口（義大利 BSI、紐西蘭 GeoNet、全球 Global CMT）、Quadtree 多解析度網格、非 Poisson 檢驗、把 GEAR1 投影到任意區域當 benchmark 的功能，以及**警報式（alarm-based）模型的評估工具**。以紐西蘭案例貫穿示範。

- **關鍵觀念與公式**：

  **1. Scoring rule 的觀念框架（本文最有價值的一節）**
  一致性檢驗與模型排名，本質上都建立在 **scoring rule** $S(P, D)$ 之上——一個把「預報 $P$」與「資料 $D$」映到分數的函數。**Proper score** 的定義：當被評估的預報最接近真實分布時，該預報得到最高（正向定義）分數。observed joint log-likelihood 就是一個 proper score。用 improper score 排名機率式地震預報會產生系統性偏誤（Serafini et al. 2022）。
  重要澄清：**分數本身不告訴你模型「對不對」，只告訴你哪個模型給觀測資料更高的 likelihood**。一致性檢驗與排名是兩件事。

  **2. 非 Poisson 一致性檢驗（pyCSEP 已內建）**
  - `binomial_number_test()`：negative-binomial N-test，變異數由與測試期等長的歷史不重疊期估計（紐西蘭 8 年用 $\sigma_{8\text{yr}} \approx 93.2$）。
  - `binary_spatial_test()` 與 `binary_conditional_likelihood_test()`：Bayona et al. (2022) 的 binary likelihood 版本。

  **3. Brier score**
  用於二元事件（例如「這格至少發生一顆超過規模門檻的地震」）的 proper score。觀測 $y \in \{0,1\}$、預報機率 $p \in [0,1]$，分數本質上是 **預報機率與觀測之間的平方差**，跨箱加總 $S_B(\mathbf p, \mathbf y) = \sum_i S_B(p_i, y_i)$。
  （註：原文此處的常數項在 pdftotext 轉檔中損毀，此處只保留觀念層面的描述，實作請回查原文。）
  **相對於 log score 的優勢**：地震預報在低活動區常出現機率為零的箱，log 類分數會給出 $-\infty$ 的無限懲罰；Brier score 基於平方差，**永遠有限**，懲罰較溫和。作者建議 **Brier score 與 log score 併看**，才能看到完整的模型表現輪廓。

  **4. Kagan information score $I_1$**
  相對於「均勻 Poisson 過程」這個參考模型的平均資訊增益，以 bit 為單位：

  $$I_1 = \frac{1}{n}\sum_{i=1}^{n}\log_2\!\frac{\lambda(x_i)}{E}$$

  其中 $\lambda(x_i)$ 是第 $i$ 顆觀測地震所在箱的預報率、$E > 0$ 是參考均勻 Poisson 過程的率。**注意 $\lambda(x_i)$ 不正規化到觀測事件數。**
  與 Rhoades IGPE 的差別：$I_1$ 是**對照一個固定的無資訊參考模型**（用來排名一群模型），IGPE 是**兩個模型的成對比較**，且 IGPE 額外含有期望總數差的修正項。
  $I_1$ 的兩個限制：（a）完全不管「沒有地震發生的格子」裡預報了多少率；（b）**只要有一顆地震落在率為零的箱，分數就掉到 $-\infty$**。常見補救是為零率箱設一個 "water level" 底線率。這個陷阱對所有含 $\log$ 的分數都適用。

  **5. 警報式（alarm-based）檢驗：從 contingency table 到 area skill score**
  作法是把連續的預報率用一個門檻 $\lambda_{\text{th}}$ 二值化：率高於門檻的格子「發布警報」。地震落在警報格內算 predicted，落在外面算 missed。於是得到 2 × 2 列聯表：

  | | 觀測 Yes | 觀測 No |
  |---|---|---|
  | **警報 Yes** | (a) hits 命中 | (b) false alarms 誤報 |
  | **警報 No** | (c) misses 漏報 | (d) correct negatives 正確不報 |

  **逐步改變 $\lambda_{\text{th}}$ 就得到一整族列聯表**——這是所有警報式檢驗的共同起點。

  - **ROC diagram**：以 hit rate $H = a/(a+c)$ 對 false alarm rate $F = b/(b+d)$ 作圖。對角線 $H = F$ 代表隨機模型（例如 SUP，spatially uniform Poisson）。曲線愈往左上角 $(0,1)$ 靠愈好；落到對角線以下就是**比亂猜還差**。
  - **concentration ROC**：把格子按預報率由大到小排序，率與面積各自正規化到總和為 1 再累積作圖，用來看「預報的活動集中度」與「觀測地震集中度」是否吻合（Kagan 2009）。
  - **ROC 的弱點**：它隱含以「地震在空間均勻分布」為參考模型，而這幾乎從不成立；且未考慮地震的空間叢集性。
  - **Molchan diagram**：對每個 $\lambda_{\text{th}}$，畫**漏報率** $\nu$ 對**警報所佔時空體積比例** $\tau$：

  $$\tau = \frac{a+b}{a+b+c+d},\qquad \nu = \frac{c}{a+c}$$

  門檻極高 → 沒有任何警報 → 點 $(\tau,\nu) = (0,1)$（左上角）；門檻極低 → 全區警報、抓到所有地震 → 點 $(1,0)$（右下角）。連接兩點的對角線 $\nu = 1-\tau$ 代表「按時空體積比例隨機發警報」的表現。**軌跡愈往左下角，模型愈好；越過對角線就輸給隨機猜測。**
  Molchan diagram 勝過 ROC 之處：**可以自由選擇參考模型**。特別是可以用該格的長期地震活動度**加權**警報體積，讓「在容易預測的地方發警報」不再佔便宜。
  - **Area skill score（Zechar & Jordan 2008, 2010）**：把 Molchan 軌跡積分成單一分數，

  $$\mathrm{AS}(\tau) = \frac{1}{\tau}\int_{0}^{\tau}\big[1 - \nu(t)\big]\,\mathrm{d}t$$

  沿整條軌跡（$\tau: 0 \to 1$）計算即得軌跡上方的面積，值域 $[0,1]$，**隨機參考模型的 AS = 0.5**。AS 的優點是**一次涵蓋所有警報門檻**，不必人為挑一個操作點。

  **6. Quadtree 網格整合**：以 $N_{\max}$（每格最多資料點數）與 $L_{\max}$（最大縮放層級）兩個參數，從目錄自動生成資料驅動的多解析度網格；紐西蘭示範中以 1985–2006 的 $M_w \ge 3.95$ 建格，$N_{\max} = 100$、$L_{\max} = 12$。
  一個誠實的警語：2010 年 Canterbury 序列的起始位置，在依 2006 年前資料建立的 Quadtree 網格上恰好是**低解析度區**——資料驅動的網格會繼承歷史資料的盲點。

  **7. GEAR1 投影功能**：使用者提供區域座標、區域平均 b 值、GEAR1 原始檔、面積檔與網格檔，即可把全球 $M_w$ 5.95+ 率投影／外推到任意區域與規模門檻，作為區域模型的 benchmark。

- **教學上可用的洞見**：
  - **likelihood 為何是共同貨幣**：所有這些工具——Poisson / binary / negative-binomial 一致性檢驗、T-test 的 IGPE、Kagan $I_1$、catalog-based pseudolikelihood——底層都是「模型賦予**實際發生的那些事件**多少機率」。likelihood 提供了跨模型、跨方法的共同尺度；Brier score 與 Molchan/AS 則是刻意跳出 likelihood 框架的補充視角（前者換成平方差、後者換成決策導向的 hit/miss 權衡）。
  - **機率式與警報式是兩種語言**：機率式預報回答「哪裡的率多高」，警報式回答「要不要發布警示」。門檻 $\lambda_{\text{th}}$ 就是兩者之間的轉譯器。**任何機率式預報都可以被二值化後用警報式方法評估**——這對政府決策情境（要不要疏散、要不要提高警戒）特別重要。
  - **$-\infty$ 陷阱**是實作課必講的一課：任何含 $\log$ 的分數，只要有一顆地震落在零率箱就全盤崩潰。這在教學上可以做成一個「故意打壞模型」的練習。
  - **不同分數懲罰不同性質**：作者明講，用一組多樣化的檢驗，目的是**搞清楚模型在哪個面向失敗**，以便改進模型；也讓不同需求的使用者（在意時間窗、規模門檻或機率水準）各取所需。**不要追求「一個總分」。**

- **與台灣的關聯**：文中未提及台灣。

---

## [2020] Rhoades Entropy — The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data

- **書目**：D. A. Rhoades, S. J. Rastin, A. Christophersen (2020), *Entropy* 22(11), 1264. DOI: 10.3390/e22111264
  （注意：檔名 `entropy-22-01264` 的 22 是**卷號**不是年份，本文出版於 2020 年 11 月 6 日。）

- **這篇在做什麼**：
  嚴格說這**不是一篇檢驗方法論文，而是 EEPAS 模型族的論文**，主題是「目錄前置時間（lead time）」對中期預報的影響。EEPAS（Every Earthquake a Precursor According to Scale）建立在 precursory scale increase（Ψ 現象）之上：大地震之前，周邊小地震的規模與發生率會上升，且遵循可預測的尺度關係。規模愈大的目標地震，前兆時間愈長——長到可能**早於目錄起始時間**，這時前兆貢獻就不完整了。作者推導出「前兆貢獻完整度」關於目標規模與 lead time 的公式，並提出兩個新版本：FLEEPAS（固定 lead time，用來檢視 lead time 的效應）與 FLCEEPAS（額外補償前兆不完整性）。套用於紐西蘭資料。

- **為何收進「預報檢驗」這一章**：
  它示範了 **information gain 作為模型比較工具的實際用法**，而且用的是連續時空點過程（非網格化）的版本，正好與 CSEP 的網格化框架形成對照。模型內容本身應歸到網站的 **EEPAS／Ψ 現象**章節（參見 reference/ 下的 `[2004] rhoades2004`、`[2007] Rhoades`、`Taiwan__[2026] Taiwan_EEPAS` 等）。

- **關鍵觀念與公式**：
  - **點過程的 log-likelihood**（連續版本，與 CSEP 的網格版本對照）：

  $$\ln L = \sum_{i=1}^{N}\ln\lambda(t_i, m_i, x_i, y_i) - E_{\lambda}(N)$$

  $\lambda$ 為率密度（rate density），$E_\lambda(N)$ 是把率密度對時間範圍 $(t_1, t_2)$、規模範圍 $(m_c, m_u)$ 與監測區 $R$ 積分所得的期望目標地震數。這裡也是「用參數擬合時，log-likelihood 就是被最佳化的目標函數」。

  - **模型 X 相對於模型 Y 的 information gain**：

  $$I(X, Y) = \frac{\ln L_X - \ln L_Y}{N}$$

  分母除以 $N$，所以單位是 **per earthquake**——與 CSEP 的 IGPE 是同一個概念的連續版本。本文中的參考模型 $Y$ 是 **SUP（Stationary Uniform Poisson）**：率密度只取決於目標地震數與 Gutenberg–Richter b 值，是「最小資訊」的基準線。

  - **Ψ 現象的量化**：累積規模異常（cumulative magnitude anomaly, cumag）

  $$C(t) = \sum_{t_s \le t_i < t}\big[M_i - (M_{\text{thres}} - 0.1)\big] - k(t - t_s)$$

  $$k = \sum_{t_s \le t_i < t_f}\big[M_i - (M_{\text{thres}} - 0.1)\big] \big/ (t_f - t_s)$$

  每顆 $M_i \ge M_{\text{thres}}$ 的地震讓 $C(t)$ 向上跳，減項是長期平均趨勢；活動低於平均時 $C(t)$ 下降。$C(t)$ 的上升段即前兆活化。

  - **混合模型與 $\phi$ 參數**：以兩個端點模型的凸線性組合補償前兆不完整，

  $$\lambda_C(t, m, x, y) = \phi\,\lambda_A(t, m, x, y) + (1-\phi)\,\lambda_B(t, m, x, y),\qquad 0 \le \phi \le 1$$

  $\phi$ 由最大概似法擬合。$\phi$ 接近 1 表示大部分前兆貢獻仍然完整。

  - **主要結果**：未重新擬合的 FLEEPAS，在 lead time 短於 11 年時 information gain 急遽下降（3 年時掉了約 0.7）；重新擬合的 FLEEPAS 只掉約 0.2。此外發現**時空權衡（space-time trade-off）**：lead time 縮短時，擬合出的時間尺度因子 $10^{a_t}$ 下降、空間尺度因子 $\sigma_A^2$ 上升——前兆訊號在時間上壓縮時，會在空間上攤開。FLCEEPAS 只擬合一個參數 $\phi$，就能達到與擬合三個參數的 FLEEPAS 相近的 information gain。

- **教學上可用的洞見**：
  - **information gain 是本章的樞紐概念**：它同時出現在網格化的 CSEP T-test（IGPE）、Kagan $I_1$ 與此處的連續點過程 $I(X,Y)$。三者都是「每顆地震的平均 log-likelihood 差」，只是參考模型與離散化方式不同。**教學時先建立這個共同結構，再談各自的變體，遠比逐一介紹省力。**
  - **參考模型（baseline）的選擇決定了數字的意義**：對 SUP 的 gain 與對 HKJ 的 gain 完全不是同一回事。看到「information gain = 0.5」時，第一個問題永遠是「相對於什麼？」
  - **每多擬合一個參數，就多欠一分懷疑**：FLCEEPAS 用 1 個參數達到 3 個參數的效果——這是模型節儉性的正面示範，也呼應了 IGPEc 的 AICc 懲罰項為何必要。
  - 誠實的方法論警示：本文的 information gain 是在**擬合資料上**計算的（模型參數對同一批紐西蘭目錄調過），因此屬於回溯性能。要主張真實預報技巧，仍須經過 CSEP 式的前瞻檢驗——這正好是本章的核心論點。

- **與台灣的關聯**：文中未提及台灣（案例為紐西蘭）。惟 EEPAS 已被套用於台灣（見 reference/ 下 `Taiwan__[2026] Taiwan_EEPAS`），該工作若要建立可信度，適用的正是本章所述的 CSEP 前瞻檢驗流程。

---

## 跨篇綜合：這個主題教什麼

### 一、建議的教學順序

1. **先講「為什麼需要檢驗」**——從 CSEP 的起源講起（pyCSEP 2022 開頭那段）。地震預測研究長期被質疑，根本問題不是想法不夠聰明，而是**結果無法被重現、模型技巧無法被獨立驗證**。CSEP 的解方是把預報變成 zero-degree-of-freedom 的可否證陳述：模型參數、預報格式、目標資料來源，全部在觀測之前釘死。
2. **建立預報的資料結構**：測試區 → 空間格 × 規模箱 → 每箱的期望數 $\lambda$。這一步決定了後面所有統計的形式。同時點出另一條路線：catalog-based forecast（一堆模擬目錄），它不需要 Poisson 假設。
3. **一次講清楚 Poisson likelihood**：$P(\omega|\lambda)$ → POLL → jPOLL。這是後面所有檢驗的原子。
4. **拆解一致性檢驗 N / M / S / cL**：教學關鍵是講清楚**每個檢驗刻意抽掉了什麼**——N-test 只留總數；S-test 正規化掉總數只留空間形狀；M-test 只留規模分布；cL-test 檢空間×規模的聯合結構。順帶說明 L-test 為何被冷落（對總數過度敏感，混淆了兩種失敗模式）。
5. **講模擬程序**：沒有解析分布，就用模擬造出分布。這是理解「分位數分數 = p-value」的必要一步，也是學生最容易卡住的地方。
6. **進入比較檢驗**：IGPE 與 T-test。強調從「模型 vs. 資料」轉成「模型 vs. 模型」的視角轉換。
7. **放寬假設**：negative binomial N-test、binary S/cL-test、catalog-based tests。此時學生已經知道 Poisson 錯在哪（變異數不足以容納叢集），新方法才會有意義。
8. **警報式檢驗**：contingency table → ROC → Molchan → area skill score。放在最後，作為「另一種語言」呈現。
9. **最後上功效與多重檢定**：Khawaja 2023 與 Bayona 2022 的 Bonferroni 段落。學生此時已經跑過檢驗，才會真正被「均勻模型通過 S-test」震撼到。
10. **pyCSEP 動手**：把上面的觀念在紐西蘭或義大利的公開資料上跑一遍。

### 二、各檢驗的互補性

一致性檢驗是**診斷矩陣**，不是及格線。它們刻意設計成各自針對預報的一個面向：

| 檢驗 | 檢什麼 | 抽掉什麼 | 典型失敗訊息 |
|---|---|---|---|
| N-test | 總地震率 | 空間、規模 | 模型高估／低估活動度 |
| M-test | 規模分布 | 空間、總數 | b 值或截止規模不對 |
| S-test | 空間圖形 | 規模、總數 | 震央落在模型沒預期的地方 |
| cL-test | 空間 × 規模聯合 | （不正規化總數） | 綜合的空間規模結構不合 |
| T-test / W-test | 相對於另一模型的資訊增益 | — | 沒比 benchmark 好 |
| Molchan / AS | 全門檻的警報表現 | — | 警報效率不如隨機 |

Bayona 2022 用模擬量化了互補性：**S-test 與 cL-test 高度相關（合理，cL 含 S），兩者都與 N-test 幾乎獨立（$R \approx 0.01$–$0.03$）**。所以「數量」與「形狀」是兩個真正獨立的資訊軸，Bonferroni 校正時有效檢驗數約為 2 而非 4。

不同分數也互補：log score 對零機率箱給無限懲罰、Brier score 永遠有限、Kagan $I_1$ 不管沒地震的格子、Molchan/AS 完全跳出 likelihood 改用決策視角。**併看，才能定位模型失敗在哪個面向。**

### 三、prospective vs. retrospective 為何是 CSEP 的核心精神

這是整章最重要的一句話，而且有硬證據：

- **Bayona 2022**：16 個 multiplicative hybrid 在回溯評估中對 HKJ 有 0.2–0.5 的顯著 information gain（而且已用 AICc 懲罰過參數數量）。在 2011–2020 的前瞻期，**沒有任何一個顯著優於 HKJ**。作者的診斷是「形成 hybrid 的擬合在時間上不穩定」——換句話說，回溯期的優勢是對特定時期地震分布的過擬合。
- **Bayona 2023**：義大利的模型排名與 Taroni et al.（2018）不同，主因只是測試期是否涵蓋 2012 Emilia 序列。**單一地震序列可以翻轉整份排名。**

由此推出的教學結論：

1. **回溯評估只是 sanity check**（Bayona 2022 引 Werner et al. 2010 的原話），永遠不足以宣稱預報技巧。
2. **前瞻性必須是制度性的**，不能靠自律。所以 CSEP 要有 testing center、要在觀測之前釘死所有規格、要 zero degree of freedom。「我保證沒偷看資料」不是科學程序。
3. **前瞻資料累積極慢**，這是這個領域無法迴避的結構性困難。CSEP 的因應是多區域並行實驗與全球實驗（Bird et al. 2015 的 GEAR1 路線）——擴大空間換取時間。
4. **可重現性套件**（reproducibility package）是前瞻精神的技術落實：把程式碼、資料、環境（Docker）、版本一起凍結。

### 四、常見誤解與警語

**誤解 1：「模型通過了所有一致性檢驗，所以它是好模型。」**
最強的反例來自 Khawaja 2023：一個宣稱「地球上每處地震機率相同」的均勻模型，在 0.1° 全球網格上**通過 S-test**。檢驗沒有拒絕，往往只代表**檢驗沒有力氣拒絕**。正確的說法是：「在此網格與此樣本數下，資料不足以拒絕該模型」，而且應該同時報告該設定下的統計功效。

**誤解 2：「likelihood 分數低就是模型爛。」**
少數叢集事件會主宰整個分數。Bayona 2022 中，2016 Hawthorne 群震與 2019 Ridgecrest 序列擠在極少數格子裡，Poisson S-test 對這幾格的重罰幾乎決定了所有模型的成敗。這正是 binary likelihood 被提出的動機——只問「這格有沒有」而不問「有幾顆」，讓 $\omega \ge 2$ 的箱不再擁有過大的話語權。

**誤解 3：「解析度愈高，檢驗愈嚴格。」**
恰好相反。格子細到每顆地震各佔一格時，空間資訊被稀釋，S-test 失去鑑別力。Khawaja 2023 估算，在慣用的全球高解析度網格上要達到有力檢驗需要 32,000 顆地震（約 300 年）；改用資料驅動的多解析度 Quadtree 網格，**8 顆就夠**。網格是研究者的自由變數，應當作實驗設計認真對待。

**誤解 4：「information gain = 0.5，所以模型好。」**
永遠要問「相對於什麼 baseline」。相對於 SUP（最小資訊模型）的 gain 與相對於 HKJ（當時最佳模型）的 gain 完全不是同一個量級的成就。同樣地，Kagan $I_1$ 的參考是均勻 Poisson 過程，與成對比較的 IGPE 不可混用。

**誤解 5：「同時跑很多檢驗比較保險。」**
跑愈多檢驗，至少一個假陽性的機率愈高。要做多重檢定校正（如 Bonferroni $\alpha/T$），而 $T$ 應該是**有效獨立檢驗數**——Bayona 2022 用模擬證明 S 與 cL 高度相關，所以 $T = 2$ 而非 4。

**誤解 6：「模型愈在地、資料解析度愈高，一定愈準。」**
Bayona 2023 在三個全世界儀器最密的地震區測試，全球模型 GEAR1 排名 1/2/3。地方資料的優勢是需要被證明的假設。

**誤解 7：「檢驗不通過就要淘汰模型。」**
CSEP 自己的立場是：**不因為某個檢驗失敗就正式 reject 模型**，而是把分位數分數當診斷指標，指出模型與資料在哪裡不合、值得科學上追究。檢驗的目的是**改進模型**，不是頒發及格證書。

**額外的實作陷阱：**
- **$\log$ 的 $-\infty$ 陷阱**：一顆地震落在零率箱，整個分數崩潰。實務上為零率箱設 "water level" 底線率。
- **對齊定義**：規模尺度（$M_w$ vs. $M_L$）、深度範圍、是否去餘震——不對齊就沒有可比性。Bayona 2023 為此逐項檢查、逐項換算，是好的示範。
- **資料驅動的網格會繼承歷史盲點**：2010 Canterbury 序列的起始位置，在依 2006 年前資料建立的 Quadtree 網格上恰好落在低解析度區。
- **回溯擬合過的模型要用 IGPEc 而非 IGPE**；前瞻期目標事件與擬合資料獨立時，才可用未懲罰的 IGPE。

### 五、與台灣的接點

本組六篇論文**均未提及台灣**（案例為加州、紐西蘭、義大利、全球）。方法論上的接點是清楚的：台灣測試區面積小、M 4.95+ 目標事件數少，正落在 Khawaja 2023 所警告的低功效區間，因此（a）任何「台灣模型通過檢驗」的宣稱都必須附上功效評估，（b）多解析度 Quadtree 網格對台灣特別值得考慮，（c）GEAR1 投影到台灣測試區可以作為現成的 benchmark。網站的台灣章節（`Taiwan__[2026] Taiwan_EEPAS` 等）若要建立可信度，正是要走這一套前瞻檢驗流程。
