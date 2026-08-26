# Ensemble／Hybrid 模型論文筆記
> 主題：混合與加乘預報模型 | 來源：reference/Ensemble/ 4 篇 | 供第二部 14 章（ensemble）等使用

---

## [2009] Mixture Models — Mixture Models for Improved Short-Term Earthquake Forecasting

- **書目**：Rhoades, D. A., & Gerstenberger, M. C. (2009). *Bulletin of the Seismological Society of America*, 99(2A), 636–646. DOI: 10.1785/0120080063

- **這篇在做什麼**：把兩個「時間尺度完全不同」的模型湊在一起看能不能更好。STEP 是短期餘震群聚模型（Omori–Utsu 衰減 + Gutenberg–Richter），主力是預報餘震；EEPAS 是中長期模型，靠 precursory scale increase 提前數月到數十年給出大地震的機率抬升。作者用加州 ANSS 目錄 1984–2004 的 152 個 $M \ge 5.0$ 事件，以最大概似擬合各種線性混合，比較各混合模型的資訊量。

- **關鍵觀念與公式**

  兩個母模型本身就已經是「混合」的結構，這點值得先講清楚：

  $$\lambda_{\mathrm{STEP}} = \max\left(\lambda_{\mathrm{CLUST}},\ \lambda_{\mathrm{STAT}}\right)$$

  $$\lambda_{\mathrm{EEPAS}} = \lambda_{\Psi} + \mu\,\lambda_{\mathrm{PPE}}$$

  符號：$\lambda_{\mathrm{CLUST}}$ 是 STEP 的時變群聚項，$\lambda_{\mathrm{STAT}}$ 是靜態背景（平滑地震度）；$\lambda_{\Psi}$ 是 EEPAS 的前兆尺度增長項，$\lambda_{\mathrm{PPE}}$（proximity to past earthquakes）是準靜態背景，$\mu \in [0,1]$ 代表「沒有前兆尺度增長的大地震比例」。EEPAS 因此本身就是「時變項 + 背景項」的加法混合。

  本文新造的混合模型有三型：

  $$\lambda_{\mathrm{SE1}} = \lambda_{\mathrm{CLUST}} + q\,\lambda_{\mathrm{EEPAS}} \qquad (0 \le q \le 1)$$

  $$\lambda_{\mathrm{SE2}} = (1-r)\,\lambda_{\mathrm{STEP}} + r\,\lambda_{\mathrm{EEPAS}} \qquad (0 \le r < 1)$$

  $$\lambda_{\mathrm{SE3}} = \lambda_{\mathrm{CLUST}} + p(m_j)\,\lambda_{\mathrm{EEPAS}}$$

  - SE2 是**凸組合（convex combination）**：權重和為 1，因此混合後的總期望地震數自動維持在區域平均地震率上，不需另外正規化。這是「凸」這個條件在預報上的實際意義——它保證 N-test 層次的守恆。
  - SE3 讓權重隨規模變化：$p(m)$ 是「該事件是主震（獨立事件）的機率」，用 logistic regression 從 Reasenberg declustering 的獨立性機率擬合：

    $$\ln\frac{p}{1-p} = a_l + b_l m,\qquad \hat a_l = -0.51,\ \hat b_l = 0.19$$

    亦即規模愈大愈可能是主震 → EEPAS 的權重愈高（M5 約 60%，M7 約 70%）。概念很漂亮，但結果顯示它相對 SE1 並沒有實質增益（$\Delta \ln L$ 只有 0.1–0.2），最後仍偏好較簡單的 SE1。這是個很好的「奧坎剃刀」教材。

  - 權重怎麼定：全部用**最大概似**，目標函數是空間–時間–規模格點上的 Poisson 對數概似

    $$\ln L = \sum_{i=1}^{N} \ln \lambda(t_{j_i}, m_{j_i}, x_{j_i}, y_{j_i}) - \sum_{j=1}^{n} \lambda(t_j, m_j, x_j, y_j)$$

    第一項是「有地震的格子」的貢獻，第二項是全區總期望數。用 Nelder–Mead simplex 最佳化，再用 AIC 扣參數：$\mathrm{AIC} = -2\ln L + 2k$，並定義每次地震的資訊率 $\Delta \mathrm{AIC}/(2N)$、平均機率增益 $G = \exp[\Delta \mathrm{AIC}/(2N)]$。

- **量化結果（相對於 stationary uniform Poisson, SUP）**

  | 模型 | 資訊率 $\Delta\mathrm{AIC}/2N$ | 機率增益 $G$ |
  |---|---|---|
  | PPE | 1.67 | 5.3 |
  | EEPAS | 2.17 | 8.8 |
  | STEP | 約 2.1–2.4 | 8–12 |
  | SE1（CLUST + EEPAS） | 約 2.9–3.2 | 18–25 |
  | **SE2（STEP ⊕ EEPAS 凸組合）** | **約 3.2–3.3** | **25–28** |

  最佳解是 $r = 0.42$：**0.58 × STEP + 0.42 × EEPAS**，相對 STEP 的平均機率增益 2.72，相對 EEPAS 也超過 2。單一模型「翻倍」等級的增益，只靠一個權重參數就拿到了。

- **教學上可用的洞見**
  - **投資組合比喻**：STEP 像高波動短線部位（地震發生後機率瞬間飆高再快速衰減），EEPAS 像長期持有（緩慢上升、更緩慢衰減）。兩者用的是同一份地震目錄，但抽取的是完全不同的訊號，所以相關性低——低相關資產組合起來，風險調整後報酬才會提升。文中 Figure 1 的「尖峰 vs 緩坡」對照圖就是最好的教學圖。
  - **「爛模型也可能有用」**：STEP 的靜態背景項整體上不如 PPE，但把它整個拿掉（SE1 vs SE2）反而損失資訊。因為它在少數格子上有其他模型沒有的強特徵。教學警語：不要只看模型的總體排名就把它踢出 ensemble。
  - **兩個模型 rate 比值跨 12 個數量級**（STEP/PPE），一半落在 0.1–10 之間。差異大 = 互補性強 = 混合的空間大。可以做成 histogram 當教學圖。
  - **警語一**：本文是**回溯（retrospective）擬合**，模型是在同一份資料上調參後再比較，作者自己說「文獻裡回溯排名沒被前瞻測試證實的例子俯拾即是」。這正好接到 2022 那篇。
  - **警語二**：STEP 的對數概似因為計算量太大而用抽樣（每 60 天取一次）估計總期望數，帶有 bootstrap 信賴區間，所以表中很多數字是區間而非點值。誠實報告計算限制，也是教材。
  - 結論最後留了一個伏筆：Aki（1981）、Utsu（1983）、Imoto（2007）的**條件獨立前兆乘法理論**暗示，如果用乘法而不是加法組合，可能可以保留完整的機率增益。這正是 2014 那篇的起點。

- **與台灣的關聯**：台灣同時有短期 ETAS/STEP 類與中期 EEPAS 類的應用（見 reference/Taiwan/ 的 Dapu ETAS 快報與 Taiwan EEPAS），本文的 SE2 凸組合是台灣做第一個 ensemble 最直接、參數最少的入手方式（只需擬合一個 $r$）。

---

## [2014] RELM II — Regional Earthquake Likelihood Models II: Information Gains of Multiplicative Hybrids

- **書目**：Rhoades, D. A., Gerstenberger, M. C., Christophersen, A., Zechar, J. D., Schorlemmer, D., Werner, M. J., & Jordan, T. H. (2014). *Bulletin of the Seismological Society of America*, 104(6), 3072–3083. DOI: 10.1785/0120140035

- **這篇在做什麼**：RELM 五年實驗（加州，$M \ge 4.95$）的一階結論是 Helmstetter et al. (2007) 的平滑地震度模型（HKJ）最好。這篇問：能不能以 HKJ 為 baseline，把其他模型當成「修正因子」乘上去，做出比 HKJ 更好的模型？Marzocchi et al. (2012) 已經證實 RELM 的**加法**組合（Bayesian model averaging）打不贏最佳單一模型，所以這裡改走乘法路線。

- **關鍵觀念與公式**

  **Multiplicative hybrid 的基本形式**：baseline 的每格期望數保留不動，其他模型（conjugate models）轉換成乘數：

  $$\lambda_H(j,k) = \lambda_1(j,k)\,\exp\left[a + \sum_{i=2}^{n_i} f_i\big(\lambda_i(j,\cdot)\big)\right]$$

  $$f_i(\lambda) = b_i\big[\log(1+\lambda)\big]^{c_i},\qquad b_i \ge 0,\ c_i > 0$$

  符號：$j$ 是空間格、$k$ 是規模格；$\lambda_1$ 是 baseline（HKJ）的格點期望數；$\lambda_i(j,\cdot) = \sum_k \lambda_i(j,k)$ 是第 $i$ 個 conjugate 模型**對規模求和後**的空間期望數；$a$ 是整體正規化參數，$b_i, c_i$ 是形狀參數。

  幾個設計要點值得在課堂上逐一拆解：

  1. **對規模求和**：隱含假設是「conjugate 模型除了 baseline 已有的規模資訊外，不再提供額外的規模資訊」。對於遵守空間不變 b 值 Gutenberg–Richter 的模型是合理的。換句話說，**規模分布由 baseline 全權負責，conjugate 只負責調整空間分布**。
  2. **只用排序、不用數值**：求和後的期望數被當成 Zechar & Jordan (2008) 意義下的 *alarm function*，重要的是空間格之間的**排序**，不是絕對值。$f_i$ 是保序（order-preserving）的單調非遞減函數，所以排序被保留，但尺度被重新學習。這一步讓不同單位、不同校正基準的模型（應變率、斷層滑移率、PI 指標）可以塞進同一個框架。
  3. **形式的彈性與節制**：$f_i$ 每個 conjugate 只有兩個參數，在 $\lambda \ll 1$（$\log(1+\lambda) \approx \lambda$）的範圍內可近似線性、上凸、下凹或常數。目標事件只有 31 個，參數再多就是過度配適。
  4. **退化保證**：當 $a = 0$ 且所有 $b_i = 0$ 時乘數恆為 1，hybrid 退化成 baseline。所以最佳化後的 hybrid **對數概似不可能低於 baseline**——這也正是為什麼一定要扣參數懲罰。

  **對照組：加法 hybrid**

  $$\lambda_H(j,k) = \sum_{i=1}^{n} a_i \lambda_i(j,k),\qquad a_i \ge 0$$

  **資訊增益的度量（有懲罰版）**：用 corrected AIC（Hurvich & Tsai, 1989）

  $$\mathrm{AICc} = -2\ln L + 2p + \frac{p+1}{N-p-1},\qquad \mathrm{IGPEc} = \frac{-\Delta}{2N}$$

  展開後可寫成可做 T 檢定的形式：

  $$\mathrm{IGPEc} = \frac{\hat N_1 - \hat N_H}{N} - \frac{1}{2N}\left[2p + \frac{p+1}{N-p-1}\right] + \frac{1}{N}\sum_{n=1}^{N}\left[\ln\lambda_H(j_n,k_n) - \ln\lambda_1(j_n,k_n)\right]$$

  符號：$N$ 是目標地震數、$p$ 是擬合參數數、$\hat N$ 是模型預期的總地震數。第一項罰「總數預報不準」，第二項罰「參數太多、資料太少」，第三項才是真正的「命中格子的對數機率增益」。這個三項分解本身就是很好的教學素材：**資訊增益 = 總量校正 + 複雜度懲罰 + 空間命中**。

- **量化結果**
  - 全加州兩模型 hybrid：HKJ ⊗ Bird & Liu Neokinema（$\Delta\ln L = 11.4$）與 HKJ ⊗ Holliday PI（$\Delta\ln L = 11.2$），IGPEc 皆約 0.25。
  - 南加州：HKJ ⊗ Shen et al. 大地測量模型，增益超過 0.5。
  - 最佳三模型 hybrid：全加州 HKJ ⊗ Neokinema ⊗ PI，IGPEc $= 0.35 \pm 0.17$（勉強達 95% 顯著）；南加州 HKJ ⊗ Shen ⊗ PI，IGPEc $= 0.79 \pm 0.27$（明確顯著）。但**第三個模型帶來的額外增益（$\Delta$IGPEc = 0.09 與 0.22）都不顯著**——邊際報酬遞減。
  - **乘法完勝加法**：全加州用同一批模型做的最佳加法 hybrid 只有 $\Delta \ln L = 4.8$，乘法是 11.4。

- **教學上可用的洞見**
  - **為什麼乘法比加法強**：加法組合本質上是加權平均，某格的期望值幾乎必然落在各成分模型該格期望值的**區間內**；乘法組合則可以跑到所有成分模型的範圍**之外**（更高或更低），前提是資料整體支持這個乘數。一句話：加法只能內插，乘法可以外推。這是本主題最核心的一句教學台詞。
  - **多元迴歸類比**：作者自己下的比喻——乘法 hybrid 就像多元迴歸，每加一個新的解釋變數（新資料、新模型、應變圖、應力變化圖、甚至一個「有沒有前兆」的二元變數），只多付兩個參數，就能解釋更多變異。而且和多元迴歸一樣：**觀測數要多、解釋變數要少**，否則不穩健。
  - **異質性帶來增益**：增益最大的組合，通常是概念或資料來源差異最大的組合。用大地測量資料（Shen、Neokinema、Ward geodetic）當 conjugate 特別有效，因為 baseline HKJ 只用地震目錄。反過來，ALM 與 Ebel 這種「粗糙版的平滑地震度」幾乎沒有增益——它們和 baseline 講的是同一件事。
  - **有效 conjugate 的判準**：不是「它自己單獨很強」，而是「它和地震發生的相關方式，是 baseline 沒有捕捉到的」。PI 模型單獨表現平平，但作為 conjugate 兩模型、三模型都上榜。
  - **警語**：作者反覆強調這是回溯擬合，「回溯測試無論統計量設計得多小心，都無法完全排除過度配適」，所有 16 個 hybrid 已送進 CSEP 做 2011–2015 的前瞻測試。這個伏筆由下一篇揭曉。

- **與台灣的關聯**：台灣的 GNSS 應變率、活動斷層滑移率、地下水／地磁異常等，都可以視為 conjugate 的候選輸入——本文框架不要求 conjugate 是一個完整的預報模型，任何格點化的量（strain map、stress-change map、二元前兆指標、警報等級的類別變數）都能塞進 $f_i$。這對本書「多種觀測資料整合」的主軸是直接可用的方法論。

---

## [2022] Prospective Test — Prospective Evaluation of Multiplicative Hybrid Earthquake Forecasting Models in California

- **書目**：Bayona, J. A., Savran, W. H., Rhoades, D. A., & Werner, M. J. (2022). *Geophysical Journal International*, 229(3), 1736–1753. DOI: 10.1093/gji/ggac018

- **這篇在做什麼**：把 2014 那 16 個乘法 hybrid 加上 6 個原始 RELM 模型，用 2011-01-01 至 2020-12-31 的 40 個 $M \ge 4.95$ 事件（含 2016 Hawthorne 群震與 2019 Ridgecrest 序列）做**真正的前瞻測試**，並額外引進一套「binary likelihood」版本的 CSEP 檢驗來降低群聚事件的影響。結論很殘酷：**沒有任何 hybrid 在前瞻測試中顯著優於 HKJ**。

- **關鍵觀念與公式**
  - hybrid 的定義與 IGPEc 公式沿用 2014（見上），此處不重複。差別在於：因為目標事件與擬合資料獨立，改用不含參數懲罰的 IGPE。
  - **Poisson vs binary 概似**：傳統 CSEP 的 N-test / S-test / cL-test 假設每格地震數獨立且服從 Poisson。但群震序列會讓少數幾格塞進 3–4 個事件，Poisson 概似對此懲罰極重。作者改用**二元概似（binary likelihood）**：只問「這格有沒有事件」，把 $\omega \ge 1$ 都算成一次命中。POLL（Poisson log-likelihood）→ BILL（binary log-likelihood）。
  - 另外用 **negative-binomial N-test** 取代 Poisson N-test，放寬地震總數的過度離散問題。
  - 模型間比較仍用 Rhoades et al. (2011) 的 paired T-test。

- **量化結果與診斷**
  - **N-test**：所有模型都**高估**了地震數，全數未通過 Poisson N-test；只有少數通過 negative-binomial 版。作者判讀為 2011–2020 是加州相對平靜的十年，而非模型全錯。
  - **S-test（空間）**：只有 KAGAN 通過 Poisson S-test；改用 binary 概似後，EBEL-C、KAGAN、SHEN 與三個 Ward 系 hybrid 通過。以 HKJ 為例，Poisson 下 7652 個空格 + 24 個單震格 + 3 個雙震格 + 2 個三震格 + 1 個四震格，其中多震格貢獻了 39% 的懲罰；改成 binary 後這些格子只貢獻 18%。**binary 概似確實有效降低了群聚敏感度，但多數模型仍然過不了空間檢驗。**
  - **T-test**：HKJ 是這十年最有資訊量的模型。回溯時增益 +0.25 / +0.25 / +0.5 的 BIRD、PI、SHEN 三個 conjugate，前瞻時對 HKJ 的 IGPE 變成 **−0.42 / −0.71 / −0.68**。符號完全翻轉。

  作者給出四個可能原因，這四點是本篇最有教學價值的部分：
  1. **時間不穩定性**：2014 的權重是用 RELM 頭五年僅 31 個事件擬合的，參數與信賴區間會隨時間大幅漂移。樣本太小 → 擬合不穩。
  2. **測試期沒有大型斷層上的地震**，這對用斷層／大地測量資料的模型不利。
  3. **conjugate 模型自身退化**：BIRD、PI、EBEL-C 相對 HKJ 的 IGPE 從 −0.70／−0.31／−1.64 掉到 −1.62／−3.38／−2.95。conjugate 變差，hybrid 自然跟著變差。
  4. **乘法結構的放大效應（最關鍵的機制解釋）**：所有 hybrid 共用 HKJ 為 baseline、規模分布相同、總量已正規化，**唯一的差別就是空間分布**。而在乘法結構中，conjugate 空間率低的格子會得到 $<1$ 的乘數，使 hybrid 在該格的率**比 baseline 更低**。一旦目標地震剛好落在這種格子裡，hybrid 的表現就會比 baseline 更糟。2006–2010 期間幾乎所有模型都輕鬆通過 S-test，所以乘法組合佔便宜；2011–2020 期間多數模型（含 HKJ、BIRD、PI）都明確未通過 S-test，乘法組合就反過來吃虧。

- **教學上可用的洞見**
  - **本主題的頭號警語**：回溯資訊增益 ≠ 前瞻資訊增益。這一組（2014 → 2022）是統計地震學裡少見的、完整走完「提出 → 送測 → 十年後開獎」流程的案例，適合當作整章的敘事骨幹：**先講乘法多漂亮，再講它為什麼在真實測試中翻車。**
  - **乘法是雙刃劍**：加法組合有「內插」的天然保護——最差不過是被權重稀釋；乘法組合的乘數可以把率壓到所有成分模型之下，錯得比任何一個成分模型都更離譜。加法輸的是天花板，乘法輸的是地板。
  - **參數少 ≠ 穩健**：每個 conjugate 只多兩個參數看起來很節制，但相對 31 個目標事件仍嫌奢侈。「有效樣本數」在地震預報裡是以**目標地震顆數**計，不是以格子數計——這點初學者最容易搞錯（本例中格子有數萬個，事件只有 31 個）。
  - **群聚讓 Poisson 概似說謊**：可以用 Ridgecrest 序列當例子，展示同一個模型在 POLL 與 BILL 下的分數差多少。這也直接鋪路到第 15 章的檢驗方法。
  - 結尾仍給出建設性方向：**相關的模型用加法組合、獨立的模型用乘法組合**，值得繼續研究。

- **與台灣的關聯**：台灣若要建立 CSEP 式的測試中心，這篇是「測試設計」的範本——特別是台灣地震序列群聚極強（車籠埔、池上、大埔），Poisson 概似的失真問題只會比加州更嚴重，binary likelihood 與 negative-binomial N-test 幾乎是必需品。

---

## [2023] Maximizing Skill — Maximizing the Forecasting Skill of an Ensemble Model

- **書目**：Herrmann, M., & Marzocchi, W. (2023). *Geophysical Journal International*, 234(1), 73–87. DOI: 10.1093/gji/ggad020

- **這篇在做什麼**：問一個看似細節、其實根本的問題——**ensemble 的權重到底該怎麼定？**過去要嘛等權重，要嘛依各模型自己的表現給權重（如 OEF-Italy 用的 Score Model Averaging, SMA）。作者主張這都不對：權重應該**直接最大化 ensemble 的表現**，而不是反映各模型的個別表現。方法是多元 logistic regression。實測對象是義大利的作業型地震預報系統 OEF-Italy。

- **關鍵觀念與公式**

  **一般化的 ensemble 記法**：$f_E = g[\{f_i, \pi_i\}]$，$g[\cdot]$ 是 ensemble 運算子，$f_i$ 是候選預報，$\pi_i$ 是權重。最常見的形式是加權平均（線性組合）：

  $$f_E(x) \equiv \bar f(x) = \sum_{i=1}^{m} f_i(x)\,\pi_i$$

  當 $\pi_i$ 被詮釋成「$f_i$ 是真模型的機率」時，這就是 **Bayesian model averaging**。文獻裡的 mixture model（Rhoades 2009/2013）、hybrid（Rhoades 2014、Bayona 2022）、multimodel ensemble、superensemble 講的都是這個框架的變形。

  **用 logistic regression 學權重**：把觀測二元化——時空格 $j$ 內有沒有至少一個目標事件 $Y_j \in \{0,1\}$：

  $$\Pr(Y_j = 1 \mid \mathbf{u}_j) = p(\mathbf{u}_j) = \frac{1}{1 + e^{-g(\mathbf{u}_j)}},\qquad g(\mathbf{u}_j) = \beta_0 + \beta_1 u_{1,j} + \cdots + \beta_m u_{m,j}$$

  關鍵細節：**自變數是各模型預報率取對數**，$u_{i,j} = \ln \phi_{i,j}$。所以這是在 log-odds 尺度上做線性組合，等價於在率的尺度上做**乘冪組合** $\prod_i \phi_i^{\beta_i}$——logistic regression 骨子裡是個乘法結構。$\beta_0$ 與 $\beta_i$ 由最大化下式求得：

  $$\ell = \sum_j \left[Y_j \ln p_j + (1-Y_j)\ln(1-p_j)\right]$$

  **為什麼不能直接把 logistic 模型當 ensemble 用**（本文最精采的轉折）：
  - $\beta_0$（baseline log-odds）反映的是各模型**絕對**表現的偏移。若過去這些模型在目標格給的率偏高（表現好），$\beta_0 > 0$；反之 $\beta_0 < 0$。拿這個 $\beta_0$ 去組合**同一批**模型，等於「讓好的更好、讓壞的更壞」，是自我強化的偏誤。
  - $\beta$ 只能用「預報視窗已結束」的過去資料估計，天生有延遲。短期模型在序列期間預報變化極快，延遲的 $\beta_0$ 會嚴重失準。
  - 目標事件極稀（本例約 0.006%），稀有事件本來就會讓截距有偏（King & Zeng, 2001）。

  **解法：只取相對資訊，丟掉截距**。把係數映射成權重：

  $$w_i = \begin{cases} e^{\beta_i} - e^{\tau}, & \beta_i > \tau \\ 0, & \text{otherwise}\end{cases},\qquad \pi_i = \frac{w_i}{\sum_j w_j}$$

  取 $\tau = 0$：$\beta_i \le 0$ 的模型（無法解釋觀測）直接給零權重。用指數映射的理由是 $e^{\beta_i}$ 恰好是該模型的 **odds ratio**（$u_i$ 每增加一單位，勝算乘以 $e^{\beta_i}$）；減去 $e^\tau$ 讓權重平滑歸零。最後用這組 $\pi_i$ 回頭做**加權平均** ensemble $\phi_E^{\mathrm{WA}}$，完全擺脫對 $\beta_0$ 的依賴。

  **評分**：相對參考模型的每事件資訊增益（IGPE, Rhoades et al. 2011）

  $$I_{A,\mathrm{ref}} = \frac{1}{N}\sum_{k=1}^{N}\left(\ln \lambda^A_{jk} - \ln \lambda^{\mathrm{ref}}_{jk}\right) - \frac{\hat N_A - \hat N_{\mathrm{ref}}}{N}$$

  並逐時間累積成 CumIGPE，用 t-test 給 95% 信賴區間。

- **實作規模與量化結果**
  - OEF-Italy：3 個候選模型（ETAS-LM、ETES-FCM、STEP-LG），$M_L \ge 3.95$，$0.1° \times 0.1°$ 共 8993 格，2005–2020 共 6227 個重疊預報視窗 → 每模型約 5600 萬個樣本，其中 99.994% 是無事件格。用 scikit-learn liblinear、正則化強度壓到 $10^{-8}$（等於不懲罰係數），並把無事件格**降採樣到 10%**（截距的偏誤事後校正）。
  - 兩種擬合方案：**#1 用全部歷史資料**（權重隨時間收斂）；**#2 只用最近一年**（權重會「遺忘」，反映近期表現）。
  - 相對 SMA ensemble 的 CumIGPE：$\phi_E^{\mathrm{WA}}$ 為 $0.064 \pm 0.017$（方案 #1）、$0.078 \pm 0.023$（方案 #2）；若改用 $M_L \ge 2.95$ 的事件來擬合權重（但預報門檻仍是 3.95），可達 $0.112 \pm 0.019$。
  - 相對**最佳單一模型 ETAS-LM**：0.002 ± 0.019（#1）、0.016 ± 0.028（#2）皆不顯著；但用 $M_L \ge 2.95$ 擬合時達 $0.050 \pm 0.019$，**顯著勝過最佳單一模型**。
  - 直接用 logistic 模型當 ensemble（$\phi_E^{\mathrm{logistic}}$）表現是**負的**（−0.06 至 −0.10），比 SMA 還差——印證了上面對 $\beta_0$ 的診斷。
  - 權重的時序本身就有解讀價值：L'Aquila（2009）與中義大利（2016）序列期間 ETAS-LM 權重飆高，Emilia（2012）期間 ETES-FCM 較佳，STEP-LG 常常拿到接近零的權重。

- **教學上可用的洞見**
  - **一句話重點**：「依個別表現給權重」與「最大化整體表現給權重」是兩件事，後者才是對的。前者（SMA）會讓表現差的模型還是佔到份額，拖累整體。
  - **遺忘機制**：只用近一年資料讓 ensemble 更「適應性」，能在序列爆發時快速把權重轉給當下最強的模型。作者實測不同視窗長度：1.5 年最佳（0.097）、2 年 0.074、半年 0.049（用 $M_L \ge 3.95$ 擬合時）。但這是可調的自由度，也是過度配適的風險來源——文中多項改良「單獨用都有效、一起用就沒效」。
  - **用小地震學權重**：預報目標是 $M \ge 3.95$，但可以用 $M \ge 2.95$ 的事件來擬合權重，樣本量大增、權重更靈敏。這是很實用的技巧（本書談 b 值、完整性規模時可以呼應）。
  - **ensemble 的天花板來自候選池**：作者自評增益「相對溫和（IGPE ~0.05）」，並直指原因是三個候選模型太像——都是統計式的群聚模型、背景率決定法也雷同。**要提升 ensemble，必須增加候選模型的多樣性**（統計 + 物理、以及在特定情境下特別強的模型）。這和 2014 那篇「概念差異愈大、增益愈大」完全一致，是兩篇跨越九年的共同結論。
  - **權重的副產品**：加權平均保留了各成分預報的**離散度**，可以用加權變異數估計知識論不確定性（epistemic uncertainty）；而乘法／collapse 型的組合直接把所有預報壓成單一分布，這個資訊就沒了。這是選擇加法而非乘法的一個常被忽略的理由。
  - **實務可靠性論證**：加權平均 ensemble「不一定大幅贏過最佳模型，但從來沒輸過」。而實驗開始時你根本不知道哪個模型會最好——所以 ensemble 是**先驗上的理性選擇**。這句話適合當整章的結論。

- **與台灣的關聯**：台灣若建置作業型地震預報（OEF），這篇提供的是「多模型如何合成單一權威預報」的完整工程流程：權重擬合、稀有事件偏誤校正、降採樣、時間視窗選擇、以及用 pyCSEP 計算 IGPE。方法已釋出為 Python class（Zenodo, doi:10.5281/zenodo.7477998），並規劃併入 pyCSEP。

---

## 跨篇綜合：這個主題教什麼

### 一、為何 ensemble 幾乎總是贏

三個層次的理由，建議按此順序講：

1. **統計層次——分散化**。不同模型抓的是資料裡不同的訊號。STEP／ETAS 抓短期餘震群聚，EEPAS 抓中期前兆尺度增長，PPE／HKJ 抓長期空間分布，大地測量模型抓應變累積。它們的預報率在同一格可以差到 12 個數量級（2009）。就像投資組合，**低相關資產的組合能提升風險調整後報酬**——這裡的「報酬」就是資訊增益。
2. **決策層次——先驗的理性**。實驗開始時沒人知道哪個模型會最好。加權平均 ensemble「不保證大贏，但從不落後」（2023），所以它是**不需要事先押注**的選擇。這對作業型預報（OEF）尤其重要：需要一個權威預報，又不想主觀挑一個模型。
3. **知識論層次——保留不確定性**。加權平均保留了各模型的離散度，可以量化 epistemic uncertainty（2023）。

**但要誠實地加上但書**：ensemble 的天花板來自候選池。三個高度相似的統計模型，最多也只換到 IGPE ~0.05（2023）；而概念差異大的組合（地震目錄 × 大地測量）在回溯上能到 0.25–0.79（2014）。**多樣性才是 ensemble 的燃料。**

### 二、混合（加法）vs 加乘（乘法）怎麼教

| | **加法 / mixture** | **乘法 / multiplicative hybrid** |
|---|---|---|
| 形式 | $\lambda_H = \sum_i \pi_i \lambda_i$ | $\lambda_H = \lambda_1 \exp[a + \sum_i f_i(\lambda_i)]$ |
| 角色 | 各模型平等，權重總和為 1 | 一個 baseline + 若干修正因子 |
| 值域 | 只能**內插**：落在各成分之間 | 可以**外推**：可高於或低於所有成分 |
| 正規化 | 凸組合自動守恆 | 需要正規化參數 $a$ |
| 規模分布 | 各模型各自負責 | 由 baseline 全權負責（conjugate 先對規模求和） |
| 參數量 | 每模型 1 個 | 每 conjugate 2 個 + 全域 1 個 |
| 失敗模式 | 被稀釋，最差不過是平庸 | 乘數 < 1 時把率壓得比誰都低，可以錯得更離譜 |
| 實測 | RELM 加法組合打不贏最佳單一模型（Marzocchi 2012）；但 OEF-Italy 加權平均打贏了（2023，關鍵在**權重怎麼學**） |回溯上大勝加法（$\Delta\ln L$ 11.4 vs 4.8, 2014）；前瞻上全數不顯著（2022） |

**教學動線建議**（也是這四篇天然的敘事順序）：

1. **2009**：先從最簡單的凸組合入手。只擬合一個權重 $r$，就拿到相對兩個母模型各超過 2 倍的機率增益。建立「互補性 → 增益」的直覺。
2. **2014**：加法把增益「稀釋」掉了（EEPAS 對 PPE 有 0.5 的優勢，但 SE1 對 SP1 只剩 0.13）。乘法理論（Aki、Utsu、Imoto 的條件獨立前兆）暗示可以保留完整增益。於是有了 multiplicative hybrid，以及「加法只能內插、乘法可以外推」這句核心台詞。
3. **2022**：開獎。前瞻測試全數失敗。教學上這比「乘法很棒」有價值得多——它教的是**方法論的誠實**：怎麼診斷失敗（四個原因）、以及乘法結構本身如何在 baseline 空間表現不佳時放大錯誤。
4. **2023**：回到加法，但把重點從「用什麼運算子」轉移到「權重怎麼學」。用 logistic regression 最大化 ensemble 本身的表現，而不是各模型的個別表現；並且示範了「該用的是係數，不是模型本身」這種微妙但關鍵的實作智慧。

一個好用的收束：**運算子（加法 vs 乘法）的選擇沒有普世答案，但「權重要用資料學、而且要用最大化整體表現的方式學」是共識。** 2022 的結尾也給出後續建議的研究方向：**相關的模型可試加法組合、獨立的模型可試乘法組合**（作者措辭是「值得進一步研究」，不是已確立的經驗法則）。

### 三、與第 15 章（模型檢驗）的關係

Ensemble 這一章和檢驗那一章是互相咬合的，四篇論文提供了三個接點：

1. **概似（likelihood）既是擬合目標，也是評分標準**。權重用最大概似擬合，資訊增益也用對數概似算。所以講 ensemble 一定要先把 $\ln L = \sum_i \ln \lambda_i - \sum_j \lambda_j$ 講透。這也解釋了為什麼 ensemble 幾乎「不可能變差」——最佳化本身保證 $\ln L$ 不下降（2014 明確指出 hybrid 退化成 baseline 的條件），所以**必須**用 AIC/AICc 扣參數，否則所有比較都是自欺。

2. **回溯 vs 前瞻的斷裂**。IGPEc（含參數懲罰）與 IGPE（無懲罰、用於獨立資料）的差別，正是這條斷裂線的形式化。2014 與 2022 這一組是最佳教材：同一批模型、同一個測試區、同一個統計量，回溯 +0.25～+0.5，前瞻 −0.42～−0.71。**沒有前瞻測試就沒有結論。**

3. **檢驗的假設會反噬**。2022 引進 binary likelihood 是因為 Poisson 概似在群震序列（Hawthorne 群震、Ridgecrest 序列）面前失真——少數幾個多震格就吃掉近 40% 的懲罰。同理，negative-binomial N-test 處理總數的過度離散。可以在第 15 章展開 N-test / S-test / cL-test / T-test 的分工，並用「同一模型在 POLL 與 BILL 下分數差多少」當作實作練習。台灣的地震序列群聚性極強，這個議題只會更尖銳。

4. **有效樣本數的陷阱**。格子有數萬個，但擬合權重時的有效樣本數是**目標地震顆數**（2014 只有 31 個，2009 是 152 個）。2022 診斷失敗的第一個原因就是「31 個事件擬合出來的權重時間上不穩定」。這一點應該在兩章都強調。

---

## 一句話速記

- **2009 Rhoades & Gerstenberger**：短期（STEP）+ 中期（EEPAS）的凸組合，只多一個權重參數，機率增益翻倍以上。
- **2014 Rhoades et al.**：baseline × 保序轉換後的修正因子 = 乘法 hybrid；回溯上大勝加法，且模型概念差異愈大增益愈大。
- **2022 Bayona et al.**：十年前瞻測試開獎——乘法 hybrid 全數未顯著勝過 HKJ baseline，回溯增益不可信，乘法結構會放大空間預報的錯誤。
- **2023 Herrmann & Marzocchi**：權重應該最大化 ensemble 的表現而非各模型的個別表現；用 logistic regression 的係數（而非 logistic 模型本身）做加權平均，在 OEF-Italy 顯著勝過最佳單一模型。
