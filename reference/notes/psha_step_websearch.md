# PSHA 與 STEP 補充素材（WebSearch）

> 來源：網路搜尋（本檔所有書目均經搜尋核實）| 供第二部 13 章（STEP/OEF）與 16 章（PSHA）使用
> 搜尋日期：2026-08-26

**核實規則**：以下每一條書目都經由 WebFetch 抓取出版社／期刊／作者官方頁面確認，或由多筆搜尋結果交叉比對。凡是「只在搜尋摘要出現、未能抓到原始頁面」的內容，都在該段落明確標註為「未核實」。

---

## A. PSHA

### Baker 白皮書

**書目（PDF 內頁「Preferred citation」逐字抄錄）**

> Baker, Jack W. (2013) *Probabilistic Seismic Hazard Analysis*. White Paper Version 2.0.1, 79 pp.

- 下載連結：<https://www.jackwbaker.com/Publications/Baker_(2013)_Intro_to_PSHA_v2.pdf>（Stanford 鏡像：`https://web.stanford.edu/~bakerjw/Publications/Baker_(2013)_Intro_to_PSHA_v2.pdf`）
- 舊版 v1.3（2008）仍可下載：`https://www.jackwbaker.com/Publications/Baker_(2008)_Intro_to_PSHA_v1_3.pdf`
- **版本號注意**：搜尋摘要常寫成「version 2.1」，但 PDF 版權頁自述為 **2.0.1**。以 PDF 內文為準。
- **文件已被教科書取代**，PDF 第一頁明寫：
  > "This document has been superseded by a new textbook: Baker, J. W., Bradley, B. A., and Stafford, P. J. (2021). *Seismic Hazard and Risk Analysis*. Cambridge University Press, Cambridge, England."
  > 配套資源網站：<http://www.pshabook.com/>
- 授權：Apache License 2.0（教學重製友善，可安心在網站上引用圖表結構，仍須註明出處）

**章節架構（目錄逐字）**

| 章 | 標題 | 起始頁 |
|---|---|---|
| 1 | Introduction | 13 |
| 2 | An overview of PSHA | 15 |
| 2.1 | Deterministic versus probabilistic approaches | 15 |
| 2.2 | Probabilistic seismic hazard analysis calculations | 18 |
| 2.3 | Example PSHA calculations | 33 |
| 3 | Extensions of PSHA | 43 |
| 3.1 | Deaggregation | 43 |
| 3.2 | Bounds on considered magnitudes and distances | 48 |
| 3.3 | Rates, probabilities and return periods | 49 |
| 3.4 | The uniform hazard spectrum | 52 |
| 3.5 | Joint distributions of two intensity measures | 52 |
| 4 | Conclusions | 59 |
| A | Review of probability | 61 |
| B | Further study | 75 |

**PSHA 框架：Baker 原文寫的是「五步驟」不是四步驟**

原文（2.2 節）逐字：

> "At its most basic level, PSHA is composed of five steps.
> 1. Identify all earthquake sources capable of producing damaging ground motions.
> 2. Characterize the distribution of earthquake magnitudes (the rates at which earthquakes of various magnitudes are expected to occur).
> 3. Characterize the distribution of source-to-site distances associated with potential earthquakes.
> 4. Predict the resulting distribution of ground motion intensity as a function of earthquake magnitude, distance, etc.
> 5. Combine uncertainties in earthquake size, location and ground motion intensity, using a calculation known as the total probability theorem."

常見的「四步驟」講法（震源描述 → 規模-頻率分布 → GMPE → 危害積分）是把上述第 1、3 步合併成「震源幾何與距離分布」。**寫章節時若要用四步驟，須說明這是常見的濃縮版，不要寫成「Baker 說有四步」。**

各步驟關鍵內容（原文重點）：

- **步驟 1 震源**：斷層源（planar surfaces，由歷史地震位置與地質證據界定）；若無法辨識個別斷層（如美東），則用**面震源**（areal source，地震可在區域內任意位置發生）。
- **步驟 2 規模分布**：Gutenberg–Richter 遞減律 `log10 λm = a − b·m`，λm 為規模大於 m 的地震發生率。
- **步驟 4 地動預測**：以對數常態分布描述，`P(PGA > x | m, r) = 1 − Φ((ln x − ln PGA_bar)/σ_lnPGA)`。文中示範用 Cornell et al. 模型，M6.5 在 3/10/30 km 的 ln PGA 平均值為 −0.5765 / −0.9788 / −1.7937，σ_lnPGA = 0.57。**這組數字很適合當教學範例**（不確定性的量級：一個標準差在 PGA 上就是約 1.8 倍）。
- **步驟 5 全機率定理積分**：即白皮書的 equation 2.25（危害積分主方程）。

**Hazard curve（危害曲線）**

- 縱軸為「年超越率」（annual rate of exceedance，通常 10⁻¹ ~ 10⁻⁵ 對數尺度），橫軸為地動強度（PGA, g）。
- 白皮書有一組很好的教學圖：同一場址、只改 `m_min`（4.5 / 5.0 / 5.5）三條 hazard curve 疊圖，示範**小規模地震截斷門檻對高強度端幾乎無影響**。文中並說明實務上通常只考慮 M > 4.5~5.0，因為被略去的小地震「不被認為有破壞結構的能力」。

**回歸期 vs 超越機率（3.3 節，教學上最容易講錯的一段）**

原文逐字重點：

> "The return period is defined as the reciprocal of the rate of occurrence. ... if a given ground motion intensity has a 0.01 annual rate of occurrence, then the return period is equal to 1/0.01 = 100 years. This does not imply that the ground motion will be exceeded exactly once every 100 years, but rather that the average (or mean) time between exceedances is 100 years. For this reason, the reciprocal of the exceedance rate is more precisely termed the **mean return period**."

Baker 甚至建議：

> "...one may avoid some confusion regarding the implied time between exceedances by simply reporting rates rather than return periods."

- 由發生率換算成時窗機率需要「地震間隔的機率分布」假設，**幾乎一律假設 Poisson**，理由有三：數學式簡單、多數情況符合觀測、更複雜的模型通常不顯著改變最終結果。
- Poisson 假設下：`P(至少一次事件於時間 t 內) = 1 − exp(−λt)`。
- 這正是「50 年 10% 超越機率 ≈ 475 年回歸期」這個工程慣用數字的來源（`λ = −ln(0.9)/50 ≈ 0.0021/yr`，倒數約 475 年）。

**Deaggregation（去聚合／解積，3.1 節）**

- 動機原文：PSHA 把所有情境加總是優點也是缺點——算完後「哪一種地震情境最可能造成 PGA > x？」答案並不直觀。
- 定義（equation 3.1）：`P(M = m | IM > x) = λ(IM > x, M = m) / λ(IM > x)`，即 Bayes 定理的應用；分母就是 PSHA 主方程，分子只是「不對 M 積分」而已。可推廣成 M、R 的聯合去聚合（equation 3.8/3.9），再加上 ε（epsilon，地動殘差標準差數）成為三維去聚合。
- **教學上非常好用的定性結論**（白皮書 Example 1）：同一場址兩條斷層，在較低的 PGA = 0.376 g 時，`P(M=6.5) = 0.77`、`P(M=7.5) = 0.23`（近而活躍的小斷層主導）；到 PGA = 1 g 時變成 0.58 / 0.42（遠而大的斷層貢獻上升）。**「危害水準越高，主控震源越大越遠」**這個現象一句話就能講清楚 deaggregation 的價值。
- 白皮書引用的 USGS 去聚合範例（Palo Alto, SA(1.0s) > 0.6385 g, Mean Return Time = 2475 yrs, Mean (R, M, ε₀) = 10.3 km, 7.65, 1.25）是很典型的去聚合圖版面，可作為畫圖參考。
- 術語小註（原文腳註）：deaggregation 與 disaggregation 兩詞並用，尚無統一；disaggregation 才是字典收錄的字，但 deaggregation 目前用得比較多。
- 出處文獻：Bazzurro & Cornell (1999)、McGuire (1995)（白皮書引用，本次未另行核實其卷期）。

**延伸概念**：3.4 節的 uniform hazard spectrum（等危害度反應譜）、3.5 節兩個強度指標的聯合分布（PGA 與 SA 的相關係數 ρ ≈ 0.7，見 Baker & Cornell 2006、Baker & Jayaram 2008），可作為進階選讀。

---

### 台灣地震模型 TEM

**TEM PSHA2015**（核實：TAO 官方頁面）

> Wang, Y.-J., C.-H. Chan, Y.-T. Lee, K.-F. Ma, J. B. H. Shyu, R.-J. Rau, and C.-T. Cheng (2016). Probabilistic Seismic Hazard Assessment for Taiwan. *Terrestrial, Atmospheric and Oceanic Sciences (TAO)*, **27**(3), 325–340. doi: 10.3319/TAO.2016.05.03.01(TEM)

模型要點（摘自官方摘要，屬公開發表的定性描述）：

- 採用 TEM 地質學家辨識的 **38 條孕震構造**（seismogenic structures）之震源參數。
- 除斷層源分類外，地震活動另分為 **淺層（shallow）、隱沒帶板塊內（subduction intraplate）、隱沒帶板塊間（subduction interplate）** 三類。
- 對地殼型與隱沒帶型地震分別採用對應的 GMPE。
- 危害圖概貌（原文結論）：**最高危害機率位於台灣西南部與東部縱谷**；西部人口稠密的直轄市中以 **台中、台南、新北** 危害最高；以 PGA 而言 **台南最高**；以擬譜加速度而言，**台南在短週期較高、台中在長週期較高**，因此原文指出「台南低樓層建物、台中高樓層建物的耐震設計特別重要」。
- 後續文獻補充（核實：Chan et al. 2020 開放版 PDF 引言逐字）：TEM PSHA2015「adopted **38 seismogenic structures** identified by Shyu et al. (2016), **28 shallow-area sources**, **4 subduction-interplate sources**, and **12 subduction-intraplate sources**」；且只評估至**工程岩盤（Vs30 = 760 m/s）**，「therefore site amplification was neglected」。
- 補充：TEM 構造資料庫在 2020 年版已擴充到 **44 條孕震構造**（PSHA2015 當時只用了其中 38 條）。

**TEM PSHA2020**（核實：SAGE / Earthquake Spectra 期刊頁面）

> Chan, C.-H., K.-F. Ma, J. B. H. Shyu, Y.-T. Lee, Y.-J. Wang, J.-C. Gao, Y.-T. Yen, and R.-J. Rau (2020). Probabilistic seismic hazard assessment for Taiwan: TEM PSHA2020. *Earthquake Spectra*, **36**(1_suppl), 137–159. doi: 10.1177/8755293020951587

相對 2015 版的更新（摘自官方摘要）：

- 更新孕震構造資料庫，納入**新辨識出的三維幾何構造**。
- 地震目錄更新至 **2016 年**。
- 新的一組 GMPE，並**加入場址放大係數**（2015 版沒有）。
- 新增**多構造同時破裂**（multi-fault rupture）的可能性。
- 對孕震構造源導入 **Brownian passage time (BPT) 模型** 以描述「斷層記憶」——**這是 16 章接到 13 章時變觀念的天然橋樑**。
- 無法歸屬到特定構造的地殼地震，同時使用**面震源**與**平滑核（smoothing kernel）**兩種模型。

**TEM PSHA2025（僅為預印本，勿當期刊論文引用）**

搜尋到 SSRN 預印本：Gao, J.-C., J.-C. Kao, C.-H. Chan, R. Y. Chuang, C. H. Chen, J. B. H. Shyu, K.-E. Ching, Y. Wang, K.-F. Ma, "Probabilistic Seismic Hazard Assessment for Taiwan: Updates and Improvements in TEM PSHA2025", SSRN abstract id 5736969。**尚未核實其正式期刊出版狀態，寫章節時只能以「最新一代模型正在發展中」帶過，不要給卷期頁碼。**

其他可用連結：TEM 官網 <https://tem.tw/>（TEM2020 專頁 <https://tem.tw/TEM2020/portfolio-overview.html>）；GEM 基金會 Taiwan Hazard 產品頁 <https://www.globalquakemodel.org/product/taiwan-hazard>。

---

### 時變危害（time-dependent hazard / OEF 接 PSHA）

1. **ICEF 報告（OEF 的定調文件，核實：Annals of Geophysics 期刊頁面）**

   > Jordan, T. H., Y.-T. Chen, P. Gasparini, R. Madariaga, I. Main, W. Marzocchi, G. Papadopoulos, G. Sobolev, K. Yamaoka, and J. Zschau (2011). Operational Earthquake Forecasting: State of Knowledge and Guidelines for Utilization. *Annals of Geophysics*, **54**(4). doi: 10.4401/ag-5350（線上出版日 2011-08-03）

   關鍵論點：報告區分 **prediction（確定性陳述）** 與 **forecast（機率性評估）**；結論是「診斷性前兆的搜尋尚未產生成功的短期預測方案」，因此主張以 OEF 作為溝通地震風險的主要工具。並明確要求 OEF 必須「與 PSHA 的長期預報一致地」提供完整的危害描述——**地動超越機率**，而不只是短期破裂機率。這句話正是 13 章接 16 章的論證主軸。

2. **PSHA 現況回顧（含時變議題）**

   > Gerstenberger, M. C., W. Marzocchi, T. Allen, M. Pagani, J. Adams, L. Danciu, et al. (2020). Probabilistic Seismic Hazard Analysis at Regional and National Scales: State of the Art and Future Challenges. *Reviews of Geophysics*, **58**, e2019RG000653. doi: 10.1029/2019RG000653
   >
   > （卷號與文章編號來自搜尋結果交叉比對，未逐字抓到期刊頁面；引用前建議再核一次。）

3. **OEF 全面回顧（reference/ 庫內已有：`[2024] OEF_Review.pdf` 極可能即為此篇）**

   > Mizrahi, L., I. Dallo, N. J. van der Elst, A. Christophersen, I. Spassiani, M. Werner, P. Iturrieta, J. A. Bayona, I. Iervolino, M. Schneider, M. T. Page, J. Zhuang, M. Herrmann, A. J. Michael, G. Falcone, W. Marzocchi, D. Rhoades, M. Gerstenberger, L. Gulia, D. Schorlemmer, J. Becker, M. Han, L. Kuratle, M. Marti, and S. Wiemer (2024). Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions. *Reviews of Geophysics*, **62**(3), e2023RG000823. doi: 10.1029/2023RG000823

4. **台灣本地的時變切入點**：TEM PSHA2020 的 BPT 斷層記憶模型（見上）即為時變危害在台灣的官方實作，不必另找外國案例。

---

## B. STEP

### 書目核實

**Gerstenberger et al. 2005（核實：Nature 官方頁面，抓到摘要全文）**

> Gerstenberger, M. C., S. Wiemer, L. M. Jones, and P. A. Reasenberg (2005). Real-time forecasts of tomorrow's earthquakes in California. *Nature*, **435**, 328–331. doi: 10.1038/nature03622（2005-05-19 出版；PMID 15902254）

摘要重點（原文逐字節錄）：

> "Our model builds upon the basic principles of this generic forecast model in two ways: it recasts the forecast in terms of the **probability of strong ground shaking**, and it combines an existing **time-independent earthquake occurrence model based on fault data and historical earthquakes** with increasingly complex models describing the **local time-dependent earthquake clustering**. The result is a **time-dependent map showing the probability of strong shaking anywhere in California within the next 24 hours**."

摘要中也點出既有「通用（generic）短期叢集模型」的兩個限制：**未針對當下序列調校**、**不含餘震可能位置的資訊**——STEP 就是為了補這兩點而生。這兩句話是 13 章切入 STEP 動機最好的引言。

原論文的即時系統網址為 `pasadena.wr.usgs.gov/step`（已停止服務，僅供歷史說明）。

**Reasenberg & Jones 1989（核實：Science 官方 DOI 頁面 + Google Scholar 條目）**

> Reasenberg, P. A., and L. M. Jones (1989). Earthquake Hazard After a Mainshock in California. *Science*, **243**(4895), 1173–1176. doi: 10.1126/science.243.4895.1173（PMID 17799897）

- 後續更正／更新：Reasenberg, P. A., and L. M. Jones (1994). Earthquake Aftershocks: Update. *Science*, **265**(5176), 1251–. doi: 10.1126/science.265.5176.1251
  （**部分核實**：標題、卷 265、期 5176、起始頁 1251 由 DOI 字串本身確定；**結束頁未核實**，science.org 頁面回 403 無法抓取。引用時只寫起始頁，或先自行確認頁碼範圍。）

**Reasenberg–Jones 公式與通用參數（核實：GNS/NZ 官方報告 PDF 原文）**

R–J 把修正大森律（Omori–Utsu）的生產力參數 K 換成主震規模的函數，得到主震 M_m 之後時刻 t、規模 M 以上的餘震發生率：

```
R(t, M) = 10^(a' + b(M_m − M)) · (t + c)^(−p)
```

- 加州通用參數（分析 **62 個加州地震序列** 得到）：`a' = −1.67, p = 1.08, b = 0.91, c = 0.05`。
- 時窗機率：`P = 1 − exp(−∫R dt)`，此式給出「一個地震在一週內被同規模或更大地震跟隨」的機率 **10.5%**；但南加州實際只有 **6.0 ± 0.5%** 的地震是前震（Jones 1985）——**這個「模型高估、觀測偏低」的對照組是講模型驗證的好例子**。
- Utsu (1969) 對日本地震用同一式，參數 `a' = −1.83, p = 1.3, b = 0.85, c = 0.3`，同樣算法得到 **4.2%**——**同一個模型換區域參數，答案差一倍以上，正好講「參數在地化」的必要性。**

（以上出自 GNS Science / EQC 報告：*Towards a New Zealand model for short-term earthquake probability*, Project No. 6OPR1B，<https://www.naturalhazards.govt.nz/assets/Publications-Resources/2311-Towards-a-New-Zealand-model-for-short-term-earthquake-probability.pdf>。報告本身未核實作者全名與正式出版年，引用時以機構報告形式標註較安全。）

---

### 模型結構

**三層元素 + AIC 權重組合（部分核實，見下方標註）**

已抓到的**逐字**描述來自 Steacy et al. (2014)：

> "The STEP model (Gerstenberger et al. 2005) is also based on superimposed Omori sequences, where every earthquake is allowed to generate its own aftershock sequence. A forecast is created using an **ensemble of three different models of increasing complexity, which are combined using weights calculated from the Akaike Information Criterion** (Burnham & Anderson 2002)."
>
> 出處：Steacy, S., M. Gerstenberger, C. Williams, D. Rhoades, and A. Christophersen (2014). A new hybrid Coulomb/statistical model for forecasting aftershock rates. *Geophysical Journal International*, **196**(2), 918–923. doi: 10.1093/gji/ggt404

三層元素的**觸發條件**（核實：GNS/NZ 報告原文，該段直接引用 Gerstenberger 2003/2004/2005）：

1. **generic（通用層）**——序列剛開始、資料不足時，直接用加州通用參數 `a' = −1.67, p = 1.08, b = 0.91`。空間上先用 Gerstenberger (2003) 的三步法：以主震為中心、半徑 `r(M) = 10^(0.59M − 2.44)` km（由 Wells & Coppersmith 1994 的次表面破裂長度尺度律導出，下限至少 5 km）的圓形區域向外遞減；有餘震資料後再據以估計斷層跡並修正餘震分布區。
2. **sequence-specific（序列特定層）**——「一旦記錄到**完整度規模以上至少 100 個餘震**，就計算該序列自身的模型參數」。
3. **spatially-varying（空間變化層）**——「對於生產力特別高的序列，參數容許**在空間上變化**」（Gerstenberger et al., 2005）。Woessner et al. (2011) 對 STEP-1 的描述同義：參數「先對整個序列估計，一旦序列地震量足以解析空間變異，就允許隨空間變化」。原始 STEP 的空間網格為 5 km × 5 km（此格距僅在搜尋摘要層級看到，未抓到原始出處，寫章節時可省略具體數字）。

另外，STEP 的預報 = 上述**時變叢集項** + 一個**時間獨立背景項**（由去除叢集後的目錄平滑而得）。Woessner et al. (2010, Ann. Geophys. **53**(3), 141–154, doi: 10.4401/ag-4812；作者序 Woessner, Christophersen, Zechar, Monelli）義大利版 STEP 的摘要即明說「兩種實作都結合一個時變貢獻與一個時間不變貢獻」，並示範用 STEP 自身的叢集演算法去叢集後再平滑建背景。

> **未核實項**：AIC 權重的**精確公式**、以及「generic / sequence-specific / spatially-varying」這組**英文標籤是否為 Gerstenberger et al. 2005 原文用語**，都沒有抓到原始出處（Nature 2005 的 Supplementary Methods 在本次搜尋中無法取得）。教學時建議寫成「三個複雜度遞增的模型元素，以 AIC 權重加權組合」（有 Steacy 2014 逐字支持），元素名稱可以中文描述其**觸發條件**（資料不足 → 通用參數；累積足夠餘震 → 序列自身參數；序列夠豐富 → 允許參數空間變化），這部分有 GNS 報告支持。

---

### 現行作業系統範例

**USGS Operational Aftershock Forecasting（OAF）**

- 產品說明頁：<https://earthquake.usgs.gov/data/oaf/overview.php>
- 軟體頁（開源，CC0，作者 Michael Barall、Nicholas J. van der Elst）：<https://www.usgs.gov/software/operational-aftershock-forecasting>
- **模型**：Reasenberg–Jones (1989, 1994) 或 ETAS (Ogata 1988)；參數有三種來源——**generic**（由過去相似序列導出）、**sequence-specific**（用本次序列資料擬合）、**Bayesian**（兩者結合）。**注意：現行 USGS 系統已不是 2005 年的 STEP，13 章要講清楚這個世代差異。**
- **呈現方式**（適合直接當教學截圖範例）：
  - 出現在個別地震事件頁上的一張「Aftershock Forecast」卡片，點入後有四個分頁：**Summary / Commentary / Forecast Table / Model Parameters**。
  - 互動圖形：使用者選規模級距（M3+, M4+, M5+, M6+, M7+）與時窗（day / week / month / year），以填色弧形顯示機率，並附餘震數量的機率分布。
  - 表格：呈現「某規模以上、某時窗內至少發生一次的機率」，或在預期事件數較多時直接給數量。
  - 文字摘要：非技術性語言的說明。
- **發布規則**：美國本土 M4+、其他州與屬地 M5+ 的地震會發布預報；系統在雲端持續監看 ComCat 目錄，多數地震約 **20 分鐘**內自動發出首報，第一年內更新 **75 次**（此兩項數字來自 USGS 頁面搜尋摘要，未逐字核實）。

**GNS Science / GeoNet（紐西蘭）**

- 總覽頁：<https://www.gns.cri.nz/our-science/natural-hazards-and-risks/earthquakes/earthquake-forecasting/>
- 預報入口：<https://www.geonet.org.nz/earthquake/forecast/>（分區：Canterbury、Central New Zealand、Kaikōura）
- **呈現方式（核實：Canterbury 分頁）**：一張表，欄位為三個規模級距 **M5.0–5.9 / M6.0–6.9 / M≥7.0**，每格給三個數字——**平均預期事件數、可能範圍（min–max）、發生一次以上的機率**。時窗為自發布日起**一年**。範例列（2025-02-01 版）：M5.0–5.9 平均 0.4、範圍 0–2、機率 29%；M6.0–6.9 平均 0.04、範圍 0–1、機率 4%。
- **語言處理值得學習**：機率旁附文字判讀，如「Within the next year, there is a 29% probability (**unlikely**) of one or more earthquakes of magnitude 5.0 to 5.9 occurring」、「It is **extremely unlikely** (less than 1%) that there will be an earthquake of magnitude 7 or greater」。同時明確標出涵蓋範圍的經緯度（171.6–173.2°E, 43.3–43.9°S）並附地圖。
- 紐西蘭自 2010 年 Darfield 地震起即公開發布地震預報；現行為混合式預報工具（HFT），結合三種時間尺度的模型：短期 ETAS 類叢集、中期（年至十年）叢集、長期（十年以上至近似時間獨立）。

---

## 教學建議

**16 章（PSHA）**

1. **主軸用 Baker 白皮書的骨架**：它是 Apache 2.0 授權、79 頁、有完整算例的教學文件，比任何教科書都適合當網站的觀念級底稿。切記標註「已被 Baker, Bradley & Stafford (2021) 教科書取代」，並給 pshabook.com 當延伸閱讀。
2. **步驟數不要寫死成「四步」**。用 Baker 的五步呈現，再加一句「文獻中常把識別震源與距離分布併為一步，故也常見四步驟的說法」。
3. **回歸期是最容易講錯的觀念**，建議照 Baker 的處理：先講「回歸期 = 發生率倒數」，立刻強調它是**平均**間隔而非週期，再用 Poisson 的 `1 − e^(−λt)` 導出「50 年 10%↔475 年」。Baker 本人建議直接報發生率避免混淆——這個「作者的忠告」本身就是很好的教學橋段。
4. **Deaggregation 用「PGA 從 0.376 g 提高到 1 g，主控震源從 M6.5 換成 M7.5」這個算例**，一張圖一句話講完概念，不必推導。
5. **台灣落地就用 TEM**：2015 版給「38 條孕震構造 + 分三類地震活動 + 西南部與縱谷危害最高、台南短週期／台中長週期」的定性描述；2020 版給「3D 幾何、多構造破裂、BPT 斷層記憶、面震源＋平滑核、加入場址放大」。**BPT 是把 16 章接回 13 章時變觀念的最短路徑。**
6. **避開的內容**：TEM PSHA2025 只有 SSRN 預印本，不要給卷期；任何具體危害圖數值（PGA 幾 g、幾年回歸期的分區數字）本次都沒有核實到，網站上只放定性描述或直接連到 TEM 官網圖資。

**13 章（STEP / OEF）**

1. **敘事線**：Reasenberg & Jones (1989) 的通用加州參數 → 它的限制（不針對當下序列、沒有位置資訊，Nature 2005 摘要親口點名）→ STEP 的三層遞增複雜度 + AIC 權重 → 輸出從「餘震數量機率」升級為「24 小時內強震動機率地圖」。
2. **兩個現成的數值對照組非常好用**：加州通用參數算出「一週內被更大地震跟隨」10.5% vs 南加州實測前震率 6.0%（模型 vs 觀測）；加州參數 10.5% vs 日本 Utsu 參數 4.2%（區域 vs 區域）。兩者都只需一行算式。
3. **AIC 權重只講觀念**：「三個複雜度遞增的模型，用 AIC 依資料量與擬合度自動決定各自權重——資料少時通用模型佔優，資料累積後序列特定與空間變化模型接手」。**不要寫出 AIC 權重的具體公式**，本次沒核實到原始定義。
4. **元素名稱處理**：中文可寫「通用層 / 序列特定層 / 空間變化層」，並用觸發條件（無資料 → ≥100 個完整度以上餘震 → 序列生產力足夠）說明切換邏輯，這部分有 GNS 報告支持。**不要宣稱這三個英文標籤出自 Gerstenberger et al. 2005 原文。**
5. **作業系統範例建議並列 USGS 與 GeoNet**：USGS 展示「規模級距 × 時窗」的互動介面與四分頁結構；GeoNet 展示「平均數 + 範圍 + 機率」三欄表格與**機率的文字判讀語彙**（unlikely / extremely unlikely）。後者對「如何把機率講給非專業者聽」這個主題特別有價值。
6. **必須說清楚的世代差異**：Nature 2005 的 STEP 即時系統（pasadena.wr.usgs.gov/step）已下線，USGS 現行 OAF 用的是 R–J 或 ETAS 加上 generic / sequence-specific / Bayesian 三種參數來源。**把 STEP 講成「USGS 現在在跑的系統」是錯的。**
7. **避開的內容**：STEP 空間網格 5 km、USGS 「20 分鐘首報、第一年更新 75 次」這些細節只在搜尋摘要層級出現，若要寫進網站請先自行點開來源確認，或直接略去。

**兩章的接點**：Jordan et al. (2011) 那句「OEF 必須與 PSHA 的長期預報一致地提供完整的危害描述——地動超越機率，而不只是短期破裂機率」是把 13 章與 16 章縫在一起的最佳引言；台灣端則用 TEM PSHA2020 的 BPT 斷層記憶當本地例證。
