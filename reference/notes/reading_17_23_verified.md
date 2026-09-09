# 第 17–23 章延伸閱讀核實稿

核實日期：2026-09-09。用途：以下各章區塊可直接附在對應教材末尾；核實紀錄不附入教材。
本次讀取 book/_conventions.md、各章 Markdown、testing/ensemble/stats/oef/taiwan/psha_step_websearch 筆記，並以原始 PDF 與出版社／機構官網交叉核對。僅引用已公開材料，不使用 embargo 研究。

<!-- CHAPTER 17 START -->
## 參考資料與延伸閱讀

- [Theory of CSEP Tests](https://docs.cseptesting.org/getting_started/theory.html) — pyCSEP 開發團隊，官方文件（免費）。先對照各檢驗的目標、模擬方式與分位數分數，再看 N、S、M 與 conditional L-test 的程式範例，可把本章公式接到實際檢驗流程。
- [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。從「觀測是否與預報一致」出發整理檢驗方法，適合先建立觀念，再回頭理解本章為何要拆開事件數、空間與規模。
- [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。本章 Poisson／二元似然、負二項計數與多重檢定的主要實例來源；閱讀時留意餘震叢集如何影響各種評分。
<!-- CHAPTER 17 END -->

<!-- CHAPTER 18 START -->
## 參考資料與延伸閱讀

- [Theory of CSEP Tests：Forecast comparison tests](https://docs.cseptesting.org/getting_started/theory.html#forecast-comparison-tests) — pyCSEP 開發團隊，官方文件（免費）。從資訊增益 IGPE 的定義與 T-test 範例開始，對照本章的基準模型、率修正與信賴區間，理解「分數較高」與「差異顯著」的區別。
- [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。接著閱讀警報式預測與誤差圖的部分，比較不同評估方法回答的問題，銜接本章的 Molchan 圖與面積技能分數。
- [Statistical power of spatial earthquake forecast tests](https://doi.org/10.1093/gji/ggad030) — Asim M. Khawaja 等（2023），*Geophysical Journal International*；[免費機構典藏全文](https://gfzpublic.gfz.de/pubman/item/item_5015770_1)。以空間檢驗說明樣本量與網格如何改變統計功效，正好延伸本章「通過檢驗不等於模型有辨識力」的討論；可接著比較等寬網格與 Quadtree 的設計。
- [Enhancing the Statistical Evaluation of Earthquake Forecasts—An Application to Italy](https://doi.org/10.1785/0220240209) — Jonas R. Brehmer、Kristof Kraus、Tilmann Gneiting、Marcus Herrmann、Warner Marzocchi（2025），*Seismological Research Letters*（出版社全文可能需訂閱；[免費作者稿](https://arxiv.org/abs/2405.10712)）。從義大利預報案例延伸本章的模型比較、可靠度圖與校準，重點是評分函數如何對應預報目標，以及如何分開檢查校準與鑑別力。
<!-- CHAPTER 18 END -->

<!-- CHAPTER 19 START -->
## 參考資料與延伸閱讀

- [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等（2024），*Reviews of Geophysics*（免費全文）。先讀模型發展與各國作業系統中組合模型的段落，了解為何需要結合不同預報，以及權重選擇如何與檢驗制度連在一起。
- [Mixture Models for Improved Short-Term Earthquake Forecasting](https://doi.org/10.1785/0120080063) — David A. Rhoades、Matthew C. Gerstenberger（2009），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[SCEC 免費摘要](https://central.scec.org/node/3964)）。本章加法組合的原始案例，將 STEP 與 EEPAS 的不同時間尺度結合；重點是凸組合如何利用互補資訊，以及回溯估出的權重仍須接受前瞻檢驗。
- [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。對照本章乘法 hybrid 的前瞻測試，觀察回溯資訊增益為何未必延續到未來；適合與前一篇的組合動機一起讀。
- [Maximizing the forecasting skill of an ensemble model](https://doi.org/10.1093/gji/ggad020) — Marcus Herrmann、Warner Marzocchi（2023），*Geophysical Journal International*（免費全文）。本章 logistic 權重學習的主要來源；重點是直接最佳化整體組合的技巧，而非只依各成分模型的單獨成績配權重。
- [Regional Earthquake Likelihood Models II: Information Gains of Multiplicative Hybrids](https://doi.org/10.1785/0120140035) — D. A. Rhoades、M. C. Gerstenberger、A. Christophersen、J. D. Zechar、D. Schorlemmer、M. J. Werner、T. H. Jordan（2014），*Bulletin of the Seismological Society of America*（出版社全文可能需訂閱；[Bristol 大學免費全文](https://research-information.bris.ac.uk/files/49915043/Rhoades_Hybrids_RELM_BSSA_2014.pdf)）。本章乘法 hybrid、保序轉換與正規化的直接來源；讀完建構方法後，對照上面的 Bayona 等人前瞻測試，區分擬合改善與未來預報表現。
<!-- CHAPTER 19 END -->

<!-- CHAPTER 20 START -->
## 參考資料與延伸閱讀

- [Related Distributions：Hazard Function 與 Survival Function](https://itl.nist.gov/div898/handbook/eda/section3/eda362.htm) — NIST／SEMATECH，*e-Handbook of Statistical Methods*（免費）。先看機率密度、存活函數與危害函數的關係，釐清本章「已經等到現在，再發生的瞬時率」與一般發生機率的差別。
- [A Brownian model for recurrent earthquakes](https://doi.org/10.1785/0120010267) — Mark V. Matthews、William L. Ellsworth、Paul A. Reasenberg（2002），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[USGS 免費摘要](https://pubs.usgs.gov/publication/70024443)）。BPT 複發模型的核心論文，從帶有布朗擾動的載入過程推到首達時間，銜接本章的物理直覺與長時間危害率。
- [Modeling the earthquake occurrence with time-dependent processes: a brief review](https://doi.org/10.1007/s11600-019-00284-4) — Ourania Mangira、Christos Kourouklas、Dimitris Chorozoglou、Aggelos Iliopoulos、Eleftheria Papadimitriou（2019），*Acta Geophysica*（全文可能需訂閱）。將本章的更新過程、複發分布與應力釋放放回時間相依模型的全貌，閱讀時比較各模型適用的時間尺度、資料與假設。
<!-- CHAPTER 20 END -->

<!-- CHAPTER 21 START -->
## 參考資料與延伸閱讀

- [Introduction to Probabilistic Seismic Hazard Analysis](https://scits.stanford.edu/sites/g/files/sbiybj22081/files/media/file/baker_2013_intro_psha_v2_0.pdf) — Jack W. Baker（2013），White Paper Version 2.0（免費教材）。本章 PSHA 推導的入門來源，依五個步驟走過震源、規模、距離、地動與積分；接著讀反聚合及回歸期，理解危害曲線背後有哪些事件在貢獻。
- [TEM 計畫概要與研究成果](https://tem.tw/TEM2020/portfolio-overview.html) — Taiwan Earthquake Model 計畫團隊，TEM2020 專頁（免費）。先看孕震構造、測地資料與危害度評估如何分工，再循頁面列出的論文了解台灣模型的輸入來源。
- [Probabilistic seismic hazard assessment for Taiwan: TEM PSHA2020](https://doi.org/10.1177/8755293020951587) — Chung-Han Chan 等（2020），*Earthquake Spectra*（出版社全文可能需訂閱；[中央大學免費全文](https://www.gep.ncu.edu.tw/storage/thesis/2020/2020%20Chung-Han%20Chan_ES2.pdf)）。本章台灣案例的正式文獻，說明孕震構造、地震目錄、地動預估式與場址效應如何納入 PSHA；可對照教材檢查每一類不確定性出現在哪一步。
<!-- CHAPTER 21 END -->

<!-- CHAPTER 22 START -->
## 參考資料與延伸閱讀

- [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS，官方餘震預報說明（免費）。先看預報表、時間窗、機率與事件數範圍的讀法，再閱讀模型參數頁籤說明，對照本章從統計模型到發布產品的流程。
- [Earthquake forecasts](https://www.geonet.org.nz/earthquake/forecast/) — GeoNet，官方預報說明（免費）。用較少公式說明機率預報的用途與限制，適合比較本章紐西蘭案例如何把數字連到不同使用者的需求。
- [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等（2024），*Reviews of Geophysics*（免費全文）。本章各國系統與機率溝通的主要回顧來源，優先讀義大利、紐西蘭與美國的制度比較及使用者共同設計；論文記錄的是出版時的狀態，後續產品可配合前兩項官網閱讀。
- [Real-time forecasts of tomorrow's earthquakes in California](https://doi.org/10.1038/nature03622) — Matthew C. Gerstenberger、Stefan Wiemer、Lucile M. Jones、Paul A. Reasenberg（2005），*Nature*（全文可能需訂閱）。本章 STEP 系統的原始論文，說明如何結合背景地震率與叢集模型，產生未來 24 小時的強震動機率圖；可與上面的 USGS OAF 說明比較不同時期的預報產品。
<!-- CHAPTER 22 END -->

<!-- CHAPTER 23 START -->
## 參考資料與延伸閱讀

- [最近地震](https://www.cwa.gov.tw/V8/C/E/index.html) — 交通部中央氣象署，官方資料頁（持續更新；免費）。從實際地震報告查看發震時間、震央、深度、規模與震度，對照本章對資料尺度與目錄用途的區分；單筆事件報告與經整理的研究目錄各有不同用途。
- [An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized Mw magnitudes](https://doi.org/10.1186/s40623-016-0414-4) — Wen-Yen Chang、Kuei-Pao Chen、Yi-Ben Tsai（2016），*Earth, Planets and Space*（免費全文）。本章規模均一化與歷史目錄的主要來源，重點是跨年代、跨測網的規模如何對齊；使用轉換結果時，一併查看文章頁連結的勘誤。
- [Fast report: performance of the ETAS model in forecasting aftershock occurrence and site-specific ground-shaking intensity for the 2025 Dapu, Taiwan, earthquake sequence](https://doi.org/10.1007/s44195-025-00097-7) — Ming-Che Hsieh、Min-Hsuan Chang、Yu-Chen Tai、Chun-Te Chen、Ting-Ying Lu（2025），*Terrestrial, Atmospheric and Oceanic Sciences*（免費全文）。沿著本章大埔案例閱讀即時目錄、ETAS 模擬與場址地動的完整流程，並留意單一序列的結果能支持哪些判斷、仍有哪些跨序列驗證待做。
<!-- CHAPTER 23 END -->

## 原始 reference 檔案與章節關係

下列 PDF 已用 pdftotext 直接讀首頁及摘要；有封面者續讀第 2–3 頁。未複製原始 PDF 到公開網站。

| 原始相對檔案 | 核實內容與引用關係 |
|---|---|
| `reference/[2022] Tests.pdf` | Bayona et al. 2022；首頁 DOI ggac018、作者、OA 與摘要；引用於 17、19。 |
| `reference/[2023] Stat_power_test.pdf` | Khawaja et al. 2023；首頁 DOI ggad030、作者與網格／功效摘要；引用於 18。 |
| `reference/Ensemble/[2009] Mixture Models for Improved Short-Term Earthquake Forecasting.pdf` | Rhoades & Gerstenberger 2009；首頁 DOI 0120080063、STEP/EEPAS 凸組合與尚待前瞻檢驗的摘要；引用於 19。 |
| `reference/Ensemble/[2023] Maximizing the forecasting skill of an ensemble model.pdf` | 第 1 頁為下載封面；第 2 頁正式首頁列 2023、ggad020 與 logistic 組合摘要；引用於 19。 |
| `reference/Ensemble/[2014] Multiplicative models.pdf` | 第 1 頁 Bristol 典藏書目，第 3 頁正式首頁；確認 Rhoades et al. 2014 乘法 hybrid 是教材公式的直接來源；引用於 19，與 Bayona 2022 的前瞻驗證搭配。 |
| `reference/[2019] Brief Revie Process.pdf` | 第 1 頁封面、第 3 頁文章摘要；Mangira et al. 2019、s11600-019-00284-4；引用於 20。第 2 頁為個人使用版權說明，不提供此本機副本的公開下載。 |
| `reference/[2024] OEF_Review.pdf` | 首頁作者、2024、2023RG000823、CC-BY 與模型／檢驗／溝通三主題；引用於 19、22。 |
| `reference/Taiwan/[2016] An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized M w magnitudes.pdf` | 首頁三位作者、s40623-016-0414-4、OA 與均一化目錄摘要；引用於 23。 |
| `reference/Taiwan/[2025] Fast_Report_Dapu_ETAS.pdf` | 首頁五位作者、2025、36:12、s44195-025-00097-7、OA 及 ETAS+GMM 摘要；引用於 23。 |
| `reference/[2024] Italy_exp.pdf` | 首頁為 arXiv:2405.10712v1（2024），舊題 Comparative evaluation of earthquake forecasting models: An application to Italy，五位作者與摘要直接核對；正式發表已改題並於 2025 出版，引用正式版本於 18。 |

Baker 2013、TEM2020 是 `reference/notes/psha_step_websearch.md` 原已使用的網路來源，本次重新開啟官網核實。CORSSA、NIST、pyCSEP、USGS、GeoNet 與 CWA 為學生補充免費入口；Matthews 2002 為教材 BPT 的核心論文，經作者所屬 USGS 出版資料庫核實。

## 網路核實記錄

- CORSSA：直接開啟上述 PDF，第 1 頁指定引用為 Zechar 2010、DOI 10.5078/corssa-77337879；本文包含一致性、模型比較與警報型預測評估。
- pyCSEP：上述官方 theory 頁及 forecast-comparison-tests 段已讀取，含各檢驗目的、公式與範例；未把爬取日期當出版年。
- Bayona 2022：[出版社文章頁](https://academic.oup.com/gji/article/229/3/1736/6510818)，作者、卷頁、DOI、CC-BY 全文均核實。
- Khawaja 2023：[GFZ 機構典藏](https://gfzpublic.gfz.de/pubman/item/item_5015770_1) 有書目、摘要與 publisher-version PDF；論文第一作者依 PDF 正式署名採 Asim M. Khawaja，避免典藏重排姓名造成歧義。
- Rhoades & Gerstenberger 2009：[SCEC 研究成果頁](https://central.scec.org/node/3964) 與本機 PDF 核對兩位作者、2009、BSSA 99(2A):636–646、DOI；僅確定免費摘要，全文保守標示可能需訂閱。
- Herrmann & Marzocchi 2023：[出版社全文](https://academic.oup.com/gji/article/234/1/73/6994524) 核對題名、兩位作者、2023、234(1):73–87、DOI ggad020 與 CC-BY。
- NIST：上述官方 handbook 頁含 h=f/S 與 survival 定義；無出版年。
- Matthews et al. 2002：[USGS 出版資料庫](https://pubs.usgs.gov/publication/70024443) 核實三位作者、2002、BSSA 92(6):2233–2250、DOI 0120010267 與 BPT 摘要；未確認免費正式全文，不標示 OA。
- Mangira et al. 2019：[Springer 正式頁](https://link.springer.com/article/10.1007/s11600-019-00284-4) 核實五位作者、2019-04-10、67:739–752、DOI；明示 subscription content，故標示全文可能需訂閱。
- Baker 2013：Stanford SCITS PDF 79 頁，前頁 Preferred citation 2013 Version 2.0；目錄與圖 2.4 明示 PSHA 五步，含 deaggregation/return periods。使用真正全文 URL，不使用僅留下新版書籍通知的舊地址。
- TEM2020：[官方計畫頁](https://tem.tw/TEM2020/portfolio-overview.html) 與 [出版社書目](https://journals.sagepub.com/doi/10.1177/8755293020951587) 核對 Chan 等、2020、36(1_suppl):137–159、DOI；中央大學作者機構 PDF 已由搜尋取得全文，教材提供該正式機構連結。
- USGS OAF、GeoNet 及 CWA：上述官方頁均成功開啟；OAF Model Parameters 明列 R–J 或 ETAS，引文沒有將 STEP 誤稱為現行 OAF。網頁以無年份／持續更新標示。
- Mizrahi et al. 2024：[出版社頁](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2023RG000823) 與 [GFZ 典藏](https://gfzpublic.gfz.de/pubman/faces/ViewItemFullPage.jsp?itemId=item_5028991_3) 核對題名、作者、2024、62(3)、e2023RG000823、CC-BY。
- Chang et al. 2016：[Springer 正式頁](https://link.springer.com/article/10.1186/s40623-016-0414-4) 核實三位作者、2016-03-17、68:45、OA，並發現有 2016-07-18 勘誤連結；導讀已提醒。
- Hsieh et al. 2025：[Springer 正式全文](https://link.springer.com/article/10.1007/s44195-025-00097-7) 核實作者、2025-04-03、36:12、DOI、CC-BY 與 ETAS+GMM 流程；不用機構年報所記 6 月作正式出版日。

- Gerstenberger et al. 2005：[Nature 正式頁](https://www.nature.com/articles/nature03622) 核實四位作者、2005-05-19、435:328–331、nature03622、STEP 摘要與訂閱標示；原先已在 psha_step_websearch.md 引用，這次納入第 22 章。

- Rhoades et al. 2014：[Bristol 作者機構典藏](https://research-information.bris.ac.uk/en/publications/regional-earthquake-likelihood-models-ii-information-gains-of-mul/) 核實七位作者、104(6):3072–3083、0120140035，頁面提供 Final published version PDF；新增於第 19 章第 5 筆。
- Brehmer et al.：[作者 GitHub](https://github.com/jbrehmer42/Earthquakes_Italy)、[更新後 arXiv 頁](https://arxiv.org/abs/2405.10712) 與[作者機構典藏](https://www.iris.unina.it/handle/11588/1005033) 共同確認本機 2024 預印本已正式發表為 2025 SRL 96(3):1966–1988、0220240209；正式題名 Enhancing the Statistical Evaluation of Earthquake Forecasts—An Application to Italy。本文的 consistent scoring／reliability／calibration 與 18.6、研究前沿直接相關，補為第 18 章第 4 筆。出版社轉址受工具限制，未臆測 OA，提供 arXiv 免費作者稿。
- 學生區塊依需求省略官網的「未標示年份」，僅在內部記錄保留來源無年分資訊。

## 限制與排除

- DOI 轉址偶有 web 工具錯誤，但已用上述出版社正式 URL 或作者機構典藏確認書目，並非把抓取錯誤當作不存在。
- 未將付費文章的本機 PDF 公開；未使用未發表 EEPAS 軟體／台灣研究，也未擴張原文的預測能力主張。
- 這是延伸閱讀核實，不是正文數學與現行國家制度的全章重新稽核；未對原有正文數字作新的正確性宣告。

- 最終 HTTP 核對發現 Zechar 的 DOI 雖有效登錄，但轉址至舊 PDF 路徑回傳 404；教材標題改連已驗證 200 的 CORSSA 官方 PDF，DOI 保留為書目文字。
