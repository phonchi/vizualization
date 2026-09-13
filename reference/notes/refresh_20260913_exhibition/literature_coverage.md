# 文獻概念覆蓋與核對深度

本表登錄reference/的59個PDF檔案，對應目前20章主線、第二部及附錄。**不是59篇全文核讀完成證明。** 原PDF與文字快取不公開；這裡只保存書目、概念位置、核對範圍與限制。年份以內文為準，未核實處明列。

採用document-folder-ingest既有快取。原索引的59個ok是命令狀態，其中ETAS/[1998] ETAS.pdf只有24個分頁字元，列為**抽取內容不可用、原文未核**。同類核心概念由Reinhart 2018與Jalilian 2019的可讀原文支撐。

## 如何讀本表

- 選段只核對所列原PDF實體頁碼，不表示全文、全部公式或數值實驗已通過。
- 首頁只辨識書目/摘要，章節是相關主題位置，不是每項細節都已教完的認證。
- 數字為顯示章號，A–G為附錄；29為台灣章；隔離不可放進公開教材。
- SHA-256、快取位置、canonical ID與閱讀深度另存同目錄coverage.json。

## 本輪重要決策與未核項

1. Hawkes已補入11/12/A，連結marked時空ETAS；分支與穩定性沿用已有推導。
2. b-positive補入7/8/B.5，區分正相鄰差、絕對差、配對與漏測假設。已發表比較依據是Tinti–Gasperini 2024 PDF pp.3–4，不能宣稱正差法普遍更佳。SeismoStats本地為預印本，b-more-positive原始方法論文尚未獨立核讀，不假裝完成。
3. Girona–Drymoni 2024（Ncom）補入14，說明森林投票比例不直接等於校準格箱RATE；原文有訓練與後段資料分工，不誣指全資料訓練。
4. Rhoades 2020 PDF p.4定義catalogue lead time為目錄起點至目標事件，與發報delay不同；已交主代理同步全站術語。
5. 每個實際檔案都保留一列。未讀篇的次要模型不逐篇新增章；使用其獨有結論前，先核讀相關原文。

## 重複與錯名

- 二進位相同2組：2007兩份CA EEPAS；2023兩份Italy EEPAS。
- 非二進位同篇3組：兩份Tsai等2012；Bayona等2022根目錄Tests/Ensemble版本；Mizrahi等2020作者稿/2021發表版。
- Calibrated ETAS實為Mizrahi等2021除叢與規模分布，不是新校準ETAS模型。
- ETAS_R首頁為2019；Ncom為Girona–Drymoni2024；Italy_exp本地為2024 arXiv評分研究；soft_intro為溝通研究；EETAS實為EEPAS回顧；entropy-22實為2020。

## 逐來源覆蓋表

| ID／來源 | 真正題名與年份 | 概念／章節 | 核對深度（原PDF頁） | 公開限制 |
|---|---|---|---|---|
| REF-01<br>ETAS/[1995] Omori.pdf | The Centenary of the Omori Formula for a Decay Law of Aftershock Activity<br>1995 | Omori；有限窗；早期漏測<br>位置：10／B | 選段 pp.1<br>Omori率及早期漏測；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-02<br>ETAS/[1998] ETAS.pdf | Space-Time Point-Process Models for Earthquake Occurrences<br>1998（書目交叉辨識；原PDF未核） | 時空ETAS；原文抽取不可用<br>位置：12／13／C | 原文未核 pp.—<br>ok只是工具狀態；24字元全為分頁。標題年份交叉辨識，原PDF未核。 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-03<br>ETAS/[1999] Review_model.pdf | Seismicity Analysis through Point-process Modeling: A Review<br>1999 | Hawkes；條件強度；更新；應力釋放<br>位置：11／12／A／C／F | 選段 pp.26,27,28,29,33<br>Hawkes/更新/應力釋放與時間變換；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-04<br>ETAS/[2002] ETAS_Declustering.pdf | Stochastic Declustering of Space-Time Earthquake Occurrences<br>2002 | 隨機除叢；親代與背景權重<br>位置：13／C | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-05<br>ETAS/[2008] ETAS using EM.pdf | Estimation of Space–Time Branching Process Models in Seismology Using an EM–Type Algorithm<br>2008 | 潛在親代；EM估計<br>位置：13／C | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-06<br>ETAS/[2018] Review_model.pdf | A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications<br>2018 | Hawkes；分支；標記與ground process<br>位置：11／12／A／C | 選段 pp.2,3,4,5,6<br>Hawkes、分支、marked/ground、概似；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-07<br>ETAS/[2023] Bayes ETAS.pdf | Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package<br>2023 | 貝氏ETAS；近似後驗<br>位置：13／A／C | 選段 pp.1,2<br>後驗模式與近似；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-08<br>Ensemble/[2009] Mixture Models for Improved Short-Term Earthquake Forecasting.pdf | Mixture Models for Improved Short-Term Earthquake Forecasting<br>2009 | STEP–EEPAS混合；時間互補<br>位置：19／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-09<br>Ensemble/[2014] Multiplicative models.pdf | Regional Earthquake Likelihood Models II: Information Gains of Multiplicative Hybrids<br>2014 | 乘法hybrid；資訊增益<br>位置：19／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-10<br>Ensemble/[2022] Testing Multiplicative models in California.pdf | Prospective evaluation of multiplicative hybrid earthquake forecasting models in California<br>2022 | 乘法hybrid；前瞻檢驗<br>位置：19／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置；同篇/版本關係：REF-41 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-11<br>Ensemble/[2023] Maximizing the forecasting skill of an ensemble model.pdf | Maximizing the forecasting skill of an ensemble model<br>2023 | logistic權重；組合分數<br>位置：19／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-12<br>Intro_to_the_theory_of_point_processes.pdf | An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, Second Edition<br>2003 | Poisson；更新；群集；條件強度<br>位置：4／11／A／C／F | 首頁 pp.1<br>只核書名版次；本輪未重讀全書或定理 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-13<br>Psi regression/[2021] Space Time Tradeoff.pdf | Space-Time Trade-Off of Precursory Seismicity in New Zealand and California Revealed by a Medium-Term Earthquake Forecasting Model<br>2021 | Ψ時空抵換；迴歸方向<br>位置：14／D | 選段 pp.1,2<br>時空抵換問題設定；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-14<br>Psi regression/[2024] psi.pdf | Algorithmic Identification of the Precursory Scale Increase Phenomenon in Earthquake Catalogs<br>2024 | 自動Ψ搜尋；累積規模異常<br>位置：14／D | 選段 pp.1,2<br>自動搜尋及C(t)定義；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-15<br>Taiwan/[2011] Omori Law Taiwan.pdf | New Empirical Tests of the Multifractal Omori Law for Taiwan<br>2012；檔名2011 | 多重分形應力活化；p與主震規模<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-16<br>Taiwan/[2012] Aftershock Statistics of the 1999 Chi–Chi, Taiwan Earthquake and the Concept of Omori Times.pdf | Aftershock Statistics of the 1999 Chi–Chi, Taiwan Earthquake and the Concept of Omori Times<br>2013卷期；2012版權 | 集集餘震；Omori times<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-17<br>Taiwan/[2012] New Empirical Tests of the Multifractal Omori Law for Taiwan.pdf | New Empirical Tests of the Multifractal Omori Law for Taiwan<br>2012 | 與REF-15同篇不同PDF<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置；同篇/版本關係：REF-15 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-18<br>Taiwan/[2013] Taiwan Bath.pdf | Maximum magnitudes in aftershock sequences in Taiwan<br>2013 | 最大餘震；Båth；GR<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-19<br>Taiwan/[2015] b-Values Observations in Taiwan A Review.pdf | b-Values Observations in Taiwan: A Review<br>2015 | 台灣b值；時空異質性<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-20<br>Taiwan/[2016] An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized M w magnitudes.pdf | An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized Mw magnitudes<br>2016 | 均一化規模；歷史目錄<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-21<br>Taiwan/[2016] Studies on Aftershocks in Taiwan A Review.pdf | Studies on Aftershocks in Taiwan: A Review<br>2016 | 台灣餘震；空間時間與機制<br>位置：29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-22<br>Taiwan/[2024] AutoBats CMT Catalog.pdf | Taiwan Regional Moment Tensor Solution: The AutoBATS CMT Catalog<br>2024-01-30講稿 | AutoBATS；CMT；規模比較<br>位置：29／G | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 講稿/報告；公開授權未全面核實，不重傳PDF |
| REF-23<br>Taiwan/[2024] Some characteristics of foreshocks and aftershocks of the 2022 ML6.8 Chihshang, Taiwan, earthquake sequence.pdf | Some characteristics of foreshocks and aftershocks of the 2022 ML 6.8 Chihshang, Taiwan, earthquake sequence<br>2024 | 池上前震與餘震<br>位置：29 | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-24<br>Taiwan/[2025] Fast_Report_Dapu_ETAS.pdf | Fast report: performance of the ETAS model in forecasting aftershock occurrence and site-specific ground-shaking intensity for the 2025 Dapu, Taiwan, earthquake sequence<br>2025 | 台灣ETAS；場址地動<br>位置：29／F | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-25<br>Taiwan/[2026] Taiwan_EEPAS.pdf | First Application of the EEPAS Earthquake Forecasting Model to Taiwan: Catalog Sensitivity, Lead-Time Compensation, and Retrospective Validation<br>2026稿件；審稿中 | 僅書目隔離<br>位置：隔離 | 隔離書目 pp.1<br>只記書目，不採用稿件研究內容 | 審稿中隔離：不得公開未發表數字結果、檢驗分數、加速倍數或內部圖表；公開程式與義大利資料、已發表設定及參數仍可使用；台灣EEPAS正文僅說在地化工作正在進行中。公開repo不等於稿件授權。 |
| REF-26<br>Taiwan/[2026] Taiwan_EEPAS_sup.pdf | Electronic Supplement for: First Application of the EEPAS Earthquake Forecasting Model to Taiwan<br>2026審稿中補充 | 僅書目隔離<br>位置：隔離 | 隔離書目 pp.1<br>只記書目，不採用稿件研究內容 | 審稿中隔離：不得公開未發表數字結果、檢驗分數、加速倍數或內部圖表；公開程式與義大利資料、已發表設定及參數仍可使用；台灣EEPAS正文僅說在地化工作正在進行中。公開repo不等於稿件授權。 |
| REF-27<br>Taiwan/臺灣地區112年中大型地震震源資訊之快速彙整與提供 .pdf | 臺灣地區112年中大型地震震源資訊之快速彙整與提供<br>112年度（2023）；正式出版日未核 | RMT；震源反演；技術報告<br>位置：29／G | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 講稿/報告；公開授權未全面核實，不重傳PDF |
| REF-28<br>Taiwan/臺灣預報發展.pdf | 臺灣地震測報的發展<br>出版年未核 | 觀測網；速報；預警發展<br>位置：21／29／G | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 講稿/報告；公開授權未全面核實，不重傳PDF |
| REF-29<br>[2004] rhoades2004.pdf | Long-range Earthquake Forecasting with Every Earthquake a Precursor According to Scale<br>2004 | EEPAS原式；三核；權重<br>位置：9／14／15／D | 選段 pp.1,4<br>三核與事件權重原式；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-30<br>[2006] Space–time ETAS models and an improved extension.pdf | Space–time ETAS models and an improved extension<br>2006 | ETAS空間核；規模尺度<br>位置：12／13／C | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-31<br>[2007] Application of the EEPAS Model to Forecasting for CA.pdf | Application of the EEPAS Model to Forecasting Earthquakes of Moderate Magnitude in Southern California<br>2007 | 跨地區與目標規模<br>位置：15／D | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-32<br>[2007] Rhoades.pdf | Application of the EEPAS Model to Forecasting Earthquakes of Moderate Magnitude in Southern California<br>2007 | 同REF-31二進位重複<br>位置：15／D | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置；同篇/版本關係：REF-31 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-33<br>[2009] Long-range earthquake forecasting allowing for aftershocks.pdf | Long-range earthquake forecasting allowing for aftershocks<br>2009 | EAS；預報主震的餘震貢獻<br>位置：15／D | 選段 pp.1,2<br>對未來主震添加餘震的模型角色；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-34<br>[2011] Japan.pdf | Application of a long-range forecasting model to earthquakes in the Japan mainland testing region<br>2011 | EEPAS日本與跨規模應用<br>位置：15／D | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-35<br>[2017] ETAS_R.pdf | ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data<br>2019；檔名2017 | marked Hawkes；條件歷史；ETAS估計<br>位置：3／11／12／13／A／C | 選段 pp.1,3,4,5<br>真年份、marked Hawkes、核與穩定性；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-36<br>[2017] SEDA.pdf | SEDA: A software package for the Statistical Earthquake Data Analysis<br>2017 | 完整度；統計地震軟體<br>位置：7／8／13 | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-37<br>[2019] Brief Revie Process.pdf | Modeling the earthquake occurrence with time-dependent processes: a brief review<br>2019 | 更新；Weibull/Lognormal/BPT<br>位置：20／F | 選段 pp.1,4<br>複發分布與更新過程；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-38<br>[2021] Estimate b.pdf | The Effect of Declustering on the Size Distribution of Mainshocks<br>2020作者稿；2021發表版 | 除叢選擇效應；b值<br>位置：7／10／29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置；同篇/版本關係：REF-43 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-39<br>[2022] A 20-Year Journey of Forecasting with the EETAS model.pdf | A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model<br>2022；檔名EETAS錯字 | EEPAS；Janus；時間補償<br>位置：14／15／19／D／E | 選段 pp.1,8<br>Janus與EEPAS變體概念表；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-40<br>[2022] Pycesp.pdf | pyCSEP: A Python Toolkit for Earthquake Forecast Developers<br>2022 | 預報格式；目錄；一致性檢驗<br>位置：3／16／17／18／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-41<br>[2022] Tests.pdf | Prospective evaluation of multiplicative hybrid earthquake forecasting models in California<br>2022 | 與REF-10同篇；二元與乘法評分<br>位置：17／19／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-42<br>[2023] Application of the EEPAS earthquake forecasting model to Italy.pdf | Application of the EEPAS earthquake forecasting model to Italy<br>2023 | 義大利規格；參數；檢驗；三核<br>位置：2／3／9／12／13／15／16／17／18／19／20／D／E | 選段 pp.1,3,4,8,9,15,16,17,18,19,20<br>規格、表3、負二項/二元檢驗、附錄核與ETAS；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-43<br>[2023] Calibrated ETAS.pdf | The Effect of Declustering on the Size Distribution of Mainshocks<br>2021；檔名Calibrated ETAS與2023皆錯 | 除叢與規模分布；非新校準ETAS<br>位置：7／10／29／B | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-44<br>[2023] Italy_EEPAS.pdf | Application of the EEPAS earthquake forecasting model to Italy<br>2023 | 同REF-42二進位重複<br>位置：2／3／9／12／13／15／16／17／18／19／20／D／E | 選段 pp.1,3,4,8,9,15,16,17,18,19,20<br>與REF-42同一PDF共享選段核對；選段核對，不是全文或數值實驗驗證；同篇/版本關係：REF-42 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-45<br>[2023] SimpleETAS.pdf | SimplETAS: A Benchmark Earthquake Forecasting Model Suitable for Operational Purposes and Seismic Hazard Analysis<br>2023線上；2024卷期 | 七形狀固定；兩振幅；MPS19<br>位置：12／13／C | 選段 pp.1,2,3<br>Hawkes家族、固定參數理由與背景；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-46<br>[2023] Stat_power_test.pdf | Statistical power of spatial earthquake forecast tests<br>2023 | S-test功效；事件數<br>位置：18／E | 選段 pp.1,2<br>功效與S-test問題設定；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-47<br>[2023] Tests.pdf | Are Regionally Calibrated Seismicity Models More Informative than Global Models? Insights from California, New Zealand, and Italy<br>2023 | 全球/區域基準；資訊比較<br>位置：18／20／E | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-48<br>[2023] The need for open, transdisciplinary.pdf | The need for open, transdisciplinary, and ethical science in seismology<br>2023 | 開放科學；研究溝通<br>位置：20 | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-49<br>[2024] Estimate b.pdf | The estimation of b-value of the frequency–magnitude distribution and of its 1σ intervals from binned magnitude data<br>2024 | 離散GR；區間；正差/絕對差/配對<br>位置：7／8／B | 選段 pp.1,2,3,4,5<br>離散估計、正差/絕對差、配對；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-50<br>[2024] Italy_exp.pdf | Comparative evaluation of earthquake forecasting models: An application to Italy<br>2024 arXiv v1；現行出版狀態未核 | 一致性分數；校準/鑑別；非EEPAS論文<br>位置：18／E | 選段 pp.1,2<br>預印本與校準/一致性分數；選段核對，不是全文或數值實驗驗證 | 本地是預印本，不聲稱已發表；不公開全文快取 |
| REF-51<br>[2024] Ncom.pdf | Abnormal low-magnitude seismicity preceding large-magnitude earthquakes<br>2024 | 隨機森林異常分類；回溯；FEM<br>位置：14 | 選段 pp.1,2,6,7,8<br>投票比例、案例、FEM及訓練分工；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-52<br>[2024] New Pycsep.pdf | New Features in the pyCSEP Toolkit for Earthquake Forecast Development and Evaluation<br>2024 | QuadTree；檢驗；目錄式預報<br>位置：3／17／E | 選段 pp.1,5<br>QuadTree與預報格式；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-53<br>[2024] OEF_Review.pdf | Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions<br>2024 | OEF；R–J/STEP/ETAS；溝通<br>位置：1／13／20／F | 選段 pp.1,8<br>OEF與STEP定位；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-54<br>[2024] soft_intro.pdf | Long-term communication of aftershock forecasts: The Canterbury earthquake sequence in New Zealand<br>2024 | 預報溝通；非軟體導論<br>位置：20 | 首頁 pp.1<br>僅首頁書目/摘要；正文未核讀，章節僅主題位置 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |
| REF-55<br>[2025] EEPAS_Software.pdf | The EEPAS Model Revisited: Statistical Formalism and a High-Performance, Reproducible Open-Source Framework<br>2025稿件；submitted to GJI | 僅書目隔離<br>位置：隔離 | 隔離書目 pp.1<br>只記書目，不採用稿件研究內容 | 審稿中隔離：不得公開未發表數字結果、檢驗分數、加速倍數或內部圖表；公開程式與義大利資料、已發表設定及參數仍可使用；台灣EEPAS正文僅說在地化工作正在進行中。公開repo不等於稿件授權。 |
| REF-56<br>[2025] SeismosStats.pdf | SeismoStats: A Python Package for Statistical Seismology<br>2025 arXiv v1；現行出版狀態未核 | b-positive/b-more-positive；完整度；軟體<br>位置：7／8／B | 選段 pp.1,15,17,18<br>預印本、b-positive與延伸名稱；選段核對，不是全文或數值實驗驗證 | 本地是預印本，不聲稱已發表；不公開全文快取 |
| REF-57<br>[2026  PyEEPAS_sup.pdf | Supplementary Materials for: PyEEPAS: An Open-Source Python Framework for EEPAS Earthquake Forecasting with Formal Statistical Foundations<br>2026-04-08審稿中補充 | 僅書目隔離<br>位置：隔離 | 隔離書目 pp.1<br>只記書目，不採用稿件研究內容 | 審稿中隔離：不得公開未發表數字結果、檢驗分數、加速倍數或內部圖表；公開程式與義大利資料、已發表設定及參數仍可使用；台灣EEPAS正文僅說在地化工作正在進行中。公開repo不等於稿件授權。 |
| REF-58<br>[2026] PyEEPAS.pdf | PyEEPAS: An Open-Source Python Framework for EEPAS Earthquake Forecasting with Formal Statistical Foundations<br>2026審稿中稿件 | 僅書目隔離；公開repo另論<br>位置：隔離 | 隔離書目 pp.1<br>只記書目，不採用稿件研究內容 | 審稿中隔離：不得公開未發表數字結果、檢驗分數、加速倍數或內部圖表；公開程式與義大利資料、已發表設定及參數仍可使用；台灣EEPAS正文僅說在地化工作正在進行中。公開repo不等於稿件授權。 |
| REF-59<br>entropy-22-01264.pdf | The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data<br>2020；22為卷號 | catalogue lead time；缺失早期前兆<br>位置：15／D | 選段 pp.1,4<br>catalogue lead time的原文定義；選段核對，不是全文或數值實驗驗證 | 僅自寫摘要與來源連結；不公開原PDF或全文快取 |

## 使用邊界

此表是可接續的來源登錄，不是模型正確性證書。原始數值與軟體行為依各章驗證紀錄查證；首頁已讀與全文已驗證有明確區別。教材或來源更新時，同步本表及JSON。
