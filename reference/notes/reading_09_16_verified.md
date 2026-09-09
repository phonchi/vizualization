# 第 09–16 章延伸閱讀：已核實素材

核實日期：2026-09-09。依 book/_conventions.md 的章節分工與公開限制選材。
以下每個 `### book/...` 與下一個同級標題之間，是可直接附入教材的區塊。
第 15、16 章以已發表的開放回顧論文作為入門讀物；不引用審稿中材料。

### book/09_forecasting_intro.py

## 參考資料與延伸閱讀

- [Can you predict earthquakes?](https://www.usgs.gov/faqs/can-you-predict-earthquakes) — USGS；免費官方說明。先讀這份短問答，釐清「預測特定地震」與「估計發生機率」的差別，再回看本章的三種用語。
- [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS；免費官方說明。用實際預報介面解釋期望次數、規模門檻與時間窗，適合對照本章「一份預報長什麼樣」。
- [Operational Earthquake Forecasting: State of Knowledge and Guidelines for Utilization](https://doi.org/10.4401/ag-5350) — Thomas H. Jordan 等，2011，*Annals of Geophysics*；[免費出版版全文](https://www.research.ed.ac.uk/files/10773224/AnnalsofGeoPhys.pdf)。這是 OEF 的重要報告，建議先讀摘要與建議，了解為什麼機率模型、公開溝通與權責制度必須一起考慮。

- [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等，2024，*Reviews of Geophysics*；免費開放全文。從第 2 節的領域導覽讀起，再比較義大利、紐西蘭與美國系統，了解本章所說的基準比較、前瞻檢驗與透明溝通如何進入實務。

### book/10_point_process.py

## 參考資料與延伸閱讀

- [Discrete Stochastic Processes, Chapter 2: Poisson Processes](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/resources/mit6_262s11_chap02/) — Robert G. Gallager／MIT OpenCourseWare，2011；免費講義。從等待時間與事件計數建立 Poisson 過程，適合在閱讀本章條件強度與概似推導前補足機率基礎。
- [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01) — Abdollah Jalilian，2019，*Journal of Statistical Software*；[免費全文與程式](https://www.jstatsoft.org/article/view/v088c01)。先讀模型介紹與估計方法，把本章的點過程語言對上實際地震模型；理解數學部分不需要先學 R。
- [Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes](https://doi.org/10.1080/01621459.1988.10478560) — Yosihiko Ogata，1988，*Journal of the American Statistical Association*；出版社全文可能需訂閱，[研究機構提供的全文](https://bemlar.ism.ac.jp/zhuang/Refs/Refs/ogata1988.pdf)可免費閱讀。重點是模型概似與時間變換後的殘差分析，對應本章「如何檢查條件強度是否漏掉結構」。

### book/11_catalog_completeness_b.py

## 參考資料與延伸閱讀

- [Estimate magnitude of completeness](https://seismostats.readthedocs.io/latest/user/estimate_mc.html) — SeismoStats 開發團隊；免費官方文件。用公式與短例子比較 MAXC、KS 與 b 值穩定法，閱讀時可對照本章對門檻選擇及方法假設的討論。
- [Estimate b-value](https://seismostats.readthedocs.io/latest/user/estimate_b.html) — SeismoStats 開發團隊；免費官方文件。整理傳統與差分式估計法，適合查閱規模離散化、完整度門檻與各方法適用條件。
- [The estimation of b-value of the frequency–magnitude distribution and of its 1σ intervals from binned magnitude data](https://doi.org/10.1093/gji/ggae159) — S. Tinti、P. Gasperini，2024，*Geophysical Journal International*；[免費出版版全文與補充資料](https://cris.unibo.it/handle/11585/980514)。深入比較分箱規模資料的估計式與不確定度，適合讀完本章連續、離散估計推導後，進一步理解哪些修正有統計依據。

- [b-Values Observations in Taiwan: A Review](https://doi.org/10.3319/TAO.2015.04.28.01%28T%29) — Jeen-Hwa Wang、Kou-Cheng Chen、Pei-Ling Leu、Jeng-Hsin Chang，2015，*Terrestrial, Atmospheric and Oceanic Sciences*；[免費全文](https://tao.cgu.org.tw/index.php/articles/archive/geophysics/item/1361-2015042801t)。整理台灣 b 值的時空變化、構造背景及估計研究，適合對照本章的台灣目錄例子；閱讀前兆相關解釋時，仍需分辨回溯觀察與前瞻檢驗。

### book/12_clustering_laws.py

## 參考資料與延伸閱讀

- [Aftershock Forecast Overview](https://earthquake.usgs.gov/data/oaf/overview.php) — USGS；免費官方說明。先看「How the Aftershock Forecasts Work」，把餘震產能、發生率衰減與規模分布三件事分開理解，再回頭閱讀本章的經驗律。
- [The Centenary of the Omori Formula for a Decay Law of Aftershock Activity](https://doi.org/10.4294/jpe1952.43.1) — Tokuji Utsu、Yosihiko Ogata、Ritsuko S. Matsu'ura，1995，*Journal of Physics of the Earth*；[免費全文](https://www.jstage.jst.go.jp/article/jpe1952/43/1/43_1_1/_article)。這篇回顧整理 Omori 公式、參數擬合與早期漏測問題，尤其適合追讀本章對 c 值與 p 值的解釋。
- [The Effect of Declustering on the Size Distribution of Mainshocks](https://doi.org/10.1785/0220200231) — Leila Mizrahi、Shyam Nandan、Stefan Wiemer，2021，*Seismological Research Letters*；出版社全文可能需訂閱，[免費作者預印本](https://arxiv.org/abs/2012.09053)。以觀測與合成目錄研究除叢如何改變規模分布，對應本章「除叢也是一種資料選擇」的討論。

- [New Empirical Tests of the Multifractal Omori Law for Taiwan](https://doi.org/10.1785/0120110237) — Ching-Yi Tsai、Guy Ouillon、Didier Sornette，2012，*Bulletin of the Seismological Society of America*；全文可能需訂閱。這是本章台灣 p 值與主震規模關係的原始來源，重點閱讀如何處理長短期目錄完整度，以及如何用不同除叢方法檢查結果。

### book/13_etas_structure.py

## 參考資料與延伸閱讀

- [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01) — Abdollah Jalilian，2019，*Journal of Statistical Software*；[免費全文與程式](https://www.jstatsoft.org/article/view/v088c01)。先讀模型結構與參數說明，把本章的背景率、時間核、空間核對照到一套完整的時空 ETAS 表達式。
- [ETAS: Epidemic-Type Aftershock Sequence](https://github.com/lmizrahi/etas) — Leila Mizrahi 等；免費作者程式庫與說明。可從目錄模擬範例追讀 ETAS 如何生成事件序列，對照本章的世代分解與 branching 模擬；程式所用參數慣例需先與本章核對。
- [Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes](https://doi.org/10.1080/01621459.1988.10478560) — Yosihiko Ogata，1988，*Journal of the American Statistical Association*；出版社全文可能需訂閱，[免費研究機構全文](https://bemlar.ism.ac.jp/zhuang/Refs/Refs/ogata1988.pdf)。回到 ETAS 的經典來源，閱讀事件歷史如何進入條件強度，以及作者如何用資料比較不同叢集模型。

- [Space–time ETAS models and an improved extension](https://doi.org/10.1016/j.tecto.2005.10.016) — Yosihiko Ogata、Jiancang Zhuang，2006，*Tectonophysics*；出版社全文可能需訂閱，[免費作者機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。這是本章時空 ETAS 結構的重要來源，說明如何把 Omori 衰減與餘震區尺度關係放進條件強度，並比較模型延伸。

### book/14_etas_estimation.py

## 參考資料與延伸閱讀

- [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01) — Abdollah Jalilian，2019，*Journal of Statistical Software*；[免費全文與程式](https://www.jstatsoft.org/article/view/v088c01)。沿著資料整理、模型擬合與結果診斷讀一遍，對照本章的背景率估計及隨機除叢；數學方法與 R 範例可以分開閱讀。
- [ETAS: Epidemic-Type Aftershock Sequence](https://github.com/lmizrahi/etas) — Leila Mizrahi 等；免費作者程式庫與說明。README 將參數估計、模擬及變動完整度的實作連回各自論文，適合延伸本章「不完整資料如何影響參數」的問題。
- [The Effect of Declustering on the Size Distribution of Mainshocks](https://doi.org/10.1785/0220200231) — Leila Mizrahi、Shyam Nandan、Stefan Wiemer，2021，*Seismological Research Letters*；出版社全文可能需訂閱，[免費作者預印本](https://arxiv.org/abs/2012.09053)。用合成目錄檢查除叢後的估計偏差，適合思考本章的背景／觸發分類結果能否直接當成物理事實。

- [SimplETAS: A Benchmark Earthquake Forecasting Model Suitable for Operational Purposes and Seismic Hazard Analysis](https://doi.org/10.1785/0220230199) — Simone Mancini、Warner Marzocchi，2023 年線上發表，*Seismological Research Letters*；全文可能需訂閱，[作者程式庫](https://github.com/smancini2/simplETAS)可免費閱讀。直接對照本章參數精簡的理由，思考固定哪些參數能減少估計困難，以及簡化模型應如何接受樣本外檢驗。

### book/15_psi_phenomenon.py

## 參考資料與延伸閱讀

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 2 節的 Ψ 現象與尺度關係，再讀限制與未解問題，掌握本章經驗觀察如何連到機率模型。
- [Algorithmic Identification of the Precursory Scale Increase Phenomenon in Earthquake Catalogs](https://doi.org/10.1785/0220240233) — Annemarie Christophersen、David A. Rhoades、Sebastian Hainzl，2024，*Seismological Research Letters*；[免費機構典藏全文](https://gfzpublic.gfz.de/rest/items/item_5029405_4/component/file_5029659/content)。這是本章自動辨識與對照實驗的已發表來源，建議比較矩形、圓形搜尋與隨機化目錄的設計；辨識到統計現象仍須與前瞻預報能力分開判斷。
- [Space–Time Trade-Off of Precursory Seismicity in New Zealand and California Revealed by a Medium-Term Earthquake Forecasting Model](https://doi.org/10.3390/app112110215) — Sepideh J. Rastin、David A. Rhoades、Annemarie Christophersen，2021，*Applied Sciences*；免費開放全文。研究時間與空間參數的取捨，適合延伸本章「同一事件的 Ψ 辨識不唯一」及參數解讀問題。

### book/16_eepas_ppe.py

## 參考資料與延伸閱讀

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 3–5 節，建立 EEPAS、模型組合與缺失前兆補償的全貌，再回到本章的正規化推導。
- [Long-range Earthquake Forecasting with Every Earthquake a Precursor According to Scale](https://doi.org/10.1007/s00024-003-2434-9) — David A. Rhoades、Frank F. Evison，2004，*Pure and Applied Geophysics*；全文可能需訂閱。這是 EEPAS 原始論文，重點是如何把尺度關係轉成每個事件對未來地震率的貢獻，而不是先判定哪個事件必然是前兆。
- [Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123) — Emanuele Biondini、David A. Rhoades、Paolo Gasperini，2023，*Geophysical Journal International*；出版社全文可能需訂閱，[免費機構典藏全文](https://www.earth-prints.org/handle/2122/17084)。將本章公式連到義大利目錄的實際應用，閱讀 PPE、ETAS 與 EEPAS 的比較時，特別留意學習期、測試期和預報時間窗。

- [The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data](https://doi.org/10.3390/e22111264) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2020，*Entropy*；免費開放全文。說明目錄開始前遺漏的事件如何影響 EEPAS，適合延伸本章的 lead time 與時間完整度補償推導。
- [Long-range earthquake forecasting allowing for aftershocks](https://doi.org/10.1111/j.1365-246X.2008.04083.x) — D. A. Rhoades，2009，*Geophysical Journal International*；[出版社網頁全文](https://academic.oup.com/gji/article/178/1/244/644120)可免費閱讀。閱讀 EEPAS 如何加入預報事件的餘震貢獻，對照本章 EAS 延伸與「降低輸入餘震權重」的不同角色。

### 核實記錄（內部，不附入教材）

- USGS 的 FAQ 與 OAF overview 均由官方頁面確認主題及免費內容；未將抓取日期冒充出版年。OAF overview 使用實際存在的 `earthquake.usgs.gov/data/oaf/overview.php`，不是不存在的 programs 路徑。
- Jordan et al.：University of Edinburgh 作者機構[書目頁](https://www.research.ed.ac.uk/en/publications/operational-earthquake-forecasting-state-of-knowledge-and-guideli/)與出版版 PDF 確認 2011、54(4)、315–391、DOI 10.4401/ag-5350、CC BY 開放全文。
- Gallager：MIT OCW 課程資源頁確認作者、Spring 2011 與 Chapter 2 Poisson Processes PDF。
- Jalilian：JSS 出版頁確認完整標題、作者、2019、Vol. 88 Code Snippet 1、DOI 10.18637/jss.v088.c01，頁面附 PDF 與程式。修正本機素材檔名的 2017 誤導。
- Ogata：統計數理研究所提供的出版版 PDF 與[作者論文表](https://www.ism.ac.jp/~ogata/papers.html)確認 1988、83(401)、9–27 與內容；DOI 10.1080/01621459.1988.10478560 亦由文獻索引對照。DOI 出版頁在研究工具內無法讀取，不據此宣稱出版社免費全文；附可讀的機構 PDF。
- SeismoStats：直接讀取官方 `estimate_mc.html`、`estimate_b.html`；它們是持續更新的網站，不用未核實的軟體論文出版年份代替文件年份。這裡只是延伸閱讀，不要求教材環境安裝新套件。
- Tinti & Gasperini：University of Bologna 作者機構[典藏頁](https://cris.unibo.it/handle/11585/980514)核對 2024、238(1)、433–458、DOI、兩位作者與 CC BY 出版版 PDF；Oxford 238(1) 卷期頁交叉確認。
- Utsu et al.：J-STAGE 出版頁確認 1995、43(1)、1–33、DOI、免費 PDF。作者解析頁偶將 Ritsuko S. Matsu'ura 拆開，採出版版的三位作者寫法。
- Mizrahi et al.：作者 [ETH/SED 論文頁](https://www.seismo.ethz.ch/de/home/bssed-imagewithtext/Leila-Mizrahi-Publications/)與作者 GitHub README 確認 2021、三位作者與 DOI；arXiv 2012.09053 為相同工作的 2020 作者預印本，因此教材使用正式論文年份 2021，免費連結明標版本。
- 作者 GitHub：直接讀取 lmizrahi/etas README，確認模擬、估計、完整度功能與兩篇 2021 論文的對應；不把 repository 當同儕審查論文。
- Rhoades et al. 2022：MDPI 原始出版搜尋內容核對三位作者、2022-09-19、12(9)、349、DOI 與全文；直接 open 偶有工具錯誤，搜尋可取得出版者全文內容。
- Christophersen et al. 2024：GFZ 作者機構全文含書目封面，確認三位作者、95(6)、3464–3481、DOI；典藏提供作者稿與附加頁，教材稱「機構典藏全文」，不稱免費出版社版。
- Rastin et al. 2021：MDPI 的 DOI 出版頁與出版社書籍收錄 PDF 確認三位作者、2021、11、10215、DOI 與 CC BY。
- Rhoades & Evison 2004：直接讀取 Springer 原始出版頁，核對兩位作者、161、47–72、DOI 與摘要；該頁提供購買／機構存取選項，標為全文可能需訂閱。
- Biondini et al. 2023：INGV Earth-prints 作者機構典藏核對已發表狀態、三位作者、234、1681–1700、DOI 與合法自存政策，並提供出版 PDF；另由 Bologna 典藏 PDF 首頁交叉確認。

本次沒有把禁止公開來源、未核實的參數或數值結果寫入章末素材。書目中的部分 DOI 在 web 工具內直接開啟失敗，已使用原始出版網站或作者機構來源核實；不將工具失敗當成 DOI 不存在。


### 原始 reference PDF 對照與抽取核實

本次對下列原始檔直接執行 `pdftotext -f 1 -l 1 -layout <PDF> -`，讀取出版首頁與摘要；未修改或再發布本地 PDF。對照不是根據檔名猜測。

| 教材章節 | 原始相對檔案 | 首頁核實與引用關係 |
|---|---|---|
| 09 | `reference/[2024] OEF_Review.pdf` | Mizrahi et al. 2024，10.1029/2023RG000823，CC BY；補入領域回顧。 |
| 10、13、14 | `reference/[2017] ETAS_R.pdf` | 首頁明寫 January 2019，Vol. 88 Code Snippet 1；教材引用 2019，涵蓋模型、估計與診斷。 |
| 11 | `reference/[2024] Estimate b.pdf` | Tinti & Gasperini 2024，10.1093/gji/ggae159；規模分箱估計與不確定度。 |
| 12、14 | `reference/[2021] Estimate b.pdf` | 首頁為 2020-12-16 投稿稿；同篇正式發表版見下一列。 |
| 12、14 | `reference/[2023] Calibrated ETAS.pdf` | 實為 Mizrahi, Nandan & Wiemer 2021，92、2333–2342，10.1785/0220200231；不是 2023 校準論文。 |
| 13 | `reference/[2006] Space–time ETAS models and an improved extension.pdf` | Ogata & Zhuang，Tectonophysics 413 (2006) 13–23，10.1016/j.tecto.2005.10.016；加入時空結構來源。 |
| 14 | `reference/[2023] SimpleETAS.pdf` | PDF cite-as 為 Mancini & Marzocchi (2023)，卷 95、38–49；教材明寫 2023 線上發表，避免混同 2024 卷期年份。 |
| 15、16 | `reference/[2022] A 20-Year Journey of Forecasting with the EETAS model.pdf` | 首頁實為 EEPAS 模型回顧，2022-09-19，10.3390/geosciences12090349；不沿用錯字 EETAS。 |
| 15 | `reference/Psi regression/[2021] Space Time Tradeoff.pdf` | Rastin, Rhoades & Christophersen，2021-10-31，10.3390/app112110215；時空取捨。 |
| 15 | `reference/Psi regression/[2024] psi.pdf` | Christophersen, Rhoades & Hainzl，2024，95、3464–3481，10.1785/0220240233；兩種演算法與隨機化對照。 |
| 16 | `reference/[2004] rhoades2004.pdf` | Rhoades & Evison，2004，161、47–72，10.1007/s00024-003-2434-9；原始 EEPAS 定義。 |
| 16 | `reference/[2023] Application of the EEPAS earthquake forecasting model to Italy.pdf` | Biondini, Rhoades & Gasperini，2023，234、1681–1700，10.1093/gji/ggad123；跨實作、基準模型及偽前瞻比較。 |

額外線上核實：2006 論文的 [ScienceDirect 出版頁](https://www.sciencedirect.com/science/article/abs/pii/S0040195105004889)確認年、作者、DOI 與摘要；simplETAS 作者程式庫及大學機構研究計畫文獻表確認標題與 2023 年引用，與本機出版 PDF 一致。OEF 回顧由 [USGS 作者機構出版頁](https://www.usgs.gov/publications/developing-testing-and-communicating-earthquake-forecasts-current-practices-and-future)及 Wiley 原始出版搜尋頁交叉核對。


### 補回既有原始來源（第 11、12、16 章）

| 教材章節 | 原始相對檔案 | PDF 首頁與線上核實 |
|---|---|---|
| 11 | `reference/Taiwan/[2015] b-Values Observations in Taiwan A Review.pdf` | 四位作者、2015、26(5)、475–492、DOI 10.3319/TAO.2015.04.28.01(T) 與 TAO 原始出版頁一致；頁面附免費 PDF。DOI 的括號以 URL percent encoding 保留。 |
| 12 | `reference/Taiwan/[2012] New Empirical Tests of the Multifractal Omori Law for Taiwan.pdf` | Tsai、Ouillon、Sornette；2012、102(5)、2128–2138、DOI 10.1785/0120110237；出版者提交的 Crossref multiple-resolution 書目與作者 arXiv 1109.5017 交叉確認。 |
| 16 | `reference/entropy-22-01264.pdf` | Rhoades、Rastin、Christophersen；2020-11-06、22、1264、DOI 10.3390/e22111264；MDPI 作者 2022 回顧文獻表交叉確認，不把卷號 22 當年份。 |
| 16 | `reference/[2009] Long-range earthquake forecasting allowing for aftershocks.pdf` | 單一作者 D. A. Rhoades，2009、178、244–256、DOI 10.1111/j.1365-246X.2008.04083.x；Oxford 出版頁確認年份與免費 HTML 內文，不把 DOI 內 2008 當出版年。 |

以上四份 PDF 均直接以 pdftotext 讀取首頁與摘要，沒有重寫原始 PDF。第 12 章 2011 作者預印本與 2012 出版版的摘要係數不同，因此學生清單只連正式 DOI，也不轉述數值；本機 2012 出版版是章內數值的適用來源。學生清單已移除未標年份與持續維護字樣；未標日期的網站不另外補造年份。
