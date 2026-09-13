# 地震統計與預報：從一份義大利實驗到臺灣觀測

2016 年 8 月 24 日凌晨，義大利中部的 Amatrice 發生 Mw 6.2 地震。接下來兩個月，
同一區域又出現三次規模 5.5 以上的地震。事後回頭看，這一串地震落在哪裡、
間隔多久，都有跡可循。真正的問題是：**在它們發生之前**，一份寫好的預報
能說出多少？

這個網站帶你親手回答這個問題。你會用一份公開的義大利地震目錄，
訂出一張「預報實驗規格卡」，做出四張不同模型的預報，
再用 2012–2021 年真正發生的目標地震去檢驗它們。
你只需要修過一門基礎統計；地震學、點過程與程式都在路上補齊。

## 第一部：一份地震預報實驗，從資料到檢驗（第 1–20 章）

第一部沿著一份實驗的順序推進。每一章只做一件事，也只在需要時才引入新名詞。

| 段落 | 章 | 你會學到 |
|---|---|---|
| 題目 | 1–3 | 預報長什麼樣、怎麼讀一份目錄、規格卡的每一欄 |
| 統計工具箱 | 4–6 | Poisson 過程與第一張預報、概似與估計、模擬式檢定與分數 |
| 讀目錄 | 7–8 | 大小地震的比例（b 值）、目錄從哪個規模開始可信（Mc） |
| 模型 | 9–15 | PPE、叢集律、Hawkes 與條件強度、ETAS、Ψ 現象、EEPAS |
| 檢驗 | 16–18 | 數量、位置與規模、相對資訊：三種問法 |
| 之後 | 19–20 | 把模型加起來；讀一篇論文、進入決策 |

同一張規格卡會一章一章填滿。模型章先讀預報形狀，第 16 章起再比較分數。
這是教學順序；實驗是否前瞻，取決於資訊隔離與設定何時固定。
本站使用修訂目錄做回溯性的擬前瞻重演，並非當年的即時發報。
從{doc}`01_forecast_question`開始。

## 第二部：臺灣地球物理觀測（第 21–29 章）

第二部回到臺灣，從地下水、地磁、地動與 GNSS 的儀器紀錄出發，
用第一部學到的統計判斷去讀真實觀測。2024 年花蓮地震的案例把不同觀測放在同一條
時間軸上；最後一章把臺灣地震目錄的統計面貌、與「義大利規格卡換成臺灣版要改哪些欄」
一起整理。已熟悉統計模型的讀者可以直接從{doc}`01_overview`進入。

## 第三部：技術附錄

主文保留動機、圖與必要公式；推導與演算法集中在附錄，正文會連到對應段落。

| 想追的內容 | 入口 |
|---|---|
| 記號表、兩套術語對照、統計工具推導、點過程 | {doc}`appendix_a_point_process` |
| HORUS 取得、完整度、叢集律 | {doc}`appendix_b_catalog` |
| ETAS 全式、simplETAS 釘參理由 | {doc}`appendix_c_etas` |
| Ψ 與 EEPAS 的核、η 與 Δ | {doc}`appendix_d_eepas` |
| 檢驗、十欄格式與重新分格、組合 | {doc}`appendix_e_testing` |
| 危害與複發模型 | {doc}`appendix_f_hazard` |
| 觀測資料處理與重現指南 | {doc}`appendix_g_observations` |

原始 notebook 保留產生每張圖的程式。閱讀正文不需要安裝環境；
想重現時，依 `README.md` 的資料準備步驟取得 HORUS 目錄與觀測資料，
並遵守資料提供者的使用規範。

## 參考資料與延伸閱讀

- Biondini, E., Rhoades, D. A. 與 Gasperini, P.（2023），[Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123)，*Geophysical Journal International* 234, 1681–1700，開放取用。本站義大利實驗的設定、參數與檢驗結果都以這篇為準；建議先讀 APPLICATION TO ITALY 一節，再對照第 2、3 章的規格卡。
- Mizrahi, L. 等（2024），[Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823)，*Reviews of Geophysics*，開放取用。先讀 Plain Language Summary，看模型開發、檢驗與溝通如何串成一個流程。
- USGS，[Can you predict earthquakes?](https://www.usgs.gov/faqs/can-you-predict-earthquakes)，免費官方說明。釐清「預測特定地震」與「機率預報」的差別。
- 蕭乃祺（2019），[臺灣地震測報的發展](https://www.ntsec.edu.tw/liveSupply/detail.aspx?a=6829&cat=6841&lid=16154&p=1&print=1)，《科學研習》58 卷 6 期，免費文章。認識第二部觀測系統的建置背景。
