# 第一部讀者審查與五項修正核對

## 範圍

讀者人設：修過基礎統計，沒有地震學背景。初次審查已通讀新第一部第 1–20 章，重點是術語首次定義、概念順序、圖文與數字的一致性、未核實斷言，以及教材是否交付已定閱讀示範。收尾只核對原先列出的五項，不擴大通讀。

這是**讀者審查**，不是全部公式、原始論文、notebook 執行或瀏覽器驗收。初次審查僅以 CSV、CELLE 與投影作一次有界的目標事件數核對；沒有執行 notebook。收尾沒有重算實驗、沒有修改章節或 toolkit、沒有啟動 Chromium。Fig. 8 的確認只包含圖片檔存在、內容可讀及章節引用，不代表整頁手機渲染已驗收。

## 五項發現與現況

| # | 原始問題 | 修正核對與現況來源 | 結果 |
|---|---|---|---|
| 1 | 第 4 章首次 Poisson 定義出現 `Lambda`、`lambda`、`rac{...}`，遺失 LaTeX 反斜線；第 12 章背景記號定義與公式不同，章末 `nu` 斷行。 | `book/04_poisson_and_sup.py:69`、`:73`、`:110` 已恢復 `\Lambda`、`\lambda`、`\frac`。`book/12_etas_structure.py:74` 定義 `\mu_0(x,y)`，與 `:78` 公式一致；`:335` 已恢復 `\nu`。 | **確認無誤：原符號與來源格式缺陷已修正。** 未以本項宣稱全書 MathJax 渲染通過。 |
| 2 | 第 3 章說四模型都更新到發報時刻，只有 EEPAS 有 50 天 delay，與 PPE 章及實作矛盾；d20 圖亦把本站 PPE 誤畫成固定形狀。 | `book/03_experiment_spec.py:111` 已區分 SUP 固定學習率與其他模型更新歷史；`:115` 明列 PPE／EEPAS 都有 50 天延遲、ETAS 沒有。d20 generator `book/_static/diagrams/sources/build_teaching_diagrams.py:151` 及生成片段明列 PPE 緩慢更新與 50 天延遲，沒有把時間尺度宣告成家族必然範圍。 | **確認無誤：共同規格與圖中的矛盾已消除。** |
| 3 | 第 3 章從 25 vs 27 的差額直接斷言兩顆地震降到門檻以下或移出 R，卻沒有逐事件版本對照。 | `book/03_experiment_spec.py:77` 保留兩份事件數；`:79–81` 明寫尚未取得逐事件對照，不能斷定是哪兩顆，也不能將差額全部歸因單一欄位。 | **確認無誤：已將未核實原因與已知計數分開。** |
| 4 | 第 2 章把年計數增加全部歸因儀器，並用曲線平坦解釋完整度，讓描述圖被誤當成原因證據。 | `book/02_reading_horus.py:118` 改為「尋找完整度的線索」；`:155–156` 明說圖不能單獨分開觀測能力與真實活動；`:163–165` 說明年代門檻來自文獻、完整目錄仍可有強烈率變化。 | **確認無誤：描述與因果證據已分開。** |
| 5 | 第 20 章只摘述 Fig. 8，沒有交付「一張來源結果圖加旁註」的閱讀示範。 | `book/20_beyond_forecast.py:83` 已嵌入 `book/_static/reading/biondini2023_fig8.png`，`:85–87` 列圖號、頁碼、來源與區間；後續保留 IGPA／IGPE 分母差異及逐句旁註。已直接開啟圖片，確認五個窗長面板、模型標籤、IGPA 軸、零線與信賴線可辨讀。 | **確認無誤：缺少來源結果圖的問題已修正。** 圖中數字未重新估讀。 |

以上五項均已關閉；未以這次有界收尾建立新問題清單。

## 初次審查的確認無誤項

- 直接讀取 `data/cache/italy/horus_clean.csv`，依地震事件、CPTI15、深度不超過 40 km、取整規模 `[5.0,7.5)`、2012–2021 年及 CELLE 測試區篩選，使用 EPSG:7794 投影核對，得到測試期 **25 顆**。Emilia 2012-05-20 至 2012-05-29 為 **7 顆**，中義大利 2016-08-24 至 2017-01-31 為 **10 顆**。
- 原始規模低於 5、取整後成為目標的三筆為：2012-01-25 的 4.98、2013-01-25 的 4.98、2016-11-01 的 4.95。與第 3 章列出的摘要一致。本項沒有核實原論文 27 顆各自的身分。
- CSEP、PPE、ETAS、EEPAS、POLL／jPOLL／BILL 有全名；概似、MLE、bootstrap、分位數與 p 值已安排在後續估計及檢驗之前。
- 第 11 章區分發報期望數與沿實際歷史的補償子；第 13 章明列估計示範與快取預報的參數及門檻差異。
- 第 18 章說明配對不等於事件獨立；第 19 章列出總期望數不變的條件，亦說明凹性不保證勝過最佳成分。
- 在新第一部來源清單中，未發現台灣資料計算殘留；台灣僅出現在第 20 章轉場規格卡。`PyEEPAS` 只命中公開 GitHub 連結。本項不代替全站 embargo 終檢。

## 圖檔的獨立靜態驗證

七張新增片段及七張同名離線預覽，正式 diagram-design `self_check.py` **14／14 通過**；七張預覽的 `verify-geometry.py` **7／7 零 finding**。靜態集合共有 **66 個唯一 ID**，ARIA 參照完整、每片段一個 SVG、沒有外部依賴。收尾已按現有圖檔刷新 SHA-256；未重建圖檔，也未修改 archify 的已交付 JSON／HTML。

收據：

- `book/_static/diagrams/standalone/checks/new_diagrams_self_check.json`
- `book/_static/diagrams/standalone/checks/new_diagrams_geometry.json`
- `book/_static/diagrams/standalone/checks/new_diagrams_static.json`

瀏覽器與全書建置由主代理另行處理；本報告不宣稱其結果。

## 收尾快照

核對時間（UTC）：2026-09-13T10:55:10.626954+00:00。檔案若後續再修改，以下雜湊可辨識本次核對的版本。

| 檔案 | SHA-256 |
|---|---|
| `book/02_reading_horus.py` | `922f238bb5eeb0d75d3f43697bb187fca3b7a766499a155991789eac793e25ee` |
| `book/03_experiment_spec.py` | `30b6bb46217515d7928afdc089bc136e9affa815bd58ff3f1f278c134e0b2edb` |
| `book/04_poisson_and_sup.py` | `4b7c182406552fe670b5819797caa199c2f9802397d0075d61321b02f3c29a35` |
| `book/12_etas_structure.py` | `e5fa235d2f649e1fc2302985a363c7d42352b318aef2692c642c0c944b40a8fc` |
| `book/20_beyond_forecast.py` | `b160a6215aa1024ccbcdad15da5d44e1fad9bd09c4fc0a03d6164f832686b1e8` |
| `book/_static/reading/biondini2023_fig8.png` | `3c58d134f1712e97620b2c5d6b8bcb62814adc587c67847b52ca5748becd4e3a` |
| `book/_static/diagrams/d20_model_family.html` | `159532a1317a5b87e99dd8b6d17b75f616d3d358cfb2b866db8c7bd1331c9b11` |
| `book/_static/diagrams/sources/build_teaching_diagrams.py` | `bdceaec18eb16dcc2a8b3497592054dc9501493a819a7574e276646c9bd5bd11` |

## 瀏覽器後續修正：觸控與 hover（2026-09-13T11:27:00.482837+00:00）

主代理在原讀者審查之後，以實際互動檢查發現 d09 第二選項雖已 checked，說明仍可能被滑鼠殘留的 `:hover` 規則遮住。這是瀏覽器後續發現，不是前述讀者審查或靜態檢查已涵蓋的結果。

主代理已將 hover 規則限定於 `@media (hover:hover) and (pointer:fine)`，移除非選中 SVG 的半透明設定，並重生十二張片段與預覽。本子工作只對重生後的檔案刷新正式靜態驗證：self-check 24／24 通過、geometry 12／12 零 finding；十二圖的 112 個 ID 唯一、ARIA 參照完整、沒有外部依賴。沒有修改 generator／章節，也沒有啟動 browser。

目前收據為 `book/_static/diagrams/standalone/checks/all_twelve_self_check.json`、`all_twelve_geometry.json`、`all_twelve_mobile_static.json`，七圖及五圖的分組收據也已同步刷新。舊静態收據保留在 `book/_static/diagrams/standalone/checks/history/20260913T112700Z` 並標示 historical_only。本報告前節的舊 generator／d20 雜湊保留為原讀者核對快照；現況雜湊以新收據為準。

**觸控 bug 是否已在實際瀏覽器消除，由主代理重跑互動及含圖頁面驗收判定；本次靜態通過不替代該結果。**
