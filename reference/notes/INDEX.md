# reference/notes/ 索引

目前教材採20章義大利主線、第二部21–29章與附錄。新版來源入口：

- [逐來源概念覆蓋與核對深度](refresh_20260913_exhibition/literature_coverage.md)
- [機器可讀覆蓋資料](refresh_20260913_exhibition/coverage.json)

新版登錄涵蓋59個實際PDF檔案，包括2組二進位重複與3組同篇不同版本。
**登錄不等於59篇全文已核讀。** 部分只核對首頁；ETAS/[1998] ETAS.pdf雖被抽取工具標成ok，快取實際只有分頁字元，原文仍未核。
原始PDF與全文文字快取位於被忽略的本地資料，不公開到網站。

**Embargo**：2025 EEPAS_Software、2026 PyEEPAS（含補充）、2026 Taiwan_EEPAS（含補充）仍按審稿中處理。
不得公開這些稿件的方法、數字、參數、圖表或研究結果；台灣EEPAS正文只能說「在地化工作正在進行中」。
舊筆記含有稿件摘要與數值，不代表它們已獲准成為教材來源。
可用的是已發表EEPAS文獻、Biondini 2023義大利設定/參數，以及明確公開的PyEEPAS GitHub程式與義大利資料；
公開repo不等於解禁稿件結果，也不授權公開EEPAS_TW內部內容。

**已核實的錯名／版本**：
- 2007兩份CA EEPAS、2023兩份Italy EEPAS為二進位相同檔。
- Calibrated ETAS實為Mizrahi等2021除叢與規模分布論文；Estimate b為其2020作者稿。
- 兩份台灣Multifractal Omori為Tsai等2012同篇不同PDF。
- 根目錄2022 Tests與Ensemble同名研究為Bayona等2022同篇版本。
- ETAS_R真正年份2019；Ncom為Girona–Drymoni2024異常低規模地震研究。
- Italy_exp本地是2024預印本評分研究；soft_intro是2024預報溝通研究。
- entropy-22真正年份2020；22是卷號。完整書目與核對範圍見新版表。

下方保留早期九份主題筆記與歷史章號，供查找素材；不作目前章序或公開授權的依據。

## 主題檔案 → 章節對應

2026-09-09 全站延伸閱讀已按目前第 01–23 章重新核對，含原始 PDF 與公開連結：
[首頁與第一部](reading_00_08_verified.md)、[第 09–16 章](reading_09_16_verified.md)、
[第 17–23 章](reading_17_23_verified.md)。以下舊筆記的供應章號保留原記錄；
定位目前教材先用新版覆蓋表，再參照 `book/_toc.yml`、`book/_conventions.md` 與上述核實稿為準。

| 檔案 | 內容 | 主要供應章節 |
|---|---|---|
| [eepas.md](eepas.md) | EEPAS 原始論文（Rhoades & Evison 2004）到 20 年回顧、CA/日本/義大利應用、EAS 餘震擴充 | 12（EEPAS 與 PPE）、14 |
| [etas.md](etas.md) | 時空 ETAS（Ogata & Zhuang 2006）、R 實作與估計、simplETAS、declustering 效應 | 11（ETAS）、10 |
| [testing.md](testing.md) | N/M/S/L 一致性檢驗、T-test、Molchan/ROC、統計功效、pyCSEP、Quadtree 網格 | 15（CSEP 檢驗） |
| [stats.md](stats.md) | b 值估計（binned 修正、b-positive、declustering 選擇效應）、Mc、時變點過程回顧、SeismoStats | 10（目錄統計進階） |
| [ensemble.md](ensemble.md) | 凸組合混合（STEP+EEPAS）、multiplicative hybrid、十年前瞻測試教訓、權重最佳化 | 14（ensemble） |
| [taiwan.md](taiwan.md) | 台灣 Mc/b 值/Omori p/Båth、均一化目錄、池上與大埔序列、CWA 速報預警史、本土化 ETAS | 10、17（台灣） |
| [oef.md](oef.md) | OEF 全景回顧（STEP 結構、義/紐/美系統、Delphi 共識）、機率溝通、時變危害、SEDA | 09（導論）、13（STEP/OEF）、16 |
| [psi.md](psi.md) | Ψ 前兆尺度增加現象、自動辨識演算法、時空抵換線、hindsight 偏誤 | 12 |
| [psha_step_websearch.md](psha_step_websearch.md) | Baker 白皮書五步驟、TEM PSHA 2015/2020、Gerstenberger 2005/R&J 1989 書目核實、GeoNet/USGS 現行系統 | 13、16（PSHA） |

## 論文 → 筆記對照

### 根目錄
| 論文 | 筆記 |
|---|---|
| [2004] rhoades2004 | eepas.md |
| [2006] Space–time ETAS models | etas.md |
| [2007] Application of the EEPAS Model to CA（=[2007] Rhoades） | eepas.md |
| [2009] Long-range forecasting allowing for aftershocks | eepas.md |
| [2011] Japan | eepas.md |
| [2017] ETAS_R（實際 2019, JSS） | etas.md |
| [2017] SEDA | oef.md |
| [2019] Brief Review Process | stats.md |
| [2021] Estimate b（Mizrahi declustering） | stats.md |
| [2022] 20-Year Journey EEPAS | eepas.md |
| [2022] Pycesp | testing.md |
| [2022] Tests（Bayona GJI） | testing.md、ensemble.md |
| [2023] EEPAS Italy（=[2023] Italy_EEPAS） | eepas.md |
| [2023] Calibrated ETAS（實為 Mizrahi 2021） | etas.md |
| [2023] SimpleETAS | etas.md |
| [2023] Stat_power_test | testing.md |
| [2023] Tests（Bayona TSR） | testing.md |
| [2023] open, transdisciplinary | oef.md |
| [2024] Estimate b（Tinti & Gasperini） | stats.md |
| [2024] Italy_exp（Brehmer 評分方法論） | eepas.md（附警示） |
| [2024] Ncom（Girona & Drymoni） | oef.md |
| [2024] New Pycsep | testing.md |
| [2024] OEF_Review（Mizrahi RoG） | oef.md |
| [2024] soft_intro（Wein 溝通研究） | oef.md |
| [2025] EEPAS_Software（submitted） | eepas.md |
| [2025] SeismosStats | stats.md |
| [2026] PyEEPAS + sup ⚠️embargo | eepas.md（僅方法背景） |
| entropy-22-01264（2020, EEPAS 前置時間） | testing.md |

### Ensemble/
四篇（2009 Mixture、2014 Multiplicative、2022 Testing CA、2023 Maximizing skill）→ ensemble.md

### Psi regression/
兩篇（2021 Space-Time Tradeoff、2024 psi）→ psi.md

### Taiwan/
全部 → taiwan.md（[2026] Taiwan_EEPAS + sup ⚠️embargo 僅方法背景；
AutoBats 為 workshop 投影片、112 年報告為政府委辦報告，引用時標明性質）
