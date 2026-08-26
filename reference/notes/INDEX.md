# reference/notes/ 索引

`reference/` 下約 40 篇論文的結構化中文筆記，按主題分成 9 個檔案，供第二部
「地震預報模型」各章撰寫使用。每篇筆記含：書目、模型／方法概述、關鍵觀念
與公式（LaTeX）、教學洞見、與台灣的關聯。

**⚠️ Embargo**：`[2026] PyEEPAS`(+sup) 與 `[2026] Taiwan_EEPAS`(+sup) 為審稿中
論文，對應筆記節僅含方法背景，其任何數字、參數、圖表、結果**不得出現於教學
網站**。`[2025] EEPAS_Software`（submitted to GJI）與 EEPAS_TW-main.zip 的程式
與結果比照辦理。

**檔案勘誤（ingest 時發現）**：
- `[2007] Rhoades.pdf` ＝ `[2007] Application...CA.pdf`（md5 相同，取一份）
- `[2023] Italy_EEPAS.pdf` ＝ `[2023] Application...Italy.pdf`（md5 相同，取一份）
- `[2023] Calibrated ETAS.pdf` 內容實為 Mizrahi et al. 2021（SRL declustering 論文），與 `[2021] Estimate b.pdf` 同篇
- `Taiwan/[2011] Omori Law Taiwan.pdf` 與 `Taiwan/[2012] New Empirical Tests...pdf` 為同一篇（Tsai, Ouillon & Sornette 2012 BSSA）
- `[2024] Italy_exp.pdf` 非 EEPAS 論文，實為 Brehmer et al. 2024 預報評分方法論
- `[2024] soft_intro.pdf` 非軟體導論，實為 Wein et al. 2024 餘震預報溝通研究
- `entropy-22-01264.pdf` 的 22 是卷號，實際年份 2020（Rhoades et al., EEPAS 前置時間）

## 主題檔案 → 章節對應

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
