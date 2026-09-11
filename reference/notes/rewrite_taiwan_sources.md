# 台灣觀測篇改寫與來源核對（2026-09-11）

範圍：01–08、23 主文與新附錄 G。舊數字章號保留供主端統一編號；章間改以文件標題連結。
未使用 embargo 材料，未執行下載、數值實驗、build 或 push。

## 連貫的主文

順序為測網與觀測機制 → 資料產品與品質 → 地下水 → 地磁 → 波形／目錄 → GNSS → 花蓮共同案例 → 證據設計 → 台灣展望。
前置統計知識改為已讀過模型篇；08 不再引出另一個「第二部」，23 改接08並收束主文。
各章不套固定「誤解清單／研究前沿／一句話」模板；技術細節另存 G。

## 核實來源與使用範圍

- 原有章末書目與導讀參考 `reading_00_08_verified.md`、`reading_verification_20260909.md`；此改寫沒有宣稱重抓全部既有文獻全文。
- USGS 地下水 FAQ 已讀官方內容：震波造成振盪、震後偏移與可能的滲透性改變。https://www.usgs.gov/faqs/how-does-earthquake-affect-groundwater-levels-and-water-quality-wells
- Sneed、Galloway、Cunningham（2003）官方出版頁核實遠近地震皆可能造成水文反應與書目：doi:10.3133/fs09603。舊 pubs 入口直接抓取403，官方出版頁可读： https://www.usgs.gov/publications/earthquakes-rattling-earths-plumbing-system
- USGS Geomagnetism Overview 已讀：對指定已報告前兆的重新檢查；不泛化成所有電磁效應不存在。https://www.usgs.gov/programs/geomagnetism/science/overview
- IGS RINEX 官方入口、NGL Plug and Play 產品頁及 Blewitt et al.（2018）Eos 頁已讀：產品、參考框架、日解／較高頻解與設備階變資料的區別。未新增當前測站數或服務延遲。
- ESA Navipedia，Sanz Subirana et al.（2011）已讀式(5)–(6)，G 的 ECEF→ENU 矩陣相符。正向大地座標公式與共變異數傳播另作幾何／線性代數說明。
- SciPy `signal.welch` 官方文件已讀：預設 density，輸出單位為輸入量平方/Hz；支持修正水位 PSD 指標名稱。
- FDSN Channel codes 舊連結本次 web 抓取失敗，沿用既有核對書目；本次不聲稱實際 HTTP 成功。頻道解讀縮限，刪除「加速度儀不會爆表」。
- Zheng et al.（2024）Nature/Scientific Reports 官方頁已讀，發表日期2024-12-28，核對 combined geodetic/seismic inversion；主文不由餘震點雲直接認定破裂面。
- `reference/Taiwan/[2015] b-Values Observations in Taiwan A Review.pdf` 首頁書目核對，歷史資料背景只保留測網／規模定義變更，刪除未核實精確現行站數。
- `reference/Taiwan/[2016] An updated and refined catalog of earthquakes in Taiwan (1900–2014) with homogenized M w magnitudes.pdf` 摘要及背景核對：跨年代規模均一化、來源優先序。正文不再抄三套可能適用範圍不同的換算式。
- `reference/Taiwan/[2013] Taiwan Bath.pdf` §2、PDF p.3 核對：當時陸上／外海完整度差異；明確區分本文資料與教材MAXC示範。
- `reference/Taiwan/[2025] Fast_Report_Dapu_ETAS.pdf` §3.2、Table2、結論核對，PDF pp.5–6。原Table2第三列為4Days，正文與結論為7days，教材依正文保留7天並註記不一致，不擅自變更機率／觀測值。
- `reference/[2022] A 20-Year Journey of Forecasting with the EETAS model.pdf` 首頁與PPE段落核對：PPE背景空間分布是過往震央鄰近性的模型。
- `reference/[2024] OEF_Review.pdf` 首頁與摘要核對：基準、前瞻、透明度、使用者共同設計；新增為08、23章末來源，不以其2024描述認證目前台灣制度。

## 重要解讀修正

地下水不普遍等於直接孔隙壓力量測；水位公分不等於地殼位移公分。
PSD峰值平方根為cm/√Hz，cm/hPa氣壓斜率不直接等於無量綱氣壓效率。
兩站地磁相減不保證去除全部外源場；低殘差也可能減掉目標。
一秒取樣Nyquist為0.5Hz；一分鐘為1/120Hz，沒有宣稱整個ULF頻帶都保留或都消失。
波形counts未經校正不稱速度；單站到時差近似距離不直接稱震央距離。
GNSS近似座標不是精密時序，日解不能辨認秒級过程，慢滑移也非只有GNSS能觀察。
花蓮未辨認出訊號限定圖、站、窗與方法；差分std使用震前10分鐘／震後5分鐘，如實說明，未聲稱已做效應檢定。
08刪除高誤報即毫無價值、疊加保證消除巧合等說法。
23刪除静態模型資訊量為零、ETAS只能震後有用、binary計分消除相依、時間獨立模型必須對除叢資料評分等敘述。

## Code cells 與快取整合

全部 code cells 原順序與數量保留，沒有搬到另一頁，沒有刪除數值計算。
以 AST 將字串常數置空後逐cell比較HEAD與新檔，程式結構一致；變動僅metadata、圖說／軸標／print文字及docstring。
下表用「該頁第幾個code cell」計數，包含setup cell，不含Markdown：

| 檔案 | code cells | 變更序號 |
|---|---:|---|
| `01_overview.py` | 5 | 無 |
| `02_download.py` | 4 | 1, 2, 3, 4 |
| `03_groundwater.py` | 9 | 3, 5, 7, 8, 9 |
| `04_geomagnetic.py` | 8 | 無 |
| `05_seismic.py` | 11 | 4, 8, 10 |
| `06_gnss.py` | 5 | 2, 3, 5 |
| `07_case_hualien2024.py` | 9 | 3 |
| `23_taiwan_outlook.py` | 6 | 2, 6 |

- 02 全四cell為登入／申請／下載清單／頻道回顯，新增remove-output及remove-input；程式原文保留。
- 03 第3cell的info/head回顯新增remove-output；第5cell圖名改半日尺度起伏；第7cell半日頻帶指標及docstring；第8cell氣壓響應斜率；第9cell水位缺測圖名。數值未動。
- 05 第4cell改未校正counts軸名，第8cell修正符號大小敘述，第10cell改主震後事件數圖名／軸名，未改計算。
- 06 第2、3、5cell原始檔名、檔頭與座標回顯新增remove-output及remove-input。
- 07 第3cell圖名改主震後事件分布，避免親代分類宣稱。
- 23 第2cell儀器沿革annotation縮限；第6cell時間尺度圖的世界排名／就緒／缺口現況字串改用途描述，數值區間原留並明說只是示意。

整合時需要對上述字串變更更新既有Plotly/text輸出；不應盲目保留舊圖標，也不必重跑下載或實驗。

## 已做與未做

已做：8份py的AST解析、code cell數量與非字串結構保持、{doc}目標存在、每頁單一章末閱讀、控制字元與git diff --check。
未做：ipynb同步、Jupyter Book建置、瀏覽器渲染、資料重算；由主端統一整合與驗證。部分既有付費文獻只沿用專案先前的書目核對與閱讀定位，沒有宣稱本輪讀到全文。
