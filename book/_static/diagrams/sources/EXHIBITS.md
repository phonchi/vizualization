# 科學展件：設計與重現

## 已接受的設計方向

模式：Redesign · Overhaul。讀者有基礎統計，沒有地震學背景；以科學展覽的寬幅舞台、敘事動畫與資料探索理解概念。保留 `.teaching-diagram`／`show_diagram(slug)`、原 12 個 slug、實驗資料與已有科學界線；新增 Hawkes 與預報實驗協定兩展件。

Design Read：資料驅動的互動科學展覽；visual variance 8、motion 6、density 5、asset dependence 7（真實義大利區域）、brand fidelity 8。使用宿主海軍藍／礦石白／青綠 token，珊瑚紅表示目標或示意事件。圖形採曲線、地圖、時間帶、直方圖、方格、樹及矩陣等不同構圖。文字沿宿主字型，圖中標籤按實際容器寬度重新排版，不把桌機畫布整張縮小到手機。

文字敘事、圖形層級與資料移動優先；不加裝飾性地震波或假傳播。全部模型參數與示例固定，滑桿只移動資料讀取位置、時間或敘事進度，不做模型重估、門檻調參或新預報實驗。

採用 web-design-engineer 的視覺與響應式流程、diagram-design 的圖型及 accessible SVG 原則。新交互控制依使用者明確要求放入全站獨立 JS；不是舊的 canonical motion widget。舊 archify JSON／deliver HTML 與 receipt 保留為歷史設計證據，不冒稱它們驗收了目前的互動展件。

## 元件契約

- 根元素：`.quake-exhibit[data-exhibit="slug"][data-kind="kind"][data-config="JSON"]`，放在既有 `.teaching-diagram` 內。
- 由 `_config.yml` 全域載入 `_static/exhibits.css` 及 `_static/exhibits.js`；片段沒有 inline script。standalone 預覽自包含 CSS／JS，不需要伺服器或網路。
- `window.EarthquakeExhibits.mount(scope=document)` 可重複呼叫；DOMContentLoaded 與 MutationObserver 自動處理後插入的片段。
- `pauseAll()` 暫停全部播放。每張圖各自保存 state、RAF、ResizeObserver、IntersectionObserver；重複 slug 的 DOM 實例仍有不同 SVG ID。
- 控制：`[data-action="play|reset|boundary"]`、`[data-control="progress|compare"]`、`[data-layer="S|R|events|contours|sources"]`。
- 每圖原生 range、select、button、checkbox 供鍵盤與觸控操作；沒有 hover 蓋過 checked 的 CSS 規則。
- 播放由使用者啟動；畫面離開 viewport／分頁隱藏會暫停，不自動重啟。reduced-motion 初始為完整靜態畫面，播放鍵改為單步前進，仍能手動讀取資料。
- `EarthquakeExhibits.render(kind,p,state,data,width)` 與 `svgHTML(...)` 是無 DOM 的純函式，用於重現靜態第一畫面及有界檢查。瀏覽器依舞台實際寬度重新計算幾何。

## 展件與操作

| slug | 圖型／真正改變的內容 | 固定資料 |
|---|---|---|
| d01_forecast_pipeline | 分支工作流程，播放時沿目錄→規格→四模型→預報→檢驗前進 | 實驗流程 |
| d02_regions_s_r | 真實地圖：S／177 格／震央圖層、歷史截止年 | CPTI15、CELLE、HORUS 128 顆 M≥5 |
| d03_timeline_windows | 40 窗 scrubber、目標點與發報／延遲時間帶 | 真實 40 窗與 25 目標 |
| d03_thresholds | 規模箱沿原始數值移動，2.45／4.95 邊界跳點 | 取整公式、門檻固定 |
| d04_density_to_counts | 36 個密度小格逐格積分、Poisson 機率柱隨累積值改變 | 固定合成密度 |
| d06_quantile_pvalue | 直方圖、階梯 CDF、含同值的上下尾／雙尾 | 20 次合成結果 `[1,2,4,6,4,2,1]` |
| d09_kernel_smoothing | 等值圈＋可移動截線＋對應一維核總和 | 兩個固定 Gaussian 核 |
| d11_hawkes_process | 時間 raster、已知歷史、各核與總率逐項疊加 | 時間 `[2,5,5.8]`、背景 0.3、核 `1.5 exp(-u/1.3)` |
| d12_etas_branching | 逐代揭露親子樹 | 固定七節點家族；不估分支比 |
| d15_eepas_kernels | 規模／時間／空間核到乘積格箱 | 固定概念形狀；不是義大利重算 |
| d16_csep_tests_flow | N／S／M／cL／L 條件矩陣＋示意抽樣點 | 固定示例；不是正式分數 |
| d19_convex_mix | 兩成分的加權柱累積 | A=(2,6,1,3)、B=(6,2,5,3)、權重各½；總數14 |
| d20_model_family | 可選焦點的模型關係圖 | 固定資料使用關係 |
| d03_forecast_protocols | 回溯／擬前瞻／真前瞻比較時序、受控延遲與部分前瞻視圖 | 規則與資料可用性的概念示意 |

## 真實義大利資料供共用

`../data/italy_exhibit_data.json`：EPSG:7794 投影公里；`cells` 為 177 列 `[x0,x1,y0,y1]`，`polygon` 為 15 個 `[x,y]` 頂點；`events` 為 128 列 `[x,y,year,mb,timeISO,cellIndex]`（`-1` 代表 R 外）；`targets` 為 25 列同格式；`windows` 為 40 列 `{index,start,end}`。

所有時間皆為帶 `Z` 的 UTC ISO，避免瀏覽器本地時區改變落窗位置。這只是從既有 CSV／MAT 的輕量匯出，未擬合模型或產生預報。來源 SHA-256 在 `exhibit_data_provenance.json`。地圖只畫真實資料中的區域與格子，沒有偽造海岸線。首頁可共用同一個 JSON。

## 重建與靜態檢查

```bash
.venv/bin/python -B book/_static/diagrams/sources/build_exhibits.py
node book/_static/diagrams/sources/emit_exhibits.cjs
node book/_static/diagrams/sources/check_exhibits.cjs
```

舊 `build_teaching_diagrams.py` 是先前窄版 radio 圖的歷史產生器，不用來重建目前展件。

目前收據：`../standalone/checks/exhibits_manifest.json`、`exhibits_static.json`、`exhibit_fragments_self_check.json`。Pure renderer 檢查 14 種圖、四尺寸、五進度及比較模式，共 520 組；另外核對真實資料量、Hawkes 背景、固定混合總數、離散上尾及 4.95 邊界。十四片段通過正式 diagram-design accessible SVG／單檔安全 self-check。

這些檢查沒有啟動瀏覽器，不保證視覺沒有碰撞；實際互動、字型、跨視窗與全書內嵌驗收由主代理集中執行。standalone 含新全域控制器，不宣稱通過舊 canonical-motion-controller 限制。舊圖的 geometry／browser pass 只適用其對應 hash，不適用目前展件。

## protocol 語意核對

依 `book/03_experiment_spec.py` 的檢驗方式段落及 `reference/notes/refresh_20260913_exhibition/literature_coverage.md` 核對：回溯可為樣本內或樣本外；擬前瞻是回溯子類。受控延遲前瞻須預先鎖定流程／規格，並隔離模型及建模者對測試資料的接觸，不是一般晚點評分。部分前瞻另以已知／未知時段視圖呈現。本站採修訂 HORUS 的事件時間重演，沒有宣稱逐版還原當時目錄 availability。此段只核對展件與已核讀主文的一致性，未新增全文文獻驗證聲明。

Hawkes 展件另核對 t=2、5、5.8 的精確左極限：歷史點、caption 與率都使用 `event_time < t`。細線只畫每顆事件的核貢獻，背景 0.3 另以虛線顯示；粗線為背景加所有核，沒有重複背景。靜態收據包含邊界歷史數與單事件的分解，並核對 SVG 細線採零基準。
