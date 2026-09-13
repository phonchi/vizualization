# 互動示意圖清單（2026-09-13 重組）

每個檔案是自包含 HTML 片段（inline SVG＋`<style>`＋少量 `<script>`，無外部資源），
在章節中以 `from gdms_toolkit.viz import show_diagram; show_diagram("d02_regions_s_r")` 嵌入。
規格：`viewBox` 設定、`max-width:100%`、字級 ≥ 12px、只用 hover 說明與逐步揭露；
顏色一律用 CSS 變數（`--dg-ink`、`--dg-muted`、`--dg-accent`、`--dg-paper`），
深色模式由 `teaching.css` 依主題屬性覆寫變數。id 加檔名前綴避免同頁衝突。

| 檔名 | 用於章 | 內容 |
|---|---|---|
| `d01_forecast_pipeline` | 1、3 | 預報實驗流程：目錄→規格→模型→預報→檢驗（archify workflow） |
| `d02_regions_s_r` | 2 | 收集區 S ⊃ 測試區 R、邊界效應（hover 各區） |
| `d03_timeline_windows` | 3 | warm-up／learning／testing、40 個滾動窗、發報時刻與資料截止（逐步揭露） |
| `d03_thresholds` | 3、8 | M_c／m_0／m_T 三門檻在規模軸上的位置 |
| `d04_density_to_counts` | 4 | 率密度→積分→格內期望數→Poisson 計數分布 |
| `d06_quantile_pvalue` | 6、16 | 統計量→模擬分布→分位數分數→p 值 |
| `d09_kernel_smoothing` | 9 | 幾個震央的核疊加成一張率面 |
| `d12_etas_branching` | 12 | 背景事件與世代分支樹、分支比 |
| `d15_eepas_kernels` | 15 | 一顆輸入事件的規模核、時間核、空間核如何形成貢獻 |
| `d16_csep_tests_flow` | 16、17 | N／S／M／cL／L 各固定什麼、模擬什麼 |
| `d19_convex_mix` | 19 | 兩張預報凸組合、期望數守恆 |
| `d20_model_family` | 1、20 | 模型家族：時間獨立（SUP、PPE）／短期（ETAS）／中期（EEPAS）／混合 |

## 續作成果與驗證

12 張列出的嵌入圖均已存在。2026-09-13 補完的 7 張採用原生 radio＋CSS 的說明切換，
不需要 JavaScript；滑鼠 hover、鍵盤和觸控皆有操作入口。兩張 archify 流程圖另保留
JSON、獨立 Viewer 與正式交付收據。來源與重建方式見 [`sources/README.md`](sources/README.md)。

新增圖的 accessible SVG／單檔安全與靜態 geometry 檢查已通過；
實際瀏覽器、手機互動及 Jupyter Book 嵌入驗收由主代理統一執行，不能以靜態檢查代替。

## 390 px 手機修正

原五圖已於手機驗收後重排為 480 viewBox／20 px SVG 標籤，與新七圖使用相同 HTML 說明與原生 radio 操作。
圖的完整內容不裁切，也不提供整圖橫向捲動。五圖的新靜態檢查收據以 `legacy_five_mobile_*` 命名；
十二圖現況雜湊在 `standalone/checks/all_twelve_mobile_static.json`。
原五圖的舊截圖與 `*.result.json` 只是先前版本紀錄，新瀏覽器驗收由主代理統一執行。
