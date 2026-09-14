# 此文件為先前圖版的歷史說明

目前的十四張互動科學展件及重建方式，請讀 [EXHIBITS.md](EXHIBITS.md)。以下文字與舊 radio 圖、舊收據相對應。

# 2026-09-13 十二張教學圖與手機重排

`build_teaching_diagrams.py` 重建十二個嵌入片段及同名的 `../standalone/` 離線預覽。
樣式沿用 `.diagram-design` 指定的 `quake-teaching`，顏色使用站台 `--dg-*` tokens。
主文嵌入格式與 `gdms_toolkit.viz.show_diagram()` 相容。

## 來源與語意

- 教學規格：`reference/notes/rewrite_20260913/CHAPTER_BRIEF.md` 與本目錄上層 `README.md`。
- `d01_forecast_pipeline`、`d16_csep_tests_flow`：archify workflow v2 JSON 為流程來源；`../standalone/*.archify.html` 由正式 `deliver` 產生，收據在 `../standalone/checks/*.delivery.json`。兩份皆通過 showcase 全 9 項檢查、0 composition errors／warnings。
- 嵌入版重排成單欄小圖，以保留 390 px 內容區的可讀性。d01 保留五步；d16 合併「選問題」到第一步，合併「算統計量」與「比較」到最後一步，另用原生選項逐一列出 N／S／M／cL／L 的條件。未更動檢驗語意。
- archify Viewer 的固定介面及 `<html lang>` 依套件回退為 English；作者文字為繁體中文。教學片段與同名離線預覽為繁體中文介面。
- `d09_kernel_smoothing` 使用兩個合成 Gaussian 空間剖面：中心為 144、304，標準差 52；實線為逐點和。只展示加法，不宣稱是義大利 PPE 數值或核形式。
- `d12_etas_branching` 是一棵合成的兩代樹；不能從這棵樹的後代數估計分支比。
- `d15_eepas_kernels` 的三核形狀為概念示意；時間形狀為正時間的對數常態示意，空間以等值圈表示。未代入任何審稿中論文數字。
- `d19_convex_mix` 用兩格合成預報 A=(2,6)、B=(6,2)，權重各 1/2，混合=(4,4)。柱高由同一尺度直接計算，總數皆 8。
- `d20_model_family` 將 SUP 的固定率與 PPE／ETAS／EEPAS 的歷史更新分開；PPE 與 EEPAS 均設 50 天延遲。

## 互動及驗證

教學版以原生 radio、label 與 CSS `:has()` 切換說明。滑鼠 hover、鍵盤 Tab／方向鍵／Space、觸控都可使用，沒有自訂 JS 或外部字型下載。SVG 始終保留完整語意；列印顯示所有說明。字型回退到系統已安裝字型。

- `../standalone/checks/new_diagrams_self_check.json`：7 片段＋7 預覽通過 diagram-design 的正式 accessible SVG／單檔安全 self-check，含檔案 SHA-256。
- `../standalone/checks/new_diagrams_geometry.json`：7 預覽的 geometry checker 均無 finding。
- `../standalone/checks/new_diagrams_static.json`：圖檔、ID、ARIA 參照、選項與外部依賴檢查。
- **瀏覽器與實際視覺驗收尚待主代理統一執行**；本子工作沒有啟動 Chromium。靜態通過不等於手機文字、hover 和實際 Jupyter Book 嵌入已驗收。

重建：

```bash
python -B book/_static/diagrams/sources/build_teaching_diagrams.py
```

archify 的已驗證 JSON 與 `.archify.html` 不由這個 generator 重寫。

## 原五圖的手機重排

主代理的 390 px 截圖顯示原五圖採 720 viewBox、12 px 字級，實際縮到約 6 px，無法閱讀。
已將 d02、d03 兩圖、d04、d06 改為 480 viewBox、20 px SVG 標籤，長說明使用 16 px HTML 與原生 radio 選項。
所有 SVG 語意保持可見，不靠縮字、裁切或橫向捲動容納內容。選項只切換下方長說明。

- d02 使用矩形集合圖，明說不是實際義大利 polygon，也未假裝畫出全部 177 格。
- d03 timeline：三段期間是用途卡，卡片寬度不表示年數；下方單窗時間軸按天數比例畫出發報前 50 天及其後 91.31 天。PPE／EEPAS 都延遲 50 天，SUP 固定率、ETAS 不採此延遲。說明保留測試期最後 0.6 天不評分。
- d03 thresholds：分開 Mc、m₀、mT 的角色。名目 2.5／5.0 分別對應原始箱界 2.45／4.95；只在相同表示法下比較差距。25 個規模箱按實際個數示意。
- d04：保留密度→積分→期望數→計數分布。明說只有固定密度時，積分才是單個密度值乘範圍；Poisson 是另加假設。
- d06：使用頻數 `[1,2,4,6,4,2,1]` 的 20 次合成模擬，觀測統計量 2，因此 q=7/20、同值=4/20、上尾=17/20。q 只有在指定下尾問題時，才是該尾 p 值的模擬近似。

新收據：`legacy_five_mobile_self_check.json`（10／10）、`legacy_five_mobile_geometry.json`（5／5 零 finding）、`all_twelve_mobile_static.json`（12 圖 SHA-256、112 個唯一 ID）。都位於 `../standalone/checks/`。
原五份 `*.result.json` 及截圖保留為歷史結果，已加 superseded 註記，不能用來驗收重排後內容。
瀏覽器驗收由主代理另行執行，本子工作未啟動 Chromium。

只重建指定五圖：

```bash
python -B book/_static/diagrams/sources/build_teaching_diagrams.py d02_regions_s_r d03_timeline_windows d03_thresholds d04_density_to_counts d06_quantile_pvalue
```

## 觸控與可讀性修正

Hover 預覽只在 `(hover:hover) and (pointer:fine)` 裝置啟用，避免觸控保留滑鼠位置時蓋住已選說明。所有 SVG 標籤保持完整不透明度；選項與說明區呈現目前焦點。
