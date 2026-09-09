# 全站延伸閱讀驗證紀錄（2026-09-09）

本次範圍為首頁與第 01–23 章，依 `book/_toc.yml`，共 24 頁、87 筆閱讀條目。
每頁 3–5 筆，至少一項免費可讀來源；論文與網站皆附中文導讀。

## 選材與來源

優先回查原始 `reference/` PDF 與筆記，再核對出版社、官方網站或作者機構。
書目及原始檔案對照見 [00–08](reading_00_08_verified.md)、
[09–16](reading_09_16_verified.md)、[17–23](reading_17_23_verified.md)。
第 8 章原有三篇綜述均保留，改為有連結與導讀的完整條目。
不公開本機 PDF，不使用專案規範禁止公開的研究。

## 確認無誤

- 全部 24 頁有單一章末「參考資料與延伸閱讀」，共 87 筆，每頁 3–5 筆。
- 全部 22 組 `.py`／`.ipynb` 的 cell 類型與文字一致。
- 相對本次起點 `820ab327f088f6b61b0b8487413e5bdd5ecb87ee`，22 份 notebook 的全部原始 cells、metadata、程式碼與輸出完全相同，只附加一個 Markdown cell。
- 建置用暫存來源逐一核對程式碼與 `book/_build/jupyter_execute/` 快取一致，帶入 155 筆既有輸出及快取的 `language_info`，停用 notebook 執行。沒有重跑資料下載或分析。
- Jupyter Book 完整建置成功，0 warnings；全部 24 頁的輸出區塊、Plotly 圖及圖片數與修改前 HTML 相同。
- Chromium 實際檢查首頁、第 3、15、23 章，1440×1000 與 390×844 共 8 個情境。新增區塊無水平溢出、外部連結存在且可鍵盤聚焦；已檢視截圖總覽。
- `git diff --check` 通過。

## 連結與限制

- 目前共有 88 個不同的外部網址。直接 HTTP 核對：62 個回傳 2xx、22 個受出版社／典藏網站的自動化防護而回傳 403、4 個在本機信任憑證驗證時失敗；沒有剩餘 404。
- 403 不作為書目無效的判準：已用出版搜尋內容、官方頁面、作者機構及原始 PDF 交叉核實。不能據此宣稱所有網址在每種網路環境都可直接開啟。
- 4 個本機憑證失敗的來源為 GDMS、氣象署 GDMS 專刊、氣象署最近地震與 TEC 花蓮事件報告。web 工具已確認官方內容或索引；未關閉 TLS 驗證。
- Zechar 2010 的 DOI 轉址指向已刪除的舊 PDF，實測 404。兩章標題均改連實測 200 的 CORSSA 官方 PDF，DOI 保留為書目文字。
- 瀏覽器仍可見 Thebe 的 `THEBE_JS_URL` 重複宣告；修改前 HTML 已有相同兩次宣告，本次未修改相關腳本。此紀錄不宣稱既有全站 JavaScript 無錯誤。
- 本次為延伸閱讀整理，不重新稽核正文的數學、觀測結果或現行制度。

## 可重現證據

本機完整證據保存於 `book/_build/reading-verification-20260909/`（建置產物，不進 main 版控）：

- `preservation.json`：逐 notebook 快取對照與保留輸出數。
- `build_readings.py`、`build-command.json`、`build.log`：停用執行的建置方式與完整紀錄。
- `links_initial.json`、`links.json`：初次與修復後的連結結果；`check_reading_links.py` 是初次檢查方法。
- `browser_readings.py`、`browser.json`、`baseline-browser.json`：瀏覽器檢查方法與結果。
- `screenshots/`、`montage.png`：8 個實際畫面與全部截圖總覽。
- `rendered/_build/html/`：本次建置的完整公開網站。

視覺處理採現有介面的延伸：保留 Jupyter Book 字型、配色、導覽及 Markdown 清單，沒有新主題、動畫或介面元件。
