# 地震統計與預報：從模型到台灣觀測

面向有基礎統計的大學生的連貫教學網站。先從地震目錄建立統計問題，再理解模型與評估，最後將相同的判斷帶到台灣物理觀測。

- 第一部：20 章，以 HORUS 義大利目錄及 177 格預報實驗貫穿統計工具、SUP／PPE／ETAS／EEPAS 與檢驗。
- 第二部：台灣地球物理觀測，9章。從觀測機制與資料品質，逐步讀地下水、地磁、地動與GNSS，接到花蓮案例與台灣展望。
- 第三部：A–G技術附錄。長推導、演算法與資料操作集中查閱，與主文往返連結。

正文以直覺、動機和圖例為主，互動科學展以地圖、時間軸、Hawkes率曲線、分支、分布及比較視圖帶入概念。網站不顯示程式輸入，原始notebook保留可重現性。
第一部採新檔名及 01–20 章號；第二部保留檔名，顯示 21–29 章。舊第一部網址由新版章序取代。

## 原始資料與程式

`book/` 是Jupyter Book來源；以Jupytext配對 `.py`／`.ipynb`，另有Markdown頁。
`gdms_toolkit/` 提供資料讀取、義大利預報快取、檢驗、視覺化與教學合成目錄。`reference/notes/` 保存文獻核對及改編紀錄。
觀測資料快取在 `data/cache/`，不納入版控；PDF原文亦不公開。使用者其他研究目錄不屬於教學網站建置範圍。

## 本機準備

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

只閱讀或建置已有輸出的教學頁不需要GDMS帳密。自行下載觀測資料時，另依附錄G準備帳號、資料與 `.env`，不要將帳密放進notebook。

## 重現義大利實驗

```bash
python scripts/fetch_italy_data.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/fetch_italy_data.py --forecasts
```

依公開 EEPAS README 指向的 Drive 資料夾取得 HORUS；兩個區域檔從固定的公開 Git commit 下載。
不需要另外複製 EEPAS 程式庫。原始檔、轉換檔與模型快取皆留在 `data/cache/`，不進版控。
HORUS 資料使用條款與來源見快取中的原始 README／DataOrigin；不要把資料和本書原始碼一同重新散布。

先準備預報，再執行章節；不要同時開多個預報計算或瀏覽器工作。
例如：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/execute_teaching_pages.py --in-process 01_forecast_question
python scripts/verify_teaching_math.py
```

`--in-process` 使用 IPython 捕捉完整輸出，不啟動 kernel socket；驗證時停用外部資料請求，並限制記憶體及每 cell 時間。

## 同步與建置

```bash
python scripts/sync_teaching_notebooks.py
python scripts/build_teaching_offline.py
```

同步工具將 `.py` 文字與標籤更新到notebook，僅承接程式來源完全相同的既有輸出，並列出需執行的頁面。
此建置命令使用已保存輸出，不執行資料取得。程式變更需先重新執行相關頁；
`scripts/execute_teaching_pages.py` 可明確選頁並停用外部請求。直接使用
`jupyter-book build book/` 則遵循原有快取執行設定，可能執行尚未快取的程式，
資料取得頁不應在無帳號或缺資料時盲目重跑。
本次改編的離線建置、圖表執行與瀏覽器驗收記錄見 `reference/notes/refresh_20260913_exhibition/`。

完整重現觀測分析需依附錄G準備公開資料；來源規範見中央氣象署GDMS，DOI：
[10.7914/SN/T5](https://doi.org/10.7914/SN/T5)。

## 發布

本機建置輸出在 `book/_build/html/`。只有明確授權後才commit、push或透過ghp-import發布；不在一般建置時自動發布。

## 展件與品質檢查

十四個互動展件共用 `book/_static/exhibits.css`／`exhibits.js`，各自有不同的幾何與操作。
重建展件與離線預覽：

```bash
python book/_static/diagrams/sources/build_teaching_diagrams.py
```

新展件來源或共用資產修改後，重新執行使用該展件的章節，讓保存的輸出與來源一致。
執行入口會先準備可寫的 Matplotlib 快取；未處理的 stderr 不可進入公開 notebook。
建置工具會檢查公開內容是否殘留警告、粗體標記或缺圖占位。

```bash
OPENBLAS_NUM_THREADS=1 python scripts/verify_teaching_math.py
python scripts/verify_teaching_presentation.py
```

瀏覽器驗收另需 Playwright 與 Chromium；本機使用既有安裝，在3GB記憶體上限內逐頁執行
`check_exhibition_browser.py`、`check_exhibition_lifecycle.py`及`check_exhibit_geometry.py`。
文獻來源、去重與核對深度見本輪 `literature_coverage.md`；59個PDF登錄不等於59篇全文核讀。

2026-09-14 的全文讀者審查與側欄修正記錄在 `reference/notes/reader_review_20260914/`。
桌面常駐章節目錄，窄螢幕保留收合操作；`scripts/check_teaching_navigation.py` 檢查章節跳轉與響應式切換。
同步建置後可用 `scripts/check_teaching_source_integrity.py` 檢查成對來源、頁面內部連結與錨點。
