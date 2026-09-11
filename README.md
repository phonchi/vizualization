# 地震統計與預報：從模型到台灣觀測

面向有基礎統計的大學生的連貫教學網站。先從地震目錄建立統計問題，再理解模型與評估，最後將相同的判斷帶到台灣物理觀測。

- 第一部：地震統計與預報，17章。新增目錄、隨機性與統計推論三章，串起點過程、ETAS、EEPAS、檢驗、複發與危害。
- 第二部：台灣地球物理觀測，9章。從觀測機制與資料品質，逐步讀地下水、地磁、地動與GNSS，接到花蓮案例與台灣展望。
- 第三部：A–G技術附錄。長推導、演算法與資料操作集中查閱，與主文往返連結。

正文以直覺、動機和圖例為主，短補充預設收合。網站不顯示程式輸入，原始notebook保留可重現性。
既有檔名保留舊網址；顯示章號依 `book/_toc.yml`，不依檔名前綴。

## 原始資料與程式

`book/` 是Jupyter Book來源；以Jupytext配對 `.py`／`.ipynb`，另有Markdown頁。
`gdms_toolkit/` 提供資料讀取、視覺化與教學合成目錄。`reference/notes/` 保存文獻核對及改編紀錄。
觀測資料快取在 `data/cache/`，不納入版控；PDF原文亦不公開。使用者其他研究目錄不屬於教學網站建置範圍。

## 本機準備

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

只閱讀或建置已有輸出的教學頁不需要GDMS帳密。自行下載觀測資料時，另依附錄G準備帳號、資料與 `.env`，不要將帳密放進notebook。

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
本次改編的離線建置、圖表執行與瀏覽器驗收記錄見 `reference/notes/rewrite_20260911/`。

完整重現觀測分析需依附錄G準備公開資料；來源規範見中央氣象署GDMS，DOI：
[10.7914/SN/T5](https://doi.org/10.7914/SN/T5)。

## 發布

本機建置輸出在 `book/_build/html/`。只有明確授權後才commit、push或透過ghp-import發布；不在一般建置時自動發布。
