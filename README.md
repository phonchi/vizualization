# 台灣地球物理觀測資料探索與地震預報（教學網站）

以中央氣象署 **GDMS（臺灣地震與地球物理資料管理系統）** 的真實資料為核心的研究生教學網站，分兩部：

- **第一部（1–8 章）：地球物理觀測資料探索**——認識台灣的地下水位、地磁、地動（地震波形／目錄）、GNSS 四類觀測，學會程式化下載與視覺化，探索它們與地震的關係（含 2024 花蓮 M7.2 案例分析）。
- **第二部（9–17 章）：地震預報模型**——觀念導向的地震預報教學：目錄統計進階（Mc、b 值、除叢）、ETAS、EEPAS 與 PPE、STEP 與作業化預報（OEF）、模型組合、CSEP 檢驗、PSHA，以及台灣的預報現況與展望。示意圖以輕量模擬與 1973–2025 台灣長期目錄產生，程式碼預設摺疊。

- 教學內容：`book/`（Jupyter Book，共 17 章）
- 下載工具：`gdms_toolkit/`（GDMS 登入、申請、取件、讀檔；`load_taiwan_catalog()` 載入 1973–2025 長期目錄）
- 文獻筆記：`reference/notes/`（第二部所依據約 40 篇論文的主題式中文筆記；PDF 原檔不進版控）
- 測站清單：`data/stations/`（7 個測網共 980 站，取自 GDMS）
- 資料快取：`data/cache/`（下載的原始資料，不進版控）

## 快速開始

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                # 安裝 gdms_toolkit
cp .env.example .env            # 填入 GDMS 帳號（免費註冊：https://gdms.cwa.gov.tw/）
```

### 建置網站

```bash
jupyter-book build book/
# 開啟 book/_build/html/index.html
```

建置時會實際執行所有 notebook；第一次需要 `data/cache/` 內的資料（依第 2 章
的申請清單下載，或直接執行 `book/02_download.ipynb`）。第二部（10、12、17 章）
另需 1973–2025 長期目錄：`reference/Taiwan/EEPAS_TW-main.zip`（不進版控）存在時
會自動解壓合併並快取為 `data/cache/catalog_1973_2025.csv`；沒有 zip 的環境請直接
提供該快取 CSV。

### 用 toolkit 下載資料（三行版）

```python
import gdms_toolkit as gt
g = gt.GDMSSession()   # 讀 .env 登入
gt.request_geophysical(g, "GW", "TUN", "2024-04-01", "2024-04-07", label="myreq")
gt.wait_and_fetch(g, "myreq")   # 等打包完成並下載到 data/cache/
```

## 章節撰寫格式

各章以 jupytext py:percent 格式撰寫（`book/0X_*.py`），修改後轉回 notebook：

```bash
jupytext --to ipynb book/0X_*.py
```

## 部署到 GitHub Pages（選用）

```bash
pip install ghp-import
ghp-import -n -p -f book/_build/html
```

## 資料引用

觀測資料版權屬中央氣象署，透過 GDMS（<https://gdms.cwa.gov.tw/>）取得，
DOI: [10.7914/SN/T5](https://doi.org/10.7914/SN/T5)。發表時請依 GDMS 規範引用。
