# 台灣地球物理觀測資料探索（教學網站）

以中央氣象署 **GDMS（臺灣地震與地球物理資料管理系統）** 的真實資料為核心的研究生教學網站：認識台灣的地下水位、地磁、地動（地震波形／目錄）、GNSS 四類觀測，學會程式化下載與視覺化，探索它們與地震的關係。

- 教學內容：`book/`（Jupyter Book，8 章，含 2024 花蓮 M7.2 地震案例分析）
- 下載工具：`gdms_toolkit/`（GDMS 登入、申請、取件、讀檔）
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
的申請清單下載，或直接執行 `book/02_download.ipynb`）。

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
