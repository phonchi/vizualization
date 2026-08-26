"""GDMS 地震目錄（需會員登入）。

對應網頁「地震目錄下載」（catalogDownload.php），資料為中央氣象署
正式地震目錄，含定位品質資訊（測站數、間隙角、誤差、品質等級）。
時間為 UTC；台灣時間 = UTC + 8。回傳的 time 欄為無時區的 UTC 時間戳。
"""

import zipfile
from pathlib import Path

import pandas as pd

from .auth import BASE, GDMSSession
from .download import CACHE_DIR

# 台灣周邊預設範圍
TAIWAN_BOX = dict(stlat="21.0", edlat="26.0", stlon="119.0", edlon="123.5")

# 1973–2025 長期目錄的來源 zip（CWA 公開目錄資料，隨 repo 提供）
_LONG_CATALOG_ZIP = (Path(__file__).resolve().parent.parent
                     / "reference" / "Taiwan" / "EEPAS_TW-main.zip")
_LONG_CATALOG_MEMBERS = [
    f"EEPAS_TW-main/data/GDMScatalog_{span}.csv"
    for span in ("1973_1990", "1991_2001", "2002_2010", "2011_2019", "2020_2025")
]


def load_taiwan_catalog(min_ml: float = 0.0) -> pd.DataFrame:
    """載入 1973–2025 台灣長期地震目錄（約 35 萬筆，M ≥ 2，UTC）。

    第一次呼叫時從 reference/Taiwan/EEPAS_TW-main.zip 解出五段 CSV
    合併後快取到 data/cache/catalog_1973_2025.csv，之後直接讀快取。
    欄位名稱對齊 gdms_earthquake_catalog（latitude / longitude / time）。
    """
    cache_csv = CACHE_DIR / "catalog_1973_2025.csv"
    if cache_csv.exists():
        cat = pd.read_csv(cache_csv, parse_dates=["time"])
    else:
        with zipfile.ZipFile(_LONG_CATALOG_ZIP) as zf:
            parts = [pd.read_csv(zf.open(m)) for m in _LONG_CATALOG_MEMBERS]
        cat = pd.concat(parts, ignore_index=True)
        cat["time"] = pd.to_datetime(cat["date"] + " " + cat["time"],
                                     format="mixed", errors="coerce")
        cat = (cat.drop(columns=["date"])
                  .rename(columns={"lat": "latitude", "lon": "longitude"})
                  .sort_values("time").reset_index(drop=True))
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cat.to_csv(cache_csv, index=False)
    if min_ml > 0:
        cat = cat.loc[cat["ML"] >= min_ml].reset_index(drop=True)
    return cat


def gdms_earthquake_catalog(gdms: GDMSSession, start: str, end: str,
                            min_ml: float = 3.0, max_ml: float = 10.0,
                            min_depth: float = 0, max_depth: float = 350,
                            box: dict | None = None) -> pd.DataFrame:
    """下載 GDMS 地震目錄，回傳 DataFrame（時間欄已合併為 datetime）。

    start / end: 'YYYY-MM-DD'
    box: dict(stlat, edlat, stlon, edlon)，預設為台灣周邊

    >>> cat = gdms_earthquake_catalog(gdms, "2024-03-01", "2024-05-01", min_ml=4)
    """
    r = gdms.post(f"{BASE}/php/dbconnect/getCatalog.php", data={
        "stdate": start, "eddate": end,
        "stML": str(min_ml), "edML": str(max_ml),
        "sttime": "00:00:00", "edtime": "23:59:59",
        "stdep": str(min_depth), "eddep": str(max_depth),
        **(box or TAIWAN_BOX),
    }, timeout=180)
    r.raise_for_status()
    data = r.json()
    if isinstance(data, dict) and "error" in data:
        raise RuntimeError(f"目錄查詢失敗：{data['error']}")
    df = pd.DataFrame(data)
    for col in ("latitude", "longitude", "depth", "ML"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["time"] = pd.to_datetime(df["date"] + " " + df["time"], errors="coerce")
    return df.drop(columns=["date", "ms"]).sort_values("time").reset_index(drop=True)
