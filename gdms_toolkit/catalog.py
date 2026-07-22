"""GDMS 地震目錄（需會員登入）。

對應網頁「地震目錄下載」（catalogDownload.php），資料為中央氣象署
正式地震目錄，含定位品質資訊（測站數、間隙角、誤差、品質等級）。
時間為 UTC；台灣時間 = UTC + 8。回傳的 time 欄為無時區的 UTC 時間戳。
"""

import pandas as pd

from .auth import BASE, GDMSSession

# 台灣周邊預設範圍
TAIWAN_BOX = dict(stlat="21.0", edlat="26.0", stlon="119.0", edlon="123.5")


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
