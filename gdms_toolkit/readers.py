"""資料讀取輔助（依 GDMS 實際下載格式撰寫）。

- read_waveform():    地震波形 miniSEED / SAC → obspy.Stream
- read_groundwater(): 地下水 tgz 內的每日 CSV（YYYYMMDD.站碼）→ pandas
- read_geomagnetic(): 地磁 IAGA-2002 秒資料（站碼YYYYMMDDdsec.sec）→ pandas

時間皆為 UTC；台灣時間 = UTC + 8（to_taipei=True 自動轉換）。
缺測值：地下水 9999、地磁 88888，讀入時轉為 NaN。
"""

import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd


def _to_taipei(df: pd.DataFrame, to_taipei: bool) -> pd.DataFrame:
    df.index = df.index.tz_localize("UTC")
    if to_taipei:
        df.index = df.index.tz_convert("Asia/Taipei")
    return df


def read_waveform(path):
    """讀取 miniSEED / SAC 波形，回傳 obspy.Stream。

    >>> st = read_waveform("data/cache/edu-wave-hualien2024.mseed")
    >>> st.plot()
    """
    from obspy import read
    return read(str(path))


def read_groundwater(tgz_path, station: str, start: str | None = None,
                     end: str | None = None, resample: str | None = "1min",
                     to_taipei: bool = False) -> pd.DataFrame:
    """從 GDMS 地下水 tgz 讀取單一測站的時間序列。

    tgz 內為每日一檔的 CSV（檔名 YYYYMMDD.站碼，1 秒取樣、86400 列），
    欄位：Time, WaterLevel(cm), Atmospheric pressure(hPa),
          Bottom Temperature, Upper Temperature

    station: 站碼（如 'TUN'）
    start/end: 'YYYY-MM-DD'，限制讀取的日期範圍（預設全部）
    resample: 重取樣間隔（預設 '1min' 平均，原始 1 秒資料量大）；None 表示不重取樣
    """
    tgz_path = Path(tgz_path)
    frames = []
    with tarfile.open(tgz_path, "r:gz") as tar:
        for m in sorted(tar.getmembers(), key=lambda m: m.name):
            name = Path(m.name).name  # YYYYMMDD.STA
            if not name.endswith(f".{station.upper()}"):
                continue
            day = name.split(".")[0]
            if start and day < start.replace("-", ""):
                continue
            if end and day > end.replace("-", ""):
                continue
            df = pd.read_csv(tar.extractfile(m), index_col="Time",
                             parse_dates=True)
            frames.append(df)
    if not frames:
        raise FileNotFoundError(f"{tgz_path.name} 內找不到 {station} 的資料")
    df = pd.concat(frames).sort_index()
    df.columns = ["water_level_cm", "pressure_hPa", "temp_bottom_C", "temp_upper_C"]
    df = df.replace({9999.00: np.nan, 88.0: np.nan})
    if resample:
        df = df.resample(resample).mean()
    return _to_taipei(df, to_taipei)


def read_geomagnetic(tgz_path, station: str, start: str | None = None,
                     end: str | None = None, resample: str | None = "1min",
                     to_taipei: bool = False) -> pd.DataFrame:
    """從 GDMS 地磁 tgz 讀取單一測站的 IAGA-2002 秒資料。

    欄位：X, Y, Z, F（nT）。另計算水平分量 H = sqrt(X^2 + Y^2)。
    station: 站碼（如 'XCG'）；resample 預設 '1min' 平均。
    """
    tgz_path = Path(tgz_path)
    frames = []
    with tarfile.open(tgz_path, "r:gz") as tar:
        for m in sorted(tar.getmembers(), key=lambda m: m.name):
            name = Path(m.name).name  # staYYYYMMDDdsec.sec
            if not (name.startswith(station.lower()) and name.endswith(".sec")):
                continue
            day = name[len(station):len(station) + 8]
            if start and day < start.replace("-", ""):
                continue
            if end and day > end.replace("-", ""):
                continue
            frames.append(_parse_iaga2002(tar.extractfile(m)))
    if not frames:
        raise FileNotFoundError(f"{tgz_path.name} 內找不到 {station} 的資料")
    df = pd.concat(frames).sort_index()
    if resample:
        df = df.resample(resample).mean()
    df["H"] = np.hypot(df["X"], df["Y"])          # 水平分量
    df["F_calc"] = np.hypot(df["H"], df["Z"])     # 由三分量計算的全磁力
    return _to_taipei(df, to_taipei)


def _parse_iaga2002(fh) -> pd.DataFrame:
    """解析 IAGA-2002 文字格式（標頭行以 | 結尾）。"""
    text = fh.read()
    if isinstance(text, bytes):
        text = text.decode("ascii", errors="replace")
    lines = text.splitlines()
    data_start = next(i for i, l in enumerate(lines) if l.startswith("DATE")) + 1
    df = pd.read_csv(
        io.StringIO("\n".join(lines[data_start:])), sep=r"\s+", header=None,
        names=["date", "time", "doy", "X", "Y", "Z", "F"],
    )
    df.index = pd.to_datetime(df["date"] + " " + df["time"])
    df = df[["X", "Y", "Z", "F"]].replace({88888.00: np.nan, 99999.00: np.nan})
    return df


def read_gdms_timeseries(path, **kw):
    """（保留的泛用介面）依副檔名自動選擇讀取器。"""
    p = Path(path)
    if p.suffix in (".mseed", ".sac"):
        return read_waveform(p)
    raise ValueError("請改用 read_groundwater() / read_geomagnetic() 並指定站碼")
