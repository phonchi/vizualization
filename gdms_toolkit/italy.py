"""義大利標準預報實驗的資料層（教學用）。

實驗設定沿用 Biondini, Rhoades & Gasperini (2023, GJI 234:1681)：
- 目錄：HORUS（Lolli et al. 2020），Mw 均一化，1960 年起。
- 收集區 S：CPTI15 polygon；測試區 R：177 個邊長 30√2 km 的方格。
- 座標：RDN2008 / Italy zone (E-N)，EPSG:7794，以公里為單位。
- 深度 ≤ 40 km；輸入門檻 m0 = 2.45（名目 2.5）；目標門檻 mT = 5.0（有效 4.95）。
- warm-up 1960–1989、learning 1990–2011、testing 2012–2021；三個月滾動窗。

原始檔取自 PyEEPAS 公開 repo 的 data/README 所指 Google Drive 資料夾，
以 scripts/fetch_italy_data.py 下載到 data/cache/italy/（不進版控）。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.io import loadmat

from .download import CACHE_DIR

ITALY_DIR = CACHE_DIR / "italy"
RAW_CATALOG = ITALY_DIR / "HORUS_Ita_Catalog.txt"
CLEAN_CATALOG = ITALY_DIR / "horus_clean.csv"
FORECAST_DIR = CACHE_DIR / "italy_forecasts"
_REGION_DIR = ITALY_DIR / "regions"
_LOCAL_REGION_DIR = Path(__file__).resolve().parent.parent / "EEPAS" / "data"
CELLS_MAT = (_REGION_DIR if (_REGION_DIR / "CELLE_ter.mat").exists() else _LOCAL_REGION_DIR) / "CELLE_ter.mat"
POLYGON_MAT = (_REGION_DIR if (_REGION_DIR / "CPTI15.mat").exists() else _LOCAL_REGION_DIR) / "CPTI15.mat"
DRIVE_FOLDER = "https://drive.google.com/drive/folders/170WNb8M8PQJDX1B2JSQYfqj4ad80TC0U"

CRS_KM = "EPSG:7794"          # RDN2008 / Italy zone (E-N)，公尺；本模組一律換算成公里
T0 = pd.Timestamp("1960-01-01")   # 目錄起點，decimal days 的零點
DAYS_PER_YEAR = 365.2425


@dataclass(frozen=True)
class ExperimentSpec:
    """規格卡：一份預報實驗要先說清楚的所有欄位。"""
    catalog: str = "HORUS（Lolli et al. 2020），Mw 均一化，1960 年起"
    depth_max_km: float = 40.0
    m0: float = 2.45          # 輸入（前兆）門檻，名目 2.5
    mT: float = 5.0           # 目標門檻，名目 5.0（有效 4.95）
    m_max: float = 7.5
    mag_bin: float = 0.1
    warmup: tuple = (1960, 1989)
    learning: tuple = (1990, 2011)
    testing: tuple = (2012, 2021)
    window_days: float = 91.31          # 三個月滾動窗
    delay_days: float = 50.0            # EEPAS：事件後 50 天內不當前兆
    b_value: float = 1.084              # Biondini 2023：讓 SUP 在學習期剛好預期 27 顆
    n_cells: int = 177
    cell_side_km: float = 30 * np.sqrt(2)
    crs: str = CRS_KM


SPEC = ExperimentSpec()
BETA = SPEC.b_value * np.log(10)


# ---------------------------------------------------------------- 座標 ----
def _transformer(inverse: bool = False):
    from pyproj import Transformer
    if inverse:
        return Transformer.from_crs(CRS_KM, "EPSG:4326", always_xy=True)
    return Transformer.from_crs("EPSG:4326", CRS_KM, always_xy=True)


def lonlat_to_km(lon, lat):
    """WGS84 經緯度 → RDN2008 Italy zone 公里座標 (x_km 東, y_km 北)。"""
    x, y = _transformer().transform(np.asarray(lon), np.asarray(lat))
    return x / 1000.0, y / 1000.0


def km_to_lonlat(x_km, y_km):
    lon, lat = _transformer(inverse=True).transform(np.asarray(x_km) * 1000.0,
                                                    np.asarray(y_km) * 1000.0)
    return lon, lat


def decimal_days(t: pd.Series | pd.Timestamp) -> np.ndarray | float:
    """自 1960-01-01 起算的 decimal days（ETAS_R 與 EEPAS 都用這種時間軸）。"""
    delta = pd.to_datetime(t) - T0
    if isinstance(delta, pd.Timedelta):
        return delta.total_seconds() / 86400.0
    return (delta.dt.total_seconds() / 86400.0).to_numpy()


def year_start_days(year: int) -> float:
    return decimal_days(pd.Timestamp(f"{year}-01-01"))


# ---------------------------------------------------------------- 區域 ----
def collection_polygon() -> pd.DataFrame:
    """收集區 S（CPTI15 polygon）的頂點，含經緯度與公里座標；最後一點重複首點。"""
    v = loadmat(POLYGON_MAT)["cpti15"]
    lon, lat = v[:, 0], v[:, 1]
    x_km, y_km = lonlat_to_km(lon, lat)
    return pd.DataFrame(dict(lon=lon, lat=lat, x_km=x_km, y_km=y_km))


def testing_cells() -> pd.DataFrame:
    """測試區 R 的 177 個方格：公里邊界、格心經緯度、四角經緯度（供畫圖）。"""
    m = loadmat(CELLS_MAT)["CELLESD"]
    df = pd.DataFrame(dict(cell=np.arange(len(m)),
                           x0=m[:, 0], x1=m[:, 1], y0=m[:, 2], y1=m[:, 3]))
    df["xc"] = (df.x0 + df.x1) / 2
    df["yc"] = (df.y0 + df.y1) / 2
    df["area_km2"] = (df.x1 - df.x0) * (df.y1 - df.y0)
    df["lon_c"], df["lat_c"] = km_to_lonlat(df.xc, df.yc)
    corners = []
    for _, r in df.iterrows():
        xs = [r.x0, r.x1, r.x1, r.x0, r.x0]
        ys = [r.y0, r.y0, r.y1, r.y1, r.y0]
        lo, la = km_to_lonlat(xs, ys)
        corners.append((lo, la))
    df["corner_lon"] = [c[0] for c in corners]
    df["corner_lat"] = [c[1] for c in corners]
    return df


def cell_of(x_km, y_km, cells: pd.DataFrame | None = None) -> np.ndarray:
    """每個點落在哪一格（0..176），不在 R 內回傳 -1。"""
    cells = testing_cells() if cells is None else cells
    x = np.asarray(x_km, float)
    y = np.asarray(y_km, float)
    out = np.full(x.shape, -1, dtype=int)
    for _, c in cells.iterrows():
        m = (x >= c.x0) & (x < c.x1) & (y >= c.y0) & (y < c.y1)
        out[m] = int(c.cell)
    return out


def in_collection_region(lon, lat) -> np.ndarray:
    """自己算一次點是否在 CPTI15 polygon 內（與 HORUS 的 Geo-CPTI15 旗標對照）。"""
    poly = collection_polygon()
    path = MplPath(np.c_[poly.lon, poly.lat])
    return path.contains_points(np.c_[np.asarray(lon), np.asarray(lat)])


# ---------------------------------------------------------------- 目錄 ----
_RAW_COLS = ["year", "month", "day", "hour", "minute", "second",
             "lat", "lon", "depth", "mw", "sig_mw",
             "geo_ita", "geo_cpti15", "ev_type", "iside"]


def load_horus(refresh: bool = False) -> pd.DataFrame:
    """讀完整 HORUS 目錄（1960 起，約 49 萬筆），不做任何實驗篩選。

    欄位：time（UTC）、lat、lon、depth（km）、mw、sig_mw、
    in_italy（Geo-Ita 旗標）、in_cpti15（Geo-CPTI15 旗標）、is_earthquake、
    mb（規模取整到 0.1，Biondini 2023 式 2）。
    """
    if CLEAN_CATALOG.exists() and not refresh:
        return pd.read_csv(CLEAN_CATALOG, parse_dates=["time"])
    if not RAW_CATALOG.exists():
        raise FileNotFoundError(
            f"找不到 {RAW_CATALOG}；請先執行 python scripts/fetch_italy_data.py")
    raw = pd.read_csv(RAW_CATALOG, sep="\t", dtype=str, keep_default_na=False)
    raw = raw.iloc[:, :len(_RAW_COLS)]
    raw.columns = _RAW_COLS
    for c in _RAW_COLS[:11]:
        raw[c] = pd.to_numeric(raw[c].str.strip(), errors="coerce")
    for c in ("geo_ita", "geo_cpti15", "ev_type", "iside"):
        raw[c] = raw[c].str.strip()
    sec = raw["second"].clip(lower=0, upper=59.999)
    time = pd.to_datetime(dict(year=raw.year.astype(int), month=raw.month.astype(int),
                               day=raw.day.astype(int)), errors="coerce")
    time = time + pd.to_timedelta(raw.hour.clip(0, 23), unit="h") \
        + pd.to_timedelta(raw.minute.clip(0, 59), unit="m") \
        + pd.to_timedelta(sec, unit="s")
    cat = pd.DataFrame(dict(
        time=time, lat=raw.lat, lon=raw.lon, depth=raw.depth,
        mw=raw.mw, sig_mw=raw.sig_mw,
        in_italy=raw.geo_ita == "*", in_cpti15=raw.geo_cpti15 == "*",
        is_earthquake=raw.ev_type != "x",
        iside=raw.iside.replace("", np.nan)))
    cat["mb"] = np.floor(cat.mw * 10 + 0.5) / 10
    cat = cat.sort_values("time").reset_index(drop=True)
    ITALY_DIR.mkdir(parents=True, exist_ok=True)
    cat.to_csv(CLEAN_CATALOG, index=False)
    return cat


def experiment_catalog(spec: ExperimentSpec = SPEC, end_year: int | None = None) -> pd.DataFrame:
    """套用規格卡的資料側篩選後的工作目錄。

    保留：地震事件（排除爆炸等）、在收集區 S 內、深度 ≤ 40 km、
    時間 ≤ testing 期末。新增 x_km/y_km、t_days（自 1960-01-01）、
    cell（0..176，R 外為 -1）、in_R、period（warmup/learning/testing）。
    """
    end_year = spec.testing[1] if end_year is None else end_year
    cat = load_horus()
    keep = (cat.is_earthquake & cat.in_cpti15 & (cat.depth <= spec.depth_max_km)
            & (cat.time < pd.Timestamp(f"{end_year + 1}-01-01")))
    cat = cat.loc[keep].copy()
    cat["x_km"], cat["y_km"] = lonlat_to_km(cat.lon, cat.lat)
    cat["t_days"] = decimal_days(cat.time)
    cat["year"] = cat.time.dt.year
    cat["cell"] = cell_of(cat.x_km, cat.y_km)
    cat["in_R"] = cat.cell >= 0
    period = np.where(cat.year < spec.learning[0], "warmup",
                      np.where(cat.year <= spec.learning[1], "learning", "testing"))
    cat["period"] = period
    return cat.reset_index(drop=True)


def target_events(cat: pd.DataFrame, period: str = "testing",
                  spec: ExperimentSpec = SPEC) -> pd.DataFrame:
    """目標地震：在 R 內、取整規模 ≥ mT、落在指定期間。"""
    sel = cat.in_R & (cat.mb >= spec.mT) & (cat.mb < spec.m_max) & (cat.period == period)
    return cat.loc[sel].reset_index(drop=True)


# ---------------------------------------------------------------- 分箱 ----
def magnitude_edges(spec: ExperimentSpec = SPEC) -> np.ndarray:
    """規模箱邊界 5.0, 5.1, …, 7.5（25 箱）。取整規模 mb 落在 [m1, m2)。"""
    n = int(round((spec.m_max - spec.mT) / spec.mag_bin))
    return np.round(spec.mT + spec.mag_bin * np.arange(n + 1), 2)


def forecast_windows(spec: ExperimentSpec = SPEC) -> pd.DataFrame:
    """測試期的滾動預報窗（decimal days 與日期），共 40 窗。"""
    t_start = year_start_days(spec.testing[0])
    total = year_start_days(spec.testing[1] + 1) - t_start
    n = int(total // spec.window_days)
    t1 = t_start + spec.window_days * np.arange(n)
    t2 = t1 + spec.window_days
    return pd.DataFrame(dict(window=np.arange(1, n + 1), t1=t1, t2=t2,
                             start=T0 + pd.to_timedelta(t1, unit="D"),
                             end=T0 + pd.to_timedelta(t2, unit="D")))


def bin_targets(targets: pd.DataFrame, windows: pd.DataFrame | None = None,
                spec: ExperimentSpec = SPEC) -> np.ndarray:
    """把目標地震數到 (窗, 格, 規模箱) 的觀測計數陣列 ω。"""
    windows = forecast_windows(spec) if windows is None else windows
    edges = magnitude_edges(spec)
    omega = np.zeros((len(windows), spec.n_cells, len(edges) - 1), dtype=int)
    k = np.searchsorted(edges, targets.mb.to_numpy(), side="right") - 1
    w = np.searchsorted(windows.t2.to_numpy(), targets.t_days.to_numpy(), side="right")
    ok = ((w >= 0) & (w < len(windows)) & (k >= 0) & (k < len(edges)-1)
          & (targets.cell.to_numpy() >= 0) & (targets.cell.to_numpy() < spec.n_cells))
    valid_window = np.minimum(w, len(windows)-1)
    ok &= targets.t_days.to_numpy() >= windows.t1.to_numpy()[valid_window]
    for wi, ci, ki in zip(w[ok], targets.cell.to_numpy()[ok], k[ok]):
        omega[wi, ci, ki] += 1
    return omega


# ---------------------------------------------------------------- 輸出 ----
def write_csep_10col(rate_by_cell_bin: np.ndarray, path: Path,
                     spec: ExperimentSpec = SPEC, depth=(0.0, 40.0)) -> Path:
    """以 177 格的經緯度外接框輸出 CSEP 十欄格式（教學示範，非重新分格）。

    欄位：LON_0 LON_1 LAT_0 LAT_1 Z_0 Z_1 MAG_0 MAG_1 RATE FLAG。
    方格在投影座標下是正方形，在經緯度下略呈梯形；正式 0.1° 重新分格見附錄 E。
    """
    cells = testing_cells()
    edges = magnitude_edges(spec)
    rows = []
    for _, c in cells.iterrows():
        lon0, lon1 = min(c.corner_lon), max(c.corner_lon)
        lat0, lat1 = min(c.corner_lat), max(c.corner_lat)
        for k in range(len(edges) - 1):
            rows.append((lon0, lon1, lat0, lat1, depth[0], depth[1],
                         edges[k], edges[k + 1], rate_by_cell_bin[int(c.cell), k], 1))
    out = pd.DataFrame(rows, columns=["LON_0", "LON_1", "LAT_0", "LAT_1", "Z_0", "Z_1",
                                      "MAG_0", "MAG_1", "RATE", "FLAG"])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False, float_format="%.6e")
    return path


def spec_card(rows: list[tuple[str, str]], title: str = "預報實驗規格卡",
              pending: str = "待填") -> str:
    """把 (欄位, 內容) 清單轉成規格卡 HTML；內容為 None 時顯示「待填」。"""
    cells_html = "".join(
        f'<div class="spec-row{" pending" if v is None else ""}">'
        f'<div class="spec-key">{k}</div><div class="spec-val">{pending if v is None else v}</div></div>'
        for k, v in rows)
    return f'<div class="spec-card"><div class="spec-title">{title}</div>{cells_html}</div>'


# ---------------------------------------------------------------- 重新分格 ----
_REGRID_CACHE = CACHE_DIR / "italy" / "regrid_01deg.npz"


def regrid_matrix_01deg(step: float = 0.1, n_sub: int = 5, refresh: bool = False):
    """177 個投影方格 → 0.1° 經緯度格的面積分配矩陣 W（n_01 × 177）。

    W[g, c] ≈ 面積(g ∩ c) / 面積(c)，用每個 0.1° 格內 n_sub×n_sub 個取樣點估計；
    對每一來源格的面積權重再正規化，rate_01 = W @ rate_177 保留總期望數。
    交集形狀仍是取樣近似，總量守恆不代表邊界精確。回傳 (W, grid)，
    grid 是每個 0.1° 格的 lon0/lon1/lat0/lat1。
    """
    if step <= 0 or n_sub < 1:
        raise ValueError("step 與 n_sub 必須為正")
    if _REGRID_CACHE.exists() and not refresh and step == 0.1 and n_sub == 5:
        z = np.load(_REGRID_CACHE)
        W = z["W"]
        if np.any(W.sum(axis=0) <= 0):
            raise ValueError("重新分格遺失來源格；請增加取樣密度")
        return W / W.sum(axis=0), pd.DataFrame(dict(lon0=z["lon0"], lon1=z["lon1"], lat0=z["lat0"], lat1=z["lat1"]))
    cells = testing_cells()
    lon_all = np.concatenate(cells.corner_lon.to_list())
    lat_all = np.concatenate(cells.corner_lat.to_list())
    lon_edges = np.arange(np.floor(lon_all.min() * 10) / 10, np.ceil(lon_all.max() * 10) / 10 + step / 2, step)
    lat_edges = np.arange(np.floor(lat_all.min() * 10) / 10, np.ceil(lat_all.max() * 10) / 10 + step / 2, step)
    g_lon0, g_lat0 = np.meshgrid(lon_edges[:-1], lat_edges[:-1], indexing="ij")
    g_lon0, g_lat0 = g_lon0.ravel(), g_lat0.ravel()
    u = (np.arange(n_sub) + 0.5) / n_sub
    W = np.zeros((g_lon0.size, len(cells)))
    for i, (lo0, la0) in enumerate(zip(g_lon0, g_lat0)):
        lo = lo0 + u * step
        la = la0 + u * step
        LO, LA = np.meshgrid(lo, la, indexing="ij")
        x, y = lonlat_to_km(LO.ravel(), LA.ravel())
        cid = cell_of(x, y, cells)
        # 該 0.1° 格的近似面積（km²）
        area_g = (111.32 * np.cos(np.deg2rad(la0 + step / 2)) * step) * (110.57 * step)
        for c in np.unique(cid[cid >= 0]):
            W[i, c] += (cid == c).mean() * area_g
    W /= cells.area_km2.to_numpy()[None, :]
    keep = W.sum(axis=1) > 0
    W = W[keep]
    if np.any(W.sum(axis=0) <= 0):
        raise ValueError("重新分格遺失來源格；請增加取樣密度")
    W /= W.sum(axis=0)
    grid = pd.DataFrame(dict(lon0=g_lon0[keep], lon1=g_lon0[keep] + step,
                             lat0=g_lat0[keep], lat1=g_lat0[keep] + step))
    _REGRID_CACHE.parent.mkdir(parents=True, exist_ok=True)
    if step == 0.1 and n_sub == 5:
        np.savez(_REGRID_CACHE, W=W, lon0=grid.lon0, lon1=grid.lon1, lat0=grid.lat0, lat1=grid.lat1)
    return W, grid


def write_csep_10col_regridded(rate_by_cell_bin: np.ndarray, path: Path,
                               spec: ExperimentSpec = SPEC, depth=(0.0, 40.0)) -> Path:
    """正式版十欄輸出：先把 177 格重新分配到 0.1° 格（保留總期望數），再逐箱寫出。"""
    W, grid = regrid_matrix_01deg()
    edges = magnitude_edges(spec)
    rate01 = W @ np.asarray(rate_by_cell_bin)            # (n_01, 25)
    rows = []
    for g, r in grid.iterrows():
        for k in range(len(edges) - 1):
            rows.append((r.lon0, r.lon1, r.lat0, r.lat1, depth[0], depth[1],
                         edges[k], edges[k + 1], rate01[g, k], 1))
    out = pd.DataFrame(rows, columns=["LON_0", "LON_1", "LAT_0", "LAT_1", "Z_0", "Z_1",
                                      "MAG_0", "MAG_1", "RATE", "FLAG"])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False, float_format="%.6e")
    return path
