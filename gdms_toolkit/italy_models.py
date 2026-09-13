"""義大利實驗的四個教學版預報模型：SUP、PPE、ETAS-lite、EEPAS。

所有模型都輸出同一種東西：每個預報窗 × 177 格 × 25 個規模箱的期望地震數，
陣列形狀 (n_windows, 177, 25)。參數全部取自 Biondini, Rhoades & Gasperini
(2023, GJI 234:1681) 表 3（main shocks + aftershocks 資料集）；本模組不做參數擬合。

教學版與論文實作的差異（正文會明說）：
- ETAS 只算「背景 + 已知歷史事件的第一代觸發」的期望率，不模擬後代觸發。
- EEPAS 採等權重版本（論文的 EEPAS-NW），不計算餘震降權 w_i。
- 格內積分：EEPAS 空間核為常態，對方格可用 erf 精確積分；PPE 用細網格取樣；
  ETAS 的尖峰空間核用固定種子的抽樣配置到方格。
"""

from __future__ import annotations

import json
import time as _time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import erf
from scipy.stats import norm

from . import italy
from .italy import SPEC, BETA, ExperimentSpec

# ---------------------------------------------------------------- 參數 ----
PARAMS = {
    "SUP": dict(b=1.084),
    "PPE": dict(a=0.62, d=30.0, s=9.0e-13, m_ref=5.0,
                source="Biondini et al. 2023, Table 3 (PPE)"),
    "EEPAS": dict(aM=1.22, bM=1.0, sM=0.25, aT=2.55, bT=0.35, sT=0.150,
                  bA=0.52, sA=1.00, mu=0.18, delay=50.0,
                  source="Biondini et al. 2023, Table 3 (EEPAS-NW, 等權重)"),
    "ETAS": dict(K=0.029, c=0.004, p=1.042, D=1.04, gamma=0.45, alpha=1.12, q=1.5,
                 nu=0.264, m_c=2.5,
                 source="Biondini et al. 2023, Table 3 (ETAS-SUP)；時間單位以天處理"),
}
MODELS = ("SUP", "PPE", "ETAS", "EEPAS")


def _gr_bin_mass(edges: np.ndarray, m_ref: float, beta: float = BETA) -> np.ndarray:
    """GR 密度 β e^{-β(m-m_ref)} 在各規模箱的質量（相對 m ≥ m_ref 的總數）。"""
    return np.exp(-beta * (edges[:-1] - m_ref)) - np.exp(-beta * (edges[1:] - m_ref))


# ---------------------------------------------------------------- SUP ----
def sup_forecast(cat: pd.DataFrame, windows: pd.DataFrame,
                 spec: ExperimentSpec = SPEC) -> np.ndarray:
    """空間均勻 Poisson：學習期 R 內 M≥mT 的平均率，均攤到面積與時間。"""
    cells = italy.testing_cells()
    edges = italy.magnitude_edges(spec)
    targets = italy.target_events(cat, "learning", spec)
    t_learn = italy.year_start_days(spec.learning[1] + 1) - italy.year_start_days(spec.learning[0])
    mu0 = len(targets) / (cells.area_km2.sum() * t_learn)      # 每 km² 每天，M≥mT
    mass = _gr_bin_mass(edges, spec.mT)
    dt = (windows.t2 - windows.t1).to_numpy()
    rate = mu0 * dt[:, None, None] * cells.area_km2.to_numpy()[None, :, None] * mass[None, None, :]
    return rate


# ---------------------------------------------------------------- PPE ----
def _cell_sample_points(cells: pd.DataFrame, n: int = 12):
    """每格 n×n 個中點取樣，回傳 (177, n*n) 的 x、y 與每點面積。"""
    u = (np.arange(n) + 0.5) / n
    xs = cells.x0.to_numpy()[:, None] + (cells.x1 - cells.x0).to_numpy()[:, None] * u[None, :]
    ys = cells.y0.to_numpy()[:, None] + (cells.y1 - cells.y0).to_numpy()[:, None] * u[None, :]
    X = np.repeat(xs, n, axis=1)
    Y = np.tile(ys, (1, n))
    dA = (cells.area_km2 / (n * n)).to_numpy()
    return X, Y, dA


def ppe_space_matrix(src: pd.DataFrame, cells: pd.DataFrame, d: float) -> np.ndarray:
    """每個來源事件的 Cauchy 核 1/(π(d²+r²)) 在每格的積分，形狀 (n_src, 177)。"""
    X, Y, dA = _cell_sample_points(cells)
    out = np.zeros((len(src), len(cells)))
    xs, ys = src.x_km.to_numpy(), src.y_km.to_numpy()
    for j in range(len(cells)):
        r2 = (X[j][None, :] - xs[:, None]) ** 2 + (Y[j][None, :] - ys[:, None]) ** 2
        out[:, j] = (1.0 / (np.pi * (d * d + r2))).sum(axis=1) * dA[j]
    return out


def ppe_forecast(cat: pd.DataFrame, windows: pd.DataFrame,
                 spec: ExperimentSpec = SPEC, params: dict | None = None) -> np.ndarray:
    """鄰近過去地震（PPE）：過去 M≥mT 事件的平滑核 × 1/(t−t0) × GR 分箱。"""
    P = PARAMS["PPE"] if params is None else params
    cells = italy.testing_cells()
    edges = italy.magnitude_edges(spec)
    src = cat.loc[cat.mb >= P["m_ref"]].reset_index(drop=True)      # S 內所有 M≥mT
    space = ppe_space_matrix(src, cells, P["d"])                    # (n_src, 177)
    w = P["a"] * (src.mb.to_numpy() - P["m_ref"])                   # 規模權重 a(m_i − m_T)
    mass = _gr_bin_mass(edges, P["m_ref"])
    delay = PARAMS["EEPAS"]["delay"]
    out = np.zeros((len(windows), len(cells), len(mass)))
    for i, (t1, t2) in enumerate(zip(windows.t1, windows.t2)):
        use = src.t_days.to_numpy() < t1 - delay
        h = (w[use, None] * space[use]).sum(axis=0) + P["s"] * cells.area_km2.to_numpy()
        out[i] = np.log(t2 / t1) * h[:, None] * mass[None, :]
    return out


# ---------------------------------------------------------------- EEPAS ----
def eepas_pieces(src: pd.DataFrame, cells: pd.DataFrame, spec: ExperimentSpec = SPEC,
                 params: dict | None = None):
    """與預報窗無關、可預先算好的三塊：η_i、規模箱質量 (n_src,25)、空間格質量 (n_src,177)。"""
    P = PARAMS["EEPAS"] if params is None else params
    m_i = src.mb.to_numpy()
    edges = italy.magnitude_edges(spec)
    # η（等權重：E[w]=1）
    eta = (P["bM"] * (1 - P["mu"])) * np.exp(-BETA * (P["aM"] + (P["bM"] - 1) * m_i
                                                     + BETA * P["sM"] ** 2 / 2))
    # 規模核 g(m|m_i)/Δ(m)，每箱 5 點中點積分
    sub = 5
    mm = (edges[:-1, None] + (np.arange(sub)[None, :] + 0.5) / sub * spec.mag_bin).ravel()
    delta = norm.cdf((mm - P["aM"] - P["bM"] * spec.m0 - P["sM"] ** 2 * BETA) / P["sM"])
    g = norm.pdf((mm[None, :] - P["aM"] - P["bM"] * m_i[:, None]) / P["sM"]) / P["sM"]
    mag = (g / delta[None, :]).reshape(len(m_i), len(edges) - 1, sub).sum(axis=2) * (spec.mag_bin / sub)
    # 空間核：圓對稱常態，σ = σ_A 10^{b_A m_i / 2}，方格可精確積分
    sig = P["sA"] * 10 ** (P["bA"] * m_i / 2)
    x, y = src.x_km.to_numpy(), src.y_km.to_numpy()
    r2 = np.sqrt(2)
    ex = (erf((cells.x1.to_numpy()[None, :] - x[:, None]) / (r2 * sig[:, None]))
          - erf((cells.x0.to_numpy()[None, :] - x[:, None]) / (r2 * sig[:, None])))
    ey = (erf((cells.y1.to_numpy()[None, :] - y[:, None]) / (r2 * sig[:, None]))
          - erf((cells.y0.to_numpy()[None, :] - y[:, None]) / (r2 * sig[:, None])))
    space = 0.25 * ex * ey
    return eta, mag, space


def eepas_time_mass(t_i: np.ndarray, m_i: np.ndarray, t1: float, t2: float,
                    params: dict | None = None) -> np.ndarray:
    """對數常態時間核在 [t1, t2] 的質量；t 為 decimal days。"""
    P = PARAMS["EEPAS"] if params is None else params
    mu_t = P["aT"] + P["bT"] * m_i
    tau_l = np.maximum(t1 - t_i, 1e-9)
    tau_u = np.maximum(t2 - t_i, 1e-9)
    return 0.5 * (erf((np.log10(tau_u) - mu_t) / (np.sqrt(2) * P["sT"]))
                  - erf((np.log10(tau_l) - mu_t) / (np.sqrt(2) * P["sT"])))


def eepas_forecast(cat: pd.DataFrame, windows: pd.DataFrame, ppe: np.ndarray | None = None,
                   spec: ExperimentSpec = SPEC, params: dict | None = None) -> np.ndarray:
    """EEPAS = μ·PPE + Σ_i η_i f_i(t) g_i(m) h_i(x,y)，來源為 S 內 m_b ≥ 2.5 的事件。"""
    P = PARAMS["EEPAS"] if params is None else params
    cells = italy.testing_cells()
    src = cat.loc[cat.mb >= np.round(spec.m0 + 0.05, 2)].reset_index(drop=True)
    eta, mag, space = eepas_pieces(src, cells, spec, P)
    t_i, m_i = src.t_days.to_numpy(), src.mb.to_numpy()
    ppe = ppe_forecast(cat, windows, spec) if ppe is None else ppe
    out = np.zeros_like(ppe)
    for i, (t1, t2) in enumerate(zip(windows.t1, windows.t2)):
        use = t_i <= t1 - P["delay"]
        tm = eepas_time_mass(t_i[use], m_i[use], t1, t2, P) * eta[use]
        out[i] = (space[use] * tm[:, None]).T @ mag[use] + P["mu"] * ppe[i]
    return out


# ---------------------------------------------------------------- ETAS ----
def etas_space_matrix(src: pd.DataFrame, cells: pd.DataFrame, params: dict | None = None,
                      n_samples: int = 400, seed: int = 20260913) -> np.ndarray:
    """尖峰空間核（式 C3）在各格的質量，用固定種子抽樣配置：形狀 (n_src, 177)。"""
    P = PARAMS["ETAS"] if params is None else params
    rng = np.random.default_rng(seed)
    scale2 = P["D"] ** 2 * np.exp(P["gamma"] * (src.mb.to_numpy() - P["m_c"]))
    out = np.zeros((len(src), len(cells)))
    u = rng.random((len(src), n_samples))
    th = rng.random((len(src), n_samples)) * 2 * np.pi
    r = np.sqrt(scale2[:, None] * (u ** (1.0 / (1.0 - P["q"])) - 1.0))
    X = src.x_km.to_numpy()[:, None] + r * np.cos(th)
    Y = src.y_km.to_numpy()[:, None] + r * np.sin(th)
    cid = italy.cell_of(X.ravel(), Y.ravel(), cells).reshape(X.shape)
    for j in range(len(cells)):
        out[:, j] = (cid == j).mean(axis=1)
    return out


def etas_forecast(cat: pd.DataFrame, windows: pd.DataFrame, sup: np.ndarray | None = None,
                  spec: ExperimentSpec = SPEC, params: dict | None = None) -> np.ndarray:
    """ETAS-lite：ν·(均勻背景) + 已知歷史事件的第一代 Omori 觸發（不含後代）。"""
    P = PARAMS["ETAS"] if params is None else params
    cells = italy.testing_cells()
    edges = italy.magnitude_edges(spec)
    src = cat.loc[cat.mb >= P["m_c"]].reset_index(drop=True)
    space = etas_space_matrix(src, cells, P)
    t_i, m_i = src.t_days.to_numpy(), src.mb.to_numpy()
    prod = P["K"] * np.exp(P["alpha"] * (m_i - P["m_c"]))
    mass = _gr_bin_mass(edges, P["m_c"])
    # 背景：M≥m_c 的均勻率 × ν，再用 GR 從 m_c 換算到目標箱
    learn = cat.loc[(cat.period == "learning") & cat.in_R & (cat.mb >= P["m_c"])]
    t_learn = italy.year_start_days(spec.learning[1] + 1) - italy.year_start_days(spec.learning[0])
    mu0 = len(learn) / (cells.area_km2.sum() * t_learn)
    dt = (windows.t2 - windows.t1).to_numpy()
    bg = P["nu"] * mu0 * dt[:, None, None] * cells.area_km2.to_numpy()[None, :, None] * mass[None, None, :]
    out = np.zeros_like(bg)
    for i, (t1, t2) in enumerate(zip(windows.t1, windows.t2)):
        use = t_i < t1
        om = ((t1 - t_i[use] + P["c"]) ** (1 - P["p"]) - (t2 - t_i[use] + P["c"]) ** (1 - P["p"])) / (P["p"] - 1)
        h = (space[use] * (prod[use] * om)[:, None]).sum(axis=0)
        out[i] = bg[i] + h[:, None] * mass[None, :]
    return out


# ---------------------------------------------------------------- 快取 ----
def learning_windows(spec: ExperimentSpec = SPEC) -> pd.DataFrame:
    """學習期也切成同樣長度的窗（用來核對模型在學習期的期望總數）。"""
    t_start = italy.year_start_days(spec.learning[0])
    total = italy.year_start_days(spec.learning[1] + 1) - t_start
    n = int(round(total / spec.window_days))
    t1 = t_start + spec.window_days * np.arange(n)
    t2 = np.minimum(t1 + spec.window_days, italy.year_start_days(spec.learning[1] + 1))
    return pd.DataFrame(dict(window=np.arange(1, n + 1), t1=t1, t2=t2,
                             start=italy.T0 + pd.to_timedelta(t1, unit="D"),
                             end=italy.T0 + pd.to_timedelta(t2, unit="D")))


def _cache_path(model: str, period: str) -> Path:
    return italy.FORECAST_DIR / f"{model}_{period}.npy"


def get_forecast(model: str, period: str = "testing", cat: pd.DataFrame | None = None,
                 rebuild: bool = False) -> np.ndarray:
    """讀取預報快取；只有明確指定 rebuild=True 才執行模型計算。"""
    path = _cache_path(model, period)
    if path.exists() and not rebuild:
        return np.load(path)
    if not rebuild:
        raise FileNotFoundError(
            f"缺少預報快取 {path}；請先執行 python scripts/fetch_italy_data.py --forecasts")
    cat = italy.experiment_catalog() if cat is None else cat
    windows = italy.forecast_windows() if period == "testing" else learning_windows()
    fn = {"SUP": sup_forecast, "PPE": ppe_forecast, "EEPAS": eepas_forecast, "ETAS": etas_forecast}[model]
    arr = fn(cat, windows)
    italy.FORECAST_DIR.mkdir(parents=True, exist_ok=True)
    np.save(path, arr)
    return arr


def build_all_forecasts(cat: pd.DataFrame | None = None, periods=("learning", "testing")) -> dict:
    """預先計算四張預報並寫快取；回傳每模型每期的期望總數（供核對論文表 5／表 7）。"""
    cat = italy.experiment_catalog() if cat is None else cat
    summary = {}
    for period in periods:
        for model in MODELS:
            t0 = _time.perf_counter()
            arr = get_forecast(model, period, cat, rebuild=True)
            summary[f"{model}_{period}"] = dict(total=float(arr.sum()),
                                                seconds=round(_time.perf_counter() - t0, 1))
            print(f"{model:6s} {period:9s} 期望總數 {arr.sum():7.2f}  ({summary[f'{model}_{period}']['seconds']} s)")
    meta = dict(built=pd.Timestamp.now().isoformat(timespec="seconds"), params=PARAMS,
                spec=SPEC.__dict__, summary=summary)
    (italy.FORECAST_DIR / "META.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2, default=str))
    return summary
