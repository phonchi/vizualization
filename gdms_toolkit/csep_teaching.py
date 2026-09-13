"""CSEP 一致性與比較檢驗的教學版實作（純 numpy／scipy，不依賴 pycsep）。

輸入慣例：預報 `rate` 與觀測 `omega` 都是形狀相同的陣列，
通常是 (n_windows, n_cells, n_mag_bins) 或已加總的 (n_cells, n_mag_bins)。
每個檢驗回傳 dict：觀測統計量、模擬分布（或解析分布）、分位數分數。
分位數分數 = 模擬值 ≤ 觀測值的比例（Zechar et al. 2010 的用法）。
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import gammaln


# ---------------------------------------------------------------- 分數 ----
def poll(rate: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """逐箱 Poisson 對數概似 POLL：−λ + ω ln λ − ln ω!（λ=0 且 ω=0 時為 0）。"""
    lam = np.asarray(rate, float)
    w = np.asarray(omega, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(w > 0, w * np.log(lam), 0.0)
    return -lam + term - gammaln(w + 1)


def jpoll(rate, omega) -> float:
    """聯合 Poisson 對數概似 jPOLL = Σ POLL。"""
    return float(poll(rate, omega).sum())


def bill(rate: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """二元對數概似 BILL：只問「箱內有沒有事件」（Bayona et al. 2022）。"""
    lam = np.asarray(rate, float)
    x = (np.asarray(omega) > 0).astype(float)
    p1 = -np.expm1(-lam)
    with np.errstate(divide="ignore"):
        return np.where(x > 0, np.log(p1), -lam)


# ---------------------------------------------------------------- N ----
def n_test(rate, omega, alpha: float = 0.05) -> dict:
    """N-test（Poisson 版）：觀測總數是否落在模型總期望數的雙尾區間內。

    回傳 δ1 = P(N ≥ n_obs)（低估的證據）、δ2 = P(N ≤ n_obs)（高估的證據）。
    """
    n_fore = float(np.sum(rate))
    n_obs = int(np.sum(omega))
    delta1 = stats.poisson.sf(n_obs - 1, n_fore)
    delta2 = stats.poisson.cdf(n_obs, n_fore)
    lo, hi = stats.poisson.ppf([alpha / 2, 1 - alpha / 2], n_fore)
    return dict(n_forecast=n_fore, n_obs=n_obs, delta1=float(delta1), delta2=float(delta2),
                interval=(int(lo), int(hi)), passed=(delta1 > alpha / 2) and (delta2 > alpha / 2))


def n_test_nb(rate, omega, variance: float, alpha: float = 0.05) -> dict:
    """負二項版 N-test：均值取模型期望數，變異數另外給（例如歷史十年計數的變異數）。"""
    mu = float(np.sum(rate))
    n_obs = int(np.sum(omega))
    if variance <= mu:
        raise ValueError("負二項需要變異數 > 均值；否則就用 Poisson 版")
    p = mu / variance                 # scipy 參數化：n = 大小, p = 成功機率
    n = mu * p / (1 - p)
    delta1 = stats.nbinom.sf(n_obs - 1, n, p)
    delta2 = stats.nbinom.cdf(n_obs, n, p)
    lo, hi = stats.nbinom.ppf([alpha / 2, 1 - alpha / 2], n, p)
    return dict(n_forecast=mu, n_obs=n_obs, variance=variance, delta1=float(delta1),
                delta2=float(delta2), interval=(int(lo), int(hi)),
                passed=(delta1 > alpha / 2) and (delta2 > alpha / 2))


# ---------------------------------------------------------------- 模擬 ----
def _sim_jpoll(lam: np.ndarray, n_events, n_sim: int, rng) -> np.ndarray:
    """模擬 n_sim 份目錄的 jPOLL，只追蹤有事件的箱（記憶體 O(n_sim × n_events)）。

    n_events 可為整數（每份目錄固定總數，cL／S／M 用）或長度 n_sim 的陣列（L-test 用）。
    jPOLL = −Σλ + Σ_事件 ln λ_箱 − Σ_箱 ln ω_箱!，最後一項只在同箱多顆時不為零。
    """
    total = lam.sum()
    if total == 0:
        counts = np.broadcast_to(np.asarray(n_events, int), (n_sim,))
        return np.where(counts == 0, 0.0, -np.inf)
    # Build the categorical CDF once; rng.choice(p=...) rebuilds it per catalogue.
    cumulative = np.cumsum(lam / total)
    cumulative[-1] = 1.0
    logl = np.log(np.where(lam > 0, lam, 1.0))
    const = -lam.sum()
    n_events = np.broadcast_to(np.asarray(n_events, int), (n_sim,))
    out = np.empty(n_sim)
    for i, n in enumerate(n_events):
        if n == 0:
            out[i] = const
            continue
        idx = np.searchsorted(cumulative, rng.random(int(n)), side="right")
        _, cnt = np.unique(idx, return_counts=True)
        out[i] = const + logl[idx].sum() - gammaln(cnt + 1).sum()
    return out


def _quantile_score(sim: np.ndarray, obs: float) -> float:
    tolerance = 32 * np.finfo(float).eps * max(1.0, abs(obs)) if np.isfinite(obs) else 0.0
    return float(np.mean(sim <= obs + tolerance))


def cl_test(rate, omega, n_sim: int = 10000, seed: int = 1, alpha: float = 0.05) -> dict:
    """條件概似檢驗 cL-test：固定總數，比較 jPOLL 與多項模擬分布（單尾，低就不合）。"""
    rng = np.random.default_rng(seed)
    lam = np.asarray(rate, float).ravel()
    w = np.asarray(omega).ravel()
    n_obs = int(w.sum())
    if lam.sum() == 0 and n_obs:
        raise ValueError("零總率無法條件在正事件數上")
    lam_c = lam * (n_obs / lam.sum()) if lam.sum() else lam.copy()
    obs = jpoll(lam_c, w)
    sim_ll = _sim_jpoll(lam_c, n_obs, n_sim, rng)
    q = _quantile_score(sim_ll, obs)
    return dict(statistic=obs, simulated=sim_ll, quantile=q, passed=q > alpha)


def l_test(rate, omega, n_sim: int = 10000, seed: int = 1, alpha: float = 0.05) -> dict:
    """L-test：不固定總數，模擬目錄的總數也服從 Poisson(Σλ)。"""
    rng = np.random.default_rng(seed)
    lam = np.asarray(rate, float).ravel()
    w = np.asarray(omega).ravel()
    obs = jpoll(lam, w)
    n_sims = rng.poisson(lam.sum(), size=n_sim)
    sim_ll = _sim_jpoll(lam, n_sims, n_sim, rng)
    q = _quantile_score(sim_ll, obs)
    return dict(statistic=obs, simulated=sim_ll, quantile=q, passed=q > alpha)


def s_test(rate, omega, n_sim: int = 10000, seed: int = 1, alpha: float = 0.05,
           mag_axis: int = -1) -> dict:
    """S-test：把規模箱加總後只看空間分布，其餘與 cL-test 相同。"""
    return cl_test(np.sum(rate, axis=mag_axis), np.sum(omega, axis=mag_axis),
                   n_sim=n_sim, seed=seed, alpha=alpha)


def m_test(rate, omega, n_sim: int = 10000, seed: int = 1, alpha: float = 0.05,
           mag_axis: int = -1) -> dict:
    """M-test：把空間（與時間）加總後只看規模分布。"""
    axes = tuple(a for a in range(np.ndim(rate)) if a != (mag_axis % np.ndim(rate)))
    return cl_test(np.sum(rate, axis=axes), np.sum(omega, axis=axes),
                   n_sim=n_sim, seed=seed, alpha=alpha)


# ---------------------------------------------------------------- 比較 ----
def information_gain(rate_a, rate_b, omega) -> dict:
    """模型 A 相對基準 B 的資訊增益：逐事件對數概似比，配對 t 區間（T-test）。

    IGPE = [ln L_A − ln L_B] / N；配對結構來自「同一顆目標地震在兩個模型下的分數」。
    回傳每顆事件的 X_A − X_B、IGPE、95% t 區間與 Poisson 部分 (N_B − N_A)/N。
    """
    la = np.asarray(rate_a, float).ravel()
    lb = np.asarray(rate_b, float).ravel()
    w = np.asarray(omega).ravel()
    idx = np.repeat(np.arange(w.size), w.astype(int))          # 每顆事件所在箱
    n = idx.size
    if n == 0:
        raise ValueError("每事件資訊增益需要至少一顆目標地震")
    with np.errstate(divide="ignore"):
        x = np.log(la[idx]) - np.log(lb[idx])                   # 逐事件對數率比
    rate_term = (lb.sum() - la.sum()) / n
    igpe = x.mean() + rate_term
    if n > 1:
        s = x.std(ddof=1)
        half = stats.t.ppf(0.975, n - 1) * s / np.sqrt(n)
    else:
        s, half = np.nan, np.nan
    return dict(igpe=float(igpe), per_event=x, rate_term=float(rate_term), n=n,
                ci=(float(igpe - half), float(igpe + half)),
                significant=bool(np.isfinite(half) and (igpe - half > 0 or igpe + half < 0)))


def power_by_simulation(rate_true, rate_model, n_sim: int = 2000, seed: int = 1,
                        alpha: float = 0.05, test=cl_test, n_inner: int = 1000) -> float:
    """統計功效示意：若真實率是 rate_true，檢驗 rate_model 被拒絕的比例。"""
    rng = np.random.default_rng(seed)
    lam_true = np.asarray(rate_true, float).ravel()
    rejected = 0
    for k in range(n_sim):
        w = rng.poisson(lam_true)
        if w.sum() == 0:
            continue
        res = test(rate_model, w.reshape(np.shape(rate_model)), n_sim=n_inner, seed=int(rng.integers(1 << 31)), alpha=alpha)
        rejected += (not res["passed"])
    return rejected / n_sim
