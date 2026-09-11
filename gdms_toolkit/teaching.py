"""Deterministic synthetic examples shared by the introductory teaching pages.

These are illustrations conditional on a fixed event count, not observed earthquakes
or an ETAS fit. Keep the seed and construction explicit for reproducibility.
"""
import numpy as np
import pandas as pd


def learning_catalog(kind="cluster", seed=20260911):
    """Return 80 marked illustrative events in [0, 120) days.

    ``poisson`` is homogeneous Poisson conditional on N=80; ``regular`` is
    jittered regular timing; ``cluster`` is 16 groups of 5 with circular edges.
    Locations and magnitudes are illustrative marks with no causal interpretation.
    """
    rng = np.random.default_rng(seed)
    n, horizon = 80, 120.0
    if kind == "poisson":
        times = rng.uniform(0, horizon, n)
    elif kind == "regular":
        times = (np.arange(n) + 0.5) * horizon / n + rng.uniform(-0.2, 0.2, n)
    elif kind == "cluster":
        centers = rng.uniform(0, horizon, 16)
        times = (np.repeat(centers, 5) + rng.normal(0, 0.65, n)) % horizon
    else:
        raise ValueError(f"Unknown illustrative catalog kind: {kind}")
    times.sort()
    return pd.DataFrame({
        "day": times,
        "x_km": rng.normal(0, 18, n),
        "y_km": rng.normal(0, 8, n),
        "m": 3.0 + rng.exponential(1 / np.log(10), n),
    })
