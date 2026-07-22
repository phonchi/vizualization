"""測站清單：讀取 data/stations/ 內建 CSV，或直接向 GDMS 更新。

測站清單端點（getStationList.php）不需登入即可存取。
"""

from pathlib import Path

import pandas as pd

from .auth import BASE, make_session
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "stations"

# 測網代碼 → 中文名稱
NETWORKS = {
    "GW": "地下水觀測網",
    "MAGNET": "地磁觀測網",
    "CWASN": "中央氣象署地震觀測網（速度型）",
    "TSMIP": "強地動觀測網（加速度型）",
    "GNSS": "GNSS 衛星定位觀測網（氣象署）",
    "GNSS_IES": "GNSS 衛星定位觀測網（中研院地球所）",
    "GNSS_ETEC": "GNSS 衛星定位觀測網（地震工程研究中心）",
}


def load_stations(network: str | None = None) -> pd.DataFrame:
    """載入內建測站清單。network=None 時合併全部測網。

    >>> gw = load_stations("GW")          # 地下水觀測井
    >>> all_st = load_stations()          # 全部 980 站
    """
    if network is not None:
        network = network.upper()
        if network not in NETWORKS:
            raise ValueError(f"未知測網 {network}，可用：{list(NETWORKS)}")
        nets = [network]
    else:
        nets = list(NETWORKS)
    dfs = [pd.read_csv(DATA_DIR / f"{n}.csv", encoding="utf-8-sig") for n in nets]
    df = pd.concat(dfs, ignore_index=True)
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    return df


def update_stations(network: str, timeout: int = 30) -> pd.DataFrame:
    """直接向 GDMS 抓最新測站清單（免登入），並回傳 DataFrame。"""
    r = make_session().post(
        f"{BASE}/php/dbconnect/getStationList.php",
        data={"network_code": network.upper()}, timeout=timeout,
    )
    r.raise_for_status()
    df = pd.DataFrame(r.json())
    df.insert(0, "network", network.upper())
    return df
