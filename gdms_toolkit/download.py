"""GDMS 資料下載（需會員登入）。

GDMS 採「申請 → 系統打包 → 會員下載清單取件」的非同步流程：

1. request_geophysical() / request_waveform() 送出申請
   （等同網頁上的下載表單，POST php/sendEqdownload.php）
2. 系統在背景打包（數分鐘），完成後出現在會員下載清單
3. list_my_downloads() 查看清單（POST php/dbconnect/getMemberDownList.php）
   fetch_download() 下載檔案到 data/cache/（本地快取，不重複下載）

實測整理的介面備註：
- 地球物理資料（GW / MAGNET / GNSS）時間格式為 'YYYY-MM-DD'
- 連續波形（CWASN / TSMIP）時間格式為 'YYYY-MM-DDTHH:MM:SS'（UTC），
  且必須指定 channel（如 'HH?'），可用 list_channels() 查詢
- GNSS 需指定 output：'o'（觀測檔）或 'n'（導航檔）
- 波形 output：mseed / fullseed / ascii1 / ascii2 / plot
"""

import re
import time
from pathlib import Path

from .auth import BASE, GDMSSession

CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "cache"


def _send(gdms: GDMSSession, data: dict) -> dict:
    r = gdms.post(f"{BASE}/php/sendEqdownload.php", data=data)
    r.raise_for_status()
    try:
        result = r.json()
    except ValueError:
        raise RuntimeError(f"GDMS 回應非 JSON（可能未登入）：{r.text[:200]}")
    if result.get("status") != 1:
        raise RuntimeError(f"申請失敗：{result}")
    print(f"已送出申請：{data['network']} {data['station']} "
          f"{data['stDatetime']} ~ {data['edDatetime']}（label={data.get('label')}）")
    return result


def request_geophysical(gdms: GDMSSession, network: str, station: str,
                        start: str, end: str, label: str = "edu",
                        output: str = "") -> dict:
    """申請地球物理資料（地下水 GW / 地磁 MAGNET / 衛星定位 GNSS*）。

    station: 站碼，多站以逗號分隔（如 'HWA,CHI'）
    start / end: 'YYYY-MM-DD'（日為單位）
    output: GNSS 需指定 'o'（RINEX 觀測檔）或 'n'（導航檔），其他留空
    """
    network = network.upper()
    if network.startswith("GNSS") and not output:
        output = "o"
    return _send(gdms, {
        "output": output, "network": network, "station": station,
        "stDatetime": start, "edDatetime": end, "label": label,
    })


def request_waveform(gdms: GDMSSession, network: str, station: str,
                     start: str, end: str, channel: str = "HH?",
                     location: str = "all", output: str = "mseed",
                     label: str = "edu") -> dict:
    """申請連續地震波形（CWASN / TSMIP）。

    start / end: 'YYYY-MM-DDTHH:MM:SS'，時間為 UTC
    channel: 頻道代碼，支援萬用字元（'HH?' 寬頻、'HN?' 強震加速度、'EH?' 短週期）
    output: mseed / fullseed / ascii1 / ascii2 / plot
    """
    return _send(gdms, {
        "output": output, "network": network.upper(), "station": station,
        "channel": channel, "location": location,
        "stDatetime": start, "edDatetime": end, "label": label,
    })


def list_channels(gdms: GDMSSession, network: str, station: str,
                  location: str = "all") -> list[str]:
    """查詢某站可用的波形頻道代碼。"""
    r = gdms.post(f"{BASE}/php/dbconnect/getOneStationChannel.php",
                  data={"nw": network.upper(), "station": station,
                        "location": location})
    r.raise_for_status()
    return r.json()


def list_my_downloads(gdms: GDMSSession) -> list[dict]:
    """會員下載清單。每筆含 label、狀態與（打包完成時的）下載網址。"""
    r = gdms.post(f"{BASE}/php/dbconnect/getMemberDownList.php")
    r.raise_for_status()
    items = []
    for row in r.json():
        status_html = row.get("show_status") or ""
        m = re.search(r'href="([^"]+)"', status_html)
        url = m.group(1) if m else None
        if url and not url.startswith("http"):
            url = f"{BASE}/{url.lstrip('./')}"
        status = ("ready" if url else
                  "processing" if "Processing" in status_html else
                  "error" if row.get("error_message") else "unavailable")
        items.append({
            "id": row["id"],
            "type": row["script"]["zh"],
            "requested_at": row["datetime"],
            "detail": row.get("para_detail", ""),
            "label": (re.search(r'"label":"([^"]*)"', row.get("para_detail") or "")
                      or [None, ""])[1],
            "status": status,
            "error": row.get("error_message"),
            "url": url,
        })
    return items


def fetch_download(gdms: GDMSSession, item: dict, dest_dir: Path | str = CACHE_DIR,
                   pause: float = 1.0) -> Path:
    """下載清單中一筆已完成（status='ready'）的項目到 dest_dir。"""
    if not item.get("url"):
        raise RuntimeError(f"此項目尚未打包完成：{item['label']}（{item['status']}）")
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    name = Path(item["url"].split("?")[0]).name or f"gdms_{item['id']}.tar"
    dest = dest_dir / name
    if dest.exists() and dest.stat().st_size > 0:
        print(f"已在快取：{dest}")
        return dest
    time.sleep(pause)  # 禮貌性間隔，避免對伺服器造成負擔
    r = gdms.get(item["url"])
    r.raise_for_status()
    dest.write_bytes(r.content)
    print(f"已下載：{dest}（{len(r.content):,} bytes）")
    return dest


def wait_and_fetch(gdms: GDMSSession, label: str, timeout_min: int = 30,
                   poll_sec: int = 60, dest_dir: Path | str = CACHE_DIR) -> list[Path]:
    """等待指定 label 的申請打包完成並全部下載（教學用的便利函式）。"""
    deadline = time.time() + timeout_min * 60
    while True:
        mine = [i for i in list_my_downloads(gdms) if i["label"] == label]
        if mine and all(i["status"] != "processing" for i in mine):
            break
        if time.time() > deadline:
            raise TimeoutError(f"等待 {label} 打包逾時（{timeout_min} 分鐘）")
        print(f"{label}: 打包中，{poll_sec} 秒後再確認…")
        time.sleep(poll_sec)
    paths = []
    for i in mine:
        if i["status"] == "ready":
            paths.append(fetch_download(gdms, i, dest_dir))
        else:
            print(f"注意：{label} 有一筆無法下載（{i['status']}：{i['error']}）")
    return paths
