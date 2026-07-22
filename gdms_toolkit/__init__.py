"""gdms_toolkit — 台灣地球物理資料教學工具箱

提供 GDMS（中央氣象署臺灣地震與地球物理資料管理系統）與公開資料源的
下載、讀取輔助函式，配合 book/ 目錄下的教學 notebook 使用。
"""

from .stations import load_stations, NETWORKS
from .auth import GDMSSession
from .download import (request_geophysical, request_waveform, list_channels,
                       list_my_downloads, fetch_download, wait_and_fetch)
from .catalog import gdms_earthquake_catalog
from .readers import (read_gdms_timeseries, read_waveform,
                      read_groundwater, read_geomagnetic)

__all__ = [
    "load_stations", "NETWORKS", "GDMSSession",
    "request_geophysical", "request_waveform", "list_channels",
    "list_my_downloads", "fetch_download", "wait_and_fetch",
    "gdms_earthquake_catalog",
    "read_gdms_timeseries", "read_waveform",
    "read_groundwater", "read_geomagnetic",
]
