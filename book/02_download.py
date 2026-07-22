# %% [markdown]
# # 2. 資料下載機制
#
# 本章說明 GDMS 的資料取得流程，以及本專案 `gdms_toolkit` 如何把整個流程
# 自動化成幾行 Python。
#
# ## 2.1 GDMS 的下載流程
#
# GDMS（<https://gdms.cwa.gov.tw/>）的資料下載採**會員制＋非同步打包**：
#
# > 免費註冊會員 → 登入 → 填下載表單（測網／測站／時間範圍）
# > → 系統背景打包（數分鐘） → 會員下載清單取得檔案
#
# 1. **註冊**：首頁「登入」→「註冊成為會員」，填 Email 即可，免費。
# 2. **填表單**：「資料下載」頁依資料類型（地震波形／地球物理）填
#    測網、測站、時間範圍送出。
# 3. **取件**：系統打包完成後，到「會員專區 → 下載清單」下載壓縮檔。
#
# ## 2.2 用 gdms_toolkit 自動化
#
# 手動點網頁適合下載一兩筆；要系統性研究就需要程式化。`gdms_toolkit`
# 模擬上述流程：
#
# | 函式 | 做什麼 |
# |---|---|
# | `GDMSSession()` | 登入（帳密讀 `.env`） |
# | `request_geophysical()` | 申請地下水／地磁／GNSS 資料 |
# | `request_waveform()` | 申請連續地震波形（需指定頻道） |
# | `list_my_downloads()` | 查看打包狀態與下載連結 |
# | `fetch_download()` / `wait_and_fetch()` | 下載到 `data/cache/`（自動快取） |
# | `gdms_earthquake_catalog()` | 直接取得地震目錄（DataFrame） |
#
# 先建立登入 session（請先把 `.env.example` 複製成 `.env` 並填入帳密）：

# %%
import gdms_toolkit as gt

gdms = gt.GDMSSession()
print("登入成功！")

# %% [markdown]
# ### 申請資料
#
# 下面示範本書各章使用的資料申請。**注意**：同樣的資料只需申請一次，
# 這裡用「快取檔案不存在才申請」的寫法，避免重複打擾伺服器。

# %%
from gdms_toolkit.download import CACHE_DIR

requests_needed = {
    # label（也是檔名）: (函式, 參數)
    "edu-gw-hualien2024": dict(
        kind="geo", network="GW", station="HWA,CHI,LIU,NAB,TUN,DON",
        start="2024-03-01", end="2024-05-01"),
    "edu-mag-hualien2024": dict(
        kind="geo", network="MAGNET", station="HLN,XCG,YLI",
        start="2024-03-25", end="2024-04-10"),   # HLN/YLI 該期間無資料，僅回傳 XCG
    "edu-mag-csg": dict(
        kind="geo", network="MAGNET", station="CSG",
        start="2024-03-25", end="2024-04-10"),
    "edu-mag-ttn": dict(
        kind="geo", network="MAGNET", station="TTN",
        start="2024-03-25", end="2024-04-10"),
    "edu-gnss-hualien2024": dict(
        kind="geo", network="GNSS", station="HUAL,SHUL",
        start="2024-03-30", end="2024-04-05"),
    "edu-wave-hualien2024": dict(
        kind="wave", network="CWASN", station="HWA,ESL", channel="HH?",
        start="2024-04-02T23:30:00", end="2024-04-03T00:30:00"),
}

for label, p in requests_needed.items():
    cached = list(CACHE_DIR.glob(f"{label}.*"))
    if cached:
        print(f"✓ {label} 已在快取：{cached[0].name}")
        continue
    if p["kind"] == "geo":
        gt.request_geophysical(gdms, p["network"], p["station"],
                               p["start"], p["end"], label=label)
    else:
        gt.request_waveform(gdms, p["network"], p["station"],
                            p["start"], p["end"], channel=p["channel"],
                            label=label)
    gt.wait_and_fetch(gdms, label)

# %% [markdown]
# ### 查看下載清單

# %%
import pandas as pd

dl = pd.DataFrame(gt.list_my_downloads(gdms))
dl[["id", "type", "label", "requested_at", "status"]].head(8)

# %% [markdown]
# ## 2.3 各類資料的檔案格式
#
# | 資料 | 打包格式 | 內容 | 讀取工具 |
# |---|---|---|---|
# | 地下水 | `.tgz` | 每日一個 CSV（`YYYYMMDD.站碼`），1 秒取樣 | `gt.read_groundwater()` |
# | 地磁 | `.tgz` | IAGA-2002 秒資料（`站碼YYYYMMDDdsec.sec`） | `gt.read_geomagnetic()` |
# | GNSS | `.tgz` | RINEX 觀測檔（`站碼DDD0.YYo.gz`） | 第 6 章 |
# | 波形 | `.mseed` | miniSEED 連續波形 | `gt.read_waveform()`（ObsPy） |
# | 地震目錄 | JSON | 事件參數＋定位品質 | `gt.gdms_earthquake_catalog()` |
#
# 缺測值慣例：地下水 `9999`、地磁 `88888`——讀取器都已自動轉為 `NaN`。
#
# ```{admonition} 禮貌原則
# :class: warning
# GDMS 是公共研究資源。批次申請時請：一次申請合理的時間範圍（而不是
# 一天一筆送幾百次）、重複使用 `data/cache/` 的快取、在迴圈中加
# `time.sleep()` 間隔。
# ```
#
# ## 2.4 疑難排解
#
# - **登入失敗**：確認 `.env` 帳密；GDMS 登入 session 約一小時後自動過期，
#   toolkit 每次執行都重新登入，通常不受影響。
# - **打包很久**：長時間範圍＋多測站的申請可能要等十幾分鐘，
#   `wait_and_fetch()` 會自動輪詢。
# - **`Err: No data available`**：該站該時段確實沒資料
#   （儀器故障或尚未釋出），換個時段或測站。
#
# ## 2.5 看懂波形頻道代碼
#
# 申請波形時要指定頻道，這個代碼不是隨便編的。國際上通用的 SEED 命名
# 規則用三個字母描述一個頻道：第一碼是取樣率與頻寬，第二碼是感測器類型，
# 第三碼是方向。先看花蓮氣象站（HWA）實際有哪些頻道：

# %%
channels = gt.list_channels(gdms, "CWASN", "HWA")
print(channels)

# %% [markdown]
# 把常見的第一、二碼對照一下，就能讀懂這些代碼在說什麼：
#
# | 代碼 | 意義 |
# |---|---|
# | `HH?` | 高取樣率（~100 Hz）寬頻地震儀，記錄地表速度 |
# | `EH?` | 高取樣率短週期地震儀 |
# | `HN?` | 高取樣率加速度儀（強震），大地震不會爆表 |
# | `?H?` | 第二碼 H＝高增益地震儀；`?N?` 則是加速度儀 |
# | `??Z / ??N / ??E` | 第三碼是方向：垂直、南北、東西 |
# | `??1 / ??2` | 有些站水平向不對正南北，改用 1、2 編號 |
#
# 選頻道其實就是在選「用哪支感測器、看多快的震動」。研究遠震或微震用
# 寬頻（`HH?`），分析強震動、做工程用途則要加速度（`HN?`）。同一個
# 地震在不同頻道上的長相差很多，第 5 章會實際比較。
#
# 想查完整規則，可參考
# [IRIS 的 SEED 頻道命名說明](https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/)。
