# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: kernelspec,jupytext
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 19. 資料產品與品質：取得資料之前先定義問題
#
# 上一章把測站放回地圖，這一章把每一筆資料放回產生它的流程。假設我們要問
# 「花蓮地震後，某口井是否留下持續的水位偏移？」只靠地震清單不夠，還需要
# 水位連續序列，以及可以解釋背景變化的氣壓、降雨和井的紀錄。若問
# 「餘震發生率是否下降」，事件目錄才是主要材料。
#
# 這個選擇會決定後面能做的推論。原始波形保留短時間震動的細節；事件目錄
# 將多站波形整理成時間、位置與規模，便於統計，但不再保留所有訊號。GNSS
# 原始觀測則還沒有直接給出位移。讀懂資料產品，是把物理問題翻譯成統計問題
# 的第一步。
#
# ## 19.1 原始觀測、處理成果與衍生指標
#
# | 想回答的問題 | 主要產品 | 還需要的資訊 |
# |---|---|---|
# | 震波何時到站 | 波形 | 時鐘、頻道、儀器響應 |
# | 某區有多少地震 | 事件目錄 | 規模尺度、完整度、資料版本 |
# | 井水位是否偏移 | 水位序列 | 氣壓、水文背景及井的紀錄 |
# | 地磁變化是否局部 | 多站地磁 | 頻帶、太空天氣、參考站品質 |
# | 測站是否移動 | 已解算座標序列 | 參考框架、誤差及裝置變更 |
#
# 連續序列和事件目錄也可以互相檢查。目錄突然變少時，波形或資料可用率能
# 幫助分辨「真的較安靜」與「當時比較難偵測」。本章原始檔保留取得教材資料
# 的程式；主文聚焦資料的意思，介面與重現細節見{doc}`附錄 G <appendix_g_observations>`。
#
# %% tags=["remove-input", "remove-output"]
import gdms_toolkit as gt

gdms = gt.GDMSSession()
print("登入成功！")

# %% [markdown]
# ## 19.2 觀測視窗要包含背景
#
# 只下載紅線附近的幾分鐘，很容易失去判斷異常的基準。檢查同震振盪需要
# 高頻資料；估計日變化需要完整的日夜；比較季節或長期趨勢則需要更長紀錄。
# 資料範圍影響檔案大小，也決定我們能否檢查其他可能的解釋。
#
# 本書以 2024 花蓮事件為共同案例，保留地下水、地磁、波形及 GNSS 的資料
# 產品。各產品的時間範圍並不完全相同，讀圖時須以圖上日期與原始 metadata
# 為準。把不同產品放上共同時間軸，只能比較它們真正重疊的部分。
#
# %% tags=["remove-input", "remove-output"]
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
# ## 19.3 缺測、延遲與「沒有事件」
#
# 資料未到站、尚未釋出和感測器記為缺值，都是觀測流程的狀態；它們與
# 「儀器有正常記錄，而且沒有看到訊號」不同。事件目錄中沒有事件，也不能
# 單靠空白判定當地完全安靜，還要問那個時段的偵測能力。
#
# 對連續資料，應保留時間戳、缺值位置與原始品質標記。內插可以協助畫圖，
# 卻不是補回一次真實觀測。若缺口剛好跨過主震，內插甚至可能畫出平滑曲線，
# 使本來未知的同震反應看起來像沒有變化。
#
# %% tags=["remove-input", "remove-output"]
import pandas as pd

dl = pd.DataFrame(gt.list_my_downloads(gdms))
dl[["id", "type", "label", "requested_at", "status"]].head(8)

# %% [markdown]
# ## 19.4 檔案格式背後的物理量
#
# 以下是本書快取資料採用的格式。格式名稱說明資料如何儲存，不保證內容已
# 校正，也不保證每個欄位都有有效資料。
#
# | 資料 | 本書使用的格式 | 解讀上的關鍵 |
# |---|---|---|
# | 地下水 | 打包的每日 CSV | 水位、氣壓與水溫分別檢查 |
# | 地磁 | IAGA-2002 | 分量方向、缺值旗標及取樣間隔 |
# | GNSS | RINEX | 偽距與載波相位，需解算才能得到座標 |
# | 地震波形 | miniSEED | 儀器響應另查，counts 不等於物理單位 |
# | 地震目錄 | 事件參數列 | 時間基準、規模尺度及定位品質 |
#
# 本書讀取器會把所使用產品的缺值旗標轉為 `NaN`。例如水位缺測時，水溫仍
# 可能正常；不能因一欄缺失就把整個測站刪掉，也不能把其他欄位的正常當成
# 水位可靠的證據。
#
# ## 19.5 時間和頻道也是資料
#
# 本書地震目錄與多數圖採 UTC；臺灣當地時間為 UTC 加八小時，因此花蓮
# 主震在當地是 4 月 3 日，在 UTC 是 4 月 2 日深夜。GNSS 的時間系統應按
# 檔頭與解算產品查核，不能看到日期就直接假定與其他資料相同。
#
# 頻道代碼則幫我們辨認感測器與方向。同一站的垂直分量和水平分量量到的
# 並不是同一條曲線；把不同方向混接成長序列，可能製造沒有物理意義的階變。
#
# %% tags=["remove-input", "remove-output"]
channels = gt.list_channels(gdms, "CWASN", "HWA")
print(channels)

# %% [markdown]
# ## 19.6 從官方入口走到可引用的資料
#
# GDMS 是本書觀測資料的官方查詢入口。實際申請時應依入口提供的產品說明
# 確認測網、測站與時間範圍，並保留申請條件、取得日期及檔案版本。官方服務
# 與已下載的教材快取是兩個層次：閱讀本書圖表不必重新申請同一批資料。
#
# SEED 頻道通常以三個字元描述頻帶／取樣、來源與方向。`HHZ` 的最後一碼
# 表示垂直分量；`HH1`、`HH2` 則須再查水平軸的實際方位。`HN?` 的來源
# 是加速度感測器，但是否飽和仍受儀器動態範圍限制，不能只憑代碼保證。
# 完整代碼表、資料申請函式與疑難排查移到{doc}`附錄 G <appendix_g_observations>`，
# 有需要重現時再回查。
#
# 接下來從{doc}`地下水位 <03_groundwater>`開始。它很適合練習這套讀法：
# 曲線看起來只是水位升降，背後卻同時包含壓力、潮汐、氣象與含水層的反應。
#
# %% [markdown]
# ## 參考資料與延伸閱讀
#
# - **操作入口・免費網站**：中央氣象署，〈[臺灣地震與地球物理資料管理系統（GDMS）](https://gdms.cwa.gov.tw/)〉。
#   檢視官方資料產品、涵蓋範圍與使用說明；重現時依官方登入流程取得資料，與本書快取分開記錄。
# - **格式查詢・免費檔案**：FDSN，〈[Channel codes](https://docs.fdsn.org/projects/source-identifiers/en/v1.0/channel-codes.html)〉。
#   查閱 band、source 與 subsource codes，對照本章的 HHZ、HH1 與 HH2。重點是頻道代碼描述什麼訊號，不要只把它當成檔名的一部分。
# - **核心論文**：Beyreuther et al.（2010），〈[ObsPy: A Python Toolbox for Seismology](https://doi.org/10.1785/gssrl.81.3.530)〉，*Seismological Research Letters*；[出版學會免費全文](https://www.seismosoc.org/Publications/SRL/SRL_81/srl_81-3_es/)。
#   瞭解波形格式、儀器 metadata 與處理工具如何銜接。論文中的程式介面屬於早期版本，實際寫程式時搭配[目前的官方教學](https://docs.obspy.org/tutorial/)。
