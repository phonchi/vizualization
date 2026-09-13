# gdms_toolkit 義大利實驗 API 摘要（自動產生）

## `gdms_toolkit.italy`

義大利標準預報實驗的資料層（教學用）。

實驗設定沿用 Biondini, Rhoades & Gasperini (2023, GJI 234:1681)：
- 目錄：HORUS（Lolli et al. 2020），Mw 均一化，1960 年起。
- 收集區 S：CPTI15 polygon；測試區 R：177 個邊長 30√2 km 的方格。
- 座標：RDN2008 / Italy zone (E-N)，EPSG:7794，以公里為單位。
- 深度 ≤ 40 km；輸入門檻 m0 = 2.45（名目 2.5）；目標門檻 mT = 5.0（有效 4.95）。
- warm-up 1960–1989、learning 1990–2011、testing 2012–2021；三個月滾動窗。

原始檔取自 PyEEPAS 公開 repo 的 data/README 所指 Google Drive 資料夾，
以 scripts/fetch_italy_data.py 下載到 data/cache/italy/（不進版控）。

- `ExperimentSpec(catalog: 'str' = 'HORUS（Lolli et al. 2020），Mw 均一化，1960 年起', depth_max_km: 'float' = 40.0, m0: 'float' = 2.45, mT: 'float' = 5.0, m_max: 'float' = 7.5, mag_bin: 'float' = 0.1, warmup: 'tuple' = (1960, 1989), learning: 'tuple' = (1990, 2011), testing: 'tuple' = (2012, 2021), window_days: 'float' = 91.31, delay_days: 'float' = 50.0, b_value: 'float' = 1.084, n_cells: 'int' = 177, cell_side_km: 'float' = np.float64(42.42640687119285), crs: 'str' = 'EPSG:7794') -> None` — 規格卡：一份預報實驗要先說清楚的所有欄位。
- `bin_targets(targets: 'pd.DataFrame', windows: 'pd.DataFrame | None' = None, spec: 'ExperimentSpec' = SPEC) -> 'np.ndarray'` — 把目標地震數到 (窗, 格, 規模箱) 的觀測計數陣列 ω。
- `cell_of(x_km, y_km, cells: 'pd.DataFrame | None' = None) -> 'np.ndarray'` — 每個點落在哪一格（0..176），不在 R 內回傳 -1。
- `collection_polygon() -> 'pd.DataFrame'` — 收集區 S（CPTI15 polygon）的頂點，含經緯度與公里座標；最後一點重複首點。
- `decimal_days(t: 'pd.Series | pd.Timestamp') -> 'np.ndarray | float'` — 自 1960-01-01 起算的 decimal days（ETAS_R 與 EEPAS 都用這種時間軸）。
- `experiment_catalog(spec: 'ExperimentSpec' = SPEC, end_year: 'int | None' = None) -> 'pd.DataFrame'` — 套用規格卡的資料側篩選後的工作目錄。
- `forecast_windows(spec: 'ExperimentSpec' = SPEC) -> 'pd.DataFrame'` — 測試期的滾動預報窗（decimal days 與日期），共 40 窗。
- `in_collection_region(lon, lat) -> 'np.ndarray'` — 自己算一次點是否在 CPTI15 polygon 內（與 HORUS 的 Geo-CPTI15 旗標對照）。
- `km_to_lonlat(x_km, y_km)` —
- `load_horus(refresh: 'bool' = False) -> 'pd.DataFrame'` — 讀完整 HORUS 目錄（1960 起，約 49 萬筆），不做任何實驗篩選。
- `lonlat_to_km(lon, lat)` — WGS84 經緯度 → RDN2008 Italy zone 公里座標 (x_km 東, y_km 北)。
- `magnitude_edges(spec: 'ExperimentSpec' = SPEC) -> 'np.ndarray'` — 規模箱邊界 5.0, 5.1, …, 7.5（25 箱）。取整規模 mb 落在 [m1, m2)。
- `regrid_matrix_01deg(step: 'float' = 0.1, n_sub: 'int' = 5, refresh: 'bool' = False)` — 177 個投影方格 → 0.1° 經緯度格的面積分配矩陣 W（n_01 × 177）。
- `spec_card(rows: 'list[tuple[str, str]]', title: 'str' = '預報實驗規格卡', pending: 'str' = '待填') -> 'str'` — 把 (欄位, 內容) 清單轉成規格卡 HTML；內容為 None 時顯示「待填」。
- `target_events(cat: 'pd.DataFrame', period: 'str' = 'testing', spec: 'ExperimentSpec' = SPEC) -> 'pd.DataFrame'` — 目標地震：在 R 內、取整規模 ≥ mT、落在指定期間。
- `testing_cells() -> 'pd.DataFrame'` — 測試區 R 的 177 個方格：公里邊界、格心經緯度、四角經緯度（供畫圖）。
- `write_csep_10col(rate_by_cell_bin: 'np.ndarray', path: 'Path', spec: 'ExperimentSpec' = SPEC, depth=(0.0, 40.0)) -> 'Path'` — 以 177 格的經緯度外接框輸出 CSEP 十欄格式（教學示範，非重新分格）。
- `write_csep_10col_regridded(rate_by_cell_bin: 'np.ndarray', path: 'Path', spec: 'ExperimentSpec' = SPEC, depth=(0.0, 40.0)) -> 'Path'` — 正式版十欄輸出：先把 177 格重新分配到 0.1° 格（保留總期望數），再逐箱寫出。
- `year_start_days(year: 'int') -> 'float'` —

## `gdms_toolkit.italy_models`

義大利實驗的四個教學版預報模型：SUP、PPE、ETAS-lite、EEPAS。

所有模型都輸出同一種東西：每個預報窗 × 177 格 × 25 個規模箱的期望地震數，
陣列形狀 (n_windows, 177, 25)。參數全部取自 Biondini, Rhoades & Gasperini
(2023, GJI 234:1681) 表 3（main shocks + aftershocks 資料集）；本模組不做參數擬合。

教學版與論文實作的差異（正文會明說）：
- ETAS 只算「背景 + 已知歷史事件的第一代觸發」的期望率，不模擬後代觸發。
- EEPAS 採等權重版本（論文的 EEPAS-NW），不計算餘震降權 w_i。
- 格內積分：EEPAS 空間核為常態，對方格可用 erf 精確積分；PPE 用細網格取樣；
  ETAS 的尖峰空間核用固定種子的抽樣配置到方格。

- `build_all_forecasts(cat: 'pd.DataFrame | None' = None, periods=('learning', 'testing')) -> 'dict'` — 預先計算四張預報並寫快取；回傳每模型每期的期望總數（供核對論文表 5／表 7）。
- `eepas_forecast(cat: 'pd.DataFrame', windows: 'pd.DataFrame', ppe: 'np.ndarray | None' = None, spec: 'ExperimentSpec' = SPEC, params: 'dict | None' = None) -> 'np.ndarray'` — EEPAS = μ·PPE + Σ_i η_i f_i(t) g_i(m) h_i(x,y)，來源為 S 內 m_b ≥ 2.5 的事件。
- `eepas_pieces(src: 'pd.DataFrame', cells: 'pd.DataFrame', spec: 'ExperimentSpec' = SPEC, params: 'dict | None' = None)` — 與預報窗無關、可預先算好的三塊：η_i、規模箱質量 (n_src,25)、空間格質量 (n_src,177)。
- `eepas_time_mass(t_i: 'np.ndarray', m_i: 'np.ndarray', t1: 'float', t2: 'float', params: 'dict | None' = None) -> 'np.ndarray'` — 對數常態時間核在 [t1, t2] 的質量；t 為 decimal days。
- `etas_forecast(cat: 'pd.DataFrame', windows: 'pd.DataFrame', sup: 'np.ndarray | None' = None, spec: 'ExperimentSpec' = SPEC, params: 'dict | None' = None) -> 'np.ndarray'` — ETAS-lite：ν·(均勻背景) + 已知歷史事件的第一代 Omori 觸發（不含後代）。
- `etas_space_matrix(src: 'pd.DataFrame', cells: 'pd.DataFrame', params: 'dict | None' = None, n_samples: 'int' = 400, seed: 'int' = 20260913) -> 'np.ndarray'` — 尖峰空間核（式 C3）在各格的質量，用固定種子抽樣配置：形狀 (n_src, 177)。
- `get_forecast(model: 'str', period: 'str' = 'testing', cat: 'pd.DataFrame | None' = None, rebuild: 'bool' = False) -> 'np.ndarray'` — 取得快取的預報陣列 (n_windows, 177, 25)；缺快取時報錯；只有 rebuild=True 才計算。
- `learning_windows(spec: 'ExperimentSpec' = SPEC) -> 'pd.DataFrame'` — 學習期也切成同樣長度的窗（用來核對模型在學習期的期望總數）。
- `ppe_forecast(cat: 'pd.DataFrame', windows: 'pd.DataFrame', spec: 'ExperimentSpec' = SPEC, params: 'dict | None' = None) -> 'np.ndarray'` — 鄰近過去地震（PPE）：過去 M≥mT 事件的平滑核 × 1/(t−t0) × GR 分箱。
- `ppe_space_matrix(src: 'pd.DataFrame', cells: 'pd.DataFrame', d: 'float') -> 'np.ndarray'` — 每個來源事件的 Cauchy 核 1/(π(d²+r²)) 在每格的積分，形狀 (n_src, 177)。
- `sup_forecast(cat: 'pd.DataFrame', windows: 'pd.DataFrame', spec: 'ExperimentSpec' = SPEC) -> 'np.ndarray'` — 空間均勻 Poisson：學習期 R 內 M≥mT 的平均率，均攤到面積與時間。

## `gdms_toolkit.csep_teaching`

CSEP 一致性與比較檢驗的教學版實作（純 numpy／scipy，不依賴 pycsep）。

輸入慣例：預報 `rate` 與觀測 `omega` 都是形狀相同的陣列，
通常是 (n_windows, n_cells, n_mag_bins) 或已加總的 (n_cells, n_mag_bins)。
每個檢驗回傳 dict：觀測統計量、模擬分布（或解析分布）、分位數分數。
分位數分數 = 模擬值 ≤ 觀測值的比例（Zechar et al. 2010 的用法）。

- `bill(rate: 'np.ndarray', omega: 'np.ndarray') -> 'np.ndarray'` — 二元對數概似 BILL：只問「箱內有沒有事件」（Bayona et al. 2022）。
- `cl_test(rate, omega, n_sim: 'int' = 10000, seed: 'int' = 1, alpha: 'float' = 0.05) -> 'dict'` — 條件概似檢驗 cL-test：固定總數，比較 jPOLL 與多項模擬分布（單尾，低就不合）。
- `information_gain(rate_a, rate_b, omega) -> 'dict'` — 模型 A 相對基準 B 的資訊增益：逐事件對數概似比，配對 t 區間（T-test）。
- `jpoll(rate, omega) -> 'float'` — 聯合 Poisson 對數概似 jPOLL = Σ POLL。
- `l_test(rate, omega, n_sim: 'int' = 10000, seed: 'int' = 1, alpha: 'float' = 0.05) -> 'dict'` — L-test：不固定總數，模擬目錄的總數也服從 Poisson(Σλ)。
- `m_test(rate, omega, n_sim: 'int' = 10000, seed: 'int' = 1, alpha: 'float' = 0.05, mag_axis: 'int' = -1) -> 'dict'` — M-test：把空間（與時間）加總後只看規模分布。
- `n_test(rate, omega, alpha: 'float' = 0.05) -> 'dict'` — N-test（Poisson 版）：觀測總數是否落在模型總期望數的雙尾區間內。
- `n_test_nb(rate, omega, variance: 'float', alpha: 'float' = 0.05) -> 'dict'` — 負二項版 N-test：均值取模型期望數，變異數另外給（例如歷史十年計數的變異數）。
- `poll(rate: 'np.ndarray', omega: 'np.ndarray') -> 'np.ndarray'` — 逐箱 Poisson 對數概似 POLL：−λ + ω ln λ − ln ω!（λ=0 且 ω=0 時為 0）。
- `power_by_simulation(rate_true, rate_model, n_sim: 'int' = 2000, seed: 'int' = 1, alpha: 'float' = 0.05, test=<function cl_test at 0x79da2bb1d760>, n_inner: 'int' = 1000) -> 'float'` — 統計功效示意：若真實率是 rate_true，檢驗 rate_model 被拒絕的比例。
- `s_test(rate, omega, n_sim: 'int' = 10000, seed: 'int' = 1, alpha: 'float' = 0.05, mag_axis: 'int' = -1) -> 'dict'` — S-test：把規模箱加總後只看空間分布，其餘與 cL-test 相同。

## `gdms_toolkit.viz`

- `plot_forecast_map(rate_by_cell, cells=None, targets=None, polygon=None, title='', log_scale=True, ...)` — 全書統一預報地圖
- `apply_layout(fig, **kw)`、`setup_plotly()`、`ACCENT`、`QUAKE_COLOR`、`SEQUENTIAL`、`PALETTE`

## 快取位置

- `data/cache/italy/horus_clean.csv`：HORUS 乾淨版
- `data/cache/italy_forecasts/{SUP,PPE,ETAS,EEPAS}_{learning,testing}.npy`：形狀 (n_windows, 177, 25)；learning 88 窗、testing 40 窗
- `data/cache/italy/regrid_01deg.npz`：177 格→0.1° 分配矩陣
