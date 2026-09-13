# 全站教學撰寫規範（內部文件，不列入 _toc.yml）

2026-09-13 版，配合「義大利標準預報實驗」主線的全面重組。撰寫或修改任何一章之前先讀這裡。
讀者版記號摘要與兩套術語對照放在附錄 A。第一部檔名前綴即顯示章號（01–20）；
第二部檔名不動，顯示章號 21–29。

---

## 0. 全書主線與讀者

- 讀者：學過一門基礎統計（機率、期望值、常態分布、簡單迴歸、假設檢定與 p 值的概念）的大學生或碩士新生，沒有地震學、點過程、程式背景。
- 可驗收的承諾：讀完能讀懂常見的網格率預報論文，辨認其假設、證據與需要補讀的部分。
- 主線句：一份義大利 2012–2021 的擬前瞻預報實驗，從讀目錄、訂規格、做四張預報（SUP、PPE、ETAS、EEPAS）到用測試期真實目標地震檢驗。
- 第一部只用義大利資料（HORUS）與合成小圖；台灣資料只出現在第二部。

## 1. 標準預報實驗規格（Biondini et al. 2023；工具箱 `gdms_toolkit.italy.SPEC`）

| 欄位 | 內容 | 首次填入 |
|---|---|---|
| 目錄與規模尺度 | HORUS（Lolli et al. 2020），Mw 均一化，1960 年起；規模取整到 0.1 | 2 |
| 深度 | ≤ 40 km | 2 |
| 收集區 S | CPTI15 polygon（`EEPAS/data/CPTI15.mat`） | 2 |
| 測試區 R | 177 個 30√2 km 方格（`CELLE_ter.mat`，EPSG:7794 公里座標） | 2 |
| 三段期間 | warm-up 1960–1989、learning 1990–2011、testing 2012–2021 | 3 |
| 發報時間軸 | 測試期切 40 個 91.31 天的窗；每窗只用發報時刻前的事件；PPE／EEPAS 另加 50 天 delay | 3 |
| 門檻 | $m_0=2.45$（輸入，名目 2.5）、$m_T=5.0$（目標，有效 4.95）、$M_c$ 由第 8 章估 | 3、8 |
| 規模箱 | 5.0–7.5，寬 0.1，共 25 箱 | 3 |
| 輸出格式 | 每窗 × 177 格 × 25 箱的期望數；CSEP 十欄示範見第 9 章、重新分格見附錄 E | 3、9 |
| 評分方式 | 一致性（N／S／M／cL／L）與比較（IGPE、配對 T）| 16–18 |

- 目標地震數：學習期 27 顆（與論文一致）；測試期在 2024 年版 HORUS 為 25 顆（論文 27 顆，目錄版本差異，第 3 章明說）。
- 論文的 177 格選區用到 2021 年前的資訊，正文要說這是「既有研究區的重演」，不能宣稱實驗設計全在測試前凍結。
- PPE 是 quasi time-dependent；SUP 使用凍結的學習期率；PPE、ETAS、EEPAS 逐窗更新可用歷史。PPE 與 EEPAS 均提前 50 天截止。
- ETAS 教學版＝背景＋已知歷史事件的第一代期望率，**明標為近似**；EEPAS 採等權重（EEPAS-NW 參數）。

## 2. 統一記號

| 概念 | 寫法 | 說明 |
|---|---|---|
| 事件 | $(t_i, x_i, y_i, m_i)$ | 時間（自 1960-01-01 起算的天數）、公里座標、取整規模 |
| 歷史 | $H_t$ | $t$ 之前的全部事件 |
| 條件強度 | $\lambda^*(t,x,y,m)\equiv\lambda(t,x,y,m\mid H_t)$ | 星號＝條件於歷史 |
| 格內期望數 | $\Lambda_{jk}$（窗、格 $j$、規模箱 $k$） | 預報陣列的一個元素 |
| 觀測數 | $\omega_{jk}$ | 同一格箱的目標地震數 |
| 三門檻 | $M_c$（目錄性質）、$m_0$（模型輸入）、$m_T$（預報目標） | 第 3 章擁有定義、第 8 章給數值 |
| GR | $b$；$\beta=b\ln 10$；$s(m)=\beta e^{-\beta(m-m_0)}$ | $\beta$ 只代表 $b\ln10$ |
| Omori 密度 | $g(\tau)=\frac{p-1}{c}\left(1+\frac{\tau}{c}\right)^{-p}$ | 需 $p>1$ |
| 產能 | $\kappa(m)=K e^{\alpha(m-m_0)}$ | Biondini 記號 |
| 分支比 | $n$ | 第 12 章 |
| PPE | $a, d, s$；$h_0(x,y)$ | 第 9 章 |
| EEPAS | $a_M,b_M,\sigma_M$／$a_T,b_T,\sigma_T$／$b_A,\sigma_A$／$\eta(m)$、$\Delta(m)$、$\mu_E$ | 混合權重寫 $\mu_E$，不寫 $\mu$ |
| Ψ | $M_P, T_P, A_P, M_m, C(t)$ | 第 14 章 |
| 檢驗 | POLL、jPOLL、BILL、$\delta_1,\delta_2$、分位數分數 $q$、IGPE | 第 16–18 章 |
| 機率 | 一律 $P(\cdot)$ | $p$ 保留給 Omori 指數與 p 值（後者寫「p 值」） |

衝突改名沿用舊規範：BPT 平均複發時間 $T_r$、aperiodicity $c_v$、Weibull 形狀 $k$、Janus 權重 $\pi_{\rm ETAS}$。

## 3. 公式所有權（其他章只引用不重推）

| 公式／主題 | 擁有章 |
|---|---|
| Poisson 分布與過程、$\Lambda=\int\lambda$、$P(N\ge1)=1-e^{-\Lambda}$、過度離散 | 4 |
| 概似、對數概似、MLE、bootstrap | 5 |
| 分位數↔p 值、對數分數、概似比 | 6 |
| GR 律、$\beta=b\ln10$、Aki 公式與分箱修正 | 7 |
| $M_c$ 三法 | 8 |
| PPE 核與正規化、GR 分配到規模箱 | 9 |
| Omori–Utsu、Båth、產能、GK 除叢 | 10 |
| $\ln L=\sum\ln\lambda^*-\int\lambda^*$、補償子 | 11 |
| ETAS 全式、分支比 | 12 |
| 低維 MLE、隨機除叢 $\rho_{ij}$／$\phi_j$、R–J 特例 | 13 |
| $C(t)$、Ψ 三迴歸 | 14 |
| EEPAS 三核、$\eta(m)$、$\Delta(m)$ | 15 |
| N-test（Poisson／負二項）、Bonferroni | 16 |
| POLL／jPOLL／BILL、S／M／cL／L | 17 |
| IGPE、配對 t 區間、功效 | 18 |
| 凸組合（守恆的是期望數）| 19 |
| 危害積分、survival／hazard 直覺 | 20（推導在附錄 F）|

## 4. 兩套術語（第 3 章對照表；正文以左欄為主，首次出現括號附右欄）

| CSEP／EEPAS（本站用語） | ETAS R 套件（Jalilian 2019） |
|---|---|
| 測試區 R（testing region） | study region |
| 收集區 S（collection／neighbourhood region） | complementary events 所在區 |
| learning period／testing period | study period |
| 目標地震（target） | target events |
| 提供歷史的事件（precursor／trigger） | complementary events |
| decimal days 自 1960-01-01 | decimal days 自 `time.begin` |

polygon 頂點順序：ETAS R 要求反時針、EEPAS 框架要求順時針，第 2 章明講。

## 5. 教學主線與句法

- 主文以核心概念完整、分節清楚為準，長推導移附錄；避免為字數刪去模型關係；圖以教學需要為準；每章末固定一段「本章填入的規格欄位」，未填欄一律寫「待填（第 N 章）」。
- **術語先定義再使用**：每個新名詞第一次出現時，同一段內用一句話定義並附英文；縮寫第一次出現時展開。
- 每段先給一個能抓住的東西（一個數字、一欄資料、一張圖、一條公式），再講判斷與限制。
- 句子控制在 30 字內；主詞明確（「這張圖」「HORUS 目錄」「PPE 模型」），少用「我們」「這」開頭；一段只做一件事。
- 否定式定義不可單獨存在（不能只說「背景率不是觸發項」，要先說它是什麼）。
- 每個模型章末放同一格式的「預報圖＋目標事件疊圖（不評分）」；檢驗章之前不看分數。
- 真實資料、模擬與示意不得混稱；圖說數字一律由程式帶入。
- 零習題、零 quiz、不新增模型調參實驗；可以播放、拖曳資料時間、切換觀測圖層與比較視圖。
- 每頁末列「參考資料與延伸閱讀」，由入門到進階，至少一項免費可讀，附中文導讀。

## 6. 規格卡、定義框與術語卡（HTML／CSS 契約）

- 規格卡：`gdms_toolkit.italy.spec_card(rows)` 產生 `<div class="spec-card">`，以 `display(HTML(...))` 輸出；或 MyST admonition ```{admonition} 規格卡``` 加 `:class: spec-card`。
- 定義框：```{admonition} 定義：條件強度``` 加 `:class: definition`。
- 術語卡（英文對照）：`:class: term-card`。
- 互動示意圖：`display(HTML(open(_static/diagrams/xxx.html).read()))`，樣板見 `book/_static/diagrams/README.md`；SVG 需 `viewBox`＋`max-width:100%`，互動只用 hover 與逐步揭露，無外部資源。
- 深色模式 selector 依 `book/_static/teaching.css` 檔頭註解。

## 7. 排版規則（違反會導致渲染失敗）

1. 粗體前後的標點放粗體外面：寫 `**內文**。`，不要 `**內文。**`。
2. 未啟用 amsmath：多行推導用 `$$\begin{aligned}…\end{aligned}$$`。
3. 數學式內不放中文；中文與行內 `$…$` 間留半形空格。
4. 表格 ≤ 5 欄。
5. 只有跨章引用的式子才編號。
6. 每行約 70–75 個全形字寬折行。

## 8. 程式與圖

- 章首固定 cell：`# %% tags=["remove-input"]` → `from gdms_toolkit.viz import setup_plotly; setup_plotly()`。
- 公開頁程式 cell 一律 `tags=["remove-input"]`；每 cell < 30 秒；亂數 `np.random.default_rng(seed)`。
- 資料一律經 `gdms_toolkit.italy`（`experiment_catalog()`、`target_events()`、`forecast_windows()`）與 `italy_models.get_forecast(model, period)`（讀 `data/cache/italy_forecasts/` 快取，不在章內重算）。檢驗用 `csep_teaching`。
- **記憶體**：任何模擬先估 n_sim × n_bins × 8 bytes；超過幾百 MB 改稀疏或分塊。章內 n_sim ≤ 3,000。
- 圖：Plotly＋`apply_layout()`；`ACCENT` 藍、`QUAKE_COLOR` 紅只給目標地震、`SEQUENTIAL="Blues"`；預報圖用 `viz.plot_forecast_map()` 統一樣式。
- 可用套件：numpy、scipy、pandas、plotly、matplotlib、pyproj。沒有 pycsep、sklearn、statsmodels、geopandas。

## 9. 內容紅線

- **Embargo**：`[2026] PyEEPAS`(+sup)、`[2025] EEPAS_Software` 稿件中的擬合結果、log-likelihood、檢驗分數、加速倍數不得出現；`[2026] Taiwan_EEPAS`(+sup) 與 `EEPAS_TW` 任何結果不得出現，台灣 EEPAS 只能寫「在地化工作正在進行中」。
- 可用：PyEEPAS 公開 GitHub（https://github.com/phonchi/EEPAS）的程式與義大利資料、Biondini 2023 已發表設定與參數、CSEP 格式與兩區域／三期術語。
- 未核實項不得寫入：STEP 的 5 km 網格、USGS「20 分鐘首報／第一年 75 次更新」出處、`generic/sequence-specific/spatially-varying` 是否為 2005 原文用語。STEP 不是 USGS 現行系統（現行為 OAF）。Baker PSHA 是五步驟。
- 引用一律已發表文獻作者與年份；找不到出處的數字不寫。

## 10. 素材對照

`reference/notes/`：`etas.md`→12、13｜`eepas.md`→9、15｜`psi.md`→14｜`testing.md`→16–18｜`ensemble.md`→19｜`oef.md`→1、20｜`stats.md`→7、8｜`psha_step_websearch.md`→20｜`taiwan.md`→29。
本輪工作紀錄：`reference/notes/rewrite_20260913/`（`toc_draft.yml`、`chapter_map.json`、驗證紀錄）。

## 11. 互動科學展與格式驗收

- 術語採 `**中文名詞**（English）`；不得將結尾括號包在粗體中又直接接中文字，避免CommonMark保留星號。
- Matplotlib在匯入前由執行入口設定可寫的`data/cache/matplotlib`；stderr保存在診斷log並使執行失敗，不進入公開輸出。
- 展件依概念採地圖、時間軸、率曲線、分支、等高線、分布與比較矩陣；共享`exhibits.css/js`，不以單一方框模板替代所有圖形。
- 網站用`museum_layout`保留原URL與全文；目錄、搜尋、明暗與閱讀位置可用鍵盤操作。
- 可重建的離線展件在`_static/diagrams/standalone/`；notebook保存內容與靜態圖，完整展件操作在網站或獨立HTML進行。
- 文獻覆蓋位置及核對深度見`reference/notes/refresh_20260913_exhibition/literature_coverage.md`，不得將索引成功當成全文核對完成。
