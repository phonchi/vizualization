# 地震預報教學網站全面重組計畫（義大利標準預報實驗主線）

## Context

使用者反映現行教學網站（`book/`，Jupyter Book；第一部 17 章＋第二部 9 章＋附錄 7 篇）「讀起來十分拗口」。
讀者是學過基礎統計、沒有地震學背景的學生；目標是讀完能看懂 CSEP 式地震預報論文，並理解「一個數字算得出來，不等於它站得住腳」。
唯讀盤點證實四類症狀：太抽象缺例子與定義；統計與地震術語先用再解釋（jPOLL／cL／SUP／CSEP 全站從未展開，「宣告區域 polygon」全站 0 次，T-test 配對結構、AIC 本體、分位數↔p 值橋接都沒有）；長句書面化、主詞不明；概念順序錯置（`foundation_inference` 排在所有使用概似／MLE／bootstrap 的章之後；第 1 章就用非齊次 Poisson、自激發、條件強度）。學生也看不到一個完整的預報實驗：基礎章用合成目錄、各模型章各自挑資料。前兩輪只做詞彙層 humanize 與圖內術語修正，未動章序，所以這次從結構、資料主線與句法一起重來。

## 已定決策（grill 結果，不再重議）

| 項目 | 決定 |
|---|---|
| 範圍 | 第一部全部重寫＋附錄重組；第二部只調銜接，另把第一部的台灣內容**全部併入** `23_taiwan_outlook` |
| 重寫深度 | 重組＋每章正文全面重寫；可用的圖與程式 cell 搬用；程式一律摺疊（`remove-input`），觀念為主，零習題 |
| 貫穿主線 | 一張「標準預報實驗規格卡」貫穿全部模型章與檢驗章；**義大利**為主線；第一部不再出現台灣資料 |
| 義大利規格 | Biondini et al. (2023)：測試區 R＝177 個 30√2 km 方格（`EEPAS/data/CELLE_ter.mat`，座標實為 RDN2008 投影公里，需反投影）、鄰域／收集區 S＝CPTI15 polygon（`EEPAS/data/CPTI15.mat`）、HORUS 目錄、深度 <40 km、m0＝2.45、mT＝5.0、學習期 1990–2011、pseudo-prospective 測試期 2012–2021、三個月滾動預報窗、規模箱 0.1（5.0–7.5，25 箱）、CSEP 十欄格式（LON_0 LON_1 LAT_0 LAT_1 Z_0 Z_1 MAG_0 MAG_1 RATE FLAG） |
| 評分格 | 教學評分用 177 格 × 25 箱（對得上 Biondini）；0.1° 格只做「匯出格式」示範 |
| 術語 | 教兩套並對照：CSEP／EEPAS 的 warm-up／learning／testing period、testing region R ⊂ collection region S；ETAS R 套件（Jalilian 2019 JSS）的 study region／study period、target vs complementary events。polygon 頂點順序差異（ETAS R 反時針、EEPAS 順時針）明講 |
| Embargo 新邊界 | 可用：義大利區域／目錄／實驗設定、Biondini 2023 已發表參數、PyEEPAS 公開 GitHub（https://github.com/phonchi/EEPAS）。**不可用**：審稿中 PyEEPAS／EEPAS_Software 稿件的數字結果、檢驗分數、加速倍數；Taiwan_EEPAS 與 EEPAS_TW 任何結果（`reference/notes/eepas.md` L238–240 一類數字禁用） |
| 資料來源與存放 | 依 `EEPAS/data/README.md`：從 Google Drive 資料夾 `170WNb8M8PQJDX1B2JSQYfqj4ad80TC0U`（已確認可公開讀取）取 `HORUS_Ita_Catalog.txt`、`HORUS_Ita_DataOrigin.txt`、`HORUS_Italy_polygon_filtered.mat`；區域檔用 `EEPAS/data/*.mat`。**只放 `data/cache/`，不進版控**，提供 `scripts/fetch_italy_data.py` 下載＋轉換＋預算四張預報快取 |
| 複發模型／PSHA | 複發模型（BPT、Weibull、hazard function）移出主文併入附錄 F；PSHA 縮成第 17 章一節 |
| 檔名 | 第一部與附錄採新檔名、新編號（`01_`–`20_` 新 slug），舊 `.py`＋`.ipynb` 刪除；第二部 `01`–`08`、`23` 檔名不動，H1 顯示章號機械式改為 21–29 |
| 章數 | 20 章（codex 批判後由 17 章拆分；使用者同意第二部章號順勢重編） |
| 視覺 | 全部做：web-design-engineer 重做主題 CSS；diagram-design（概念圖）／archify（流程圖）產出**可互動**（hover 說明、逐步揭露）SVG 嵌入頁面；Plotly 圖依 dataviz 規範與示意圖共用同一套 token（以 ACCENT 藍／QUAKE 紅為基礎） |
| 驗收與發布 | 離線 build 無警告 → Playwright 桌機＋手機截圖 → 「只修過基礎統計、不懂地震」人設的讀者 subagent 通讀 → codex 交叉審查 → embargo grep 終檢 → 數學重算；**通過後 commit main 並 ghp-import 發布 gh-pages**（使用者已授權） |

### 我自行採用的預設（核准時如不同意請說）
- `pyproj` 加入允許套件（venv 已有 3.7.2，用於 CELLE 反投影與 RDN2008↔WGS84）。
- 共用模型與檢驗程式放進 `gdms_toolkit`，章內只顯示簡化片段（仍摺疊）。
- 互動示意圖用 `IPython.display.HTML` 輸出（myst-nb 保證渲染），不用 `{raw}` 指令。

## 新第一部章序（20 章；設計子代理 17 章版經 codex 唯讀批判後拆分）

全書主線句：一份義大利 2012–2021 的擬前瞻預報實驗，從讀目錄、訂規格、做四張預報（SUP、PPE、ETAS、EEPAS）到用測試期真實目標地震檢驗。可驗收的承諾改為「能讀常見的網格率預報論文，辨認其假設、證據與需要補讀的部分」。每章 2,000–3,000 字、圖 ≤7 張、每個 cell <30 秒；每章末固定「本章填入的規格欄位」段，未填欄寫「待填（第 N 章）」。每個模型章末放同一格式的「預報圖＋目標事件疊圖（不評分）」。第二部 H1 顯示章號機械式改為 21–29（檔名不動）。

| # | 檔名 | 標題（面向學生） | 新定義概念（統計 ｜ 地震） | 資料 | 規格卡 | 搬用舊資產 |
|---|---|---|---|---|---|---|
| 1 | `01_forecast_question.py` | 地震可以「預報」嗎？一份預報長什麼樣 | 期望數 vs 機率（口語）、基準模型 ｜ 規模／震央／深度、餘震序列、CSEP 全名、前瞻／擬前瞻、OEF 一句話 | HORUS M≥5 地圖；第 15 章 EEPAS 預報圖當終點預覽（圖說當場解釋顏色＝格內期望數、點＝之後才發生的目標地震；EEPAS 只是名字）；合成三種時間線 | 只出現空白骨架 | `09` 間隔／每日計數雙圖改用 2016 Amatrice；`foundation_randomness` 時間線 |
| 2 | `02_reading_horus.py` | 讀一份地震目錄：HORUS 與義大利 | 截斷樣本、直方圖／經驗分布 ｜ Mw 均一化、規模的十分位取整（2.45／4.95 是箱邊界）、深度截斷與固定深度假值、經緯度 vs 投影公里、宣告區域 polygon、S ⊃ R 與邊界效應、CPTI15 vs 儀器目錄、**漏報與完整度一句話定義**（估計留第 8 章）、1960 起早期完整門檻較高 | HORUS 全目錄、CELLE、CPTI15 | **章末立卡**：目錄來源／規模尺度、深度、S、R（註明 177 格為既有研究區重演，選區含 2021 年前資訊）、目錄時間範圍 | `foundation_catalog` 三種讀法圖；`23` 年計數雙圖結構 |
| 3 | `03_experiment_spec.py` | 標準預報實驗規格：把題目一次說清楚 | 學習／測試分離、資料窺探（提前解釋檢驗不是洩漏，依結果改模型才是）、**可追蹤的發報時間軸**（步長、窗是否重疊、資料截止、測試期已發生事件可否進入下一窗歷史、參數凍結≠歷史凍結）、格與箱＝題目解析度 ｜ warm-up／learning／testing、target vs complementary（只教集合位置，估計角色留第 13 章）、Jalilian 對照表、m0／mT／Mc 三門檻專節（角色）、十欄格式（欄位意義，數值意義待第 4 章）、擬前瞻 vs 真前瞻 | 學習期＋測試期子集；目標事件表（由載入器數出） | 填三期、窗定義、m0、mT、規模箱、輸出格式 | 新寫；`17` 12.1 概念 |
| 4 | `04_poisson_and_sup.py` | 統計工具箱 I：計數、Poisson 過程，與第一張預報 SUP | Poisson 分布、齊次／非齊次 Poisson 過程、率 λ 與 Λ＝∫λ（密度→積分→格內期望數→計數分布的橋接）、P(N≥1)＝1−e^{−Λ}、過度離散、蒙地卡羅與種子、經驗分位數 ｜ SUP（Spatially Uniform Poisson）當 Poisson 的具體成果、GR 分配到規模箱先用給定 b | HORUS M≥5 年計數；2012 Q1 窗 | SUP 卡（第一次整卡填滿） | `foundation_randomness` 3.1–3.2；`10` 非齊次示意 |
| 5 | `05_likelihood_estimation.py` | 統計工具箱 II：概似與估計 | 概似／對數概似（Poisson 率；固定事件數下的時間密度概似，供第 10 章）、MLE 與不確定性、抽樣分布 vs 標準誤 vs 預測不確定性、bootstrap ｜ 學習期 M≥5 年率 MLE | HORUS M≥5 年計數 | 不填 | `foundation_inference` 7.2–7.3 |
| 6 | `06_simulation_tests_scores.py` | 統計工具箱 III：模擬、檢定與分數 | 檢定重述（統計量→模擬分布→分位數→p 值）、虛無假設與「未拒絕≠證實」、對數分數、概似比、資訊增益一句話 ｜ 測試期計數是否與學習期率相容（未命名的 N-test）。AIC／貝氏／KL 進附錄 A | HORUS M≥5 年計數 | 填「評分方式」骨架 | `foundation_inference` 7.4–7.5 |
| 7 | `07_gr_bvalue.py` | 大小地震的比例：GR 律與 b 值 | 指數分布 MLE（Aki）、分箱修正、偏差、bootstrap 區間 ｜ GR、a／b、β＝b ln10、b 隨時空變化一圖 | HORUS（S 內、<40 km） | 填 b（自估 vs Biondini 1.084） | `11`：`b_exact`、`b_trim` |
| 8 | `08_completeness.py` | 從哪個規模開始相信目錄：完整度 Mc | 估計偏差與選擇效應、移動窗 ｜ Mc 一法詳講（MaxC＋修正）、其餘（b 穩定、KS）比較表、Mc 時空變化、早期餘震不完整、warm-up 期不完整對 EEPAS 長時間核的影響 | HORUS | 填 Mc、確認 m0＝2.45 ≥ Mc（用箱邊界） | `11`：`maxc`、`mc_bstability`、`mc_ks` |
| 9 | `09_ppe_forecast.py` | 鄰近過去地震：PPE 預報 | 核與帶寬、空間密度正規化、177 格內積分近似 ｜ PPE（Proximity to Past Earthquakes，Biondini 參數）、quasi time-dependent（本站用凍結學習期的固定版並標明）、預報圖固定樣式、十欄匯出（177→0.1° 是重新分格：守總期望數、交集面積、遮罩） | 學習期 m≥2.45；2012 Q1 窗 | PPE 卡 | `16` 11.6 PPE 核平滑（改 177 格） |
| 10 | `10_clustering_laws.py` | 地震會互相引發：Omori、Båth 與產能 | 冪律衰減與正規化、有限窗 MLE、模型 vs 不完整混淆 ｜ Omori–Utsu、Båth、產能 α、主震／餘震／群震、GK 除叢一段 | L'Aquila 2009、Umbria–Marche 1997（學習期） | 不填；埋負二項伏筆 | `12`：`omori_rate`、`omori_integral`、`omori_mle`、`gk_window` |
| 11 | `11_conditional_intensity.py` | 讓率隨歷史更新：條件強度與點過程概似 | λ*、H_t、補償子、Σln λ*−∫λ*（本章擁有）、自激發 vs 非齊次 Poisson、「給定發報當下歷史的未來分布」vs「沿實際歷史的條件概似」 ｜ 把 Omori 核掛到每顆地震上。稀疏化模擬、殘差／時間變換進附錄 A | 合成為主；L'Aquila 補償子檢查一張 | 不填 | `10`：`simulate_branching`、`compensator` |
| 12 | `12_etas_structure.py` | ETAS：把叢集寫成分支過程 | 分支過程、世代、分支比與臨界 ｜ ETAS 五組成與全式、simplETAS「七釘二估」、α＝β 需 Mmax | 合成模擬 | 預告 ETAS 卡差異欄 | `13`：`A_from_n`、`simulate_etas`、`simulate_family` |
| 13 | `13_etas_italy_forecast.py` | ETAS 在義大利：估兩個參數、出一張三個月預報 | 固定形狀的低維 MLE、不可辨識性、target／complementary 在概似中的角色 ｜ 本站近似＝背景＋已知歷史的第一代期望率（**明標為近似，非 ETAS 家族必然**）、200 條模擬目錄對照、R–J 特例與隨機除叢摺疊、背景 μ(x,y) 用學習期平滑地震度取代 MPS19 | 學習期 m≥3.95（估）；m≥2.45 歷史（預報） | **ETAS 卡**，標與主線不同欄 | `14`：`fit_mu_A_alpha`（改估 ν、A）、`omori_matrices`、`rj_prob` |
| 14 | `14_psi_precursory_scale.py` | 大地震之前，小地震先變多變大？Ψ 現象 | C(t)、三迴歸與回歸稀釋、固定搜尋規則 vs 事後選擇、時空抵換 ｜ Ψ、M_P／T_P／A_P／M_m、回溯示意 | 合成為主；L'Aquila 2009 回溯圖標「回溯示意」 | 不填 | `15`：`cumag`、`scan_max_z`、三迴歸圖 |
| 15 | `15_eepas_italy_forecast.py` | EEPAS：每顆地震都是它尺度上的前兆 | 對數轉換與對數常態、聯合／條件密度、可分離乘積核與正規化、混合權重 μ_E ｜ **用一顆輸入事件示範三核如何形成貢獻，再疊加＋背景**；η(m)、Δ(m) 只講用途（推導進附錄 D）；等權（除叢權重摺疊）；三階段擬合觀念；前置時間 | 1960 起 m≥2.45（warm-up，標明早期不完整）；同一窗 | **EEPAS 卡**（Biondini 參數表） | `16`：`g_mag`、`f_time`、`sigma_space`、`completeness` |
| 16 | `16_test_number.py` | 檢驗 I：數量對不對？N-test | 離散尾端機率、雙尾 δ₁／δ₂、負二項與過度離散、多重檢定與 Bonferroni（在此定義）｜ **預先指定**逐窗累積圖與序列貢獻分解（同一序列多顆目標不是獨立證據）、Poisson 版先跑再換負二項（不可看結果才換）、「通過」≠「率抓對」 | 四張預報 × 40 窗（快取）；測試期目標事件 | 填「評分方式：一致性（數量）」 | `17`：負二項圖 |
| 17 | `17_test_space_magnitude.py` | 檢驗 II：位置與規模對不對？S、M、cL、L | multinomial、邊際化與條件化、POLL／jPOLL、條件多項模擬、L vs cL、BILL 一段 ｜ 哪些量被固定、哪些被模擬；等面積 SUP 的條件 S-test 為何退化 | 同 16 | 填「評分方式：一致性（位置／規模）」，整卡完成 | `17`：`jpoll`；`18`：`joint_poll` |
| 18 | `18_test_comparison.py` | 檢驗 III：哪張預報比較有資訊？ | 配對、共變異與樣本相依、配對 t 區間（T-test 配對結構）、IGPE／IGPEc、功效模擬、參數不確定性 vs 地震隨機性 ｜ 基準選擇改變結論、勝過 SUP≠實用。Molchan、可靠度、Brier 進附錄 E | 同 16 | 填「評分方式：比較」與基準 | `18`：`igpe`、`power` |
| 19 | `19_ensembles.py` | 把模型加起來：凸組合 | 凸組合守恆的是**期望數**（需正規化條件；概似不守恆）、凹性、權重不可用最終測試資料選 ｜ ETAS＋EEPAS 時間尺度互補、Bayona 2022 前瞻退步教訓。乘法 hybrid、logistic 權重進附錄 E | 第 13、15 章預報做凸組合 | 混合卡（多「成分與權重」列） | `19`：`loglik`、`fit_beta`、`to_weights` |
| 20 | `20_beyond_forecast.py` | 預報之後：讀一篇論文、進入決策，與回到儀器 | **有旁註的原文閱讀示範**（Biondini 2023 一段方法、一張結果圖、一句結論：找假設、證據、限制）、事件瞬時發生率 vs 地動危害、survival／hazard 直覺（距上次事件多久會改變條件機率，Poisson 不是共同前提）、回歸期 ｜ OEF-Italy／USGS OAF、GMPE 一句話、PSHA 濃縮一節（五步驟進附錄 F）、預報溝通 | 合成危害曲線一張；OAF 示意一張 | 無；章末「義大利卡換台灣卡要改哪些欄」（全標待定）→ 前指 `01_overview` | `21`：`hazard`；`22`：`p_year`；`20` 的 hazard function 直覺段 |

### Codex 唯讀批判已併入的修正（job `795cb536-203f-4097-a2e0-9eb8923aca75`）
拆章（第 5／6／14–15）；SUP 併入第 4 章；第 2 章一句話定義完整度；第 3 章可追蹤發報時間軸；規模箱取整與邊界；177 格選區含事後資訊須明說；warm-up 不完整；PPE quasi time-dependent；177→0.1° 是重新分格；ETAS 近似明標；2016 序列預先指定分析；混合模型守恆期望數而非概似；第 20 章加原文閱讀示範與 survival 直覺；否決倒敘案的理由改為「依測試結果改模型才是問題」。

### 未定義概念認領表（驗收用，讀者 subagent 對照）

| 概念 | 章 | 概念 | 章 |
|---|---|---|---|
| CSEP 全名、前瞻／擬前瞻 | 1, 3 | 分位數↔p 值橋接、虛無假設、對數分數 | 6 |
| 宣告區域 polygon、S ⊃ R、邊界效應、完整度一句話、規模箱邊界 | 2 | 概似、MLE、bootstrap、標準誤 vs 預測不確定性 | 5 |
| depth cutoff、規模尺度均一化 | 2 | GR、b 值 | 7 |
| target／learning／testing、Jalilian 對照、發報時間軸 | 3 | Mc 估計、偏差 | 8 |
| m0／mT／Mc 三門檻 | 3（角色）、8（數值） | 核 kernel、PPE 全名、quasi time-dependent、重新分格 | 9 |
| 期望數 vs 機率、密度→積分→期望數、過度離散、蒙地卡羅、SUP 全名 | 4 | 條件強度、補償子、自激發 | 11 |
| Omori、Båth、產能、除叢 | 10 | 分支比、ETAS 五組成 | 12 |
| Ψ、C(t)、時空抵換 | 14 | EEPAS 三核、η／Δ 用途、μ_E | 15 |
| 負二項 N-test、δ₁／δ₂、Bonferroni（定義與使用） | 16 | POLL／jPOLL／BILL、L vs cL、M-test、multinomial | 17 |
| T-test 配對結構、IGPE、功效 | 18 | 凸組合期望數守恆 | 19 |
| survival／hazard 直覺、回歸期、OEF／OAF | 20 | AIC／KL／貝氏／稀疏化／殘差／Molchan／Brier／乘法 hybrid | 附錄 A／E |

### 附錄對應

| 新附錄 | 內容 | 來源 |
|---|---|---|
| A 記號、統計工具與點過程 | 讀者版記號表、兩套術語對照、工具箱推導、AIC／KL／貝氏、概似／時間變換／殘差／稀疏化 | a ＋ `foundation_inference` 推導 |
| B 目錄、完整度與叢集律 | 新增 B.0「HORUS 取得、快取與 CELLE 反投影」 | b |
| C ETAS | 加 simplETAS 釘參理由、α＝β 的 Mmax 條件 | c |
| D Ψ 與 EEPAS | 三階段擬合、Δ／η 推導 | d |
| E 檢驗與組合 | 加十欄格式、177→0.1° 重新分格規則、教學版檢驗程式說明、Molchan／可靠度／Brier、乘法 hybrid 與 logistic 權重 | e ＋ 舊 18／19 部分 |
| F 危害與複發 | 吸收舊 `20_recurrence_models` 文字與公式、PSHA 五步驟（附錄為 .md，圖改為預先輸出 PNG 或放棄） | f ＋ 20 ＋ 21 |
| G 觀測 | 不動 | g |

（備選「先檢驗後模型」倒敘式已評估並否決：模型章邊做邊看測試期分數違反擬前瞻紀律，且只有 SUP 時 S-test 無意義。）

## 第二部銜接與台灣內容併入

- `23_taiwan_outlook.py`（新顯示 29 章）吸收：舊 `11`／`12` 的台灣 Mc／b 值時空變化、Omori／Båth；舊 `09` 的 2024 花蓮序列間隔／每日計數圖；舊 `16` 的台灣 PPE 示意（改寫為「沿用第 9 章義大利 PPE 核形式，參數為台灣示意」）。結構：「用第 3 章規格卡逐欄改成台灣版（目錄 CWA、ML、深度、R、期間皆待定）→ 台灣目錄的統計面貌 → 大埔／花蓮案例 → 展望」。**Taiwan EEPAS 仍只寫「在地化工作進行中」**。此章將超過 15–20 分鐘常態，使用者已知悉接受。
- 第二部九個檔的 H1 章號 18–26 → 21–29（含 `08_explore_ideas.md`）；`01_overview.py` 首段一句改為「前一部用義大利實驗建立了…」；第二部所有 `{doc}` 與「第 N 章」字樣指向新章。
- `00_intro.md`：入口改 `01_forecast_question`，附錄表與導讀更新。

## 資料與工具層（Phase 0，必須先做並通過閘門）

1. `scripts/fetch_italy_data.py`：從 Drive 下載三個檔到 `data/cache/italy/`（需 `gdown` 或直接檔案 ID；記錄下載日期與 `HORUS_Ita_DataOrigin.txt` 的授權說明）；轉 CSV；**預先計算四張預報 × 40 窗存 `data/cache/italy_forecasts/`**（章內只載入，確保每 cell <30 秒）。
2. `gdms_toolkit/italy.py`：`load_horus()`、`testing_cells()`（177 格反投影為經緯度，用 CPTI15 十五點 lon/lat↔x_km/y_km 往返驗證，判定 EPSG 6875 vs 7794）、`collection_polygon()`、`SPEC` 字典、`target_events()`、`grid_01deg_in_R()`、`write_csep_10col()`、`spec_card()`（輸出 HTML 規格卡）。
3. `gdms_toolkit/italy_models.py`（SUP、PPE、ETAS-lite、EEPAS 率函式，事件×格心權重矩陣實作；`get_forecast(model, window)` 快取）、`gdms_toolkit/csep_teaching.py`（`n_test`、`s_test`、`m_test`、`cl_test`、`l_test`、`t_test_igpe`、負二項版）、`gdms_toolkit/viz.py` 加 `plot_forecast_map()` 統一樣式與新 token。
4. **閘門（任何一項不符即回頭改規格卡與第 3 章）**：原始 HORUS 起始年是否到 1960（否則 warm-up 改寫）；CELLE 反投影誤差；所有 Biondini 參數（EEPAS a_M…μ_E、PPE a／d、b＝1.084、m0／mT／delay／p／c）逐一對照 `reference/[2023] Italy_EEPAS.pdf` 原文而非筆記或 JSON；目標事件實際顆數（「27」只是預期）；學習期 m≥2.45／m≥3.95 事件數與各模型計算時間。
5. 重寫 `book/_conventions.md`：記號表、公式所有權改新章號、兩套術語對照、規格卡模板（MyST admonition `:class: spec-card`）、§6 紅線改為新 embargo 邊界（保留 STEP／AIC 權重等未核實項警告）；新 `book/_toc.yml`；`reference/notes/rewrite_20260911/chapter_map.json` 標記已被取代並新增新映射。

## 視覺層

- **一頁 spike 先行**：在測試頁確認 `display(HTML())` 內的 inline `<svg>`＋`<script>` 在 build 後可執行、`_ext/teaching_render.py` 的 inline-script 去重不會吃掉它、390 px 可用；再擴到各章。
- 主題（web-design-engineer，Redesign · Preserve 模式）：重寫 `book/_static/teaching.css`（字型與行長、定義框、`spec-card`、術語卡、章首導覽、深色模式），必要時 `book/_templates/` 覆寫；保留 sphinx-book-theme 骨架。
- 示意圖（diagram-design 概念圖、archify 流程圖）：先以站台 token（ACCENT `#2a78d6`／QUAKE `#e34948`／`Blues`）完成 diagram-design 的 style-guide onboarding 並存 profile；輸出到 `book/_static/diagrams/`。清單（每章 1–2 張）：預報實驗流程（archify workflow）、S ⊃ R 兩區域、發報時間軸（三段期＋滾動窗＋資料截止）、規格卡骨架、m0／mT／Mc 三門檻、密度→積分→格內期望數→計數分布橋接、模型家族層級、ETAS 分支樹、EEPAS 一顆事件的三核貢獻、CSEP 檢驗流程（哪些量固定、哪些模擬）、凸組合。互動限 hover 說明與逐步揭露，無外部依賴。
- 圖表：`gdms_toolkit.viz` 的 PALETTE／layout 與示意圖 token 對齊；寫任何圖前依 dataviz 規範。

## 執行順序

- Phase 0（單一 agent）：資料與工具層 1–5、視覺 spike、token／profile、規格卡模板。
- Phase 1（四個 subagent 並行，各組不寫同一檔）：A 組 1–6 章｜B 組 7–11 章（需 PPE 函式）｜C 組 12–15 章（需 ETAS／EEPAS 函式）｜D 組 16–20 章（需 `get_forecast` 與檢驗函式）。每章寫完自跑 nbclient 執行與計時。
- Phase 2（單一 agent）：附錄 A–F、`00_intro.md`、第二部銜接與第 23 章併入、每頁「參考資料與延伸閱讀」（沿用 `reference/notes/reading_*_verified.md` 區塊重新對應新章，並補義大利來源：Biondini 2023、Jalilian 2019、Savran 2022／Graham 2024、Mancini & Marzocchi 2023、HORUS／CPTI15 出處）、全站 `{doc}` 與章號校正、刪舊檔、README 更新。
- Codex 觸點（`claude-codex-collab`，記錄精確 job ID）：(a) 大綱唯讀批判已完成（job `795cb536-203f-4097-a2e0-9eb8923aca75`，log `~/.local/share/claude-codex-collab/runs/795cb536-…`）；(b) Phase 1 後以 `--resume-job` 續問同一 session 做讀者視角審查；(c) 最終 embargo／數學稽核。

## 驗證

1. `python scripts/sync_teaching_notebooks.py` → `python scripts/build_teaching_offline.py`（`-W --keep-going` 無警告）；`.py` 與 `.ipynb` 都進版控。
2. `python scripts/check_teaching_browser.py`（桌機 1440×900、手機 390×844），含示意圖互動檢查與截圖 montage。
3. 讀者 subagent（人設：只修過基礎統計、完全不懂地震）逐章回報「哪個詞沒定義就用」「哪句讀不懂」，對照認領表；codex 交叉審查。
4. `grep -rin "PyEEPAS\|Taiwan_EEPAS\|EEPAS_TW" book/` 只允許 GitHub 連結命中；稿件數字零命中。
5. `scripts/verify_teaching_math.py` 擴充到新章公式（Omori 積分、EEPAS 核正規化、jPOLL、IGPE）。
6. 通過後 commit main、`ghp-import` 發布 gh-pages，線上全部頁面 HTTP 200 核對。
7. 完成後更新記憶：`embargo-under-review-papers`（新邊界）、`teaching-site-style`（工具箱內隱藏實作取代「不手刻完整模型實作」）。
