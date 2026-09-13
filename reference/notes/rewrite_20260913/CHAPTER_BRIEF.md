# 第一部各章撰寫簡報（給撰寫子代理；必讀 `book/_conventions.md` 與 `TOOLKIT_API.md`）

> 此文件保留啟動時的撰寫簡報。現況修正與驗證以 `continuation/implementation_notes.md`、`continuation/math_checks.json`、正文及 `_conventions.md` 為準：PPE 也延遲 50 天；學習末窗已裁齊；25／27 顆差異未逐事件證實；檢驗數字使用實算。

## 共同規則（摘要，細節在 `_conventions.md`）
- 讀者只修過基礎統計、不懂地震；**每個新名詞第一次出現時，同一段就定義並附英文**；縮寫先展開。
- 每段先給一個能抓住的東西（數字、欄位、圖、公式），再講判斷。句子 ≤ 30 字，主詞明確，一段一件事。
- 每章 2,000–3,000 字、≤ 7 張圖、每 cell < 30 秒；亂數固定 seed；模擬 n_sim ≤ 3,000，先估記憶體。
- 檔案：`book/NN_slug.py`（jupytext py:percent，無 front-matter 以外的 header 沿用舊檔樣板），章首固定 `setup_plotly()` cell；程式 cell 一律 `# %% tags=["remove-input"]`；每圖以裸 `fig` 結尾。
- 章末固定兩段：「本章填入的規格欄位」（用 `italy.spec_card(rows)` 輸出，未填欄給 `None`）與前指下一章的 `{doc}`；最後是「參考資料與延伸閱讀」（沿用 `reference/notes/reading_*_verified.md` 可對應的區塊，並補義大利來源）。
- 資料只經 `gdms_toolkit.italy`／`italy_models.get_forecast()`（讀快取）／`csep_teaching`；預報圖用 `viz.plot_forecast_map()`。**不在章內重算四張預報**。
- 模型章末放「預報圖＋測試期目標地震疊圖（圖說明寫：只疊圖不評分，評分在第 16–18 章）」。
- Embargo：不引用審稿中 PyEEPAS／EEPAS_Software 稿件的數字；台灣 EEPAS 不提；PyEEPAS 只給 GitHub 連結。
- 互動示意圖：從 `book/_static/diagrams/` 讀現成 HTML 以 `display(HTML(...))` 嵌入（Phase 0 會提供清單與檔名）；沒有對應檔就不要自己畫 SVG，改用 Plotly 或文字。
- 寫完自跑：`ulimit -v 4000000 && .venv/bin/python3 scripts/execute_teaching_pages.py <file>`（或 nbclient），確認無錯、每 cell 時間，並用 `grep -c "^# " ` 粗估字數。

## 已知資料事實（可直接寫進正文）
- HORUS 2024-06 版：493,418 筆，1960-01-03 起；欄位 Year…Se、Lat、Lon、Depth、Mw、sigMw、Geo-Ita、Geo-CPTI15、Ev.type、ISIDe。
- 實驗目錄（S 內、深度 ≤ 40、排除非地震、≤ 2021）：435,237 筆；R 內 415,412 筆。
- 目標地震（R 內、取整 Mw ≥ 5.0）：學習期 27 顆（論文亦 27）；測試期 25 顆（論文 27；目錄版本差異）。測試期有事件的窗 12 個。
- 完整度（Lolli 2020，Biondini 引用）：M≥4.0 自 1960、≥3.0 自 1981、≥2.5 自 1990、≥2.1 自 2003、≥1.8 自 2005。
- 四張預報期望總數（本站教學版）：學習期 SUP 26.95、PPE 26.95、EEPAS 27.82、ETAS 13.14（論文表 5：27、27、27.7、27.5；ETAS 差異來自只算第一代）；測試期 SUP 12.25、PPE 14.01、EEPAS 14.44、ETAS 7.50（論文表 7 三個月窗：12.27、14.03、14.41、8.36）。
- 檢驗（測試期 25 顆）：Poisson N-test 四模型皆低估（δ1 ≤ 0.007）；負二項（σ²=67.76）四模型皆過；IGPE 對 SUP：PPE 0.40、EEPAS 0.70、ETAS 1.05；S-test SUP 分位數 0.000。
- 參數（Biondini 2023 表 3，M+A）：b=1.084；PPE a=0.62、d=30 km、s=9e-13；EEPAS-NW aM=1.22、bM=1、σM=0.25、aT=2.55、bT=0.35、σT=0.150、bA=0.52、σA=1.00、μ=0.18、delay 50 天；ETAS-SUP K=0.029、c=0.004、p=1.042、D=1.04、γ=0.45、α=1.12、q=1.5、ν=0.264。
- 座標：EPSG:7794（RDN2008 Italy zone E-N）公里；`CELLE_ter.mat` 為投影公里，README 寫成經緯度是錯的。polygon 頂點順序：ETAS R 反時針、EEPAS 順時針。

## 各章規格
（欄位：檔名｜標題｜學習成果｜新定義概念｜資料與圖｜規格卡｜可搬舊資產）

1. `01_forecast_question.py`｜地震可以「預報」嗎？一份預報長什麼樣｜能分辨預測／預報／預警，看懂一張網格預報圖｜期望數 vs 機率（口語）、基準模型；規模／震央／深度、餘震序列、CSEP 全名、前瞻／擬前瞻、OEF 一句話｜HORUS M≥5 地圖；EEPAS 十年合計預報圖＋目標事件（終點預覽，圖說當場解釋顏色與點；EEPAS 只是名字）；2016 Amatrice 序列間隔／每日計數雙圖；合成三種時間線（`gdms_toolkit.teaching.learning_catalog`）｜只出現空白骨架卡｜舊 `09` 雙圖、`foundation_randomness` 時間線。
2. `02_reading_horus.py`｜讀一份地震目錄：HORUS 與義大利｜能逐欄讀 HORUS，說出為何統一規模、限制深度、分 S 與 R｜截斷樣本、直方圖；Mw 均一化、規模取整（2.45／4.95 是箱邊界）、深度截斷與固定深度假值、經緯度 vs 投影公里、宣告區域 polygon、S ⊃ R 與邊界效應、CPTI15 vs 儀器目錄、**漏報與完整度一句話**、1960 起早期完整門檻較高｜欄位表（程式輸出）、年計數 vs 最小規模、深度直方圖、S／R 地圖（`plot_forecast_map` 全 0 或 `Scattergeo`）、規模–時間圖、旗標 vs 自算 polygon 對照｜章末立卡：目錄、深度、S、R（註明選區含 2021 前資訊）、時間範圍｜舊 `foundation_catalog`、`23` 年計數。
3. `03_experiment_spec.py`｜標準預報實驗規格：把題目一次說清楚｜能把論文方法段逐項對回規格卡、兩套術語互譯｜學習／測試分離與資料窺探（檢驗不是洩漏，依結果改模型才是）、**發報時間軸**（40 窗、91.31 天、發報時刻、資料截止、delay、參數凍結≠歷史凍結）、格與箱＝解析度；三期、target vs complementary（只教集合位置）、Jalilian 對照表、m0／mT／Mc 三門檻角色、十欄格式欄位意義、擬前瞻 vs 真前瞻｜時間軸圖（互動示意圖 `timeline`）、目標事件地圖與表（25 顆，含 2016 序列標示）、十欄格式前 10 列（`write_csep_10col` 於 SUP 十年合計）｜填三期、窗、m0、mT、規模箱、輸出格式｜新寫。
4. `04_poisson_and_sup.py`｜統計工具箱 I：計數、Poisson 過程，與第一張預報 SUP｜能分辨 Poisson 分布／過程／率，親手做出 SUP｜Poisson 分布與變異數＝均值、齊次／非齊次過程、率 λ 與期望數 Λ=∫λ（密度→積分→格內期望數→計數分布橋接）、P(N≥1)、過度離散、蒙地卡羅與 seed、經驗分位數；SUP 全名、GR 分配到規模箱（b 先給定 1.084，第 7 章再估）｜R 內 M≥5 年計數 1960–2021 條圖＋Poisson 理論；模擬 1,000 次年計數；SUP 一窗預報圖＋十年合計圖｜SUP 卡（第一次整卡填滿）｜舊 `foundation_randomness` 3.1–3.2、`10` 非齊次示意。
5. `05_likelihood_estimation.py`｜統計工具箱 II：概似與估計｜能讀懂 MLE、標準誤、bootstrap 區間｜概似／對數概似（Poisson 率；固定事件數下的時間密度概似，供第 10 章）、MLE、抽樣分布 vs 標準誤 vs 預測不確定性、bootstrap｜學習期 M≥5 年率概似曲線、bootstrap 直方圖｜不填｜舊 `foundation_inference` 7.2–7.3。
6. `06_simulation_tests_scores.py`｜統計工具箱 III：模擬、檢定與分數｜能把「分位數分數」接回 p 值，讀懂對數分數｜檢定重述（統計量→模擬分布→分位數→p 值）、虛無假設與「未拒絕≠證實」、對數分數、概似比、資訊增益一句話（AIC／貝氏／KL 進附錄 A）｜測試期計數 25 vs 學習期率的模擬分布（未命名的 N-test）、對數分數比較兩個率｜填評分骨架｜舊 `foundation_inference` 7.4–7.5。
7. `07_gr_bvalue.py`｜大小地震的比例：GR 律與 b 值｜能估 b 與區間並解讀｜指數分布 MLE（Aki）、分箱修正、偏差、bootstrap 區間；GR、a／b、β=b ln10｜HORUS 學習期 R 內 M≥2.5 規模–頻率圖、b 隨時間一圖｜填 b（自估 vs 1.084）｜舊 `11`：`b_exact`、`b_trim`。
8. `08_completeness.py`｜從哪個規模開始相信目錄：完整度 Mc｜能估 Mc 並確認 m0 ≥ Mc｜估計偏差與選擇效應、移動窗；Mc 一法詳講（MaxC＋修正）、其餘（b 穩定、KS）比較表、Mc 時空變化、早期餘震不完整、warm-up 不完整對 EEPAS 長時間核的影響｜Mc 隨年、Mc 地圖（177 格）、L'Aquila 序列早期不完整｜填 Mc、確認 m0＝2.45 ≥ Mc｜舊 `11`：`maxc`、`mc_bstability`、`mc_ks`。
9. `09_ppe_forecast.py`｜鄰近過去地震：PPE 預報｜能親手做出 PPE 並匯出十欄格式｜核與帶寬、空間密度正規化、格內積分近似；PPE 全名與參數、quasi time-dependent（本站凍結版）、預報圖樣式、177→0.1° 是重新分格｜128 顆 M≥5 來源事件地圖、核形狀圖、PPE 一窗預報圖＋十年合計圖、十欄輸出對照（`write_csep_10col_regridded`）｜PPE 卡｜舊 `16` 11.6。
10. `10_clustering_laws.py`｜地震會互相引發：Omori、Båth 與產能｜能從 L'Aquila 2009 估 Omori 參數｜冪律衰減與正規化、有限窗 MLE、模型 vs 不完整混淆；Omori–Utsu、Båth、產能 α、主震／餘震／群震、GK 除叢一段｜L'Aquila 2009 與 Umbria–Marche 1997 序列率衰減、Omori 擬合、產能 vs 規模｜不填｜舊 `12`。
11. `11_conditional_intensity.py`｜讓率隨歷史更新：條件強度與點過程概似｜能讀 λ*(t|H_t) 與概似兩項｜λ*、H_t、補償子、Σln λ*−∫λ*、自激發 vs 非齊次 Poisson、「發報當下的未來分布」vs「沿實際歷史的條件概似」（稀疏化、殘差進附錄 A）｜合成自激發示意、L'Aquila 補償子檢查｜不填｜舊 `10`：`simulate_branching`、`compensator`。
12. `12_etas_structure.py`｜ETAS：把叢集寫成分支過程｜能說出 ETAS 五組成與分支比｜分支過程、世代、分支比與臨界；ETAS 全式、simplETAS「七釘二估」、α=β 需 Mmax｜合成模擬世代分解、參數效應｜預告 ETAS 卡差異欄｜舊 `13`。
13. `13_etas_italy_forecast.py`｜ETAS 在義大利：估兩個參數、出一張三個月預報｜能解讀 ETAS 三個月預報與其近似｜固定形狀的低維 MLE、不可辨識性、target／complementary 在概似中的角色；本站近似＝背景＋第一代（明標）、200 條模擬目錄對照、R–J 與隨機除叢摺疊｜學習期 M≥4.0 事件 MLE 估 ν、K 示意（`cKDTree` 限距，<30 秒）、ETAS 2016 Q3 預報圖＋十年合計、模擬對照｜ETAS 卡（標差異欄）｜舊 `14`。
14. `14_psi_precursory_scale.py`｜大地震之前，小地震先變多變大？Ψ 現象｜能解釋 Ψ 三尺度關係與回溯 vs 前瞻的落差｜C(t)、三迴歸與回歸稀釋、固定搜尋規則 vs 事後選擇、時空抵換；Ψ、M_P／T_P／A_P／M_m｜合成示意為主；L'Aquila 2009 回溯圖標「回溯示意」｜不填｜舊 `15`。
15. `15_eepas_italy_forecast.py`｜EEPAS：每顆地震都是它尺度上的前兆｜能讀三核公式並解讀 EEPAS 預報｜對數轉換與對數常態、聯合／條件密度、可分離乘積核、μ_E；**一顆輸入事件的三核貢獻示範**，再疊加＋背景；η／Δ 用途；等權；三階段擬合觀念；前置時間｜單一事件三核圖（互動示意圖 `eepas-kernels`）、EEPAS vs PPE 同一窗、十年合計＋目標事件｜EEPAS 卡｜舊 `16`：`g_mag`、`f_time`、`sigma_space`、`completeness`。
16. `16_test_number.py`｜檢驗 I：數量對不對？N-test｜能算 δ1／δ2 與負二項版｜離散尾端機率、雙尾、負二項與過度離散、Bonferroni（定義）；預先指定逐窗累積圖與序列分解、Poisson 先跑再換負二項、「通過」≠「率抓對」｜四模型逐窗累積期望 vs 觀測、N-test 區間圖（Poisson／NB）｜填一致性（數量）｜舊 `17` 負二項圖。
17. `17_test_space_magnitude.py`｜檢驗 II：位置與規模對不對？S、M、cL、L｜能解釋哪些量固定、哪些模擬｜multinomial、邊際化與條件化、POLL／jPOLL／BILL、條件多項模擬、L vs cL、等面積 SUP 的 S-test 退化｜四模型的 S／M／cL 模擬分布與觀測值｜填一致性（位置／規模），整卡完成｜舊 `17`、`18`。
18. `18_test_comparison.py`｜檢驗 III：哪張預報比較有資訊？｜能算 IGPE 與配對 t 區間並判讀｜配對、共變異與相依、配對 t 區間、IGPE、功效模擬、參數不確定性 vs 地震隨機性；基準選擇改變結論、勝過 SUP≠實用（Molchan、Brier 進附錄 E）｜逐事件對數率比、IGPE 區間圖、換基準（PPE）後的圖、功效模擬｜填比較與基準｜舊 `18`：`igpe`、`power`。
19. `19_ensembles.py`｜把模型加起來：凸組合｜能解釋為何混合常勝過單一模型、前瞻退步｜凸組合守恆的是期望數（需正規化條件）、凹性、權重不可用測試資料選；ETAS＋EEPAS 互補、Bayona 2022 教訓（乘法 hybrid 進附錄 E）｜ETAS＋EEPAS 凸組合掃 π、學習期選權重 vs 測試期表現｜混合卡｜舊 `19`。
20. `20_beyond_forecast.py`｜預報之後：讀一篇論文、進入決策，與回到儀器｜能用旁註讀法讀一段方法、一張圖、一句結論｜有旁註的原文閱讀示範（Biondini 2023 一段方法、一張結果圖、一句結論）、事件瞬時發生率 vs 地動危害、survival／hazard 直覺、回歸期；OEF-Italy／USGS OAF、GMPE 一句話、PSHA 濃縮一節、預報溝通｜合成危害曲線、OAF 示意｜無；章末「義大利卡換台灣卡要改哪些欄」全標待定→前指 `01_overview`｜舊 `21`、`22`、`20`。
