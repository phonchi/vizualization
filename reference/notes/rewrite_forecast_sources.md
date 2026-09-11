# 預報後半部教學改編：來源與驗算紀錄

日期：2026-09-11。範圍：原 15–22 章與附錄 D/E/F。章號由整站重編流程統一處理。

## 改編方式

逐一重寫原 markdown cell，保留原繪圖資料與大部分 code cell，依閱讀順序串起 Ψ → EEPAS → 一致性 → 比較 → 組合 → 複發 → PSHA → 作業 → 台灣觀測。刪除每章重複的參數表、誤解清單、前沿與內嵌附錄固定骨架；將主要代數重新推導成三個獨立附錄。原網址與 eq label 保留，label 僅在主文定義一次。

此改編不是原論文的逐字摘要。刪去未核實、過度概括與跟主線無關的研究數字，仍保留既有圖；圖文已逐圖對照，區分合成、歷史案例與台灣公開目錄示意。

## 直接檢查的來源

- `reference/Psi regression/[2024] psi.pdf`，印刷頁 3465，式 (1)–(2)：cumag 使用 `M_i - (m_c - 0.1)`，不是 `M_i - m_c - 0.1`。本文定義 minimum 為 onset，並未證明它最大化未加權兩段率差。
- `reference/Intro_to_the_theory_of_point_processes.pdf`：版權／書名核對為 Daley–Vere-Jones，Volume I，第二版；目錄確認第 4 章 Renewal Processes，4.1 Basic Properties、4.2 Stationarity and Recurrence Times，7.2 Conditional Intensities, Likelihoods, and Compensators。依其主題補入抽樣起點、剩餘等待與條件強度/補償子區別。
- `reference/[2023] Application of the EEPAS earthquake forecasting model to Italy.pdf`：原文與附錄核對三核、不同規模門檻、權重與測試安排。長代數另獨立重算，不沿用可能損壞的 PDF 文字公式。
- `reference/[2011] Japan.pdf`：首頁 DOI 10.5047/eps.2010.08.002，摘要明載 2000–2009 retrospective fitting；正文確認圖中 0.24、0.42、1.02 資訊差及 0.74、0.61 移植增益，N=1040/396/148。補入主文與附錄 D 延伸閱讀，正文標明不是事後完成的前瞻成績。
- `reference/Ensemble/[2009] Mixture Models for Improved Short-Term Earthquake Forecasting.pdf`：首頁摘要確認 STEP/EEPAS 混合、加州 1984–2004、最適權重是擬合結果，將來需實時測試。主文不再把無顯著差異當成可證的內點條件。
- `reference/ETAS/[2023] Bayes ETAS.pdf`：首頁確認 Naylor、Serafini、Lindgren、Main，2023，DOI 10.3389/fams.2023.1126759；補入模型內參數不確定性與模型間差異的閱讀連結。
- 既有筆記 `psi.md`、`eepas.md`、`ensemble.md`、`testing.md`、`reading_17_23_verified.md`、`psha_step_websearch.md` 作索引；正文只沿用適用的已發表材料。
- 即時開啟確認 pyCSEP `https://docs.cseptesting.org/getting_started/theory.html`、USGS `https://earthquake.usgs.gov/data/oaf/overview.php`、Springer `https://link.springer.com/book/10.1007/b97277`。官方文档用來核對檢驗目標、產品入口與教材版本，不將固定出版案例宣稱為最新完整國家系統現況。

未使用審稿中的 PyEEPAS、Taiwan_EEPAS 或 EEPAS_Software。

## 數學修正

1. cumag 符號與 code 改為 `me - mc + 0.1`。修正後每次合格事件為正跳升。
2. cumag 最小值只最大化帶長度權重的兩段率差：`C(u)=[(u-ts)(tf-u)/(tf-ts)](r_minus-r_plus)`；不能宣稱最大未加權率差、最低前段平均率或真實變點。
3. 回歸方向差異只在明確誤差模型下解釋；兩斜率不能唯一識別真斜率與量測誤差。RMA 不作普遍正確估计。負共變異數也不保證任意第三變數的 R² 必然提高。
4. EEPAS 的時間核正規化、三個摘要與 eta/Delta 有限門檻重算；原對數常態 mean/median 數值錯誤修為 sigma=0.20 時 1.111864、sigma=0.81 時 5.693287。c>0 的 Omori 在零點有限，不發散；核交點不是機制界線。
5. Poisson 网格平均不決定計數分布；每事件 log score 仍有 log(箱體積)，不因除 N 消失。
6. Bonferroni 採實際檢驗數、不要求獨立，不能只憑相關就用更小有效數替代。
7. 完整同質 Poisson 分支群集的 Fano 比 1/(1-n)^2 有模型及完整窗條件；不能泛化為所有有限窗、規模混合產能 ETAS。
8. IGPE 的 N=0 邊界、相依事件區間、AICc 適用條件與有限樣本限制明寫。Molchan AS 的隨機基準為 tau/2，全區時才是 0.5。
9. 組合概似凹性不保證嚴格凹或內點；需事件處至少有不同率及兩端導數條件；未顯著不等於樣本增益精確為零。凸組合總量相同不保證 N-test 通過。
10. 模型率平均與完整 Poisson 分布混合分開；後者通常也不等於負二項。logistic 移除截距/負係數不是定理。
11. 更新過程條件下一事件與右設限；BPT 只是理想化載入模型，不宣稱剛破裂區域安全。
12. PSHA 的平均超越率與至少一次機率分開；非齊次 Poisson 公式要求確定率和 Poisson 假設，不能塞入平均 ETAS 強度。Cox 的 E(exp(-A)) 也不能不加条件套到自激發。
13. PSHA code 改最低規模門檻時，使用同一有限最大規模的尾積分比，確保未正規化規模率不變；舊 code 用無上界 GR 倍率後又有限上界正規化，兩者不一致。

## 已做驗證

- 八個 Python 檔 AST parse：確認無誤。
- 本次章與附錄的 `{doc}` 目標均存在，無非預期控制字元：確認無誤。
- 對數常態比值獨立計算如上。
- eta/Delta 有限門檻數值積分（beta=ln10,m0=3,aM=1.2,bM=.7,sigma=.4,m=6）與解析值相對差 2.22e-16：確認無誤。
- PSHA 改門檻後在同一規模位置的「率×密度」相等：確認無誤。
- 未執行整章 notebook 或建置／瀏覽器驗收；交主代理整站處理。未逐項重新推算所有沿用文獻圖數值；日本增益圖已核對原 PDF，其餘既有來源有保留。

## Code 變更，需更新輸出

- `15_psi_phenomenon.py`：cumag 的規模增量正負號修正，影響其後依賴 onset 的圖。
- `18_testing_comparison.py`：只改第一圖題，將「跨零就是平手」改成「尚無法判定增益方向」。
- `21_psha.py`：source_nu 新增 m_max，使用有限上界尾積分比；所有呼叫傳入對應震源 mx。
- `22_operational_systems.py`：只改 Kaikōura 圖題，標歷史情境教學重繪，移除「必有一個會發生」過度概括。
- 16、17、19、20 的 code cell 未改；16 的台灣圖依賴既有快取，主代理應保留快取避免網路重抓。
