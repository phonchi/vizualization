# 第二部撰寫規範（內部文件，不列入 _toc.yml）

本檔是第二部（09–23 章）的**單一事實來源**：記號、公式所有權、風格與排版規則。
撰寫或修改任何一章之前先讀這裡。讀者版的記號摘要放在 10.8 節。

---

## 1. 統一記號

### 1.1 核心

| 概念 | 本部寫法 | 說明 |
|---|---|---|
| 事件 | $(t_i, x_i, y_i, m_i)$ | 時間、位置、規模 |
| 歷史 | $H_t$ | $t$ 之前的全部事件 |
| 條件強度 | $\lambda^*(t,x,y,m) \equiv \lambda(t,x,y,m \mid H_t)$ | 星號一律代表「條件於歷史」 |
| 期望數 | $\Lambda$ | 積分後的量（$\Lambda(t)=\int_0^t\lambda^*$、網格 $\Lambda_{jk}$） |
| 觀測數 | $\omega_{jk}$ | 空間格 $j$ × 規模箱 $k$ |
| 目錄完整度 | $M_c$ | 目錄的性質（估計量） |
| 模型輸入門檻 | $m_0$ | 模型設定，通常 $m_0 \ge M_c$ |
| 預報目標門檻 | $m_T$ | 目標地震集的下限 |
| GR 斜率 | $b$；$\beta = b\ln 10$ | $\beta$ **只**代表 $b\ln 10$ |
| 規模密度 | $s(m)=\beta e^{-\beta(m-m_0)}$ | 積分為 1 |
| Omori 未正規化 | $n(t)=K(t+c)^{-p}$ | $K$ 帶量綱 |
| Omori 密度 | $g(t)=\frac{p-1}{c}\left(1+\frac{t}{c}\right)^{-p}$ | 積分為 1，需 $p>1$ |
| 產能 | $\kappa(m)=A\,e^{\alpha(m-m_0)}$ | $A$ 搭配正規化核；**不可與 Ogata 的 $K$ 直接比較** |
| 空間核 | $f(x,y;m)=\frac{q-1}{\pi D e^{\gamma(m-m_0)}}\left[1+\frac{r^2}{D e^{\gamma(m-m_0)}}\right]^{-q}$ | |
| ETAS 背景率 | $\mu(x,y)$ | 恆帶引數 |
| 分支比 | $n = \dfrac{A\beta}{\beta-\alpha}$ | 需 $\alpha<\beta$ |
| 觸發／背景機率 | $\rho_{ij}$ / $\phi_j$ | 隨機除叢 |

### 1.2 必須改名的衝突符號（文獻沿用會撞名，本部一律改寫）

| 文獻用法 | 本部寫法 |
|---|---|
| EEPAS 混合權重 $\mu$ | $\mu_E$ |
| BPT 平均複發時間 $\mu$ | $T_r$ |
| BPT aperiodicity $\alpha$ | $c_v$ |
| Weibull 形狀參數 $\beta$（或 $b$） | $k$（尺度用 $\theta$） |
| 應力釋放模型的 $a, b, c$ | $a_s, b_s, c_s$ |
| Janus 混合權重 $q$ | $\pi_{\rm ETAS}$ |
| EAS 餘震規模指數 $\alpha$ | $\beta_A$ |
| 各種機率 $p$ | 一律 $P(\cdot)$（$p$ 保留給 Omori 指數） |
| 乘法 hybrid 的 $b_i, c_i$ | $u_i, v_i$ |

### 1.3 專屬符號

- **Ψ**：$M_P, T_P, A_P, M_m, M^-, C(t), r$（率比）、$Z$
- **EEPAS**：$a_M, b_M, \sigma_M$／$a_T, b_T, \sigma_T$／$b_A, \sigma_A$／$\eta(m), \Delta(m), w_i, \mu_E$
- **檢驗**：IGPE、IGPEc、IGPA、$I_1$、$AS(\tau)$、$\tau$、$\nu$（Molchan 漏報率，只在 18.4 用）、$S_B$（Brier）
- **組合**：$\pi_i$（權重，$\sum\pi_i=1$）、$\lambda_1$（baseline）、$f_i$（保序轉換）、$a_0$（正規化）
- **危害**：$\lambda_{IM}(x)$（年超越率）、$\sigma_{\ln IM}$、$\varepsilon$、$T_R$（回歸期）

### 1.4 規模尺度
目錄規模一律標明尺度（$M_L$／$M_w$／$M_D$）；模型式中的 $m$ 表示「已宣告尺度的規模」。
跨目錄／跨模型比較前先出示轉換式（11.1 擁有）。

---

## 2. 公式所有權（防重複；其他章一律引用不重推）

| 公式／主題 | 擁有章節 |
|---|---|
| $\ln L=\sum\ln\lambda^*-\int\lambda^*$、隨機時間變換、殘差、反函數法抽樣 | 10 |
| GR 律、$\beta=b\ln10$、b 值 MLE 家族、$M_c$ 三法 | 11 |
| Omori–Utsu 與正規化、Utsu–Seki、Båth 三估計法、除叢與 $m_x$ | 12 |
| 分支比 $n$、世代分解 | 13 |
| $\rho_{ij}$／$\phi_j$、R–J 特例、simplETAS 釘死理由 | 14 |
| $C(t)$、Ψ 三迴歸、時空抵換與回歸稀釋 | 15 |
| EEPAS 三核、$\eta(m)$、$\Delta(m)$、PPE／SUP、常態卷積 | 16 |
| POLL／BILL／負二項、N/M/S/cL、模擬程序 | 17 |
| IGPE／IGPEc／$I_1$／Molchan／$AS$／功效／可靠度／成本–損失 | 18 |
| 凸組合守恆、乘法 hybrid、logistic 權重 | 19 |
| hazard function、BPT 首達、應力釋放、Zöller 橋 | 20 |
| 危害積分、截斷 GR、反聚合、NHPP | 21 |

---

## 3. 章節骨架（15 章一致）

```
# N. 章名
（開場：從上一章接續的動機，2–3 段散文；一句「這章會用到前面哪一條式子」）
## N.1 … N.k    正文：先動機 → 再數學 → 再限制
## 參數與典型值（表）      ← 有參數的章
## 常見誤解與陷阱
## 研究前沿與未解問題
## 附錄：本章推導細節      ← 長代數放這裡，正文只留結果與它改變了什麼
（結尾：反思式散文，以 {doc} 前指下一章）
```

零習題、零 quiz。admonition 每章 ≤1。

---

## 4. 排版規則（**違反會導致渲染失敗**）

1. **粗體**：`**…**` 的**開頭前一字元**與**結尾後一字元**若為標點，容易觸發 CommonMark
   flanking 失敗。規則：標點一律放在粗體**外面**——寫 `**內文**。接著`，不要寫
   `**內文。**接著`；`「**內文**」` 不要寫 `**「內文」**`。
2. **數學**：本專案未啟用 `amsmath`。多行推導一律用
   `$$\begin{aligned} … \end{aligned}$$`，不要用 `\begin{align}`。
3. 數學式內**不放中文**（`\text{中文}` 字型可能缺字）；要標註就寫在行文或圖說。
4. 中文與行內 `$…$` 之間留一個半形空格。
5. 表格欄數 ≤5（窄螢幕會爆版）；更多欄就拆兩張表。
6. 需要跨章引用的式子才編號（`$$…$$ (eq:label)` ＋ `{eq}`），全書約 20 條，其餘不編。
7. 每行約 70–75 個全形字寬折行（沿用既有風格）。

---

## 5. 程式與圖規格

- 章首樣板 cell（`# %% tags=["remove-input"]`）一律：
  ```python
  from gdms_toolkit.viz import setup_plotly
  setup_plotly()
  ```
- 所有繪圖 cell 加 `tags=["hide-input"]`；圖以 Plotly 為主，套 `apply_layout()`，
  用 `ACCENT`／`QUAKE_COLOR`（紅只給事件標記）／`SEQUENTIAL="Blues"`／`PALETTE`。
- **每章圖 ≤6–8 張**；每個 cell 執行 **< 30 秒**；亂數一律 `np.random.default_rng(seed)`。
- 資料：`gdms_toolkit.load_taiwan_catalog()`（1973–2025，35 萬筆，M≥2）、
  `data/cache/catalog_2024spring.csv`。**禁止對 35 萬筆做 $O(N^2)$**；
  最近鄰用 `scipy.spatial.cKDTree` 或取十年子集；網格統計用 `np.histogram2d`。
- 圖說中的數字一律由程式以 f-string 帶入，不要手寫（避免圖文不一致）。
- 可用套件：numpy、scipy、pandas、plotly、matplotlib、obspy、folium。
  **沒有** sklearn／statsmodels／pycsep／geopandas。

---

## 6. 內容紅線

- **Embargo（審稿中，任何數字／參數／結果／圖表一律不得出現）**：
  `[2026] PyEEPAS`(+sup)、`[2026] Taiwan_EEPAS`(+sup)、`[2025] EEPAS_Software`。
  Ψ 自動辨識演算法只引 **Christophersen et al. (2024, SRL)**；台灣 EEPAS 只能寫
  「在地化工作正在進行中」，不得提及軟體、開源化、參數或結果。
  （`EEPAS_TW-main.zip` 內的 **CWA 公開目錄資料**可用，其 EEPAS 結果不可用。）
- **未核實項不得寫入**（見 `reference/notes/psha_step_websearch.md`）：STEP 的 5 km 網格、
  USGS「20 分鐘首報／第一年更新 75 次」的出處、AIC 權重的具體公式、
  `generic/sequence-specific/spatially-varying` 是否為 2005 原文用語、TEM PSHA2025 卷期。
- Baker 是**五**步驟（常見四步驟是合併版，要註明）；**STEP 不是 USGS 現行系統**
  （現行為 OAF，引擎是 R–J 或 ETAS）。
- 引用一律用已發表文獻的作者與年份；找不到出處的數字不要寫。

---

## 7. 素材對照（`reference/notes/`）

`etas.md`→13,14｜`stats.md`→11,20｜`taiwan.md`→11,12,14,23｜`psi.md`→15｜
`eepas.md`→16｜`testing.md`→17,18｜`ensemble.md`→19｜`oef.md`→09,22｜
`psha_step_websearch.md`→21,22｜`INDEX.md`（索引與檔案勘誤）
