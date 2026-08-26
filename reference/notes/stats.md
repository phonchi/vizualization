# b 值與統計基礎論文筆記

> 主題：GR 律、b 值估計、點過程基礎 | 來源：reference/ 下 4 篇 | 供第二部 10 章（目錄統計進階）、09 章（導論）等使用

---

## [2021] Declustering 與 b 值 — The Effect of Declustering on the Size Distribution of Mainshocks

- **書目**：Leila Mizrahi、Shyam Nandan、Stefan Wiemer（Swiss Seismological Service, ETH Zurich）。全文標頭僅寫「Manuscript submitted to *Seismological Research Letters*」，日期 2020-12-16。**期刊卷期與 DOI 文中未見**（PDF 是投稿版手稿）。

- **這篇在做什麼**：問一個很少人認真問的問題——「主震目錄的 b 值到底是不是真的？」作者拿 1980 年起的加州 ANSS/ComCat 目錄，套上五類常用的 declustering 演算法（Reasenberg、三種 window method、Zaliapin 最近鄰、兩種 ETAS 隨機去叢），比較去叢前後的規模分布。接著再用**已知 b 值**的 ETAS 合成目錄重跑同一套流程。結論是：去叢後 b 值下降最多達 30%，而且合成資料（設計上主震與餘震共用同一個 b）也會出現同樣的下降——所以這個下降主要是演算法的產物，不是主震的物理性質。

- **關鍵觀念與公式**

  - **兩條 GR 直線不能同時成立**。若全目錄與主震目錄都被硬套上線性關係：

    $$\log_{10} N(m) = a - b\,m, \qquad \log_{10} N_{\text{main}}(m) = a_{\text{main}} - b_{\text{main}}\,m$$

    且 $b_{\text{main}} \neq b$，兩條線必在某個規模相交：

    $$m_{+} = \frac{a - a_{\text{main}}}{b - b_{\text{main}}}$$

    當 $b_{\text{main}} < b$（實際觀測到的情形），在 $m > m_{+}$ 之上模型會預測「主震數多於總地震數」——邏輯上不可能。加州資料算出的 $m_{+}$ 落在 **6.9 到 8.8** 之間，正好是工程上最關心的規模區間。

  - **為什麼會偏**。多數 declustering 演算法把「叢集中最大的事件」定義為主震。小地震比較不容易成為某個叢集的最大者，因此被移除的比例高於大地震，相對頻率往大規模端傾斜，斜率（b）自然變小。這是**選擇效應**，不是物理。理論上 Lombardi (2003)、Zhuang & Ogata (2006) 早就指出主震規模分布不是 Pareto 分布，全目錄 GR 律只在 $m \to \infty$ 的漸近意義下對主震成立。

  - **關鍵對照組：ETAS-Background**。作者刻意做兩種 ETAS 去叢：
    - `ETAS-Main`：叢集中最大事件 = 主震（沿用傳統定義）→ b 值明顯下降。
    - `ETAS-Background`：非被觸發的事件 = 主震（ETAS 原生定義，背景事件可以很小）→ **b 值與全目錄無顯著差異**。

    這組對照非常乾淨地把「b 值下降」歸因到「取最大值」這個動作本身。

  - **量級**：加州目錄全目錄 $b = 1.01$（$m_c = 3.6$，$\Delta M = 0.2$ binning，用 Clauset et al. (2009) 的方法聯合估 $m_c$ 與 $b$）；去叢後 b 值散布在 **0.73–1.00**，中間沒有明顯間隙。主震率在最「積極」與最「保守」的演算法之間差了 **6.1 倍**。

  - **合成資料的效應比真實資料更強**。原因推測是：所有 declustering 演算法都假設餘震空間分布是等向的——這對合成 ETAS 目錄成立，對真實地震不成立。所以合成目錄的叢集偵測更「成功」，小事件被剔得更乾淨。

- **教學上可用的洞見**

  - 「主震的 b 值比較低」這句在教科書層級流傳很廣的敘述，**至少有很大一部分是方法造成的**。Gulia et al. (2018) 說餘震序列 b 值平均比主震高 20%——這篇的證據顯示其中相當比例可以用去叢演算法解釋掉。教學時值得當成「觀測 ≠ 物理」的標準案例。
  - **警語**：報告 b 值時，一定要同時報告「用了哪個 declustering 方法、什麼參數、$m_c$ 取多少」。這三者任一改變，b 值就會動。作者的敏感度分析顯示，b 值下降在任何合理的 $m_c$ 選擇下都會出現，但**下降幅度是各演算法的特徵值**。
  - **PSHA 的三重缺陷**（可直接當講義的 bullet）：(1) 主震定義本身無法驗證；(2) 忽略餘震會低估危害度；(3) 被壓低的 b 值會在 $m > m_{+}$ 之上**高估**危害度。作者的說法很直白：「兩個錯誤不會變成一個對的——只有在 $m_{+}$ 這一個規模上才剛好抵銷。」
  - **替代路線**：直接用 ETAS 做地震率估計（Field et al. 2015；Nandan et al. 2019a）。ETAS 只需要全目錄的 GR 律，不必對「被任意挑出來的大事件」假設任何規模分布。這是本篇通往預報模型章節的自然接點。

- **與台灣的關聯**：無直接台灣案例（研究區為加州）。台灣的 b 值與餘震統計綜述另見 `reference/Taiwan/` 群組的 [2015] b-Values in Taiwan、[2016] Aftershocks in Taiwan 兩篇。

---

## [2024] Binned b 值估計 — The estimation of b-value of the frequency–magnitude distribution and of its 1σ intervals from binned magnitude data

- **書目**：S. Tinti（Università di Bologna）、P. Gasperini（Bologna + INGV Bologna）。*Geophysical Journal International* **238**(1), 433–458, 2024。DOI: `10.1093/gji/ggae159`。Open Access。程式碼：`https://github.com/pgaspy/b-value-testing`。

- **這篇在做什麼**：把 b 值估計從「continuous 近似」正式搬到「discrete（binned）」的理論地基上，並用大規模模擬（每組 K = 10 000 個合成樣本）把各家估計式排一次名。它同時是對 van der Elst (2021) b-positive 方法的第一份系統性檢驗——結論部分同意（用規模「差」比用規模本身穩健）、部分反駁（「非得用正差不可」缺乏證據）。附錄 A–I 提供離散指數分布與離散 Laplace 分布的完整動差推導，以及前人沒給的 1σ 區間公式。

- **關鍵觀念與公式**

  - **Aki (1965) MLE 與它的兩個偏差來源**。連續指數分布假設下：

    $$b = \frac{1}{\ln(10)\,(\bar{M} - M_c)}, \qquad \sigma_b = \frac{b}{\sqrt{N}}$$

    兩個偏差方向相反：
    - **Binning 讓 b 被高估**（規模被四捨五入到 0.1 格，$\bar{M}$ 被系統性壓低）。
    - **Incompleteness 讓 b 被低估**（$m_c$ 取太低，小事件缺漏，分布尾巴被削平）。

    最惡毒的一點在結論裡：**這兩個偏差有時會剛好互相抵銷**，讓未修正的 Aki 公式「看起來運作良好」——這是假象，不是驗證。

  - **Shi & Bolt (1982) 標準差**（比 Aki 的 $b/\sqrt{N}$ 誠實，因為它用資料的實際離散度）：

    $$\sigma_b = \ln(10)\, b^2 \sqrt{\frac{\sum_{i=1}^{N}\left(M_i - \bar{M}\right)^2}{N\,(N-1)}}$$

  - **Utsu (1966) 修正**：把 $M_c$ 往下挪半個 bin，$\delta = \Delta M / 2$（例如 0.05）：

    $$b = \frac{1}{\ln(10)\,(\bar{M} - M_c + \delta)}$$

    本篇附錄 F 證明：**Utsu 修正正好是精確式展開到二階的截斷**。所以 bin 小的時候它很好用，bin 大（$2\delta = 0.5$，例如由巨觀震度轉換而來的規模）時它會**系統性低估** b，模擬中連檢定門檻都過不了。

  - **離散情形的精確式**（Guttorp & Hopkins 1986 = Tinti & Mulargia 1987 = van der Elst 2021 的 coth 形式，三者代數上等價）：

    $$b = \frac{1}{2\delta \ln(10)} \ln\!\left(1 + \frac{2\delta}{\bar{M} - M_c}\right) = \frac{1}{\delta\ln(10)}\,\coth^{-1}\!\left(\frac{\bar{M} - M_c + \delta}{\delta}\right)$$

    這是 binned 目錄的正確預設選項。教學上可以直接告訴學生：**Aki 公式不要用在 binned 資料上**，模擬中它在 $N = 1000$、$b = 1$ 時給出 $\bar{b} = 1.126$，效能指標 $p = 0.0006$（遠低於 $\alpha = 0.05$）。

  - **規模差（magnitude differences）路線**。若規模服從離散指數分布，兩個規模之差服從**離散 Laplace 分布**。用絕對差的平均 $\overline{|\Delta M|}$：

    $$b = \frac{1}{2\delta\ln(10)}\,\mathrm{csch}^{-1}\!\left(\frac{\overline{|\Delta M|}}{2\delta}\right) = \frac{1}{2\delta\ln(10)}\ln\!\left[\frac{2\delta + \sqrt{4\delta^2 + \overline{|\Delta M|}^2}}{\overline{|\Delta M|}}\right]$$

    若做 **trimming**（丟掉 $|\Delta M| < \Delta M'_c$ 的差，最基本是丟掉零差，即 $\Delta M'_c = 2\delta$），估計式退化成與 binned 精確式同形，只是把 $\bar{M}, M_c$ 換成 $\overline{|\Delta M|}, \Delta M'_c$：

    $$b = \frac{1}{2\delta\ln(10)}\ln\!\left[\frac{\overline{|\Delta M|} - \Delta M'_c + 2\delta}{\overline{|\Delta M|} - \Delta M'_c}\right]$$

  - **差要怎麼取？獨立性的代價**。兩種取法：

    $$\text{(A)}\ \Delta M_i = M_{i+1} - M_i,\ i = 1,\dots,N-1 \qquad \text{(B)}\ \Delta M_i = M_{2i} - M_{2i-1},\ i = 1,\dots,N/2$$

    (A) 資料量最大但相鄰差共用同一個 $M_i$，**引入相關性**；(B) 保證獨立但資料減半。模擬結果很具體：用**絕對差**時，(A) 的理論 $\sigma$ 比實際散布小了約 22–23%，等於有效樣本數只剩 $N_e \approx 0.67N$——所以絕對差**必須用 (B)**。但用**單號差（只取正差或只取負差）**時，(A) 的相關性影響消失（$\bar{\sigma}_K / S_K \approx 1$），因此單號差**應該用 (A)** 以免資料被砍到四分之一。這是實作上很容易踩到的坑，值得寫進講義。

  - **1σ 區間（本篇新增）**。令 $c = 10^{2\delta \tilde{b}}$（$\tilde{b}$ 為估計值），則區間端點為：

    $$b_1 = \frac{1}{2\delta\ln(10)}\ln\!\left[\frac{c + c/\sqrt{N}}{1 + c/\sqrt{N}}\right], \qquad b_2 = \frac{1}{2\delta\ln(10)}\ln\!\left[\frac{c - c/\sqrt{N}}{1 - c/\sqrt{N}}\right]$$

    $$\sigma = \frac{\sigma_1 + \sigma_2}{2} = \frac{b_2 - b_1}{2}$$

    適用於 $N \gtrsim 30$–40（讓 $\bar{M}$ 或 $\overline{|\Delta M|}$ 近似常態）。注意 $\sigma_2 > \sigma_1$ 是系統性的——**b 值的分布本來就不對稱**，教學上不該把 $b \pm \sigma$ 講成對稱誤差棒。

  - **效能指標 $p$**（本篇自創，Student's t 的無母數版本）。給 K 個估計值 $\tilde{b}_i$、樣本均值 $\bar{b}_K$、真值 $b$：

    $$p = \begin{cases} L^{+}/K^{+} & \text{若 } \bar{b}_K < b \\ L^{-}/K^{-} & \text{若 } \bar{b}_K > b \end{cases}$$

    其中 $K^{\pm}$ 是大於／小於 $\bar{b}_K$ 的個數，$L^{\pm}$ 是大於／小於 $b$ 的個數。$p \in [0,1]$，當作虛無假設檢定的 p 值用，$p < 0.05$ 即判定估計式不可接受。

  - **完整目錄的結論**：精確式（規模）、絕對差、trimmed 差三條路線**表現相當**（$N=1000$、$2\delta=0.1$、$b=1$ 時 $\bar{b}$ 都在 1.001–1.002）。唯一出局的是未修正的 Aki 式。

  - **不完整目錄的結論**：規模的直方圖嚴重偏離指數，但**規模差的直方圖仍近似指數**——這是差分方法穩健性的視覺化根據（論文 Fig. 1 vs Fig. 2）。量化的門檻很實用：
    - 用**規模**的估計式（精確式）：需要 $M_c \geq M_{\text{maxc}} + 0.2$ 才給出正確 b 值。
    - 用**規模差**的估計式：$M_c \geq M_{\text{maxc}}$ 就夠了。

    前者正好驗證了文獻慣例（Wiemer & Wyss 2000；Mignan & Woessner 2012）「maximum curvature 加 0.2」的合理性。

  - **對 b-positive 的反駁**。van der Elst (2021) 主張只有正差可用。本篇在 (i) 合成不完整目錄、(ii) 完整度隨時間變化的合成餘震序列（$m = 5.6$ 主震、Omori-Utsu $p_O = 1$、$c_O = 0.01$）、(iii) 2016/10/30 義大利 Norcia $M_w$ 6.6 真實餘震序列（Horus 目錄，前 1000 個事件、50 km 內、17 小時內）三個層次上，都找不到「正差優於絕對差或負差」的證據。Norcia 序列三種取法給出 $b = 1.04 \pm 0.05$、$1.03 \pm 0.05$、$1.01 \pm 0.05$，彼此無法區分。**真正有效的是 trimming**：把 $\Delta M'_c$ 從 $2\delta$ 一路加到 $10\delta$，三種取法的偏差都單調縮小到幾乎可忽略，而且付出的樣本數代價遠比提高 $M_c$ 溫和。

- **教學上可用的洞見**

  - **樣本數不是萬靈丹**。論文報告了一個「表面上的悖論」：在完整度隨時間變化的餘震序列上，$N \gtrsim 4000$ 的大樣本估計**比小樣本更差**。原因是大樣本裡混進更多時間不完整的早期資料，偏差不會被平均掉、只會被固化。加大 trimming 後這個現象才消失。這是很好的反直覺教材：**系統偏差不隨 $\sqrt{N}$ 縮小**。
  - **決策樹**（可直接畫成講義流程圖）：完整目錄 + 小 bin → 精確式即可；bin 大於 0.1 → 一定要用精確式，不能用 Utsu；短期不完整（餘震序列早期）→ 用差分法 + trimming，並且把 trimming 門檻當成可調參數去測敏感度。
  - **最後的警語（論文自己寫的）**：以上一切只在「規模分布真的是指數」的前提下成立。Herrmann & Marzocchi (2020) 分析南加州與義大利中部的高解析目錄後指出，這個前提未必成立——如果不成立，**b 值這個概念本身就失去意義**。這句話很適合當作 10 章的收尾。

- **與台灣的關聯**：無直接台灣案例（實測資料為義大利 Norcia 序列）。但方法可直接搬到台灣 CWA 目錄——台灣目錄 $M_L$ 通常記到 0.1，即 $\delta = 0.05$，屬於「精確式與 Utsu 差異不大、但仍應用精確式」的情境。

---

## [2019] 時變過程回顧 — Modeling the earthquake occurrence with time-dependent processes: a brief review

- **書目**：Ourania Mangira、Christos Kourouklas、Dimitris Chorozoglou、Aggelos Iliopoulos、Eleftheria Papadimitriou（Aristotle University of Thessaloniki）。*Acta Geophysica*，2019。DOI: `10.1007/s11600-019-00284-4`。Review Article。

- **這篇在做什麼**：一篇緊湊的地圖式回顧，把「時間相依的地震發生模型」攤開成兩大家族——**fault-based**（單一斷層或斷層段的複發時間）與 **seismicity-based**（區域內所有震源的統計性質）——再依時間尺度（長期 vs 短期）細分。對教學最有用的地方是它把 renewal process、stress release、ETAS、EEPAS 放在同一張座標上比較，讓學生看得到 ETAS 不是憑空冒出來的，而是這條光譜上的一個位置。

- **關鍵觀念與公式**

  - **條件機率（renewal model 的核心）**。已知距上次事件已過 $T$ 年，未來 $\Delta T$ 內發生的機率：

    $$P(T \leq t \leq T + \Delta T \mid t > T) = \frac{\int_{T}^{T+\Delta T} f(t)\,\mathrm{d}t}{\int_{T}^{\infty} f(t)\,\mathrm{d}t}$$

    整個 fault-based 家族的差異，全在於**選哪個 $f(t)$**。

  - **三大複發時間分布**：

    - **Weibull**：$f(t \mid a, \beta) = \dfrac{\beta}{a}\left(\dfrac{t}{a}\right)^{\beta - 1}\exp\left[-\left(\dfrac{t}{a}\right)^{\beta}\right]$

      形狀參數 $\beta$ 的物理意義很漂亮：$\beta = 1$ 退化成指數分布（無記憶、Poisson）；$\beta < 1$ 代表**短期叢集**；$\beta > 1$ 代表**準週期**（fault-based 的典型情形）。**注意：原文用 $b$ 表示這個形狀參數，與 GR 的 b 值同名不同物**——教學時務必改用別的符號，這是很容易讓學生混淆的地雷。

    - **Lognormal**：$f(t \mid \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\,\sigma t}\exp\left[-\dfrac{(\ln t - \mu)^2}{2\sigma^2}\right]$，重尾分布，對超過平均複發期的長間隔給出偏高的機率。

    - **Brownian Passage Time（BPT）**：$f(t \mid \mu, \alpha) = \sqrt{\dfrac{\mu}{2\pi\alpha^2 t^3}}\exp\left[-\dfrac{(t-\mu)^2}{2\mu\alpha^2 t}\right]$

      源自 **Brownian Relaxation Oscillator**：應力 $X(t) = \lambda t + \sigma W(t)$（$\lambda$ 為定值加載率，$W(t)$ 為標準布朗運動），事件發生於 $X$ 首次觸及門檻。$\alpha$ 是 **aperiodicity**（等同變異係數）：$\alpha \to 0$ 趨近完全週期，$\alpha \to \infty$ 趨近完全隨機。實務研究多落在 $0.3 \leq \alpha \leq 0.7$。

  - **BPT 與 b 值的橋樑（Zöller et al. 2008）**。假設斷層上的中小地震（卸載）與斷層外地震（加載）效應大致相抵，則 aperiodicity 可以由儀器目錄的 b 值推出：

    $$\alpha = \sqrt{\frac{b}{3 - b}}, \qquad 0 < b < 3$$

    這是本組四篇之間最直接的連結：**b 值不只是規模分布的斜率，它可以回頭餵給複發時間模型**。$b = 1$ 給出 $\alpha = \sqrt{1/2} \approx 0.71$，恰好在文獻慣用範圍的上緣。

  - **hazard function 的形狀決定一切**（這是選分布時真正該看的東西，不是 pdf）：
    - **BPT**：事件剛發生後 hazard 極低 → 隨時間上升 → 在平均複發期附近達最大 → 之後漸近趨於 $1/(2\mu\alpha^2)$（有限值）。物理上最合理，也是 UCERF2/UCERF3 採用的模型。
    - **Weibull**（$\beta > 1$）：hazard 單調遞增，永不飽和。
    - **Lognormal**：hazard 上升到最大後**漸近趨於零**——意思是「等愈久愈不會發生」，這對地震顯然不合物理。

    Polidoro et al. (2013) 的觀察值得引用：當「距上次事件的時間」約為複發期的一半時，各模型結果差不多；**時間拖得愈久，選錯分布的代價愈大**。

  - **Stress Release Model（Vere-Jones 1978）——連結物理與點過程**。應力水準 $X(t) = X(0) + \rho t - S(t)$（$\rho$ 加載率，$S(t)$ 為累積釋放量），條件強度取指數形式：

    $$\lambda^{*}(t) = \exp\{a + b[t - c\,S(t)]\}$$

    **Linked SRM**（Liu et al. 1998）加入子區域之間的應力轉移：

    $$\lambda_i^{*}(t) = \exp\left\{a_i + b_i\left[t - \sum_j c_{ij} S(t, j)\right]\right\}$$

    這是一個 **self-correcting**（自我修正）點過程：事件發生後強度下降。與 ETAS 的 **self-exciting**（自我激發，事件後強度上升）**方向相反**——這組對比是講點過程時最有力的一刀。

  - **ETAS（Ogata 1988, 1998）——self-exciting / Hawkes 家族**。時間型：

    $$\lambda(t, m) = \beta e^{-\beta(m - m_0)}\left\{\mu + A\sum_{i:\,t_i < t} e^{\alpha(m_i - m_0)} f(t - t_i)\right\}$$

    其中 $f(t) = \dfrac{p-1}{c}\left(1 + \dfrac{t}{c}\right)^{-p}$ 即 Omori–Utsu 律的機率密度，$\mu$ 是背景率，$\alpha$ 控制產能（大地震生多少子代），而 **$\beta = b\ln 10$ 直接把 GR 律嵌進條件強度裡**。時空型再加入背景空間密度 $h(x)$ 與觸發空間核 $g(x - x_i)$。

    ETAS 的哲學價值（本篇講得很清楚）：**它不需要事先區分主震／前震／餘震**。既然這個分類本身是任意的（見 [2021] 那篇的證據），把它從模型裡拿掉就是一種進步。

  - **EEPAS（Evison & Rhoades 2004）**：

    $$\lambda^{*}(t, m, x) = \mu\,\lambda_0(t, m, x) + \sum_{t_i < t} w_i\,\eta(m_i)\,r(M \mid M_i)\,f(t - t_i \mid M_i)\,g(x - x_i \mid M_i)$$

    形式上和 ETAS 幾乎一樣，**但 $f$ 和 $g$ 不是 Omori 律**——它們來自 precursory scale increase（前兆尺度增長）的經驗尺度關係，描述的是**大震前數月到數十年**的中小地震活動增強，而不是大震後的餘震衰減。所以：**ETAS 往後看（觸發），EEPAS 往前看（前兆）**，兩者在時間尺度與物理假設上是互補而非競爭。權重 $w_i$ 通常設為 1，但可由 ETAS 隨機去叢得到。

  - **PPE（Proximity to Past Earthquakes）** 是 EEPAS 常用的參考率密度（虛無假設）：

    $$\lambda_{\text{PPE}}(t, m, x, y) = g_0(m)\,h(t, x, y), \quad g_0(m) = \beta e^{-\beta(m - m_c)}$$
    $$h(t,x,y) = \frac{1}{t - t_0}\sum_{t_i < t} h_{\text{eq}}(i), \qquad h_{\text{eq}}(i) = \frac{a\,(m_i - m_c)}{\pi (d^2 + \Delta_i^2)} + s$$

    $d$ 是平滑距離、$s$ 是空間均勻的背景率、$\Delta_i$ 是到第 $i$ 個歷史地震的距離。同樣地 $\beta = b\ln 10$。

  - **其他家族（點到為止即可）**：Accelerating Moment Release（Benioff 應變 $\varepsilon(t) = \sum\sqrt{E_i} = A - B(t_f - t)^m$，$0.1 \leq m \leq 0.5$）、hidden Markov / semi-Markov model（隱藏狀態 = 應力水準）、Markovian Arrival Process、rate-and-state 約束的叢集模型（Dieterich 1994）、double branching model。

- **教學上可用的洞見**

  - 這篇最適合當 **09 章導論**的骨架：先建立「Poisson = 無記憶 = 虛無假設」，再沿著兩個軸把模型攤開——**時間尺度**（長期複發 ↔ 短期觸發）與 **記憶方向**（self-correcting ↔ self-exciting）。ETAS 位在「短期 + self-exciting」，renewal/BPT 位在「長期 + self-correcting」，EEPAS 位在「中長期 + 前兆型」。
  - **警語**：fault-based 模型的資料量通常極小（一條斷層的歷史破裂記錄可能只有 3–5 次），此時分布選擇幾乎完全由假設而非資料決定。UCERF3 用 5760 個 logic-tree 分支來表達這種認知不確定性，本身就說明了問題有多嚴重。
  - **CSEP 的角色**：本篇結尾強調，模型選擇不能靠作者自評，必須交給 prospective、預先約定資料集與統計檢定的第三方測試平台。這是通往 pyCSEP／預報檢定章節的接點。

- **與台灣的關聯**：無直接台灣案例（實例多為希臘、日本、義大利、加州）。

---

## [2025] SeismoStats — SeismoStats: A Python Package for Statistical Seismology

- **書目**：Aron Mirwald、Nicolas Schmid、Leila Mizrahi、Marta Han、Alicia Rohnacher、Vanille A. Ritz、Stefan Wiemer（Swiss Seismological Service, ETH Zurich）。arXiv:`2511.04521v1` [physics.geo-ph]，2025-11-06。**期刊 DOI 文中未見**（預印本）。文件：`https://seismostats.readthedocs.io/`；程式碼：`https://github.com/swiss-seismological-service/SeismoStats`。

- **這篇在做什麼**：一份軟體論文，目標很明確——**取代 ZMAP**。ZMAP 綁 MATLAB、已停止維護，且沒有實作近年（van der Elst 2021、van der Elst & Page 2023、Lippiello & Petrillo 2024）在 a/b 值估計上的進展。SeismoStats 用 Python 重新實作 GR 律三參數（$m_c$、a 值、b 值）的估計，並附上目錄下載（FDSN Web Services、QuakeML）、繪圖與工具函式。對教學來說，這篇的價值在於它把「該教哪些方法」整理成一份**已實作、可跑、有測試的清單**。

- **關鍵觀念與公式**

  - **GR 律的寫法**（注意這裡的 a 值是相對 $m_c$ 定義的）：

    $$N(m) = 10^{\,a - b(m - m_c)}$$

  - **架構**：核心是 `Catalog` 類別（pandas `DataFrame` 的子類別，至少需要 `magnitude` 欄），加上三個子套件 `analysis` / `plots` / `utils`。a 值與 b 值估計器實作成**類別**（因為共用輸入介面：magnitudes、$m_c$、bin size $\Delta m$，共用屬性 `n`、`std`、`magnitudes`）；$m_c$ 估計則是**函式**（共通點太少），統一回傳 `(mc, dict)` 的 tuple 以便互換。

  - **三種 $m_c$ 估計法**：

    - **Maximum Curvature（MAXC）**：取非累積 FMD 的峰值再加修正項。

      $$m_c = \arg\max_{m_i} N(m_i) + \delta, \qquad \delta = 0.2 \ \text{(預設)}$$

      注意套件刻意區分兩個「bin」：`delta_m`（$\Delta m$，資料本身的離散化格距）與 `fmd_bin`（$\Delta m^{*}$，計算 $m_c$ 時的直方圖寬度）。兩者可以不同（例如 $\Delta m = 0.01$ 但 $\Delta m^{*} = 0.1$）——這是實作上很常見的混淆點。

    - **b 值穩定度（Cao & Gao 2002；Woessner & Wiemer 2005）**：利用「不完整目錄的 b 值被低估、抬高截切門檻後 b 值會上升」這個性質，取 b 值開始穩定的那一點：

      $$m_c = \min\left\{ m_i \ \middle|\ \left|\frac{1}{K}\sum_{k=1}^{K} b(m_i + k\Delta m^{*}) - b(m_i)\right| < \sigma_{b(m_i)} \right\}$$

      $L = K \cdot \Delta m^{*}$ 是判定穩定所用的規模範圍，預設 0.5。

    - **KS 距離（goodness-of-fit 路線；Clauset et al. 2009；Mizrahi et al. 2021）**：比較觀測 CDF 與理論指數 CDF：

      $$m_c = \min\left\{ m_i \ \middle|\ p(D_{\text{KS}}^{i}) \geq p_{\text{th}} \right\}, \qquad p_{\text{th}} = 0.1 \ \text{(預設)}$$

      $p(D_{\text{KS}})$ 由模擬（預設 10 000 次）估得。**陷阱**：KS 檢定對樣本數敏感——目錄愈大，愈微小的偏離都會被判為顯著，導致 $m_c$ 被推得過高。對超大目錄這個方法可能不實用。

  - **四種 b 值估計器**：

    - `ClassicBValueEstimator`——離散精確 MLE（Tinti & Mulargia 1987；Marzocchi & Sandri 2003）：

      $$\hat{b} = \frac{1}{\ln(10)\,\Delta m}\ln p, \qquad p = 1 + \frac{\Delta m}{\frac{1}{n}\sum_{i=1}^{n}(m_i - m_c)}$$

      $\Delta m \to 0$ 時退化成 Aki 的 $\hat{b} = \log e \big/ \left[\frac{1}{n}\sum (m_i - m_c)\right]$。

    - `UtsuBValueEstimator`——近似式 $\hat{b} = \log e \big/ \left[\frac{1}{n}\sum(m_i - m_c) + \Delta m / 2\right]$。

    - `BPositiveBValueEstimator`——**b-positive**（van der Elst 2021）。動機是 **STAI（short-term aftershock incompleteness，短期餘震不完整）**：大地震後波形重疊，小地震偵測不到，等於 $m_c$ 在事件後短暫暴增。b-positive 的核心假設是「**每一刻的完整度至少是上一個記錄到的規模再加一點餘裕 $\delta m_c$**」。作法：取連續規模差 $m_i - m_{i-1}$，只保留大於 $\delta m_c$ 的正差，再套用經典 MLE——理論依據是兩個指數變數之差服從 Laplace 分布。**限制**：它只解決 STAI，解決不了「網路偵測能力造成的一般不完整」（Lippiello & Petrillo 2024）；補救方式是先截到 $m_c$ 之上再取差，或加大 $\delta m_c$。

    - `BMorePositiveBValueEstimator`——**b-more-positive**（Lippiello & Petrillo 2024）。對每個 $m_i$，往後找**第一個** $\geq m_i + \delta m_c$ 的事件來取差，因此差的數量比 b-positive 多。但作者的合成測試發現**實際標準差大於式(6)的理論值**，所以這個估計器改用 bootstrap 算誤差。

  - **標準差（Shi & Bolt 1982）**，套件對所有方法（b-more-positive 除外）統一使用：

    $$\sigma(\hat{b}) = \frac{\ln(10)\,b^2}{\sqrt{n_m - 1}}\sqrt{\mathrm{var}(m)}$$

    **作者自己標註的妥協**：這個式子是在連續規模假設下推導的，套件卻在離散規模上照用——理由是「沒有考慮離散化的替代方案，而且大家都這樣用」。對照 [2024] Tinti & Gasperini 那篇提供的離散 1σ 公式（$b_1, b_2$），這正是一個「文獻已經解決、但軟體還沒跟上」的缺口，很適合當作業題目。

  - **a 值估計**：$a = \log N(m_c)$。兩個常見修正：換算到參考規模 $a_{m_{\text{ref}}} = a - b(m_{\text{ref}} - m_c)$；以及用 scaling factor 換算成單位時間（率）或單位體積。另有 **a-positive** 與 **a-more-positive**（van der Elst & Page 2023），概念是估計「被用到的 inter-event 時間佔總觀測時間的比例」：

    $$a^{+} = \log n^{+} - \log\frac{\sum_{i=1}^{n^{+}} \Delta t_i}{T}$$

    a-more-positive 還要依 GR 律縮放時間差 $\tau_i = \Delta t_i \, 10^{-b(m_i + \delta m_c)}$，並把「後面沒出現更大事件」的開放區間 $T_j = (T - t_j)10^{-b(m_j + \delta m_c)}$ 一併算進去以免偏差。**a 值估計器沒有 `std` 屬性**——因為多數文獻根本沒給 a 值的誤差估計。

  - **檢驗工具**（這三個很值得單獨教）：
    - **Shi & Bolt 標準差**——b 值的不確定度。
    - **Lilliefors 檢定**（Lilliefors 1969；Herrmann & Marzocchi 2021）——檢定樣本是否真的服從指數分布，也就是**檢查 GR 律這個前提本身**。
    - **b-significance 方法**（Mirwald et al. 2024）——檢定 b 值的變化是否**顯著**，而不是把任何波動都當成訊號。

- **教學上可用的洞見**

  - 這是 10 章最好的**動手工具**：從 FDSN 下載目錄 → `plot_fmd` / `plot_cum_fmd` 看 FMD → 三種 $m_c$ 方法畫在同一張圖上（論文 Fig. 3 就是這樣做的）→ 比較四種 b 值估計器 → 用 `plot_b_series_constant_nm` 畫 b 值時間序列，classic 與 b-positive 疊在一起看差異。整條 pipeline 可以直接搬成一節課的 notebook。
  - **「b 值在變」的舉證責任**。套件把 Lilliefors 檢定與 b-significance 檢定放在同一層級，等於在說：宣稱 b 值有時空變化之前，你得先證明 (a) 分布真的是指數，(b) 變化幅度超過估計不確定度。論文 Pitfalls 一節直接引用 Marzocchi et al. (2020)〈How to be fooled searching for significant variations of the b-value〉——這個標題本身就是最好的講義標語。
  - **軟體不會幫你檢查資料**（作者原話）。已知的資料品質地雷：$m_c$ 的**時空變化**（套件目前所有 $m_c$ 方法都假設 $m_c$ 不隨時間變）、採石場爆破等人為事件、目錄裡**混雜多種規模型別**（$M_L$ / $M_w$ 轉換並不直接，見 Deichmann 2017）。
  - 生態系定位（可畫成一張圖）：ObsPy 管波形、SeisComP 管即時、SeisBench 管機器學習、pyCSEP 管預報檢定、OpenQuake 管危害度——**SeismoStats 補的是「目錄的統計分析」這一塊**。

- **與台灣的關聯**：無直接台灣案例（範例資料為 2024 年瑞士目錄）。但 FDSN 下載介面理論上可直接接台灣的資料服務，`Catalog` 類別只要求 `magnitude` 欄，把 CWA 目錄讀成 DataFrame 即可套用全部功能。

---

## 跨篇綜合：這個主題教什麼

### 一、b 值估計的教學順序

建議按「**每一步都在修前一步的一個錯**」來鋪陳，讓學生看到方法演進的動機：

1. **GR 律與 b 值的意義**（[2025]）。$N(m) = 10^{a - b(m - m_c)}$，b 值 ≈ 1 代表規模每降 1 級、次數增為 10 倍。先建立直覺：**b 低 = 大地震佔比高**。
2. **Aki (1965) MLE**（[2024]）。$b = 1/[\ln(10)(\bar{M} - M_c)]$，推導一次，讓學生知道它從「連續指數分布」來。
3. **第一個問題：$m_c$ 從哪來？**（[2025]）。介紹 MAXC + 0.2、b 值穩定度、KS 距離三種路線，並強調 **$m_c$ 與 b 值是耦合問題**——不知道 $m_c$ 就估不了 b，反之亦然。實作上畫「$b(m_c)$ 曲線」是最直觀的教法。
4. **第二個問題：規模是離散的**（[2024]）。Utsu 半格修正 → 精確式。點明 Utsu 只是精確式的二階截斷，bin 大於 0.1 時會失效。
5. **不確定度**（[2024] + [2025]）。Aki 的 $b/\sqrt{N}$ → Shi & Bolt → 離散情形的非對稱 1σ 區間。強調 **b 值的誤差棒本質上不對稱**。
6. **第三個問題：短期不完整**（[2025] 的方法 + [2024] 的檢驗）。STAI → b-positive / b-more-positive → trimming。
7. **回頭質疑前提**（[2024] 結尾 + [2025] 的 Lilliefors）。如果規模分布根本不是指數，前面六步全部作廢。

### 二、哪些「常識」已被近年文獻推翻或大幅削弱

| 流傳的說法 | 近年文獻的修正 | 出處 |
|---|---|---|
| 「主震的 b 值天生比較低」 | 大部分是 declustering 演算法「取叢集最大值」造成的選擇效應；用 ETAS-Background 定義主震時 b 值不變 | [2021] |
| 「binning 造成的偏差很小，可以忽略」 | 模擬證實未修正的 Aki 式在 $N=1000$、$\Delta M=0.1$ 就已顯著高估（$\bar{b} = 1.126$，$p = 0.0006$）；明確反駁 Marzocchi et al. (2020) 的「可忽略」說法 | [2024] |
| 「Utsu 修正夠用了」 | 只是二階截斷；$2\delta = 0.5$ 時系統性低估到不可接受 | [2024] |
| 「b-positive 必須用正差」 | 合成資料與 Norcia 真實序列都顯示正差／絕對差／負差表現相當；**真正有效的是 trimming 門檻，不是取正號** | [2024] vs [2025] |
| 「樣本愈大估計愈準」 | 完整度隨時間變化時，大樣本反而更差——系統偏差不隨 $\sqrt{N}$ 縮小 | [2024] |
| 「b 值時空變化就是應力變化的訊號」 | 必須先過 Lilliefors（分布真的是指數嗎）與 b-significance（變化超過不確定度嗎）兩關；$m_c$ 選擇、規模型別混用、爆破事件都會製造假訊號 | [2025] |
| 「lognormal 適合描述複發時間」 | 它的 hazard function 最終趨於零（等愈久愈不會發生），物理上不合理；BPT 的 hazard 趨於有限正值才合理 | [2019] |

**最需要反覆強調的一句**：b 值對 $m_c$ 的選擇**極度敏感**。$m_c$ 取太低 → b 被低估；取太高 → 樣本銳減、誤差爆炸。任何 b 值報告若沒附上 $m_c$ 的決定方式與敏感度分析，都不能當成結論。

### 三、與預報模型（ETAS / EEPAS）的銜接點

- **b 值直接住在 ETAS 裡**。ETAS 條件強度中的 $\beta = b\ln 10$——b 值估錯，觸發率與規模分布一起錯。這是「為什麼要花一整章講 b 值」最有力的答案。
- **b 值也能餵給複發模型**。Zöller et al. (2008) 的 $\alpha = \sqrt{b/(3-b)}$ 把儀器目錄的 b 值接到 BPT 的 aperiodicity，是長期危害度與短期統計之間少見的量化橋樑（[2019]）。
- **declustering 是可以繞開的**。[2021] 論證主震目錄的 b 值不可信；[2019] 指出 ETAS 的設計哲學正是**不需要區分主震／餘震**。兩篇合起來是一個完整的論證：與其在任意的分類上做統計，不如用能同時描述背景與觸發的點過程模型。這條線直通 ETAS 章節。
- **self-exciting vs self-correcting 是理解點過程的主軸**（[2019]）。Stress release（事件後強度下降）與 ETAS（事件後強度上升）在同一個 $\lambda^{*}(t)$ 框架下方向相反；Poisson 是兩者的共同退化情形。用這一組對比開場，比直接寫 Hawkes 過程定義好懂得多。
- **ETAS 與 EEPAS 是互補而非競爭**（[2019]）。條件強度形式幾乎相同，差別在 $f$、$g$ 的來源：ETAS 用 Omori 律往後看觸發（時間尺度：天到月），EEPAS 用 precursory scale increase 往前看前兆（時間尺度：月到數十年）。教學時把兩者的 $\lambda^{*}$ 並排寫出來，差異一目了然。
- **軟體鏈路**（[2025]）。SeismoStats 負責目錄統計（$m_c$、a、b）→ 這些參數是 ETAS/EEPAS 的輸入 → pyCSEP 負責預報結果的統計檢定。作者也明說改善與 pyCSEP 的互通性是後續目標。這條 pipeline 可以當作整個第二部的骨架圖。
- **模型選擇必須是 prospective 的**（[2019]）。CSEP 的設計原則——預先約定資料集與檢定方法、由第三方執行——是本組筆記通往預報檢定章節的最後一塊拼圖。
