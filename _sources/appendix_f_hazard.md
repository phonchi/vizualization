
# 附錄 F：複發時間與地震危害積分

本附錄需要條件機率、積分與基本隨機過程。對應 {doc}`預報之後 <20_beyond_forecast>`。本附錄接收複發模型與 PSHA 的完整說明；先從等待時間讀起，再將事件率轉成場址地動。

## F.0 從斷層等待時間走向場址地動

一段斷層已經多久沒有大破裂，可能是下一次破裂的資訊。
要使用它，先定義哪些破裂算同一類事件，以及古地震年代的不確定性。
這與 ETAS 描述區域內一串觸發活動是不同分析物件。
同一斷層段剛破裂，不表示周邊區域沒有餘震危害。

更新過程（renewal process）把相鄰指定事件的間隔視為獨立同分布。
事件發生後，模型中的年齡歸零。
指數間隔具有無記憶性，等待多久都不改變下一段相同長度的條件機率。
Weibull、對數常態與 BPT 則讓等待時間影響未來。

| 分布 | 主要參數 | 等待很久時的危害率 |
|---|---|---|
| Weibull | 形狀 k、尺度 a | k>1 持續上升；k=1 固定；k<1 下降 |
| 對數常態 | 對數平均與標準差 | 最終趨近零 |
| BPT | 平均間隔 Tr、變異係數 cv | 趨近正的有限常數 |

分布相似的中央部分，可能配上不同尾端。
資料只有少數完整間隔時，不能只看最佳擬合線就認定尾部機制。
右設限指觀測結束時下一次破裂還沒來；那段等待仍包含資訊，須放入概似。
年代誤差、漏掉古事件與起始時刻的選擇，也須保留在解讀中。

BPT 是 Brownian passage time，從平均累積並帶布朗擾動的載入量出發。
首次到達門檻就記一次破裂。用平均間隔 $T_r$ 與變異係數 $c_v$ 表示，

$$f(t)=\sqrt{\frac{T_r}{2\pi c_v^2t^3}}
\exp\left[-\frac{(t-T_r)^2}{2c_v^2T_rt}\right],\qquad t>0.$$ (eq:bpt-pdf)

這是理想化機制，不是直接量得全部斷層應力。F.3 推導其首達時間。

### PSHA 的五個步驟

機率式地震危害分析（probabilistic seismic hazard analysis，PSHA）
將各種地震情境按發生率加權，回答場址地動超過某門檻的頻率。
地動指標記為 $IM$，例如尖峰地表加速度；它與地震規模不同。
依 Baker 教材，工作分為五步：

1. 辨認可能影響場址的震源及各自事件率。
2. 描述震源的規模分布，指定上下限。
3. 描述震源與場址的距離及破裂幾何。
4. 用地動模型給出各情境下的地動分布。
5. 對情境與隨機變異積分，再加總震源。

地動模型（ground-motion model）在指定規模、距離與場址條件下，
給出地動的中央趨勢和散布。一次近處中型地震與一次遠處大震都可能造成明顯地動。
因此不能只挑最大規模，也不能只挑最常見情境。

$$\lambda_{IM}(x)=\sum_s\nu_s\int\!\!\int
P(IM>x\mid m,r,s)f_{M,R\mid s}(m,r)\,dm\,dr.$$ (eq:psha)

內部積分是一次事件的超越機率，乘年事件率後得到年超越率。
危害曲線逐一改變 $x$；不同門檻可能由不同規模和距離主導。
F.5 的反聚合把各種情境的貢獻拆開。它解釋危害來源，不預言下一場地震。

使用上下截斷的 GR 密度時，正規化形式為

$$f_M(m)=\frac{\beta e^{-\beta(m-m_{\min})}}
{1-e^{-\beta(m_{\max}-m_{\min})}},\quad m_{\min}\le m<m_{\max}.$$ (eq:trunc-gr)

在獨立 Poisson 超越事件的假設下，期間 $T$ 的超越機率為

$$P(N_x(T)\ge1)=1-\exp\left[-\int_0^T\nu_x(t)\,dt\right].$$ (eq:nhpp)

固定年率下的回歸期 $1/\nu_x$ 是平均時間尺度，不是等間隔日程。
例如年超越率 0.002 對應 500 年回歸期，50 年超越機率為 $1-e^{-0.1}$，約 9.5%。
未來群集會改變事件率時，需先處理未來歷史；F.6 說明為何不能直接套用零事件公式。

危害描述搖晃，風險還需要建物脆弱度、暴露量與損失。
本附錄的數值只是教學例子；工程設計須使用對應場址、法規和正式模型。

## F.1 存活函數、危害率與右設限

由 $S=1-F$、$S'=-f$，

$$h(t)=\frac{f(t)}{S(t)}=-\frac{d}{dt}\ln S(t),\qquad
S(t)=\exp\left[-\int_0^th(u)\,du\right].$$

所以已經等待 $T$ 後的下一次事件機率為

$$1-\frac{S(T+\Delta)}{S(T)}
=1-\exp\left[-\int_T^{T+\Delta}h(u)\,du\right].$$

這是一次間隔的條件存活關係，不要求整段更新過程具有獨立增量。事件發生後年齡會重設，不能繼續沿原來的 $h(u)$ 來描述第二次破裂。

若觀測從一次已知事件開始，完整間隔為 $x_1,\ldots,x_n$，末端還等待了 $c$ 而沒有事件，概似為

$$L(\theta)=\prod_{i=1}^nf_\theta(x_i)\,S_\theta(c).$$

最後一項是右設限資訊。若一開始已知年齡為 $a$，第一個等待 $u$ 的密度為 $f(a+u)/S(a)$；若年齡未知，就需要另外指定抽樣或先驗條件。平穩更新過程在任意時刻看到的剩餘等待，通常不等於原始間隔分布，這是長間隔較容易被抽到的結果。

## F.2 Weibull 與對數常態

Weibull 尺度 $\theta>0$、形狀 $k>0$：

$$F(t)=1-e^{-(t/\theta)^k},\quad
h(t)=\frac{k}{\theta}(t/\theta)^{k-1},\quad
\mathbb E[T]=\theta\Gamma(1+1/k).$$

由換元 $u=(t/\theta)^k$ 可得 $\mathbb E[T^r]=\theta^r\Gamma(1+r/k)$，因此

$$c_v^2=\frac{\Gamma(1+2/k)}{\Gamma(1+1/k)^2}-1.$$

對數常態若 $\ln T\sim N(\mu,\sigma^2)$，則 $\mathbb E[T]=e^{\mu+\sigma^2/2}$、$c_v^2=e^{\sigma^2}-1$。已知平均 $T_r$ 和 $c_v$ 時，$\sigma^2=\ln(1+c_v^2)$、$\mu=\ln T_r-\sigma^2/2$。

其危害函數是 $\phi(z)/[t\sigma\bar\Phi(z)]$，$z=(\ln t-\mu)/\sigma$。由 Mills 比 $\phi(z)/\bar\Phi(z)\sim z$，長時間 $h(t)\sim(\ln t-\mu)/(\sigma^2t)\to0$。這是尾部條件化的結果，不是已觀測到的斷層卸載機制。

## F.3 BPT 的首達時間

令載入 $X(t)=\rho t+\sigma W_t$，$X(0)=0$，門檻 $a>0$，$\rho,\sigma>0$；首次到達時間 $\tau=\inf\{t:X(t)=a\}$。無漂移情形，由反射原理

$$P(\tau\le t)=2\left[1-\Phi\left(\frac{a}{\sigma\sqrt t}\right)\right],$$

微分得到 $f_0(t)=a(\sigma\sqrt{2\pi t^3})^{-1}e^{-a^2/(2\sigma^2t)}$。

加入常漂移時，對已在 $t$ 首次碰到 $a$ 的路徑，漂移改變測度的因子為 $e^{\rho a/\sigma^2-\rho^2t/(2\sigma^2)}$。這一步使用布朗運動的指數鞅；乘回並配方得

$$\begin{aligned}
f(t)&=\frac{a}{\sigma\sqrt{2\pi t^3}}
\exp\left[-\frac{a^2-2\rho at+\rho^2t^2}{2\sigma^2t}\right]\\
&=\frac{a}{\sigma\sqrt{2\pi t^3}}
\exp\left[-\frac{(a-\rho t)^2}{2\sigma^2t}\right].
\end{aligned}$$

其 Laplace 變換為

$$\mathbb E[e^{-s\tau}]=\exp\left[\frac{a}{\sigma^2}
\left(\rho-\sqrt{\rho^2+2\sigma^2s}\right)\right].$$

在 $s=0$ 微分得 $\mathbb E[\tau]=a/\rho$、$\mathrm{Var}(\tau)=a\sigma^2/\rho^3$。設 $T_r=a/\rho$、$c_v^2=\sigma^2/(a\rho)$，即得到 {eq}`eq:bpt-pdf`。

常用 inverse Gaussian 形狀尺度參數為平均 $T_r$、形狀 $T_r/c_v^2$。SciPy 的 `invgauss(mu, scale)` 則用 `mu=c_v**2`、`scale=T_r/c_v**2`；這是參數化對照，不代表名稱相同的參數可跨軟體直接複製。

由密度可寫 $\ln f(t)=\mathrm{const}-(3/2)\ln t-t/(2c_v^2T_r)-T_r/(2c_v^2t)$。尾部 $S(t)\sim2c_v^2T_rf(t)$，故

$$\lim_{t\to\infty}h(t)=\frac1{2c_v^2T_r}.$$

離散步長模擬布朗路徑可能錯過步間首達，須做步長收斂或使用適當橋接修正；離散路徑只能提供近似。

## F.4 截斷規模分布與事件率

令 $d=m_{\max}-m_{\min}>0$，積分 $e^{-\beta(m-m_{\min})}$ 後正規化，即得 {eq}`eq:trunc-gr`。其累積分布為

$$F_M(m)=\frac{1-e^{-\beta(m-m_{\min})}}{1-e^{-\beta d}}.$$

取 $U\sim\mathrm{Unif}(0,1)$，反解得

$$M=m_{\min}-\beta^{-1}\ln[1-U(1-e^{-\beta d})].$$

分部積分可得

$$\mathbb E[M]=m_{\min}+\frac1\beta-\frac{d}{e^{\beta d}-1}.$$

$d\to\infty$ 時回到 $m_{\min}+1/\beta$。如果原發生率是門檻 $m_a$ 以上事件的年率 $\nu_a$，改成較高門檻 $m_b$ 時，應乘以原分布的 $P(M\ge m_b)$；不能只重設密度而保持率不變。

## F.5 全機率積分與反聚合

給定震源 $s$ 的規模距離聯合密度 $f_{M,R\mid s}$，一次事件超越門檻機率為

$$q_s(x)=\int\!\!\int P(IM>x\mid m,r,s)f_{M,R\mid s}(m,r)\,dm\,dr.$$

乘年率並加總即為 {eq}`eq:psha`。對數常態地動假設下，可再展開標準化殘差 $\varepsilon$：

$$\lambda_{IM}(x)=\sum_s\nu_s\int\!\!\int\!\!\int_{\varepsilon_0(x,m,r,s)}^\infty
\phi(e)f_{M,R\mid s}(m,r)\,de\,dm\,dr,$$

其中 $\varepsilon_0=(\ln x-\mu_{\ln IM})/\sigma_{\ln IM}$。規模距離若相依，保留聯合密度；若截斷殘差分布，須對其重新正規化。

令 $c_{jk}(x)$ 是規模距離格對此積分的非負貢獻。若 $\lambda_{IM}(x)>0$，反聚合權重 $w_{jk}=c_{jk}/\lambda_{IM}(x)$ 相加為一。它給超越事件的相對來源，不是未來第一場地震的預報。

## F.6 何時可以用 Poisson 零事件公式

假設源事件為 Poisson，且給定事件標記後，各次是否超越門檻依指定獨立標記機制產生。獨立 thinning 使超越事件也是 Poisson；多個獨立 Poisson 震源相加仍為 Poisson。

如果超越率為確定函數 $\nu_x(t)$，將時間切成小段，各段零事件機率為 $1-\nu_x(t_i)\Delta t+o(\Delta t)$。由獨立增量相乘、取對數並令步長趨零，

$$P(N_x(T)=0)=\exp\left[-\int_0^T\nu_x(t)\,dt\right],$$

即 {eq}`eq:nhpp`。對固定率，回歸期是 $1/\nu_x$，仍不表示等間隔。

若強度為外生隨機過程且條件下為 Poisson（Cox），則是 $P(N=0)=\mathbb E[e^{-A}]$，$A=\int\nu_x(t)dt$，而不是 $e^{-\mathbb E[A]}$。例如等機率選 $A=0$ 或 $2$，均值同為一，零事件機率卻是 $(1+e^{-2})/2$，不等於 $e^{-1}$。

ETAS 的強度由事件本身改變，不能直接當成外生 Cox 路徑使用這個條件混合公式。發布時應傳遞所有未來歷史，透過相應點過程方法或完整模擬求至少一次超越；若每條模擬還生成地動，須保留模型要求的標記相依性。

更新過程的下一次破裂可用 F.1 的條件存活公式，但某次未超越地動門檻的破裂仍會重設斷層年齡。故不能直接把更新 hazard 乘超越比例，再不加說明地套非齊次 Poisson 公式。

## F.7 載入與釋放：自我修正和應力釋放模型

一段斷層長期承受載入，事件又可能釋放部分累積量。這提供與Hawkes正向觸發不同的建模方向。**自我修正過程**（self-correcting process）的一個簡單形式是

$$\lambda^*(t)=\exp\{a+bt-cN(t^-)\},\qquad b>0,\ c>0.$$

時間單位固定後，$a$ 控制參考對數率，$b$ 控制無事件時的增加速度，$c$ 控制一次事件的降低幅度。$N(t^-)$ 是時刻 $t$ 前的累積事件數；此例無標記，每次事件具有相同降低幅度。式子保持率為正，並非將負的Hawkes核直接加到背景。

例如取 $c=\ln2$，每次事件之後的率乘 $1/2$。若 $b=\ln2/10$（每天），連續十天無事件，率就乘二。這兩個數字只是公式示意，沒有用義大利地震擬合；它們描述隨機發生率的變化，不代表十天後必然發生事件。

### 讓每次釋放量由標記決定

**應力釋放模型**（stress-release model）用一個隱含狀態表示載入及地震釋放。可考慮

$$X(t^-)=X_0+\rho t-\sum_{t_i<t}\Delta_i,\qquad
\lambda^*(t)=\lambda_{\rm ref}\exp\{\gamma X(t^-)\},$$

其中 $\rho>0$ 是載入率，$\Delta_i>0$ 是第 $i$ 次釋放量，$\gamma>0$ 控制狀態對率的影響；需使指數中的量無量綱。若 $\Delta_i$ 全取同一常數，便回到上面的計數自我修正形式。若 $\Delta_i$ 依事件規模決定，就須另定規模—釋放關係和標記分布。

$X$ 是模型狀態，沒有應力觀測與單位校準時，不能把它直接當成已量得的斷層應力。Ogata（1999）以這類模型說明條件強度能表示負向回饋，也指出將真實應力資料引入模型是另一項工作。不同stress-release版本不必使用相同指數連結或釋放律。

### 為什麼它不是一般的更新過程？

更新過程在每次指定事件後把年齡設為零，下一次間隔從同一分布重新開始。此處事件後的 $X$ 通常仍保留過去累積狀態，不會回到同一固定值，所以一般不具有iid等待間隔。只有另加適當重設條件，才能得到相應更新模型。

同一斷層段的大破裂可以採載入—釋放的負回饋；其周圍餘震可同時具有正向觸發。兩種敘事的事件集合與空間尺度不同，並不必然互相矛盾。實際混合兩種機制時，仍須指定完整強度並檢驗，不能只依物理語言選擇模型。

附錄A.11並列Poisson、Cox、renewal、Hawkes與self-correcting。這些模型的零事件機率都要依其條件歷史處理；本節的事件後率下降，不會自動使整段事件數服從Poisson。
## 參考資料與延伸閱讀

- Ogata（1999），[Seismicity Analysis through Point-process Modeling: A Review](https://doi.org/10.1007/s000240050275)。刊頁496討論應力狀態與資料，刊頁497–499並列更新、自我修正與Hawkes過程；用來對照不同歷史回饋。

- [Related Distributions：Hazard Function 與 Survival Function](https://itl.nist.gov/div898/handbook/eda/section3/eda362.htm) — NIST／SEMATECH，*e-Handbook of Statistical Methods*（免費）。先看機率密度、存活函數與危害函數的關係，釐清本章「已經等到現在，再發生的瞬時率」與一般發生機率的差別。

- [A Brownian model for recurrent earthquakes](https://doi.org/10.1785/0120010267) — Mark V. Matthews、William L. Ellsworth、Paul A. Reasenberg（2002），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[USGS 免費摘要](https://pubs.usgs.gov/publication/70024443)）。BPT 複發模型的核心論文，從帶有布朗擾動的載入過程推到首達時間，銜接本章的物理直覺與長時間危害率。

- [Modeling the earthquake occurrence with time-dependent processes: a brief review](https://doi.org/10.1007/s11600-019-00284-4) — Ourania Mangira、Christos Kourouklas、Dimitris Chorozoglou、Aggelos Iliopoulos、Eleftheria Papadimitriou（2019），*Acta Geophysica*（全文可能需訂閱）。將本章的更新過程、複發分布與應力釋放放回時間相依模型的全貌，閱讀時比較各模型適用的時間尺度、資料與假設。

- [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。

- [Introduction to Probabilistic Seismic Hazard Analysis](https://scits.stanford.edu/sites/g/files/sbiybj22081/files/media/file/baker_2013_intro_psha_v2_0.pdf) — Jack W. Baker（2013），White Paper Version 2.0（免費教材）。本章 PSHA 推導的入門來源，依五個步驟走過震源、規模、距離、地動與積分；接著讀反聚合及回歸期，理解危害曲線背後有哪些事件在貢獻。

- [TEM 計畫概要與研究成果](https://tem.tw/TEM2020/portfolio-overview.html) — Taiwan Earthquake Model 計畫團隊，TEM2020 專頁（免費）。先看孕震構造、測地資料與危害度評估如何分工，再循頁面列出的論文了解臺灣模型的輸入來源。

- [Probabilistic seismic hazard assessment for Taiwan: TEM PSHA2020](https://doi.org/10.1177/8755293020951587) — Chung-Han Chan 等（2020），*Earthquake Spectra*（出版社全文可能需訂閱；[中央大學免費全文](https://www.gep.ncu.edu.tw/storage/thesis/2020/2020%20Chung-Han%20Chan_ES2.pdf)）。本章臺灣案例的正式文獻，說明孕震構造、地震目錄、地動預估式與場址效應如何納入 PSHA；可對照教材檢查每一類不確定性出現在哪一步。
