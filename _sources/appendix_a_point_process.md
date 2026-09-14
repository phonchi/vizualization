# 附錄 A：記號、統計工具與點過程

**對應主文章**：
第 3 章 {doc}`03_experiment_spec`；
第 4 章 {doc}`04_poisson_and_sup`；
第 5 章 {doc}`05_likelihood_estimation`；
第 6 章 {doc}`06_simulation_tests_scores`；
第 11 章 {doc}`11_conditional_intensity`。

本附錄提供記號回查與工具推導。
主文先說用途，附錄再列適用條件。
先備知識為基礎機率、微分與積分。

## A.1 讀者版記號表

條件強度（conditional intensity）是歷史下的瞬時率。
率密度還須指定面積與規模單位。
格內期望數則已完成這些積分。

| 概念 | 記號 | 定義或單位 |
|---|---|---|
| 事件 | $(t_i,x_i,y_i,m_i)$ | 時間、位置與取整規模 |
| 時間 | $t_i$ | 自 1960-01-01 起算的天數 |
| 位置 | $(x_i,y_i)$ | 投影公里座標 |
| 歷史 | $H_t$ | $t$ 之前可取得的事件 |
| 條件強度 | $\lambda^*(t,x,y,m)$ | $\lambda(t,x,y,m\mid H_t)$ |
| 格內期望數 | $\Lambda_{jk}$ | 指定窗、格 $j$、規模箱 $k$ |
| 觀測數 | $\omega_{jk}$ | 同窗、同格箱的事件數 |
| 完整度 | $M_c$ | 目錄可充分記錄的規模下限 |
| 輸入門檻 | $m_0$ | 模型納入歷史事件的下限 |
| 目標門檻 | $m_T$ | 預報與檢驗的規模下限 |
| 機率 | $P(\cdot)$ | 一律以大寫 $P$ 表示 |

Gutenberg–Richter 律簡稱 GR 律。
它描述規模增加時的事件數衰減。
核函數（kernel）分配單個來源的貢獻。
產能（productivity）控制其總量。

| 概念 | 記號 | 定義或提醒 |
|---|---|---|
| GR 斜率 | $b,\beta$ | $\beta=b\ln10$ |
| GR 密度 | $s(m)$ | $\beta e^{-\beta(m-m_0)}$ |
| Omori 時間密度 | $g(\tau)$ | $(p-1)c^{-1}(1+\tau/c)^{-p}$ |
| 產能 | $\kappa(m)$ | $K e^{\alpha(m-m_0)}$ |
| 分支比 | $n$ | 每事件的平均直接後代數 |
| PPE 參數 | $a,d,s$ | 權重、平滑距離與底率 |
| PPE 空間項 | $h_0(x,y)$ | 過去震央的平滑貢獻 |
| EEPAS 規模核 | $a_M,b_M,\sigma_M$ | 截距、斜率與散布 |
| EEPAS 時間核 | $a_T,b_T,\sigma_T$ | 對數時間的截距、斜率與散布 |
| EEPAS 空間核 | $b_A,\sigma_A$ | 規模斜率與長度尺度 |
| EEPAS 修正 | $\eta(m),\Delta(m)$ | 產能正規化與門檻保留比例 |
| EEPAS 混合權重 | $\mu_E$ | 背景成分的係數 |
| Ψ 尺度 | $M_P,T_P,A_P,M_m$ | 前兆規模、時間、面積與主震規模 |
| 累積規模異常 | $C(t)$ | 去除全窗平均趨勢的累積量 |

PPE 是 Proximity to Past Earthquakes。
它以過去震央的鄰近程度分配率。
EEPAS 的全名如下。
Every Earthquake a Precursor According to Scale。
模型讓每個事件依規模貢獻前兆訊息。
Ψ 指前兆尺度增加現象。
英文為 precursory scale increase。
其量測定義見 {doc}`14_psi_precursory_scale`。

| 概念 | 記號 | 定義或提醒 |
|---|---|---|
| 逐箱 Poisson 對數概似 | POLL | Poisson log-likelihood |
| 聯合 Poisson 對數概似 | jPOLL | joint Poisson log-likelihood |
| 二元對數概似 | BILL | binary log-likelihood |
| 數量檢驗尾機率 | $\delta_1,\delta_2$ | 含等號的上尾與下尾 |
| 分位數分數 | $q$ | 模擬統計量不大於觀測的比例 |
| 每事件資訊增益 | IGPE | information gain per earthquake |
| BPT 平均複發時間 | $T_r$ | 與時間軸同單位 |
| 間隔變異係數 | $c_v$ | 標準差除以平均間隔 |
| Weibull 形狀 | $k$ | 不使用 $\beta$ |
| Janus 組合權重 | $\pi_{\rm ETAS}$ | ETAS 成分的權重 |

BPT 指 Brownian passage time。
它以隨機載入的首達時間描述複發。
Janus 指結合不同時間尺度的組合模型。
相關推導見 {doc}`appendix_f_hazard`。
組合方法見 {doc}`19_ensembles`。

$\beta$ 只表示 $b\ln10$。
$p$ 用於 Omori 指數；檢定時寫「p 值」。
ETAS 空間尾指數也沿用文獻的 $q$。
它與分位數分數須依上下文區分。
表中的 $g$ 要求 $p>1,c>0$。
不同核正規化下的 $K$ 不可直接互換。
具體換算見 {doc}`appendix_c_etas`。

## A.2 兩套術語對照

CSEP 是地震可預報性協作研究計畫。
英文為 Collaboratory for the Study of Earthquake Predictability。
本站採其預報實驗用語。
ETAS 是類流行病餘震序列模型。
英文為 Epidemic-Type Aftershock Sequence。
下表對照 Jalilian（2019）的 R 套件用語。

| 本站用語 | 英文 | ETAS R 對照 |
|---|---|---|
| 測試區 $R$ | testing region | study region |
| 收集區 $S$ | collection／neighbourhood region | complementary events 所在區 |
| 學習期、測試期 | learning／testing period | 依該次分析指定 study period |
| 目標地震 | target events | target events |
| 提供歷史的事件 | precursor／trigger | complementary events 的對應角色 |
| 十進位天數 | decimal days，自 1960-01-01 | 自 `time.begin` 起算 |

收集區提供可能影響測試區的事件。
目標事件才進入指定評估窗的求和。
目標事件發生後，也可成為後續歷史。
「前兆」與「觸發源」仍屬不同模型概念。
術語對照不代表已辨認物理因果。

區域邊界稱多邊形（polygon）。
ETAS R 要求頂點逆時針排列。
EEPAS 框架要求順時針排列。
兩者的時間原點也必須明確換算。

## A.3 Poisson 過程與計數

Poisson 過程以独立增量描述事件。
此處的率 $\lambda(t)$ 是確定函數。
令 $N(B)$ 為時間區間 $B$ 的事件數。
將區間切細，每段多事件機率可忽略。
各段獨立的 Bernoulli 計數取極限，得到

$$\begin{aligned}
\Lambda(B)&=\int_B\lambda(t)\,dt,\\
P(N(B)=n)&=e^{-\Lambda(B)}\frac{\Lambda(B)^n}{n!},\\
\mathbb E[N(B)]&=\operatorname{Var}(N(B))=\Lambda(B),\\
P(N(B)\ge1)&=1-e^{-\Lambda(B)}.
\end{aligned}$$

齊次（homogeneous）表示率固定。
此時長度 $T$ 的期望數為 $\lambda T$。
非齊次表示率可隨時間改變。
非齊次本身不代表事件互相觸發。

過度離散（overdispersion）指變異數大於平均。
叢集可造成這種現象。
混合不同年份的率也可能如此。
不能只憑計數變異就判定觸發機制。

## A.4 概似、MLE 與 bootstrap

概似（likelihood）比較参数對資料的支持。
它是固定資料後的參數函數。
最大概似估計簡稱 MLE。
英文為 maximum likelihood estimation。

固定觀察期長 $T$，觀測數為 $N$。
齊次 Poisson 的對數概似為

$$\begin{aligned}
\ell(\lambda)&=N\ln(\lambda T)-\lambda T-\ln(N!),\\
\ell'(\lambda)&=N/\lambda-T,\\
\widehat\lambda&=N/T.
\end{aligned}$$

$N>0$ 時，二階導數為負。
$N=0$ 時，最大值位於零率邊界。
資訊量（information）描述概似的曲率。
其期望值為 $I(\lambda)=T/\lambda$。
標準誤（standard error）近似為

$$\operatorname{SE}(\widehat\lambda)\simeq\frac{\sqrt N}{T}.$$

零事件時，此近似不能表示沒有不確定性。
應回到完整計數分布建立區間。
概似歸一化也不會自動變成參數機率。

Bootstrap 是重抽樣估計法。
它用重抽的資料近似估計量變動。
普通版本從原樣本有放回抽樣。
每份樣本重新估參數，再取經驗分位數。
參數式版本則從擬合模型產生資料。

地震序列常有相依性。
逐事件重抽可能破壞叢集結構。
區塊重抽樣保留一段內的相依。
模型式模擬則依賴生成模型正確。
所報區間須註明採用哪種抽樣單位。

## A.5 分位數、p 值與對數分數

統計量（statistic）是資料的一個摘要。
令觀測值為 $s_{\rm obs}$。
虛無模型產生 $B$ 個模擬值 $s_r$。
分位數分數估計其下尾位置：

$$\widehat q=\frac1B\sum_{r=1}^B\mathbf1\{s_r\le s_{\rm obs}\}.$$

若低分才是不相容證據，便採下尾檢定。
連續分布下，其 p 值為 $q$。
若高分才是證據，則採上尾。
離散分布須保留等號機率。

$$\begin{aligned}
P(S\ge s)&=1-F_S(s^-),\\
P(S\le s)&=F_S(s),\\
p_{\rm two}&=\min\{1,2\min[P(S\le s),P(S\ge s)]\}.
\end{aligned}$$

最後一式是等尾雙尾慣例。
它不是每種檢定唯一的雙尾定義。
離散資料不可一律使用 $1-q$ 作上尾。
模擬零次落入尾端，也不表示機率為零。

對數分數（log score）獎勵觀測處的機率。
給觀測 $y$ 的分數為 $\ln P(y)$。
連續觀測改用相同測度下的密度。
兩模型相減得到對數概似比。
網格、單位與目標不同時，不可直接比較。

## A.6 AIC、KL 與貝氏推論

KL 散度比較兩個分布的差異。
全名為 Kullback–Leibler divergence。
令真實分布為 $P$，候選分布為 $Q$。
在離散且支撐相容的情況下，

$$\begin{aligned}
D_{\rm KL}(P\Vert Q)&=\sum_zP(z)\ln\frac{P(z)}{Q(z)}\ge0,\\
\mathbb E_P[\ln Q(Z)]&=\mathbb E_P[\ln P(Z)]-D_{\rm KL}(P\Vert Q).
\end{aligned}$$

提高期望對數分數，相當於減少 KL 損失。
有限測試資料只能估計這個期望。
它不能直接揭露未知真實分布。

AIC 是赤池資訊準則。
英文為 Akaike information criterion。
令自由參數數目為 $d$，則

$$\operatorname{AIC}=-2\ell(\widehat\theta)+2d.$$

懲罰項近似修正訓練分數的樂觀偏差。
較小的 AIC 較受此準則支持。
推導需要正則性與大樣本近似。
相依資料、邊界參數須另查條件。
AIC 權重也不是自動得到的後驗機率。

貝氏推論（Bayesian inference）結合先驗與資料。
先驗密度記為 $\pi(\theta)$。
後驗密度記為 $\pi(\theta\mid\mathcal D)$。

$$\pi(\theta\mid\mathcal D)=
\frac{L(\theta;\mathcal D)\pi(\theta)}
{\int L(u;\mathcal D)\pi(u)\,du}.$$

以形狀、率參數定義 Gamma 先驗。
令 $\lambda\sim\operatorname{Gamma}(a,r)$。
Poisson 資料給出

$$\lambda\mid N\sim\operatorname{Gamma}(a+N,r+T).$$

這是共軛（conjugate）更新。
意思是後驗仍屬同一分布族。
未來長度 $h$ 的計數變異數為

$$\operatorname{Var}(N_{\rm future}\mid\mathcal D)
=h\mathbb E[\lambda\mid\mathcal D]
+h^2\operatorname{Var}(\lambda\mid\mathcal D).$$

第一項來自事件本身的波動。
第二項來自未知參數。
ETAS 的後驗預報還須整合未來觸發。
模型未描述的漏測也不會自動消失。

## A.7 下一事件與點過程概似

存活函數（survival function）表示尚未發生。
令 $S_0(t)=P(T_1>t\mid H_{t_0})$。
沿沒有新事件的路徑，危害率記為 $\lambda_0(t)$。
小時間段的乘法關係給出

$$\begin{aligned}
S_0(t+h)&=S_0(t)[1-\lambda_0(t)h+o(h)],\\
S'_0(t)&=-\lambda_0(t)S_0(t),\\
S_0(t)&=\exp\left[-\int_{t_0}^t\lambda_0(u)\,du\right],\\
f_0(t)&=\lambda_0(t)S_0(t).
\end{aligned}$$

此式要求相應密度存在且率局部可積。
若還有未知外部狀態，須一併平均。
不能任取一條外部路徑代入。

連乘每次等待密度，再乘最後空白，得到

$$\begin{aligned}
L(\theta)&=\prod_{i=1}^N\lambda^*_\theta(t_i)
\exp\left[-\int_0^T\lambda^*_\theta(u)\,du\right],\\
\ell(\theta)&=\sum_{i=1}^N\ln\lambda^*_\theta(t_i)
-\int_0^T\lambda^*_\theta(u)\,du.
\end{aligned}$$

這是事件時間的密度概似。
精確實數時間本身沒有正機率。
過程須簡單、非爆炸且符合可積條件。
非爆炸指有限期間不產生無限事件。

令評估窗為 $[a,b]\times R\times[m_T,m_u)$。
加上位置與規模後，

$$\ell(\theta)=\sum_{i\in\mathcal I}\ln\lambda^*_\theta(t_i,x_i,y_i,m_i)
-\int_a^b\int_R\int_{m_T}^{m_u}\lambda^*_\theta(t,x,y,m)\,dm\,dx\,dy\,dt.$$

$\mathcal I$ 只含評估窗內的目標事件。
$S$ 內及窗前事件可提供歷史。
估計輸入事件模型時，須另換評估門檻。
求和與積分必須使用同一觀察範圍。

若強度分解為 $\lambda^*_{ST}s(m)$，規模項可分開。
分開估參數還要求參數空間獨立。
例如 $\alpha=\beta$ 會把兩項重新連結。
當次規模獨立，也不代表歷史規模無作用。

## A.8 補償子、殘差與時間變換

補償子（compensator）是沿歷史累積的率。
本節只看時間，定義 $A(t)=\int_0^t\lambda^*(u)\,du$。
鞅（martingale）的未來增量條件期望為零。
適當条件下，$N(t)-A(t)$ 是鞅。
因此

$$\begin{aligned}
\mathbb E[N(t)]&=\mathbb E[A(t)],\\
\mathbb E[N(T)-N(t_0)\mid H_{t_0}]
&=\mathbb E\left[\int_{t_0}^T\lambda^*(u)\,du\mid H_{t_0}\right].
\end{aligned}$$

外側期望包含所有可能的未來歷史。
目前歷史凍結後的積分通常只是近似。
尤其小事件可能先出現，再觸發目標事件。

殘差（residual）比較觀測與模型預期。
可預測權重 $v$ 只使用事件前的資訊。
相應加權殘差為

$$R_v=\sum_i v(t_i)-\int v(t)\lambda^*(t)\,dt.$$

適當可積條件下，$\mathbb E[R_v]=0$。
倒數率權重會放大低率處的誤差。
事後家族標籤不可當成可預測資訊。

時間變換（time rescaling）改用累積率計時。
若 $A$ 連續且終將趨於無限，令 $\tau_i=A(t_i)$。
正確強度下，

$$\begin{aligned}
P(\tau_{i+1}-\tau_i>s\mid H_{t_i})&=e^{-s},\\
U_i&=1-e^{-(\tau_i-\tau_{i-1})}.
\end{aligned}$$

變換間隔應獨立且服從單位指數分布。
$U_i$ 應獨立且服從均勻分布。
率為零的平臺須用廣義反函數處理。
有限窗另有末端截尾。

Q–Q 圖對照兩分布的相同分位數。
KS 統計量是累積分布的最大距離。
兩者只檢查邊際形狀仍不夠。
間隔相依性也須檢查。
同資料估參數後，應模擬並重新擬合。
如此才能校準診斷的參考分布。

## A.9 稀疏化與分支模擬

稀疏化（thinning）先提候選，再決定保留。
它需要候選區間內有效的上界 $\bar\lambda$。

1. 以 $\bar\lambda$ 產生候選等待時間。
2. 以 $\lambda^*(s)/\bar\lambda$ 保留候選。
3. 接受後更新歷史與上界。
4. 拒絕後只推進時間。
5. 到觀察窗末端時停止。

比值大於一表示上界失效。
直接截成一不能修復抽樣機制。
固定背景與遞減 Omori 核可用當下率作界。
延遲峰值或上升背景須另找上界。

分支模擬（branching simulation）先抽背景。
每個事件再生成直接後代。
後代仍能繼續繁衍。
它需要模型具有相應 Poisson 分支表示。

Omori 核的反函數抽樣為

$$\tau=c[(1-U)^{-1/(p-1)}-1],\qquad U\sim\operatorname{Unif}(0,1).$$

這裡要求 $p>1$。
對面積尺度 $v>0$ 的空間核，

$$\begin{aligned}
h(r)&=\frac{q-1}{\pi v}(1+r^2/v)^{-q},\qquad q>1,\\
P(R\le r)&=1-(1+r^2/v)^{1-q},\\
R&=\sqrt{v[(1-U)^{-1/(q-1)}-1]}.
\end{aligned}$$

角度均勻抽於 $[0,2\pi)$。
推導使用極座標面積元 $r\,dr\,d\vartheta$。
$\sqrt v$ 才是長度尺度。
長尾樣本不能直接壓到窗邊。
區外事件若能觸發區內後代，仍須保留。
空歷史起步也不等於平穩目錄。

## A.10 固定總數與事件對

固定 $N$ 個均勻時間點，再排序。
包含兩端的間隔共有 $N+1$ 個。
任一間隔 $W$ 滿足

$$P(W>w\mid N)=\left(1-\frac wT\right)^N,\qquad 0\le w\le T.$$

平均間隔為 $T/(N+1)$。
這些間隔互相相依。
它們不是未固定總數的獨立指數間隔。

時間切成 $J$ 個等長箱後，計數服從多項分布。
若樣本變異數分母為 $J-1$，則

$$\mathbb E[s_N^2]=N/J.$$

Fano 比（Fano factor）是變異數除以平均。
上述固定總數設定的期望基準仍為一。
單次比值偏離一，不等於顯著叢集。

一階矩描述平均計數。
二階階乘矩描述不同事件對。
對區間 $B,C$，

$$M^{[2]}(B\times C)=\mathbb E\left[
\sum_{i\ne j}\mathbf1\{t_i\in B\}\mathbf1\{t_j\in C\}\right].$$

排除自身配對後，才能比較事件對結構。
相同平均率仍可有不同叢集程度。
非均勻背景與邊界也會改變近鄰數。

## A.11 依條件強度辨認模型家族

第11章以Hawkes作為事件回饋的主例。下表比較的是生成機制，而不是單看圖形是否成群。非齊次Poisson的高率時段、Cox的共同環境與Hawkes的觸發，都可能產生視覺上密集的活動。

| 模型 | 率使用哪些資訊？ | 事件出現後的作用 |
|---|---|---|
| 齊次／非齊次 Poisson | 確定常數或確定時間函數 $\lambda(t)$ | 不改變已指定的率 |
| Cox | 外生隨機強度路徑；給定路徑為 Poisson | 更新對未知路徑的推測，不是事件回饋生成該路徑 |
| 更新過程（renewal） | 自上次指定事件的年齡 $a$；率為 $f(a)/S(a)$ | 年齡歸零；率如何跳變由間隔分布決定 |
| 線性 Hawkes | 背景加所有已發生事件的非負因果核 | 加入新的未來觸發貢獻 |
| 自我修正（self-correcting） | 載入與累積事件造成的負向回饋 | 在典型形式中降低率，平靜時再逐漸上升 |

這五類不必構成互斥且窮盡的分類。齊次Poisson也是指數間隔的更新過程；更一般模型還能混合外生變動與事件回饋。比較表用來找主要假設，不作所有隨機過程的分類定理。

### 標準線性 Hawkes 的積分形式

以 $N(du)$ 表示時間 $u$ 的事件計數測度。有限起點模型可寫成

$$\lambda^*(t)=\mu+\int_{[0,t)}\phi(t-u)\,N(du)
=\mu+\sum_{0\le t_i<t}\phi(t-t_i).$$

若有更早歷史，另加入那些事件的核貢獻；平穩版本則考慮無限過去。$\phi$ 取非負、因果且可積的函數，$n=\int_0^\infty\phi(u)\,du$ 是未標記模型的平均直接後代數。對標準固定正背景模型，$n<1$ 給出有限均值的平穩版本。

「線性」是對事件計數測度的可加作用，並不要求 $\phi(u)$ 對時間線性。非線性連結、負核抑制或多類型互激發需要另立條件，不能照搬這個單一非負核的分支論證。第12章及附錄C已給出規模依賴產能後的對應平均。

### 從標記強度到 ground process

**標記**（mark）是事件附帶量，例如規模。**基礎過程**（ground process）保留時間與位置，去掉標記。令 $z=(t,x,y)$，帶標記的強度可分解為

$$\lambda^*(z,m)=\lambda_g^*(z)\,f(m\mid z,H_t),\qquad
\lambda_g^*(z)=\int \lambda^*(z,m)\,dm.$$

此處 $\lambda_g^*$ 仍以可用的完整歷史為條件，可以依賴過去規模；「去掉標記」不是把歷史規模的資訊一律丟棄。若只觀測ground process，就要改成相應較少資訊下的條件強度。

**不可由過去預知的標記**（unpredictable marks）描述新標記的分布不依過去事件及其標記。標準ETAS用固定 $s(m)$ 生成新規模，同時讓過去 $m_i$ 決定產能和空間尺度。附錄A.7的概似分解正是這個設定的例子。

**獨立標記**（independent marks）是另一個較強、不同的條件：給定全部事件時間與位置後，各標記互相獨立，分布只依各自位置時間。不能從標準ETAS的「新規模按固定GR抽樣」直接推出此條件，因為稍後的事件位置與數量含有先前規模造成的訊息。

例如，一顆規模較大的事件通常容許更多後代。看到它之後一整團事件，會改變對其規模的條件推測；這不違反事件生成當下的GR抽樣假設。Reinhart（2018）第2.3節的區分，正是避免把這兩種獨立混用。
## 參考資料與延伸閱讀

- Gallager（2011），[Poisson Processes](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/resources/mit6_262s11_chap02/)。免費講義。先讀等待時間與計數的關係。
- Reinhart（2018），[自激發時空點過程綜述](https://doi.org/10.1214/17-STS629)。[免費作者稿](https://arxiv.org/abs/1708.02647)。串起估計、分支與殘差。
- Daley 與 Vere-Jones（2003），[點過程理論，卷一](https://doi.org/10.1007/b97277)。需館藏或訂閱。回查事件對、標記與時間變換。
- Ogata（1999），[地震點過程建模綜述](https://doi.org/10.1007/s000240050275)。全文可能需訂閱。對照地震概似與模擬。
- Jalilian（2019），[ETAS R 套件](https://doi.org/10.18637/jss.v088.c01)。免費全文。核對兩區域與歷史事件用語。
- Naylor 等（2023），[ETAS.inlabru 貝氏建模](https://doi.org/10.3389/fams.2023.1126759)。免費全文。比較參數不確定性與事件波動。
- Biondini 等（2023），[EEPAS 在義大利的應用](https://doi.org/10.1093/gji/ggad123)。[免費機構版本](https://www.earth-prints.org/handle/2122/17084)。把記號對回義大利實驗。
- Savran 等（2022），[pyCSEP 工具介紹](https://doi.org/10.21105/joss.03658)。免費全文。連結預報格式與統計檢驗。
