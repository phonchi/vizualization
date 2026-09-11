# 附錄 A：點過程的機率、概似與模擬

這裡接續{doc}`10_point_process`。正文先建立率如何隨歷史更新的直覺，本附錄
補上機率與概似的推導、時間變換的條件，以及生成圖表所用的抽樣原理。
只想沿主文學習，可以回到{doc}`11_catalog_completeness_b`。

## 從下一事件的危害率得到存活函式

給定時間 $t_0$ 的歷史 $H_{t_0}$，令 $T_1$ 是下一事件時間，
$S_0(t)=P(T_1>t\mid H_{t_0})$。在直到 $t$ 都沒有新事件的條件下，
令 $\lambda_0(t)$ 為相應的危害率。對小量 $h>0$，

$$S_0(t+h)=S_0(t)[1-\lambda_0(t)h+o(h)].$$

假設所需的絕對連續性與局部可積性成立，取極限得到
$S'_0(t)=-\lambda_0(t)S_0(t)$。以 $S_0(t_0)=1$ 積分：

$$\begin{aligned}
\frac{\mathrm d}{\mathrm dt}\log S_0(t)&=-\lambda_0(t),\\
S_0(t)&=\exp\left[-\int_{t_0}^t\lambda_0(u)\,\mathrm du\right],\\
f_0(t)&=\lambda_0(t)S_0(t).
\end{aligned}$$

對只含可觀測事件歷史與已知協變數的模型，$\lambda_0$ 可由「無新事件」的
路徑求得。含未觀測外部隨機狀態時，必須使用相對於所給資訊的正確危害率，
不能將某一條任意外部狀態路徑代入，當成條件於較少資訊的存活機率。

這是下一個任意事件的結果。若目標是 $m\ge m_T$，其他規模的事件仍可
先發生並影響後續目標事件，因此一般不能將目前歷史凍結、乘上一個規模比例
後就當成完整的目標預報。在點過程模型下，常以給定目前歷史的未來目錄模擬
計算目標事件出現比例。這個問題也不同於把一條已觀測的補償子稱作未來期望數。

## 將事件密度與最後的空白連乘

對 $0<t_1<\cdots<t_N\le T$，條件於觀察開始前的歷史，使用每次等待的
條件密度，最後再乘上 $(t_N,T]$ 無事件的條件機率。沿已觀測歷史評估強度後，

$$\begin{aligned}
L(\theta)&=\left\{\prod_{i=1}^N\lambda^*_\theta(t_i)
\exp\left[-\int_{t_{i-1}}^{t_i}\lambda^*_\theta(u)\,\mathrm du\right]\right\}
\exp\left[-\int_{t_N}^{T}\lambda^*_\theta(u)\,\mathrm du\right]\\
&=\left\{\prod_{i=1}^N\lambda^*_\theta(t_i)\right\}
\exp\left[-\int_0^T\lambda^*_\theta(u)\,\mathrm du\right].
\end{aligned}$$

取對數即得到正文的{eq}`eq:pp-loglik`。這是相對於事件時間的密度形式，
不是某一串精確實數時間發生的正機率。它適用於具有相應條件密度、強度局部
可積且不爆發的常見簡單點過程；不能不加條件地推廣到所有隨機點集合。

空間與規模加入後，觀察區為 $S$、輸入規模下限為 $m_0$：

$$\ell(\theta)=\sum_i\log\lambda^*_\theta(t_i,x_i,y_i,m_i)
-\int_0^T\int_S\int_{m_0}^{\infty}\lambda^*_\theta(t,x,y,m)
\,\mathrm dm\,\mathrm dx\,\mathrm dy\,\mathrm dt.$$

求和只包含評估窗內事件，歷史則可能需要窗前或區域外的補充事件。
比較概似時須固定資料、觀察範圍、規模門檻與座標測度。每事件正規化不會
自動消除不同單位、網格、目標或背景選擇的差異。

## 標記分解不等於所有參數都能分開估

若 $\lambda^*(t,x,y,m)=\lambda^*_{\theta}(t,x,y)s_\beta(m)$，且
$s_\beta$ 在規模範圍積分為一，則

$$\ell=\sum_i\log s_\beta(m_i)+\sum_i\log\lambda^*_{\theta}(t_i,x_i,y_i)
-\int_0^T\int_S\lambda^*_{\theta}(t,x,y)\,\mathrm dx\,\mathrm dy\,\mathrm dt.$$

要分別最佳化，還需要兩部分的參數空間彼此獨立。例如將觸發產能指數設為
$\alpha=\beta$，會使時空項也依賴 $\beta$，此時不能只最大化規模項。
當次規模按固定密度抽取，也不表示過去規模不影響未來時間或位置。

## 補償子、殘差與時間變換

定義 $\Lambda(t)=\int_0^t\lambda^*(u)\,\mathrm du$。在常見可積性條件下，
$N(t)-\Lambda(t)$ 是鞅，因而 $E[N(t)]=E[\Lambda(t)]$。
對歷史相依模型，$\Lambda(t)$ 通常隨觀測歷史變化。給定 $H_{t_0}$ 的未來
期望數則滿足

$$E[N(T)-N(t_0)\mid H_{t_0}]
=E\left[\int_{t_0}^T\lambda^*(u)\,\mathrm du\mid H_{t_0}\right],$$

不能省略外側對未來歷史的期望。

對可預測權重 $h$（只依賴事件之前可取得的資訊），在適當可積性下，

$$E\left[\sum_i h(t_i,x_i,y_i)\right]
=E\left[\int h(t,x,y)\lambda^*(t,x,y)\,\mathrm dt\,\mathrm dx\,\mathrm dy\right].$$

因此兩者之差形成零平均的加權殘差。$h=1$、$1/\lambda^*$ 或
$1/\sqrt{\lambda^*}$ 強調不同區域，但率接近零時後兩者可能數值不穩定，
需要檢查可積性與估計誤差。用事後知道的家族分類當作可預測權重，不符合上述條件。

當累積強度連續、可作相應時間變換且終將趨於無限時，令
$\tau_i=\Lambda(t_i)$。下一事件存活公式給出

$$P(\tau_{i+1}-\tau_i>s\mid H_{t_i})=e^{-s}.$$

條件分布不依賴歷史，逐次套用可得獨立單位指數間隔。若強度有零值區段，
應使用廣義反函式處理平臺；若只觀察有限期間，則有邊界與截尾效應。
再由機率積分變換，得到正文{eq}`eq:time-rescale`。
用同一資料估參數後，不能假裝強度是事先已知的真值；可採參數模擬並重新擬合，
校準診斷統計量的參考分布。

## 稀疏化與分支模擬

稀疏化需要每一步都有對候選區間有效的上界 $\bar\lambda$。以該率提出
候選時間 $s$，再以 $\lambda^*(s)/\bar\lambda$ 接受。接受後更新事件歷史，
拒絕則不加入歷史。若上界不再有效，必須縮短區間或更新上界；將大於一的比值
直接截為一，不能修復錯誤的提案機制。

正文的正規化Omori核在事件之間遞減，背景率又固定，因此事件後的當下強度
是直到下一事件之前的有效上界。若核有延遲峰值，或背景率會上升，就需要
另找對候選區間有效的界，不能直接複製正文的選法。

具有相應Poisson分支表示的線性自激發模型，可以先抽背景事件，再抽每個事件的
後代數、延遲及位置。初始歷史、空間邊界與窗外後代處理必須和預定模型一致。
正文的兩種方法比較使用空歷史；若目標是平穩目錄，需處理觀察前歷史，而不是
直接假定空歷史等同穩態。

對無限時間可正規化的Omori核，$p>1,c>0$ 時，

$$g(u)=\frac{p-1}{c}(1+u/c)^{-p},\qquad
G(u)=1-(1+u/c)^{1-p}.$$

令 $U\sim U(0,1)$，反解 $G(u)=U$ 得
$u=c[(1-U)^{-1/(p-1)}-1]$。這個 $p>1$ 的限制來自無限時間正規化；
有限窗的未正規化衰減率可以有 $p\le1$，見{doc}`appendix_b_catalog`。

對平面上各向同性的核
$f(x,y)=(q-1)(\pi\sigma)^{-1}(1+r^2/\sigma)^{-q}$，$q>1,\sigma>0$，
極座標面積元素是 $r\,\mathrm dr\,\mathrm d\vartheta$，因此

$$P(R\le r)=\int_0^r2\pi u f(u)\,\mathrm du
=1-(1+r^2/\sigma)^{1-q}.$$

半徑可抽為 $R=\sqrt{\sigma[(1-U)^{-1/(q-1)}-1]}$，角度均勻抽於
$[0,2\pi)$。$\sigma$ 是面積尺度，$\sqrt{\sigma}$ 才是長度尺度。
計算時需避免端點造成溢位，並按預定觀察窗處理長尾樣本；直接截短延遲再保留
為窗內事件，會改變原模型。

## 描述平均與描述事件對

### 固定事件總數後的間隔與計數

{doc}`foundation_randomness`的均勻序列是將 $N$ 個獨立均勻點放在 $[0,T]$
後排序。包含兩端在內的 $N+1$ 個間隔可交換；任一間隔 $W$ 的存活機率為

$$P(W>w\mid N)=\left(1-\frac{w}{T}\right)^N,\qquad 0\le w\le T.$$

以內部某個間隔為例，要求它大於 $w$，再把它右側所有點向左平移 $w$，
可將允許區域對應到長度 $T-w$ 的有序點區域。$N$ 維體積比就是上式。
因此平均間隔為 $T/(N+1)$；指數參考的平均 $T/N$ 是未固定總數的另一種設定。
正文只取內部完整間隔，但其邊際分布仍是上式，指數不會因取了 $N-1$ 個間隔
而改成 $N-1$。這些間隔也不是互相獨立的。

若把時間平均切成 $k$ 箱，固定總數下的箱計數服從多項分布，各箱平均為 $N/k$。
以分母 $k-1$ 計算樣本變異數 $s^2$ 時，

$$E[s^2]=\frac{k(N/k)(1-1/k)}{k-1}=\frac Nk,
\qquad E\!\left[\frac{s^2}{N/k}\right]=1.$$

所以固定總數不改變這個Fano統計量的期望基準，但單一實現與少量箱數仍會
造成波動；這個期望等式本身不是顯著性檢定。

### 一階與二階資訊

對時間區間 $B$，一階矩測度 $M(B)=E[N(B)]$ 描述平均計數；
二階階乘矩測度描述不同事件對，例如

$$M^{[2]}(B\times C)=E\!\left[\sum_{i\ne j}
\mathbf1_{\{t_i\in B\}}\mathbf1_{\{t_j\in C\}}\right].$$

排除 $i=j$ 可避免把每個事件與自身配對。兩份目錄可有相同的一階結構，
但不同的事件對結構。非均勻背景與邊界會影響比較，因此簡單比較近鄰數量
不會自動辨認物理觸發。這是{doc}`foundation_randomness`二階統計的形式化版本。

## 符號回查

| 符號 | 意義 | 使用時留意 |
|---|---|---|
| $H_t$ | $t$ 以前的事件歷史 | 不含未來事件 |
| $\lambda^*$ | 條件率或率密度 | 交代時間、空間與規模單位 |
| $\Lambda(t)$ | 沿歷史累積的條件率 | 一般是隨機補償子 |
| $\Lambda_{jk}$ | 預報網格期望數 | 依已宣告預報條件定義 |
| $\omega_{jk}$ | 網格觀測數 | 和預報使用同樣目標 |
| $M_c,m_0,m_T$ | 完整度、輸入門檻、目標門檻 | 三者目的不同 |
| $b,\beta=b\ln10$ | 規模分布斜率 | $\beta$ 不作其他模型的形狀參數 |
| $\mu(x,y)$ | ETAS背景率 | 是模型成分，不等同物理分類 |
| $A,K$ | 正規化／未正規化核的產能係數 | 不直接跨參數化比較 |

回到{doc}`10_point_process`，或接著讀{doc}`11_catalog_completeness_b`。

## 參考資料與延伸閱讀

- Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I，第二版](https://doi.org/10.1007/b97277)。第5.4節支撐事件對，第6.3–6.4節處理群集與標記，第7.2–7.6節處理本附錄的強度、概似、時間變換和模擬；全文需訂閱或館藏權限。
- Ogata, Y.（1999），[Seismicity Analysis through Point-process Modeling: A Review](https://doi.org/10.1007/s000240050275)。附錄A整理地震點過程的概似、診斷與模擬，可與一般教材對照；出版社全文可能需訂閱。
- Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)；[免費作者預印本](https://arxiv.org/abs/1708.02647)。第2–3節說明分支表示、估計與殘差工具，適合作為教材與應用之間的橋梁。
- Jalilian, A.（2019），[ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01)，[免費全文](https://www.jstatsoft.org/article/view/v088c01)。對照時間核、空間核及研究區邊界的具體定義。
