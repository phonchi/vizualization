
# 附錄 D：Ψ 的量測與 EEPAS 的核函數

本附錄需要條件機率、簡單線性迴歸與一變數積分。它補足 {doc}`Ψ <14_psi_precursory_scale>`及 {doc}`EEPAS <15_eepas_italy_forecast>`的計算；讀完每節可回到對應主文繼續閱讀。

## D.1 累積規模異常的性質與限制

令 $q_i=M_i-(m_c-0.1)>0$，$A(t)=\sum_{t_s<t_i\le t}q_i$，$D=t_f-t_s$。由第 14 章累積規模異常的定義，$k=A(t_f)/D$，故

$$C(t_s)=0,\qquad C(t_f)=A(t_f)-\frac{A(t_f)}D D=0.$$

事件間斜率為 $-k$，在事件處跳升 $q_i$。取內部切點 $u$，記切點前後的平均累積速率為

$$r_-(u)=\frac{A(u)}{u-t_s},\qquad r_+(u)=\frac{A(t_f)-A(u)}{t_f-u}.$$

直接代入可得

$$C(u)=\frac{(u-t_s)(t_f-u)}{D}\,[r_-(u)-r_+(u)].$$

因此最小值對應的是帶時間窗長度權重的率差，不是未加權率差最大的切點，也不保證前段平均率在所有切點中最低。若最小值為負，可推出該切點前段低於全窗平均、後段高於全窗平均；這項代數事實本身不是變點檢定。

在跳躍曲線上，最小值可能位於事件前的左極限。計算時須同時保留事件前、事件後的值，不能只用事件後累積量。更換分析窗會更換 $k$；搜尋很多窗後選最小者，也需要以同樣搜尋程序校準。

## D.2 測量誤差與回歸方向

設 $X=X^*+\epsilon$、$Y=a+bX^*+\delta$，且 $X^*,\epsilon,\delta$ 互不相關、誤差均值為零。令 $v=\mathrm{Var}(X^*)>0$，則

$$\mathrm{Cov}(X,Y)=bv,\quad
\mathrm{Var}(X)=v+\sigma_\epsilon^2,\quad
\mathrm{Var}(Y)=b^2v+\sigma_\delta^2.$$

所以正向最小平方斜率的大樣本極限是

$$b_{Y|X}=\frac{bv}{v+\sigma_\epsilon^2}
=\frac{b}{1+\sigma_\epsilon^2/v}.$$

反向斜率為 $b_{X|Y}=bv/(b^2v+\sigma_\delta^2)$；若 $b\ne0$，換回原座標得到

$$\frac1{b_{X|Y}}=b\left(1+\frac{\sigma_\delta^2}{b^2v}\right).$$

兩個方向一般不重合。若兩方向的誤差變異數相同，令
$s^2=\sigma_\epsilon^2/v=\sigma_\delta^2/v$，表示誤差變異數相對於真實訊號變異數的比例。
在 $b=-1$ 的例子中，兩支為 $-1/(1+s^2)$ 與 $-(1+s^2)$。這是有條件的誤差模型解釋，不能用兩條實測斜率就唯一辨認真斜率及所有誤差變異數。

由樣本平方和也有恆等式

$$\hat b_{Y|X}\hat b_{X|Y}
=\frac{S_{XY}^2}{S_{XX}S_{YY}}=\hat r^2.$$

帶主震固定截距的回歸，可先做組內去中心化，再使用相同關係。整體模型的 $R^2$ 包含截距解釋的組間差異，不等於組內相關係數平方。

reduced major axis 斜率為 $\mathrm{sign}(S_{XY})\sqrt{S_{YY}/S_{XX}}$，等於兩方向斜率在原座標下的帶符號幾何平均；它不是在所有誤差模型下都一致的真斜率估計。

## D.3 為什麼分群前後方向可以不同

令 $G$ 表示主震群組。全共變異數公式為

$$\mathrm{Cov}(X,Y)=\mathbb E[\mathrm{Cov}(X,Y\mid G)]
+\mathrm{Cov}(\mathbb E[X\mid G],\mathbb E[Y\mid G]).$$

第一項可為負，第二項可為正；混合後接近零，不表示兩層都沒有關係。要比較各主震，可以在組內模型 $Y_{ij}=a_i+bX_{ij}+e_{ij}$ 中選共同參考 $X_0$，以 $a_i+bX_0$ 比較，而不是讓每組不同的 $X$ 分布混入組間差異。

此外 $\mathrm{Var}(X+Y)=\mathrm{Var}(X)+\mathrm{Var}(Y)+2\mathrm{Cov}(X,Y)$。負共變異數可降低和的變異數，卻不能無條件保證和對第三變數的 $R^2$ 必然提升，更不能由此推導物理守恆。

## D.4 三個正規化核

令 $\tau=t-t_i>0$，$u=\log_{10}\tau$，$\mu_{T,i}=a_T+b_Tm_i$。若 $U\sim N(\mu_{T,i},\sigma_T^2)$，因 $du/d\tau=1/(\tau\ln10)$，

$$f(\tau\mid m_i)=\frac{1}{\tau\sigma_T\ln10\sqrt{2\pi}}
\exp\left[-\frac{(\log_{10}\tau-\mu_{T,i})^2}{2\sigma_T^2}\right].$$

令 $f=0$ 於 $\tau\le0$；換回 $u$ 後積分等於常態密度積分，故正規化為一。

規模核與空間核可寫成

$$g(m\mid m_i)=\frac{1}{\sigma_M\sqrt{2\pi}}
\exp\left[-\frac{(m-a_M-b_Mm_i)^2}{2\sigma_M^2}\right],$$

$$h(x,y\mid x_i,y_i,m_i)=\frac{1}{2\pi v_i}
\exp\left[-\frac{(x-x_i)^2+(y-y_i)^2}{2v_i}\right],\qquad
v_i=\sigma_A^2 10^{b_Am_i}.$$

$g$ 對整條實線正規化，$h$ 對平面正規化；實際預報區域或規模範圍只會收到其中一部分，不能把有限區域積分直接當成一。空間座標須為相同長度單位，不能將經緯度差直接視為公里。

三核相乘表達給定輸入事件後的可分離假設。面積尺度回歸只決定尺度如何隨規模變動，不唯一決定圓對稱常態形狀。若以圓形常態內含質量 $q$ 的面積表示，

$$P(R\le r)=1-e^{-r^2/(2v_i)},\qquad
A_q=-2\pi v_i\ln(1-q).$$

因此把面積截距換成 $\sigma_A$ 之前，需指定面積涵蓋的機率及幾何定義。

## D.5 對數常態的三個時間摘要

由 $\tau=e^{U\ln10}$ 及常態動差生成函數，

$$\mathrm{median}(\tau)=10^{\mu_T},\qquad
\mathbb E[\tau]=10^{\mu_T}e^{\sigma_T^2(\ln10)^2/2}.$$

微分 $\ln f(\tau)$：

$$\frac{d\ln f}{d\tau}=-\frac1\tau
-\frac{\log_{10}\tau-\mu_T}{\sigma_T^2\tau\ln10}=0,$$

得 $\mathrm{mode}(\tau)=10^{\mu_T-\sigma_T^2\ln10}$。$\sigma_T>0$ 時眾數小於中位數、小於平均數。當 $\sigma_T=0.20$，平均與中位數比約 $1.112$；當 $\sigma_T=0.81$，約 $5.69$。比較文獻的「等待時間」前須先確認使用哪個摘要。

## D.6 規模權重與有限輸入門檻

以 $z$ 表示輸入規模，假設輸入長期率密度為 $C e^{-\beta z}$，平均事件權重 $\bar w>0$ 不隨規模變化，且 $b_M>0$。採用

$$\eta(z)=\frac{b_M(1-\mu_E)}{\bar w}
\exp\left[-\beta\{a_M+(b_M-1)z\}-\frac12\beta^2\sigma_M^2\right].$$

令 $u=a_M+b_Mz$，並記 $\phi_\sigma$ 為中心零的常態密度。若先將積分下限形式上延伸至 $-\infty$，

$$\begin{aligned}
\int e^{-\beta z}\eta(z)\bar w\,g(m\mid z)\,dz
&=(1-\mu_E)e^{-\beta^2\sigma_M^2/2}
\int e^{-\beta u}\phi_{\sigma_M}(m-u)\,du\\
&=(1-\mu_E)e^{-\beta m}.
\end{aligned}$$

最後一步用換元 $v=m-u$，其積分為

$$e^{-\beta m}\int e^{\beta v}\phi_{\sigma_M}(v)\,dv
=e^{-\beta m+\beta^2\sigma_M^2/2}.$$

延伸下限是分析基準，不是實際存在無限多低規模事件的可用目錄。回到 $z\ge m_0$，以配方

$$-\beta u-\frac{(m-u)^2}{2\sigma_M^2}
=-\frac{[u-(m-\beta\sigma_M^2)]^2}{2\sigma_M^2}
-\beta m+\frac12\beta^2\sigma_M^2$$

可得保留比例

$$\Delta(m)=\Phi\left(\frac{m-a_M-b_Mm_0-\beta\sigma_M^2}{\sigma_M}\right).$$

有限門檻下的長期平均規模率為理想值乘 $\Delta(m)$。以 $1/\Delta(m)$ 補償或使用其他補償版本，都須說明假設；$\Delta$ 很小時會放大模型誤差。若平均權重隨規模改變，前面的常數 $\bar w$ 不能照搬。

## D.7 空間卷積與時間截斷

若未知主震位置 $X\sim N(x_i,\sigma_1^2)$，相對位移 $E\sim N(0,\sigma_2^2)$ 且獨立，餘震位置 $Z=X+E$，其密度為

$$\int\phi_{\sigma_2}(z-x)\phi_{\sigma_1}(x-x_i)\,dx
=\phi_{\sqrt{\sigma_1^2+\sigma_2^2}}(z-x_i).$$

完整配方可令 $v_*=(\sigma_1^{-2}+\sigma_2^{-2})^{-1}$，兩核指數之和為

$$-\frac{[x-v_*(z/\sigma_2^2+x_i/\sigma_1^2)]^2}{2v_*}
-\frac{(z-x_i)^2}{2(\sigma_1^2+\sigma_2^2)}.$$

對 $x$ 積分給 $\sqrt{2\pi v_*}$，與兩核常數合併即得上式。兩軸獨立時逐軸套用；一般多變數常態則將共變異數矩陣相加。

對時間核，給定某輸入事件，可用累積分布 $F_i$ 算其在預報窗 $[T_0,T_1]$ 的比例：

$$F_i(T_1-t_i)-F_i(T_0-t_i),$$

其中 $F_i(u)=0$ 於 $u\le0$。目錄開始前事件的貢獻仍未知，不能只靠這個有限窗積分恢復；需指定過去活動模型或使用文獻的完整度補償方法。

### D.7.1 目標地震收到多少前兆貢獻？

**時間完整度**（temporal completeness）衡量可用歷史保留了多少預期前兆貢獻。
它不是目錄中「有多少事件被偵測」的直接比例。
單顆事件的 $F_i$ 只能回答該事件的有限窗積分。
要問規模 $m$ 的目標還缺多少，必須把可能的輸入規模一起平均。

Rhoades 等（2020）§2.3–2.4、式 10–17 以平穩歷史及 GR 規模分布為前提。
以下用與本書核記號一致的積分寫出同一結構。
先以共同平均空間因子處理區域活動；空間邊界造成的損失需另外檢查。
比例的分母須為正，且來源權重與積分都須有限。
令輸入規模為 $z$、長期輸入率密度為 $r_{\rm in}(z)$，時間核累積分布為 $F_z$。
以 $g(m\mid z)$ 表示規模核，並令

$$W(z;m)=\eta(z)g(m\mid z)r_{\rm in}(z),\qquad
r_{\rm in}(z)\propto10^{-bz}.$$

共同的規模門檻補償若只依目標 $m$ 變化，會在下面的比例中消去。
等權、$b_M=1$ 時，$\eta(z)$ 也是常數；若放寬這些條件，就須在分子、分母保留相同權重。

**目錄前置時間**（catalogue lead time）從目錄起點量到目標地震的時刻。
令目錄起點為 $t_{\rm cat}$、目標時刻為 $t$，則 $L_{\rm cat}=t-t_{\rm cat}$。
它不是發報到目標的距離，也不是發報到資料截止的間隔。
Rhoades 等（2020）緒論的定義及 §2.3 式 10–13 都以目錄來源時間來計算前兆貢獻。
下面改用明確的兩個端點表示其時間積分，避免不同的 $L$ 記法混用。

**輸入延遲**（input delay）記為 $d$，由發報時刻 $t_f$ 回推資料截止 $t_f-d$。
**歷史時間差**（time lag）記為 $T_{\rm lag}=t-(t_f-d)$。
它從資料截止量到目標，因此窗內的目標越晚，time lag 越長。
**預報期限**（forecast horizon）是發報到預報窗結束的長度 $\Delta t$。
每個目標距發報的 $t-t_f$ 位於該期限內；這又是另一個時間量。

| 量 | 兩個端點 | 本站的意義 |
|---|---|---|
| $L_{\rm cat}$ | 目錄起點 → 目標時刻 | 1960 年起的歷史到該目標有多長 |
| $\Delta t$ | 發報時刻 → 預報窗末 | 通常為 91.31 天 |
| $d$ | 資料截止 → 發報時刻 | PPE／EEPAS 設 50 天 |
| $T_{\rm lag}$ | 資料截止 → 目標時刻 | $(t-t_f)+d$，隨窗內目標改變 |

當 $L_{\rm cat}\ge T_{\rm lag}\ge0$，來源時間區間是
$[t-L_{\rm cat},t-T_{\rm lag}]=[t_{\rm cat},t_f-d]$。
對每種輸入規模平均後，預期保留的貢獻為

$$c(T_{\rm lag},L_{\rm cat},m)=\int_{m_0}^{m_u^{\rm in}}W(z;m)
\left[F_z(L_{\rm cat})-F_z(T_{\rm lag})\right]\,dz.$$

若 $L_{\rm cat}<T_{\rm lag}$，沒有可用歷史，本式改為零而非負積分。
$m_u^{\rm in}$ 是此補償模型指定的輸入上限，不可悄悄當成本站輸出箱上限。
此處 $c$ 的第二引數是最早來源距目標的時間差；若文獻寫成上限 $T+L$，須先換算這個端點。
只考慮目錄起點截尾、令 $T_{\rm lag}=0$ 時，原文記為 $p(L,m)$。
這裡寫 $p_{\rm comp}$，避免混同 Omori 指數：

$$p_{\rm comp}(L_{\rm cat},m)=\frac{c(0,L_{\rm cat},m)}{c(0,\infty,m)}\in[0,1].$$

因 $F_z(L_{\rm cat})$ 隨 $L_{\rm cat}$ 增加，這個比例也不減；$L_{\rm cat}\to\infty$ 時趨近 1。
在原文 EEPAS 尺度關係下，大目標通常對應更久以前的輸入，因此同一 $L_{\rm cat}$ 的完整度較低。
這個規模方向取決於核與尺度設定，不能推成任意模型的定理。
$M_c$、$\Delta(m)$ 與 $p_{\rm comp}(L_{\rm cat},m)$ 分別處理目錄偵測、輸入規模截尾、時間截尾。

本站 50 天 delay 不會把所有目標的 time lag 固定成 50 天，也不會讓時間核往後平移。
固定 1960 年目錄起點時，catalogue lead time 隨目標日期增長。
**固定前置時間 EEPAS**（Fixed Lead Time EEPAS，FLEEPAS）則讓歷史起點隨目標時間移動，
使每個目標都採同長的 catalogue lead time。
它不能與本書固定目錄起點、逐窗發報的設定混稱。
本節補償例先處理 $T_{\rm lag}=0$ 的目錄起點損失；若同時補最近資料截尾，須使用兩端點比例。

### D.7.2 補到平滑背景，或保留時變形狀？

將有限歷史的時變部分記為 $v(t,m,x,y)$，原率為 $\mu_E\lambda_0+v$。
在上述平均模型下，若完整目標規模率為 $r(m)$，則

$$\overline v=(1-\mu_E)p_{\rm comp}(L_{\rm cat},m)r(m).$$

Rhoades 等（2020）式 15–16 提出兩個補償端點。
省略共同引數，令 $p_c=p_{\rm comp}(L_{\rm cat},m)>0$：

$$\lambda_{\rm smooth}=\left[\mu_E+(1-\mu_E)(1-p_c)\right]\lambda_0+v,$$

$$\lambda_{\rm signal}=\mu_E\lambda_0+\frac{v}{p_c}.$$

第一式把缺少的平均量補到背景圖，第二式放大現有時變形狀。
只要 $\overline{\lambda_0}=r(m)$，兩式的長期平均都回到 $r(m)$。
但有限窗中的空間與時間圖不相同。
在兩者間作凸組合，可由學習資料選擇分配方式；這便是補償版本的一個核心選擇。
$p_c$ 很小時，放大操作也會放大模型誤差；$p_c=0$ 時第二式無定義。
平均量補回來不代表找回遺失事件的真實位置。
本站四張預報快取沒有加入這套時間補償，不能把本節公式當成已實作的模型。

### D.7.3 EAS：預報主震之後，還可能跟著餘震

**EEPAS 的餘震延伸**（EAS）處理尚未發生的預報主震可能帶來的餘震。
Rhoades（2009）式 10–16 將主震率 $\lambda_M$ 與餘震率 $\lambda_A$ 分開，再加背景：

$$\lambda_{\rm EAS}=\mu_E\lambda_0+\lambda_M+\lambda_A.$$

令 $u$、$\mathbf r$、$m'$ 是未來主震的時間、位置與規模。
令 $k_A(t,m,\mathbf x\mid u,m',\mathbf r)$ 為該主震帶來的餘震期望率核。
它包含指定餘震產能，不必是積分為一的機率密度。
在發報後 $t_f<u<t$ 的可能主震上平均，結構可寫成

$$\lambda_A(t,m,\mathbf x)=\int_{t_f}^{t}\int\!\!\int
\lambda_M(u,m',\mathbf r)\,
 k_A(t,m,\mathbf x\mid u,m',\mathbf r)\,d\mathbf r\,dm'\,du.$$

這是對「可能發生的主震」積分，不是等主震真的出現才更新。
原文利用餘震時間尺度較短的近似，簡化時間卷積；空間部分則出現 D.7 的常態卷積。
該近似需尺度分離，不是任意短窗都成立，也不能據此宣稱已涵蓋無限後代。

| 做法 | 調整哪一部分 | 不能混同的事 |
|---|---|---|
| EEPAS-NW | 已知輸入事件等權 | 不代表每顆都是已確認的物理前兆 |
| EEPAS-W | 以事件權重降低部分已知餘震的輸入貢獻 | 沒有因此補上未來主震的餘震 |
| EAS | 對預報主震整合其可能餘震率 | 不是只把既有來源的權重改掉 |
| 本站 ETAS 第一代近似 | 計算已知歷史的直接觸發 | 未含窗內新事件繼續繁殖 |

原文 EEPAS-W 的事件權重來自其採用的除叢／餘震模型，不能只憑時間距離任意指定。
本站採 NW，沒有估那些權重，也沒有新增 EAS 快取或成績。
本節讓讀者辨認模型結構，數值來源仍限已發表研究。

## D.8 義大利的三階段擬合

Biondini 等（2023）第 3 節將 $b_M$ 全程固定為一，其餘參數分成三次最佳化。以下 $\mu_E$ 對應原文的 $\mu$。

| 階段 | 固定參數 | 估計參數 |
|---|---|---|
| 第一次 | $b_T=0.40,b_A=0.35,\sigma_M=0.32,\sigma_T=0.23$ | $a_T,a_M,\sigma_A,\mu_E$ |
| 第二次 | 第一次得到的 $a_T,a_M,\sigma_A$ | $b_T,b_A,\sigma_M,\sigma_T,\mu_E$ |
| 第三次 | $b_M=1$ | 上述八個參數一起估計 |

第三次以先前結果為起始值；$\mu_E$ 是三次都重新估計的參數。這套安排提供可解釋的起點，但不是全域最佳保證。第 15 章直接使用論文表 3 的 EEPAS-NW 已發表結果，不在教學頁執行最佳化。

## D.9 發報截止與格箱積分

本站時間以天、位置以公里計。設發報時間為 $t_0$、delay 為 $d=50$ 天，該窗來源集合固定為 $I_0=\{i:t_i\le t_0-d,m_i\ge m_0\}$。在窗 $[t_0,t_1)$、空間格 $R_j$、規模箱 $B_k$ 中，時變部分的期望數為

$$\Lambda^{\mathrm{var}}_{jk}=
\sum_{i\in I_0}\eta(m_i)
\left[\int_{t_0}^{t_1}f_i(t)\,dt\right]
\left[\int_{B_k}\frac{g_i(m)}{\Delta(m)}\,dm\right]
\left[\int_{R_j}h_i(x,y)\,dx\,dy\right].$$

最後再加 $\mu_E\Lambda^{\mathrm{PPE}}_{jk}$。可分離假設使三個低維積分相乘；$g_i/\Delta$ 是補償後的規模貢獻，不能再當成總積分為一的原始密度。

對方格 $[x_0,x_1]\times[y_0,y_1]$，空間積分等於兩個常態累積分布差的乘積。本站時間積分也使用累積分布差，規模箱使用五個中點近似。模型總數的驗證必須涵蓋這些有限區域積分，不能只檢查無限範圍核正規化。

完整歷史的可用起點、發報截止與 50 天 delay 是三個不同設定。規模門檻補償也不會恢復目錄起點以前的事件；時間完整度延伸須另外指定模型。暖機資料早期漏測的影響，見第 8 章及 Rhoades 等（2020）。

## 參考資料與延伸閱讀

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 2 節的 Ψ 現象與尺度關係，再讀限制與未解問題，掌握本章經驗觀察如何連到機率模型。

- [Algorithmic Identification of the Precursory Scale Increase Phenomenon in Earthquake Catalogs](https://doi.org/10.1785/0220240233) — Annemarie Christophersen、David A. Rhoades、Sebastian Hainzl，2024，*Seismological Research Letters*；[免費機構典藏全文](https://gfzpublic.gfz.de/rest/items/item_5029405_4/component/file_5029659/content)。這是本章自動辨識與對照實驗的已發表來源，建議比較矩形、圓形搜尋與隨機化目錄的設計；辨識到統計現象仍須與前瞻預報能力分開判斷。

- [Space–Time Trade-Off of Precursory Seismicity in New Zealand and California Revealed by a Medium-Term Earthquake Forecasting Model](https://doi.org/10.3390/app112110215) — Sepideh J. Rastin、David A. Rhoades、Annemarie Christophersen，2021，*Applied Sciences*；免費開放全文。研究時間與空間參數的取捨，適合延伸本章「同一事件的 Ψ 辨識不唯一」及參數解讀問題。

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 3–5 節，建立 EEPAS、模型組合與缺失前兆補償的全貌，再回到本章的正規化推導。

- [Long-range Earthquake Forecasting with Every Earthquake a Precursor According to Scale](https://doi.org/10.1007/s00024-003-2434-9) — David A. Rhoades、Frank F. Evison，2004，*Pure and Applied Geophysics*；全文可能需訂閱。這是 EEPAS 原始論文，重點是如何把尺度關係轉成每個事件對未來地震率的貢獻，而不是先判定哪個事件必然是前兆。

- [Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123) — Emanuele Biondini、David A. Rhoades、Paolo Gasperini，2023，*Geophysical Journal International*；出版社全文可能需訂閱，[免費機構典藏全文](https://www.earth-prints.org/handle/2122/17084)。將本章公式連到義大利目錄的實際應用，閱讀 PPE、ETAS 與 EEPAS 的比較時，特別留意學習期、測試期和預報時間窗。

- [The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data](https://doi.org/10.3390/e22111264) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2020，*Entropy*；免費開放全文。說明目錄開始前遺漏的事件如何影響 EEPAS，適合延伸本章的 lead time 與時間完整度補償推導。

- [Long-range earthquake forecasting allowing for aftershocks](https://doi.org/10.1111/j.1365-246X.2008.04083.x) — D. A. Rhoades，2009，*Geophysical Journal International*；[出版社網頁全文](https://academic.oup.com/gji/article/178/1/244/644120)可免費閱讀。閱讀 EEPAS 如何加入預報事件的餘震貢獻，對照本章 EAS 延伸與「降低輸入餘震權重」的不同角色。

- [Application of a long-range forecasting model to earthquakes in the Japan mainland testing region](https://doi.org/10.5047/eps.2010.08.002) — David A. Rhoades（2011），Earth, Planets and Space（免費全文）。檢視不同目標規模的擬合差異與已發表比較的回溯條件。
