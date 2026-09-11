
# 附錄 D：Ψ 的量測與 EEPAS 的核函數

本附錄需要條件機率、簡單線性迴歸與一變數積分。它補足 {doc}`Ψ <15_psi_phenomenon>`及 {doc}`EEPAS <16_eepas_ppe>`的計算；讀完每節可回到對應主文繼續閱讀。

## D.1 累積規模異常的性質與限制

令 $q_i=M_i-(m_c-0.1)>0$，$A(t)=\sum_{t_s<t_i\le t}q_i$，$D=t_f-t_s$。由主文 {eq}`eq:cumag`，$k=A(t_f)/D$，故

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

$g$ 對整條實線正規化，$h$ 對平面正規化；實際預報區域或規模範圍只會收到其中一部分，不能把有限區域積分直接當成一。空間座標須為相同長度單位，不能將經緯度差直接視為公裡。

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

有限門檻下的時變長期規模率為理想值乘 $\Delta(m)$。以 $1/\Delta(m)$ 補償或使用其他補償版本，都須說明假設；$\Delta$ 很小時會放大模型誤差。若平均權重隨規模改變，前面的常數 $\bar w$ 不能照搬。

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

## 參考資料與延伸閱讀

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 2 節的 Ψ 現象與尺度關係，再讀限制與未解問題，掌握本章經驗觀察如何連到機率模型。

- [Algorithmic Identification of the Precursory Scale Increase Phenomenon in Earthquake Catalogs](https://doi.org/10.1785/0220240233) — Annemarie Christophersen、David A. Rhoades、Sebastian Hainzl，2024，*Seismological Research Letters*；[免費機構典藏全文](https://gfzpublic.gfz.de/rest/items/item_5029405_4/component/file_5029659/content)。這是本章自動辨識與對照實驗的已發表來源，建議比較矩形、圓形搜尋與隨機化目錄的設計；辨識到統計現象仍須與前瞻預報能力分開判斷。

- [Space–Time Trade-Off of Precursory Seismicity in New Zealand and California Revealed by a Medium-Term Earthquake Forecasting Model](https://doi.org/10.3390/app112110215) — Sepideh J. Rastin、David A. Rhoades、Annemarie Christophersen，2021，*Applied Sciences*；免費開放全文。研究時間與空間參數的取捨，適合延伸本章「同一事件的 Ψ 辨識不唯一」及參數解讀問題。

- [A 20-Year Journey of Forecasting with the “Every Earthquake a Precursor According to Scale” Model](https://doi.org/10.3390/geosciences12090349) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2022，*Geosciences*；免費開放全文。先讀第 3–5 節，建立 EEPAS、模型組合與缺失前兆補償的全貌，再回到本章的正規化推導。

- [Long-range Earthquake Forecasting with Every Earthquake a Precursor According to Scale](https://doi.org/10.1007/s00024-003-2434-9) — David A. Rhoades、Frank F. Evison，2004，*Pure and Applied Geophysics*；全文可能需訂閱。這是 EEPAS 原始論文，重點是如何把尺度關係轉成每個事件對未來地震率的貢獻，而不是先判定哪個事件必然是前兆。

- [Application of the EEPAS earthquake forecasting model to Italy](https://doi.org/10.1093/gji/ggad123) — Emanuele Biondini、David A. Rhoades、Paolo Gasperini，2023，*Geophysical Journal International*；出版社全文可能需訂閱，[免費機構典藏全文](https://www.earth-prints.org/handle/2122/17084)。將本章公式連到義大利目錄的實際應用，閱讀 PPE、ETAS 與 EEPAS 的比較時，特別留意學習期、測試期和預報時間窗。

- [The Effect of Catalogue Lead Time on Medium-Term Earthquake Forecasting with Application to New Zealand Data](https://doi.org/10.3390/e22111264) — David A. Rhoades、Sepideh J. Rastin、Annemarie Christophersen，2020，*Entropy*；免費開放全文。說明目錄開始前遺漏的事件如何影響 EEPAS，適合延伸本章的 lead time 與時間完整度補償推導。

- [Long-range earthquake forecasting allowing for aftershocks](https://doi.org/10.1111/j.1365-246X.2008.04083.x) — D. A. Rhoades，2009，*Geophysical Journal International*；[出版社網頁全文](https://academic.oup.com/gji/article/178/1/244/644120)可免費閱讀。閱讀 EEPAS 如何加入預報事件的餘震貢獻，對照本章 EAS 延伸與「降低輸入餘震權重」的不同角色。

- [Application of a long-range forecasting model to earthquakes in the Japan mainland testing region](https://doi.org/10.5047/eps.2010.08.002) — David A. Rhoades（2011），Earth, Planets and Space（免費全文）。檢視不同目標規模的擬合差異與主文增益圖的回溯條件。
