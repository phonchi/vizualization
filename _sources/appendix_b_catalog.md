# 附錄 B：目錄統計與叢集律的推導

本附錄接續 {doc}`11_catalog_completeness_b` 與 {doc}`12_clustering_laws`。先備知識為機率密度、期望值、微分與最大概似；不需要先讀完本附錄才能繼續主文。以下推導會在分布的支撐集中明確交代「資料的選取方式」，因為門檻、格距與時間窗都是估計問題的一部分。

## B.1 從 GR 累積數到規模密度

固定分析區域與期間，設門檻以上的期望累積數為 $N(\ge m)=10^{a-bm}$，$m\ge m_0$。除以門檻以上的總量，得到

$$S(m)=\frac{N(\ge m)}{N(\ge m_0)}=10^{-b(m-m_0)}=e^{-\beta(m-m_0)},\qquad\beta=b\ln10.$$

對 $m$ 取負導數，即為主文的 {eq}`eq:gr-density`。這是未截斷 GR 的機率表示；若設上限 $M_{\max}$，須另除以 $1-e^{-\beta(M_{\max}-m_0)}$，其影響見附錄 C。

若兩種規模尺度有確定關係 $M_2=a_sM_1+d_s$、$a_s>0$，則

$$a_1-b_1M_1=\left(a_1+\frac{b_1d_s}{a_s}\right)-\frac{b_1}{a_s}M_2.$$

所以 $b_2=b_1/a_s$。這只處理代數轉換，不包含迴歸係數、規模量測與轉換殘差的不確定性。

## B.2 連續規模的概似與有限樣本偏差

設 $X_i=M_i-m_0$ 是 $N$ 個獨立、完整的指數樣本，則

$$\ell(\beta)=N\ln\beta-\beta\sum_iX_i,\qquad
\frac{\partial\ell}{\partial\beta}=\frac N\beta-\sum_iX_i.$$

因此 $\widehat\beta=N/\sum_iX_i=1/\bar X$，$\widehat b=\widehat\beta/\ln10$。資訊量為 $N/\beta^2$，給出大樣本近似 $\operatorname{SE}(\widehat b)\simeq b/\sqrt N$。

有限樣本時，$S_N=\sum_iX_i\sim\operatorname{Gamma}(N,\text{rate}=\beta)$，而不是未調整尺度的 $\bar X$ 服從這個 Gamma 分布。對 $N>1$，

$$\mathbb E[S_N^{-1}]=\frac\beta{N-1},\qquad
\mathbb E[\widehat\beta]=\frac{N}{N-1}\beta.$$

因此 $(N-1)\widehat\beta/N$ 在這個固定門檻的理想模型下無偏。若門檻由同一資料選出、規模離散或存在漏測，不能不加條件地沿用這個修正。

對一般樣本標準差 $s_M$，delta method 給出 Shi–Bolt 型近似

$$\operatorname{SE}(\widehat b)\simeq\ln10\,\widehat b^2\frac{s_M}{\sqrt N}.$$

理想指數模型下 $s_M\to1/(b\ln10)$，便回到 $b/\sqrt N$。這些公式只反映所假設抽樣結構下的隨機誤差。

## B.3 離散 GR 與精確最大概似

規模格點為 $m_0+k\Delta M$，$k=0,1,\ldots$。在固定格距、門檻以上完整的離散 GR 模型中，

$$P(K=k)=(1-r)r^k,\qquad r=e^{-\beta\Delta M}.$$

給定 $N$ 筆格點資料，

$$\ell(r)=N\ln(1-r)+\left(\sum_i k_i\right)\ln r.$$

令導數為零，得 $\widehat r=\bar k/(1+\bar k)$。寫 $u=\bar m-m_0=\Delta M\bar k$，便得到主文 {eq}`eq:b-exact`：

$$\widehat b=\frac{1}{\Delta M\ln10}\ln\left(1+\frac{\Delta M}{u}\right).$$

當 $u=0$，最優值在 $r=0$，對應 $b\to\infty$；資料無法給出有限斜率。若格距不固定或包含不同四捨五入規則，需使用相應的區間機率重新建立概似。

用半格 $\delta=\Delta M/2$ 表示，與文獻中的反雙曲餘切寫法相同，因為

$$\operatorname{arccoth}x=\tfrac12\ln\frac{x+1}{x-1},\qquad
\frac1{\delta\ln10}\operatorname{arccoth}\frac{u+\delta}{\delta}
=\frac1{2\delta\ln10}\ln\left(1+\frac{2\delta}{u}\right).$$

這裡的「精確」指所假設離散模型的最大概似解，不代表估計量有限樣本無偏，也不代表目錄必然符合模型。

## B.4 半格近似與區間為何不對稱

令 $\delta=\Delta M/2$。在 $2\delta/u<1$ 的展開範圍內，

$$\begin{aligned}
 b_{\rm exact}&=\frac1{\ln10}\left(\frac1u-\frac\delta{u^2}+\frac{4\delta^2}{3u^3}+\cdots\right),\\
 b_{\rm Utsu}&=\frac1{\ln10}\frac1{u+\delta}
 =\frac1{\ln10}\left(\frac1u-\frac\delta{u^2}+\frac{\delta^2}{u^3}+\cdots\right).
\end{aligned}$$

差異從 $\delta^2/(3u^3\ln10)$ 開始。這是兩種公式在相同 $u$ 上的差異，不應說成近似式在每種抽樣情境都必然低估真實 $b$；估計的抽樣偏差與資料偏差還要分開。

令 $f(u)=\ln(1+2\delta/u)/(2\delta\ln10)$，則

$$f'(u)=-\frac1{\ln10\,u(u+2\delta)},\qquad
f''(u)=\frac{2u+2\delta}{\ln10\,u^2(u+2\delta)^2}>0.$$

因此將 $u$ 上的對稱區間轉換為 $b$ 時，方向反轉且上側常較長。Tinti 與 Gasperini（2024）討論的離散近似區間，可用 $C=1+2\delta/u$ 寫成

$$\begin{aligned}
b_1&=\frac1{2\delta\ln10}\ln\frac{C+\sqrt{C/N}}{1+\sqrt{C/N}},\\
b_2&=\frac1{2\delta\ln10}\ln\frac{C-\sqrt{C/N}}{1-\sqrt{C/N}}.
\end{aligned}$$

這個近似需要平均值的分布足以使用相應漸近展開，且 $N>C$ 使端點有意義；它不是任意小樣本的精確涵蓋區間。正文的 bootstrap 圖比較區間形狀，並未替其在所有目錄上的涵蓋率做驗證。

## B.5 規模差、修剪與相關性

若 $X,Y$ 獨立且同為 $\operatorname{Exp}(\beta)$，對 $d\ge0$，

$$f_{X-Y}(d)=\int_0^\infty\beta e^{-\beta(y+d)}\beta e^{-\beta y}\,dy
=\frac\beta2e^{-\beta d}.$$

由對稱性，$f_{X-Y}(d)=\beta e^{-\beta|d|}/2$。條件在 $X-Y>d_0$，超額差值 $X-Y-d_0$ 仍服從 $\operatorname{Exp}(\beta)$。這是差分式估計的理想起點；實際的漏測與配對選擇是否保留這個分布，需要另行分析。

離散情形中，若 $K_1,K_2$ 獨立服從 $(1-r)r^k$，

$$P(K_1-K_2=j)=\frac{1-r}{1+r}r^{|j|}.$$

零差值有自己的機率質量。若只保留非零絕對差，$P(|j|=k\mid|j|\ge1)=(1-r)r^{k-1}$，$k\ge1$，回到移位幾何分布。

即使規模本身獨立，相鄰差 $D_i=X_{i+1}-X_i$ 也共用觀測：

$$\operatorname{Cov}(D_i,D_{i+1})=-\frac1{\beta^2},\qquad
\operatorname{Var}(D_i)=\frac2{\beta^2}.$$

符號差的相關係數為 $-1/2$。**取絕對值後不能沿用負相關結論**。令 $A_i=|X_{i+1}-X_i|$，利用條件獨立與 $\mathbb E[|X-y|]=y-1/\beta+2e^{-\beta y}/\beta$，可得

$$\operatorname{Cov}(A_i,A_{i+1})=\frac1{3\beta^2},\qquad
\operatorname{Var}(A_i)=\frac1{\beta^2}.$$

所以相鄰絕對差在理想模型下是正相關。選取正差、修剪或採用其他配對後，相關性還會改變；不能把某一篇模擬中的有效樣本數比例當成通用修正係數。

## B.6 Omori 核：有限窗與無限時間

設 $c>0$、$0\le S<T<\infty$。對有限的 $p$，

$$I(S,T;c,p)=\int_S^T(t+c)^{-p}\,dt
=\begin{cases}
\dfrac{(T+c)^{1-p}-(S+c)^{1-p}}{1-p},&p\ne1,\\
\ln\dfrac{T+c}{S+c},&p=1.
\end{cases}$$

兩種寫法在 $p\to1$ 時連續。只有把上限推到無限遠時，才需要 $p>1$，並得到 $I(0,\infty;c,p)=c^{1-p}/(p-1)$。除以這個量，即得 {eq}`eq:omori-density`。

若指定觸發期限 $T_c$，對所有有限 $p$ 都可定義

$$g_{T_c}(t)=\frac{(t+c)^{-p}}{I(0,T_c;c,p)}\mathbf1_{[0,T_c]}(t).$$

這是改變觸發時間模型，不只是裁切觀測窗。若保持未正規化係數 $K$ 不變，期限內產能為 $K I(0,T_c;c,p)$，會依期限而變；若改為保持正規化產能 $A$ 不變，則改變的是後代的時間分配。這兩種操作需分清楚。

## B.7 Omori 的剖面概似

對固定 $[S,T]$ 中的 $N$ 個事件，非齊次 Poisson 假設下

$$\ell(K,c,p)=N\ln K-p\sum_i\ln(t_i+c)-K I(S,T;c,p).$$

若 $N>0$，令 $\partial\ell/\partial K=0$，可得 $\widehat K(c,p)=N/I$，代回後為

$$\ell_{\rm prof}(c,p)=N\ln\frac N{I(S,T;c,p)}-p\sum_i\ln(t_i+c)-N.$$

原本三個參數的搜尋降為兩個。若用負剖面概似的 Hessian 估計標準誤，需確認最優點在內部且資訊矩陣正定。對負變異數取絕對值再開根號，沒有統計依據；這種情況應標示區間無法由局部二次近似得到。$N=0$ 時不能使用上式的同一數值流程。

## B.8 空間核的尺度與正規化

寫 $\sigma(m)=D e^{\gamma(m-m_0)}$，採用正文的各向同性核。極座標下面積元為 $r\,dr\,d\varphi$，令 $v=1+r^2/\sigma$，可得

$$\int_{\mathbb R^2}f\,dx\,dy=(q-1)\int_1^\infty v^{-q}\,dv=1,\qquad q>1.$$

若面積尺度依 $10^{m-m_0}$ 變化，$\gamma=\ln10$；相應長度尺度 $\sqrt\sigma$ 的規模指數為 $\gamma/2$。經驗擬合的差異不能全歸因於記號：還可能包含區域、幾何、定位誤差與模型形狀的差異。

## B.9 最大值與最大規模差

設 $X_1,\ldots,X_N$ 為固定數量的獨立 $\operatorname{Exp}(\beta)$ 樣本，最大值的累積分布為

$$P(X_{\max}\le x)=(1-e^{-\beta x})^N,\qquad x\ge0.$$

由大到小排序 $X_{(1)}>\cdots>X_{(N)}$，定義 $D_k=X_{(k)}-X_{(k+1)}$ 及 $D_N=X_{(N)}$。利用指數分布的無記憶性，最小值先以率 $N\beta$ 出現，剩餘 $N-1$ 個扣掉它後仍是獨立指數，遞迴得到 $D_k\sim\operatorname{Exp}(k\beta)$ 且各間距獨立。因此

$$\begin{aligned}
\mathbb E[X_{(1)}]&=\frac1\beta\sum_{k=1}^N\frac1k,\\
\operatorname{Var}(X_{(1)})&=\frac1{\beta^2}\sum_{k=1}^N\frac1{k^2},\\
\mathbb E[X_{(1)}-X_{(2)}]&=\frac1\beta\quad(N\ge2).
\end{aligned}$$

當 $N\to\infty$，最大值的標準差趨於 $\pi/(\sqrt6\beta)$，不會隨樣本數消失。最大與次大之差仍有整個指數分布，故 $1/\beta$ 是期望，並非硬下界。

若把某個主震規模固定、再獨立抽 $N$ 個餘震規模，這個公式可以計算最大餘震的參考分布；若再要求餘震不得超過主震，便改變了抽樣條件，必須重新截斷分布。依最大值選群、依規模決定家族大小的情況，也不能直接使用固定 $N$ 的獨立樣本結論。

## B.10 除叢後的 GR 交點

設兩條擬合直線為 $a-bm$ 與 $a_{\rm main}-b_{\rm main}m$。相減令零即得 {eq}`eq:mx`，而累積數比值為

$$r(m)=10^{a_{\rm main}-a+(b-b_{\rm main})m}.$$

若 $b_{\rm main}<b$，比值隨規模增加，通過交點後超過 1。子目錄的觀測計數不可能超過全目錄，因此不相容的是兩條擬合線的共同外推，而非真實資料違反集合關係。由此不能直接推斷完整 PSHA 的高估或低估，因為危害度還包含位置、地動模型與其他權重。

## 參考資料與延伸閱讀

- Tinti, S. 與 Gasperini, P.（2024），[The estimation of b-value of the frequency–magnitude distribution and of its 1σ intervals from binned magnitude data](https://doi.org/10.1093/gji/ggae159)。離散最大概似、區間與差分方法的主要對照來源；[免費機構版本](https://cris.unibo.it/handle/11585/980514)。
- SeismoStats 團隊，[Estimate b-value](https://seismostats.readthedocs.io/latest/user/estimate_b.html) 與 [Estimate magnitude of completeness](https://seismostats.readthedocs.io/latest/user/estimate_mc.html)。免費文件，可將本附錄的門檻、格距與估計前提對照到實際資料格式。
- Utsu, T., Ogata, Y. 與 Matsu'ura, R. S.（1995），[The Centenary of the Omori Formula for a Decay Law of Aftershock Activity](https://doi.org/10.4294/jpe1952.43.1)。[免費全文](https://www.jstage.jst.go.jp/article/jpe1952/43/1/43_1_1/_article)，對照有限觀測窗、$c,p$ 與估計問題。
- Ogata, Y. 與 Zhuang, J.（2006），[Space–time ETAS models and an improved extension](https://doi.org/10.1016/j.tecto.2005.10.016)。[免費機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)，對照空間核和產能的分離；出版社可能需訂閱。
- Chan, C.-H. 與 Wu, Y.-M.（2013），[Maximum magnitudes in aftershock sequences in Taiwan](https://doi.org/10.1016/j.jseaes.2013.05.006)。比較最大餘震的不同統計摘要及序列選取；出版社全文可能需訂閱。
- Mizrahi, L., Nandan, S. 與 Wiemer, S.（2021），[The Effect of Declustering on the Size Distribution of Mainshocks](https://doi.org/10.1785/0220200231)。[免費作者預印本](https://arxiv.org/abs/2012.09053)，用於理解選群如何改變規模分布及外推。

返回 {doc}`11_catalog_completeness_b` 或 {doc}`12_clustering_laws`。
