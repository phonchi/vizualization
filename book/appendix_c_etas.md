# 附錄 C：ETAS 的分支、概似與推論

本附錄接續 {doc}`13_etas_structure` 與 {doc}`14_etas_estimation`。需要機率密度、期望、概似與基本矩陣微分；涉及機率生成函數的段落會先定義符號。正文先解釋各個結果的直覺，這裡集中列出推導與適用條件。

## C.1 正規化與分支比

使用 {eq}`eq:etas-intensity`，時間、空間及新事件規模密度都正規化時，規模 $m$ 的事件有平均 $\kappa(m)=Ae^{\alpha(m-m_0)}$ 個直接後代。無上界 GR 下，

$$\begin{aligned}
n&=\int_{m_0}^{\infty}\kappa(m)s(m)\,dm\\
 &=A\beta\int_0^{\infty}e^{-(\beta-\alpha)u}\,du
 =\frac{A\beta}{\beta-\alpha},\qquad \alpha<\beta.
\end{aligned}$$

這給出正文 {eq}`eq:branching-ratio` 的有限均值條件。$\alpha\ge\beta$ 時平均直接後代數無限，不能使用有限均值的平穩分支公式；它不同於「有限歷史上的每一項條件率都無限」。

若 $\Delta=M_{\max}-m_0>0$，截斷規模密度為

$$s_T(m)=\frac{\beta e^{-\beta(m-m_0)}}{1-e^{-\beta\Delta}},\qquad m_0\le m\le M_{\max}.$$

因此

$$n_T=\frac{A\beta}{1-e^{-\beta\Delta}}
\begin{cases}
\dfrac{1-e^{-(\beta-\alpha)\Delta}}{\beta-\alpha},&\alpha\ne\beta,\\
\Delta,&\alpha=\beta.
\end{cases}$$

兩種寫法在 $\alpha\to\beta$ 時相接，且對有限 $\alpha$ 都有有限值。$\alpha<\beta$ 且 $\Delta\to\infty$ 時回到未截斷公式。有限 $n_T$ 仍可能大於 1；規模截斷與次臨界是不同條件。

### 改變輸入門檻並不是單純刪除資料

若假設同一個直接觸發機制可向下延伸，將門檻降低 $d>0$，規模尾端依 GR 增加 $e^{\beta d}$ 倍，而親代產能相對新門檻需再調整 $e^{-\alpha d}$。在未截斷、相同形狀的外推假設下，$A'=Ae^{(\beta-\alpha)d}$，從而 $n'=ne^{(\beta-\alpha)d}$。

這是指定模型外推，不是對任意重估 ETAS 的恆等式。若只是把小事件從觀測目錄刪掉，未觀測中間世代可能改變有效觸發核；固定上限時也須重新計算規模正規化。因此不能說 $\alpha=\beta$ 時分支比在所有門檻設定下都不變。

## C.2 世代期望、平穩率與局部有限性

若每個成員的規模獨立抽自相同 $s(m)$，且 $n<\infty$，第 $k$ 代平均事件數為 $n^k$。完整家族包括始祖，故 $n<1$ 時

$$\mathbb E[Y]=\sum_{k=0}^{\infty}n^k=\frac1{1-n}.$$

在有限均值的平穩時間 ETAS 中，對條件率取期望可得

$$\bar\lambda=\mu+n\bar\lambda,\qquad\bar\lambda=\frac\mu{1-n}.$$

平均觸發率為 $n\bar\lambda$，其佔比為 $n$。有限窗實現、沒有窗前歷史的起始模擬，以及有空間外流的區域都不必滿足同樣的實現比例。

$n>1$ 時不存在這個正背景、有限均值的平穩解，但不能據此說每份有限目錄的概似都發散。若背景在有限窗可積、已知歷史有限、$c>0$ 且核局部可積，觀測歷史代入率後的有限窗積分可以有限。過程是否非爆炸還需要另行驗證其生成條件。例如有限規模上限使產能有界、$c>0$ 使時間核局部有界時，可以用線性出生率控制有限期間的增長；此時超臨界代表長期增長，並不等於有限時間內出現無限事件。

## C.3 直接後代超過親代的機率

給定始祖規模 $m_1$，其直接後代數 $K\sim\operatorname{Poisson}(\kappa(m_1))$。若每個後代獨立以機率 $\theta$ 超過某個目標規模，其保留數的生成函數為

$$\mathbb E[z^{K_>}]=\mathbb E[(1-\theta+\theta z)^K]
=\exp\{\kappa(m_1)\theta(z-1)\}.$$

因此 $K_>\sim\operatorname{Poisson}(\kappa(m_1)\theta)$，至少一個的機率為 $1-e^{-\kappa(m_1)\theta}$。未截斷 GR 且目標為親代本身時，$\theta=e^{-\beta(m_1-m_0)}$。

若限定直接後代在時間區間 $[a,b]$、空間區域 $\mathcal R$，再乘上對應核積分 $G$、$F$ 即可。這是單一親代的直接後代計數；多世代的完整目標事件機率不能只乘一個固定的平均放大倍率後套入 Poisson 公式。

## C.4 平均場家族的生成函數

這一節改用一個較簡單的模型：每個事件的後代數均為 $\operatorname{Poisson}(n)$，與規模無關。它是正文灰色參考線的模型，不是規模依賴產能的 ETAS。

令 $H(z)=\mathbb E[z^Y]$，始祖貢獻一個 $z$，每個後代再開啟一個同分布子家族，因此

$$H(z)=z\,\mathbb E[H(z)^K]=z e^{n(H(z)-1)}.$$

$n<1$ 時 $H'(1)=1/(1-n)$。由 Lagrange 反演，其家族大小是 Borel 分布：

$$P(Y=k)=\frac{e^{-nk}(nk)^{k-1}}{k!},\qquad k=1,2,\ldots.$$

令單個後代超過始祖規模的機率為 $\theta$。在這個規模與繁衍獨立的平均場模型中，沒有任何後代超過始祖的機率為 $H(1-\theta)/(1-\theta)$。設 $w=H(1-\theta)$，則

$$w=(1-\theta)e^{n(w-1)},\qquad
P_>=1-\frac{w}{1-\theta}.$$

固定 $n<1$、$\theta\to0$ 時，$P_>\simeq n\theta/(1-n)$。不能在這個近似下再任意令 $n\to1$；機率不會發散。

對連續規模分布抽出的始祖，$\theta$ 在 $(0,1)$ 均勻分布。當 $n\to1$，整體平均 $P_>$ 也不等於 1：單成員家族的機率為 $e^{-n}$，它們不可能出現更大後代，因此平均機率至多 $1-e^{-n}$。ETAS 的規模會影響家族大小，不能用此 Borel 模型作為其精確家族分布。

## C.5 分支法模擬與觀測窗

指定背景模型、規模分布、正規化觸發核及觀測範圍後，先生成背景事件；對每個待處理事件 $i$，抽取 $K_i\sim\operatorname{Poisson}(\kappa(m_i))$，再為各後代抽取相對時間、位置與規模，加入待處理佇列。處理到佇列為空或已宣告的計算限制為止。

Omori 時間密度的累積分布為 $G(t)=1-(1+t/c)^{1-p}$，反函數抽樣可寫成

$$\Delta t=c\left[(1-U)^{-1/(p-1)}-1\right],\qquad U\sim U(0,1),\quad p>1.$$

截斷規模的反函數為

$$M=m_0-\frac1\beta\ln\left[1-U(1-e^{-\beta\Delta})\right].$$

只有落在窗內的後代才列入輸出，但若允許窗外事件的後代回到區內，就不能把所有窗外生成事件直接丟棄。從有限起點開始、沒有暖機的模擬，也不是長期平穩樣本。家族大小上限是計算截斷，不是物理參數；應記錄是否觸及上限，而不能把未完成家族當成完整結果。

## C.6 有限時間與空間的概似積分

設評分窗為 $[S,T]\times\mathcal R$，使用對新事件規模積分後的 $\lambda^*_{ST}$。忽略與參數無關常數時，

$$\ell=\sum_{j:(t_j,z_j)\in[S,T]\times\mathcal R}\ln\lambda^*_{ST}(t_j,z_j)
-\int_S^T\int_{\mathcal R}\lambda^*_{ST}(t,z)\,dz\,dt.$$

窗前及鄰近區域的已知事件可以作為條件歷史，但不參與窗內事件求和。對 $t_i<T$，定義

$$a_i=\max(S,t_i)-t_i,\qquad b_i=T-t_i,$$

則正規化 Omori 的時間佔比為

$$G_i=(1+a_i/c)^{1-p}-(1+b_i/c)^{1-p}.$$

只有 $t_i\ge S$ 時才退化成 $1-(1+(T-t_i)/c)^{1-p}$。較早的事件不能漏掉第一個下限項。另令

$$F_i=\int_{\mathcal R}f(z-z_i;m_i)\,dz,$$

固定背景下的積分為

$$(T-S)\int_{\mathcal R}\mu(z)\,dz+\sum_{i:t_i<T}\kappa(m_i)G_iF_i.$$

$\kappa(m_i)G_iF_i$ 是沿實際歷史對單個來源累積的窗內貢獻；整個補償子通常是依歷史變動的隨機量，不能無條件稱為預報起點已知的未來期望數。

空間核若以來源為圓心，半徑 $R$ 的圓內積分為

$$F(R)=1-(1+R^2/\sigma)^{1-q},\qquad\sigma=De^{\gamma(m_i-m_0)}.$$

若區域相對來源是星形、每條射線自 0 到 $R(\varphi)$ 都在區域內，可再平均方位角得到 $F$。任意多邊形或區外來源可能有多段進出邊界，不能只放一個 $R(\varphi)$；須將每段徑向積分相加，或採其他可靠的區域積分法。核在 $\sigma>0$ 時原點有限，不需宣稱 Jacobian 消除了不存在的奇異點。

## C.7 未知親代與 EM

令潛在標籤 $Z_j=0$ 表示背景，$Z_j=i$ 表示事件 $i$ 的直接後代。將事件 $j$ 附近的率貢獻除以總率，便得到 {eq}`eq:rho-phi`。共同的新事件規模密度在分子分母消去。

給定參數 $\theta^{(k)}$ 的 E 步計算 $\phi_j^{(k)},\rho_{ij}^{(k)}$。M 步最大化期望完整資料對數概似，其結構為

$$\begin{aligned}
Q(\theta\mid\theta^{(k)})={}&\sum_j\phi_j^{(k)}\ln\mu(z_j)\\
&+\sum_j\sum_{i:t_i<t_j}\rho_{ij}^{(k)}\ln\{\kappa(m_i)g(t_j-t_i)f(z_j-z_i;m_i)\}\\
&-\int_S^T\int_{\mathcal R}\lambda^*_{ST,\theta}(t,z)\,dz\,dt.
\end{aligned}$$

這裡省略可分離的規模項，並假設分支模型的條件歷史處理一致。精確 E 步與使 $Q$ 不下降的 M 步，可保證觀測概似不下降；不保證全域最大。Veen 與 Schoenberg（2008）的 EM 型方法應連同其邊界與計算近似閱讀。

Zhuang et al.（2002）的背景估計包含加權變頻寬平滑。若任意加入頻寬重選、截斷或固定次數迭代，就不能直接引用精確 EM 的單調保證。正文地圖固定外部觸發參數，只示範背景權重更新，並未最大化上述完整 $Q$。

## C.8 參數補償與資訊矩陣

先考慮一個可分離的簡化模型：$\ell(\eta)=\sum_i\ell_i(u_i)$，$u_i=\eta_1+\eta_2\Delta m_i$。例如已知配對的空間核，可用 $\eta_1=\ln D$、$\eta_2=\gamma$。令 $w_i=-\mathbb E[\partial^2\ell_i/\partial u_i^2]>0$、$W=\sum_iw_i$，則

$$I=\begin{pmatrix}\sum_iw_i&\sum_iw_i\Delta m_i\\
\sum_iw_i\Delta m_i&\sum_iw_i\Delta m_i^2\end{pmatrix},\qquad
\det I=W^2\operatorname{Var}_w(\Delta m).$$

因此 $I^{-1}_{22}=1/[W\operatorname{Var}_w(\Delta m)]$。當規模完全相同，參數只能透過一個線性組合出現；規模跨度與總資訊量都會影響估計精度。真實 ETAS 的親代未知、混合強度與其他參數會引入更多交叉項，不能將此簡化矩陣宣稱為完整 ETAS 的精確資訊矩陣。

對負對數概似 $\xi$ 使用 DFP 時，令 $s_k=\theta_{k+1}-\theta_k$、$y_k=\nabla\xi(\theta_{k+1})-\nabla\xi(\theta_k)$，反 Hessian 近似更新為

$$H_{k+1}=H_k+\frac{s_ks_k^{\mathsf T}}{s_k^{\mathsf T}y_k}
-\frac{H_ky_ky_k^{\mathsf T}H_k}{y_k^{\mathsf T}H_ky_k}.$$

正定性需要合適初始矩陣與曲率條件。最佳化器內部的近似矩陣不一定足以作為統計共變異數；最後仍須檢查目標函數、座標轉換、參數邊界與局部二次近似。

## C.9 R–J 的受限對照與全世代平均

只保留一個指定主震 $M_m$，忽略背景與新事件再觸發，並採未截斷 GR，對目標規模 $M\ge m_0$ 積分，可得

$$R(t,M)=A(p-1)c^{p-1}e^{\alpha(M_m-m_0)}e^{-\beta(M-m_0)}(t+c)^{-p}.$$

若 $\alpha=\beta$，便是正文 {eq}`eq:rj-rate` 的形式，且在相同時間單位下

$$10^{a'}=A(p-1)c^{p-1}.$$

這個代數對照要求 $p>1$ 以使用上述正規化核，並只描述指定來源的直接後代。R–J 自身在有限窗可使用其他 $p$；實際對整段序列擬合得到的 $a',p$ 也不必等於直接觸發核的參數。

若恢復完整分支，給定始祖規模，其全世代後代期望數在次臨界、無窗截斷下為 $\kappa(M_m)/(1-n)$。時間率則是卷積之和：

$$r_{\rm all}(t\mid M_m)=\kappa(M_m)\{g(t)+n(g*g)(t)+n^2(g*g*g)(t)+\cdots\}.$$

每個卷積積分為 1，所以總積分具有 $1/(1-n)$ 放大因子；但時間形狀通常不再是原來的 $g$。因此不能用這個總數關係，把實際 R–J 產能反推為唯一 ETAS 分支比。

## C.10 後驗預報與診斷的條件

貝氏參數後驗為 $p(\theta\mid\mathcal D)\propto L(\theta;\mathcal D)p(\theta)$。對未來事件集合 $B$，後驗預報將參數不確定性平均：

$$P(N(B)\ge1\mid\mathcal D)=\int P_\theta(N(B)\ge1\mid H_{t_0})p(\theta\mid\mathcal D)\,d\theta.$$

內層是給定參數的條件預報，仍須包含模型允許的未來級聯。後驗只反映指定先驗、概似和觀測模型，沒有建模的漏測、機制變化不會自動被納入。

時間變換診斷使用補償子 $\tau_j=\int_0^{t_j}\lambda^*(u)du$。在正確強度與相應正則條件下，變換間隔應如單位指數分布；同資料估參數、選模型後，原始 KS 標準分布通常不再是正確校準，須透過模擬重估等方式處理。只檢查邊際分布也不足以確認間隔獨立。

## 參考資料與延伸閱讀

- Jalilian, A.（2019），[ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://doi.org/10.18637/jss.v088.c01)。[免費全文](https://www.jstatsoft.org/article/view/v088c01)，對照條件歷史、空間積分及估計流程。
- Reinhart, A.（2018），[A Review of Self-Exciting Spatio-Temporal Point Processes and Their Applications](https://doi.org/10.1214/17-STS629)。[免費作者預印本](https://arxiv.org/abs/1708.02647)，串起分支表示、邊界、估計與診斷。
- Daley, D. J. 與 Vere-Jones, D.（2003），[An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, Second Edition](https://doi.org/10.1007/b97277)。第 6 章群集與標記模型、第 7 章條件強度與概似為本附錄的理論背景；出版社全文可能需訂閱。
- Ogata, Y.（1988），[Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes](https://doi.org/10.1080/01621459.1988.10478560)。[免費研究機構全文](https://bemlar.ism.ac.jp/zhuang/Refs/Refs/ogata1988.pdf)，閱讀時間 ETAS 與模型殘差。
- Ogata, Y. 與 Zhuang, J.（2006），[Space–time ETAS models and an improved extension](https://doi.org/10.1016/j.tecto.2005.10.016)。[免費機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)，對照產能與空間尺度的分離；出版社可能需訂閱。
- Zhuang, J., Ogata, Y. 與 Vere-Jones, D.（2002），[Stochastic Declustering of Space-Time Earthquake Occurrences](https://doi.org/10.1198/016214502760046925)。背景加權估計及隨機親代的主要來源；出版社全文可能需訂閱。
- Veen, A. 與 Schoenberg, F. P.（2008），[Estimation of Space–Time Branching Process Models in Seismology Using an EM–Type Algorithm](https://doi.org/10.1198/016214508000000148)。未知親代如何形成不完整資料問題，以及 EM 型估計的計算方式；出版社全文可能需訂閱。
- Naylor, M., Serafini, F., Lindgren, F. 與 Main, I. G.（2023），[Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759)。免費開放全文，用合成案例理解聯合後驗與漏測影響。
- Reasenberg, P. A. 與 Jones, L. M.（1989），[Earthquake Hazard After a Mainshock in California](https://doi.org/10.1126/science.243.4895.1173)。R–J 餘震預報的原始來源；出版社全文可能需訂閱。

- [Seismicity Analysis through Point-process Modeling: A Review](https://doi.org/10.1007/s000240050275) — Yosihiko Ogata，1999，*Pure and Applied Geophysics*，155:471–507。從條件強度一路連到概似、ETAS 與殘差，適合把本章放回統計地震學的發展脈絡；出版社全文可能需訂閱。

返回 {doc}`13_etas_structure` 或 {doc}`14_etas_estimation`。
