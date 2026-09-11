
# 附錄 E：預報檢驗與模型組合的推導

本附錄需要機率分布、微分與期望值。對應 {doc}`一致性檢驗 <17_testing_consistency>`、{doc}`模型比較 <18_testing_comparison>`及 {doc}`模型組合 <19_ensembles>`，各節推導後可回到主文的圖例。

## E.1 從 Poisson 計數到條件形狀

獨立箱計數 $\omega_j\sim\mathrm{Poisson}(\Lambda_j)$ 的乘積機率取對數，即得 {eq}`eq:jpoll`。總數 $N=\sum_j\omega_j$ 的均值為 $\Lambda=\sum_j\Lambda_j$。條件在總數 $N=n$ 上：

$$P(\boldsymbol\omega\mid N=n)
=\frac{n!}{\prod_j\omega_j!}\prod_j\left(\frac{\Lambda_j}{\Lambda}\right)^{\omega_j},
\qquad \sum_j\omega_j=n.$$

這是多項分布，率的整體尺度消去，只剩形狀。S 與 M 分別邊際化到空間和規模；cL 保留聯合形狀。模擬時固定 $n$，按正規化率分配事件，再對模擬與觀測使用相同統計量。

若在固定 Poisson 點過程下使用很細的等體積箱，每箱最多一個事件，$\Lambda_j\approx\lambda(x_j)\Delta V$，可得

$$\ln L_{\rm grid}\approx\sum_{i=1}^N\ln\lambda(x_i)-\int\lambda(x)\,dx
+N\ln\Delta V.$$

同網格兩模型的最後一項抵消，但除以 $N$ 仍留下 $\ln\Delta V$。因此每事件分數也不能直接消除跨網格尺度。一般歷史相依強度的連續概似使用順序條件化，不能由一組無條件獨立空間箱原封不動推出。

## E.2 叢集的變異數：一個明確的模型

令背景叢集數 $B\sim\mathrm{Poisson}(\mu_B)$，各叢集總大小 $G_i\ge1$ 獨立同分布，且與 $B$ 獨立。若 $N=\sum_{i=1}^B G_i$，全期望與全變異數給出

$$\mathbb E[N]=\mu_B\mathbb E[G],\qquad
\mathrm{Var}(N)=\mu_B\mathbb E[G^2].$$

所以 $\mathrm{Var}(N)/\mathbb E[N]=\mathbb E[G^2]/\mathbb E[G]\ge1$，只要 $P(G>1)>0$ 就嚴格大於一。這是完整獨立叢集模型的結果；有限時空窗截斷叢集時需要重新計算。

本小節以 $g=\mathbb E[G]$、$v=\mathrm{Var}(G)$ 表示完整家族大小的平均與變異數。
再假設每事件直接後代數 $K\sim\mathrm{Poisson}(n)$、$n<1$，利用 $G=1+\sum_{i=1}^{K}G_i$：

$$g=1+ng\Rightarrow g=\frac1{1-n},\qquad
v=nv+ng^2\Rightarrow v=\frac{n}{(1-n)^3}.$$

故完整叢集計數的 Fano 比為 $1/(1-n)^2$。真實帶規模產能的 ETAS 直接後代數可能是混合 Poisson，並非每事件同均值的 Poisson；不能只憑分支比就將此數值套到任何有限窗 ETAS 計數。

## E.3 負二項與二元分數

採負二項參數 $r>0,0<q<1$。這一節的 $\mu$、$v$ 分別表示計數 $N$ 的平均與變異數，

$$P(N=n)=\frac{\Gamma(n+r)}{\Gamma(r)n!}q^r(1-q)^n,
\quad \mu=\frac{r(1-q)}q,\quad v=\frac{r(1-q)}{q^2}.$$

若 $v>\mu>0$，則 $q=\mu/v$、$r=\mu^2/(v-\mu)$；$v\downarrow\mu$ 時趨近 Poisson。這只指定邊際分布，沒有自動指定箱間聯合分布。

若 Poisson 計數均值為 $\Lambda$，至少一次機率 $P=1-e^{-\Lambda}$。二元分數 BILL 為 $y\ln P+(1-y)\ln(1-P)$。$n=0$ 時與 POLL 都為 $-\Lambda$；$n\ge1$ 時兩者差為

$$\begin{aligned}
\mathrm{POLL}-\mathrm{BILL}
&=-\Lambda+n\ln\Lambda-\ln(n!)-\ln(1-e^{-\Lambda})\\
&=(n-1)\ln\Lambda-\ln(n!)-\frac\Lambda2+O(\Lambda^2).
\end{aligned}$$

最後一步使用 $(1-e^{-\Lambda})/\Lambda=1-\Lambda/2+O(\Lambda^2)$。因此稀少單事件箱的差異小，重複事件箱的差異可能很大。一般計數分布的至少一次機率應直接用 $1-P(N=0)$，不能只看均值。

## E.4 模擬分位數的誤差

固定觀測統計量 $s_{\rm obs}$，獨立模擬 $R$ 次，令 $I_r=1\{s_r\le s_{\rm obs}\}$，$\hat\gamma=R^{-1}\sum_r I_r$。則

$$R\hat\gamma\sim\mathrm{Binomial}(R,\gamma),\qquad
\mathrm{SE}(\hat\gamma)=\sqrt{\frac{\gamma(1-\gamma)}R}.$$

在 $\gamma=0.05$ 時，$R=10^4$ 對應標準誤約 $0.00218$。這只是模擬誤差，沒有包括模型參數不確定性。尾端或零次計數時應使用合適的二項區間；正式 Monte Carlo 檢定可預先採含觀測秩的加一修正，並處理 ties。

多重檢驗若採 Bonferroni，以檢驗家族內實際數目 $J$ 使用門檻 $\alpha/J$，不要求獨立。把相關檢驗直接折成較小的「有效數目」不保證同樣控制；需另行校準聯合拒絕機率。

## E.5 資訊增益與區間

同一目錄的階乘項相消，網格概似差為

$$\ln L_Z-\ln L_1=(\Lambda_1-\Lambda_Z)
+\sum_{j,k}\omega_{jk}\ln\frac{\Lambda^Z_{jk}}{\Lambda^1_{jk}}.$$

每個事件所在箱重複一次，把後項改寫為 $\sum_nd_n$，除以 $N>0$ 即得 {eq}`eq:igpe`。若 $N=0$，仍有總分差，但每事件資訊增益未定義。

樣本變異數可寫成

$$s_d^2=\frac{\sum_n(d_n-\bar d)^2}{N-1}
=\frac{\sum_nd_n^2}{N-1}-\frac{(\sum_nd_n)^2}{N(N-1)}.$$

將總量差視為固定平移、事件貢獻近似獨立時，常見區間為 $\mathrm{IGPE}\pm t_{N-1,1-\alpha/2}s_d/\sqrt N$。事件相關會改變均值變異數：

$$\mathrm{Var}(\bar d)=\frac1{N^2}\left[\sum_n\mathrm{Var}(d_n)
+2\sum_{i<j}\mathrm{Cov}(d_i,d_j)\right].$$

所以叢集中的區間可能需要區塊或模型式校準。未拒絕零增益不等於分數恰好相同，也不等於通過等價性檢定。

AICc 型修正採 $\mathrm{AICc}=-2\ln L+2k+2k(k+1)/(N-k-1)$，兩模型相減、除以 $2N$ 即得 {eq}`eq:igpec`。分母要求 $N>k+1$；此為有條件近似，非所有相依點過程的通用修正。

## E.6 警報面積與校準

令 $\nu(\tau)$ 是警報比例 $\tau$ 下的漏報率，常用面積摘要為

$$\mathrm{AS}(\tau)=\frac1\tau\int_0^\tau[1-\nu(u)]\,du,\quad \tau>0.$$

若在所選參考測度下隨機警報，$\nu=1-\tau$，故 $\mathrm{AS}(\tau)=\tau/2$；只有完整範圍 $\tau=1$ 時基準為 $1/2$。等面積與背景率加權使用不同參考測度。

Brier 分解在預報實際取離散值 $P_k$ 的分組下精確成立。第 $k$ 組大小 $n_k$、觀測比例 $\bar y_k$，$\bar y=\sum_kn_k\bar y_k/N$，先展開

$$\frac1{n_k}\sum_{i\in k}(P_k-y_i)^2
=(P_k-\bar y_k)^2+\bar y_k(1-\bar y_k).$$

交叉項因組內殘差和為零而消去。再用全變異數，得

$$S_B=\underbrace{\sum_k\frac{n_k}N(P_k-\bar y_k)^2}_{\mathrm{REL}}
-\underbrace{\sum_k\frac{n_k}N(\bar y_k-\bar y)^2}_{\mathrm{RES}}
+\underbrace{\bar y(1-\bar y)}_{\mathrm{UNC}}.$$

若只是將不同連續預報任意分箱，以箱平均代入不一定保持精確分解。對真實機率 $Q$，期望單次 Brier 分數為 $Q(P-1)^2+(1-Q)P^2=(P-Q)^2+Q(1-Q)$，唯一最小值在 $P=Q$，這就是嚴格適當性。

## E.7 成本與損失

假設行動成本 $C$、不行動且事件發生損失 $L$、行動可完全避免該損失，$0<\alpha=C/L<1$。已知事件機率 $P$ 時，比較 $C$ 與 $PL$，得到門檻 $P>\alpha$。

若事件比例為 $s$，警報比例 $\tau$，漏報比例 $\nu$，每次平均支出除以 $L$ 為 $e_f=\alpha\tau+s\nu$。只用長期比例的最佳基準為 $e_c=\min(\alpha,s)$；完美預報為 $e_p=\alpha s$。因此相對價值

$$V=\frac{\min(\alpha,s)-(\alpha\tau+s\nu)}{\min(\alpha,s)-\alpha s},\qquad 0<s<1.$$

此式依賴相同成本、完全防護及可及時行動等假設，並不推出所有預報的價值峰值必在 $\alpha=s$。實際行動門檻改變時，$\tau,\nu$ 也會改變。

## E.8 凸組合的凹性與端點

固定率函數 $\lambda_1,\lambda_2$，令 $\lambda_r=(1-r)\lambda_1+r\lambda_2$，在事件處率皆正，總積分為 $\Lambda_r=(1-r)\Lambda_1+r\Lambda_2$。對數概似的導數是

$$\ell'(r)=\sum_n\frac{\lambda_2(n)-\lambda_1(n)}{\lambda_r(n)}-(\Lambda_2-\Lambda_1),$$

$$\ell''(r)=-\sum_n\frac{[\lambda_2(n)-\lambda_1(n)]^2}{\lambda_r(n)^2}\le0.$$

至少一個觀測事件處兩率不同時才嚴格凹。兩函數可在未觀測處不同，卻在所有事件處相同，故函數不相同不足以保證嚴格性。

嚴格凹情形有唯一最大值；在 $(0,1)$ 內的條件為 $\ell'(0)>0$ 且 $\ell'(1)<0$。若 $\Lambda_1=\Lambda_2$，令 $\rho_n=\lambda_2(n)/\lambda_1(n)$，條件簡化成

$$\overline\rho>1,\qquad\overline{\rho^{-1}}>1.$$

若樣本平均 $\overline{\ln\rho}=0$ 恰成立，且事件處比值非常數，AM–GM 可推出兩條件。但統計檢驗未顯著不代表這個精確等式。反例是所有事件處 $\rho_n=2$ 且總量相同，此時概似隨 $r$ 增加，最佳值在端點 $r=1$。

Jensen 另給出 $\ell(r)\ge(1-r)\ell(0)+r\ell(1)$；這是低於最佳成分與否都可能的下界，不是保證超過最佳成分。

## E.9 乘法修正的 KL 表達

將基準與修正後空間分布正規化為 $w_j,q_j$，且 $q_j=w_jm_j$、$\sum_jq_j=1$。任意真實分布 $p$ 下，

$$\sum_jp_j\ln m_j
=\sum_jp_j\ln\frac{q_j}{w_j}
=D_{\rm KL}(p\|w)-D_{\rm KL}(p\|q).$$

若 $p=w$，增益為 $-D_{\rm KL}(w\|q)\le0$，表示對本來正確的基準做不必要修正會受損；不表示任何乘法修正都會系統性失敗。

加法若保留 $q_j\ge\pi_1w_j$、$\pi_1>0$，則

$$D_{\rm KL}(p\|q)\le D_{\rm KL}(p\|w)+\ln(1/\pi_1).$$

它限制相對基準的損失，但不保證候選基準本身足夠好。各式需在所比較支撐上有正率或採適當無窮大慣例。

## E.10 混合分布、平均率與 logistic

先抽模型指標 $I$，$P(I=i)=\pi_i$，再抽 $N\mid I=i\sim\mathrm{Poisson}(\Lambda_i)$。全變異數給出

$$\mathrm{Var}(N)=\sum_i\pi_i\Lambda_i+\sum_i\pi_i(\Lambda_i-\bar\Lambda)^2,
\qquad\bar\Lambda=\sum_i\pi_i\Lambda_i.$$

相較之下，$N\sim\mathrm{Poisson}(\bar\Lambda)$ 只有第一項。混合 Poisson 通常不是負二項；以負二項近似還需額外理由。

若 logistic 迴歸為 $\ln[P/(1-P)]=b_0+\sum_ib_i\ln\phi_i$，直接取指數得到

$$\frac{P}{1-P}=e^{b_0}\prod_i\phi_i^{b_i}.$$

這是勝算的乘冪關係，不是點過程率的等式。由係數再建立非負率權重，需要特定轉換與獨立檢驗；截去負係數或移除截距是方法設計，並非 logistic 理論必然要求。

## 參考資料與延伸閱讀

- [Theory of CSEP Tests](https://docs.cseptesting.org/getting_started/theory.html) — pyCSEP 開發團隊，官方文件（免費）。先對照各檢驗的目標、模擬方式與分位數分數，再看 N、S、M 與 conditional L-test 的程式範例，可把本章公式接到實際檢驗流程。

- [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。從「觀測是否與預報一致」出發整理檢驗方法，適合先建立觀念，再回頭理解本章為何要拆開事件數、空間與規模。

- [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。本章 Poisson／二元似然、負二項計數與多重檢定的主要實例來源；閱讀時留意餘震叢集如何影響各種評分。

- [An Introduction to the Theory of Point Processes, Volume I: Elementary Theory and Methods, 2nd ed.](https://doi.org/10.1007/b97277) — D. J. Daley、D. Vere-Jones（2003），Springer（全文可能需訂閱）。第 2 章建立 Poisson 過程，第 4 章說明更新過程與等待時間，第 6–7 章連結群集、條件強度、概似與補償子；可按本章主題選讀。

- [Theory of CSEP Tests：Forecast comparison tests](https://docs.cseptesting.org/getting_started/theory.html#forecast-comparison-tests) — pyCSEP 開發團隊，官方文件（免費）。從資訊增益 IGPE 的定義與 T-test 範例開始，對照本章的基準模型、率修正與信賴區間，理解「分數較高」與「差異顯著」的區別。

- [Evaluating earthquake predictions and earthquake forecasts: a guide for students and new researchers](https://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/zechar.pdf_2063069264.pdf) — J. Douglas Zechar（2010），CORSSA 免費教材（DOI：10.5078/corssa-77337879）。接著閱讀警報式預測與誤差圖的部分，比較不同評估方法回答的問題，銜接本章的 Molchan 圖與面積技能分數。

- [Statistical power of spatial earthquake forecast tests](https://doi.org/10.1093/gji/ggad030) — Asim M. Khawaja 等（2023），*Geophysical Journal International*；[免費機構典藏全文](https://gfzpublic.gfz.de/pubman/item/item_5015770_1)。以空間檢驗說明樣本量與網格如何改變統計功效，正好延伸本章「透過檢驗不等於模型有辨識力」的討論；可接著比較等寬網格與 Quadtree 的設計。

- [Enhancing the Statistical Evaluation of Earthquake Forecasts—An Application to Italy](https://doi.org/10.1785/0220240209) — Jonas R. Brehmer、Kristof Kraus、Tilmann Gneiting、Marcus Herrmann、Warner Marzocchi（2025），*Seismological Research Letters*（出版社全文可能需訂閱；[免費作者稿](https://arxiv.org/abs/2405.10712)）。從義大利預報案例延伸本章的模型比較、可靠度圖與校準，重點是評分函數如何對應預報目標，以及如何分開檢查校準與鑑別力。

- [Developing, Testing, and Communicating Earthquake Forecasts: Current Practices and Future Directions](https://doi.org/10.1029/2023RG000823) — Leila Mizrahi 等（2024），*Reviews of Geophysics*（免費全文）。先讀模型發展與各國作業系統中組合模型的段落，瞭解為何需要結合不同預報，以及權重選擇如何與檢驗制度連在一起。

- [Mixture Models for Improved Short-Term Earthquake Forecasting](https://doi.org/10.1785/0120080063) — David A. Rhoades、Matthew C. Gerstenberger（2009），*Bulletin of the Seismological Society of America*（全文可能需訂閱；[SCEC 免費摘要](https://central.scec.org/node/3964)）。本章加法組合的原始案例，將 STEP 與 EEPAS 的不同時間尺度結合；重點是凸組合如何利用互補資訊，以及回溯估出的權重仍須接受前瞻檢驗。

- [Prospective evaluation of multiplicative hybrid earthquake forecasting models in California](https://doi.org/10.1093/gji/ggac018) — J. A. Bayona、W. H. Savran、D. A. Rhoades、M. J. Werner（2022），*Geophysical Journal International*（免費全文）。對照本章乘法 hybrid 的前瞻測試，觀察回溯資訊增益為何未必延續到未來；適合與前一篇的組合動機一起讀。

- [Maximizing the forecasting skill of an ensemble model](https://doi.org/10.1093/gji/ggad020) — Marcus Herrmann、Warner Marzocchi（2023），*Geophysical Journal International*（免費全文）。本章 logistic 權重學習的主要來源；重點是直接最佳化整體組合的技巧，而非只依各成分模型的單獨成績配權重。

- [Regional Earthquake Likelihood Models II: Information Gains of Multiplicative Hybrids](https://doi.org/10.1785/0120140035) — D. A. Rhoades、M. C. Gerstenberger、A. Christophersen、J. D. Zechar、D. Schorlemmer、M. J. Werner、T. H. Jordan（2014），*Bulletin of the Seismological Society of America*（出版社全文可能需訂閱；[Bristol 大學免費全文](https://research-information.bris.ac.uk/files/49915043/Rhoades_Hybrids_RELM_BSSA_2014.pdf)）。本章乘法 hybrid、保序轉換與正規化的直接來源；讀完建構方法後，對照上面的 Bayona 等人前瞻測試，區分擬合改善與未來預報表現。

- [Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package](https://doi.org/10.3389/fams.2023.1126759) — Mark Naylor、Francesco Serafini、Finn Lindgren、Ian G. Main（2023），Frontiers in Applied Mathematics and Statistics（免費全文）。對照模型間差異與單一模型內的參數不確定性；貝氏 ETAS 的後驗分布並不等同於模型組合權重。
