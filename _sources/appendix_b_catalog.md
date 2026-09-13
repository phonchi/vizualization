# 附錄 B：目錄、完整度與叢集律

**對應主文章**：
第 2 章 {doc}`02_reading_horus`；
第 3 章 {doc}`03_experiment_spec`；
第 7 章 {doc}`07_gr_bvalue`；
第 8 章 {doc}`08_completeness`；
第 10 章 {doc}`10_clustering_laws`。

本附錄補上資料處理與估計推導。
門檻、格距與觀察窗都是模型條件。
先備知識為機率密度、微分與概似。

## B.0 HORUS 取得、快取與 CELLE 反投影

HORUS 是義大利均一化儀器地震目錄。
英文為 Homogenized Instrumental Seismic Catalogue。
均一化指把不同規模資料轉到共同尺度。
本實驗使用矩規模 $M_w$。
它以地震矩衡量地震大小。

資料入口是 [PyEEPAS 公開程式庫](https://github.com/phonchi/EEPAS)。
其 [data/README](https://github.com/phonchi/EEPAS/blob/main/data/README.md)
指向 [Google Drive 資料夾](https://drive.google.com/drive/folders/170WNb8M8PQJDX1B2JSQYfqj4ad80TC0U)。
公開程式是資料入口，不是目錄的原始作者。
目錄來源應引用 Lolli 等（2020）。

`scripts/fetch_italy_data.py` 負責下載。
它也建立清理後的文字快取。
快取（cache）是供重複讀取的處理結果。
下表說明各檔案用途。

| 位置 | 用途 |
|---|---|
| `data/cache/italy/HORUS_Ita_Catalog.txt` | 原始目錄 |
| `data/cache/italy/HORUS_Ita_DataOrigin.txt` | 資料來源說明 |
| `data/cache/italy/README.txt` | 目錄版本與欄位說明 |
| `data/cache/italy/PROVENANCE.txt` | 來源網址與下載日期 |
| `data/cache/italy/horus_clean.csv` | 清理後的快取 |
| `EEPAS/data/CPTI15.mat` | 收集區多邊形 |
| `EEPAS/data/CELLE_ter.mat` | 177 個測試方格 |

`load_horus()` 讀完整目錄。
`experiment_catalog()` 套用實驗篩選。
它保留收集區內、深度不超過 40 km 的地震。
事件類型旗標用來排除非地震事件。
收集區目前依 `Geo-CPTI15` 旗標篩選。
`in_collection_region()` 可另做幾何對照。
兩種方法不應未經核對就視為相同。

測試區 $R$ 用於計數與評分。
收集區 $S$ 提供區內外的歷史事件。
`target_events()` 再選指定期間的目標事件。
原始規模、取整規模與事件識別碼都應保留。

### 規模取整與三個門檻

取整（rounding）把規模映到固定格點。
依 Biondini（2023）式 2，規則為

$$m_{\rm binned}=\frac{\operatorname{int}(10m_{\rm raw}+0.5)}{10}.$$

文字式寫作 `int(m×10+0.5)/10`。
本式用於本實驗的正規模資料。
程式欄位 `mw` 保留原值，`mb` 保存取整值。
名目 5.0 的有效下緣是 4.95。
名目 2.5 的有效下緣是 2.45。

$m_0=2.45$ 是模型輸入下緣。
$m_T=5.0$ 是取整後的目標門檻。
$M_c$ 是目錄完整度，須由資料評估。
比較門檻前，須先統一格點與箱邊界。

2024 年版 HORUS 的測試期目標為 25 顆。
已發表論文的同期間為 27 顆。
學習期則同為 27 顆。
這些是專案已記錄的版本核對結果。
本附錄沒有重新計數。
差異不能解讀為模型多漏報了兩顆。

### CELLE 座標與反投影

反投影（inverse projection）把平面座標轉回經緯度。
`CELLE_ter.mat` 的 `CELLESD` 是投影公里座標。
前四欄為 $x_0,x_1,y_0,y_1$。
方格邊長為 $30\sqrt2$ km。
它們不是經緯度邊界。

本站採 EPSG:7794。
名稱為 RDN2008 / Italy zone (E-N)。
EPSG 原生長度單位為公尺。
檔案卻以公里保存，因此要先乘 1,000。
再轉成 EPSG:4326 的經緯度。
經度在前，緯度在後。

`km_to_lonlat()` 執行上述換算。
`lonlat_to_km()` 則先投影，再除以 1,000。
`testing_cells()` 回傳格心與角點經緯度。
`cell_of()` 仍以投影邊界判定格號。
其規則為左閉右開、下閉上開。

公開 README 的 CELLE 經緯度說明有誤。
本文依專案已核定座標與資料層說明。
經緯度外接框也不等於原投影方格。
重新分格見 {doc}`appendix_e_testing`。

177 格的選區曾使用截至 2021 年的資訊。
因此本實驗是既有研究區的重演。
不能稱所有區域設計都在測試前凍結。

## B.1 從 GR 累積數到規模密度

GR 是 Gutenberg–Richter 規模頻率律。
令固定期間的期望累積數為 $10^{a-bm}$。
除以門檻以上總數，得到存活機率。
存活機率表示規模仍超過指定值。

$$\begin{aligned}
P(M\ge m\mid M\ge m_0)&=10^{-b(m-m_0)}=e^{-\beta(m-m_0)},\\
s(m)&=\beta e^{-\beta(m-m_0)},\qquad \beta=b\ln10.
\end{aligned}$$

密度由存活機率取負導數得到。
規模上限有限時，還須重新正規化。
詳見 {doc}`appendix_c_etas`。

若兩尺度滿足 $M_2=a_sM_1+d_s$，且 $a_s>0$，

$$a_1-b_1M_1=a_1+\frac{b_1d_s}{a_s}-\frac{b_1}{a_s}M_2.$$

因此 $b_2=b_1/a_s$。
這只是確定轉換的代數結果。
實際換算还須納入量測與迴歸誤差。

## B.2 連續規模的概似與有限樣本偏差

令 $X_i=m_i-m_0$ 為獨立指數樣本。
完整度與門檻先固定，樣本數為 $N$。
最大概似估計（MLE）給出

$$\begin{aligned}
\ell(\beta)&=N\ln\beta-\beta\sum_iX_i,\\
\widehat\beta&=\frac{N}{\sum_iX_i},\qquad
\widehat b=\frac{\widehat\beta}{\ln10}.
\end{aligned}$$

令 $S_N=\sum_iX_i$。
它服從形狀 $N$、率 $\beta$ 的 Gamma 分布。
對 $N>1$，積分可得

$$\mathbb E[S_N^{-1}]=\frac\beta{N-1},\qquad
\mathbb E[\widehat\beta]=\frac N{N-1}\beta.$$

因此 $(N-1)\widehat\beta/N$ 在此模型下無偏。
門檻由同資料挑選時，不能直接沿用。
離散規模與漏測也會改變結論。

標準誤描述估計量的抽樣散布。
令規模樣本標準差為 $s_M$。
以一階誤差傳播近似，

$$\operatorname{SE}(\widehat b)\simeq
\ln10\,\widehat b^2\frac{s_M}{\sqrt N}\simeq\frac b{\sqrt N}.$$

最後一步要求理想指數模型成立。
它不包含目錄系統誤差。

## B.3 離散 GR 與精確最大概似

令最低格點為 $m_\ell$，格距為 $\Delta m$。
$m_\ell$ 與有效下緣 $m_0$ 須分開。
規模格點為 $m_\ell+j\Delta m$。
離散 GR 等同幾何分布：

$$P(J=j)=(1-r)r^j,\qquad r=e^{-\beta\Delta m}.$$

令 $u=\bar m-m_\ell$，最大化概似得到

$$\begin{aligned}
\ell(r)&=N\ln(1-r)+\sum_i j_i\ln r,\\
\widehat r&=\frac{\bar j}{1+\bar j},\\
\widehat b&=\frac1{\Delta m\ln10}\ln\left(1+\frac{\Delta m}{u}\right).
\end{aligned}$$

「精確」指這個離散模型的最大值。
它不保證有限樣本無偏。
$u=0$ 時，只能得到 $b\to\infty$ 的邊界解。
混合不同取整規則時，須重建區間概似。

## B.4 半格近似與不對稱區間

令 $\delta=\Delta m/2$。
Utsu 半格近似使用 $1/[(u+\delta)\ln10]$。
展開兩種公式，得到

$$\begin{aligned}
b_{\rm exact}&=\frac1{\ln10}\left(\frac1u-\frac\delta{u^2}+\frac{4\delta^2}{3u^3}+\cdots\right),\\
b_{\rm Utsu}&=\frac1{\ln10}\left(\frac1u-\frac\delta{u^2}+\frac{\delta^2}{u^3}+\cdots\right).
\end{aligned}$$

級數要求 $2\delta/u<1$。
兩式差異從二階項開始。
此差異不等同估計量對真值的偏差。

$b(u)$ 遞減且向上凸。
因此對稱的 $u$ 區間會轉成不對稱的 $b$ 區間。
Tinti 與 Gasperini（2024）給出近似端點。
令 $C=1+2\delta/u$，則

$$\begin{aligned}
b_1&=\frac1{2\delta\ln10}\ln\frac{C+\sqrt{C/N}}{1+\sqrt{C/N}},\\
b_2&=\frac1{2\delta\ln10}\ln\frac{C-\sqrt{C/N}}{1-\sqrt{C/N}}.
\end{aligned}$$

端點要求 $N>C$，並依賴漸近近似。
小樣本不能宣稱精確涵蓋率。
Bootstrap 也須配合資料相依性。
其定義見 {doc}`appendix_a_point_process`。

## B.5 規模差、修剪與相關性

修剪（trimming）是按預定規則刪除部分樣本。
規模差方法常先保留超過某門檻的差值。
若 $X,Y$ 獨立且同為指數分布，

$$f_{X-Y}(d)=\frac\beta2e^{-\beta|d|}.$$

這是 Laplace 雙指數分布。
條件於 $X-Y>d_0\ge0$，超額仍為指數分布。
漏測與配對選擇可能破壞這個性質。

離散樣本的差值則滿足

$$P(J_1-J_2=j)=\frac{1-r}{1+r}r^{|j|}.$$

零差值具有獨立的機率質量。
刪去它後，非零絕對差是移位幾何分布。

相鄰差共用一個觀測，因而相依。
令 $D_i=X_{i+1}-X_i$，則

$$\operatorname{Cov}(D_i,D_{i+1})=-\beta^{-2},\qquad
\operatorname{Var}(D_i)=2\beta^{-2}.$$

符號差的相關係數是 $-1/2$。
絕對差 $A_i=|D_i|$ 卻滿足

$$\operatorname{Cov}(A_i,A_{i+1})=\frac1{3\beta^2},\qquad
\operatorname{Var}(A_i)=\beta^{-2}.$$

其相關係數為正的 $1/3$。
取絕對值後，不能沿用負相關結論。
進一步修剪也會改變有效抽樣結構。

## B.6 Omori 核：有限窗與無限時間

Omori–Utsu 律描述餘震率隨時間衰減。
令 $\tau$ 為距主震的時間，$c>0$。
有限窗 $[a,b]$ 的積分為

$$I(a,b;c,p)=\int_a^b(\tau+c)^{-p}\,d\tau
=\begin{cases}
\dfrac{(b+c)^{1-p}-(a+c)^{1-p}}{1-p},&p\ne1,\\
\ln\dfrac{b+c}{a+c},&p=1.
\end{cases}$$

兩式在 $p\to1$ 時相接。
有限窗允许 $p\le1$。
無限時間正規化才要求 $p>1$。
此時積分為 $c^{1-p}/(p-1)$。
除以此值便得到主文的密度 $g$。

若直接設定觸發期限 $T_c$，則

$$g_{T_c}(\tau)=\frac{(\tau+c)^{-p}}{I(0,T_c;c,p)}
\mathbf1\{0\le\tau\le T_c\}.$$

這已改變觸發模型。
它不同於只裁切觀察窗。
固定未正規化係數時，總產能會隨期限改變。

## B.7 Omori 的剖面概似

剖面概似（profile likelihood）先消去部分參數。
令未正規化率為 $K_O(\tau+c)^{-p}$。
對固定窗內的 $N>0$ 個事件，

$$\begin{aligned}
\ell(K_O,c,p)&=N\ln K_O-p\sum_i\ln(\tau_i+c)-K_OI,\\
\widehat K_O(c,p)&=N/I,\\
\ell_{\rm prof}(c,p)&=N\ln(N/I)-p\sum_i\ln(\tau_i+c)-N.
\end{aligned}$$

原本三維搜尋降為兩維。
此推導假設事件是非齊次 Poisson。
多世代觸發與漏測可能讓此假設失效。

Hessian 是目標函數的二階導數矩陣。
用它估區間時，最優點須適合局部二次近似。
矩陣非正定時，區間可能無法辨認。
對負變異數取絕對值沒有統計依據。

## B.8 空間核的尺度與正規化

各向同性（isotropic）表示不依方位角改變。
令面積尺度為 $v(m)>0$，尾指數為 $q>1$。

$$h(r;m)=\frac{q-1}{\pi v(m)}[1+r^2/v(m)]^{-q}.$$

令 $z=1+r^2/v$，極座標積分給出

$$\int_{\mathbb R^2}h\,dx\,dy=(q-1)\int_1^\infty z^{-q}\,dz=1.$$

若 $v(m)=v_0e^{\gamma(m-m_0)}$，
長度尺度的指數便是 $\gamma/2$。
有些文獻的 $D$ 是面積，有些是長度。
Biondini 的形式使用 $D^2$ 作面積尺度。
比較前須先確認參數化。

## B.9 最大值、Båth 律與最大規模差

Båth 律描述主震與最大餘震的平均規模差。
它是序列統計规律，不是逐序列硬下界。
序列選取與餘震門檻都會影響結果。

固定 $N$ 個獨立指數超額規模，則

$$P(X_{\max}\le x)=(1-e^{-\beta x})^N.$$

由大到小排序，令第 $j$ 個間距為 $D_j$。
指數無記憶性給出 $D_j\sim\operatorname{Exp}(j\beta)$。
各排序間距互相獨立，因此

$$\begin{aligned}
\mathbb E[X_{(1)}]&=\frac1\beta\sum_{j=1}^N\frac1j,\\
\operatorname{Var}(X_{(1)})&=\frac1{\beta^2}\sum_{j=1}^N\frac1{j^2},\\
\mathbb E[X_{(1)}-X_{(2)}]&=1/\beta.
\end{aligned}$$

最後一式要求 $N\ge2$。
最大與次大之差仍有完整分布。
它不能直接當作 Båth 律的推導。
固定主震後再抽餘震，是另一個抽樣問題。
若限制餘震不超過主震，還須截斷規模分布。

## B.10 除叢後的 GR 交點

除叢（declustering）依規則分離事件群。
Gardner–Knopoff 方法使用時空窗。
被保留的子目錄會受選群規則影響。

令全目錄與子目錄擬合線分別為
$a-bm$ 與 $a_d-b_dm$。
若 $b\ne b_d$，交點及比值為

$$m_x=\frac{a-a_d}{b-b_d},\qquad
r(m)=10^{a_d-a+(b-b_d)m}.$$

$b_d<b$ 時，外推比值終將超過一。
子目錄觀測數卻不可能超過全目錄。
問題在兩條擬合線的共同外推。
它不能單獨判定危害估計的偏差方向。

## B.11 完整度的三種估法

完整度（completeness）指事件被充分記錄的範圍。
$M_c$ 是其中的最低規模。
估法須指定期間、區域與規模格距。

| 方法 | 操作 | 主要限制 |
|---|---|---|
| 最大曲率 MaxC | 找非累積規模直方圖的峰值 | 漸進漏測時可能偏低 |
| $b$ 值穩定法 | 提高門檻，尋找斜率穩定區 | 大門檻樣本少也會看似穩定 |
| KS 配適法 | 比較候選門檻以上的 GR 累積分布 | 估參數後須重新校準檢定 |

MaxC 全名為 maximum curvature。
其簡化實作以直方圖峰值定位轉折。
經驗修正量須依資料校準。
不能把固定加值當成普遍定律。

KS 指 Kolmogorov–Smirnov 方法。
其統計量為

$$D=\sup_m|\widehat F(m)-F_{\widehat b}(m)|.$$

每個候選門檻都重新估 $b$。
參數由同資料估出，標準 KS 臨界值未必適用。
模擬也應重做選門檻與估計步驟。

HORUS 的歷史完整度會變動。
Lolli（2020）的區域性摘要如下。

| 起始年 | 名目完整規模 |
|---|---|
| 1960 | 4.0 |
| 1981 | 3.0 |
| 1990 | 2.5 |
| 2003 | 2.1 |
| 2005 | 1.8 |

區域摘要不保證每格、每時段都完整。
L'Aquila 2009 序列可用來檢查短期漏測。
強震後的小事件可能被重疊波形掩蓋。
Omori 的 $c,p$ 因而會受觀測機制影響。
暖機期較高的 $M_c$ 也會減少已知歷史。
EEPAS 的長時間核尤其需要注意此點。

## 參考資料與延伸閱讀

- Lolli 等（2020），*The homogenized instrumental seismic catalog (HORUS) of Italy from 1960 to present*。刊於 SRL 91，3208–3222。先讀規模均一化與完整度。
- Biondini 等（2023），[EEPAS 在義大利的應用](https://doi.org/10.1093/gji/ggad123)。[免費機構版本](https://www.earth-prints.org/handle/2122/17084)。核對資料選取與式 2。
- Tinti 與 Gasperini（2024），[分箱規模的 b 值與區間](https://doi.org/10.1093/gji/ggae159)。[免費機構版本](https://cris.unibo.it/handle/11585/980514)。對照離散估計與修剪。
- SeismoStats 團隊，[完整度估計文件](https://seismostats.readthedocs.io/latest/user/estimate_mc.html)。免費文件。比較三種估法的資料條件。
- Utsu、Ogata 與 Matsu'ura（1995），[Omori 公式百年回顧](https://doi.org/10.4294/jpe1952.43.1)。[免費全文](https://www.jstage.jst.go.jp/article/jpe1952/43/1/43_1_1/_article)。閱讀有限窗與漏測問題。
- Ogata 與 Zhuang（2006），[時空 ETAS 的延伸](https://doi.org/10.1016/j.tecto.2005.10.016)。[免費機構全文](https://bemlar.ism.ac.jp/zhuang/pubs/ogata2006tectno.pdf)。對照空間核與產能。
- Mizrahi 等（2021），[除叢對規模分布的影響](https://doi.org/10.1785/0220200231)。[免費作者稿](https://arxiv.org/abs/2012.09053)。理解選群如何改變外推。
