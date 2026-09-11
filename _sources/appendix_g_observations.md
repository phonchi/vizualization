# 附錄 G：觀測資料、處理與幾何細節

本附錄承接臺灣觀測篇。主文保留觀測機制與讀圖所需的條件，這裡集中
說明資料取得介面、取樣、單位與座標幾何。程式儲存在教材原始檔，閱讀
下列說明不需要先登入資料服務或執行下載。

## G.1 資料產品與可重現的取得流程

GDMS 提供不同類型的地震與地球物理資料。重現本書圖表時，應先使用
對應的教材快取，並核對測站、日期、取樣與規模門檻；改下載另一版本的
資料，不一定得到完全相同的事件數或時間序列。

本專案的 `gdms_toolkit` 將登入、申請、查詢打包結果與下載分開。
它是教材輔助工具，不應被誤認為官方保證穩定的公開 API。介面行為以
專案程式及官方產品說明為準，服務的認證或欄位調整可能需要更新工具。

| 介面 | 用途 | 對應原始檔 |
|---|---|---|
| `GDMSSession` | 建立已登入的工作階段 | `gdms_toolkit/auth.py` |
| `request_geophysical` | 申請地球物理資料 | `gdms_toolkit/download.py` |
| `request_waveform`、`list_channels` | 波形申請與頻道查詢 | `gdms_toolkit/download.py` |
| `list_my_downloads`、`wait_and_fetch` | 查詢打包與取得結果 | `gdms_toolkit/download.py` |
| `read_groundwater`、`read_geomagnetic`、`read_waveform` | 解析觀測產品 | `gdms_toolkit/readers.py` |

專案需要 Python 3.10 以上；本機安裝依 repository 的環境說明。登入設定
由本機 `.env` 讀入；帳密與 cookie 都是個人憑證，不能寫進公開教材或
輸出紀錄。遇到驗證碼應使用官方互動登入流程，登入失敗與沒有觀測資料
是不同問題。這裡不假定工作階段一定在某個固定分鐘數到期。

申請後尚未完成打包、帳號看不到某資料，以及該站該期間沒有釋出資料，
需要分別處理。保留原始申請條件與狀態即可協助追查，不能只由「無資料」
回覆就斷定感測器壞了。重複分析應重用已取得的檔案，避免反覆申請。

## G.2 缺值、時間與降取樣

本書使用的地下水與地磁產品有各自的缺值旗標；讀取器把指定旗標轉為
`NaN`，再依時間範圍取資料。遇到不同版本的產品，仍應查核旗標、單位
與欄位順序。保留原始檔與可用率，才能分辨是資料處理改變了結果，還是
觀測本身不同。

UTC 與臺灣當地時間相差八小時。GNSS 原始觀測及解算產品另可能使用
GNSS 時間系統，須按檔頭確認。若把兩種時間系統直接合併，即使曲線形狀
都合理，也可能將到時或短暫反應錯開。

設取樣間隔為 $\Delta t$，取樣頻率為 $f_s=1/\Delta t$。Nyquist
頻率為 $f_s/2$，所以一秒資料的上限為 0.5 Hz，一分鐘資料約為
0.0083 Hz。這只是取樣界線，實際可用頻帶還受儀器、抗混疊與噪訊影響。

一分鐘平均會先把一段時間混合成一個值，與單純刪去多餘點不同。理想矩形
平均窗的頻率響應具有 sinc 形狀，會衰減部分頻率；再降低取樣率時，
高頻成分若未被足夠抑制，仍可能混疊。為畫長時間趨勢而平均的資料，
不能再被當成完整保留短暫峰值的原始資料。

單純每隔幾點取一點的顯示抽樣，也不等同經過適當抗混疊的研究用降取樣。
花蓮案例的波形瀏覽圖使用前者，精確峰值與到時應回到原始波形處理。

## G.3 地下水的水頭、頻譜與氣壓斜率

在靜水近似下，水頭可寫為

$$h=z+\frac{P_w}{\rho_w g},$$

其中 $z$ 為高度，$P_w$ 為相應壓力，$\rho_w$ 為水密度。
井水位能否代表含水層壓力水頭，取決於水力連通、井內儲水與排水條件。
這個式子不足以把所有觀測井都視為直接孔隙壓力計，也沒有給出水位到地殼
應變的通用換算。

### 頻譜峰值不是水位振幅

SciPy 的 Welch 估計預設 `scaling="density"`。水位單位若是 cm，
頻率單位是 Hz，功率譜密度的單位就是 cm²/Hz。因此對某頻帶內的
最大 PSD 開根號，得到 cm/√Hz，並不是 cm 振幅。這正是地下水章
既有 `m2_amplitude` 輔助函式所計算的量；函式名稱沿用，正文稱為
半日頻帶峰值指標。

若 $S_h(f)$ 是單邊功率譜密度，頻帶變異量可由

$$V_{[f_1,f_2]}\approx\int_{f_1}^{f_2} S_h(f)\,\mathrm{d}f$$

估計。只有在適當分離背景、頻帶包含完整峰值且訊號近似單一正弦等
條件下，才可將 $\sqrt{2V_{[f_1,f_2]}}$ 對應到該正弦振幅。
正式潮汐分析還需分開接近的頻率及其相位，不能用一個峰值代替。

### 氣壓響應斜率與氣壓效率

主文示範對水位與氣壓取差分後作簡單迴歸，輸出斜率 $s$ 的單位為
cm/hPa。若採理想化線性反應的約定

$$\mathrm{d}h=-\frac{BE}{\rho_w g}\,\mathrm{d}P_a,$$

則 $BE=-\rho_w g\,\mathrm{d}h/\mathrm{d}P_a$ 是無量綱量。
由 cm/hPa 換到 m/Pa，斜率需乘 $10^{-4}$；再乘水的重量密度與
負號，才符合這個定義。因而不能把以 cm/hPa 輸出的數字直接叫作
無量綱效率。

差分不會保證消除潮汐、抽水與氣壓的延遲反應。原始示範先刪除缺值再
差分，缺口兩端也可能成為相鄰有效紀錄；正式分析應另檢查時間間隔，
避免把不等長差分混作同一種變化率。這個示範用來說明迴歸概念，
不構成完整水文地質參數辨識。

## G.4 波形頻道與物理校正

SEED 頻道代碼通常包含頻帶、來源及方向。下表協助閱讀本書用到的
代碼；完整範圍仍以 FDSN 規範與個別站 metadata 為準。

| 代碼例子 | 可以先辨認的資訊 | 不能只由代碼保證的資訊 |
|---|---|---|
| `HHZ` | H 頻帶、高增益地震儀、垂直分量 | 是否飽和、實際靈敏度 |
| `HH1`、`HH2` | 編號水平分量 | 是否對準正北、正東 |
| `HN?` | 加速度感測器來源 | 大震時一定不超出動態範圍 |

原始 counts 反映數位化後的儀器輸出。去均值、去趨勢或帶通濾波，不會
自動完成儀器響應移除。比較不同儀器的實際速度、加速度或地動峰值前，
需要正確的響應資訊、單位以及適合的頻帶處理。

P、S 到時差與震源距離的簡單關係，也只能在已指定傳播速度等近似下
使用。在均勻介質、同一路徑距離 $r$ 的近似下，

$$\Delta t_{SP}=r\left(\frac1{v_S}-\frac1{v_P}\right),
\qquad r=\frac{v_Pv_S}{v_P-v_S}\Delta t_{SP}.$$

這裡的 $r$ 是傳播路徑的近似距離，不能不計深度就稱為震央距離。
多站定位還需速度模型與到時誤差，正式做法不由單一固定乘數取代。

## G.5 RINEX、座標轉換與地表位移

RINEX 依版本定義標頭、觀測量與曆元格式。以教材中的 RINEX 2 檔名
`hual0930.24o.gz` 為例，它包含站碼、年積日與年份等資訊；年積日必須
按該年是否閏年解讀。標頭中的 `APPROX POSITION XYZ` 是地心直角
近似座標，不是整天各曆元的精密位置。

### 從地心直角座標到大地座標

設橢球長半軸為 $a$、第一偏心率平方為 $e^2$，大地緯度為 $\varphi$、
經度為 $\ell$、橢球高為 $h$，令

$$\nu(\varphi)=\frac{a}{\sqrt{1-e^2\sin^2\varphi}}.$$

正向關係為

$$\begin{aligned}
X&=(\nu+h)\cos\varphi\cos\ell,\\
Y&=(\nu+h)\cos\varphi\sin\ell,\\
Z&=[\nu(1-e^2)+h]\sin\varphi.
\end{aligned}$$

反算時，經度由 $\operatorname{atan2}(Y,X)$ 得到；緯度與橢球高可
迭代求解。教材保留的短函式適合本例測站位置查核，不宣稱處理極點等
所有數值邊界。橢球高也不等於相對平均海面的正高，兩者需要大地水準面
資訊才能互換。

### 從座標差到東、北、上

若兩個時刻的精密座標已在同一框架中，先取差
$\Delta\mathbf x=(\Delta X,\Delta Y,\Delta Z)^T$，再在參考位置
以旋轉矩陣換成局部座標下的位移：

$$\begin{pmatrix}\Delta E \\ \Delta N \\ \Delta U\end{pmatrix}
=\begin{pmatrix}
-\sin\ell&\cos\ell&0\\
-\sin\varphi\cos\ell&-\sin\varphi\sin\ell&\cos\varphi\\
\cos\varphi\cos\ell&\cos\varphi\sin\ell&\sin\varphi
\end{pmatrix}\Delta\mathbf x.$$

這是幾何表示的轉換，並沒有消除座標誤差。對固定旋轉矩陣 $R$，
共變異數的線性傳播為 $\Sigma_{ENU}=R\Sigma_{XYZ}R^T$。
日解、較高頻解與長期速度產品有不同時間解析度，讀圖時必須使用
對應的不確定性，不能用每日階變宣稱解析了秒級運動。

## G.6 花蓮與大埔圖表的解讀範圍

花蓮水位示範的波動指標是相鄰一秒讀值差分的標準差，震前使用十分鐘、
震後使用五分鐘。它不是對所有同震反應的充分統計量，也沒有在本書
另行校準最小可檢測幅度。若要報告檢測能力，應事先指定訊號形狀及
假警報控制，再以合適背景時窗或模擬評估。

臺灣展望的大埔機率圖，來源為 Hsieh et al.（2025）§3.2 與 Table 2。
原表第三列印為「4 Days」，但該節正文的觀測敘述與結論採七天，並列
一、三、七、十天的模擬視窗。教材依正文保留七天標示與原有機率數值，
明確記錄這個原文不一致；未自行重新模擬或改寫論文數值。

這張圖中的不同時窗共享發報時刻，且彼此包含；不是四個獨立實驗。
機率和實際次數也不共用量尺，不能由雙軸曲線的相對高度判斷預報
偏高或偏低。計數評估需完整計數分布，二元評估則需將觀測轉成
是否至少一次，並保留跨窗、跨事件的相依性。

## 參考資料與延伸閱讀

- 中央氣象署，〈[臺灣地震與地球物理資料管理系統](https://gdms.cwa.gov.tw/)〉，官方資料入口。查核產品、時間範圍與申請說明，與教材快取分開看；程式介面是本專案的輔助工具。
- FDSN，〈[Channel codes](https://docs.fdsn.org/projects/source-identifiers/en/v1.0/channel-codes.html)〉，免費規範。核對頻帶、來源與方向代碼，並搭配測站響應資料讀取波形。
- IGS／RTCM RINEX Working Group，〈[RINEX](https://igs.org/wg/rinex/)〉，免費格式文件。按原始檔版本選規格，特別區分標頭近似座標、時間系統與逐曆元觀測量。
- Sanz Subirana、Juan Zornoza、Hernández-Pajares（2011），〈[Transformations between ECEF and ENU coordinates](https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates)〉，ESA Navipedia，免費教材。對照局部座標旋轉矩陣與橢球緯度的幾何意義；轉換不會提升定位精度。
- SciPy 開發團隊，〈[welch](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html)〉，免費官方文件。重點檢視 `density` 與 `spectrum` 的單位差異，理解為何 PSD 峰值開根號不直接是水位振幅。
- ObsPy 開發團隊，〈[Tutorial](https://docs.obspy.org/tutorial/)〉，免費教材。從波形讀取、濾波與響應移除的區別，回查本書原始檔所做的處理。
- Roeloffs, E. A.（1988），〈[Hydrologic precursors to earthquakes: A review](https://doi.org/10.1007/BF00878996)〉，*Pure and Applied Geophysics*，全文需訂閱。延伸水位、含水層反應與背景因素的關係；可先讀 USGS 的[免費地下水效應說明](https://www.usgs.gov/faqs/how-does-earthquake-affect-groundwater-levels-and-water-quality-wells)。
- Hsieh、Chang、Tai、Chen、Lu（2025），〈[Fast report: performance of the ETAS model in forecasting aftershock occurrence and site-specific ground-shaking intensity for the 2025 Dapu, Taiwan, earthquake sequence](https://doi.org/10.1007/s44195-025-00097-7)〉，*Terrestrial, Atmospheric and Oceanic Sciences*，免費全文。原始預報與觀測數值出自 §3.2、Table 2，閱讀時注意本附錄指出的視窗標示不一致。
