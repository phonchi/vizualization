# 2026-09-11 核心目錄／ETAS 教學改編核對

## 範圍與交付

編修 `book/11_catalog_completeness_b.py`、`12_clustering_laws.py`、`13_etas_structure.py`、`14_etas_estimation.py`；新增 `appendix_b_catalog.md`、`appendix_c_etas.md`。正文改為順讀式教學，移除固定研究綜述骨架；必要定義保留正文，長代數集中附錄。既有 28 個 code cells 全數保留，僅有下述兩處計算修正。未同步 notebook、未建置、未發布；交由主代理整合。

## 實讀來源

- `reference/ETAS/[1995] Omori.pdf`：Utsu、Ogata、Matsu'ura 1995，doi:10.4294/jpe1952.43.1。全文文字擷取；重點為有限資料窗的概似、早期漏測與 c/p 解讀。用於12與附錄B。
- `reference/ETAS/[1999] Review_model.pdf`：Ogata 1999，*Seismicity Analysis through Point-process Modeling: A Review*，155:471–507。實讀条件強度、概似、ETAS與殘差相關段落，供架構核對。
- `reference/ETAS/[2002] ETAS_Declustering.pdf`：Zhuang、Ogata、Vere-Jones，97:369–380，doi:10.1198/016214502760046925。實讀首頁、模型背景與權重估計，補入12、14及附錄C。
- `reference/ETAS/[2008] ETAS using EM.pdf`：Veen、Schoenberg，103:614–624，doi:10.1198/016214508000000148。實讀未知親代、incomplete data與EM型估計，補入14及附錄C；區分精確EM保證和任意背景平滑迭代。
- `reference/ETAS/[2018] Review_model.pdf`：Reinhart，33:299–318，doi:10.1214/17-STS629。實讀分支表示、邊界、概似及診斷，補入13與附錄C。
- `reference/ETAS/[2023] Bayes ETAS.pdf`：Naylor、Serafini、Lindgren、Main，doi:10.3389/fams.2023.1126759。實讀聯合後驗、率相關漏測及合成案例，補入14及附錄C。
- `reference/ETAS/[1998] ETAS.pdf`：`pdftotext` 僅取出24字並有字型警告，未宣稱直接核實其正文。時空結構改由2006來源、1999回顧與Jalilian2019支撐；未把1998細節新增為已驗證結論。
- `reference/Intro_to_the_theory_of_point_processes.pdf`：確認是Daley–Vere-Jones 2003第二版Volume I；讀取版權頁、目錄、6章群集及7章概似的相关文字，作為13/14與附錄C的概念來源。
- `reference/[2024] Estimate b.pdf`：Tinti–Gasperini，doi:10.1093/gji/ggae159。實讀式(2)–(20)及配對差分段落。發現原稿式18/19轉錄平方根錯誤，見下。
- `reference/[2006] Space–time ETAS models and an improved extension.pdf` 與 `reference/notes/etas.md`：核對核參數慣例。`reference/notes/stats.md` 用作線索，不沿用其未限定的數學結論。
- `reference/Taiwan/臺灣地區112年中大型地震震源資訊之快速彙整與提供 .pdf`：文字擷取確認子計畫輸出參數列（0.5098、0.6188、0.0031、1.1733、1.0616、0.00005、1.5934、0.6786）。14地圖保留既有外部參數實驗，明確標示固定參數、短窗、經緯度平面距離、未完整邊界校正的演示範圍。

網頁核對僅補正式來源連結：CWA官方報告、Chan–Wu2013出版社／機構書目、Helmstetter2005出版社／作者機構、van der Elst2021出版社。各正文頁與附錄延伸閱讀已補實際使用之主要來源，沒有引用 embargo 材料。

## 數學與範圍修正

1. 有限窗、c>0的Omori積分對有限p均有限；p>1只對本書無限時間正規化核必要。
2. 面積尺度 D exp(gamma Delta m) 搭配 r²，面積∝10^m 時gamma=ln10，長度指數才是ln10/2；不能把所有實測差異歸因於慣例。
3. Båth約1.2是經驗平均；1/beta是iid指數最大與次大差的期望，非硬下界；固定N極值不能直接套到規模依賴ETAS家族。
4. 除叢GR交點代表兩條線的外推不相容，不可據此直接判定完整PSHA偏差方向。較平滑累積曲線也不證明Poisson獨立性。
5. alpha<beta是無上界GR的有限平均產能條件，不是所有有限歷史條件強度的可定義條件。超臨界n>1不自動導致有限窗概似發散。
6. 完整家族1/(1-n)與平穩觸發率占比n，不能當作有限窗實現恆等式。有限規模截斷時需重算分支比；alpha=beta不保證各種門檻下n不變。
7. 平均場家族大於始祖的機率在n→1不必趨1；孤立家族即提供上界1-exp(-n)。家族模擬含計算上限，保留限制。
8. 窗前觸發源的G_i需保留下限項，不是所有事件均用1-(1+(T-ti)/c)^(1-p)。多邊形徑向積分須考慮星形條件或多段進出區域。
9. rho/phi分母是對新規模積分後的時空率；權重是模型內潛在分支機率，不是已證明物理因果。固定外部參數的背景權重圖不是完整MLE/EM。
10. 空間資訊矩陣的加權變異公式只在所聲明的可分離／已知配對簡化模型成立；不宣稱完整ETAS精確資訊矩陣。
11. R–J只在保留單一指定觸發源、忽略再觸發／背景、alpha=beta等限制下有同形率；實際dressed時間核包含卷積，不能由1/(1-n)唯一換算全部參數。
12. 移除兩份重疊台灣估計被稱獨立bootstrap、p<1幾乎必定錯設、親代機率直接代表物理來源等舊稿過度結論。
13. 相鄰帶符號差相關係數-1/2，絕對差不能沿用負相關；iid Exp(beta)下相鄰絕對差共變異為1/(3 beta²)，附錄提供條件計算。
14. Tinti–Gasperini區間項為sqrt(C/N)，原稿與程式錯寫C/sqrt(N)。幾何分布Var(X)=u(u+Delta M)=u²C，故對u作delta區間後也獨立導出sqrt(C/N)，端點條件N>C。

## Code cells 實際修改與重跑需求

- 11章原第5號（0-based，包括setup）code cell：`sN = c_val / np.sqrt(n_sub)` 改為 `sN = np.sqrt(c_val / n_sub)`。需重跑bootstrap區間圖；前置資料／函數維持原設定。
- 12章原第1號code cell：Omori標準誤不再對Hessian逆矩陣對角值取abs；要求best.success、非邊界解、Hessian正定，否則sig_p=np.nan。需重跑首張擬合圖及呼叫omori_mle的相依圖。點估計演算法不變。
- 13、14章所有code cells與HEAD相同，seed、資料、圖表計算未變。13模擬上限、14固定參數演示保留並於正文明示限制。
- 主代理需同步remove-input及ipynb、章序與顯示章號，再完成建置與實際渲染。

## 已執行最小驗證

- 四份.py以ast.parse確認語法有效；與HEAD比較，每章均保留7個code cells，僅上面11/12各一個改變。
- 原有8個公式label在四個主文中各保留一次；附錄使用引用，沒有另複製同名定義。
- scipy.quad獨立數值核對：有限Omori p=0.7、1、1.3；截斷分支alpha=1、beta、3；窗前來源時間積分；幾何區間兩側與delta映射；相鄰絕對差共變異。全部絕對誤差<1e-9，最大1.24e-14。
- 未執行整章notebook、Jupyter Book建置或瀏覽器驗收，不能據此宣稱視覺通過。
