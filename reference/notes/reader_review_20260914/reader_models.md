# 讀者完整順讀：模型、檢驗與推導

日期：2026-09-14（Asia/Taipei）。Git root：/home/phonchi/vizualization；閱讀時HEAD：95391eada6e9ed6fe6cc7280beb84634ad7b5b37。本報告依未提交工作目錄快照，不將HEAD視為全部閱讀內容。

## 範圍與做法

依序完整閱讀第11–20章.py，包括文字、全部公式、程式cell、圖說、規格卡、延伸閱讀；完整閱讀附錄A/C/D/E/F。再逐本解析10份ipynb的所有輸出：表格、Markdown、規格卡、Plotly內嵌JSON、座標資料、圖例、圖說與模擬樣本。沒有用關鍵字掃描代替全文閱讀。

本輪所有10份.py與ipynb的cell來源在修改前一致。未發現error輸出。未開瀏覽器，不能由本報告宣稱實際排版、觸控互動或SVG視覺驗收通過。沒有重新執行整章、重建預報、重跑大型ETAS或蒙地卡羅；小型數值核對只載入既有資料/快取與已保存模擬樣本。

## 已核實的三項問題與替換文案

### 1. 第14章把去趨勢曲線說成水平階梯（內容錯誤）

位置：book/14_psi_precursory_scale.py:328–330；對照同章:59–78、cumulative_anomaly程式及附錄D:12。

Before：
「累積規模曲線的水平段表示沒有新的入選事件。向上的階梯來自新增事件，其高度取決於採用的規模權重。」

本頁实际畫的是C(t)=A(t)-k(t-t_s)，不是原始A(t)。無事件時C'(t)=-k；向上跳躍才是新增規模超額。直接讀取ipynb cell4的322個C(t)座標，在所有非跳躍線段重算斜率，範圍為−0.088210643941866至−0.088210643941626，並非零。

After：
「沒有新的入選事件時，原始累積量A(t)保持不變，但圖中的C(t)會以斜率−k下降。新增事件使C(t)向上跳躍，跳幅是該事件的規模超額。先區分原始累積量與扣除平均趨勢後的曲線，再與時間、位置及搜尋條件對照。」

### 2. 第15章同頁保留兩套互相衝突的lead time定義（術語錯誤）

位置：book/15_eepas_italy_forecast.py:111；對照同頁:178–183、appendix_d_eepas.md:179–207。

Before：
「前置時間（lead time）則描述預報起點與可用歷史截止之間的距離。」

這其實是前述input delay；後段又將catalogue lead time定義為目錄起點到目標事件。Rhoades等2020原始PDF entropy-22-01264.pdf實體p.4明定後者；p.7式10需另處理資料截止至目標的time lag。

After：
「目錄前置時間（catalogue lead time）從目錄起點量到目標地震。它描述可用歷史向前涵蓋多久；發報時刻與資料截止之間的50天則是輸入延遲（input delay）。例如目前資料即時可用，也不會補回1960年以前未收錄的前兆貢獻。」

### 3. 附錄C空間積分與資訊矩陣例子沿用另一套D單位（公式/跨章不一致）

位置：appendix_c_etas.md:124、146；對照12_etas_structure.py:106–112、13_etas_italy_forecast.py:132及附錄C對照表:205前後。

本站D為公里，空間面積尺度為D² exp(γΔm)。但C.6寫σ=D exp(γΔm)，使R²/σ不再無量綱；C.8又以η₁=lnD作對數面積尺度。這是舊D代表面積的慣例殘留，並非兩種均可直接帶入本站數字。

Before：σ=D exp(γΔm)；η₁=lnD。
After：σ=D² exp(γΔm)，並說明σ在該式是面積尺度；η₁=2lnD、η₂=γ，說明u_i為面積尺度的對數（固定公里單位）。

獨立數值核對：設定D=2 km、R=2 km、γΔm=0、q=1.5。由主文面密度作2πr dr積分得到0.2928932188134524；附錄舊式得到0.42264973081037427。改為D²後與積分相等。

## 讀者承接觀察

11→12的Hawkes/marked ETAS已連得上，12→13也明說估計子目錄和正式快取不同。14的Ψ回溯選擇與15的逐事件三核有明確區隔；16→17拆總數/分配，18→19拆比較/權重選取，20再接地動而非直接稱建築風險，這些承接確認無誤。

第12章輸出n=1.579後，接著展示n=0.55的合成家族。兩數都正確、來源也有標示；若進一步潤稿，可增加一句「前者的無限域參數推算屬超臨界，以下另用次臨界合成設定講家族平均」，讓初學者少回頭一次。這是閱讀補強建議，不列為數值錯誤，也未據此更改模型。

附錄的證明深度高於只修過基礎統計的讀者；章首已說先備知識。沒有要求將完整EM、Borel生成函數、布朗運動變換證明搬回主文。未發現需要再增加模型章的斷裂。

## 獨立重算與確認無誤

資料使用現有HORUS實驗資料與npy快取，未調用預報重建。分數以獨立寫出的Poisson log-pmf和SciPy分布函數核對，不只比對舊報告。

- 實驗目錄435237筆；學習期目標27顆、測試期25顆；2016年6顆；12個測試窗有事件。與章內表格一致。
- 學習窗現已截止2012-01-01 00:00；測試40窗截止2021-12-31 09:36，附錄E已說明末端0.6天不評分。
- 四張測試預報總數：SUP 12.2492242022、PPE 14.0099613734、EEPAS 14.4371098671、ETAS 7.5021832615。各章地圖customdata加總與文字顯示一致。
- Poisson上尾P(N≥25)：SUP 0.000908043477、PPE 0.005062956062、EEPAS 0.007220790442、ETAS 3.769364399e−7。四者在未校正雙尾α=0.05均拒絕，與第16章一致。
- Poisson中央95%分位點區間依序[6,20]、[7,22]、[8,22]、[3,13]；外部變異67.76的NB區間依序[1,33]、[2,34]、[3,34]、[0,30]。與已保存圖線一致；沒有聲稱變異67.76由本站資料重新估得。
- 相對SUP的IGPE：PPE 0.3969582617、EEPAS 0.6986716123、ETAS 1.0543740328。逐事件項與總率修正直接相加，與第18章表格及圖點一致。
- 對SUP近似95%配對t區間：PPE [0.0389374957,0.7549790278]；EEPAS [0.2942345429,1.1031086817]；ETAS [0.2434735480,1.8652745175]。這只是既定獨立事件近似，正文已說未校正序列相依。
- 第19章51個候選以學習期分數重新計算，最大值在π_ETAS=0.76；測試混合總數9.1665656468，相對EEPAS IGPE=0.4173916996。與保存結果一致；未用測試期選權重。
- 第17章不重跑模擬。從保存的每組1000個模擬分數重新算觀測分數及低尾秩：S=(0,.432,.138,1)，M=(.240,.240,.107,.240)，cL=(.110,.174,.149,.778)，模型順序SUP/PPE/EEPAS/ETAS；與圖標/表一致。L只有表格、未保存每個模擬分數，未獨立重演L模擬。
- Hawkes核積分1.95；ETAS正規化K=0.8706872560、n=1.5793850312。原始參數和兩套K慣例的換算一致。
- 13章合成兩份200筆計數陣列之和1850、2676，平均9.25、13.38，最大20、39；直方圖0–100未切掉任何樣本。
- 14章合成C(t)最低點605.4181308天，與圖說605.4一致；起點600是生成設定，正文已區別。
- 附錄D對數常態mean/median比：σ=.20時1.1118640845；σ=.81時5.6932865084。與1.112、5.69一致。
- 附錄F的500年回歸期、50年超越機率重算為0.09516258196；Cox零事件例為0.56766764162，而固定均值Poisson為0.36787944117。
- 第20章OAF合成至少一顆機率為8.4137%、19.4783%、28.0263%、40.9953%，顯示8.4%、19.5%、28.0%、41.0%正確。
- 附錄E的AICc修正和條件形狀增益已分開命名，先前同標籤錯誤已不存在。
- 附錄D的時間補償兩端點、EAS/NW/W區別及附錄E乘法正規化的代數確認無誤。

## 限制

本輪不是重做全部研究或原論文的數值復現。未重新估13章兩個參數，只核對模型/資料集合/程式與保存結果的對應；未重建EEPAS/ETAS快取或執行完整級聯。未重新外網逐項查全部DOI，亦未驗收第20章原圖的瀏覽器顯示。

小型唯讀Python匯入套件時，matplotlib自動建立暫存字型設定目錄；沒有更改資料、book程式或快取預報。正式證據留於本報告。

## 全文閱讀清單與修改前快照

- book/11_conditional_intensity.py：全文順讀／輸出解析；SHA-256 36a0fd817d640eccb6c967f981c5fd90fbfd5b3f114bd9e6022d643355f14c4e
- book/11_conditional_intensity.ipynb：全文順讀／輸出解析；SHA-256 6113aa2f15cfb984aaeed17c95a831298f4273f55536c7aa5a12bd6fc9d173c5
- book/12_etas_structure.py：全文順讀／輸出解析；SHA-256 09013f62e436f312e23ff27b0c18de29b008c3804deb0e123f4e8defdd04c14c
- book/12_etas_structure.ipynb：全文順讀／輸出解析；SHA-256 630c1f45fe461fae1a9e7ec0a791fee94153c4939e9e994abf1fb1c73f5cfb88
- book/13_etas_italy_forecast.py：全文順讀／輸出解析；SHA-256 8680d5c1580d1cc7c00c0d36f75e126e6327b9f831961c3cc8b60243872c5df0
- book/13_etas_italy_forecast.ipynb：全文順讀／輸出解析；SHA-256 fd7422df1daae280e0dc72c33ac55d72100d6e31db5adae783bdca04d97e7ff7
- book/14_psi_precursory_scale.py：全文順讀／輸出解析；SHA-256 2568e98c09e4f586b2d4d3cab0748137758a822f41cf69c50f056df7249dcced
- book/14_psi_precursory_scale.ipynb：全文順讀／輸出解析；SHA-256 a0d44437ab33c3e74e46ea06928fb558ea10dd5c945251b240c272a9193cdfb4
- book/15_eepas_italy_forecast.py：全文順讀／輸出解析；SHA-256 c62f260fbc8aa9d183693ebe44ce30b787fe15e2bcb20dc27bf631616274109c
- book/15_eepas_italy_forecast.ipynb：全文順讀／輸出解析；SHA-256 13e03348c4e10d5d4e6939c49cfb86cd29525ab92ec3a2f8468a233a82f43c73
- book/16_test_number.py：全文順讀／輸出解析；SHA-256 7207fa6ce803f43bcd9945ac45e62f0cd047df4cf2d14dcd9f4c51b0618f072a
- book/16_test_number.ipynb：全文順讀／輸出解析；SHA-256 ce66642c92d68de73c6944a5a42887507451fd39157ec7c8e5903b58a1446f2e
- book/17_test_space_magnitude.py：全文順讀／輸出解析；SHA-256 0b5c63623f3b510948b775173bf0cb776c964b925e532921eeb2e7318f598eba
- book/17_test_space_magnitude.ipynb：全文順讀／輸出解析；SHA-256 c3ed38e35535933f50f936ec3d71975287685020823798053362ac21584690b6
- book/18_test_comparison.py：全文順讀／輸出解析；SHA-256 23001cf355b84ce81abaaf5aa4f9a7ebda08ae89ed2cc9abc9861d52869281b4
- book/18_test_comparison.ipynb：全文順讀／輸出解析；SHA-256 346b5c79ba0c503953fcdc7027b3e9a117bcdf08f5cb7804bc7b14a22fbe46f4
- book/19_ensembles.py：全文順讀／輸出解析；SHA-256 74a25d36abfff7f484d17ed909da9c6758efe4caa2c633e621a7c8a3411dae2f
- book/19_ensembles.ipynb：全文順讀／輸出解析；SHA-256 7365c94ed3c19e4e223ad26fe008e31c449cea677ce20cd97545e1201fa7cbd6
- book/20_beyond_forecast.py：全文順讀／輸出解析；SHA-256 0d00410035d34b398b3acc23373cc0d04c79a2e418608e167fcf858bfd4e7076
- book/20_beyond_forecast.ipynb：全文順讀／輸出解析；SHA-256 87ca9411891b3a0588feeabcd3ee3ecfc90310a7d95e538bf6bf9579709cd199
- book/appendix_a_point_process.md：全文順讀／輸出解析；SHA-256 7796290c64ef69d0a0396c57062366e76fff7f2bef952f222bbab30d62445300
- book/appendix_c_etas.md：全文順讀／輸出解析；SHA-256 fad3699d665422bfc367bad22adc83ae79dbcfab40a62ceb92797528ae04cc36
- book/appendix_d_eepas.md：全文順讀／輸出解析；SHA-256 59993672edee17dcbe182e4244a1f3e71d14864cbd214aff0ed6bdb99a2bbc2b
- book/appendix_e_testing.md：全文順讀／輸出解析；SHA-256 f1867a3bbba45467e988ce52f7dd197c3421c54784ec41319d3bbddfdc861555
- book/appendix_f_hazard.md：全文順讀／輸出解析；SHA-256 d558f5a052ea3a1393762fc0935001a35588578eac0a66dcf1aab42acc8cddeb

## 修正狀態

報告建立時三項修正已獲主代理授權；只修.py/.md文字與公式，不改code cell、ipynb或共用程式。下方完成後補記實際狀態。

修正完成：book/14_psi_precursory_scale.py、book/15_eepas_italy_forecast.py、book/appendix_c_etas.md。兩份.py全部code cell來源逐一對比完全不變，AST通過；未修改ipynb，由主代理統一同步。C.6與C.8一起統一D的長度慣例。
