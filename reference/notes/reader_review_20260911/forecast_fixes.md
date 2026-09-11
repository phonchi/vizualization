# 讀者回饋修正：預報後半部

日期：2026-09-11。範圍僅原 15–18 四章，沿用已读 speak-human-tw 修飾新增教學文字。沒有更改其他章、同步 notebook、執行分析、建置或發布。

## 修正內容

1. Ψ 第一圖之前定義 M^- 為先前期最大三個規模的平均，r 為兩期事件率之比；說明 M.U. 是規模單位、累積量不是地震能量。保留原图與數值。
2. Ψ 多重辨識圖之前分開七個合成點與愛琴海的真實兩點，說明不合併估斜率；直接連結已核實的出版來源。
3. EEPAS 三事件時間空間疊加圖及圖說從 11.7 移至 11.3 的正規化例子之後、11.4 之前。此處 B_GR 已定義，且 imports、LN10、DAY_YR 已可用。原 code cell 逐字保留，只有位置改變；11.7 留未知中間位置的卷積直覺與附錄 D 連結。
4. EEPAS 日本成績圖之前說明 I 為每事件平均對數分數，取指數的差以 1 為持平；大於一不代表各窗機率同倍放大，連往預報比較章。
5. 二元圖之前將 POLL/BILL 與兩種單格對數分數對上；縱軸是兩種計分方式之差，不是模型間資訊增益。
6. 決策圖之前說明 V 的零、一及負值基準；可靠度圖與決策圖兩處『氣候基準』改為『長期基準』，沒有改數值。

## 愛琴海兩點的直接來源

- 已讀 `reference/notes/psi.md`，該筆記沒有這兩點的直接數值，不能用筆記替代原文確認。
- 直接讀 `reference/[2022] A 20-Year Journey of Forecasting with the EETAS model.pdf`。首頁確認 Rhoades、Rastin、Christophersen，2022，Geosciences 12, 349，DOI 10.3390/geosciences12090349。
- PDF 第 2 頁正文／第 3 頁 Figure 1：2020-10-30 愛琴海 M6.7，T_P=10,220 天、A_P=3,203 km²。
- PDF 第 19 頁 §6.3／Figure 12：同事件另一辨識，T_P=6,392 天、A_P=8,091 km²；正文明說可對照 Figure 1。
- 同篇已在原章延伸閱讀中；本次只補圖前直接歸屬與 DOI。未將兩點誤配給 2024 自動辨識論文。

## Code cell 與輸出

- 15、17：只改 markdown，全部 code cells 逐字不變。
- 16：只移動一個既有 code cell，全部 code cells 的內容集合逐字不變。依賴順序靜態檢查通過；可按原始 code hash 重用輸出，不需為圖重算。
- 18：兩個 code cells 各只改圖中 annotation 的『氣候』→『長期』。其他公式與數值不變。若不更新既有 Plotly JSON 的對應 annotation，需重跑該頁以更新這兩圖。
- 四檔 AST parse、EEPAS equation label 唯一性、B_GR/LN10/DAY_YR 在移動圖之前可用：確認無誤。沒有执行 notebook 或build，實際渲染由主代理統一驗證。
