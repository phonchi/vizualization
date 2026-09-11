# 連貫教學改編：交付與驗收

完成日期：2026-09-11。Git root：`/home/phonchi/vizualization`。
基底 HEAD：`f970df6f708da5b0623a56711213899c1e1024ed`。此為提交前驗收；使用者已授權審查修正完成後直接推送與更新網站。

## 交付內容

首頁、26章主文與A–G七組附錄，共34個閱讀頁面。
第一部17章由目錄、隨機性與估計鋪陳地震預報；第二部9章把模型判斷帶回台灣觀測；
第三部集中技術細節。既有頁面檔名與網址保留，章號依新閱讀順序排列。

新增 `foundation_catalog`、`foundation_randomness`、`foundation_inference`。
各章全面改寫主文、圖前問題與圖後解釋，取消重複的研究綜述骨架；短補充收合、
長推導移附錄。網站不顯示程式輸入，配對notebook保留程式與必要圖表輸出。
操作性帳號／下載回覆也從可下載notebook的輸出移除。

本機成果：`book/_build/html/00_intro.html`；來源入口：`book/00_intro.md`。
章序對照見 [chapter_map.json](chapter_map.json)。

## 已完成的驗證

| 項目 | 實際結果與證據 |
|---|---|
| 全站建置 | 確認無誤；34頁建置成功，`-W --keep-going` 無警告；[最終建置log](build_after_humanize.log) |
| Notebook配對 | 確認無誤；25組 `.py`／`.ipynb` 內容與標籤一致，共155個保存輸出，沒有error輸出；[配對結果](pair_validation.json) |
| 程式執行 | 確認無誤；16個新增或程式有變動的頁面離線執行成功；其餘僅承接相同程式來源的輸出；[執行摘要](execution_summary.json) |
| 數學獨立重算 | 確認無誤；13項針對有限窗積分、正規化、估計、預測變異與共用資料的檢查通過；[數值結果](math_checks.json) |
| 來源交叉引用 | 確認無誤；110個章節連結、19個公式引用，無遺失或重複公式label；[結構檢查](structure_checks.json) |
| HTML站內連結 | 確認無誤；2,329個本地檔案／錨點連結無缺漏；[結果](html_links.json) |
| 瀏覽器 | 確認無誤；Chromium桌機1440×900、手機390×844，全34頁共68個情境；無頁面例外、可見程式輸入、公式渲染錯誤或整頁橫向溢出；[結果](browser/results.json) |
| 收合與手機圖表 | 確認無誤；教學收合開關正常，9個代表圖表實測局部水平捲動；地圖版權標示獨立保留，未當成教學收合 |
| 原始碼差異 | `git diff --check` 通過；來源hash見 [manifest](source_manifest.json)，工作目錄見 [Git快照](git_status.txt) |

完整過程（包含已修復的失敗嘗試）存於 [run.log](run.log)。
51張代表截圖與 [全項目montage](browser/montage.png) 保存於 `browser/`。
實際執行使用本地資料快取並停用外部資料請求；觀測tgz解析較耗時，沒有把它宣稱為30秒內的教學模擬。

## 主要科學修正

- 區分沿無新事件路徑的下一事件存活機率，與會經過較小中介事件的目標規模預報。
- 補償子一般依賴隨機歷史；超臨界不等於有限觀察窗概似發散；Omori的 $p>1$ 限制指無限時間正規化。
- 修正離散 $b$ 值區間中的平方根、空間面積／長度尺度、Båth期望與硬下界混淆，以及無效Hessian標準誤。
- 修正Ψ累積規模權重的符號、起點左極限的事件分組，並取消時間完整度必然隨lag單調下降的宣稱。
- 修正Poisson日計數漏掉零事件日、截斷直方圖的正規化、KS統計量、PSHA有限規模門檻的率換算。
- 區分數位counts與校正後波形、PSD峰值與振幅、氣壓斜率與無量綱效率，以及統計觸發權重與物理因果。

這些是有針對性的重算與來源核對，並非宣稱所有現存數值程式已通過完整研究級驗證。
文獻對照：[核心模型](../rewrite_core_sources.md)、[預報與檢驗](../rewrite_forecast_sources.md)、[台灣觀測](../rewrite_taiwan_sources.md)。
Daley–Vere-Jones教材核實為2003年第二版Volume I；新增ETAS七篇均納入文獻盤點。
1998掃描版由早期唯讀文獻代理以記憶體OCR核對，核心編修代理另以相關原文交叉確認，來源紀錄保留各自查核限制。
大埔論文表格的4／7天標示矛盾已在相關教學頁與附錄說明；未自行改造原文機率數值。

## 呈現修正與範圍限制

Sphinx將相同inline script分別以 `None` 與空字串註冊，造成Thebe重複宣告。
新增專案內extension僅去除相同inline資產，不改套件環境；[根因與修正](render_fix.md)。
新序列圖補上刻線stroke。窄螢幕保留圖表可讀寬度並在圖框內捲動，避免把科學座標與圖例壓縮。

使用者明確授權後，Claude審查已完成（exit code 0）；建議逐項核實並修正，審查內固定N間隔公式的指數錯誤亦已更正。詳見 [Claude核實紀錄](claude_status.md)。

最後使用 speak-human-tw 完成215處校訂，全34頁保真比對通過；詳見 [校訂摘要](humanize_summary.md)。

未修改使用者的 `EEPAS/`、`psi-modeling-explained/`、既有私人講義或PNG。
未啟用付費服務；PDF原文與未發表材料未加入公開網站。
