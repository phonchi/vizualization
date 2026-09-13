# 本機驗收完成

- 第一部義大利20章、第二部台灣9章、附錄7篇及首頁，共37頁。
- 28對Jupytext來源確認一致；20章及台灣整合章共21頁已執行。最慢cell 1.507秒；文字修改保留相同程式的輸出。
- 81項獨立數值／邊界檢查通過；2,650個內部連結含錨點通過；離線建置零警告。
- 148個桌機／手機、明／暗頁面組合通過；每頁HTML雜湊已對回現況。
- 12張互動圖，24組觸控／鍵盤／減少動態偏好與12張桌機hover測試通過；實際主題按鈕通過。
- 已檢視手機與桌機總覽、明暗示意圖與代表性圖表；archify兩張獨立流程圖的四尺寸／兩主題自動證據與實際圖像檢閱分開保存。

主要證據：`math_checks.json`、`structure_checks.json`、`internal_links.json`、`browser_acceptance.json`、`browser/results.json`、`interactions/`、`reader_review.md`及八張`*-all-pages.png`／`*-all-diagrams.png`。
所有逐頁完整截圖仍在專案內`browser/`；版控保存總覽與機器可讀結果，避免重複存放數百張中間截圖。

發布後的線上位元組核對已完成，41／41檔一致，詳見`publication.json`及`COMPLETION.md`。
