# 互動科學展改版驗收

保留第一部義大利、第二部臺灣、20章預報主線與37個閱讀頁的結構。新版由上一輪 `f7781f2` 接續，計畫與文獻範圍見 `scope.json`、`literature_coverage.md`。

## 完成範圍

- 以展覽式首頁、章節導讀、可收合目錄、閱讀進度及深淺主題重新呈現全站。
- 14組具有不同幾何、比較方式與操作的展件，含真實義大利格網／目標事件回放、Hawkes疊加、ETAS分支、預報評估時間軸與混合預報。固定合成例皆標示，不將動畫當成新模型實驗。
- Matplotlib匯入前檢查可寫快取；未處理的stderr禁止保存；空的新執行輸出不再被舊警告覆蓋。修正Markdown來源而非隱藏格式錯誤。
- 補齊Hawkes、marked point process、self-correcting／stress release、SVP、b-positive、EEPAS完整度修正、conjugate multiplier及評分分布等核心概念。明確區分回溯、擬前瞻、真正前瞻與部分前瞻，並分開time-dependent、rolling、realtime、operational等軸。
- 59份PDF建立索引；58份可抽取文字，不宣稱逐篇全文核讀。重要概念的來源與已核對段落見文獻對照表；1998 ETAS原件文字不可用及未出版資料限制已記錄。

## 已完成驗證與證據

- `execution/`：21個指定程式頁成功執行，無stderr外洩。
- `math/math_checks.json`：131項有限範圍的獨立計算檢查。
- `presentation_regressions.json`：8項格式／快取根因回歸檢查。
- `rendered_content.json`：37頁無警告、粗體殘碼或缺圖占位；離線建置無警告。
- `browser_final/summary.json`：37頁×桌機／手機×深／淺色，共148個組合通過；記錄操作、公式、圖片、版面及JS錯誤檢查。
- `lifecycle_checks.json`：14展件的鍵盤、reduced-motion、適用展件離屏暫停，及首頁資料回放、目錄、搜尋、切換主題通過。
- `svg_text_geometry.json`：84種圖形狀態文字邊界檢查通過；此幾何檢查不能代替實際截圖。
- `source_integrity.json`：28對Jupytext逐cell來源一致、AST與控制字元檢查通過；2993個內部引用、390個唯一錨點與80個本機資產確認無誤。
- `accepted_build_hashes.json`：瀏覽器驗收後的HTML仍與本機建置逐位元組相符。
- 四種裝置／主題的全頁與全展件montage、`comparison/`保留前後對照。完整原始截圖留在本專案，版控保存總覽與機器報告。

取景補充：element screenshot自動捲動會把固定導覽收入展件裁切。最終展件截圖使用文件座標裁切；保留實際網站導覽，沒有用CSS隱藏元件。`browser_final/capture_checks.json`確認未聚焦時skip link離開可視區；鍵盤與互動結果仍以完整148組驗收記錄為準。

發布的HTTP狀態與雜湊驗證另見發布後產生的 `publication.json`，不把本機建置通過等同已上線。
