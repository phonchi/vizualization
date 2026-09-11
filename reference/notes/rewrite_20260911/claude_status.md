# Claude 審查完成

使用者已明確授權指定稿件外部唯讀審查，並要求審查修正後直接推送。

- Job：`4b392d0f-a99a-479c-9ce3-f226b9f58d88`
- Child session：`7dedc168-ec46-441e-a977-b684505bd670`
- 終端狀態：`completed`，exit code 0，cleanup confirmed。
- 範圍：指定八份教學稿；沒有讀取原始PDF或帳密。
- 審查原文：[claude_review.md](claude_review.md)；執行身分：[claude_job.json](claude_job.json)。

## 核實與處理

1. 將存活曲線的「驗證」改成時間變換抽樣，明說不是獨立驗證。
2. 間隔CDF新增Poisson指數參考；附錄補固定N的精確關係。審查原文所寫指數N−1不正確，獨立推導應為N，未照抄。
3. Fano圖加入1的基準，說清楚固定總數、樣本變異數分母k−1時的期望仍為1。
4. 補明Omori核在事件間遞減，才能以當下率作上界。
5. 補明獨立疏化先作用於期望數，再換算至少一次機率。
6. 明說合成位置、規模與時間獨立，空間中央較密來自指定密度，不代表觸發。
7. 說清後驗與標準化概似可能外觀相近，但機率意義不同。
8. 更正參考線用語、補附錄連結；核實現有3006筆目錄皆在宣告時間窗，並讓兩圖明確使用同一篩選範圍。

新增核對：[claude_math_checks.json](claude_math_checks.json)。
重跑09導論、foundation_randomness、10點過程三頁；[執行log](execution_claude_revision.log)。
此後完成speak-human-tw校訂，程式與數學皆保真；沒有宣稱Claude再次審過最後文字潤稿。
