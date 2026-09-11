# 讀者通讀後的觀測與讀圖修正

範圍：03_groundwater.py、04_geomagnetic.py、21_psha.py、appendix_g_observations.md。沿用已讀 speak-human-tw 人工修順語句，未機械替換全站用語。

## 修正逐項紀錄

1. `book/03_groundwater.py` 原始第 25 行：先以水在孔隙中流動、水柱被壓力托高的圖像，介紹含水層、承壓與壓力水頭；原有條件與科學限制保留。
2. `book/03_groundwater.py` 原始第 69 行：將圖前說明對齊實際兩列圖，不再承諾圖中含水溫。
3. `book/03_groundwater.py` 原始第 84 行：圖後的縱軸說明同步限定水位、氣壓。
4. `book/03_groundwater.py` 原始第 212 行：將全欄缺值印象改為部分時段缺測；區分讀取器非缺值、儀器品質與真正可用觀測。
5. `book/03_groundwater.py` 原始第 219 行：僅改既有計數輸出的名稱，保留notna計算與全部資料數值。
6. `book/03_groundwater.py` 原始第 225 行：圖題不再讓人以為水位整欄缺值；水溫曲線與資料不變。
7. `book/04_geomagnetic.py` 原始第 138 行：把資料是否足以比較放在殘差大小之前，不預先假定TTN必定能算。
8. `book/04_geomagnetic.py` 原始第 146 行：依有限配對數與std有限性分支輸出人話；CSG仍用去平均後std，TTN不再直接印nan。
9. `book/04_geomagnetic.py` 原始第 150 行：明確說本次TTN比較不成立，保留一般殘差解讀但書，不捏造精確缺測原因。
10. `book/04_geomagnetic.py` 原始第 158 行：將ULF全名、中文意義與實際頻帶需求放在第一張相關圖前。
11. `book/04_geomagnetic.py` 原始第 178 行：補上Nyquist直覺定義與抗混疊條件，保留原一秒／一分鐘數值。
12. `book/21_psha.py` 原始第 75 行：在第一張圖前介紹PGA與g，並區分教學門檻與實測值，不更動危害公式。
13. `book/appendix_g_observations.md` 原始第 5 行：附錄開頭加入觀測篇入口，提供明確回程。
14. `book/appendix_g_observations.md` 原始第 34 行：補對應主文的{doc}回程，不猜測HTML錨點。
15. `book/appendix_g_observations.md` 原始第 57 行：補對應主文的{doc}回程，不猜測HTML錨點。
16. `book/appendix_g_observations.md` 原始第 101 行：補對應主文的{doc}回程，不猜測HTML錨點。
17. `book/appendix_g_observations.md` 原始第 125 行：補對應主文的{doc}回程，不猜測HTML錨點。
18. `book/appendix_g_observations.md` 原始第 170 行：補對應主文的{doc}回程，不猜測HTML錨點。
19. `book/appendix_g_observations.md` 原始第 187 行：補對應主文的{doc}回程，不猜測HTML錨點。

## Code 變動與重跑需求

- `03_groundwater.py`：HWA水溫圖的同一cell（從 `hwa = gt.read_groundwater(...)` 開始）僅改兩個print標籤及Plotly圖題；`notna().sum()`、10分鐘平均、資料與水溫曲線不變。需更新該cell的stream與圖表輸出。沒有加入時間範圍推測，也沒有新增分析。
- `04_geomagnetic.py`：TTN參考站比較cell（從 `mag_ttn = ...` 開始）以有限差分配對數檢查至少兩筆，再檢查去平均std是否有限。有效時照常顯示CSG/TTN數值及配對數；不足時印資料不足／無法比較，不顯示nan當作可排名數字。需要重跑此cell更新stream。
- `21_psha.py` 僅新增圖前定義，所有code cells逐字未改。
- 附錄G只有說明與已確認目標存在的{doc}連結；未新增數值公式或程式。

## 保留事項與限制

- 全部既有數學公式逐字一致。ULF／Nyquist既有一秒0.5Hz、一分鐘約0.0083Hz的數值保留。PGA的0.5g例子明確為教學門檻。
- HWA非缺值不宣稱品質正常；TTN不能比較的精確資料缺失原因未在本子任務判定。
- code cells數量與順序不變；只改上述兩個cell。
- 已做AST解析、公式保留與{doc}目標存在的靜態檢查。未同步notebook、未執行cell、未build、未push。新程式分支的實際輸出由主端重跑確認。
