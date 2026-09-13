# 互動科學展件

目前為 14 個科學展件，保留原十二個 slug，新增 `d11_hawkes_process` 與 `d03_forecast_protocols`。

主文沿用 `show_diagram("slug")`，全域載入 `exhibits.css`／`exhibits.js`。各展件直接改變圖形與資料閱讀狀態，不使用先前的固定窄圖與 radio 說明切換。

來源、14種操作、元件契約、實際義大利共用 JSON、重建命令與驗證限制，見 [sources/EXHIBITS.md](sources/EXHIBITS.md)。

每個展件都有 `standalone/<slug>.html` 自包含離線預覽。真實資料與固定合成示意各自明標。舊來源／收據留作歷史，不代表目前展件已通過瀏覽器驗收。
