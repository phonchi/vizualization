# 互動示意圖嵌入 spike 結論（2026-09-13）

- `IPython.display.HTML` 內的 inline SVG＋`<script>` 在 Jupyter Book 1.0.4 建置後可執行；`{raw} html` 也可；`_ext/teaching_render.py` 不干擾。
- 同頁多份 widget 以不同 id 前綴互不干擾；Plotly、MathJax 3 共存正常；admonition 自訂 class（`spec-card`）保留。
- 390 px：SVG `max-width:100%` 不溢出；Plotly 760 px 由 `teaching.css` 既有規則橫向捲動。
- 深色模式 selector：`html[data-theme=dark]`（切換 JS 寫 `dataset.theme`；`dataset.mode` 是偏好值，不用來寫樣式）。
- 深色下 pydata 主題會給 `.cell_output .text_html` 鋪淺色底；覆寫規則（specificity 必須 ≥ (0,6,3)）：
  `html[data-theme=dark] .bd-content div.cell_output div.output.text_html:has(.teaching-diagram){background-color:transparent;color:inherit;padding:0}`
  （`spec-card` 等 HTML 輸出同樣需要。）
- 改 `_static/teaching.css` 後增量 build 不會複製新 CSS，需清 `_build` 重建。
- Console 既有雜訊：plotly 6.9 `notebook_connected` 的 `<script type="module">import "https://cdn.plot.ly/plotly-3.7.0.min"</script>` 403（可在 `setup_plotly()` 剝掉）；`Got invalid theme mode` 為主題首次造訪雜訊。
- tooltip 只綁 mouseenter 時觸控無法看到，需要 tap／focus fallback。
- 工具：playwright 在 miniconda python；Chromium 需 `systemd-run --user --scope -p MemoryMax=3G`，不能用 `ulimit -v`。
- 樣板與最小書：scratchpad `spike_template.html`、`spike_book/`。
