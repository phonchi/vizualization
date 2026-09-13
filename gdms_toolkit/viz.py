"""共用視覺化設定：固定調色盤與 Plotly 版面。

配色規則：每個測網有固定顏色（色彩跟著實體走），所有章節一致。
調色盤經過色覺辨認度驗證（相鄰色 CVD ΔE ≥ 8）。
"""
import re


# 分類調色盤（依固定順序指派，不循環）
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

# 測網固定用色（全書一致）
NETWORK_COLORS = {
    "GW": "#2a78d6",        # 地下水 — 藍
    "MAGNET": "#4a3aa7",    # 地磁 — 紫
    "CWASN": "#eb6834",     # 地震觀測網 — 橘
    "TSMIP": "#e34948",     # 強震觀測網 — 紅
    "GNSS": "#1baf7a",      # GNSS（氣象署）— 青
    "GNSS_IES": "#008300",  # GNSS（中研院）— 綠
    "GNSS_ETEC": "#e87ba4", # GNSS（震工中心）— 洋紅
}

# 單一序列／量值用
SEQUENTIAL = "Blues"
ACCENT = "#2a78d6"
QUAKE_COLOR = "#e34948"   # 地震事件標記（紅，保留給事件，不作序列色）

PLOTLY_LAYOUT = dict(
    template="plotly_white",
    font=dict(size=13),
    colorway=PALETTE,
    margin=dict(l=60, r=20, t=50, b=40),
    hovermode="x unified",
)


def apply_layout(fig, **kw):
    """套用全書一致的 Plotly 版面設定。"""
    fig.update_layout(**{**PLOTLY_LAYOUT, **kw})
    return fig


_MATHJAX2_TAG = re.compile(
    r'<script src="https://cdnjs\.cloudflare\.com/[^"]*mathjax[^"]*"></script>',
    re.IGNORECASE)
# plotly 6 的 notebook_connected renderer 另外注入一支 ES module import（cdn.plot.ly 回 403），
# 圖本身由 classic script tag 載入，這支只會在 console 留下錯誤，一併剝掉。
_PLOTLY_MODULE_TAG = re.compile(
    r'<script type="module">\s*import\s+"https://cdn\.plot\.ly/[^"]*"\s*;?\s*</script>',
    re.IGNORECASE)


def setup_plotly():
    """設定各章共用的 Plotly 輸出方式。

    Plotly 的 notebook renderer 會在每張圖的 HTML 裡硬塞一支 MathJax 2，
    與 Jupyter Book 的 MathJax 3 相衝，導致整頁數學式停在原始 LaTeX。
    這裡把那支 script 從輸出中移掉（圖本身不受影響）。
    """
    import plotly.io as pio
    import plotly.io._base_renderers as _br

    if not getattr(_br.NotebookRenderer, "_gdms_no_mathjax", False):
        _orig = _br.NotebookRenderer.to_mimebundle

        def to_mimebundle(self, fig_dict):
            bundle = _orig(self, fig_dict)
            html = bundle.get("text/html")
            if html:
                html = _MATHJAX2_TAG.sub("", html)
                bundle["text/html"] = _PLOTLY_MODULE_TAG.sub("", html)
            return bundle

        _br.NotebookRenderer.to_mimebundle = to_mimebundle
        _br.NotebookRenderer._gdms_no_mathjax = True

    pio.renderers.default = "notebook_connected"


# ---------------------------------------------------------------- 預報地圖 ----
def plot_forecast_map(rate_by_cell, cells=None, targets=None, polygon=None,
                      title="", log_scale=True, zmin=None, zmax=None, height=520,
                      colorbar_title="每格期望數"):
    """全書統一的網格預報圖：177 格上色、S 邊界、目標地震紅點（只疊圖，不評分）。

    rate_by_cell：長度 177 的期望數；cells／polygon 省略時自動從 gdms_toolkit.italy 取。
    targets：含 lon、lat、mb 欄位的 DataFrame（通常是測試期目標地震）。
    """
    import numpy as np
    import plotly.graph_objects as go

    from . import italy

    cells = italy.testing_cells() if cells is None else cells
    polygon = italy.collection_polygon() if polygon is None else polygon
    z = np.asarray(rate_by_cell, float)
    features = []
    for _, c in cells.iterrows():
        ring = [[float(lo), float(la)] for lo, la in zip(c.corner_lon, c.corner_lat)]
        features.append(dict(type="Feature", id=int(c.cell),
                             geometry=dict(type="Polygon", coordinates=[ring])))
    geojson = dict(type="FeatureCollection", features=features)
    if log_scale:
        zz = np.log10(np.maximum(z, 1e-6))
        cb = dict(title=colorbar_title, tickvals=[-4, -3, -2, -1, 0],
                  ticktext=["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "1"])
    else:
        zz = z
        cb = dict(title=colorbar_title)
    fig = go.Figure()
    fig.add_trace(go.Choropleth(
        geojson=geojson, locations=cells.cell, z=zz, colorscale=SEQUENTIAL,
        zmin=zmin, zmax=zmax, marker_line_color="rgba(36,48,64,0.35)", marker_line_width=0.5,
        colorbar=cb, customdata=np.c_[z],
        hovertemplate="格 %{location}<br>期望數 %{customdata[0]:.3g}<extra></extra>",
        name="預報"))
    fig.add_trace(go.Scattergeo(lon=polygon.lon, lat=polygon.lat, mode="lines",
                                line=dict(color="rgba(36,48,64,0.7)", width=1.5),
                                name="收集區 S", hoverinfo="skip"))
    if targets is not None and len(targets):
        fig.add_trace(go.Scattergeo(
            lon=targets.lon, lat=targets.lat, mode="markers",
            marker=dict(color=QUAKE_COLOR, size=4 + 3 * (targets.mb - 5.0) * 2,
                        line=dict(color="white", width=0.8)),
            text=[f"Mw {m:.1f}" for m in targets.mb], name="目標地震（事後疊圖）",
            hovertemplate="%{text}<extra></extra>"))
    fig.update_geos(fitbounds="locations", visible=False, showcountries=True,
                    countrycolor="rgba(120,130,150,0.5)", showcoastlines=True,
                    coastlinecolor="rgba(120,130,150,0.5)", showland=True,
                    landcolor="rgba(238,242,247,0.6)", projection_type="mercator")
    apply_layout(fig, title=title, height=height, hovermode="closest",
                 margin=dict(l=10, r=10, t=50, b=10),
                 legend=dict(orientation="h", y=-0.02))
    return fig


# ---------------------------------------------------------------- 示意圖 ----
_DIAGRAM_DIR = None


def show_diagram(name: str, caption: str | None = None):
    """嵌入 book/_static/diagrams/<name>.html；缺圖即失敗，避免發布占位說明。

    互動展件隨輸出帶入本機 CSS／JS，讓 notebook 也保留靜態樣式；
    瀏覽器單例保護避免與網站共用資產重複註冊。
    """
    from pathlib import Path

    from IPython.display import HTML, display

    global _DIAGRAM_DIR
    if _DIAGRAM_DIR is None:
        _DIAGRAM_DIR = Path(__file__).resolve().parent.parent / "book" / "_static" / "diagrams"
    path = _DIAGRAM_DIR / f"{name}.html"
    if path.exists():
        html = path.read_text(encoding="utf-8")
    else:
        raise FileNotFoundError(f"缺少教學示意圖：{path}")
    if 'class="quake-exhibit"' in html:
        static = _DIAGRAM_DIR.parent
        css = (static / 'exhibits.css').read_text(encoding='utf-8')
        js = (static / 'exhibits.js').read_text(encoding='utf-8')
        html = '<style>' + css + '</style>' + html + '<script>' + js + '</script>'
    if caption:
        html += f'<p class="diagram-caption">{caption}</p>'
    display(HTML(f'<div class="teaching-diagram" data-diagram="{name}">{html}</div>'))
