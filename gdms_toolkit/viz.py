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
                bundle["text/html"] = _MATHJAX2_TAG.sub("", html)
            return bundle

        _br.NotebookRenderer.to_mimebundle = to_mimebundle
        _br.NotebookRenderer._gdms_no_mathjax = True

    pio.renderers.default = "notebook_connected"
