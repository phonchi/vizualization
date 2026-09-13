# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: kernelspec,jupytext
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # N. 章名（面向學生的問句或陳述）
#
# 開場：承接前一章尚未解決的問題（一段），本章要回答什麼（一段）。
# 第一段就給一個能抓住的東西：一個數字、一欄資料或一張圖。

# %% tags=["remove-input"]
from gdms_toolkit.viz import setup_plotly
setup_plotly()

# %% tags=["remove-input"]
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML, display

from gdms_toolkit import italy
from gdms_toolkit.viz import ACCENT, PALETTE, QUAKE_COLOR, SEQUENTIAL, apply_layout, show_diagram

# %% [markdown]
# ## N.1 小節標題
#
# ```{admonition} 定義：術語（English term）
# :class: definition
# 一句話定義；符號；單位。
# ```

# %% tags=["remove-input"]
show_diagram("d02_regions_s_r", caption="圖說：一句話說明讀者要看什麼。")

# %% [markdown]
# ## N.k 本章填入的規格欄位

# %% tags=["remove-input"]
display(HTML(italy.spec_card([
    ("目錄與規模尺度", "HORUS，Mw 均一化"),
    ("深度", "≤ 40 km"),
    ("學習期／測試期", None),
])))

# %% [markdown]
# 收尾一段，前指下一章：{doc}`下一章標題 <NN_next_slug>`。
#
# ## 參考資料與延伸閱讀
#
# - 作者（年份），[標題](連結)。中文導讀一句。
