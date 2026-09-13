#!/usr/bin/env python3
"""把 book/_static/diagrams/<name>.html 嵌入片段包成 diagram-design 風格的獨立 HTML。

用法：python3 build_standalone.py d02_regions_s_r:nested d04_density_to_counts:process ...
獨立檔與片段共用同一份 SVG；只在這裡加 Google Fonts、頁首、深色 token 與版權列。
（d01 的獨立檔由 archify deliver 產生，不經這裡。）
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DIAG = HERE.parent

TEMPLATE = """<!DOCTYPE html>
<html lang="zh-Hant">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <link href="https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&family=Geist:wght@400;500;600&family=Geist+Mono:wght@400;500;600&family=Noto+Sans+TC:wght@400;500;600&family=Noto+Serif+TC:wght@400&display=swap" rel="stylesheet">
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    :root {{
      --dg-paper:#fbfcfd; --dg-paper-2:#eef2f7; --dg-ink:#243040; --dg-muted:#56657a; --dg-soft:#7d8797;
      --dg-rule:rgba(36,48,64,.12); --dg-accent:#2a78d6; --dg-accent-tint:rgba(42,120,214,.08); --dg-quake:#e34948;
      --font-sans: 'Geist', 'Noto Sans TC', 'PingFang TC', 'Microsoft JhengHei', system-ui, sans-serif;
      --font-serif: 'Instrument Serif', 'Noto Serif TC', serif;
      --font-mono: 'Geist Mono', ui-monospace, monospace;
    }}
    @media (prefers-color-scheme: dark) {{
      :root:not([data-theme="light"]) {{
        --dg-paper:#1e232b; --dg-paper-2:#2a313c; --dg-ink:#e9eef4; --dg-muted:#b6c0cd; --dg-soft:#8d98a8;
        --dg-rule:rgba(233,238,244,.12); --dg-accent:#6aa5ec; --dg-accent-tint:rgba(106,165,236,.12); --dg-quake:#f07a79;
      }}
    }}
    :root[data-theme="dark"] {{
      --dg-paper:#1e232b; --dg-paper-2:#2a313c; --dg-ink:#e9eef4; --dg-muted:#b6c0cd; --dg-soft:#8d98a8;
      --dg-rule:rgba(233,238,244,.12); --dg-accent:#6aa5ec; --dg-accent-tint:rgba(106,165,236,.12); --dg-quake:#f07a79;
    }}
    body {{ font-family: var(--font-sans); background: var(--dg-paper); color: var(--dg-ink);
      min-height: 100vh; display: flex; align-items: center; justify-content: center; padding: 3rem 1rem; }}
    .frame {{ max-width: 960px; width: 100%; }}
    .eyebrow {{ font-family: var(--font-mono); font-size: 0.66rem; font-weight: 500; letter-spacing: 0.18em;
      text-transform: uppercase; color: var(--dg-muted); margin-bottom: 0.5rem; }}
    h1 {{ font-family: var(--font-serif); font-size: clamp(1.5rem, 2.4vw + 0.75rem, 2rem); font-weight: 400;
      letter-spacing: -0.02em; line-height: 1.15; color: var(--dg-ink); margin-bottom: 1.5rem; }}
    .teaching-diagram {{ max-width: 100%; }}
    footer {{ margin-top: 1.5rem; padding-top: 0.75rem; border-top: 1px solid var(--dg-rule);
      font-family: var(--font-mono); font-size: 0.66rem; color: var(--dg-soft); letter-spacing: 0.06em; }}
  </style>
</head>
<body>
  <div class="frame">
    <p class="eyebrow">{eyebrow} · Diagram Design · quake-teaching</p>
    <h1>{title}</h1>
    <div class="teaching-diagram" data-diagram="{name}">
{fragment}
    </div>
    <footer>{name}.html · 與 book/_static/diagrams/{name}.html 共用同一份 SVG；顏色全走 --dg-* 變數</footer>
  </div>
</body>
</html>
"""


def main(specs: list[str]) -> None:
    for spec in specs:
        name, _, kind = spec.partition(":")
        frag = (DIAG / f"{name}.html").read_text(encoding="utf-8")
        m = re.search(r"<title id=\"[^\"]+\">([^<]+)</title>", frag)
        title = m.group(1) if m else name
        out = HERE / f"{name}.html"
        out.write_text(TEMPLATE.format(title=title, eyebrow=(kind or "diagram").capitalize(), name=name, fragment=frag), encoding="utf-8")
        print("wrote", out, f"{out.stat().st_size/1024:.1f} KB")


if __name__ == "__main__":
    main(sys.argv[1:])
