#!/usr/bin/env python3
"""驗收 book/_static/diagrams/<name>.html 嵌入片段。

用法：ulimit -v 3000000; python3 check_fragment.py d02_regions_s_r [more names]

每個片段：包進模擬 show_diagram() 的最小 HTML 殼（<div class="teaching-diagram">），
以 chromium headless 開啟 1440×900 與 390×844，檢查
  - 無水平溢出（scrollWidth <= innerWidth）
  - console 無 error／pageerror
  - 每個 .dg-hot[data-for] 元素 hover 與 focus 時，對應說明元素可見；移開後隱藏
  - 每個 radio 步驟控制：點選第 k 步後，[data-step<=k] 可見、[data-step>k] 隱藏
  - 深色變數覆寫（html[data-theme=dark] 層級）後截圖，並確認 svg 內沒有寫死的白／黑色
截圖與 <name>.result.json 寫在本目錄。一次只開一個 chromium。
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
DIAG = HERE.parent.parent  # book/_static/diagrams

DARK = {
    "--dg-paper": "#1e232b", "--dg-paper-2": "#2a313c", "--dg-ink": "#e9eef4",
    "--dg-muted": "#b6c0cd", "--dg-soft": "#8d98a8", "--dg-rule": "rgba(233,238,244,.12)",
    "--dg-accent": "#6aa5ec", "--dg-accent-tint": "rgba(106,165,236,.12)", "--dg-quake": "#f07a79",
}

SHELL = """<!doctype html><html lang="zh-Hant"{theme_attr}><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
html,body{{margin:0}}
body{{padding:16px;background:var(--dg-paper,#fbfcfd);color:var(--dg-ink,#243040);
 font-family:system-ui,"Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif}}
.teaching-diagram{{max-width:100%}}
{dark_css}
</style></head><body>
<div class="teaching-diagram" data-diagram="{name}">{fragment}</div>
</body></html>"""


def build_shell(name: str, fragment: str, dark: bool) -> str:
    dark_css = ""
    theme_attr = ""
    if dark:
        theme_attr = ' data-theme="dark"'
        decl = ";".join(f"{k}:{v}" for k, v in DARK.items())
        dark_css = f'html[data-theme="dark"]{{{decl}}}'
    return SHELL.format(name=name, fragment=fragment, dark_css=dark_css, theme_attr=theme_attr)


JS_VISIBLE = """(el) => { if (!el) return false; const cs = getComputedStyle(el);
  if (cs.visibility === 'hidden' || cs.display === 'none') return false;
  if (parseFloat(cs.opacity) === 0) return false; return true; }"""


def check_one(browser, name: str) -> dict:
    frag_path = DIAG / f"{name}.html"
    fragment = frag_path.read_text(encoding="utf-8")
    result: dict = {"name": name, "fragment_kb": round(frag_path.stat().st_size / 1024, 1),
                    "viewports": {}, "hover": [], "steps": [], "errors": [], "static": {}}

    # 靜態檢查：外部資源、未加前綴的 id、寫死顏色
    result["static"]["external_refs"] = re.findall(r'(?:src|href)=["\'](https?:)?//[^"\']+', fragment)
    ids = re.findall(r'\bid="([^"]+)"', fragment)
    result["static"]["unprefixed_ids"] = [i for i in ids if not i.startswith(f"dg-{name.split('_')[0]}")]
    result["static"]["hardcoded_hex"] = sorted(set(re.findall(r'#(?:fff|ffffff|000|000000)\b', fragment, re.I)))
    result["static"]["has_viewbox"] = 'viewBox="' in fragment
    result["static"]["has_html_tags"] = bool(re.search(r"<(html|head|body)\b", fragment, re.I))
    result["static"]["has_script"] = "<script" in fragment
    sizes = [float(s) for s in re.findall(r'font-size:\s*([0-9.]+)px', fragment)]
    result["static"]["min_font_px"] = min(sizes) if sizes else None

    for dark in (False, True):
        shell = build_shell(name, fragment, dark)
        shell_path = HERE / f"_shell_{name}{'_dark' if dark else ''}.html"
        shell_path.write_text(shell, encoding="utf-8")
        vps = [(1440, 900), (390, 844)] if not dark else [(1440, 900)]
        for (w, h) in vps:
            ctx = browser.new_context(viewport={"width": w, "height": h}, device_scale_factor=1)
            page = ctx.new_page()
            errs: list[str] = []
            page.on("console", lambda m: errs.append(f"console.{m.type}: {m.text}") if m.type == "error" else None)
            page.on("pageerror", lambda e: errs.append(f"pageerror: {e}"))
            page.goto(shell_path.as_uri())
            page.wait_for_timeout(300)
            page.evaluate("document.fonts && document.fonts.ready")
            key = f"{'dark_' if dark else ''}{w}x{h}"
            metrics = page.evaluate("""() => ({
                scrollWidth: document.documentElement.scrollWidth, innerWidth: window.innerWidth,
                svgWidth: (document.querySelector('.teaching-diagram svg')||{}).getBoundingClientRect ?
                    document.querySelector('.teaching-diagram svg').getBoundingClientRect().width : null,
                scrollHeight: document.documentElement.scrollHeight})""")
            metrics["no_h_overflow"] = metrics["scrollWidth"] <= metrics["innerWidth"]
            if not dark and w == 1440:
                # 每個 <text> 的 bbox 必須落在 viewBox 內（右緣留 8px）；隱藏中的說明也一併量
                result["static"]["text_overflow"] = page.evaluate("""() => {
                  const svg = document.querySelector('.teaching-diagram svg'); if (!svg) return ['no svg'];
                  const vb = svg.viewBox.baseVal; const out = [];
                  svg.querySelectorAll('text').forEach(t => { const b = t.getBBox();
                    if (b.x < 0 || b.x + b.width > vb.width - 8 || b.y + b.height > vb.height)
                      out.push(t.textContent.trim().slice(0, 30) + ' → right=' + Math.round(b.x + b.width)); });
                  return out; }""")
            shot = HERE / f"{name}_{key}.png"
            page.screenshot(path=str(shot), full_page=True)
            metrics["screenshot"] = shot.name
            if dark:
                # 檢查深色下 svg 文字實際顏色不是深色 ink（反色是否生效）
                metrics["sample_text_fill"] = page.evaluate(
                    "() => { const t=document.querySelector('.teaching-diagram svg text'); return t?getComputedStyle(t).fill:null }")
                metrics["svg_bg_fill"] = page.evaluate(
                    "() => { const r=document.querySelector('.teaching-diagram svg rect'); return r?getComputedStyle(r).fill:null }")
            result["viewports"][key] = metrics

            if not dark and w == 1440:
                # hover / focus 檢查
                hots = page.query_selector_all(".teaching-diagram .dg-hot[data-for]")
                for i, hot in enumerate(hots):
                    tip_id = hot.get_attribute("data-for")
                    tip = page.query_selector(f"#{tip_id}")
                    entry = {"hot": hot.get_attribute("id") or f"hot{i}", "tip": tip_id}
                    entry["hidden_before"] = not page.evaluate(JS_VISIBLE, tip)
                    target = hot.query_selector("[data-hover-target]") or hot.query_selector("text") or hot
                    target.hover(timeout=5000)
                    page.wait_for_timeout(300)
                    entry["visible_on_hover"] = page.evaluate(JS_VISIBLE, tip)
                    if i == 0:
                        page.screenshot(path=str(HERE / f"{name}_hover.png"), full_page=True)
                    page.mouse.move(0, 0)
                    page.wait_for_timeout(300)
                    hot.focus()
                    page.wait_for_timeout(300)
                    entry["visible_on_focus"] = page.evaluate(JS_VISIBLE, tip)
                    page.evaluate("document.activeElement && document.activeElement.blur()")
                    page.mouse.move(0, 0)
                    page.wait_for_timeout(300)
                    entry["hidden_after"] = not page.evaluate(JS_VISIBLE, tip)
                    entry["ok"] = all([entry["hidden_before"], entry["visible_on_hover"],
                                       entry["visible_on_focus"], entry["hidden_after"]])
                    result["hover"].append(entry)
                # 逐步揭露檢查
                radios = page.query_selector_all(".teaching-diagram input[type=radio]")
                for r in radios:
                    val = r.get_attribute("value")
                    rid = r.get_attribute("id")
                    page.click(f'label[for="{rid}"]') if rid else r.check()
                    page.wait_for_timeout(300)
                    states = page.evaluate("""(k) => Array.from(document.querySelectorAll('.teaching-diagram [data-step]')).map(el => {
                        const cs = getComputedStyle(el);
                        const vis = !(cs.visibility==='hidden'||cs.display==='none'||parseFloat(cs.opacity)===0);
                        return {step: parseInt(el.dataset.step), visible: vis, expect: parseInt(el.dataset.step) <= k};
                      })""", int(val))
                    bad = [s for s in states if s["visible"] != s["expect"]]
                    result["steps"].append({"radio": rid, "value": val, "n_marked": len(states), "mismatch": bad, "ok": not bad})
                    page.screenshot(path=str(HERE / f"{name}_step{val}.png"), full_page=True)
            result["errors"].extend(f"[{key}] {e}" for e in errs)
            ctx.close()
        shell_path.unlink(missing_ok=True)

    result["ok"] = (
        all(v["no_h_overflow"] for v in result["viewports"].values())
        and not result["errors"]
        and all(h["ok"] for h in result["hover"])
        and all(s["ok"] for s in result["steps"])
        and not result["static"]["external_refs"]
        and not result["static"]["has_html_tags"]
        and result["static"]["has_viewbox"]
        and not result["static"]["unprefixed_ids"]
        and not result["static"].get("text_overflow")
    )
    (HERE / f"{name}.result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main(names: list[str]) -> int:
    rc = 0
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True, args=["--disable-gpu", "--no-sandbox"])
        for n in names:
            r = check_one(browser, n)
            flag = "OK " if r["ok"] else "FAIL"
            print(f"[{flag}] {n}: {r['fragment_kb']} KB, viewports={{{', '.join(k+':'+('ok' if v['no_h_overflow'] else 'OVERFLOW') for k,v in r['viewports'].items())}}}, "
                  f"hover={sum(h['ok'] for h in r['hover'])}/{len(r['hover'])}, steps={sum(s['ok'] for s in r['steps'])}/{len(r['steps'])}, "
                  f"errors={len(r['errors'])}, static={ {k:v for k,v in r['static'].items() if v} }")
            rc |= 0 if r["ok"] else 1
        browser.close()
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
