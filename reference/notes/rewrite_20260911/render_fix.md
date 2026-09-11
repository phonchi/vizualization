# 重複 inline JavaScript：根因與局部修正

日期：2026-09-11。症狀：乾淨建置的每頁仍有兩次 `const THEBE_JS_URL`，浏览器宣告衝突；`togglebuttonSelector` 也重複。不是 notebook cache 或 doctree 重用所致。

## 原因證據

本機套件皆在 `.venv/lib/python3.13/site-packages/`：

1. `sphinx/application.py:1069–1074`，`app.add_js_file` 同時送入 registry 及已建立的 builder。
2. `sphinx_thebe/__init__.py:82` 在 `env-before-read-docs` 回呼加入 filename=None 的 inline const；其註冊在 235–238 行。`sphinx_togglebutton/__init__.py:44` 同樣在 builder-inited 加入 filename=None 的 selector。
3. `sphinx/builders/html/__init__.py:516–517` 在 prepare_writing 再次取 registry 並以 `js_filename or ''` 加入；此時相同 inline 的 filename 變成空字串。
4. `sphinx/builders/html/_assets.py:94–103` 的 `_JavaScript.__eq__` 同時比較 filename、priority 及 attributes；None 與 '' 不相等，後續 dict 去重也不能消除。
5. 唯讀實驗構造 `_JavaScript(None, body=...)` 與 `_JavaScript('', body=...)`，相等比較 False，dict 去重長度仍為 2。
6. `get_final_config` 實際產生的 extension 名單沒有重複。toggleHintShow/Hide 在 config-inited、builder 建立前加入，因此沒有同樣的雙份 late registration 症狀。

## 修正

新增 `book/_ext/teaching_render.py`，透過公開 `html-page-context` hook，在模板呈現之前去掉完全相同的 inline asset。正規化 None/空字串；鍵同時保留全部 attributes（含 body/type 等）與 priority；保留第一次出現位置。外連腳本及不同屬性的 inline 不動，不修改套件或所有 extension 註冊表。

由主代理在 `sphinx.local_extensions` 設定 `teaching_render: _ext`。注意 Jupyter Book 把這個路徑相對於「實際傳入的 config YAML 所在目錄」解析；離線腳本把 YAML 寫到 book/_build 時，須在寫出前把 local_extensions 路徑轉成以原始 book 目錄為基準的絕對路徑。一般直接 `jupyter-book build book` 則使用原 `_config.yml`，不需另處理。

## 驗證與限制

以真實 Sphinx `_JavaScript` 在記憶體驗證：None/空字串相同 inline 合併；不同 body、type 或 priority 保留；外連腳本保留；重複執行結果不變。未在此子任務建置、變更 config 或處理 dropdown 語法；由主代理配置後重建並檢查瀏覽器，雙層 toggle 是否另有獨立問題尚待驗收。
