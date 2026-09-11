"""Pair teaching sources, preserving outputs only for identical code sources.

Run with the project's Python. This does not execute notebooks or contact GDMS.
"""
from pathlib import Path
from collections import defaultdict, deque
import copy
import json
import os
import jupytext
import nbformat

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "book"


def sync():
    report = []
    for source in sorted(BOOK.glob("*.py")):
        if source.name.startswith("_"):
            continue
        paired = source.with_suffix(".ipynb")
        prior_ids = defaultdict(deque)
        if paired.exists():
            for old_cell in nbformat.read(paired, as_version=4).cells:
                if old_cell.get("id"):
                    prior_ids[(old_cell.cell_type, old_cell.source)].append(old_cell.id)
        cached = BOOK / "_build/jupyter_execute" / paired.name
        prior_outputs = {}
        for candidate in (cached, paired):
            if candidate.exists():
                prior = nbformat.read(candidate, as_version=4)
                for cell in prior.cells:
                    if cell.cell_type == "code" and cell.get("outputs"):
                        prior_outputs[cell.source.strip()] = cell
        notebook = jupytext.read(source)
        restored, missing = 0, []
        for i, cell in enumerate(notebook.cells):
            identities = prior_ids[(cell.cell_type, cell.source)]
            if identities:
                cell.id = identities.popleft()
            if cell.cell_type != "code":
                continue
            tags = list(cell.metadata.get("tags", []))
            tags = [tag for tag in tags if tag != "hide-input"]
            if "remove-input" not in tags:
                tags.append("remove-input")
            cell.metadata["tags"] = tags
            old = prior_outputs.get(cell.source.strip())
            if old is not None:
                cell.outputs = copy.deepcopy(old.outputs)
                cell.execution_count = old.get("execution_count")
                restored += 1
            elif "remove-output" not in tags:
                # Import/setup cells need no visible output. Any code change is
                # separately detected against cached code below.
                missing.append(i)
            if "remove-output" in tags:
                # Operational account/download responses do not belong in the
                # downloadable teaching notebook either. Source code remains.
                cell.outputs = []
        notebook.metadata["kernelspec"] = {
            "display_name": "Python 3", "language": "python", "name": "python3"
        }
        notebook.metadata["jupytext"] = {
            "formats": "ipynb,py:percent", "cell_metadata_filter": "tags,-all",
            "notebook_metadata_filter": "kernelspec,jupytext",
        }
        changed = True
        if cached.exists():
            old_codes = [c.source.strip() for c in nbformat.read(cached,4).cells if c.cell_type == "code"]
            new_codes = [c.source.strip() for c in notebook.cells if c.cell_type == "code"]
            changed = old_codes != new_codes
        nbformat.write(notebook, paired)
        jupytext.write(notebook, source, fmt="py:percent")
        report.append({"page": source.stem, "restored_output_cells": restored,
                       "code_changed_vs_build": changed, "cells_without_cached_output": missing})
    target = ROOT / os.environ.get("TEACHING_REPORT_DIR", "reference/notes/rewrite_20260911") / "notebook_sync.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    sync()
