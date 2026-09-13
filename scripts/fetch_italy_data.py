#!/usr/bin/env python3
"""下載義大利預報實驗資料到 data/cache/italy/ 並建立乾淨快取。

來源：PyEEPAS 公開 repo（https://github.com/phonchi/EEPAS）data/README.md 所指的
Google Drive 資料夾（HORUS 目錄；Lolli et al. 2020）。區域檔 CELLE_ter.mat、
CPTI15.mat 直接讀 EEPAS/data/。

用法：
    python scripts/fetch_italy_data.py            # 下載 + 轉換快取
    python scripts/fetch_italy_data.py --forecasts  # 另外預先計算四張預報 × 40 窗
"""

from _teaching_runtime import prepare_plotting_environment
prepare_plotting_environment()

import argparse
import datetime as dt
import subprocess
import sys
import hashlib
import json
import shutil
import requests
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gdms_toolkit import italy  # noqa: E402

NEEDED = ["HORUS_Ita_Catalog.txt", "HORUS_Ita_DataOrigin.txt", "README.txt"]
REGION_COMMIT = "eafe229c23c4c788c5474a7e88e7e3fce07f07c6"


def prepare_regions():
    """Obtain the two public region files without requiring an untracked EEPAS checkout."""
    target = italy.ITALY_DIR / "regions"
    target.mkdir(parents=True, exist_ok=True)
    provenance = {}
    for name, attr in [("CELLE_ter.mat", "CELLS_MAT"), ("CPTI15.mat", "POLYGON_MAT")]:
        path = target / name
        url = f"https://raw.githubusercontent.com/phonchi/EEPAS/{REGION_COMMIT}/data/{name}"
        if not path.exists():
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            path.write_bytes(response.content)
        from scipy.io import loadmat
        loadmat(path)
        setattr(italy, attr, path)
        provenance[name] = dict(url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (target / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2)+"\n")


def download():
    italy.ITALY_DIR.mkdir(parents=True, exist_ok=True)
    missing = [f for f in NEEDED if not (italy.ITALY_DIR / f).exists()]
    if missing:
        print("下載 Google Drive 資料夾（需 gdown）…", flush=True)
        subprocess.run([sys.executable, "-m", "gdown", "--folder", italy.DRIVE_FOLDER,
                        "-O", str(italy.ITALY_DIR)], check=True)
    still = [f for f in NEEDED if not (italy.ITALY_DIR / f).exists()]
    if still:
        raise SystemExit(f"下載後仍缺檔案：{still}")
    stamp = italy.ITALY_DIR / "PROVENANCE.txt"
    if not stamp.exists():
        stamp.write_text(
            "來源：{}\n下載日期：{}\n目錄：HORUS (Lolli, Randazzo, Vannucci & Gasperini), "
            "詳見同資料夾 README.txt；資料依原提供者條款使用，不隨本 repo 散布。\n"
            "區域檔：EEPAS/data/CELLE_ter.mat、CPTI15.mat（Biondini et al. 2023）。\n".format(
                italy.DRIVE_FOLDER, dt.date.today().isoformat()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forecasts", action="store_true", help="預先計算並快取四張預報")
    ap.add_argument("--refresh", action="store_true", help="重建乾淨快取")
    args = ap.parse_args()
    download()
    prepare_regions()
    cat = italy.load_horus(refresh=args.refresh)
    print(f"HORUS 乾淨快取：{len(cat):,} 筆，{cat.time.min():%Y-%m-%d} – {cat.time.max():%Y-%m-%d}")
    exp = italy.experiment_catalog()
    print(f"實驗目錄（S 內、深度≤40、≤2021）：{len(exp):,} 筆；R 內 {exp.in_R.sum():,} 筆")
    for p in ("learning", "testing"):
        print(f"  {p} 目標地震（R 內、mb≥{italy.SPEC.mT}）：{len(italy.target_events(exp, p))}")
    if args.forecasts:
        from gdms_toolkit import italy_models
        italy_models.build_all_forecasts(exp)


if __name__ == "__main__":
    main()
