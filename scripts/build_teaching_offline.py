"""Build the book using saved notebook outputs, without data/API execution."""
from pathlib import Path
import subprocess
import sys
import yaml

ROOT=Path(__file__).resolve().parents[1]
BOOK=ROOT/'book'
config=yaml.safe_load((BOOK/'_config.yml').read_text())
config['execute']['execute_notebooks']='off'
for name,path in config.get('sphinx',{}).get('local_extensions',{}).items():
    config['sphinx']['local_extensions'][name]=str((BOOK/path).resolve())
settings=config.get('sphinx',{}).get('config',{})
if 'html_static_path' in settings:
    settings['html_static_path']=[str((BOOK/path).resolve()) for path in settings['html_static_path']]
override=BOOK/'_build/offline_config.yml'
override.parent.mkdir(parents=True,exist_ok=True)
override.write_text(yaml.safe_dump(config,allow_unicode=True,sort_keys=False))
raise SystemExit(subprocess.call([
    str(Path(sys.executable).parent/'jupyter-book'),'build',str(BOOK),
    '--config',str(override),'--all','-W','--keep-going',
],cwd=ROOT))
