"""Process-local plotting setup for the teaching execution/build commands."""
from pathlib import Path
import os
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def prepare_plotting_environment():
    """Choose a genuinely writable cache before any Matplotlib import occurs."""
    configured = os.environ.get('MPLCONFIGDIR')
    candidates = [Path(configured).expanduser()] if configured else []
    candidates.append(ROOT / 'data' / 'cache' / 'matplotlib')
    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryFile(dir=candidate):
                pass
        except OSError:
            continue
        os.environ['MPLCONFIGDIR'] = str(candidate.resolve())
        os.environ.setdefault('MPLBACKEND', 'Agg')
        return candidate
    raise RuntimeError('No writable Matplotlib cache; prepare data/cache/matplotlib before execution')
