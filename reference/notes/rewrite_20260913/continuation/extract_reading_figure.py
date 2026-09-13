"""Extract the user-supplied published Fig. 8 for the annotated reading example."""
from pathlib import Path
import pymupdf as fitz
import hashlib,json
source=Path('reference/[2023] Italy_EEPAS.pdf')
target=Path('book/_static/reading/biondini2023_fig8.png');target.parent.mkdir(exist_ok=True)
doc=fitz.open(source)
clip=fitz.Rect(102,39,474,214)
doc[11].get_pixmap(matrix=fitz.Matrix(3,3),clip=clip).save(target)
report=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),page=12,printed_page=1692,figure=8,clip=list(clip),scale=3,output=str(target),output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),citation='Biondini, Rhoades & Gasperini (2023), GJI 234, doi:10.1093/gji/ggad123')
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
