"""Build contact sheets covering every reading page and every embedded diagram."""
from pathlib import Path
from PIL import Image,ImageDraw
import yaml
root=Path(__file__).parent
browser=root/'browser'
toc=yaml.safe_load(Path('book/_toc.yml').read_text())
pages=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']]
diagrams=sorted(p.stem for p in Path('book/_static/diagrams').glob('d*.html'))

def sheet(files,out,width,height,cols):
    canvas=Image.new('RGB',(width*cols,height*((len(files)+cols-1)//cols)),'#e5e7eb')
    draw=ImageDraw.Draw(canvas)
    for i,p in enumerate(files):
        with Image.open(p) as original:
            im=original.convert('RGB');im.thumbnail((width-12,height-32))
            x=(i%cols)*width;y=(i//cols)*height
            canvas.paste(im,(x+(width-im.width)//2,y+26))
            draw.text((x+6,y+7),p.stem[:52],fill='#243040')
    canvas.save(out)

for device in ['desktop-light','mobile-light','desktop-dark','mobile-dark']:
    w,h,cols=(400,284,3) if device.startswith('desktop') else (195,456,6)
    sheet([browser/f'{device}-{p}-top.png' for p in pages],root/f'{device}-all-pages.png',w,h,cols)
    sheet([browser/f'{device}-{d}.png' for d in diagrams],root/f'{device}-all-diagrams.png',384,610,3)
