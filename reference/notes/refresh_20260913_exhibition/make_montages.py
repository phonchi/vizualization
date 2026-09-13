from pathlib import Path
from PIL import Image,ImageDraw
import yaml,json
root=Path(__file__).parent;images=root/'browser_final'
toc=yaml.safe_load(Path('book/_toc.yml').read_text());pages=[toc['root']]+[c['file'] for p in toc['parts'] for c in p['chapters']]
slugs=sorted(p.stem for p in Path('book/_static/diagrams').glob('d*.html'))

def montage(files,dest,size,columns):
    width,height=size;canvas=Image.new('RGB',(width*columns,height*((len(files)+columns-1)//columns)),'#e9e8e1');draw=ImageDraw.Draw(canvas)
    for i,p in enumerate(files):
        with Image.open(p) as original:
            image=original.convert('RGB');image.thumbnail((width-12,height-30))
            x=i%columns*width;y=i//columns*height;canvas.paste(image,(x+(width-image.width)//2,y+25));draw.text((x+6,y+7),p.stem[:54],fill='#18363d')
    canvas.save(dest)
for device in ['desktop-light','mobile-light','desktop-dark','mobile-dark']:
    montage([images/f'{device}-{p}-top.png' for p in pages],root/f'{device}-all-pages.png',(400,310) if device.startswith('desktop') else (195,452),3 if device.startswith('desktop') else 6)
    montage([images/f'{device}-{p}.png' for p in slugs],root/f'{device}-all-exhibits.png',(390,660),3)
comparison=root/'comparison';comparison.mkdir(exist_ok=True)
old=Path('reference/notes/rewrite_20260913/continuation/browser')
montage([old/'desktop-light-00_intro-top.png',images/'desktop-light-00_intro-top.png'],comparison/'homepage-before-after.png',(720,530),2)
montage([old/'mobile-light-d12_etas_branching.png',images/'mobile-light-d12_etas_branching.png'],comparison/'branching-before-after.png',(390,740),2)
(root/'montage_inventory.json').write_text(json.dumps({'reading_pages':pages,'exhibits':slugs,'full_screenshots':str(images),'before_after':str(comparison)},ensure_ascii=False,indent=2)+'\n')
