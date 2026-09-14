from pathlib import Path
from PIL import Image,ImageDraw
import yaml
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
T=yaml.safe_load((ROOT/'book/_toc.yml').read_text());pages=[T['root']]+[c['file'] for p in T['parts'] for c in p['chapters']]
def montage(files,target,w,h,columns):
    canvas=Image.new('RGB',(w*columns,h*((len(files)+columns-1)//columns)),'#e9e8e1');draw=ImageDraw.Draw(canvas)
    for i,f in enumerate(files):
        im=Image.open(f).convert('RGB');im.thumbnail((w-12,h-30));x=i%columns*w;y=i//columns*h
        canvas.paste(im,(x+(w-im.width)//2,y+25));draw.text((x+6,y+7),f.stem[:55],fill='#15333c')
    canvas.save(target)
for device in ['desktop-light','mobile-light','desktop-dark','mobile-dark']:
    w,h,cols=(400,310,3) if device.startswith('desktop') else (195,452,6)
    montage([OUT/'browser'/f'{device}-{p}-top.png' for p in pages],OUT/f'{device}-all-pages.png',w,h,cols)
old=ROOT/'reference/notes/refresh_20260913_exhibition/browser_final/desktop-light-11_conditional_intensity-top.png'
new=OUT/'browser/desktop-light-11_conditional_intensity-top.png'
montage([old,new],OUT/'sidebar-before-after.png',720,530,2)
