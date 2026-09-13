"""Normalize source-level emphasis while protecting code and math fences."""
from pathlib import Path
import argparse
import json
import re
import jupytext
from markdown_it import MarkdownIt

ROOT=Path(__file__).resolve().parents[1]
PARSER=MarkdownIt('commonmark')
PAIR=re.compile(r'\*\*([^*\n]+?)\*\*')
DEFINITION=re.compile(r'^(.*?)（([^（）]*[A-Za-z][^（）]*)）$')


def normalize(text):
    result=[];protected_fence=False;marker=''
    for line in text.splitlines(keepends=True):
        stripped=line.lstrip()
        opening=re.match(r'(`{3,}|~{3,})(.*)',stripped)
        if opening:
            token,info=opening.groups()
            if marker and token[0]==marker[0] and len(token)>=len(marker) and not info.strip():
                marker='';protected_fence=False
            elif not marker:
                marker=token
                protected_fence=not info.startswith(('{admonition','{note','{tip','{important','{warning','{seealso','{dropdown}'))
            result.append(line);continue
        if protected_fence:result.append(line);continue
        segments=re.split(r'(`+[^`]*`+)',line)
        for i in range(0,len(segments),2):
            def definition(match):
                m=DEFINITION.match(match[1])
                return f'**{m[1]}**（{m[2]}）' if m and m[1] and '$' not in m[1] else match[0]
            segments[i]=PAIR.sub(definition,segments[i])
            # Repair only delimiter pairs that CommonMark still leaves as literal text.
            bad={m[0] for m in PAIR.finditer(segments[i]) if m[0] in PARSER.renderInline(segments[i])}
            for token in bad:segments[i]=segments[i].replace(token,' '+token+' ')
        result.append(''.join(segments))
    return ''.join(result)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('pages',nargs='*');args=parser.parse_args()
    paths=[ROOT/'book'/p for p in args.pages] if args.pages else [*sorted((ROOT/'book').glob('*.py')),*sorted((ROOT/'book').glob('*.md'))]
    records=[]
    for path in paths:
        before=path.read_text()
        if path.suffix=='.py':
            nb=jupytext.read(path)
            changed=0
            for cell in nb.cells:
                if cell.cell_type=='markdown':
                    new=normalize(cell.source)
                    if new!=cell.source:changed+=1;cell.source=new
            if changed:jupytext.write(nb,path,fmt='py:percent')
        else:
            after=normalize(before);changed=int(after!=before)
            if changed:path.write_text(after)
        if changed:records.append({'file':str(path.relative_to(ROOT)),'markdown_cells_changed':changed})
    out=ROOT/'reference/notes/refresh_20260913_exhibition'
    out.mkdir(parents=True,exist_ok=True)
    (out/'format_repairs.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(records,ensure_ascii=False))


if __name__=='__main__':main()
