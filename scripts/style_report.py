"""Apply explicit, accessible table formatting after Quarto DOCX rendering."""
import sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as E

W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
E.register_namespace('w',W)
def tag(s): return '{'+W+'}'+s
def put(parent,name,**attrs):
    node=parent.find(tag(name))
    if node is None: node=E.Element(tag(name));parent.append(node)
    if name in ('rPr','pPr','tblPr','tcPr','trPr'):
        parent.remove(node);parent.insert(0,node)
    for key,value in attrs.items(): node.set(tag(key),str(value))
    return node

path=Path(sys.argv[1])
with zipfile.ZipFile(path) as z: parts={name:z.read(name) for name in z.namelist()}
root=E.fromstring(parts['word/document.xml'])
for table in root.iter(tag('tbl')):
    rows=table.findall(tag('tr'));nc=len(rows[0].findall(tag('tc')))
    header=' '.join(t.text or '' for t in rows[0].iter(tag('t')))
    if 'Assembly accession' in header:widths=[1050,2550,1200,1900,2660]
    else:widths={2:[4000,5360],3:[3120,3120,3120],4:[2250,2500,2500,1200],5:[2250,1777,1777,1778,1778]}[nc]
    props=put(table,'tblPr');put(props,'tblW',w=sum(widths),type='dxa');put(props,'tblLayout',type='fixed')
    borders=put(props,'tblBorders')
    for edge in ('top','left','bottom','right','insideH','insideV'): put(borders,edge,val='single',sz=4,color='D9D9D9')
    grid=put(table,'tblGrid');grid.clear()
    for width in widths: E.SubElement(grid,tag('gridCol'),{tag('w'):str(width)})
    for i,row in enumerate(rows):
        rp=put(row,'trPr');put(rp,'cantSplit')
        if i==0: put(rp,'tblHeader')
        for j,cell in enumerate(row.findall(tag('tc'))):
            cp=put(cell,'tcPr');put(cp,'tcW',w=widths[j],type='dxa');put(cp,'vAlign',val='center')
            put(cp,'shd',fill='243746' if i==0 else ('F2F5F7' if i%2==0 else 'FFFFFF'),val='clear')
            margins=put(cp,'tcMar')
            for side in ('top','bottom','left','right'):put(margins,side,w=80,type='dxa')
            for p in cell.iter(tag('p')):
                pp=put(p,'pPr');put(pp,'jc',val='left' if j==0 else 'center');put(pp,'spacing',before=0,after=0,line=240,lineRule='auto')
                if len(rows)<=14:put(pp,'keepNext',val=1 if i<len(rows)-1 else 0)
                for run in p.findall(tag('r')):
                    rpr=put(run,'rPr');put(rpr,'color',val='FFFFFF' if i==0 else '000000');put(rpr,'sz',val=20)
                    if i==0:put(rpr,'b')
parts['word/document.xml']=E.tostring(root,encoding='utf-8',xml_declaration=True)
tmp=path.with_suffix('.styled.docx')
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
    for name,data in parts.items():z.writestr(name,data)
tmp.replace(path)
print('Styled result tables with fixed widths, repeat headers and visible borders.')
