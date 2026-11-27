"""Check epidemiological reporting and preservation using the standard library."""
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import hashlib,json,re
def main():
    root=Path(__file__).resolve().parents[1]
    e=root/'epidemiology';v=json.loads((e/'validation.json').read_text());m=json.loads((e/'original_results_manifest.json').read_text())
    assert v['status']=='passed'
    for name,record in m['files'].items():
        b=(root/name).read_bytes()
        if hashlib.sha256(b).hexdigest()!=record['sha256']:
            assert hashlib.sha256(b.replace(b'\r\n',b'\n')).hexdigest()==record.get('canonical_sha256'),name
    assert v['original_scientific_files_preserved']==len(m['files'])
    readme=(root/'README.md').read_text(encoding='utf-8')
    for section in ['## Abstract','## Introduction','### Aim and objectives','## Methods','## Results','## Discussion','## Conclusion']:
        assert section in readme,section
    assert '@@' not in readme
    images=re.findall(r'!\[.*?\]\(([^)]+)\)',readme)
    assert len(images)==v['figures']
    for path in re.findall(r'\]\(([^)]+)\)',readme):
        if '://' not in path and not path.startswith('#'):
            assert (root/path.split('#')[0]).exists(),path
    word=re.search(r'\[Word report\]\(([^)]+)\)',readme).group(1)
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    with ZipFile(root/word) as z:
        doc=ET.fromstring(z.read('word/document.xml'));text=' '.join(doc.itertext())
        assert len(doc.findall('.//w:tbl',ns))==v['tables']
        assert len(doc.findall('.//w:drawing',ns))==v['figures']
        assert 'Vusumuzi Mabasa' in z.read('docProps/core.xml').decode()
        assert 'simulat' in text.lower() or 'synthetic' in text.lower()
        assert '@@' not in text
        assert all(edge.get('{'+ns['w']+'}val') in ('none','nil') for b in doc.findall('.//w:pBdr',ns) for edge in b)
        styles=ET.fromstring(z.read('word/styles.xml'))
        title=[s for s in styles.findall('w:style',ns) if s.get('{'+ns['w']+'}styleId')=='Title'][0]
        assert not title.findall('.//w:pBdr',ns)
    print(json.dumps(dict(status='passed',original_scientific_files=len(m['files']),figures=v['figures'],native_word_tables=v['tables'],word=word)))

if __name__ == "__main__":
    main()
