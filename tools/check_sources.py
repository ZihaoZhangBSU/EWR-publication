"""Check syntax and source-manifest hashes without running model/data processing."""
from pathlib import Path
import ast, hashlib, json

ROOT=Path(__file__).resolve().parents[1]
def main():
    count=0
    for folder in ['analysis','workflow/src','workflow/tests','tools','tests']:
        for path in (ROOT/folder).rglob('*'):
            if path.suffix=='.py':code=[path.read_text(encoding='utf-8-sig')]
            elif path.suffix=='.ipynb':
                notebook=json.loads(path.read_text(encoding='utf-8'))
                code=[''.join(c['source']) for c in notebook['cells'] if c['cell_type']=='code']
            else:continue
            for i,source in enumerate(code):compile(source,f'{path.name}:cell{i}','exec')
            count+=1
    manifest=json.loads((ROOT/'docs/source_manifest.json').read_text())
    for record in manifest:
        actual=hashlib.sha256((ROOT/record['file']).read_bytes()).hexdigest()
        if actual!=record['packaged_sha256']:raise ValueError('Source manifest mismatch: '+record['file'])
    print(f'PASS: syntax for {count} scripts/notebooks; {len(manifest)} source fingerprints.')
if __name__=='__main__':main()
