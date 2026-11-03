"""Verify the portable release manifest without importing model dependencies.
Editor: Wenjie Liu; Date: 2026-10-31.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();root=args.root.resolve()
    manifest=json.loads((root/'docs/release/RELEASE_MANIFEST.json').read_text(encoding='utf-8'))
    failures=[]
    for name,expected in manifest['files_sha256'].items():
        path=root/name
        if not path.resolve().is_relative_to(root) or not path.is_file():
            failures.append({'file':name,'problem':'missing or outside root'});continue
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=expected:failures.append({'file':name,'problem':'SHA256 mismatch'})
    print(json.dumps({'status':'FAIL' if failures else 'PASS','release':manifest['release'],
        'checked_files':len(manifest['files_sha256']),'failures':failures},ensure_ascii=False,indent=2))
    return 1 if failures else 0


if __name__=='__main__':sys.exit(main())
