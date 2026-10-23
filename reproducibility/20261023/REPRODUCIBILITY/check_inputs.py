"""Verify full SHA-256 identities of the two public Pantheon+ source products."""
from pathlib import Path
import hashlib, json, os, sys
ROOT = Path(__file__).resolve().parent
EXPECTED = {
    'Pantheon+SH0ES.dat': '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8',
    'Pantheon+SH0ES_STAT+SYS.cov': 'abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc',
}
def validate(root=None):
    root = Path(root or os.environ.get('PANTHEON_DATA_DIR', ROOT/'inputs'))
    errors=[]; observed={}
    for name, expected in EXPECTED.items():
        path=root/name
        if not path.is_file():
            errors.append(f'missing: {path}')
            continue
        with path.open('rb') as handle:
            sha=hashlib.file_digest(handle,'sha256').hexdigest()
        observed[name]=sha
        if sha!=expected:
            errors.append(f'SHA-256 mismatch: {name}: got {sha}, expected {expected}')
    return {'ok':not errors,'directory':str(root),'files_checked':len(observed),'sha256':observed,'errors':errors}
if __name__ == '__main__':
    result=validate();print(json.dumps(result,indent=2));sys.exit(0 if result['ok'] else 1)
