"""Fetch immutable upstream inputs; verify before writing; never overwrite mismatches."""
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import quote
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def verify(data,entry):
    if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:
        raise ValueError('Upstream input SHA256/size mismatch: '+entry['local_name'])
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if blob!=entry['official_tree_blob_sha1']:raise ValueError('Git blob mismatch')
def main():
    manifest=json.loads((ROOT/'inputs/OFFICIAL_BLOB_VERIFICATION.json').read_text())
    for entry in manifest['files']:
        dest=ROOT/'inputs'/entry['local_name']
        if dest.exists():verify(dest.read_bytes(),entry);continue
        url='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/'+manifest['commit']+'/'+quote(entry['repository_path'],safe='/')
        with urlopen(url,timeout=180) as response:data=response.read()
        verify(data,entry)
        with dest.open('xb') as stream:stream.write(data)
        print('Verified',entry['local_name'])
if __name__=='__main__':main()
