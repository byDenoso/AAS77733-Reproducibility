from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
v=json.load(open(r/'inputs/OFFICIAL_BLOB_VERIFICATION.json'))
for f in v['files']:
 b=(r/'inputs'/f['local_name']).read_bytes();assert hashlib.sha256(b).hexdigest()==f['sha256'];assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==f['official_tree_blob_sha1']
print('PASS: all four inputs match the pinned official repository blobs and SHA256.')
