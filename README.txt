AAS77733 REPRODUCIBILITY PACKAGE

Manuscript:
A global test of an apparent mid-redshift residual feature in Pantheon+ Type Ia supernovae

This repository isolates the frozen code paths used for the revised AAS77733 analysis.

Primary publication source:
byDenoso/Pantheon@1e944e34f30d3d4948f60a6f045e5557d733c8b6

Additional R1 hardening source:
byDenoso/Pantheon@3aaec79fd6e6b43931a29d03afa4676c2163af6d

Reference runtime:
Python 3.13
numpy==2.2.6
scipy==1.15.3

Publication battery:
python -m unittest discover -s science/aas77733r1 -p 'test*.py' -v
python -m science.aas77733r1.run --mocks 10000 --seed 77733 --cache-dir .cache/aas77733 --output artifacts/aas77733-r1-publication.json

Hardening battery:
python -m unittest discover -s science/aas77733r1_hardening -p 'test*.py' -v
python -m science.aas77733r1_hardening.run --mocks 5000 --seed 77733 --cache-dir .cache/aas77733 --output artifacts/aas77733-r1-hardening-lm.json

Third-party datasets are not redistributed. Exact public inputs are pinned by commit and Git blob hash in science/aas77733/manifest.json.
