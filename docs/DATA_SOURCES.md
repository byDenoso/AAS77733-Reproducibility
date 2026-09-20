# Data sources

Third-party Pantheon+ and DES-SN5YR data products are not vendored here.

Exact public inputs are pinned in `science/aas77733/manifest.json` and, for additional hardening inputs, in `science/aas77733r1_hardening/hardening.py`.

Principal upstream repositories:

- PantheonPlusSH0ES/DataRelease @ `c447f0fea703fcd0fff57de5000947b5ca81286b`
- des-science/DES-SN5YR @ `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`

The data loader verifies Git blob SHA-1 values before caching downloads.
