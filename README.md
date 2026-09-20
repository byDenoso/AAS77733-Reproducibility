# AAS77733 Reproducibility Package

Reproducibility repository for the revised manuscript:

**A global test of an apparent mid-redshift residual feature in Pantheon+ Type Ia supernovae**

This repository isolates the code paths and tests used to support the revised AAS77733 analysis from the larger `byDenoso/Pantheon` development repository.

## Scope

The package contains:

- the publication battery `G0–G19`;
- the R1 hardening checks `H0–H6`;
- immutable upstream data manifests;
- exact source provenance;
- GitHub Actions workflows used for the publication and hardening batteries;
- documentation describing the relation between tests, data inputs, and the manuscript.

The scientific conclusion supported by the revised analysis is a **non-detection of a robust localized residual feature** under the tested statistic and null hierarchy. This repository is intended to reproduce the tests, not to introduce a new model or claim.

## Frozen upstream provenance

Primary R1 implementation:

`byDenoso/Pantheon@1e944e34f30d3d4948f60a6f045e5557d733c8b6`

Additional R1 hardening implementation:

`byDenoso/Pantheon@3aaec79fd6e6b43931a29d03afa4676c2163af6d`

The hardening commit is one commit ahead of the primary frozen baseline and adds the `aas77733r1_hardening` layer without rewriting the frozen G0–G19 implementation.

## Quick start

Python 3.13 is the reference runtime used by the publication workflows.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

python -m unittest science.aas77733r1.test_contract
python -m unittest science.aas77733r1_hardening.test_hardening
```

For the full publication battery, use the commands encoded in the workflows under `.github/workflows/`.

## Data

The repository does not vendor third-party Pantheon+ or DES-SN5YR data products. The exact public upstream repositories, commits, paths, and Git blob hashes are pinned in:

`science/aas77733/manifest.json`

The data loader verifies immutable Git blob hashes before caching downloaded inputs.

## Important numerical checkpoints

The revised manuscript reports, among other checks:

- official Pantheon+ likelihood: Ωm = 0.331576, χ² = 1402.919;
- fiducial non-calibrator sample with zHD > 0.0233: N = 1365;
- blind-scan maximum: Δχ² = 2.224 at zt = 0.03;
- global covariance-null calibration: p = 0.8378 from 10,000 mocks;
- historical zt = 0.36 diagnostic: Δχ² = 1.841.

These are checkpoints for reproducibility, not independent claims of significance.

## Repository status

This repository was assembled from the exact frozen source files referenced by the revised manuscript. A successful file-level provenance check does **not** by itself imply that every scientific workflow has been executed in this repository. Workflow/run status should be inspected separately.

## Citation and archival release

A DOI-bearing archival release, e.g. Zenodo, should be created from the final accepted reproducibility snapshot. The manuscript data/software availability section should then cite that DOI.

## License

No new software license is asserted here until the license status of the upstream code is explicitly resolved.
