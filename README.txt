AAS77733 REPRODUCIBILITY PACKAGE

Manuscript:
A global test of an apparent mid-redshift residual feature in Pantheon+ Type Ia supernovae

Purpose:
This repository is the frozen reproducibility package for the revised AAS77733 analysis. It supports a non-detection / identifiability result for the historical z~0.4 residual candidate.

Reference runtime:
Python 3.13
numpy==2.2.6
scipy==1.15.3

Frozen implementation provenance:
Primary publication source:
byDenoso/Pantheon@1e944e34f30d3d4948f60a6f045e5557d733c8b6

Additional R1 hardening source:
byDenoso/Pantheon@3aaec79fd6e6b43931a29d03afa4676c2163af6d

Public data provenance:
PantheonPlusSH0ES/DataRelease@c447f0fea703fcd0fff57de5000947b5ca81286b
  Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat
  Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov
  Pantheon+_Data/4_DISTANCES_AND_COVAR/README
  Pantheon+_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon+_only_cosmosis_likelihood.py

DES external comparison:
des-science/DES-SN5YR@c9a4fcafc4cbd19bd750dee47fc76194a45c181f

Third-party datasets are not redistributed. Exact paths, commits and Git blob hashes are pinned in science/aas77733/manifest.json and the accompanying provenance files.

Publication battery:
python -m unittest discover -s science/aas77733r1 -p 'test*.py' -v
python -m science.aas77733r1.run --mocks 10000 --seed 77733 --cache-dir .cache/aas77733 --output artifacts/aas77733-r1-publication.json

Hardening battery:
python -m unittest discover -s science/aas77733r1_hardening -p 'test*.py' -v
python -m science.aas77733r1_hardening.run --mocks 5000 --seed 77733 --cache-dir .cache/aas77733 --output artifacts/aas77733-r1-hardening-lm.json

Headline checkpoints:
Official Pantheon+ likelihood reproduction:
  N = 1590
  Omega_m = 0.331576
  chi2 = 1402.919

Fiducial sample:
  N = 1365
  zHD > 0.0233

Blind scan:
  Delta chi2 = 2.224
  z_t = 0.03

Historical candidate:
  z_t = 0.36
  Delta chi2 = 1.841

Global covariance null:
  mocks = 10000
  p_global = 0.8378
  q95(T) = 9.007

Classification:
NON_DETECTION

Repository layout:
README.txt                         human-readable archive guide
requirements.txt                   pinned numerical dependencies
CITATION.cff                       citation metadata
science/aas77733/manifest.json     immutable public data provenance
science/aas77733r1/                G0-G19 publication battery
science/aas77733r1_hardening/      H0-H6 R1 hardening layer
docs/DATA_SOURCES.md               public dataset locations
docs/PROVENANCE.md                 frozen implementation provenance
docs/COVARIANCE_SEMANTICS.md       covariance convention trace
docs/REFEREE_TRACEABILITY.md       reviewer concern -> test mapping
docs/RESPONSE_TO_REFEREE_R1.md     working point-by-point response
docs/MANUSCRIPT_R2_PATCH.md         final manuscript patch record
artifacts/expected_headline_results.json
                                   manuscript numerical checkpoints

Archival release / DOI:
The final journal-upload snapshot must be deposited in a DOI-issuing repository such as Zenodo and submitted to the AAS Journals Community. The DOI must then be inserted into the manuscript using the AASTeX \dataset[]{} macro.

This repository does not invent or predeclare a DOI. Use the metadata template in ZENODO_METADATA.json for the archival deposit.

AAS manuscript ID:
AAS77733

Author:
Dener Pereira
Instituto de Fisica, Universidade Federal do Rio de Janeiro (UFRJ), Rio de Janeiro, Brazil
