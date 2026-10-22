# Independent Pantheon+ mean-step replication, 21 October 2026

Companion numerical implementation for **No Globally Significant Redshift Step in the Pantheon+ Type Ia Supernova Hubble Diagram**. This layer supersedes the older R1 headline calibration for the current manuscript; old implementations remain intact for historical provenance. It is not a journal submission or a DOI archive.

## Scientific scope

The 62-pivot mean-step statistic is calibrated with a refitted background under a conditional Gaussian plug-in null. The separate joint-background alternative, all 1,090 admissible partitions, redshift cuts, survey intercepts, sparse high-redshift bins and injection power are explicitly distinguished.

- 1,365 selected light curves: zHD > 0.0233 and IS_CALIBRATOR = 0.
- Primary observed T = 2.224483643; p = 0.653834617 (10,000 baseline draws).
- Joint observed T = 3.771316803; p = 0.619638036.
- All-partition p = 0.428257174.
- Exploratory maximum over 11 simultaneous GLS bins: p = 0.012098790. This is not a global correction for all prior data inspection.
- One-bin global/survey-intercept p = 0.018298170 / 0.065967016; nuisance freedom changes uncertainty and calibration, not evidence that it explains away the approximately +0.239 mag offset.
- Historical weak-step injection detection power: 8.45%. Non-detection does not exclude weak steps.

No discovery claim is made. The design is retrospective, not prospective preregistration. Calibration is conditional on the specified plug-in Gaussian null, not a uniform-coverage guarantee. Native covariance rounding asymmetry is symmetrized as (C + C.T)/2; display-error variances are not added to the official full covariance. Extra-diagonal runs are stress tests only.

## Environment and inputs

Reference environment: CPython 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Install this layer's requirements in a separate environment, not the older R1 root requirements. No cosmology backend, credentials or paid service is required.

```sh
python3.12 -m venv .venv-pantheon
. .venv-pantheon/bin/activate
python -m pip install -r reproducibility/20261021/requirements.txt
python reproducibility/20261021/code/fetch_inputs.py
python reproducibility/20261021/code/verify_inputs.py
python -m unittest discover -s reproducibility/20261021/tests -v
```

Inputs are fetched from PantheonPlusSH0ES/DataRelease at immutable commit c447f0fea703fcd0fff57de5000947b5ca81286b. Four exact paths, sizes, SHA-256 and Git blob hashes are pinned in inputs/OFFICIAL_BLOB_VERIFICATION.json. Third-party data and code are not redistributed here; retain upstream attribution and applicable terms. A mismatch fails closed and is never silently replaced.

## Reproduction without overwriting evidence

Choose a destination that does not exist:

```sh
# Observed-only numerical regression, no Monte Carlo rerun:
python reproducibility/20261021/code/reproduce.py /tmp/pantheon-observed --mode observed
# Regenerate figures using fresh observed quantities and frozen calibration summaries:
python reproducibility/20261021/code/reproduce.py /tmp/pantheon-figures --mode figures
# Explicitly rerun the complete simulation battery:
python reproducibility/20261021/code/reproduce.py /tmp/pantheon-full --mode full
```

Full mode generates 20,000 main null draws plus 40 pilot draws, five paired injection contrasts reusing 2,000 noise vectors (10,000 fits, not 10,000 new independent vectors), 40 direct full-vector checks, summaries and figure arrays/PDFs. Seeds and controls are encoded in code/independent_step.py. BLAS runs with one thread. Runtime is hardware-dependent. The historical 10,000-draw comparison bank and editorial package are not dependencies of this implementation.

Expected JSON files preserve numerical reference outputs from the completed 21 October run. They are regression fixtures, not evidence that a new run occurred. Figures mode deliberately reuses these frozen summaries and must not be described as fresh Monte Carlo validation. Full mode produces fresh simulated statistics. Floating-point differences across platforms are possible; observed regression tolerance is 1e-5.

The immutable scientific engine, direct validator and plotting source match the recorded 21 October implementation (see SOURCE_PROVENANCE.json). New packaging scripts isolate public reproduction from editorial material. No private correspondence, reviewer reports, personal data, unrelated models or credentials are included. No new license is asserted; the repository's LICENSE_PENDING.md remains applicable.
