# AAS77733 latest source-only reproduction snapshot

This version contains the original Python scripts from the supplied latest pre-submission package for **No Significant Residual Step across 62 Tested Pantheon+ Redshift Partitions**. Numerical algorithms have not been edited here. The manifest binds the original source bytes to this snapshot. Documentation is adapted for this code-only distribution.

## Scope and status

This is source code, not the full archival package. Stored arrays, figures, CSV products, journal source, third-party AAS classes and editorial correspondence are deliberately excluded. No new science run or Monte Carlo replay was performed for this source import. Syntax inspection is not scientific validation.

Zenodo DOI **10.5281/zenodo.23306992** is reserved; the deposit was not published at the time this snapshot was prepared. Do not describe it as a publicly available archival release. No new software or data license is assigned by this update; the repository's existing license-pending status is preserved.

## Inputs and environment

Use Python 3.13 and the package versions in `requirements.txt`. Obtain the two original Pantheon+ files from the pinned official upstream revision in `inputs/README.txt`, and put them in this layer's `inputs/` directory. Check their SHA-256 identities before running.

The commit `c447f0fea703fcd0fff57de5000947b5ca81286b` identifies **PantheonPlusSH0ES/DataRelease**, not this analysis repository. Some scripts accept environment overrides, but several require local `inputs/`; therefore use the local directory for a complete run.

From this `REPRODUCIBILITY` directory, in a separate working copy:

```
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python check_inputs.py
OPENBLAS_NUM_THREADS=1 python run_all.py
OPENBLAS_NUM_THREADS=1 python analysis/equal_count_bins.py
```

These commands perform real numerical calculations, overwrite locally generated outputs, and may take substantial time. They were not executed as part of this GitHub preparation. The equal-count script retains its supplied fixed fitted background; it is a descriptive display, not a new global test.

Additional diagnostic scripts are supplied separately: `analysis/survey_deletion_influence.py`, `analysis/diagnostics.py` and `analysis/parameter_projection.py`. They are not executed by `run_all.py`. Run them only after the baseline outputs exist and with an appropriate scientific purpose.

## Archived versus regenerated verification

`verify_reproduction.py` is retained unchanged for users who also have the complete journal/archive assets. It requires the sibling `JOURNAL_UPLOAD`, generated/archived outputs and `analysis/r1_equal_count_reference.csv`, none of which is included in this source-only layer. It is **not** a standalone test for this checkout. Its `archived_integrity_only` result does not mean Monte Carlo was rerun.

`sync_journal_assets.py` is likewise optional and requires generated figures. In the supplied version it copies Figures 1-6 and numeric macros, not Supplementary Figure S1. Do not treat it as synchronization of every manuscript asset.

## Data redistribution boundary

The local numerical pipeline produces `result_arrays.npz`, which includes selected observed magnitudes, redshifts and covariance, in addition to calculated quantities. Such arrays are not included here. They must not be committed or redistributed without checking the applicable upstream rights and approved publication scope. The local ignore rules reduce accidental inclusion but do not replace this check.

The covariance-scale routine estimates a scalar from the observed sample and holds it fixed in matched simulations. Its result remains a sensitivity diagnostic, not a fully calibrated treatment of covariance-estimation uncertainty.
