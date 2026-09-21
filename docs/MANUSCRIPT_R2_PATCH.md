# AAS77733 — Manuscript R2 patch set

This file records the remaining text/figure edits to apply to the current 6-page revised manuscript before resubmission.

## 1. Abstract provenance patch

Replace:

> I test a previously reported exploratory split near z~0.4 in Pantheon+ without fixing either the existence or the location of a transition.

With:

> I test a candidate split near z~0.4 identified in an earlier exploratory analysis, treating both its existence and location as unconfirmed.

Reason: removes the implication that a published Pantheon+ paper previously reported the feature.

## 2. Covariance-validation patch

Replace the first paragraph of Section III with:

> The covariance interpretation was checked before evaluating the historical candidate. Reproducing the public Pantheon+ likelihood gives Omega_m=0.331576 and chi^2=1402.919 for the 1590-row public selection. The Pantheon+ DataRelease identifies `Pantheon+SH0ES_STAT+SYS.cov` as the statistical+systematic covariance product, and the released CosmoSIS Pantheon+ likelihood passes that file directly as its default covariance input (PantheonPlusSH0ES/DataRelease, commit c447f0fea703fcd0fff57de5000947b5ca81286b, `Pantheon+_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon+_only_cosmosis_likelihood.py`). The accompanying DataRelease README states that the STATONLY matrix contains the statistical distance uncertainties, including duplicate-SN off-diagonal intrinsic-scatter terms, while the STAT+SYS matrix additionally contains the systematic covariance between SNe and Cepheid-host covariance. We therefore use the released STAT+SYS matrix directly in the fiducial likelihood convention rather than adding a second distance-modulus variance diagonal.

Keep the existing sensitivity paragraph immediately after it:

> The alternative prescription remains useful as a stress test. In a matched covariance-sensitivity configuration, the direct released covariance gives Delta chi^2=4.349, while adding the additional distance-modulus variance diagonal reduces the value to 1.376. The convention changes the numerical score but not the direction of the inference: the more conservative prescription weakens the apparent split.

## 3. New observational figure panel

Add a descriptive panel before the current discovery-statistic figure:

**Panel:** Hubble residual vs redshift for the fiducial non-calibrator sample (N=1365, z_HD>0.0233).

Required content:
- individual residual points with low opacity;
- binned residual means with uncertainties as a visual guide;
- vertical marker at z=0.36 labelled "historical exploratory candidate";
- optional horizontal mean levels on either side of z=0.36;
- no spline/smoothed curve;
- no significance annotation derived from this visual panel.

Caption text:

> **Observed Hubble residuals around the historical candidate.** Individual Pantheon+ Hubble residuals for the fiducial non-calibrator sample are shown with low opacity, with binned means and their uncertainties overlaid only as a visual guide. The vertical line marks the historical exploratory candidate at z=0.36. No visually sharp discontinuity is present; the separation is small compared with the object-to-object scatter. This panel is descriptive only. Statistical evidence is quantified by the covariance-aware blind scan and its globally calibrated null distribution, not by visual inspection of the binned residuals.

## 4. Data and software availability replacement

Replace the current section with:

> **Data and software availability.** Pantheon+ distance and covariance products are publicly available from the Pantheon+SH0ES DataRelease (Scolnic et al. 2022; Brout et al. 2022). The exact data table and covariance products used here are pinned to `PantheonPlusSH0ES/DataRelease@c447f0fea703fcd0fff57de5000947b5ca81286b`, with paths and Git blob hashes recorded in the reproducibility package. DES-Dovekie Hubble-diagram products are likewise pinned to `des-science/DES-SN5YR@c9a4fcafc4cbd19bd750dee47fc76194a45c181f`. Revision-specific code, frozen source provenance, dependency versions, test manifests, and rerun commands are available at `https://github.com/byDenoso/AAS77733-Reproducibility`. A DOI-bearing archival snapshot of the final resubmission package will be cited here with the AASTeX `\\dataset[]{}` macro once the DOI deposit has been created. No private data are required for the headline results.

After Zenodo deposition, replace the final sentence about the future DOI with the actual `\\dataset[<label>]{<DOI>}` citation. Do not leave future-tense DOI wording in the submitted version.

## 5. Figure ordering

Recommended final order:
1. Observed Hubble residuals around z=0.36 + blind-scan landmarks/covariance sensitivity.
2. Null hierarchy and independent sensitivity forks.
3. Pivot identifiability and survey-model sensitivity.
4. Observed-scale power and DES-Dovekie external fork.

This directly answers the referee's complaint that the claimed feature was never shown in the data.

## 6. Stop rule

No additional scientific battery is required unless one of the above edits reveals an inconsistency with the frozen headline checkpoints. The remaining publication blockers are documentary/figure-level:
- render the residual panel from the frozen fiducial sample;
- create the DOI-bearing archival snapshot;
- insert the DOI in Data Availability;
- compile the revised source;
- perform final text/figure/citation QA.
