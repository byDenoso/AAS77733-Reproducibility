# AAS77733 — Response to Reviewer 1 (working resubmission text)

Manuscript: **Global calibration and identifiability of an apparent mid-redshift residual feature in Pantheon+ Type Ia supernovae**

We thank the referee for the detailed report. The manuscript has been substantially rewritten rather than locally patched. The revised paper is shorter, removes the previous "inference-audit framework" terminology, introduces Pantheon+ and the SN Ia context explicitly, defines the candidate split and the global statistic before use, moves interpretation and limitations into a single Discussion, removes the previous appendices, expands the relevant literature, and narrows the scientific claim to a falsification/identifiability result.

The central conclusion has also changed materially. The revised analysis does **not** claim evidence for a physical transition near z~0.4. The historical candidate is retained only as provenance for why the audit was performed. In the revised fiducial analysis the historical location is weak, the blind maximum is elsewhere, the global probability is large, and the preferred pivot is unstable.

## 1. Existence and provenance of the z~0.4 feature

**Referee concern:** the manuscript did not show the feature clearly, did not establish where it came from in the literature or data, and could be read as if a published Pantheon+ feature had already been established.

**Response:** We agree that the original wording was ambiguous. The revised manuscript now states explicitly that the target arose from an earlier exploratory, post-hoc internal analysis and is not a previously published Pantheon+ detection. The Abstract wording is changed to:

> "I test a candidate split near z~0.4 identified in an earlier exploratory analysis, treating both its existence and location as unconfirmed."

The Introduction now states that the earlier exploratory score is retained only as provenance for the candidate and does not define the confirmatory test. The blind analysis does not fix z_t~0.36. For the fiducial non-calibrator sample (N=1365, z_HD>0.0233), the blind maximum is Delta chi^2=2.224 at z_t=0.03, the historical z_t=0.36 location gives Delta chi^2=1.841, and the maximum in 0.34<=z_t<=0.45 is Delta chi^2=2.037 at z_t=0.43. The globally calibrated probability is p_global=0.8378.

A revised observational panel is also required in the resubmission figure set: Hubble residual versus redshift, with low-opacity individual points, binned means with uncertainties, and the historical z=0.36 candidate marked as a provenance reference rather than as a detected transition. The caption will state explicitly that this panel is descriptive and is not the significance statistic.

## 2. Pantheon+ context and terminology

**Referee concern:** Pantheon+ was not adequately introduced; SN Ia/cosmology language was sparse; several statistical terms were undefined or field-inappropriate.

**Response:** The revised Introduction now begins from SN Ia cosmology rather than from a statistical framework. Pantheon+ is described as 1701 light curves of 1550 spectroscopically confirmed SNe Ia, with citations to Scolnic et al. (2022) and Brout et al. (2022). SNe Ia are introduced as standardizable candles and as a principal low-redshift distance-redshift probe. The previous "inference-audit framework" and "detector" wording has been removed. Generalized least squares is written out before the abbreviation is used, z_t is defined as a candidate split redshift, and the hard step is described as a minimal diagnostic of a mean residual change, not as a physical phase-transition model.

The previous broad references to voids, dust, and external tracers that lacked a direct role in the main claim have been removed from the main analysis.

## 3. Low-redshift cut and calibrator exclusion

**Referee concern:** z_HD>0.01 was weakly motivated and calibrator removal was unexplained.

**Response:** The fiducial cut is now z_HD>0.0233, matching the modern Hubble-flow motivation highlighted by the referee. The text explains that nearby SNe are more sensitive to peculiar velocities and frame corrections, and the minimum-redshift choice is varied explicitly as a sensitivity branch rather than selected after looking at the feature.

Cepheid calibrators are removed from the fiducial residual sample to form a homogeneous Hubble-flow diagnostic sample. The revised text explicitly states that this does not imply that calibrator SNe are defective; they belong to the distance-ladder branch because they carry independent host-distance information.

## 4. Step model, statistic, and look-elsewhere effect

**Referee concern:** the step template, z_t, the statistic, and the look-elsewhere correction were insufficiently motivated or defined.

**Response:** The revised Methods now defines the residual vector, the intercept-only reference model, the step vector s_i(z_t), the generalized least-squares estimator, and the local improvement Delta chi^2(z_t). The discovery statistic is defined as

T = max_{z_t in Z} Delta chi^2(z_t),

over the declared scan. The paper explicitly states that the hard step is a minimal change-point diagnostic and not a physical transition model.

Because the pivot is selected from the data, the significance statement is based on the null distribution of the maximum T from the complete scan rather than on the local value at a selected z_t. The fiducial global calibration uses 10,000 covariance realizations and gives p_global=0.8378; the 95% null quantile is T=9.007.

## 5. Critical Concern 1 — covariance semantics

**Referee concern:** the released Pantheon+ STAT+SYS matrix might require an additional diagonal variance term, C_tot=C_stat+sys+diag(sigma_mu^2), and therefore the original analysis might have used the covariance incorrectly.

**Response:** We treated this as a decisive methodological issue and rebuilt the covariance validation before interpreting the historical candidate.

The revised analysis reproduces the public Pantheon+ likelihood for its 1590-row selection, obtaining

Omega_m = 0.331576,
chi^2 = 1402.919.

The official Pantheon+ DataRelease identifies `Pantheon+SH0ES_STAT+SYS.cov` as the statistical+systematic covariance product. More importantly, the released CosmoSIS likelihood uses that file directly as its default covariance input:

`Pantheon+_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon+_only_cosmosis_likelihood.py`

at frozen upstream commit

`PantheonPlusSH0ES/DataRelease@c447f0fea703fcd0fff57de5000947b5ca81286b`.

The corresponding DataRelease README states that the STATONLY matrix contains the statistical distance uncertainties (including duplicate-SN off-diagonal intrinsic-scatter terms) and that STAT+SYS additionally contains the systematic covariance between SNe and Cepheid-host covariance. We therefore use the released STAT+SYS matrix directly in the fiducial likelihood convention.

We nevertheless implemented the referee's proposed additional-diagonal prescription as a matched sensitivity branch rather than dismissing it. In that branch the split gain decreases from Delta chi^2=4.349 under the released covariance to 1.376 after adding the additional distance-modulus variance diagonal. Thus the alternative convention weakens the candidate further and does not restore a physical-transition interpretation.

The revised manuscript now cites the exact public likelihood implementation and the DataRelease covariance README, rather than relying only on the phrase "inspection and numerical reproduction."

## 6. Structured nulls and survey organization

**Referee concern:** shuffle/null terminology was unclear and the paper did not explain what the tests answered.

**Response:** Each null family is now tied to a specific question. The baseline global null draws Gaussian residual realizations from the Pantheon+ covariance and repeats the full blind scan. Additional nulls preserve progressively more catalogue organization: intra-survey shuffles, redshift shuffles, and survey-by-redshift block permutations. They are not treated as independent experiments to be combined.

The empirical probabilities are 0.9023 for intra-survey shuffles, 0.9307 for redshift shuffles, and 0.9908-0.9960 for the survey-by-redshift block families. The result is therefore not merely "non-significant"; maxima at least as large as observed are routine once survey/redshift structure is preserved.

## 7. Results presentation and figure captions

**Referee concern:** figures/tables were hard to interpret and captions were insufficient.

**Response:** The revised manuscript reduces the result set to one headline table and four figures with self-contained captions. The old redundant tables and appendices have been removed. The revised figures state the tested quantity, the numerical landmarks, and the intended inference directly in the caption.

The remaining addition before resubmission is the descriptive Hubble-residual-versus-redshift panel requested above, which will make the small visual scale of the historical candidate explicit rather than forcing the reader to infer it from the scan statistic alone.

## 8. Results / Interpretation / Limitations structure

**Referee concern:** Results, Interpretation, and Limitations were fragmented.

**Response:** These have been reorganized into Results followed by a single Discussion with subsections for the artifact-class interpretation, relation to current SN-systematics work, covariance sensitivity, unresolved questions, and limitations. The previous appendices have been deleted.

## 9. Literature context

**Referee concern:** the paper was disconnected from the Pantheon+ and broader SN Ia literature.

**Response:** The revised Introduction and Discussion now explicitly contextualize the analysis with the Pantheon+ cosmology and redshift/peculiar-velocity papers, earlier SN Ia internal-robustness work, recent Pantheon+ anisotropy/alternative-background analyses, Pantheon+ vs DES-SN5YR comparisons, and Dovekie/DES calibration work. The manuscript also states what these papers do *not* establish: none provides prior evidence for a localized physical change point at z~0.4.

## 10. Critical Concern 2 — writing and accessibility

**Referee concern:** the manuscript was overcomplicated, jargon-heavy, and insufficiently connected to SN Ia cosmology.

**Response:** The paper has been substantially rewritten and shortened. The statistical workflow is now introduced only after the physical data context. Terms that carried project-specific meanings have been removed. The revised paper avoids the prior "framework/package/detector/layer" vocabulary in the scientific narrative and keeps provenance details in the reproducibility material rather than in the main prose.

## 11. Critical Concern 3 — evidence for the claimed feature

**Referee concern:** the manuscript lacked evidence that a z~0.4 feature exists and could be interpreted as building a statistical analysis around an unsupported anomaly.

**Response:** The revised paper agrees with the referee's substantive concern. The current full-covariance, blind, globally calibrated analysis does not recover a significant z~0.4 feature. The paper is now framed as the falsification of a historical exploratory candidate, not as evidence for a new Pantheon+ anomaly. The historical location gives Delta chi^2=1.841, p_global for the full blind maximum is 0.8378, the pivot is broad and unstable, and the DES-Dovekie branch does not recover a stable counterpart.

The conclusion is therefore bounded: the tested data/statistic combination does not identify a physical transition epoch near z~0.4. It also does not exclude arbitrarily weak physical structure below the resolving power of the statistic.

## 12. Data Editor review

We have replaced the previous future-tense Data Availability language with a public reproducibility repository:

https://github.com/byDenoso/AAS77733-Reproducibility

The repository contains README.txt, requirements.txt, pinned upstream repositories and commits, immutable Git blob hashes for third-party inputs, source-provenance records, the G0-G19 publication battery, the H0-H6 hardening layer, test-to-referee traceability, and headline numerical checkpoints. Third-party Pantheon+ and DES-SN5YR data are not redistributed; the exact public input locations and hashes are recorded.

Before final resubmission, the final tagged package will be archived in a DOI-issuing repository and submitted to the AAS Journals Zenodo community. The DOI will be inserted in the manuscript using the AASTeX `\\dataset[]{}` macro. No DOI is claimed in the current working text until that deposit exists.

