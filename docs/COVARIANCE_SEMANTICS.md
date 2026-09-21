# Covariance semantics — official Pantheon+ source trace

Frozen upstream commit:
`PantheonPlusSH0ES/DataRelease@c447f0fea703fcd0fff57de5000947b5ca81286b`

Official covariance README:
`Pantheon+_Data/4_DISTANCES_AND_COVAR/README`

The README describes:
- `Pantheon+SH0ES_STATONLY.cov` as the statistical covariance product;
- `Pantheon+SH0ES_STAT+SYS.cov` as the statistical+systematic covariance product;
- STATONLY includes the statistical distance uncertainties and duplicate-SN off-diagonal intrinsic-scatter terms;
- STAT+SYS additionally includes covariance between SNe and Cepheid-host covariance from systematic uncertainties.

Released CosmoSIS likelihood:
`Pantheon+_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon+_only_cosmosis_likelihood.py`

The released likelihood sets `Pantheon+SH0ES_STAT+SYS.cov` as the default covariance input and loads that covariance in `build_covariance()`.

AAS77733 numerical reproduction:
- N = 1590
- Omega_m = 0.331576
- chi^2 = 1402.919

Matched referee-prescription sensitivity:
- released covariance directly: Delta chi^2 = 4.349
- released covariance + extra distance-modulus variance diagonal: Delta chi^2 = 1.376

Interpretive boundary:
The public-likelihood convention is the fiducial branch. The extra-diagonal construction is retained as a stress test. The stress test weakens the historical candidate and therefore does not rescue a physical-transition claim.
