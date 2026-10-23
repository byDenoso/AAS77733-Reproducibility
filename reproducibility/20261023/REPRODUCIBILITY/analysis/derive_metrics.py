"""Derive all secondary manuscript/figure quantities from computed artifacts."""
import json, os
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular, cho_factor, cho_solve
from scipy.integrate import cumulative_trapezoid, quad

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / 'analysis'
a = np.load(A / 'result_arrays.npz')
data = np.genfromtxt(ROOT / 'inputs' / 'Pantheon+SH0ES.dat', names=True, dtype=None, encoding='utf8')
idx = a['selected_idx']
z = a['z']; zh = a['zh']; C = a['C']; residual = a['residual']

# Unique supernovae.
unique_count = int(len(set(np.asarray(data['CID'][idx]).tolist())))

# Covariance-weighted descriptive bin means used in Figure 1.
be = np.array([.0233,.05,.1,.2,.3,.4,.5,.65,.8,1,1.3,1.8,2.5])
bix = np.digitize(z, be, right=True) - 1
X = np.column_stack([bix == j for j in range(len(be)-1)]).astype(float)
cf = cho_factor(C, lower=True, check_finite=False)
ones = np.ones(len(z)); ci1 = cho_solve(cf, ones, check_finite=False)
intercept = float((ci1 @ residual) / (ci1 @ ones))
r0 = residual - intercept
V = np.linalg.inv(X.T @ cho_solve(cf, X, check_finite=False))
means = V @ (X.T @ cho_solve(cf, r0, check_finite=False))
errs = np.sqrt(np.diag(V)); counts = X.sum(axis=0).astype(int)

# Covariance-metric overlap between Omega_m response and step at z=0.36.
L = cholesky(C, lower=True, check_finite=False)
a1 = solve_triangular(L, np.ones(len(z)), lower=True, check_finite=False)
a1 = a1 / np.linalg.norm(a1)
def white(v):
    q = solve_triangular(L, v, lower=True, check_finite=False)
    return q - a1 * (a1 @ q)
def mu(omega):
    g = np.linspace(0, max(z)*1.001, 8192)
    v = 1/np.sqrt(omega*(1+g)**3 + 1-omega)
    integ = cumulative_trapezoid(v, g, initial=0)
    return 5*np.log10((299792.458/70)*(1+zh)*np.interp(z,g,integ)) + 25
om = float(a['om_hat'])
u = white((mu(om+1e-5)-mu(om-1e-5))/(2e-5))
s36 = white((z > .36).astype(float))
rho = float((u @ s36) / np.sqrt((u @ u)*(s36 @ s36)))
se_inflation = float(1/np.sqrt(1-rho*rho))

# Read primary and joint output records.
records = [json.loads(line) for line in (A/'console.jsonl').read_text().splitlines() if line.strip().startswith('{')]
rec = {r['stage']:r for r in records}
joint_om = float(rec['joint_background_alternative']['best_omega_m'])
step43_idx = int(np.argmin(abs(a['pivots']-.43)))
step43_amp = float(a['step_amp'][step43_idx])

# Analytic response from baseline Omega_m to the jointly optimized alternative.
grid = np.linspace(.005, 1.2, 900)
def distmod_grid(omega):
    E = np.sqrt(omega*(1+grid)**3 + 1-omega)
    d = cumulative_trapezoid(1/E, grid, initial=0)
    d += quad(lambda x:1/np.sqrt(omega*(1+x)**3+1-omega), 0, .005, epsabs=1e-13)[0]
    return 5*np.log10((299792.458/70)*(1+grid)*d) + 25
m0 = distmod_grid(om); m1 = distmod_grid(joint_om)
shift = -(m1-m0); shift -= np.interp(.08, grid, shift)
shift43 = float(np.interp(.43, grid, shift))
distance_fraction = float((np.log(10)/5)*step43_amp)

# Deliberate covariance-addition stress test.
var_records = [json.loads(line) for line in (A/'variance_sensitivity.jsonl').read_text().splitlines() if line.strip().startswith('{')]
var = {r['stage']:r for r in var_records}
var_obs = var['observed']; var_null = var['null_calibration']['refit_background_each_mock_then_fixed_scan']

out = {
    'unique_supernovae': unique_count,
    'display_bins': {
        'edges': be.tolist(), 'counts': counts.tolist(), 'means_mag': means.tolist(), 'sigma_mag': errs.tolist(),
        'bin_03_04': {'mean_mag':float(means[4]), 'sigma_mag':float(errs[4]), 'count':int(counts[4])},
        'bin_04_05': {'mean_mag':float(means[5]), 'sigma_mag':float(errs[5]), 'count':int(counts[5])},
        'bin_10_13': {'mean_mag':float(means[9]), 'sigma_mag':float(errs[9]), 'count':int(counts[9])},
    },
    'geometry': {
        'omega_baseline': om, 'omega_joint': joint_om, 'rho_step036': rho,
        'profile_se_inflation': se_inflation, 'step43_amplitude_mag': step43_amp,
        'residual_shift_at_043_mag': shift43, 'step43_distance_fraction': distance_fraction,
        'anchor_redshift': 0.08,
    },
    'variance_addition': {
        'omega_m': float(var_obs['omega_m']), 'chi2': float(var_obs['chi2']),
        'T': float(var_obs['T']), 'pivot': float(var_obs['best_z']),
        'p': float(var_null['p_plus_one']), 'se': float(var_null['mc_se']), 'q95': float(var_null['q95']),
    }
}
(A/'derived_summary.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out))
