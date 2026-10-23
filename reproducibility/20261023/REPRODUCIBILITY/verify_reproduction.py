"""Submission/reproduction integrity checks that do not rerun the Monte Carlo calculations."""
from pathlib import Path
import json, hashlib, re, sys
ROOT=Path(__file__).resolve().parent; A=ROOT/'analysis'; J=ROOT.parent/'JOURNAL_UPLOAD'
errors=[]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
required=[A/'console.jsonl',A/'result_arrays.npz',A/'survey20_summary.json',A/'survey9_summary.json',A/'baseline_summary.json',A/'derived_summary.json',A/'direct_gls_check.json',A/'covariance_scale_sensitivity.json',ROOT/'generated_results.tex']
for p in required:
    if not p.exists() or not p.stat().st_size: errors.append(f'missing/empty: {p}')
# Generated assets must be byte-identical to the journal copy.
for name in ['fig1_residuals.pdf','fig2_scan_and_null.pdf','fig3_power.pdf','fig4_survey_robustness.pdf','fig5_background_response.pdf','fig6_survey_composition.pdf','fig7_equal_count_bins.pdf']:
    a=ROOT/'figures'/name; b=J/name
    if not a.exists() or not b.exists() or sha(a)!=sha(b): errors.append(f'journal figure not synchronized: {name}')
if sha(ROOT/'generated_results.tex')!=sha(J/'results.tex'): errors.append('results.tex not synchronized')
# Direct GLS check must agree at tight numerical precision.
c=json.loads((A/'direct_gls_check.json').read_text())
if abs(c['step_vs_reduced_difference'])>1e-9 or abs(c['bin_vs_reduced_difference'])>1e-9: errors.append('direct GLS crosscheck exceeds 1e-9')
# Figure-build scripts must not embed the reported scientific result literals.
pat=re.compile(r'0\.652|0\.606|0\.0187|0\.0700|8\.47|2\.224|6\.225|0\.625|0\.843|0\.9716|0\.239|0\.077|3\.492|3\.597|9\.323|7\.065')
for p in A.glob('*.py'):
    if pat.search(p.read_text()): errors.append(f'hard-coded reported result in {p.name}')
# AAS revision copies require line numbers in the source used for the editorial PDF.
tex=(J/'manuscript.tex').read_text()
if 'linenumbers' not in tex.splitlines()[0]: errors.append('line numbers missing in referee source')

if 'equal-occupancy' not in tex or not (A/'equal_count_bins.py').exists(): errors.append('R1 complementary visual analysis absent')
if 'tab:influence' not in tex or not (A/'survey_deletion_influence.py').exists(): errors.append('survey influence table or reproduction code missing')
# Compare newly regenerated GLS quantile-bin output against exact R1 ledger CSV.
import csv, math
bin_ref=A/'r1_equal_count_reference.csv'
bin_new=ROOT/'figures'/'fig7_equal_count_values.csv'
if not bin_ref.is_file() or not bin_new.is_file():
    errors.append('missing R1 equal-count source/result CSV')
else:
    with bin_ref.open(newline='') as fa, bin_new.open(newline='') as fb:
        ar=list(csv.DictReader(fa)); br=list(csv.DictReader(fb))
    if len(ar)!=12 or len(br)!=12: errors.append('invalid equal-count bin length')
    elif any(set(a)!=set(b) or any(not math.isclose(float(a[k]),float(b[k]),rel_tol=0,abs_tol=1e-12) for k in a) for a,b in zip(ar,br)):
        errors.append('R1 recovered equal-count values differ')
from check_inputs import validate
inputs=validate()
if any('SHA-256 mismatch' in e for e in inputs['errors']): errors.extend(inputs['errors'])
print(json.dumps({'ok':not errors,'archived_integrity_only':True,'source_inputs_present_and_verified':inputs['ok'],'inputs_required_for_recompute':True,'errors':errors,'input_status':inputs['errors']},indent=2))
if errors: sys.exit(1)
