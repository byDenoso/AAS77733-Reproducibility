"""Recompute AAS77733 numerical outputs and regenerate every manuscript asset.

Run from the REPRODUCIBILITY directory. This script never edits the scientific
manuscript text. It regenerates analysis arrays, summaries, figures, and the
LaTeX macro file that supplies numerical values to the manuscript.
"""
from pathlib import Path
import subprocess, sys, os
ROOT=Path(__file__).resolve().parent
A=ROOT/'analysis'; F=ROOT/'figures'
env=os.environ.copy(); env['OPENBLAS_NUM_THREADS']='1'
steps=[
    ([sys.executable,'check_inputs.py'], A/'input_hash_check_console.json'),
    ([sys.executable,'analysis/reproduce_calibration.py'], A/'console.jsonl'),
    ([sys.executable,'analysis/variance_addition_sensitivity.py'], A/'variance_sensitivity.jsonl'),
    ([sys.executable,'analysis/survey_bin_recalibration.py'], A/'survey_recalibration_console.txt'),
    ([sys.executable,'analysis/derive_metrics.py'], A/'derive_metrics_console.json'),
    ([sys.executable,'analysis/covariance_scale_sensitivity.py'], A/'covariance_scale_sensitivity_console.json'),
    ([sys.executable,'analysis/make_figures.py'], A/'figure_build_console.txt'),
    ([sys.executable,'analysis/plot_survey_diagnostics.py'], A/'figure6_build_console.txt'),
    ([sys.executable,'analysis/direct_gls_crosscheck.py'], A/'direct_gls_console.json'),
    ([sys.executable,'analysis/make_macros.py'], A/'macros_console.txt'),
]
for cmd,log in steps:
    print('RUN', ' '.join(cmd), flush=True)
    with log.open('w') as fh:
        subprocess.run(cmd,cwd=ROOT,env=env,stdout=fh,stderr=subprocess.STDOUT,check=True)
print('COMPLETE', ROOT/'generated_results.tex')
