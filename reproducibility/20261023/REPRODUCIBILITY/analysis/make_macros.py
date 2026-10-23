"""Generate manuscript macros only from recomputed structured results."""
from pathlib import Path
import json, numpy as np
ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'analysis'
rows=[json.loads(line) for line in (A/'console.jsonl').read_text().splitlines() if line.strip().startswith('{')]
t={x['stage']:x for x in rows}
a=np.load(A/'result_arrays.npz')
d=json.loads((A/'derived_summary.json').read_text())
s9=json.loads((A/'survey9_summary.json').read_text()); s20=json.loads((A/'survey20_summary.json').read_text()); sb=json.loads((A/'baseline_summary.json').read_text())
s=t['observed']; n=t['null_calibration']; inj=t['refitted_injection']; j=t['joint_background_alternative']
fmtp=lambda x:f"{x:.4f}" if x<0.1 else f"{x:.3f}"
vals={
'Nrows':f"{s['N']:,}", 'Nunique':f"{d['unique_supernovae']:,}", 'Nsplit':str(len(a['pivots'])),
'OmegaFit':f"{s['omega_m']:.4f}", 'ChiBaseline':f"{s['chi2']:.2f}",
'Tmaximum':f"{s['T']:.3f}",'BestPivot':f"{s['best_z']:.2f}",'BestAmplitude':f"{s['step_amp']:+.4f}",
'Ttarget':f"{s['target036']:.3f}",'Twindow':f"{s['target_window_T']:.3f}",
'WindowPivot':'0.43','WindowAmplitude':f"{float(a['step_amp'][np.argmin(abs(a['pivots']-.43))]):+.4f}",
'Pglobal':f"{n['refit_background_each_mock_then_fixed_scan']['p_plus_one']:.3f}",
'PglobalError':f"{n['refit_background_each_mock_then_fixed_scan']['mc_se']:.3f}",
'PrimaryExceedances':str(n['refit_background_each_mock_then_fixed_scan']['exceedances']),
'Tthreshold':f"{n['refit_background_each_mock_then_fixed_scan']['q95']:.3f}",
'Pfixed':f"{n['original_fixed_background']['p_plus_one']:.3f}", 'PfixedError':f"{n['original_fixed_background']['mc_se']:.3f}",
'TfixedThreshold':f"{n['original_fixed_background']['q95']:.3f}",
'Pjoint':f"{j['calibration']['p_plus_one']:.3f}",'PjointError':f"{j['calibration']['mc_se']:.3f}",
'Tjoint':f"{j['observed_T']:.3f}",'JointPivot':f"{j['best_z']:.2f}",'JointOmega':f"{j['best_omega_m']:.4f}",
'TjointThreshold':f"{j['calibration']['q95']:.3f}",
'MidPower':f"{100*inj['results'][1]['detection_fraction']:.2f}",'MidLoc':f"{100*inj['results'][1]['localization_fraction_pm003']:.2f}",
'LowPower':f"{100*inj['results'][0]['detection_fraction']:.2f}",'LowLoc':f"{100*inj['results'][0]['localization_fraction_pm003']:.2f}",
'DisplayMeanA':f"{d['display_bins']['bin_03_04']['mean_mag']:.4f}",'DisplaySigmaA':f"{d['display_bins']['bin_03_04']['sigma_mag']:.4f}",
'DisplayMeanB':f"{d['display_bins']['bin_04_05']['mean_mag']:.4f}",'DisplaySigmaB':f"{d['display_bins']['bin_04_05']['sigma_mag']:.4f}",
'HighDisplayMean':f"{d['display_bins']['bin_10_13']['mean_mag']:.3f}",'HighDisplaySigma':f"{d['display_bins']['bin_10_13']['sigma_mag']:.3f}",'HighDisplayCount':str(d['display_bins']['bin_10_13']['count']),
'RhoStep':f"{d['geometry']['rho_step036']:.3f}",'ProfileInflation':f"{d['geometry']['profile_se_inflation']:.2f}",
'GeometryShift':f"{d['geometry']['residual_shift_at_043_mag']:.4f}",'DistancePct':f"{100*d['geometry']['step43_distance_fraction']:.2f}",
'SurveyTwentyOmega':f"{s20['best_omega_m']:.4f}",'SurveyTwentyChi':f"{s20['no_step_chi2']:.2f}",'SurveyTwentyT':f"{s20['step_search']['T']:.3f}",'SurveyTwentyPivot':f"{s20['step_search']['pivot']:.2f}",'SurveyTwentyP':f"{s20['step_search']['p']:.3f}",'SurveyTwentyPErr':f"{s20['step_search']['se']:.3f}",'SurveyTwentyQ':f"{s20['step_search']['q95']:.3f}",
'SurveyNineT':f"{s9['step_search']['T']:.3f}",'SurveyNinePivot':f"{s9['step_search']['pivot']:.2f}",'SurveyNineP':f"{s9['step_search']['p']:.3f}",'SurveyNinePErr':f"{s9['step_search']['se']:.3f}",
'BinOneT':f"{sb['bin_family']['T']:.3f}",'BinOneP':f"{sb['bin_family']['p']:.4f}",'BinOnePErr':f"{sb['bin_family']['se']:.4f}",'BinOneAmp':f"{sb['bin_family']['high_bin_amplitude_mag']:.3f}",'BinOneAmpErr':f"{sb['bin_family']['high_bin_conditional_sigma_mag']:.3f}",
'BinNineT':f"{s9['bin_family']['T']:.3f}",'BinNineP':f"{s9['bin_family']['p']:.4f}",'BinNinePErr':f"{s9['bin_family']['se']:.4f}",
'BinTwentyT':f"{s20['bin_family']['T']:.3f}",'BinTwentyP':f"{s20['bin_family']['p']:.4f}",'BinTwentyPErr':f"{s20['bin_family']['se']:.4f}",'BinTwentyAmp':f"{s20['bin_family']['high_bin_amplitude_mag']:.3f}",'BinTwentyAmpErr':f"{s20['bin_family']['high_bin_conditional_sigma_mag']:.3f}",
'VarianceOmega':f"{d['variance_addition']['omega_m']:.4f}",'VarianceChi':f"{d['variance_addition']['chi2']:.2f}",'VarianceT':f"{d['variance_addition']['T']:.3f}",'VariancePivot':f"{d['variance_addition']['pivot']:.2f}",'VarianceP':f"{d['variance_addition']['p']:.4f}",'VariancePErr':f"{d['variance_addition']['se']:.4f}",
}
out=ROOT/'generated_results.tex'; out.write_text(''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in vals.items()))
print(out.read_text())
