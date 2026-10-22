"""Collect current numerical outputs, without private editorial provenance."""
import json
from independent_step import ROOT,save_json
out=ROOT/'outputs'
obs=json.loads((out/'observed_all.json').read_text())
names=['baseline','cut001','cut0023','cut0024','cut01','survey']
null={n:json.loads((out/n/'summary.json').read_text()) for n in names}
b=obs['baseline'];s=obs['survey'];i=9
high={'simultaneous_GLS':b['simultaneous_bins']['amplitudes'][i],'simultaneous_sigma':b['simultaneous_bins']['sigma'][i], 'one_bin_amplitude':b['bins']['amplitudes'][i],'one_bin_sigma':b['bins']['sigma'][i], 'survey_one_bin_amplitude':s['bins']['amplitudes'][i],'survey_one_bin_sigma':s['bins']['sigma'][i]}
save_json(ROOT/'RESULT.json',{'observed':obs,'null_calibrations':null,'injection_results':json.loads((out/'injection.json').read_text()),'direct_validation':json.loads((out/'direct40_validation.json').read_text()),'diagnostics':{'high_z_bin':high}})
