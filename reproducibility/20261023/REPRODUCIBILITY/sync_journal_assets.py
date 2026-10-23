"""Copy regenerated numeric macros/figures into the sibling journal source and verify SHA-256 identity."""
from pathlib import Path
import hashlib, json, shutil
ROOT=Path(__file__).resolve().parent
J=ROOT.parent/'JOURNAL_UPLOAD'
F=ROOT/'figures'
assets=['fig1_residuals.pdf','fig2_scan_and_null.pdf','fig3_power.pdf','fig4_survey_robustness.pdf','fig5_background_response.pdf','fig6_survey_composition.pdf']
J.mkdir(parents=True,exist_ok=True)
shutil.copy2(ROOT/'generated_results.tex',J/'results.tex')
for x in assets: shutil.copy2(F/x,J/x)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest={'results.tex':{'generated':sha(ROOT/'generated_results.tex'),'journal':sha(J/'results.tex')},'figures':{}}
for x in assets: manifest['figures'][x]={'generated':sha(F/x),'journal':sha(J/x)}
manifest['all_match']=all(v['generated']==v['journal'] for v in manifest['figures'].values()) and manifest['results.tex']['generated']==manifest['results.tex']['journal']
(ROOT/'journal_asset_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
if not manifest['all_match']: raise SystemExit('asset mismatch')
