"""Run in a fresh destination, preserving the source and frozen evidence."""
from pathlib import Path
import argparse,shutil,subprocess,sys,os

def main():
    p=argparse.ArgumentParser();p.add_argument('destination',type=Path)
    p.add_argument('--mode',choices=['observed','figures','full'],default='observed');a=p.parse_args()
    source=Path(__file__).resolve().parents[1];dest=a.destination.resolve()
    if dest.exists():raise SystemExit('Destination must not exist.')
    dest.mkdir(parents=True)
    for name in ['code','inputs','expected','tests']:
        shutil.copytree(source/name,dest/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    (dest/'outputs').mkdir();(dest/'figures').mkdir()
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(script,*args):subprocess.run([sys.executable,str(dest/'code'/script),*args],cwd=dest,env=env,check=True)
    run('fetch_inputs.py');run('verify_inputs.py');run('independent_step.py','observed')
    if a.mode=='full':
        for mode in ['pilot','baseline','controls','injection']:run('independent_step.py',mode)
        run('validate_full_draws.py');run('summarize.py');run('make_figures.py')
    if a.mode=='figures':
        import json
        o=dest/'outputs';expected=dest/'expected'
        for n,v in json.loads((expected/'null_calibrations.json').read_text()).items():
            (o/n).mkdir();(o/n/'summary.json').write_text(json.dumps(v)+'\n')
        for name in ['injection.json','direct40_validation.json']:shutil.copyfile(expected/name,o/name)
        run('summarize.py');run('make_figures.py')
    print('Reproduction directory:',dest)
if __name__=='__main__':main()
