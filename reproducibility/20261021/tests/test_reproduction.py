from pathlib import Path
import sys,json,unittest,hashlib,tempfile,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from independent_step import Engine,mc_summary,vector_minimize
from fetch_inputs import verify
class Contracts(unittest.TestCase):
    def test_mc_plus_one(self):
        r=mc_summary(np.array([0.,1.,2.]),1.);self.assertEqual(r['hits'],2);self.assertEqual(r['p_plus_one'],.75)
    def test_input_corruption_rejected(self):
        b=b'abc';e={'bytes':3,'sha256':hashlib.sha256(b).hexdigest(),'official_tree_blob_sha1':hashlib.sha1(b'blob 3\0'+b).hexdigest(),'local_name':'test'}
        verify(b,e)
        with self.assertRaises(ValueError):verify(b'abd',e)
    def test_optimizer_includes_endpoints(self):
        targets=np.array([-2.,.2,2.]);x,_,_=vector_minimize(lambda a:(a-targets)**2,(3,))
        np.testing.assert_allclose(x,[-1,.2,1],atol=1e-7)
    def test_existing_destination_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=subprocess.run([sys.executable,str(ROOT/'code/reproduce.py'),d],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0);self.assertIn('must not exist',p.stderr)
    def test_observed_baseline(self):
        if not (ROOT/'inputs/Pantheon+SH0ES.dat').exists():self.skipTest('Run fetch_inputs.py for data-dependent regression')
        e=Engine(all_partitions=True);a=e.observed();b=json.loads((ROOT/'expected/observed_all.json').read_text())['baseline']
        self.assertEqual(a['n'],1365);self.assertEqual(a['n_all_partitions'],1090)
        for k in ['omega','primary_T','joint_T','all_partitions_T','simultaneous_bin_T']:
            self.assertAlmostEqual(a[k],b[k],delta=1e-5,msg=k)
        np.testing.assert_allclose(a['scan']['T'],b['scan']['T'],atol=1e-5,rtol=1e-6)
    def test_frozen_calibration_counts(self):
        d=json.loads((ROOT/'expected/null_calibrations.json').read_text())
        self.assertEqual(sum(x['primary_T']['N'] for x in d.values()),20000)
        self.assertAlmostEqual(d['baseline']['primary_T']['p_plus_one'],.6538346165383462)
if __name__=='__main__':unittest.main()
