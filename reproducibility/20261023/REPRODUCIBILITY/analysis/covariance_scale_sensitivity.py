"""Matched global scalar-covariance sensitivity for the frozen Pantheon+ step test.

C' = a C with a = chi2(no step)/(N-2) taken from the *observed* sample.
The relative data weighting and fitted best Omega_m are unchanged, while the
observed statistic is T/a. The Gaussian response manifold and search scores
are rescaled appropriately before refitting Omega_m within each null draw.
The same archived standard-Gaussian score realizations are reused to enable an
exact paired comparison at fixed mock noise; they are NOT new observations.
This is a conditional plug-in stress test. Its p-value does not account for
estimating a from the data and does not validate structured covariance changes.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
A=Path(__file__).resolve().parent
q=np.load(A/'result_arrays.npz',allow_pickle=False)
N=len(q['z']);degrees=N-2
alpha=float(q['best_chi'])/degrees
v=q['injection_noise'];omega=q['om_grid']; rank=int(q['rank'])
B=q['B']/np.sqrt(alpha); D=q['D']/alpha; cross=q['cross']/np.sqrt(alpha)
noiseB=v[:,:rank]@B;noiseS=v[:,rank:]
loss=D[None,:]-2*noiseB
pick=np.argmin(loss,axis=1);mid=np.clip(pick,1,len(omega)-2);row=np.arange(len(v))
L=loss[row,mid-1];M=loss[row,mid];R=loss[row,mid+1]
den=L-2*M+R
frac=np.divide(.5*(L-R),den,out=np.zeros(len(v)),where=np.abs(den)>1e-14)
frac=np.clip(frac,-1,1);frac[(pick==0)|(pick==len(omega)-1)]=0
om=omega[mid]+frac*(omega[1]-omega[0]);om[pick==0]=omega[0];om[pick==len(omega)-1]=omega[-1]
scores=noiseS-CubicSpline(omega,cross,axis=0)(om)
null_T=np.max(scores*scores,axis=1)
obs=float(np.max(q['scan_T']))/alpha
exceed=int((null_T>=obs).sum())
p=(exceed+1)/(len(null_T)+1)
se=float(np.sqrt((exceed/len(v))*(1-exceed/len(v))/len(v)))
result={
  'sensitivity_only':True,
  'covariance_model':'a_times_released_STAT+SYS',
  'scale_estimated_from_observed_data':True,
  'reestimate_scale_inside_each_null':False,
  'reused_archived_standard_gaussian_scores':True,
  'N':N,'nominal_dof':degrees,'chi2_no_step':float(q['best_chi']),
  'scale_a':alpha,'observed_T':obs,'mocks':len(v),
  'exceedances':exceed,'p_plus_one':p,'mc_se':se,
  'null_q95':float(np.quantile(null_T,.95)),
  'omega_mock_min':float(np.min(om)),'omega_mock_max':float(np.max(om)),
  'scope':'Global scalar rescaling cannot validate survey/redshift covariance structure or calibrated uncertainty in a.'
}
(A/'covariance_scale_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
