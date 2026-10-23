"""Full-covariance Pantheon+ tests with profiled survey zero points and bin trials.

Reuses verified primary R4 arrays, and reads the original public table.
Every simulated null repeats the background refit, using a joint Gaussian
representation of its exact relevant covariance-whitened score directions.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from pathlib import Path
from scipy.linalg import cholesky,solve_triangular
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
a=np.load(ROOT/'analysis'/'result_arrays.npz')
d=np.genfromtxt(Path(os.environ.get('PANTHEON_DATA_DIR',ROOT/'inputs'))/'Pantheon+SH0ES.dat',names=True,dtype=None,encoding='utf-8')
ids=np.asarray(d['IDSURVEY'][a['selected_idx']],dtype=int)
n=len(ids);z=a['z'];piv=a['pivots'];oms=a['om_grid'];U0=a['U'];B0=a['B']
L=cholesky(a['C'],lower=True,check_finite=False)
w=lambda x:solve_triangular(L,x,lower=True,check_finite=False)
one=w(np.ones(n));one/=np.linalg.norm(one)
def demean(x):
 if x.ndim==1:return x-one*(one@x)
 return x-one[:,None]*(one@x)[None,:]
rb=demean(w(a['residual']))
S0=a['Q']*a['Cholesky_norm'][None,:]
edges=np.array([.0233,.05,.1,.2,.3,.4,.5,.65,.8,1,1.3,1.8,2.5])
bi=np.digitize(z,edges,right=True)-1
Z0=demean(w(np.column_stack([bi==j for j in range(len(edges)-1)]).astype(float)))

def summary(obs,sims):
 count=int((sims>=obs).sum());p=(count+1)/(len(sims)+1);phat=count/len(sims)
 return dict(T=float(obs),exceedances=count,p=float(p),se=float(np.sqrt(phat*(1-phat)/len(sims))),q95=float(np.quantile(sims,.95)))

def analyze(label,X,seed,mocks=10000):
 if X.shape[1]:
  WX=demean(w(X));E,Rx=np.linalg.qr(WX,mode='reduced')
  if np.min(np.abs(np.diag(Rx)))<1e-8:raise RuntimeError('rank deficient nuisance model')
 else:E=np.zeros((n,0))
 proj=lambda v:v-E@(E.T@v)
 Um,R=np.linalg.qr(proj(U0),mode='reduced');B=R@B0
 br=Um.T@proj(rb)
 cs=CubicSpline(oms,B,axis=1)
 fit=minimize_scalar(lambda x:float(np.sum((cs(x)-br)**2)),bounds=(oms[0],oms[-1]),method='bounded',options={'xatol':1e-11})
 best=float(fit.x);fitres=proj(rb)-Um@cs(best)
 def components(W):
  Wp=proj(W);norm=np.linalg.norm(Wp,axis=0)
  if norm.min()<1e-6:raise RuntimeError('step/bin projected to zero')
  return Wp/norm[None,:],norm
 Qs,ns=components(S0);Qz,nz=components(Z0)
 obs_s=fitres@Qs;obs_z=fitres@Qz
 deltas_s=obs_s**2;deltas_z=obs_z**2
 Brel=B-cs(best)[:,None];D=np.sum(Brel*Brel,axis=0)
 cross_s=Brel.T@(Um.T@Qs);cross_z=Brel.T@(Um.T@Qz)
 K=np.column_stack((Um,Qs,Qz));G=K.T@K;G=(G+G.T)/2
 eig,V=np.linalg.eigh(G)
 if eig.min() < -1e-8:raise RuntimeError('negative covariance eigenvalue')
 F=V*np.sqrt(np.maximum(eig,0))
 rng=np.random.default_rng(seed);q=rng.standard_normal((mocks,K.shape[1]))@F.T
 ek=q[:,:len(B)];es=q[:,len(B):len(B)+len(piv)];ez=q[:,len(B)+len(piv):]
 cost=D[None,:]-2*ek@Brel
 j=np.argmin(cost,axis=1);jj=np.clip(j,1,len(oms)-2);ii=np.arange(mocks)
 cm,c0,cp=cost[ii,jj-1],cost[ii,jj],cost[ii,jj+1]
 frac=np.divide(.5*(cm-cp),cm-2*c0+cp,out=np.zeros(mocks),where=np.abs(cm-2*c0+cp)>1e-12)
 frac=np.clip(frac,-1,1);frac[(j==0)|(j==len(oms)-1)]=0
 omega_mock=oms[jj]+frac*(oms[1]-oms[0]);omega_mock[j==0]=oms[0];omega_mock[j==len(oms)-1]=oms[-1]
 null_s=(es-CubicSpline(oms,cross_s,axis=0)(omega_mock))**2
 null_z=(ez-CubicSpline(oms,cross_z,axis=0)(omega_mock))**2
 im=int(np.argmax(deltas_s));iz=int(np.argmax(deltas_z));hi=9
 result={
 'label':label,'n':n,'n_survey_ids':len(np.unique(ids)),'n_extra_intercepts':X.shape[1],
 'best_omega_m':best,'no_step_chi2':float(fitres@fitres),'seed':seed,'mocks':mocks,
 'step_search':{**summary(deltas_s.max(),null_s.max(axis=1)), 'pivot':float(piv[im]),'amplitude_mag':float(obs_s[im]/ns[im]),
 'local_chi2_at_036':float(deltas_s[np.argmin(abs(piv-.36))]),'local_chi2_at_043':float(deltas_s[np.argmin(abs(piv-.43))])},
 'bin_family':{**summary(deltas_z.max(),null_z.max(axis=1)),'max_bin':[float(edges[iz]),float(edges[iz+1])],
 'high_bin_local_chi2':float(deltas_z[hi]),'high_bin_amplitude_mag':float(obs_z[hi]/nz[hi]),
 'high_bin_conditional_sigma_mag':float(1/nz[hi]),
 'high_bin_pointwise_p':summary(deltas_z[hi],null_z[:,hi])['p'],
 'high_bin_family_p':summary(deltas_z[hi],null_z.max(axis=1))['p']},
 'checks':{'gram_eig_min':float(eig.min()),'model_survey_projection':float(np.max(np.abs(Um.T@E))) if E.shape[1] else 0,
 'mock_omega_min':float(min(omega_mock)),'mock_omega_max':float(max(omega_mock)),'boundary_count':int(((j==0)|(j==len(oms)-1)).sum())}}
 (ROOT/'analysis'/f'{label}_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 np.savez_compressed(ROOT/'analysis'/f'{label}_results.npz',observed_step=deltas_s,observed_step_amp=obs_s/ns,
  observed_bin=deltas_z,observed_bin_amp=obs_z/nz,null_step=null_s.max(axis=1),null_bin=null_z.max(axis=1),
  bin_edges=edges,pivots=piv,survey_ids=np.unique(ids),omega_mock=omega_mock)
 print(json.dumps(result),flush=True)
 return result

if __name__=='__main__':
 vals,ct=np.unique(ids,return_counts=True)
 allx=np.column_stack([ids==v for v in vals[1:]]).astype(float)
 a=analyze('survey20',allx,2026102201)
 b=analyze('baseline',np.zeros((n,0)),2026102202)
 main=vals[np.argsort(ct)[-8:]];group=np.array([int(x) if x in main else -1 for x in ids])
 grp=np.unique(group);gx=np.column_stack([group==g for g in grp[1:]]).astype(float)
 c=analyze('survey9',gx,2026102203)
 print('FINISHED',json.dumps({'full_step':a['step_search'],'pooled_step':c['step_search'],'baseline_bins':b['bin_family'],'full_bins':a['bin_family']}),flush=True)
