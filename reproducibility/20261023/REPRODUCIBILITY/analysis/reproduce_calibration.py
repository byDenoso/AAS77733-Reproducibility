import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
import numpy as np, scipy, json,zipfile,io,hashlib,time
from scipy.linalg import cholesky,solve_triangular
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from pathlib import Path
base=Path(os.environ.get("PANTHEON_DATA_DIR",Path(__file__).resolve().parents[1]/"inputs"))
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from check_inputs import validate
check=validate(base)
if not check['ok']: raise RuntimeError('Pantheon+ input identity verification failed: '+str(check['errors']))
db=(base/"Pantheon+SH0ES.dat").read_bytes();cb=(base/"Pantheon+SH0ES_STAT+SYS.cov").read_bytes()
d=np.genfromtxt(io.BytesIO(db),names=True,dtype=None,encoding="utf8")
flat=np.fromstring(cb.decode(),sep=" "); n=int(flat[0]); C=flat[1:].reshape(n,n); C=(C+C.T)/2
mask=(d["IS_CALIBRATOR"]==0)&(d["zHD"]>.0233);idx=np.flatnonzero(mask)
z=d["zHD"][idx];zh=d["zHEL"][idx]; y=d["m_b_corr"][idx];C=C[np.ix_(idx,idx)]
L=cholesky(C,lower=True); a=solve_triangular(L,np.ones(len(z)),lower=True);a=a/np.linalg.norm(a)
def white(v):
 w=solve_triangular(L,v,lower=True);return w-a*(a@w) if w.ndim==1 else w-a[:,None]*(a@w)[None,:]
yz=white(y)
def mu(om):
 g=np.linspace(0,max(z)*1.001,8192);v=1/np.sqrt(om*(1+g)**3+1-om);integ=cumulative_trapezoid(v,g,initial=0)
 return 5*np.log10((299792.458/70)*(1+zh)*np.interp(z,g,integ))+25
def fun(om):r=yz-white(mu(om));return r@r
fit=minimize_scalar(fun,bounds=(.15,.5),method="bounded",options={"xatol":1e-11})
best=fit.x; rb=yz-white(mu(best)); chi=rb@rb
ps=np.round(np.arange(.01,1.50001,.01),10);ps=np.array([p for p in ps if (z<=p).sum()>=100 and (z>p).sum()>=100])
S=white(np.column_stack([z>p for p in ps]).astype(float)); norms=np.linalg.norm(S,axis=0);Q=S/norms
score=rb@Q; T=score**2; print(json.dumps({"stage":"observed","N":len(z),"omega_m":best,"chi2":chi,"pivot_min":min(ps),"pivot_max":max(ps),"T":max(T),"best_z":ps[T.argmax()],"step_amp":score[T.argmax()]/norms[T.argmax()],"target036":T[np.argmin(abs(ps-.36))],"target_window_T":max(T[(ps>=.34)&(ps<=.45)]),"input_sha256":{"table":hashlib.sha256(db).hexdigest(),"covariance":hashlib.sha256(cb).hexdigest()}}),flush=True)
oms=np.linspace(.15,.5,701)
M=white(np.column_stack([mu(o) for o in oms]))-white(mu(best))[:,None]
U,sv,V=np.linalg.svd(M,full_matrices=False);rank=int((sv>1e-8).sum());U=U[:,:rank];B=U.T@M
D=np.sum(M*M,axis=0);cross=M.T@Q
nsim=10000;rng=np.random.default_rng(20261018);A=np.column_stack([U,Q])
# Joint Gaussian scores are sampled from their exact covariance, not a data-vector approximation.
G=A.T@A;ev,evec=np.linalg.eigh((G+G.T)/2);factor=evec*np.sqrt(np.clip(ev,0,None))
v=rng.standard_normal((nsim,len(G)))@factor.T
eB=v[:,:rank]@B;eQ=v[:,rank:]
cost=D[None,:]-2*eB
j=np.argmin(cost,axis=1);edge=(j==0)|(j==len(oms)-1);jj=np.clip(j,1,len(oms)-2);ar=np.arange(nsim)
cm=cost[ar,jj-1];c0=cost[ar,jj];cp=cost[ar,jj+1]
frac=np.clip(.5*(cm-cp)/(cm-2*c0+cp),-1,1);frac[edge]=0
oh=oms[jj]+frac*(oms[1]-oms[0])
sp=CubicSpline(oms,cross,axis=0)
refitscore=eQ-sp(oh)
Tm=np.max(refitscore**2,axis=1);fixed=np.max(eQ**2,axis=1)
def summ(v,t):return {"exceedances":int((v>=t).sum()),"p_plus_one":float((1+(v>=t).sum())/(1+len(v))),"mc_se":float(np.sqrt(np.mean(v>=t)*(1-np.mean(v>=t))/len(v))),"q95":float(np.quantile(v,.95))}
print(json.dumps({"stage":"null_calibration","seed":20261018,"mocks":nsim,"model_svd_rank":rank,"max_svd_reconstruction_error":float(np.max(abs(U@B-M))),"omega_grid_step":.0005,"omega_boundary_mocks":int(edge.sum()),"original_fixed_background":summ(fixed,max(T)),"refit_background_each_mock_then_fixed_scan":summ(Tm,max(T)),"omega_mock_range":[float(min(oh)),float(max(oh))]}),flush=True)
# Validate grid-plus-quadratic null refits using direct continuous minimization on 30 preselected simulations in reduced model coordinates.
bs=CubicSpline(oms,B,axis=1);ds=CubicSpline(oms,D)
errors=[]
for k in range(30):
 f=lambda o:float(ds(o)-2*v[k,:rank]@bs(o))
 opt=minimize_scalar(f,bounds=(.15,.5),method="bounded",options={"xatol":1e-11})
 sc=eQ[k]-sp(opt.x);errors.append((abs(opt.x-oh[k]),abs(max(sc*sc)-Tm[k])))
print(json.dumps({"stage":"validation","n":30,"max_omega_interp_error":max(x[0] for x in errors),"max_T_interp_error":max(x[1] for x in errors),"versions":{"numpy":np.__version__,"scipy":scipy.__version__}}),flush=True)


# Observed-amplitude injection calibrated with the refitted-background null.
def refit_scan_mean(amp,piv):
 ss=white((z>piv).astype(float));mean=amp*ss
 vb=v[:,:rank]+mean@U;vq=eQ+mean@Q
 cc=D[None,:]-2*(vb@B);j=np.argmin(cc,axis=1);jj=np.clip(j,1,len(oms)-2);ar=np.arange(nsim)
 cm=cc[ar,jj-1];c0=cc[ar,jj];cp=cc[ar,jj+1];frac=np.clip(.5*(cm-cp)/(cm-2*c0+cp),-1,1)
 oo=oms[jj]+frac*.0005; scores=vq-sp(oo);mx=np.max(scores*scores,axis=1);pbest=ps[np.argmax(scores*scores,axis=1)]
 return {"amplitude_mag":float(amp),"pivot":float(piv),"detection_fraction":float(np.mean(mx>=np.quantile(Tm,.95))),"localization_fraction_pm003":float(np.mean(abs(pbest-piv)<=.0300001)),"MC_N":nsim}
print(json.dumps({"stage":"refitted_injection","results":[refit_scan_mean(score[T.argmax()]/norms[T.argmax()],ps[T.argmax()]),refit_scan_mean(score[np.argmin(abs(ps-.43))]/norms[np.argmin(abs(ps-.43))],.43)]}),flush=True)
# Supplementary joint background-plus-step likelihood ratio; observed optimizer independent per split.
robs=rb@M;obscost=np.sum(M*M,axis=0)-2*robs;cs=CubicSpline(oms,obscost)
profiles=[]
for k in range(len(ps)):
 objective=lambda o:float(cs(o)-(score[k]-sp(o)[k])**2)
 oo=minimize_scalar(objective,bounds=(.15,.5),method="bounded",options={"xatol":1e-10})
 profiles.append((-oo.fun,oo.x))
jointobs=max(profiles);jointTs=[]
for start in range(0,nsim,100):
 c=cost[start:start+100];eq=eQ[start:start+100]
 vals=c[:,:,None]-(eq[:,None,:]-cross[None,:,:])**2
 jointTs.extend(np.min(c,axis=1)-np.min(vals,axis=(1,2)))
jointTs=np.array(jointTs)
print(json.dumps({"stage":"joint_background_alternative","observed_T":float(jointobs[0]),"best_z":float(ps[np.argmax([p[0] for p in profiles])]),"best_omega_m":float(jointobs[1]),"calibration":summ(jointTs,jointobs[0]),"grid_omega_step":.0005,"claim":"SUPPLEMENTARY_NEW_STATISTIC_NOT_ORIGINAL_FROZEN_SCAN"}),flush=True)


from pathlib import Path
np.savez_compressed(str(Path(__file__).with_name("result_arrays.npz")), z=z, zh=zh, y=y, pivots=ps, score=score, step_amp=score/norms, scan_T=T, null_T=Tm, null_fixed_T=fixed, joint_null_T=jointTs, rank=rank, om_grid=oms, om_hat=best, best_chi=chi, injection_noise=v, covariance_projected_variance=np.diag(C), residual=y-mu(best), Cholesky_norm=norms, U=U, B=B, Q=Q, D=D, cross=cross, C=C, selected_idx=idx)
print(json.dumps({"stage":"files","result_arrays":"analysis/result_arrays.npz","status":"computed"}),flush=True)
