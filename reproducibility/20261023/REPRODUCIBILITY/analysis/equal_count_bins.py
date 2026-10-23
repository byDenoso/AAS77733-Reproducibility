# Derived from the R1 verified full-covariance equal-count bin method; no new global hypothesis test.
"""Descriptive Pantheon+ covariance-aware residual-bin figure.
Reuses the previously fitted background (not a new search or null test).
Outputs bin data and PDF; error bars are conditional on the fitted background.
"""
import os, io, json, hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pathlib import Path
P=Path(__file__).resolve().parents[1]
DAT=Path(os.environ.get('PANTHEON_DAT', str(P/'inputs'/'Pantheon+SH0ES.dat')))
COV=Path(os.environ.get('PANTHEON_COV', str(P/'inputs'/'Pantheon+SH0ES_STAT+SYS.cov')))
db=DAT.read_bytes()
cb=COV.read_bytes()
assert hashlib.sha256(db).hexdigest()=='1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8'
assert hashlib.sha256(cb).hexdigest()=='abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc'
d=np.genfromtxt(io.BytesIO(db),names=True,dtype=None,encoding='utf8')
flat=np.fromstring(cb.decode(),sep=' ')
n=int(flat[0]); C=flat[1:].reshape(n,n); C=(C+C.T)/2
idx=np.flatnonzero((d['IS_CALIBRATOR']==0)&(d['zHD']>0.0233));assert len(idx)==1365
z=d['zHD'][idx]; zh=d['zHEL'][idx];y=d['m_b_corr'][idx];C=C[np.ix_(idx,idx)]
Omega=0.3276969243165704
g=np.linspace(0,max(z)*1.001,8192)
integ=cumulative_trapezoid(1/np.sqrt(Omega*(1+g)**3+1-Omega),g,initial=0)
mu=5*np.log10((299792.458/70)*(1+zh)*np.interp(z,g,integ))+25
r=y-mu
L=cholesky(C,lower=True)
yw=solve_triangular(L,r,lower=True)
a=solve_triangular(L,np.ones(len(z)),lower=True)
Mhat=(a@yw)/(a@a)
resid=r-Mhat
# Equal-number redshift bins; GLS profile of all bin levels jointly with C.
q=np.quantile(z,np.linspace(0,1,13))
assert (np.diff(q)>0).all()
bi=np.searchsorted(q[1:-1],z,side='right')
B=(bi[:,None]==np.arange(12)[None,:]).astype(float)
W=solve_triangular(L,B,lower=True)
means=np.linalg.solve(W.T@W,W.T@solve_triangular(L,resid,lower=True))
covbin=np.linalg.inv(W.T@W)
errs=np.sqrt(np.diag(covbin))
median=np.array([np.median(z[bi==k]) for k in range(12)])
count=np.array([(bi==k).sum() for k in range(12)])
assert sum(count)==1365
records=[{'bin':k+1,'z_min':float(q[k]),'z_max':float(q[k+1]),'z_median':float(median[k]),'n_light_curves':int(count[k]),'gls_residual_mag':float(means[k]),'sigma_conditional_mag':float(errs[k])} for k in range(12)]
(P/'figures').mkdir(exist_ok=True)
(P/'analysis'/'equal_count_bin_summary.json').write_text(json.dumps({'description':'Full-covariance GLS descriptive redshift bin offsets relative to a fitted flat-LambdaCDM background; uncertainty conditional on fixed Omega_m and M','Omega_m':Omega,'Mhat':float(Mhat),'sample':1365,'binning':'12 equal-number redshift bins','bin':records},indent=2))
with open(P/'figures'/'fig7_equal_count_values.csv','w') as out:
 out.write('bin,z_min,z_max,z_median,n_light_curves,gls_residual_mag,sigma_conditional_mag\n')
 for x in records: out.write(','.join(str(x[k]) for k in ['bin','z_min','z_max','z_median','n_light_curves','gls_residual_mag','sigma_conditional_mag'])+'\n')
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.labelsize':10,'axes.titlesize':10,'pdf.fonttype':42,'ps.fonttype':42,'savefig.pad_inches':0.04})
fig,ax=plt.subplots(figsize=(6.8,2.65))
ax.axhline(0,color='#686868',lw=.9)
ax.axvline(.36,color='#b66b37',lw=1,linestyle=(0,(5,3)),label=r'Historical candidate $z_t=0.36$')
ax.errorbar(median,means,yerr=errs,fmt='o',color='#24547e',ecolor='#24547e',elinewidth=1.25,capsize=3,markersize=5.0,label='Binned GLS residuals')
ax.set_xscale('log')
ax.set_xlim(.021,max(2.05, max(z)*1.1))
ax.set_xticks([.03,.05,.1,.2,.36,.6,1.,2.])
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xlabel(r'Hubble-diagram redshift $z_{\mathrm{HD}}$ (log scale)')
ax.set_ylabel('Residual bin offset (mag)')
ax.grid(axis='y',color='#d2d6da',lw=.5,alpha=.75)
ax.legend(loc='upper right',fontsize=8,frameon=False)
ax.set_axisbelow(True)
fig.tight_layout(pad=.4)
fig.savefig(P/'figures'/'fig7_equal_count_bins.pdf',bbox_inches='tight')
fig.savefig(P/'figures'/'fig7_equal_count_bins.png',dpi=190,bbox_inches='tight')
print(json.dumps({'range_z':[float(min(z)),float(max(z))],'fit_intercept':float(Mhat),'bins':records},indent=2))
