"""Independent Pantheon+ step audit. No archived scientific implementation imported.
GL64 distances; full-covariance whitening; QR nuisance projection; Chebyshev
interpolation of the entire whitened model (no SVD rank truncation); continuous
bounded minimization with explicit endpoints and a global-bracketing mesh.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import json, hashlib, time, argparse
import numpy as np
from numpy.polynomial.chebyshev import chebval,chebfit,chebpts1
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky,solve_triangular,qr
from scipy.optimize import minimize_scalar
from scipy.stats import beta,chi2,ks_2samp
ROOT=Path(__file__).resolve().parents[1]


def save_json(path,obj):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(path)

def interval(hits,n):
    return [0. if hits==0 else float(beta.ppf(.025,hits,n-hits+1)),1. if hits==n else float(beta.ppf(.975,hits+1,n-hits))]

def mc_summary(v,t):
    hits=int(np.count_nonzero(v>=t));n=len(v)
    return {'N':n,'observed':float(t),'hits':hits,'p_plus_one':(hits+1)/(n+1),'binomial_95':interval(hits,n),'q95':float(np.quantile(v,.95)),'mean':float(np.mean(v))}


def vector_minimize(fun,shape,grid_n=65,iterations=48):
    """Bounded global mesh plus continuous golden refinement, including endpoints.
    Any multiple mesh-detected minima cause an exception rather than silent selection.
    Arrays of independent 1D objectives are evaluated in parallel, not approximated.
    """
    grid=np.linspace(-1,1,grid_n)
    vals=np.stack([fun(np.full(shape,x)) for x in grid],axis=-1)
    if not np.isfinite(vals).all():raise FloatingPointError('Nonfinite minimization objective')
    local=(vals[...,1:-1]<vals[...,:-2])&(vals[...,1:-1]<vals[...,2:])
    nlocal=np.sum(local,axis=-1)
    if np.any(nlocal>1):raise RuntimeError('Multiple local minima detected; explicit refinement required')
    j=np.argmin(vals,axis=-1)
    a=grid[np.maximum(j-1,0)].copy();b=grid[np.minimum(j+1,grid_n-1)].copy()
    g=(np.sqrt(5)-1)/2;c=b-g*(b-a);d=a+g*(b-a);fc=fun(c);fd=fun(d)
    for _ in range(iterations):
        take=fc<fd
        a=np.where(take,a,c);b=np.where(take,d,b)
        c=b-g*(b-a);d=a+g*(b-a);fc=fun(c);fd=fun(d)
    x=(a+b)/2;f=fun(x)
    for ep in [-1.,1.]:
        fe=fun(np.full(shape,ep));use=fe<f;x=np.where(use,ep,x);f=np.where(use,fe,f)
    if np.any(f>np.min(vals,axis=-1)+1e-8):raise RuntimeError('Refined minimum worse than mesh')
    return x,f,{'boundary_count':int(np.count_nonzero(np.abs(x)>1-1e-8)),'max_mesh_local_minima':int(np.max(nlocal))}


class Engine:
    def __init__(self,cut=.0233,survey=False,bounds=(.15,.50),extra_diag=False,calibrators=False,all_partitions=False):
        self.cut=cut;self.survey=survey;self.bounds=bounds
        raw=np.genfromtxt(ROOT/'inputs/Pantheon+SH0ES.dat',names=True,dtype=None,encoding='utf8')
        flat=np.fromstring((ROOT/'inputs/Pantheon+SH0ES_STAT+SYS.cov').read_text(),sep=' ')
        n=int(flat[0]);C=flat[1:].reshape(n,n)
        self.raw_asym=float(np.max(np.abs(C-C.T)))
        if self.raw_asym>3.1e-8:raise ValueError('Unexpected native covariance asymmetry')
        C=(C+C.T)/2 # explicit rounding-level symmetrization, identical to R3 contract
        self.idx=np.flatnonzero((raw['zHD']>cut)&((raw['IS_CALIBRATOR']==0)|calibrators))
        self.data=raw[self.idx];self.z=self.data['zHD'];self.zh=self.data['zHEL'];self.y=self.data['m_b_corr'];self.n=len(self.y)
        self.C=C[np.ix_(self.idx,self.idx)].copy()
        if extra_diag:self.C+=np.diag(self.data['m_b_corr_err_DIAG']**2)
        self.L=cholesky(self.C,lower=True)
        self.groups=np.unique(self.data['IDSURVEY'])
        self.X=(self.data['IDSURVEY'][:,None]==self.groups).astype(float) if survey else np.ones((self.n,1))
        wx=self.white(self.X);self.U,self.R=qr(wx,mode='economic')
        if np.linalg.matrix_rank(wx)!=self.X.shape[1]:raise ValueError('Rank deficient nuisance design')
        self.reference_om=np.mean(bounds);self.reference_mu=self.mu(self.reference_om)
        self.yw=self.project(self.white(self.y-self.reference_mu))
        self.degree=48 if bounds[1]-bounds[0]>.5 else 24;nodes=chebpts1(self.degree+1)
        dm=np.column_stack([self.mu(self.om(x))-self.reference_mu for x in nodes]);wm=self.project(self.white(dm))
        self.coeff=chebfit(nodes,wm.T,self.degree).T # N x degree+1
        dn=chebpts1(2*self.degree+1);wv=chebval(dn,self.coeff.T)
        self.dcoeff=chebfit(dn,np.sum(wv*wv,axis=0),2*self.degree)
        piv=np.round(np.arange(.01,1.50001,.01),10)
        self.pivots=np.array([p for p in piv if np.count_nonzero(self.z<=p)>=100 and np.count_nonzero(self.z>p)>=100])
        self.edges=np.array([.0233,.05,.1,.2,.3,.4,.5,.65,.8,1.,1.3,2.5])
        steps=(self.z[:,None]>self.pivots).astype(float)
        bins=np.column_stack([(self.z>=a)&(self.z<b) for a,b in zip(self.edges[:-1],self.edges[1:])]).astype(float)
        self.nstep=len(self.pivots);self.nbins=len(self.edges)-1
        self.bin_operator=None
        if not survey and cut==.0233:
            bw=self.white(bins);bf=bw.T@bw
            if np.linalg.matrix_rank(bf)!=len(bf):raise ValueError("Unidentifiable simultaneous bins")
            op=np.linalg.solve(bf,bw.T)
            self.bin_operator=self.project(op.T).T
            self.bin_sigma=np.linalg.norm(self.bin_operator,axis=1)
            self.bin_model_coeff=self.bin_operator@self.coeff
        self.templates=np.column_stack([steps,bins]);self.setup_templates()
        self.all_pivots=None
        if all_partitions:
            self.all_pivots=np.array([p for p in np.unique(self.z) if np.count_nonzero(self.z<=p)>=100 and np.count_nonzero(self.z>p)>=100])
            wa=self.project(self.white((self.z[:,None]>self.all_pivots).astype(float)))
            self.all_norm=np.linalg.norm(wa,axis=0);self.all_Q=wa/self.all_norm
            self.all_cross=self.coeff.T@self.all_Q
        self.fit_x,self.fit_cost,_=self.fit_projected(self.yw[None,:],joint=False)
        self.fit_x=float(self.fit_x[0]);self.fit_om=float(self.om(self.fit_x));self.fit_cost=float(self.fit_cost[0])
        self.fit_mu=self.mu(self.fit_om);self.fit_beta=np.linalg.solve(self.R,self.U.T@self.white(self.y-self.fit_mu))
        self.fit_mean=self.fit_mu+self.X@self.fit_beta
        self.fit_resid=self.project(self.white(self.y-self.fit_mu));self.fit_chi2=float(self.fit_resid@self.fit_resid)
        self.interpolation_validation=self.validate_interpolation()

    def white(self,v):return solve_triangular(self.L,v,lower=True,check_finite=True)
    def project(self,v):return v-self.U@(self.U.T@v)
    def om(self,x):return self.bounds[0]+(np.asarray(x)+1)*(self.bounds[1]-self.bounds[0])/2
    def x(self,om):return 2*(np.asarray(om)-self.bounds[0])/(self.bounds[1]-self.bounds[0])-1
    def mu(self,om,order=64):
        q,w=leggauss(order);zz=self.z[:,None]*(q+1)/2
        integral=self.z/2*np.sum(w/np.sqrt(om*(1+zz)**3+1-om),axis=1)
        return 5*np.log10(299792.458/70*(1+self.zh)*integral)+25
    def setup_templates(self):
        ws=self.project(self.white(self.templates));self.norms=np.linalg.norm(ws,axis=0)
        self.valid=self.norms>1e-8
        if not self.valid[:self.nstep].all():raise ValueError('Unidentifiable step template')
        self.Q=np.zeros_like(ws);self.Q[:,self.valid]=ws[:,self.valid]/self.norms[self.valid]
        self.cross=self.coeff.T@self.Q
    def validate_interpolation(self):
        err=[];whiteerr=[];quad=[]
        for om in np.linspace(*self.bounds,41):
            direct=self.project(self.white(self.mu(om)-self.reference_mu));interp=chebval(self.x(om),self.coeff.T)
            whiteerr.append(np.max(np.abs(direct-interp)));quad.append(np.max(np.abs(self.mu(om)-self.mu(om,128))))
        out={'max_projected_whitened_model_error':float(max(whiteerr)),'max_GL64_vs128_mag':float(max(quad)),'nodes':self.degree+1,'validation_points':41,'model_rank_truncation':False}
        if max(whiteerr)>1e-7 or max(quad)>1e-7:raise ValueError('Interpolation/quadrature validation failed: '+repr(out))
        return out
    def coefficients_cost(self,yw):
        ec=yw@self.coeff;c=np.broadcast_to(self.dcoeff,(len(yw),len(self.dcoeff))).copy();c[:,:self.degree+1]-=2*ec
        return c
    def fit_projected(self,yw,joint=False):
        c=self.coefficients_cost(yw)
        if not joint:
            f=lambda x:chebval(x,c.T,tensor=False)
            return vector_minimize(f,(len(yw),))
        e=yw@self.Q[:,:self.nstep];cross=self.cross[:,:self.nstep]
        f=lambda x:chebval(x,c.T[:,:,None],tensor=False)-(e-chebval(x,cross[:,None,:],tensor=False))**2
        return vector_minimize(f,(len(yw),self.nstep))
    def evaluate(self,yw,joint=True,all_partitions=False):
        x,c,flags=self.fit_projected(yw)
        scores=yw@self.Q-chebval(x,self.cross).T
        primary=np.max(scores[:,:self.nstep]**2,axis=1)
        bt=scores[:,self.nstep:]**2;bt[:,~self.valid[self.nstep:]]=0
        out={'omega':self.om(x),'primary_T':primary,'primary_pivot':self.pivots[np.argmax(scores[:,:self.nstep]**2,axis=1)],'bin_T':np.max(bt,axis=1),'combined_T':np.maximum(primary,np.max(bt,axis=1)),'chi2':np.sum(yw*yw,axis=1)+c,'boundary':np.abs(x)>1-1e-8}
        if self.bin_operator is not None:
            bm=yw@self.bin_operator.T-chebval(x,self.bin_model_coeff.T).T
            out['simultaneous_bin_T']=np.max((bm/self.bin_sigma)**2,axis=1)
        if joint:
            xj,cj,jflags=self.fit_projected(yw,joint=True);j=np.argmin(cj,axis=1)
            out.update(joint_T=c-np.min(cj,axis=1),joint_omega=self.om(xj[np.arange(len(xj)),j]),joint_pivot=self.pivots[j],joint_boundary=np.abs(xj[np.arange(len(xj)),j])>1-1e-8)
        if all_partitions:
            sa=yw@self.all_Q-chebval(x,self.all_cross).T
            out.update(all_partitions_T=np.max(sa*sa,axis=1),all_partitions_pivot=self.all_pivots[np.argmax(sa*sa,axis=1)])
        return out
    def observed(self):
        out=self.evaluate(self.yw[None,:],all_partitions=self.all_pivots is not None)
        out={k:(bool(v[0]) if v.dtype==bool else float(v[0])) for k,v in out.items()}
        scores=self.fit_resid@self.Q;amps=np.divide(scores,self.norms,out=np.zeros_like(scores),where=self.valid)
        direct=minimize_scalar(lambda o:float(np.sum(self.project(self.white(self.y-self.mu(o)))**2)),bounds=self.bounds,method='bounded',options={'xatol':1e-12})
        if not direct.success:raise RuntimeError('Direct observed optimizer failed')
        derivative=self.project(self.white((self.mu(self.fit_om+1e-5)-self.mu(self.fit_om-1e-5))/2e-5))
        info=float(derivative@derivative)
        dof=self.n-self.X.shape[1]-1
        out.update(n=self.n,unique_SN=len(np.unique(self.data['CID'])),n_nuisance=self.X.shape[1],survey_groups=self.groups.tolist(),fit_intercepts=self.fit_beta.tolist(),fit_chi2=self.fit_chi2,nominal_dof=dof,nominal_chi2_upper_tail=float(chi2.sf(self.fit_chi2,dof)),nominal_chi2_lower_tail=float(chi2.cdf(self.fit_chi2,dof)),direct_optimizer={'omega':float(direct.x),'chi2':float(direct.fun),'success':bool(direct.success)},interpolation=self.interpolation_validation,norms=self.norms.tolist(),Omega_information_after_linear_nuisance=info,raw_covariance_asymmetry=self.raw_asym,cut=self.cut,bounds=list(self.bounds),survey_nuisance=self.survey,scan={'pivots':self.pivots.tolist(),'T':(scores[:self.nstep]**2).tolist(),'amplitudes':amps[:self.nstep].tolist()},bins={'edges':self.edges.tolist(),'counts':self.templates[:,self.nstep:].sum(axis=0).astype(int).tolist(),'T':(scores[self.nstep:]**2).tolist(),'amplitudes':amps[self.nstep:].tolist(),'sigma':np.divide(1.,self.norms[self.nstep:],out=np.zeros(self.nbins),where=self.valid[self.nstep:]).tolist(),'valid':self.valid[self.nstep:].tolist()})
        if self.bin_operator is not None:
            out['simultaneous_bins']={'amplitudes':(self.bin_operator@self.fit_resid).tolist(),'sigma':self.bin_sigma.tolist(),'covariance':(self.bin_operator@self.bin_operator.T).tolist(),'definition':'Same simultaneous-bin GLS estimator centered by fitted common intercept as R3 figure; covariance conditional on Omega; mock statistic refits Omega'}
        if self.all_pivots is not None:out['n_all_partitions']=len(self.all_pivots)
        return out


def run_mock(name,n,seed,cut=.0233,survey=False,batch=250,pilot=False):
    started=time.monotonic();e=Engine(cut=cut,survey=survey,all_partitions=(name=='baseline'))
    obs=e.observed();outdir=ROOT/'outputs'/name;outdir.mkdir(exist_ok=True)
    save_json(outdir/'observed.json',obs)
    rng=np.random.default_rng(seed);parts={};model_at_null=chebval(e.x(e.fit_om),e.coeff.T)
    noise_first=[]
    for start in range(0,n,batch):
        eps=rng.standard_normal((min(batch,n-start),e.n))
        if start<2000:noise_first.append(eps[:min(len(eps),2000-start)])
        yw=e.project(eps.T).T+model_at_null
        result=e.evaluate(yw,all_partitions=name=='baseline')
        for k,v in result.items():parts.setdefault(k,[]).append(v)
        arrays={k:np.concatenate(v) for k,v in parts.items()}
        np.savez_compressed(outdir/'draws.npz',**arrays,seed=seed)
        save_json(outdir/'PROGRESS.json',{'status':'RUNNING' if start+len(eps)<n else 'COMPLETED','done':start+len(eps),'target':n,'seed':seed,'elapsed_seconds':time.monotonic()-started})
        print(json.dumps({'contract':name,'done':start+len(eps),'target':n,'seconds':round(time.monotonic()-started,2)}),flush=True)
    summary={k:mc_summary(arrays[k],obs[k]) for k in ['primary_T','joint_T','bin_T','combined_T']+(['all_partitions_T'] if name=='baseline' else [])+(['simultaneous_bin_T'] if e.bin_operator is not None else [])}
    summary['Omega_boundary_count']=int(arrays['boundary'].sum());summary['joint_Omega_boundary_count']=int(arrays['joint_boundary'].sum())
    summary['gof_upper']=mc_summary(arrays['chi2'],obs['fit_chi2']);summary['gof_lower']=mc_summary(-arrays['chi2'],-obs['fit_chi2'])
    summary['seed']=seed;summary['runtime_seconds']=time.monotonic()-started
    save_json(outdir/'summary.json',summary)
    if name=='baseline':np.savez_compressed(outdir/'first2000_noise.npz',eps=np.concatenate(noise_first),seed=seed,indices=np.arange(min(n,2000)))
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['observed','pilot','baseline','controls','injection']);args=p.parse_args()
    if args.mode=='observed':
        results=json.load(open(ROOT/'outputs/observed_all.json')) if (ROOT/'outputs/observed_all.json').exists() else {}
        for name,kwargs in [('baseline',{'all_partitions':True}),('cut001',{'cut':.01}),('cut0023',{'cut':.023}),('cut0024',{'cut':.024}),('cut005',{'cut':.05}),('cut01',{'cut':.1}),('survey',{'survey':True}),('expanded',{'bounds':(.05,.95)}),('extra_diag',{'extra_diag':True}),('include_cal001',{'cut':.01,'calibrators':True}),('include_cal00233',{'calibrators':True})]:
            if name in results:continue
            e=Engine(**kwargs);results[name]=e.observed();print(name,results[name]['omega'],results[name]['primary_T'],results[name]['joint_T'],flush=True)
            save_json(ROOT/'outputs/observed_all.json',results)
    elif args.mode=='pilot':run_mock('pilot',40,2026102199,batch=40)
    elif args.mode=='baseline':run_mock('baseline',10000,2026102101)
    elif args.mode=='controls':
        for name,cut,seed in [('cut001',.01,2026102102),('cut0023',.023,2026102103),('cut0024',.024,2026102104),('cut01',.1,2026102105)]:run_mock(name,2000,seed,cut=cut)
        run_mock('survey',2000,2026102106,survey=True)
    elif args.mode=='injection':
        e=Engine();eps=np.load(ROOT/'outputs/baseline/first2000_noise.npz')['eps'];baseline=json.load(open(ROOT/'outputs/baseline/summary.json'));thr=baseline['primary_T']['q95'];out=[]
        for pivot,amp in [(.03,.02330572628987407),(.43,-.01777699228724288),(.43,-.02),(.43,-.05),(.43,-.10)]:
            mu=chebval(e.x(e.fit_om),e.coeff.T)+e.project(e.white(amp*(e.z>pivot)))
            yw=e.project(eps.T).T+mu;res=e.evaluate(yw,joint=False);hits=int((res['primary_T']>=thr).sum());near=int((np.abs(res['primary_pivot']-pivot)<=.0300001).sum())
            out.append({'pivot':pivot,'amplitude':amp,'N':len(eps),'power_hits':hits,'power':hits/len(eps),'power_binomial_95':interval(hits,len(eps)),'threshold':thr,'localization_hits_pm003':near,'localization_fraction':near/len(eps),'Omega_boundary_count':int(res['boundary'].sum())})
            np.savez_compressed(ROOT/f'outputs/injection_z{pivot}_A{amp}.npz',**res)
        save_json(ROOT/'outputs/injection.json',out)
if __name__=='__main__':main()
