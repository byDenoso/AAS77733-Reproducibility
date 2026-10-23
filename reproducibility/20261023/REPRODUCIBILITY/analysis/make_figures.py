"""Generate Figures 1--5 only from observed inputs and recomputed analysis artifacts."""
import os, json
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import cho_factor,cho_solve
from scipy.integrate import cumulative_trapezoid, quad
ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'analysis'; F=ROOT/'figures'; F.mkdir(exist_ok=True)
a=np.load(A/'result_arrays.npz')
obs=np.genfromtxt(ROOT/'inputs'/'Pantheon+SH0ES.dat',names=True,dtype=None,encoding='utf8')
records=[json.loads(line) for line in (A/'console.jsonl').read_text().splitlines() if line.strip().startswith('{')]
rec={x['stage']:x for x in records}
der=json.loads((A/'derived_summary.json').read_text())
idx=a['selected_idx']; z=a['z']; y=a['y']; C=a['C']; N=len(z)
sl=cho_factor(C,lower=True,check_finite=False)
ones=np.ones(N); Ci1=cho_solve(sl,ones,check_finite=False); resid=a['residual']; intercept=(Ci1@resid)/(Ci1@ones); r=resid-intercept
cgray='#424B5B'; blue='#28618A'; verm='#B0433A'; gray='#7F8C93'; brown='#9A6A3A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'savefig.bbox':'tight','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def save(name):
    plt.savefig(F/(name+'.pdf'),bbox_inches='tight',pad_inches=.07)
    plt.savefig(F/(name+'.png'),dpi=170,bbox_inches='tight',pad_inches=.07)
    plt.close()
def csv(name,header,arr,fmt='%.8g'):
    np.savetxt(F/name,np.asarray(arr),delimiter=',',header=header,comments='',fmt=fmt)

# Figure 1: observed residuals and covariance-weighted descriptive bins.
be=np.asarray(der['display_bins']['edges']); ix=np.digitize(z,be,right=True)-1
means=np.asarray(der['display_bins']['means_mag']); errs=np.asarray(der['display_bins']['sigma_mag']); counts=np.asarray(der['display_bins']['counts']); mids=np.sqrt(be[:-1]*be[1:])
fig,(ax,bx)=plt.subplots(2,1,figsize=(7.0,4.25),sharex=True,gridspec_kw={'height_ratios':[1.45,1],'hspace':.08})
ax.scatter(z,r,s=7,c=cgray,alpha=.24,rasterized=True,linewidths=0); ax.axhline(0,c='black',lw=.7); ax.axvspan(.34,.45,color=blue,alpha=.085); ax.set_ylim(-.77,.77); ax.set_ylabel('Hubble residual (mag)')
bx.errorbar(mids,means,yerr=errs,fmt='o',ms=3.2,c=blue,lw=1.0,capsize=2); bx.axhline(0,c='black',lw=.7); bx.axvspan(.34,.45,color=blue,alpha=.085)
bx.set_xscale('log'); bx.set_xlim(.023,2.55); bx.set_ylim(-.13,.37); bx.set_xlabel(r'Hubble-diagram redshift $z_{\rm HD}$'); bx.set_ylabel('GLS bin mean (mag)')
bx.set_xticks([.03,.05,.1,.2,.4,.8,1.5,2.3]); bx.set_xticklabels(['.03','.05','.1','.2','.4','.8','1.5','2.3'])
hi=9; bx.annotate(f'{counts[hi]} light curves',xy=(mids[hi],means[hi]),xytext=(.65,.30),fontsize=7,arrowprops={'arrowstyle':'-','color':gray,'lw':.7})
save('fig1_residuals')
csv('fig1_observed_rows.csv','zHD,residual_mag,display_diagonal_sigma_mag',np.column_stack([z,r,np.sqrt(np.diag(C))]))
csv('fig1_covariance_bin_means.csv','z_bin_low,z_bin_high,row_count,gls_mean_mag,gls_sigma_mag',np.column_stack([be[:-1],be[1:],counts,means,errs]))

# Figure 2: actual primary scan and matched null distribution.
ps=a['pivots']; ts=a['scan_T']; null=a['null_T']; qt=float(np.quantile(null,.95)); tobs=float(ts.max())
fig,(ax,bx)=plt.subplots(1,2,figsize=(7.0,2.9),gridspec_kw={'width_ratios':[1.1,1],'wspace':.30})
ax.plot(ps,ts,c=blue,lw=1.5,label=r'Observed $\Delta\chi^2(z_t)$'); ax.axhline(qt,c=verm,linestyle='--',lw=1.0,label='Global 95% threshold'); ax.axvspan(.34,.45,color=blue,alpha=.085); ax.plot([ps[np.argmax(ts)]],[tobs],'o',c=blue,ms=3.7)
ax.set_xlabel(r'Trial redshift $z_t$'); ax.set_ylabel(r'$\Delta\chi^2$'); ax.set_xlim(.02,.65); ax.set_ylim(0,max(8,qt*1.18)); ax.legend(loc='upper right',frameon=False)
bins=np.arange(0,max(24,float(null.max())+1),.55); bx.hist(null,bins=bins,density=True,color=gray,alpha=.78,edgecolor='white',linewidth=.25)
bx.axvline(tobs,c=blue,lw=1.5,label=f'Observed {tobs:.3f}'); bx.axvline(qt,c=verm,lw=1.1,ls='--',label=f'95th percentile {qt:.3f}')
bx.set_xlim(0,max(17,float(null.max())*1.05)); bx.set_xlabel(r'$T_{\max}$ under fitted-background null'); bx.set_ylabel('Probability density'); bx.legend(frameon=False,loc='upper right')
save('fig2_scan_and_null')
csv('fig2_scan.csv','candidate_split_redshift,observed_improvement_chi2,estimated_amplitude_mag',np.column_stack([ps,ts,a['step_amp']]))
csv('fig2_null_values.csv','matched_global_max_delta_chi2',null[:,None])

# Figure 3: injection/recovery fractions from the recomputed primary calibration.
inj=rec['refitted_injection']['results']; piv=np.array([x['pivot'] for x in inj]); amp=np.array([x['amplitude_mag'] for x in inj]); power=np.array([x['detection_fraction'] for x in inj]); loc=np.array([x['localization_fraction_pm003'] for x in inj]); nm=np.array([x['MC_N'] for x in inj])
labels=[f'$z_t={p:.2f}$\n$A={aa:+.4f}$ mag' for p,aa in zip(piv,amp)]; pos=np.arange(len(inj)); fig,ax=plt.subplots(figsize=(3.35,2.75)); w=.31
ax.bar(pos-w/2,100*power,width=w,color=blue,label='Threshold crossing'); ax.bar(pos+w/2,100*loc,width=w,color=gray,label=r'Location within $\pm0.03$')
for k,v in enumerate(power): ax.text(pos[k]-w/2,100*v+1.2,f'{100*v:.1f}%',ha='center',fontsize=7)
for k,v in enumerate(loc): ax.text(pos[k]+w/2,100*v+1.2,f'{100*v:.1f}%',ha='center',fontsize=7)
ax.set_xticks(pos,labels); ax.set_ylim(0,max(65,100*max(loc)*1.2)); ax.set_ylabel(f'Fraction of {int(nm[0]):,} injections (%)'); ax.legend(frameon=False,loc='upper right',fontsize=6.3)
save('fig3_power')
csv('fig3_injections.csv','pivot_redshift,amplitude_mag,detection_fraction,localization_fraction,n_mocks',np.column_stack([piv,amp,power,loc,nm]))

# Figure 4: exact flat-LambdaCDM response, with comparison omega taken from the joint fit output.
grid=np.linspace(.005,1.2,900)
def distmod(om):
    E=np.sqrt(om*(1+grid)**3+1-om); d=cumulative_trapezoid(1/E,grid,initial=0); d+=quad(lambda x:1/np.sqrt(om*(1+x)**3+1-om),0,.005,epsabs=1e-13)[0]
    return 5*np.log10((299792.458/70)*(1+grid)*d)+25
om0=der['geometry']['omega_baseline']; om1=der['geometry']['omega_joint']; step43=der['geometry']['step43_amplitude_mag']; anchor=der['geometry']['anchor_redshift']
ms=distmod(om0); ml=distmod(om1); shift=-(ml-ms); shift-=np.interp(anchor,grid,shift)
fig,ax=plt.subplots(figsize=(3.35,2.7)); ax.plot(grid,shift,c=blue,lw=1.8,label=rf'$\Omega_m: {om0:.4f}\to{om1:.4f}$'); ax.axvspan(.34,.45,color=blue,alpha=.09); ax.axhline(step43,c=verm,ls='--',lw=.9,label=rf'Fitted $z_t=0.43$: {step43:+.4f} mag'); ax.axhline(0,c=gray,lw=.7)
ax.set_xlim(.05,1.2); ax.set_ylim(min(-.068,float(shift.min())*1.08),.01); ax.set_xlabel('Redshift $z$'); ax.set_ylabel('Change in residual (mag)'); ax.legend(frameon=False,loc='lower left',fontsize=6.1)
save('fig5_background_response'); csv('fig5_flatlcdm_response.csv','redshift,delta_observed_residual_mag',np.column_stack([grid,shift]))

# Figure 5: observed survey-by-redshift structure.
sv=np.asarray(obs['IDSURVEY'][idx],dtype=int); vals,cts=np.unique(sv,return_counts=True); top=list(vals[np.argsort(cts)[-8:][::-1]])
name_map={1:'SDSS',15:'PS1MD',10:'DES',4:'SNLS',150:'FOUND',5:'CSP',57:'LOSS2',64:'CfA3K'}
# Preserve field-recognizable names when available; otherwise show survey ID.
names=[name_map.get(int(x),f'ID {int(x)}') for x in top]+['Other surveys']
edges=np.array([.0233,.10,.20,.30,.40,.50,.65,.8,1.0,1.3,2.5]); sbin=np.digitize(z,edges,right=True)-1
rows=[[np.sum((sv==num)&(sbin==i)) for i in range(len(edges)-1)] for num in top]; rows.append([np.sum((~np.isin(sv,top))&(sbin==i)) for i in range(len(edges)-1)]); cm=np.array(rows,dtype=int)
fig,ax=plt.subplots(figsize=(7.0,2.8)); im=ax.imshow(np.log1p(cm),aspect='auto',cmap='Blues',interpolation='nearest')
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        if cm[i,j]>0: ax.text(j,i,str(cm[i,j]),ha='center',va='center',fontsize=6,color=('white' if cm[i,j]>65 else cgray))
ax.set_yticks(np.arange(len(names)),names); ax.set_xticks(np.arange(len(edges)-1),[f'{edges[i]:g}-{edges[i+1]:g}' for i in range(len(edges)-1)],rotation=30,ha='right'); ax.set_xlabel('Hubble-diagram redshift interval'); ax.set_ylabel('Survey origin'); fig.colorbar(im,ax=ax,label=r'$\ln(1+N)$',fraction=.019,pad=.012)
save('fig6_survey_composition'); np.savetxt(F/'fig6_counts.csv',cm,delimiter=',',fmt='%d',header=','.join(f'{edges[i]:g}-{edges[i+1]:g}' for i in range(len(edges)-1)),comments='')
print('FIGURES',[(x.name,x.stat().st_size) for x in sorted(F.glob('fig*.pdf'))])
