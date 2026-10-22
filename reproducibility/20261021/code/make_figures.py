import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MPLCONFIGDIR']='/tmp/pantheon_step_mpl'
os.environ['XDG_CACHE_HOME']='/tmp/pantheon_step_cache'
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from independent_step import Engine,ROOT
R=json.load(open(ROOT/'RESULT.json'));e=Engine(all_partitions=True);o=R['observed']['baseline'];s=R['null_calibrations']['baseline'];f=ROOT/'figures';f.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.4,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':True,'axes.spines.right':True,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True})
blue='#29618A';orange='#BC6B2F';gray='#72797E'
fig,(ax,bx)=plt.subplots(2,1,figsize=(3.35,4.1),sharex=True,gridspec_kw={'hspace':.1})
ax.scatter(e.z,e.y-e.fit_mean,s=2.7,color=gray,alpha=.22,linewidths=0)
bins=o['simultaneous_bins'];cent=np.sqrt(e.edges[:-1]*e.edges[1:]);bx.errorbar(cent,bins['amplitudes'],yerr=bins['sigma'],fmt='o',ms=3.2,capsize=2,color=blue,elinewidth=.8)
for a in [ax,bx]:
 a.axhline(0,color=gray,ls='--',lw=.7);a.axvspan(.34,.45,color=orange,alpha=.1);a.axvspan(1,1.3,color=blue,alpha=.1);a.set_xscale('log');a.set_xlim(.023,2.4)
ax.set_ylabel('SN residual (mag)');ax.set_ylim(-.88,.88)
bx.set_ylabel('Simultaneous GLS bin (mag)');bx.set_ylim(-.17,.36);bx.set_xlabel(r'Hubble-diagram redshift $z_{\mathrm{HD}}$');bx.set_xticks([.03,.1,.4,1,2]);bx.set_xticklabels(['.03','.1','.4','1','2'])
bx.annotate('9 SNe',xy=(cent[9],bins['amplitudes'][9]),xytext=(.44,.30),arrowprops={'arrowstyle':'-','lw':.6},fontsize=7)
fig.subplots_adjust(left=.2,right=.97,bottom=.12,top=.98);fig.savefig(f/'fig0_observed_residuals.pdf',metadata={'CreationDate':None,'ModDate':None,'Creator':None});plt.close(fig)
fig,ax=plt.subplots(figsize=(3.35,2.6));zsort=np.argsort(e.z);anchor=np.argmin(abs(e.z-.08))
for om,c,ls in [(.315,gray,':'),(.3,blue,'--'),(.29,orange,'-')]:
 d=e.fit_mu-e.mu(om);d-=d[anchor];ax.plot(e.z[zsort],d[zsort],c=c,ls=ls,lw=1.1,label=r'$\Omega_m='+str(om)+'$')
ax.axhline(0,color=gray,lw=.6);ax.axvspan(.34,.45,color=orange,alpha=.1);ax.set(xlim=(.02,.65),ylim=(-.04,.01),xlabel=r'$z_{\mathrm{HD}}$',ylabel='Relative residual shift (mag)');ax.legend(frameon=False,loc='lower left');fig.subplots_adjust(left=.21,right=.97,bottom=.2,top=.97);fig.savefig(f/'fig1_background_projection.pdf',metadata={'CreationDate':None,'ModDate':None,'Creator':None});plt.close(fig)
fig,ax=plt.subplots(figsize=(3.35,2.8));score=e.fit_resid@e.all_Q
ax.plot(e.all_pivots,score**2,c=gray,lw=.65,alpha=.8,label='All 1,090 partitions')
ax.plot(o['scan']['pivots'],o['scan']['T'],'o-',ms=2,lw=.9,c=blue,label='Original 62 pivots')
ax.axhline(s['primary_T']['q95'],c=blue,ls='--',lw=.8,label='62-pivot 95% threshold')
ax.axhline(s['all_partitions_T']['q95'],c=gray,ls=':',lw=.9,label='All-partition 95% threshold')
ax.axvspan(.34,.45,color=orange,alpha=.12);ax.set(xlim=(.02,.66),ylim=(0,9.5),xlabel='Trial split $z_t$',ylabel=r'$\chi^2$ improvement');ax.legend(frameon=True,facecolor='white',framealpha=1,edgecolor='none',loc='upper right');fig.subplots_adjust(left=.17,right=.97,bottom=.19,top=.97);fig.savefig(f/'fig2_observed_full_scan.pdf',metadata={'CreationDate':None,'ModDate':None,'Creator':None});plt.close(fig)
fig,ax=plt.subplots(figsize=(3.35,2.65));inj=[x for x in R['injection_results'] if x['pivot']==.43];a=np.array([-x['amplitude'] for x in inj]);y=np.array([x['power'] for x in inj]);ci=np.array([x['power_binomial_95'] for x in inj]).T
ax.errorbar(a,y,yerr=np.vstack([y-ci[0],ci[1]-y]),fmt='o-',c=blue,capsize=2,lw=1,ms=3)
ax.axhline(.05,c=gray,lw=.7,ls='--');ax.set(xlim=(0,.108),ylim=(0,1.04),xlabel='Injected downward step magnitude (mag)',ylabel='Detection probability');ax.text(.005,.89,'$z_t=0.43$',fontsize=8);fig.subplots_adjust(left=.18,right=.97,bottom=.21,top=.97);fig.savefig(f/'fig3_detection_power.pdf',metadata={'CreationDate':None,'ModDate':None,'Creator':None});plt.close(fig)
fig,ax=plt.subplots(figsize=(3.35,2.65));h=R['diagnostics']['high_z_bin'];vals=[h['simultaneous_GLS'],h['one_bin_amplitude'],h['survey_one_bin_amplitude']];err=[h['simultaneous_sigma'],h['one_bin_sigma'],h['survey_one_bin_sigma']]
ax.errorbar(vals,[2,1,0],xerr=err,fmt='o',color=blue,capsize=3);ax.axvline(0,c=gray,lw=.7,ls='--');ax.set_yticks([2,1,0]);ax.set_yticklabels(['Simultaneous bins','One-bin / global M','One-bin / survey M']);ax.set(xlim=(-.02,.39),ylim=(-.6,2.6),xlabel='High-redshift residual offset (mag)');fig.subplots_adjust(left=.4,right=.97,bottom=.22,top=.97);fig.savefig(f/'fig4_highz_nuisance.pdf',metadata={'CreationDate':None,'ModDate':None,'Creator':None});plt.close(fig)
np.savez_compressed(f/'figure_data.npz',z=e.z,residual=e.y-e.fit_mean,bin_edges=e.edges,bin_amplitudes=np.array(bins['amplitudes']),bin_sigma=np.array(bins['sigma']),all_partition_z=e.all_pivots,all_partition_T=score**2,grid_z=np.array(o['scan']['pivots']),grid_T=np.array(o['scan']['T']),q95_grid=s['primary_T']['q95'],q95_all=s['all_partitions_T']['q95'])
