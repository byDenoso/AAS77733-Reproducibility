"""Figure 6 from recomputed survey-intercept simulations and summaries."""
import json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; F=ROOT/'figures'; A=ROOT/'analysis'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
b=np.load(A/'baseline_results.npz'); n9=np.load(A/'survey9_results.npz'); n20=np.load(A/'survey20_results.npz')
sb=json.loads((A/'baseline_summary.json').read_text()); s9=json.loads((A/'survey9_summary.json').read_text()); s20=json.loads((A/'survey20_summary.json').read_text())
blue='#326788'; gray='#667784'; brown='#a9663d'; red='#9a3838'
fig,(ax,bx)=plt.subplots(1,2,figsize=(7.1,2.9),gridspec_kw={'wspace':.38,'width_ratios':[1.3,1]})
x=b['pivots']; ax.plot(x,b['observed_step'],c=blue,lw=1.3,label='Common intercept'); ax.plot(x,n9['observed_step'],c=gray,lw=1.3,ls='--',label='9 survey groups'); ax.plot(x,n20['observed_step'],c=brown,lw=1.3,ls='-.',label='20 survey offsets'); ax.axvspan(.34,.45,color=gray,alpha=.1)
ax.set_xlim(.025,.645); ax.set_ylim(0,max(5.6,float(max(b['observed_step'].max(),n9['observed_step'].max(),n20['observed_step'].max()))*1.35)); ax.set_xlabel(r'Candidate split $z_t$'); ax.set_ylabel(r'Observed improvement $\Delta\chi^2$'); ax.legend(fontsize=6.6,frameon=False,loc='upper right')
summary=[sb,s9,s20]; ys=[2,1,0]; colors=[blue,gray,brown]; labels=['One intercept','9 survey groups','20 survey offsets']
ps=[x['bin_family']['p'] for x in summary]; ses=[x['bin_family']['se'] for x in summary]
for i,(xval,se,c) in enumerate(zip(ps,ses,colors)):
    bx.errorbar(xval,ys[i],xerr=se,fmt='o',c=c,ms=5,capsize=2,lw=1.2); bx.text(xval+.003,ys[i]+.11,f'{xval:.3f}',fontsize=7,c=c)
bx.axvline(.05,c=red,lw=1,ls='--'); bx.set_xlim(0,max(.09,max(ps)+.02)); bx.set_ylim(-.5,2.65); bx.set_yticks(ys,labels); bx.set_xlabel('Twelve-bin family probability'); bx.grid(axis='x',lw=.3,alpha=.2)
fig.subplots_adjust(left=.09,right=.97,bottom=.22,top=.95); fig.savefig(F/'fig4_survey_robustness.pdf',bbox_inches='tight',pad_inches=.08); fig.savefig(F/'fig4_survey_robustness.png',dpi=200,bbox_inches='tight',pad_inches=.08); plt.close(fig)
np.savetxt(F/'fig4_survey_results.csv',np.column_stack([b['pivots'],b['observed_step'],n9['observed_step'],n20['observed_step']]),delimiter=',',header='candidate_pivot,common_intercept_delta_chi2,nine_survey_groups_delta_chi2,twenty_survey_offsets_delta_chi2',comments='')
print('FIG6',F/'fig4_survey_robustness.pdf')
