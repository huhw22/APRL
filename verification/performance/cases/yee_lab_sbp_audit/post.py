import json,os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import require,save
R=Path(__file__).resolve().parents[1]
def read(p):return np.genfromtxt(p,delimiter=',',names=True,dtype=np.longdouble)
def main():
    require((R/'CALCULATION_COMPLETE.json').exists(),'Calculation incomplete');p=json.loads((R/'PLAN.json').read_text());d=R/'yee-oblique';require((d/'complete').exists(),'Case incomplete')
    a=read(d/'history.csv');s=read(d/'sbp-state.csv');f=read(d/'sbp-flux.csv');base=read(Path(p['baseline'])/'history.csv');n=len(a)-1
    require(np.array_equal(a['step'],np.arange(n+1)) and np.array_equal(s['step'],a['step']) and np.array_equal(f['step'],np.arange(n)),'Step coverage')
    require(all(np.isfinite(z[k]).all() for z in [a,s,f] for k in z.dtype.names),'Nonfinite data')
    meta=json.loads((d/'frame_parameters.json').read_text());G=np.longdouble(meta['gamma']);b=np.longdouble(meta['beta']);scale=a['lab_charge_J'][0];u0=s['U_sync_J'][0]
    cum=lambda x:np.r_[np.longdouble(0),np.cumsum(x,dtype=np.longdouble)]
    def stock(x):return x-x[0]
    maxrel=lambda x:float(np.max(abs(x))/abs(scale))
    scalar_err={k:float(np.max(abs(a[k]-base[k]))/(np.max(abs(base[k]))+np.finfo(np.longdouble).tiny)) for k in a.dtype.names if k not in ['lab_balance_pct','boost_balance_pct','divB_fixed_scale']}
    require(max(scalar_err.values())<1e-11,'Old propagation scalar history changed')
    require(maxrel(a['Uprime_J']-s['U_sync_J'])<1e-11 and maxrel(a['cPzprime_J']-s['cP_sync_J'])<1e-11,'Independent stock definitions mismatch')
    stockerr=s['U_sync_J']-s['H_cross_J']-s['B_time_variance_J'];require(maxrel(stockerr)<1e-11,'Time-variance identity')
    RE=stock(s['H_cross_J'])+cum(f['Q_H_J']);RP=stock(s['cP_sync_J'])+cum(f['Q_P_J']);gauss=cum(f['Gauss_E_average_J']+f['Gauss_B_J']);canonical=G*(RE+b*RP);canonical_with_gauss=canonical+G*b*gauss
    require(maxrel(RE)<1e-9 and maxrel(RP+gauss)<1e-9 and maxrel(canonical_with_gauss)<1e-9,'Canonical Yee budget fails')
    physical=stock(a['lab_charge_J'])+a['lab_outward_J'];pieces=dict(B_stock_time_variance=G*stock(s['B_time_variance_J']),energy_flux_spatial_pairing=G*cum(f['old_Q_S_J']-f['Q_H_J']),stress_spatial_pairing=G*b*cum(f['old_Q_cT_J']-f['Q_P_at_mean_J']),stress_time_pairing=G*b*cum(f['Q_P_at_mean_J']-f['Q_P_J']))
    direct_time=G*b*cum(f['stress_time_variance_direct_J']);require(maxrel(direct_time-pieces['stress_time_pairing'])<1e-11,'Independent stress time identity')
    prediction=canonical+sum(pieces.values());reconstruction_error=physical-prediction;require(maxrel(reconstruction_error)<1e-10,'Independent decomposition fails')
    require(maxrel(a['boost_outward_J']-cum(f['old_Q_S_J']))<1e-11 and maxrel(a['lab_outward_J']-G*cum(f['old_Q_S_J']+b*f['old_Q_cT_J']))<1e-11,'Old measured ports mismatch')
    omitted=max(pieces,key=lambda k:float(np.max(abs(pieces[k]))));negative=physical-(prediction-pieces[omitted]);require(maxrel(negative)>1e-6,'Omission negative control failed')
    channels={'original_direct_lab_residual':physical,'Yee_canonical_lab_residual':canonical,'Yee_canonical_with_Gauss_check':canonical_with_gauss,**pieces,'independent_sum':prediction,'reconstruction_error':reconstruction_error}
    data=np.column_stack([a['step'],a['time_box_s']]+[100*v/scale for v in channels.values()]);np.savetxt(d/'decomposition_percent.csv',data,delimiter=',',header='step,time_box_s,'+','.join(k+'_percent' for k in channels),comments='')
    np.savetxt(d/'decomposition_J.csv',np.column_stack([a['step'],a['time_box_s']]+list(channels.values())),delimiter=',',header='step,time_box_s,'+','.join(k+'_J' for k in channels),comments='')
    answer=dict(complete=True,steps=n,original_lab_final_percent=float(100*physical[-1]/scale),channels_final_percent={k:float(100*v[-1]/scale) for k,v in channels.items()},max_abs_canonical_lab_percent=float(100*maxrel(canonical)),max_abs_canonical_with_gauss_percent=float(100*maxrel(canonical_with_gauss)),max_reconstruction_error_percent=float(100*maxrel(reconstruction_error)),max_stock_identity_error_relative=maxrel(stockerr),same_propagation_scalar_relative_errors=scalar_err,negative_control=dict(omitted=omitted,max_error_percent=float(100*maxrel(negative)),detected=True),normalization=dict(denominator_J=float(scale),definition='fixed initial native time-centered lab time-translation charge, shared by every displayed channel'),production_archive_unchanged=True,physical_continuum_accuracy_accepted=False,Cowan_tested=False,full_19174_root_cause_resolved=False)
    save(R/'RESULT.json',answer)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained');t=a['time_box_s']*1e12
    axes[0].plot(t,100*physical/scale,label='Original direct balance');axes[0].plot(t,100*prediction/scale,'--',label='Independent reconstruction');axes[0].plot(t,100*canonical/scale,label='Matched Yee discrete balance')
    labels=['Magnetic stock time variance','Energy flux spatial pairing','Stress spatial pairing','Stress time pairing']
    for (k,v),label in zip(pieces.items(),labels):axes[1].plot(t,100*v/scale,label=label)
    for ax in axes:ax.set(xlabel='Boost time (ps)',ylabel='% of fixed initial lab field charge');ax.axhline(0,color='black',lw=.5);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_title('Original residual and independent reconstruction');axes[1].set_title('Independently measured definition differences');fig.savefig(R/'decomposition.png',dpi=180);fig.savefig(R/'decomposition.pdf');plt.close(fig)
    print(json.dumps(answer,indent=2),flush=True)
if __name__=='__main__':main()
