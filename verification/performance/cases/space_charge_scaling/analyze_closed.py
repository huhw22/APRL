#!/usr/bin/env python3
"""Independently check tiny mean gains and interpret a completed charge scan."""
import csv,hashlib,json,sys
from decimal import Decimal,localcontext
from pathlib import Path
import h5py,numpy as np
from scipy.constants import m_e,c,e
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_csv(path,rows):
    with Path(path).open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def records(path):
    with h5py.File(path) as h:
        g=h['particle_plane'];assert int(g['complete'][()])==1
        a=g['records'][:];assert len(a)==int(g['committed_records'][()])==32768
    return a[np.argsort(a['source_id'])]

def decimal_mean_delta(initial,final,weights):
    # Direct LAB four-energy from stored LAB momenta; no inverse/forward boost.
    with localcontext() as ctx:
        ctx.prec=70;D=Decimal.from_float
        mc2=D(float(m_e))*D(float(c))**2/D(float(e));num=Decimal(0);den=Decimal(0)
        for ua,ub,w in zip(initial,final,weights):
            aa=sum((D(float(v))**2 for v in ua),Decimal(1))
            bb=sum((D(float(v))**2 for v in ub),Decimal(1))
            dg=(bb-aa)/(bb.sqrt()+aa.sqrt());q=D(float(abs(w)))
            num+=q*mc2*dg;den+=q
        return float(num/den)

def main(root):
    root=Path(root).resolve();result=json.loads((root/'review/RESULT.json').read_text())
    assert result['all_cases_closed'] and not result['missing']
    frozen={}
    for name in ['FROZEN.json','POSTPROCESS_FROZEN.json','THEORY_SOURCE.json']:
        for p,h in json.loads((root/name).read_text()).items():
            assert digest(p)==h,p;frozen[p]=h
    out=root/'analysis_final';out.mkdir(exist_ok=True)
    last=result['end_plane'];summary=[];precision=[];integrity=[]
    for row in last:
        d=Path(row['source_directory']);a=records(d/'detectors/entry.h5');b=records(d/'detectors/probe-5.h5')
        assert np.array_equal(a['source_id'],b['source_id'])
        exact=decimal_mean_delta(a['proper_velocity'],b['proper_velocity'],a['charge_C'])
        difference=exact-row['PIC_mean_lab_gain_eV'];assert abs(difference)<1e-5,(d,difference)
        precision.append(dict(charge_factor=row['charge_factor'],PIC_mean_lab_gain_decimal70_eV=exact,
                              prior_inverse_boost_longdouble_eV=row['PIC_mean_lab_gain_eV'],difference_eV=difference,
                              bias_eV=row['signed_fixed_plane_bias_eV']))
        original=json.loads((d/'analysis/RESULT.json').read_text());norm=original['relative_energy_budget']
        gain=original['planes'][-1]
        item=dict(charge_factor=row['charge_factor'],theory_gain_eV=row['theory_fixed_plane_mean_lab_gain_eV'],
                  PIC_gain_eV=row['PIC_mean_lab_gain_eV'],signed_bias_eV=row['signed_fixed_plane_bias_eV'],
                  signed_bias_pct=row['signed_fixed_plane_bias_pct'],bias_div_charge_factor_eV=row['signed_fixed_plane_bias_eV']/row['charge_factor'],
                  theory_div_charge_factor_squared_eV=row['theory_fixed_plane_mean_lab_gain_eV']/row['charge_factor']**2,
                  gain_std_relative_error_pct=100*(row['PIC_std_lab_gain_eV']/row['theory_std_lab_gain_eV']-1),
                  individual_gain_error_rms_div_theory_std_pct=100*row['individual_gain_error_rms_eV']/row['theory_std_lab_gain_eV'],
                  rest_gain_std_relative_error_pct=100*(gain['PIC_std_rest_gain_eV']/gain['reference_std_rest_gain_eV']-1),
                  maximum_ledger_gap_pct_fixed_D=row['maximum_ledger_gap_pct_fixed_D'],
                  ledger_fixed_D_J=norm['fixed_denominator_J'],ledger_final_gap_J=original['new_total_final_residual_J'],
                  ledger_maximum_absolute_gap_J=norm['fixed_denominator_J']*norm['maximum_absolute_total_residual_pct_fixed_scale']/100,
                  field_operator_maximum_residual_J=original['new_field_operator_maximum_residual_J'],
                  particle_grid_final_gap_J=original['particle_grid_exchange_final_J'])
        summary.append(item)
        with h5py.File(d/'ledger/energy-ledger.h5') as h:
            g=h['energy_ledger'];assert int(g['complete'][()])==1
            x=g['records'][:];assert len(x)==int(g['committed_records'][()])==1339
            assert all(np.isfinite(x[n]).all() for n in x.dtype.names if x.dtype[n].kind=='f')
        integrity.append(dict(charge_factor=row['charge_factor'],directory=str(d),ledger_records=len(x),
                              unique_particle_IDs=len(np.unique(a['source_id'])),last_step=int(x['step'][-1])))
    with (root/'review/all_planes.csv').open() as h:allrows=list(csv.DictReader(h))
    local=[]
    for plane in ['probe-1','probe-2','probe-3','probe-4','probe-5']:
        a=next(r for r in allrows if r['plane']==plane and float(r['charge_factor'])==.001)
        b=next(r for r in allrows if r['plane']==plane and float(r['charge_factor'])==.01)
        ea,eb=float(a['signed_fixed_plane_bias_eV']),float(b['signed_fixed_plane_bias_eV'])
        assert ea>0 and eb>0
        local.append(dict(plane=plane,z_m=float(a['z_m']),bias_q0p001_eV=ea,bias_q0p01_eV=eb,
                          local_bias_Q_exponent=float(np.log(eb/ea)/np.log(10)),
                          local_signal_Q_exponent=float(np.log(float(b['theory_fixed_plane_mean_lab_gain_eV'])/float(a['theory_fixed_plane_mean_lab_gain_eV']))/np.log(10))))
    write_csv(out/'precision_check.csv',precision);write_csv(out/'summary.csv',summary);write_csv(out/'weak_exponents_by_plane.csv',local)
    x=np.array([r['charge_factor'] for r in summary]);fig,ax=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    ax[0,0].plot(x,[r['signed_bias_eV'] for r in summary],'o-',color='#a64660')
    ax[0,0].set(xscale='log',yscale='symlog',ylabel='Signed mean bias (eV/electron)',title='Bias changes sign at stronger charge')
    ax[0,0].set_ylim(-1000,30)
    ax[0,0].axhline(0,color='0.5',lw=.7)
    ax[0,1].plot(x,[r['bias_div_charge_factor_eV'] for r in summary],'o-',color='#4b78a4')
    ax[0,1].set(xscale='log',ylabel='Mean bias / charge factor (eV)',title='Weak-charge bias is approximately linear')
    ax[0,1].axhline(0,color='0.5',lw=.7)
    for q in [.001,.01,.1,1.]:
        a=[r for r in allrows if float(r['charge_factor'])==q]
        ax[1,0].plot([float(r['z_m']) for r in a],[float(r['signed_fixed_plane_bias_eV'])/q for r in a],'.-',label=f'Q/Q0={q:g}')
    ax[1,0].set(xlabel='Lab plane z (m)',ylabel='Mean bias / charge factor (eV)',title='The trend is present at all five planes')
    ax[1,0].legend(fontsize=8)
    ax[1,1].semilogx(x,[r['gain_std_relative_error_pct'] for r in summary],'o-',label='Lab gain std')
    ax[1,1].semilogx(x,[r['rest_gain_std_relative_error_pct'] for r in summary],'s--',label='Rest-frame gain std')
    ax[1,1].axhline(0,color='0.5',lw=.7)
    ax[1,1].set(ylabel='PIC / theory std - 1 (%)',title='Gain-spread differences persist at weak charge',ylim=(-2,14))
    ax[1,1].legend(fontsize=8)
    for a in ax.ravel():a.grid(alpha=.22)
    for a in [ax[0,0],ax[0,1],ax[1,1]]:a.set_xlabel('Charge / original charge')
    fig.savefig(out/'scaling_diagnostics.png',dpi=180);fig.savefig(out/'scaling_diagnostics.pdf');plt.close(fig)
    save(out/'VALIDATION.json',dict(checks_pass=True,unique_frozen_files=len(frozen),closed_case_integrity=integrity,
         decimal70_maximum_mean_difference_eV=max(abs(r['difference_eV']) for r in precision),precision=precision,
         scope='Independent direct LAB high-precision energy reconstruction; not a new convergence test or physical acceptance.'))
    save(out/'CONCLUSIONS.json',dict(summary=summary,weak_exponents_by_plane=local,
         conclusion='Absolute bias decreases with charge; no fixed keV floor in this measured range. Fractional mean-gain bias is not constant. At the weakest pair the bias is approximately Q^0.94 at the last plane, versus signal Q^2.',
         caveats=['Do not fit a single global power law across the sign change.','The lowest mean gain is a cancellation of much larger signed individual changes.','Charge scaling at fixed spherical geometry is not an error bound for changing the aspect ratio.','Ledger gaps do not vanish merely because the field operator audit closes.']))
    print(json.dumps({'precision':precision,'summary':summary,'weak_exponents':local},indent=2))
if __name__=='__main__':main(sys.argv[1])
