#!/usr/bin/env python3
"""Read closed charge scans; compare absolute and fractional physical errors."""
import csv,hashlib,json,sys
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main(root):
    root=Path(root).resolve();sys.path.insert(0,str(root/'scripts'))
    from analyze import read_closed,rest_u,kinetic_ev,mean_std
    from model import G_BEAM,B_BEAM,c,m_e,e
    plan=json.loads((root/'plan.json').read_text());out=root/'review';out.mkdir(exist_ok=True)
    for p,h in json.loads((root/'FROZEN.json').read_text()).items():assert digest(p)==h,p
    post=root/'POSTPROCESS_FROZEN.json'
    if post.exists():
        for p,h in json.loads(post.read_text()).items():assert digest(p)==h,p
    cases=[(1.,Path(plan['baseline']))]+[(v['factor'],Path(v['directory'])) for v in plan['cases'].values()]
    missing=[];rows=[];LD=np.longdouble;mc2=LD(m_e)*LD(c)**2/LD(e)
    def lab_delta(a,b):
        aa=np.sum(a*a,axis=1);bb=np.sum(b*b,axis=1)
        dg=(bb-aa)/(np.sqrt(1+bb)+np.sqrt(1+aa))
        return LD(G_BEAM)*mc2*(dg+LD(B_BEAM)*(b[:,2]-a[:,2]))
    for factor,d in cases:
        if factor!=1.:
            marker=d/'COMPLETION.json'
            if not marker.exists() or not json.loads(marker.read_text()).get('checks_pass'):
                missing.append(str(d));continue
        result_path=d/'analysis/RESULT.json'
        if not result_path.exists():missing.append(str(d));continue
        prior=json.loads(result_path.read_text());meta=json.loads((d/'case.json').read_text())
        expected=np.load(d/'reference_input/expected_events.npz')
        entry=read_closed(d/'detectors/entry.h5','particle_plane');entry=entry[np.argsort(entry['source_id'])]
        entry_u=rest_u(entry['proper_velocity']);ref_entry=expected['entry__rest_u'].astype(LD)
        w=np.abs(entry['charge_C']).astype(LD);w/=np.sum(w)
        assert len(entry)==32768 and np.array_equal(entry['source_id'],np.arange(1,32769))
        with h5py.File(d/'ledger/energy-ledger.h5') as h:assert int(h['energy_ledger/complete'][()])==1
        old_by_plane={r['plane']:r for r in prior['planes']}
        for plane in meta['planes'][1:]:
            name=plane['name'];a=read_closed(d/'detectors'/(name+'.h5'),'particle_plane');a=a[np.argsort(a['source_id'])]
            assert np.array_equal(a['source_id'],entry['source_id'])
            actual=lab_delta(entry_u,rest_u(a['proper_velocity']))
            ref=lab_delta(ref_entry,expected[name+'__rest_u'].astype(LD))
            am,astd=mean_std(actual,w);rm,rstd=mean_std(ref,w);bias=am-rm
            check=old_by_plane[name]
            assert abs(am-check['PIC_mean_lab_gain_eV'])<1e-3,(factor,name,am,check['PIC_mean_lab_gain_eV'])
            row=dict(charge_factor=factor,electrons=meta['model']['electrons'],charge_pC=meta['model']['electrons']*e*1e12,
                     plane=name,z_m=plane['z'],PIC_mean_lab_gain_eV=am,theory_fixed_plane_mean_lab_gain_eV=rm,
                     signed_fixed_plane_bias_eV=bias,absolute_fixed_plane_bias_eV=abs(bias),
                     signed_fixed_plane_bias_pct=100*bias/abs(rm) if rm else None,
                     theory_at_actual_event_mean_lab_gain_eV=check['reference_mean_lab_gain_eV'],
                     signed_actual_event_bias_eV=am-check['reference_mean_lab_gain_eV'],
                     PIC_std_lab_gain_eV=astd,theory_std_lab_gain_eV=rstd,
                     individual_gain_error_rms_eV=float(np.sqrt(np.sum(w*(actual-ref)**2))),
                     rms_x_relative_error_pct=100*(check['PIC_sigma_x_m']/check['reference_sigma_x_m']-1),
                     maximum_ledger_gap_pct_fixed_D=prior['relative_energy_budget']['maximum_absolute_total_residual_pct_fixed_scale'],
                     source_directory=str(d))
            rows.append(row)
    save(out/'STATUS.json',dict(all_cases_closed=not missing,missing=missing,accepted_physics=False))
    if not rows:raise RuntimeError('No closed science case')
    with (out/'all_planes.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    last=sorted([r for r in rows if r['plane']=='probe-5'],key=lambda r:r['charge_factor'])
    with (out/'end_plane.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=last[0]);writer.writeheader();writer.writerows(last)
    exponents=[]
    for a,b in zip(last,last[1:]):
        value=dict(lower_factor=a['charge_factor'],upper_factor=b['charge_factor'])
        for key,label in [('absolute_fixed_plane_bias_eV','absolute_bias'),('theory_fixed_plane_mean_lab_gain_eV','physical_signal')]:
            aa,bb=abs(a[key]),abs(b[key]);value[label+'_local_Q_exponent']=float(np.log(bb/aa)/np.log(b['charge_factor']/a['charge_factor'])) if aa>0 and bb>0 else None
        value['bias_changes_sign']=a['signed_fixed_plane_bias_eV']*b['signed_fixed_plane_bias_eV']<0
        exponents.append(value)
    save(out/'RESULT.json',dict(all_cases_closed=not missing,missing=missing,physics_accepted=False,
         primary_scope='Same fixed lab input and output planes; each theory uses its own exact crossing events. Cold initial rest radius/velocities fixed across charge family.',
         end_plane=last,local_exponents=exponents,
         limitations=['A local exponent can be unstable near cancellation or a sign change.','Charge scaling at fixed spherical geometry is not aspect-ratio scaling.','No extrapolated production or radiation error bound.']))
    x=np.array([r['charge_factor'] for r in last]);fig,ax=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    ax[0].loglog(x,[abs(r['theory_fixed_plane_mean_lab_gain_eV']) for r in last],'o-',label='Theory physical gain')
    ax[0].loglog(x,[abs(r['PIC_mean_lab_gain_eV']) for r in last],'s--',label='PIC gain magnitude')
    ax[0].set(ylabel='Mean gain magnitude (eV/electron)',title='Physical signal')
    ax[1].loglog(x,[r['absolute_fixed_plane_bias_eV'] for r in last],'o-',label='Fixed-plane bias')
    ax[1].loglog(x,[abs(r['signed_actual_event_bias_eV']) for r in last],'x--',label='Actual-event diagnostic')
    ax[1].set(ylabel='Absolute bias (eV/electron)',title='Absolute error')
    ax[2].loglog(x,[abs(r['signed_fixed_plane_bias_pct']) for r in last],'o-',color='#b34562')
    ax[2].set(ylabel='Absolute bias / theory gain (%)',title='Fractional error')
    for a in ax:a.set_xlabel('Charge / original charge');a.grid(alpha=.2,which='both')
    ax[0].legend(fontsize=8);ax[1].legend(fontsize=8)
    fig.suptitle('Cold sphere at fixed size, mesh and 1.3262 m plane'+(' — incomplete' if missing else ''))
    fig.savefig(out/'charge_scaling.png',dpi=180);fig.savefig(out/'charge_scaling.pdf');plt.close(fig)
    print(json.dumps({'closed':not missing,'end_plane':last,'local_exponents':exponents},indent=2))
    if missing:raise RuntimeError('Some cases incomplete; partial results retained')
if __name__=='__main__':main(sys.argv[1])
