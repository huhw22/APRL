"""Freeze one historical same-input six-segment Puffin control, no resampling."""
import math
import shutil
import sys
from pathlib import Path
import h5py
import numpy as np
import yaml
from scipy.fft import next_fast_len

W=Path('/public/home/wfollowm/APRL_Workspace')
P=W/'results/puffin_eleven_four_slip_20261113_r2'
sys.path.insert(0,str(P/'scripts'))
from common import C,E,read,save,need,sha
from puffin_units import scaling


def deck(root,directory,plan,scale,nz,preflight=False):
    periods=1 if preflight else plan['module_periods'];steps=80
    (directory/'beam.in').write_text("&NBLIST\nnbeams=1\ndtype='h5'\n/\n&BH5LIST\ndist_f='../particles.h5'\n/\n")
    (directory/'matched.latt').write_text(f"UN 'aprl-sparc11' {periods} 1.0 0.0 {steps} 0.0 1.0 0.0 0.0\n")
    (directory/'writes.in').write_text('WR 1\n'+('' if preflight else ''.join(f"WR {x['step']}\n" for x in plan['outputs'])))
    evolve='.false.' if preflight else '.true.'
    (directory/'main.in').write_text(f'''&MDATA
qScaled=.true.
qOneD=.false.
qFieldEvolve={evolve}
qWriteA={evolve}
qElectronsEvolve=.true.
qElectronFieldCoupling=.true.
qDiffraction=.true.
qFocussing=.false.
qUndEnds=.true.
qFMesh_G=.false.
q_noise=.false.
qFilter=.true.
qMeasure=.false.
meshType=0
beam_file='beam.in'
lattFile='matched.latt'
wr_file='writes.in'
iNumNodesX=81
iNumNodesY=81
nodesPerLambdar=21
sFModelLengthX={.0032/scale['transverse_scale_m']:.17g}
sFModelLengthY={.0032/scale['transverse_scale_m']:.17g}
sFModelLengthZ2={(nz-1-1e-7)*scale['lambda_reference_m']/20/scale['cooperation_length_m']:.17g}
iRedNodesX=49
iRedNodesY=49
sFiltFrac=.3
sDiffFrac=.25
sBeta=1
sElectronThreshold=0
srho=.005
saw={plan['K']:.17g}
sgamma_r={plan['gamma_ref']:.17g}
lambda_w=.028
sux=0
suy=1
zundType='aprl-sparc11'
stepsPerPeriod=80
nPeriods={periods}
sZ0=0
iWriteNthSteps={periods*steps+1}
iWriteIntNthSteps=80
sRedistLen=.5
iRedistStp=20
ioutInfo=2
/
''')


def main(root):
    root=Path(root).resolve();root.mkdir(exist_ok=False)
    template=Path(__file__).resolve().parent;shutil.copytree(template,root/'scripts')
    for name in ['common.py','binding.py','puffin_units.py']:
        shutil.copy2(P/'scripts'/name,root/'scripts'/name)
    for name in ['inputs','import-check','main','preparation','reference']: (root/name).mkdir()
    old=W/'results/sparc_six_full_ledger_20261112_r3';aprl=read(old/'PLAN.json');prior=read(P/'PLAN.json')
    source=Path(aprl['input_file']);need(sha(source)==aprl['input_sha256'],'Historical APRL input changed')
    binary=Path(prior['binary']);oldfreeze=read(P/'FROZEN_PHYSICS.json')
    need(sha(binary)==oldfreeze[str(binary)],'Validated Puffin binary changed')
    need(read(P/'FIELD_VALIDATED.json')['passed'] and read(P/'SMOKE_VALIDATED.json')['passed'],'Missing solver/magnet qualification')
    magnets=yaml.safe_load((old/'energy-binomial3/run.yaml').read_text())['sources']['magnetic_elements']
    need(magnets==prior['lattice']['magnets'][:11],'Six-segment magnetic lattice differs from qualified prefix')
    with h5py.File(source) as f:a=f['particles/records'][:]
    need(len(a)==aprl['input_particles'] and np.isfinite(a['macro_weight']).all(),'Input population/weights invalid')
    need(np.all(a['macro_weight']>0),'Nonpositive macro weight')
    n=aprl['input_represented_electrons'];weight=a['macro_weight']*(n/np.sum(a['macro_weight'],dtype=np.longdouble))
    weight=np.asarray(weight,dtype=float)
    # APRL's reader also uses electrons * input_weight / sum(input_weights).
    # This is the physical charge conversion, not a fitted field normalization.
    K=aprl['physics']['K'];gamma_ref=aprl['physics']['gamma'];scale=scaling(gamma_ref,K,.028)
    u=a['proper_velocity'];gamma=np.sqrt(1+np.sum(u*u,axis=1));drift=.01
    arrival=a['arrival_time_offset_s']+drift*gamma/(C*u[:,2])
    xy=a['plane_position_m']+drift*u[:,:2]/u[:,2,None]
    lam=aprl['physics']['wavelength_m'];anchor=float(arrival.min()-128*lam/C)
    native=np.column_stack((xy/scale['transverse_scale_m'],
        scale['code_c_m_per_s']*(arrival-anchor)/scale['cooperation_length_m'],
        u[:,0]/K,-u[:,1]/K,gamma/gamma_ref,weight*E/(scale['qe_C']*scale['npkbar'])))
    with h5py.File(root/'particles.h5','x') as f:
        f['electrons']=native;f['source_id']=a['source_id'];f['physical_electrons_per_record']=weight
        f.attrs.update(physical_time_origin_s=anchor,physical_entrance_m=.01,
            scope='Historical correlated input preserved solely for comparison with APRL jobs19743/19744')
    back=np.column_stack((native[:,:2]*scale['transverse_scale_m'],native[:,2]*scale['cooperation_length_m']/scale['code_c_m_per_s']+anchor,
        native[:,3]*K,-native[:,4]*K,native[:,5]*gamma_ref))
    expected=np.column_stack((xy,arrival,u[:,:2],gamma));error=np.max(abs(back-expected),axis=0)
    need(np.all(error<np.array([1e-15,1e-15,1e-24,1e-12,1e-12,1e-10])),'Event mapping roundtrip failed')
    charge=float(np.sum(native[:,6].astype(np.longdouble))*scale['npkbar']*scale['qe_C'])
    need(abs(charge/(n*E)-1)<2e-12,'Physical charge changed')
    span=float(np.ptp(arrival)*scale['code_c_m_per_s']/scale['lambda_reference_m'])
    nz=128*next_fast_len(math.ceil(((span+600+256)*20+1)/128))
    outputs=prior['outputs'][:6];periods=math.ceil((outputs[-1]['actual_lab_z_m']-.01)/.028)
    plan=dict(root=str(root),binary=str(binary),binary_sha256=sha(binary),validated_Puffin_parent=str(P),
        APRL_parent=str(old),APRL_jobs=[19743,19744],input_file=str(source),input_sha256=aprl['input_sha256'],
        sample_records=len(a),physical_electrons=n,source_charge_C=n*E,duration_s=aprl['physics']['duration_s'],
        segments=6,periods_per_segment=77,module_periods=periods,steps_per_period=80,outputs=outputs,
        module_final_lab_z_m=.01+periods*.028,K=K,gamma_ref=gamma_ref,wavelength_m=lam,
        dx_dy_m=4e-5,transverse_nodes=81,field_nodes_z2=nz,intervals_per_wavelength=20,
        apertures_half_width_m=[.00052,.00056],band_half_fraction=.05,filter_minimum_omega_over_reference=.3,
        main_nodes=1,main_ranks=128,preflight_timeout_s=900,physics_timeout_s=22*3600,analysis_timeout_s=3600,
        historical_correlated_loading=True,new_independent_beam=False,no_resampling=True,no_extra_noise=True,
        common_drift_m=.01,drift_mapping='Exact ballistic transport of original crossing events; Puffin omits PIC self-fields in this drift.',
        magnetic_lattice=magnets,
        magnetic_scope='Qualified aprl-sparc11 binary stopped after its exact first six magnets; one continuous module, no intermediate correctTrans.',
        output_scope='First physical dump per step, before final correctTrans; actual segment planes retained.',
        energy_scope='Full-cohort weighted lab particle changes; all retained forward/paraxial Puffin field, plus Rect/Hann resonant bands; no claim to PIC near-field worldtube identity.',
        physical_acceptance=False)
    save(root/'PLAN.json',plan)
    physical_gamma_mean=float(np.sum(weight.astype(np.longdouble)*gamma.astype(np.longdouble))/np.sum(weight,dtype=np.longdouble))
    info=dict(passed=True,records=len(a),physical_electrons=n,total_charge_C=charge,scale=scale,
        gamma_ref=gamma_ref,K=K,arrival_origin_s=anchor,physical_gamma_initial_mean=physical_gamma_mean,
        initial_mean_gamma_scaled=float(np.sum(native[:,6].astype(np.longdouble)*native[:,5].astype(np.longdouble))/np.sum(native[:,6],dtype=np.longdouble)),
        roundtrip_max_abs_errors=error.tolist(),native_sha256=sha(root/'particles.h5'),
        input_sha256=aprl['input_sha256'],macro_weight_min=float(weight.min()),macro_weight_max=float(weight.max()),
        loading_preserved=True,historical_transverse_duplicates_present=len(np.unique(a['plane_position_m'],axis=0))<len(a),
        unique_source_IDs=len(np.unique(a['source_id']))==len(a),field_nodes_z2=nz,
        minimum_padding_optical_wavelengths=128,additional_slippage_padding_wavelengths=600,
        physical_bunch_duration_s=aprl['physics']['duration_s'],particle_output_has_persistent_ID=False)
    need(info['unique_source_IDs'],'Duplicate original particle IDs')
    save(root/'INPUT_VALIDATED.json',info)
    review=old/'review-20261113-r1'
    shutil.copy2(review/'RESULT.json',root/'reference/RESULT.json')
    for name in ['energy-binomial3','energy-momentum']:
        dest=root/'reference/lab-v2'/name;dest.mkdir(parents=True)
        shutil.copy2(review/'lab-v2'/name/'lab_energy_absolute.csv',dest/'lab_energy_absolute.csv')
    for name in ['FIELD_VALIDATED.json','SMOKE_VALIDATED.json','SOURCE_METHOD_AUDIT.json']:
        shutil.copy2(P/name,root/('QUALIFIED_'+name))
    for name in ['import-check','main']:deck(root,root/name,plan,scale,nz,preflight=name=='import-check')
    refs=[source,binary,P/'FIELD_VALIDATED.json',P/'SMOKE_VALIDATED.json',P/'SOURCE_METHOD_AUDIT.json',
        P/'source/puffin/lib/undulator/bfields.f90']
    files=[root/'PLAN.json',root/'particles.h5',root/'INPUT_VALIDATED.json',*refs,
        *[p for p in (root/'scripts').iterdir() if p.is_file()],
        *[p for p in (root/'reference').rglob('*') if p.is_file()],
        *[p for name in ['import-check','main'] for p in (root/name).iterdir() if p.is_file()]]
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})
    print('PREPARED',root,len(a),charge,nz,periods,flush=True)


if __name__=='__main__':main(sys.argv[1])
