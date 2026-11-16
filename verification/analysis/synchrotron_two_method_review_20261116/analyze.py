"""Read-only review of eight completed self-consistent bend/control cases."""
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c,e,h,hbar,mu_0,m_e
from report_discrete import load,calculate
from relative_ledger import construct,write_tables
import whole
import lab

LD=np.longdouble
SHAPES={'F3':'linear-binomial3','MM':'quintic-compensated'}

def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def sha(p):
    hsh=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):hsh.update(b)
    return hsh.hexdigest()
def table(p,data):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with Path(p).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(list(data));w.writerows(zip(*data.values()))
def ident(p):
    s=Path(p).stat();return dict(inode=s.st_ino,size=s.st_size,mtime_ns=s.st_mtime_ns)
def intensity(E,B):return e/(h*mu_0*c)*(abs(E[...,:2])**2+c*c*abs(B[...,:2])**2).sum(axis=-1)
def poynting(E,B):return 2*e/(h*mu_0)*np.real(E[...,0]*B[...,1].conj()-E[...,1]*B[...,0].conj())

def plane(d,label):
    with h5py.File(d/'detectors'/(label+'-substeps.h5')) as f:
        g=f['particle_substep_plane'];assert int(g['complete'][()])==1
        a=g['records'][:];assert len(a)==int(g['committed_records'][()])==100000
    assert np.all(a['crossing_direction']==1) and len(np.unique(a['source_id']))==100000
    return a[np.argsort(a['source_id'])]

def particle_stats(d,out):
    a0=plane(d,'entry');ids=a0['source_id'];K0=a0['kinetic_lab_J'].astype(LD)
    data=[]
    for label,z in [('entry',0.),('audit-base',.02),('magnet-entry',.05),('bend-exit',.0524),
                    ('z0060-particles',.06),('z0100-particles',.1),('z0200-particles',.2),('exit',.3524)]:
        a=plane(d,label);assert np.array_equal(a['source_id'],ids)
        loss=(K0-a['kinetic_lab_J'].astype(LD))/LD(e)
        u=a['proper_velocity'].astype(LD);angle=np.arctan2(u[:,0],u[:,2])
        data.append(dict(plane=label,z_lab_m=z,particles=len(a),charge_C=float(abs(a['charge_C'].astype(LD)).sum()),
            mean_loss_eV=float(loss.mean()),std_loss_eV=float(loss.std()),
            projected_energy_spread_eV=float(a['kinetic_lab_J'].astype(LD).std()/LD(e)),
            mean_x_original_m=float(a['position_m'][:,0].astype(LD).mean()+LD('.00012')),
            sigma_x_m=float(a['position_m'][:,0].astype(LD).std()),sigma_y_m=float(a['position_m'][:,1].astype(LD).std()),
            mean_x_angle_rad=float(angle.mean()),rms_x_angle_rad=float(angle.std())))
    table(out/'particle_statistics.csv',{k:[r[k] for r in data] for k in data[0]})
    return data

def budget(raw,out,case,binary):
    d=raw/case['directory_relative'];o=out/'cases'/case['tag'];o.mkdir(parents=True)
    rows,sums,ranks,steps=load(d/'ledger',d/'ledger/energy-ledger.h5')
    assert ranks==16 and steps==case['steps']
    ix,q,summary=calculate(rows,sums,steps);absolute,pct,norm=construct(rows,sums,ix,q)
    write_tables(o/'boost',ix,absolute,pct,norm)
    table(o/'boost/complete_components.csv',dict(time_box_s=rows['time_box_s'],**absolute))
    a=lambda k:rows[k].astype(LD)
    faces=a('outward_field_energy_face_J')
    phys=dict(particle_change_J=a('delta_particle_kinetic_box_J'),field_change_J=a('delta_field_energy_J'),
        front_outward_minus_injection_J=faces[:,4]-a('moving_inserted_field_J')-a('inlet_overwrite_work_J'),
        back_outward_plus_cropped_J=faces[:,5]+a('moving_cropped_field_J'),side_outward_J=faces[:,:4].sum(axis=1),
        moving_interface_J=-a('moving_remap_interface_correction_J'),
        external_and_constraint_work_signed_J=-a('prescribed_work_J'))
    phys['physical_component_sum_J']=sum(phys.values())
    table(o/'boost/physical_components.csv',dict(time_box_s=a('time_box_s'),**phys))
    values,wr=whole.one(raw,out/'whole',case['directory_relative'])
    table(o/'boost/lab_slice_components.csv',dict(time_box_s=a('time_box_s'),**values))
    _,_,lr=lab.analyze_case(raw,out/'lab',case['directory_relative'],binary)
    lab_source=out/'lab'/case['directory_relative']/'lab_energy_absolute.csv'
    (o/'lab').mkdir();(o/'lab/complete_components.csv').write_bytes(lab_source.read_bytes())
    stats=particle_stats(d,o)
    result=dict(**case,boost=dict(final_discrete_sum_J=float(absolute['total_residual_J'][-1]),
        max_abs_discrete_sum_J=float(max(abs(absolute['total_residual_J']))),
        discrete_sum_range_J=float(np.ptp(absolute['total_residual_J'])),normalization=norm,
        final_physical_component_sum_J=float(phys['physical_component_sum_J'][-1]),**summary),
        lab_slice=wr,fixed_lab=lr,particle_statistics=stats)
    save(o/'RESULT.json',result);print('CASE_DONE',case['tag'],json.dumps(lr['final_components_J']),flush=True)
    return result

def radiation(raw,out,ws,cases):
    root=out/'radiation';root.mkdir();records=[];spectra={};maps={};reference_checks=[]
    library=ws/'results/synchrotron_sampling_grid_controls_20261110_r1/transverse_single_frequency_8eV_modes_r1'
    source=library/'transverse_maps_8eV.npz';audit=read(library/'SOURCE_AUDIT.json')
    assert sha(source)==audit['outputs_sha256'][source.name]
    dense=ws/'results/synchrotron_sample_matched_lw_20261110_r2/sample_matched_spectra.csv'
    refspec=np.genfromtxt(dense,delimiter=',',names=True)
    energies=refspec['photon_eV'];assert len(energies)==321
    for seed in [20261111,20261112]:
        old=ws/f'results/synchrotron_independent_electron_100k_seed{seed}_20261110_r1'
        oldtimes={}
        for kind in ['bend','straight']:
            with h5py.File(old/kind/'detectors/z03524.h5') as f:oldtimes[kind]=f['field_plane/time_s'][:]
        for tag,shape in SHAPES.items():
            bykind={};axis={};shared=None
            for kind in ['bend','straight']:
                d=raw/shape/str(seed)/kind
                with h5py.File(d/'detectors/z03524.h5') as f:
                    g=f['field_plane'];assert int(g['complete'][()])==1
                    t=g['time_s'][:];assert np.array_equal(t,oldtimes[kind]),(tag,seed,kind,'record gate changed')
                    nx,ny=int(g.attrs['nx']),int(g.attrs['ny']);dx,dy=float(g.attrs['dx_m']),float(g.attrs['dy_m'])
                    x=float(g.attrs['x_first_m'])+np.arange(nx)*dx+.00012
                    y=float(g.attrs['y_first_m'])+np.arange(ny)*dy
                    dt=(t[-1]-t[0])/(len(t)-1)
                    assert np.max(abs(t-t[0]-np.arange(len(t))*dt))/dt<1e-4
                    if shared is None:shared=(t,x,y)
                    else:assert all(np.array_equal(a,b) for a,b in zip(shared,(t,x,y)))
                    ix=np.sort(np.argsort(abs(x))[:2]);iy=np.sort(np.argsort(abs(y))[:2])
                    assert np.max(abs(x[ix]))<=dx/2+1e-18
                    ds=g['fields'];patch=ds[:,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
                    phase=np.exp(1j*energies[:,None]*e/hbar*(t-t[0])[None,:])*dt
                    axis[kind]=(phase@patch['electric_V_per_m'].mean(axis=(1,2)),phase@patch['magnetic_T'].mean(axis=(1,2)))
                    ef=np.zeros((ny,nx,3),complex);bf=np.zeros_like(ef)
                    ph=np.exp(1j*8*e/hbar*(t-t[0]))*dt
                    for start in range(0,len(t),32):
                        block=ds[start:start+32]
                        ef+=np.einsum('t,tyxd->yxd',ph[start:start+32],block['electric_V_per_m'],optimize=True)
                        bf+=np.einsum('t,tyxd->yxd',ph[start:start+32],block['magnetic_T'],optimize=True)
                    bykind[kind]=(ef,bf)
                print('SPECTRUM_READ',tag,seed,kind,flush=True)
            paired=poynting(axis['bend'][0]-axis['straight'][0],axis['bend'][1]-axis['straight'][1])
            reference=refspec[f'100k_seed{str(seed)[-2:]}_sample_LW_J_per_eV_m2']
            spectra[f'{tag}_{seed}']=dict(photon_eV=energies,PIC_paired_J_per_eV_m2=paired,
                prescribed_orbit_LW_J_per_eV_m2=reference)
            table(root/f'{tag}_{seed}_axis_spectrum.csv',spectra[f'{tag}_{seed}'])
            with np.load(source) as z:
                assert np.max(abs(z['x_m']-x))<1e-18 and np.max(abs(z['y_m']-y))<1e-18
                assert np.isclose(float(z['dt_s']),dt,rtol=1e-13)
                E=z[f'seed{seed}_bend_LW_E'];B=z[f'seed{seed}_bend_LW_B']
                Ae=z[f'seed{seed}_acceleration_LW_E'];Ab=z[f'seed{seed}_acceleration_LW_B']
            pic=intensity(*bykind['bend']);full=intensity(E,B);acc=intensity(Ae,Ab)
            region=np.zeros_like(pic,dtype=bool);region[10:-10,10:-10]=True
            xx,yy=np.meshgrid(x,y);ap=(xx*xx+yy*yy<=(75e-6)**2)&region
            maps[f'{tag}_{seed}']=dict(x=x,y=y,PIC=pic,LW_full=full,LW_acc=acc,aperture=ap)
            np.savez_compressed(root/f'{tag}_{seed}_maps.npz',**maps[f'{tag}_{seed}'])
            table(root/f'{tag}_{seed}_map.csv',dict(x_original_m=xx.ravel(),y_m=yy.ravel(),
                PIC_intensity_J_per_eV_m2=pic.ravel(),prescribed_full_LW_intensity_J_per_eV_m2=full.ravel(),
                prescribed_acceleration_LW_intensity_J_per_eV_m2=acc.ravel(),in_optical_aperture=ap.ravel()))
            rows=dict(method=tag,seed=seed,photon_eV=8.,PIC_over_full_LW_full_plane=float(pic[region].sum()/full[region].sum()),
                PIC_over_acceleration_LW_aperture=float(pic[ap].sum()/acc[ap].sum()),
                PIC_full_plane_J_per_eV=float(pic[region].sum()*dx*dy),LW_full_plane_J_per_eV=float(full[region].sum()*dx*dy),
                PIC_aperture_J_per_eV=float(pic[ap].sum()*dx*dy),LW_acc_aperture_J_per_eV=float(acc[ap].sum()*dx*dy),
                absolute_map_relative_L2=float(np.linalg.norm(pic[region]-full[region])/np.linalg.norm(full[region])),
                map_cosine=float(np.vdot(pic[region],full[region]).real/(np.linalg.norm(pic[region])*np.linalg.norm(full[region]))))
            for lo,hi in [(4,12),(8,12)]:
                mask=(energies>=lo)&(energies<=hi)
                rows[f'axis_PIC_over_prescribed_LW_{lo}_{hi}eV']=float(np.trapz(paired[mask],energies[mask])/np.trapz(reference[mask],energies[mask]))
            records.append(rows)
            reference_checks.append(dict(method=tag,seed=seed,full_gate_equal_to_original=True,
                original_xy_restored=True,original_input_sha256=read(old/'campaign.json')['input_sha256'],
                theory_is_prescribed_external_orbit_not_new_self_consistent_worldlines=True))
    table(root/'comparison_metrics.csv',{k:[r[k] for r in records] for k in records[0]})
    save(root/'RESULT.json',dict(complete=True,rows=records,reference_checks=reference_checks,
        library_sha256=sha(source),source_library=str(source),source_axis_reference=str(dense),source_axis_reference_sha256=sha(dense),
        map_observable='positive symmetric transverse E/B spectral intensity, not full-frequency signed lab ledger',
        axis_observable='bend-minus-straight complex Fourier field before signed spectral flux; native Yee half-step convention matches old theory',
        no_near_field_subtraction_performed=True,amplitude_fitting=False,
        reference_scope='Exact original sampled input, same finite gate and original coordinates, prescribed magnetic orbit. New self-consistent orbit effects remain part of PIC/reference differences.',
        physical_accuracy_accepted=False))
    return records

def main(raw,out,ws):
    assert os.environ.get('SLURM_JOB_ID'),'Use the compute queue'
    raw,out,ws=map(lambda p:Path(p).resolve(),[raw,out,ws]);out.mkdir(exist_ok=True)
    for file in ['FROZEN.json','PROTECTED.json']:
        for p,digest in read(raw/file).items():assert sha(p)==digest,p
    cases=[]
    for tag,shape in SHAPES.items():
        assert read(raw/f'COMPLETE_{shape}.json')['complete']
        for seed in [20261111,20261112]:
            for kind in ['bend','straight']:
                d=raw/shape/str(seed)/kind;rec=read(d/'COMPLETE.json')
                assert rec['complete'] and read(d/'LAUNCH.json')['exit_code']==0
                assert 'phase=mpi-finalize-complete' in (d/'aprl.log').read_text()
                cases.append(dict(tag=f'{tag}_{seed}_{kind}',method=tag,shape=shape,seed=seed,kind=kind,
                    directory_relative=str(d.relative_to(raw)),steps=rec['actual_steps']))
    protect={str(p):ident(p) for case in cases for p in (raw/case['directory_relative']).rglob('*') if p.is_file()}
    save(out/'PROTECTED_OUTPUTS.json',protect)
    binary=out/'worldtube_integrator'
    subprocess.run(['c++','-std=c++17','-O3','-fopenmp',str(Path(__file__).with_name('worldtube.cpp')),'-o',str(binary)],check=True)
    (out/'geometry').mkdir();lab.geometry_test(out/'geometry',binary)
    results=[budget(raw,out,x,binary) for x in cases]
    metrics=radiation(raw,out,ws,cases)
    for p,before in protect.items():assert ident(p)==before,p
    for p,digest in read(raw/'FROZEN.json').items():assert sha(p)==digest,p
    save(out/'RESULT.json',dict(complete=True,read_only=True,physical_jobs=[19814,19815],cases=results,
        radiation=metrics,original_file_identities_unchanged=len(protect),physics_accepted=False,no_new_physics=True))
    print('ANALYSIS_COMPLETE',flush=True)

if __name__=='__main__':main(*sys.argv[1:])
