"""One targeted GENESIS representation control of the same closed w64 sample.

Move optical bucket boundaries by half a wavelength. Repack, do not resample;
compensate theta by the known common pi. Physical arrival events remain fixed.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c
from scipy.integrate import simpson

from analyze import read,save,need,table,LD,TWO,W


def spectral_parts(root):
    analysis=root/'analysis';rows=[];aperture=[]
    for plane in ['u07','u77']:
        for program in ['PIC','GENESIS']:
            with np.load(analysis/f'{plane}-{program}.npz') as z:
                x=z['x_m'];y=z['y_m'];frequency=z['frequency_Hz'];values=z['forward_x_spectrum']
            dx=x[1]-x[0];dy=y[1]-y[0]
            wx=np.maximum(0,np.minimum(x+dx/2,.00056)-np.maximum(x-dx/2,-.00056))/dx
            wy=np.maximum(0,np.minimum(y+dy/2,.00056)-np.maximum(y-dy/2,-.00056))/dy
            weights=wy[:,None]*wx[None,:]*dx*dy
            spectrum=(abs(values)**2*weights).sum(axis=(1,2));mid=len(frequency)//2
            for width in [.005,.01,.02,.05]:
                use=abs(frequency/frequency[mid]-1)<=width+1e-14
                rows.append(dict(plane=plane,program=program,half_fractional_bandwidth=width,
                    band_forward_x_J=float(simpson(spectrum[use],x=frequency[use])),
                    exact_carrier_spectral_density_J_per_Hz=float(spectrum[mid])))
    table(root/'spectral_parts.csv',rows)
    return rows


def prepare(root):
    old=W/'results/sparc_random_startup_20261112_r2';d=old/'w64-a'
    need(read(d/'COMPLETE.json')['complete'],'Parent not complete')
    info=read(d/'INPUT_VALIDATED.json');mapping=read(d/'NATIVE_MAPPING.json');plan=read(old/'PLAN.json')
    def sha(path):
        value=hashlib.sha256()
        with Path(path).open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):value.update(block)
        return value.hexdigest()
    need(sha(plan['genesis_binary'])==plan['genesis_binary_sha256'],'GENESIS program changed')
    need(sha(info['input'])==info['input_sha256'],'Physical input changed')
    need(sha(d/'common.par.h5')==mapping['native_sha256'],'Baseline native input changed')
    control=root/'half-bucket';control.mkdir(exist_ok=False)
    with h5py.File(info['input'],'r') as f:
        need(int(f['generation/complete'][()])==1,'Input incomplete');a=f['particles/records'][:]
    lam=LD(plan['wavelength_m']);coordinate=-LD(c)*a['arrival_time_offset_s'].astype(LD)/lam
    offset=LD('.5');integer=np.floor(coordinate-offset).astype(np.int64)
    theta=np.asarray(TWO*(coordinate-offset-integer.astype(LD)),float)
    theta=np.minimum(theta,np.nextafter(float(TWO),0.))
    sort=np.argsort(integer,kind='stable');bins=integer[sort];origin=int(integer.min())-256
    ns=int(np.ceil((int(integer.max())-origin+1+256)/32))*32
    with h5py.File(d/'common.par.h5','r') as parent:
        need(ns==int(parent['slicecount'][0]),'Changed total native window slice count')
    gamma=np.sqrt(1+np.sum(a['proper_velocity']**2,axis=1));indices=[]
    count=0;max_time_error=0.;max_phase_error=0.;charge=0.
    with h5py.File(control/'common.par.h5','x') as f:
        for key,value in [('refposition',float((LD(origin)+offset)*lam)),('slicelength',float(lam)),
            ('slicespacing',float(lam)),('slicecount',ns),('beamletsize',1),('one4one',0)]:f[key]=[value]
        ref=LD(f['refposition'][0])
        for j in range(ns):
            lo,hi=np.searchsorted(bins,[origin+j,origin+j+1]);rows=sort[lo:hi]
            group=f.create_group(f'slice{j+1:06d}');current=len(rows)*info['actual_weight']*1.602176634e-19*c/float(lam)
            group['current']=[current];charge+=current*float(lam)/c
            for key,values in [('x',a['plane_position_m'][:,0]),('y',a['plane_position_m'][:,1]),
                ('px',a['proper_velocity'][:,0]),('py',a['proper_velocity'][:,1]),('gamma',gamma),('theta',theta)]:
                group[key]=values[rows] if len(rows) else [plan['source_config']['phase_space']['gamma_mean'] if key=='gamma' else 0.]
            group['source_id']=a['source_id'][rows] if len(rows) else np.array([0],dtype='u8')
            if len(rows):
                indices.append(j+1);count+=len(rows)
                recovered=-(ref+LD(j)*lam+group['theta'][:].astype(LD)/TWO*lam)/LD(c)
                max_time_error=max(max_time_error,float(np.max(abs(recovered-a['arrival_time_offset_s'][rows].astype(LD)))))
                original_phase=np.exp(-1j*np.asarray(TWO*LD(c)*a['arrival_time_offset_s'][rows].astype(LD)/lam,dtype=float))
                restored_phase=np.exp(1j*(theta[rows]+np.pi))
                max_phase_error=max(max_phase_error,float(np.max(abs(original_phase-restored_phase))))
    need(count==len(a) and max_time_error<1e-23 and max_phase_error<5e-11,'Physical crossing events changed')
    need(abs(charge/(info['physical_electrons']*1.602176634e-19)-1)<2e-12,'Charge changed')
    for name in ['reference.lat','reference.in']:
        (control/name).write_bytes((d/'genesis-off'/name).read_bytes())
    validation=dict(passed=True,physical_input=str(Path(info['input']).resolve()),records=len(a),
        source_ids_and_5D_phase_space_preserved=True,max_reconstructed_time_error_s=max_time_error,
        max_restored_carrier_phase_error=max_phase_error,charge_C=charge,
        only_representation_change='Optical bucket origin +lambda/2; relative theta shifted -pi with physical-time reconstruction retained.',
        origin_shift_wavelength=.5,active_groups=indices,native_slices=ns,one4one=False,nbins=1,shotnoise=False,
        fixed_seed=2026111201,new_seeds=0,new_pic_runs=0)
    save(control/'INPUT_VALIDATED.json',validation)
    return control,indices,plan


def output_compare(root,control,plan):
    frequency=np.load(root/'analysis/initial_sources.npz')['frequency_Hz'];f0=c/plan['wavelength_m'];result=[]
    for plane,dump in [('u07','reference.29.fld.h5'),('u77','reference.309.fld.h5')]:
        with h5py.File(control/dump,'r') as f:
            ns=int(f['slicecount'][0]);ng=int(f['gridpoints'][0]);dx=float(f['gridsize'][0]);dt=float(f['slicespacing'][0])/c
            x=(np.arange(ng)-(ng-1)/2)*dx;ix=np.flatnonzero(abs(x)<.0009);v=np.empty((ns,len(ix),len(ix)),complex)
            for j in range(ns):
                group=f[f'slice{j+1:06d}'];a=group['field-real'][:]+1j*group['field-imag'][:]
                v[j]=a.reshape(ng,ng)[np.ix_(ix,ix)].conj()/dx
        F=(np.exp(2j*np.pi*(frequency-f0)[:,None]*np.arange(ns)*dt)*dt@v.reshape(ns,-1)).reshape(len(frequency),len(ix),len(ix))
        w=np.maximum(0,np.minimum(x[ix]+dx/2,.00056)-np.maximum(x[ix]-dx/2,-.00056))/dx
        area=w[:,None]*w[None,:]*dx*dx;energy=float(simpson((abs(F)**2*area).sum(axis=(1,2)),x=frequency))
        with np.load(root/'analysis'/f'{plane}-GENESIS.npz') as f:
            old=f['forward_x_spectrum'];need(np.allclose(f['x_m'],x[ix],rtol=0,atol=1e-15),'Spatial grid changed')
        before=float(simpson((abs(old)**2*area).sum(axis=(1,2)),x=frequency));mid=len(frequency)//2
        difference=float(np.sqrt(np.sum(abs(F[mid]+old[mid])**2*area)/np.sum(abs(old[mid])**2*area)))
        result.append(dict(plane=plane,old_band_J=before,half_bucket_band_J=energy,
            half_over_old_band_energy=energy/before,carrier_field_difference_after_known_pi_rotation=difference,
            carrier_spectral_power_ratio=float(np.sum(abs(F[mid])**2*area)/np.sum(abs(old[mid])**2*area))))
        np.savez_compressed(control/(plane+'-spectrum.npz'),frequency_Hz=frequency,x_m=x[ix],
            band_intensity_J_per_m2=simpson(abs(F)**2,x=frequency,axis=0))
    table(control/'comparison.csv',result)
    save(control/'RESULT.json',dict(complete=True,comparisons=result,
        scope='Same particles, seed, charge, 5D phase space and physical arrival events; one native bucket-origin change. This tests representation sensitivity, not selection of a preferred origin or physical validation.'))


def main(root):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use a bounded allocation');root=Path(root)
    spectral_parts(root);control,indices,plan=prepare(root)
    command=['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1','--ntasks=32','--ntasks-per-node=32',
        '--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (control/'binding.log').open('x') as stream:
        subprocess.run(command+[str(W/'build/cluster-tools/binding_probe'),'32','1'],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=120)
    save(control/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],time_unix=time.time(),binary=plan['genesis_binary']))
    with (control/'genesis.log').open('x') as stream:
        run=subprocess.run(['timeout','--signal=TERM','--kill-after=30','900',*command,plan['genesis_binary'],'reference.in'],
            cwd=control,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=940)
    need(run.returncode==0 and 'End of Track' in (control/'genesis.log').read_text(),'Control did not finish')
    with h5py.File(control/'common.par.h5','r') as a,h5py.File(control/'initial.par.h5','r') as b:
        for i in indices:
            key=f'slice{i:06d}'
            for field in ['current','x','y','px','py','gamma','theta']:
                need(np.array_equal(a[key][field][:],b[key][field][:]),'Control import changed data')
    output_compare(root,control,plan)
    save(control/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],new_seeds=0,new_pic_runs=0,native_import_bitwise=True))


if __name__=='__main__':main(sys.argv[1])
