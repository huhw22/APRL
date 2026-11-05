"""Closed-only force-interval accounting for selected complete phase families.

The binary records contain full converged active intervals, not their internal
particle substeps. Force and work budgets apply to the full intervals exactly;
spatial profiles use their endpoint chords. Existing substep-plane records remain
the authority for fixed-plane crossing times.
"""
import argparse
import csv
import ctypes
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import h5py
import numpy as np
import yaml

C = np.longdouble(299792458)
QE = np.longdouble('1.602176634e-19')
LD = np.longdouble
PI = np.arccos(LD(-1))
HEADER = np.dtype([('magic','S8'),('version','<u4'),('recordBytes','<u4'),
    ('headerBytes','<u4'),('endianMarker','<u4'),('rank','<u8'),('thread','<u8'),
    ('selectedIds','<u8'),('reserved0','<u8'),('reserved1','<u8')])
RECORD = np.dtype([('step','<u8'),('sourceId','<u8'),('runtimeId','<u8'),
    ('substeps','<u4'),('flags','<u4'),('time','<f8'),('dt','<f8'),
    ('charge','<f8'),('mass','<f8'),('weight','<f8'),
    ('positionBefore','<f8',(3,)),('positionAfter','<f8',(3,)),
    ('properBefore','<f8',(3,)),('properAfter','<f8',(3,)),
    ('gridElectric','<f8',(3,)),('gridMagnetic','<f8',(3,)),
    ('gammaBefore','<f8'),('gammaAfter','<f8'),
    ('prescribedImpulseElectric','<f8',(3,)),
    ('prescribedImpulseMagnetic','<f8',(3,)),('prescribedDeltaGamma','<f8')])
assert HEADER.itemsize == 64 and RECORD.itemsize == 288


def need(value, text):
    if not value:
        raise RuntimeError(text)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x,np.generic) else str(x))+'\n')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def identity(path):
    p = Path(path)
    s = p.stat()
    return dict(path=str(p.resolve()), bytes=s.st_size, mtime_ns=s.st_mtime_ns,
                sha256=digest(p))


def fnv_helper(output):
    """Check every record byte in C++; the helper never changes input files."""
    tools = output/'checksum-tools'
    tools.mkdir()
    source = tools/'fnv.cpp'
    source.write_text(r'''#include <cstdio>
#include <cstdint>
extern "C" std::uint64_t records_fnv(const char* name,int* error) {
  *error=0;FILE* f=std::fopen(name,"rb");
  if(!f){*error=1;return 0;}if(std::fseek(f,64,SEEK_SET)){*error=2;std::fclose(f);return 0;}
  std::uint64_t h=14695981039346656037ULL;unsigned char b[65536];std::size_t n;
  while((n=std::fread(b,1,sizeof(b),f)))for(std::size_t i=0;i<n;++i){h^=b[i];h*=1099511628211ULL;}
  if(std::ferror(f))*error=3;if(std::fclose(f))*error=4;return h;
}
''')
    library = tools/'fnv.so'
    command = [os.environ.get('CXX','c++'),'-O3','-std=c++11','-shared','-fPIC',
               str(source),'-o',str(library)]
    subprocess.run(command, check=True, capture_output=True, timeout=60)
    lib = ctypes.CDLL(str(library.resolve()))
    method = lib.records_fnv
    method.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int)]
    method.restype = ctypes.c_uint64
    def run(path):
        status = ctypes.c_int()
        result = int(method(os.fsencode(path),ctypes.byref(status)))
        need(status.value == 0,'Cannot read complete checksum stream: '+str(path))
        return result
    return run, dict(method='All record bytes, C++ FNV-1a64; header independently checked',
        source_sha256=digest(source),library_sha256=digest(library),command=command)


def gate_case(case, expected_ranks=None):
    marker = case/'PIC_COMPLETE.json'
    need(marker.exists(),'Missing final PIC_COMPLETE.json; no trace or H5 read')
    closed = json.loads(marker.read_text())
    need(closed.get('analysis_complete') is True and
         closed.get('run',{}).get('exit_code') == 0,
         'PIC did not complete its final normal-stop checks')
    ledger_path = case/'ledger/energy-ledger.h5'
    with h5py.File(ledger_path,'r') as f:
        group=f['energy_ledger']
        need(int(group['complete'][()]) == 1,'Energy ledger not closed')
        ranks=int(group.attrs['mpi_size'])
    if expected_ranks is not None:
        need(ranks == expected_ranks,'Unexpected MPI size')
    return ranks


def load_records(directory, ranks, output):
    directory = Path(directory).resolve()
    expected = [f'force-rank-{r:05d}.json' for r in range(ranks)]
    need(sorted(p.name for p in directory.glob('force-rank-*.json')) == expected,
         'Missing or extra rank completion manifest')
    manifests=[]; selected=None; steps=None; paths=[]; total=0
    for rank,name in enumerate(expected):
        manifest_path=directory/name
        m=json.loads(manifest_path.read_text())
        need(m.get('schema') == 'APRL_particle_force_interval_v1' and
             m.get('frame') == 'boost' and m.get('rank') == rank and
             m.get('complete') is True,'Invalid or unclosed rank manifest')
        need(m['record_bytes'] == 288 and m['header_bytes'] == 64 and
             m['endianness_marker'] == 0x01020304,'Trace ABI differs')
        marker=manifest_path.with_suffix('.complete').read_text().split()
        need(marker == ['1',str(m['step']),str(m['records'])],
             'Rank completion marker differs from manifest')
        ids=np.asarray(m['selected_ids'],dtype=np.uint64)
        need(len(ids)>0 and ids[0]>0 and np.all(np.diff(ids)>0),
             'Selected IDs must be unique, nonzero and sorted')
        if selected is None:
            selected=ids;steps=int(m['step'])
        need(np.array_equal(ids,selected) and int(m['step']) == steps,
             'Ranks disagree on selected IDs or completed step')
        count=0
        for entry in m['files']:
            p=Path(entry['path']).resolve()
            need(p.parent == directory,'Probe stream lies outside selected directory')
            need(p.name == f'force-rank-{rank:05d}-thread-{entry["thread"]:03d}.bin',
                 'Probe stream filename differs from rank/thread')
            need(p.stat().st_size == 64+288*entry['records'],'Truncated trace stream')
            paths.append((rank,p,entry));count+=int(entry['records'])
        need(count == m['records'],'Manifest record counts disagree')
        total+=count;manifests.append(dict(path=str(manifest_path),sha256=digest(manifest_path)))
    families,counts=np.unique((selected-1)//16,return_counts=True)
    need(np.all(counts == 16),'Selected cohort contains incomplete phase families')
    need(len(paths) == len(set(str(p) for _,p,_ in paths)),'Duplicate stream in manifests')
    need(set(directory.glob('force-rank-*-thread-*.bin')) == set(p for _,p,_ in paths),
         'Unlisted or missing binary trace stream')
    checksum,checksum_meta=fnv_helper(output)
    records=np.empty(total,dtype=RECORD);offset=0;files=[]
    for rank,path,entry in paths:
        before=identity(path)
        with path.open('rb') as stream:
            h=np.frombuffer(stream.read(64),dtype=HEADER,count=1)[0]
            need(h['magic'] == b'APRLFP1' and h['version'] == 1 and
                 h['recordBytes'] == 288 and h['headerBytes'] == 64 and
                 h['endianMarker'] == 0x01020304 and h['rank'] == rank and
                 h['thread'] == entry['thread'] and h['selectedIds'] == len(selected) and
                 h['reserved0'] == 0 and h['reserved1'] == 0,'Invalid binary trace header')
            data=np.fromfile(stream,dtype=RECORD)
        need(len(data) == entry['records'],'Record count differs')
        need(checksum(path) == int(entry['records_fnv1a64'],16),'Record FNV checksum differs')
        if len(data):
            need(int(data['step'][0]) == entry['first_step'] and
                 int(data['step'][-1]) == entry['last_step'],'Stream step endpoints differ')
            need(np.all(data['step']>=1) and np.all(data['step']<=steps),
                 'Trace step outside completed run')
            need(np.all(data['flags'] == 1) and np.all(data['substeps']>0),
                 'Missing device-work validation or substep count')
            for key in RECORD.names:
                if RECORD[key].base.kind == 'f':
                    need(np.isfinite(data[key]).all(),'Nonfinite trace field: '+key)
            need(np.all(data['dt']>0) and np.all(data['mass']>0) and
                 np.all(data['weight']>0),'Invalid trace interval or macro particle')
        records[offset:offset+len(data)]=data;offset+=len(data)
        need(before == identity(path),'Closed trace changed while reading')
        files.append(dict(**before,rank=rank,thread=entry['thread'],records=len(data),
                          verified_fnv1a64=entry['records_fnv1a64']))
    need(total>0,'No selected active intervals recorded')
    order=np.lexsort((records['step'],records['time'],records['sourceId']))
    records=records[order]
    need(np.array_equal(np.unique(records['sourceId']),selected),'Selected source IDs missing from trace')
    duplicate=(records['sourceId'][1:]==records['sourceId'][:-1]) & (records['time'][1:]==records['time'][:-1])
    need(not np.any(duplicate),'Duplicate interval/source time across rank streams')
    return records,selected,dict(ranks=ranks,completed_steps=steps,selected_ids=len(selected),
        complete_families=len(families),records=total,manifests=manifests,files=files,
        checksum=checksum_meta)


def frame_from(case, geometry_path, energy_path=None):
    cfg=yaml.safe_load((case/'run.yaml').read_text())
    geometry=json.loads(Path(geometry_path).read_text())
    if energy_path is None:
        energy_path=case/'detectors/entry-transport-energy.h5'
    with h5py.File(energy_path,'r') as f:
        g=f['energy_plane'];need(int(g['complete'][()])==1,'Entry energy plane not closed')
        gamma=LD(g.attrs['boost_gamma']);beta=LD(g.attrs['boost_beta']);origin=LD(g.attrs['time_lab_origin_s'])
    need(abs(gamma-LD(geometry['boost_gamma']))<LD('1e-12')*gamma,'Boost differs from frozen geometry')
    if 'predicted_frame_time_origin_lab_s' in geometry:
        need(abs(origin-LD(geometry['predicted_frame_time_origin_lab_s']))<
             LD(geometry['dt_box_s'])/gamma*LD('1e-6'),'Frame time origin differs')
    result=dict(gamma=gamma,beta=beta,time_origin=origin,
        z_origin=LD(cfg['beam']['reference']['activation_plane_z']),
        translation=LD(geometry['predicted_box_translation_m']))
    return result,cfg,geometry


def positions_lab(records, frame):
    g,b,t0,z0,tr=(frame[k] for k in ['gamma','beta','time_origin','z_origin','translation'])
    t=records['time'].astype(LD);h=records['dt'].astype(LD)
    p0=records['positionBefore'].astype(LD);p1=records['positionAfter'].astype(LD)
    zbefore=z0+g*(p0[:,2]-tr+b*C*t)
    zafter=z0+g*(p1[:,2]-tr+b*C*(t+h))
    # Stable retarded times avoid subtracting nanosecond lab t and z/c.
    taubefore=t0-z0/C+g*(1-b)*(t-(p0[:,2]-tr)/C)
    tauafter=t0-z0/C+g*(1-b)*(t+h-(p1[:,2]-tr)/C)
    return zbefore,zafter,taubefore,tauafter


def budgets(records, frame):
    r=records;u0=r['properBefore'].astype(LD);u1=r['properAfter'].astype(LD)
    g0=np.sqrt(1+np.sum(u0*u0,axis=1));g1=np.sqrt(1+np.sum(u1*u1,axis=1))
    actual_dg=np.sum((u1-u0)*(u1+u0),axis=1)/(g0+g1)
    dx=r['positionAfter'].astype(LD)-r['positionBefore'].astype(LD)
    E=r['gridElectric'].astype(LD);B=r['gridMagnetic'].astype(LD)
    factor=r['charge'].astype(LD)/(r['mass'].astype(LD)*C)
    du_e=factor[:,None]*r['dt'][:,None].astype(LD)*E
    du_b=factor[:,None]*np.cross(dx,B)
    du_de=r['prescribedImpulseElectric'].astype(LD)
    du_db=r['prescribedImpulseMagnetic'].astype(LD)
    dg_grid=factor/C*np.sum(E*dx,axis=1)
    dg_dev=r['prescribedDeltaGamma'].astype(LD)
    du_res=u1-u0-du_e-du_b-du_de-du_db
    dg_res=actual_dg-dg_grid-dg_dev
    g,b=frame['gamma'],frame['beta']
    dz_lab=g*(dx[:,2]+b*C*r['dt'].astype(LD))
    lab_ez=factor/C*E[:,2]*dz_lab
    lab_et=factor/C*g*((E[:,0]+b*C*B[:,1])*dx[:,0]+(E[:,1]-b*C*B[:,0])*dx[:,1])
    lab_grid=g*(dg_grid+b*(du_e[:,2]+du_b[:,2]))
    lab_device=g*(dg_dev+b*(du_de[:,2]+du_db[:,2]))
    lab_actual=g*(actual_dg+b*(u1[:,2]-u0[:,2]))
    rest_eV=r['mass'].astype(LD)*C*C/(r['weight'].astype(LD)*QE)
    channels={
        'lab_grid_Ez_work_eV':lab_ez*rest_eV,
        'lab_grid_Eperp_work_eV':lab_et*rest_eV,
        'lab_grid_total_work_eV':lab_grid*rest_eV,
        'lab_device_work_eV':lab_device*rest_eV,
        'lab_actual_delta_energy_eV':lab_actual*rest_eV,
        'box_grid_electric_delta_uz':du_e[:,2],
        'box_grid_magnetic_delta_uz':du_b[:,2],
        'box_device_electric_delta_uz':du_de[:,2],
        'box_device_magnetic_delta_uz':du_db[:,2],
        'box_actual_delta_uz':u1[:,2]-u0[:,2]}
    errors=dict(max_abs_impulse_residual=float(np.max(abs(du_res))),
        max_abs_gamma_residual=float(np.max(abs(dg_res))),
        max_abs_saved_gamma_onshell_difference=float(max(np.max(abs(g0-r['gammaBefore'])),np.max(abs(g1-r['gammaAfter'])))),
        max_abs_lab_grid_work_split_residual=float(np.max(abs(lab_grid-lab_ez-lab_et))),
        max_abs_static_device_lab_delta_gamma=float(np.max(abs(lab_device))))
    need(errors['max_abs_impulse_residual']<2e-10,'Full-interval force budget does not match actual endpoints')
    need(errors['max_abs_gamma_residual']<2e-10,'Full-interval energy budget does not match actual endpoints')
    need(errors['max_abs_lab_grid_work_split_residual']<2e-10,'Lab grid work does not match longitudinal/transverse split')
    need(errors['max_abs_static_device_lab_delta_gamma']<2e-10,'Static prescribed lab magnet supplies nonzero lab work')
    return channels,errors


def interp(x, xp, fp):
    """Long-double interpolation retains the small retarded-time differences."""
    idx=np.clip(np.searchsorted(xp,x,side='right')-1,0,len(xp)-2)
    f=(x-xp[idx])/(xp[idx+1]-xp[idx])
    return fp[idx]+f*(fp[idx+1]-fp[idx])


def profiles(records, selected, frame, z, wavelength):
    edges=np.r_[0,np.flatnonzero(np.diff(records['sourceId']))+1,len(records)]
    values={};initial_tau=[];weights=[];errors={};continuity={
        'max_time_gap_s':0.,'max_position_gap_m':0.,'max_proper_velocity_gap':0.,
        'max_start_lab_z_m':0.,'minimum_last_lab_z_m':float('inf'),
        'intervals_crossing_profile_endpoint':0,'discarded_intervals_after_endpoint':0,
        'maximum_lab_interval_dz_m':0.}
    tau=np.zeros((len(selected),len(z)),dtype=LD)
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        r=records[a:b]
        need(len(r)>=2,'Selected ID has fewer than two intervals')
        need(np.all(np.diff(r['step'])==1),'Missing field-step interval for selected ID')
        need(len(np.unique(r['runtimeId']))==1,'Selected ID changed runtime identity')
        for key in ['mass','charge','weight']:
            need(np.all(r[key]==r[key][0]),'Macro particle properties changed: '+key)
        dtgap=r['time'][1:].astype(LD)-(r['time'][:-1].astype(LD)+r['dt'][:-1].astype(LD))
        pgap=r['positionBefore'][1:].astype(LD)-r['positionAfter'][:-1].astype(LD)
        ugap=r['properBefore'][1:].astype(LD)-r['properAfter'][:-1].astype(LD)
        timegap=float(np.max(abs(dtgap)));positiongap=float(np.max(abs(pgap)));propergap=float(np.max(abs(ugap)))
        need(timegap<5e-23,'Discontinuous selected-ID time intervals')
        need(positiongap<2e-12,'Discontinuous selected-ID positions')
        need(propergap<2e-12,'Discontinuous selected-ID proper velocities')
        continuity['max_time_gap_s']=max(continuity['max_time_gap_s'],timegap)
        continuity['max_position_gap_m']=max(continuity['max_position_gap_m'],positiongap)
        continuity['max_proper_velocity_gap']=max(continuity['max_proper_velocity_gap'],propergap)
        za,zb,ta,tb=positions_lab(r,frame)
        need(np.all(zb>za) and np.all(np.diff(zb)>0),'Nonforward or nonmonotone lab trajectory')
        need(abs(za[0]-z[0])<LD('1e-8'),'Selected ID does not begin at activation/entry plane')
        need(zb[-1]>=z[-1],'Selected ID never reaches requested profile exit')
        continuity['max_start_lab_z_m']=max(continuity['max_start_lab_z_m'],float(abs(za[0])))
        continuity['minimum_last_lab_z_m']=min(continuity['minimum_last_lab_z_m'],float(zb[-1]))
        continuity['maximum_lab_interval_dz_m']=max(continuity['maximum_lab_interval_dz_m'],float(np.max(zb-za)))
        continuity['intervals_crossing_profile_endpoint']+=int(np.sum((za<z[-1])&(zb>z[-1])))
        continuity['discarded_intervals_after_endpoint']+=int(np.sum(za>=z[-1]))
        channels,budget_errors=budgets(r,frame)
        for key,error in budget_errors.items():errors[key]=max(errors.get(key,0.),error)
        total_uz=sum(np.sum(channels[key],dtype=LD) for key in [
            'box_grid_electric_delta_uz','box_grid_magnetic_delta_uz',
            'box_device_electric_delta_uz','box_device_magnetic_delta_uz'])
        total_lab_work=np.sum(channels['lab_grid_total_work_eV']+channels['lab_device_work_eV'],dtype=LD)
        ugap_total=abs(total_uz-(LD(r['properAfter'][-1,2])-LD(r['properBefore'][0,2])))
        workgap_total=abs(total_lab_work-np.sum(channels['lab_actual_delta_energy_eV'],dtype=LD))
        errors['max_abs_whole_history_delta_uz_residual']=max(errors.get('max_abs_whole_history_delta_uz_residual',0.),float(ugap_total))
        errors['max_abs_whole_history_lab_work_residual_eV_per_electron']=max(errors.get('max_abs_whole_history_lab_work_residual_eV_per_electron',0.),float(workgap_total))
        knots=np.r_[z[0],zb]
        initial_tau.append(ta[0]);weights.append(r['weight'][0])
        tau[i]=interp(z,knots,np.r_[ta[0],tb])
        for key,delta in channels.items():
            if key not in values:values[key]=np.zeros_like(tau)
            values[key][i]=interp(z,knots,np.r_[LD(0),np.cumsum(delta,dtype=LD)])
    initial_tau=np.array(initial_tau,dtype=LD);weights=np.array(weights,dtype=LD)
    omega=2*PI*C/LD(wavelength)
    phase_in=omega*(initial_tau-np.median(initial_tau))
    phase_delta=omega*(tau-initial_tau[:,None])
    phase=phase_in[:,None]+phase_delta
    return values,phase_in,phase_delta,phase,tau,weights,dict(
        numerical_identity_checks=errors,continuity=continuity,
        profile_method='Cumulative full-interval endpoint chords, linearly clipped/interpolated in lab z; no extrapolation',
        cutoff_scope='Only lab z in requested [entry,exit] is reported. All stored intervals are checked, including post-exit intervals.',
        limitation='Grid fields are constant over a recorded interval. Device impulses/work are true substep sums, but splitting a full interval at a profile z is a chord approximation; actual substep-plane times remain authoritative.')


def within_family(array, ids, weights):
    result=np.array(array,dtype=LD,copy=True)
    for fam in np.unique((ids-1)//16):
        mask=(ids-1)//16==fam
        result[mask]-=np.sum(weights[mask,None]*result[mask],axis=0)/np.sum(weights[mask])
    return result


def harmonics(array, ids, weights, phase_in):
    centered=within_family(array,ids,weights)
    family_values=[];family_weights=[]
    basis=np.exp(1j*phase_in.astype(np.clongdouble))
    for fam in np.unique((ids-1)//16):
        mask=(ids-1)//16==fam;w=weights[mask]
        family_values.append(np.sum(w[:,None]*centered[mask]*basis[mask,None],axis=0)/sum(w))
        family_weights.append(sum(w))
    h=np.array(family_values);fw=np.array(family_weights)
    global_h=np.sum(fw[:,None]*h,axis=0)/sum(fw)
    rms=np.sqrt(np.sum(fw[:,None]*abs(h)**2,axis=0)/sum(fw))
    meanabs=np.sum(fw[:,None]*abs(h),axis=0)/sum(fw)
    coherence=np.divide(abs(global_h),meanabs,out=np.zeros_like(meanabs),where=meanabs>0)
    return dict(real=np.asarray(global_h.real,float),imag=np.asarray(global_h.imag,float),
        absolute=np.asarray(abs(global_h),float),family_rms=np.asarray(rms,float),
        phase_alignment=np.asarray(coherence,float)),centered


def core_selection(case, ids):
    plan_path=case/'run-plan.json'
    if not plan_path.exists():return np.ones(len(ids),dtype=bool),dict(scope='All selected complete families')
    plan=json.loads(plan_path.read_text())
    if 'core_first_slice_inclusive' not in plan:return np.ones(len(ids),dtype=bool),dict(scope='All selected complete families')
    lo=int(plan['core_first_slice_inclusive'])*64;hi=int(plan['core_end_slice_exclusive'])*64
    family=(ids-1)//16
    mask=(family>=lo)&(family<hi)
    need(np.any(mask),'No selected families in input central-slice range')
    return mask,dict(scope='Central input slice range; front selected families are startup-only',
        core_first_slice_inclusive=lo//64,core_end_slice_exclusive=hi//64,
        primary_ids=int(np.sum(mask)),startup_only_ids=int(np.sum(~mask)),
        primary_families=int(len(np.unique(family[mask]))),
        startup_only_families=int(len(np.unique(family[~mask]))))


def truth_planes(case, ids, tau, z, weights, wavelength):
    result=[];omega=2*PI*C/LD(wavelength)
    for path in sorted((case/'detectors').glob('*-particles-substeps.h5')):
        with h5py.File(path,'r') as f:
            group=f['particle_substep_plane']
            need(int(group['complete'][()])==1,'Actual substep plane not complete: '+str(path))
            data=group['records'][()]
            if 'committed_records' in group:
                need(int(group['committed_records'][()])==len(data),'Substep plane count differs')
            plane_z=float(group.attrs['plane_z_m'])
        if plane_z<z[0]-1e-12 or plane_z>z[-1]+1e-12:continue
        subset=data[np.isin(data['source_id'],ids)]
        if 'crossing_direction' in subset.dtype.names:subset=subset[subset['crossing_direction']>0]
        order=np.argsort(subset['source_id']);subset=subset[order]
        need(np.array_equal(subset['source_id'],ids),'Actual plane missing or duplicating selected cohort')
        true_tau=subset['time_s'].astype(LD)-LD(plane_z)/C
        approx=np.array([interp(np.array([plane_z],dtype=LD),z,row)[0] for row in tau])
        reference=np.median(true_tau)
        truebf=np.sum(weights*np.exp(1j*omega*(true_tau-reference)))/sum(weights)
        approxbf=np.sum(weights*np.exp(1j*omega*(approx-reference)))/sum(weights)
        difference=approx-true_tau
        result.append(dict(path=str(path),z_lab_m=plane_z,selected_ids=len(ids),
            actual_substep_bf_absolute=float(abs(truebf)),full_interval_chord_bf_absolute=float(abs(approxbf)),
            bf_absolute_difference=float(abs(approxbf-truebf)),
            crossing_time_difference_max_s=float(np.max(abs(difference))),
            crossing_time_difference_rms_s=float(np.sqrt(np.sum(weights*difference**2)/sum(weights))),
            phase_difference_max_rad=float(omega*np.max(abs(difference)))))
    return result


def scalar_at_end(metrics):
    return {key:{name:float(values[-1]) for name,values in metric.items()}
            for key,metric in metrics.items()}


def draw(output,z,metrics,bunching,core):
    os.environ.setdefault('MPLCONFIGDIR',str(output/'matplotlib-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#1f77b4','#d17a22','#5b3d8a','#38855b']
    channels=['lab_grid_Ez_work_eV','lab_grid_Eperp_work_eV','lab_actual_delta_energy_eV','lab_device_work_eV']
    labels=['Grid longitudinal E','Grid transverse E','Actual particle energy','Prescribed magnet']
    for key,label,color in zip(channels,labels,colors):
        axs[0,0].plot(z,metrics[key]['absolute'],label=label,color=color)
        axs[0,1].plot(z,metrics[key]['phase_alignment'],label=label,color=color)
    axs[0,0].set_ylabel('Within-family first harmonic |h₁| (eV/electron)')
    axs[0,1].set_ylabel('Family phase alignment |Σ W h₁| / Σ W |h₁|')
    axs[0,1].set_ylim(-.03,1.03)
    axs[1,0].plot(z,bunching,label='Selected-cohort actual BF',color='#222222')
    axs[1,0].plot(z,metrics['actual_phase_change_rad']['absolute'],label='Within-family phase-change |h₁| (rad)',color='#5b3d8a')
    axs[1,0].set_ylabel('BF / phase-change harmonic (rad)')
    momentum=['box_grid_electric_delta_uz','box_grid_magnetic_delta_uz',
              'box_device_magnetic_delta_uz','box_actual_delta_uz']
    for key,label,color in zip(momentum,['Grid electric','Grid magnetic','Device magnetic','Actual Δu_z'],colors):
        axs[1,1].plot(z,metrics[key]['absolute'],label=label,color=color)
    axs[1,1].set_ylabel('Within-family first harmonic |h₁(Δu_z)|')
    for ax in axs.flat:
        ax.set_xlabel('Lab propagation position z (m)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Actual converged-interval force/work diagnostics — '+str(core.get('primary_families','all'))+' central families\nHarmonic association does not establish a causal counterfactual',fontsize=12)
    fig.savefig(output/'force_phase_profiles.png',dpi=180)
    fig.savefig(output/'force_phase_profiles.pdf');plt.close(fig)


def analyze(case,output,geometry_path,probe_directory=None,energy_plane=None,
            wavelength=4.913172351123507e-7,exit_z=None,expected_ranks=None,points=161):
    case=Path(case).resolve();output=Path(output).resolve()
    need(not output.exists(),'Refusing to overwrite analysis output')
    ranks=gate_case(case,expected_ranks)
    # All rank completion markers are checked before opening binary intervals.
    output.mkdir(parents=True)
    records,ids,trace=load_records(probe_directory or case/'force-probe',ranks,output)
    frame,cfg,geometry=frame_from(case,geometry_path,energy_plane)
    endpoint=float(geometry['last_magnet_z_m'] if exit_z is None else exit_z)
    need(endpoint>0 and points>=3 and wavelength>0,'Invalid propagation range, points or wavelength')
    z=np.linspace(0,endpoint,points,dtype=LD)
    # Ensure actual comparison planes and magnet entrance are interpolation knots.
    extra=[geometry.get('first_magnet_z_m',0),endpoint/2]
    for path in sorted((case/'detectors').glob('*-particles-substeps.h5')):
        with h5py.File(path,'r') as f:
            group=f['particle_substep_plane'];need(int(group['complete'][()])==1,'Unclosed actual plane')
            extra.append(float(group.attrs['plane_z_m']))
    z=np.unique(np.r_[z,np.array([v for v in extra if 0<=v<=endpoint],dtype=LD)])
    values,phase_in,phase_delta,phase,tau,weights,checks=profiles(records,ids,frame,z,wavelength)
    mask,selection=core_selection(case,ids);primary_ids=ids[mask];w=weights[mask]
    metrics={};centered={};primary_values={key:array[mask] for key,array in values.items()}
    primary_values['actual_phase_change_rad']=phase_delta[mask]
    for key,array in primary_values.items():
        metrics[key],centered[key]=harmonics(array,primary_ids,w,phase_in[mask])
    bunching=np.asarray(abs(np.sum(w[:,None]*np.exp(1j*phase[mask]),axis=0)/sum(w)),float)
    truth=truth_planes(case,primary_ids,tau[mask],z,w,wavelength)
    correlation={}
    p=centered['actual_phase_change_rad'][:,-1]
    for key,array in centered.items():
        v=array[:,-1];pv=np.sum(w*p*v);den=np.sqrt(np.sum(w*p*p)*np.sum(w*v*v))
        correlation[key]=None if den==0 else float(pv/den)
    with (output/'harmonics.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['z_lab_m','channel','real','imag','absolute','family_rms','phase_alignment'])
        for key,metric in metrics.items():
            for j,zj in enumerate(z):writer.writerow([float(zj),key]+[metric[n][j] for n in ['real','imag','absolute','family_rms','phase_alignment']])
    np.savez_compressed(output/'profiles.npz',z_lab_m=np.asarray(z,float),source_id=ids,
        primary_mask=mask,weight=np.asarray(weights,float),phase_in_rad=np.asarray(phase_in,float),
        actual_phase_change_rad=np.asarray(phase_delta,float),retarded_time_s=np.asarray(tau,float),
        selected_primary_bf=bunching,**{key:np.asarray(array,float) for key,array in values.items()})
    draw(output,np.asarray(z,float),metrics,bunching,selection)
    result=dict(schema='APRL_force_interval_analysis_v1',analysis_complete=True,
        scope='Observed converged-interval identities and within-family harmonics; no causal force-removal prediction',
        case=str(case),trace=trace,frame={key:float(val) for key,val in frame.items()},
        geometry_identity=identity(geometry_path),reader_identity=identity(__file__),
        selection=selection,wavelength_m=wavelength,exit_z_lab_m=endpoint,
        checks=checks,actual_substep_plane_comparisons=truth,
        final_harmonics=scalar_at_end(metrics),
        final_within_family_phase_correlation=correlation,
        selected_primary_bf_entry=float(bunching[0]),selected_primary_bf_exit=float(bunching[-1]),
        harmonic_definition='h_f = weighted mean[(A_i-weighted family mean A) exp(+i omega tau_entry_i)]; global is family-charge-weighted mean. Reported |h| is a Fourier coefficient, not 2|h| sinusoid peak.',
        initial_state='Each particle uses its actual first active interval at lab z=0; work/energy/momentum deltas are measured from that state. Every harmonic is centered within its original complete 16-clone family.',
        evidence_limits=['Static lab prescribed-magnet work is checked independently through gamma_b*(delta_gamma_device+beta_b*delta_u_z_device).',
            'Full-interval chord crossing times are compared with actual internal-substep plane times; they do not replace those truth records.',
            'Longitudinal/transverse grid E work is reconstructed by transforming the actual grid fields and actual lab displacement; it is not inferred from a final energy residual.',
            'Work/phase correlation and coherent first harmonics identify candidate driving channels; field self-consistency means they do not establish which field component first caused the bunching.'])
    save(output/'RESULT.json',result)
    (output/'README.md').write_text('# Converged-interval force and phase accounting\n\n'
        'Only normally closed PIC output and complete traces from every MPI rank were read. Every record byte was checked against its recorded FNV-1a checksum, and all stream identities were checked before and after reading.\n\n'
        f'The primary curves contain {selection.get("primary_ids",len(ids))} selected particles in complete 16-particle phase families. All {len(ids)} selected IDs are checked for uninterrupted time, position and momentum histories across rank/thread migration. Front startup samples are excluded from the primary curves when the run plan defines a central range.\n\n'
        'The grid electric and magnetic impulses, prescribed-device impulses and electric work are computed independently and compared with the recorded actual endpoints. The lab grid work is also independently split into longitudinal and transverse electric work after Lorentz transformation. The static lab device must supply approximately zero lab work; its boost electric and magnetic impulses need not vanish.\n\n'
        'Profiles use full-interval endpoint chords and clip the interval crossing the requested last plane. Device contributions within such a partial interval are linearly apportioned, which is an interpolation approximation. The recorded actual internal-substep crossing times remain authoritative; their differences are listed in RESULT.json. No predictor trajectory is used.\n\n'
        'Each harmonic removes the weighted mean within its original 16-particle family and uses the actual entry arrival phase. The plotted coefficient is |h₁|, half the peak amplitude of a pure sinusoid. Alignment measures the agreement of harmonic phase between families. These observations test an association; they are not a field-removal counterfactual or proof of a unique cause.\n\n'
        '![Force/work and phase profiles](force_phase_profiles.png)\n\n'
        'The full numerical identities, input hashes, cohort selection, endpoint values and actual substep-plane comparisons are in [RESULT.json](RESULT.json); [harmonics.csv](harmonics.csv) and [profiles.npz](profiles.npz) retain the plotted data.\n')
    return result


def fixture_fnv(data):
    value=14695981039346656037
    for byte in data:value=((value^byte)*1099511628211)&((1<<64)-1)
    return format(value,'x')


def synthetic_case(case,kind):
    """Independent on-shell electric and Lorentz-transformed magnetic fixtures."""
    case.mkdir(parents=True);(case/'force-probe').mkdir();(case/'ledger').mkdir();(case/'detectors').mkdir()
    gb=LD(2);beta=np.sqrt(LD(3))/2;me=LD('9.1093837139e-31');wave=LD('4.913172351123507e-7')
    frame=dict(gamma=gb,beta=beta,time_origin=LD(0),z_origin=LD(0),translation=LD(0))
    allrecords=[];endpoints=[]
    for j in range(16):
        weight=LD(100+j);q=-QE*weight;mass=me*weight
        tl=LD(j)*wave/(16*C);xlab=np.array([LD(j)*LD('1e-7'),0,0],dtype=LD)
        ulab=np.array([LD('.2') if kind=='magnetic' else 0,0,3],dtype=LD)
        glab=np.sqrt(1+sum(ulab**2))
        t=gb*(tl-beta*xlab[2]/C)
        x=np.array([xlab[0],0,gb*(xlab[2]-beta*C*tl)],dtype=LD)
        u=np.array([ulab[0],0,gb*(ulab[2]-beta*glab)],dtype=LD)
        for step in [1,2]:
            r=np.zeros((),dtype=RECORD);r['sourceId']=j+1;r['runtimeId']=100+j;r['step']=step
            r['substeps']=1;r['flags']=1;r['time']=t;r['charge']=q;r['mass']=mass;r['weight']=weight
            r['positionBefore']=x;r['properBefore']=u;r['gammaBefore']=np.sqrt(1+sum(u*u))
            if kind=='electric':
                h=LD('3e-11');E=np.array([0,0,LD(200)*(1+LD('.1')*np.sin(2*PI*j/16))],dtype=LD)
                u1=u+q*h*E/(mass*C);g0=np.sqrt(1+sum(u*u));g1=np.sqrt(1+sum(u1*u1))
                dx=C*h*(u+u1)/(g0+g1);r['gridElectric']=E
            else:
                h_lab=LD('1e-10');angle=LD('.003')*(1+LD('.1')*np.sin(2*PI*j/16))
                ulab1=np.array([ulab[0]*np.cos(angle)+ulab[2]*np.sin(angle),0,
                    -ulab[0]*np.sin(angle)+ulab[2]*np.cos(angle)],dtype=LD)
                dxlab=C*h_lab*(ulab+ulab1)/(2*glab)
                h=gb*(h_lab-beta*dxlab[2]/C)
                dx=np.array([dxlab[0],0,gb*(dxlab[2]-beta*C*h_lab)],dtype=LD)
                u1=np.array([ulab1[0],0,gb*(ulab1[2]-beta*glab)],dtype=LD)
                by=-2*mass*glab*np.tan(angle/2)/(q*h_lab)
                ee=np.array([-gb*beta*C*by,0,0],dtype=LD);bb=np.array([0,gb*by,0],dtype=LD)
                r['prescribedImpulseElectric']=q*h*ee/(mass*C)
                r['prescribedImpulseMagnetic']=q*np.cross(dx,bb)/(mass*C)
                r['prescribedDeltaGamma']=q*sum(ee*dx)/(mass*C*C)
                ulab=ulab1
            r['dt']=h;r['positionAfter']=x+dx;r['properAfter']=u1
            r['gammaAfter']=np.sqrt(1+sum(u1*u1));allrecords.append(r)
            t+=h;x+=dx;u=u1
        endpoints.append(float(positions_lab(np.array([r],dtype=RECORD),frame)[1][0]))
    endpoint=.9*min(endpoints)
    geometry=dict(boost_gamma=2.,boost_beta=float(beta),dt_box_s=3e-11,
        predicted_frame_time_origin_lab_s=0.,predicted_box_translation_m=0.,
        first_magnet_z_m=.01,last_magnet_z_m=endpoint)
    save(case/'GEOMETRY.json',geometry)
    (case/'run.yaml').write_text(yaml.safe_dump({'beam':{'reference':{'activation_plane_z':0.}}}))
    save(case/'PIC_COMPLETE.json',dict(analysis_complete=True,run=dict(exit_code=0)))
    for name,group_name,attrs in [('ledger/energy-ledger.h5','energy_ledger',{'mpi_size':2}),
        ('detectors/entry-transport-energy.h5','energy_plane',
        {'boost_gamma':2.,'boost_beta':float(beta),'time_lab_origin_s':0.,'z_lab_m':0.})]:
        with h5py.File(case/name,'w') as f:
            g=f.create_group(group_name);g.create_dataset('complete',data=1)
            for key,val in attrs.items():g.attrs[key]=val
    for rank in [0,1]:
        data=np.array([r for r in allrecords if r['step']==rank+1],dtype=RECORD)
        h=np.zeros((),dtype=HEADER);h['magic']=b'APRLFP1';h['version']=1;h['recordBytes']=288
        h['headerBytes']=64;h['endianMarker']=0x01020304;h['rank']=rank;h['selectedIds']=16
        path=case/f'force-probe/force-rank-{rank:05d}-thread-000.bin'
        path.write_bytes(h.tobytes()+data.tobytes())
        m=dict(schema='APRL_particle_force_interval_v1',frame='boost',rank=rank,complete=True,
            step=2,record_bytes=288,header_bytes=64,endianness_marker=0x01020304,
            selected_ids=list(range(1,17)),records=16,
            files=[dict(path=str(path.resolve()),thread=0,records=16,first_step=rank+1,
                last_step=rank+1,records_fnv1a64=fixture_fnv(data.tobytes()))])
        save(case/f'force-probe/force-rank-{rank:05d}.json',m)
        (case/f'force-probe/force-rank-{rank:05d}.complete').write_text('1 2 16\n')
    return endpoint


def self_test(output):
    output=Path(output).resolve();need(not output.exists(),'Refusing to overwrite self-test')
    output.mkdir(parents=True);results={}
    for kind in ['electric','magnetic']:
        case=output/(kind+'-case');endpoint=synthetic_case(case,kind)
        result=analyze(case,output/(kind+'-analysis'),case/'GEOMETRY.json',points=17)
        errors=result['checks']['numerical_identity_checks']
        need(errors['max_abs_impulse_residual']<2e-14,'Synthetic impulse identity failed')
        need(errors['max_abs_gamma_residual']<2e-14,'Synthetic work identity failed')
        need(errors['max_abs_lab_grid_work_split_residual']<2e-14,'Synthetic lab split failed')
        if kind=='magnetic':
            need(errors['max_abs_static_device_lab_delta_gamma']<2e-14,'Static lab magnet control failed')
            need(result['final_harmonics']['lab_grid_total_work_eV']['absolute']==0.,'No-grid fixture obtained grid work')
        else:
            need(result['final_harmonics']['lab_grid_Ez_work_eV']['absolute']>.01,'Electric fixture lacks phase modulation')
            need(result['final_harmonics']['lab_grid_Eperp_work_eV']['absolute']==0.,'Pure Ez fixture obtained transverse work')
        results[kind]=dict(passed=True,endpoint_z_m=endpoint,errors=errors,
            continuity_across_rank_migration=result['checks']['continuity'])
    # Gating and corruption controls use only our tiny disposable synthetic files.
    case=output/'electric-case';marker=case/'PIC_COMPLETE.json';text=marker.read_text()
    marker.unlink()
    try:
        try:gate_case(case);raise AssertionError('Unclosed run was accepted')
        except RuntimeError:pass
    finally:marker.write_text(text)
    manifest=case/'force-probe/force-rank-00001.json';text=manifest.read_text();manifest.unlink()
    try:
        tmp=output/'missing-rank';tmp.mkdir()
        try:load_records(case/'force-probe',2,tmp);raise AssertionError('Missing rank was accepted')
        except RuntimeError:pass
    finally:manifest.write_text(text)
    stream=case/'force-probe/force-rank-00001-thread-000.bin';raw=stream.read_bytes()
    altered=bytearray(raw);altered[-1]^=1;stream.write_bytes(altered)
    try:
        tmp=output/'corrupted-record';tmp.mkdir()
        try:load_records(case/'force-probe',2,tmp);raise AssertionError('Corrupted checksum was accepted')
        except RuntimeError as error:need('checksum' in str(error),'Corruption failed for an unrelated reason')
    finally:stream.write_bytes(raw)
    results['negative_controls']={'unclosed_pic_rejected':True,'missing_rank_rejected':True,'record_checksum_corruption_rejected':True}
    save(output/'SELF_TEST.json',dict(passed=True,tests=results,
        meaning='Reader numerical identities and closed-file checks only; no claim of PIC physics acceptance'))
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',nargs='?');parser.add_argument('output',nargs='?')
    parser.add_argument('--geometry');parser.add_argument('--probe-directory');parser.add_argument('--energy-plane')
    parser.add_argument('--wavelength-m',type=float,default=4.913172351123507e-7)
    parser.add_argument('--exit-z-m',type=float);parser.add_argument('--expected-ranks',type=int)
    parser.add_argument('--points',type=int,default=161);parser.add_argument('--self-test',metavar='OUTPUT')
    args=parser.parse_args()
    if args.self_test:
        need(not args.case and not args.output,'Self-test takes no case/output positional paths')
        self_test(args.self_test);print('Reader self-test PASS: '+args.self_test);return
    need(args.case and args.output and args.geometry,'Case, new output directory and --geometry are required')
    output_new=not Path(args.output).exists()
    try:
        result=analyze(args.case,args.output,args.geometry,args.probe_directory,args.energy_plane,
            args.wavelength_m,args.exit_z_m,args.expected_ranks,args.points)
    except (RuntimeError,OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        if output_new:
            Path(args.output).mkdir(parents=True,exist_ok=True)
            save(Path(args.output)/'RESULT.json',dict(analysis_complete=False,case=str(Path(args.case).resolve()),
                error=str(error),scope='No successful analysis or physics conclusion is reported'))
        print('Force-probe analysis incomplete: '+str(error),file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps(dict(analysis_complete=result['analysis_complete'],output=str(Path(args.output).resolve()),
        selected_ids=result['trace']['selected_ids'],primary_cohort=result['selection'],
        checks=result['checks']['numerical_identity_checks']),indent=2))


if __name__=='__main__':main()
