"""Closed raw-plane fields sampled at actual substep particle crossings.

This is bilinear cell-centre x/y plus linear recorded-time interpolation. It
is NOT the conservative source-adjoint mean field used by the particle push.
"""
import argparse
import json
from pathlib import Path
import time

import h5py
import numpy as np
from scipy.constants import c,e,m_e

from bf_particle_phase import (save,mean,rms,corr,familymean,within,phase,
                               harmonic_fit,ramp_fit,read)

LD=np.longdouble


def stat(path):
    p=Path(path);s=p.stat()
    return dict(path=str(p.resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns)


def interpolate_block(block,lower_index,ix,iy,fx,fy,ft):
    """Pure interpolator; block is [time,y,x,component]."""
    v0=(1-fx)[:,None]*block[lower_index,iy,ix]+fx[:,None]*block[lower_index,iy,ix+1]
    v1=(1-fx)[:,None]*block[lower_index,iy+1,ix]+fx[:,None]*block[lower_index,iy+1,ix+1]
    first=(1-fy)[:,None]*v0+fy[:,None]*v1
    v0=(1-fx)[:,None]*block[lower_index+1,iy,ix]+fx[:,None]*block[lower_index+1,iy,ix+1]
    v1=(1-fx)[:,None]*block[lower_index+1,iy+1,ix]+fx[:,None]*block[lower_index+1,iy+1,ix+1]
    last=(1-fy)[:,None]*v0+fy[:,None]*v1
    return (1-ft)[:,None]*first+ft[:,None]*last


def sample(path,events,batch=32):
    before=stat(path);N=len(events);out=np.empty((N,6));started=time.monotonic()
    with h5py.File(path,'r') as f:
        g=f['field_plane'];assert int(g['complete'][()])==1
        times=g['time_s'][:];ds=g['fields'];assert len(times)==len(ds)==int(g['committed_samples'][()])
        assert np.all(np.diff(times)>0)
        event_times=events['time_s'].astype(LD);t=times.astype(LD)
        assert np.all(event_times>=t[0]) and np.all(event_times<=t[-1]),'Time extrapolation forbidden'
        lo=np.searchsorted(t,event_times,side='right')-1;lo=np.clip(lo,0,len(t)-2)
        ft=np.asarray((event_times-t[lo])/(t[lo+1]-t[lo]),float)
        x=(events['position_m'][:,0]-float(g.attrs['x_first_m']))/float(g.attrs['dx_m'])
        y=(events['position_m'][:,1]-float(g.attrs['y_first_m']))/float(g.attrs['dy_m'])
        ix=np.floor(x).astype(int);iy=np.floor(y).astype(int);fx=x-ix;fy=y-iy
        assert np.all(ix>=0) and np.all(ix+1<ds.shape[2]) and np.all(iy>=0) and np.all(iy+1<ds.shape[1]),'XY extrapolation forbidden'
        for first in range(int(lo.min()),int(lo.max())+1,batch):
            last=min(first+batch,int(lo.max())+1);m=(lo>=first)&(lo<last)
            if not m.any():continue
            v=ds[first:last+1]
            assert np.isfinite(v['electric_V_per_m']).all() and np.isfinite(v['magnetic_T']).all()
            index=lo[m]-first
            for offset,key in [(0,'electric_V_per_m'),(3,'magnetic_T')]:
                out[m,offset:offset+3]=interpolate_block(v[key],index,ix[m],iy[m],fx[m],fy[m],ft[m])
        metadata=dict(raw_plane=before,recorded_frames=len(t),read_frame_range_inclusive=[int(lo.min()),int(lo.max()+1)],
            events=N,xy_fraction_range=[float(min(fx.min(),fy.min())),float(max(fx.max(),fy.max()))],
            time_fraction_range=[float(ft.min()),float(ft.max())],seconds=time.monotonic()-started,
            magnetic_time_stagger=g.attrs['magnetic_time_stagger'].decode() if isinstance(g.attrs['magnetic_time_stagger'],bytes) else str(g.attrs['magnetic_time_stagger']))
    assert before==stat(path),'Closed raw source changed'
    return out,metadata


def moment(value,w,phi,group):
    v=np.asarray(value,dtype=LD);z=mean(v*phi,w)
    bins=np.asarray(group-group.min(),dtype=np.int64);W=np.bincount(bins,weights=np.asarray(w,float));good=W>0
    H=(np.bincount(bins,weights=np.asarray(w*v*phi.real,float))+1j*np.bincount(bins,weights=np.asarray(w*v*phi.imag,float)))[good]/W[good]
    ma=np.sum(W[good]*abs(H))/W.sum()
    return dict(global_real=float(z.real),global_imag=float(z.imag),global_abs=float(abs(z)),
        global_phase_rad=float(np.angle(z)),local_10slice_RMS=float(np.sqrt(np.sum(W[good]*abs(H)**2)/W.sum())),
        phase_alignment=float(abs(z)/ma) if ma else 0.)


def summarize(value,w,phi,relative_phi,r,group):
    v=np.asarray(value,dtype=LD);vw=within(v,w)
    return dict(mean=float(mean(v,w)),RMS=rms(v,w),within_family_RMS=rms(vw,w),
        first_harmonic_raw=moment(v,w,phi,group),first_harmonic_within=moment(vw,w,phi,group),
        ramp=ramp_fit(vw,w,r),first_harmonic_fit=harmonic_fit(vw,w,phi),
        family_relative_first_harmonic_fit=harmonic_fit(vw,w,relative_phi))


def self_test():
    # A nonseparable-in-position coefficient and nonuniform time fractions
    # expose x/y swaps, incorrect corner weights, and time-bracket mistakes.
    yy,xx=np.meshgrid(np.arange(4),np.arange(5),indexing='ij')
    block=np.empty((3,4,5,3))
    for n in range(3):
        for d in range(3):block[n,:,:,d]=2+d+(3+d)*xx+(5-d)*yy+(7+d)*n
    ix=np.array([0,2,3]);iy=np.array([0,1,2]);fx=np.array([0.,.3,.9]);fy=np.array([.4,.7,.1]);ft=np.array([.2,.8,1.]);lo=np.array([0,1,1])
    got=interpolate_block(block,lo,ix,iy,fx,fy,ft)
    expected=np.column_stack([2+d+(3+d)*(ix+fx)+(5-d)*(iy+fy)+(7+d)*(lo+ft) for d in range(3)])
    assert np.max(abs(got-expected))<2e-14
    return dict(passed=True,max_absolute_error=float(np.max(abs(got-expected))),scope='Interpolation algebra only; not an actual adjoint-force validation')


def run(workspace,output):
    workspace=Path(workspace).resolve();root=workspace/'results/energy_critical_control_20261105';campaign=root/'campaign-r4';out=Path(output).resolve();out.mkdir(exist_ok=False)
    complete=json.loads((campaign/'on/PIC_COMPLETE.json').read_text());assert complete['analysis_complete'] and complete['run']['exit_code']==0
    sel=json.loads((campaign/'INPUT_SELECTION.json').read_text());wl=4.913172351123507e-7
    with h5py.File(campaign/'inputs/particles.h5','r') as f:a=f['particles/records'][:]
    ss=(a['source_id']-1)//1024;mask=(ss>=sel['core_first_slice_inclusive'])&(ss<sel['core_end_slice_exclusive']);a=a[mask]
    ids=a['source_id'];w=a['macro_weight'].astype(LD);assert len(a)==78848
    t0=a['arrival_time_offset_s'].astype(LD);phi=phase(t0,wl);r=within((t0-t0.min())/LD(wl/c),w)
    relative_phi=phase(t0-t0.reshape(-1,16).min(axis=1).repeat(16),wl)
    order=np.argsort(t0.reshape(-1,16),axis=1);rank=np.empty((len(a)//16,16),int);np.put_along_axis(rank,order,np.tile(np.arange(16),(len(a)//16,1)),axis=1);rank=rank.ravel()
    group=(ids-1)//10240;events={};fields={};metadata={};planes={};profiles=[]
    for plane in ['entry','und-entry']:
        p=read(campaign/'on/detectors'/(plane+'-particles-substeps.h5'));index=np.searchsorted(p['source_id'],ids);p=p[index];assert np.array_equal(p['source_id'],ids) and np.array_equal(p['weight'],w)
        events[plane]=p;field,metadata[plane]=sample(campaign/'on/detectors'/(plane+'.h5'),p);fields[plane]=field
        Ex,Ey,Ez,Bx,By,Bz=field.T
        values=dict(Ez_V_per_m=Ez,Ex_V_per_m=Ex,Ey_V_per_m=Ey,
            forward_x=(Ex+c*By)/2,forward_y=(Ey-c*Bx)/2,
            backward_x=(Ex-c*By)/2,backward_y=(Ey+c*Bx)/2)
        planes[plane]={name:summarize(v,w,phi,relative_phi,r,group) for name,v in values.items()}
        for j in range(16):
            m=rank==j;ww=w[m]
            row=dict(plane=plane,temporal_rank=j,relative_time_cycles=float(mean(r[m],ww)))
            for name,v in values.items():row[name]=float(mean(within(np.asarray(v,dtype=LD),w)[m],ww))
            profiles.append(row)
    before,after=events['entry'],events['und-entry'];q=before['charge_C'].astype(LD)/(w*LD(e));assert np.max(abs(q+1))<1e-12
    displacement=after['position_m'].astype(LD)-before['position_m'].astype(LD)
    Emean=(fields['entry'][:,:3].astype(LD)+fields['und-entry'][:,:3].astype(LD))/2
    proxy_components=q[:,None]*Emean*displacement
    dK=(after['kinetic_lab_J'].astype(LD)-before['kinetic_lab_J'].astype(LD))/(w*LD(e));dKwithin=within(dK,w)
    proxy=proxy_components.sum(axis=1);proxywithin=within(proxy,w);ezproxy=within(proxy_components[:,2],w)
    work=dict(actual_scalar_delta_K=summarize(dK,w,phi,relative_phi,r,group),
        endpoint_trapezoid_E_work=summarize(proxy,w,phi,relative_phi,r,group),
        endpoint_trapezoid_Ez_work=summarize(proxy_components[:,2],w,phi,relative_phi,r,group),
        within_actual_vs_proxy_correlation=corr(dKwithin,proxywithin,w),
        within_actual_vs_Ez_proxy_correlation=corr(dKwithin,ezproxy,w),
        within_difference_RMS_eV=rms(dKwithin-proxywithin,w),
        limitation='Two-face endpoint trapezoid proxy, not the pusher work history or the matched orbit-adjoint mean electric field; does not certify causality.')
    for j in range(16):
        m=rank==j;ww=w[m]
        work.setdefault('head_tail_profile',[]).append(dict(temporal_rank=j,relative_time_cycles=float(mean(r[m],ww)),
            actual_scalar_delta_K_eV=float(mean(dKwithin[m],ww)),
            endpoint_trapezoid_E_work_eV=float(mean(proxywithin[m],ww)),
            endpoint_trapezoid_Ez_work_eV=float(mean(ezproxy[m],ww))))
    np.savez_compressed(out/'entry_field_samples.npz',source_id=ids,weight=np.asarray(w,float),initial_arrival_time_s=np.asarray(t0,float),
        entry_fields_lab=fields['entry'],und_entry_fields_lab=fields['und-entry'],actual_scalar_delta_K_lab_eV=np.asarray(dK,float))
    result=dict(complete=True,new_simulations=0,case='on',records=len(a),core_slices=[654,731],
        interpolation='Bilinear raw cell-centre xy and linear recorded lab time. This is a diagnostic proxy, NOT the source-adjoint orbit-mean force.',
        phase_definition='exp(+i 2pi c original_arrival_time/lambda), no waveform fit to final BF; relative phase measured from each family earliest actual input arrival.',
        sampling=metadata,planes=planes,head_tail_profiles=profiles,drift_work_proxy=work,self_test=self_test(),
        caveats=['Ez_lab=Ez_boost for the longitudinal Lorentz axis, but raw spatial/time sampling differs from the particle coupling gather.',
                 'Forward/backward E±cB combinations are directional algebraic proxies, not unique radiation or near-field separation.',
                 'Native stored Yee half-time B retained; no magnetic phase/time correction.',
                 'First harmonic of particle-phase-conditioned field samples is not a full time-frequency filter.',
                 'Field/energy endpoint associations do not establish what force caused later BF growth.'])
    save(out/'ENTRY_FIELD_REVIEW.json',result)
    lines=['# 入射粒子过面场值的只读核对','','仅七点核已闭合工况，中央77切片。双线性xy和线性记录时间插值，**不是实际推进的伴随轨道平均取场**。', '',
      '|面|Ez族内rms/V m⁻¹|Ez头尾斜率/V m⁻¹ cycle⁻¹|Ez一阶拟合峰值/V m⁻¹|','|---|---:|---:|---:|']
    for name,row in planes.items():
        v=row['Ez_V_per_m'];lines.append(f"|{name}|{v['within_family_RMS']:.9g}|{v['ramp']['slope_per_optical_cycle']:.9g}|{v['first_harmonic_fit']['fitted_peak']:.9g}|")
    lines+=['',f"入射到磁入口的实际族内能量rms={work['actual_scalar_delta_K']['within_family_RMS']:.9g} eV；两端E场梯形功代理rms={work['endpoint_trapezoid_E_work']['within_family_RMS']:.9g} eV，相关系数={work['within_actual_vs_proxy_correlation']:.9g}。",'',
      '这只检验已有入场场值是否带有相同头尾/一阶调制，不能替代实际沿轨迹做功记录。E±cB不是近场/辐射分离，且原生B半时间层保持。完整复谐波、局部分组、端点功代理和头尾表在ENTRY_FIELD_REVIEW.json；采样值在entry_field_samples.npz。']
    (out/'ENTRY_FIELD_REVIEW.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({p:planes[p]['Ez_V_per_m'] for p in planes},indent=2),flush=True)
    print(json.dumps(work,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('workspace',type=Path,nargs='?');p.add_argument('output',type=Path,nargs='?');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:print(json.dumps(self_test()))
    else:assert a.workspace is not None and a.output is not None;run(a.workspace,a.output)
