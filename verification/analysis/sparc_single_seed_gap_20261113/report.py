"""Summarize the completed one-sample bucket-origin control, without tracking.

Only completed native particle dumps and small analysis caches are read here.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze import W, read, save, need, table


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def cohort_bunching(root):
    parent=W/'results/sparc_random_startup_20261112_r2'
    case=parent/'w64-a';control=root/'half-bucket'
    need(read(control/'COMPLETE.json')['complete'],'Control is not complete')
    prefix=np.asarray(read(case/'COHORT.json')['bin_prefix'],dtype=np.uint64)
    groups=read(control/'INPUT_VALIDATED.json')['active_groups']
    ids=[];phases=[[],[]]
    with h5py.File(control/'common.par.h5','r') as source, \
         h5py.File(control/'initial.par.h5','r') as initial, \
         h5py.File(control/'exit.par.h5','r') as final:
        for index in groups:
            key=f'slice{index:06d}';source_ids=source[key]['source_id'][:]
            need(len(source_ids)==len(final[key]['theta']),'Population changed')
            ids.append(source_ids)
            for j,handle in enumerate((initial,final)):
                phases[j].append(np.exp(1j*handle[key]['theta'][:].astype(np.longdouble)))
    ids=np.concatenate(ids)
    need(len(np.unique(ids))==len(ids),'Repeated source ID')
    bins=np.searchsorted(prefix,ids-1,side='right')-1
    use=(bins>=77)&(bins<154);group=(bins[use]-77)//10
    counts=np.bincount(group,minlength=8).astype(np.longdouble)
    measured=[]
    for label,phase in zip(('entry','u77'),phases):
        phase=np.concatenate(phase)[use];b=np.zeros(8,np.clongdouble)
        np.add.at(b,group,phase)
        rms=np.sqrt((abs(b)**2/counts).sum()/counts.sum())
        measured.append(dict(program='GENESIS half-bucket',plane=label,local_rms_b1=float(rms),
                             cohort='original generation bins [77,154), 10-bin groups, last group 7'))
    baseline=next(x for x in rows(parent/'reviewed-r1/ten_cell_summary.csv') if x['case']=='w64-a')
    for name,prefix_key in [('PIC','PIC'),('GENESIS original','GENESIS-off')]:
        for plane,key in [('entry','initial'),('u77','exit')]:
            measured.append(dict(program=name,plane=plane,
                local_rms_b1=float(baseline[f'{prefix_key}_{key}_local_b1']),cohort=measured[0]['cohort']))
    need(abs(measured[0]['local_rms_b1']/float(baseline['GENESIS-off_initial_local_b1'])-1)<1e-9,
         'Initial physical-cohort BF differs')
    table(root/'matched_cohort_bunching.csv',measured)
    return measured


def summarize(root):
    parent=W/'results/sparc_random_startup_20261112_r2';case=parent/'w64-a'
    plan=read(parent/'PLAN.json');control=root/'half-bucket'
    need(read(root/'analysis/COMPLETE.json')['complete'],'Source analysis incomplete')
    comparisons=read(control/'RESULT.json')['comparisons']
    field=read(root/'analysis/FIELD_RESULT.json')['measurements']
    parts=rows(root/'spectral_parts.csv');summary=[]
    for value in comparisons:
        plane=value['plane']
        pic=next(x for x in field if x['plane']==plane and x['program']=='PIC' and x['half_width_m']==.00056)
        p=next(x for x in parts if x['plane']==plane and x['program']=='PIC')
        g=next(x for x in parts if x['plane']==plane and x['program']=='GENESIS')
        carrier_ratio=float(p['exact_carrier_spectral_density_J_per_Hz'])/float(g['exact_carrier_spectral_density_J_per_Hz'])
        summary.append(dict(plane=plane,PIC_band_J=pic['band_flux_J'],
            GENESIS_original_band_J=value['old_band_J'],GENESIS_half_bucket_band_J=value['half_bucket_band_J'],
            PIC_over_original=pic['band_flux_J']/value['old_band_J'],
            PIC_over_half_bucket=pic['band_flux_J']/value['half_bucket_band_J'],
            half_bucket_over_original=value['half_over_old_band_energy'],
            PIC_over_original_carrier_spectral_power=carrier_ratio,
            half_bucket_over_original_carrier_spectral_power=value['carrier_spectral_power_ratio'],
            PIC_over_half_bucket_carrier_spectral_power=carrier_ratio/value['carrier_spectral_power_ratio']))
    table(root/'comparison_summary.csv',summary)
    bunching=cohort_bunching(root)
    identities={}
    info=read(case/'INPUT_VALIDATED.json');mapping=read(case/'NATIVE_MAPPING.json')
    for path,expected in [(info['input'],info['input_sha256']),
                          (case/'common.par.h5',mapping['native_sha256']),
                          (plan['genesis_binary'],plan['genesis_binary_sha256']),
                          (plan['pic_binary'],plan['pic_binary_sha256'])]:
        identities[str(path)]=sha(path);need(identities[str(path)]==expected,'Original identity changed')
    for manifest in ('SCRIPT_HASHES.json','CONTROL_PLAN.json'):
        data=read(root/manifest)
        for path,expected in data.get('script_sha256',data).items():
            need(sha(path)==expected,'Frozen analysis/control script changed')
    save(root/'VERIFICATION.json',dict(passed=True,original_identity_sha256=identities,
        frozen_scripts_checked=True,matched_initial_cohort_bunching=True,
        native_import_bitwise=read(control/'COMPLETE.json')['native_import_bitwise'],
        source_carrier_difference=read(root/'analysis/SOURCE_RESULT.json')['carrier_B_difference'],
        physical_event_time_error_s=read(control/'INPUT_VALIDATED.json')['max_reconstructed_time_error_s'],
        no_running_H5_read=True,new_physical_seeds=0,new_PIC_runs=0,new_GENESIS_runs=1,
        physical_acceptance=False))
    save(root/'RESULT.json',dict(complete=True,summary=summary,bunching=bunching,
        cause_established='Finite-band GENESIS source/output is sensitive to optical bucket origin for the same physical particles and hard-truncated pulse.',
        not_established='The half-bucket origin is not an exact solution or preferred reference; remaining spectral, transverse and dynamical differences are unresolved.',
        physical_acceptance=False))
    plot(root,summary,bunching)


def plot(root,summary,bunching):
    with np.load(root/'analysis/initial_sources.npz') as f:
        freq=f['frequency_Hz'];offset=100*(freq/freq[len(freq)//2]-1)
        exact=abs(f['exact'][:,0])**2;projected=abs(f['projected'][:,0])**2
        mean=abs(f['analytic_mean'])**2;mean_projected=abs(f['analytic_projected_mean'])**2
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    a=axes[0]
    a.plot(offset,projected*1e6,label='Native bucket source',color='#c46b10',alpha=.85,lw=1.2)
    a.plot(offset,exact*1e6,label='Exact particle source',color='#1769a5',lw=1.2)
    a.plot(offset,mean_projected*1e6,color='#c46b10',ls=':',lw=1,label='Mean profile: bucket')
    a.plot(offset,mean*1e6,color='#1769a5',ls=':',lw=1,label='Mean profile: exact')
    a.set(xlabel='Frequency offset from carrier (%)',ylabel=r'Initial $|B|^2$ ($10^{-6}$)',
          title='Same sample: finite-band source differs',xlim=(-5,5))
    a.text(.03,.98,'Carrier B matches to 1.1e-16\nExact/native band factor = 0.352',
           transform=a.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    a.legend(loc='upper right',fontsize=7.7);a.grid(alpha=.17)
    colors=['#1769a5','#c46b10','#32936f'];labels=['PIC','GENESIS original','GENESIS half-bucket']
    a=axes[1];x=np.arange(2);width=.25
    for j,(key,label,color) in enumerate(zip(['PIC_band_J','GENESIS_original_band_J','GENESIS_half_bucket_band_J'],labels,colors)):
        vals=[r[key]*1e12 for r in summary]
        bar=a.bar(x+(j-1)*width,vals,width,label=label,color=color)
        a.bar_label(bar,fmt='%.1f',fontsize=8,padding=3)
    a.set(xticks=x,xticklabels=['7 periods','77 periods'],ylabel='Band energy (pJ)',
          title='Same physical beam, different buckets',ylim=(0,440))
    a.legend(fontsize=8);a.grid(axis='y',alpha=.17);a.set_axisbelow(True)
    a=axes[2]
    final=[next(r['local_rms_b1'] for r in bunching if r['program']==name and r['plane']=='u77') for name in labels]
    bar=a.bar(np.arange(3),np.array(final)*1e3,color=colors,width=.65)
    a.bar_label(bar,fmt='%.3f',fontsize=9,padding=3)
    initial=next(r['local_rms_b1'] for r in bunching if r['plane']=='entry')*1e3
    a.axhline(initial,color='#555555',ls='--',lw=1,label='Same initial cohort BF')
    a.set(xticks=np.arange(3),xticklabels=['PIC','GENESIS\noriginal','GENESIS\nhalf-bucket'],
          ylabel=r'Local RMS bunching ($10^{-3}$)',title='Exit bunching: same original particle IDs',ylim=(0,max(final)*1.2e3))
    a.legend(fontsize=8,loc='lower left');a.grid(axis='y',alpha=.17);a.set_axisbelow(True)
    fig.suptitle('w64-a: one physical sample, no resampling or amplitude matching',fontsize=14)
    fig.supxlabel('Energy: whole recorded pulse, Rect ±5%, aperture ±560 µm. BF: central original 77 bins, groups of 10.\n'
                  'Changing bucket origin is a diagnostic; the half-bucket result is not a corrected reference.',fontsize=9)
    fig.savefig(root/'bucket_origin_comparison.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':summarize(Path(sys.argv[1]))
