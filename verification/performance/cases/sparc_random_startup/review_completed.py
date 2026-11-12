"""Review closed independent-loading runs without changing their frozen scripts."""
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

from common import load, require, save, sha

C=np.longdouble(299792458)
PI=np.arccos(np.longdouble(-1))
E=np.longdouble('1.602176634e-19')


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def grouped(ids,weight,phase,prefix,first,last,width=1):
    bins=np.searchsorted(prefix,ids-1,side='right')-1
    use=(bins>=first)&(bins<last)
    index=(bins[use]-first)//width
    size=(last-first+width-1)//width
    q=np.zeros(size,np.longdouble);b=np.zeros(size,np.clongdouble)
    np.add.at(q,index,weight[use]);np.add.at(b,index,weight[use]*phase[use])
    good=q>0;total=q.sum();mean=b.sum()/total
    rms2=(abs(b[good])**2/q[good]).sum()/total
    return dict(count=int(use.sum()),charge_e=float(total),global_b1=float(abs(mean)),
                local_rms_b1=float(np.sqrt(rms2)),common_phase_fraction=float(abs(mean)**2/rms2))


def main(root,out):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Use a bounded read-only allocation')
    root=Path(root);out=Path(out);out.mkdir(exist_ok=False)
    plan=load(root/'PLAN.json');pic=[];gen=[];bands=[];summary=[];identities={};cache_error=0.
    checks=[]
    for file in [root/'FROZEN_PREPARATION.json',*[Path(x['directory'])/'FROZEN.json' for x in plan['cases']]]:
        for path,digest in load(file).items():
            resolved=str(Path(path).resolve())
            if resolved not in identities:identities[resolved]=sha(resolved)
            require(identities[resolved]==digest,'Frozen identity differs: '+path)
    for entry in plan['cases']:
        name=entry['name'];d=Path(entry['directory']);require(load(d/'COMPLETE.json')['complete'],'Case incomplete')
        info=load(d/'INPUT_VALIDATED.json');prefix=np.asarray(load(d/'COHORT.json')['bin_prefix'],dtype=np.uint64)
        with h5py.File(info['input'],'r') as f:
            require(int(f['generation/complete'][()])==1,'Input incomplete')
            original=f['particles/records'][:]
        source_ids=original['source_id'];cache=list(csv.DictReader((d/'analysis/pic_bunching.csv').open()))
        for plane,z in plan['particle_planes']:
            with h5py.File(d/'detectors'/(plane+'-particles.h5'),'r') as f:
                g=f['particle_plane'];require(int(g['complete'][()])==1,'Particle plane incomplete')
                a=g['records'][:]
            require(np.array_equal(np.sort(a['source_id']),source_ids),'Missing/duplicated actual source IDs')
            require(np.all(a['crossing_direction']==1),'Non-forward crossing')
            weight=abs(a['charge_C'].astype(np.longdouble))/E
            require(np.max(abs(weight/info['actual_weight']-1))<1e-12,'Charge weight differs')
            phase=np.exp(2j*PI*(C*a['time_s'].astype(np.longdouble)-np.longdouble(z))/np.longdouble(plan['wavelength_m']))
            for scope,lo,hi in [('all',0,231),('core',77,154)]:
                for width in (1,10):
                    row=dict(case=name,program='PIC',plane=plane,z_m=z,scope=scope,group_cells=width,
                             **grouped(a['source_id'],weight,phase,prefix,lo,hi,width))
                    pic.append(row)
                    if width==1:
                        old=next(x for x in cache if x['scope']==scope and x['plane']==plane)
                        for field in ('global_b1','local_rms_b1'):
                            cache_error=max(cache_error,abs(row[field]-float(old[field])))
            checks.append(dict(case=name,plane=plane,all_input_ids_present=True,actual_physical_weights=True))
        mapping=load(d/'NATIVE_MAPPING.json')
        for mode in ('off','longrange'):
            gd=d/('genesis-'+mode);require(load(gd/'COMPLETE.json')['complete'],'GENESIS incomplete')
            with h5py.File(d/'common.par.h5','r') as src,h5py.File(gd/'initial.par.h5','r') as initial,h5py.File(gd/'exit.par.h5','r') as final:
                ids=[];p0=[];p1=[]
                for index in mapping['active_groups']:
                    key=f'slice{index:06d}';g=src[key]
                    for field in ('current','x','y','px','py','gamma','theta'):
                        require(np.array_equal(g[field][:],initial[key][field][:]),'Native import differs')
                    source=g['source_id'][:]
                    require(len(final[key]['theta'])==len(source),'GENESIS changed slice populations')
                    ids.append(source)
                    p0.append(np.exp(1j*initial[key]['theta'][:].astype(np.longdouble)))
                    p1.append(np.exp(1j*final[key]['theta'][:].astype(np.longdouble)))
                ids=np.concatenate(ids);p0=np.concatenate(p0);p1=np.concatenate(p1)
                require(np.array_equal(np.sort(ids),source_ids),'Native source-ID audit incomplete')
            weight=np.full(len(ids),info['actual_weight'],dtype=np.longdouble)
            for plane,z,phase in [('entry',0.,p0),('u77',2.166,p1)]:
                for scope,lo,hi in [('all',0,231),('core',77,154)]:
                    for width in (1,10):
                        gen.append(dict(case=name,program='GENESIS-'+mode,plane=plane,z_m=z,scope=scope,group_cells=width,
                            **grouped(ids,weight,phase,prefix,lo,hi,width)))
            cached=list(csv.DictReader((d/'analysis/genesis_bunching.csv').open()))
            last=[x for x in cached if x['program']=='GENESIS-'+mode][-1]
            fresh=next(x for x in gen if x['case']==name and x['program']=='GENESIS-'+mode and x['scope']=='all' and x['plane']=='u77' and x['group_cells']==1)
            require(abs(fresh['global_b1']-float(last['global_b1']))<2e-10,'Final GENESIS particle BF differs from diagnostics')
            with h5py.File(gd/'reference.out.h5','r') as f:
                z=f['Lattice/zplot'][:];ds=float(f['Global/lambdaref'][0])*float(f['Global/sample'][0])
                energy=f['Field/Global/energy'][:];power=f['Field/power'][:]
                require(np.allclose(power.sum(axis=1)*ds/float(C),energy,rtol=2e-11,atol=1e-24),'GENESIS power/energy mismatch')
                dumpfiles=sorted(gd.glob('reference.*.fld.h5'),key=lambda p:int(p.name.split('.')[1]))
                for dump,(_,expected) in zip(dumpfiles,plan['particle_planes'][2:]):
                    step=int(dump.name.split('.')[1])
                    require(abs(float(z[step])-expected)<1e-8,'Wrong Genesis field dump coordinate')
        energies=load(d/'analysis/pic_spectra.json')
        reference=list(csv.DictReader((d/'analysis/genesis_band_energy.csv').open()))
        for x in energies:
            for mode in ('off','longrange'):
                ref=next(y for y in reference if y['program']=='GENESIS-'+mode and y['plane']==x['plane'])
                bands.append(dict(case=name,plane=x['plane'],z_m=x['z_m'],reference='GENESIS-'+mode,
                    PIC_rect_J=x['rect_band_J'],PIC_hann_J=x['hann_band_J'],GENESIS_rect_J=float(ref['rect_band_J']),
                    PIC_over_GENESIS_rect=x['rect_band_J']/float(ref['rect_band_J'])))
        select=lambda source,plane:next(x for x in source if x['case']==name and x['scope']=='core' and x['group_cells']==1 and x['plane']==plane)
        start=select(pic,'und-entry');end=select(pic,'u77')
        row=dict(case=name,actual_weight=info['actual_weight'],PIC_core_entry_local_b1=start['local_rms_b1'],
            PIC_core_exit_local_b1=end['local_rms_b1'],PIC_core_local_change_pct=100*(end['local_rms_b1']/start['local_rms_b1']-1),
            PIC_core_entry_common_phase_pct=100*start['common_phase_fraction'],PIC_core_exit_common_phase_pct=100*end['common_phase_fraction'])
        for mode in ('off','longrange'):
            subset=[x for x in gen if x['program']=='GENESIS-'+mode]
            a=select(subset,'entry');b=select(subset,'u77')
            row['GENESIS_'+mode+'_core_local_change_pct']=100*(b['local_rms_b1']/a['local_rms_b1']-1)
            row['GENESIS_'+mode+'_exit_local_b1']=b['local_rms_b1']
        s=next(x for x in bands if x['case']==name and x['plane']=='u77' and x['reference']=='GENESIS-off')
        row.update(PIC_exit_rect_J=s['PIC_rect_J'],PIC_exit_hann_J=s['PIC_hann_J'],GENESIS_off_exit_rect_J=s['GENESIS_rect_J'],PIC_over_GENESIS_off_rect=s['PIC_over_GENESIS_rect'])
        summary.append(row)
    require(cache_error<2e-10,'PIC cache BF mismatch')
    table(out/'pic_fixed_cohorts.csv',pic);table(out/'genesis_matched_endpoints.csv',gen)
    table(out/'band_energy.csv',bands);table(out/'summary.csv',summary)
    save(out/'REVIEW.json',dict(complete=True,physical_acceptance=False,summary=summary,particle_checks=checks,
        frozen_unique_files=len(identities),max_independent_PIC_BF_difference=cache_error,
        genesis_id_contract='Source IDs reconstructed from native slice and unchanged particle-array index; one4one=false makes Beam::sort() a no-op; tracking updates particle entries without reordering. All active populations, initial arrays and final global BF independently checked.',
        spectrum_scope='Reuse completed matched-aperture ±5% Rect caches and their Parseval checks; PIC Hann is not compared with GENESIS Rect.',
        limitations=['Weighted random sources, two seeds, cropped finite pulse, one segment; no one-electron/full-FEL convergence.',
                     'Raw PIC total fields; no near-field separation. Residual absolute-band discrepancy remains.',
                     'Historical clone campaign had different source length/grid/noise structure; not a single-variable causal power subtraction.']))
    save(out/'VERIFICATION.json',dict(passed=True,frozen_files=identities,particle_checks=len(checks),
        complete_source_id_coverage=True,physical_charge_checked=True,initial_genesis_import_bitwise=True,
        genesis_end_global_bf_checked=True,field_dump_z_checked=True,genesis_power_energy_checked=True,
        source_scope='Read-only closed outputs; no new simulation',job=os.environ['SLURM_JOB_ID']))
    ten_cell_summary(out,pic,gen,plan)
    plot(out,pic,gen,bands,summary,plan)
    print(json.dumps(summary,indent=2),flush=True)


def ten_cell_summary(out,pic,gen,plan):
    result=[]
    for entry in plan['cases']:
        name=entry['name'];row=dict(case=name)
        for source,program in [(pic,'PIC'),(gen,'GENESIS-off'),(gen,'GENESIS-longrange')]:
            subset=[x for x in source if x['case']==name and x['program']==program and x['scope']=='core' and x['group_cells']==10]
            a=next(x for x in subset if x['plane']=='entry');b=next(x for x in subset if x['plane']=='u77')
            row[program+'_initial_local_b1']=a['local_rms_b1'];row[program+'_exit_local_b1']=b['local_rms_b1']
            row[program+'_local_change_pct']=100*(b['local_rms_b1']/a['local_rms_b1']-1)
        row['PIC_over_GENESIS_off_exit_local_b1']=row['PIC_exit_local_b1']/row['GENESIS-off_exit_local_b1']
        result.append(row)
    table(out/'ten_cell_summary.csv',result)
    return result


def plot(out,pic,gen,bands,summary,plan):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    colors={'PIC':'#167887','GENESIS-off':'#b66324','GENESIS-longrange':'#626975'}
    fig,axes=plt.subplots(2,3,figsize=(15,8),sharex=True)
    fig.subplots_adjust(top=.83,bottom=.16,left=.075,right=.985,wspace=.26,hspace=.34)
    for col,entry in enumerate(plan['cases']):
        name=entry['name'];a=axes[0,col];b=axes[1,col]
        rows=[x for x in pic if x['case']==name and x['scope']=='core' and x['group_cells']==10]
        start=next(x for x in rows if x['plane']=='entry');end=next(x for x in rows if x['plane']=='u77')
        a.plot([(x['z_m']-.01)/.028 for x in rows],[x['local_rms_b1'] for x in rows],'-o',color=colors['PIC'],label='PIC fixed core')
        for mode,marker in [('off','s'),('longrange','^')]:
            rows=[x for x in gen if x['case']==name and x['scope']=='core' and x['group_cells']==10 and x['program']=='GENESIS-'+mode]
            a.plot([(x['z_m']-.01)/.028 for x in rows],[x['local_rms_b1'] for x in rows],linestyle='None',marker=marker,
                   markerfacecolor='none',markersize=8,color=colors['GENESIS-'+mode],label='GENESIS '+mode+' (endpoints)')
        a.set_title(name+' | fixed central cohort',pad=12);a.set_ylabel('Local RMS bunching, 10-cell groups')
        a.set_ylim(0,.013);a.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        change=100*(end['local_rms_b1']/start['local_rms_b1']-1);a.text(.04,.08,f'PIC input → exit: {change:+.2f}%',transform=a.transAxes)
        for program,style in [('PIC','-o'),('GENESIS-off','--s'),('GENESIS-longrange',':^')]:
            selected=[x for x in bands if x['case']==name and x['reference']==('GENESIS-longrange' if program=='GENESIS-longrange' else 'GENESIS-off')]
            b.plot([(x['z_m']-.01)/.028 for x in selected],
                   [x['PIC_rect_J' if program=='PIC' else 'GENESIS_rect_J']*1e9 for x in selected],style,color=colors[program],label=program)
        b.set_ylabel('Whole-pulse ±5% band energy (nJ)');b.set_xlabel('Undulator periods after magnetic entrance')
        b.set_ylim(0,1.15);b.set_xlim(-1,79)
        for ax in (a,b):ax.set_xticks([0,19,37,57,77]);ax.grid(alpha=.18)
    handles=[Line2D([],[],color=colors[p],linestyle=style,marker=marker,label=p)
             for p,style,marker in [('PIC','-','o'),('GENESIS-off','--','s'),('GENESIS-longrange',':','^')]]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.948),ncol=3,frameon=False,fontsize=11)
    fig.suptitle('Independent 6D random loading: first-segment growth',fontsize=16,y=.986)
    fig.text(.5,.022,'Top: fixed central77-cell cohort, grouped by10 (last group7); GENESIS markers show initial/final dumps only.\nBottom: whole-pulse energy, same ±560µm aperture and ±5% Rect band; no amplitude matching. Weighted entrance noise and raw PIC bound fields are retained.',ha='center',fontsize=9)
    fig.savefig(out/'startup_comparison.png',dpi=160);plt.close(fig)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
