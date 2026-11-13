"""Full laboratory worldtube: seven planes, four sides and two null-time caps."""
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from worldtube_reference import Tube

LD=np.longdouble


def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def table(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def config(path,m,states,zs,lo,hi):
    values=[len(states),int(m['nz']),m['gamma'],m['c_m_per_s'],m['dz_box_m'],
            m['frame_lab_z_origin_m'],m['frame_box_z_origin_m'],lo,hi,len(zs),*zs]
    path.write_text(' '.join(str(x) for x in values)+'\n')


def geometry_test(root,binary):
    d=root/'geometry-test';d.mkdir();g=LD('929.1305492858543');c=LD(299792458)
    b=np.sqrt((g-1)*(g+1))/g;eta=1/(g*(1+b));dt=LD('4.1839819784122236e-15');dz=c*dt;nt,nz=97,129
    m=dict(gamma=g,c_m_per_s=c,dz_box_m=dz,frame_lab_t_origin_s=LD(0),
           frame_lab_z_origin_m=LD(0),frame_box_z_origin_m=LD(0),nz=nz)
    st=np.zeros(nt,dtype=[('time_box_s',LD),('origin_z_box_m',LD)])
    st['time_box_s']=np.arange(nt,dtype=LD)*dt;st['origin_z_box_m']=-64*dz-np.arange(nt,dtype=LD)*dz
    q=np.zeros((nt,nz,6),dtype=LD)
    for n in range(nt):
        box=st['origin_z_box_m'][n]+(np.arange(nz,dtype=LD)+LD('.5'))*dz
        tau=eta*(st['time_box_s'][n]-box/c);z=g*(box+b*c*st['time_box_s'][n]);T=LD('1e-15')
        q[n,:,0]=2*tau/T;q[n,:,1]=c*LD('.3')*tau/T-z/T
        q[n,:,2:4]=LD('.2')/T;q[n,:,4:6]=-LD('.7')/T
    q=np.asarray(q,dtype='<f8');q.tofile(d/'field_lines.bin')
    with (d/'states.csv').open('w') as f:
        f.write('step,time_box_s,origin_z_box_m\n')
        for i in range(nt):f.write(f"{i},{st['time_box_s'][i]},{st['origin_z_box_m'][i]}\n")
    zs=np.array(['-.02','-.01','0','.01','.02'],dtype=LD);tube=Tube(m,st,q.astype(LD));lo,hi=tube.gate(zs[0],zs[-1])
    config(d/'config.txt',m,st,zs,lo,hi)
    subprocess.run([str(binary),str(d),str(d/'config.txt'),str(d/'fields.csv')],check=True)
    rows=np.genfromtxt(d/'fields.csv',delimiter=',',names=True,dtype=LD);maximum=LD(0)
    for j,z in enumerate(zs):
        ref=tube.balance(zs[0],z,lo,hi)
        pairs={'field_entrance_signed_J':'field_in','field_exit_signed_J':'field_out',
               'time_cap_early_signed_J':'cap_early','time_cap_late_signed_J':'cap_late',
               'transverse_xminus_signed_J':'xminus','transverse_xplus_signed_J':'xplus',
               'transverse_yminus_signed_J':'yminus','transverse_yplus_signed_J':'yplus'}
        scale=max(abs(ref['field_in']),abs(ref['field_out']),abs(ref['cap_change']),abs(ref['side_sum']),LD('1e-100'))
        maximum=max(maximum,*[abs(rows[name][j]-ref[key])/scale for name,key in pairs.items()])
    last=rows[-1];field_names=[x for x in rows.dtype.names if x.endswith('_J')]
    net=sum(last[k] for k in field_names);scale=max(abs(last[k]) for k in field_names)
    wrong=net-2*(last['time_cap_late_signed_J']+last['time_cap_early_signed_J'])
    assert maximum<LD('1e-10') and abs(net)/scale<LD('1e-10') and abs(wrong)/scale>1
    answer=dict(passed=True,python_frozen_geometry_max_difference=float(maximum),
        manufactured_complete_sum_relative=float(abs(net)/scale),wrong_time_cap_sign_rejected=True,
        scope='Independent conserved manufactured current; ports, sides and time caps all nonzero. Numerical integration test, no new PIC.')
    save(root/'GEOMETRY_CHECK.json',answer)


def particle_plane(d,name,half):
    with h5py.File(d/'detectors'/(name+'-particles-substeps.h5')) as f:
        g=f['particle_substep_plane'];assert int(g['complete'][()])==1
        a=g['records'][:];assert int(g['committed_records'][()])==len(a);z=LD(g.attrs['plane_z_m'])
    assert np.isfinite(a['kinetic_lab_J']).all() and np.isfinite(a['time_s']).all()
    inside=(a['position_m'][:,0]>=-half)&(a['position_m'][:,0]<half)&(a['position_m'][:,1]>=-half)&(a['position_m'][:,1]<half)
    return z,a,inside


def analyze_case(parent,out,case,binary):
    d=parent/case;o=out/case;o.mkdir(parents=True,exist_ok=True)
    assert read(d/'PIC_COMPLETE.json')['complete']
    plan=read(parent/'PLAN.json');expected=plan['steps']
    assert (d/'worldtube/COMPLETE').read_text().split()==['1',str(expected)]
    m=json.loads((d/'worldtube/METADATA.json').read_text(),parse_float=LD)
    states=np.genfromtxt(d/'worldtube/states.csv',delimiter=',',names=True,dtype=LD)
    assert len(states)==expected and np.array_equal(states['step'],np.arange(expected))
    assert (d/'worldtube/field_lines.bin').stat().st_size==len(states)*int(m['nz'])*6*8
    # Geometry setup needs only array shape; never copies the 14GB field map.
    dummy=type('ShapeOnly',(),{'shape':(len(states),int(m['nz']),6)})()
    tube=Tube(m,states,dummy)
    half=min(m['physical_half_width_x_m'],m['physical_half_width_y_m'])
    labels=['audit-entry',*[x[0] for x in plan['worldtube_exit_planes']]]
    zs=[LD(str(plan['worldtube_entrance_z_m'])),*[LD(str(x[1])) for x in plan['worldtube_exit_planes']]]
    lo,hi=tube.gate(zs[0],zs[-1]);assert hi>lo
    source=np.max(states['maximum_frozen_z_lab_m'][states['frozen_count']>0])+9*tube.g*tube.dz
    lower=tube.zorig+tube.g*(tube.o-tube.boxorig+tube.b*tube.c*tube.t);upper=lower+tube.g*tube.dz*tube.nz
    assert source<zs[0] and np.max(lower)+9*tube.g*tube.dz<zs[0] and np.min(upper)-9*tube.g*tube.dz>zs[-1]
    cfg=o/'tube-config.txt';config(cfg,m,states,zs,lo,hi)
    # Reuse only an exactly matching completed first-pass surface integration.
    prior=out.parent/'lab'/case
    reuse=(prior/'field_surfaces.csv').is_file() and (prior/'tube-config.txt').read_bytes()==cfg.read_bytes()
    if reuse:shutil.copy2(prior/'field_surfaces.csv',o/'field_surfaces.csv')
    else:subprocess.run([str(binary),str(d/'worldtube'),str(cfg),str(o/'field_surfaces.csv')],check=True)
    fields=np.genfromtxt(o/'field_surfaces.csv',delimiter=',',names=True,dtype=LD)
    z0,entry,entry_inside=particle_plane(d,labels[0],half)
    assert len(entry)==plan['input_particles'] and len(np.unique(entry['source_id']))==len(entry)
    order=np.argsort(entry['source_id']);reference_ids=entry['source_id'][order];reference_K=entry['kinetic_lab_J'][order].astype(LD)
    index=lambda ids:np.searchsorted(reference_ids,ids)
    def locate(ids):
        k=index(ids);assert np.all(k<len(reference_ids)) and np.array_equal(reference_ids[k],ids);return k
    entry_index=locate(entry['source_id']);entry_sign=entry['crossing_direction'].astype(np.int64)
    assert np.all(entry_sign==1)
    kin=np.sum(entry['kinetic_lab_J'][entry_inside].astype(LD)*entry_sign[entry_inside],dtype=LD)
    entry_tau=entry['time_s'].astype(LD)-z0/tube.c-tube.tau_ref
    tolerance=LD('2e-21');assert np.min(entry_tau)>=lo-tolerance and np.max(entry_tau)<=hi+tolerance
    side_parts=[]
    for rank in range(16):
        path=d/'worldtube'/f'particle-sides-rank-{rank:05d}.csv'
        a=np.genfromtxt(path,delimiter=',',names=True,dtype=None,encoding='utf8')
        if a.size:side_parts.append(np.atleast_1d(a))
    side=np.concatenate(side_parts);side_index=locate(side['source_id']);side_direction=side['outward_direction'].astype(np.int64)
    assert np.all(abs(side_direction)==1)
    relevant=(side['z_lab_m']>=zs[0])&(side['z_lab_m']<=zs[-1])
    tau=side['time_lab_s'].astype(LD)-side['z_lab_m'].astype(LD)/tube.c-tube.tau_ref
    assert np.all((tau[relevant]>=lo-tolerance)&(tau[relevant]<=hi+tolerance))
    rows=[];coverage=[]
    for j,(name,z) in enumerate(zip(labels,zs)):
        if j==0:a=entry;inside=entry_inside;actual=z0
        else:actual,a,inside=particle_plane(d,name,half)
        assert abs(actual-z)<LD('1e-12') and np.all(a['crossing_direction']==1)
        ai=locate(a['source_id']);pt=a['time_s'].astype(LD)-actual/tube.c-tube.tau_ref
        in_gate=(pt>=lo-tolerance)&(pt<=hi+tolerance)
        outside_gate=int(np.count_nonzero(inside&~in_gate))
        inside=inside&in_gate
        side_use=relevant&(side['z_lab_m']<=z)
        counts=np.zeros(len(entry),dtype=np.int16)
        np.add.at(counts,entry_index[entry_inside],-entry_sign[entry_inside])
        np.add.at(counts,ai[inside],a['crossing_direction'][inside])
        np.add.at(counts,side_index[side_use],side_direction[side_use])
        unmatched=int(np.count_nonzero(counts))
        # The first five tubes have no occupied particle caps. The run ends
        # before the entire beam crosses p6: retain its missing late-cap term.
        assert not unmatched or (j==len(zs)-1 and np.all((counts==0)|(counts==-1))), (case,name,unmatched)
        kout=np.sum(a['kinetic_lab_J'][inside].astype(LD),dtype=LD)
        lateral=np.sum(side['kinetic_lab_J'][side_use].astype(LD)*side_direction[side_use],dtype=LD)
        # With zero signed crossing count for each ID, subtracting its fixed
        # entry energy from every event is an exact identity and avoids loss
        # of precision in the ~15mJ incoming/outgoing totals.
        stable=np.sum((a['kinetic_lab_J'][inside].astype(LD)-reference_K[ai[inside]]),dtype=LD)
        stable+=np.sum((side['kinetic_lab_J'][side_use].astype(LD)-reference_K[side_index[side_use]])*side_direction[side_use],dtype=LD)
        # Entry differences are identically zero, including outside-aperture
        # particles which enter later through a side.
        direct=kout-kin+lateral
        baseline_count_term=np.sum(reference_K*counts,dtype=LD)
        assert abs(stable+baseline_count_term-direct)<LD('1e-17')
        f=fields[j];terms={k:LD(f[k]) for k in fields.dtype.names if k.endswith('_J')}
        gain=sum(terms.values());ports=terms['field_exit_signed_J']+terms['field_entrance_signed_J']
        caps=terms['time_cap_late_signed_J']+terms['time_cap_early_signed_J']
        transverse=sum(terms[k] for k in terms if k.startswith('transverse_'))
        row=dict(z_lab_m=z,particle_entrance_signed_J=-kin,particle_exit_signed_J=kout,
            particle_ports_net_J=kout-kin,particle_lateral_signed_J=lateral,
            particle_net_transport_J=(stable if not unmatched else direct),
            particle_energy_change_with_entry_baseline_J=stable,
            particle_time_caps_J=(LD(0) if not unmatched else None),**terms,
            field_ports_net_J=ports,field_time_caps_net_J=caps,field_transverse_net_J=transverse,
            independent_field_gain_J=gain,
            physical_lab_sum_J=(stable+gain if not unmatched else None),
            measured_open_sum_excluding_particle_caps_J=direct+gain,
            measured_particle_cap_count=unmatched,particle_caps_known=(unmatched==0),
            incomplete_without_particle_sides_J=kout-kin+gain,
            maximum_coverage_relative_error=LD(f['max_coverage_relative_error']))
        rows.append(row)
        coverage.append(dict(plane=name,full_crossings=len(a),aperture_crossings=int(inside.sum()),
            side_events=int(side_use.sum()),unmatched_particle_ID_counts=unmatched,
            signed_ID_counts_sum=int(counts.sum()),aperture_port_events_outside_field_gate=outside_gate,
            direct_vs_baseline_corrected_sum_difference_J=float(direct-stable-baseline_count_term),
            occupied_late_particle_cap=bool(unmatched),particle_cap_energy_measured=not bool(unmatched)))
        print('LAB_BALANCE',case,name,float(direct),float(gain),
            float(stable+gain) if not unmatched else 'UNKNOWN_PARTICLE_CAP',unmatched,flush=True)
    table(o/'lab_energy_absolute.csv',rows)
    valid=[r for r in rows if r['particle_caps_known']]
    D=max(abs(r['particle_net_transport_J']) for r in valid);assert D>0
    relative=[{k:(100*v/D if k.endswith('_J') and v is not None else v) for k,v in row.items()} for row in rows]
    table(o/'lab_energy_percent.csv',relative)
    summary=dict(analysis_complete=True,all_six_physical_budgets_complete=False,
        frame='laboratory fixed-space worldtube with common retarded-time caps',
        common_relative_tau_s=[float(lo),float(hi)],tau_reference_s=float(tube.tau_ref),
        particle_time_caps_empty_through_plane='p5',particle_side_transport_measured=True,
        final_particle_cap_energy_J=None,final_particle_cap_energy_missing=True,
        external_electric_work_J=0.,source_exclusion_bound_m=float(source),observer_half_width_m=float(half),
        coverage=coverage,normalization_denominator_J=float(D),normalization_definition='max absolute full-gate net particle transport through fully measured p0-p5 tubes, including measured sides',
        final_components_J={k:(float(v) if v is not None else None) for k,v in rows[-1].items() if k.endswith('_J')},
        last_complete_plane='p5',last_complete_residual_J=float(valid[-1]['physical_lab_sum_J']),
        last_complete_residual_percent=float(100*valid[-1]['physical_lab_sum_J']/D),
        maximum_abs_complete_residual_percent=float(max(abs(r['physical_lab_sum_J']) for r in valid)*100/D),
        maximum_surface_coverage_error=float(max(r['maximum_coverage_relative_error'] for r in rows)),
        numerical_Cowan_terms_added=False,residual_derived_missing_terms=False,physical_acceptance=False,
        source_scope='Original correlated phase-clone SPARC beam; valid software/energy diagnostic, not independent-random FEL benchmark.')
    save(o/'LAB_RESULT.json',summary)
    return rows,relative,summary


def plots(out,cases):
    titles=['Energy matching + three-point kernel','Energy matching + momentum matching']
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for col,((rows,relative,summary),title) in enumerate(zip(cases,titles)):
        valid=[r for r in relative if r['particle_caps_known']];z=[r['z_lab_m'] for r in valid]
        for key,label,style in [('particle_net_transport_J','Particle net incl. sides','-'),
            ('field_ports_net_J','Field exit - entrance','--'),('field_transverse_net_J','Four field sides','-.'),
            ('field_time_caps_net_J','Late - early field cap',':'),('physical_lab_sum_J','Full physical sum','-')]:
            axes[0,col].plot(z,[r[key] for r in valid],style,marker='o',label=label)
        axes[0,col].set(title=title,ylabel='% of fixed net-particle scale');axes[0,col].legend(fontsize=8)
        axes[1,col].plot(z,[r['physical_lab_sum_J'] for r in valid],'o-',label='Complete physical sum')
        axes[1,col].plot(z,[r['incomplete_without_particle_sides_J'] for r in valid],'s--',label='If particle sides omitted')
        axes[1,col].set(xlabel='Laboratory exit plane z (m)',ylabel='% of same fixed scale');axes[1,col].legend(fontsize=8)
    for ax in axes.flat:
        ax.axhline(0,color='grey',lw=.6);ax.grid(alpha=.2);ax.set_xlim(0,15.3)
        ax.axvspan(12.8,15.3,color='grey',alpha=.12);ax.text(.91,.93,'p6: particle\ncap unmeasured',transform=ax.transAxes,ha='center',va='top',fontsize=8)
    fig.suptitle('Laboratory physical energy budget | all frequencies, full common field gate\nPorts + four sides + null-time caps; complete particle coverage through p5 only')
    fig.savefig(out/'lab_energy.png',dpi=175);fig.savefig(out/'lab_energy.pdf');plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    for col,((rows,_,_),title) in enumerate(zip(cases,titles)):
        z=[r['z_lab_m'] for r in rows]
        groups=[['field_entrance_signed_J','field_exit_signed_J'],
            ['time_cap_early_signed_J','time_cap_late_signed_J','transverse_xminus_signed_J','transverse_xplus_signed_J','transverse_yminus_signed_J','transverse_yplus_signed_J'],
            ['particle_ports_net_J','particle_lateral_signed_J','particle_net_transport_J','physical_lab_sum_J']]
        for row,keys in enumerate(groups):
            for key in keys:
                use=[r for r in rows if r[key] is not None and (row<2 or r['particle_caps_known'])]
                axes[row,col].plot([r['z_lab_m'] for r in use],[r[key]*1e6 for r in use],marker='o',label=key.removesuffix('_J'))
            axes[row,col].set(xlabel='Laboratory z (m)',ylabel='Signed energy (µJ)');axes[row,col].legend(fontsize=7)
            axes[row,col].axhline(0,color='grey',lw=.6);axes[row,col].grid(alpha=.2)
        axes[0,col].set_title(title)
    fig.suptitle('Individual measured laboratory channels | outward signs retained')
    fig.savefig(out/'lab_surface_terms.png',dpi=170);fig.savefig(out/'lab_surface_terms.pdf');plt.close(fig)


def main(parent,out):
    parent,out=Path(parent),Path(out);out.mkdir(exist_ok=True);binary=out/'worldtube_integrator'
    subprocess.run(['c++','-std=c++17','-O3','-fopenmp',str(Path(__file__).with_name('worldtube.cpp')),'-o',str(binary)],check=True)
    geometry_test(out,binary)
    answers=[]
    for case in ['energy-binomial3','energy-momentum']:answers.append(analyze_case(parent,out,case,binary))
    plots(out,answers)
    save(out/'LAB_SUMMARY.json',dict(analysis_complete=True,all_six_physical_budgets_complete=False,
        cases={name:x[2] for name,x in zip(['energy-binomial3','energy-momentum'],answers)}))


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
