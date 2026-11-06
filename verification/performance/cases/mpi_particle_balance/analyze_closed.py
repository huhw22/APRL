"""Closed scientific output comparisons for changed integer MPI ownership."""
import json
from pathlib import Path

import h5py
import numpy as np
import yaml

from common import require
from original_closed import verify,bitwise
import report_discrete


def verify_case(case):
    if case.get('bounded_duration_prefix') is True:
        from prefix_closed import verify_prefix_case
        return verify_prefix_case(case)
    d=Path(case['directory']);require((d/'PIC_COMPLETE.json').exists(),'No closed PIC result')
    card=yaml.safe_load((d/'run.yaml').read_text());cfg=card['detectors']
    expected={x['name']+'.h5' for x in cfg['field_planes']}|{x['name']+'-energy.h5' for x in cfg.get('energy_planes',[])}
    expected|={x['name']+'.h5' for x in cfg['particle_planes']}
    expected|={x['name']+'-substeps.h5' for x in cfg['particle_planes']}
    require({p.name for p in (d/'detectors').glob('*.h5')}==expected,'Missing/extra detector outputs')
    scientific=verify(d)
    require_all=case.get('require_all_particles_through_planes',True)
    counts={}
    if not require_all:
        require(case.get('input_source_ids_contiguous') is True,
            'Partial-prefix validation needs an independently verified contiguous source-ID range')
        first,last=case['expected_source_id_first'],case['expected_source_id_last']
        require(last-first+1==case['expected_particles'],'Partial-prefix source-ID range disagrees with input count')
    for name in expected:
        if name.endswith('-energy.h5') or name in {x['name']+'.h5' for x in cfg['field_planes']}:continue
        with h5py.File(d/'detectors'/name,'r') as h:
            a=h[next(iter(h))]['records'][:]
            counts[name]=len(a)
            if require_all:
                require(len(a)==case['expected_particles'],'Missing complete source through particle plane')
            else:
                require(len(a)<=case['expected_particles'],'Partial plane exceeds complete source count')
                require(np.all((a['source_id']>=first)&(a['source_id']<=last)),
                    'Partial plane contains an ID outside the accepted source')
            require(len(np.unique(a['particle_id']))==len(a) and len(np.unique(a['source_id']))==len(a),'Repeated particle identity')
    rows,sums,ranks,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5')
    require(ranks==case['ranks'] and steps==case['expected_steps'],'Wrong rank count or physical duration')
    _,_,ledger=report_discrete.calculate(rows,sums,steps)
    histories=sorted((d/'reference').glob('source-coupling-rank-*.csv'))
    require(len(histories)==ranks,'Missing coupling history')
    max_force=max_work=0.
    for p in histories:
        require(p.with_suffix('.complete').read_text().split()==['1',str(steps)],'Unclosed coupling history')
        a=np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True))
        require(np.array_equal(a['step'],np.arange(1,steps+1)),'Incomplete coupling steps')
        require(all(np.isfinite(a[n]).all() for n in a.dtype.names),'Nonfinite coupling history')
        max_force=max(max_force,float(np.max(a['self_force_defect'])));max_work=max(max_work,float(np.max(a['relative_fixed_point_work_defect'])))
    require(max_force<1e-13 and max_work<1e-8,'Fixed-point iteration failed')
    log=(d/'aprl.log').read_text();events=log.count('[particle-balance] step=')
    require(not case.get('balance_enabled',case.get('mode') in ['on','balanced']) or events>0,'No actual repartition occurred')
    return dict(passed=True,program_integrity_pass=True,physical_acceptance=False,steps=steps,ranks=ranks,
        scientific_files=scientific,strict_ledger=ledger,max_force_defect=max_force,max_work_defect=max_work,
        balance_events=events,particle_plane_counts=counts,
        complete_source_through_all_planes=require_all,
        verification_scope=case.get('verification_scope','complete configured detector/source comparison'))


def load_ledger(directory):
    directory=Path(directory)
    if (directory/'PREFIX_COMPLETE.json').exists():
        case=json.loads((directory/'run-plan.json').read_text())
        require(case.get('bounded_duration_prefix') is True,'Unexpected prefix marker on a full-run card')
        marker=json.loads((directory/'PREFIX_COMPLETE.json').read_text())
        require(marker.get('prefix_complete') is True and marker.get('physical_complete') is False,
            'Invalid closed-prefix marker')
        from prefix_closed import load_prefix_ledger
        return load_prefix_ledger(directory,case)
    return report_discrete.load(directory/'ledger',directory/'ledger/energy-ledger.h5')


def numerical(first,second,expected_particles):
    first,second=Path(first),Path(second);metrics={}
    for p in sorted((first/'detectors').glob('*.h5')):
        with h5py.File(p,'r') as h,h5py.File(second/'detectors'/p.name,'r') as z:
            a=h[next(iter(h))];b=z[next(iter(z))];fields={}
            for key in ['fields','records','time_s','coordinates']:
                if key not in a:continue
                x=a[key][:];y=b[key][:]
                require(x.shape==y.shape and x.dtype==y.dtype,'Detector layout changed')
                if x.dtype.names:
                    particle='source_id' in x.dtype.names
                    if particle:x=np.sort(x,order='source_id');y=np.sort(y,order='source_id')
                    for n in x.dtype.names:
                        if n.endswith('_low'):continue
                        if x[n].dtype.kind not in 'fc':
                            require(np.array_equal(x[n],y[n]),'Changed integer/particle identity: '+n);continue
                        xx=x[n].astype(np.longdouble);yy=y[n].astype(np.longdouble)
                        if n+'_low' in x.dtype.names:xx=xx+x[n+'_low'];yy=yy+y[n+'_low']
                        difference=float(np.max(abs(xx-yy),initial=0));scale=float(np.max(abs(xx),initial=0))
                        if particle and n in ['weight','charge_C','mass_kg']:
                            require(np.array_equal(xx,yy),'Particle weights/species changed')
                        else:require(difference<=1e-10*max(scale,1e-30),'Scientific output differs: '+p.name+'/'+n)
                        fields[key+'/'+n]=dict(max_abs_difference=difference,baseline_scale=scale)
                else:
                    difference=float(np.max(abs(x-y),initial=0));scale=float(np.max(abs(x),initial=0))
                    require(difference<=1e-10*max(scale,1e-30),'Field/time output differs')
                    fields[key]=dict(max_abs_difference=difference,baseline_scale=scale)
            metrics[p.name]=fields
    r0,s0,_,steps=load_ledger(first)
    r1,s1,_,steps1=load_ledger(second)
    require(steps==steps1 and np.array_equal(r0['step'],r1['step']),'Ledger time axes differ')
    D=float(np.max(abs(np.asarray(r0['delta_particle_kinetic_box_J'],dtype=np.longdouble))))
    require(D>0,'This particle test needs a nonzero common particle-change scale')
    _,q0,l0=report_discrete.calculate(r0,s0,steps);_,q1,l1=report_discrete.calculate(r1,s1,steps)
    gaps={n:float(np.max(abs(q0[n]-q1[n]))/D) for n in ['discrete_total_residual_J','particle_grid_coupling_gap_J']}
    dk=float(np.max(abs(np.asarray(r0['delta_particle_kinetic_box_J'],dtype=np.longdouble)-np.asarray(r1['delta_particle_kinetic_box_J'],dtype=np.longdouble)))/D)
    require(max([dk]+list(gaps.values()))<=1e-8,'Ledger/particle changes differ on the common exchange scale')
    return dict(detectors=metrics,common_denominator_J=D,max_delta_K_difference_over_D=dk,
        max_ledger_difference_over_D=gaps,baseline_ledger=l0,candidate_ledger=l1)


def compare_cases(plan):
    groups={}
    for case in plan['cases']:groups.setdefault(case['family'],{})[case['mode']]=case
    result={}
    for name,g in groups.items():
        old=g.get('old');off=g.get('off',g.get('new-off'));on=g.get('on',g.get('new-on'))
        require(off and on,'Missing static/balanced comparison')
        item={}
        if old:item['off_science_bitwise']=bitwise(old['directory'],off['directory'])
        item['balanced_science']=numerical(off['directory'],on['directory'],on['expected_particles'])
        result[name]=item
    return result
