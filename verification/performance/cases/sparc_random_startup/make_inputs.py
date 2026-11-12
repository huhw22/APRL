"""Sample independent physical events on Slurm and prepare exact common inputs."""
import copy
import math
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

from common import W, frozen, load, native, require, save, sha
from source_math import C, E, DTYPE, plan_bins, sample_records, validate_config
from source_verify import analytic_phase_space, recover_parameter_plane, distribution_report
from phase_map import phase_coordinates, gamma_from_u
from binding import binding


def choose_rows(n, count, seed):
    return np.sort(np.random.Generator(np.random.Philox(seed)).choice(n,count,replace=False))


def generate(root, plan):
    for seed in sorted({x['physical_seed'] for x in plan['cases']}):
        cfg=copy.deepcopy(plan['source_config']);cfg['seed']=seed
        cfg=validate_config(cfg);bins=plan_bins(cfg);n=cfg['electrons']
        entries=[x for x in plan['cases'] if x['physical_seed']==seed]
        count={x['name']:round(n/x['target_weight']) for x in entries}
        largest=max(count,key=count.get)
        rows={largest:choose_rows(n,count[largest],seed+1000)}
        for name in count:
            if name!=largest:
                rows[name]=np.sort(rows[largest][choose_rows(count[largest],count[name],seed+2000)])
        arrays={name:np.empty(size,dtype=DTYPE) for name,size in count.items()}
        cursor=0
        for j,nbin in enumerate(bins['counts']):
            end=cursor+int(nbin)
            for start in range(cursor,end,262144):
                stop=min(end,start+262144)
                block,_=sample_records(cfg,start,stop-start,bins['cdf_values'][j],bins['cdf_values'][j+1])
                for name,index in rows.items():
                    lo,hi=np.searchsorted(index,[start,stop])
                    if hi>lo: arrays[name][lo:hi]=block[index[lo:hi]-start]
            cursor=end
            if j%32==0: print('independent source',seed,'bin',j,'/',len(bins['counts']),flush=True)
        require(cursor==n,'Conceptual source truncated')
        prefix=np.r_[0,np.cumsum(bins['counts'])]
        for entry in entries:
            a=arrays[entry['name']]
            # v4 weights are relative. Integer unit values sum exactly in the
            # existing reader; beam.input.electrons supplies the physical N.
            a['macro_weight']=1.
            source_bins=np.searchsorted(prefix,a['source_id']-1,side='right')-1
            require(len(np.unique(a['plane_position_m'],axis=0))==len(a),'Copied transverse positions')
            v,t,_=recover_parameter_plane(a,cfg['phase_space'])
            values=np.column_stack([v,t/np.longdouble(cfg['profile']['full_duration_s'])])
            values-=analytic_phase_space(cfg)[0]
            recovered=np.searchsorted(bins['edges'].astype(np.longdouble),t,side='right')-1
            require(np.array_equal(recovered,source_bins),'Generation-bin recovery differs')
            counts=np.bincount(source_bins,minlength=len(bins['counts']))
            statistics=distribution_report(cfg,len(a),values.sum(axis=0),values.T@values,
                                          counts,counts,bins['edges'])
            path=root/'inputs'/(entry['name']+'.h5')
            with h5py.File(path,'x') as f:
                g=f.create_group('particles');g.attrs.update(format_version=4,input_plane_z_m=0.,
                    electron_count=n,macro_records=len(a),scope='Independent physical source uniformly thinned; no clones',
                    macro_weight=n/len(a),input_relative_weight=1.,physical_seed=seed)
                g.create_dataset('records',data=a,chunks=(min(65536,len(a)),))
                g=f.create_group('generation');g['complete']=np.uint8(1)
                g['committed_records']=np.uint64(len(a));g['physical_bin_counts']=bins['counts']
                g['selected_bin_counts']=counts;g['bin_edges_time_s']=bins['edges']
                g['config_json']=__import__('json').dumps(cfg)
            info=dict(input=str(path),input_sha256=sha(path),macro_records=len(a),
                physical_electrons=n,actual_weight=n/len(a),physical_seed=seed,
                transverse_unique_positions=len(a),zero_transverse_duplicates=True,
                source_ids_unchanged=True,input_relative_weights_all_one=True,phase_space_statistics=statistics,
                noise_expectation='sum(w^2)/sum(w)^2=1/Nmacro; not physical one-electron shot noise',
                profile_warning='Central231-wavelength hard crop of the original full parabola; core77 wavelengths reported separately.')
            save(Path(entry['directory'])/'INPUT_VALIDATED.json',info)
            prepare_case(root,plan,entry,a,prefix)


def geometry(a,plan):
    G=plan['boost_gamma'];bb=math.sqrt(1-1/G**2)
    dz=G*(1+bb)*plan['wavelength_m']/20;dt=dz/C
    u=a['proper_velocity'];gamma=np.sqrt(1+(u*u).sum(axis=1));beta=u[:,2]/gamma
    w=a['macro_weight'].astype(np.longdouble);t=a['arrival_time_offset_s'].astype(np.longdouble)
    mb=(w*beta).sum()/w.sum();mbt=(w*beta*t).sum()/w.sum()
    light=(np.hypot(1.,u[:,2])+u[:,2])/(G+math.sqrt(G-1)*math.sqrt(G+1))
    hu=.5*(light-1/light);held=hu/np.hypot(1.,hu)
    center=float(C*(mb*t.max()-mbt)+2*dz/G)
    ref=float((np.longdouble(center)/C+mbt)/mb)
    snap=beta*C*(ref-a['arrival_time_offset_s']);inp=-snap/(C*beta)
    origin=float(np.min(inp-2*dz/(G*C*(bb+held))))
    z=-G*C*(bb+held)*(inp-origin);translation=-.5*float(z.min()+z.max())
    last=.01+77*.028
    span=float(t.max()-t.min())+77*plan['wavelength_m']/C+(1-bb)*last/C+4e-15
    steps=math.ceil(span*G/dt);active=math.ceil(last/(G*dz))+64+128+2
    extent=max((1+held.max())*(steps+2)+66-translation/dz,translation/dz-64+active+2)
    full=math.ceil(2*extent/256)*256
    stop=G*G*(center+bb*C*origin)+math.sqrt(G*G-1)*dz*steps
    require(held.min()>0,'Held motion must be forward')
    return dict(center=center,reference_time_s=ref,origin_time_s=origin,translation_m=translation,
        dz=dz,dt=dt,steps=steps,active_nz=active,reference_nz=full,stop=stop,
        input_span_s=float(t.max()-t.min()),grid_entry_fraction=float((translation/dz+full/2)%1))


def prepare_case(root,plan,entry,a,prefix):
    d=Path(entry['directory']);info=load(d/'INPUT_VALIDATED.json');geo=geometry(a,plan)
    cfg=yaml.safe_load((Path(plan['prior_campaign'])/'energy-binomial3/run.yaml').read_text())
    cfg['output']['manifest']=str(d/'run-manifest.yaml')
    cfg['mesh'].update(cells=[plan['transverse_cells'],plan['transverse_cells'],geo['reference_nz']],duration=(geo['steps']+.5)*geo['dt'])
    require(float(np.abs(a['plane_position_m']).max())+2*cfg['mesh']['cell_size'][0]
            < (plan['transverse_cells']/2-6)*cfg['mesh']['cell_size'][0],
            'Random transverse tails need a larger physical mesh; do not clip particles')
    cfg['mesh']['moving_window'].update(cells=geo['active_nz'],reference_directory=str(d/'reference'),particle_collar_cells=8)
    cfg['beam']['reference']['initial_center_z']=geo['center']
    cfg['beam']['input'].update(file=info['input'],electrons=info['physical_electrons'])
    cfg['sources']['magnetic_elements']=cfg['sources']['magnetic_elements'][:1]
    cfg['detectors'].update(directory=str(d/'detectors'),substep_particle_planes=True,
        particle_planes=[dict(name=n+'-particles',z=z,buffer_records=32768,compression=0) for n,z in plan['particle_planes']],
        field_planes=[dict(name=n,z=z,rhythm=geo['dt']/plan['boost_gamma'],buffer_samples=8,compression=1,
                          particle_background={'enabled':False}) for n,z in plan['particle_planes'] if n.startswith('u') and n!='und-entry'],
        energy_planes=[])
    cfg['energy_ledger'].update(directory=str(d/'ledger'),sample_interval_steps=10,discrete_operator_audit=True)
    cfg['runtime']['resource_monitor']['progress_interval_steps']=50
    cfg['stop']=dict(mode='reference-center-z',z=geo['stop'])
    (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
    save(d/'GEOMETRY.json',geo)
    save(d/'COHORT.json',dict(core_source_id_first=int(prefix[77])+1,
        core_source_id_end_exclusive=int(prefix[154])+1,bin_prefix=prefix.tolist(),
        generation_bins=[77,154],source_ids_stable=True))
    native_path=d/'common.par.h5'
    integer,theta,rounded=phase_coordinates(a['arrival_time_offset_s'],plan['wavelength_m'])
    lo=int(integer.min())-256;ns=math.ceil((int(integer.max())-lo+1+256)/32)*32
    sort=np.argsort(integer,kind='stable');sorted_bins=integer[sort]
    g=gamma_from_u(a['proper_velocity']);index=[]
    with h5py.File(native_path,'x') as f:
        for key,value in [('refposition',lo*plan['wavelength_m']),('slicelength',plan['wavelength_m']),
            ('slicespacing',plan['wavelength_m']),('slicecount',ns),('beamletsize',1),('one4one',0)]:
            f[key]=[value]
        for j in range(ns):
            left,right=np.searchsorted(sorted_bins,[lo+j,lo+j+1]);rows=sort[left:right]
            group=f.create_group(f'slice{j+1:06d}')
            group['current']=[len(rows)*info['actual_weight']*E*C/plan['wavelength_m']]
            # A single zero-current record makes empty non-one4one slices safe.
            for key,values in [('x',a['plane_position_m'][:,0]),('y',a['plane_position_m'][:,1]),
                ('px',a['proper_velocity'][:,0]),('py',a['proper_velocity'][:,1]),('gamma',g),('theta',theta)]:
                group[key]=values[rows] if len(rows) else [plan['source_config']['phase_space']['gamma_mean'] if key=='gamma' else 0.]
            group['source_id']=a['source_id'][rows] if len(rows) else np.array([0],dtype='u8')
            if len(rows): index.append(j+1)
    save(d/'NATIVE_MAPPING.json',dict(native=str(native_path),native_sha256=sha(native_path),
        active_groups=index,slices=ns,integer_origin=lo,rounded_phase_values=rounded,
        one4one=False,nbins=1,sorting_enabled=False,shotnoise=False,
        scope='Exact actual event phases and 6D states; variable slice populations, constant macro-weight encoded in slice current. No re-uniformization or cloning.'))
    # Same quarter-period averaged end ramps as the accepted first segment.
    aw=plan['K']/math.sqrt(2);definitions=['PREP: DRIFT={l=.01};','M: MARKER={dumpfield=1};']
    ramps=[]
    for j in range(4):
        name=f'R{j}';ramps.append(name)
        definitions.append(f'{name}: UNDULATOR={{lambdau=.028,nwig=.25,aw={aw*.5*(1-math.cos(math.pi*(j+.5)/4)):.17g}}};')
    line=['PREP',*ramps]
    for label,nwig in [('A',6),('B',12),('C',18),('D',20),('E',19)]:
        definitions.append(f'{label}: UNDULATOR={{lambdau=.028,nwig={nwig},aw={aw:.17g}}};')
        line.append(label)
        if label!='E':line.append('M')
    line.extend([*ramps[::-1],'M'])
    definitions.append('FEL: LINE={'+','.join(line)+'};')
    for mode in ('off','longrange'):
        gd=d/('genesis-'+mode)
        (gd/'common.par.h5').symlink_to(native_path)
        (gd/'reference.lat').write_text('\n'.join(definitions)+'\n')
        setup=f'''&setup
rootname=reference
lattice=reference.lat
beamline=FEL
lambda0={plan['wavelength_m']:.17g}
gamma0={plan['source_config']['phase_space']['gamma_mean']:.17g}
delz=.007
shotnoise=false
one4one=false
npart=1024
nbins=1
beam_global_stat=true
field_global_stat=true
&end
&importbeam
file=common.par.h5
time=true
&end
&write
beam=imported
&end
'''
        (gd/'import.in').write_text(setup)
        track='''&field
power=0
dgrid=.0012
ngrid=161
&end
'''
        if mode=='longrange':track+='''&efield
longrange=true
rmax=2
nz=1
nphi=0
ngrid=100
&end
'''
        track+='''&track
output_step=1
&end
&write
beam=exit
field=exit
&end
'''
        (gd/'reference.in').write_text(setup.replace('beam=imported','beam=initial')+track)
    files=[Path(info['input']),d/'run.yaml',d/'COHORT.json',d/'GEOMETRY.json',d/'NATIVE_MAPPING.json',native_path]
    files.extend(p for p in d.glob('genesis-*/*') if p.is_file())
    save(d/'FROZEN.json',{str(p):sha(p) for p in files})


def main(root):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm for generation/import checks')
    root=Path(root);frozen(root/'FROZEN_PREPARATION.json');plan=load(root/'PLAN.json')
    from run_case import launch, srun, verify_import
    env=os.environ.copy();env['OMP_NUM_THREADS']='1'
    launch(srun(32,1)+[str(W/'build/cluster-tools/binding_probe'),'32','1'],
           root/'preparation',env,120,'binding.log')
    save(root/'preparation/BINDING.json',binding(root/'preparation/binding.log',1,32,1))
    generate(root,plan)
    for entry in plan['cases']:
        d=Path(entry['directory']);gd=d/'genesis-off'
        launch(srun(32,1)+[plan['genesis_binary'],'import.in'],gd,env,300,'import.log')
        verify_import(d,gd/'imported.par.h5')
        save(d/'IMPORT_VERIFIED.json',dict(passed=True,native_import_preserved=True))
        frozen(d/'FROZEN.json')
    frozen(root/'FROZEN_PREPARATION.json')
    save(root/'INPUTS_COMPLETE.json',dict(complete=True,cases=[x['name'] for x in plan['cases']],
        physical_acceptance=False,job=os.environ['SLURM_JOB_ID']))


if __name__=='__main__': main(sys.argv[1])
