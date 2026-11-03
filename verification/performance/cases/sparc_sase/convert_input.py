"""Convert closed native GENESIS shot-noise dump to common crossing events."""
import argparse,copy,math
from pathlib import Path
import h5py,numpy as np,yaml
from common import *

DTYPE=np.dtype([('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),
 ('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')])

def main(root):
    root=Path(root);plan=load(root/'plan.json');currents=slice_profile();n=len(currents)
    records=root/'particles.h5';reference=root/'genesis/common.par.h5'
    scores={str(m):[] for m in [1,2,3]};quadrants=[[] for _ in range(4)]
    maxphase=0.;sumw=np.longdouble(0.);sumwb=np.longdouble(0.);sumwbt=np.longdouble(0.)
    tmin=1.;tmax=-1.;umin=1.;umax=-1.;xmax=0.;mean_u=np.zeros(3,dtype=np.longdouble)
    beta_b=np.sqrt(1-1/GBOOST**2)
    with h5py.File(root/'native-noise/initial.par.h5','r') as src,h5py.File(records,'x') as out,h5py.File(reference,'x') as gen:
        pad=512
        count=int(src['slicecount'][0]);need(count>=n,'Native loading truncated')
        for key,val in [('refposition',-pad*LAMBDA),('slicelength',LAMBDA),
         ('slicespacing',LAMBDA),('slicecount',count+2*pad),('beamletsize',NBINS),('one4one',0)]:
            gen[key]=np.array([val],dtype='i4' if key in ['slicecount','beamletsize','one4one'] else 'f8')
        group=out.create_group('particles');group.attrs['format_version']=4
        group.attrs['scope']='Native GENESIS physical shot noise, common lab-plane crossing events; finite-band SASE benchmark'
        group.attrs['input_plane_z_m']=0.
        dataset=group.create_dataset('records',shape=(n*NPART,),dtype=DTYPE,chunks=(65536,))
        template=None;cursor=0;cache=[]
        for i in range(count):
            old=src[f'slice{i+1:06d}'];new=gen.create_group(f'slice{i+pad+1:06d}')
            vals={key:old[key][:] for key in ['x','y','px','py','gamma','theta']}
            uz=np.sqrt(vals['gamma']**2-1-vals['px']**2-vals['py']**2)
            # The native dump specifies the published entrance state. Move it
            # upstream through an independently defined magnet-free drift.
            desired_arrival=(i+np.mod(-vals['theta'],2*np.pi)/(2*np.pi))*LAMBDA/c
            lag=(1+vals['px']**2+vals['py']**2)/(uz*(vals['gamma']+uz))
            reference_lag=1/np.sqrt(1-1/GAMMA**2)-1
            vals['x']=vals['x']-PREP*vals['px']/uz
            vals['y']=vals['y']-PREP*vals['py']/uz
            vals['theta']=vals['theta']+2*np.pi/LAMBDA*PREP*(lag-reference_lag)
            cur=float(currents[i]) if i<n else 0.
            new['current']=[cur]
            for key,val in vals.items():
                need(np.isfinite(val).all(),'Nonfinite native input')
                new[key]=val
            if template is None:template=new
            if i>=n:continue
            need(len(vals['x'])==NPART,'Unexpected native particle count')
            weight=cur*LAMBDA/(e*c*NPART);ne=weight*NPART
            theta=vals['theta'];arrival=desired_arrival-PREP/c*(lag-reference_lag)
            r=np.zeros(NPART,dtype=DTYPE)
            r['plane_position_m']=np.c_[vals['x'],vals['y']]
            r['arrival_time_offset_s']=arrival
            r['proper_velocity'][:,0]=vals['px'];r['proper_velocity'][:,1]=vals['py']
            r['proper_velocity'][:,2]=uz
            r['source_id']=np.arange(i*NPART+1,(i+1)*NPART+1);r['macro_weight']=weight
            cache.append(r)
            if len(cache)==64 or i==n-1:
                batch=np.concatenate(cache);dataset[cursor:cursor+len(batch)]=batch;cursor+=len(batch);cache=[]
            b=uz/vals['gamma']
            sumw+=np.longdouble(weight)*NPART;sumwb+=np.longdouble(weight)*np.sum(b.astype(np.longdouble))
            sumwbt+=np.longdouble(weight)*np.sum(b.astype(np.longdouble)*arrival.astype(np.longdouble))
            mean_u+=np.longdouble(weight)*np.sum(r['proper_velocity'].astype(np.longdouble),axis=0)
            # Stable boost light-cone momenta, retaining transverse proper velocity.
            transverse=1+vals['px']**2+vals['py']**2
            light=(vals['gamma']+uz)/(GBOOST+np.sqrt(GBOOST**2-1))
            up=(light-transverse/light)/2;gp=(light+transverse/light)/2;held=up/gp
            umin=min(umin,float(held.min()));umax=max(umax,float(held.max()))
            tmin=min(tmin,float(arrival.min()));tmax=max(tmax,float(arrival.max()))
            xmax=max(xmax,float(np.abs(r['plane_position_m']).max()))
            reconstructed=np.exp(-2j*np.pi*c*arrival/LAMBDA)
            maxphase=max(maxphase,float(np.max(np.abs(reconstructed-np.exp(1j*theta)))))
            if cur>.2*IP:
                for m in [1,2,3]:scores[str(m)].append(float(ne*abs(np.exp(1j*m*theta).mean())**2))
                for q in range(4):
                    mask=(vals['x']>=0)==bool(q&1);mask&=(vals['y']>=0)==bool(q&2)
                    nq=int(mask.sum())
                    quadrants[q].append(float(ne*nq/NPART*abs(np.exp(1j*theta[mask]).mean())**2))
        for j in [*range(pad),*range(pad+count,pad*2+count)]:
            g=gen.create_group(f'slice{j+1:06d}');g['current']=[0.]
            for key in ['x','y','px','py','gamma','theta']:g[key]=template[key]
        need(cursor==n*NPART,'Incomplete crossing input')
    stats={m:dict(samples=len(v),mean=float(np.mean(v)),standard_error=float(np.std(v,ddof=1)/np.sqrt(len(v)))) for m,v in scores.items()}
    qstats=[float(np.mean(v)) for v in quadrants]
    need(all(abs(v['mean']-1)<.10 for v in stats.values()),'Physical harmonic noise variance failed')
    need(all(abs(v-1)<.15 for v in qstats),'Transverse shot noise normalization failed')
    need(maxphase<1e-10,'Optical phase changed during lab-plane conversion')
    need(umin>0,'Current injection model requires forward held motion')
    save(root/'NOISE_VALIDATION.json',dict(harmonics_Ne_times_abs_b_squared=stats,
      quadrants_mean_Ne_q_times_abs_b_q_squared=qstats,max_optical_phase_error=maxphase,
      input_particles=n*NPART,input_electrons=float(sumw),charge_C=float(sumw)*e,
      physical_charge_relative_error=float(sumw*e/Q-1),pass_noise_variance=True,
      method='Unchanged GENESIS Fawley loader, imported without adding a second noise realization',
      preparation_drift_m=PREP,
      entrance_mapping='Exact ballistic inverse through 1 cm, retaining the desired entrance arrival phases and transverse phase space in the magnet-only limit',
      scope='Physical noise in registered near-resonant band and harmonics 1-3; not all-frequency point-electron noise'))
    mean_beta=sumwb/sumw;mean_bt=sumwbt/sumw
    dx,dy=DX,DY
    dz=GBOOST*(1+beta_b)*LAMBDA/20
    dt=dz/c
    center=float(c*(mean_beta*np.longdouble(tmax)-mean_bt)+2*dz/GBOOST)
    ref=float((np.longdouble(center)/c+mean_bt)/mean_beta)
    t0=float(tmin-ref-2*dz/(GBOOST*c*(beta_b+umin)))
    z0max=-GBOOST*c*(beta_b+umax)*(tmin-ref-t0)
    z0min=-GBOOST*c*(beta_b+umax)*(tmax-ref-t0)
    translation=-.5*(z0min+z0max)
    lattice=load(root/'genesis/lattice.json');last=lattice['field_far_z_m']
    span=tmax-tmin+6*77*LAMBDA/c+(1-beta_b)*last/c+.3/(2*GAMMA**2*c)+4e-15
    steps=int(np.ceil(span*GBOOST/dt));duration=(steps+.5)*dt
    active=int(np.ceil(last/(GBOOST*dz)))+64+128+2
    extent=max((1+umax)*(steps+2)+66-translation/dz,translation/dz-64+active+2)
    full=int(np.ceil(2*extent/256))*256
    stop=GBOOST**2*(center+beta_b*c*t0)+np.sqrt(GBOOST**2-1)*dz*steps
    base=yaml.safe_load((WS/'results/source_coupling_full_comparison_20261101/free_drift_0p5m/run.yaml').read_text())
    planes=[(f'p{i+1}',PREP+i*(77*LU+.4)+77*LU) for i in range(6)]+[('far',last)]
    cases=[]
    for name,coupling,shape in [('boris','boris','linear'),('energy','conservative-source','linear'),('momentum','conservative-source','quintic-compensated')]:
        d=root/name;d.mkdir()
        cfg=copy.deepcopy(base)
        cfg['output']['manifest']=str(d/'run-manifest.yaml')
        cfg['mesh'].update(cells=[NX,NY,full],cell_size=[dx,dy,dz],center=[0.,0.,0.],
          duration=duration,boost_gamma=GBOOST,particle_coupling=coupling,longitudinal_shape=shape,
          particle_steps_per_undulator_period=80)
        cfg['mesh']['moving_window'].update(cells=active,reference_directory=str(d/'reference'))
        cfg['beam']['reference'].update(input_plane_z=0.,initial_center_z=center,activation_plane_z=0.)
        cfg['beam']['input'].update(file=str(records),electrons=float(sumw))
        cfg['sources']['magnetic_elements']=lattice['magnets']
        cfg['detectors'].update(directory=str(d/'detectors'),
          field_planes=[dict(name=n,z=z,rhythm=dt/GBOOST,buffer_samples=8,compression=1,
                             particle_background={'enabled':False}) for n,z in planes],
          particle_planes=[dict(name=n+'-particles',z=z,buffer_records=32768,compression=0) for n,z in [('entry',0.),*planes]],
          energy_planes=[])
        cfg['radiation_resolution']['maximum_photon_energy_eV']=3.0
        cfg['energy_ledger'].update(enabled=True,directory=str(d/'ledger'),sample_interval_steps=10,discrete_operator_audit=True)
        cfg['trajectory']['enabled']=False
        cfg['stop']={'mode':'reference-center-z','z':stop}
        (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
        cases.append(dict(name=name,directory=str(d)))
    # Preflight preserves full beam, full reference and active field geometry.
    # Only stop time/output path changes, so no field-sized optimistic proxy.
    d=root/'timing';d.mkdir()
    cfg=yaml.safe_load((root/'boris/run.yaml').read_text())
    def paths(x):
        if isinstance(x,dict):return {k:paths(v) for k,v in x.items()}
        if isinstance(x,list):return [paths(v) for v in x]
        return x.replace(str(root/'boris'),str(d)) if isinstance(x,str) else x
    cfg=paths(cfg);cfg['mesh']['duration']=120.5*dt
    cfg['stop']={'mode':'reference-center-z','z':stop}
    cfg['runtime']['resource_monitor']['progress_interval_steps']=20
    (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
    plan.update(cases=cases,geometry=dict(nx=NX,ny=NY,reference_nz=full,active_nz=active,
      dx_m=dx,dy_m=dy,dz_m=dz,dt_s=dt,steps=steps,stop_reference_z_m=stop,
      reference_time_offset_s=ref,frame_time_origin_s=t0,held_beta_z_range=[umin,umax],
      field_planes=planes,uncompressed_EB_upper_GiB_per_case=len(planes)*(steps+1)*NX*NY*48/2**30,
      physical_half_width_m=(NX/2-10)*dx),
      input_electrons=float(sumw),input_particles=n*NPART,max_transverse_coordinate_m=xmax)
    save(root/'plan.json',plan)
    files=[records,reference,root/'plan.json',*list((root/'genesis').glob('*.lat')),
          *[Path(q['directory'])/'run.yaml' for q in cases],d/'run.yaml']
    save(root/'FROZEN_CARDS.json',{str(p):sha(p) for p in files})
    print(json.dumps(plan['geometry'],indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);main(ap.parse_args().root)
