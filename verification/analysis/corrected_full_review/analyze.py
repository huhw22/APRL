#!/usr/bin/env python3
"""Read a closed full run: exact identity coverage, particles, ledger and FFT.

Original run outputs and its failed legacy completion record remain immutable.
All large reads run inside one compute allocation. Only an ID-assumption error
in the old checker is permitted; no other incomplete or failed run is accepted.
"""
import argparse,csv,hashlib,importlib.util,json,os,re,sys,time
from decimal import Decimal,localcontext
from pathlib import Path
import h5py,numpy as np,yaml
from scipy.constants import c,e,m_e,h

LD=np.longdouble
MC2=LD(m_e)*LD(c)**2/LD(e)
HERE=Path(__file__).resolve().parent

def save(path,data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:float(v))+'\n')

def sha(path):
    d=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):d.update(b)
    return d.hexdigest()

def load_module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def gfactor(u):return np.sqrt(1+np.sum(u*u,axis=1))

def add(s,key,a):s[key]=s.get(key,LD(0))+np.sum(a,dtype=LD)

def finite(a):
    for n in a.dtype.names:
        if a[n].dtype.kind in 'fc':assert np.isfinite(a[n]).all(),n

def plane_summary(s,cx,cy,anchor,wavelengths,z,records,source_ids):
    W=s['w'];mg=s['g']/W;vg=max(LD(0),s['g2']/W-mg*mg)
    ms=s['s']/W;vs=max(LD(0),s['s2']/W-ms*ms);cov=s['sg']/W-ms*mg
    m=s['loss_g']/W;v=max(LD(0),s['loss_g2']/W-m*m)
    mb=s['loss_box_g']/W;vb=max(LD(0),s['loss_box_g2']/W-mb*mb)
    return {'plane_z_m':z,'records':records,'source_ID_min':int(source_ids.min()),'source_ID_max':int(source_ids.max()),
      'unique_particle_IDs':records,'source_IDs_match_original_input':True,'all_values_finite':True,
      'represented_electrons':float(W),'mean_gamma':float(anchor+mg),'std_gamma':float(np.sqrt(vg)),
      'projected_std_keV':float(MC2*np.sqrt(vg)/1000),'relative_projected_spread':float(np.sqrt(vg)/(anchor+mg)),
      'linear_detrended_std_keV':float(MC2*np.sqrt(max(LD(0),vg-cov*cov/vs))/1000),
      'rms_x_um':float(np.sqrt(max(LD(0),s['x2']/W-(s['x']/W)**2))*1e6),
      'rms_y_um':float(np.sqrt(max(LD(0),s['y2']/W-(s['y']/W)**2))*1e6),
      'rms_arrival_length_um':float(np.sqrt(vs)*1e6),'physical_aperture_fraction':float(s['physical']/W),
      'mean_loss_gamma':float(m),'std_loss_gamma':float(np.sqrt(v)),
      'std_over_mean_loss':float(np.sqrt(v)/m) if m!=0 else None,
      'mean_loss_eV':float(m*MC2),'std_loss_eV':float(np.sqrt(v)*MC2),
      'mean_loss_box_eV':float(mb*MC2),'std_loss_box_eV':float(np.sqrt(vb)*MC2),
      'total_loss_lab_J':float(m*LD(m_e)*LD(c)**2*W),
      'bunching':{f'{lam*1e9:.9g}_nm':float(np.hypot(cx[j],cy[j])/W) for j,lam in enumerate(wavelengths)}}

def scan_particles(root,out,plan,allow_exit_subset=False):
    n=plan['records'];files=[root/'detectors'/f'{v}-particles.h5' for v in ('entry','exit')]
    with h5py.File(plan['input_file']) as f:
        d=f['particles/records'];assert len(d)==n
        original=d[:]
    assert np.all(original['macro_weight']==1)
    source_ids=original['source_id'];sorted_ids=np.sort(source_ids)
    assert np.all(sorted_ids[1:]>sorted_ids[:-1]);del sorted_ids
    entry_u=np.empty((n,3),dtype=np.float64);entry_q=np.empty(n);entry_m=np.empty(n)
    entry_t=np.empty(n);entry_xy=np.empty((n,2));loss=np.full(n,np.nan);summaries=[];size=131072
    beta_gamma=LD(np.sqrt(plan['geometry']['gamma_boost']**2-1))
    G=LD(plan['geometry']['gamma_boost']);beta=beta_gamma/G
    with h5py.File(files[0]) as f:
        first=f['particle_plane/records'][0]
        tref=LD(first['time_s'])-LD(first['position_m'][2])/LD(c)
        anchor=gfactor(np.asarray(first['proper_velocity'],dtype=LD)[None,:])[0]
    wavelengths=(13.5e-9,h*c/e/91.5)
    for which,path in enumerate(files):
        seen=np.zeros(n,dtype=bool);s={};cx=np.zeros(2,dtype=LD);cy=cx.copy();max_input_u=0.;max_input_xy=0.
        with h5py.File(path) as f:
            group=f['particle_plane'];data=group['records']
            count=len(data)
            assert int(group['complete'][()])==1 and count==int(group['committed_records'][()])
            assert count==n if which==0 or not allow_exit_subset else 0<count<=n
            z=float(group.attrs['plane_z_m'])
            for start in range(0,count,size):
                a=data[start:min(count,start+size)];finite(a)
                raw=a['particle_id'];assert np.all((raw>=1)&(raw<=n));idx=raw.astype(np.int64)-1
                assert np.unique(idx).size==len(idx) and not np.any(seen[idx]);seen[idx]=True
                assert np.array_equal(a['source_id'],source_ids[idx])
                assert np.all(a['crossing_direction']==1)
                assert np.max(abs(a['position_m'][:,2]-z))<1e-12
                w=-np.asarray(a['charge_C'],dtype=LD)/LD(e)
                assert np.all(abs(w-1)<1e-12) and np.all(a['mass_kg']>0)
                assert np.all(abs(a['mass_kg']/m_e-1)<1e-12)
                u=np.asarray(a['proper_velocity'],dtype=LD);g=gfactor(u)
                if which==0:
                    entry_u[idx]=a['proper_velocity'];entry_q[idx]=a['charge_C'];entry_m[idx]=a['mass_kg'];entry_t[idx]=a['time_s'];entry_xy[idx]=a['position_m'][:,:2]
                    dl=np.zeros(len(a),dtype=LD);db=dl.copy()
                    max_input_u=max(max_input_u,float(np.max(abs(a['proper_velocity']-original['proper_velocity'][idx]))))
                    max_input_xy=max(max_input_xy,float(np.max(abs(a['position_m'][:,:2]-original['plane_position_m'][idx,:2]))))
                else:
                    assert np.array_equal(a['charge_C'],entry_q[idx]) and np.array_equal(a['mass_kg'],entry_m[idx])
                    u0=np.asarray(entry_u[idx],dtype=LD);g0=gfactor(u0)
                    # Stable per-ID gamma(entry)-gamma(exit), avoiding subtraction
                    # of ~600 MeV scalar energies in ordinary double precision.
                    dl=np.sum((u0-u)*(u0+u),axis=1)/(g0+g)
                    db=G*dl-beta_gamma*(u0[:,2]-u[:,2])
                    loss[idx]=np.asarray(dl*MC2,dtype=float)
                ss=LD(c)*(np.asarray(a['time_s'],dtype=LD)-LD(z)/LD(c)-tref)
                gg=g-anchor
                for k,v in {'w':w,'g':w*gg,'g2':w*gg*gg,'s':w*ss,'s2':w*ss*ss,'sg':w*ss*gg,
                    'x':w*a['position_m'][:,0],'x2':w*a['position_m'][:,0]**2,
                    'y':w*a['position_m'][:,1],'y2':w*a['position_m'][:,1]**2,
                    'loss_g':w*dl,'loss_g2':w*dl*dl,'loss_box_g':w*db,'loss_box_g2':w*db*db,
                    'physical':w*((abs(a['position_m'][:,0])<=.00035)&(abs(a['position_m'][:,1])<=.00035))}.items():add(s,k,v)
                dt=np.asarray(np.asarray(a['time_s'],dtype=LD)-LD(z)/LD(c)-tref,dtype=float)
                for j,lam in enumerate(wavelengths):
                    phase=2*np.pi*c*dt/lam;cx[j]+=np.sum(w*np.cos(phase),dtype=LD);cy[j]+=np.sum(w*np.sin(phase),dtype=LD)
                if start//size%64==0:print(f'particles {which}: {min(count,start+size)}/{count}',flush=True)
        assert np.count_nonzero(seen)==count
        if which==0 or not allow_exit_subset:assert np.all(seen)
        summ=plane_summary(s,cx,cy,anchor,wavelengths,z,count,source_ids)
        if which==0:summ.update(entry_vs_input_max_abs_u=max_input_u,entry_vs_input_max_abs_xy_m=max_input_xy)
        if which==0:entry_accumulator=(s.copy(),cx.copy(),cy.copy());entry_all=summ.copy()
        summaries.append(summ)
    missing=np.flatnonzero(~seen);missing_totals={k:'0' for k in ['particle_kinetic_lab_J','particle_cpz_lab_J','particle_kinetic_box_J','particle_cpz_box_J']}
    if len(missing):
        ss0,cx0,cy0=entry_accumulator;u=np.asarray(entry_u[missing],dtype=LD);g=gfactor(u)
        w=-np.asarray(entry_q[missing],dtype=LD)/LD(e);gg=g-anchor
        xx,yy=entry_xy[missing,0],entry_xy[missing,1]
        ss=LD(c)*(np.asarray(entry_t[missing],dtype=LD)-LD(summaries[0]['plane_z_m'])/LD(c)-tref)
        for key,value in {'w':w,'g':w*gg,'g2':w*gg*gg,'s':w*ss,'s2':w*ss*ss,'sg':w*ss*gg,
          'x':w*xx,'x2':w*xx*xx,'y':w*yy,'y2':w*yy*yy,
          'physical':w*((abs(xx)<=.00035)&(abs(yy)<=.00035))}.items():ss0[key]-=np.sum(value,dtype=LD)
        dt=np.asarray(np.asarray(entry_t[missing],dtype=LD)-LD(summaries[0]['plane_z_m'])/LD(c)-tref,dtype=float)
        for j,lam in enumerate(wavelengths):
            phase=2*np.pi*c*dt/lam;cx0[j]-=np.sum(w*np.cos(phase),dtype=LD);cy0[j]-=np.sum(w*np.sin(phase),dtype=LD)
        summaries[0]=plane_summary(ss0,cx0,cy0,anchor,wavelengths,summaries[0]['plane_z_m'],n-len(missing),source_ids[seen])
        summaries[0]['cohort']='Entry subset whose IDs also reached the exit plane'
        summaries[1]['source_ID_min']=summaries[0]['source_ID_min'];summaries[1]['source_ID_max']=summaries[0]['source_ID_max']
        plus=(g+u[:,2])/(G+beta_gamma);transverse=1+np.sum(u[:,:2]**2,axis=1)
        uzb=(plus-transverse/plus)/2;gb=(plus+transverse/plus)/2
        masses=np.asarray(entry_m[missing],dtype=LD);unit=masses*LD(c)**2
        for key,value in {'particle_kinetic_lab_J':unit*np.sum(u*u,axis=1)/(g+1),
          'particle_cpz_lab_J':unit*u[:,2],'particle_kinetic_box_J':unit*(gb-1),'particle_cpz_box_J':unit*uzb}.items():missing_totals[key]=str(np.sum(value,dtype=LD))
        with (out/'missing_exit_particles.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['particle_id','source_id','entry_x_m','entry_y_m','entry_gamma','ballistic_exit_x_m','ballistic_exit_y_m'])
            distance=plan['field_plane_z_m']-entry_all['plane_z_m']
            for j,index in enumerate(missing):writer.writerow([int(index+1),int(source_ids[index]),entry_xy[index,0],entry_xy[index,1],str(g[j]),
              float(entry_xy[index,0]+distance*entry_u[index,0]/entry_u[index,2]),float(entry_xy[index,1]+distance*entry_u[index,1]/entry_u[index,2])])
    actual_loss=loss[seen];assert np.isfinite(actual_loss).all()
    quantiles=np.quantile(actual_loss,[0,.001,.01,.1,.5,.9,.99,.999,1])
    lo,hi=quantiles[1],quantiles[-2]
    if lo==hi:lo-=max(1e-9,abs(lo)*1e-12);hi+=max(1e-9,abs(hi)*1e-12)
    edges=np.linspace(lo,hi,241);counts,_=np.histogram(actual_loss,edges)
    hist={'edges_loss_eV':edges,'count':counts}
    np.savez_compressed(out/'particle_distribution.npz',**hist,quantiles_eV=quantiles)
    save(out/'particle_statistics.json',{'planes':summaries,'entry_all_input_cohort':entry_all,
      'available_counts':{'entry':n,'exit':int(seen.sum()),'missing_exit':len(missing)},
      'missing_entry_totals_J':missing_totals,'full_exit_coverage':len(missing)==0,
      'loss_quantile_probabilities':[0,.001,.01,.1,.5,.9,.99,.999,1],
      'loss_quantiles_eV':quantiles.tolist(),'histogram_outside_display_count':int(len(actual_loss)-counts.sum()),
      'definitions':{'loss':'gamma_entry-gamma_exit; positive is loss; exact particle_id plus original source_id',
       'weights':'abs(charge_C)/e, unit physical electrons here; do not multiply relative weight twice',
       'spread':'Charge-weighted population std; linear detrending is a fitted chirp subtraction, not slice spread',
       'bunching':'Magnitude of weighted Fourier sum of lab arrival times at fixed wavelength; retarded-time origin affects phase only'}})
    return summaries

def read_energy(root,out,n,coverage=None):
    planes=[];series={}
    for tag in ('entry','exit'):
        path=root/'detectors'/f'{tag}-transport-energy.h5'
        with h5py.File(path) as f:
            g=f['energy_plane'];a=g['records'][:];finite(a)
            assert int(g['complete'][()])==1 and int(g['committed_records'][()])==len(a)
            expected=n if coverage is None else coverage['available_counts'][tag]
            last=a[-1];assert int(last['particle_downstream_count'])==expected and int(last['particle_upstream_count'])==0
            sums={}
            with localcontext() as ctx:
                ctx.prec=70
                for key in ('particle_kinetic_lab_J','particle_cpz_lab_J','particle_kinetic_box_J','particle_cpz_box_J'):
                    sums[key]=str(Decimal.from_float(float(last[key]))+Decimal.from_float(float(last[key+'_low'])))
            planes.append({'name':tag,'records':len(a),'particle_count':expected,'totals_decimal_J':sums,
              'physical_field_energy_lab_J':float(last['physical_field_energy_lab_J']),
              'full_field_energy_lab_J':float(last['field_energy_lab_J']),
              'field_gap_count':int(last['field_gap_count']),'outside_field_window_steps':int(last['outside_field_window_steps']),
              'physical_field_power_last_W':float(last['physical_field_power_lab_W'])})
            assert planes[-1]['field_gap_count']==planes[-1]['outside_field_window_steps']==0
            series[tag+'_time_lab_s']=a['time_lab_s'];series[tag+'_physical_power_W']=a['physical_field_power_lab_W']
    with localcontext() as ctx:
        ctx.prec=70
        raw_changes={k:str(Decimal(planes[1]['totals_decimal_J'][k])-Decimal(planes[0]['totals_decimal_J'][k])) for k in planes[0]['totals_decimal_J']}
        matched_entry={k:str(Decimal(v)-Decimal('0' if coverage is None else coverage['missing_entry_totals_J'][k])) for k,v in planes[0]['totals_decimal_J'].items()}
        changes={k:str(Decimal(planes[1]['totals_decimal_J'][k])-Decimal(matched_entry[k])) for k in matched_entry}
    G=LD(yaml.safe_load((root/'run.yaml').read_text())['mesh']['boost_gamma']);gb=np.sqrt((G-1)*(G+1))
    box=LD(changes['particle_kinetic_box_J']);cpz=LD(changes['particle_cpz_box_J']);lab=LD(changes['particle_kinetic_lab_J'])
    matched_n=n if coverage is None else coverage['available_counts']['exit']
    data={'planes':planes,'exit_minus_entry_decimal_J':changes,'matched_entry_totals_decimal_J':matched_entry,
      'unmatched_all_entry_to_subset_exit_delta_J':raw_changes,'matched_cohort_count':matched_n,
      'missing_entry_correction_scope':'Subtract missing IDs entry raw four-momenta from scalar entry totals; no exit energies invented.' if coverage is not None else 'Full coverage; no correction',
      'lab_change_from_boost_kinetic_J':float(G*box),'lab_change_from_boost_momentum_J':float(gb*cpz),
      'Lorentz_identity_difference_J':float(lab-G*box-gb*cpz),
      'scalar_endpoint_mean_loss_eV':float(-lab/(LD(matched_n)*LD(e))),
      'physical_field_exit_minus_entry_J':planes[1]['physical_field_energy_lab_J']-planes[0]['physical_field_energy_lab_J'],
      'scope':'Two independent fixed lab planes. Their energy difference lacks side flux and in-volume storage for a closed spacetime control volume.'}
    save(out/'energy_planes.json',data);np.savez_compressed(out/'energy_plane_series.npz',**series)
    return data

def read_ledger(root,out,title='19098 | corrected full-input free drift, 0.5 m'):
    sys.path.insert(0,str(HERE));import report_discrete,relative_ledger
    rows,sums,ranks,steps=report_discrete.load(root/'ledger',root/'ledger/energy-ledger.h5')
    index,q,summary=report_discrete.calculate(rows,sums,steps)
    absval,pct,meta=relative_ledger.construct(rows,sums,index,q)
    ledout=out/'ledger';ledout.mkdir();relative_ledger.write_tables(ledout,index,absval,pct,meta)
    relative_ledger.plot(ledout,index,pct,meta,title)
    final={name:str(rows[-1][name]) for name in rows.dtype.names}
    info={'ranks':ranks,'steps':steps,'normalization':meta,'summary':summary,'final_native_ledger_record':final}
    save(out/'ledger_summary.json',info)
    return info

def spectrum(root,out,simulation_job=19098,completion_scope='Closed solver and full source identity checked in the independent review; legacy batch exit 1 came from a sparse source-ID range assumption.'):
    m=load_module('fft_kernel',HERE/'analyze_spectra.py');target=out/'spectrum'
    m.analyze(root,root/'detectors/exit-field.h5',target,1,8,256)
    p=target/'summary.json';s=json.loads(p.read_text())
    s['notes']=[t.replace('and the known MPI-interface diagnostic defect','; the MPI diagnostic stencil is repaired in this run') for t in s['notes']]
    s.update(simulation_job=simulation_job,slab_cells=None,sampling_comb_spacing_eV=None,
      completion_scope=completion_scope,
      moving_window=True,batch_exit_code=1,simulation_step_exit_code=0)
    save(p,s)
    checker=load_module('independent_fft_review',HERE/'independent_review.py')
    review=checker.review_case(target);save(out/'spectrum_independent_review.json',review)
    return s

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();assert os.environ.get('SLURM_JOB_ID'),'Large scans require a compute allocation'
    root=a.run.resolve();out=a.out.resolve();assert not out.exists();out.mkdir()
    plan=json.loads((root/'run-plan.json').read_text());prior=json.loads((root/'completion-check.json').read_text())
    assert prior['error']=='Particle IDs out of range' and prior['field']['complete'] and prior['field']['all_values_finite']
    log=(root/'aprl.log').read_text();assert all(t in log for t in ('phase=solve-complete','phase=mpi-finalize-complete','phase=detector-close-complete'))
    assert 'Simulation failure' not in log
    for name,value in json.loads((root/'frozen.json').read_text()).items():assert sha(Path(name))==value,name
    before={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [root/'completion-check.json',root/'aprl.log',*(root/'detectors').glob('*.h5'),root/'ledger/energy-ledger.h5']}
    print('Start complete input-row / original source-ID validation',flush=True)
    particles=scan_particles(root,out,plan);print('Particle statistics complete',flush=True)
    energy=read_energy(root,out,plan['records']);ledger=read_ledger(root,out)
    save(out/'CORE_RESULT.json',{'particles':particles,'energy_planes':energy,'ledger':ledger})
    print('Particle/energy/ledger results committed; starting FFT',flush=True)
    spec=spectrum(root,out)
    assert all([Path(p).stat().st_size,Path(p).stat().st_mtime_ns]==value for p,value in before.items())
    save(out/'COMPLETE.json',{'passed':True,'solver_complete':True,'original_batch_exit_code':1,
      'legacy_checker_error':'Assumed arbitrary source_id was sequential 1..N; corrected analysis checks unique particle_id and original SDDS source_id.',
      'particles_per_plane':plan['records'],'field_samples':spec['samples'],'source_stat_unchanged':before,
      'new_pic_simulations':0,'simulation_job':19098,'analysis_job':os.environ['SLURM_JOB_ID']})
    print('COMPLETE: closed full drift reviewed; original outputs unchanged',flush=True)

if __name__=='__main__':main()
