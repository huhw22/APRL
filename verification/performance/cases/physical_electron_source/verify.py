"""Independent streaming checks of a unit-electron crossing-event master.

The distribution checks reconstruct the parameter-plane events and compare
their moments with analytic expectations.  Genesis is checked by independently
recovering arrival times and looking up the immutable master by source ID.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

LD = np.longdouble
C = LD(299792458)
QE = LD('1.602176634e-19')
DTYPE = np.dtype([('plane_position_m','<f8',(2,)),
    ('arrival_time_offset_s','<f8'),('proper_velocity','<f8',(3,)),
    ('source_id','<u8'),('macro_weight','<f8')])


def need(value, message):
    if not value:
        raise RuntimeError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def identity(path, checksum=False):
    p=Path(path).resolve();s=p.stat()
    result=dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
    if checksum:result['sha256']=sha(p)
    return result


def scalar(dataset):
    value=np.asarray(dataset[()])
    need(value.size==1,'Expected one scalar: '+dataset.name)
    return value.reshape(-1)[0].item()


def read_config(f):
    if 'generation/config_json' not in f:return None
    text=f['generation/config_json'][()]
    if isinstance(text,bytes):text=text.decode()
    return json.loads(str(text))


def master_gate(f):
    need('generation/complete' in f and scalar(f['generation/complete'])==1,
         'Master is not complete; no particle data read')
    need('particles/records' in f,'Missing master records')
    group=f['particles'];data=group['records']
    need(int(group.attrs['format_version'])==4,'Expected APRL v4 crossing events')
    need(data.dtype==DTYPE,'Unexpected master record ABI')
    need(len(data)>0,'Empty master')
    need(int(group.attrs['electron_count'])==len(data),'Electron count differs from rows')
    need('generation/committed_records' in f and scalar(f['generation/committed_records'])==len(data),
         'Master committed-record count differs')
    need(np.isfinite(float(group.attrs['input_plane_z_m'])),'Invalid entrance plane')
    return data,read_config(f)


def recover_parameter_plane(records, phase_space):
    u=records['proper_velocity'].astype(LD)
    # Direct long-double norm, independent of the converter's scaled norm.
    gamma=np.sqrt(1+np.sum(u*u,axis=1))
    distance=LD(phase_space['parameter_plane_z_m'])-LD(phase_space['output_plane_z_m'])
    xy=records['plane_position_m'].astype(LD)+distance*u[:,:2]/u[:,2,None]
    time=records['arrival_time_offset_s'].astype(LD)+distance*gamma/(C*u[:,2])
    variables=np.column_stack([xy,u[:,:2],gamma])
    return variables,time,gamma


def analytic_phase_space(config):
    p=config['phase_space'];sx=LD(p['sigma_x_m']);sy=LD(p['sigma_y_m'])
    ex=LD(p['emittance_n_x_m']);ey=LD(p['emittance_n_y_m'])
    ax=LD(p['alpha_x']);ay=LD(p['alpha_y']);sg=LD(p['gamma_sigma'])
    profile=config['profile'];lo=LD(profile['window_start_s'])/LD(profile['full_duration_s']);hi=LD(profile['window_stop_s'])/LD(profile['full_duration_s'])
    norm=(3*hi*hi-2*hi**3)-(3*lo*lo-2*lo**3)
    def moment(power):return 6*((hi**(power+2)-lo**(power+2))/(power+2)-(hi**(power+3)-lo**(power+3))/(power+3))/norm
    rt=moment(1);vt=moment(2)-rt*rt
    mean=np.array([0,0,0,0,p['gamma_mean'],rt],dtype=LD)
    cov=np.diag(np.array([sx*sx,sy*sy,(ex/sx)**2*(1+ax*ax),
        (ey/sy)**2*(1+ay*ay),sg*sg,vt],dtype=LD))
    cov[0,2]=cov[2,0]=-ex*ax;cov[1,3]=cov[3,1]=-ey*ay
    fourth=moment(4)-4*rt*moment(3)+6*rt*rt*moment(2)-3*rt**4
    return mean,cov,fourth


def profile_bin_probabilities(config,edges):
    """Integrate the quadratic current independently, without inverse CDF."""
    T=LD(config['profile']['full_duration_s'])
    x=np.asarray(edges,dtype=LD)/T
    primitive=3*x*x-2*x*x*x
    probs=np.diff(primitive)
    return probs/sum(probs)


def distribution_report(config,n,sums,products,observed_bins,expected_bins,edges):
    target,targetcov,time_fourth=analytic_phase_space(config)
    centered_mean=sums/n;mean=target+centered_mean;cov=products/n-np.outer(centered_mean,centered_mean)
    scale=np.sqrt(np.diag(targetcov));mean_z=np.zeros(6,dtype=LD)
    active=scale>0;mean_z[active]=(mean[active]-target[active])/(scale[active]/np.sqrt(LD(n)))
    stderr=np.sqrt((targetcov*targetcov+np.outer(np.diag(targetcov),np.diag(targetcov)))/LD(n))
    stderr[5,5]=np.sqrt(max(LD(0),time_fourth-targetcov[5,5]**2)/LD(n))
    covariance_z=np.divide(cov-targetcov,stderr,out=np.zeros_like(cov),where=stderr>0)
    fixed=scale==0
    tolerance=256*np.finfo(float).eps*np.maximum(1,abs(target))
    need(np.all(abs(mean[fixed]-target[fixed])<=tolerance[fixed]),'A deterministic phase-space variable changed')
    for index in np.flatnonzero(fixed):
        bound=tolerance[index]*(scale+tolerance)
        need(np.all(abs(cov[index])<=bound),'A deterministic variable acquired a covariance')
    probs=profile_bin_probabilities(config,edges)
    expect=LD(n)*probs
    count_variance=LD(n)*probs*(1-probs)
    count_z=np.divide(observed_bins-expect,np.sqrt(count_variance),
        out=np.zeros_like(expect),where=count_variance>0)
    pearson=float(np.sum((observed_bins-expect)**2/expect))
    result=dict(parameter_plane_z_m=config['phase_space']['parameter_plane_z_m'],
        variable_order=['x_m','y_m','ux','uy','gamma','parameter_time_over_full_duration'],sample_mean=np.asarray(mean,float).tolist(),
        theoretical_mean=np.asarray(target,float).tolist(),sample_covariance=np.asarray(cov,float).tolist(),
        theoretical_covariance=np.asarray(targetcov,float).tolist(),
        mean_standard_errors=np.asarray(mean_z,float).tolist(),
        maximum_covariance_standard_errors=float(np.max(abs(covariance_z))),
        maximum_bin_count_standard_errors=float(np.max(abs(count_z))),
        count_pearson=pearson,count_pearson_degrees_of_freedom=len(probs)-1,
        generation_bin_counts_recovered_exactly=bool(np.array_equal(observed_bins,expected_bins)),
        scope='Full6D Gaussian transverse/energy plus analytic parabolic time moments and cross-covariances; fixed-N multinomial counts. Not complete physical acceptance')
    need(np.max(abs(mean_z))<9,'Parameter-plane mean differs by more than nine standard errors')
    need(np.max(abs(covariance_z))<10,'Parameter-plane covariance differs by more than ten standard errors')
    need(np.max(abs(count_z))<10,'A generation-bin count differs by more than ten standard errors')
    need(np.array_equal(observed_bins,expected_bins),'Recovered generation-bin occupancy differs')
    return result


def verify_master(master,chunk_records=65536,check_distribution=True):
    before=identity(master);n=0;sums=np.zeros(6,dtype=LD);products=np.zeros((6,6),dtype=LD)
    max_gamma=0.;minimum_uz=float('inf');tmin=float('inf');tmax=-float('inf')
    duplicate_adjacent=0;previous=None;observed_bins=None;payload=hashlib.sha256()
    with h5py.File(master,'r') as f:
        data,config=master_gate(f);N=len(data);plane=float(f['particles'].attrs['input_plane_z_m'])
        edges=f['generation/bin_edges_time_s'][:].astype(LD) if 'generation/bin_edges_time_s' in f else None
        expected_bins=f['generation/bin_counts'][:] if 'generation/bin_counts' in f else None
        shift=analytic_phase_space(config)[0] if config is not None else None
        if check_distribution:
            need(config is not None and edges is not None and expected_bins is not None,
                 'Distribution validation needs saved configuration and generation bins')
            observed_bins=np.zeros(len(edges)-1,dtype=np.int64)
        for start in range(0,N,chunk_records):
            a=data[start:min(N,start+chunk_records)]
            payload.update(a.tobytes(order='C'))
            need(np.array_equal(a['source_id'],np.arange(start+1,start+len(a)+1,dtype=np.uint64)),
                 'Master source IDs must be contiguous one-based record IDs')
            need(np.all(a['macro_weight']==1),'Master contains a non-unit macro weight')
            for key in a.dtype.names:
                if a.dtype[key].base.kind=='f':need(np.isfinite(a[key]).all(),'Nonfinite master field: '+key)
            need(np.all(a['proper_velocity'][:,2]>0),'Nonforward master particle')
            minimum_uz=min(minimum_uz,float(a['proper_velocity'][:,2].min()))
            gamma=np.sqrt(1+np.sum(a['proper_velocity'].astype(LD)**2,axis=1));max_gamma=max(max_gamma,float(gamma.max()))
            tmin=min(tmin,float(a['arrival_time_offset_s'].min()));tmax=max(tmax,float(a['arrival_time_offset_s'].max()))
            xy=a['plane_position_m']
            duplicate_adjacent+=int(np.sum(np.all(xy[1:]==xy[:-1],axis=1)))
            if previous is not None:duplicate_adjacent+=int(np.array_equal(previous,xy[0]))
            previous=xy[-1].copy()
            if config is not None:
                values,tau,_=recover_parameter_plane(a,config['phase_space'])
                values=np.column_stack([values,tau/LD(config['profile']['full_duration_s'])])
                values-=shift
                sums+=np.sum(values,axis=0,dtype=LD);products+=values.T@values
                if check_distribution:
                    lo,hi=LD(config['profile']['window_start_s']),LD(config['profile']['window_stop_s'])
                    tolerance=LD('2e-25')
                    need(np.all(tau>=lo-tolerance) and np.all(tau<=hi+tolerance),
                         'Recovered time lies outside configured current support')
                    bins=np.searchsorted(edges,tau,side='right')-1
                    bins=np.where((bins==len(edges)-1)&(abs(tau-edges[-1])<=tolerance),len(edges)-2,bins)
                    need(np.all((bins>=0)&(bins<len(observed_bins))),'Event outside generation bins')
                    observed_bins+=np.bincount(bins.astype(np.int64),minlength=len(observed_bins))
            n+=len(a)
        need(n==N,'Master scan truncated')
        claimed=f['generation'].attrs.get('record_payload_sha256','')
        if isinstance(claimed,bytes):claimed=claimed.decode()
        if claimed:need(payload.hexdigest()==str(claimed),'Master record-payload checksum differs')
        result=dict(passed=True,records=N,physical_electrons=N,macro_weights='All exactly one',
            actual_charge_C=float(LD(N)*QE),input_plane_z_m=plane,minimum_uz=minimum_uz,
            maximum_reconstructed_gamma=max_gamma,arrival_time_range_s=[tmin,tmax],
            adjacent_exact_xy_duplicates=duplicate_adjacent,
            duplicate_scope='Adjacent source records only; full pair-correlation tests are a separate small-sample check',
            master_identity=before,record_payload_sha256=payload.hexdigest())
        if check_distribution:result['distribution']=distribution_report(config,N,sums,products,observed_bins,expected_bins,edges)
    need(before==identity(master),'Closed master changed during verification')
    return result


def reference_rows(data,ids,chunk_records):
    """Bounded-memory gather from a contiguous one-based master ID contract."""
    result=np.empty(len(ids),dtype=DTYPE)
    block=(ids-1)//chunk_records
    for key in np.unique(block):
        mask=block==key;start=int(key)*chunk_records
        values=data[start:min(start+chunk_records,len(data))]
        result[mask]=values[(ids[mask]-1-start).astype(np.int64)]
    return result


def verify_genesis(master,genesis,chunk_records=65536):
    source_before=identity(master);before=identity(genesis)
    max_error=dict(x_m=0.,y_m=0.,ux=0.,uy=0.,gamma=0.,time_s=0.,phase_rad=0.,slice_current_A=0.)
    count=0;nonempty=0;empty=0
    with h5py.File(master,'r') as m,h5py.File(genesis,'r') as f:
        data,_=master_gate(m);N=len(data)
        need('conversion/complete' in f and scalar(f['conversion/complete'])==1,'Genesis conversion is incomplete')
        need(int(scalar(f['one4one']))==1 and int(scalar(f['beamletsize']))==1,'Genesis is not a unit-electron source')
        wavelength=LD(scalar(f['slicelength']));spacing=LD(scalar(f['slicespacing']));origin=LD(scalar(f['refposition']))
        need(wavelength>0 and spacing==wavelength,'First converter contract requires sample=1')
        total_slices=int(scalar(f['slicecount']));anchor=LD(f.attrs['time_anchor_s'])
        need(int(f.attrs['electron_count'])==N,'Genesis electron count differs')
        claimed=f.attrs.get('source_master_sha256','');claimed=claimed.decode() if isinstance(claimed,bytes) else str(claimed)
        need(claimed==sha(master),'Genesis master SHA256 differs')
        if 'conversion/slice_counts' in f:declared=f['conversion/slice_counts'][:]
        else:declared=None
        need(total_slices>0,'No Genesis slices')
        id_seen=np.zeros(N,dtype=np.bool_) # One byte per electron; explicit audit memory.
        for j in range(total_slices):
            key=f'slice{j+1:06d}';need(key in f,'Missing Genesis slice '+key);group=f[key]
            need(all(k in group for k in ['x','y','px','py','gamma','theta','source_id','current']),
                 'Native slice lacks coordinate or audit ID')
            n=len(group['x']);need(all(len(group[k])==n for k in ['y','px','py','gamma','theta','source_id']),
                'Native slice arrays differ in length')
            if declared is not None:need(int(declared[j])==n,'Declared slice population differs')
            if n:nonempty+=1
            else:empty+=1
            target_current=LD(n)*QE*C/spacing
            current=LD(scalar(group['current']));current_error=float(abs(current-target_current))
            max_error['slice_current_A']=max(max_error['slice_current_A'],current_error)
            need(current_error<float(abs(target_current)*LD('3e-15')+LD('1e-22')),'Slice current is not actual unit-electron count')
            previous=0
            for start in range(0,n,chunk_records):
                stop=min(n,start+chunk_records);ids=group['source_id'][start:stop].astype(np.uint64)
                need(np.all((ids>=1)&(ids<=N)),'Unknown Genesis source ID')
                need(np.all(ids[1:]>ids[:-1]) and (not len(ids) or int(ids[0])>previous),'Slice IDs not stable increasing order')
                if len(ids):previous=int(ids[-1])
                need(not np.any(id_seen[ids-1]),'Duplicate Genesis source ID');id_seen[ids-1]=True
                a=reference_rows(data,ids,chunk_records);u=a['proper_velocity'].astype(LD)
                gamma=np.sqrt(1+np.sum(u*u,axis=1));theta=group['theta'][start:stop].astype(LD)
                values={k:group[k][start:stop] for k in ['x','y','px','py','gamma']}
                need(all(np.isfinite(v).all() for v in values.values()) and np.isfinite(theta).all(),'Nonfinite native coordinate')
                need(np.all((theta>=0)&(theta<2*np.arccos(LD(-1))+LD('1e-14'))),'Slice phase outside sample=1 interval')
                s=origin+LD(j)*spacing+theta*wavelength/(2*np.arccos(LD(-1)))
                time=anchor-s/C
                delta=time-a['arrival_time_offset_s'].astype(LD)
                errors=dict(x_m=np.max(abs(values['x'].astype(LD)-a['plane_position_m'][:,0])),
                    y_m=np.max(abs(values['y'].astype(LD)-a['plane_position_m'][:,1])),
                    ux=np.max(abs(values['px'].astype(LD)-u[:,0])),uy=np.max(abs(values['py'].astype(LD)-u[:,1])),
                    gamma=np.max(abs(values['gamma'].astype(LD)-gamma)),time_s=np.max(abs(delta)),
                    phase_rad=np.max(abs(delta*2*np.arccos(LD(-1))*C/wavelength)))
                for name,error in errors.items():max_error[name]=max(max_error[name],float(error))
                need(errors['x_m']==0 and errors['y_m']==0 and errors['ux']==0 and errors['uy']==0,
                     'Genesis changed transverse master coordinates')
                need(errors['gamma']<LD('4e-13')*max(LD(1),np.max(gamma)),'Native gamma differs from stored master mass shell')
                need(errors['phase_rad']<LD('3e-9'),'Native conversion changed optical arrival phase')
                count+=len(ids)
        need(count==N and np.all(id_seen),'Genesis conversion lost master particles')
        result=dict(passed=True,records=count,slices=total_slices,nonempty_slices=nonempty,empty_slices=empty,
            exact_unit_electrons=True,maximum_errors=max_error,source_identity=source_before,genesis_identity=before,
            audit_bitmap_bytes=N,scope='Deterministic initial conversion only; actual Genesis import/dump tested separately')
    need(before==identity(genesis) and source_before==identity(master),'Closed input changed during conversion verification')
    return result


def verify(master,genesis=None,chunk_records=65536,check_distribution=True):
    need(chunk_records>0,'Chunk size must be positive')
    result=dict(passed=True,master=verify_master(master,chunk_records,check_distribution))
    if genesis is not None:result['genesis']=verify_genesis(master,genesis,chunk_records)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('master',type=Path)
    parser.add_argument('--genesis',type=Path);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--chunk-records',type=int,default=65536);parser.add_argument('--structural-only',action='store_true')
    args=parser.parse_args();need(not args.output.exists(),'Refusing to overwrite verification output')
    try:result=verify(args.master,args.genesis,args.chunk_records,not args.structural_only)
    except Exception as error:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(dict(passed=False,error=str(error)),indent=2)+'\n');raise
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=True,records=result['master']['records'],output=str(args.output)),indent=2))

if __name__=='__main__':main()
