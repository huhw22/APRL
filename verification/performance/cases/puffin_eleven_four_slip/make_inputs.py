"""Generate only selected records of a conceptual independent physical source."""
import concurrent.futures
import math
import os
from pathlib import Path
import h5py
import numpy as np
from scipy.fft import next_fast_len
from common import C,E,read,save,need,sha
from source_math import validate_config,plan_bins,sample_records,DTYPE
from source_verify import analytic_phase_space,recover_parameter_plane,distribution_report
from puffin_units import scaling


def generate(root):
    plan=read(root/'PLAN.json');cfg=validate_config(plan['source_config']);bins=plan_bins(cfg)
    n=cfg['electrons'];count=plan['sample_records']
    index=np.sort(np.random.Generator(np.random.Philox(cfg['seed']+1000)).choice(n,count,replace=False))
    a=np.empty(count,dtype=DTYPE);prefix=np.r_[0,np.cumsum(bins['counts'])]
    tasks=[]
    for j in range(len(bins['counts'])):
        for start in range(int(prefix[j]),int(prefix[j+1]),262144):
            stop=min(start+262144,int(prefix[j+1]));lo,hi=np.searchsorted(index,[start,stop])
            if hi>lo: tasks.append((j,start,stop,int(lo),int(hi)))
    def block(task):
        j,start,stop,lo,hi=task
        values,_=sample_records(cfg,start,stop-start,bins['cdf_values'][j],bins['cdf_values'][j+1])
        return lo,hi,values[index[lo:hi]-start]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for k,(lo,hi,values) in enumerate(pool.map(block,tasks)):
            a[lo:hi]=values
            if k%128==0: print('SOURCE_BLOCK',k,len(tasks),hi,count,flush=True)
    need(np.array_equal(a['source_id'],index.astype('u8')+1),'Selected physical IDs differ')
    # A unique x key is stronger than the exact (x,y) clone rejection and is
    # also used to verify actual native import after Puffin redistribution.
    need(len(np.unique(a['plane_position_m'][:,0]))==count,'Repeated transverse x / cloned state')
    selected_bins=np.searchsorted(prefix,a['source_id']-1,side='right')-1
    counts=np.bincount(selected_bins,minlength=len(bins['counts']))
    mean=analytic_phase_space(cfg)[0]
    total=np.zeros_like(mean);second=np.zeros((len(mean),len(mean)),dtype=np.longdouble)
    for lo in range(0,count,262144):
        sub=a[lo:lo+262144];v,t,_=recover_parameter_plane(sub,cfg['phase_space'])
        values=np.column_stack([v,t/np.longdouble(cfg['profile']['full_duration_s'])])-mean
        recovered=np.searchsorted(bins['edges'].astype(np.longdouble),t,side='right')-1
        need(np.array_equal(recovered,selected_bins[lo:lo+len(sub)]),'Source bin roundtrip failed')
        total+=values.sum(axis=0);second+=values.T@values
    report=distribution_report(cfg,count,total,second,counts,counts,bins['edges'])
    path=root/'inputs/independent-w64.h5'
    with h5py.File(path,'x') as f:
        g=f.create_group('particles');g.attrs.update(format_version=4,input_plane_z_m=0.,electron_count=n,
            macro_records=count,macro_weight=n/count,input_relative_weight=1.,physical_seed=cfg['seed'],
            scope='Full parabolic independent fixed-N physical source, uniform thinning, no copied transverse states')
        g.create_dataset('records',data=a,chunks=(65536,))
        g=f.create_group('generation');g['complete']=np.uint8(1);g['committed_records']=np.uint64(count)
        g['physical_bin_counts']=bins['counts'];g['selected_bin_counts']=counts;g['bin_edges_time_s']=bins['edges']
        g['config_json']=__import__('json').dumps(cfg)
    weight=n/count;phase=cfg['phase_space'];scale=scaling(phase['gamma_mean'],plan['K'],.028)
    u=a['proper_velocity'];gamma=np.sqrt(1+np.sum(u*u,axis=1));drift=.01
    arrival=a['arrival_time_offset_s']+drift*gamma/(C*u[:,2])
    xy=a['plane_position_m']+drift*u[:,:2]/u[:,2,None]
    cycle=plan['wavelength_m']/C;anchor=float(arrival.min()-128*cycle)
    native=np.column_stack((xy[:,0]/scale['transverse_scale_m'],xy[:,1]/scale['transverse_scale_m'],
        scale['code_c_m_per_s']*(arrival-anchor)/scale['cooperation_length_m'],u[:,0]/plan['K'],
        -u[:,1]/plan['K'],gamma/phase['gamma_mean'],np.full(count,weight*E/(scale['qe_C']*scale['npkbar']))))
    with h5py.File(root/'particles.h5','x') as f:
        f['electrons']=native;f['source_id']=a['source_id'];f.attrs['physical_time_origin_s']=anchor
        f.attrs['physical_entrance_m']=.01
    recovered=np.column_stack((native[:,:2]*scale['transverse_scale_m'],
        native[:,2]*scale['cooperation_length_m']/scale['code_c_m_per_s']+anchor,
        native[:,3]*plan['K'],-native[:,4]*plan['K'],native[:,5]*phase['gamma_mean']))
    expected=np.column_stack((xy,arrival,u[:,:2],gamma));error=np.max(abs(recovered-expected),axis=0)
    need(np.all(error<np.array([1e-15,1e-15,1e-24,1e-12,1e-12,1e-10])),'Puffin event mapping failed')
    charge=float(native[:,6].sum()*scale['npkbar']*scale['qe_C'])
    need(abs(charge/(n*E)-1)<2e-12,'Physical charge changed')
    span=float(np.ptp(arrival)*scale['code_c_m_per_s']/scale['lambda_reference_m'])
    # 1100 extra optical wavelengths cover nominal847 plus inter-segment
    # drift slip, with128 wavelengths on each side. Padding is not bunch length.
    needed=int(math.ceil((span+1100+256)*20+1));nz=256*next_fast_len(math.ceil(needed/256))
    info=dict(complete=True,input=str(path),input_sha256=sha(path),native_sha256=sha(root/'particles.h5'),
        physical_electrons=n,records=count,actual_weight=weight,total_charge_C=charge,scale=scale,
        arrival_origin_s=anchor,source_seed=cfg['seed'],zero_transverse_duplicates=True,
        roundtrip_max_abs_errors=error.tolist(),phase_space_statistics=report,field_nodes_z2=nz,
        field_wavelength_span=(nz-1)/20,physical_bunch_duration_s=cfg['profile']['full_duration_s'],
        independent_fixed_N_source=True,no_extra_noise=True,no_clones=True,full_parabolic_support=True,
        particle_output_has_persistent_ID=False,APRL_input_ready=True)
    save(root/'INPUT_VALIDATED.json',info)
    for directory,evolve in [(root/'import-check',False),(root/'main',True)]:
        periods=990 if evolve else 1;steps=80
        (directory/'beam.in').write_text("&NBLIST\nnbeams=1\ndtype='h5'\n/\n&BH5LIST\ndist_f='../particles.h5'\n/\n")
        (directory/'matched.latt').write_text(f"UN 'aprl-sparc11' {periods} 1.0 0.0 {steps} 0.0 1.0 0.0 0.0\n")
        (directory/'writes.in').write_text('WR 1\n'+''.join(f"WR {x['step']}\n" for x in plan['outputs']) if evolve else 'WR 1\n')
        toggle='.true.' if evolve else '.false.'
        deck=f'''&MDATA
qScaled=.true.
qOneD=.false.
qFieldEvolve={toggle}
qWriteA={toggle}
qElectronsEvolve=.true.
qElectronFieldCoupling=.true.
qDiffraction=.true.
qFocussing=.false.
qUndEnds=.true.
qFMesh_G=.false.
q_noise=.false.
qFilter=.true.
qMeasure=.false.
meshType=0
beam_file='beam.in'
lattFile='matched.latt'
wr_file='writes.in'
iNumNodesX=81
iNumNodesY=81
nodesPerLambdar=21
sFModelLengthX={.0032/scale['transverse_scale_m']:.17g}
sFModelLengthY={.0032/scale['transverse_scale_m']:.17g}
sFModelLengthZ2={(nz-1-1e-7)*scale['lambda_reference_m']/20/scale['cooperation_length_m']:.17g}
iRedNodesX=49
iRedNodesY=49
sFiltFrac=.3
sDiffFrac=.25
sBeta=1
sElectronThreshold=0
srho=.005
saw={plan['K']:.17g}
sgamma_r={phase['gamma_mean']:.17g}
lambda_w=.028
sux=0
suy=1
zundType='aprl-sparc11'
stepsPerPeriod={steps}
nPeriods={periods}
sZ0=0
iWriteNthSteps={periods*steps+1}
iWriteIntNthSteps={steps}
sRedistLen=.5
iRedistStp=20
ioutInfo=2
/
'''
        (directory/'main.in').write_text(deck)
    # The closed native check uses the same full field extent as production.
    print('INPUT_READY',count,n,weight,nz,flush=True)
    return info
