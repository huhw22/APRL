"""Build a bounded Gaussian finite-dipole test; no kernel/default edits."""
import copy,hashlib,json,math,sys
from pathlib import Path
import h5py,numpy as np,yaml
from scipy.special import roots_hermitenorm
from theory import constants,c,e
ROOT=Path(sys.argv[1]).resolve();WS=ROOT.parents[2]
G=929.1305492858543;DZ=1.2543262415359034e-6;DT=DZ/c
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,message):
    if not ok:raise RuntimeError(message)
def native(v):
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:native(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [native(x) for x in v]
    return v
def save(p,v):Path(p).write_text(json.dumps(native(v),indent=2,allow_nan=False)+'\n')
def main():
    binary=WS/'build/longitudinal-main-flow-r1-20261103/aprl'
    assert sha(binary)=='c693d067153934ea5b683305bff5d5be27284cd065396aa0cd11181c392720d0'
    cfg=yaml.safe_load((WS/'results/source_coupling_full_comparison_20261101/free_drift_0p5m/run.yaml').read_text())
    p=constants(length=.0024,screen=.3024);sx=sy=0.;sz=0.;N=1.
    nodes=[];weights=[]
    for order in [1,1,1]:
        x,w=roots_hermitenorm(order);nodes.append(x);weights.append(w/np.sqrt(2*np.pi))
    xx,yy,zz=np.meshgrid(*nodes,indexing='ij');wx,wy,wz=np.meshgrid(*weights,indexing='ij')
    w=(wx*wy*wz).ravel();w/=w.sum();count=len(w)
    dtype=[('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')]
    a=np.zeros(count,dtype=dtype);a['plane_position_m']=np.c_[xx.ravel()*sx,yy.ravel()*sy]
    a['arrival_time_offset_s']=zz.ravel()*sz/(p['beta']*c);a['proper_velocity'][:,2]=np.sqrt((p['gamma']-1)*(p['gamma']+1))
    a['source_id']=np.arange(1,count+1);a['macro_weight']=w*count
    ROOT.mkdir(parents=True,exist_ok=False)
    directory=ROOT/'template';directory.mkdir(exist_ok=False)
    inputfile=directory/'particles.h5'
    with h5py.File(inputfile,'x') as f:
        g=f.create_group('particles');g.attrs['format_version']=4
        g.attrs['scope']='Product Gauss-Hermite weighted Gaussian continuum; cold parallel input, fixed lab crossing events'
        g.create_dataset('records',data=a)
    gb=np.sqrt((G-1)*(G+1));bb=gb/G
    u=np.longdouble(a['proper_velocity'][0,2]);plus=(np.sqrt(1+u*u)+u)/(np.longdouble(G)+np.longdouble(gb));uz=(plus-1/plus)/2
    held=float(uz/np.sqrt(1+uz*uz));beta=np.longdouble(p['beta']);arrival=a['arrival_time_offset_s'].astype(np.longdouble)
    mean=np.sum(w.astype(np.longdouble)*arrival)
    center=float(max(np.longdouble(1e-7),beta*c*(np.max(arrival)-mean)+256*np.finfo(float).eps))
    ref=float(np.longdouble(center)/(beta*c)+mean)
    t0=float(np.min(arrival-ref-2*DZ/(G*c*(bb+held))))
    source=-np.longdouble(G)*c*(bb+held)*(arrival-ref-t0);translation=float(-.5*(source.min()+source.max()))
    last=.05+p['screen'];steps=math.ceil(5e-15*G/DT);duration=(steps+.5)*DT
    active=math.ceil(last/(G*DZ))+64+128+2
    need=max((1+held)*(steps+2)+66-translation/DZ,translation/DZ-64+active+2)
    full=math.ceil(2*need/256)*256;origin=-.5*full*DZ;offset=math.floor((translation-origin)/DZ)-64
    assert offset-(1+held)*(steps+2)-1>0 and offset+active<full
    initial_reference=G*G*(center+bb*c*t0)
    stop=initial_reference+gb*DZ*steps
    cfg['output']['manifest']=str(directory/'run-manifest.yaml')
    cfg['mesh'].update(cells=[96,80,full],center=[.12e-3,0.,0.],duration=duration,particle_steps_per_undulator_period=1024)
    cfg['mesh']['moving_window'].update(cells=active,reference_directory=str(directory/'reference'))
    cfg['runtime']['resource_monitor']['progress_interval_steps']=50
    cfg['beam']['reference'].update(input_plane_z=0.,initial_center_z=center,activation_plane_z=0.)
    cfg['beam']['input'].update(file=str(inputfile),electrons=N)
    cfg['sources']['magnetic_elements']=[
      {'type':'uniform-dipole','entrance_z':.05,'length':p['length'],'field_T':p['field_T'],'polarization_angle_rad':np.pi/2},
      {'type':'planar-undulator','entrance_z':.02,'period':.009,'periods':1,'strength_parameter':0.,'gaussian_fringe':False}]
    # A zero-strength element activates EXISTING Boris subcycling, since the
    # current program derives its sampling only from undulator periods. It has
    # exactly zero field and adds no physical magnet or interaction extent.
    cfg['detectors'].update(directory=str(directory/'detectors'),field_planes=[{'name':'screen','z':last,'rhythm':DT/G,'buffer_samples':8,'compression':0,'particle_background':{'enabled':False}}],
      particle_planes=[{'name':name,'z':z,'buffer_records':8192,'compression':0} for name,z in [('entry',0.),('bend-exit',.05+p['length']),('exit',last)]],
      energy_planes=[{'name':name,'z':z,'rhythm':DT/G,'buffer_samples':64,'compression':0} for name,z in [('entry-energy',0.),('exit-energy',last)]])
    cfg['radiation_resolution']['maximum_photon_energy_eV']=200.
    cfg['energy_ledger'].update(directory=str(directory/'ledger'),sample_interval_steps=10,enabled=True,discrete_operator_audit=True)
    cfg['trajectory']['enabled']=False
    cfg['stop']={'mode':'reference-center-z','z':stop}
    (directory/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
    first=t0+bb*last/c;magtime=-ref+.05/(p['beta']*c);first_u=first-magtime-p['screen']/c
    # Product quadrature must reproduce the analytic Gaussian form factor in
    # the registered observation band, before seeing any PIC output.
    tests=[]
    for ev in [8.,10.,12.]:
        for nx in [-.0002,0.,.0002]:
            omega=ev*e/1.0545718176461565e-34;k=omega/c
            phase=omega*a['arrival_time_offset_s']-k*nx*a['plane_position_m'][:,0]
            finite=np.sum(w*np.exp(1j*phase));exact=np.exp(-.5*((omega*sz/(p['beta']*c))**2+(k*nx*sx)**2))
            tests.append(float(abs(finite-exact)))
    assert max(tests)<1e-8
    plan={'binary':str(binary),'binary_sha256':sha(binary),'binding_probe':str(WS/'build/cluster-tools/binding_probe'),
      'case_directory':str(directory),'source_core':'a157f1c',
      'macroparticles':count,'represented_electrons':N,'input_sha256':sha(inputfile),
      'parameters':p,'sigma_x_m':sx,'sigma_y_m':sy,'sigma_z_m':sz,
      'magnet_entrance_z_m':.05,'screen_lab_z_m':last,'reference_time_offset_s':ref,
      'reference_frame_time_origin_s':t0,'screen_first_relative_u_s':first_u,
      'registered_band_eV':[50.,100.],'registered_aperture_m':[-10e-6,10e-6,-10e-6,10e-6],
      'comparison':'PIC raw same-record signed field spectrum vs independent finite-arc acceleration and complete LW fields; no fit, no subtraction',
      'geometry':{'reference_nz':full,'active_nz':active,'nx':96,'ny':80,'field_steps':steps,'mesh_duration_s':duration,'stop_reference_lab_z_m':stop},
      'form_factor_quadrature_max_absolute_error':max(tests),
      'zero_field_sampling_element':'K=0; only activates existing prescribed-device Boris subcycling; no kernel change',
      'timeout_seconds':2400,'active_budget_seconds':2600,'physics_accepted':False}
    save(ROOT/'PLAN.json',plan)
    make_methods(plan,cfg,directory)
def make_methods(common,config,template):
    import shutil
    from copy import deepcopy
    scripts=ROOT/'scripts'
    shutil.copytree(Path(__file__).parent,scripts,ignore=shutil.ignore_patterns('__pycache__'))
    cases=[]
    for label,coupling,shape in [('boris','boris','linear'),('energy','conservative-source','linear'),('momentum','conservative-source','quintic-compensated')]:
        d=ROOT/label;d.mkdir()
        def rewrite(x):
            if isinstance(x,dict):return {k:rewrite(v) for k,v in x.items()}
            if isinstance(x,list):return [rewrite(v) for v in x]
            return x.replace(str(template),str(d)) if isinstance(x,str) else x
        cfg=rewrite(config);cfg['beam']=deepcopy(config['beam'])
        cfg['mesh']['particle_coupling']=coupling;cfg['mesh']['longitudinal_shape']=shape
        (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
        plan=dict(common,method=label,case_directory=str(d),requested_conservation='Energy iteration and approximate longitudinal matching; not exact lab four-momentum closure')
        save(d/'PLAN.json',plan)
        cases.append(dict(name=label,directory=str(d),timeout_seconds=2400))
    save(ROOT/'campaign.json',dict(common=common,cases=cases,binary=common['binary'],binding_probe=common['binding_probe'],scope='Identical 0.8mrad finite bend in production gamma/boost/dx/dz; continuum coherent radiation only'))
    files=[ROOT/'campaign.json',Path(common['binary']),template/'particles.h5',*scripts.glob('*.py')]
    for case in cases:files.extend([Path(case['directory'])/'run.yaml',Path(case['directory'])/'PLAN.json'])
    save(ROOT/'FROZEN.json',{str(p):sha(p) for p in files})
    print(ROOT)
if __name__=='__main__':main()
