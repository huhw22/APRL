"""Prebunched, finite-pulse FEL test and imported-particle GENESIS references.

All APRL methods share exactly the same crossing events. GENESIS receives
the same quiet transverse distribution and optical phases, with an independent
magnet-only entrance-fringe transport to its period-averaged start plane.
"""
import argparse
import copy
import json
import math
import shutil
from pathlib import Path
import h5py
import numpy as np
import yaml
from scipy.constants import c,e,m_e
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import j0,j1
from prepare import save,sha,require

G=929.1305492858543
GAMMA=1174.171072271572
DZ=1.2543262415359034e-6
DT=DZ/c
PERIOD=.0234
K=1.092714657005
ENTRANCE=.05
NCYCLE=128
NPERIOD=100
CURRENT=1000.

def entrance_reference(x,y):
    """SI Lorentz equations parametrized by z; magnetic gamma is exact."""
    ku=2*np.pi/PERIOD;B0=K*m_e*c*ku/e
    extent=np.sqrt(-2*np.log(1e-9))/ku
    n=len(x)
    initial=np.array([x,y,np.zeros(n),np.zeros(n),np.zeros(n)]).ravel()
    def rhs(z,state):
        xx,yy,ux,uy,tau=state.reshape(5,n)
        uz=np.sqrt(GAMMA*GAMMA-1-ux*ux-uy*uy)
        distance=ENTRANCE-z
        if distance>=extent:
            a=b=0.
        else:
            q=distance/extent
            taper=1-10*q**3+15*q**4-6*q**5
            derivative=(-30*q*q+60*q**3-30*q**4)/extent
            gaussian=np.exp(-.5*(ku*distance)**2)
            a=gaussian*(derivative-ku*ku*distance*taper)/ku
            b=gaussian*taper
        By=B0*np.cosh(ku*yy)*a
        Bz=B0*np.sinh(ku*yy)*b
        f=e/(m_e*c)
        return np.array([ux/uz,uy/uz,f*(By-uy/uz*Bz),f*ux/uz*Bz,
                         (1+ux*ux+uy*uy)/(c*uz*(GAMMA+uz))]).ravel()
    sol=solve_ivp(rhs,(0,ENTRANCE),initial,method='DOP853',rtol=2e-11,
                  atol=1e-19,max_step=PERIOD/128)
    require(sol.success,'Independent entrance orbit failed')
    xx,yy,ux,uy,tau=sol.y[:,-1].reshape(5,n)
    # Remove known quiver momentum at sin-field central entrance.
    canonical_x=ux+K*np.cosh(2*np.pi*yy/PERIOD)
    return xx,yy,canonical_x,uy,tau

def imported_beam(path,sample,current,wavelength,phase,transport):
    xx,yy,px,py,tau=transport
    ntrans=len(xx);nphase=len(phase)
    grid=np.meshgrid(np.arange(nphase),np.arange(ntrans),indexing='ij')
    ip,it=map(np.ravel,grid)
    omega=2*np.pi*c/wavelength
    # A common carrier phase is arbitrary because the reference starts at E=0.
    theta=-phase[ip]-omega*(tau[it]-np.mean(tau))
    values=dict(x=xx[it],y=yy[it],px=px[it],py=py[it],theta=theta,
                gamma=np.full(len(it),GAMMA))
    pre=128//sample;active=NCYCLE//sample;post=128//sample
    with h5py.File(path,'x') as f:
        for key,value,dtype in [('refposition',-128*wavelength,'f8'),
           ('slicelength',wavelength,'f8'),('slicespacing',sample*wavelength,'f8'),
           ('slicecount',pre+active+post,'i4'),('beamletsize',16,'i4'),('one4one',0,'i4')]:
            f.create_dataset(key,data=np.array([value],dtype=dtype))
        for index in range(pre+active+post):
            g=f.create_group(f'slice{index+1:06d}')
            g.create_dataset('current',data=[current if pre<=index<pre+active else 0.])
            for key,value in values.items():g.create_dataset(key,data=value)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('workspace',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--reuse-references',type=Path)
    a=ap.parse_args();w=a.workspace.resolve();root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    scripts=root/'scripts';shutil.copytree(Path(__file__).parent,scripts,ignore=shutil.ignore_patterns('__pycache__'))
    binary=w/'build/longitudinal-main-flow-r1-20261103/aprl'
    require(sha(binary)=='c693d067153934ea5b683305bff5d5be27284cd065396aa0cd11181c392720d0','Unexpected APRL binary')
    wavelength=PERIOD*(1+K*K/2)/(2*GAMMA*GAMMA)
    sigma=20e-6
    rr=sigma*np.sqrt(-2*np.log(1-(np.arange(8)+.5)/8))
    rr*=sigma/np.sqrt(np.mean(rr*rr)/2)
    angle=2*np.pi*np.arange(8)[None,:]/8+np.arange(8)[:,None]*np.pi/8
    xx=(rr[:,None]*np.cos(angle)).ravel();yy=(rr[:,None]*np.sin(angle)).ravel()
    bparam=brentq(lambda q:j1(2*q)-.03,.02,.05,xtol=1e-15)
    phase=(np.arange(16)+.5)*2*np.pi/16
    phase=phase-2*bparam*np.sin(phase)
    cc,pp,tt=np.meshgrid(np.arange(NCYCLE),np.arange(16),np.arange(64),indexing='ij')
    cc,pp,tt=map(np.ravel,(cc,pp,tt))
    arrival=(cc+phase[pp]/(2*np.pi))*wavelength/c
    electrons=CURRENT*NCYCLE*wavelength/(e*c)
    dtype=[('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),
           ('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')]
    records=np.zeros(len(cc),dtype=dtype)
    records['plane_position_m']=np.c_[xx[tt],yy[tt]]
    records['arrival_time_offset_s']=arrival
    records['proper_velocity'][:,2]=np.sqrt((GAMMA-1)*(GAMMA+1))
    records['source_id']=np.arange(1,len(records)+1)
    records['macro_weight']=1.
    inputfile=root/'particles.h5'
    if a.reuse_references:
        previous=a.reuse_references.resolve()
        require((previous/'EXIT_CODE').exists(),'Previous allocation is not closed')
        with h5py.File(previous/'particles.h5','r') as f:
            require(np.array_equal(f['particles/records'][:],records),'Physical input changed')
        inputfile.hardlink_to(previous/'particles.h5')
    else:
        with h5py.File(inputfile,'x') as f:
            g=f.create_group('particles');g.attrs['format_version']=4
            g.attrs['scope']='Cold quiet finite-pulse coherent FEL test; weighted continuum, not SASE shot-noise validation'
            g.create_dataset('records',data=records)
    gb=np.sqrt((G-1)*(G+1));bb=gb/G;beta=np.sqrt(1-1/GAMMA**2)
    u=np.longdouble(records['proper_velocity'][0,2])
    plus=(np.sqrt(1+u*u)+u)/(np.longdouble(G)+np.longdouble(gb))
    uz=(plus-1/plus)/2;held=float(uz/np.sqrt(1+uz*uz))
    mean=np.mean(arrival.astype(np.longdouble))
    center=float(max(np.longdouble(1e-7),np.longdouble(beta)*c*(np.max(arrival)-mean)+256*np.finfo(float).eps))
    ref=float(np.longdouble(center)/(beta*c)+mean)
    t0=float(np.min(arrival.astype(np.longdouble)-ref-2*DZ/(G*c*(bb+held))))
    source=-np.longdouble(G)*c*(bb+held)*(arrival.astype(np.longdouble)-ref-t0)
    translation=float(-.5*(source.min()+source.max()))
    last=ENTRANCE+NPERIOD*PERIOD+.3
    pulse=NCYCLE*wavelength/c
    # First lab samples are early by (1-beta_b) z/c. Include the full finite
    # bunch, undulator slippage, free drift lag, fringe lag and a 2 fs margin.
    span=pulse+NPERIOD*wavelength/c+(1-bb)*last/c+.35/(2*GAMMA**2*c)+2e-15
    steps=math.ceil(span*G/DT);duration=(steps+.5)*DT
    active=math.ceil(last/(G*DZ))+64+128+2
    need=max((1+held)*(steps+2)+66-translation/DZ,translation/DZ-64+active+2)
    full=math.ceil(2*need/256)*256
    origin=-.5*full*DZ;offset=math.floor((translation-origin)/DZ)-64
    require(offset-(1+held)*(steps+2)-9>0 and offset+active+9<full,'Reference does not cover all field queries')
    initial_reference=G*G*(center+bb*c*t0);stop=initial_reference+gb*DZ*steps
    base=yaml.safe_load((w/'results/source_coupling_full_comparison_20261101/free_drift_0p5m/run.yaml').read_text())
    planes=[('p'+str(n),ENTRANCE+n*PERIOD) for n in [20,40,60,80,100]]+[('far',last)]
    cases=[]
    for name,coupling,shape,feedback in [
      ('boris','boris','linear',True),('energy','conservative-source','linear',True),
      ('momentum','conservative-source','quintic-compensated',True),
      ('prescribed','boris','linear',False)]:
        d=root/name;d.mkdir()
        cfg=copy.deepcopy(base)
        cfg['output']['manifest']=str(d/'run-manifest.yaml')
        cfg['mesh'].update(cells=[80,80,full],center=[0.,0.,0.],duration=duration,
                          particle_coupling=coupling,longitudinal_shape=shape,
                          self_field_feedback=feedback,particle_steps_per_undulator_period=80)
        cfg['mesh']['moving_window'].update(cells=active,reference_directory=str(d/'reference'))
        cfg['beam']['reference'].update(input_plane_z=0.,initial_center_z=center,activation_plane_z=0.)
        cfg['beam']['input'].update(file=str(inputfile),electrons=electrons)
        cfg['sources']['magnetic_elements']=[{'type':'planar-undulator','entrance_z':ENTRANCE,
          'period':PERIOD,'periods':NPERIOD,'strength_parameter':K,
          'polarization_angle_rad':np.pi/2,'gaussian_fringe':True,'fringe_relative_cutoff':1e-9}]
        cfg['detectors'].update(directory=str(d/'detectors'),
          field_planes=[dict(name=n,z=z,rhythm=DT/G,buffer_samples=8,compression=0,
                             particle_background={'enabled':False}) for n,z in planes],
          particle_planes=[dict(name=n if n=='entry' else n+'-particles',z=z,
                               buffer_records=8192,compression=0) for n,z in [('entry',0.),*planes]],
          energy_planes=[])
        cfg['radiation_resolution']['maximum_photon_energy_eV']=100.
        cfg['runtime']['resource_monitor']['progress_interval_steps']=50
        cfg['energy_ledger'].update(enabled=True,directory=str(d/'ledger'),
                                   sample_interval_steps=10,discrete_operator_audit=True)
        cfg['trajectory']['enabled']=False
        cfg['stop']={'mode':'reference-center-z','z':stop}
        detector_names=[p['name'] for key in ['field_planes','particle_planes','energy_planes']
                        for p in cfg['detectors'][key]]
        require(len(detector_names)==len(set(detector_names)),'Detector names collide')
        # yaml safe_dump rejects NumPy scalar parameters.
        from prepare import native
        (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
        limits={'boris':1800,'energy':3600,'momentum':4800,'prescribed':1800}
        cases.append(dict(name=name,directory=str(d),timeout_seconds=limits[name],feedback=feedback))
    transport=entrance_reference(xx,yy)
    np.savez(root/'entrance_reference.npz',x=transport[0],y=transport[1],px=transport[2],
             py=transport[3],tau_s=transport[4],initial_x=xx,initial_y=yy)
    genesis=Path('/public/software/Genesis-1.3-Version4-4.6.4/build/genesis4')
    refs=[]
    variants=[('strong',2,241,4,CURRENT,True),('weak',2,241,4,1.,True),
              ('time-refined',1,241,4,CURRENT,True),
              ('space-refined',2,481,4,CURRENT,True),
              ('step-refined',2,241,8,CURRENT,True),
              ('without-longitudinal-SC',2,241,4,CURRENT,False)]
    for name,sample,ngrid,nstep,current,sc in variants:
        if a.reuse_references:
            previous_plan=json.loads((previous/'campaign.json').read_text())
            q=next(q for q in previous_plan['genesis_references'] if q['name']==name)
            status=json.loads((Path(q['directory'])/'genesis.log.launch.json').read_text())
            require(status.get('exit_code')==0,'Incomplete reference process')
            with h5py.File(Path(q['directory'])/'reference.out.h5','r') as f:
                require(abs(f['Lattice/z'][-1]+f['Lattice/dz'][-1]-(NPERIOD*PERIOD+.3))<1e-9,
                        'Incomplete reference lattice')
            require(sha(Path(q['directory'])/'reference.in')==
                    json.loads((previous/'FROZEN.json').read_text())[str(Path(q['directory'])/'reference.in')],
                    'Reference input changed')
            refs.append(dict(q,reuse_closed=True))
            continue
        d=root/'genesis'/name;d.mkdir(parents=True)
        imported_beam(d/'beam.par.h5',sample,current,wavelength,phase,transport)
        fielddump=20*nstep if name in ['strong','weak','space-refined'] else 0
        deck=f"""&setup
rootname=reference
lattice=reference.lat
beamline=FEL
lambda0={wavelength:.17g}
gamma0={GAMMA:.17g}
delz={PERIOD/nstep:.17g}
shotnoise=false
npart=1024
nbins=16
beam_global_stat=true
field_global_stat=true
&end
&importbeam
file=beam.par.h5
time=true
&end
&field
power=0
dgrid=.0006
ngrid={ngrid}
waist_size=.00004
&end
"""
        if sc:deck+='&efield\nlongrange=true\nrmax=2\nnz=1\nnphi=0\nngrid=100\n&end\n'
        deck+=f'&track\noutput_step={nstep}\nfield_dump_step={fielddump}\n&end\n&write\nfield=far\n&end\n'
        (d/'reference.in').write_text(deck)
        (d/'reference.lat').write_text(
          f'UND: UNDULATOR={{ lambdau={PERIOD},nwig={NPERIOD},aw={K/math.sqrt(2):.17g},helical=false }};\n'
          'D: DRIFT={l=.015};\nFEL: LINE={UND,20*D};\n')
        refs.append(dict(name=name,directory=str(d),ranks=16,current_A=current,
                         sample=sample,ngrid=ngrid,substeps_per_period=nstep,timeout_seconds=600))
    ku=2*np.pi/PERIOD;zeta=K*K/(4+2*K*K);JJ=j0(zeta)-j1(zeta)
    IA=4*np.pi*8.8541878128e-12*m_e*c**3/e
    rho=(CURRENT/IA*(K*JJ)**2/(16*GAMMA**3*sigma*sigma*ku*ku))**(1/3)
    require(abs(np.mean(np.exp(1j*phase))-.03)<1e-12,'Incorrect input bunching')
    require(abs(np.mean(xx*xx)/sigma**2-1)<1e-13 and abs(np.mean(yy*yy)/sigma**2-1)<1e-13,'Incorrect transverse RMS')
    plan=dict(workspace=str(w),binary=str(binary),binary_sha256=sha(binary),
      binding_probe=str(w/'build/cluster-tools/binding_probe'),cases=cases,
      genesis_binary=str(genesis),genesis_sha256=sha(genesis),genesis_references=refs,
      parameters=dict(gamma=GAMMA,boost_gamma=G,undulator_K=K,period_m=PERIOD,
        periods=NPERIOD,entrance_z_m=ENTRANCE,resonant_wavelength_m=wavelength,
        peak_current_A=CURRENT,represented_electrons=electrons,macroparticles=len(records),
        sigma_x_m=sigma,sigma_y_m=sigma,input_bunching=.03,pulse_cycles=NCYCLE,
        rho_1D=rho,gain_length_1D_m=PERIOD/(4*np.pi*np.sqrt(3)*rho)),
      geometry=dict(nx=80,ny=80,reference_nz=full,active_nz=active,field_steps=steps,
        expected_record_s=steps*DT/G,reference_time_offset_s=ref,
        frame_time_origin_s=t0,stop_reference_z_m=stop,field_planes=planes),
      primary_radiation=dict(aperture='same PIC pixels with radius 90-130 um',band_eV_halfwidth=3,
        common_retarded_time_window=True,windows=['rectangular','energy-normalized Hann'],
        raw_total_fields=True,near_field_subtraction=False),
      limitations=['Coherent continuum and seeded feedback; not physical SASE shot noise',
        'GENESIS is period averaged; entrance-fringe radiation and transverse space charge are not included in the direct reference',
        'Finite imported quiet transverse distribution, not an exact continuum Gaussian',
        'Use measured reference convergence and model limits; do not claim acceptance from inter-method agreement alone'])
    save(root/'campaign.json',plan)
    if a.reuse_references:
        save(root/'REPAIR.json',dict(previous_campaign=str(previous),previous_job=19281,
          reason='Unique names required across detector namespaces; full GENESIS lattice end is checked independently of diagnostic output stride',
          previous_PIC_started=False,physical_input_sha256=sha(inputfile),
          kernels_changed=False,reference_reruns=0))
    files=[root/'campaign.json',binary,genesis,inputfile,root/'entrance_reference.npz',*scripts.glob('*.py')]
    files.extend(Path(q['directory'])/'run.yaml' for q in cases)
    files.extend(Path(q['directory'])/n for q in refs for n in ['reference.in','reference.lat','beam.par.h5'])
    if a.reuse_references:
        files.extend([previous/'GENESIS_STATUS.json',previous/'EXIT_CODE'])
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})
    print(root)

if __name__=='__main__':main()
