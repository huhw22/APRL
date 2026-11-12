"""Prepare a frozen Puffin source copy and the existing w64-a physical events.

No noise generation, sampling, quiet loading, or moment matching is performed.
"""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.fft import next_fast_len

W=Path('/public/home/wfollowm/APRL_Workspace')
UPSTREAM=Path('/public/home/wfollowm/Puffin')
C=299792458.0
E=1.602176634e-19


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def read(path):return json.loads(Path(path).read_text())


def sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):value.update(block)
    return value.hexdigest()


def scaling(gamma,aw,period,rho=.005):
    # The upstream constants omit _wp on literal RHS: reproduce the actual
    # default-real -> wp promotion, rather than using newer SI constants here.
    pi=float(np.float32(np.pi));qe=float(np.float32(1.60217653e-19))
    eps=float(np.float32(8.854187817e-12));me=float(np.float32(9.1093826e-31));code_c=float(np.float32(C))
    beta=np.sqrt(gamma*gamma-1-aw*aw/2)/gamma
    eta=(1-beta)/beta;lg=period/(4*pi*rho);lr=period*eta;lc=lr/(4*pi*rho)
    kappa=aw/(2*rho*gamma)
    npkbar=lg*lc**2*eps*me/qe**2*gamma**3*rho**3*(4*code_c*2*pi/period/aw)**2
    ef=gamma*me*code_c*code_c/(qe*kappa*lg)
    return dict(rho=rho,gamma_ref=gamma,aw_peak=aw,period_m=period,pi=pi,qe_C=qe,eps0=eps,me_kg=me,
        eta=eta,gain_length_m=lg,lambda_reference_m=lr,cooperation_length_m=lc,npkbar=npkbar,code_c_m_per_s=code_c,
        kappa=kappa,field_scale_V_per_m=ef,transverse_scale_m=float(np.sqrt(lg*lc)))


PATCH='''
  ! APRL-compatible external magnetic field; the Puffin solver is unchanged.
  ! This explicit type is for the existing 77-period, one-period-ramp case.
  subroutine getAPRL77(sy, sz, bxj, byj, bzj, frame)
    real(kind=wp), contiguous, intent(in) :: sy(:)
    real(kind=wp), intent(in) :: sz
    real(kind=wp), contiguous, intent(out) :: bxj(:), byj(:), bzj(:)
    type(tFELFrame), intent(in) :: frame
    real(kind=wp) :: z, length, ku, distance, envelope, derivative, ph, ft, fl
    z=sz*frame%gain_length
    length=77.0_wp*frame%lambda_w
    ku=2.0_wp*acos(-1.0_wp)/frame%lambda_w
    ft=0.0_wp
    fl=0.0_wp
    if (z >= 0.0_wp .and. z <= length) then
      envelope=1.0_wp
      derivative=0.0_wp
      distance=min(z,length-z)
      if (distance < frame%lambda_w) then
        ph=acos(-1.0_wp)*distance/frame%lambda_w
        envelope=0.5_wp*(1.0_wp-cos(ph))
        derivative=0.5_wp*acos(-1.0_wp)/frame%lambda_w*sin(ph)
        if (z >= frame%lambda_w) derivative=-derivative
      end if
      ft=envelope*sin(ku*z)-derivative/ku*cos(ku*z)
      fl=envelope*cos(ku*z)
    end if
!$OMP WORKSHARE
    bxj=0.0_wp
    byj=cosh(ku*sy*sqrt(frame%gain_length*frame%cooperation_length))*ft
    bzj=sinh(ku*sy*sqrt(frame%gain_length*frame%cooperation_length))*fl
!$OMP END WORKSHARE
  end subroutine getAPRL77

'''


def source_copy(root):
    clean=subprocess.check_output(['git','status','--porcelain'],cwd=UPSTREAM,text=True)
    if clean:raise RuntimeError('Puffin upstream checkout is not clean')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=UPSTREAM,text=True).strip()
    files=subprocess.check_output(['git','ls-files','-z'],cwd=UPSTREAM).decode().split('\0')
    target=root/'source';target.mkdir()
    hashes={}
    for name in filter(None,files):
        src=UPSTREAM/name;dst=target/name;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(src,dst);hashes[name]=sha(src)
    path=target/'puffin/lib/undulator/bfields.f90';old=path.read_text()
    marker='  call getBXfield(sx, sy, sz, bxj, und, frame)'
    addition="""  if (und%undulator_type == 'aprl-ramp77') then
    call getAPRL77(sy, sz, bxj, byj, bzj, frame)
    return
  end if

"""
    if old.count(marker)!=1:raise RuntimeError('External field insertion point changed')
    modified=old.replace(marker,addition+marker,1)
    end='end module';position=modified.lower().rfind(end)
    if position<0:raise RuntimeError('Magnetic module end missing')
    modified=modified[:position]+PATCH+modified[position:];path.write_text(modified)
    probe=Path(__file__).with_name('field_probe.f90')
    shutil.copyfile(probe,target/'field_probe.f90')
    cmake=target/'CMakeLists.txt'
    cmake.write_text(cmake.read_text()+'''\n# Isolated external-magnet verification executable.\nadd_executable(field_probe field_probe.f90)\ntarget_include_directories(field_probe PRIVATE ${CMAKE_BINARY_DIR}/puffin ${HDF5_INCLUDE_DIRS})\ntarget_link_libraries(field_probe puffin_mpi)\n''')
    save(root/'SOURCE.json',dict(upstream=str(UPSTREAM),commit=commit,original_tracked_sha256=hashes,
        changed_file='puffin/lib/undulator/bfields.f90',changed_sha256=sha(path),
        only_change='Add explicit aprl-ramp77 external magnet option; all existing paths and field/particle solvers retained.',
        build_only_changed_file='CMakeLists.txt',added_file='field_probe.f90',
        source_copy=str(target.resolve())))


def input_file(root):
    if (root/'FROZEN.json').exists():raise RuntimeError('Refuse to rewrite frozen inputs')
    parent=W/'results/sparc_random_startup_20261112_r2';plan=read(parent/'PLAN.json')
    info=read(parent/'w64-a/INPUT_VALIDATED.json');original=Path(info['input'])
    if sha(original)!=info['input_sha256']:raise RuntimeError('Physical input changed')
    with h5py.File(original,'r') as f:
        if int(f['generation/complete'][()])!=1:raise RuntimeError('Physical input incomplete')
        a=f['particles/records'][:]
    phase=plan['source_config']['phase_space'];gamma_ref=phase['gamma_mean'];aw=2.060074987322757
    scale=scaling(gamma_ref,aw,.028);lg=scale['gain_length_m'];lc=scale['cooperation_length_m']
    u=a['proper_velocity'];gamma=np.sqrt(1+np.sum(u*u,axis=1));drift=.01
    arrival=a['arrival_time_offset_s']+drift*gamma/(C*u[:,2])
    xy=a['plane_position_m'][:,:2]+drift*u[:,:2]/u[:,2,None]
    anchor=float(arrival.min()-48*plan['wavelength_m']/C)
    z2=scale['code_c_m_per_s']*(arrival-anchor)/lc;trans=scale['transverse_scale_m']
    # Puffin complex transverse momentum is px - i py, scaled by peak K.
    weight=info['actual_weight']*E/(scale['qe_C']*scale['npkbar'])
    data=np.column_stack((xy[:,0]/trans,xy[:,1]/trans,z2,u[:,0]/aw,-u[:,1]/aw,
                          gamma/gamma_ref,np.full(len(a),weight)))
    dest=root/'particles.h5'
    with h5py.File(dest,'x') as f:
        f['electrons']=data;f['source_id']=a['source_id'];f.attrs['units']='Puffin native scaled'
        f.attrs['physical_time_origin_s']=anchor;f.attrs['initial_lab_plane_m']=drift
    reconstructed=np.column_stack((data[:,:2]*trans,data[:,2]*lc/scale['code_c_m_per_s']+anchor,
                                   data[:,3]*aw,-data[:,4]*aw,data[:,5]*gamma_ref))
    expected=np.column_stack((xy,arrival,u[:,:2],gamma))
    errors=np.max(abs(reconstructed-expected),axis=0)
    if not np.all(errors<np.array([1e-15,1e-15,1e-24,1e-12,1e-12,1e-10])):
        raise RuntimeError('Native mapping roundtrip failed')
    charge=float(data[:,6].sum()*scale['npkbar']*scale['qe_C'])
    expected_charge=info['physical_electrons']*E
    if abs(charge/expected_charge-1)>2e-12:raise RuntimeError('Total charge changed')
    lam=scale['lambda_reference_m'];span=float(np.ptp(arrival)*scale['code_c_m_per_s']/lam)
    waves=int(np.ceil(span))+96
    nz_fast=int(next_fast_len(waves*20+1))
    native=dict(case='w64-a',physical_seed=2026111201,records=len(a),original_input=str(original),
        original_sha256=info['input_sha256'],native_sha256=sha(dest),scale=scale,arrival_origin_s=anchor,
        physical_entrance_m=drift,original_actual_weight=info['actual_weight'],total_charge_C=charge,
        roundtrip_max_abs_errors=errors.tolist(),roundtrip_columns=['x_m','y_m','time_s','ux','uy','gamma'],
        minimum_field_waves=waves,full_field_waves=(nz_fast-1)/20,field_nodes_z2=nz_fast,
        transverse_nodes=81,transverse_spacing_m=.00004,
        nodes_per_wavelength=21,intervals_per_wavelength=20,steps_per_period=80,
        upfront_map='Exact ballistic transport of original crossing events over common 1cm field-free drift.',
        no_extra_noise=True,no_resampling=True,no_clones=True,physical_weight_preserved=True,
        limits=['No space charge or velocity near field in Puffin.',
                'Puffin may redistribute arrays; output does not carry original source IDs. Fixed-cohort BF is not claimed.'])
    save(root/'INPUT_VALIDATED.json',native)
    for name,nperiods,steps,evolve in [('import-check',1,4,False),('short-u07',7,80,True)]:
        d=root/name;d.mkdir(exist_ok=True)
        # The no-lattice path allocates zundtype_arr without assigning the
        # requested type. A smooth custom magnet also must bypass matchIn/Out.
        (d/'matched.latt').write_text(f"UN 'aprl-ramp77' {nperiods} 1.0 0.0 {steps} 0.0 1.0 0.0 0.0\n")
        (d/'beam.in').write_text("&NBLIST\nnbeams=1\ndtype='h5'\n/\n&BH5LIST\ndist_f='../particles.h5'\n/\n")
        deck=f"""&MDATA
qScaled=.true.
qOneD=.false.
qFieldEvolve={'.true.' if evolve else '.false.'}
qWriteA={'.true.' if evolve else '.false.'}
qElectronsEvolve={'.true.' if evolve else '.false.'}
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
iNumNodesX=81
iNumNodesY=81
nodesPerLambdar=21
sFModelLengthX={.0032/trans:.17g}
sFModelLengthY={.0032/trans:.17g}
sFModelLengthZ2={(nz_fast-1-1e-7)*lam/20/lc:.17g}
iRedNodesX=49
iRedNodesY=49
sFiltFrac=.3
sDiffFrac=.25
sBeta=1
sElectronThreshold=0
srho={scale['rho']:.17g}
saw={aw:.17g}
sgamma_r={gamma_ref:.17g}
lambda_w=.028
sux=0
suy=1
zundType='aprl-ramp77'
stepsPerPeriod={steps}
nPeriods={nperiods}
sZ0=0
iWriteNthSteps={nperiods*steps}
iWriteIntNthSteps={steps}
sRedistLen=.5
iRedistStp=20
ioutInfo=2
/
"""
        (d/'main.in').write_text(deck)
    print(json.dumps({k:v for k,v in native.items() if k not in ('limits','scale')},indent=2))


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve();source_copy(root);input_file(root)
