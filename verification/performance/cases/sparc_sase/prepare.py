"""Prepare immutable input-generation and preflight allocation; no simulation."""
import argparse,shutil
from pathlib import Path
import h5py,numpy as np
from common import *
import common

def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path)
    ap.add_argument('--reuse-input-and-tests',type=Path)
    ap.add_argument('--reuse-tests',type=Path)
    ap.add_argument('--slippage-multiple',type=float,choices=[2.,3.])
    ap.add_argument('--transverse-cell-um',type=float)
    ap.add_argument('--transverse-cells',type=int)
    ap.add_argument('--nodes',type=int,choices=[1,2]);a=ap.parse_args()
    short=a.slippage_multiple is not None
    settings=dict(slippage_multiple=a.slippage_multiple,
      transverse_cell_um=a.transverse_cell_um or (40 if short else 20),
      transverse_cells=a.transverse_cells or (80 if short else 128),
      nodes=a.nodes or (1 if short else 2))
    common.configure(settings)
    globals().update({k:getattr(common,k) for k in ['T','Q','IP','DX','DY','NX','NY','EXEC_NODES']})
    root=a.output.resolve();root.mkdir(exist_ok=False)
    stage_scripts(root,settings)
    native=root/'native-noise';native.mkdir()
    currents=slice_profile();n=len(currents);rounded=int(np.ceil(n/32)*32)
    padded=np.pad(currents,(0,rounded-n),constant_values=1e-12)
    with h5py.File(native/'current.h5','x') as f:
        f['s']=np.arange(rounded)*LAMBDA
        f['current']=np.maximum(padded,1e-12)
        f['physical_current']=np.pad(currents,(0,rounded-n))
    (native/'initial.lat').write_text('D: DRIFT={l=.001};\nFEL: LINE={D};\n')
    # Native, unchanged GENESIS loading supplies the physical Fawley noise.
    bx=132e-6**2*GAMMA/2.5e-6;by=75e-6**2*GAMMA/2.9e-6
    (native/'initial.in').write_text(f'''&setup
rootname=initial
lattice=initial.lat
beamline=FEL
lambda0={LAMBDA:.17g}
gamma0={GAMMA:.17g}
seed=20261103
shotnoise=true
npart={NPART}
nbins={NBINS}
&end
&time
s0=0
slen={rounded*LAMBDA:.17g}
sample=1
time=true
&end
&profile_file
label=current
xdata=current.h5/s
ydata=current.h5/current
&end
&beam
gamma={GAMMA:.17g}
current=@current
delgam={GAMMA*.0002:.17g}
ex=2.5e-6
ey=2.9e-6
betax={bx:.17g}
betay={by:.17g}
alphax=.938
alphay=-.705
&end
&write
beam=initial
&end
''')
    main_lattice(root/'genesis')
    binary=WS/'build/sparc-sase-r1-20261103/aprl'
    need(binary.exists(),'Missing new source binary')
    plan=dict(workspace=str(WS),root=str(root),binary=str(binary),binary_sha256=sha(binary),
      genesis=str(GENESIS),genesis_sha256=sha(GENESIS),
      seed=20261103,NPART=NPART,NBINS=NBINS,active_slices=n,generated_slices=rounded,
      physics=dict(kinetic_energy_MeV=151.9,gamma=GAMMA,boost_gamma=GBOOST,
       period_m=LU,K=K,wavelength_m=LAMBDA,charge_C=Q,duration_s=T,peak_current_A=IP,
       bunch_length_m=c*T,total_nominal_slippage_m=TOTAL_SLIP,slippage_multiple=a.slippage_multiple,
       reference='Published long-bunch SPARC' if not short else 'Short-bunch derivative preserving peak current, emittances, energy spread and magnetic lattice'),
      execution=dict(nodes=EXEC_NODES,ranks=16*EXEC_NODES,threads_per_rank=8),
      launch_stage='Input, noise audit, required regressions and bounded full-geometry timing',
      production_methods=['boris','energy','momentum'],first_campaign_noise_seeds=1,
      limitations=['Single noise seed cannot give a 5 percent SASE statistical acceptance',
       'Reconstructed common magnetic ends and phase matching are explicitly recorded',
       'Raw total EB; no automatic near-field subtraction'])
    if a.reuse_tests:
        previous=a.reuse_tests.resolve()
        need(load(previous/'required.log.launch.json')['exit_code']==0,'Previous required tests did not pass')
        need(load(previous/'plan.json')['binary_sha256']==plan['binary_sha256'],'Tests used a different binary')
        plan['reuse_tests']=str(previous)
    if a.reuse_input_and_tests:
        previous=a.reuse_input_and_tests.resolve()
        need(load(previous/'required.log.launch.json')['exit_code']==0,'Previous required tests did not pass')
        need(load(previous/'native-noise/initial.log.launch.json')['exit_code']==0,'Native input generation incomplete')
        need(load(previous/'plan.json')['binary_sha256']==plan['binary_sha256'],'Previous tests used different binary')
        need((native/'initial.in').read_bytes()==(previous/'native-noise/initial.in').read_bytes(),'Native input deck changed')
        with h5py.File(native/'current.h5','r') as f,h5py.File(previous/'native-noise/current.h5','r') as g:
            need(all(np.array_equal(f[k][:],g[k][:]) for k in f),'Current profile changed')
        (native/'initial.par.h5').hardlink_to(previous/'native-noise/initial.par.h5')
        plan['reuse_input_and_tests']=str(previous)
        plan['reused_native_input_sha256']=sha(native/'initial.par.h5')
    save(root/'plan.json',plan)
    save(root/'preparation_plan.json',plan)
    files=[root/'preparation_plan.json',root/'benchmark.json',binary,GENESIS,*[p for p in (root/'scripts').iterdir() if p.is_file()],
           *native.glob('*'),root/'genesis/reference.lat',root/'genesis/lattice.json']
    save(root/'PREPARED.json',{str(p):sha(p) for p in files})
    print(root)

if __name__=='__main__':main()
