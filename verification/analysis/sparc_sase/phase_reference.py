"""Derive nominal phase closure from the actual GENESIS ODE, not FEL power."""
import sys,json,re,shutil,hashlib
from pathlib import Path
import h5py,numpy as np

def wrap(x):return float(np.angle(np.exp(1j*x)))
def phase_rows(directory,plan):
    with h5py.File(Path(directory)/'reference.out.h5','r') as f:
        z=f['Lattice/z'][:];dz=f['Lattice/dz'][:];aw=f['Lattice/aw'][:]
        ku=f['Lattice/ku'][:];ps=f['Lattice/phaseshift'][:]
    gamma=plan['physics']['gamma'];ks=2*np.pi/plan['physics']['wavelength_m']
    uz=np.sqrt(gamma*gamma-1-aw*aw)
    rate=np.where(ku==0,ks/(2*gamma*gamma),ku)-ks*(1+aw*aw)/(uz*(gamma+uz))
    result=[]
    for j in range(5):
        start=.01+j*(77*.028+.4);end=start+77*.028+.4
        m=(z>=start-1e-8)&(z<end-1e-8)
        natural=float(np.sum(rate[m]*dz[m]))
        result.append(dict(gap=j,natural_phase_rad=natural,target_phi_rad=wrap(natural),
          observed_autophase_rad=float(ps[m].sum()),entrance_residual_rad=wrap(natural+ps[m].sum())))
    return result

def prepare(review,repo):
    review=Path(review).resolve();repo=Path(repo).resolve()
    root=review.parent/'compact-production-r1';plan=json.loads((root/'plan.json').read_text())
    old=review/'genesis-aligned-r1';out=review/'genesis-phase-r1';out.mkdir(exist_ok=False)
    rows=phase_rows(old,plan);text=(old/'reference.lat').read_text()
    for r in rows:
        pattern=r'(P'+str(r['gap'])+r': PHASESHIFTER=\{l=\.242,phi=)[^}]+'
        text,n=re.subn(pattern,lambda m:m[1]+format(r['target_phi_rad'],'.17g'),text)
        assert n==1
    (out/'reference.lat').write_text(text)
    shutil.copy2(old/'reference.in',out/'reference.in')
    (out/'common.par.h5').symlink_to(old/'common.par.h5')
    tools=review/'phase-tools';tools.mkdir(exist_ok=False)
    src=repo/'verification/analysis/sparc_sase'
    for name in ['phase_reference.py','phase_reference.sbatch','aligned_reference.py','analyze_explicit.py','plot_radiation.py']:
        shutil.copy2(src/name,tools/name)
    proof=dict(old_nominal_phases=rows,change='Only five finite GENESIS phase-shifter phi values; optical carrier fallback in BeamSolver drifts included',
      physical_radiation_not_used_to_choose_phi=True,input_particle_file_unchanged=True,
      PIC_not_rerun=True,old_reference_retained=True)
    (out/'PHASE_DERIVATION.json').write_text(json.dumps(proof,indent=2)+'\n')
    files=[*tools.iterdir(),out/'reference.lat',out/'reference.in',out/'common.par.h5']
    (review/'PHASE_FROZEN.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
    print(json.dumps(rows,indent=2))

def verify(review):
    review=Path(review).resolve();root=review.parent/'compact-production-r1'
    plan=json.loads((root/'plan.json').read_text());ref=review/'genesis-phase-r1'
    rows=phase_rows(ref,plan)
    assert max(abs(r['entrance_residual_rad']) for r in rows)<1e-8
    with h5py.File(review/'genesis-aligned-r1/reference.out.h5','r') as a,h5py.File(ref/'reference.out.h5','r') as b:
        for k in ['z','dz','aw','ku','qf','slippage']:
            assert np.array_equal(a['Lattice/'+k][:],b['Lattice/'+k][:]),k
    data=json.loads((review/'PHASE_FROZEN.json').read_text())
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in data.items())
    (ref/'PHASE_CLOSURE_VERIFIED.json').write_text(json.dumps(dict(pass_nominal_closure=True,rows=rows),indent=2)+'\n')
    print('Nominal phase closure verified',flush=True)
if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare(*sys.argv[2:4])
    else:verify(sys.argv[2])
