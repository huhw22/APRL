"""Radiation reference with G4 head-index direction matched to PIC arrivals."""
import sys,json,shutil,hashlib,time
from pathlib import Path
import h5py,numpy as np
from scipy.constants import c

def prepare(root,out):
    root=Path(root).resolve();out=Path(out).resolve();out.mkdir(exist_ok=False)
    plan=json.loads((root/'plan.json').read_text());n=plan['active_slices'];pad=512
    old=root/'genesis/common.par.h5';new=out/'common.par.h5';shutil.copy2(old,new)
    with h5py.File(new,'r+') as f:
        for i in range(n):f.move(f'slice{pad+i+1:06d}',f'temporary{i:06d}')
        for i in range(n):f.move(f'temporary{i:06d}',f'slice{pad+n-i:06d}')
    with h5py.File(old,'r') as a,h5py.File(new,'r') as b:
        for i in range(n):
            ga=a[f'slice{pad+i+1:06d}'];gb=b[f'slice{pad+n-i:06d}']
            assert set(ga)==set(gb)
            for key in ga:assert np.array_equal(ga[key][:],gb[key][:])
    (out/'reference.lat').symlink_to(root/'genesis/reference.lat')
    shutil.copy2(root/'genesis/reference.in',out/'reference.in')
    proof=dict(active_slices=n,physical_particles_unchanged=True,
      all_active_slice_datasets_bitwise_verified=True,phase_not_conjugated=True,
      input_change='Reverse active slice groups only; same optical phases, transverse phase space, per-particle energy and charge',
      source_direction_proof='Control::applySlippage sends old last slice to next rank and sets first=old last; physical light advances toward larger index, so fixed-plane arrival order is reversed',
      original_input_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
      aligned_input_sha256=hashlib.sha256(new.read_bytes()).hexdigest())
    (out/'INPUT_DIRECTION_PROOF.json').write_text(json.dumps(proof,indent=2)+'\n')

def analyze(root,review,reference,repo,output_name='radiation_aligned_r1'):
    root=Path(root);review=Path(review);reference=Path(reference)
    assert (reference/'TRACK_COMPLETE').read_text().strip()=='0'
    # Only committed cache outputs from the already closed PIC case are used.
    for _ in range(40):
        if (review/'pic.json').exists():break
        time.sleep(5)
    assert (review/'pic.json').exists(),'PIC radiation reduction not available'
    plan=json.loads((root/'plan.json').read_text());pic=json.loads((review/'pic.json').read_text())
    out=review/output_name;out.mkdir(exist_ok=False)
    for p in pic:
        for suffix in ['-spectra.npz','-summary.json']:(out/(p['plane']+suffix)).hardlink_to(review/(p['plane']+suffix))
    shutil.copy2(review/'pic.json',out/'pic.json')
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    from analyze_explicit import genesis
    gen=genesis(root,out,plan,pic,reference)
    from plot_radiation import main
    main(out,root)
    (out/'COMPLETE.json').write_text(json.dumps(dict(complete=True,physical_acceptance=False,
      scope='Radiation power only; corrected G4 slice direction; existing PIC spectra reused'),indent=2)+'\n')
if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare(sys.argv[2],sys.argv[3])
    else:analyze(*sys.argv[2:])
