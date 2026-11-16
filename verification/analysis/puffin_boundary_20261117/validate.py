"""Independent checks of the calculation and the portable explanation."""
import json
import re
import sys
from pathlib import Path
import numpy as np
from prepare import sha,save

def main(root):
    root=Path(root);r=json.loads((root/'RESULT.json').read_text())
    n=np.load(root/'data/operator_illustration.npz')
    rows=np.genfromtxt(root/'data/prescribed_position_scan.csv',names=True,delimiter=',')
    a=n['input_A'];b=n['output_A'];u=np.sum(abs(a)**2,axis=(1,2));v=np.sum(abs(b)**2,axis=(1,2))
    err=100*np.sqrt(np.sum(abs(b-a)**2,axis=(1,2))/u)
    assert np.allclose(err,rows['complex_field_L2_difference_percent'],rtol=1e-13,atol=1e-13)
    assert np.allclose(100*(v/u-1),rows['energy_proxy_change_percent'],rtol=1e-13,atol=1e-12)
    delta=b-a;identity=np.sum(2*np.real(a.conj()*delta)+abs(delta)**2,axis=(1,2))
    identity_error=float(np.max(abs((v-u)-identity)/u));assert identity_error<1e-14
    # Independent explicit DFT matrices, rather than a second FFT calculation.
    size=81;j=np.arange(size)
    f=np.exp(-2j*np.pi*j[:,None]*j[None,:]/size)
    g=r['grid'];k=2*np.pi*np.fft.fftfreq(size,d=g['dx_bar'])
    d=np.exp(-g['h_bar']*(abs(k[:,None])+abs(k[None,:]))/np.sqrt(2*g['carrier_kz2_bar']))
    selected=128;masked=n['mask']*a[selected]
    explicit=(1-n['mask'])*a[selected]+f.conj()@(d*(f@masked@f.T))@f.conj().T/size**2
    dft_error=float(np.linalg.norm(explicit-b[selected])/np.linalg.norm(b[selected]));assert dft_error<1e-12
    native=np.genfromtxt(root/'data/native_tail_masks.csv',names=True,delimiter=',')
    native=native[native['nz']==96000]['mask_native']
    rank=np.genfromtxt(root/'data/native_tail_actual_rank.csv',names=True,delimiter=',')['mask_native']
    rank_error=float(np.max(abs(native-rank)));assert rank_error<1e-7
    assert r['validations']['native_transverse_vs_independent_formula_max']<1e-13
    assert r['phase_only_example']['complex_L2_difference']>0.999999999999
    assert abs(r['phase_only_example']['norm2_ratio']-1)<1e-14
    for p,h in r['sources'].items():assert sha(p)==h,p
    links=[]
    for p in root.rglob('*.md'):
        if 'verification' in p.parts:continue
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('http:','https:','#')):continue
            if (p.parent/target).resolve()==(root/'VALIDATION.json').resolve():
                links.append([str(p.relative_to(root)),target]);continue
            assert (p.parent/target).exists(),(p,target)
            links.append([str(p.relative_to(root)),target])
    rendering=json.loads((root/'RENDER_VALIDATION.json').read_text())
    for p,h in rendering['sources'].items():assert sha(root/p)==h,p
    for v in rendering['figures']+rendering['videos']:assert sha(root/v['path'])==v['sha256']
    save(root/'VALIDATION.json',dict(complete=True,
       native_mask_called_from_existing_production_library=True,
       native_transverse_vs_literal_independent_formula_max=r['validations']['native_transverse_vs_independent_formula_max'],
       actual_128_rank_last_slab_vs_37_node_function_evaluation_max=rank_error,
       independent_explicit_DFT_vs_FFT_relative_L2=dft_error,
       all_141_field_and_norm2_table_rows_recomputed=True,
       energy_interference_identity_max_relative=identity_error,
       phase_only_equal_norm_example_checked=True,
       source_and_binary_sha256_unchanged=True,local_markdown_links=len(links),
       rendered_source_hashes_verified=True,
       full_Puffin_reflection_or_FEL_impact_validated=False,
       no_scientific_H5_read=True,no_particle_or_along_z_evolution=True))
    print((root/'VALIDATION.json').read_text())

if __name__=='__main__':main(sys.argv[1])
