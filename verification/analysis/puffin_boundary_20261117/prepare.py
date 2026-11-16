"""Source-backed boundary calculations; no particle or propagation run."""
import csv
import hashlib
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path
import numpy as np

W = Path('/public/home/wfollowm/APRL_Workspace')
P = W/'results/puffin_eleven_four_slip_20261113_r2'
B = W/'build/puffin-eleven-four-slip-20261113-r2/puffin'
E = W/'results/fel_independent_w355_20261114_r1/eleven-single-node-r1'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(1048576),b''): h.update(x)
    return h.hexdigest()

def save(p,x):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def table(p,rows):
    with Path(p).open('w',newline='') as f:
        z=csv.DictWriter(f,fieldnames=list(rows[0]));z.writeheader();z.writerows(rows)

def prescribed_case(x0,x,y,mask,d):
    X,Y=np.meshgrid(x,y)
    # Independent prescribed test fields, with intensity sigma = 80 um.
    a=np.exp(-((X-x0)**2+Y**2)/(4*(80e-6)**2))*np.exp(1j*20000*X)
    out=(1-mask)*a+np.fft.ifft2(d*np.fft.fft2(mask*a))
    u=float(np.sum(abs(a)**2));v=float(np.sum(abs(out)**2))
    return a,out,dict(x0_um=x0*1e6,input_norm2=u,output_norm2=v,
        energy_proxy_change_percent=100*(v/u-1),
        complex_field_L2_difference_percent=100*float(np.linalg.norm(out-a)/np.linalg.norm(a)),
        boundary_input_norm2_fraction=float(np.sum(abs(mask*a)**2)/u))

def main(out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    assert not (out/'RESULT.json').exists(), 'Refuse to overwrite a completed calculation'
    data=out/'data';data.mkdir(exist_ok=True);work=out/'verification';work.mkdir(exist_ok=True)
    template=Path(__file__).resolve().parent
    scale=json.loads((E/'INPUT_VALIDATED.json').read_text())['scale']
    paths=[P/'source/puffin/lib/diffraction/masks.f90',
           P/'source/puffin/lib/diffraction/diffraction.f90',
           P/'source/puffin/lib/diffraction/transforms.f90',
           P/'source/puffin/lib/setup/setup_calcs.f90',
           P/'source/puffin/lib/common/puffin_constants.f90',
           P/'source/puffin/lib/utilities/base_functions.f90',
           P/'SOURCE_METHOD_AUDIT.json',E/'INPUT_VALIDATED.json',E/'main/main.in',
           B/'libpuffin_mpi.a', B/'puffin', B/'CMakeFiles/puffin.dir/link.txt']
    before={str(p):sha(p) for p in paths}
    frozen=json.loads((P/'FROZEN_PREPARATION.json').read_text())
    for p in paths[:6]: assert frozen[str(p)]==before[str(p)]
    assert before[str(B/'puffin')]=='c4ecf60611e86628f7ce77b90c752f49dddf6daedb8ca2877ee1783ad5c5c54e'
    dz2=4*scale['pi']*scale['rho']/20
    dxbar=40e-6/scale['transverse_scale_m']
    # Reuse the exact library and link dependencies, replacing only its main.
    original=shlex.split((B/'CMakeFiles/puffin.dir/link.txt').read_text())
    deps=original[next(i for i,s in enumerate(original) if s.startswith('-Wl,-rpath')):]
    deps=[str(B/'libpuffin_mpi.a') if s=='libpuffin_mpi.a' else s for s in deps]
    cmd=[original[0],'-ffree-line-length-none','-fopenmp','-O0','-I'+str(B),
         str(template/'mask_values.f90'),'-o',str(work/'mask_values')]+deps
    log=subprocess.run(cmd,capture_output=True,text=True,check=False)
    (work/'compile.log').write_text(log.stdout+log.stderr)
    assert log.returncode==0,log.stderr
    subprocess.run([str(work/'mask_values'),str(data),repr(dz2),repr(dxbar)],check=True)
    native=np.genfromtxt(data/'native_transverse_mask.csv',names=True,delimiter=',')
    mask=native['mask_native'].reshape(81,81)
    # Independent evaluation of the transverse formula, preserving native pi.
    one=np.zeros(81);one[:16]=np.cos(scale['pi']*np.arange(16)/(2*15))**2
    one[-16:]=one[:16][::-1]
    # Preserve the actual Maskf4 formula, including the native promoted pi.
    envelope=np.ones(81)
    envelope[:16]=np.cos(scale['pi']*(np.arange(16)-15)/(2*15))**2
    envelope[-16:]=np.cos(scale['pi']*np.arange(16)/(2*15))**2
    expected=one[None,:]+envelope[None,:]*one[:,None]
    transverse_error=float(np.max(abs(mask-expected)))
    assert transverse_error<1e-13
    tails=np.genfromtxt(data/'native_tail_masks.csv',names=True,delimiter=',')
    tailrows=[];tailcheck=0.
    for row in tails:
        n=int(row['nz']);idx=int(row['node_zero_based']);x=idx*dz2;lb=36*dz2
        actual=float(row['mask_native'])
        bad=np.cos(scale['pi']*(x-(n-dz2-lb/2))/lb)**2
        tailcheck=max(tailcheck,abs(actual-bad))
        geometry=np.cos(scale['pi']*(x-((n-1)*dz2-lb/2))/lb)**2
        tailrows.append(dict(nz=n,node_zero_based=idx,
             distance_from_tail_um=(idx-(n-1))*dz2*scale['cooperation_length_m']*1e6,
             mask_native=actual,geometric_parentheses_reference=float(geometry)))
    assert tailcheck<1e-7
    table(data/'tail_mask_comparison.csv',tailrows)
    h=.25*4*scale['pi']*scale['rho']
    # GetKValues uses its own double-precision pi, not puffin_constants.pi.
    kz=2*math.pi/(20*dz2)
    modes=[]
    for k in [0,10000,20000,30000,40000]:
        d=math.exp(-h*abs(k)*scale['transverse_scale_m']/math.sqrt(2*abs(kz)))
        modes.append(dict(kx_per_m=k,ky_per_m=0,D_amplitude=d,D2_norm2=d*d,
            approximate_lab_angle_mrad=k*scale['lambda_reference_m']/(2*math.pi)*1000))
    table(data/'modal_factors.csv',modes)
    x=(np.arange(81)-40)*40e-6;y=x.copy()
    kk=2*np.pi*np.fft.fftfreq(81,d=dxbar)
    d=np.exp(-h*(abs(kk[None,:])+abs(kk[:,None]))/math.sqrt(2*abs(kz)))
    rows=[];ain=[];aout=[]
    for x0 in np.linspace(0,1.4e-3,141):
        a,b,row=prescribed_case(x0,x,y,mask,d);rows.append(row);ain.append(a);aout.append(b)
    table(data/'prescribed_position_scan.csv',rows)
    np.savez_compressed(data/'operator_illustration.npz',x_m=x,y_m=y,mask=mask,
         x0_um=[r['x0_um'] for r in rows],input_A=np.array(ain),output_A=np.array(aout))
    selected=[rows[i] for i in [0,96,128,140]]
    samesum=np.linalg.norm(ain[0])**2
    shifted=np.asarray(ain[0])*np.exp(1j*np.pi/3)
    phase=dict(phase_radians=math.pi/3,
       norm2_ratio=float(np.linalg.norm(shifted)**2/samesum),
       complex_L2_difference=float(np.linalg.norm(shifted-ain[0])/np.linalg.norm(ain[0])),
       scope='Exact mathematical phase example, not a measured Puffin field error.')
    samples=[r for r in tailrows if r['nz']==96000]
    result=dict(complete=True,date_local='2026-11-17',
        scope='Native mask-function evaluation and prescribed single-operator formula illustration; no Puffin particle or along-z evolution.',
        production_case='19799, eleven segments, w about 355',
        grid=dict(nx=81,ny=81,nz2=96000,dx_m=40e-6,nbx=16,nby=16,nbz2=37,
             native_pi=scale['pi'],dz2_bar=dz2,dx_bar=dxbar,h_bar=h,carrier_kz2_bar=kz,
             step_physical_m=h*scale['gain_length_m'],beta=1,
             node_extent_um=[-1600,1600],transverse_fft_period_um=3240,
             transverse_mask_support_inner_um=1000),
        modal_factors=modes, prescribed_cases=selected,phase_only_example=phase,
        tail_mask=dict(native_first=samples[0]['mask_native'],native_last=samples[-1]['mask_native'],
          max_native_vs_parentheses_difference=max(abs(r['mask_native']-r['geometric_parentheses_reference']) for r in samples),
          physical_width_um=36*dz2*scale['cooperation_length_m']*1e6,
          total_bar_length_native=96000-dz2,total_bar_length_geometric=95999*dz2,
          production_impact_measured=False,reference_is_accepted_fix=False),
        validations=dict(native_transverse_vs_independent_formula_max=transverse_error,
          native_tail_vs_independent_formula_max=tailcheck,
          existing_source_and_binary_hashes_unchanged=all(sha(Path(k))==v for k,v in before.items())),
        sources=before,scientific_h5_read=False,new_physical_simulation=False,
        puffin_reflection_coefficient_measured=False,
        pulse_test=dict(intensity_sigma_um=80,kx_per_m=20000,ky_central_per_m=0,
            single_frequency=True,longitudinal_boundary_omitted=True,
            x0_is_prescribed_parameter_not_propagation_coordinate=True))
    save(out/'RESULT.json',result)
    save(work/'NATIVE_READER.json',dict(command=cmd,executable_sha256=sha(work/'mask_values'),
          calls=['getMask','getZ2Mask'],no_MPI_initialization=True,no_physical_program_rebuild=True))
    # Human-readable, distributable source excerpts; full source and binaries stay out of ZIP.
    evidence=out/'evidence';evidence.mkdir(exist_ok=True)
    snippets=[('masks.f90',P/'source/puffin/lib/diffraction/masks.f90',[(135,189),(260,435)]),
              ('diffraction.f90',P/'source/puffin/lib/diffraction/diffraction.f90',[(128,171),(431,522)]),
              ('setup_calcs.f90',P/'source/puffin/lib/setup/setup_calcs.f90',[(237,247),(846,858)])]
    for name,path,ranges in snippets:
        lines=path.read_text().splitlines()
        text='Source: '+str(path)+'\nSHA256: '+sha(path)+'\n\n'
        for lo,hi in ranges:
            text+='\n'.join(f'{i}: {lines[i-1]}' for i in range(lo,min(hi,len(lines))+1))+'\n\n'
        (evidence/(name+'.txt')).write_text(text)
    save(evidence/'SOURCE_LINKS.json',dict(official_documentation='https://ukfels.github.io/puffinDocs/',
       official_masks='https://github.com/UKFELs/Puffin/blob/master/puffin/lib/diffraction/masks.f90',
       official_diffraction='https://github.com/UKFELs/Puffin/blob/master/puffin/lib/diffraction/diffraction.f90',
       authority='Numerical identities use the local frozen production library; current public pages are supporting source context.'))
    print(json.dumps({k:result[k] for k in ['modal_factors','prescribed_cases','tail_mask','validations']},ensure_ascii=False,indent=2))

if __name__=='__main__': main(sys.argv[1])
