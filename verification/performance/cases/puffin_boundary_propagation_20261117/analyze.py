"""Quantify absorption, wrap-around, reverse-x modes and central-field changes."""
import csv
import json
import sys
from pathlib import Path
import numpy as np
from prepare import save,sha

def table(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main(root):
    root=Path(root);p=json.loads((root/'PLAN.json').read_text())
    assert json.loads((root/'PHYSICS_COMPLETE.json').read_text())['complete']
    out=root/'analysis';out.mkdir(exist_ok=True)
    loaded={};protected={}
    for c in p['cases']:
        d=root/'cases'/c['name'];comp=json.loads((d/'COMPLETE.json').read_text());m=json.loads((d/'META.json').read_text())
        hist=np.genfromtxt(d/'history.csv',names=True,delimiter=',')
        a=np.fromfile(d/'field_xy.bin',dtype='<c16').reshape(comp['frames'],81,81)
        line=np.fromfile(d/'line_x.bin',dtype='<f8').reshape(c['steps']+1,81)
        assert np.isfinite(a).all() and np.isfinite(hist['total_norm2_m2']).all()
        loaded[c['name']]=dict(meta=m,hist=hist,a=a,line=line)
        for f in ['field_xy.bin','line_x.bin','history.csv','META.json','COMPLETE.json','LAUNCH.json']:
            path=d/f;st=path.stat();protected[str(path)]=dict(size=st.st_size,inode=st.st_ino,mtime_ns=st.st_mtime_ns)
    names=['small-native','small-no-absorption','wide-reference','wider-reference']
    ref=loaded['wide-reference'];wider=loaded['wider-reference']
    x=(np.arange(81)-40)*40e-6;dx=40e-6
    wt=np.ones(27);wt[[0,-1]]=.5;weights=wt[:,None]*wt[None,:]
    def roi(a):return a[:,27:54,27:54]
    def norm(a):return np.sum(abs(roi(a))**2*weights,axis=(1,2))*dx**2
    u0=float(ref['hist']['total_norm2_m2'][0]);zs=ref['hist']['z_m'][::5]
    for name in names:
        q=loaded[name]
        assert np.array_equal(q['hist']['z_m'],ref['hist']['z_m'])
        assert np.array_equal(q['a'][0],ref['a'][0]),name
        assert np.allclose(norm(q['a']),q['hist']['roi_norm2_m2'][::5],rtol=2e-14,atol=1e-24)
    referr=np.sqrt(norm(ref['a']-wider['a'])/u0)
    assert float(referr.max())<1e-6,('reference convergence',referr.max())
    noabs=loaded['small-no-absorption']['hist']['total_norm2_m2']
    assert float(np.max(abs(noabs/noabs[0]-1)))<1e-10
    # Independent analytic paraxial Gaussian, at the same retarded coordinate.
    X,Y=np.meshgrid(x,x);sig=200e-6;x0=-300e-6;kx=30000.
    theta=ref['meta']['theta_x_rad'];k=kx/theta
    theory=[]
    for z in zs:
        q=1+1j*z/(2*k*sig**2)
        theory.append(np.exp(1j*kx*X-1j*kx*kx*z/(2*k))/q*
                      np.exp(-((X-x0-theta*z)**2+Y*Y)/(4*sig**2*q)))
    theory=np.array(theory);theoryerr=np.sqrt(norm(ref['a']-theory)/u0)
    assert float(theoryerr.max())<1e-6,('analytic reference',theoryerr.max())
    records=[];metrics={};rows=[]
    for name in names:
        q=loaded[name];d=q['a']-ref['a'];err2=norm(d)/u0
        measured=norm(q['a']);reference=norm(ref['a'])
        cross=2*np.real(np.sum(roi(ref['a']).conj()*roi(d)*weights,axis=(1,2)))*dx**2
        closure=np.max(abs((measured-reference)-cross-norm(d))/u0)
        assert closure<1e-13
        neg=[]
        if name.startswith('small-'):
            # Use the complete native small-domain Fourier basis, not a cropped reference spectrum.
            f=np.fft.fft2(q['a'],axes=(1,2));kxgrid=np.fft.fftfreq(81,d=dx)
            neg=np.sum(abs(f[:,:,kxgrid<0])**2,axis=(1,2))*dx**2/(81*81*u0)
        else:neg=np.full(len(zs),np.nan)
        for i,z in enumerate(zs):
            rows.append(dict(case=name,step=i*5,z_m=z,total_norm2_over_initial=q['hist']['total_norm2_m2'][i*5]/u0,
                 roi_norm2_over_initial=measured[i]/u0,roi_reference_norm2_over_initial=reference[i]/u0,
                 roi_signed_norm2_difference_over_initial=(measured[i]-reference[i])/u0,
                 roi_error_field_norm2_over_initial=err2[i],roi_field_L2_difference_over_initial=err2[i]**.5,
                 negative_kx_norm2_over_initial=neg[i] if np.isfinite(neg[i]) else ''))
        late=zs>=1.3
        metrics[name]=dict(final_total_norm2_over_initial=float(q['hist']['total_norm2_m2'][-1]/u0),
             net_norm2_removed_fraction=float(1-q['hist']['total_norm2_m2'][-1]/u0),
             max_central_field_L2_change_over_initial=float(np.sqrt(err2).max()),
             max_central_error_norm2_over_initial=float(err2.max()),
             late_max_central_norm2_over_initial=float(np.max(measured[late])/u0),
             late_max_reference_central_norm2_over_initial=float(np.max(reference[late])/u0),
             interference_identity_max_relative=float(closure))
        if name.startswith('small-'):
            metrics[name]['negative_kx_initial_fraction']=float(neg[0])
            metrics[name]['negative_kx_max_fraction']=float(neg.max())
            metrics[name]['negative_kx_final_fraction']=float(neg[-1])
    table(out/'comparison.csv',rows)
    history=[]
    for name in names:
        for v in loaded[name]['hist']:
            history.append(dict(case=name,step=int(v['step']),z_m=v['z_m'],total_norm2_over_initial=v['total_norm2_m2']/u0,
                           roi_norm2_over_initial=v['roi_norm2_m2']/u0))
    table(out/'history.csv',history)
    np.savez_compressed(out/'visual_arrays.npz',x_m=x,z_m=ref['hist']['z_m'],snapshot_z_m=zs,
       reference=ref['a'],native=loaded['small-native']['a'],no_absorption=loaded['small-no-absorption']['a'],
       reference_line=ref['line'],native_line=loaded['small-native']['line'],no_absorption_line=loaded['small-no-absorption']['line'])
    result=dict(complete=True,native_frozen_DiffractionStep_invoked=True,production_FEL=False,particles=False,
       boundary_scope='Transverse only; native one-wavelength periodic longitudinal branch disables the longitudinal end mask.',
       initial_conditions=dict(sigma_intensity_um=200,x0_um=-300,kx_per_m=kx,theta_x_mrad=theta*1000),
       propagation_m=float(zs[-1]),metrics=metrics,
       nominal_geometrical_markers_m=dict(mask_start_center=(.001-x0)/theta,
           rightmost_node_center=(.0016-x0)/theta,
           wrap_reentry_at_left_ROI_center=(.00324-.00052-x0)/theta),
       validation=dict(plane_mode=json.loads((root/'cases/plane-preflight/PLANE_CHECK.json').read_text()),
         wide_vs_wider_central_field_max_relative_to_input=float(referr.max()),
         wide_vs_analytic_central_field_max_relative_to_input=float(theoryerr.max()),
         no_absorption_global_norm_drift=float(np.max(abs(noabs/noabs[0]-1))),
         all_initial_crop_samples_bitwise_equal=True),
       interpretation='Field feedback means changes of the interior radiation field, not electron feedback. Negative-kx norm is a transverse modal diagnostic, not a measured backwards-z reflection coefficient.',
       source_sha256_unchanged=all(sha(path)==h for path,h in p['source_hashes'].items()))
    assert result['source_sha256_unchanged']
    save(out/'RESULT.json',result);save(out/'PROTECTED_OUTPUTS.json',protected)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main(sys.argv[1])
