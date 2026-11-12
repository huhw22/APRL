"""Figures and small summaries from closed source/field/propagation caches."""
import csv
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import read,save


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))


def main(root):
    root=Path(root);control=root/'physical-route-control-r1'
    if not read(control/'COMPLETE.json')['complete']:raise RuntimeError('Control incomplete')
    stages=read(control/'analysis/MEASUREMENT_RESULT.json')['stages']
    times=read(control/'phase-analysis/TIMING.json');interp=read(control/'phase-analysis/INTERPOLATION.json')
    sources=rows(control/'analysis/source_summary.csv')
    pick=lambda program:next(x for x in sources if x['plane']=='u07' and x['program']==program)
    p=pick('Puffin');a=pick('PIC')
    summary=dict(complete=True,new_seeds=0,new_PIC_runs=0,one_corrected_Puffin_run=19770,
        original_Puffin_route_not_a_matched_benchmark=True,
        source_CIC_Jx_band_norm_ratio=float(p['CIC_Jx_band_norm'])/float(a['CIC_Jx_band_norm']),
        source_CIC_inverse_uz_band_norm_ratio=float(p['CIC_quiver_density_band_norm'])/float(a['CIC_quiver_density_band_norm']),
        field_measurement_stages=stages,timing=times,interpolation=interp,physical_acceptance=False,
        conclusions=['Old single-undulator route did not preserve the requested custom type; explicit lattice and model-ends flag are required.',
                     'Native field-node power and native interpolation at the same physical points are different observables on these coarse meshes.',
                     'Physical source maps nearly coincide; off-axis Cowan propagation differs and was verified against the real C++ field kernels.',
                     'No unique percentage attribution for the remaining coupled phase/intensity difference has been established.'])
    save(root/'RESULT.json',summary)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    old=read(root/'analysis/MEASUREMENT_RESULT.json')['stages'][0]['band_energy_J']
    values=[old,stages[0]['band_energy_J'],interp['bilinear_native_Puffin_energy_J'],interp['PIC_energy_J']]
    bars=axes[0].bar(['Old Puffin\nroute','Matched\nnodes','Matched\ninterpolation','PIC\npixels'],np.array(values)*1e12,
                     color=['#999999','#32936f','#7cbfa4','#1769a5'])
    axes[0].bar_label(bars,fmt='%.2f',fontsize=8,padding=3)
    axes[0].set(ylabel='Rect ±5% energy (pJ)',title='Route and readout are separate stages',ylim=(0,38))
    old_times=read(root/'phase-analysis/TIMING.json')
    values=1e15*np.array([times['on_axis_external_field_exact_reference_delay_s'],times['delays_s']['PIC'],
                         old_times['delays_s']['Puffin'],times['delays_s']['Puffin']])
    bars=axes[1].bar(['On-axis\nreference','PIC','Old Puffin','Correct-route\nPuffin'],values,color=['#555555','#1769a5','#999999','#32936f'])
    axes[1].bar_label(bars,fmt='%.3f',fontsize=8,padding=3)
    axes[1].set(ylabel='Mean radiation slippage (fs)',title='Actual particle arrival time',ylim=(10.5,11.75))
    disp=rows(root/'analysis/vacuum_dispersion.csv')
    for axis,color in [('x','#1769a5'),('y','#c46b10'),('diagonal','#8c63a9')]:
        rr=[x for x in disp if x['axis']==axis and float(x['k_transverse_per_m'])<=40000]
        axes[2].plot([float(x['k_transverse_per_m'])/1e4 for x in rr],
                     [float(x['PIC_minus_exact_phase_7period_rad']) for x in rr],label='Cowan '+axis,color=color)
    rr=[x for x in disp if x['axis']=='y' and float(x['k_transverse_per_m'])<=40000]
    axes[2].plot([float(x['k_transverse_per_m'])/1e4 for x in rr],
                 [float(x['Puffin_paraxial_minus_exact_phase_7period_rad']) for x in rr],color='#32936f',ls='--',label='Puffin paraxial')
    axes[2].set(xlabel=r'$k_\perp$ ($10^4$ m$^{-1}$)',ylabel='Phase error after 0.196 m (rad)',title='Homogeneous vacuum propagation')
    axes[2].legend(fontsize=8);axes[2].grid(alpha=.2)
    fig.suptitle('One w64-a sample: current source, readout and propagation audit',fontsize=14)
    fig.supxlabel('No amplitude fitting. Native interpolation changes the observable; the on-axis reference is not an ensemble prediction.',fontsize=9)
    fig.savefig(root/'stage_diagnosis.png',dpi=180);plt.close(fig)
    with np.load(control/'analysis/source-PIC-u07.npz') as f:
        xx=f['xy_m'];s1=f['CIC_maps_C'][64,:,:,2];dd=xx[1]-xx[0]
    with np.load(control/'analysis/source-Puffin-u07.npz') as f:s2=f['CIC_maps_C'][64,:,:,2]
    with np.load(control/'analysis/fields.npz') as f:
        freq=f['frequency_Hz'];pic=f['PIC'];xp=f['PIC_x_m'];yp=f['PIC_y_m'];puf=f['Puffin'];x=f['x_m'];y=f['y_m']
    with np.load(control/'phase-analysis/bilinear_field.npz') as f:pb=f['Puffin_bilinear']
    fluence1=simpson(np.sum(abs(pic)**2,axis=3),x=freq,axis=0)*1e6
    fluence2=simpson(np.sum(abs(pb)**2,axis=3),x=freq,axis=0)*1e6
    use=np.flatnonzero(abs(xp)<.00056-1e-12);uy=np.flatnonzero(abs(yp)<.00056-1e-12)
    fluences=[v[np.ix_(uy,use)] for v in [fluence1,fluence2]]
    fig,ax=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    limit=max(float(abs(s1).max()),float(abs(s2).max()))/dd**2*1e9
    for i,(source,name) in enumerate([(s1,'PIC particle current'),(s2,'Puffin particle current')]):
        im=ax[0,i].imshow(abs(source)/dd**2*1e9,origin='lower',extent=[(xx[0]-dd/2)*1e6,(xx[-1]+dd/2)*1e6]*2,
                          vmin=0,vmax=limit,cmap='magma')
        ax[0,i].set(xlim=(-560,560),ylim=(-560,560),xlabel='x (µm)',ylabel='y (µm)',title=name+' at carrier')
    fig.colorbar(im,ax=ax[0,:],label=r'$|\widetilde J_x|$ (nC m$^{-2}$), same audit CIC basis',shrink=.85)
    limit=max(v.max() for v in fluences)
    for i,(v,name) in enumerate(zip(fluences,['PIC detector','Puffin native interpolation'])):
        im=ax[1,i].imshow(v,origin='lower',extent=[-560,560,-560,560],vmin=0,vmax=limit,cmap='magma')
        ax[1,i].set(xlabel='x (µm)',ylabel='y (µm)',title=name+' | ±5% band')
    fig.colorbar(im,ax=ax[1,:],label='Band fluence (µJ m$^{-2}$)',shrink=.85)
    fig.suptitle('Same physical sample: source agreement does not guarantee the same propagated field',fontsize=12)
    fig.supxlabel('Top: common source representation, not radiation intensity. Bottom: native field evaluation at identical physical pixels.',fontsize=9)
    fig.savefig(root/'transverse_source_and_field.png',dpi=180);plt.close(fig)
    # Whole stored transverse aperture, not the cropped comparison aperture.
    ft=np.fft.fft2(puf,axes=(1,2));dd=x[1]-x[0];kk=2*np.pi*np.fft.fftfreq(len(x),dd)
    kt=np.hypot(kk[:,None],kk[None,:]);modal=simpson(np.sum(abs(ft)**2,axis=3),x=freq,axis=0)
    save(root/'ANGULAR_CONTENT.json',dict(scope='Whole stored Puffin transverse field, homogeneous Fourier diagnostic; not aperture-limited radiated power.',
        energy_fraction_above_k={str(k):float(modal[kt>k].sum()/modal.sum()) for k in [10000,20000,30000,40000]}))


if __name__=='__main__':main(sys.argv[1])
