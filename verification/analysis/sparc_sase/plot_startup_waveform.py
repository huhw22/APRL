"""Display observed energy-case phase waveform; no new dynamics or momentum data."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def main(top, out):
    top,out=Path(top),Path(out);out.mkdir(exist_ok=True)
    old=top/'energy_startup_origin_20261105/r2/analysis'
    independent=top/'momentum_startup_assessment_20261105/independent'
    fits=json.loads((independent/'energy_waveform_check.json').read_text())
    rows=list(csv.DictReader((independent/'negativecontrol_profiles.csv').open()))
    selected={name:[r for r in rows if float(r['origin_shift_rad'])==0 and r['profile']==name]
              for name in ['all_particle','within_family']}
    original=list(csv.DictReader((old/'phase_bins.csv').open()))
    theta=np.array([float(r['phase_rad']) for r in original])
    mean=np.array([float(r['phase_within_mean_rad']) for r in original])
    a0,a=fits['phase_within_mean_rad']['linear_with_wrap_at_zero']['coefficients']
    h=np.array(fits['phase_within_mean_rad']['first_harmonic']['coefficients'])
    ramp=a0+a*(theta-np.pi)
    sine=h[0]+h[1]*np.cos(theta)+h[2]*np.sin(theta)
    counter=fits['phase_ramp_counterfactual']
    actual=json.loads((old/'RESULT.json').read_text())['p1']['bunching']['actual']
    values=[counter['reference']['local10_rms'],
            counter['reference_plus_fitted_periodic_ramp']['local10_rms'],actual['local_rms']]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.6,4.9),layout='constrained')
    axes[0].plot(theta,mean,'o',ms=3.5,color='#24648c',label='Observed within-family mean')
    full=selected['all_particle']
    axes[0].plot([float(r['phase_mid_rad']) for r in full],
                 [float(r['mean_phase_residual_rad']) for r in full],'-',lw=1.3,
                 color='#555555',label='Full residual (no family centering)')
    axes[0].plot(theta,ramp,'--',color='#b18428',lw=1.8,label='Periodic-ramp profile fit')
    axes[0].plot(theta,sine,':',color='#8e4b75',lw=1.8,label='Single-harmonic profile fit')
    axes[0].set_xlabel('Undulator-entrance optical phase (rad)')
    axes[0].set_ylabel('Mean arrival-phase residual (rad)')
    axes[0].set_title('(a) Observed phase waveform')
    axes[0].set_xticks([0,np.pi,2*np.pi],['0',r'$\pi$',r'$2\pi$'])
    axes[0].legend(fontsize=8.5)
    axes[1].plot(range(3),np.array(values)*1e3,'o',ms=8,color='#24648c',ls='none')
    axes[1].set_xticks(range(3),['No-feedback\nreference','Observed-ramp\nreconstruction','Existing PIC\nrecords'])
    axes[1].set_ylabel(r'Local RMS BF ($10^{-3}$)')
    axes[1].set_ylim(0,6.1)
    axes[1].set_xlim(-.3,2.3)
    for x,y in enumerate(values):
        axes[1].annotate(f'{y*1e3:.3f}',(x,y*1e3),xytext=(0,10),textcoords='offset points',ha='center')
    axes[1].set_title('(b) Phase-shape reconstruction; no force prediction')
    for ax in axes:ax.grid(alpha=.18)
    fig.suptitle('Existing energy-coupled data: optical-period ramp as a startup clue',fontsize=13)
    fig.savefig(out/'observed_phase_ramp.png',dpi=220,bbox_inches='tight');plt.close(fig)
    (out/'README.md').write_text('# 已有能量版：周期相位斜坡线索\n\n'
        '![周期相位斜坡](observed_phase_ramp.png)\n\n'
        '左图为已有记录的条件均值与拟合；右图把从结果拟合的斜坡施加于同一个无反馈参考。'
        '它复现了BF的形状和量级，但不是独立增长预测，不证明哪一种场或数值环节产生斜坡。'
        '未读取动量版H5，也未运行新仿真。\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('top');p.add_argument('out')
    a=p.parse_args();main(a.top,a.out)
