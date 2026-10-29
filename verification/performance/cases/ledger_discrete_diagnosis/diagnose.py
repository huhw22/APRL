#!/usr/bin/env python3
"""Identify independently defined spatial, flux and time terms in the ledger.

No PIC run, field change, FFT, fitted compensation or production modification.
The three-term identity is for the particle-free axial pulse where Cowan
reduces to the ordinary 1-D staggered update. It is not a 3-D Cowan invariant.
"""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import random
import h5py
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-ledger-diagnosis-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD=np.longdouble


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def exact_rational_checks():
    """Independent exact SBP algebra, with arbitrary E/B and fixed E endpoints."""
    rng=random.Random(20261029);omission_failures=0
    for case in range(64):
        n=2+case%9
        e=[F(rng.randrange(-7,8),5) for _ in range(n+1)]
        b=[F(rng.randrange(-7,8),7) for _ in range(n)]
        bp=[b[k]-(e[k+1]-e[k]) for k in range(n)]
        ep=e[:]
        for k in range(1,n):ep[k]=e[k]-(bp[k]-bp[k-1])
        def uc(x,y):return sum((((x[k]+x[k+1])/2)**2+y[k]**2)/2 for k in range(n))
        def un(x,y):return sum(((x[k]**2+x[k+1]**2)/2+y[k]**2)/2 for k in range(n))
        def p(x,y):return (x[-2]+x[-1])/2*y[-1]-(x[0]+x[1])/2*y[0]
        spatial=(uc(ep,bp)-un(ep,bp))-(uc(e,b)-un(e,b))
        variance=lambda x:-sum((x[k+1]-x[k])**2 for k in range(n))/8
        assert uc(e,b)-un(e,b)==variance(e)
        assert spatial==variance(ep)-variance(e)
        q=(p(e,b)+p(ep,bp))/2
        q_sbp=e[-1]*bp[-1]-e[0]*bp[0]
        time=(sum((ep[k]-e[k])**2 for k in range(1,n))-
              sum((bp[k]-b[k])**2 for k in range(n)))/2
        old=uc(ep,bp)-uc(e,b)+q
        assert old==spatial+(q-q_sbp)+time
        if time!=0:omission_failures+=int(old!=spatial+(q-q_sbp))
    assert omission_failures>0
    return {'pass':True,'exact_fraction_cases':64,
            'cases_that_fail_if_time_term_is_omitted':omission_failures,
            'no_floating_point_tolerance_used':True}


def pulse_terms(start,width):
    k=np.arange(65,dtype=LD)
    eps=LD('8.8541878128e-12');mu=LD('1.25663706212e-6');c=LD('299792458')
    dz=LD('1.2543262415359034e-6');dt=dz/c;area=(16*LD('5e-6'))**2
    def pulse(x):
        phase=(x-start)/width
        return np.where((phase>0)&(phase<1),LD('1e6')*np.sin(np.arccos(LD(-1))*phase)**2,LD(0))
    def uc(e,b):return area*dz/2*np.sum(eps*((e[:-1]+e[1:])/2)**2+b*b/mu)
    def un(e,b):return area*dz/2*np.sum(eps*(e[:-1]**2+e[1:]**2)/2+b*b/mu)
    def p(e,b):return area/mu*((e[-2]+e[-1])/2*b[-1]-(e[0]+e[1])/2*b[0])
    values=[]
    for n in range(28):
        e=pulse(k-2*n);b=pulse(k[:-1]+1-2*n)/c
        ep=pulse(k-1-2*n);ep[0]=e[0];ep[-1]=e[-1]
        bp=pulse(k[:-1]-2*n)/c
        spatial=(uc(ep,bp)-un(ep,bp))-(uc(e,b)-un(e,b))
        variance=lambda x:-area*dz*eps/8*np.sum(np.diff(x)**2)
        assert abs(spatial-(variance(ep)-variance(e)))<LD('1e-28')
        q=dt/2*(p(e,b)+p(ep,bp))
        q_sbp=area*dt/mu*(e[-1]*bp[-1]-e[0]*bp[0])
        time=area*dz/2*(eps*np.sum((ep[1:-1]-e[1:-1])**2)-np.sum((bp-b)**2)/mu)
        residual=uc(ep,bp)-uc(e,b)+q
        values.append((spatial,q-q_sbp,time,residual,q,q_sbp))
    return np.array(values,dtype=LD)


def real_case_spatial_check(prior):
    case=prior/'stage-replay-r2/free-drift'
    names=('center_energy_J','native_electric_energy_J','native_magnetic_energy_J')
    sums={n:np.zeros((257,8),dtype=LD) for n in names}
    paths=sorted((case/'field-stages').glob('*.csv'));assert len(paths)==16
    power_names=[f'power_{a}{side}_W' for a in 'xyz' for side in ('minus','plus')]
    powers=np.zeros((257,8,6),dtype=LD)
    for path in paths:
        with path.open() as f:rows=list(csv.DictReader(f))
        assert len(rows)==257*8
        for n in names:sums[n]+=np.array([LD(x[n]) for x in rows]).reshape(257,8)
        for a,n in enumerate(power_names):powers[:,:,a]+=np.array([LD(x[n]) for x in rows]).reshape(257,8)
    with h5py.File(case/'ledger/energy-ledger.h5') as f:
        g=f['energy_ledger'];r=g['records'][:]
        assert int(g['complete'][()])==1 and int(g['committed_records'][()])==258
    assert np.all(r['boundary_source_work_J']==0) and np.all(r['prescribed_boris_work_J']==0)
    c=sums['center_energy_J'];native=sums['native_electric_energy_J']+sums['native_magnetic_energy_J']
    dc=c[:,6]-c[:,0]-(c[:,5]-c[:,4])
    dn=native[:,6]-native[:,0]-(native[:,5]-native[:,4])
    dt=LD(str(r['time_box_s'][1]))
    q=dt/2*np.sum(powers[:,0,:]+powers[:,6,:]-(powers[:,5,:]-powers[:,4,:]),axis=1)
    source=np.diff(np.array(r['particle_grid_work_J'],dtype=LD)-np.array(r['frozen_constraint_work_J'],dtype=LD)-np.array(r['carrier_constraint_work_J'],dtype=LD))
    rc=np.r_[LD(0),np.cumsum(dc+q+source)];rn=np.r_[LD(0),np.cumsum(dn+q+source)]
    recorded=np.array(r['field_exchange_residual_J'],dtype=LD)
    assert np.max(abs(rc-recorded))<LD('1e-29')
    effects=np.r_[LD(0),np.cumsum(dc-dn)]
    result={}
    for label,i,j in [('whole',0,257),('released_no_exit',61,180),('zero_current_tail',244,257)]:
        result[label]={'start_step':i,'end_step':j,'reported_center_residual_J':float(recorded[j]-recorded[i]),
          'independently_reconstructed_center_residual_J':float(rc[j]-rc[i]),
          'native_energy_residual_with_original_flux_J':float(rn[j]-rn[i]),
          'spatial_energy_definition_effect_J':float(effects[j]-effects[i]),
          'spatial_effect_over_reported_residual':float((effects[j]-effects[i])/(recorded[j]-recorded[i]))}
    return {'checks_pass':True,'phases':result,
            'maximum_ledger_reconstruction_error_J':float(np.max(abs(rc-recorded))),
            'limit':'Native-energy variant also measures its own inlet/remap changes; it is not claimed to be an exact Cowan invariant or a corrected full 3-D ledger.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();ws=a.workspace.resolve();out=a.output.resolve();out.mkdir(exist_ok=True,parents=True)
    prior=ws/'results/energy_residual_origin_20261029';vac=prior/'vacuum-r1'
    assert json.loads((vac/'VACUUM_REVIEW.json').read_text())['checks_pass']
    rational=exact_rational_checks();results=[];frozen={}
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for start,width in ((16,8),(48,8),(8,16),(8,32)):
        path=vac/f'pulse-start-{start}-width-{width}.csv';frozen[str(path)]=sha(path)
        with path.open() as f:trace=list(csv.DictReader(f))
        actual=np.array([LD(x['residual_J']) for x in trace]);initial=LD(trace[0]['field_J'])
        terms=pulse_terms(start,width);cumulative=np.vstack([np.zeros(6,dtype=LD),np.cumsum(terms,axis=0)])
        corrected=actual-np.sum(cumulative[:,:3],axis=1)
        assert np.max(abs(corrected))/initial<LD('1e-10')
        record={'start_cell':start,'width_cells':width,'initial_field_J':float(initial),
          'reported_residual_J':float(actual[-1]),'spatial_term_J':float(cumulative[-1,0]),
          'flux_term_J':float(cumulative[-1,1]),'time_update_term_J':float(cumulative[-1,2]),
          'residual_after_independent_three_term_subtraction_J':float(corrected[-1]),
          'maximum_remaining_over_initial_field':float(np.max(abs(corrected))/initial),
          'percent_initial':{n:float(100*x/initial) for n,x in zip(('spatial','flux','time_update','reported','remaining'),
             [*cumulative[-1,:3],actual[-1],corrected[-1]])}}
        results.append(record)
        with (out/f'pulse-{start}-{width}-decomposition.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['step','actual_residual_J','spatial_J','flux_J','time_update_J','remaining_J'])
            for n in range(29):writer.writerow([n,str(actual[n])]+[str(x) for x in cumulative[n,:3]]+[str(corrected[n])])
        if (start,width)==(16,8):
            for index,label,color in [(0,'Spatial energy definition','#377eb8'),(1,'Boundary flux definition','#4daf4a'),(2,'Discrete time/update term','#ff7f00')]:
                axes[0].plot(np.arange(29),100*cumulative[:,index]/initial,label=label,color=color)
            axes[0].plot(np.arange(29),100*actual/initial,label='Reported residual',color='black',linestyle='--')
            axes[0].plot(np.arange(29),100*corrected/initial,label='After independent terms',color='#984ea3')
            bars=[record['percent_initial'][n] for n in ('spatial','flux','time_update','reported','remaining')]
            axes[1].bar(range(5),bars,color=['#377eb8','#4daf4a','#ff7f00','#555555','#984ea3'])
            axes[1].set_xticks(range(5));axes[1].set_xticklabels(['Spatial','Flux','Time/update','Reported','Remaining'],rotation=20)
            for x,v in enumerate(bars):axes[1].annotate(f'{v:+.3f}%',(x,v),xytext=(0,5 if v>=0 else -13),textcoords='offset points',ha='center',fontsize=9)
    real=real_case_spatial_check(prior)
    result={'checks_pass':True,'frame':'boost','new_pic_jobs':0,'rational_identity_checks':rational,
       'pulse_results':results,'real_case_spatial_energy_check':real,
       'identity':'R_center = Delta(Ucenter-Unative) + (Qcenter_trapezoid-Qnative_sbp) + T; inlet/remap use their separately measured changes',
       'time_update_definition':'T = epsilon0*dV/2 sum(Enew-Eold)^2 - dV/(2*mu0) sum(Bnew-Bold)^2, with native quadrature and actually held E endpoints',
       'native_boundary_flux_1d':'dt*area/mu0*(Eold[N]*Bnew[N-1]-Eold[0]*Bnew[0])',
       'limits':['The three-term identity is verified for the source-free 1-D axial case. 3-D Cowan needs its modified spatial operator and boundary handling included.',
         'T is a diagnostic term from the discrete update, not new physical work or fitted compensation. It interacts with the held sacrificial end plane.',
         'Real-case spatial effect is measured. Its time/flux/Cowan terms are not all available from the saved scalar traces.',
         'The former real-case 32% residual has not been wholly explained, and no propagation kernel was changed.']}
    (out/'DIAGNOSIS.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    frozen[str(vac/'VACUUM_REVIEW.json')]=sha(vac/'VACUUM_REVIEW.json')
    for path in (prior/'stage-replay-r2/free-drift/field-stages').glob('*.csv'):frozen[str(path)]=sha(path)
    path=prior/'stage-replay-r2/free-drift/ledger/energy-ledger.h5';frozen[str(path)]=sha(path)
    frozen[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
    (out/'INPUT_SHA256.json').write_text(json.dumps(frozen,indent=2,sort_keys=True)+'\n')
    axes[0].set_xlabel('Moving-window step');axes[0].set_ylabel('% of initial field energy')
    axes[0].set_title('Known correct 8-cell vacuum pulse');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('% of initial field energy');axes[1].set_title('Final additive decomposition')
    axes[1].set_ylim(min(bars)-.7,max(bars)+.5)
    for ax in axes:ax.axhline(0,color='#999999',linewidth=.7);ax.grid(axis='y',alpha=.2)
    fig.savefig(out/'ledger_term_decomposition.png',dpi=180);fig.savefig(out/'ledger_term_decomposition.pdf');plt.close(fig)
    print(json.dumps({'rational_checks':rational,'first_pulse':results[0],'real_case':real},indent=2))


if __name__=='__main__':main()
