"""Native energy/momentum field identities and their Lorentz combination.

These are charges on boost-time slices, not same-lab-time volume energies.
They are kept separate from the fixed lab worldtube physical balance.
"""
import concurrent.futures
import csv
import json
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
LD=np.longdouble
SOURCES=['magnetic_supply','magnetic_cpml','magnetic_sync','electric_hold',
         'electric_cpml','electric_inlet','electric_sync','reference_current']


def load_rank(path):
    with path.open() as f:names=f.readline().strip().split(',')
    a=np.loadtxt(path,delimiter=',',skiprows=1,usecols=[0]+list(range(2,len(names))),dtype=LD)
    n=26923;assert a.shape==(3*n,len(names)-1)
    assert np.array_equal(a[:,0],np.repeat(np.arange(1,n+1),3))
    with path.open() as f:
        next(f)
        for i,line in enumerate(f):assert line.split(',',2)[1]==['begin','fixed_end','remap'][i%3]
    numeric=names[2:];stages={s:{k:a[j::3,i+1] for i,k in enumerate(numeric)} for j,s in enumerate(['begin','fixed_end','remap'])}
    assert all(np.isfinite(v).all() for s in stages.values() for v in s.values())
    start,end,remap=[stages[s] for s in ['begin','fixed_end','remap']]
    sourceH=sum(end[s+'_U_J'] for s in SOURCES)+remap['remap_H_J']
    sourceP=sum(end[s+'_P_J'] for s in SOURCES)+remap['remap_P_J']
    dh=remap['H_ref_J']-start['H_ref_J'];dp=remap['cP_ref_J']-start['cP_ref_J']
    rh=dh+end['Q_H_J']+end['WJ_particle_J']-sourceH-end['C_H_J']
    rp=dp+end['Q_P_J']-end['PJ_particle_J']-sourceP-end['C_P_J']+end['GaussE_J']+end['GaussB_J']
    scale=abs(dh)+abs(dp)+sum(abs(v) for v in end.values())+abs(sourceH)+abs(sourceP)
    bound=4096*LD(np.finfo(float).eps)*np.maximum(scale,LD('1e-100'))
    assert np.all(abs(rh-end['field_R_H_J'])<=bound) and np.all(abs(rp-end['field_R_P_J'])<=bound)
    # A rank is an open subdomain: the stored stencil identity contains
    # interface exchanges that cancel only after summing all ranks.
    # Preserve these values and test the global identity below.
    return stages,float(np.max(abs(rh))),float(np.max(abs(rp)))


def one(parent,out,case):
    d=parent/case;o=out/case;o.mkdir(parents=True,exist_ok=True)
    paths=sorted((d/'lab-audit').glob('whole-cowan-rank-*.csv'));assert len(paths)==16
    for rank,p in enumerate(paths):assert p.name==f'whole-cowan-rank-{rank:05d}.csv' and p.with_suffix('.complete').read_text().split()==['1','26923']
    stages={};errH=errP=0.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for result,h,p in pool.map(load_rank,paths):
            errH=max(errH,h);errP=max(errP,p)
            for stage,values in result.items():
                if stage not in stages:stages[stage]={k:v.copy() for k,v in values.items()}
                else:
                    for key,value in values.items():stages[stage][key]+=value
    start,end,remap=[stages[x] for x in ['begin','fixed_end','remap']]
    with h5py.File(d/'ledger/energy-ledger.h5') as f:rows=f['energy_ledger/records'][:]
    G=LD(rows['boost_gamma'][0]);b=LD(rows['boost_beta'][0])
    dh=remap['H_ref_J']-start['H_ref_J'];dp=remap['cP_ref_J']-start['cP_ref_J']
    initial=G*(start['H_ref_J'][0]+b*start['cP_ref_J'][0])
    delta_field=G*(dh+b*dp)
    outward=G*(end['Q_H_J']+b*end['Q_P_J'])
    cowan=G*(end['C_H_J']+b*end['C_P_J'])
    gauss=G*b*(end['GaussE_J']+end['GaussB_J'])
    current=G*(end['WJ_particle_J']-b*end['PJ_particle_J'])
    sources={s:G*(end[s+'_U_J']+b*end[s+'_P_J']) for s in SOURCES}
    moving=G*(remap['remap_H_J']+b*remap['remap_P_J'])
    field_identity=delta_field+outward+current-sum(sources.values())-moving-cowan+gauss
    native_identity=G*(end['field_R_H_J']+b*end['field_R_P_J'])
    global_scale=abs(dh)+abs(dp)+sum(abs(v) for v in end.values())
    global_bound=8192*LD(np.finfo(float).eps)*np.maximum(global_scale,LD('1e-100'))
    global_rh=dh+end['Q_H_J']+end['WJ_particle_J']-sum(end[s+'_U_J'] for s in SOURCES)-remap['remap_H_J']-end['C_H_J']
    global_rp=dp+end['Q_P_J']-end['PJ_particle_J']-sum(end[s+'_P_J'] for s in SOURCES)-remap['remap_P_J']-end['C_P_J']+end['GaussE_J']+end['GaussB_J']
    assert np.all(abs(global_rh)<=global_bound) and np.all(abs(global_rp)<=global_bound)
    cumulative=lambda v:np.r_[LD(0),np.cumsum(v,dtype=LD)]
    particle=rows['delta_particle_kinetic_lab_J'].astype(LD)
    values=dict(particle_kinetic_change_lab_J=particle,field_Lorentz_charge_change_J=cumulative(delta_field),
        outward_energy_momentum_exchange_J=cumulative(outward),Cowan_energy_momentum_channel_J=cumulative(cowan),
        Gauss_energy_momentum_channel_J=cumulative(gauss),deposited_current_Lorentz_work_J=cumulative(current),
        remap_Lorentz_exchange_J=cumulative(moving),field_Lorentz_identity_residual_J=cumulative(field_identity),
        stored_field_identity_residual_J=cumulative(native_identity))
    for s,v in sources.items():values[s+'_Lorentz_exchange_J']=cumulative(v)
    gap=particle-cumulative(current)
    values['particle_vs_deposited_current_Lorentz_gap_J']=gap
    # This full-window quantity includes frozen/reference sources and is not
    # the physical fixed-lab worldtube closure residual.
    values['full_window_Lorentz_accounting_difference_J']=gap+cumulative(field_identity)
    prediction=(-values['field_Lorentz_charge_change_J']-values['outward_energy_momentum_exchange_J']
        +sum(cumulative(v) for v in sources.values())+values['remap_Lorentz_exchange_J']
        +values['Cowan_energy_momentum_channel_J']-values['Gauss_energy_momentum_channel_J'])
    values['field_boundary_source_prediction_lab_J']=prediction
    assert np.max(abs(particle-prediction-values['full_window_Lorentz_accounting_difference_J']))<LD('1e-14')
    max_continuity=max(float(np.max(abs(start['H_ref_J'][1:]-remap['H_ref_J'][:-1]))),
                       float(np.max(abs(start['cP_ref_J'][1:]-remap['cP_ref_J'][:-1]))))
    np.savez_compressed(o/'whole_cache.npz',**values,**{s+'_'+k:v for s,x in stages.items() for k,v in x.items()})
    with (o/'whole_lab_slice.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['step']+list(values))
        for i in range(len(rows)):writer.writerow([i]+[str(v[i]) for v in values.values()])
    result=dict(complete=True,source='All16rank complete whole-cowan ledgers; every stage and step checked.',
        scope='Lorentz energy-momentum charge on boost-time slices, not a same-laboratory-time volume integral or the lab worldtube.',
        maximum_per_rank_energy_identity_J=errH,maximum_per_rank_momentum_identity_J=errP,
        global_field_identities_passed=True,rank_local_identities_are_open_subdomains=True,
        max_interstep_H_or_cP_continuity_J=max_continuity,
        max_cumulative_Lorentz_field_identity_J=float(np.max(abs(values['field_Lorentz_identity_residual_J']))),
        max_identity_reconstruction_vs_stored_J=float(np.max(abs(cumulative(field_identity)-cumulative(native_identity)))),
        final_components_J={k:float(v[-1]) for k,v in values.items()},physical_acceptance=False,
        missing_finer_particle_force_split='Actual particle substep force-geometry trace was not enabled in these production jobs; no external/grid electric/magnetic impulse split is inferred from residuals.')
    (o/'WHOLE_RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return values,result


def main(parent,out):
    parent,out=Path(parent),Path(out);out.mkdir(exist_ok=True);answers=[]
    fig,axes=plt.subplots(2,2,figsize=(14,9),sharex=True,layout='constrained')
    for col,case in enumerate(['energy-binomial3','energy-momentum']):
        values,result=one(parent,out,case);answers.append(result);x=np.arange(len(next(iter(values.values()))))
        for key in ['field_Lorentz_charge_change_J','outward_energy_momentum_exchange_J','Cowan_energy_momentum_channel_J',
                    'Gauss_energy_momentum_channel_J','remap_Lorentz_exchange_J','deposited_current_Lorentz_work_J']:
            axes[0,col].plot(x,values[key]*1e6,label=key.removesuffix('_J'),lw=1)
        for s in SOURCES:axes[1,col].plot(x,values[s+'_Lorentz_exchange_J']*1e6,label=s,lw=1)
        axes[0,col].set_title(case)
        for row in range(2):axes[row,col].set(xlabel='Maxwell step',ylabel='Signed cumulative channel (µJ)');axes[row,col].legend(fontsize=6);axes[row,col].grid(alpha=.2)
    fig.suptitle('Lorentz combination of complete numerical field energy and momentum ledgers\nBoost-time slices; separate geometry from fixed-lab physical worldtube')
    fig.savefig(out/'lab_slice_numerical.png',dpi=170);fig.savefig(out/'lab_slice_numerical.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8),sharex=True,layout='constrained')
    for col,case in enumerate(['energy-binomial3','energy-momentum']):
        with np.load(out/case/'whole_cache.npz') as cache:
            x=np.arange(len(cache['particle_kinetic_change_lab_J']))
            for key,label,style in [('particle_kinetic_change_lab_J','Particle kinetic change (same cohort)','-'),
                ('field_boundary_source_prediction_lab_J','Field/boundary/source prediction','--')]:
                axes[0,col].plot(x,cache[key]*1e6,style,label=label)
            for key,label in [('full_window_Lorentz_accounting_difference_J','Complete numerical difference'),
                ('particle_vs_deposited_current_Lorentz_gap_J','Particle vs deposited-current gap'),
                ('field_Lorentz_identity_residual_J','Field-operator arithmetic difference')]:
                axes[1,col].plot(x,cache[key]*1e6,label=label)
        axes[0,col].set_title(case)
        for row in range(2):
            axes[row,col].set(xlabel='Maxwell step',ylabel='Signed energy (µJ)');axes[row,col].legend(fontsize=8);axes[row,col].grid(alpha=.2)
            axes[row,col].axhline(0,color='grey',lw=.6)
    fig.suptitle('Laboratory particle-energy change and full numerical energy/momentum budget\nCharges evaluated on boost-time slices; not a fixed-lab worldtube closure test')
    fig.savefig(out/'lab_slice_balance.png',dpi=175);fig.savefig(out/'lab_slice_balance.pdf');plt.close(fig)
    (out/'WHOLE_SUMMARY.json').write_text(json.dumps(dict(complete=True,cases=answers),indent=2,allow_nan=False)+'\n')
    print(json.dumps(answers,indent=2),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
