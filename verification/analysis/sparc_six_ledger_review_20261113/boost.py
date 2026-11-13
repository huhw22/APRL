"""Review both completed six-segment boost ledgers with validated readers."""
import concurrent.futures
import json
import os
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_discrete import load,calculate
from relative_ledger import construct,write_tables

LD=np.longdouble


def one(args):
    parent,out,case=args;d=Path(parent)/case;o=Path(out)/case;o.mkdir(parents=True)
    done=json.loads((d/'PIC_COMPLETE.json').read_text());assert done['complete'] and done['run']['exit_code']==0
    rows,sums,ranks,steps=load(d/'ledger',d/'ledger/energy-ledger.h5')
    assert steps==26923 and len(rows)==26924 and ranks==16
    for rank in range(ranks):
        for prefix in ['lab-field','paired-momentum','whole-cowan']:
            assert (d/'lab-audit'/f'{prefix}-rank-{rank:05d}.complete').read_text().split()==['1',str(steps)]
    assert (d/'worldtube/COMPLETE').read_text().split()==['1',str(steps)]
    index,q,summary=calculate(rows,sums,steps)
    absolute,pct,normalization=construct(rows,sums,index,q);write_tables(o,index,absolute,pct,normalization)
    g=LD(rows['boost_gamma'][0]);b=LD(rows['boost_beta'][0])
    transform=rows['delta_particle_kinetic_lab_J'].astype(LD)-g*(rows['delta_particle_kinetic_box_J'].astype(LD)+b*rows['delta_particle_cPz_box_J'].astype(LD))
    np.savez_compressed(o/'boost_cache.npz',records=rows,**{'sum_'+k:v for k,v in sums.items()},**q)
    report=dict(case=case,complete=True,steps=steps,mpi_ranks=ranks,normalization=normalization,
        final_boost_particle_change_J=float(rows['delta_particle_kinetic_box_J'][-1]),
        final_lab_particle_change_full_cohort_J=float(rows['delta_particle_kinetic_lab_J'][-1]),
        max_particle_Lorentz_identity_error_J=float(np.max(abs(transform))),
        physical_acceptance=False,**summary)
    (o/'BOOST_RESULT.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


def main(parent,out):
    out=Path(out);out.mkdir(exist_ok=True)
    cases=['energy-binomial3','energy-momentum']
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        answers=list(pool.map(one,[(str(parent),str(out),c) for c in cases]))
    fig,axes=plt.subplots(3,2,figsize=(14,11),sharex=True,layout='constrained')
    titles=['Energy matching + three-point kernel','Energy matching + momentum matching']
    for col,(case,title,answer) in enumerate(zip(cases,titles,answers)):
        a=np.genfromtxt(out/case/'energy_budget.csv',delimiter=',',names=True,dtype=LD);x=a['step']
        for n,label,style in [('particle_change_pct','Particle kinetic change','-'),
             ('field_boundary_source_prediction_pct','Independent accounting prediction','--')]:axes[0,col].plot(x,a[n],style,label=label)
        for n,label in [('native_field_change_pct','Native field change'),('outward_curl_exchange_pct','Outward curl exchange'),
             ('prescribed_and_constraint_source_work_pct','External + constrained-source work'),
             ('inlet_energy_added_pct','Inlet overwrite'),('remap_energy_added_pct','Moving-window remap'),
             ('time_pairing_channel_pct','Time/current pairing'),('Cowan_channel_pct','Cowan operator'),
             ('physical_CPML_channel_pct','CPML')]:axes[1,col].plot(x,a[n],label=label,lw=1.)
        for n,label,style in [('total_residual_pct','Total discrete difference','-'),
             ('particle_grid_exchange_difference_pct','Particle/current exchange gap','--'),
             ('field_operator_residual_pct','Field operator arithmetic difference',':')]:axes[2,col].plot(x,a[n],style,label=label)
        axes[0,col].set_title(title)
        for row in range(3):
            axes[row,col].axhline(0,color='grey',lw=.6);axes[row,col].grid(alpha=.2)
            axes[row,col].set_ylabel('% of fixed particle-change scale');axes[row,col].legend(fontsize=7)
        axes[2,col].set_xlabel('Maxwell step')
    fig.suptitle('Boost-frame complete ledger | full interval and all recorded channels\nEach column uses D=max|particle kinetic change|; numerical channels are not radiation')
    fig.savefig(out/'boost_energy.png',dpi=170);fig.savefig(out/'boost_energy.pdf');plt.close(fig)
    (out/'BOOST_SUMMARY.json').write_text(json.dumps(dict(complete=True,cases=answers),indent=2,allow_nan=False)+'\n')
    print(json.dumps(answers,indent=2),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
