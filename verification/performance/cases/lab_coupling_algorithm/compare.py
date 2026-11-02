"""Compare closed baseline/candidate runs using the unchanged budget formula."""
import argparse,json
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def particle(path):
    with h5py.File(path) as f:
        g=f['particle_plane'];assert g['complete'][()]==1;a=g['records'][:];assert len(a)==g['committed_records'][()];a.sort(order='source_id');return a
def main():
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);a=p.parse_args();w=a.workspace;root=w/'results/lab_coupling_algorithm_20261102';base=w/'results/lab_energy_balance_20261102/coupled_sbp_r2';new=root/'r3';assert (new/'EXIT_CODE').read_text().strip()=='0';out={};fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    A=json.loads((base/'RESULT.json').read_text())['cases'];B=json.loads((new/'RESULT.json').read_text())['cases'];norm=json.loads((base/'REVIEW.json').read_text())
    for ax,name in zip(axes,['undulator-on','drift-on']):
        x=np.genfromtxt(base/name/'whole_steps_J.csv',names=True,delimiter=',',dtype=np.longdouble);y=np.genfromtxt(new/name/'whole_steps_J.csv',names=True,delimiter=',',dtype=np.longdouble);assert np.array_equal(x['step'],y['step']) and np.array_equal(x['clean'],y['clean']);mask=x['clean']==1;D=norm[name]['normalization']['denominator_J'];xx=np.r_[x['step'][mask][0]-1,x['step'][mask]]
        for value,label in [(x['total_numerical_discrepancy'],'Original total discrepancy'),(y['total_numerical_discrepancy'],'Candidate total discrepancy'),(x['magnetic_impulse_gap'],'Original magnetic pairing gap')]:ax.plot(xx,np.r_[0,100*np.cumsum(value[mask])/D],label=label)
        ax.set(title=name,xlabel='Field step (all particles active)',ylabel='% of fixed BASELINE particle-change scale');ax.axhline(0,color='black',lw=.5);ax.grid(alpha=.2)
        changes={}
        for plane in ['release','exit']:
            old=particle(base/name/'detectors'/(plane+'.h5'));candidate=particle(new/name/'detectors'/(plane+'.h5'));assert np.array_equal(old['source_id'],candidate['source_id']) and np.array_equal(old['charge_C'],candidate['charge_C']) and np.array_equal(old['mass_kg'],candidate['mass_kg']);weight=abs(old['charge_C']).astype(np.longdouble);weight/=sum(weight)
            def stats(a):
                z=a['position_m'].astype(np.longdouble);mean=np.sum(weight[:,None]*z,axis=0);rms=np.sqrt(np.sum(weight[:,None]*(z-mean)**2,axis=0));return dict(mean_position_m=list(map(float,mean)),rms_position_m=list(map(float,rms)))
            changes[plane]=dict(records=len(old),original=stats(old),candidate=stats(candidate),max_position_difference_m=float(np.max(abs(candidate['position_m']-old['position_m']))),max_proper_velocity_difference=float(np.max(abs(candidate['proper_velocity']-old['proper_velocity']))),all_records_bitwise_equal=bool(np.array_equal(old,candidate)))
        av=A[name]['clean_interval']['eV_per_electron'];bv=B[name]['clean_interval']['eV_per_electron'];oldplane=json.loads((base/name/'PAIRED_RESULT.json').read_text());newplane=json.loads((new/name/'PAIRED_RESULT.json').read_text());out[name]=dict(unchanged_clean_interval=True,original=av,candidate=bv,reference_total_discrepancy_change_percent=float(100*(abs(bv['total_numerical_discrepancy'])/abs(av['total_numerical_discrepancy'])-1)),fixed_plane_gain_eV={'original':oldplane['fixed_plane_particle_gain_eV'],'candidate':newplane['fixed_plane_particle_gain_eV']},particle_planes=changes,normalization=norm[name]['normalization'])
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside upper center',ncol=3,fontsize=8);fig.savefig(root/'magnetic_pairing_comparison.png',dpi=180);fig.savefig(root/'magnetic_pairing_comparison.pdf');plt.close(fig);(root/'COMPARISON.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({n:{k:v for k,v in q.items() if k not in ['original','candidate','particle_planes']} for n,q in out.items()},indent=2))
if __name__=='__main__':main()
