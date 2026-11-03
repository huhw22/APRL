"""Same-ID arrival-time phase differences; no fitted clock or energy proxy."""
import argparse,itertools,json
from pathlib import Path
import numpy as np
from scipy.constants import c
from analyze_fel import read_particles
from prepare import save

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);a=ap.parse_args()
    root=a.root.resolve();plan=json.loads((root/'campaign.json').read_text())
    n=plan['parameters']['macroparticles'];omega=2*np.pi*c/plan['parameters']['resonant_wavelength_m']
    available=[q for q in plan['cases'] if (Path(q['directory'])/'PIC_COMPLETE.json').exists()]
    rows=[]
    for name,z in [('entry',0.),*plan['geometry']['field_planes']]:
        filename='entry.h5' if name=='entry' else name+'-particles.h5'
        records={q['name']:read_particles(Path(q['directory'])/'detectors'/filename,n) for q in available}
        for left,right in itertools.combinations(records,2):
            aa,bb=records[left],records[right]
            weights=abs(aa['charge_C']);weights/=weights.sum()
            phase=-omega*(aa['time_s'].astype(np.longdouble)-bb['time_s'].astype(np.longdouble))
            mean=float(np.sum(weights*phase));std=float(np.sqrt(np.sum(weights*(phase-mean)**2)))
            rows.append(dict(left=left,right=right,plane=name,z_m=z,
               mean_phase_rad=mean,std_phase_rad=std,
               rms_phase_rad=float(np.sqrt(np.sum(weights*phase*phase))),
               p99_abs_phase_rad=float(np.quantile(abs(phase),.99)),max_abs_phase_rad=float(max(abs(phase))),
               cohort=n,definition='-omega_r*(t_left-t_right), original common input clock'))
    save(root/'PARTICLE_PHASE_COMPARISON.json',rows)
    print(json.dumps([q for q in rows if q['plane']=='far'],indent=2))

if __name__=='__main__':main()
