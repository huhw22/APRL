#!/usr/bin/env python3
"""Compare actual held inlet charge with the rigidly translated reference.

Uses invariant input-plane events; does not read running PIC output. Scalar
rho marginal validates the reconstructed initial held coordinates first.
"""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.constants import c,e
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'space_charge_physics'))
from model import Sphere,DZ,G_BOOST
LD=np.longdouble
def main(root):
    plan=json.loads((root/'plan.json').read_text());original=Path(plan['original']);cfg=json.loads((original/'case.json').read_text());old=json.loads((original/'analysis/RESULT.json').read_text())
    ex=np.load(original/'reference_input/expected_events.npz');model=Sphere('converging');w=model.weights.astype(LD);G=LD(G_BOOST);beta=np.sqrt(1-G**-2)
    held=np.tanh(np.arcsinh(ex['entry__lab_u'][:,2].astype(LD))-np.arccosh(G))
    tbox=G*(ex['entry__lab_time_s'].astype(LD)-LD(old['entry_time_offset_s'])-LD(cfg['geometry']['initial_frame_lab_time_s']))
    z0=LD(cfg['geometry']['box_translation_m'])-LD(c)*(beta+held)*tbox
    coord=(z0-LD(cfg['geometry']['origin_reference_m']))/LD(DZ);q=-LD(model.electrons)*LD(e)*w
    nz=cfg['geometry']['reference_cells'];lo=np.floor(coord).astype(int);f=coord-lo;prof=np.zeros(nz+1,dtype=LD)
    np.add.at(prof,lo,q*(1-f));np.add.at(prof,lo+1,q*f)
    data=np.memmap(root/'base/reference/rho.bin',dtype='float64',mode='r',shape=(3327,139,139))
    marginal=np.asarray(data[1400:1750]).sum(axis=(1,2),dtype=LD)*LD('5e-6')**2*LD(DZ)
    rel=float(np.linalg.norm(np.asarray(marginal-prof[1401:1751],float))/np.linalg.norm(np.asarray(marginal,float)));assert rel<1e-6
    bt=LD(np.mean(held));bw=np.sum(w*held);rows=[]
    def rigid(k,b):
        k=float(k);i=int(np.floor(k));f=k-i;return (1-f)*prof[i]+f*prof[i+1]
    for guard in (64,128,256):
        initial=cfg['geometry']['initial_offset']-(guard-64)
        for step in range(150):
            k=initial-step+1;actual=np.sum(q*np.maximum(0,1-abs(coord+held*step-k)))
            exact_common=np.sum(q*np.maximum(0,1-abs(coord+bt*step-k)))
            rows.append(dict(guard=guard,step=step,actual_held_charge_C=float(actual),common_velocity_exact_CIC_charge_C=float(exact_common),reference_rigid_charge_C=float(rigid(k-bt*step,bt)),
                weighted_velocity_reference_charge_C=float(rigid(k-bw*step,bw))))
    with (root/'held_inlet_charge.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    result=dict(initial_marginal_relative_L2=rel,held_beta_charge_weighted=float(bw),reference_translation_beta_unweighted=float(bt),
        held_beta_min=float(min(held)),held_beta_max=float(max(held)),scope='Charge on new inlet charge node k=1; actual held trajectories vs scalar reference rigid advection. This is a density mismatch measurement, not by itself a Gauss-defect or force attribution.')
    for guard in (64,128,256):
        a=[r for r in rows if r['guard']==guard];x=np.array([r['actual_held_charge_C'] for r in a]);y=np.array([r['reference_rigid_charge_C'] for r in a]);z=np.array([r['weighted_velocity_reference_charge_C'] for r in a]);equal=np.array([r['common_velocity_exact_CIC_charge_C'] for r in a])
        scale=max(abs(x))
        result[str(guard)]=dict(actual_peak_charge_C=float(scale),rigid_mismatch_L1_C=float(sum(abs(x-y))),
            interpolation_only_mismatch_L1_C=float(sum(abs(y-equal))),individual_velocity_mismatch_L1_C=float(sum(abs(equal-x))),
            charge_weighted_rigid_mismatch_L1_C=float(sum(abs(x-z))),maximum_density_mismatch_over_actual_peak=float(max(abs(x-y))/scale) if scale else None,
            last_nonzero_actual_node_step=max([r['step'] for r in a if r['actual_held_charge_C']!=0],default=None))
    (root/'HELD_PROFILE.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root.resolve())
