#!/usr/bin/env python3
"""Bound local Boris/time-stagger effects under independent continuum fields."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'longitudinal_momentum_origin'))
from analyze import Sphere,GR,BR,G_BOOST,DT,c,e,boris,LD
from force_reference import DenseSphereForce

def main(root):
    plan=json.loads((root/'plan.json').read_text());original=Path(plan['original']);cfg=json.loads((original/'case.json').read_text())
    result=json.loads((original/'analysis/RESULT.json').read_text());model=Sphere('converging');force=DenseSphereForce(model)
    shell=np.repeat(np.arange(model.nr),model.nmu);mu=np.tile(model.mu,model.nr);direction=np.c_[np.sqrt(1-mu*mu),np.zeros(len(mu)),mu]
    w=model.weights.reshape(model.nr,model.nmu,model.nphi).sum(axis=2).ravel()
    origin=LD(cfg['geometry']['initial_frame_lab_time_s']);offset=LD(result['entry_time_offset_s']);gb=LD(cfg['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    def state(step):
        rhs=(model.entry_center_time+(LD(step)*DT+G_BOOST*(origin+offset))/GR)/model.t0;tau=np.full(len(shell),float(rhs))
        for _ in range(8):
            r,p=model.shell_state(shell,tau);tau=np.asarray(rhs-BR*r*model.sigma*mu/(c*model.t0),dtype=float)
        r,p=model.shell_state(shell,tau);rest=(r*model.sigma)[:,None]*direction;u=(np.sqrt(model.eta)*p)[:,None]*direction
        u[:,2]=GR*(u[:,2]+BR*np.sqrt(1+np.sum(u*u,axis=1)))
        pos=rest.astype(LD);pos[:,2]=(pos[:,2]+gb*bb*c*(origin+offset))/GR+BR*c*LD(step)*DT
        return pos,u
    def fields(pos,time):
        tau=(model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/model.t0
        rest=pos.copy();rest[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
        E=force.field(rest,tau).astype(LD);E[:,:2]*=GR
        return E,np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
    rows=[]
    for step in (140,200,300,400,480,800,1000):
        pos,u=state(step);_,end=state(step+1);E,B=fields(pos,LD(step)*DT);_,Bh=fields(pos,(LD(step)+.5)*DT)
        same=boris(u,E,B,np.full(len(u),DT));stagger=boris(u,E,Bh,np.full(len(u),DT))
        rows.append(dict(step=step,mean_exact_duz=float(np.sum(w*(end[:,2]-u[:,2]))),
            local_error_same_time_B=float(np.sum(w*(same[:,2]-end[:,2]))),local_error_half_time_B=float(np.sum(w*(stagger[:,2]-end[:,2])))))
    (root/'IDEAL_CHECK.json').write_text(json.dumps(dict(scope='local map only, not a long-time PIC convergence proof',rows=rows),indent=2)+'\n')
    print(json.dumps(rows,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root.resolve())
