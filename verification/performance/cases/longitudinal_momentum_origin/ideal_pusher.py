#!/usr/bin/env python3
"""Single-step truncation of the existing Boris map under independent ideal fields.

This local consistency check is not a proof of long-time PIC convergence.
It uses no field grid, deposition, CPML, moving inlet, or fitted force.
"""
import argparse,json
from pathlib import Path
import numpy as np
from analyze import Sphere,GR,BR,G_BOOST,DT,c,e,epsilon_0,boris,LD

def main(original,out):
    config=json.loads((original/'case.json').read_text());result=json.loads((original/'analysis/RESULT.json').read_text())
    model=Sphere('uniform');shell=np.repeat(np.arange(model.nr),model.nmu);mu=np.tile(model.mu,model.nr)
    direction=np.c_[np.sqrt(1-mu*mu),np.zeros(len(mu)),mu]
    w=model.weights.reshape(model.nr,model.nmu,model.nphi).sum(axis=2).ravel()
    origin=LD(config['geometry']['initial_frame_lab_time_s']);offset=LD(result['entry_time_offset_s'])
    def state(step):
        rhs=(model.entry_center_time+(LD(step)*DT+G_BOOST*(origin+offset))/GR)/model.t0
        tau=np.full(len(shell),float(rhs))
        for _ in range(8):
            r,p=model.shell_state(shell,tau);tau=np.asarray(rhs-BR*r*model.sigma*mu/(c*model.t0),dtype=float)
        r,p=model.shell_state(shell,tau);u=(np.sqrt(model.eta)*p)[:,None]*direction
        g=np.sqrt(1+np.sum(u*u,axis=1));u[:,2]=GR*(u[:,2]+BR*g)
        return tau,r,u
    rows=[]
    for step in (140,200,300,400,480,800):
        tau,r,u=state(step);_,_,end=state(step+1)
        estar=(-model.electrons*e/(4*np.pi*epsilon_0)*model.q[shell]/(r*model.sigma)**2)[:,None]*direction
        E=estar.astype(LD);E[:,:2]*=GR;B=np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
        pushed=boris(u,E,B,np.full(len(u),DT))
        rows.append(dict(step=step,mean_local_Boris_uz_error=float(np.sum(w*(pushed[:,2]-end[:,2]))),
            mean_exact_uz_increment=float(np.sum(w*(end[:,2]-u[:,2])))))
    out.write_text(json.dumps(dict(scope='Local ideal-field truncation only; no PIC accuracy acceptance',rows=rows),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('original',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();main(a.original,a.output)
