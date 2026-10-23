#!/usr/bin/env python3
"""Stream ballistic-reference coverage using field_reconstruction's projection/grid rules."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

C=299792458.0
PADDINGS=((1,1,1),(2,1,1),(1,2,2),(2,2,2),(1,3,3),(2,3,3),(1,4,4))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--field',type=Path,required=True)
    p.add_argument('--particles',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with h5py.File(a.field,'r') as f:
        g=f['field_plane']; assert int(g['complete'][()])==1
        nt=int(g['committed_samples'][()]); ny,nx=g['fields'].shape[1:]
        times=g['time_s'][:nt]; assert np.all(np.diff(times)>0)
        dt=(times[-1]-times[0])/(nt-1)
        dx,dy=float(g.attrs['dx_m']),float(g.attrs['dy_m'])
        x0,y0=float(g.attrs['x_first_m']),float(g.attrs['y_first_m'])
        zfield=float(g.attrs['plane_z_m'])
    grids=[]
    for pt,py,px in PADDINGS:
        shape=(nt*pt,ny*py,nx*px)
        first=np.array([times[0]-((shape[0]-nt)//2)*dt,
                        y0-((shape[1]-ny)//2)*dy,x0-((shape[2]-nx)//2)*dx])
        spacing=np.array([dt,dy,dx])
        nreal=int(np.prod(shape));ncomplex=shape[0]*shape[1]*(shape[2]//2+1)
        grids.append(dict(padding=[pt,py,px],shape=list(shape),first=first,spacing=spacing,
                          last=first+(np.array(shape)-1)*spacing,outside_charge=np.longdouble(0),
                          outside_count=0,outside_time_count=0,outside_xy_count=0,
                          explicit_allocation_GiB=(4*nreal*8+5*ncomplex*16)/2**30,
                          main_transform_resident_GiB=max(3*nreal*8+3*ncomplex*16,nreal*8+5*ncomplex*16)/2**30))
    sums=np.zeros(5,dtype=np.longdouble); extrema=np.array([[np.inf]*3,[-np.inf]*3])
    with h5py.File(a.particles,'r') as f:
        g=f['ballistic_reference'];assert int(g['complete'][()])==1
        n=int(g['committed_records'][()]);ds=g['records'];assert ds.shape==(n,)
        zref=float(g.attrs['plane_z_m']);seen=np.zeros(n,dtype=bool)
        assert zfield>=zref
        for i in range(0,n,65536):
            r=ds[i:min(n,i+65536)]
            for key in ('time_s','position_m','proper_velocity','charge_C'):
                assert np.isfinite(r[key]).all(),key
            ids=r['particle_id'].astype(np.int64)-1
            assert np.all((ids>=0)&(ids<n))
            assert np.unique(ids).size==len(ids) and not seen[ids].any()
            seen[ids]=True
            assert np.all(r['crossing_direction']==1)
            u=r['proper_velocity']; gamma=np.sqrt(1+np.einsum('ij,ij->i',u,u))
            velocity=C*u/gamma[:,None];assert np.all(velocity[:,2]>0)
            flight=(zfield-zref)/velocity[:,2]
            projected=np.column_stack((r['time_s']+flight,
                                       r['position_m'][:,1]+velocity[:,1]*flight,
                                       r['position_m'][:,0]+velocity[:,0]*flight))
            assert np.isfinite(projected).all()
            w=np.abs(r['charge_C']);assert np.all(w>0)
            sums+=np.array([np.sum(w,dtype=np.longdouble),np.sum(w*gamma,dtype=np.longdouble),
                            np.sum(w*gamma*gamma,dtype=np.longdouble),
                            np.sum(w*np.sum((velocity[:,:2]/C)**2,axis=1),dtype=np.longdouble),
                            np.sum(w/(gamma**2*velocity[:,2]**2),dtype=np.longdouble)])
            extrema[0]=np.minimum(extrema[0],projected.min(axis=0));extrema[1]=np.maximum(extrema[1],projected.max(axis=0))
            for grid in grids:
                # Match interpolationWeights' coordinate test, including cell centres.
                coords=(projected-grid['first'])/grid['spacing']
                bad=(coords<0)|(coords>np.array(grid['shape'])-1)
                outside=bad.any(axis=1)
                grid['outside_charge']+=np.sum(w[outside],dtype=np.longdouble)
                grid['outside_count']+=int(outside.sum())
                grid['outside_time_count']+=int(bad[:,0].sum())
                grid['outside_xy_count']+=int(bad[:,1:].any(axis=1).sum())
            if i%(65536*64)==0: print(f'projected {min(i+65536,n)}/{n}',flush=True)
        assert seen.all()
    mean=sums[1]/sums[0]; rel=float(np.sqrt(max(np.longdouble(0),sums[2]/sums[0]-mean**2))/mean)
    for grid in grids:
        grid['outside_charge_fraction']=float(grid.pop('outside_charge')/sums[0])
        grid['passes_original_outside_limit']=grid['outside_charge_fraction']<=0.01
        for key in ('first','last','spacing'):grid[key]=grid[key].tolist()
    report=dict(records=n,unique_complete_ID_cohort=True,field_z_m=zfield,reference_z_m=zref,
                mean_gamma=float(mean),relative_gamma_spread=rel,
                rms_transverse_beta=float(np.sqrt(sums[3]/sums[0])),
                original_limits=dict(maximum_relative_gamma_spread=.05,maximum_rms_transverse_beta=.05,maximum_outside_charge_fraction=.01),
                projected_time_y_x_min=extrema[0].tolist(),projected_time_y_x_max=extrema[1].tolist(),
                grids=grids,note='Coverage calculation only; larger padding cases have not been reconstructed. Memory excludes FFTW plan temporaries and page cache.')
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
