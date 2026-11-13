"""Build external-field option, validate it, generate/import one random beam."""
import math
import os
import re
import sys
from pathlib import Path
import numpy as np
from common import W,read,save,sha,need,frozen,srun,launch
from binding import binding


def reference_points(root,plan):
    rows=[];ku=2*np.pi/.028;length=77*.028
    for j,phase in enumerate(plan['lattice']['segment_phases_rad']):
        start=j*(length+.4)
        zs=np.unique(np.r_[np.linspace(0,length,401),np.linspace(0,.028,101),
                           length-np.linspace(0,.028,101)])
        for local in zs:
            d=min(local,length-local);h=1.;hp=0.
            if d<.028:
                h=.5*(1-np.cos(np.pi*d/.028));hp=.5*np.pi/.028*np.sin(np.pi*d/.028)*(1 if local<.028 else -1)
            z=start+local;ph=ku*local+phase
            for y in [-.0008,-.000075,0.,.000075,.0008]:
                rows.append([z,.00011,y,0.,np.cosh(ku*y)*(h*np.sin(ph)-hp/ku*np.cos(ph)),np.sinh(ku*y)*h*np.cos(ph)])
        if j<10:
            for local in [.001,.05,.105-1e-7,.105+1e-7,.13,.158-1e-7,.158+1e-7,.2,.399]:
                z=start+length+local
                for x,y in [(.0001,.0002),(-.0001,.0002),(.0004,-.0003)]:
                    g=-9/.78796 if .105<=local<=.158 else 0.
                    rows.append([z,x,y,g*y,g*x,0.])
    path=root/'preparation/field-reference.txt';np.savetxt(path,np.array(rows),fmt='%.17e')
    return path,len(rows)


def main(root):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm');root=Path(root);plan=read(root/'PLAN.json')
    frozen(root/'FROZEN_PREPARATION.json')
    launch(['cmake','-S',str(root/'source'),'-B',plan['build'],'-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_Fortran_COMPILER=/usr/bin/mpif90','-DCMAKE_Fortran_FLAGS=-ffree-line-length-none',
        '-DFFTW3_ROOT='+str(W/'build/puffin-deps-20261113'),'-DHDF5_PREFER_PARALLEL=ON'],
        root,root/'preparation/configure.log',300)
    launch(['cmake','--build',plan['build'],'-j','24'],root,root/'preparation/build.log',900)
    # The field reference is independent Python evaluation of the original
    # frozen APRL magnet list, in normalized B/B0 units.
    route=root/'preparation/probe.latt';route.write_text("UN 'aprl-sparc11' 990 1.0 0.0 80 0.0 1.0 0.0 0.0\n")
    points,n=reference_points(root,plan)
    launch(srun(1,1)+[plan['field_probe'],str(route),str(points)],root,root/'preparation/field-probe.log',120)
    log=(root/'preparation/field-probe.log').read_text();need('FIELD_ROUTE_PASS' in log,'Real field route probe failed')
    save(root/'FIELD_VALIDATED.json',dict(passed=True,points=n,actual_route_validated=True,
        type=plan['external_type'],only_external_field_extended=True,log=log))
    launch(srun(1,2)+[plan['binary'],'inputs/3D/clara_test.in'],root/'smoke-matched',root/'preparation/smoke.log',240)
    from legacy_validate import smoke,native
    smoke(root)
    launch(srun(1,128)+[str(W/'build/cluster-tools/binding_probe'),'128','1'],root,root/'preparation/binding.log',120)
    save(root/'preparation/BINDING.json',binding(root/'preparation/binding.log',1,128,1))
    from make_inputs import generate
    generate(root)
    launch(srun(1,128)+[plan['binary'],'main.in'],root/'import-check',root/'preparation/import.log',600,True)
    need('Puffin simulation completed' in (root/'preparation/import.log').read_text(),'Native preflight did not close')
    native(root)
    import h5py
    from common import dump_at
    with h5py.File(root/'particles.h5') as f: initial=float(f['electrons'][:,3].mean())
    first=dump_at(root/'import-check','electrons',1)
    with h5py.File(first) as f:
        observed=(float(f['electrons'][:,3].mean())-initial)*plan['K']
    expected=-plan['K']*.5*(1-math.cos(math.pi/80))*math.cos(2*math.pi/80)
    need(abs(observed-expected)<1e-5,'Actual first-step smooth magnetic trajectory failed')
    save(root/'FIRST_STEP_VALIDATED.json',dict(passed=True,physical_step_m=.028/80,
        observed_mean_delta_ux=observed,on_axis_expected_delta_ux=expected,
        field_feedback_disabled=True,all_original_events_retained=True))
    frozen(root/'FROZEN_PREPARATION.json')
    files=[root/'PLAN.json',root/'particles.h5',root/'inputs/independent-w64.h5',root/'INPUT_VALIDATED.json',
           root/'NATIVE_VALIDATED.json',root/'FIELD_VALIDATED.json',root/'SMOKE_VALIDATED.json',root/'FIRST_STEP_VALIDATED.json',
           Path(plan['binary']),Path(plan['field_probe']),
           *[p for p in (root/'scripts').iterdir() if p.is_file()],
           *[p for p in (root/'main').iterdir() if p.is_file()]]
    save(root/'FROZEN_PHYSICS.json',{str(p):sha(p) for p in files})
    save(root/'PREPARATION_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],
        particle_count=plan['sample_records'],source_seed=plan['source_seed'],APRL_input_saved=True,
        actual_native_import_bitwise=True,field_route_and_stock_path_validated=True,physical_acceptance=False))


if __name__=='__main__': main(sys.argv[1])
