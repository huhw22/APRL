#!/usr/bin/env python3
"""Small analytic checks that an exit subset cannot bias its entry baseline."""
import argparse,json,tempfile
from pathlib import Path
import h5py,numpy as np
from scipy.constants import c,e,m_e
import analyze


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    results={}
    with tempfile.TemporaryDirectory(prefix='cohort-analysis-only-') as tmp:
        root=Path(tmp);(root/'detectors').mkdir()
        dtype=np.dtype([('particle_id','u8'),('source_id','u8'),('time_s','f8'),('position_m','f8',(3,)),
          ('proper_velocity','f8',(3,)),('charge_C','f8'),('mass_kg','f8'),('weight','f8'),('crossing_direction','i4')])
        a=np.zeros(3,dtype=dtype);a['particle_id']=[1,2,3];a['source_id']=[90001,7,330000]
        a['time_s']=[-1e-15,0,1e-15];a['position_m'][:,2]=-.05;a['position_m'][:,0]=[1e-5,2e-5,3e-5]
        gamma=np.array([1000.,1100.,1200.]);a['proper_velocity'][:,2]=np.sqrt(gamma*gamma-1)
        a['charge_C']=-e;a['mass_kg']=m_e;a['weight']=1;a['crossing_direction']=1
        itype=np.dtype([('plane_position_m','f8',(3,)),('proper_velocity','f8',(3,)),('source_id','u8'),('macro_weight','f8')])
        inputs=np.zeros(3,dtype=itype);inputs['plane_position_m']=a['position_m'];inputs['proper_velocity']=a['proper_velocity'];inputs['source_id']=a['source_id'];inputs['macro_weight']=1
        with h5py.File(root/'input.h5','w') as h:h.create_group('particles')['records']=inputs
        plan={'records':3,'input_file':str(root/'input.h5'),'geometry':{'gamma_boost':929.1305492858543},'field_plane_z_m':4.5}
        b=a[[2,1]].copy();b['position_m'][:,2]=4.5;b['time_s']+=(4.55/c)
        target_gamma=np.array([1199.997,1099.999]);b['proper_velocity'][:,2]=np.sqrt(target_gamma*target_gamma-1)
        def write(tag,values,z):
            with h5py.File(root/'detectors'/f'{tag}-particles.h5','w') as h:
                g=h.create_group('particle_plane');g['records']=values;g['complete']=1;g['committed_records']=len(values);g.attrs['plane_z_m']=z
        write('entry',a,-.05);write('exit',b,4.5)
        out=root/'valid';out.mkdir();summaries=analyze.scan_particles(root,out,plan,allow_exit_subset=True)
        detail=json.loads((out/'particle_statistics.json').read_text())
        assert detail['entry_all_input_cohort']['records']==3 and detail['available_counts']['exit']==2
        assert abs(summaries[0]['mean_gamma']-1150)<1e-10
        assert abs(summaries[1]['mean_loss_gamma']-.002)<1e-10
        assert abs(summaries[1]['std_loss_gamma']-.001)<1e-10
        assert not detail['full_exit_coverage'];results['exact_same_cohort_mean_std']=True
        strict=root/'strict';strict.mkdir()
        try:analyze.scan_particles(root,strict,plan);raise RuntimeError('Strict scan accepted absent ID')
        except AssertionError:results['strict_full_coverage_still_rejects']=True
        for name,changed in [('duplicate',b.copy()),('source_mismatch',b.copy())]:
            if name=='duplicate':changed[1]=changed[0]
            else:changed['source_id'][0]+=1
            write('exit',changed,4.5);target=root/name;target.mkdir()
            try:analyze.scan_particles(root,target,plan,allow_exit_subset=True);raise RuntimeError('Subset accepted '+name)
            except AssertionError:results[name+'_rejected']=True
        full=a.copy();full['position_m'][:,2]=4.5;full['time_s']+=4.55/c
        full_gamma=gamma-np.array([.002,.001,.003]);full['proper_velocity'][:,2]=np.sqrt(full_gamma*full_gamma-1)
        write('exit',full,4.5);target=root/'full';target.mkdir()
        expected=analyze.scan_particles(root,target,plan)
        assert abs(expected[0]['mean_gamma']-1100)<1e-10
        assert abs(expected[1]['mean_loss_gamma']-.002)<1e-10
        assert abs(expected[1]['std_loss_gamma']-np.sqrt(2/3)*.001)<1e-10
        results['default_full_cohort_still_accepts']=True
    args.output.write_text(json.dumps({'checks_pass':True,'checks':results,'new_pic_simulations':0,'scope':'Analysis cohort arithmetic only'},indent=2)+'\n')
    print(json.dumps(results))


if __name__=='__main__':main()
