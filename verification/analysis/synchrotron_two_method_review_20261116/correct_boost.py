"""Retain separate actual source works in the physical boost display."""
import csv,json,shutil,sys
from pathlib import Path
import h5py
import numpy as np
LD=np.longdouble
raw,out=map(lambda p:Path(p).resolve(),sys.argv[1:])
r=json.loads((out/'RESULT.json').read_text());hist=out/'prepublication_boost_display';hist.mkdir()
shutil.copyfile(out/'RESULT.json',hist/'RESULT.json')
checks=[]
for case in r['cases']:
    d=raw/case['directory_relative'];o=out/'cases'/case['tag'];target=o/'boost/physical_components.csv'
    shutil.copyfile(target,hist/(case['tag']+'.csv'))
    with h5py.File(d/'ledger/energy-ledger.h5') as f:a=f['energy_ledger/records'][:]
    v=lambda k:a[k].astype(LD)
    faces=v('outward_field_energy_face_J')
    q=dict(particle_change_J=v('delta_particle_kinetic_box_J'),field_change_J=v('delta_field_energy_J'),
        front_outward_minus_injection_J=faces[:,4]-v('moving_inserted_field_J')-v('inlet_overwrite_work_J'),
        back_outward_plus_cropped_J=faces[:,5]+v('moving_cropped_field_J'),side_outward_J=faces[:,:4].sum(axis=1),
        moving_interface_J=-v('moving_remap_interface_correction_J'),
        external_and_constraint_work_signed_J=-v('prescribed_boris_work_J')-v('boundary_source_work_J')-v('frozen_constraint_work_J')-v('carrier_constraint_work_J'))
    q['physical_component_sum_J']=sum(q.values())
    scale=sum(abs(x) for x in q.values());bound=2048*LD(np.finfo(float).eps)*np.maximum(scale,LD('1e-100'))
    diff=abs(q['physical_component_sum_J']-v('balance_residual_boris_J'));assert np.all(diff<=bound)
    q.update(external_prescribed_work_signed_J=-v('prescribed_boris_work_J'),boundary_source_work_signed_J=-v('boundary_source_work_J'),
        frozen_constraint_work_signed_J=-v('frozen_constraint_work_J'),carrier_constraint_work_signed_J=-v('carrier_constraint_work_J'))
    values=dict(time_box_s=v('time_box_s'),**q)
    with target.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(values);w.writerows(zip(*values.values()))
    case['boost']['final_physical_component_sum_J']=float(q['physical_component_sum_J'][-1])
    (o/'RESULT.json').write_text(json.dumps(case,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    checks.append(dict(case=case['tag'],max_native_monitor_difference_J=float(diff.max())))
(out/'RESULT.json').write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
(out/'BOOST_DISPLAY_REVIEW.json').write_text(json.dumps(dict(passed=True,cases=checks,
    separate_measured_Boris_external_and_boundary_frozen_carrier_works_retained=True,
    original_first_pass_saved=True,raw_physical_outputs_unchanged=True,no_physical_rerun=True),indent=2)+'\n')
print('BOOST_DISPLAY_CORRECTED',checks,flush=True)
