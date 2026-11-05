from pathlib import Path
import json,sys,tempfile,csv
import numpy as np,h5py,yaml
sys.path.insert(0,'/public/home/wfollowm/APRL_Workspace/source/aprl/verification/performance/cases/energy_critical_control')
import mode_review as m
root=Path(tempfile.mkdtemp(prefix='aprl-mode-review-synthetic-'))
gamma=2.;beta=np.sqrt(1-1/gamma**2);dz=1e-5;dt=dz/299792458;nt=6;t0=-2e-14
geom=dict(boost_gamma=gamma,dz_m=dz,dt_box_s=dt,predicted_frame_time_origin_lab_s=t0,mesh_center_z_m=0.,reference_nz=128,active_nz=64,predicted_box_translation_m=.25*dz,predicted_entry_fraction=.25,steps=nt)
m.save(root/'GEOMETRY.json',geom)
rng=np.random.default_rng(4);pat=rng.normal(size=(40,40))
for mode,amp,cos in [('off',1.,-.9),('on',.01,.99)]:
 d=root/mode;d.mkdir();(d/'detectors').mkdir();(d/'reference').mkdir();(d/'ledger').mkdir()
 m.save(d/'RESULT.json',dict(program_integrity_pass=True,stage='full',steps=nt,ranks=1,synthetic_fixture=True))
 m.save(d/'PIC_COMPLETE.json',dict(analysis_complete=True,stage='full',synthetic_fixture=True))
 m.save(d/'run-plan.json',dict(expected_steps=nt,ranks=1))
 (d/'run.yaml').write_text(yaml.safe_dump({'beam':{'reference':{'activation_plane_z':0.}}}))
 (d/'aprl.log').write_text('Moving window: reference_nz=128 field_nz=64 initial_offset=64\n')
 (d/'reference/window-rank-0.csv').write_text('step,origin_box_z_m\n1,'+str(-dz)+'\n')
 with h5py.File(d/'detectors/entry-transport-energy.h5','w') as f:
  g=f.create_group('energy_plane');g.create_dataset('complete',data=1)
  g.attrs.update(time_lab_origin_s=t0,boost_gamma=gamma,boost_beta=beta)
 keys=['step','native_begin_J','native_end_J','remap_J','inlet_overwrite_J','magnetic_cpml_J','time_electric_J','time_magnetic_J']
 p=d/'ledger/discrete-field-rank-00000.csv'
 with p.open('w') as f:
  w=csv.writer(f);w.writerow(keys)
  for n in range(1,nt+1):w.writerow([n,1,1.6,.1,.2,.3,2*(1-cos),0])
 p.with_suffix('.complete').write_text('1 '+str(nt)+'\n')
 for plane,z in [('entry',0.),('und-entry',.01),('p1',2.166)]:
  dtype=np.dtype([('electric_V_per_m','f8',(40,40,3)),('magnetic_T','f8',(40,40,3))])
  data=np.zeros(nt+1,dtype=dtype)
  for n in range(nt+1):
   frac=(.25+z/(gamma*dz)+(1-beta)*n)%1;e=1-2*frac;h=-2*min(frac,1-frac)
   ev=amp*pat*e;hv=amp*pat*h
   data[n]['electric_V_per_m'][...,0]=gamma*(ev+beta*hv)
   data[n]['magnetic_T'][...,1]=gamma*(hv+beta*ev)/299792458
  with h5py.File(d/'detectors'/(plane+'.h5'),'w') as f:
   g=f.create_group('field_plane');g.create_dataset('complete',data=1);g.create_dataset('committed_samples',data=nt+1)
   g.create_dataset('time_s',data=t0+beta*z/299792458+np.arange(nt+1)*dt/gamma);g.create_dataset('fields',data=data);g.attrs['plane_z_m']=z
# Completion guards must reject before selected field data is opened.
marker=root/'on/PIC_COMPLETE.json'; original_marker=marker.read_text(); marker.unlink()
try:m.compare(root)
except RuntimeError as e:assert 'Missing final PIC completion marker' in str(e)
else:raise AssertionError('Missing final PIC marker accepted')
marker.write_text(original_marker)
with h5py.File(root/'off/detectors/entry-transport-energy.h5','r+') as f:f['energy_plane/complete'][()]=0
try:m.compare(root)
except RuntimeError as e:assert 'Unclosed energy plane' in str(e)
else:raise AssertionError('Incomplete energy plane accepted')
with h5py.File(root/'off/detectors/entry-transport-energy.h5','r+') as f:f['energy_plane/complete'][()]=1
out=m.compare(root);r=json.loads((out/'RESULT.json').read_text())
ratios=[x['on_over_off_projection_rms'] for x in r['comparisons'] if x['component']=='x']
assert max(abs(np.asarray(ratios)-.01))<1e-12
res=[x['signature_residual_amplitude'] for x in r['cases']['off']['rows'] if x['component']=='x']
assert max(res)<1e-8,max(res)
assert abs(r['cases']['off']['min_field_overlap']+.9)<1e-12
assert abs(r['cases']['on']['min_field_overlap']-.99)<1e-12
m.save(Path('/public/home/wfollowm/APRL_Workspace/results/energy_critical_control_20261105/design/MODE_REVIEW_IO_VALIDATION.json'),dict(passed=True,synthetic_fixture=str(root),scope='Manufactured analysis IO only; no scientific/PIC result',unclosed_plane_rejected=True,missing_final_pic_marker_rejected=True,expected_projection_ratio=.01,max_projection_ratio_error=float(max(abs(np.asarray(ratios)-.01))),max_exact_signature_residual=float(max(res)),off_overlap=-.9,on_overlap=.99,figure=str(out/'mode_comparison.png')))
print(root)
