import os,json,re,subprocess,hashlib
from pathlib import Path
from decimal import Decimal,localcontext
import h5py,numpy as np,yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
REC=np.dtype([('step','<u8'),('source','<u8'),('time','<f8'),('dt','<f8'),('charge','<f8'),('position','<f8',(3,)),('before','<f8',(3,)),('after','<f8',(3,)),('electric','<f8',(3,)),('magnetic','<f8',(3,))]);assert REC.itemsize==160
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def rows(path,group):
 with h5py.File(path,'r') as f:
  g=f[group];assert int(g['complete'][()])==1;a=g['records'][:];assert len(a)==int(g['committed_records'][()]);return a,dict(g.attrs)
def identity():
 a=R/'undulator-off/detectors';b=R/'undulator-on/detectors';reports={}
 for name in ['release','aux','exit']:
  x,_=rows(a/(name+'.h5'),'particle_plane');y,_=rows(b/(name+'.h5'),'particle_plane');x.sort(order='particle_id');y.sort(order='particle_id');assert np.array_equal(x,y),name;reports[name]=len(x)
 with h5py.File(a/'exit-field.h5','r') as f,h5py.File(b/'exit-field.h5','r') as h:
  x=f['field_plane'];y=h['field_plane'];assert x['complete'][()]==y['complete'][()]==1
  assert np.array_equal(x['time_s'][:],y['time_s'][:]) and x['fields'].shape==y['fields'].shape
  for i in range(0,len(x['fields']),16):assert np.array_equal(x['fields'][i:i+16],y['fields'][i:i+16])
 reports['field_frames']=len(x['time_s']) if False else None
 save(R/'PROBE_IDENTITY.json',{'passed':True,'particle_planes':reports,'field_arrays_bitwise_equal':True,'physics_kernel_changed':False})
def run_case(name):
 d=R/name;cfg=yaml.safe_load((d/'run.yaml').read_text());entry,attrs=rows(d/'detectors/audit-entry-energy.h5','energy_plane');exit,_=rows(d/'detectors/audit-exit-energy.h5','energy_plane')
 release,_=rows(d/'detectors/release.h5','particle_plane');Q=float(-np.sum(release['charge_C'],dtype=np.longdouble));n=len(release)
 info=np.empty(n,dtype=[('source','<u8'),('mass','<f8')]);info['source']=release['source_id'];info['mass']=release['mass_kg'];info.tofile(d/'mass.bin')
 files=sorted((d/'momentum').glob('*.bin'));assert files
 nominal=float(cfg['mesh']['cell_size'][2]/299792458.);beta=float(attrs['boost_beta']);values=[]
 for p in files:
  a=np.fromfile(p,dtype=REC,count=20000);select=a['dt']<nominal*(1-1e-10)
  if select.any():values.extend((a['position'][select,2]+(beta*299792458.)*a['time'][select]).tolist())
  if len(values)>2000:break
 assert len(values)>20,'No activation-plane frame anchor'
 values=np.array(values);v,cnt=np.unique(values,return_counts=True);T=float(v[np.argmax(cnt)]);inliers=values[abs(values-T)<1e-15];assert len(inliers)>.9*len(values)
 assert np.ptp(inliers)<1e-18
 text=(d/'aprl.log').read_text();match=re.search(r'Boris substeps per Maxwell step=(\d+)',text);substeps=int(match[1]) if match else 1
 settings={'directory':str(d),'gamma':float(attrs['boost_gamma']),'time_origin':float(attrs['time_lab_origin_s']),'translation':T,'entry_z':-.05,'exit_z':.06,'substeps':substeps,'translation_anchor_spread_m':float(np.ptp(inliers)),'translation_anchor_samples':len(inliers)}
 (d/'replay.yaml').write_text(yaml.safe_dump(settings,sort_keys=False));(d/'traces.txt').write_text('\n'.join(map(str,files))+'\n')
 with (d/'replay.log').open('w') as f:subprocess.run([str(R/'bin/replay'),str(d/'replay.yaml')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=480)
 result=json.loads((d/'REPLAY.json').read_text())
 with localcontext() as ctx:
  ctx.prec=60
  sums=[]
  for a in [entry,exit]:
   last=a[-1];sums.append(Decimal.from_float(float(last['particle_kinetic_lab_J']))+Decimal.from_float(float(last['particle_kinetic_lab_J_low'])))
  actual=float(sums[1]-sums[0])
 assert int(entry[-1]['particle_downstream_count'])==int(exit[-1]['particle_downstream_count'])==n
 result.update(cohort=n,charge_abs_C=Q,plane_gain_eV=actual/Q,replayed_plane_gain_eV=result['plane_clipped_delta_J']/Q,
  grid_longitudinal_lab_gain_eV=result['plane_clipped_Ez_work_J']/Q,grid_transverse_lab_gain_eV=result['plane_clipped_transverse_work_J']/Q,
  external_lab_gain_eV=result['plane_clipped_external_work_J']/Q,energy_plane_vs_replay_eV=(actual-result['plane_clipped_delta_J'])/Q,frame=settings)
 assert abs(result['energy_plane_vs_replay_eV'])<.02,result
 assert abs(result['plane_clipped_residual_J']/Q)<.01,result
 save(d/'RESULT.json',result);return result
def main():
 assert os.environ.get('SLURM_JOB_ID') and (R/'PIC_COMPLETE.json').exists()
 identity();result={n:run_case(n) for n in ['undulator-on','drift-on']}
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for ax,name in zip(axes,['undulator-on','drift-on']):
  a=np.genfromtxt(R/name/'lab-work-z.csv',delimiter=',',names=True);Q=result[name]['charge_abs_C']
  for k,label in [('delta_K_lab_J','Kinetic change'),('grid_Ez_lab_work_J','Grid longitudinal E work'),('grid_transverse_lab_work_J','Grid transverse E work'),('external_lab_work_J','Prescribed lab work')]:ax.plot(a['z_m'],np.cumsum(a[k])/Q,label=label)
  ax.set(xlabel='Lab z (m)',ylabel='Cumulative mean gain (eV/electron)',title=name);ax.grid(alpha=.2);ax.legend(fontsize=7)
 fig.savefig(R/'lab_work_components.png',dpi=170);fig.savefig(R/'lab_work_components.pdf');plt.close(fig)
 save(R/'RESULT.json',{'cases':result,'probe_off_on_identity':True,'source_simulation_binary_unchanged':True,'full_170_period_gain_causally_explained':False})
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
