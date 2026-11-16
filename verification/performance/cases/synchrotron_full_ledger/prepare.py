"""Prepare the two original 100k synchrotron samples with self-consistent feedback."""
import copy,json,math,shutil,sys
from pathlib import Path
import h5py,yaml
from common import sha,save,duration_steps

W=Path('/public/home/wfollowm/APRL_Workspace');T=Path(__file__).resolve().parent
def main(root):
 root=Path(root).resolve();root.mkdir(exist_ok=False);shutil.copytree(T,root/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
 proof=W/'results/sparc_six_full_ledger_20261112_r3';gate=json.loads((proof/'PREFLIGHT_COMPLETE.json').read_text())
 assert gate['complete'] and sha(gate['binary'])==gate['binary_sha256']
 (root/'bin').mkdir();binary=root/'bin/aprl';shutil.copy2(gate['binary'],binary);assert sha(binary)==gate['binary_sha256']
 metadata=root/'metadata';metadata.mkdir()
 for name in ['PLAN.json','PREFLIGHT_COMPLETE.json','STARTUP_CONFIRMED.json','FROZEN.json']:shutil.copy2(proof/name,metadata/('reused_program_'+name))
 cases=[];inputs={};protected={};audit=[];first_card=None
 for shape in ['linear-binomial3','quintic-compensated']:
  for seed in [20261111,20261112]:
   previous=W/f'results/synchrotron_independent_electron_100k_seed{seed}_20261110_r1'
   campaign=json.loads((previous/'campaign.json').read_text());input_file=Path(campaign['input'])
   assert sha(input_file)==campaign['input_sha256']
   with h5py.File(input_file,'r') as f:assert len(f['particles/records'])==100000
   inputs[str(input_file)]=campaign['input_sha256']
   for kind in ['bend','straight']:
    source=previous/kind;old=yaml.safe_load((source/'run.yaml').read_text());original_plan=json.loads((source/'PLAN.json').read_text())
    d=root/shape/str(seed)/kind;d.mkdir(parents=True)
    for sub in ['detectors','ledger','worldtube','lab-audit']:(d/sub).mkdir()
    card=copy.deepcopy(old);card['output']['manifest']=str(d/'run-manifest.yaml')
    card['mesh'].update(longitudinal_shape=shape,self_field_feedback=True,particle_coupling='conservative-source')
    card['mesh'].pop('one_way_source_diagnostic',None)
    card['mesh']['moving_window']['reference_directory']=str(d/'reference')
    # The proven side-particle observer uses a rectangle centred at x=y=0.
    # Translate grid and particle coordinates together; this preserves every
    # relative source/grid coordinate and the spatially uniform dipole field.
    shift=card['mesh']['center'][0];assert shift==.00012 and card['mesh']['center'][1]==0
    card['mesh']['center'][0]=0.;card['beam']['input']['position_offset'][0]-=shift
    assert all(m['type']=='uniform-dipole' or (m['type']=='planar-undulator' and m['strength_parameter']==0.) for m in card['sources']['magnetic_elements'])
    assert card['beam']['input']['position_offset'][0]-card['mesh']['center'][0]==old['beam']['input']['position_offset'][0]-old['mesh']['center'][0]
    for key in ['cell_size','cells','duration','boost_gamma','particle_steps_per_undulator_period','maximum_particle_substeps']:
     assert card['mesh'][key]==old['mesh'][key]
    assert card['boundary']==old['boundary'] and card['sources']==old['sources'] and card['stop']==old['stop']
    assert card['beam']['input']['file']==old['beam']['input']['file'] and card['beam']['input']['electrons']==100000.
    planes=[('entry',0.),('audit-base',.02),('magnet-entry',.05),('bend-exit',.0524),('z0060-particles',.06),('z0100-particles',.1),('z0200-particles',.2),('exit',.3524)]
    card['detectors'].update(directory=str(d/'detectors'),substep_particle_planes=True,
       particle_planes=[dict(name=name,z=z,buffer_records=8192,compression=0) for name,z in planes])
    gamma=card['mesh']['boost_gamma'];dt=card['mesh']['cell_size'][2]/299792458.
    card['detectors']['energy_planes']=[dict(name=name+'-transport',z=z,rhythm=dt/gamma,buffer_samples=64,compression=0) for name,z in planes]
    assert card['detectors']['field_planes']==old['detectors']['field_planes']
    card['energy_ledger'].update(enabled=True,directory=str(d/'ledger'),sample_interval_steps=1,discrete_operator_audit=True)
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    case=dict(name=f'{shape}-seed{seed}-{kind}',shape=shape,seed=seed,kind=kind,directory=str(d),input=str(input_file),input_sha256=campaign['input_sha256'],
       particle_count=100000,legacy_expected_steps=original_plan['expected_steps'],duration_step_budget=duration_steps(card)[0],timeout_seconds=5400,original_case=str(source),lab_base_z_m=.02,
       restore_original_lab_x_by_adding_m=shift,feedback=True)
    save(d/'CASE.json',case);cases.append(case)
    if first_card is None:first_card=card
    protected[str(source/'run.yaml')]=sha(source/'run.yaml')
    audit.append(dict(case=case['name'],input_byte_identical=True,geometry_and_magnet_and_stop_preserved=True,
       feedback_changed_from=False,feedback_changed_to=True,source_shape=shape,
       transverse_coordinate_origin_shift_m=-shift,relative_grid_particle_coordinates_unchanged=True,
       original_output_x_recovered_by_adding_m=shift,all_frequency_ledgers=True,no_record_windowing=True))
 preflight=[]
 gamma=first_card['mesh']['boost_gamma'];dt=first_card['mesh']['cell_size'][2]/299792458.
 for name,shape,observer in [('three-point-off','linear-binomial3',False),('three-point-on','linear-binomial3',True),('momentum-on','quintic-compensated',True)]:
  d=root/'preflight'/name;d.mkdir(parents=True)
  for sub in ['detectors','ledger','worldtube','lab-audit']:(d/sub).mkdir()
  base=Path(cases[0]['directory']);card=yaml.safe_load((base/'run.yaml').read_text().replace(str(base),str(d)))
  card['mesh']['longitudinal_shape']=shape;card['mesh']['duration']=40.5*dt
  (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
  preflight.append(dict(name=name,shape=shape,observer=observer,directory=str(d),expected_steps=duration_steps(card)[0],
   duration_step_budget=duration_steps(card)[0],particle_count=100000,timeout_seconds=600))
 plan=dict(complete_preparation=True,binary=str(binary),binary_sha256=sha(binary),binding_probe=str(W/'build/cluster-tools/binding_probe'),
   reused_observer_validation=str(proof/'PREFLIGHT_COMPLETE.json'),mpi=16,threads=8,exclusive_cpus_per_job=128,
   cases=cases,preflight=preflight,input_files=inputs,self_field_feedback=True,
   scope='Two current source shapes; original two 100k samples, bend plus matching straight control. Self-consistent feedback explicitly authorized by user.',
   ledger_scope=['native boost ledger each step','strict discrete Cowan energy ledger','lab energy/momentum/stress','whole Cowan exchanges',
     'full-frequency fixed laboratory energy planes','four transverse field faces','two null-time caps','committed-substep side-particle kinetic transport','substep scalar lab-plane particle energy'],
   lab_planes_m=[0.,.02,.05,.0524,.06,.1,.2,.3524],worldtube_z_interval_m=[.02,.3524],
   physical_half_width_m=[(160/2-11)*5e-6,(128/2-11)*5e-6],
   source_support_upper_bound_m=9*gamma*first_card['mesh']['cell_size'][2],
   worldtube_source_exclusion='0.02m starts downstream of frozen-source support; transverse surface-current sheets are one cell outside observer faces; devices are static lab magnetic fields.',
   no_CXX_changes=True,no_build=True,no_default_changes=True,new_physics_results_accepted=False,
   queued_short_preflight_required=True,preflight_duration_in_dt=40.5,
   preflight_expected_steps=duration_steps(card)[0],
   step_validation='Preflight follows binary64 repeated-addition duration bound; formal runs require physical stop and complete mutually consistent records.',
   prior_failed_campaign=str(W/'results/synchrotron_full_ledger_20261112_r1'))
 assert plan['source_support_upper_bound_m']<.02
 save(root/'PLAN.json',plan);save(root/'CARD_AUDIT.json',dict(complete=True,cases=audit))
 protected[str(W/'cluster.env')]=sha(W/'cluster.env');save(root/'PROTECTED.json',protected)
 paths=[root/'PLAN.json',binary,*[Path(c['directory'])/'run.yaml' for c in cases+preflight],*[p for p in (root/'scripts').iterdir() if p.is_file()]]
 save(root/'FROZEN.json',{str(p):sha(p) for p in paths}|inputs)
 print(json.dumps(dict(cases=len(cases),preflight=len(preflight),binary_sha256=plan['binary_sha256'],source_support_upper_bound_m=plan['source_support_upper_bound_m']),indent=2))
if __name__=='__main__':main(sys.argv[1])
