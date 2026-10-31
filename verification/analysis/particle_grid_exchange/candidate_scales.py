#!/usr/bin/env python3
"""Inspectable scale estimates; no assertion that fitted density is a PIC peak."""
import hashlib,json,math,re,sys
from pathlib import Path
import h5py,numpy as np,yaml
ws=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();L=np.longdouble
statistics=ws/'results/corrected_full_comparison_20261030/undulator_review_19099_r1/science/particle_statistics.json'
phase_path=ws/'results/particle_grid_exchange_diagnosis_20261031/r1/SCALES.json'
full=ws/'results/corrected_full_comparison_20261030/r3_2nodes/undulator_4p5m'
p=json.loads(statistics.read_text())['entry_all_input_cohort'];phase=json.loads(phase_path.read_text())['plane_wave_phase_scale']
with h5py.File(full/'ledger/energy-ledger.h5') as f:
 g=f['energy_ledger'];assert int(g['complete'][()])==1
 assert len(g['records'])==int(g['committed_records'][()]);mean_gamma=g['records']['mean_gamma_box']
card=yaml.safe_load((full/'run.yaml').read_text());G=L(str(card['mesh']['boost_gamma']));lab=L(str(p['mean_gamma']))
b=L(np.sqrt(1-1/lab**2));B=L(np.sqrt(1-1/G**2));factor=1/(G*(1-b*B))
sx=L(str(p['rms_x_um']))*L('1e-6');sy=L(str(p['rms_y_um']))*L('1e-6');sz=b*L(str(p['rms_arrival_length_um']))*L('1e-6')*factor
n=L(str(p['represented_electrons']))/((2*L(str(math.pi)))**L('1.5')*sx*sy*sz)
e=L('1.602176634e-19');eps=L('8.8541878128e-12');mass=L('9.1093837015e-31');dt=L(str(phase['field_dt_s']))
wp=np.sqrt(n*e*e/(eps*mass));N=L(str(phase['Maxwell_steps_per_boost_cycle']))
line=next(x for x in (full/'aprl.log').read_text().splitlines() if '[resource] phase=pre-run' in x and 'grid_cells=' in x)
cells=list(map(int,re.search(r'grid_cells=(\d+)x(\d+)x(\d+)',line).groups()));nodes=math.prod(x+1 for x in cells)
data={'source_statistics':str(statistics),'source_statistics_sha256':hashlib.sha256(statistics.read_bytes()).hexdigest(),
 'input_source_cohort':'all 41732775 input particles, not only matched exits',
 'computational_frame':{'boost_gamma':G,'full_19099_recorded_mean_gamma_min':np.min(mean_gamma),'full_19099_recorded_mean_gamma_max':np.max(mean_gamma),'warning':'These are per-record means, not individual velocity extrema or mean beta; individual extrema from short-case observers are separate.'},
 'literal_ECSIM_nodal_storage_estimate':{'source_geometry_line':line,'global_active_cells':cells,'unique_global_nodes':nodes,'mass_blocks_per_node':14,'entries_per_block':9,'bytes_per_entry':8,'mass_storage_GiB':nodes*14*9*8/2**30,'one_nodal_three_component_work_vector_GiB':nodes*3*8/2**30,'scope':'Literal dense block layout in the cited implementation, not an algorithmic lower bound; MPI ghosts, padding and solver vectors excluded.'},
 'uncorrelated_rigid_Gaussian_fit_estimate':{'assumptions':['One constant speed derived from mean lab gamma','Arrival length is c times rms arrival time; convert to simultaneous boost length','Uncorrelated Gaussian spatial fit','No microbunch peaks or PIC-deposited density inferred'],'simultaneous_boost_sigma_z_m':sz,'sigma_z_formula':'beta_lab * c_sigma_arrival_time / [Gamma * (1-beta_lab*beta_boost)]','peak_density_fit_m_minus3':n,'classical_omega_p_dt':wp*dt,'unmagnetized_source_response_dt2_omega_p2_over4':(wp*dt)**2/4,'interpretation':'Small fitted source response motivates measuring a source-only solve; not a bound on actual density, its matrix spectrum or iteration cost. Does not bound the full Maxwell-curl linear solve.'},
 'one_dimensional_standard_difference_CN_illustration':{'steps_per_forward_boost_photon_cycle':N,'CN_phase_ratio_at_cdt_equals_dz':np.arctan(np.sin(np.pi/N))/(np.pi/N),'relative_phase_error_percent':100*(np.arctan(np.sin(np.pi/N))/(np.pi/N)-1),'scope':'Derived 1D illustration with same spatial difference, not measured 3D Cowan error.'}}
(out/'CANDIDATE_SCALES.json').write_text(json.dumps(data,indent=2,default=lambda x:float(x))+'\n')
print(json.dumps(data,indent=2,default=lambda x:float(x)))
