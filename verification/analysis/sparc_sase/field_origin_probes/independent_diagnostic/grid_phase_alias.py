from pathlib import Path
import json,numpy as np,h5py,struct
r=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103');out=r/'momentum_field_origin_20261105/independent_diagnostic'
case=r/'compact-production-r1/momentum'
meta=dict(line.split(' ',1) for line in (case/'reference/metadata.txt').read_text().splitlines() if ' ' in line)
def val(s):return struct.unpack('>d',bytes.fromhex(s))[0]
dz=val(meta['spacing_bits'].split()[2]);dt=val(meta['dt_bits']);c=299792458.
w=np.genfromtxt(case/'reference/window-rank-0.csv',delimiter=',',names=True)
w0=w['origin_box_z_m'][0]+w['step'][0]*dz
with h5py.File(case/'detectors/entry-transport-energy.h5') as f:
 g=f['energy_plane'];a=g['records'][:];G=float(g.attrs['boost_gamma']);beta=float(g.attrs['boost_beta']);timeorigin=float(g.attrs['time_lab_origin_s'])
steps=np.rint(a['time_box_s']/dt).astype(int)
assert np.max(abs(a['time_box_s']/dt-steps))<1e-8
# Absolute z_box offset cancels from differences; reconstruct phase using z_box origin=0.
coordinate=-beta*c*a['time_box_s']/dz-(w0-steps*dz)/dz
expected=(steps-steps[0])*(1-beta*c*dt/dz)
observed=coordinate-coordinate[0]
windowerr=np.max(abs(w['origin_box_z_m']-(w0-w['step']*dz)))
record=dict(gamma=G,beta=beta,one_minus_beta=1-beta,dz_m=dz,dt_box_s=dt,c_dt_over_dz=c*dt/dz,dt_lab_s=dt/G,first_step=int(steps[0]),last_step=int(steps[-1]),number_step_intervals=int(steps[-1]-steps[0]),expected_coordinate_increment_per_step_cells=1-beta*c*dt/dz,expected_total_cell_phase_drift=float(expected[-1]),observed_total_cell_phase_drift=float(observed[-1]),maximum_coordinate_affine_error_cells=float(np.max(abs(observed-expected))),window_origin_reconstruction_error_m=float(windowerr),drift_over_255_intervals_cells=255*(1-beta),nyquist_phase_drift_over_255_intervals_rad=np.pi*255*(1-beta),note='Absolute phase offset omits constant boxZOrigin; phase increments and total drift independent of it. This proves sampling geometry only, not a fitted measured mode occupancy.')
(out/'GRID_PHASE_ALIAS.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
