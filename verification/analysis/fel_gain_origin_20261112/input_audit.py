"""Inspect full input mapping and conditional 6D clone structure."""
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c, e

from prepare import save


def audit(root, p):
    n = p['active_slices']
    wavelength = p['original_plan']['physics']['wavelength_m']
    with h5py.File(p['exact_input'],'r') as f:
        records = f['particles/records'][:]
    assert len(records) == n*1024
    xyz = records['plane_position_m'].reshape(n,64,16,2)
    u = records['proper_velocity'].reshape(n,64,16,3)
    transverse_duplicate = np.all(xyz == xyz[:,:,:1,:],axis=(2,3))
    momentum_duplicate = np.all(u == u[:,:,:1,:],axis=(2,3))
    theta = -2*np.pi*c*records['arrival_time_offset_s']/wavelength
    theta = theta.reshape(n,64,16)
    weights = records['macro_weight'].reshape(n,64,16)
    ne_family = weights.sum(axis=2)
    rows = []
    for harmonic in [1,2,3,4,7,8,9,10,15,16,17,20,32]:
        phase = np.exp(1j*harmonic*theta)
        family_b = (weights*phase).sum(axis=2)/ne_family
        slice_b = (weights*phase).sum(axis=(1,2))/weights.sum(axis=(1,2))
        eligible = weights.sum(axis=(1,2)) > .2*weights.sum(axis=(1,2)).max()
        rows.append(dict(harmonic=harmonic,
            family_charge_weighted_b2=float(np.sum(ne_family[eligible]*abs(family_b[eligible])**2)/ne_family[eligible].sum()),
            family_Poisson_b2=float(np.count_nonzero(ne_family[eligible])/ne_family[eligible].sum()),
            slice_mean_Ne_b2=float(np.mean(weights.sum(axis=(1,2))[eligible]*abs(slice_b[eligible])**2)),
            harmonic_scope='Conditional-family variance is a diagnostic of joint correlations, not a target radiation mode power.'))
    phase_error = 0.
    state_error = 0.
    current_error = 0.
    with h5py.File(p['genesis_input'],'r') as f:
        assert int(f['one4one'][0]) == 0 and int(f['beamletsize'][0]) == 16
        for index in range(n):
            group = f[f'slice{512+n-index:06d}']
            r = records[index*1024:(index+1)*1024]
            ph = np.exp(-2j*np.pi*c*r['arrival_time_offset_s']/wavelength)
            phase_error = max(phase_error,float(np.max(abs(ph-np.exp(1j*group['theta'][:])))))
            for axis,key in [(0,'x'),(1,'y')]:
                state_error = max(state_error,float(np.max(abs(r['plane_position_m'][:,axis]-group[key][:]))))
            for axis,key in [(0,'px'),(1,'py')]:
                state_error = max(state_error,float(np.max(abs(r['proper_velocity'][:,axis]-group[key][:]))))
            gamma = np.sqrt(1+(r['proper_velocity']**2).sum(axis=1))
            state_error = max(state_error,float(np.max(abs(gamma-group['gamma'][:]))))
            expected = r['macro_weight'].sum()*e*c/wavelength
            current_error = max(current_error,float(abs(expected-group['current'][0])/max(expected,1e-300)))
    assert phase_error < 1e-9 and state_error < 1e-10 and current_error < 1e-12
    record = dict(complete=True, particles=len(records), clone_families=n*64,
        clones_per_family=16, independent_transverse_states_per_optical_cell=64,
        exact_16fold_transverse_duplicate_fraction=float(transverse_duplicate.mean()),
        exact_16fold_momentum_duplicate_fraction=float(momentum_duplicate.mean()),
        phase_mapping_max_complex_error=phase_error,
        state_mapping_max_absolute_error_m_or_dimensionless=state_error,
        optical_cell_current_max_relative_error=current_error,
        represented_electrons=float(weights.sum()),
        source_harmonics=rows,
        interpretation='First-three-harmonic noise tests do not certify the joint transverse/phase covariance or physical short-range two-particle correlations.',
        historical_causal_control='/public/home/wfollowm/APRL_Workspace/results/energy_critical_control_20261105/bf_origin_review/loading_correlation_test_20261105/completed-review/README.md')
    save(root/'INPUT_AUDIT.json',record)
    with (root/'input_harmonics.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({k:record[k] for k in ['particles','clone_families','exact_16fold_transverse_duplicate_fraction','phase_mapping_max_complex_error']}),flush=True)
