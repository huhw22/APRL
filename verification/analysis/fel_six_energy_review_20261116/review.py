"""Read-only full-ledger review of completed independent-w355 momentum case."""
import json
import sys
import subprocess
from pathlib import Path
import numpy as np
from report_discrete import load, calculate
from relative_ledger import construct, write_tables
import whole
import lab

LD = np.longdouble


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main(parent, out):
    parent, out = Path(parent).resolve(), Path(out).resolve()
    case = 'energy-momentum'
    d = parent / case
    assert json.loads((d / 'aprl.log.launch.json').read_text())['exit_code'] == 0
    assert 'phase=mpi-finalize-complete' in (d / 'aprl.log').read_text()
    identities = json.loads((out / 'PROTECTED_OUTPUTS.json').read_text())
    def unchanged():
        for name, before in identities.items():
            s = Path(name).stat()
            assert before == dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino), name
    unchanged()
    rows, sums, ranks, steps = load(d / 'ledger', d / 'ledger/energy-ledger.h5')
    assert ranks == 16 and steps == 28199 and len(rows) == 28200
    index, quantities, summary = calculate(rows, sums, steps)
    absolute, pct, normalization = construct(rows, sums, index, quantities)
    boost = out / 'boost'
    write_tables(boost, index, absolute, pct, normalization)
    np.savez_compressed(boost / 'cache.npz', records=rows, **absolute)
    gamma, beta = LD(rows['boost_gamma'][0]), LD(rows['boost_beta'][0])
    lorentz = rows['delta_particle_kinetic_lab_J'].astype(LD) - gamma * (
        rows['delta_particle_kinetic_box_J'].astype(LD)
        + beta * rows['delta_particle_cPz_box_J'].astype(LD))
    boost_result = dict(complete=True, steps=steps,
        final_particle_change_box_J=float(rows['delta_particle_kinetic_box_J'][-1]),
        final_particle_change_lab_full_cohort_J=float(rows['delta_particle_kinetic_lab_J'][-1]),
        max_particle_Lorentz_identity_error_J=float(np.max(abs(lorentz))),
        normalization=normalization, **summary)
    save(boost / 'RESULT.json', boost_result)
    print('BOOST', json.dumps(boost_result), flush=True)

    values, numerical = whole.one(parent, out / 'whole', case)
    # Independent bridge between the full Lorentz gap and measured particle /
    # deposited-current energy and longitudinal momentum components.
    current_energy = np.r_[LD(0), np.cumsum(
        np.load(out / 'whole' / case / 'whole_cache.npz')['fixed_end_WJ_particle_J'], dtype=LD)]
    current_momentum = np.r_[LD(0), np.cumsum(
        np.load(out / 'whole' / case / 'whole_cache.npz')['fixed_end_PJ_particle_J'], dtype=LD)]
    energy_gap = gamma * (rows['delta_particle_kinetic_box_J'].astype(LD) - current_energy)
    momentum_gap = gamma * beta * (rows['delta_particle_cPz_box_J'].astype(LD) + current_momentum)
    error = np.max(abs(energy_gap + momentum_gap - values['particle_vs_deposited_current_Lorentz_gap_J']))
    assert error < LD('1e-13')
    numerical['independent_particle_energy_momentum_bridge_max_J'] = float(error)
    numerical['particle_energy_part_final_J'] = float(energy_gap[-1])
    numerical['particle_longitudinal_momentum_part_final_J'] = float(momentum_gap[-1])
    save(out / 'whole' / case / 'WHOLE_RESULT.json', numerical)
    print('WHOLE', json.dumps(numerical), flush=True)

    laboratory = out / 'lab'
    laboratory.mkdir()
    binary = laboratory / 'worldtube_integrator'
    subprocess.run(['c++', '-std=c++17', '-O3', '-fopenmp',
                    str(Path(__file__).with_name('worldtube.cpp')), '-o', str(binary)], check=True)
    lab.geometry_test(laboratory, binary)
    _, _, physical = lab.analyze_case(parent, laboratory, case, binary)
    unchanged()
    result = dict(complete=True, physics_job='19796_1', case=case,
        independent_random_source=True, no_new_physics=True,
        boost=boost_result, lab_on_boost_slices=numerical, fixed_lab_worldtube=physical,
        original_output_identities_unchanged=len(identities))
    save(out / 'RESULT.json', result)
    print('REVIEW_COMPLETE', flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])
