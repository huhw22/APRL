"""Require a closed feedback run before starting the long radiation calculation.

This only screens gross early field growth; it is not a convergence test or a
radiation-energy measurement. The native force/work tolerances remain active.
"""
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

case = Path(sys.argv[1]).resolve()
gate = case / 'stability-gate'
steps = json.loads((case / 'PARAMETERS.json').read_text())['stability_gate_steps']
assert int((gate / 'EXIT_CODE').read_text()) != 0
log = (gate / 'aprl.log').read_text()
assert 'mesh.duration was exhausted before the configured physical stop;' in log
assert f'Moving window completed: steps={steps} particle_exits=0 cropped_carriers=0' in log
paths = sorted((gate / 'ledger').glob('discrete-field-rank-*.csv'))
assert len(paths) == 16
arrays = []
for path in paths:
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert [int(r['step']) for r in rows] == list(range(1, steps + 1)), path
    assert all(math.isfinite(float(v)) for r in rows for v in r.values()), path
    arrays.append(rows)
initial = math.fsum(float(rows[0]['native_begin_J']) for rows in arrays)
energy = np.array([math.fsum(float(rows[i]['native_end_J']) for rows in arrays)
                   for i in range(steps)])
assert initial > 0 and np.all(energy > 0)
growth = float(energy.max() / initial)
tail_ratio = float(np.median(energy[-32:] / energy[-33:-1]))
# Deliberately loose alarms for the catastrophic inlet mode seen in the old
# runs. Physical beam entry can increase the field norm by many orders.
assert growth < 1.e10, ('gross field-norm growth', growth)
assert tail_ratio < 1.2, ('sustained rapid field-norm growth', tail_ratio)
coupling = sorted((gate / 'reference').glob('source-coupling-rank-*.csv'))
assert len(coupling) == 16
max_force = max_work = 0.
max_trials = 0
for path in coupling:
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert [int(r['step']) for r in rows] == list(range(1, steps + 1)), path
    for row in rows:
        force = float(row['self_force_defect'])
        work = float(row['relative_fixed_point_work_defect'])
        assert math.isfinite(force) and force <= 1.e-13
        assert math.isfinite(work) and work <= 1.e-8
        max_force = max(max_force, force)
        max_work = max(max_work, work)
        max_trials = max(max_trials, int(row['iterations']))
result = dict(stability_gate_passed=True, recorded_steps=steps,
              expected_duration_guard=True, physical_run_complete=False,
              initial_native_field_norm_J=initial,
              last_native_field_norm_J=float(energy[-1]),
              peak_field_norm_over_initial=growth,
              median_last_32_step_norm_energy_ratio=tail_ratio,
              maximum_self_force_defect=max_force,
              maximum_relative_fixed_point_work_defect=max_work,
              maximum_coupling_trials=max_trials,
              scope='Closed early self-consistent field-growth screen; no complete-run or mesh-convergence claim. Native staggered field norm is not radiation energy.')
(gate / 'PASS.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2), flush=True)
