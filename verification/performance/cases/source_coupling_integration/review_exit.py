"""Review completed exit-repair short cards without accepting full physics."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import h5py
import numpy as np


def compare_file(a, b):
    differences = {}
    with h5py.File(a, 'r') as left, h5py.File(b, 'r') as right:
        name = next(iter(left))
        x, y = left[name], right[name]
        assert x['complete'][()] == y['complete'][()] == 1
        def visit(key, ds):
            if not isinstance(ds, h5py.Dataset) or key not in ['records', 'fields', 'time_s', 'coordinates']:
                return
            other = y[key]
            assert ds.shape == other.shape and ds.dtype == other.dtype
            for begin in range(0, len(ds), 8 if key == 'fields' else 65536):
                a = ds[begin:begin+(8 if key == 'fields' else 65536)]
                b = other[begin:begin+(8 if key == 'fields' else 65536)]
                for field in ds.dtype.names or [None]:
                    if field and 'seconds' in field:
                        continue
                    aa, bb = (a[field], b[field]) if field else (a, b)
                    if np.ascontiguousarray(aa).tobytes() == np.ascontiguousarray(bb).tobytes():
                        continue
                    tag = key+('/'+field if field else '')
                    value = float(np.max(abs(np.asarray(aa, np.longdouble)-np.asarray(bb, np.longdouble))))
                    differences[tag] = max(differences.get(tag, 0), value)
        x.visititems(visit)
    return differences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('--single', default='preflight-r1')
    parser.add_argument('--dual', default='two-nodes-r1')
    parser.add_argument('--out', default='review-r1')
    args = parser.parse_args()
    root = args.root.resolve()
    single = root/args.single
    dual = root/args.dual
    assert json.loads((single/'COMPLETE.json').read_text())['program_preflight_pass']
    assert json.loads((dual/'COMPLETE.json').read_text())['program_preflight_pass']
    sys.path.insert(0, str(single/'scripts'))
    import report_discrete
    out = root/args.out
    out.mkdir(exist_ok=False)
    report = {'full_physics_accepted': False, 'cases': {}, 'old_new_detector_comparison': {},
              'single_dual_comparison': {}, 'source_snapshot_matches': {}}
    plan = json.loads((single/'plan.json').read_text())
    binary = Path(plan['binary'])
    report['binary_sha256'] = hashlib.sha256(binary.read_bytes()).hexdigest()
    assert report['binary_sha256'] == plan['binary_sha256']
    build_plan = json.loads((Path(plan['build_record'])/'plan.json').read_text())
    snapshot = Path(build_plan['source'])/'src'
    live = Path('/public/home/wfollowm/APRL_Workspace/source/aprl/src')
    for path in snapshot.iterdir():
        if path.suffix in ['.cpp', '.h']:
            report['source_snapshot_matches'][path.name] = path.read_bytes() == (live/path.name).read_bytes()
    assert all(report['source_snapshot_matches'].values())
    for batch in [single, dual]:
        for path in sorted(batch.glob('*/RESULT.json')):
            case = json.loads(path.read_text())
            if 'strict_ledger' not in case:
                continue
            value = case['strict_ledger']
            report['cases'][batch.name+'/'+path.parent.name] = value
    old = Path('/public/home/wfollowm/APRL_Workspace/results/source_coupling_integration_20261031/preflight-r4')
    histories = {}
    for kind in ['free-drift', 'undulator']:
        case = kind+'-16-on'
        d = single/case
        old_differences = {str(p.relative_to(d)): compare_file(old/case/p.relative_to(d), p)
                           for p in sorted((d/'detectors').glob('*.h5'))}
        report['old_new_detector_comparison'][kind] = old_differences
        assert not any(old_differences.values()), old_differences
        cross = {str(p.relative_to(d)): compare_file(p, dual/case/p.relative_to(d))
                 for p in [*sorted((d/'detectors').glob('*.h5')),
                           *sorted((d/'trajectory').glob('*.h5')), d/'ledger/energy-ledger.h5']}
        report['single_dual_comparison'][kind] = cross
        assert not any(v for k, v in cross.items() if not k.startswith('ledger/')), cross
        histories[kind] = {}
        for name, folder in [('before', old/case), ('after', d)]:
            rows, sums, ranks, steps = report_discrete.load(folder/'ledger', folder/'ledger/energy-ledger.h5')
            index, curves, ledger = report_discrete.calculate(rows, sums, steps)
            k = np.asarray(rows['delta_particle_kinetic_box_J'], np.longdouble)
            residual = curves['discrete_total_residual_J']
            scale = float(np.max(abs(k)))
            report['cases'][kind+'/'+name] = dict(
                final_residual_J=float(residual[-1]), maximum_absolute_residual_J=float(np.max(abs(residual))),
                particle_change_final_J=float(k[-1]), peak_particle_change_J=scale,
                final_residual_over_abs_final_particle_change_percent=float(100*residual[-1]/abs(k[-1])),
                max_abs_residual_over_peak_particle_change_percent=float(100*np.max(abs(residual))/scale))
            histories[kind][name] = (index, residual, k)
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-source-exit-review-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for ax, kind in zip(axes, histories):
        denom = max(float(np.max(abs(x[2]))) for x in histories[kind].values())
        for name, color in [('before', '#64748b'), ('after', '#c45b28')]:
            step, residual, _ = histories[kind][name]
            ax.plot(step, 100*residual/denom, label=name, color=color)
        ax.set(title=kind, xlabel='Field step', ylabel='Signed residual / fixed common D (%)')
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle('Longitudinal exit repair: discrete energy balance')
    fig.savefig(out/'tail_balance.png', dpi=170)
    plt.close(fig)
    (out/'REVIEW.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report['cases'].items() if '/before' in k or '/after' in k}, indent=2))
    print('single-dual non-timing ledger differences', {k: v['ledger/energy-ledger.h5'] for k, v in report['single_dual_comparison'].items()})


if __name__ == '__main__':
    main()
