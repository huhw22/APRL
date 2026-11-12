"""Freeze an independent-6D first-segment campaign; no simulation here."""
import argparse
import copy
import shutil
import sys
from pathlib import Path

from common import W, load, save, sha, require


def main(root):
    root = Path(root).resolve()
    root.mkdir(exist_ok=False)
    template = Path(__file__).resolve().parent
    shutil.copytree(template, root/'scripts', ignore=shutil.ignore_patterns('__pycache__'))
    cases = template.parent
    helpers = {'source_math.py':cases/'physical_electron_source/source_math.py',
               'source_verify.py':cases/'physical_electron_source/verify.py',
               'phase_map.py':cases/'physical_electron_source/genesis.py',
               'loading_guard.py':cases/'sparc_sase/loading_guard.py',
               'binding.py':cases/'source_coupling_integration/binding.py',
               'analyze_closed.py':cases/'source_coupling_integration/analyze_closed.py'}
    for name, path in helpers.items(): shutil.copy2(path, root/'scripts'/name)
    sys.path.insert(0, str(cases/'physical_electron_source'))
    from source_math import nominal_electrons, C
    old = W/'results/sparc_sase_11segment_2slip_20261111_r1'
    prior = load(old/'plan.json')
    cfg = load(cases/'physical_electron_source/example_37periods.json')
    wavelength = prior['physics']['wavelength_m']
    cycle = wavelength/C
    first = (prior['active_slices']-231)//2
    cfg['profile'].update(full_duration_s=prior['physics']['duration_s'],
        window_start_s=first*cycle, window_stop_s=(first+231)*cycle,
        bin_width_s=cycle, peak_current_A=prior['physics']['peak_current_A'])
    cfg['electrons'] = round(nominal_electrons(cfg))
    cfg['seed'] = 2026111201
    binary = Path(prior['binary'])
    require(sha(binary) == prior['binary_sha256'], 'Prior binary differs')
    genesis = Path(load(W/'results/fel_gain_origin_20261112_r1/PLAN.json')['cases'][0]['binary'])
    for name in ('inputs', 'preparation'): (root/name).mkdir()
    entries = [dict(name='w512-a',target_weight=512,physical_seed=2026111201),
               dict(name='w512-b',target_weight=512,physical_seed=2026111202),
               dict(name='w64-a',target_weight=64,physical_seed=2026111201)]
    for entry in entries:
        directory = root/entry['name'];directory.mkdir()
        for folder in ('detectors','ledger','genesis-off','genesis-longrange'):
            (directory/folder).mkdir()
        entry['directory'] = str(directory)
        entry['timeout_seconds'] = 4*3600
    plan = dict(root=str(root),workspace=str(W),source_config=cfg,cases=entries,
        pic_binary=str(binary),pic_binary_sha256=sha(binary),
        genesis_binary=str(genesis),genesis_binary_sha256=sha(genesis),
        prior_campaign=str(old),boost_gamma=prior['physics']['boost_gamma'],
        wavelength_m=wavelength,period_m=.028,K=prior['physics']['K'],
        periods=77,source_optical_cells=231,core_generation_bins=[77,154],
        particle_planes=[['entry',0.],['und-entry',.01],
                         *[[f'u{n:02d}',.01+n*.028] for n in (7,19,37,57,77)]],
        weighted_sampling='Uniform fixed-size sampling without replacement from a conceptual independent fixed-N physical source; w512-a is nested in w64-a. Only selected records are stored.',
        source_generation='Unchanged validated Philox/source_math.py, separate physical and thinning seeds; no copied transverse states, phase adjustment or moment fitting.',
        execution=dict(pic_ranks=16,pic_threads=8,genesis_ranks=32,genesis_threads=1),
        scope='First-segment loading/growth check; weighted entrance noise is retained, no amplitude normalization, no saturation claim.')
    save(root/'PLAN.json',plan)
    paths = [root/'PLAN.json',binary,genesis,
             old/'energy-binomial3/run.yaml', *[p for p in (root/'scripts').iterdir() if p.is_file()]]
    save(root/'FROZEN_PREPARATION.json',{str(p):sha(p) for p in paths})
    print(root)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    main(parser.parse_args().output)
