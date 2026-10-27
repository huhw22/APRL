#!/usr/bin/env python3
"""Review bounded window-state/overhead probes, keeping physical scopes separate."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import math


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cost',type=Path,required=True)
    parser.add_argument('--wave',type=Path)
    parser.add_argument('--controls',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require((args.cost/'COMPLETED.json').is_file(),'Cost campaign incomplete')
    plan=json.loads((args.cost/'plan.json').read_text())
    for path,expected in json.loads((args.cost/'frozen.json').read_text()).items():
        require(digest(Path(path))==expected,'Frozen input changed: '+path)
    verification=json.loads((args.cost/'verify/metrics.json').read_text())
    require(verification['completed'] and verification['bitwise_failures']==0,'State migration failed')
    coverage=(verification['tagged_shift_cases'],verification['particle_schedules'])
    require(coverage in ((126,4),(252,8)),'Missing migration coverage')
    if coverage==(252,8):
        policies=verification.get('history_policies',[])
        require(len(policies)==2 and {p['policy'] for p in policies}==
                {'world-position','boundary-relative-z'},'Missing history-policy coverage')
        require(all(p['tagged_shift_cases']==126 and p['particle_schedules']==4 for p in policies),
                'Incomplete per-policy migration coverage')
        require(all(verification.get(key) is True for key in
                    ('zero_shift_policy_identity','excited_z_history_policy_difference',
                     'next_step_policy_field_difference')),'History-policy distinction check failed')
    result={'format':'aprl-moving-window-review-v1','source_base':plan['source_base'],
            'binary_sha256':plan['binary_sha256'],'state_migration':'PASS',
            'verification':verification,'costs':[],
            'scope':'Unchanged kernels, measured outer-layer migration; not production speedup'}
    for width in (64,128,256):
        directory=args.cost/f'bench{width}'
        m=json.loads((directory/'metrics.json').read_text())
        require(m['completed'] and m['numerical_status']=='PASS','Benchmark gate failed')
        with (directory/'timing.csv').open() as f:
            rows=list(csv.DictReader(f))
        require(len(rows)==9,'Incomplete rotated repeats')
        for row in rows:
            for key,value in row.items():
                if key!='mode':require(math.isfinite(float(value)),'Nonfinite benchmark value')
            require(int(row['width'])==width and int(row['steps'])==64,'Changed benchmark dimensions')
            require(abs(float(row['total_s'])-float(row['kernel_s'])-float(row['move_s']))<.002,'Timing components inconsistent')
        groups={mode:[r for r in rows if r['mode']==mode] for mode in ('fixed','compact','snapshot')}
        for group in groups.values():
            require(sorted(int(r['repeat']) for r in group)==[0,1,2],'Missing repeat or mode')
        modes={}
        for mode,group in groups.items():
            modes[mode]={key:statistics.median(float(r[key]) for r in group) for key in group[0] if key!='mode'}
            modes[mode]['total_s_min']=min(float(r['total_s']) for r in group)
            modes[mode]['total_s_max']=max(float(r['total_s']) for r in group)
        paired=[]
        for repeat in range(3):
            f=next(r for r in groups['fixed'] if int(r['repeat'])==repeat)
            c=next(r for r in groups['compact'] if int(r['repeat'])==repeat)
            paired.append(float(c['total_s'])/float(f['total_s'])-1)
        result['costs'].append({'width':width,'modes':modes,
            'same_width_total_overhead_fraction':modes['compact']['total_s']/modes['fixed']['total_s']-1,
            'same_width_native_kernel_difference_fraction':modes['compact']['kernel_s']/modes['fixed']['kernel_s']-1,
            'movement_over_native_kernel_fraction':modes['compact']['move_s']/modes['compact']['kernel_s'],
            'paired_total_overhead_fractions':paired,
            'preflight_bitwise_steps':m['bitwise_preflight_steps']})
    if args.wave:
        require((args.wave/'COMPLETED.json').is_file(),'Wave campaign incomplete')
        waveplan=json.loads((args.wave/'plan.json').read_text())
        require(waveplan['binary_sha256']==plan['binary_sha256'],'Different executable in wave test')
        directory=args.wave/'wave';m=json.loads((directory/'metrics.json').read_text())
        with (directory/'wave-plane.csv').open() as f:rows=list(csv.DictReader(f))
        require(len(rows)==177 and [int(r['step']) for r in rows]==list(range(177)),'Wave frames missing/repeated')
        for row in rows:
            require(all(math.isfinite(float(v)) for v in row.values()),'Nonfinite wave sample')
            n=int(row['step'])
            require(int(row['window_origin_global_k'])==-((4*n+4)//5),'Wrong integer moving origin')
        for prefix,metric in [('E','E_relative_l2_error'),('cB','B_relative_l2_error')]:
            ref=sum(float(r[prefix+'_reference_rms_V_m'])**2 for r in rows)
            err=sum(float(r[prefix+'_difference_rms_V_m'])**2 for r in rows)
            measured=math.sqrt(err/ref)
            require(abs(measured-m[metric])<1e-12*max(1,measured),'Wave norm independent recomputation failed')
        require(m['initial_overlap_bitwise_failures']==0 and m['finite'],'Wave initial/finite gate failed')
        result['vacuum_wave']=m
    if args.controls:
        require((args.controls/'COMPLETED.json').is_file(),'CPML controls incomplete')
        for path,expected in json.loads((args.controls/'frozen.json').read_text()).items():
            require(digest(Path(path))==expected,'Frozen control changed: '+path)
        common_reference=None;control_results={}
        for name in ('wave','wave-static','wave20','wave-static20'):
            directory=args.controls/name;m=json.loads((directory/'metrics.json').read_text())
            with (directory/'wave-plane.csv').open() as f:rows=list(csv.DictReader(f))
            require(len(rows)==177 and [int(row['step']) for row in rows]==list(range(177)),'Control frames missing')
            reference=[tuple(row[k] for k in row if 'reference' in k or k in ('step','time_E_box_s','time_B_box_s','detector_z_box_m')) for row in rows]
            if common_reference is None:common_reference=reference
            require(reference==common_reference,'Common large-domain trace changed')
            for prefix,metric in [('E','E_relative_l2_error'),('cB','B_relative_l2_error')]:
                measured=math.sqrt(sum(float(row[prefix+'_difference_rms_V_m'])**2 for row in rows)/
                                   sum(float(row[prefix+'_reference_rms_V_m'])**2 for row in rows))
                require(abs(measured-m[metric])<1e-12,'Independent control norm differs')
            control_results[name]={'metrics':m,
                'first_100_frames_max_difference_rms':max(float(row['combined_difference_rms_V_m']) for row in rows[:100])}
        if args.wave:
            require((args.controls/'wave/wave-plane.csv').read_bytes()==(args.wave/'wave/wave-plane.csv').read_bytes(),
                    'Default wave changed after adding controls')
        result['cpml_controls']={'cases':control_results,'common_reference_identical':True,
            'original_wave_reproduced_bitwise':bool(args.wave),
            'binary_sha256':json.loads((args.controls/'plan.json').read_text())['binary_sha256'],
            'interpretation':'Finite-time late responses at different boundary trajectories; not isolated reflection coefficients'}
    result['limitations']=[
        'Moving CPML changes physical support; state remap equality alone does not accept boundary physics.',
        'The 256-particle equivalence fixture covers direct core calls, not complete Simulation held/MPI/checkpoint flows.',
        'Each benchmark has eight bitwise preflight steps; timing uses 64 steps after four untimed warm-up steps.',
        'Synthetic inflow cost is included; production reference-field generation/queries and MPI are not included.',
        'One CPU host and single-rank timings are not a distributed production measurement.',
        'Different-width cost observations are not same-physics implementation speedups.'
    ]
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__':
    main()
