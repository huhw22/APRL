#!/usr/bin/env python3
"""Candidate-only, bounded gnode1 workload/layout experiments; no core changes."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import time

import yaml

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
BINARY_SHA = '56cc1e0931103687dcc440c3ad390a923b2e74f5da0d68ef1806e3992ccbfeac'
CORE = 'ac58e67e33d6b4ecfd12c3d00f4a9947d0e8500d'
LAYOUTS = [(16, 1, 'compact'), (8, 1, 'compact'), (4, 1, 'compact'),
           (2, 1, 'compact'), (1, 1, 'compact'), (8, 2, 'compact'),
           (4, 4, 'compact'), (2, 8, 'compact'), (1, 16, 'compact'),
           (8, 1, 'spread')]
PROFILE_LAYOUTS = [(16, 1, 'compact'), (4, 4, 'compact'), (1, 16, 'compact')]
STAGES = ('field_total', 'particle_prepare', 'particle_deposit', 'particle_migration_total',
          'diagnostics', 'control', 'field_communication', 'particle_communication')
CASE_SPECS = [('grid016', 16, 16384), ('grid032', 32, 16384),
              ('grid064', 64, 16384), ('grid128', 128, 16384),
              ('particles001024', 32, 1024), ('particles004096', 32, 4096),
              ('particles065536', 32, 65536), ('particles262144', 32, 262144)]
SEED = 20261025
RUN_LIMIT_S = 180
TOTAL_LIMIT_S = 7200


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def layout_name(layout):
    ranks, threads, placement = layout
    return f'r{ranks}-t{threads}-{placement}'


def mpi_command(layout, executable, arguments):
    ranks, threads, placement = layout
    mapping = 'ppr:1:numa' if placement == 'spread' else (
        'node:PE=16' if threads == 16 else f'slot:PE={threads}')
    command = ['mpirun', '--host', 'gnode1:16', '-np', str(ranks),
               '--map-by', mapping, '--bind-to', 'core', '--nooversubscribe',
               '--report-bindings']
    if placement == 'spread': command += ['--rank-by', 'numa']
    return command + [str(executable), *map(str, arguments)]


def environment(threads, profile):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
               OMP_PROC_BIND='close', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', APRL_PROFILE_INTERVAL=str(profile))
    return env


def prepare(workspace, root):
    root.mkdir(parents=True, exist_ok=False)
    scripts = root / 'scripts'; scripts.mkdir()
    cards = root / 'cards'; cards.mkdir()
    here = Path(__file__).resolve().parent
    for name in ('run.py', 'run.sh', 'README.md'):
        shutil.copy2(here / name, scripts / name)
    case_root = workspace / 'source/aprl/verification/performance/cases'
    shutil.copy2(case_root / 'generic_particles/run.py', scripts / 'generic_inspector.py')
    shutil.copy2(workspace / 'source/aprl/verification/performance/run.py', scripts / 'basic_inspector.py')
    base = yaml.safe_load((case_root / 'cpu_calibration/cards/balanced.yaml').read_text())
    base['beam']['input']['sigma_position'] = [.25, .25, .25]
    base['runtime']['resource_monitor']['progress_interval_steps'] = 20
    base['runtime']['shutdown_watchdog'] = dict(enabled=True, timeout_seconds=600,
                                               finalize_timeout_seconds=60)
    cases = []
    for name, transverse, particles in CASE_SPECS:
        card = copy.deepcopy(base)
        card['mesh']['cells'] = [transverse, transverse, 128]
        card['beam']['input']['macroparticles'] = particles
        file = cards / (name + '.yaml'); file.write_text(yaml.safe_dump(card, sort_keys=False))
        cases.append(dict(name=name, card=str(file), sha256=digest(file), cells=card['mesh']['cells'],
                          particles=particles, grid_cells=transverse**2 * 128,
                          represented_electrons=card['beam']['input']['electrons']))
    binary = workspace / 'build/candidate/aprl'
    probe = workspace / 'build/node-survey/binding_probe'
    if digest(binary) != BINARY_SHA: raise ValueError('Accepted candidate binary differs')
    if not probe.is_file(): raise ValueError('Prepared standalone binding probe missing')
    plan = dict(format='aprl-memory-scaling-v1', workspace=str(workspace), root=str(root),
        prepared_at=datetime.now(timezone.utc).isoformat(),
        source_revision=subprocess.check_output(['git', '-C', str(workspace/'source/aprl'), 'rev-parse', 'HEAD'], text=True).strip(),
        binary=str(binary), binary_sha256=BINARY_SHA, compiled_core_revision=CORE,
        binding_probe=str(probe), binding_probe_sha256=digest(probe), cases=cases,
        layouts=[dict(ranks=r, threads=t, placement=p, cpus=r*t) for r,t,p in LAYOUTS],
        repeats=3, order_seed=SEED, run_timeout_s=RUN_LIMIT_S, total_timeout_s=TOTAL_LIMIT_S,
        expected_host='gnode1', physical_core_budget=16, smt_workers=False,
        pilot_cases=['grid016', 'grid128'], pilot_layouts=[[16,1,'compact'], [1,1,'compact'], [1,16,'compact']],
        profile_layouts=PROFILE_LAYOUTS, profile_interval=20,
        primary_throughput_interval=[20,80], minimum_completed_steps=81,
        step_count_policy='Learn each independent case from its first completed 16x1 control; every layout/repeat must match that case.',
        comparison=dict(same_rank='strict bitwise against candidate t1 control',
                        cross_rank='record existing rtol 1e-9/scaled_atol 1e-12, do not gate or relax'),
        limits=['Every layout within a case uses byte-identical physics and output cadence.',
                'Different cases intentionally change work and may change boundary/particle-discretization physics.',
                'The Gaussian is narrow in z; load imbalance and serial deposition/ledger are confounders.',
                'Particle preparation uses 2048-particle batches and stays serial below 32 particles; a small-workload crossover is not by itself evidence of cache saturation.',
                'Compact lower-core layouts use fewer NUMA domains; 8x1-spread is a placement sensitivity check.',
                'RSS is not cache working set. No hardware cache-miss/bandwidth counters are asserted.',
                'Only profiling runs sample /proc numa_maps; throughput has no process-sampling overhead.',
                'No old executable, production input, GPU, network setting, or system setting is changed.'])
    save(root / 'plan.json', plan)
    frozen = [root/'plan.json', *sorted(cards.iterdir()), *sorted(scripts.iterdir()), binary, probe]
    save(root/'frozen.json', {str(f): digest(f) for f in frozen})
    print(json.dumps(dict(prepared=str(root), pilot_command=f'bash {scripts / "run.sh"} --root {root} --stage pilot',
        matrix_command=f'bash {scripts / "run.sh"} --root {root} --stage matrix',
        all_command=f'bash {scripts / "run.sh"} --root {root} --stage all'), indent=2))


def verify_frozen(root):
    for filename, expected in json.loads((root/'frozen.json').read_text()).items():
        if digest(filename) != expected: raise ValueError('Frozen artifact changed: ' + filename)


def topology():
    result = subprocess.run(['lscpu', '--parse=CPU,CORE,SOCKET,NODE,ONLINE'], text=True, capture_output=True, check=True)
    rows = [line.split(',') for line in result.stdout.splitlines() if line and not line.startswith('#')]
    rows = [r for r in rows if r[4].lower() in ('y','yes','1')]
    cores = {(int(r[2]),int(r[1])) for r in rows}
    if len(cores) != 16: raise ValueError('Host does not match reviewed 16-physical-core topology')
    return result.stdout, {int(r[0]): dict(core=int(r[1]),socket=int(r[2]),numa=int(r[3])) for r in rows}


def numa_snapshot(config_path):
    rows = []
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit(): continue
        try:
            argv = (directory/'cmdline').read_bytes().split(b'\0')
            if str(config_path).encode() not in argv or not argv[0].endswith(b'/aprl'): continue
            status = (directory/'status').read_text()
            maps = (directory/'numa_maps').read_text()
            nodes = {}
            for node, pages in re.findall(r'\bN(\d+)=(\d+)', maps):
                nodes[node] = nodes.get(node,0) + int(pages)
            env = (directory/'environ').read_bytes().split(b'\0')
            rank = next((x.split(b'=',1)[1].decode() for x in env if x.startswith(b'OMPI_COMM_WORLD_RANK=')), None)
            affinity = next((x.split(':',1)[1].strip() for x in status.splitlines() if x.startswith('Cpus_allowed_list:')), None)
            rows.append(dict(pid=int(directory.name), rank=rank, mapped_pages_by_numa_node=nodes,
                             cpus_allowed_list=affinity, maps=maps))
        except (OSError, ValueError): continue
    return dict(captured=bool(rows), processes=rows,
        limitation='One /proc mapping snapshot; includes shared/runtime pages, not measured bandwidth or remote-access counts.')


class Campaign:
    def __init__(self, root):
        self.root = root
        self.plan = json.loads((root/'plan.json').read_text())
        self.generic = module('memory_generic', root/'scripts/generic_inspector.py')
        self.inspector = module('memory_basic', root/'scripts/basic_inspector.py')
        report_file = root/'report.json'
        self.report = json.loads(report_file.read_text()) if report_file.exists() else dict(
            format='aprl-memory-scaling-result-v1', status='prepared', pass_=False,
            started_epoch_s=time.time(), deadline_epoch_s=time.time()+TOTAL_LIMIT_S,
            host=platform.node(), runs=[], comparisons=[], references={}, binding=[],
            case_steps={}, notes=self.plan['limits'], measured_hardware_counters=False)
        self.cases = {c['name']: c for c in self.plan['cases']}
        self.case_info = {name:self.generic.read_case(Path(c['card']))[1] for name,c in self.cases.items()}

    def checkpoint(self):
        self.report['updated_at'] = datetime.now(timezone.utc).isoformat()
        save(self.root/'report.json', self.report)

    def remaining(self):
        remaining = self.report['deadline_epoch_s'] - time.time()
        if remaining < 10: raise ValueError('Two-hour campaign limit reached; no new launch')
        return remaining

    def execute(self, command, directory, threads, profile, filename='run.log'):
        limit = min(RUN_LIMIT_S, max(1,int(self.remaining()-5)))
        command = ['timeout', '--signal=TERM', '--kill-after=10', str(limit), *command]
        env = environment(threads, profile)
        save(directory/'launch.json', dict(command=command, cwd=str(directory),
            timeout_s=limit, environment={key:env.get(key) for key in
            ('OMP_NUM_THREADS','OMP_PLACES','OMP_PROC_BIND','APRL_PROFILE_INTERVAL',
             'OMPI_MCA_btl','OMPI_MCA_pml','UCX_TLS','LD_LIBRARY_PATH')}))
        start = time.monotonic(); snapshot = None
        with (directory/filename).open('w') as output:
            process = subprocess.Popen(command, cwd=directory, env=env, stdout=output, stderr=subprocess.STDOUT)
            if profile:
                # Optional single sample only in profile trials, once stepping starts.
                while process.poll() is None:
                    text = (directory/filename).read_text(errors='replace')
                    if 'phase=progress step=20 ' in text:
                        snapshot = numa_snapshot(directory/'run.yaml')
                        break
                    time.sleep(.1)
            code = process.wait()
        if profile: save(directory/'numa-snapshot.json', snapshot or dict(captured=False, reason='No live step-20 snapshot available'))
        return code, time.monotonic()-start

    def binding(self):
        if (self.root/'BINDING_PASS').exists(): return
        text, cpu = topology(); (self.root/'cpu-topology.csv').write_text(text)
        self.report['topology'] = cpu
        self.report['environment'] = dict(loadavg_start=os.getloadavg(), python=sys.executable,
            perf_event_paranoid=Path('/proc/sys/kernel/perf_event_paranoid').read_text().strip(),
            mpirun=subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0])
        for layout in LAYOUTS:
            directory = self.root/('binding-'+layout_name(layout)); directory.mkdir()
            code, wall = self.execute(mpi_command(layout,self.plan['binding_probe'],layout[:2]),directory,layout[1],0)
            log=(directory/'run.log').read_text(errors='replace')
            if code or '[binding] result=PASS failures=0' not in log:
                raise ValueError('Binding failed: '+layout_name(layout))
            matches=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=gnode1 cpu=(\d+) package=(\d+) core=(\d+)',log)
            if len(matches)!=layout[0]*layout[1]: raise ValueError('Missing binding workers')
            workers=[dict(rank=int(r),thread=int(t),cpu=int(c),socket=int(s),core=int(k),numa=cpu[int(c)]['numa']) for r,t,c,s,k in matches]
            nodes=sorted({x['numa'] for x in workers})
            if layout[2]=='spread' and len(nodes)!=8: raise ValueError('Spread layout did not cover all 8 NUMA domains')
            record=dict(layout=layout_name(layout),workers=workers,numa_domains=nodes,wall_s=wall)
            self.report['binding'].append(record); save(directory/'binding.json',record); self.checkpoint()
        (self.root/'BINDING_PASS').write_text('All requested layouts use distinct physical cores.\n')

    def reference_arrays(self, name):
        run=next(r for r in self.report['runs'] if r['name']==name)
        return self.generic.inspect(self.inspector,run,self.case_info[run['case']])[1]

    def compare(self, run, arrays, reference_name, strict):
        reference=self.reference_arrays(reference_name)
        result=self.generic.compare(self.inspector,reference,arrays,strict)
        suffix='same-rank' if strict else 'cross-rank'
        file=Path(run['directory'])/(suffix+'-comparison.json'); save(file,result)
        self.report['comparisons'].append(dict(run=run['name'],reference=reference_name,
            mode=suffix,pass_=result['pass'],details=str(file),gating=strict))
        if strict and not result['pass']: raise ValueError(run['name']+': same-rank bitwise comparison failed')

    def one(self, case_name, layout, phase, repeat=0):
        name=f'{phase}-{case_name}-{layout_name(layout)}-rep{repeat}'
        if any(r['name']==name for r in self.report['runs']): raise ValueError('Refusing duplicate run '+name)
        directory=self.root/name; directory.mkdir()
        case=self.cases[case_name]; shutil.copy2(case['card'],directory/'run.yaml')
        ranks,threads,placement=layout; interval=20 if phase=='profile' else 0
        row=dict(name=name,case=case_name,directory=str(directory),label='candidate',
            ranks=ranks,threads=threads,placement=placement,phase=phase,repeat=repeat,
            cpus=ranks*threads,profile_interval=interval,valid=False,status='running',
            configured_grid_cells=case['grid_cells'],configured_particles=case['particles'],
            slab_cells_min=case['cells'][2]//ranks,slab_cells_max=math.ceil(case['cells'][2]/ranks),
            grid_cells_per_rank=case['grid_cells']/ranks,
            average_particles_per_rank_not_balance_measure=case['particles']/ranks,
            loadavg_before=os.getloadavg())
        self.report['runs'].append(row); self.checkpoint(); print('Running',name,flush=True)
        code,wall=self.execute(mpi_command(layout,self.plan['binary'],[directory/'run.yaml']),directory,threads,interval)
        row.update(returncode=code,process_wall_s=wall,loadavg_after=os.getloadavg()); self.checkpoint()
        if code: raise ValueError(name+f': process exit {code}; see run.log')
        if digest(directory/'run.yaml')!=case['sha256']: raise ValueError('Run card changed')
        manifest=yaml.safe_load((directory/'run-manifest.yaml').read_text())
        if manifest['source_revision']!=CORE or manifest['mpi_size']!=ranks: raise ValueError('Compute identity mismatch')
        if manifest['resolved_mesh']['cells']!=case['cells']: raise ValueError('Actual mesh differs from case')
        summary,arrays=self.generic.inspect(self.inspector,row,self.case_info[case_name])
        row.update(summary)
        progress={}
        for line in (directory/'run.log').read_text(errors='replace').splitlines():
            if not line.startswith('[resource] phase=progress '): continue
            fields=dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s,]+)',line))
            step=int(fields['step'])
            if step in progress: raise ValueError('Duplicate progress step')
            progress[step]=float(fields['elapsed_s'])
        counts=self.report.setdefault('case_steps',{})
        if case_name not in counts:
            if (ranks,threads,placement)!=(16,1,'compact'):
                raise ValueError('First physical-stop reference for each case must be compact 16x1')
            if summary['steps']<81: raise ValueError('Physical endpoint does not include the common 20-80 interval')
            counts[case_name]=summary['steps']
        if summary['steps']!=counts[case_name] or not {20,80}.issubset(progress):
            raise ValueError('Case endpoint or common 20/80 progress markers differ')
        steady=progress[80]-progress[20]
        if not math.isfinite(steady) or steady<=0: raise ValueError('Invalid steady timing')
        row.update(steady_steps=60,steady_wall_s=steady,steady_seconds_per_step=steady/60,
            steady_grid_cell_updates_per_s=case['grid_cells']*60/steady,
            steady_particle_updates_per_s=case['particles']*60/steady,
            steady_core_seconds=steady*ranks*threads,
            peak_rss_rank_max_mib=float(summary['final_resource'][-1]['peak_rss_rank_max_mib']),
            peak_rss_total_mib=float(summary['final_resource'][-1]['peak_rss_total_mib']))
        if phase=='profile':
            stage_rows=[p for p in summary['performance'] if p.get('stage') in STAGES and 20<int(p['step'])<=80]
            workload=[p for p in summary['performance'] if p.get('phase')=='workload' and 20<int(p['step'])<=80]
            if {p['stage'] for p in stage_rows}!=set(STAGES) or not workload:
                raise ValueError('Profiling did not record every stage and particle-load statistics')
            typed=[]
            for stage in STAGES:
                data=[p for p in stage_rows if p['stage']==stage]
                nsteps=sum(int(p['steps']) for p in data)
                if nsteps!=60: raise ValueError('Profile stage does not cover exactly steps 20-80')
                typed.append(dict(stage=stage,steps=nsteps,
                    rank_mean_seconds_per_step=sum(float(p['mean_s']) for p in data)/nsteps,
                    rank_max_seconds_per_step=sum(float(p['max_s']) for p in data)/nsteps,
                    maximum_rank_ids=sorted({int(p['max_rank']) for p in data})))
            load=[dict(step=int(p['step']),particles_total=int(p['particles_total']),
                particles_min=int(p['particles_min']),particles_max=int(p['particles_max']),
                particles_mean=int(p['particles_total'])/ranks,
                active_total=int(p['active_total']),frozen_total=int(p['frozen_total'])) for p in workload]
            if any(p['particles_total']!=case['particles'] for p in load):
                raise ValueError('Profile workload particle count differs from case')
            row['profile_summary']=dict(stages=typed,workload=load,
                particle_max_to_mean_max=max(p['particles_max']/p['particles_mean'] for p in load),
                any_rank_empty=any(p['particles_min']==0 for p in load),
                interpretation='Communication stages are nested and include MPI waits. Stage maxima are not additive elapsed time. RSS is not cache working set.')
            save(directory/'profile-summary.json',row['profile_summary'])
        refs=self.report['references'].setdefault(case_name,{})
        same=refs.get(str(ranks))
        if same: self.compare(row,arrays,same,True)
        elif threads!=1: raise ValueError('Missing same-rank single-thread control')
        else: refs[str(ranks)]=name
        if ranks!=16 and '16' in refs: self.compare(row,arrays,refs['16'],False)
        row.update(valid=True,status='passed'); self.checkpoint()

    def pilot(self):
        if (self.root/'PILOT_PASS').exists(): return
        self.binding()
        for case in self.plan['pilot_cases']:
            for layout in self.plan['pilot_layouts']: self.one(case,tuple(layout),'pilot')
        self.report['pilot_pass']=True; self.report['status']='pilot_passed'; self.checkpoint()
        (self.root/'PILOT_PASS').write_text('Candidate-only smallest/largest pilot validated; no matrix run implied.\n')

    def matrix(self):
        if (self.root/'PASS').exists(): return
        if not (self.root/'PILOT_PASS').exists(): raise ValueError('Pilot must pass before the full matrix')
        for case in self.cases:
            for layout in LAYOUTS: self.one(case,layout,'warmup')
        configurations=[(case,layout) for case in self.cases for layout in LAYOUTS]
        for repeat in range(1,4):
            order=configurations.copy(); random.Random(SEED+repeat).shuffle(order)
            for case,layout in order: self.one(case,layout,'throughput',repeat)
        for case in self.cases:
            for layout in PROFILE_LAYOUTS: self.one(case,layout,'profile')
        self.summarize(); self.report.update(pass_=True,status='passed'); self.checkpoint()
        (self.root/'PASS').write_text('All throughput/profile runs complete and same-rank bitwise checks passed; cross-rank comparisons are reported separately.\n')

    def summarize(self):
        groups={}
        for row in self.report['runs']:
            if row.get('valid') and row['phase']=='throughput':
                groups.setdefault((row['case'],row['ranks'],row['threads'],row['placement']),[]).append(row)
        metrics=('steady_seconds_per_step','steady_core_seconds','loop_wall_s','reported_wall_s','process_wall_s',
                 'peak_rss_rank_max_mib','peak_rss_total_mib','steady_grid_cell_updates_per_s','steady_particle_updates_per_s')
        rows=[]
        for (case,ranks,threads,placement),selection in sorted(groups.items()):
            record=dict(case=case,ranks=ranks,threads=threads,placement=placement,cpus=ranks*threads,repeats=len(selection),
                grid_cells=self.cases[case]['grid_cells'],particles=self.cases[case]['particles'])
            for metric in metrics:
                values=[r[metric] for r in selection]
                for stat,fn in (('median',statistics.median),('min',min),('max',max)):
                    record[metric+'_'+stat]=fn(values)
            rows.append(record)
        self.report['timing_summary']=rows
        cross=[c for c in self.report['comparisons'] if not c['gating']]
        self.report['cross_rank_consistency']=dict(all_within_existing_tolerance=all(c['pass_'] for c in cross),
            compared=len(cross),failed=sum(not c['pass_'] for c in cross),gating=False)
        if rows:
            with (self.root/'summary.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        profiles=[]
        for row in self.report['runs']:
            if row.get('valid') and row['phase']=='profile':
                for stage in row['profile_summary']['stages']:
                    profiles.append(dict(case=row['case'],ranks=row['ranks'],threads=row['threads'],placement=row['placement'],
                        stage=stage['stage'],rank_mean_seconds_per_step=stage['rank_mean_seconds_per_step'],
                        rank_max_seconds_per_step=stage['rank_max_seconds_per_step'],
                        particle_max_to_mean_max=row['profile_summary']['particle_max_to_mean_max'],
                        any_rank_empty=row['profile_summary']['any_rank_empty']))
        if profiles:
            with (self.root/'profile-summary.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(profiles[0])); writer.writeheader(); writer.writerows(profiles)


def self_test(workspace):
    """Metadata and validation checks only; never starts MPI or APRL."""
    import tempfile
    with tempfile.TemporaryDirectory(prefix='aprl-memory-scaling-selftest-') as temporary:
        root=Path(temporary)/'prepared'
        prepare(workspace,root); verify_frozen(root); campaign=Campaign(root)
        if len(campaign.cases)!=8 or len(LAYOUTS)!=10: raise AssertionError('Matrix dimensions changed')
        original=next(iter(campaign.case_info.values()))['physical_card']
        for name,info in campaign.case_info.items():
            card=copy.deepcopy(info['physical_card'])
            reference=copy.deepcopy(original)
            for item in (card,reference):
                item['mesh']['cells']=[0,0,0]
                item['beam']['input']['macroparticles']=0
            if card!=reference: raise AssertionError('Non-workload physics changed between cases')
            if info['represented_electrons']!=100000: raise AssertionError('Charge is not fixed')
        for layout in LAYOUTS:
            if layout[0]*layout[1]>16: raise AssertionError('CPU budget exceeded')
        if 'node:PE=16' not in mpi_command((1,16,'compact'),'probe',[]):
            raise AssertionError('Full-node thread team mapping lost')
        if 'ppr:1:numa' not in mpi_command((8,1,'spread'),'probe',[]):
            raise AssertionError('NUMA-sensitive mapping lost')
        sequences=[]
        for rep in range(1,4):
            sequence=[(case,layout) for case in campaign.cases for layout in LAYOUTS]
            random.Random(SEED+rep).shuffle(sequence); sequences.append(sequence)
            if len(sequence)!=80 or len(set(sequence))!=80: raise AssertionError('Matrix rows duplicated')
        if sequences[0]==sequences[1]: raise AssertionError('Repeats not reshuffled')
        frozen_file=next((root/'cards').iterdir()); frozen_file.write_text(frozen_file.read_text()+'\n')
        try: verify_frozen(root)
        except ValueError: pass
        else: raise AssertionError('Changed frozen card accepted')
        print('Self-test PASS: case isolation, fixed charge, layouts, repeat ordering, frozen-card rejection; no compute launched.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,default=WORKSPACE)
    parser.add_argument('--root',type=Path)
    parser.add_argument('--stage',choices=('prepare','pilot','matrix','all'))
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test: self_test(args.workspace.resolve()); return 0
    if args.root is None or args.stage is None: parser.error('--root and --stage are required')
    root=args.root.resolve()
    if args.stage=='prepare': prepare(args.workspace.resolve(),root); return 0
    if platform.node()!='gnode1' or os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS')!='1' or os.environ.get('SLURM_JOB_ID'):
        parser.error('Execution requires explicitly authorized gnode1 CPU-only environment without Slurm variables')
    verify_frozen(root)
    campaign=Campaign(root)
    if (root/'PASS').exists() or (args.stage=='pilot' and (root/'PILOT_PASS').exists()):
        print('Requested stage already completed; no duplicate launches.',flush=True); return 0
    if campaign.report['status']=='failed': raise ValueError('Failed campaign cannot be silently resumed; use a new reviewed directory')
    if not campaign.report['runs'] and os.getloadavg()[0]>2: raise ValueError('gnode1 is busy at campaign start')
    try:
        campaign.report['status']='running'; campaign.checkpoint()
        if args.stage in ('pilot','all'): campaign.pilot()
        if args.stage in ('matrix','all'): campaign.matrix()
        verify_frozen(root)
        return 0
    except Exception as error:
        campaign.report.update(status='failed',pass_=False,error=str(error)); campaign.summarize(); campaign.checkpoint()
        (root/'FAILED').write_text(str(error)+'\n'); print('FAILED:',error,flush=True); return 1
    finally:
        campaign.checkpoint()


if __name__=='__main__':
    raise SystemExit(main())
