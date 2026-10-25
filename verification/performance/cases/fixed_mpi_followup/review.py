#!/usr/bin/env python3
"""Independently review a completed fixed_mpi_followup campaign; no launches."""
from pathlib import Path
import argparse, csv, hashlib, json, math, platform, re, statistics
from collections import Counter
from datetime import datetime, timezone
import h5py, numpy as np, yaml
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,required=True,help='Completed campaign directory')
parser.add_argument('--superseded-root',type=Path,help='Optional prior directory explicitly marked NOT_LAUNCHED')
args=parser.parse_args()
ROOT=args.root.resolve()
P=json.loads((ROOT/'plan.json').read_text()); R=json.loads((ROOT/'report.json').read_text())
def check(ok,msg):
    if not ok: raise AssertionError(msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def equal(a,b):return math.isclose(float(a),float(b),rel_tol=1e-13,abs_tol=1e-15)
def fields(line):return dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s,]+)',line))
frozen=json.loads((ROOT/'frozen.json').read_text())
for path,want in frozen.items():check(sha(path)==want,'Changed frozen file: '+path)
check(R['status']=='passed' and R['pass_'] and (ROOT/'PASS').is_file(),'Incomplete campaign')
check(len(R['runs'])==48 and len(R['comparisons'])==45,'Wrong counts')
check(Counter(r['phase'] for r in R['runs'])=={'pilot':6,'warmup':6,'throughput':30,'profile':6},'Phase coverage')
check(len({r['name'] for r in R['runs']})==48,'Duplicate runs')
check(P['layout']==dict(ranks=4,threads=4,placement='compact',cpus=16) and R['host']=='gnode1','Wrong layout/host')
check(P['binaries']['control']['build']['selected']==P['binaries']['trial']['build']['selected'],'Different toolchains')
for name,b in P['binaries'].items():
    flags=Path(b['path']).parent/'CMakeFiles/direct_eb_core.dir/flags.make'
    rev=re.findall(r'-DAPRL_SOURCE_REVISION=\\"([0-9a-f]{40})\\"',flags.read_text())
    check(rev==[b['compiled_revision']] and sha(b['path'])==b['sha256'],'Binary/flags mismatch: '+name)
workers=json.loads((ROOT/'binding-r4-t4-compact/binding.json').read_text())['workers']
check(len(workers)==16 and len({(w['socket'],w['core']) for w in workers})==16,'Physical binding failed')
check(len({w['numa'] for w in workers})==8 and all(sum(w['rank']==rank for w in workers)==4 for rank in range(4)),'Binding topology')
C={c['name']:c for c in P['cases']}
PRODUCTS=(('ledger/energy-ledger.h5','energy_ledger','records','committed_records'),('detectors/particles.h5','particle_plane','records','committed_records'),('detectors/field.h5','field_plane','fields','committed_samples'))
def numeric(a,prefix):
    if a.dtype.names:
        for name in a.dtype.names:yield from numeric(a[name],prefix+'.'+name)
    elif np.issubdtype(a.dtype,np.number):yield prefix,a
def science(directory,case):
    arrays={};count=C[case]['particles']
    for filename,gn,dn,cn in PRODUCTS:
        with h5py.File(directory/filename) as f:
            g=f[gn];a=g[dn][:]
            check(int(g['complete'][()])==1 and len(a)>0 and int(g[cn][()])==len(a),'Incomplete '+str(directory/filename))
            if gn=='particle_plane':
                a=a[np.argsort(a['source_id'],kind='stable')]
                check(np.array_equal(a['source_id'],np.arange(1,count+1,dtype=a['source_id'].dtype)) and int(g.attrs['mpi_size'])==4,'Particle cohort/MPI')
            arrays.update(numeric(a,gn))
            if gn=='field_plane':
                t=g['time_s'][:];check(len(t)==len(a) and np.all(np.diff(t)>0),'Field times')
                arrays[gn+'.time_s']=t
    for name,a in arrays.items():check(np.all(np.isfinite(a)),'Nonfinite '+name)
    for name,want in [('active_macroparticles',count),('removed_macroparticles',0),('closure_valid',1),('boundary_injection_entries',0)]:
        check(np.all(arrays['energy_ledger.'+name]==want),'Ledger '+name)
    s=arrays['energy_ledger.step'];check(s[0]==0 and s[-1]>80 and np.all(np.diff(s.astype(np.int64))>=0),'Ledger step sequence')
    return arrays
refs={case:science(ROOT/name,case) for case,name in R['references'].items()}
steps={case:int(a['energy_ledger.step'][-1]) for case,a in refs.items()}
check(steps==R['case_steps'],'Case endpoint discrepancy')
comp={c['run']:c for c in R['comparisons']};raw=[];prof=[];scientific=[];evidence=[]
for run in R['runs']:
    d=Path(run['directory']);case=run['case'];b=P['binaries'][run['label']]
    check(run['valid'] and run['status']=='passed' and run['returncode']==0,'Bad run '+run['name'])
    check((run['ranks'],run['threads'],run['cpus'])==(4,4,16),'Run layout')
    check(sha(d/'run.yaml')==C[case]['sha256']==sha(C[case]['card'])==sha(Path(P['source_campaign'])/'cards'/Path(C[case]['card']).name),'Cards differ')
    m=yaml.safe_load((d/'run-manifest.yaml').read_text())
    check(m['source_revision']==b['compiled_revision']==run['compiled_revision'] and run['binary_sha256']==b['sha256'],'Manifest revision')
    check(m['mpi_size']==4 and m['resolved_mesh']['cells']==C[case]['cells'],'Manifest geometry/layout')
    launch=json.loads((d/'launch.json').read_text());cmd=launch['command'];env=launch['environment']
    check(cmd[cmd.index('-np')+1]=='4' and cmd[cmd.index('--map-by')+1]=='slot:PE=4' and b['path'] in cmd,'Actual launch')
    check(env['OMP_NUM_THREADS']=='4' and env['OMP_PLACES']=='cores' and env['OMP_PROC_BIND']=='close','OpenMP environment')
    check(env['LD_LIBRARY_PATH'].split(':')[0]==P['runtime_libraries']['directory'],'Runtime library path')
    check(env['APRL_PROFILE_INTERVAL']==('20' if run['phase']=='profile' else '0'),'Profile environment')
    log=(d/'run.log').read_text();lines=log.splitlines()
    check('Configured stop reached:' in log and 'duration was exhausted' not in log and 'stopped after a complete field step' not in log,'Abnormal physical stop')
    check(re.search(r'Compute threads per MPI rank:\s*4;',log),'Reported thread count')
    q=[fields(l) for l in lines if l.startswith('[resource] phase=progress ')]
    check(len(q)==len({x['step'] for x in q}),'Duplicate progress markers')
    q={int(x['step']):float(x['elapsed_s']) for x in q};check(20 in q and 80 in q,'Timing markers')
    seconds=(q[80]-q[20])/60;check(seconds>0 and equal(seconds,run['steady_seconds_per_step']),'Raw/report timing mismatch')
    final=[fields(l) for l in lines if l.startswith('[resource] phase=final ')]
    check(len(final)==1 and int(final[0]['step'])==steps[case]==run['steps'],'Final step mismatch')
    total=float(final[0]['peak_rss_total_mib']);rank=float(final[0]['peak_rss_rank_max_mib'])
    check(equal(total,run['peak_rss_total_mib']) and equal(rank,run['peak_rss_rank_max_mib']),'RSS mismatch')
    if run['phase']=='throughput':
        check(not any(l.startswith('[performance]') for l in lines) and not (d/'numa-snapshot.json').exists(),'Profile overhead in throughput')
        raw.append(dict(name=run['name'],case=case,variant=run['label'],repeat=run['repeat'],step20_s=q[20],step80_s=q[80],seconds_per_step=seconds,peak_rss_total_mib=total,peak_rss_rank_max_mib=rank,run_log_sha256=sha(d/'run.log')))
    if run['phase']=='profile':
        ps=[fields(l) for l in lines if l.startswith('[performance]')]
        ps=[p for p in ps if 'stage' in p and 20<int(p['step'])<=80]
        for stage in sorted({p['stage'] for p in ps}):
            selected=[p for p in ps if p['stage']==stage]
            check([int(p['step']) for p in selected]==[40,60,80] and sum(int(p['steps']) for p in selected)==60,'Profile coverage')
            prof.append(dict(case=case,variant=run['label'],stage=stage,rank_mean_seconds_per_step=sum(float(p['mean_s']) for p in selected)/60,rank_max_seconds_per_step=sum(float(p['max_s']) for p in selected)/60))
    arrays=science(d,case);check(int(arrays['energy_ledger.step'][-1])==steps[case],'H5 endpoint')
    if run['name'] in comp:
        item=comp[run['name']];saved=json.loads(Path(item['details']).read_text())
        check(item['gating'] and item['pass_'] and item['mode']=='same-rank' and item['reference']==R['references'][case],'Comparison provenance')
        reference=refs[case];check(set(arrays)==set(reference)==set(saved['fields']) and saved['pass'],'Comparison schema')
        for key,a in arrays.items():
            z=reference[key];check(a.dtype==z.dtype and a.shape==z.shape and np.ascontiguousarray(a).tobytes()==np.ascontiguousarray(z).tobytes(),'Independent bitwise mismatch: '+run['name']+' '+key)
            c=saved['fields'][key];check(c['pass'] and (c.get('exact',False) or c.get('bitwise_equal',False)),'Stored comparison differs')
        scientific.append(dict(run=run['name'],reference=item['reference'],numeric_fields=len(arrays),bitwise_equal=True))
    evidence.append(dict(run=run['name'],manifest_sha256=sha(d/'run-manifest.yaml'),card_sha256=C[case]['sha256'],normal_stop=True,steps=steps[case],products_complete=3))
summary={(r['case'],r['variant']):r for r in csv.DictReader((ROOT/'summary.csv').open())}
comparison={r['case']:r for r in csv.DictReader((ROOT/'comparison.csv').open())}
profile={(r['case'],r['variant'],r['stage']):r for r in csv.DictReader((ROOT/'profile-summary.csv').open())}
check(len(prof)==len(profile)==48,'Profile total coverage')
for p in prof:
    for key in ('rank_mean_seconds_per_step','rank_max_seconds_per_step'):check(equal(p[key],profile[(p['case'],p['variant'],p['stage'])][key]),'Raw profile differs from CSV')
results=[]
for case in C:
    groups={v:sorted([r for r in raw if r['case']==case and r['variant']==v],key=lambda r:r['repeat']) for v in ('control','trial')}
    for v,values in groups.items():
        check([r['repeat'] for r in values]==[1,2,3,4,5],'Repeat coverage')
        for metric,key in [('steady_seconds_per_step','seconds_per_step'),('peak_rss_total_mib','peak_rss_total_mib'),('peak_rss_rank_max_mib','peak_rss_rank_max_mib')]:
            for label,fn in [('median',statistics.median),('min',min),('max',max)]:check(equal(fn(r[key] for r in values),summary[(case,v)][metric+'_'+label]),'Raw summary differs from CSV')
    cs=[r['seconds_per_step'] for r in groups['control']];ts=[r['seconds_per_step'] for r in groups['trial']]
    ratios=[c/t for c,t in zip(cs,ts)];c=statistics.median(cs);t=statistics.median(ts)
    result=dict(case=case,paired_repeats=5,control_seconds_per_step_median=c,trial_seconds_per_step_median=t,control_over_trial_medians=c/t,trial_time_reduction_fraction=1-t/c,paired_speedup_median=statistics.median(ratios),paired_speedup_min=min(ratios),paired_speedup_max=max(ratios),trial_faster_pairs=sum(t<c for c,t in zip(cs,ts)),control_seconds_per_step_range=[min(cs),max(cs)],trial_seconds_per_step_range=[min(ts),max(ts)],paired_speedups=ratios)
    for key in comparison[case]:
        if key!='case':check(equal(result[key],comparison[case][key]),'Raw comparison differs from CSV')
    result['peak_rss_total_mib']={v:dict(median=statistics.median(r['peak_rss_total_mib'] for r in rows),minimum=min(r['peak_rss_total_mib'] for r in rows),maximum=max(r['peak_rss_total_mib'] for r in rows)) for v,rows in groups.items()}
    result['peak_rss_rank_max_mib_median']={v:statistics.median(r['peak_rss_rank_max_mib'] for r in rows) for v,rows in groups.items()}
    results.append(result)
superseded=None
if args.superseded_root is not None:
    old=args.superseded_root.resolve();note=json.loads((old/'NOT_LAUNCHED.json').read_text())
    check(note['status']=='superseded-before-launch' and (old/note['replacement']).resolve()==ROOT,'Old campaign identity')
    check(not (old/'report.json').exists() and not list(old.glob('*/run.log')),'Old campaign launched unexpectedly')
    superseded=dict(path=str(old),status=note['status'],no_report_or_run_logs=True,reason=note['reason'],effect_on_campaign='None; all reviewed manifests and actual compile flags match the correct full trial revision.')
review=dict(format='aprl-independent-fixed-mpi-followup-review-v1',status='PASS',reviewed_at=datetime.now(timezone.utc).isoformat(),review_host=platform.node(),campaign=str(ROOT),compute_launched=False,method='Independent raw timing/profile parser and HDF5 numeric-field byte comparison; particles sorted by source_id; all frozen hashes rechecked.',review_script=str(Path(__file__).resolve()),review_script_sha256=sha(__file__),frozen_files_checked=len(frozen),run_count=48,throughput_runs_recomputed=30,profile_runs_recomputed=6,profile_rows_recomputed=48,bitwise_comparisons_recomputed=45,numeric_fields_per_comparison=sorted({r['numeric_fields'] for r in scientific}),h5_files_checked=144,case_final_steps=steps,required_tests_passed=72,binaries=P['binaries'],binding=dict(physical_cores=16,ranks=4,threads=4,numa_domains=8,pass_=True),throughput=results,raw_throughput=raw,profile=prof,scientific_comparisons=scientific,runs=evidence,superseded_campaign=superseded,limitations=['Short fixed-4x4 gnode1 generated-gaussian workloads, not a full-scale two-node production guarantee.','The existing numeric-dataset comparison was independently recomputed; no new all-attributes equivalence claim.','Five paired repeats report observed variability; small differences or mixed signs are not stable speedup evidence.','Single separate profile per variant/case; MPI stages include waits, and nested stages or different-rank maxima cannot be added.','No precision, mesh, physical endpoint, cadence, production executable or frozen campaign modification.'])
text=[ROOT.name+'：独立验收 PASS','', '固定gnode1、4 MPI × 4 OpenMP。48次运行完整；30份原始throughput日志按(elapsed80-elapsed20)/60复算，与summary/comparison一致。', '72 required证据冻结保持；同卡各版均到正常物理终点：'+', '.join(f'{case} {count}步' for case,count in steps.items())+'。', f'144个H5重新检查完整提交、有限值、粒子守恒；45组共{sum(r["numeric_fields"] for r in scientific)}个数值字段独立逐位重比全部通过。', f'全部{len(frozen)}个冻结文件SHA保持。三张原卡逐字节相同；二进制、实际compile flags及48份manifest版本、4×4绑定一致。','', '关闭profiling、每卡5次配对：']
for r in results:
    text.append(f"- {r['case']}: {r['control_seconds_per_step_median']*1000:.6f} -> {r['trial_seconds_per_step_median']*1000:.6f} ms/步；中位耗时降低{r['trial_time_reduction_fraction']*100:.4f}%；配对加速中位{r['paired_speedup_median']:.6f}，范围{r['paired_speedup_min']:.6f}–{r['paired_speedup_max']:.6f}；{r['trial_faster_pairs']}/5对更快。")
text+=['', 'throughput总峰值RSS中位（control -> trial，MiB）：'+', '.join(f"{r['case']} {r['peak_rss_total_mib']['control']['median']:.3f} -> {r['peak_rss_total_mib']['trial']['median']:.3f}" for r in results)+'。', '', '独立profile原始日志复算一致。以下为rank-mean毫秒/步（control -> trial）：']
by_profile={(p['case'],p['variant'],p['stage']):p for p in prof}
for case in C:
    descriptions=[]
    for stage in ('field_total','particle_deposit','diagnostics','field_communication'):
        control=by_profile[(case,'control',stage)]['rank_mean_seconds_per_step']*1000
        trial=by_profile[(case,'trial',stage)]['rank_mean_seconds_per_step']*1000
        descriptions.append(f'{stage} {control:.6f} -> {trial:.6f}')
    text.append('- '+case+': '+', '.join(descriptions))
text+=['profile每版每卡仅一次；不与关闭profile计时混合、不相加rank maxima或嵌套阶段。','','程序身份：']
for variant,binary in P['binaries'].items():
    text.append(f"- {variant}: {binary['path']}; compiled {binary['compiled_revision']}; SHA256 {binary['sha256']}")
text.append('两版工具链和只读FFTW依赖一致。')
if superseded is not None:
    text+=['', '旧目录'+superseded['path']+'在启动前废止，NOT_LAUNCHED保留且无report/run.log。当前目录实际编译flags和全部manifest正确，旧元数据纠正不影响本次结果。']
text+=['', '边界：仅当前短卡同4×4代码等价和性能证据，不代表完整两节点大工况、跨rank一致性或辐射物理收敛。微小差异及配对正负混合不作为稳定收益。此次仅本地读文件并写复核结果，没有访问测试机或启动新计算。']
(ROOT/'REVIEWED_COMPLETION.json').write_text(json.dumps(review,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
(ROOT/'RESULT_REVIEW.txt').write_text('\n'.join(text)+'\n')
print(json.dumps(dict(status='PASS',frozen_files=len(frozen),fields=sum(r['numeric_fields'] for r in scientific),timing=results),indent=2))
