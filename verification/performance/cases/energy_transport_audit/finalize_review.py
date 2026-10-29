#!/usr/bin/env python3
"""Summarize completed, independently reviewed energy diagnostic acceptance."""
import argparse
import hashlib
import json
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--unit', type=Path, required=True)
a=p.parse_args(); root=a.root.resolve(); out=a.output.resolve(); unit=a.unit.resolve()
load=lambda path:json.loads(path.read_text())
stages={s:load(root/(s.upper()+'_REVIEW.json')) for s in ('science','mpi','performance')}
assert all(s['checks_pass'] for s in stages.values())
assert [len(stages[s]['cases']) for s in ('science','mpi','performance')]==[8,4,24]
assert unit.joinpath('PASS').read_text().strip()=='PASS' and unit.joinpath('EXIT_CODE').read_text().strip()=='0'
assert '100% tests passed, 0 tests failed out of 77' in unit.joinpath('required.log').read_text()
source=load(out/'SOURCE_BUILD_R2.json'); binary=Path(source['binary'])
assert hashlib.sha256(binary.read_bytes()).hexdigest()==source['binary_sha256']
assert all(s['binary_sha256']['new']==source['binary_sha256'] for s in stages.values())
science=stages['science']; science_rows={}; transport={}; perf=stages['performance']['performance']
for case in ('undulator','free-drift'):
    e=science['cases'][case+'-r1-new-on']['energy']
    science_rows[case]={'imbalance':e['physical_imbalance'],'final':e['final_components'],
        'scope':'The peak particle-work error is in the independent review; final components are copied without recomputing tiny differences from rounded inventories.'}
    # Do not substitute rounded serialized totals for the high+low reader's delta.
    transport[case]=science['cases'][case+'-r1-new-audit-on']['energy_planes']
result={'status':'FUNCTIONAL_ACCEPTANCE_PASS_PHYSICAL_CLOSURE_NOT_ACCEPTED',
 'core_commit':source['source_head'], 'binary':str(binary),'binary_sha256':source['binary_sha256'],
 'compiled_identity_note':source['note'], 'campaign':str(root),'unit':str(unit),
 'validation':{'required_tests':77,'reader_tests':15,'science_cases':8,'science_raw_bitwise_pairs':6,
               'mpi_cases':4,'mpi_raw_bitwise_pairs':2,'timing_cases':24,
               'explicit_modules':'pending statistics; 2/4-rank moving-window and energy-plane modules',
               'independent_field_decimal_checks':88},
 'boost':science_rows, 'lab_transport':transport, 'performance':perf,
 'figures':str(out/'figures-final-r3'),
 'figure_validation':'Three PNGs visually reviewed; PDFs share the export source, not independently rendered.',
 'scope':{'physics_unchanged':True,'pure_near_field_separated':False,'lab_control_volume_closed':False,
          'production_cost_established':False,'production_rerun':False,'default_switched':False,
          'pushed':False,'old_checkpoint_abi_compatible':False},
 'earlier_attempts':{'19011':'Old76 passed; new Boris fixture violated cdt=dz; fixed fixture only.',
 '19012':'Required/modules passed; before cancellation and false-partial corrections.',
 '19013':'Two closed magnetic cases; new work algebra caught false partial and history-cancellation issues; retained FAIL.',
 '19015':'Final core77required/pending/MPI modules/12reader PASS.',
 '19016':'Three complete magnetic cases; audit config had duplicate detector names. Reused only the three complete cases with62 H5 SHA checks.',
 '19017':'Final scientific acceptance PASS; five new cases plus three validated reused cases.',
 '19018':'24 profile-off timing cases PASS.', '19019':'Two-node4 science cases PASS.'}}
(out/'FINAL.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
lines=['两套能量诊断已实现并完成小工况功能验收；物理收支仍有未闭合项。',
 '', '实现与范围',
 'boost v5按实际Boris两次电踢独立记grid/prescribed功，磁旋转及舍入差单列；完整步不做时间截断，真正退场才按子步前缀记账。各功用Kahan长双精度累加，绝对历史功独立记录。',
 '补全退出粒子两系动能和cPz，稳定差量在long double中相减后写入，直接检验dK_lab=Gamma*(dK_box+beta*d(cPz_box))。保留旧外功和普通退场口径，新退场插值差单列；跨rank同一步激活并退出的旧漏计也补上。',
 '独立energy_planes输出固定lab面粒子能量/动量/电荷/计数以及总自洽E/B有符号能流，full与physical横向孔径分开。没有raw粒子面或场面也可使用；默认空列表关闭。粒子采用实际Boris子步端点标量四动量插值，避免插值三动量后重算gamma的假能损。',
 '场按每个Maxwell步积分，rhythm只控制写入；保留原Yee E/B时间层。J是累计过面能量，J/m是平面场能密度积分，不能当封闭体积存量。无Hann/扣底/辐射分离。',
 '', 'boost短工况（有磁一周期/无磁；各31675粒子、257步；不是旧全量222keV结果）',
 '工况 | 旧最大R[J] | 新最大R[J] | 新终态R[J] | 终态R/abs(净boost粒子动能变化)']
for c,label in [('undulator','有磁'),('free-drift','无磁')]:
    v=science_rows[c]['imbalance'];b=v['boris'];lines.append(f"{label} | {v['legacy']['max_abs_J']:.12g} | {b['max_abs_J']:.12g} | {b['final_J']:.12g} | {100*b['final_over_abs_net_particle_change']:.8g}%")
lines+=['', '有磁旧4.93pJ峰主要是旧外功估计与实际Boris更新不一致；新最大残差约3.03e-5pJ。粒子功残差终值有磁-1.179e-26J、无磁-1.007e-28J，远小于场收支侧残差。无磁仍约净粒子能损32.03%，不能用初始总动能很大掩盖，也不能宣称真实守恒或近场分离已验收。推进/供场/Cowan/Boris状态更新/沉积/CPML均保持。',
 '', '固定lab面（入口-0.05m、出口0.06m；以下为完整横向网格、有符号累计量；单位pJ）',
 '工况 | 粒子K_out-K_in | 场过面能量in | 场过面能量out | (粒子+场)out-in']
for c,label in [('undulator','有磁'),('free-drift','无磁')]:
    t=transport[c];ps=t['planes'];lines.append(f"{label} | {-t['full_particle_kinetic_loss_lab_J']*1e12:.9g} | {ps['entry']['signed_field_energy_lab_J']*1e12:.9g} | {ps['exit']['signed_field_energy_lab_J']*1e12:.9g} | {-t['full_particle_plus_field_entry_minus_exit_lab_J']*1e12:.9g}")
lines+=['两面均31675向下游/0反向，物理孔径也覆盖同数量；旧逐ID面覆盖另行核验，新累计文件本身不存ID，不能只凭总数推广为任意情况下的同群体证明。每面258次场采样、窗外0、缺口0；所存时间段完整不等于无限时间完整脉冲。',
 '上述两面输运差不是lab守恒残差：还缺该时空体积的侧向流和两端时间截面的场存量等。boost总账本追踪到固定trusted cutoff约0.206363825678513m，0.06m面只是中途面，不能强令两者dK相等。入口/出口时间窗也有lab传播延迟。',
 '', '计时（node3，16MPI×8线程，profile关闭，40–240步；3轮全部保留，百分比为配对中位[最小,最大]）',
 '工况 | 新账本/旧账本 | 能量面/新账本无能量面 | 全诊断/关闭']
for c,label in [('undulator','有磁'),('free-drift','无磁')]:
    comp=perf[c]['comparisons'];vals=[]
    for name in ('increment_on_new_vs_old','new_planes_on_vs_no_planes','new_all_diagnostics_on_vs_off'):
        v=comp[name];ps=v['paired_percent_changes'];vals.append(f"{v['median_paired_percent_change']:+.3f}% [{min(ps):+.3f},{max(ps):+.3f}]")
    lines.append(label+' | '+' | '.join(vals))
lines+=['有磁波动很大，配对中位与中位耗时之比甚至可异号；负值不能解释为诊断加速，也不能声称严格零开销。本轮没有old-off对照，无法估计新代码关闭诊断相对旧代码关闭诊断的精确成本。无磁新增能量面三次约4.18–7.43%，仍仅是此短卡证据。不得外推大工况成本。',
 '', '验收与文件',
 '19015:77required、pending统计、2/4rank移动窗与能量面模块PASS；最终reader15正负控PASS；独立Decimal88项最大相对差4.503e-15。',
 '19017:8科学组/6原场与粒子/轨迹逐位比较PASS；19019:双节点4组/2逐位比较PASS；19018:24计时组/18逐位比较PASS。3组科学输出来自同程序已完成r2结果并重新只读验收，旧62H5哈希保持。',
 'SOURCE_BUILD_R2.json固定核心8f2a8ab与binary SHA；二进制编译标签是前HEAD-dirty，逐文件哈希对应已提交核心，不冒称重编HEAD。',
 '图：figures-final-r3/{boost_energy_accounting,lab_plane_transport,diagnostic_performance}.png（PDF同源，三PNG目视通过）；相同目录有完整CSV和来源manifest。',
 '新程序: '+str(binary),
 '使用说明: source/aprl/verification/performance/cases/energy_transport_audit/USAGE.txt',
 '新面文件: campaign-r3/science/{undulator,free-drift}-r1-new-audit-on/detectors/{audit-entry,audit-exit}-energy.h5',
 '总账本: 同case下ledger/energy-ledger.h5；旧/新数据、所有失败尝试、独立审查和完整JSON/CSV均保留。',
 '快照native ABI升v2，旧格式明确拒绝；energy_planes与checkpoint/restart拒绝。本轮未开发快照分支、未重跑全量、未切换默认、未推远程。']
(out/'RESULT_REVIEW.txt').write_text('\n'.join(lines)+'\n')
print(out/'FINAL.json'); print(out/'RESULT_REVIEW.txt')
