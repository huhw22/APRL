"""Prepare review-only deletion allowlists; this module never deletes files."""
import collections
import csv
import json
from pathlib import Path
import sqlite3
import sys


GROUPS=[
 ('A01','旧SPARC六段三算法的克隆输入输出','已撤销独立随机束流基准资格',
  ['sparc_sase_comparison_20261103/compact-production-r1'],
  '同一静态输入1419264条全部属于横向共线复制，每1024条仅64个横向状态；不能作修复后FEL基准。保留全部账本、图表、配置、日志和原输入。'),
 ('A02','旧SPARC长束克隆加载的有界计时输出','计时及历史加载记录',
  ['sparc_sase_comparison_20261103/campaign-r4/timing'],
  '7916544条静态输入全部共线复制；本工况本来也只有有界计时，不能当完整FEL结果。保留耗时与完整性记录。'),
 ('A03','旧十一段两倍滑移束流的克隆加载输出','已撤销独立随机束流基准资格',
  ['sparc_sase_11segment_2slip_20261111_r1/energy-binomial3',
   'sparc_sase_11segment_2slip_20261111_r1/energy-momentum',
   'sparc_sase_11segment_2slip_20261111_r1/genesis'],
  '1734656条输入全部共线复制。保留历史响应及其账本，不当作独立随机源的增益／饱和。'),
 ('A04','旧六段19743和19744的场探测器与参考缓存','保留完整账本及粒子面，原场可退役',
  ['sparc_six_full_ledger_20261112_r3/energy-binomial3/detectors',
   'sparc_six_full_ledger_20261112_r3/energy-binomial3/reference',
   'sparc_six_full_ledger_20261112_r3/energy-momentum/detectors',
   'sparc_six_full_ledger_20261112_r3/energy-momentum/reference'],
  '同一旧克隆输入。这里只列场探测器与可再生参考缓存，保留粒子面、energy_planes、ledger、lab-audit及worldtube。删后不能重新作任意原场谱。世界管原始积分另列B05，不默认混删。'),
 ('A05','旧六段Puffin克隆同输入比较的大输出','相关加载历史响应',
  ['puffin_six_same_ledger_20261113_r1/main','puffin_six_same_ledger_20261113_r1/import-check'],
  '19783与19743/44比较有历史意义，但不能当新独立随机束流基准。保留完整功率/BF/能损CSV、资格纠正和原输入。'),
 ('A06','旧小横向窗口扫描的克隆加载输出','已被独立源56格方案替代',
  ['sparc_sase_transverse_compact_test_20261110_r1',
   'sparc_sase_transverse_compact_test_20261110_r2',
   'sparc_sase_transverse_compact_test_20261111_r3'],
  '来自同一1419264条克隆输入；34/38格另有CIC支撑碰Dirichlet边界的初始化拒绝。原40格结论不能直接推广到独立高斯尾部。'),
 ('A07','Puffin实际磁路由错误的短测原输出','已撤销磁场匹配资格',
  ['puffin_same_sample_20261113_r1/main','puffin_same_sample_20261113_r1/import-check',
   'puffin_same_sample_20261113_r2/main','puffin_same_sample_20261113_r2/import-check',
   'puffin_same_sample_20261113_r3/main','puffin_same_sample_20261113_r3/import-check'],
  '旧31.47pJ路径没有真正采用匹配端场，另触发人工入口匹配；以19770显式lattice及实际首步核验取代。stock smoke、source、scripts、输入与导入初始dump保留。'),
 ('A08','同步辐射稀疏横向节点的早期原输出','未收敛源／网格负控',
  ['synchrotron_shape_validation_20261110','synchrotron_shape_validation_20261110_r2',
   'synchrotron_shape_validation_20261110_r3','synchrotron_shape_feedback_control_20261110_r1',
   'synchrotron_shape_feedback_control_20261110_r2','synchrotron_shape_feedback_attribution_20261110_r1',
   'synchrotron_axis_notch_propagation_20261110_r1','synchrotron_axis_notch_transverse_order12_20261110_r1'],
  '8×8/12×12横向节点沿纵向复制的求积源产生未收敛高k及干涉深零；不是已接受的两组100k独立样本基准。保留所有诊断图表与输入。'),
 ('B01','已取消w64十一段19775的原输出','源正确但任务未完成，非物理错误',
  ['puffin_eleven_four_slip_20261113_r2/main','puffin_eleven_four_slip_20261113_r2/import-check'],
  '用户为降低成本改w355；19775停于18439/79200步。可保留性能记录后删未完成输出，但它不是坏随机源。原输入、source、软件资格和实际初始导入dump保留。'),
 ('B02','准备失败且已替换的少量原输出','未形成有效物理结果',
  ['sparc_random_startup_20261112_r1','sparc_six_full_ledger_20261112_r1',
   'sparc_six_full_ledger_20261112_r2','puffin_eleven_four_slip_20261113_r1',
   'synchrotron_axis_notch_current_filter_20261110_r1'],
  '初始化尾部／输入权重求和、准备维数或不支持的滤波组合拒绝；不等于完整物理结果。保留失败原因、输入、源码和配置。'),
 ('B03','被干净重跑替代的早期权重扫描输出','重复或准备失败，非所有大权重都错误',
  ['macro_weight_radiation_20261109/f3_oneway_large_weights_r1',
   'macro_weight_radiation_20261109/f3_oneway_scan_r1',
   'macro_weight_radiation_20261109/genesis_macro_weight_r1',
   'macro_weight_radiation_20261109/genesis_macro_weight_r2',
   'macro_weight_radiation_20261109/genesis_macro_weight_r3'],
  '已有独立r2/r4干净重跑。大权重噪声本来就是测试目标，不因此判新扫描错误；有效细化及新扫描不列候选。'),
 ('B04','粒子撞边界以去近场的负控原输出','方法失败，原仿真并非无效',
  ['synchrotron_bend_boundary_exit_20261110_r1'],
  '吸收／截断电流在边界引入亮斑，不能制造无近场。保留失败方法图表；已接受bend-only与两组100k原输出不在候选。'),
 ('B05','旧六段账本的worldtube原始积分序列','诊断有效，建议保留；可单独选择退役',
  ['sparc_six_full_ledger_20261112_r3/energy-binomial3/worldtube',
   'sparc_six_full_ledger_20261112_r3/energy-momentum/worldtube'],
  '这不是错误仿真。只有确认不再重算旧固定lab管任意时空边界积分时才选择；删除大bin后，已完成的端侧面／时间盖CSV和图仍在，但不能从原时序重算。账本H5、逐rank审计CSV和粒子面继续保留。'),
]


def under(path,prefix):return path==prefix or path.startswith(prefix+'/')


def collect_protection(workspace,db):
    r=workspace/'results';current=r/'fel_independent_w355_20261114_r1'
    prefixes={str(workspace/'inputs'):'原始粒子输入与变体',str(current):'当前运行／排队工况整体保护'}
    accepted=['accepted_physics_transfer_20261112','energy_ledger_test_record_20261030',
        'physics_extension_20261112_r1','synchrotron_full_ledger_20261112_r1',
        'synchrotron_independent_electron_20261110_r1','synchrotron_independent_electron_grid2x_20261110_r1',
        'synchrotron_independent_electron_100k_seed20261111_20261110_r1',
        'synchrotron_independent_electron_100k_seed20261112_20261110_r1',
        'synchrotron_sample_matched_lw_20261110_r1','synchrotron_sample_matched_lw_20261110_r2',
        'sparc_random_startup_20261112_r2','sparc_single_seed_gap_20261113_r1',
        'sparc_startup_space_charge_theory_20261113_r2','puffin_aprl_source_propagation_20261113_r1',
        'binomial3_fel_20261106/full-r1']
    prefixes.update({str(r/p):'已接受或仍有用途的物理／理论／独立输入基准' for p in accepted})
    exact={};manifests=[]
    for p in current.rglob('*FROZEN*.json'):
        if 'scripts' in p.relative_to(current).parts:continue
        try:a=json.loads(p.read_text())
        except (OSError,ValueError):continue
        if isinstance(a,dict):
            for key in a:
                q=Path(key)
                if q.is_absolute():exact[str(q)]='当前工况冻结清单直接依赖'
            manifests.append(str(p))
    # Preserve input paths explicitly named by historical plans/metadata.
    # This does not treat every historical output in FROZEN as a live dependency.
    names={'PLAN.json','plan.json','INPUT_READY.json','INPUT_VALIDATED.json','MASTER_READY.json','INDEX.json','PREPARED.json','campaign.json'}
    def visit(value,key=''):
        if isinstance(value,dict):
            for k,v in value.items():visit(v,k)
        elif isinstance(value,list):
            for v in value:visit(v,key)
        elif isinstance(value,str) and ('input' in key.lower() or key.lower() in ['master','master_file','source_file']):
            q=Path(value)
            if q.is_absolute() and q.suffix.lower() in ['.h5','.hdf5','.bin','.npy','.npz'] and q.exists():exact[str(q)]='历史记录声明的原始／共同输入'
    for path,size in db.execute("SELECT path,size FROM files WHERE kind='file' AND path LIKE ?",(str(r)+'/%',)):
        p=Path(path);parts=p.relative_to(r).parts
        if p.name not in names or size>4*1024**2 or len(parts)>7 or any(x in ['source','scripts','frozen_source','build-r1'] for x in parts):continue
        try:visit(json.loads(p.read_text()))
        except (OSError,ValueError):pass
    for path,reason in list(exact.items()):
        exact[str(Path(path).resolve())]=reason+'（实际目标）'
    return prefixes,exact,manifests


def main(workspace,output):
    workspace=Path(workspace).resolve();output=Path(output).resolve();db=sqlite3.connect(output/'inventory.sqlite')
    prefixes,exact,manifests=collect_protection(workspace,db)
    rules=dict(prefixes=prefixes,exact=exact,current_frozen_manifests=manifests,
        implicit_keep=['source/scripts/build trees','all original inputs','all ledger/lab-audit files','all plots/CSV/JSON/cards/logs','native initial electron dumps'],
        no_delete_authorized=True)
    (output/'PROTECTION_RULES.json').write_text(json.dumps(rules,indent=2,ensure_ascii=False)+'\n')
    groups=[];candidate=[];excluded=[];candidate_by_inode=collections.defaultdict(list);inode_paths=collections.defaultdict(list)
    all_files=list(db.execute("SELECT path,dev,ino,nlink,size,allocated,mtime_ns,kind FROM files"))
    protected_inodes=set()
    for row in all_files:
        if row[-1]=='file':
            inode_paths[row[1],row[2]].append(row[0])
            if row[0] in exact or any(under(row[0],p) for p in prefixes):protected_inodes.add((row[1],row[2]))
    for ident,title,scope,roots,reason in GROUPS:
        roots=[str(workspace/'results'/p) for p in roots]
        g=dict(id=ident,title=title,category='A_已撤销或不适合作物理基准' if ident.startswith('A') else 'B_另行选择的取消重复或方法负控',
            qualification=scope,roots=roots,reason=reason,status='候选，未批准删除',
            files=0,apparent_bytes=0,unique_allocated_bytes=0,exclusive_reclaim_bytes=0,excluded_files=0,excluded_allocated_bytes=0)
        unique={}
        for row in all_files:
            path,dev,ino,nlink,size,allocated,mtime,kind=row
            if kind!='file' or not any(under(path,p) for p in roots):continue
            p=Path(path);parts=p.relative_to(workspace/'results').parts
            if p.suffix.lower() not in ['.h5','.hdf5','.bin','.raw']:continue
            why=exact.get(path)
            if why is None:why=next((reason for prefix,reason in prefixes.items() if under(path,prefix)),None)
            if why is None and (dev,ino) in protected_inodes:why='与受保护依赖或输入共享inode'
            if why is None and any(x in ['inputs','source','src','scripts','build','frozen_source','ledger','lab-audit','smoke-stock','smoke-matched'] or x.startswith('build-') for x in parts):why='输入／源码／程序／账本保留'
            if why is None and (p.name=='particles.h5' or p.name.endswith('.par.h5') or p.name.endswith('_electrons_0.h5') or p.name.startswith('input')):why='原始或实际初始导入粒子保留'
            if why is None and ident=='A04' and ('particles' in p.name or p.name.endswith('-energy.h5')):why='旧完整账本所需粒子面和energy_planes保留'
            if why is None and any('analysis' in x or x.startswith('review') or x.startswith('diagnosis') for x in parts):why='结论／分析缓存保留'
            if why:
                excluded.append(dict(group=ident,path=path,reason=why,allocated_bytes=allocated));g['excluded_files']+=1;g['excluded_allocated_bytes']+=allocated;continue
            item=dict(group=ident,path=path,dev=dev,inode=ino,nlink=nlink,size_bytes=size,allocated_bytes=allocated,mtime_ns=mtime)
            candidate.append(item);candidate_by_inode[dev,ino].append(item);unique[dev,ino]=allocated
            g['files']+=1;g['apparent_bytes']+=size
        g['unique_allocated_bytes']=sum(unique.values());groups.append(g)
    counts=collections.Counter(x['path'] for x in candidate)
    assert not [p for p,n in counts.items() if n!=1],'Overlapping candidate groups'
    joint=0;shared=[];lookup={x['id']:x for x in groups}
    for key,rows in candidate_by_inode.items():
        all_paths=inode_paths[key];nlink=max(x['nlink'] for x in rows);allocated=max(x['allocated_bytes'] for x in rows)
        selected={x['path'] for x in rows};ids={x['group'] for x in rows}
        complete=len(selected)==nlink and set(all_paths)==selected
        if complete:
            joint+=allocated
            if len(ids)==1:lookup[next(iter(ids))]['exclusive_reclaim_bytes']+=allocated
        if nlink>1:
            shared.append(dict(dev=key[0],inode=key[1],nlink=nlink,candidate_paths=sorted(selected),
                other_observed_paths=sorted(set(all_paths)-selected),unobserved_links=max(0,nlink-len(all_paths)),
                allocation_reclaimed_if_all_candidates_selected=complete,allocated_bytes=allocated))
    a=sum(x['exclusive_reclaim_bytes'] for x in groups if x['id'].startswith('A'));b=sum(x['exclusive_reclaim_bytes'] for x in groups if x['id'].startswith('B'))
    result=dict(status='REVIEW_ONLY_NO_DELETION',workspace=str(workspace),groups=groups,
        candidate_files=len(candidate),candidate_paths_unique=True,all_groups_joint_reclaim_bytes=joint,
        A_individual_reclaim_sum_bytes=a,B_individual_reclaim_sum_bytes=b,
        additional_reclaim_from_cross_group_hardlinks_bytes=joint-a-b,
        byte_definition='allocated blocks512; only count an inode when every hardlink is selected and no protected/unknown link remains',
        original_inputs_preserved=True,ledger_and_lab_audit_directories_preserved=True,
        old_six_particle_and_energy_planes_preserved=True,raw_worldtube_retirement_is_separate_optional_B05=True,
        active_dependency_intersection=0,
        original_paths_unchanged=True,files_deleted=0,hardlink_details=shared)
    (output/'CANDIDATES.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    fields=['group','path','dev','inode','nlink','size_bytes','allocated_bytes','mtime_ns']
    with (output/'CANDIDATE_FILES.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(candidate)
    with (output/'EXCLUDED_PROTECTED_FILES.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['group','path','reason','allocated_bytes']);w.writeheader();w.writerows(excluded)
    lists=output/'file_lists';lists.mkdir(exist_ok=True)
    for g in groups:
        with (lists/(g['id']+'.jsonl')).open('w') as f:
            for item in candidate:
                if item['group']==g['id']:f.write(json.dumps(item,ensure_ascii=False)+'\n')
    for g in groups:print(g['id'],round(g['exclusive_reclaim_bytes']/2**30,3),'GiB',g['files'],'files',g['title'],flush=True)
    print('JOINT',round(joint/2**30,3),'GiB','A',round(a/2**30,3),'B',round(b/2**30,3),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
