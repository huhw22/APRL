"""Archive compact evidence and write the review; no deletion or raw movement."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys


def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def archive(workspace,output):
    required=[
      'sparc_input_loading_regression_20261112_r1/README.md',
      'sparc_input_loading_regression_20261112_r1/current-crossings.json',
      'puffin_six_same_ledger_20261113_r1/loading-audit-20261114-r1/README.md',
      'puffin_six_same_ledger_20261113_r1/loading-audit-20261114-r1/QUALIFICATION.json',
      'sparc_six_full_ledger_20261112_r3/review-20261113-r1/README.md',
      'sparc_six_full_ledger_20261112_r3/review-20261113-r1/RESULT.json',
      'sparc_six_full_ledger_20261112_r3/review-20261113-r1/both_frames.png',
      'sparc_six_full_ledger_20261112_r3/review-20261113-r1/lab_physical_absolute.png',
      'puffin_six_same_ledger_20261113_r1/review-20261114-r1/README.md',
      'puffin_six_same_ledger_20261113_r1/growth-20261114-r1/README.md',
      'puffin_six_same_ledger_20261113_r1/growth-20261114-r1/data/growth.csv',
      'puffin_six_same_ledger_20261113_r1/growth-20261114-r1/data/power_bf_growth.png',
      'puffin_same_sample_20261113_r3/README.md',
      'puffin_same_sample_20261113_r3/analysis/RESULT.json',
      'puffin_same_sample_20261113_r3/analysis/comparison.png',
      'puffin_aprl_source_propagation_20261113_r1/README.md',
      'synchrotron_axis_notch_diagnosis_20261110_r1/RESULT_REVIEW.md',
      'synchrotron_axis_notch_diagnosis_20261110_r1/high_k_diagnosis/high_k_decomposition.png',
      'synchrotron_axis_notch_transverse_order12_20261110_r1/comparison/transverse_sampling_comparison.png',
      'synchrotron_shape_validation_20261110_r3/RESULT_REVIEW.md',
      'accepted_physics_transfer_20261112/03_synchrotron/README.md',
      'accepted_physics_transfer_20261112/03_synchrotron/FAILED_NEAR_FIELD_SEPARATION_METHODS.md',
      'synchrotron_bend_boundary_exit_20261110_r1/analysis/bend_boundary_exit_8eV.png',
      'binomial3_fel_20261106/resonance_origin_20261108/source_generation_audit/README.md',
      'binomial3_fel_20261106/resonance_origin_20261108/source_generation_audit/METHOD_AUDIT.md',
      'binomial3_fel_20261106/resonance_origin_20261110/genesis_particle_state_control_r1/README.md',
      'binomial3_fel_20261106/resonance_origin_20261110/genesis_particle_state_control_r1/genesis_particle_state_comparison.png',
      'macro_weight_radiation_20261109/README.md',
      'puffin_eleven_four_slip_20261113_r2/PLAN.json',
      'puffin_eleven_four_slip_20261113_r2/NATIVE_VALIDATED.json',
      'puffin_eleven_four_slip_20261113_r2/main/puffin.log',
      'puffin_eleven_four_slip_20261113_r1/PREPARATION_STOPPED.json',
      'sparc_random_startup_20261112_r1/README.md',
      'sparc_random_startup_20261112_r1/SUPERSEDED.json',
      'ledger_discrete_diagnosis_20261029/RESULT_REVIEW.txt',
    ]
    entries=[]
    for name in required:
        source=workspace/'results'/name
        if not source.is_file():raise FileNotFoundError(source)
        assert source.suffix.lower() not in ['.h5','.hdf5','.bin','.raw']
        assert source.stat().st_size<32*1024**2,source
        destination=output/'archive/results'/name;destination.parent.mkdir(parents=True,exist_ok=True)
        before=source.stat();shutil.copy2(source,destination)
        digest=sha(source);assert digest==sha(destination)
        after=source.stat();assert (before.st_ino,before.st_size,before.st_mtime_ns)==(after.st_ino,after.st_size,after.st_mtime_ns)
        entries.append(dict(source=str(source),archive=str(destination.relative_to(output)),size=after.st_size,
            sha256=digest,source_inode=after.st_ino,source_mtime_ns=after.st_mtime_ns))
    (output/'ARCHIVE_MANIFEST.json').write_text(json.dumps(dict(complete=True,files=entries,
        payload_bytes=sum(x['size'] for x in entries),copy_verified=True,originals_unchanged=True,
        h5_or_raw_archived=False,archived_original_markdown_links_preserved=True),indent=2)+'\n')
    return entries


def main(workspace,output):
    workspace=Path(workspace).resolve();output=Path(output).resolve()
    manifest=archive(workspace,output);plan=read(output/'CANDIDATES.json');groups=plan['groups']
    byid={g['id']:g for g in groups};now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    def gb(value):return f'{value/2**30:.2f}'
    def evidence(path,label):return f'[{label}](archive/results/{path})'
    load=evidence('sparc_input_loading_regression_20261112_r1/README.md','加载回归原因')
    load2=evidence('puffin_six_same_ledger_20261113_r1/loading-audit-20261114-r1/README.md','旧六段资格复核')
    ledger=evidence('sparc_six_full_ledger_20261112_r3/review-20261113-r1/README.md','完整账本结论')
    puffin=evidence('puffin_aprl_source_propagation_20261113_r1/README.md','实际磁路由纠正与源传播诊断')
    near=evidence('accepted_physics_transfer_20261112/03_synchrotron/FAILED_NEAR_FIELD_SEPARATION_METHODS.md','五种近场分离方法的完整记录')
    source=evidence('binomial3_fel_20261106/resonance_origin_20261108/source_generation_audit/README.md','硬截断源的全量审计')
    state=evidence('binomial3_fel_20261106/resonance_origin_20261110/genesis_particle_state_control_r1/README.md','Genesis粒子状态保持对照')
    notch=evidence('synchrotron_axis_notch_diagnosis_20261110_r1/RESULT_REVIEW.md','深零与高k分解')
    correct=evidence('accepted_physics_transfer_20261112/03_synchrotron/README.md','已接受的同步辐射范围')
    first=[g['id'] for g in groups if g['id'].startswith('A')]+['B01','B03']
    first_bytes=sum(byid[i]['exclusive_reclaim_bytes'] for i in first)
    body=f'''# 错误结果大合集

旧结果中已经确认的问题主要有三类：把纵向复制加载当作独立随机束流、实际执行的外磁场与目标卡不一致、把未收敛或含近场的观测量当作纯辐射基准。它们需要退役的是特定物理结论，并不意味着所有原始场或账本都是无意义数据。

本次整理了 **{plan['candidate_files']} 个原始输出文件**，全部候选同时删除最多可释放约 **{gb(plan['all_groups_joint_reclaim_bytes'])} GiB**。其中 A01–A08 为已撤销或不适合作当前物理基准的大输出，合计 **{gb(plan['A_individual_reclaim_sum_bytes'])} GiB**；B类另行说明取消、重复和有效负控的区别。原始输入、图表、CSV、配置、日志和代码保留。**尚未删除或移动任何原结果。**

可先考虑 A01–A08、B01、B03，合计约 **{gb(first_bytes)} GiB**。B05涉及旧账本的原始世界管积分，建议先保留。完整路径、影响和逐文件清单见[删除候选确认表](DELETE_CANDIDATES.md)。

## 纵向复制加载为何让旧FEL基准失去资格

GENESIS原生quiet/Fawley输入每个光学片放1024个宏粒子，但只有64组独立横向及动量状态，每组沿纵向复制16份。把这种输入直接送入三维自洽PIC，会引入人为横纵向联合相关。全束BF、几个谐波BF和两程序入口数组一致，都不能排除这种相关。错误是将这类源宣称为独立随机六维电子束；纵向自场本身并不是应被删除的物理效应。

这次对三份已闭合历史输入逐片重查：

| 历史输入 | 记录数 | 检查光片数 | 每片不同横向位置 | 共线复制比例 |
|---|---:|---:|---:|---:|
| 旧短束六段／小窗口输入 |1,419,264|1,386|64|100%|
| 旧长束有界计时输入 |7,916,544|7,731|64|100%|
| 旧两倍滑移十一段／19743/44/19783共同输入 |1,734,656|1,694|64|100%|

完整结果在[本次历史输入审计](HISTORICAL_LOADING_AUDIT.json)。后一次回归不是独立生成器重新坏了，而是SPARC专用prepare→native-noise→convert_input绕过了已经验证的独立固定N采样流程，随后又复用旧文件。{load}；{load2}。

这些旧曲线可以说明相关加载的历史响应和算法敏感性；不能证明修复后的APRL仍有同样快的FEL启动，也不能与不同束长、权重和观测群的独立随机首段曲线直接相减作因果归因。此前新的独立w64首段，固定中央ID群BF出口PIC约0.004448、GENESIS约0.004746；正确磁路由的7周期Puffin全束BF与PIC相差约0.32%，这些有效对照继续保留。

经验：输入准入必须记录生成链和实际文件身份，并检查横向重复、联合六维统计及真实导入数组；零重复本身仍不是全部统计独立性的证明。重新使用旧输入时，必须重新核对其资格。

## Puffin磁场函数正确不等于实际运行磁场正确

旧短测的外磁函数探针通过、七列入口导入逐位相同，却没有保证正式推进真的使用该外磁选项。无显式lattice的单段路径未正确设置段磁场类型，且另一开关触发了不适用于自带平滑入口的人工matchIn/Out。旧31.47pJ与PIC约24.22pJ的“磁场已匹配但仍差30%”结论因此撤销。

修正控制19770使用显式UN磁路由及原一周期平滑端部，实际首物理步轨迹也通过。修正后的原生节点结果约29.29pJ，而在PIC像素处双线性求场后约22.56pJ；后者不能直接命名为唯一正确能量。{puffin}。

经验：既检查静态场公式，也检查经过初始化、端部处理后的真实场调用路径和首步轨迹。还要选末尾correctTrans之前的物理粒子dump，不能把收尾重定心后的粒子当作同一出口事件。

## 稀疏横向节点造成的谱深零不是辐射消失

早期同步辐射用8×8或12×12 Gauss–Hermite横向求积节点，并沿纵向节点复制。它可以是确定性积分源，却不能被不加区别地当成独立随机电子样本。稀疏源沉积到PIC网格后产生未收敛高横向波数，随后传播相位与低角场干涉，某些轴上点便接近零。

9.825eV处，8×8原始轴上值约0.0827nJ/(eV·m²)，并非全平面能量也为零。换12×12后深零消失，却变成约1087.98nJ/(eV·m²)的强峰；所以“图不再有零点”不构成正确性验收。{notch}。

后续两组100k独立样本与同粒子LW相干和的4–12eV积分比分别为0.93895和0.98751。它们与连续高斯期望谱不同，主要不能再一概归为PIC传播错误。有限样本必须和同一实际样本比较；余下振荡仍包含粒子沉积及网格离散影响。{correct}。

经验：同时检查轴谱、全平面能量、高k含量和同样本理论。增加横向节点数、缩小网格或改变滤波后某个点变好，都不能单独替代收敛和绝对量对照。

## 近场分离失败不等于bend-only原场失败

以下方法在当前工况不能可靠恢复全平面纯加速度辐射：

| 方法 | 问题及保留的教训 |
|---|---|
| bend减straight复场后再画强度 | 两种轨道的速度场中心不同，差场中会同时出现两个位置；不能解释成一束粒子分裂。|
| 事后减解析速度场 | 连续LW与离散PIC的沉积、交错、边界、时间记录和Green函数不一致，大背景相减会放大小复场误差。|
| 总强度减速度场强度 | 漏掉辐射与速度场的干涉交叉项，物理量定义不成立。|
| 让粒子撞横向边界 | 突然终止电流和吸收瞬态产生亮斑，不能冒充近场自然离开。|
| 小角度截断 | 变平滑的同时删除大量真实理论辐射，必须报告接受度及保留比例。|
| 单平面前向真空模投影 | 当前可表示速度场也主要落在前向传播锥内，单平面不足以唯一辨别源区近场。|

事后速度场相减的两个样本，全平面PIC/纯加速度LW积分比分别约9.825和8.087，不能作为纯辐射图。正确主结论继续采用bend-only原始PIC与完整LW比较，远离束流的明确孔径才与纯加速度LW比较；孔径只收集部分辐射，不能冒称全空间能量。{near}。

经验：先定义可独立计算、可测量的物理量，再看图；不能用更平滑、能量更接近或双峰消失来选择一种近场扣除方法。B04只把撞边界负控的大输出列为可选，不把正确bend-only结果列入删除。

## Genesis one-for-one的真正差别在推进中的重分片

“开启one-for-one会在入口重新生成全部六维状态”曾是一个猜测，后来实际核查否定了它：旧导入的x/y/px/py/gamma/theta/current与原H5逐位一致。明确存在的是推进时按theta跨片迁移、相位回绕和按新人数重算片电流。

固定11,438,389个入口事件，只关闭这条重排路径，目标带Rect能量由约0.0784pJ变为22.89pJ，Hann由约0.0798pJ变为15.82pJ。但关闭后仍有电流斜率带来的确定性片内源，不能据更大能量就宣布更物理。{state}。

经验：区分入口加载、运行中重分片、源表示和实际电流。对照必须检查真实入口及沿程处理；不能把模式名当成状态保持保证，也不能用事后幅度归一化掩盖源定义的差别。

## 独立随机源也可能有很强的硬截断边缘辐射

旧21λ短源来自长抛物线高流强中央区域，两端在约53.26A处直接到零。实际全量输入未发现共线复制，去平均包络后的随机涨落也处在正常散粒幅度；邻频大形状因子来自截断包络。随机生成并不会消除确定性的头尾旁瓣。{source}。

需要分清生成时截断真实源、分析时加Hann、以及两程序重建片电流：它们分别改变源、观测窗和数值源表示。只在精确载频BF较小，不能说明整个目标频带都由散粒噪声主导。

经验：生成前检查全目标频带的N|F(k)|²及端部电流支持。独立one-for-one测试保留为边界与相干响应的证据，不能因为有硬边缘就一概称为坏随机源。本次未把约313GiB的原完整37周期one-for-one数据列入候选。

## 统计窗口和横向插值不能被当成增益修正

短束团的头尾启动、滑移和记录截止可能占据很大比例；除以完整记录时长所得平均功率，不能直接与论文中央稳态瞬时功率或饱和功率比较。Hann频带能量除以哪一段物理窗口，必须在比较前固定；加窗能量也不能作为未加窗全量账本的收支项。后半段取窗仍需报告实际覆盖，不能因为曲线不成平台就宣称找到了饱和。

Puffin的两个横向实场分量共用节点，并无Yee六分量半格错位；粒子位置的三线性插值是另一件事。APRL探测器对交错场先作空间平均再计算强度，会与先平方后积分得到不同量。对同一Puffin场作四点平均时，约70%的强度下降来自相位差项，这不等于PIC实际漏记了同样多能量。不同观测算子不能排成一串可相加的能量修正，也不能用总量接近宣布复场一致。{puffin}。

经验：增益、饱和、局部BF、全束BF及场强度必须各自固定窗口、粒子群、孔径和取样位置；统计量定义的改变不能当作物理增益改变。本次保留历史增长CSV及其资格说明，不将这些统计问题自动扩大为原场全部无效。

## 全量账本需要正确时空区域和全部端盖

旧cell-centre能量、混合时间层和梯形Poynting不能自动满足离散场算子的严格恒等式；正确真空传播也能在这种账本定义下出现差额。这个诊断错误不能直接证明场推进物理错误，更不能以初始总束流能量作分母掩盖小净能损。相关推导和成功账本验收继续保留。

旧克隆六段19743/44的boost离散总差约6pJ，而lab能量—动量账本仍差约3.16/2.88µJ；不能声称动量匹配已证明lab物理闭合。第六段还缺约1.67万个晚粒子时间端盖ID的实际动能，未知量不能填零或用差额反推。{ledger}。

经验：场算子账闭合、粒子功配对和固定lab世界管是不同核对。完整端面、侧面、时间盖及真实粒子事件必须对应同一窗口。A04保留账本H5、逐rank审计CSV、粒子面和energy_planes；世界管大bin另列B05，默认建议保留。

## 不把预期噪声和失败作业标签当作删除依据

提高宏权重后散粒噪声增强，是不少扫描有意要观察的现象；w=64或w=512本身不是错误标签。有效权重曲线、真实抽样输入和最新独立源对照继续保留。B03只列已经有干净重跑替代的旧批次。

Slurm的FAILED也不等于物理数据全坏：旧Puffin六段的物理推进完成，失败只在后处理舍入容差；经独立核查的结果仍有资格。反过来，COMPLETED也不保证实际磁场、加载或观测定义正确。准备阶段的维数、输入求和、尾部越界或不支持组合拒绝，往往根本没有产生完整物理输出。

已取消19775的独立w64源是正确生成的，只因计算成本而由w355替代；B01属于未完成输出退役，不是物理错误。原两节点19797未启动，几乎没有科学数据可清理；当前单节点19799和其他三个运行／排队任务全部保护。

## 清理时保留哪些证据

每个候选只列明确的H5/bin/raw大输出，原路径不搬动，原始输入和代码不在删除清单。全部PNG/PDF/CSV/JSON、配置和日志继续留在原处，并将上述关键结论、图和审计复制成[轻量证据快照](ARCHIVE_MANIFEST.json)。快照是逐字节核验副本；原文保留历史链接，其所指大数据可能在后续获准删除，本文的关键证据入口已复制到合集。

当前19795、19796_0、19796_1、19799对应的整个工况目录及其冻结依赖受保护；有效两组100k同步辐射、独立首段、正确Puffin19770、已认可物理转移包、空间电荷与账本验收不在候选。较大的166.84GiB corrected_full_comparison等目录，没有已撤销依据，本次也不因占空间大而列为错误。

空间按真实分配块统计，考虑硬链接；删除一个别名仍被其他路径引用时，不计为释放。清单是本次快照，执行前还须重新检查文件身份和当时的作业依赖。计量明细见[盘点记录](INVENTORY.json)、[候选机器清单](CANDIDATES.json)、[全部结果目录处置表](DIRECTORY_DISPOSITION.csv)及[保护规则](PROTECTION_RULES.json)。

整理时间：{now}。状态：待用户逐项确认，删除数为0。
'''
    (output/'README.md').write_text(body)
    table=['# 可以删除的结果确认表','',
      '**本表尚未执行任何删除。** 每个编号只对应其JSONL中的具体文件，不代表删除整目录。原始输入、图表、配置、日志、代码和规定保留的账本不随组删除。','',
      f'所有候选上限 **{gb(plan["all_groups_joint_reclaim_bytes"])} GiB**；建议先考虑 A01–A08、B01、B03，约 **{gb(first_bytes)} GiB**。B05建议保留以便重算旧账本。','',
      '| 编号 | 可以释放 GiB | 文件数 | 内容 | 资格与影响 |','|---|---:|---:|---|---|']
    for g in groups:
        table.append(f'| [{g["id"]}](file_lists/{g["id"]}.jsonl) | {gb(g["exclusive_reclaim_bytes"])} | {g["files"]} | {g["title"]} | {g["qualification"]} |')
    table+=['','GiB=2³⁰字节；不是文件名义大小之和。上表按每组可独立释放的inode计量，组合硬链接增量为'+gb(plan['additional_reclaim_from_cross_group_hardlinks_bytes'])+' GiB。B02当前为0，没有实际大文件候选。','',
      '## 各编号的实际路径和保留内容','']
    for g in groups:
        table += [f'### {g["id"]} {g["title"]}','',g['reason'],'',f'可释放 **{gb(g["exclusive_reclaim_bytes"])} GiB**，{g["files"]}个文件。完整文件清单：[JSONL](file_lists/{g["id"]}.jsonl)。','', '原目录范围：','']
        table += ['- `'+p+'`' for p in g['roots']]
        table+=['']
    table+=['## 默认保留且不在上述清单','',
      '- 当前19795、19796_0、19796_1、19799及它们的全部输入、程序和冻结依赖。',
      '- `inputs/`、历史计划明确声明的共同输入、native初始导入dump、所有source/scripts/build和原配置。',
      '- 全部ledger、lab-audit目录。19743/44的粒子面和energy_planes也保留；worldtube原积分只有B05单独选择才会删除。',
      '- 已接受两组100k同步辐射及其同样本LW、正确Puffin19770、独立随机首段、空间电荷／真空脉冲与账本成功验收。',
      '- 原完整37周期独立one-for-one数据约313GiB，以及其他没有明确退役依据的大型正确结果。','',
      '## 确认方式','',
      '可以按编号回复，例如“删除A01、A02、A03、A05、B01，保留其余”。未明确选中的编号继续保留。确认后还会按当时的运行依赖和文件inode/大小/mtime重新核对，再仅删除获准清单中的文件。','',
      '[错误结果大合集](README.md) · [逐文件CSV](CANDIDATE_FILES.csv) · [受保护排除项](EXCLUDED_PROTECTED_FILES.csv) · [机器清单](CANDIDATES.json)','']
    (output/'DELETE_CANDIDATES.md').write_text('\n'.join(table))
    with (output/'DECISIONS.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id','title','reclaim_GiB','decision','status'])
        for g in groups:w.writerow([g['id'],g['title'],gb(g['exclusive_reclaim_bytes']),'','未确认'])
    summary=dict(complete=True,review_only=True,files_deleted=0,raw_results_moved=0,
        candidate_files=plan['candidate_files'],all_candidate_reclaim_GiB=plan['all_groups_joint_reclaim_bytes']/2**30,
        first_pass_suggestion_groups=first,first_pass_reclaim_GiB=first_bytes/2**30,
        suggested_groups_are_not_authorization=True,archive_files=len(manifest),archive_bytes=sum(x['size'] for x in manifest),
        user_confirmation_pending=True,current_jobs_protected=[19795,'19796_0','19796_1',19799])
    (output/'DELIVERY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(summary,ensure_ascii=False),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
