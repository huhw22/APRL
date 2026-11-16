"""Publish completed bend results into a new version of the reviewed delivery."""
import csv,hashlib,json,math,re,shutil,subprocess,sys,zipfile
from decimal import Decimal,localcontext
from pathlib import Path,PurePosixPath
import numpy as np

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()
def save(p,v):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def table(p,v):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(v[0]));w.writeheader();w.writerows(v)
def identity(p):
    s=Path(p).stat();return dict(inode=s.st_ino,size=s.st_size,mtime_ns=s.st_mtime_ns)

def main(ws,analysis,output,archive):
    ws,analysis,root,archive=map(lambda p:Path(p).resolve(),[ws,analysis,output,archive])
    raw=ws/'results/synchrotron_full_ledger_20261116_r1';result=read(analysis/'RESULT.json')
    assert result['complete'] and read(analysis/'BOOST_DISPLAY_REVIEW.json')['passed']
    assert not root.exists() and not archive.exists();root.mkdir()
    for path,before in read(analysis/'PROTECTED_OUTPUTS.json').items():assert identity(path)==before,path
    for name in ['FROZEN.json','PROTECTED.json']:
        for p,expected in read(raw/name).items():assert sha(p)==expected,p
    for p,expected in read(analysis/'FROZEN_ANALYSIS.json').items():assert sha(p)==expected,p
    old=ws/'results/physics_delivery_with_ppt_guide_20261116_r1.zip'
    assert sha(old)=='b63403f9315b6823a8d554db8d0cda1598f4923bd9785f47cee48d4239ce8b5a'
    incoming=[]
    with zipfile.ZipFile(old) as z:
        assert z.testzip() is None
        prefix='physics_delivery_with_ppt_guide_20261116_r1/'
        for n in z.namelist():
            assert n.startswith(prefix)
            rel=PurePosixPath(n[len(prefix):]);assert not rel.is_absolute() and '..' not in rel.parts
            target=root.joinpath(*rel.parts);target.parent.mkdir(parents=True,exist_ok=True)
            with z.open(n) as src,target.open('xb') as dest:shutil.copyfileobj(src,dest)
            incoming.append(dict(member=n,destination=str(rel),sha256=sha(target)))
    previous=root/'metadata/previous_delivery';previous.mkdir()
    for name in ['README.md','PPT制作指导与结果总结.md','SHA256SUMS.txt']:
        shutil.copyfile(root/name,previous/(name+'.txt' if name.endswith('.md') else name))
    shutil.copyfile(root/'metadata/VALIDATION.json',previous/'VALIDATION.json')
    shutil.copyfile(root/'metadata/SOURCE_COMMIT.json',previous/'SOURCE_COMMIT.json')
    current=root/'06_synchrotron_current';current.mkdir()
    copied={}
    def copy(src,dest):
        p=current/dest;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,p)
        assert sha(src)==sha(p);copied[str(src)]=sha(src)
    for name in ['RESULT.json','PLOTS.json','BOOST_DISPLAY_REVIEW.json','SOURCE_REUSE.json','JOB_STATUS.txt']:
        copy(analysis/name,Path(name))
    for directory in ['cases','figures','radiation']:
        for p in sorted((analysis/directory).rglob('*')):
            if p.is_file() and p.suffix in ['.csv','.json','.png','.pdf']:
                copy(p,p.relative_to(analysis))
    for name in ['JOB_LINKS.json','COMPLETE_linear-binomial3.json','COMPLETE_quintic-compensated.json','STATIC_VALIDATION.json','PREFLIGHT_COMPLETE.json']:
        copy(raw/name,Path('records')/name)
    summary=[];orbits=[];worst=Decimal(0)
    with localcontext() as ctx:
        ctx.prec=60
        for case in result['cases']:
            p=current/'cases'/case['tag'];lab=rows(p/'lab/complete_components.csv')
            for row in lab:
                total=sum(Decimal(row[k]) for k in ['particle_net_transport_J','particle_time_caps_J','field_ports_net_J','field_transverse_net_J','field_time_caps_net_J'])
                error=abs(total-Decimal(row['physical_lab_sum_J']));worst=max(worst,error)
                assert error<Decimal('1e-25')
            val=case['fixed_lab'];last=val['final_components_J']
            sums=[float(x['physical_lab_sum_J']) for x in lab]
            entry=case['particle_statistics'][0];exit=case['particle_statistics'][-1]
            theta=math.asin(.0024/3.) if case['kind']=='bend' else 0.
            expected_x=entry['mean_x_original_m']+(3*(1-math.cos(theta))+(.3524-.0524)*math.tan(theta) if theta else 0.)
            orbit=dict(case=case['tag'],mean_exit_x_deviation_from_external_orbit_m=exit['mean_x_original_m']-expected_x,
                mean_exit_angle_deviation_rad=exit['mean_x_angle_rad']-theta,exit_rms_angle_rad=exit['rms_x_angle_rad'],
                scope='Mean transverse orbit only; no bound on complete per-electron radiation phase or exact self-consistent LW.')
            orbits.append(orbit)
            summary.append(dict(case=case['tag'],method=case['method'],seed=case['seed'],kind=case['kind'],steps=case['steps'],
                boost_discrete_final_sum_J=case['boost']['final_discrete_sum_J'],
                boost_physical_final_sum_J=case['boost']['final_physical_component_sum_J'],
                lab_particle_change_J=last['particle_net_transport_J'],lab_field_ports_J=last['field_ports_net_J'],
                lab_field_sides_J=last['field_transverse_net_J'],lab_field_caps_J=last['field_time_caps_net_J'],
                lab_field_total_J=last['independent_field_gain_J'],lab_final_sum_J=last['physical_lab_sum_J'],
                lab_final_percent_fixed_D=val['last_complete_residual_percent'],lab_sum_range_J=max(sums)-min(sums),
                lab_fixed_D_J=val['normalization_denominator_J'],mean_particle_loss_eV=exit['mean_loss_eV'],
                all_lab_budgets_complete=val['all_physical_budgets_complete']))
    table(current/'SUMMARY.csv',summary);table(current/'ORBIT_MEAN_REVIEW.csv',orbits)
    lines=['# 自洽弯铁同步辐射 三点核与五次补偿','',
        '当前主结果来自19814/19815，两项COMPLETED 0:0，用时2:07:53/2:12:44。两份原100k独立单位电子样本，每方法包含弯转和匹配直线，共八例。实际1079/1084步，全部输出关闭；所有粒子面均100000条，世界管和16rank审计完整。只读后处理19821完成40秒，没有物理重跑。','',
        '## 完整能量账本','',
        '| 方法与样本 | 粒子净能变 fJ | 场输运与存量 fJ | 完整lab加和 fJ | 加和/固定D % |','|---|---:|---:|---:|---:|']
    for x in summary:
        if x['kind']=='bend':lines.append(f"| {x['method']} {x['seed']} | {x['lab_particle_change_J']*1e15:+.4f} | {x['lab_field_total_J']*1e15:+.4f} | {x['lab_final_sum_J']*1e15:+.4f} | {x['lab_final_percent_fixed_D']:+.5f} |")
    lines+=['',
        'lab完整加和包含粒子、纵向入口负/出口正、四侧净流、晚减早时间端面存量。56个管预算均完整，逐ID/电荷检查通过，粒子侧流和时间端盖在这些记录中确认为零，未以残差补项。几何为新坐标中心半宽345/265µm，起点z=.02m，终点.3524m；迟滞时间端盖C=积分(u−Sz/c)dV。','',
        '直线对照的末面残差绝对值约0.16–0.34fJ，但粒子净能变只有约8–9fJ，固定尺度百分比可达数个百分点。全部八例的终值、加和max−min、固定分母见SUMMARY.csv，不能仅按百分比跨工况排序。','',
        '计算系分别保留物理分项与独立离散通道。物理图中的前端注入、后端裁出、四侧、移窗接口、外场与约束源功全部计入；另表逐项列出Boris外磁功、边界源功、冻结与载流约束功。boost_internal_balance图用已独立记录的外磁功作代数消元，以免它掩盖较小的内场/边界交换，不改变加和值，也不从残差反推外功。完整离散加和末值约1e−24至2e−23J；近零离散加和不代表连续物理分量之和恰为零。该弯铁结果也不能代替六段FEL中尚非零的26.51%/37.56%账本。','',
        '## LW参考与不同观测量','',
        '| 方法与样本 | 8eV全平面PIC/完整LW | 75µm孔径PIC/加速度LW | 轴上8–12eV PIC/LW |','|---|---:|---:|---:|']
    for x in result['radiation']:lines.append(f"| {x['method']} {x['seed']} | {x['PIC_over_full_LW_full_plane']:.6f} | {x['PIC_over_acceleration_LW_aperture']:.6f} | {x['axis_PIC_over_prescribed_LW_8_12eV']:.6f} |")
    lines+=['',
        '全平面为8eV去掉CPML后的正定横向E/B谱场强度积分，不是全频有符号能量账本。孔径半径75µm、中心原lab辐射轴，比较加速度LW；轴上频谱是弯转减直线的复场诊断，不称作唯一分离的纯辐射。三者不得合并成一个精度指标。没有幅度拟合，没有从新PIC中扣除解析近场。','',
        'LW参考使用已验证的相同两份真实100k输入，逐电子横向位置、到达时间和权重相干求和，速度项与辐射项均保留；新旧记录时间数组逐点相同。坐标恢复x_original=x_output+120µm后才比较，原辐射轴在输出x=−120µm。','',
        '该参考仍采用规定外磁轨迹，未声称已重建新自洽的每电子完整世界线。新输出末面平均横向轨迹相对该模型偏移约34–36nm，平均角差约0.11–0.12µrad；ORBIT_MEAN_REVIEW.csv只是均值核对，不给完整相位误差界。PIC/reference差异包含反馈轨迹与离散误差，不能全部称作场求解精度。','',
        '## 使用图片与数据','',
        '- [八例账本汇总](SUMMARY.csv)','- [光学对照指标](radiation/comparison_metrics.csv)',
        '- [平均轨迹近似检查](ORBIT_MEAN_REVIEW.csv)',
        '- [三点核弯铁lab账本样本一](figures/zh-CN/F3_20261111_bend_lab_balance.png)',
        '- [五次补偿弯铁lab账本样本一](figures/zh-CN/MM_20261111_bend_lab_balance.png)',
        '- [计算系内场交换：已扣独立记录的外磁功](figures/zh-CN/F3_20261111_bend_boost_internal_balance.png)',
        '- [样本一轴上频谱](figures/zh-CN/20261111_axis_spectrum.png)',
        '- [样本一三点核总场](figures/zh-CN/20261111_F3_8eV.png)',
        '- [样本一完整LW总场](figures/zh-CN/20261111_LW_full_8eV.png)',
        '- [LW完整公式推导](../01_reviewed_physics/docs/03_LW_VELOCITY_RADIATION_AND_NEAR_FIELD.md)','',
        'figures/zh-CN与figures/en各48张独立图，PNG/PDF同源；cases目录保留两样本、两方法、弯铁/直线的所有时间与位置曲线CSV。所有图都只有一个数据坐标区，场图另外保留色标。无直方图、无联立图、无下标角标参数框。旧七点无反馈图仅在历史部分。']
    (current/'README.md').write_text('\n'.join(lines)+'\n')
    template=Path(__file__).resolve().parent
    guide=template.parent/'physics_delivery_ppt_guide_20261116/PPT_GUIDE.md'
    shutil.copyfile(guide,root/'PPT制作指导与结果总结.md')
    assets=rows(root/'PPT_ASSET_MAP.csv')
    assets=[a for a in assets if a['module']!='弯铁新两方法']
    fields=list(assets[0])
    for lang,name,title in [('zh-CN','F3_20261111_bend_boost_internal_balance','三点核计算系内场交换（已扣记录的外磁功）'),
        ('zh-CN','MM_20261111_bend_boost_internal_balance','五次补偿计算系内场交换（已扣记录的外磁功）'),
        ('zh-CN','F3_20261111_bend_lab_balance','三点核弯铁lab账本'),
        ('zh-CN','MM_20261111_bend_lab_balance','五次补偿弯铁lab账本'),
        ('zh-CN','20261111_axis_spectrum','弯铁样本一轴上频谱'),('zh-CN','20261112_axis_spectrum','弯铁样本二轴上频谱'),
        ('zh-CN','20261111_F3_8eV','8eV三点核总场'),('zh-CN','20261111_MM_8eV','8eV五次补偿总场'),
        ('zh-CN','20261111_LW_full_8eV','8eV规定轨迹完整LW'),('zh-CN','20261111_LW_acc_8eV','8eV规定轨迹加速度LW')]:
        assets.append(dict(module='当前自洽弯铁',title=title,preferred_asset=f'06_synchrotron_current/figures/{lang}/{name}.png',
            numerical_or_derivation_source='06_synchrotron_current/SUMMARY.csv' if 'balance' in name else '06_synchrotron_current/radiation/comparison_metrics.csv',
            english_asset=f'06_synchrotron_current/figures/en/{name}.png',PPT_use='当前单图主素材；LW为规定外磁轨迹参考，完整分项与指标口径见README'))
    for a in assets:
        for k in ['preferred_asset','numerical_or_derivation_source','english_asset']:
            if a[k]:assert (root/a[k]).is_file(),a[k]
    table(root/'PPT_ASSET_MAP.csv',assets)
    (root/'README.md').write_text('# 物理结果与 PPT 制作指导 更新交付\n\n'
        '当前弯铁三点核/五次补偿自洽结果已完成并整理。先读[PPT制作指导与结果总结](PPT制作指导与结果总结.md)，按[素材映射表](PPT_ASSET_MAP.csv)取用。\n\n'
        '- [当前弯铁两方法结果、LW和完整能量账本](06_synchrotron_current/README.md)\n'
        '- [真空、短漂移、球形空间电荷与历史LW](01_reviewed_physics/README.md)\n'
        '- [独立六段FEL增长与两方法账本](02_FEL_six/README.md)\n'
        '- [十一段Puffin饱和](03_FEL_eleven/README.md)\n'
        '- [谢明与空间电荷纯公式推导](04_formula_derivations/FEL_THEORY.md)\n\n'
        '05为重提时的历史提交记录；06是当前有效弯铁展示入口，原七点无反馈结果保留为历史。此前总包中的“新弯铁待补”现由06完成结果取代。\n\n'
        '新弯铁四例lab末面残差约0.26%–0.36%，非零值照实保留；LW采用同输入同记录窗的规定外磁轨迹参考，不等于新自洽世界线的精确LW。六段FEL原26.51%/37.56%残差不受此结果替代。完整真实束团PIC/ELEGANT仍仅有理论/统计定义随包，未虚构数值完成。\n\n'
        '仅含图片、PDF、视频、公式与结果CSV/JSON/原时序NPZ，无代码、可执行程序或科学H5。中英文图保留，没有生成实际PPT。旧来源与科学输出字节保持，主导航和指导均使用包内路径。\n')
    # All imported scientific assets remain byte-for-byte; only current navigation is revised.
    changed={'README.md','PPT制作指导与结果总结.md','PPT_ASSET_MAP.csv','SHA256SUMS.txt','metadata/VALIDATION.json','metadata/SOURCE_COMMIT.json'}
    for x in incoming:
        if x['destination'] not in changed:assert sha(root/x['destination'])==x['sha256']
    links=0
    for p in root.rglob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            target=target.split('#')[0]
            if not target or target.startswith(('http://','https://','mailto:')):continue
            dest=(p.parent/target).resolve();assert dest.is_relative_to(root) and dest.exists(),(str(p),target)
            links+=1
    plots=read(analysis/'PLOTS.json')['figures'];assert len(plots)==96 and all(p['data_axes']==1 for p in plots)
    for p,expected in copied.items():assert sha(p)==expected,p
    for p,before in read(analysis/'PROTECTED_OUTPUTS.json').items():assert identity(p)==before,p
    commit=subprocess.check_output(['git','-C',str(template.parents[2]),'rev-parse','HEAD'],text=True).strip()
    save(root/'metadata/SOURCE_COMMIT.json',dict(commit=commit,pushed=False,templates_sha256={p.name:sha(p) for p in template.iterdir() if p.is_file()}))
    save(root/'metadata/BEND_SOURCE_OUTPUT_SHA256.json',copied)
    save(root/'metadata/VALIDATION.json',dict(passed=True,closed_physical_cases=8,read_only_postprocess_job=19821,
        physical_boost_source_review_job=19822,original_output_identities_unchanged=len(read(analysis/'PROTECTED_OUTPUTS.json')),
        new_single_data_axis_figures=96,PNG_visual_review_passed=True,PDF_separately_rendered=False,
        incoming_scientific_members_unchanged=True,local_links_checked=links,asset_map_rows=len(assets),
        independent_decimal_lab_sum_max_error_J=str(worst),source_hashes_verified=True,
        exact_old_record_time_arrays_verified=True,original_coordinates_restored=True,
        LW_reference_is_prescribed_orbit_not_exact_new_worldlines=True,
        synchrotron_new_results_pending=False,no_physical_rerun=True,actual_PPT_created=False))
    allowed={'.md','.txt','.csv','.json','.png','.pdf','.mp4','.npz'}
    files=sorted(p for p in root.rglob('*') if p.is_file() and p!=root/'SHA256SUMS.txt')
    assert all(p.suffix.lower() in allowed or p.name=='SHA256SUMS' for p in files)
    (root/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(root)}\n' for p in files));files.append(root/'SHA256SUMS.txt')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,str(Path(root.name)/p.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for p in files:assert hashlib.sha256(z.read(str(Path(root.name)/p.relative_to(root)))).hexdigest()==sha(p)
    package=dict(complete=True,path=str(archive),archive=str(archive),directory=str(root),bytes=archive.stat().st_size,
        sha256=sha(archive),files=len(files),source_commit=commit,synchrotron_new_results_pending=False,
        physical_jobs=[19814,19815],guide=str(root/'PPT制作指导与结果总结.md'),ZIP_CRC_and_all_member_SHA256_passed=True)
    save(root/'PACKAGE.json',package);save(ws/'results/LATEST_PPT_DELIVERY.json',package)
    small=ws/'results/synchrotron_two_method_results_20261116_r1.zip'
    small_files=[]
    with zipfile.ZipFile(small,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(current.rglob('*')):
            if p.is_file():
                member=str(Path('synchrotron_two_method_results/06_synchrotron_current')/p.relative_to(current))
                z.write(p,member);small_files.append((p,member))
        theory=root/'01_reviewed_physics/docs/03_LW_VELOCITY_RADIATION_AND_NEAR_FIELD.md'
        member='synchrotron_two_method_results/01_reviewed_physics/docs/03_LW_VELOCITY_RADIATION_AND_NEAR_FIELD.md'
        z.write(theory,member);small_files.append((theory,member))
        z.writestr('synchrotron_two_method_results/SHA256SUMS.txt',''.join(f'{sha(p)}  {n.split("/",1)[1]}\n' for p,n in small_files))
    with zipfile.ZipFile(small) as z:
        assert z.testzip() is None
        for p,n in small_files:assert hashlib.sha256(z.read(n)).hexdigest()==sha(p)
    package['synchrotron_only_archive']=str(small);package['synchrotron_only_bytes']=small.stat().st_size
    package['synchrotron_only_sha256']=sha(small)
    save(root/'PACKAGE.json',package);save(ws/'results/LATEST_PPT_DELIVERY.json',package)
    print(json.dumps(package,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main(*sys.argv[1:])
