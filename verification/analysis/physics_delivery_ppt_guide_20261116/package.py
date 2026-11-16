"""Consolidate reviewed results and formulas with a portable PPT brief."""
import csv
import hashlib
import json
import re
import shutil
import stat
import subprocess
import sys
import zipfile
from pathlib import Path, PurePosixPath


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path,value):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def identity(path):
    s=Path(path).stat()
    return dict(inode=s.st_ino,size=s.st_size,mtime_ns=s.st_mtime_ns)


ALLOWED={'.png','.pdf','.md','.json','.txt','.csv','.npz','.mp4'}


def allowed(path):
    return path.suffix.lower() in ALLOWED or path.name=='SHA256SUMS'


def main(workspace,output,archive):
    ws,root,archive=map(lambda p:Path(p).resolve(),[workspace,output,archive])
    assert not root.exists() and not archive.exists()
    root.mkdir(parents=True)
    template=Path(__file__).resolve().parent
    protected={}
    mappings=[]
    packages=[]

    def protect(path):
        path=Path(path)
        protected[str(path)]=dict(sha256=sha(path),**identity(path))

    def add_file(source,dest):
        source=Path(source)
        assert source.is_file() and not source.is_symlink() and allowed(source)
        protect(source)
        path=root/dest
        path.parent.mkdir(parents=True,exist_ok=True)
        assert not path.exists()
        shutil.copyfile(source,path)
        assert sha(path)==protected[str(source)]['sha256']
        mappings.append(dict(source=str(source),source_member=None,
            destination=str(path.relative_to(root)),sha256=sha(path),bytes=path.stat().st_size))

    specs=[('physics_front_complete_bilingual_reviewed_20261112_r3.zip',
            '35c31ea5cbcca4a008b1e284d2aa4a9cd6f12dd2aaa11148c8f238fd2577580a',
            'physics_front_complete_bilingual','01_reviewed_physics'),
           ('fel_six_energy_and_theory_bilingual_20261116_r1.zip',
            '35bbc131b92b1d237aea58d050ebf7a2d7b6481d0e2f73f0e3cad592c31e9e48',
            'fel_six_energy_and_theory_20261116_r1','02_FEL_six')]
    for filename,expected,oldroot,newroot in specs:
        source=ws/'results'/filename
        protect(source)
        assert protected[str(source)]['sha256']==expected,filename
        with zipfile.ZipFile(source) as z:
            assert z.testzip() is None
            ns=[n for n in z.namelist() if not n.endswith('/')]
            assert len(ns)==len(set(ns))
            checked=0
            # Only top-level current manifests; nested ones are historical records.
            for member in ns:
                p=PurePosixPath(member)
                if len(p.parts)==2 and 'SHA256SUM' in p.name:
                    for line in z.read(member).decode().splitlines():
                        m=re.fullmatch(r'([0-9a-f]{64})  (.+)',line)
                        assert m,(member,line)
                        target=str(PurePosixPath(oldroot)/m[2])
                        assert target in ns,target
                        assert hashlib.sha256(z.read(target)).hexdigest()==m[1],target
                        checked+=1
            assert checked>0
            for member in ns:
                p=PurePosixPath(member)
                assert not p.is_absolute() and '..' not in p.parts and p.parts[0]==oldroot
                assert len(p.parts)>1
                assert not stat.S_ISLNK(z.getinfo(member).external_attr>>16)
                rel=Path(newroot).joinpath(*p.parts[1:])
                assert allowed(rel),member
                dest=root/rel
                assert dest.resolve().is_relative_to(root) and not dest.exists()
                dest.parent.mkdir(parents=True,exist_ok=True)
                h=hashlib.sha256()
                with z.open(member) as incoming,dest.open('xb') as outgoing:
                    for block in iter(lambda:incoming.read(4194304),b''):
                        outgoing.write(block)
                        h.update(block)
                assert dest.stat().st_size==z.getinfo(member).file_size
                mappings.append(dict(source=str(source),source_member=member,
                    destination=str(rel),sha256=h.hexdigest(),bytes=dest.stat().st_size))
        packages.append(dict(archive=str(source),sha256=expected,member_count=len(ns),
            original_top_level_manifest_checks=checked,destination=newroot))

    eleven=ws/'results/fel_independent_w355_20261114_r1/eleven-single-node-r1/saturation-20261116-r2'
    for name in ['README.md','SATURATION.json','VALIDATION.json','SOURCE_REVIEW.json','JOB_STATUS.txt']:
        add_file(eleven/name,Path('03_FEL_eleven')/name)
    for name in ['RESULT.json','gain.csv','native_history.csv','NATIVE_PEAK_CANDIDATES.json',
                 'power_growth.png','power_growth.pdf','bunching_growth.png','bunching_growth.pdf',
                 'saturation.png','saturation.pdf']:
        add_file(eleven/'data'/name,Path('03_FEL_eleven/data')/name)
    formulas=ws/'results/fel_theory_formulas_20261116_r1/THEORY.md'
    assert sha(formulas)=='45b8b15de6aca71718c3915f35057fe0c5c2de8dd320439da479017b235cf975'
    add_file(formulas,Path('04_formula_derivations/FEL_THEORY.md'))
    add_file(template/'PPT_GUIDE.md',Path('PPT制作指导与结果总结.md'))
    synch=ws/'results/synchrotron_full_ledger_20261116_r1'
    for name in ['README.md','JOB_LINKS.json','SUBMISSION.json','STATIC_VALIDATION.json',
                 'REGRESSION.json','SOURCE_COMMIT.json','POSTPROCESSING_PLAN.md']:
        add_file(synch/name,Path('05_synchrotron_rerun_status')/name)
    if (synch/'PREFLIGHT_COMPLETE.json').is_file():
        add_file(synch/'PREFLIGHT_COMPLETE.json',Path('05_synchrotron_rerun_status/PREFLIGHT_COMPLETE.json'))
    jobs=read(synch/'JOB_LINKS.json')
    assert (jobs['preflight'],jobs['three_point'],jobs['longitudinal_momentum_matching'])==(19813,19814,19815)

    assets=[]

    def asset(module,title,primary,data,note,english=''):
        for path in [primary,data,english]:
            if path:
                assert (root/path).is_file(),path
        assets.append(dict(module=module,title=title,preferred_asset=primary,
            numerical_or_derivation_source=data,english_asset=english,PPT_use=note))

    front='01_reviewed_physics'
    f6='02_FEL_six'
    asset('真空','实际传播',f'{front}/zh-CN/01_vacuum_and_short_drift/01_vacuum_propagation.png',
          f'{front}/01_vacuum_and_short_drift/vacuum_energy_and_error.csv',
          '来源整图；正式页用实际传播帧，理论误差单轴另页',
          f'{front}/01_vacuum_and_short_drift/01_vacuum_propagation.png')
    asset('真空','球壳传播与吸收视频',f'{front}/01_vacuum_and_short_drift/vacuum_shell_absorption.mp4',
          f'{front}/01_vacuum_and_short_drift/vacuum_time_fields.npz','播放视频或取带真实时间的场帧')
    asset('吸收层','吸收与反射视频',f'{front}/01_vacuum_and_short_drift/cpml_vs_reflecting_wall.mp4',
          f'{front}/01_vacuum_and_short_drift/cpml_probe.csv','同一入射；反射率用探针窗实测')
    asset('吸收层','反射对照',f'{front}/zh-CN/01_vacuum_and_short_drift/03_cpml_absorption_control.png',
          f'{front}/01_vacuum_and_short_drift/cpml_probe.csv','旧整图按探针时序/反射统计拆开',
          f'{front}/01_vacuum_and_short_drift/03_cpml_absorption_control.png')
    for shape in ['linear-binomial3','quintic-compensated']:
        for frame,figure in [('boost','04_free_drift_boost_energy'),('lab','05_free_drift_lab_energy_components')]:
            filename=f'boost_{shape}.csv' if frame=='boost' else f'lab_complete_{shape}.csv'
            asset('短漂移',f'{shape} {frame} 完整账本',
                  f'{front}/zh-CN/01_vacuum_and_short_drift/{figure}.png',
                  f'{front}/01_vacuum_and_short_drift/{filename}',
                  '从CSV取该方法分项与SUM；时间/位置横轴分开，侧面与端口分开',
                  f'{front}/01_vacuum_and_short_drift/{figure}.png')
    for case in ['uniform','converging','gaussian']:
        for method in ['linear-binomial3','quintic-compensated']:
            for name,table_name in [('space_charge_theory_comparison','theory_comparison'),
                                    ('boost_energy_ledger','boost_energy_budget'),('lab_energy_ledger','lab_energy_budget')]:
                asset('球形空间电荷',f'{case} {method} {name}',
                      f'{front}/zh-CN/02_space_charge/{case}/{name}.png',
                      f'{front}/02_space_charge/{case}/{method}/{table_name}.csv',
                      '旧整图包含多分量/方法，正式PPT按单物理量单轴制作',
                      f'{front}/02_space_charge/{case}/{name}.png')
    for name,title in [('correct_100k_sample_matched','同粒子LW频谱'),
                       ('correct_transverse_bend_single_frequency_8eV','全平面总场/完整LW'),
                       ('correct_transverse_bend_aperture_8eV','离束孔径/加速度LW')]:
        asset('有限弯转LW',title,f'{front}/zh-CN/04_synchrotron/accepted_LW_records/{name}.png',
              f'{front}/04_synchrotron/accepted_LW_records/RESULT_REVIEW.md',
              '历史附件；不能用作19814/19815的自洽两方法结果，主展示待新结果闭合后替换',
              f'{front}/04_synchrotron/accepted_LW_records/{name}.png')
    asset('完整束团','ELEGANT内场与注入理论',f'{front}/docs/02_FULL_BUNCH_SELF_FIELD_AND_ELEGANT.md',
          f'{front}/docs/02_FULL_BUNCH_SELF_FIELD_AND_ELEGANT.md',
          '只有理论/统计定义；未收录独立完整PIC/ELEGANT数值交付')
    for name,title in [('power','六段功率'),('bunching','六段局部RMS BF'),('particle_energy_loss','六段平均能损')]:
        asset('六段FEL',title,f'{f6}/02_growth_comparison/figures/zh-CN/{name}.png',
              f'{f6}/02_growth_comparison/data/comparison_by_segment.csv','已是独立单轴图，可直接使用',
              f'{f6}/02_growth_comparison/figures/en/{name}.png')
    for tag in ['F3','MM']:
        for name,title,data_name in [(f'{tag}_boost_balance','计算系完整分项',f'{tag}_boost_components.csv'),
                                     (f'{tag}_fixed_lab_balance','lab完整分项',f'{tag}_lab_position_curves.csv')]:
            asset('六段账本',f'{tag} {title}',f'{f6}/01_energy_curves/figures/zh-CN/{name}.png',
                  f'{f6}/01_energy_curves/data/{data_name}','已是独立单轴图；非零SUM保留',
                  f'{f6}/01_energy_curves/figures/en/{name}.png')
    for name,title in [('boost_discrete_sum','计算系离散加和'),('fixed_lab_sum','固定lab完整加和'),
                       ('energy_pairing','电功配对差'),('longitudinal_current_pairing','纵向配对差')]:
        asset('六段账本与定位',title,f'{f6}/01_energy_curves/figures/zh-CN/{name}.png',
              f'{f6}/01_energy_curves/README.md','两方法同物理量单轴对照；参考系不同的残差不能相减归因',
              f'{f6}/01_energy_curves/figures/en/{name}.png')
    asset('谢明与空间电荷','首选完整公式推导','04_formula_derivations/FEL_THEORY.md',
          '04_formula_derivations/FEL_THEORY.md','主讲用公式与假设；不增加无助于推导的理论附图')
    asset('弯铁新两方法','自洽重跑与完整账本更新计划','05_synchrotron_rerun_status/README.md',
          '05_synchrotron_rerun_status/JOB_LINKS.json','新主结果待闭合；旧七点图仅作历史，不预填新数据')
    for name,title in [('power_growth','十一段功率'),('bunching_growth','十一段BF')]:
        asset('十一段Puffin',title,f'03_FEL_eleven/data/{name}.png',
              '03_FEL_eleven/data/gain.csv','现有独立图为英文轴标；中文页题，重画时保持CSV和窗口定义',
              f'03_FEL_eleven/data/{name}.png')
    with (root/'PPT_ASSET_MAP.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(assets[0]))
        w.writeheader();w.writerows(assets)

    (root/'README.md').write_text(
        '# 物理结果与 PPT 制作指导总交付\n\n'
        '先读 [PPT 制作指导与结果总结](PPT制作指导与结果总结.md)，再按 '
        '[素材映射表](PPT_ASSET_MAP.csv) 选择图、数据和推导。'
        '本包用于后续制作 PPT，没有生成实际幻灯片。\n\n'
        '- [已复核前半部分：真空、短漂移、三球、已有 LW 同步辐射](01_reviewed_physics/README.md)\n'
        '- [独立随机六段 FEL：三方法增长和两方法能量账本](02_FEL_six/README.md)\n'
        '- [独立随机十一段 Puffin：饱和位置与功率/BF](03_FEL_eleven/README.md)\n'
        '- [首选纯公式推导：谢明与空间电荷恢复项](04_formula_derivations/FEL_THEORY.md)\n\n'
        '- [弯铁新两方法自洽重跑记录与更新计划](05_synchrotron_rerun_status/README.md)\n\n'
        '包含中英文图片、PDF、视频、CSV/JSON、用于传播视频的时序 NPZ，以及公式和结果记录。'
        '没有源代码、可执行程序或科学 H5。原交付与原数据不改写；三类已有理论附图仍在素材中，'
        'PPT 主讲按用户要求优先用公式。\n\n'
        '前半部分图文为 2026-11-12 的已验快照，六段/十一段及新推导为 2026-11-16 的交付。'
        '子包旧报告中“19731继续执行”等是原日期的历史描述，不是本轮队列查询。'
        '完整真实束团漂移/PIC–ELEGANT 仅收录理论与统计口径，没有将其完整数值对照标成已交付。'
        '旧报告可能引用原工作区或包外历史附件，现行主导航、制作指导与素材表使用包内可用路径。\n\n'
        '弯铁新三点核/五次补偿自洽任务19814/19815已提交，依赖预检19813通过；'
        '当前主结果待闭合复核，不预填谱/LW/能量收支数值，旧七点无反馈图只作历史附件。'
        '关键边界：六段动量匹配实际是五次补偿，'
        '固定 lab 管加和仍为 26.51% / 37.56%；BF 抑制仅有量级理论支持，未唯一归因。'
        '具体数值、符号、统计窗和讲述方式见制作指导。\n')

    # Verify key statements against exactly the records shipped here.
    six=read(root/'02_FEL_six/01_energy_curves/RESULT.json')['methods']
    for tag,expected in [('F3',26.5066555479),('MM',37.563021907)]:
        assert abs(six[tag]['final_fixed_lab_sum_percent_fixed_D']-expected)<1e-7
    endpoint=rows(root/'02_FEL_six/02_growth_comparison/data/comparison_by_segment.csv')[-1]
    assert abs(float(endpoint['Puffin_power_W'])/1e6-8.145483287762678)<1e-12
    sphere=rows(root/'01_reviewed_physics/SPHERE_SUMMARY.csv')
    assert len(sphere)==6 and all(abs(float(r['sigma_x_error_percent']))<.271 for r in sphere)
    sat=read(root/'03_FEL_eleven/SATURATION.json')
    assert abs(sat['first_prominent_native_peak_lab_z_m']-19.358)<1e-9
    guide=(root/'PPT制作指导与结果总结.md').read_text()
    assert '不使用此前' in guide and '5.17%' in guide and '5.16%' in guide
    assert '七点核' in guide and '未收录' in guide

    good,bad_legacy=[],[]
    for path in sorted(root.rglob('*.md')):
        strict=path.parent==root
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',path.read_text()):
            target=target.split('#',1)[0]
            if not target or target.startswith(('http://','https://','mailto:')):
                continue
            dest=(path.parent/target).resolve()
            record=dict(document=str(path.relative_to(root)),target=target)
            if dest.is_relative_to(root) and dest.exists():
                good.append(record)
            else:
                assert not strict,record
                bad_legacy.append(record)
    assert len(good)>50
    for entry in mappings:
        assert sha(root/entry['destination'])==entry['sha256']
    for path,before in protected.items():
        assert identity(path)=={k:before[k] for k in ['inode','size','mtime_ns']},path
        assert sha(path)==before['sha256'],path

    commit=subprocess.check_output(['git','-C',str(template.parents[2]),'rev-parse','HEAD'],text=True).strip()
    save(root/'metadata/SOURCE_COMMIT.json',dict(commit=commit,pushed=False,
        templates_sha256={p.name:sha(p) for p in template.iterdir() if p.is_file()}))
    save(root/'metadata/SOURCE_PACKAGES.json',packages)
    save(root/'metadata/SOURCE_FILE_IDENTITIES_AND_SHA256.json',protected)
    save(root/'metadata/INCOMING_MEMBER_MAP.json',mappings)
    save(root/'metadata/LEGACY_DOCUMENT_LINKS.json',dict(
        explanation='Unchanged historical source references; use current root navigation and PPT asset map.',
        missing_or_external_local_references=bad_legacy))
    save(root/'metadata/VALIDATION.json',dict(passed=True,incoming_members_byte_preserved=len(mappings),
        input_archives_and_source_files_unchanged=len(protected),asset_map_rows=len(assets),
        asset_paths_exist=True,current_primary_markdown_links_valid=True,
        valid_local_markdown_links=len(good),legacy_report_local_reference_count=len(bad_legacy),
        key_numbers_checked_against_shipped_results=True,formulas_included=True,
        no_new_physics=True,no_scientific_H5_read=True,no_new_figures=True,
        synchrotron_new_results_pending=True,synchrotron_rerun_jobs=[19813,19814,19815],
        actual_PPT_created=False,content_only_pictures_videos_records_and_documents=True))
    allfiles=sorted(p for p in root.rglob('*') if p.is_file())
    assert all(allowed(p) for p in allfiles)
    assert not any(p.is_symlink() for p in allfiles)
    (root/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(root)}\n' for p in allfiles))
    allfiles.append(root/'SHA256SUMS.txt')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in allfiles:
            z.write(path,str(Path(root.name)/path.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(allfiles)
        for path in allfiles:
            assert hashlib.sha256(z.read(str(Path(root.name)/path.relative_to(root)))).hexdigest()==sha(path)
    package=dict(complete=True,path=str(archive),archive=str(archive),directory=str(root),
        guide=str(root/'PPT制作指导与结果总结.md'),bytes=archive.stat().st_size,
        sha256=sha(archive),files=len(allfiles),source_commit=commit,
        scope='Reviewed vacuum/short drift/spheres + historical LW + independent six-segment FEL + eleven-segment Puffin + formulas + PPT guidance + submitted synchrotron rerun status; new bend results pending; full-bunch ELEGANT numeric delivery not included.',
        synchrotron_new_results_pending=True,synchrotron_rerun_jobs=[19813,19814,19815],
        ZIP_CRC_and_all_member_SHA256_passed=True,actual_PPT_created=False)
    save(root/'PACKAGE.json',package)
    save(ws/'results/LATEST_PPT_DELIVERY.json',package)
    print(json.dumps(package,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main(*sys.argv[1:])
