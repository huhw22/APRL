"""Version the presentation revision while preserving all scientific records."""
import csv
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from urllib.parse import unquote

from render import save, sha

OLD_NAME='physics_delivery_with_ppt_guide_20261116_r2'
NEW_NAME='physics_delivery_with_ppt_guide_20261117_r1'


def prepare(ws):
    ws=Path(ws).resolve();old=ws/'results'/OLD_NAME;root=ws/'results'/NEW_NAME
    assert not root.exists()
    expected='00cddf17d3d890d2d7beba38f43684334d1e707d2b38a3540d3878d70ecc03ba'
    assert sha(old.with_suffix('.zip'))==expected
    protected={str(p.relative_to(old)):sha(p) for p in sorted(old.rglob('*')) if p.is_file()}
    shutil.copytree(old,root)
    assert all(sha(root/rel)==h for rel,h in protected.items())
    rev=root/'07_readable_figures';rev.mkdir()
    save(rev/'PREVIOUS_PACKAGE.json',dict(directory=str(old),archive_sha256=expected,files=protected))
    # Old manifests remain historical evidence, the new root manifest is created last.
    previous=root/'metadata/before_readability_revision';previous.mkdir()
    for rel in ['README.md','PPT制作指导与结果总结.md','PPT_ASSET_MAP.csv','SHA256SUMS.txt','PACKAGE.json']:
        if (root/rel).exists():shutil.copyfile(root/rel,previous/Path(rel).name)
    save(rev/'PREPARED.json',dict(root=str(root),copied_files=len(protected),no_physics_rerun=True))
    print(root)


def finish(ws):
    ws=Path(ws).resolve();root=ws/'results'/NEW_NAME;rev=root/'07_readable_figures'
    previous=json.loads((rev/'PREVIOUS_PACKAGE.json').read_text());old=Path(previous['directory'])
    assert json.loads((rev/'RENDER_COMPLETE.json').read_text())['complete']
    assert json.loads((rev/'00_vacuum/MOVIE_VALIDATION.json').read_text())['complete']
    assert json.loads((rev/'INDEPENDENT_VALIDATION.json').read_text())['passed']
    assert json.loads((rev/'VISUAL_REVIEW.json').read_text())['passed']
    assert json.loads((rev/'POLISH_VALIDATION.json').read_text())['passed']
    template=Path(__file__).resolve().parent
    shutil.copyfile(template/'能量账本与图像阅读说明.md',rev/'能量账本与图像阅读说明.md')
    shutil.copyfile(template/'README.md',rev/'README.md')
    books=json.loads((rev/'LEDGERS.json').read_text())
    plots=json.loads((rev/'PLOTS.json').read_text())['figures']
    # A dedicated entry table prevents accidental use of old crowded multi-panel figures.
    rows=[];index=['# 清晰版图片与视频入口\n',
        '本目录是当前展示入口。每份账本按“主要分量 / 小分量与加和”成对显示，另附加和独立图。原始 CSV 和公式仍在原模块中，数值没有改动。\n',
        '[先读账本与图像说明](能量账本与图像阅读说明.md)\n',
        'F3 = 能量守恒加三点核；MM = 能量守恒加动量匹配（五次补偿）。\n',
        '快速选图：[短漂移 lab 三点核](01_short_drift/zh-CN/F3_lab_pair.png) · [短漂移 lab 动量匹配](01_short_drift/zh-CN/MM_lab_pair.png) · [均匀球 RMS 与偏差](02_space_charge/uniform/zh-CN/theory_size.png) · [六段 FEL 动量匹配 lab](03_FEL_six/zh-CN/MM_lab_pair.png) · [新弯铁三点核 lab](04_synchrotron/zh-CN/F3_20261111_bend_lab_pair.png)\n',
        '## 真空传播\n',
        '[中文四画面视频](00_vacuum/zh-CN/vacuum_shell_absorption.mp4) · [English](00_vacuum/en/vacuum_shell_absorption.mp4)\n',
        '四个画面分别是 xz 仿真 Ey、xz 理论 Ey、xy 仿真 Ey、场能量时序。也提供四个独立视频。\n']
    for name,label in [('01_xz_simulation','xz 平面：仿真电场 Ey'),('02_xz_theory','xz 平面：理论电场 Ey'),
                       ('03_xy_simulation','xy 平面：仿真电场 Ey'),('04_field_energy','场能量随时间变化'),
                       ('cpml_vs_reflecting_wall','沿 x 的场剖面：吸收层与反射壁')]:
        index.append(f'- [{label}](00_vacuum/zh-CN/{name}.mp4) · [English](00_vacuum/en/{name}.mp4)\n')
        rows.append(dict(module='真空',title=label,preferred_asset=f'07_readable_figures/00_vacuum/zh-CN/{name}.mp4',
                    numerical_or_derivation_source='01_reviewed_physics/01_vacuum_and_short_drift/vacuum_time_fields.npz' if name!='cpml_vs_reflecting_wall' else '01_reviewed_physics/01_vacuum_and_short_drift/cpml_time_fields.npz',
                    english_asset=f'07_readable_figures/00_vacuum/en/{name}.mp4',PPT_use='真实保存时间步；固定色标；保留物理正负'))
    headings={'01_short_drift':'短漂移','02_space_charge/uniform':'冷均匀球','02_space_charge/converging':'会聚均匀球',
              '02_space_charge/gaussian':'冷高斯球','03_FEL_six':'六段 FEL','04_synchrotron':'新弯铁及直线控制'}
    for cat,label in headings.items():
        index.append('\n## '+label+'\n')
        if cat.startswith('02_space_charge'):
            for name,metric in [('size','RMS 尺寸与尺寸偏差'),('mean_rest_gain','静止系平均能量变化与偏差'),
                                ('std_rest_gain','静止系能量变化标准差与偏差'),('mean_lab_gain','实验室系平均能量变化与偏差')]:
                index.append(f'- [{metric}]({cat}/zh-CN/theory_{name}.png) · [English]({cat}/en/theory_{name}.png)\n')
                rows.append(dict(module=label,title=metric,preferred_asset=f'07_readable_figures/{cat}/zh-CN/theory_{name}.png',
                    numerical_or_derivation_source=f'01_reviewed_physics/{cat}/linear-binomial3/theory_comparison.csv',
                    english_asset=f'07_readable_figures/{cat}/en/theory_{name}.png',PPT_use='一套物理理论；上图实际尺寸/能量，下图各自事件匹配的偏差'))
        for b in books:
            if b['category']!=cat:continue
            n=b['name'];p=f'{cat}/zh-CN/{n}';en=f'{cat}/en/{n}'
            index.append(f'- {n}：[完整分项成对图]({p}_pair.png) · [加和独立图]({p}_sum.png) · [English]({en}_pair.png)\n')
            for suffix,usage in [('pair','主要分量与小分量并排；读各轴绝对能量单位'),('sum','单独展示原始非零加和值；与同域同参考系比较')]:
                rows.append(dict(module=label,title=n+' '+suffix,preferred_asset=f'07_readable_figures/{p}_{suffix}.png',
                    numerical_or_derivation_source=b['source'][0],english_asset=f'07_readable_figures/{en}_{suffix}.png',PPT_use=usage))
    (rev/'展示入口.md').write_text('\n'.join(index)+'\n')
    # Record original transfer fields as history; route all energy/theory/video rows to the new figures.
    with (root/'PPT_ASSET_MAP.csv').open() as f:oldrows=list(csv.DictReader(f))
    replace_modules={'短漂移','空间电荷','球形空间电荷'}
    kept=[]
    for row in oldrows:
        path=row['preferred_asset'];title=row['title']
        replaced=['ledger','boost_balance','fixed_lab_balance','lab_balance','boost_internal_balance',
                  '04_free_drift_boost_energy','05_free_drift_lab_energy_components',
                  '06_lab_signed_surfaces','07_lab_sum_comparison','space_charge_theory_comparison']
        if (any(token in path for token in replaced) or path.endswith('.mp4') or row['module'] in replace_modules):continue
        kept.append(row)
    with (root/'PPT_ASSET_MAP.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows+kept)
    readme=(root/'README.md').read_text()
    readme=readme.replace('# 物理结果与 PPT 制作指导 更新交付','# 物理结果与 PPT 制作指导：展示修订版')
    (root/'README.md').write_text(readme.split('\n',1)[0]+'\n\n**展示优先使用 [清晰版图片与视频入口](07_readable_figures/展示入口.md)。**\n\n'
        '[账本逐项解释、xy 弱场原因及共同球壳理论说明](07_readable_figures/能量账本与图像阅读说明.md)。\n\n'+readme.split('\n',1)[1]+
        '\n旧目录的多面板图是上一版可追溯材料；本版 PPT 选图表已改指向 07 的清晰版。此次仅重绘与改说明，没有重新仿真。\n')
    guide=root/'PPT制作指导与结果总结.md';s=guide.read_text()
    override='''
## 当前展示修订（2026-11-17）

制作 PPT 时，账本、球形膨胀理论对照和真空视频统一从 [新版展示入口](07_readable_figures/展示入口.md) 取图，或按更新的 PPT_ASSET_MAP.csv 选图。原章节里的旧图路径保留来源追溯，已由本入口取代。

每种方法、每个参考系单独成页。页内并排放“主要分量”和“小分量与加和”，两侧各自标清绝对能量单位，不逐线归一化；需要强调闭合时另用加和独立图。不要再把全部曲线压在同一个大纵轴上。粒子、纵向过面、侧流、时间端面分别讲清，符号保持原记录；数值通道和外源功不能混称辐射。

球形测试只用一套相对论球壳理论，按两种方法各粒子的实际过面时刻取样；黑线和灰带展示该模型的两组取样差，差值图分别与自己的事件取样比较。RMS 尺寸近乎重合就是接近，用下方 nm 差值图区分，不平移或夸大主曲线。能量变化的标准差不写成投影能散。

视频标题写明“xz 平面：仿真电场 Ey”“xz 平面：理论电场 Ey”“xy 平面：仿真电场 Ey”和“场能量随时间的变化”。展示时说明：波前是球形，电 Hertz 矢量沿 y；xy 图沿 y 的弱场由偶极角分布与 Ey 投影决定，不表示波没有向该方向传播。吸收后仍继续传播的解析图是无限真空参考。

逐项定义和公式见 [账本与图像阅读说明](07_readable_figures/能量账本与图像阅读说明.md)。六段 FEL 的非闭合结论与新弯铁的小残差结论均保持原值，不由展示改版改变。

'''
    guide.write_text(s.split('\n',1)[0]+'\n'+override+s.split('\n',1)[1])
    # Put a current-navigation notice on the affected module documents.
    for rel in ['01_reviewed_physics/README.md','01_reviewed_physics/01_vacuum_and_short_drift/README.md',
                '01_reviewed_physics/02_space_charge/uniform/README.md','01_reviewed_physics/02_space_charge/converging/README.md',
                '01_reviewed_physics/02_space_charge/gaussian/README.md','02_FEL_six/README.md',
                '02_FEL_six/01_energy_curves/README.md','06_synchrotron_current/README.md']:
        path=root/rel
        if path.exists():
            import os
            link=os.path.relpath(rev/'展示入口.md',path.parent)
            text=path.read_text();path.write_text(f'> 当前展示请使用 [2026-11-17 清晰版图片与视频]({link})；下方原图及说明作为上版来源保留。\n\n'+text)
    # Keep source figures intact in their original modules; the original package is also untouched.
    changed=[]
    for rel,expected in previous['files'].items():
        assert sha(old/rel)==expected,('old delivery changed',rel)
        if sha(root/rel)!=expected:changed.append(rel)
    assert all(Path(rel).suffix in ['.md','.csv','.txt'] for rel in changed),changed
    # All CSVs other than the asset-selection index, all NPZs and result JSONs must remain exact.
    scientific=[rel for rel in previous['files'] if Path(rel).suffix in ['.csv','.npz','.json'] and rel not in ['PPT_ASSET_MAP.csv','PACKAGE.json']]
    for rel in scientific:assert sha(root/rel)==previous['files'][rel],rel
    # Verify all generated plot paths and all newly authored relative document links.
    for item in plots:assert (rev/item['path']).is_file()
    broken=[];checked=0
    for path in [root/'README.md',guide,rev/'展示入口.md',rev/'README.md',rev/'能量账本与图像阅读说明.md']:
        for target in re.findall(r'(?<!!)\[[^\]]*\]\(([^\n]+?)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            target=unquote(target.split('#')[0].strip('<>'))
            checked+=1
            if not (path.parent/target).is_file():broken.append([str(path.relative_to(root)),target])
    assert not broken,broken
    for row in rows+kept:
        for key in ['preferred_asset','numerical_or_derivation_source','english_asset']:
            if row[key]:assert (root/row[key]).exists(),(key,row[key])
    source_hashes={p.name:sha(p) for p in template.iterdir() if p.is_file()}
    save(rev/'SOURCE_SCRIPTS.json',source_hashes)
    save(rev/'VALIDATION.json',dict(passed=True,changed_presentation_members_before_root_manifests=changed,
         regenerated_root_manifest_and_package_sidecar=True,scientific_records_unchanged=len(scientific),
         previous_package_all_members_unchanged=True,new_document_links_checked=checked,new_plot_count=len(plots),
         ledgers=len(books),new_physics_jobs=0,source_scripts=source_hashes))
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=template.parents[2],text=True).strip()
    save(rev/'SOURCE_COMMIT.json',dict(commit=commit,source_directory=str(template),source_scripts=source_hashes))
    # PACKAGE.json is a sidecar, not self-hashed.
    (root/'PACKAGE.json').unlink(missing_ok=True)
    sums=root/'SHA256SUMS.txt'
    files=sorted(p for p in root.rglob('*') if p.is_file() and p!=sums)
    sums.write_text(''.join(sha(p)+'  '+str(p.relative_to(root))+'\n' for p in files))
    files.append(sums)
    archive=root.with_suffix('.zip');assert not archive.exists()
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in sorted(files):z.write(p,root.name+'/'+str(p.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for p in files:
            import hashlib
            assert hashlib.sha256(z.read(root.name+'/'+str(p.relative_to(root)))).hexdigest()==sha(p)
    info=dict(complete=True,path=str(archive),archive=str(archive),directory=str(root),bytes=archive.stat().st_size,
              sha256=sha(archive),files=len(files),guide=str(guide),display_index=str(rev/'展示入口.md'),
              no_new_physics=True,source_commit=commit,ZIP_CRC_and_all_member_SHA256_passed=True)
    save(root/'PACKAGE.json',info)
    save(ws/'results/LATEST_PPT_DELIVERY.json',info)
    print(json.dumps(info,ensure_ascii=False,indent=2))


if __name__=='__main__':
    {'prepare':prepare,'finish':finish}[sys.argv[1]](sys.argv[2])
