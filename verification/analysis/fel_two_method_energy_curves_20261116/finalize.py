"""Validate and package only bilingual figures and result records."""
import csv
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path
import numpy as np

LD=np.longdouble


def read(p):return json.loads(Path(p).read_text())


def save(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''):h.update(block)
    return h.hexdigest()


def main(root,archive):
    root,archive=Path(root).resolve(),Path(archive).resolve()
    result=read(root/'RESULT.json');methods=result['methods']
    plots=read(root/'PLOTS.json')['figures']
    assert len(plots)==20 and all(p['axes']==1 for p in plots)
    protected={};source_files=read(root/'SOURCE_DATA_SHA256.json')
    for path,expected in source_files.items():assert sha(Path(path))==expected,path
    for source in [root/'three-point',root.parent/'momentum-ledger-20261116-r1']:
        for path,before in read(source/'PROTECTED_OUTPUTS.json').items():
            s=Path(path).stat()
            assert before==dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino),path
            protected[path]=before
    for path,expected in read(root/'three-point/FROZEN_ANALYSIS.json').items():assert sha(Path(path))==expected,path
    frozen={}
    for p in [root.parent/'six/FROZEN_PHYSICS.json',root.parent/'six-aprl/FROZEN.json']:frozen.update(read(p))
    for path,expected in frozen.items():assert sha(Path(path))==expected,path
    checks={}
    for tag in ['F3','MM']:
        with (root/'data'/(tag+'_time_curves.csv')).open() as f:
            rows=list(csv.DictReader(f))
        assert len(rows)==28200
        split=max(abs(LD(r['energy_pairing_J'])+LD(r['longitudinal_current_pairing_J'])-LD(r['lab_on_boost_slice_sum_J'])) for r in rows)
        assert split<LD('2e-18')
        with (root/'data'/(tag+'_lab_position_curves.csv')).open() as f:
            lab=list(csv.DictReader(f))
        worst=LD(0)
        for row in lab:
            if not row['physical_lab_sum_J']:continue
            total=sum(LD(row[k]) for k in ['particle_net_transport_J','particle_time_caps_J',
                'field_ports_net_J','field_transverse_net_J','field_time_caps_net_J'])
            worst=max(worst,abs(total-LD(row['physical_lab_sum_J'])))
        assert worst<LD('1e-18')
        checks[tag]=dict(time_records=len(rows),lab_planes=len(lab),split_max_J=float(split),fixed_lab_sum_max_J=float(worst))
    # Preserve the earlier momentum case numerical values and geometry.
    old=read(root.parent/'momentum-ledger-20261116-r1/RESULT.json')
    new=read(root/'three-point/RESULT.json')
    for key in ['common_relative_tau_s','tau_reference_s','observer_half_width_m']:
        assert new['fixed_lab_worldtube'][key]==old['fixed_lab_worldtube'][key],key
    assert methods['MM']['final_lab_on_boost_slice_sum_J']==old['lab_on_boost_slices']['final_components_J']['full_window_Lorentz_accounting_difference_J']
    assert methods['MM']['final_fixed_lab_sum_J']==old['fixed_lab_worldtube']['last_complete_residual_J']
    specs={'energy_pairing':('energy_pairing_J',1e9),
        'longitudinal_current_pairing':('longitudinal_current_pairing_J',1e6),
        'longitudinal_force_pairing':('longitudinal_force_pairing_J',1e6),
        'lab_slice_sum':('lab_on_boost_slice_sum_J',1e6),
        'boost_discrete_sum':('boost_discrete_sum_J',1e18)}
    for tag_index,tag in enumerate(['F3','MM']):
        saved=np.genfromtxt(root/'data'/(tag+'_time_curves.csv'),delimiter=',',names=True,dtype=LD)
        for language in ['zh-CN','en']:
            for name,(key,scale) in specs.items():
                p=next(p for p in plots if p['name']==name and p['language']==language)
                assert np.array_equal(p['series'][tag_index]['y'],np.asarray(saved[key]*scale,float)),(tag,name)
                assert np.array_equal(p['series'][tag_index]['x'],np.asarray(saved['time_box_s'],float)*1e12),(tag,name)
    fields=[('final_energy_pairing_J',1e9,'电功配对差 nJ'),
        ('final_longitudinal_current_pairing_J',1e6,'粒子/电流纵向配对差 µJ'),
        ('final_longitudinal_force_pairing_J',1e6,'计入原生电荷电力后的纵向差 µJ'),
        ('final_lab_on_boost_slice_sum_J',1e6,'计算系等时面上的lab账本加和 µJ'),
        ('final_boost_discrete_sum_J',1e18,'计算系完整离散账本加和 aJ'),
        ('final_fixed_lab_particle_net_J',1e6,'固定lab管粒子净能变 µJ'),
        ('final_fixed_lab_sum_J',1e6,'固定lab管物理收支加和 µJ'),
        ('final_fixed_lab_sum_percent_fixed_D',1,'固定lab加和 / 同方法固定D %')]
    lines=['# 三点核与动量匹配：能量账本变化曲线','',
        '2026-11-16。同一独立w355六段束流：三点核19796_0、动量匹配19796_1。每张图都是独立单轴图，叠加两种方法；保留中英文PNG/PDF与完整CSV。动量匹配沿用用户简称，实际为quintic-compensated高阶补偿方法，不能据此认为严格总纵向动量守恒已经成立。','',
        '## 末值对照','', '| 量 | 三点核 | 动量匹配（五次补偿） |','|---|---:|---:|']
    for key,scale,label in fields:
        values=['未知' if methods[tag][key] is None else f'{methods[tag][key]*scale:+.9g}' for tag in ['F3','MM']]
        lines.append(f'| {label} | {values[0]} | {values[1]} |')
    lines+=['','两种方法各自使用整个记录中完整lab管粒子净输运绝对值的最大值作为固定D；曲线均按绝对单位直接比较，没有单独缩放曲线去匹配幅度，也不按每个时刻的粒子能变做除法。','',
        '## 直接比较的变化曲线','',
        '### 电功配对差','', '![电功配对差](figures/zh-CN/energy_pairing.png)','',
        '### 纵向冲量配对差','', '![纵向冲量配对差](figures/zh-CN/longitudinal_current_pairing.png)','',
        '### 计入原生电荷电力后的纵向差','',
        '![原生电荷电力后的纵向差](figures/zh-CN/longitudinal_force_pairing.png)','',
        '### 计算系等时面上的实验室系账本加和','',
        '![lab组合加和](figures/zh-CN/lab_slice_sum.png)','',
        '### 计算系完整离散账本加和','',
        '![计算系离散加和](figures/zh-CN/boost_discrete_sum.png)','',
        '### 固定实验室管的物理收支加和','',
        '![固定lab管加和](figures/zh-CN/fixed_lab_sum.png)','',
        '## 各方法的分项收支单图','',
        '- [三点核：计算系收支](figures/zh-CN/F3_boost_balance.png)',
        '- [动量匹配：计算系收支](figures/zh-CN/MM_boost_balance.png)',
        '- [三点核：固定lab管收支](figures/zh-CN/F3_fixed_lab_balance.png)',
        '- [动量匹配：固定lab管收支](figures/zh-CN/MM_fixed_lab_balance.png)','',
        '固定lab管中入口取负号、出口取正号，四侧能流与通过端口分开，时间端面存量单独保留；完整加和包括粒子端口、侧流、粒子时间端盖以及场的两端、四侧和两时间端盖。若粒子端盖不可确定，对应物理加和留空，不填零。','',
        '## 统计与公式','',
        '计算系相关曲线横轴为真实模拟时间（ps），全部28200个记录点保留。固定实验室管曲线横轴为7个实际探测面位置（入口加6段出口）；点之间的线只连接已有采样，不表示连续空间采样。','',
        '两种方法使用相同输入、波荡器和网格。固定管为z=0.03m至各段出口，横向±840µm，完整公共迟滞时间门和全频总场；不加Hann、不去掉近场或启动/尾场。它与功率/BF对照采用的后半高流区窗口、±520µm光学孔径是不同统计。','',
        '图中的电功与纵向差按此前审计口径定义：','',
        '    E差 = Gamma × [粒子DeltaKprime − 外磁场boost功 − 累计WJ]',
        '    P电流差 = Gamma beta × [粒子Delta(cPzprime) + 外磁场boost功/beta + 累计PJ]',
        '    P力差 = P电流差 − Gamma beta × 累计GaussE','',
        'E差与P电流差的加和等于计算系等时面上的lab配对差；完整lab离散账本还保留纯场算子恒等式舍入余差，因此这里直接沿用完整账本记录作图。P力差另计入原生电荷电力，不能把E差加P力差当成同一个完整lab加和。','',
        '外磁场冲量由静态lab磁场不做功的关系作代数消元，并非新测得的逐粒子力。当前没有完整实际子步Ez、v叉B、外磁冲量轨迹，不能把这些全局曲线精确拆分为某条推进语句的误差。','',
        '固定lab管加和与计算系等时面的lab组合使用不同的时空域和场定义，不能相减来编造误差分摊。计算系完整离散能量账本保留时间配对、Cowan、边界和吸收层通道，lab能量—动量组合另有Gauss项；计算系离散能量加和接近零不等于连续场的lab物理收支已经闭合。','',
        '## 数据与复核','',
        '- [三点核：逐时刻变化数据](data/F3_time_curves.csv)',
        '- [动量匹配：逐时刻变化数据](data/MM_time_curves.csv)',
        '- [三点核：逐lab面全部分项](data/F3_lab_position_curves.csv)',
        '- [动量匹配：逐lab面全部分项](data/MM_lab_position_curves.csv)',
        '- [三点核：计算系完整分项](data/F3_boost_components.csv)',
        '- [动量匹配：计算系完整分项](data/MM_boost_components.csv)',
        '- [三点核：lab能量动量组合分项](data/F3_lab_slice_components.csv)',
        '- [动量匹配：lab能量动量组合分项](data/MM_lab_slice_components.csv)',
        '- [完整数值摘要](RESULT.json)',
        '- [复核记录](VALIDATION.json)','',
        '同名英文图位于figures/en。动量匹配使用之前已经核对的闭合缓存；三点核只读后处理19812补算，没有新的物理仿真。原输入、程序、卡、输出与此前结果保持。']
    (root/'README.md').write_text('\n'.join(lines)+'\n')
    save(root/'VALIDATION.json',dict(passed=True,independent_methods=2,time_series_records_each=28200,
        single_axis_figures=20,PNG_visual_review_passed=True,PDF_separately_rendered=False,
        preserved_original_momentum_final_values=True,source_data_SHA256_verified=len(source_files),
        original_output_identities_unchanged=len(protected),production_frozen_SHA256_verified=len(frozen),
        source_F3_frozen_scripts_unchanged=True,numerical_checks=checks,no_new_physics=True))
    for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',(root/'README.md').read_text()):assert (root/target).exists(),target
    commit=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[3]),'rev-parse','HEAD'],text=True).strip()
    save(root/'SOURCE_COMMIT.json',dict(commit=commit,pushed=False,
        templates_sha256={p.name:sha(p) for p in Path(__file__).parent.iterdir() if p.is_file()}))
    files=sorted([p for p in root.rglob('*') if p.is_file() and p.relative_to(root).parts[0] in ['data','figures']]+
        [root/n for n in ['README.md','RESULT.json','VALIDATION.json','SOURCE_COMMIT.json','SOURCE_DATA_SHA256.json','PLOTS.json']])
    assert not any(p.suffix in ['.py','.h5','.npz','.bin','.cpp','.sbatch'] for p in files)
    (root/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(root)}\n' for p in files))
    files.append(root/'SHA256SUMS')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,str(Path(root.name)/p.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for p in files:assert hashlib.sha256(z.read(str(Path(root.name)/p.relative_to(root)))).hexdigest()==sha(p)
    package=dict(complete=True,directory=str(root),archive=str(archive),bytes=archive.stat().st_size,
        sha256=sha(archive),source_commit=commit,archive_CRC_and_member_SHA256_passed=True)
    save(root/'PACKAGE.json',package)
    save(archive.parent/'LATEST_SIX_ENERGY_CURVES.json',package)
    print(json.dumps(package,indent=2),flush=True)


if __name__=='__main__':main(*sys.argv[1:])
