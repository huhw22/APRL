#!/usr/bin/env python3
"""Closed-output review; physics comparisons and program checks stay separate."""
import argparse,hashlib,json,re
from pathlib import Path
import h5py,numpy as np

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source_work(d):
    assert read(d/'COMPLETION.json')['checks_pass']
    with h5py.File(d/'ledger/energy-ledger.h5') as h:
        g=next(iter(h));assert int(h[g+'/complete'][()])==1
        r=h[g+'/records'][:];ledger=np.longdouble(r['boundary_source_work_J'][-1]);steps=int(r['step'][-1])
    values=[]
    for p in sorted((d/'reference').glob('transverse-source-work-rank-*.csv')):
        a=np.genfromtxt(p,names=True,delimiter=',',dtype=np.longdouble)
        assert len(a)==steps and np.isfinite(a['work_J']).all()
        assert abs(np.sum(a['work_J'])-a['cumulative_work_J'][-1])<1e-24
        values.append(a['cumulative_work_J'][-1])
    total=sum(values,np.longdouble(0));assert len(values)==16 and abs(total-ledger)<1e-24
    return dict(steps=steps,source_work_J=float(total),ledger_source_work_J=float(ledger),difference_J=float(total-ledger))
def main(root):
    plan=read(root/'plan.json');pilot=read(root/'surface_500/momentum-analysis/ANALYSIS.json')
    cases={};work={'surface_500':source_work(root/'surface_500')}
    for k in ('uniform','converging','gaussian'):
        d=root/(k+'_full');assert read(d/'COMPLETION.json')['checks_pass']
        new=read(d/'analysis/RESULT.json');old=read(Path(plan['cases'][k+'_full']['original'])/'analysis/RESULT.json')
        n,o=new['planes'][-1],old['planes'][-1];t=n['reference_mean_lab_gain_eV']
        cases[k]=dict(original_lab_mean_gain_keV=o['PIC_mean_lab_gain_eV']/1e3,
            new_lab_mean_gain_keV=n['PIC_mean_lab_gain_eV']/1e3,theory_lab_mean_gain_keV=t/1e3,
            original_lab_error_keV=(o['PIC_mean_lab_gain_eV']-t)/1e3,new_lab_error_keV=(n['PIC_mean_lab_gain_eV']-t)/1e3,
            new_rest_mean_gain_eV=n['PIC_mean_rest_gain_eV'],theory_rest_mean_gain_eV=n['reference_mean_rest_gain_eV'],
            new_rest_std_gain_eV=n['PIC_std_rest_gain_eV'],theory_rest_std_gain_eV=n['reference_std_rest_gain_eV'],
            new_sigma_x_relative_error=n['PIC_sigma_x_m']/n['reference_sigma_x_m']-1.,
            new_sigma_z_rest_relative_error=n['PIC_sigma_z_rest_m']/n['reference_sigma_z_rest_m']-1.,
            common_boost_slice_final_clean=read(d/'momentum-analysis/ANALYSIS.json')['final_clean_comparison'],
            ledger=read(d/'LEDGER_CHECK.json'))
        work[k]=source_work(d)
    frozen=read(root/'FROZEN.json');assert all(sha(p)==h for p,h in frozen.items())
    assert read(root/'preflight/COMPLETION.json')['checks_pass']
    assert read(root/'two_nodes_200/COMPLETION.json')['checks_pass']
    result=dict(status='CLOSED_SMALL_SPHERE_BOUNDARY_COMPARISON',production_spectrum_accepted=False,
        pilot=pilot['final_clean_comparison'],cases=cases,source_work=work,default_unchanged=True,two_node_scientific_prefix_bitwise=True,
        limitations=['Finite Dirichlet reference with prescribed numerical boundary current; not exact infinite-space fields.',
          'Cowan-z, no current filter, transverse-only CPML, same 16 MPI x8 tested.',
          'Rest energy spread, force RMS and finite-mesh discrepancies remain measured quantities; mean correction does not prove arbitrary precision.',
          'No production ELEGANT bunch, undulator spectrum or general oblique-wave acceptance.'])
    (root/'FINAL.json').write_text(json.dumps(result,indent=2)+'\n')
    rows=['| 束团 | 理论 lab 平均ΔK (keV) | 原结果 | 边界匹配后 | 新偏差 (keV) |','|---|---:|---:|---:|---:|']
    rest=['| 束团 | rest 平均ΔK：PIC / 理论 (eV) | rest std(ΔK)：PIC / 理论 (eV) | 横向RMS尺寸偏差 |','|---|---:|---:|---:|']
    for k,v in cases.items():
        rows.append(f"| {k} | {v['theory_lab_mean_gain_keV']:.3f} | {v['original_lab_mean_gain_keV']:.3f} | {v['new_lab_mean_gain_keV']:.3f} | {v['new_lab_error_keV']:+.3f} |")
        rest.append(f"| {k} | {v['new_rest_mean_gain_eV']:.3f} / {v['theory_rest_mean_gain_eV']:.3f} | {v['new_rest_std_gain_eV']:.3f} / {v['theory_rest_std_gain_eV']:.3f} | {100*v['new_sigma_x_relative_error']:+.4f}% |")
    text='''# 横向参考场与开放边界匹配：短束团物理对照

本轮针对参考场横向截断激发的纵向冲量，采用预设等效表面电流。参考场仍在原CPML内缘终止，
只在四个界面的Jz增加参考B的横向离散curl/mu0；沿z按原参考供场速度平移。源只依赖参考场和时间，
不读取实际总场、不钳制E/B。内部Cowan、Boris、粒子沉积、窗口、分辨率和原CPML更新保持。
该源相当于维护有限参考场边缘的数值边界装置，做功明确计账；不是新增物理电子束。

## 固定lab探测面与独立空间电荷理论

三组原输入事件、32768带权宏粒子/1e8电子、释放位置/动量、全部6个粒子面与理论到达事件保持。
理论为独立相对论球壳ODE，再按每个粒子的固定lab面到达事件变换；没有把固定lab面当作同时rest截面。
下表为入口至各组原末面、逐粒子差后按权重统计，不是初始总束流能量相减。

'''+ '\n'.join(rows)+ '\n\n'+ '\n'.join(rest)+'''

![三类束团的平均能量和膨胀对照](figures/full_space_charge_comparison.png)

## 同一boost时刻的完整束团

只在全部释放、首个粒子尚未流出的区间比较，避免不同取样面的混合。
uniform500步原lab能量偏差−98.630keV，新方案−2.060keV；rest平均uz偏差由约−1.64e−4降至−3.15e−6。
480步粒子所采样平均Ez相对球壳参考偏差由+207826降至−2000V/m；局部场RMS仍约5.49e5V/m，
不能把平均偏差改善误写成每个粒子的场已精确一致。具体晚期数据在FINAL.json与各组momentum-analysis/。

![完整束团共同boost时间对照](figures/pilot_boundary_comparison.png)

## 边界能量交换

等效源的−dt*dV*sum(w*Js*(Ez_old+Ez_new)/2)并入原boundary_source_work_J，另有16份逐rank源做功CSV。
500步源做功−0.410287pJ，与主账本相差约1.36e−29J。全部完整短工况的独立合计见FINAL.json.source_work。
严格离散场算子残差、粒子与网格交换剩余项分别保留；账本闭合不能替代上面的动力学理论对照。

## 适用范围

这是三类短球形束团上的边界修复对照，尚未验收真实大束团、波荡器能谱或任意波入射。
原参考仍为有限Dirichlet/common-beta近似，剩余能散和格点误差不能宣称消失。
只验证Cowan-z、无电流滤波、横向四面CPML配置，不自动启用生产或更改cluster.env/default。
关闭选项无表面文件/循环/通信。新增源是预设加法，所以两个实际场解相减时源精确相消，传播扰动仍由原线性场/CPML算子推进。
既有普通波验证见[原Cowan横向CPML脉冲](../reference_cpml_match_20261030/wave/WAVE_CHECK.json)，
该单一正入射测试不代表普遍波精度证明。

程序与多节点检查、身份、耗时另见[PROGRAM_VALIDATION.md](PROGRAM_VALIDATION.md)，不当作物理精度证明。
源码方法记录在verification/performance/cases/reference_cpml_match/SURFACE_CURRENT.md。
'''
    (root/'README.md').write_text(text)
    program='''# 程序、身份与耗时记录

19051：77required全部通过；关闭功能100步，18个H5科学records与原已验版本逐位相同。
19052：uniform500步，原严格能量账本和只读Boris冲量探针开启。
19053：原核首步无粒子接口检查，界面Ez约0.25V/m；无源物理段最大偏差2.28V/m，原界面114903V/m。
19054/19055/19057：uniform/converging/gaussian原完整短卡，分别1338/1337/1294步。
19056：总MPI/线程固定16x8，两个节点各8rank；200步18个H5科学records与单节点前缀逐位相同。
全部128不同物理核绑定、正常完整关闭，原MPI-free600/60秒watchdog与作业内BTL保持。

单节点500步，原base与新surface均node8：推进335.052894/335.387312秒，即0.670106/0.670775秒每步，差+0.100%。
这是一次含重诊断的配对样本，不能作为稳定性能验收或外推生产成本。launch总时344.043/346.206秒。
表面文件14909552字节；两个32层缓存，开启账本另保存本地表面J/E。活动场网格不增大，逐步新增成本O(Nactive_z*(Nx+Ny))。
准备参考标量仍走原快速生成器；没有加入全空间PIC等待稳定。关闭功能不创建新文件/通信或源循环。

FROZEN.json记录原输入、脚本、卡、binary身份；SOURCE_BUILD.json记录compiled-parent-dirty的实际C++/header SHA。
最终提交与独立binary不以HEAD字符串混称，SOURCE_FINAL.json给出逐文件对应。旧原始结果、A方案未通过记录保持。
新的正入射原CPML测试反射幅值比0.0004611，范围仅限该制造脉冲；附在程序证据，未提升为束团物理验收。
'''
    (root/'PROGRAM_VALIDATION.md').write_text(program)
    print(json.dumps({k:{n:v[n] for n in ('new_lab_error_keV','new_sigma_x_relative_error')} for k,v in cases.items()},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root.resolve())
