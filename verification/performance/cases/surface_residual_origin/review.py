#!/usr/bin/env python3
"""Report completed causal controls, preserving all original outputs."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-residual-review-mpl');os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def figure(fig,p):fig.savefig(p.with_suffix('.png'),dpi=180,bbox_inches='tight');fig.savefig(p.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig)
def csv(p):return np.genfromtxt(p,names=True,delimiter=',')
def main(root,follow,probe):
    out=root/'figures';out.mkdir(exist_ok=True);plan=read(root/'plan.json');prior=Path(plan['prior'])
    cases={**{k:root/k for k in ('base','inlet_guard_256','transverse_240')},**{k:follow/k for k in ('inlet_guard_128','inlet_guard_256_long')}}
    stats={}
    sys.path.insert(0,str(root/'scripts'));import report_discrete
    sys.path.insert(0,str(prior/'scripts'));from relative_ledger import construct
    for name,d in cases.items():
        assert read(d/'COMPLETION.json')['checks_pass'];an=read(d/'analysis/ANALYSIS.json')
        rows,sums,_,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');idx,quant,budget=report_discrete.calculate(rows,sums,steps)
        _,_,relative=construct(rows,sums,idx,quant)
        stats[name]=dict(**an['final_clean_comparison'],all_released_step=an['all_released_step'],first_exit_step=an['first_exit_step'],relative_ledger=relative,
            maximum_field_operator_residual_J=budget['max_abs_cumulative_field_operator_residual_J'])
    force=read(root/'FORCE_REVIEW.json')
    fig,ax=plt.subplots(1,2,figsize=(12,4.3))
    allcurves=[('Original 64 cells',prior/'converging_full/momentum-analysis/common_boost_time.csv','#66788d'),
        ('128 cells (diagnostic control)',follow/'inlet_guard_128/analysis/common_boost_time.csv','#c5802d'),
        ('256 cells, continued to 1000',follow/'inlet_guard_256_long/analysis/common_boost_time.csv','#148572'),
        ('Wider xy, same 64-cell inlet',root/'transverse_240/analysis/common_boost_time.csv','#9366b3')]
    for label,p,color in allcurves:
        a=csv(p);ax[0].plot(a['step'],a['lab_energy_error_from_same_box_slice_eV']/1000,label=label,color=color)
    for name,label,color in [('base','64 cells','#66788d'),('inlet_guard_256','256 cells','#148572'),('transverse_240','Wider xy, 64 cells','#9366b3')]:
        a=csv(root/name/'force-review/force.csv');ax[1].plot(a['step'],a['mean_Ez_error_V_per_m']/1000,label=label,color=color)
    ax[0].set_ylabel('Lab mean-energy error on same boost slice (keV)');ax[0].set_xlabel('Field step: full cohort before first exit');ax[0].legend(fontsize=7)
    ax[1].set_ylabel('Mean sampled $E_z$ minus spherical reference (kV/m)');ax[1].set_xlabel('Field step');ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2);a.axhline(0,color='k',lw=.6)
    fig.suptitle('Converging sphere: inlet sensitivity dominates the remaining drift');fig.tight_layout();figure(fig,out/'causal_controls')
    d=probe/'probe_300';assert read(d/'COMPLETION.json')['checks_pass']
    a=np.genfromtxt(d/'diagnostic/inlet-constraint.csv',names=True,delimiter=',',dtype=None,encoding=None)
    begin=a[a['stage']=='begin'];end=a[a['stage']=='electric_complete'];birth=begin[begin['k']==1];ref=csv(root/'reference_gauss.csv')
    retained_end=end[(end['k']<=7)&(end['step']<299)];retained_begin=begin[(begin['k']>=2)&(begin['step']>=1)]
    assert np.array_equal(retained_end['step']+1,retained_begin['step']) and np.array_equal(retained_end['k']+1,retained_begin['k'])
    gauss=dict(continuity_max_C_per_m3=float(max(end['continuity_max_C_per_m3'])),ampere_divergence_max_C_per_m3=float(max(end['ampere_gauss_max_C_per_m3'])),
        ordinary_step_signed_gauss_change_max_C=float(max(abs(begin['gauss_signed_C']-end['gauss_signed_C']))),
        retained_remap_signed_gauss_change_max_C=float(max(abs(retained_end['gauss_signed_C']-retained_begin['gauss_signed_C']))),
        initial_reference_signed_gauss_max_C=float(max(abs(ref['initial_gauss_signed_C']))),
        translated_reference_gauss_max_C=float(max(abs(ref['reference_divE_charge_C']-ref['rigid_rho_charge_C']))),
        birth_signed_gauss_max_C=float(max(abs(birth['gauss_signed_C']))),
        birth_gauss_L1_sum_steps_1_149_C=float(sum(birth['gauss_L1_C'][1:150])),birth_gauss_signed_sum_steps_1_149_C=float(sum(birth['gauss_signed_C'][1:150])),
        scope='Native charge nodes; x/y central +/-120um at 5um spacing, k=1..8; no total-domain net-charge assertion.')
    fig,ax=plt.subplots(1,2,figsize=(12,4.3))
    ax[0].plot(birth['step'][:150],birth['gauss_signed_C'][:150]*1e12,label='New inlet node: signed Gauss defect',color='#b56b4b')
    ax[0].plot(ref['step'],ref['initial_gauss_signed_C']*1e12,label='Initial reference at same world nodes',color='#66788d');ax[0].set_xlabel('Field step');ax[0].set_ylabel(r'Integrated $\epsilon_0\nabla_h\cdot E-\rho$ (pC)');ax[0].legend(fontsize=7)
    ax[1].semilogy(birth['step'],np.maximum(birth['gauss_rms_C_per_m3'],1e-30),label='RMS Gauss defect on new inlet node',color='#b56b4b')
    for key,label,color in [('continuity_max_C_per_m3','Max charge-continuity defect','#148572'),('ampere_gauss_max_C_per_m3','Max Ampere-divergence update defect','#66788d')]:
        vals=np.array([max(end[end['step']==n][key]) for n in birth['step']]);ax[1].semilogy(birth['step'],np.where(vals>0,vals,np.nan),label=label+' (zeros omitted)',color=color)
    ax[1].set_ylabel('Charge-density defect (C/m$^3$)');ax[1].set_xlabel('Field step');ax[1].legend(fontsize=7)
    for z in ax:z.grid(alpha=.2)
    fig.suptitle('Defect enters through newly supplied cells; ordinary updates preserve the constraint');fig.tight_layout();figure(fig,out/'inlet_constraint')
    h=csv(root/'held_inlet_charge.csv');fig,ax=plt.subplots(1,2,figsize=(12,4.2))
    x=h[h['guard']==64];ax[0].plot(x['step'],x['actual_held_charge_C']*1e12,label='Actual individual held velocities',color='#148572');ax[0].plot(x['step'],x['reference_rigid_charge_C']*1e12,label='Rigid reference advection',color='#b56b4b',ls='--')
    ax[0].plot(x['step'],x['common_velocity_exact_CIC_charge_C']*1e12,label='Common velocity, direct CIC',color='#66788d',ls=':',lw=1)
    for guard,color in [(64,'#66788d'),(128,'#c5802d'),(256,'#148572')]:
        x=h[h['guard']==guard];ax[1].plot(x['step'],(x['reference_rigid_charge_C']-x['actual_held_charge_C'])*1e12,label=str(guard)+' cells',color=color)
    ax[0].set_ylabel('Signed charge on new inlet node (pC)');ax[1].set_ylabel('Reference charge minus actual charge (pC)')
    for z in ax:z.set_xlabel('Field step');z.grid(alpha=.2);z.legend(fontsize=7)
    fig.suptitle('The prescribed field and the supplied particles imply different inlet charge histories');fig.tight_layout();figure(fig,out/'held_charge_mismatch')
    conclusion=dict(status='MAIN_ACCUMULATING_RESIDUAL_LOCALIZED_TO_INLET_FIELD_CHARGE_INCOMPATIBILITY',cases=stats,force=force,gauss=gauss,held_profile=read(root/'HELD_PROFILE.json'),
        limitations=['The two inlet operations (Ex/Ey overwrite and six-component new-layer supply) are jointly implicated, not separately ablated.',
          'No independent full-volume force reconstruction from the Gauss defect; do not equate its L1 norm with net added particle charge.',
          'A roughly 3.5keV remaining offset and the uniform-sphere energy-spread discrepancy are not resolved.',
          'No production fix or guard/default change; all controls are short converging spheres.'])
    save(root/'FINAL.json',conclusion)
    table=['| 工况 | 左侧格数 | 500步平均lab等效能量偏差 |','|---|---:|---:|']
    for name,guard,label in [('base',64,'原窗口'),('transverse_240',64,'仅横向边界移远200um'),('inlet_guard_128',128,'左侧距离对照'),('inlet_guard_256',256,'左侧距离对照')]:
        table.append(f"| {label} | {guard} | {stats[name]['lab_energy_error_from_same_box_slice_eV']/1000:+.5f} keV/电子 |")
    energy=['| 工况 | 最大总账差额 / 固定D |','|---|---:|']
    for name,st in stats.items():energy.append(f"| {name} | {st['relative_ledger']['maximum_absolute_total_residual_pct_fixed_scale']:.6f}% |")
    text='''# 聚焦束团剩余纵向动量积累：左侧入窗场与电荷不兼容

本轮主要积累项已定位到左侧参考供场/新格点注入与实际held粒子电荷的离散匹配。
这是边界操作引入的新约束误差；本轮没有修改内部Cowan/Boris/电流沉积、参考算法或生产设置。
以下为同一原输入32768带权宏粒子、1e8电子、16MPI×8、同dx/dz/dt的短converging对照。
仍开启原boost能量账本和严格离散场账本。所有比较限于全部释放、首个粒子尚未流出的完整cohort。

## 距离控制与实际力

'''+ '\n'.join(table)+'''

256格继续至1000步仍仅−3.52340keV，原64格相同时刻−20.61398keV，已检查范围内改善不只是延后同样大小的漂移。
加宽横向仅小幅改善；左侧128/256格明显降低积累。这些增大窗口仅作定位，不能据此将大guard写入生产卡。

![因果对照与纵向电场](figures/causal_controls.png)

480步真实粒子所采样的平均Ez偏差：base +7391.47、仅横向加宽+7233.48、左侧256格−2757.67V/m。
base额外电冲量−1.81435e−8，磁旋转差约+9.20e−11；即主要差在实际电场而非磁旋转。
实际E/B离线重放Boris的最大动量差8.31e−17；独立理想球形场下局部Boris误差约1e−10/步。
这不能单独证明长时间PIC收敛，但不支持将主要积累直接归为Boris代数错误。
参考力按实际PIC位置/rest时刻计算，用更密集的球壳q(r^3)查表；256/512网格的平均Ez差<0.482V/m，远小于上面的数千V/m。
B同时刻/半步时刻口径均保留，修正后结论不变。力参考仍是球形连续模型，不是对实际扰动后的PIC密度重新求Poisson。

## 新误差产生于左侧补层，内部推进保留它

定义G=epsilon0*div_h(E)-rho_CIC，rho由原程序自带的CIC电荷诊断沉积。
只读记录新左侧第1个保留电荷节点及邻近8层，横向监测±120um，独立打开原参考场作初始同位置对照。

* 原静态参考同位置的积分G最大6.80e−25C；按刚性速度平移后的参考对其自身rho也约5.94e−25C。
* 实际新补入节点的积分G可达4.24e−14C；因此不是原初场已存在同等大小误差。
* 普通推进中电荷连续性最大缺陷8.11e−14C/m³，Ampere散度更新最大缺陷1.67e−14C/m³。
* 同一节点一次推进前后积分G变化<2.37e−29C；只看保留节点，平移后变化<3.98e−28C。

![约束误差在哪里进入](figures/inlet_constraint.png)

因此核内的电荷连续性和Ampere散度关系保持，但新插入节点的初始G没有与实际入粒子匹配。
已有G不会由普通curl推进自动消掉；它作为额外纵向场的约束误差保留，粒子经过时获得额外电冲量。
从离散恒等式可直接看出：epsilon0*Delta(div E)=-dt*div J、Delta rho=-dt*div J，所以保留节点Delta G≈0。
这些是额外Gauss误差的来源证据，不是仅凭能量账本差额推测。

监测区域前149个新层中G的绝对值积分累计约2.975pC，而带符号合计仅−0.001854pC，正负误差大幅抵消。
这不是新增2.975pC真实粒子或整域净电荷；它是约束误差绝对量。净值很小仍可有明显空间场结构。

## 参考供场与实际输入的具体不兼容

movingWindowFieldValue按一个translationBeta刚性平移最初参考场，shift()把全部6个参考分量填入新左层；
afterElectricUpdate又把左端Ex/Ey写回参考值。实际held粒子仍各按自己的frozenBetaZ运动并作原电流沉积。
本工况held beta范围0.196684–0.257262，而参考translation beta为0.22920347。

从原入射事件独立重建held初始坐标，得到rho.bin纵向电荷边缘分布相对L2差2.69e−9，验证了坐标/权重口径。
在新入窗电荷节点，刚性参考rho与真实held轨迹rho的最大差达到真实峰值的12.47%；128格对照的误差绝对量约降一个数量级，256格没有粒子穿过该入口。
将参考速度改成电荷加权均值0.22949550，在这个离线密度比较中仅小幅改善，不能靠换一个平均速度就宣称解决。

还存在格点插值的非交换误差：把已经沉积的参考rho作线性平移，不等于把粒子逐个平移后重新CIC沉积。
额外离线控制将所有粒子赋同一速度再直接CIC，64格入口两种口径的L1差仍有0.2696pC；
共同速度直接CIC与真实逐粒子速度CIC的L1差为0.5814pC。两项有正负抵消，L1不能相加，也不能当作20keV误差的能量分摊。

![参考电荷历史与实际held电荷历史](figures/held_charge_mismatch.png)

主要失配是准稳态供场与逐粒子源在左侧的离散兼容性。边界供场可以近似，但不能持续向新格点写入与实际源不一致的约束模式。
本轮尚未分别消融Ex/Ey覆盖和新层Ez填入的各自贡献，也未从全部三维G独立重构20keV的力；不能把两者中的某一个孤立称为唯一代码错误。

## 能量账本与范围

D沿用全程固定max|同cohort Delta K_boost|，不除以初始总束流动能。

'''+ '\n'.join(energy)+'''

账本保持闭合并不消除这类场—电荷约束模式；因此需同时检验约束和粒子动力学。
剩余约3.5keV及uniform能散偏差未在本轮解决。下一步应在原紧凑窗口内修正入窗兼容性并对照这些量，不能用增加大guard代替修复。

程序与可复现实验身份见[PROGRAM_VALIDATION.md](PROGRAM_VALIDATION.md)，原数据、失败启动记录和生产默认保持。
'''
    (root/'README.md').write_text(text)
    (root/'PROGRAM_VALIDATION.md').write_text('''# 程序验证与实验身份

实际物理对照19061/62/63（原/左256/横240，500步），19064（左128，500步）、19065（左256，1000步）均COMPLETED0:0。
均为原已验build/reference-surface-current-20261030/aprl，核心5c4e4b7，SHA见plan/FROZEN。base500步所有非计时科学账本字段与旧完整converging前缀逐位相同。
19066使用独立build/inlet-constraint-probe-20261030/aprl，仅增加只读rho/散度观察，不写E/B/J/粒子/CPML历史；77required通过，开启300步的18个H5科学records与旧程序前缀逐位。
19067独立原参考只读Gauss检查正常完成，不推进粒子。源码与实际二进制映射见SOURCE_FINAL.json。

最初19058/59/60因运行脚本少复制review_discrete依赖，在进入PIC/绑定前退出；补齐导入预检后新建r2，没有改物理卡/容差/原输出，不计为物理失败或成功测试。
新诊断编译初次常数名eps0写错，在编译阶段已纠正为epsilon0，未生成被用于模拟的错误程序。

每项绑定128不同物理核，16MPI×8，不固定节点；原600/60秒MPI-free watchdog及作业内BTL保持。所有原能量账本/严格审计继续开启。
新的只读探针有额外开销，不作为性能优化，不切生产默认。大guard与加宽横向均只为因果定位，不能写入生产模板。
FROZEN.json、两后续根目录FROZEN、300步探针FROZEN均保留。后处理源码快照、PNG/PDF、独立密度检查、Gauss CSV齐。
''')
    print(json.dumps(dict(cases={k:v['lab_energy_error_from_same_box_slice_eV'] for k,v in stats.items()},gauss=gauss),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('followup',type=Path);p.add_argument('probe',type=Path);a=p.parse_args();main(a.root.resolve(),a.followup.resolve(),a.probe.resolve())
