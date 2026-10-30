#!/usr/bin/env python3
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-charge-inlet-mpl');os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(p.read_text())
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def csv(p):return np.genfromtxt(p,names=True,delimiter=',')
def figure(fig,p):fig.savefig(p.with_suffix('.png'),dpi=180,bbox_inches='tight');fig.savefig(p.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig)
def main(root,aux,audit):
    plan=read(root/'plan.json');b=root.parent;prior=b/'reference_surface_current_20261030';out=root/'figures';out.mkdir(exist_ok=True)
    pilot=root/'pilot_500';assert read(pilot/'COMPLETION.json')['checks_pass']
    assert read(aux/'guard32_500/COMPLETION.json')['checks_pass'] and read(aux/'guard32_500/ENTRY_PHASE_SPACE_BITWISE.json')['pass_all']
    assert read(aux/'two_nodes_200/COMPLETION.json')['checks_pass']
    old=np.genfromtxt(b/'surface_residual_inlet_probe_20261030/probe_300/diagnostic/inlet-constraint.csv',names=True,delimiter=',',dtype=None,encoding=None)
    new=np.genfromtxt(pilot/'diagnostic/inlet-constraint.csv',names=True,delimiter=',',dtype=None,encoding=None)
    sel=lambda a:a[(a['stage']=='begin')&(a['k']==1)&(a['step']>0)&(a['step']<300)]
    x,y=sel(old),sel(new);assert np.array_equal(x['step'],y['step'])
    fig,ax=plt.subplots(1,2,figsize=(12,4.3))
    for p,label,c in [(b/'surface_residual_origin_20261030_r2/base/analysis/common_boost_time.csv','Original 64-cell inlet','#66788d'),
        (pilot/'momentum-analysis/common_boost_time.csv','Charge-matched 64-cell inlet','#148572')]:
        a=csv(p);ax[0].plot(a['step'],a['lab_energy_error_from_same_box_slice_eV']/1e3,label=label,color=c)
    if (aux/'guard32_500/COMPLETION.json').exists():
        a=csv(aux/'guard32_500/momentum-analysis/common_boost_time.csv');ax[0].plot(a['step'],a['lab_energy_error_from_same_box_slice_eV']/1e3,label='Charge-matched 32-cell inlet',color='#b97d32')
    ax[0].set_xlabel('Field step: complete released cohort');ax[0].set_ylabel('Same boost-slice lab mean-energy error (keV)');ax[0].axhline(0,color='k',lw=.6)
    for a,label,c in [(x,'Original supplied layer','#66788d'),(y,'Charge-matched supplied layer','#148572')]:ax[1].semilogy(a['step'],a['gauss_rms_C_per_m3'],label=label,color=c)
    ax[1].set_xlabel('Field step');ax[1].set_ylabel(r'RMS $\epsilon_0\nabla_h\cdot E-\rho$ (C/m$^3$)')
    ax[1].set_title('New inlet node k=1; central +/-120 um aperture',fontsize=9)
    for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('Trusted entry events fixed; upstream continuation is an auxiliary model');fig.tight_layout();figure(fig,out/'inlet_completion')
    summary={};rows=['| 束团 | 旧方案平均lab ΔK偏差 | 新方案偏差 | 新横向RMS尺寸偏差 |','|---|---:|---:|---:|'];rest=[];energy=[]
    fig,axes=plt.subplots(3,3,figsize=(13,10));colors={'old':'#66788d','new':'#148572'}
    for row,k in enumerate(('uniform','converging','gaussian')):
        d=root/(k+'_full')
        if (d/'COMPLETION.json').exists():
            assert read(d/'COMPLETION.json')['checks_pass'];entry=read(d/'ENTRY_PHASE_SPACE_BITWISE.json');assert entry['pass_all'];ledger=read(d/'LEDGER_CHECK.json')
        else:
            closure=read(d/'REVIEWED_CLOSED_OUTPUT.json');assert closure['actual_source_states_bitwise'] and closure['pic_exit_code']==0
            entry=read(d/'REVIEWED_ENTRY_DIAGNOSTIC.json');ledger=read(d/'REVIEWED_LEDGER_CHECK.json')
            assert read(audit/(k+'_SOURCE_CHECK.json'))['pass_all']
        r=read(d/'analysis/RESULT.json');o=read(prior/(k+'_full')/'analysis/RESULT.json');n,p=r['planes'][-1],o['planes'][-1]
        s=dict(old_lab_error_keV=(p['PIC_mean_lab_gain_eV']-p['reference_mean_lab_gain_eV'])/1000,new_lab_error_keV=(n['PIC_mean_lab_gain_eV']-n['reference_mean_lab_gain_eV'])/1000,
            new_lab_gain_keV=n['PIC_mean_lab_gain_eV']/1000,theory_at_new_event_lab_gain_keV=n['reference_mean_lab_gain_eV']/1000,
            new_rest_mean_gain_eV=n['PIC_mean_rest_gain_eV'],theory_rest_mean_gain_eV=n['reference_mean_rest_gain_eV'],new_rest_std_gain_eV=n['PIC_std_rest_gain_eV'],theory_rest_std_gain_eV=n['reference_std_rest_gain_eV'],
            sigma_x_relative_error=n['PIC_sigma_x_m']/n['reference_sigma_x_m']-1,sigma_z_rest_relative_error=n['PIC_sigma_z_rest_m']/n['reference_sigma_z_rest_m']-1,
            entry_plane_diagnostic=entry,actual_source_states_bitwise=True,relative_energy_budget=r['relative_energy_budget'],ledger=ledger,
            common_slice=read(d/'momentum-analysis/ANALYSIS.json')['final_clean_comparison'])
        summary[k]=s;rows.append(f"| {k} | {s['old_lab_error_keV']:+.3f} keV | {s['new_lab_error_keV']:+.3f} keV | {100*s['sigma_x_relative_error']:+.4f}% |")
        rest.append(f"| {k} | {s['new_rest_mean_gain_eV']:.3f} / {s['theory_rest_mean_gain_eV']:.3f} | {s['new_rest_std_gain_eV']:.3f} / {s['theory_rest_std_gain_eV']:.3f} |")
        rb=s['relative_energy_budget'];energy.append(f"| {k} | {rb['maximum_absolute_total_residual_pct_fixed_scale']:.6f}% | [图]({k}_full/analysis/energy_budget.png) |")
        for col,(key,scale,label) in enumerate([('mean_lab_gain_eV',.001,'Mean lab gain (keV)'),('mean_rest_gain_eV',1,'Mean rest-frame gain (eV)'),('sigma_x_m',1e6,'Transverse RMS size (um)')]):
            a=axes[row,col];z=[v['lab_z_m'] for v in r['planes']]
            a.plot(z,[v['reference_'+key]*scale for v in r['planes']],color='k',label='Theory at new PIC event times')
            a.plot(z,[v['PIC_'+key]*scale for v in o['planes']],'s--',color=colors['old'],label='Old inlet')
            a.plot(z,[v['PIC_'+key]*scale for v in r['planes']],'o-',color=colors['new'],label='Charge-matched inlet');a.set_title(k);a.set_ylabel(label);a.grid(alpha=.2)
            if row==2:a.set_xlabel('Fixed lab detector z (m)')
    axes[0,0].legend(fontsize=7);fig.suptitle('Independent spherical dynamics after the immutable input plane');fig.tight_layout();figure(fig,out/'sphere_physics')
    layer=csv(pilot/'reference/inlet-charge-match.csv');assert np.isfinite(layer['native_layer_energy_change_J']).all()
    preparation=read(pilot/'reference/preparation.json')
    def timing(path):
        import re
        line=next(x for x in path.read_text().splitlines() if 'phase=progress step=500 ' in x)
        return {n:float(re.search(r'\b'+n+r'=([0-9.eE+\-]+)',line).group(1)) for n in ('elapsed_s','current_rss_total_mib')}
    old_timing=timing(b/'surface_residual_origin_20261030_r2/base/aprl.log');new_timing=timing(pilot/'aprl.log')
    cost=dict(background_preparation_seconds=preparation['preparationSeconds'],same_500_step_old=old_timing,same_500_step_new=new_timing,
        observed_progress_seconds_relative_change=new_timing['elapsed_s']/old_timing['elapsed_s']-1,
        scope='One unpaired 500-step sample with diagnostic probes; not stable throughput or full production overhead. Preparation includes the existing scalar background solve and IO, not a full PIC relaxation.')
    result=dict(status='REVIEWED_SMALL_CASE_COMPARISON',production_spectrum_accepted=False,contract=read(root/'TRUSTED_INPUT_CONTRACT.json'),cases=summary,
        pilot=read(pilot/'momentum-analysis/ANALYSIS.json')['final_clean_comparison'],
        independent_new_layer_gauss_rms_max_C_per_m3=float(max(y['gauss_rms_C_per_m3'])),independent_old_layer_gauss_rms_max_C_per_m3=float(max(x['gauss_rms_C_per_m3'])),
        independent_gauss_probe_scope='begin-step new k=1 charge node, central |x|,|y|<=120 um ROI; steps 1..299, not a whole-aperture norm',
        pilot_native_completion_energy_J=float(np.sum(layer['native_layer_energy_change_J'].astype(np.longdouble))),
        two_node_prefix=read(aux/'two_nodes_200/COMPLETION.json'),guard32=read(aux/'guard32_500/momentum-analysis/ANALYSIS.json')['final_clean_comparison'],
        source_state_audits={k:read(audit/(k+'_SOURCE_CHECK.json')) for k in ('uniform','gaussian')},cost=cost)
    save(root/'FINAL.json',result)
    import csv as csv_module
    with (root/'comparison.csv').open('w') as f:
        names=['case','old_lab_error_keV','new_lab_error_keV','new_rest_mean_gain_eV','theory_rest_mean_gain_eV','new_rest_std_gain_eV','theory_rest_std_gain_eV','sigma_x_relative_error']
        writer=csv_module.DictWriter(f,fieldnames=names);writer.writeheader()
        for k,s in summary.items():writer.writerow(dict(case=k,**{n:s[n] for n in names[1:]}))
    text='''# 以可信输入面为约束的快速入窗初场

## 物理约束

唯一固定的是用户距波荡器前0.05m探测面上的逐粒子相空间/到达事件/权重。已知之前0.5m是漂移，
但不把held或其它延拓当作可信的完整历史，不反演唯一真实上游场。初场采用准稳态数值补全，
不额外指定上游入射辐射。具体输入契约见[TRUSTED_INPUT_CONTRACT.json](TRUSTED_INPUT_CONTRACT.json)。

本轮人工球形束团测试把输入面设为z=0且没有波荡器；它检查输入边界保持和之后的空间电荷动力学，
不冒称已做真实0.05m波荡器生产工况。

## 实现

新增`mesh.moving_window.inlet_charge_match: true`，默认false。复用已验快速背景场生成器只是实现选择；
新完成器仅接受数值入口中的粒子状态和邻接场，可以搭配其它准稳态生成器和上游延拓。
在新层的原生法向E位置，用第一个保留电荷节点的真实CIC rho和相邻已推进场完成离散Gauss关系。
只改新层Ez的横向物理内部节点，保留原64/128guard、网格/步长、输入粒子、B和CPML历史、
Cowan/Boris/沉积与内部已推进场。近似背景继续给切向场；不要求它复现held的精确物理历史。

pending只查看本步应入窗的事件，不提前admit、不扫描整个等待束团。额外内存为横向rho平面，
工作限于新层、活动窗粒子扫描及本步待入窗事件，无新增全空间迭代/PIC等待稳定或0.5m预跑。
完成器在shift内部测量inserted/remap之前执行，新增能量已计入原边界remap，不能重复加进source work。

## 可信输入保持与500步对照

输入H5和释放算法保持。converging完整束团、500步pilot和32格对照的entry面全部32768条记录逐位相同。
uniform和gaussian的入口面诊断，分别有48和136个横向动量分量相差最多4.34e-19（1 ULP），其它字段逐位相同。
入口端点转换略低于平面时，过面采样器会在第一个受力段插值极小的比例；其舍入依赖入场后的力。
为区分它与真实输入变化，另用150/270步开关回放，直接记录全部32768粒子第一次受力前的状态：
source ID、步数、时刻、子步长、电荷、位置和动量均逐位相同。诊断回放的所有科学H5前缀也与原完整记录逐位相同。
因此没有通过调整可信输入来改善输出；不能据此把那两个入口面诊断记录说成逐位相同。
详细证据见[入口源状态核验](../'''+audit.name+'''/uniform_SOURCE_CHECK.json)和[高斯核验](../'''+audit.name+'''/gaussian_SOURCE_CHECK.json)。

原64格窗口的converging500步同boost完整cohort lab等效平均能量误差−17.072→−3.371keV/电子。
独立原CIC诊断显示新左层k=1、x/y各±120um中心区域的Gauss误差降到舍入水平，
而非只用完成器内部自检；图中的RMS限于该区域，不冒称全孔径独立扫描。

![入口相容性与能量变化](figures/inlet_completion.png)

## 原三类完整短束团

以下平均能量偏差均为PIC减去对应实际过面rest时刻的独立球壳世界线参考；新旧到达时间略变，理论分别取值。
各粒子Delta K相对于自己的入口事件，不是不同同时性的两份总束能量相减。

'''+ '\n'.join(rows)+'''

| 束团 | rest平均Delta K：PIC / 理论 (eV) | rest std(Delta K)：PIC / 理论 (eV) |
|---|---:|---:|
'''+ '\n'.join(rest)+'''

![空间电荷动力学对照](figures/sphere_physics.png)

改进并非每项都单调：冷均匀球末面的平均lab能量偏差由+0.790变为−1.048keV，绝对值略增；
它的共同boost末个完整cohort时刻偏差约−4.039keV。聚焦和高斯的均值明显改善，
但不能把“新层Gauss约束满足”当作所有动力学误差消失。32格准备区500步聚焦偏差约−5.069keV，
大于64格的−3.371keV；本轮保留原64格，不据此自动缩窗。

## 相对能量账本

保持原D=max|同cohort Delta K_boost|为全程固定分母，包括保留/移除粒子，不除以初始总束流动能。
新层电场完成已进入remap能量。严格场算子闭合与总粒子—网格剩余项分开保留，不据闭合声称物理全部收敛。

| 工况 | 最大累计总差额 / D | 相对账本图 |
|---|---:|---|
'''+ '\n'.join(energy)+'''

## 验收范围

这是相空间固定、上游只作数值补全的首版薄层实现。没有验证所有上游补全、真实大束团/有磁能谱，
剩余能散和均值差必须按FINAL.json逐项看，不能因Gauss误差小而宣称所有物理正确。
冷均匀球rest能散仍约高9.6%，本轮没有消除它。该方案满足的是廉价准稳态入口的初步对照，不是唯一真实上游场的重建。
不擅改内核或大guard，不切生产default/cluster.env，不自动提交全量生产或推远程。
双节点/关闭逐位、单位检查、性能观察与源码身份另见[PROGRAM_VALIDATION.md](PROGRAM_VALIDATION.md)。
'''
    (root/'README.md').write_text(text)
    (root/'PROGRAM_VALIDATION.md').write_text('''# 程序验证与范围

19068：78required通过（含新层Gauss制造场、保留场/输入源不变、pending visitDue纯读且与真正admit一致），关闭100步18H5科学records与旧binary逐位。
19069：原64格converging500步，独立只读Gauss/真实pusher E/B冲量、原双账本；entry全部字段逐位保持。
19071：converging完整短束团，保留全部6面；同原物理输入/终点。
19070/19073：uniform/gaussian的PIC正常退出、输出闭合且各6面32768粒子齐；随后包装脚本的入口面逐位断言失败，原Slurm FAILED与FAILED.json保留。
另读闭合H5逐文件核验complete/committed/finite及账本，写REVIEWED_CLOSED_OUTPUT/REVIEWED_ENTRY_DIAGNOSTIC/REVIEWED_LEDGER_CHECK；没有篡改原COMPLETION或放宽断言。
19075/19076：150/270步开关回放，新增首次受力前只读记录，原科学H5前缀逐位；全部32768输入源状态逐位不变，入口面采样差单列。19075另过78required。
19072：两节点总16MPI×8，200步18H5科学records与单节点前缀逐位；每节点8rank，128实际物理核而非256实际线程。
19074：左32格、活动Nz减32/右端保持，500步，用于数值准备区敏感性，非自动生产选择。

同MPI/同工具链/128不同物理核，原600/60秒MPI-free watchdog和作业内BTL保持。没有改系统网络/默认程序。
辅助prepare首次把未冻结的__pycache__目录当文件计算哈希失败，发生在提交前，过滤目录后完成清单；未改变物理卡或任何已冻结结果。
编译首次unit使用未定义FieldVector等号，已在编译阶段改逐分量检查，没有用失败产物启动。
SOURCE_BUILD/SOURCE_FINAL给出实际C++/header与compiled-parent-dirty映射，不冒称编译于后来的提交。
首个source-audit准备目录仅生成未提交；在提交前另建r2，删除短审计不需要的下游探测面，避免旧下游面约束影响短回放。所有冻结目录保留。
''')
    with (root/'PROGRAM_VALIDATION.md').open('a') as f:
        f.write(f"\n背景准备 {cost['background_preparation_seconds']:.3f} s；同500步progress耗时 {old_timing['elapsed_s']:.3f} → {new_timing['elapsed_s']:.3f} s（{100*cost['observed_progress_seconds_relative_change']:+.3f}%），仅一个含重诊断样本，不能视为稳定开销或外推大工况。新增rho平面内存按(Nx+1)(Ny+1)double，pending不回到全队列逐步扫描。\n")
    print(json.dumps({k:{n:v[n] for n in ('old_lab_error_keV','new_lab_error_keV','sigma_x_relative_error')} for k,v in summary.items()},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('aux',type=Path);p.add_argument('audit',type=Path);a=p.parse_args();main(a.root.resolve(),a.aux.resolve(),a.audit.resolve())
