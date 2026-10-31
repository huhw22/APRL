#!/usr/bin/env python3
"""Render the completed read-only decomposition and unchanged-pusher probes."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--workspace',type=Path,required=True);a=p.parse_args()
r=a.root.resolve();ws=a.workspace;cases={j:json.loads((r/f'time_pairing/{j}/RESULT.json').read_text()) for j in [19099,19098]}
probe=json.loads((r.parent/'external_probe_r1/RESULT.json').read_text());assert probe['integrity_passed']
figs=r/'figures';figs.mkdir(exist_ok=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
colors=['#8546a0','#cc6727','#327e74'];names=['particle_minus_old_mesh_J','old_minus_mid_mesh_J','particle_minus_mid_mesh_J']
labels=['Particle − old-field mesh work','Old-field − midpoint-field mesh work','Actual particle − midpoint mesh work']
normalizations={}
for col,job,oldreview in [(0,19099,'undulator_review_19099_r1'),(1,19098,'free_drift_review_19098_r1')]:
    data=np.genfromtxt(r/f'time_pairing/{job}/decomposition.csv',delimiter=',',names=True)
    norm=json.loads((ws/f'results/corrected_full_comparison_20261030/{oldreview}/science/ledger_summary.json').read_text())['normalization']
    D=norm['fixed_denominator_J'];terminal=abs(norm['particle_change_final_J']);normalizations[str(job)]={'peak_denominator_J':D,'terminal_denominator_J':terminal}
    ax=axes[0,col]
    for name,color,label in zip(names,colors,labels):ax.plot(data['step'],100*data[name]/D,color=color,label=label,lw=1)
    ax.set(title=f'{job} | '+('undulator' if job==19099 else 'free drift'),xlabel='Maxwell step',ylabel='Percent of peak |ΔK_boost|')
    ax.legend(fontsize=7,frameon=False,loc='best')
    ax=axes[1,col];vals=[100*cases[job]['final'][n]/terminal for n in names]
    ax.bar(['Remaining old-field term','Time-pairing term','Actual total gap'],vals,color=colors,width=.6)
    ax.set(ylabel='Percent of terminal |ΔK_boost|',title='Same terms relative to terminal net particle change')
    bottom,top=min(vals),max(vals);span=top-bottom
    ax.set_ylim(bottom-.15*span,top+.15*span)
    for i,v in enumerate(vals):ax.annotate(f'{v:+.5g}%',(i,v),ha='center',va='bottom' if v>=0 else 'top',xytext=(0,4 if v>=0 else -4),textcoords='offset points',fontsize=9)
for ax in axes.flat:ax.axhline(0,color='#777777',lw=.6);ax.grid(axis='y',alpha=.2)
fig.suptitle('Unchanged trajectories | direct saved work decomposition | boost frame\nOne fixed denominator for all curves within each panel; no physical correction applied',fontsize=11)
fig.savefig(figs/'time_pairing_decomposition.png',dpi=180);fig.savefig(figs/'time_pairing_decomposition.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for item,color in zip(probe['input_samples'],['#657b9e','#8546a0','#327e74']):
    rows=[v for v in probe['results'] if v['magnet'] and v['sample_index']==item['row']]
    n=[v['particle_substeps'] for v in rows]
    axes[0].plot(n,[v['max_sampled_absolute_lab_energy_change_eV']/1000 for v in rows],'o-',color=color,label=f'Input row {item["row"]}')
    axes[1].plot(n,[v['lab_energy_change_eV'] for v in rows],'o-',color=color)
for ax in axes:
    ax.set_xscale('log');ax.set_xticks([1,5,10,20],labels=['1','5 (production)','10','20']);ax.set_xlabel('Particle substeps per unchanged Maxwell-size interval');ax.grid(alpha=.2)
axes[0].set_yscale('log');axes[0].set(ylabel='Maximum sampled |ΔK_lab| [keV]',title='Inside the static lab magnet');axes[0].legend(fontsize=8,frameon=False)
axes[1].set_yscale('symlog',linthresh=.001);axes[1].set(ylabel='Signed final ΔK_lab [eV]',title='After the magnet at z ≥ 4.5 m');axes[1].axhline(0,color='#777777',lw=.6)
fig.suptitle('Three actual input particles | unchanged external-only pusher and 170-period undulator\nExact reference: static lab B does no work; no Maxwell evolution or self-consistent feedback',fontsize=10)
fig.savefig(figs/'external_magnet_probe.png',dpi=180);fig.savefig(figs/'external_magnet_probe.pdf');plt.close(fig)
G=929.1305492858543;dt=1.2543262415359034e-6/299792458;energy_lab=91.49184169394452
doppler=1/(G+np.sqrt(G*G-1));energy_boost=energy_lab*doppler;samples=4.135667696e-15/(energy_boost*dt)
phase={'scope':'On-axis forward vacuum plane-wave scale estimate; not a measured phase error or attribution of the ledger gap.',
       'lab_photon_energy_eV':energy_lab,'boost_photon_energy_eV':energy_boost,'field_dt_s':dt,'Maxwell_steps_per_boost_cycle':samples,'half_field_step_phase_deg':180/samples}
(r/'SCALES.json').write_text(json.dumps({'normalization':normalizations,'plane_wave_phase_scale':phase},indent=2)+'\n')
m=cases[19099]['final'];f=cases[19098]['final'];prod=[v for v in probe['results'] if v['magnet'] and v['particle_substeps']==5]
lines=['# 场—粒子能量交换：时间配对核查','',
  '**现有数据支持优先查时间配对，但不支持只改一个时间索引就宣布解决，更不足以判定Cowan必须替换。**本轮保持原推进和物理设置，做了已保存功的只读分解，以及链接原推进库的纯外磁三粒子核查；没有新的全程序PIC或内核修改。','',
  '## 已完成的定量分解','',
  '定义P为实际Boris网格场功扣除已记录的边界／frozen／carrier约束源功；M_old和M_mid为场推进实际采用的同一份J，分别与E_old及(E_old+E_new)/2相乘后的时空累计。于是恒等地有：','',
  '$$R=P-M_{mid}=(P-M_{old})+(M_{old}-M_{mid}).$$','',
  '第一项仍含取场空间权重、粒子速度与位移、轨迹／子步及可能的边界取样差异；第二项单独量出原轨迹上的电场时间配对差。这里没有假装第一项已完全归因于空间插值。','',
  '| 有磁19099，boost系 | nJ |','|---|---:|']
for k,label in [('particle_grid_work_corrected_J','粒子网格功P，已扣约束源'),('mesh_old_work_J','沉积电流与旧场的功M_old'),('mesh_mid_work_J','沉积电流与中点场的功M_mid'),('particle_minus_old_mesh_J','余项P−M_old'),('old_minus_mid_mesh_J','时间配对项M_old−M_mid'),('particle_minus_mid_mesh_J','原实际交换差R')]:lines.append(f'| {label} | {m[k]*1e9:+.9f} |')
lines+=['',f'时间配对项为原差额绝对值的−131.19%，反号余项为+31.19%；两者部分抵消。这是原轨迹的代数分解，不是“换成中点推进后就会精确减少131%”的预测。真正改变取场将同时改变轨迹、电流及外场功，必须实际验收。','',
  f'自由漂移19098中，对应两项为 **{f["particle_minus_old_mesh_J"]*1e15:+.6f} fJ** 与 **{f["old_minus_mid_mesh_J"]*1e15:+.6f} fJ**，相消后才得到 **{f["particle_minus_mid_mesh_J"]*1e15:+.6f} fJ**。这说明不能仅凭自由漂移总余差小就认定每一项耦合都已足够准确。','',
  '![原轨迹上的时间配对与余项](figures/time_pairing_decomposition.png)','',
  '[有磁数据](time_pairing/19099/decomposition.csv) · [自由漂移数据](time_pairing/19098/decomposition.csv) · [分解完成记录](time_pairing/COMPLETE.json)','',
  '图上半部沿用此前每工况D=max|同cohort粒子ΔK_boost|；下半部另以终态净变化的绝对值作固定分母展示同三项。每个面板内所有曲线共用分母；两种尺度明确区分，不覆盖原账本或丢失J原量。','',
  '## 代码核查：当前真实的推进顺序','',
  '1. 先更新B，再建立粒子插值halo；E尚未做本步Ampere更新。',
  '2. 每个活动粒子在本步起始位置取一次网格E/B，在整个Maxwell步内保持这个取样。',
  '3. 有磁生产组有5个Boris子步；各子步重新计算外加波荡器场，但不重新计算自洽网格场。',
  '4. 推进结束后，电流按粒子路径沉积；再更新E。原严格场功自然采用更新前后E的平均值。','',
  '场在Yee时间格上交错本身是正常的；必须核对粒子动量、位置、E与B各自的时间标签及入射初始化，不能看到半步差就直接认定错误。当前顺序确实给出了需要检查的接口：粒子为何能用旧取样，而场能实际变化要求与中点功配对；B的更新和粒子子步时刻也应一并核对。','',
  '空间上，当前取场是各原生分量的三线性点值插值，电流是逐格轨迹积分。路径的功与Boris电踢的功还依赖不同的速度／位移表达；这与Maxwell用哪种色散优化算子是两个可区分的问题。','',
  f'按生产G和dt估算，91.49184 eV的同向轴上真空平面波在boost系约{energy_boost:.6f} eV，每周期只有约{samples:.3f}个Maxwell步，半个场步相当于约{180/samples:.3f}°相位。5个外场粒子子步不会自动把这个自洽场采样也加密5倍。这只是尺度估计，未证明真实差额来自该相位。','',
  '## 纯外磁核查：三条真实输入，原核保持','',
  '从原始输入取第0、中间和最后一条记录，保持各自输入位置／动量；使用原库、原boost、原170周期波荡器、原Maxwell尺寸dt，只在独立小程序中关闭网格反馈。没有Maxwell或Cowan推进。lab静态磁场的理论参照为dK_lab/dt=q v·(v×B)=0。表内是出口减入口的能量变化，正值表示增能。','',
  '| 输入行 | 5子步：出口ΔK_lab，eV | 5子步：途中采样最大绝对偏差，keV |','|---|---:|---:|']
for row in prod:lines.append(f'| {row["sample_index"]} | {row["lab_energy_change_eV"]:+.9f} | {row["max_sampled_absolute_lab_energy_change_eV"]/1000:.6f} |')
lines+=['',
  '5→10→20个子步时，途中偏差约按每次四分之一缩小，表现出这三条轨迹上的二阶趋势。12条无磁零场控制均保持常量到数值精度。三条纯外磁轨迹的出口误差不足以单独解释它们若出现221.8 keV量级的变化；但它们不是完整束流，也没有自洽反馈，不能据此排除Boris中途误差经耦合影响结果。','',
  '![原核纯外磁时间精度](figures/external_magnet_probe.png)','',
  '[原核核查数据及输入记录](../external_probe_r1/RESULT.json) · [编译和运行记录](../external_probe_r1/slurm-19108.out)','',
  '## 下一步应如何兼顾速度','',
  '优先级是把同一段物理轨迹中的“粒子电踢功—取场值乘位移—轨迹沉积旧场功—实际中点场功”四项分开量测。这样才能决定需要修的是时间取样、子步路径／速度，还是空间权重；这一阶段只加可关闭的诊断，不改实际推进。','',
  '之后先做粒子端的时间居中／一次预测原型，并保留单次全网格推进。若要改实际取场或轨迹规则，须先列明时间标签、入射动量如何保持、MPI边界和额外数据需求，按既有约定确认后再改。只增加外磁子步不能当作自洽耦合修正；完整迭代和全网格dt减半放在后面评估。',
  '全局减半dz和dt通常同时增加纵格和步数，场成本可能接近四倍；而粒子端局部额外取样／预测有希望更便宜，但活动粒子数及额外MPI需求必须实测，不能承诺零成本。','',
  '若时间标签和路径配对厘清后仍有boost下E/B抵消问题，可再比较Vay或Higuera–Cary等显式粒子推进，保留Cowan场核。更换粒子推进不自动得到电流沉积的能量一致性，不能用它跳过交换功、动量及谱的验证。本轮没有采用这些方法。','',
  '验收同时保留有磁和漂移，检查lab粒子平均／std能损及ΔPz、boost各功项与余差、窄带场谱／BF、同资源每步时长和总节点时。只让一个残差靠抵消变小，或只重标定能量，并不能达到目标。','',
  '## 理论来源和复核边界','',
  '[WarpX官方PIC说明](https://warpx.readthedocs.io/en/22.11/theory/picsar_theory.html)给出粒子／场时间交错及Boris、Vay的区别；它支持核查时间标签及boost下受力抵消的必要性，不证明本程序某一行已错。',
  '[Ricketson与Hu：显式能量守恒PIC](https://arxiv.org/html/2411.09605v1)区分时间与空间离散的能量条件；该文自身把相对论推广列为后续工作，不能直接移植为本程序的现成解法。这里没有引入其粒子能量校正。','',
  '19107只读分解3秒；19108独立原核检查及编译2秒，均在node3 COMPLETED0:0。两次都不是新的全量PIC。完整代码、冻结身份、源文件stat及验收范围保留；分析通过不等于已经修好耦合。']
(r/'README.md').write_text('\n'.join(lines)+'\n')
(r/'FINAL.json').write_text(json.dumps({'status':'DIAGNOSED_LARGE_TIME_PAIRING_WITH_NONZERO_REMAINDER','new_full_PIC_simulations':0,'physics_method_changes':False,
  'time_pairing_results':cases,'external_only_probe':probe,'phase_scale':phase,'physical_coupling_fixed':False},indent=2)+'\n')
print(json.dumps({'boost_wave_steps_per_cycle':samples,'half_step_phase_degrees':180/samples,'magnet_decomposition_nJ':[m[n]*1e9 for n in names]},indent=2))
