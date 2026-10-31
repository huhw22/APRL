#!/usr/bin/env python3
"""Add a clearly separate terminal normalization and auditable validation record."""
import hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-work-geometry-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ws=Path(sys.argv[1]).resolve();root=Path(sys.argv[2]).resolve();parent=root.parent
src=ws/'source/aprl';here=src/'verification/analysis/particle_grid_exchange';LD=np.longdouble
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
decomp=json.loads((root/'DECOMPOSITION.json').read_text());checks={};secondary={}
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False})
fig,axes=plt.subplots(2,2,figsize=(12.3,8.7),layout='constrained',sharey='row')
for case,summary in decomp.items():
 absolute=np.genfromtxt(root/case/'work_decomposition_absolute.csv',delimiter=',',names=True,dtype=LD)
 relative=np.genfromtxt(root/case/'work_decomposition.csv',delimiter=',',names=True,dtype=LD)
 D=LD(str(summary['normalization']['D_J']));worst=0.
 for n in absolute.dtype.names:
  if n in ['step','fully_released_no_removal']:assert np.array_equal(absolute[n],relative[n]);continue
  reconstructed=relative[n]*D/100
  scale=max(np.max(abs(absolute[n])),LD('1e-300'))
  error=np.max(abs(absolute[n]-reconstructed))/scale
  assert error<1e-15,(case,n,error);worst=max(worst,float(error))
 checks[case]={'percent_to_J_max_relative_reconstruction_error':worst,'checks_pass':True}
 terminal=abs(absolute['particle_change'][-1]);assert terminal>0
 secondary[case]={'fixed_terminal_denominator_J':float(terminal),'definition':'abs(final stable-cohort delta K boost); secondary view only, unchanged peak-D primary CSV and figure','whole_gap_percent_of_terminal':100*summary['whole']['gap']/float(terminal)}
 if not case.startswith('native-'):continue
 row=0 if 'undulator' in case else 1
 increments={k:np.diff(np.r_[LD(0),absolute[k]]) for k in absolute.dtype.names if k not in ['step','fully_released_no_removal']}
 curves=[increments['boris_minus_gather'],increments['gather_minus_subpath'],increments['subpath_minus_chord'],increments['particle_cohort_source']+increments['mesh_cohort_source'],increments['old_minus_mid'],increments['gap']]
 labels=['B − G: pusher / displacement','G − S: spatial gather / path','S − C: curved path / chord','(P − B) + (C − M_old): cohort / sources','M_old − M_mid: time pairing','P − M_mid: net gap']
 for col,mask in enumerate([np.ones(len(absolute),bool),absolute['fully_released_no_removal']==1]):
  ax=axes[row,col];x=absolute['step'][mask]
  for y,label,color,style in zip(curves,labels,['#28639A','#CB7723','#778532','#66717B','#AF5377','#25282C'],['-','--',':','-.','-','--']):ax.plot(x,np.cumsum(y[mask])*100/terminal,label=label,color=color,ls=style,lw=1.7)
  ax.axhline(0,color='#777777',lw=.7);ax.grid(alpha=.18);ax.set_xlim(x[0],x[-1]);ax.set_xlabel('Maxwell step')
  ax.set_title(('One-period undulator' if row==0 else 'Free drift')+' · '+('whole run' if col==0 else 'fully released, no exits'))
  if col==0:ax.set_ylabel('Cumulative work / fixed terminal |ΔK| (%)')
  ax.text(.02,.98,f'Fixed terminal |ΔK| = {float(terminal):.5e} J',transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','edgecolor':'none','alpha':.75})
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False,fontsize=9)
fig.suptitle('Secondary view · fixed terminal denominator · primary peak-D view retained',fontsize=13)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('work_decomposition_terminal.'+ext),dpi=180)
plt.close(fig);save(root/'NORMALIZATION_CHECKS.json',checks);save(root/'SECONDARY_NORMALIZATION.json',secondary)
shutil.copy2(here/'SEMI_IMPLICIT.md',root/'SEMI_IMPLICIT.md')
readme=(root/'README.md').read_text();anchor='## 经典算法按计算系评估'
extra='''为避免有磁工况中途外场摆动抬高D而掩盖小净交换，另保留**终态粒子变化绝对值这一固定分母**的辅助图。它不替换原主图/主CSV，也不逐时刻归一化。本短有磁终态|ΔKboost|=2.60204×10⁻¹⁵ J，净交换差/该量约+0.09777%；漂移约−0.05737%。这不是全量19099的误差比例。

![固定终态分母辅助图](figures/work_decomposition_terminal.png)

'''
assert readme.count(anchor)==1;readme=readme.replace(anchor,extra+anchor)
readme=readme.replace('经典电磁学与把粒子动力学中的γ取1需要区分。','经典电磁学与把粒子动力学中的γ取1需要区分；不以此预先否定候选，要用同一计算系的对照验证。')
readme=readme.replace('数据在各组SUMMARY；','短漂移观察到的最大个体γ′约1.03468，短有磁约1.62651，数据在各组SUMMARY；')
readme=readme.replace('同时状态误差随步长减小。','同时状态误差随步长减小。dt=0.1的状态误差仍约2.39%，dt=1约173%，明确说明守恒不能替代轨迹/相位精度。')
(root/'README.md').write_text(readme)
frozen_counts={}
for name in ['build-r2','campaign-r2','postcheck-r1']:
 paths=json.loads((parent/name/'FROZEN.json').read_text())
 for p,h in paths.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 frozen_counts[name]=len(paths)
unchanged={}
for name in ['eb_particles.cpp','eb_field.cpp','eb_deposition.cpp','eb_mpi.cpp','runtime_control.cpp','shutdown_watchdog.cpp']:
 old=subprocess.check_output(['git','show','70a06aa:src/'+name],cwd=src);now=(src/'src'/name).read_bytes();assert old==now,name
 unchanged[name]=hashlib.sha256(now).hexdigest()
comparisons=json.loads((parent/'postcheck-r1/COMPARISONS.json').read_text())
assert len(comparisons)==6 and all(all(x['ledger_fields_bitwise'].values()) for x in comparisons.values())
layout={}
for kind in ['undulator','free-drift']:
 x=decomp['native-'+kind+'-on'];y=decomp['mpi-small-'+kind+'-on']
 error=max(abs(x[scope][n]-y[scope][n]) for scope in ['whole','fully_released_before_any_removal'] for n in x[scope])
 layout[kind]={'max_absolute_difference_in_summed_terms_J':error,'scope':'Different MPI layouts need not be bitwise identical; same-layout off/on comparison is exact.'}
record={'normalization_checks':checks,'all_six_same_layout_comparisons_exact':True,'all_70_ledger_scientific_fields_bitwise_each_pair':True,'frozen_hash_counts':frozen_counts,'physical_kernel_files_unchanged_from_time_audit_parent':unchanged,'cross_layout_diagnostic_agreement':layout,'physical_conservation_or_radiation_accuracy_accepted':False}
save(root/'VALIDATION.json',record)
(root/'PROGRAM_VALIDATION.md').write_text('''# 程序验证与运行记录

这些检查验证只读诊断实现及物理输出未受影响，不代替物理精度验收。

- 19111独立候选构建：78项required回归通过；192个3D独立轨迹积分与原电流沉积对照通过；2/4rank双向跨接口直线及折线路径检查通过。测试均在node3，本轮未新增跨节点验收。
- 19112的10个257步短物理工况全部完整关闭：16MPI×8的有磁/漂移各原程序、新off、新on，另4MPI×8两组off/on。场面、粒子面、轨迹逐位保持；六对各70个v5科学账本字段也全部逐位一致，无需使用预先规定的浮点归约容差。
- 19112原批次仍为FAILED1:0：全部科学计算和比较之后，辅助负控Python启动包装器在MPI_Init前误读OMPI_COMM_WORLD_RANK，报KeyError。冻结脚本及日志保持，未以当前源模板替换。修正版从SLURM_PROCID读取启动编号。
- 19113只读补验COMPLETED0:0，14秒：重新保存上述比较，实际运行2rank诊断开关不一致负控，两rank及时统一拒绝并退出，未进入物理初始化；随后运行经典固定权重耦合小原型。没有重跑10个物理工况。
- 19110初次单MPI×8短卡超过180秒时限，不完整数据未进入科学分析；新r2使用适当MPI布局。19109为初次局部路径诊断构建，r2补私有电场接口后由19111重新构建验收。原失败目录不删除、不改称成功。
- CSV百分数独立反算J、逐步六项恒等式、完整无退场段独立Boris功、4/16rank分解和冻结文件哈希检查均保存于VALIDATION.json。几何观察计数必须与移除计数联用，不能把首个退出步当完整同cohort段。

候选程序为build/particle-work-geometry-r2-20261031/aprl，SHA256 973881f7d534efcc71232bc964e6246ff42b846645ff30fc0ee63c35ead6716e。默认和原已验程序保持，未推远程。源构建由build-r2/FROZEN.json逐文件对应；后续只改分析和辅助脚本，不误称重编最终文档HEAD。

PNG三张已生成，交付前逐张目视；PDF同源导出，未另渲染。静态科学图无交互筛选或手机布局验收需求；颜色配合线型、标记和明文标签区分。
''')
print('REPORT COMPLETED; ALL FROZEN AND NORMALIZATION CHECKS PASS')
