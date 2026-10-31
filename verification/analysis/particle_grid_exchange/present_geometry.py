#!/usr/bin/env python3
"""Present only closed geometry diagnostics; keep raw J and fixed-D percentages."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);p.add_argument('postcheck',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
assert json.loads((a.postcheck/'COMPLETE.json').read_text())['checks_passed']
root=a.output.resolve();root.mkdir(exist_ok=False);(root/'figures').mkdir()
plan=json.loads((a.campaign/'plan.json').read_text());N=plan['particles'];LD=np.longdouble
terms=['boris_minus_gather','gather_minus_subpath','subpath_minus_chord','particle_cohort_source','mesh_cohort_source','old_minus_mid']
names=['B − G','G − S','S − C','P − B','C − M_old','M_old − M_mid']
colors=['#28639A','#CB7723','#778532','#66717B','#A1A6AC','#AF5377']
source_hashes={};stats={};series={}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,indent=2,default=lambda x:float(x) if isinstance(x,np.generic) else str(x))+'\n')
def intervals(x):
 if len(x)==0:return []
 edges=np.r_[0,np.where(np.diff(x)!=1)[0]+1,len(x)]
 return [[int(x[edges[i]]),int(x[edges[i+1]-1])] for i in range(len(edges)-1)]
for case in plan['cases']:
 if case['mode']!='on':continue
 directory=Path(case['directory']);result=json.loads((directory/'RESULT.json').read_text());assert result['integrity_passed']
 files=sorted((directory/'geometry').glob('*.csv'));assert len(files)==case['ranks'];acc=None
 for path in files:
  source_hashes[str(path)]=sha(path);source_hashes[str(path.with_suffix('.complete'))]=sha(path.with_suffix('.complete'))
  v=np.genfromtxt(path,delimiter=',',names=True,dtype=LD)
  if acc is None:acc={n:v[n].copy() for n in v.dtype.names}
  else:
   for n in acc:
    if n=='step':assert np.array_equal(acc[n],v[n])
    elif n=='max_observed_gamma_box':acc[n]=np.maximum(acc[n],v[n])
    else:acc[n]+=v[n]
 path=directory/'mesh_work.npz';source_hashes[str(path)]=sha(path)
 with np.load(path) as z:m={k:np.asarray(z[k],dtype=LD) for k in z.files}
 with h5py.File(directory/'ledger/energy-ledger.h5') as f:removed=f['energy_ledger/records']['removed_macroparticles']
 steps=acc['step'].astype(int);assert np.array_equal(steps,m['step'][1:])
 # A numerical full push may later be clipped at an exit. Counts alone do not
 # certify that cohort: exclude the first removal step and every later step.
 clean=(acc['eligible_pushes']==N)&(acc['observed_pushes']==N)&(acc['excluded_release_pushes']==0)&(acc['excluded_boundary_interface_pushes']==0)&(removed[1:]==0)&(removed[:-1]==0)
 assert clean.any()
 B=acc['eligible_boris_grid_J'];G=acc['eligible_gather_displacement_J'];S=acc['eligible_subpath_old_mesh_J'];C=acc['eligible_chord_old_mesh_J']
 rawP=np.diff(m['particle_grid_work_J']);P=rawP-np.diff(m['constraints_J']);old=np.diff(m['mesh_old_J']);mid=np.diff(m['mesh_mid_J'])
 values=dict(zip(terms,[B-G,G-S,S-C,P-B,C-old,old-mid]));R=P-mid
 identity=sum(values.values())-R
 bound=256*np.finfo(LD).eps*max(np.max(abs(x)) for x in [B,G,S,C,P,old,mid])
 assert np.max(abs(identity))<=bound
 # All selected physical particles independently reproduce the raw Boris work.
 raw_bound=512*np.finfo(float).eps*np.maximum(abs(m['particle_grid_work_J'][1:])+abs(m['particle_grid_work_J'][:-1])+abs(B),LD('1e-100'))
 assert np.all(abs(rawP-B)[clean]<=raw_bound[clean])
 assert abs(R.sum()-result['strict_summary']['final']['particle_grid_coupling_gap_J'])<=512*np.finfo(float).eps*max(np.max(abs(m['particle_grid_work_J'])),np.max(abs(m['mesh_mid_J'])))
 D=np.max(abs(m['particle_change_J']));assert D>0
 out=root/case['name'];out.mkdir()
 normalization={'D_J':D,'definition':'max(abs(delta_particle_kinetic_box_J)) over the entire completed stable-cohort history, including removed events','common_to_all_terms_and_both_whole_clean_views':True,'not_initial_total_particle_energy':True}
 save(out/'normalization.json',normalization)
 columns={'step':steps,'fully_released_no_removal':clean.astype(int),**{k:np.cumsum(v,dtype=LD) for k,v in values.items()},'gap':np.cumsum(R,dtype=LD),'particle_change':m['particle_change_J'][1:]}
 for filename,scale in [('work_decomposition_absolute.csv',LD(1)),('work_decomposition.csv',LD(100)/D)]:
  with (out/filename).open('w') as f:
   writer=csv.writer(f);writer.writerow(list(columns))
   for i in range(len(steps)):writer.writerow([int(columns[k][i]) if k in ['step','fully_released_no_removal'] else format(columns[k][i]*scale,'.20e') for k in columns])
 summary={'normalization':normalization,'clean_intervals_steps':intervals(steps[clean]),'clean_steps':int(clean.sum()),'source_particles':N,
  'observed_numerical_pushes':int(acc['observed_pushes'].sum()),'eligible_native_support_push_fraction':acc['eligible_pushes'].sum()/acc['observed_pushes'].sum(),
  'maximum_observed_individual_gamma_box':np.max(acc['max_observed_gamma_box']),
  'maximum_identity_error_J':np.max(abs(identity)),'clean_raw_boris_agreement_max_J':np.max(abs(rawP-B)[clean]),
  'excluded_absolute_boris_work_J':acc['excluded_abs_push_boris_grid_J'].sum(),
  'whole':{**{k:np.sum(v,dtype=LD) for k,v in values.items()},'gap':R.sum()},
  'fully_released_before_any_removal':{**{k:np.sum(v[clean],dtype=LD) for k,v in values.items()},'gap':R[clean].sum()}}
 save(out/'SUMMARY.json',summary);stats[case['name']]=summary;series[case['name']]=(steps,values,R,clean,D)

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'font.family':'DejaVu Sans','axes.unicode_minus':False})
fig,axes=plt.subplots(2,2,figsize=(12.3,8.7),layout='constrained',sharey='row')
for row,kind in enumerate(['undulator','free-drift']):
 key='native-'+kind+'-on';steps,v,R,clean,D=series[key]
 for col,mask in enumerate([np.ones(len(steps),bool),clean]):
  ax=axes[row,col];x=steps[mask];scale=100/D
  curves=[v[terms[0]],v[terms[1]],v[terms[2]],v[terms[3]]+v[terms[4]],v[terms[5]],R]
  labels=['B − G: pusher / displacement','G − S: spatial gather / path','S − C: curved path / chord','(P − B) + (C − M_old): cohort / sources','M_old − M_mid: time pairing','P − M_mid: net gap']
  for y,label,color,style in zip(curves,labels,[colors[0],colors[1],colors[2],colors[3],colors[5],'#25282C'],['-','--',':','-.','-','--']):
   ax.plot(x,np.cumsum(y[mask],dtype=LD)*scale,label=label,color=color,ls=style,lw=1.7)
  ax.axhline(0,color='#777777',lw=.7);ax.set_xlim(x[0],x[-1]);ax.set_xlabel('Maxwell step')
  ax.set_title(('One-period undulator' if row==0 else 'Free drift')+' · '+('whole run' if col==0 else f'all {N:,} released, no exits'))
  if col==0:ax.set_ylabel('Cumulative signed work / fixed D (%)')
  ax.text(.02,.98,f'D = {float(D):.5e} J'+('\nclean interval only; same whole-run D' if col else ''),transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','edgecolor':'none','alpha':.75})
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False,fontsize=9)
fig.suptitle('Exchange-work decomposition · 16 MPI × 8 threads · unchanged physical evolution',fontsize=13)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('work_decomposition.'+ext),dpi=180)
plt.close(fig)

toy=json.loads((a.postcheck/'linear-response/RESULT.json').read_text());save(root/'LINEAR_RESPONSE.json',toy)
fig,axs=plt.subplots(1,2,figsize=(10.6,4.6),layout='constrained')
labels={'explicit_old_field':'Explicit old field','one_predictor_correction':'One predictor correction','linear_implicit_response':'Linear implicit response'}
for mode,color,marker in zip(labels,[colors[0],colors[1],colors[5]],['o','s','^']):
 rows=sorted([x for x in toy['results'] if x['mode']==mode],key=lambda x:x['dt_normalized'])
 for ax,metric in zip(axs,['maximum_absolute_relative_energy_error','final_state_relative_error_vs_exact_linear_ODE']):
  ax.loglog([r['dt_normalized'] for r in rows],[max(r[metric],1e-16) for r in rows],label=labels[mode],color=color,marker=marker,lw=1.7)
for ax in axs:ax.set_xlabel('Normalized time step');ax.set_xticks([.1,.5,1],['0.1','0.5','1']);ax.legend(frameon=False,fontsize=9)
axs[0].set_ylabel('Maximum |ΔH| / H0');axs[1].set_ylabel('Final state error / exact state norm')
axs[0].set_title('Exchange-energy error');axs[1].set_title('Independent exact linear-ODE reference')
fig.suptitle('Classical fixed-weight coupling toy · no spatial motion, Maxwell curl or CPML',fontsize=12)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('classical_linear_response.'+ext),dpi=180)
plt.close(fig)

times={}
for case in plan['cases']:
 result=json.loads((Path(case['directory'])/'RESULT.json').read_text());t=result['progress_elapsed_s'];seconds=(t['240']-t['40'])/200
 times[case['name']]={'progress_40_to_240_seconds_per_step':seconds,'single_run_no_confidence_interval':True}
for kind in ['undulator','free-drift']:
 b=times[f'native-{kind}-baseline']['progress_40_to_240_seconds_per_step'];off=times[f'native-{kind}-off']['progress_40_to_240_seconds_per_step'];on=times[f'native-{kind}-on']['progress_40_to_240_seconds_per_step']
 times[kind+'_native_pair_percent']={'candidate_off_vs_accepted':100*(off/b-1),'diagnostic_on_vs_off':100*(on/off-1)}
save(root/'TIMING.json',times);save(root/'DECOMPOSITION.json',stats);save(root/'SOURCE_HASHES.json',source_hashes)
table=['|工况|完整无退场区间|B−G / aJ|G−S / aJ|S−C / aJ|P−B / aJ|C−Mold / aJ|Mold−Mmid / aJ|净差 / aJ|','|---|---|---:|---:|---:|---:|---:|---:|---:|']
for key,s in stats.items():
 z=s['fully_released_before_any_removal'];table.append('|'+key+'|'+str(s['clean_intervals_steps'])+'|'+'|'.join(f'{float(z[k]*1e18):+.6g}' for k in terms+['gap'])+'|')
timing_table=['|16MPI×8|已验程序 s/步|新程序关闭诊断 s/步|开启诊断 s/步|关闭/原增量|开启/关闭增量|','|---|---:|---:|---:|---:|---:|']
for kind in ['undulator','free-drift']:
 vals=[times[f'native-{kind}-{mode}']['progress_40_to_240_seconds_per_step'] for mode in ['baseline','off','on']];pct=times[kind+'_native_pair_percent'];timing_table.append('|'+kind+'|'+'|'.join(f'{x:.6f}' for x in vals)+f'|{pct["candidate_off_vs_accepted"]:+.2f}%|{pct["diagnostic_on_vs_off"]:+.2f}%|')
readme='''# 场—粒子交换功分解与经典半隐式探索

10个短工况均完整关闭。此轮增加默认关闭的只读观察，原Cowan、Boris、电流沉积、CPML和供场的数值更新保持；尚未启用新的耦合方法。

## 交换差不只来自时间层

定义同一步中：B为选中粒子的原Boris网格功；G为原先保持的取场值沿真实子步位移的功；S为用原生沉积对应空间权重沿子步路径积分的旧场功；C为同样权重沿整步直线弦的旧场功。P为全账本粒子网格功减已记录的约束源功，Mold/Mmid为同一份实际电流分别与旧场/新旧平均场配对。

独立逐项量测满足：

\[
 P-M_{mid}=(B-G)+(G-S)+(S-C)+(P-B)+(C-M_{old})+(M_{old}-M_{mid}).
\]

前三项分别区分Boris电踢与实际位移、原点值取场与轨迹权重、子步折线与整步沉积弦。P−B和C−Mold显式保留全局/子集、终止裁切与约束源口径；不冒充纯空间误差，也不用残差反推。

几何CSV记录的是原数值推送尝试。其空间支撑/释放筛选**不能独立证明没有后续物理退场裁切**，故下表另与v5的累计移除数联用：所有粒子已释放、全部观察完整、当步与前步均无移除。这样避开首个退出事件。完整记录也保留，但它含释放/退出与源项口径。

'''+ '\n'.join(table)+'''

1 aJ=10⁻¹⁸ J，均为boost系。短有磁工况中空间取场／路径项与时间项都有显著贡献并部分抵消，不能据此认为一次中点时间预测即可解决。原生轨迹积分使用独立Gauss求积，并已与原沉积在随机3D场及跨MPI接口上对照。

**只有一周期的短有磁工况，不能把这些占比或符号直接外推到19099的170周期全量。** 这是在原轨迹上解释交换通道，不是修改方法后的因果结果，更不是221.8 keV能损已被解释或修复。

![相对交换功分解](figures/work_decomposition.png)

图及主CSV统一用100×带符号能量/D，D为各工况完整稳定cohort历史的max|ΔKboost|，包括移除事件；同一工况的完整记录和中段使用相同D。原J数据另存。不同工况分母不同，绝对大小请查原J表。相对净差很小不代替对小净交换的验收。

## 经典算法按计算系评估

算法筛选不以lab高γ为排除条件；经典线性响应已进入独立探索。经典电磁学与把粒子动力学中的γ取1需要区分。新诊断记录个体计算系γ，数据在各组SUMMARY；它不是全束速度分布，不能把一个最大值或平均值当整个束团的速度。

固定CIC权重的独立经典子系统采用mv²/2，与矩阵指数精确解比较旧场显式、一次预测校正和直接线性半隐式响应。线性响应在该子系统维持交换能量到舍入量级，同时状态误差随步长减小。它**没有空间运动、Maxwell curl、CPML和电荷连续性验收**；固定矩阵可复用分解，其耗时不能当全PIC成本。

![经典耦合探索](figures/classical_linear_response.png)

半隐式论文、公式、内存和色散成本的研究见[方法分析](SEMI_IMPLICIT.md)。优先考虑保留原波传播、只求解电场—电流—粒子交换的结构；适配Cowan仍待实现和验收，当前默认方法不变。

## 只读观察的耗时与验收范围

'''+ '\n'.join(timing_table)+'''

每项仅同一分配中的一次运行，取40→240步；全部保留原严格账本。波动没有重复样本/置信区间，不据此证明关闭诊断严格零成本，也不外推全量。开启诊断额外交换私有电场接口并做轨迹求积。4MPI×8是接口/布局验证，不是等资源速度比较。

原软件回归、默认/开关逐位比较、冻结校验、启动负控以及失败批次的准确状态见[程序验收记录](PROGRAM_VALIDATION.md)，与物理可用性分开。尚未证明新的耦合算法能保持全量速度或消除能量偏差。
'''
(root/'README.md').write_text(readme)
shutil.copy2(Path(__file__).parent/'SEMI_IMPLICIT.md',root/'SEMI_IMPLICIT.md')
save(root/'ANALYSIS.json',{'checks_passed':True,'cases':list(stats),'physical_coupling_corrected':False,'identity_and_clean_raw_boris_checks_pass':True,'figures_png_rendered_pdf_same_source_not_separately_inspected':True})
print(root)
