"""Review closed short controls and static provenance; never read running H5."""
import argparse,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-cowan-matched-mpl')
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import sha,save,require
def closed_physics(r):
 out={};root=r/'physics_r1'
 for name in ['uniform','converging','gaussian']:
  if not (root/(name+'_COMPLETE.json')).exists():continue
  require((root/(name+'_EXIT_CODE')).read_text().strip()=='0','Physical driver failed '+name);d=json.loads((root/name/'RESULT.json').read_text());s=d['physics_summary'];last=s['planes'][-1];old=Path(json.loads((root/'plan.json').read_text())['cases'][name]['original']);prior=json.loads((old/'analysis/RESULT.json').read_text())['planes'][-1];require(all(d['original_first_force_state_bitwise'].values()),'Release state changed');require(d['execution_integrity_pass'],'Case integrity failed');out[name]=dict(execution_integrity_pass=True,release_state_bitwise_identical=True,particles_per_plane=s['particles_per_plane'],steps=s['steps'],last_plane=last,previous_last_plane=prior,last_mean_rest_gain_error_pct=100*(last['PIC_mean_rest_gain_eV']-last['reference_mean_rest_gain_eV'])/abs(last['reference_mean_rest_gain_eV']),last_mean_lab_gain_error_eV=last['PIC_mean_lab_gain_eV']-last['reference_mean_lab_gain_eV'],last_std_rest_gain_error_pct=100*(last['PIC_std_rest_gain_eV']/last['reference_std_rest_gain_eV']-1),previous_std_rest_gain_error_pct=100*(prior['PIC_std_rest_gain_eV']/prior['reference_std_rest_gain_eV']-1),last_sigma_x_error_pct=100*(last['PIC_sigma_x_m']/last['reference_sigma_x_m']-1),old_unweighted_ledger_maximum_residual_percent_D=s['relative_energy_budget']['maximum_absolute_total_residual_pct_fixed_scale'],old_unweighted_ledger_final_residual_J=s['new_total_final_residual_J'],independently_measured_particle_grid_work_difference_J=s['particle_grid_exchange_final_J'],field_operator_residual_max_J=s['new_field_operator_maximum_residual_J'],physical_accuracy_accepted=False)
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();r=w/'results/cowan_matched_coupling_20261102';b=r/'build_r3';z=json.loads((b/'RESULT.json').read_text());checks={}
 for name in ['build_r1','build_r2','build_r3','physics_r1','mpi_r1','flow_r1']:
  d=r/name;f=d/'FROZEN.json'
  if f.exists():checks[name+'_frozen']=all(sha(k)==v for k,v in json.loads(f.read_text()).items());require(checks[name+'_frozen'],'Frozen mismatch '+name)
 original=w/'results/lab_energy_balance_20261102/coupled_sbp_r2/source/src';new=r/'build_r2/source/src';unchanged=[]
 for f in original.iterdir():
  if f.suffix in ['.cpp','.h'] and f.name not in ['source_coupling.cpp','source_coupling_math.h']:
   require(sha(f)==sha(new/f.name),'Unapproved core change '+f.name);unchanged.append(f.name)
 checks['other_core_files_identical']=unchanged
 doc_record=w/'results/lab_energy_balance_20261102/yee_sbp_r1/PREEXISTING_DOC_HASHES.json';checks['preexisting_docs_unchanged']=all(sha(w/'source/aprl'/k)==v for k,v in json.loads(doc_record.read_text()).items());require(checks['preexisting_docs_unchanged'],'Preexisting documents changed')
 checks['production_program_unchanged']=sha(w/'build/source-coupling-exit-r3-20261031/aprl')=='865ebef5b79b6b8670d40c81080f15f916ab1938033afcabd050bcd3c8198752';require(checks['production_program_unchanged'],'Production program changed')
 cards={};plan=json.loads((r/'physics_r1/plan.json').read_text())
 for name,case in plan['cases'].items():
  d=Path(case['directory']);old=Path(case['original']);now=yaml.safe_load((d/'run.yaml').read_text().replace(str(d),str(old)));before=yaml.safe_load((old/'run.yaml').read_text());require(now==before,'Scientific card changed '+name);cards[name]=dict(scientific_card_identical_except_outputs=True,original=str(old),expected_events_sha256=sha(d/'reference_input/expected_events.npz'),steps=case['steps'])
 physics=closed_physics(r)
 fig,axes=plt.subplots(1,3,figsize=(13.8,4.2),constrained_layout=True);metrics=[];Gamma=np.longdouble('929.1305492858543');beta=np.sqrt(1-1/Gamma**2)
 for ax,item in zip(axes,[x for x in z['periodic'] if x['wave']=='oblique']):
  q=item['quiet_refinement'];data=np.genfromtxt(b/f'periodic-oblique-{q}.csv',delimiter=',',names=True,dtype=np.longdouble);D=np.max(np.abs(data['particle_delta_Kprime_J']));require(D>0,'Zero normalization');step=data['step'];rk=data['matched_energy_residual_J'];rp=data['particle_delta_cPzprime_J']+data['matched_field_cPz_change_J'];lab=data['matched_lab_residual_J'];require(np.max(np.abs(lab/Gamma-rk-beta*rp))<1e-14*D,'Lorentz residual reconstruction')
  for y,label,color in [(rk,'Matched boost energy','#0072B2'),(rp,'Longitudinal cP sum','#D55E00'),(lab/Gamma,'Lab combination / Gamma','#009E73'),(data['particle_delta_Kprime_J']+data['native_field_change_J'],'Unweighted synchronized energy','#CC79A7')]:ax.plot(step,100*y/D,label=label,color=color,linewidth=1.3)
  ax.set(title=f'{item["particles"]:,} macro-particles',xlabel='Maxwell step');ax.set_yscale('symlog',linthresh=1e-7);ax.grid(alpha=.2);metrics.append(dict(wave='oblique',particles=item['particles'],quiet_refinement=q,normalization_D_boost_particle_J=float(D),normalization='fixed max abs particle boost kinetic change over complete record',matched_boost_residual_max_percent_D=float(100*np.max(np.abs(rk))/D),momentum_residual_max_percent_D=float(100*np.max(np.abs(rp))/D),lab_residual_max_percent_particle_lab_change=100*item['lab_residual_over_particle_change_scale'],electric_impulse_gap_J=item['electric_impulse_gap_J'],magnetic_impulse_gap_J=item['magnetic_impulse_gap_J'],physical_accuracy_accepted=False))
 axes[0].set_ylabel("Signed % of fixed max |particle Delta K'|");handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,fontsize=9,loc='lower center',bbox_to_anchor=(.5,-.13),ncol=2);fig.suptitle('Closed oblique-wave / pair-plasma control: same grid, particle sampling changes',fontsize=12);out=r/'figures';out.mkdir(exist_ok=True);fig.savefig(out/'periodic_relative.png',dpi=170,bbox_inches='tight');fig.savefig(out/'periodic_relative.pdf',bbox_inches='tight');plt.close(fig)
 mpi=None
 if (r/'mpi_r1/EXIT_CODE').exists() and (r/'mpi_r1/EXIT_CODE').read_text().strip()=='0':
  mpi=json.loads((r/'mpi_r1/checks.log').read_text());require(mpi['ranks']==4 and mpi['weighted_work_relative_max']<5e-12 and mpi['weighted_magnetic_impulse_relative_max']<5e-12,'Cross-node identity failed')
 result=dict(status='CANDIDATE_IMPLEMENTED_PHYSICAL_ACCEPTANCE_PENDING',validation=checks,space_charge_cards=cards,closed_space_charge_results=physics,periodic=metrics,cross_node_implementation_check=mpi,physical_acceptance=False,production_replaced=False,all_periodic_loadings_accepted=False,known_failure='The 192-macro-particle loading has large electric-impulse and lab residuals; increasing sampling at fixed grid is not a general convergence proof',full_moving_boundary_matched_metric_ledger_implemented=False,production_cost_measured=False);save(r/'REVIEW.json',result)
 lines=['# Cowan 配套自场候选：2026-11-02','', '已实现局部自电场/自磁场取样候选。原 Cowan、粒子动量推进、沉积、CPML、参考供场及输入释放相空间保持；正式程序与默认入口未替换。**尚未宣布物理验收通过。**','', '## 已闭合短控制','', '程序验证：1/2/4 rank 的真实沉积—取场功/磁冲量配对通过，最大相对误差分别约 7.7e-15 / 1.13e-14。跨节点结果见 REVIEW.json；它属于程序通信验证。','', '物理控制：相同生产网格间距、dt 和 boost 的闭合周期带电粒子盒，160 步。正负粒子初始重合，非零电流、场和粒子能量交换。此处是验证控制，不能直接代表实际电子束。轴向对称例与具有 Ez 的斜向例分开。','', '|斜向例宏粒子数|配套 boost 总能差峰值 / D|lab 组合差峰值 / 实际粒子 lab 能量变化尺度|','|---:|---:|---:|']
 for x in metrics:lines.append(f'|{x["particles"]:,}|{x["matched_boost_residual_max_percent_D"]:.6g}%|{x["lab_residual_max_percent_particle_lab_change"]:.6g}%|')
 lines+=['', 'D=max|完整记录粒子 boost 动能变化|，每例固定一个分母。lab 列使用该系实际粒子能量变化的峰值为分母。各图面板纵轴范围不同，图内所有曲线共同除以该例固定 D；lab 组合先除以 Gamma 再与 boost 量比较。保留带符号原始 J 数据，不以初始总束流动能稀释差额。','', '![相对控制](figures/periodic_relative.png)','', '粗加载例的明显差额保留，并已独立定位到电力冲量配对。加密粒子加载后小得多；这不等于网格/时间收敛，更不保证误差已经低于目标辐射能损。磁冲量差和场动量恒等式另测，不能把前者闭合作为总动量闭合。紫线为未加权、时间同步的原生能量收支，特意保留以显示不同能量定义的差别。','', '## 正在进行的原物理对照','', '|作业|工况|日志|','|---:|---|---|','|19239|冷均匀球|[aprl.log](physics_r1/uniform/aprl.log)|','|19240|会聚均匀球|[aprl.log](physics_r1/converging/aprl.log)|','|19241|冷高斯|[aprl.log](physics_r1/gaussian/aprl.log)|','', '三组沿用原输入、电荷、理论事件、网格和停止条件；单节点16 MPI×8线程，每例2700秒硬限、作业50分钟硬限，不固定节点。仍开原严格账本并自动在闭合后与理论比较。`physics_r1/plan.json` 和 `FROZEN.json` 保存来源及身份。运行中的 H5 不读取；启动状态仅来自日志、静态卡和绑定记录。','', '新度量在完整移动边界上的账本尚未接入；本轮不能用周期泛函闭合替代 CPML、移窗和供场的全程验收。全量4.5m未提交，成本配对未测。','', '## 证据和实现','', '* [短控制结果](build_r3/RESULT.json)、[逐步CSV所在目录](build_r3/)、[源码构建身份](build_r2/SOURCE_BUILD.json)。','* [空间电荷计划](physics_r1/plan.json)、[静态复核](REVIEW.json)。','* [实现说明](../../source/aprl/verification/performance/cases/cowan_matched_coupling/README.md)、[算子推导与限制](../../source/aprl/verification/performance/cases/cowan_matched_coupling/METHOD.md)。','', '旧证据、生产程序和四个此前已修改的理论文档保持。build_r1/r2 的粗加载及后续细化证据均保留；build_r3 增加斜向控制。PNG需目视核查，PDF与PNG同源导出。']
 if len(physics)==3:
  lines=[x.replace('## 正在进行的原物理对照','## 原空间电荷物理对照（均已正常完成）') for x in lines]
  lines+=['','## 空间电荷结果与尚未闭合的口径','', '|工况|末端静止系平均能量变化 / 理论|均值误差|能散误差（旧→新）|末端横向rms束宽误差|原未加权账本差额峰值 / D|','|---|---:|---:|---:|---:|---:|']
  for name,v in physics.items():
   last=v['last_plane'];lines.append(f'|{name}|{last["PIC_mean_rest_gain_eV"]:.4f} / {last["reference_mean_rest_gain_eV"]:.4f} eV|{v["last_mean_rest_gain_error_pct"]:+.4f}%|{v["previous_std_rest_gain_error_pct"]:+.3f}% → {v["last_std_rest_gain_error_pct"]:+.3f}%|{v["last_sigma_x_error_pct"]:+.5f}%|{v["old_unweighted_ledger_maximum_residual_percent_D"]:.6f}%|')
  lines+=['','三组32768个粒子的首次释放状态逐字段逐位相同，全部粒子完整过面；原膨胀、会聚和回弹过程保留。理论为原独立球壳世界线，在各粒子对应事件上比较；旧版与新版到达事件略有差别时，查询理论的时刻也略有差别，没有重新拟合理论。结果支持这次取场改动保留基本动力学趋势；均匀球能散相对连续理论的偏差仍约9.15%（旧版10.46%），并非所有统计量都在1%内，不构成辐射精度验收。','', '**原未加权账本出现新增的功差，不能把它写成能量已经守恒。** 它独立记录 raw J·E，而粒子现在接收轨道伴随的 J·ME E。该差与独立测得的场粒子功差几乎一致，场算子自身残差仍约1e-24 J。配套 H_C 在周期盒已闭合，但完整移窗、供场和CPML的 H_C 账本尚未接入；在补齐并核对前保留全部原差额，不改称外场功或真实辐射。','', '* [冷均匀球理论对照](physics_r1/uniform/analysis/comparison.png)','* [会聚均匀球理论对照](physics_r1/converging/analysis/comparison.png)','* [冷高斯理论对照](physics_r1/gaussian/analysis/comparison.png)','', '## 双节点短全流程复核','', '[flow_r1](flow_r1/) 作业19245：短波荡器和自由漂移各257步，两节点、总16 MPI×8线程，保留原分片数；每例600秒限制、总25分钟。确认收敛和通信、读取原参考口径的分项账本，**不是**宣称新配套度量的完整边界账本已完成。日志为 [有磁](flow_r1/undulator-on/aprl.log)、[无磁](flow_r1/drift-on/aprl.log)。']
 (r/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps({k:result[k] for k in ['status','physical_acceptance','cross_node_implementation_check']},ensure_ascii=False))
if __name__=='__main__':main()
