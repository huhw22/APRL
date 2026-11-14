"""Read-only provenance audit: independent short runs vs cloned six-segment run."""
import csv
import sys
from pathlib import Path
import h5py
import numpy as np
from common import read,save,sha,need


def inspect(path,expected):
    need(sha(path)==expected,'Original input hash changed')
    with h5py.File(path,'r') as f:a=f['particles/records'][:]
    xy=a['plane_position_m'];_,inverse,count=np.unique(xy,axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(inverse);same=inverse[order][1:]==inverse[order][:-1]
    u=a['proper_velocity'][order]
    copied_momentum=bool(np.array_equal(u[1:][same],u[:-1][same])) if same.any() else False
    return dict(path=str(path),sha256=expected,records=len(a),distinct_transverse_xy=len(count),
        minimum_xy_multiplicity=int(count.min()),maximum_xy_multiplicity=int(count.max()),
        fraction_in_repeated_xy_columns=float(np.mean(count[inverse]>1)),
        same_momenta_inside_repeated_columns=copied_momentum,
        recorded_arrival_span_s=float(np.ptp(a['arrival_time_offset_s'])))


def main(root,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    w=root.parents[1];current=read(root/'PLAN.json')
    short=w/'results/sparc_random_startup_20261112_r2';si=read(short/'w64-a/INPUT_VALIDATED.json')
    old=inspect(Path(current['input_file']),current['input_sha256'])
    new=inspect(Path(si['input']),si['input_sha256'])
    need(old['sha256']!=new['sha256'],'Expected distinct campaigns')
    need(old['fraction_in_repeated_xy_columns']==1 and new['fraction_in_repeated_xy_columns']==0,'Loading evidence differs from campaign classification')
    with (short/'reviewed-r1/ten_cell_summary.csv').open() as f:gen=next(x for x in csv.DictReader(f) if x['case']=='w64-a')
    source=w/'results/puffin_aprl_source_propagation_20261113_r1/physical-route-control-r1/analysis/source_summary.csv'
    with source.open() as f:rr=list(csv.DictReader(f))
    bf={x['program']:float(x['global_bf']) for x in rr if x['plane']=='u07'}
    growth=read(root/'growth-20261114-r1/data/RESULT.json')
    first={x['method']:x for x in growth['rows'] if x['segment']==1}
    answer=dict(complete=True,no_new_simulation=True,current_six_segment_input=old,
        previously_verified_independent_input=new,previous_GENESIS_w64_77periods=gen,
        previous_corrected_Puffin_7periods_global_BF=bf,
        previous_corrected_Puffin_BF_ratio_to_PIC=bf['Puffin']/bf['PIC'],
        current_first_segment_global_BF={k:v['global_b1'] for k,v in first.items()},
        current_first_segment_local10lambda_RMS_BF={k:v['central_10lambda_rms_b1'] for k,v in first.items()},
        current_curves_validate_corrected_independent_loading=False,
        conclusion='The six-segment historical same-input comparison reused the pre-fix phase-cloned input. It does not contradict the independent-random short-run BF result.',
        causal_scope='Loading/provenance mismatch is established. Length, charge, weights, grid and BF definitions also differ; no matched six-segment independent-vs-cloned control exists here to attribute the entire discrepancy quantitatively.')
    save(out/'RESULT.json',answer)
    lines=['# 六段增长图的输入资格复核','',
        '**当前六段图没有使用此前恢复的独立随机束流，而是19743/19744旧账本的冻结相关加载输入。它不能作为“修正后仍然异常快增长”的证据。**','',
        '| 项目 | 此前独立随机首段对照 | 这次六段对照 |','| --- | --- | --- |',
        f"| 原始记录数 | {new['records']} | {old['records']} |",
        f"| 不同横向(x,y)数 | {new['distinct_transverse_xy']} | {old['distinct_transverse_xy']} |",
        f"| 每个横向位置重复次数 | {new['minimum_xy_multiplicity']}…{new['maximum_xy_multiplicity']} | {old['minimum_xy_multiplicity']}…{old['maximum_xy_multiplicity']} |",
        f"| 属于重复横向列的粒子比例 | {100*new['fraction_in_repeated_xy_columns']:.1f}% | {100*old['fraction_in_repeated_xy_columns']:.1f}% |",
        f"| 输入实际到达时间跨度 | {new['recorded_arrival_span_s']*1e12:.6f} ps | {old['recorded_arrival_span_s']*1e12:.6f} ps |",'',
        '旧输入重复列内的动量也逐位相同，沿纵向分布成共线副本；新输入没有精确重复的横向位置。两文件SHA256分别从各自冻结输入清单独立核验。','',
        f"此前独立随机w64-a、77周期固定中央ID群10λ RMS BF：PIC输入{float(gen['PIC_initial_local_b1']):.9f}，出口{float(gen['PIC_exit_local_b1']):.9f}；GENESIS SC-off出口{float(gen['GENESIS-off_exit_local_b1']):.9f}，PIC比GENESIS低约{100*(1-float(gen['PIC_over_GENESIS_off_exit_local_b1'])):.2f}%。当时的结论是旧倍数级异常快启动没有重现，并非所有长度和所有观测量已经严格相等。",'',
        f"后来正确磁路由的同一独立样本7周期Puffin控制：全束global BF PIC={bf['PIC']:.10f}、Puffin={bf['Puffin']:.10f}，相差{100*(bf['Puffin']/bf['PIC']-1):.3f}%。此处使用19770正确磁场控制，不复用已撤销匹配资格的31.47pJ旧磁路由作为基准。",'',
        '这次为了对19743/19744的完整账本作同输入程序比较，Puffin19783也导入了旧克隆束团。因此当前图中的APRL快增长是旧加载条件下的历史响应，不是独立随机加载修复之后的新复测。入口BF一致只约束一个纵向复矩，不能排除高维横纵向相关；相同旧输入分别进入含三维自场的PIC和不含纵向自场的Puffin，沿程仍可分开。','',
        '统计口径也不同：此前GENESIS主对照是固定中央ID群10λ RMS；六段图是后半中央Eulerian窗10λ RMS。当前全束global BF也存在早期差异，因此不能仅靠改BF分组消除六段图里的差别。束长、总电荷、宏权重、网格和传播长度同时不同，不能把两轮曲线直接相减，宣称已经定量分出克隆/自场的贡献。','',
        '本轮只确认输入与结论适用范围，没有重跑物理、换种子、重新生成输入或查询其它运行任务。要验证修正后六段增长，需要新的独立随机六段匹配对照；当前历史账本数据不能代替它。']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    print(answer['conclusion']);print('OLD',old);print('INDEPENDENT',new)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
