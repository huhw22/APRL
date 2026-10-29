#!/usr/bin/env python3
"""Render only cases with closed ledger markers; never open running case H5."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-relative-ledger-mpl')
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
try:
    import report_discrete
except ModuleNotFoundError as error:
    if error.name!='report_discrete':raise
    sys.path.insert(0,str(Path(__file__).resolve().parent.parents[1]/'performance/cases/discrete_energy_ledger'))
    import report as report_discrete
from relative_ledger import construct,write_tables,plot

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--source-job-ended',action='store_true');a=p.parse_args()
    plan=json.loads(a.plan.read_text());root=Path(plan['output']);reports={}
    for case in plan['cases']:
        name=case['name'];marker=Path(case['closed_marker']);ledger=Path(case['ledger'])
        if not marker.exists():
            status='NO_CLOSED_RESULT_AFTER_JOB_END' if a.source_job_ended else 'WAITING_FOR_CLOSED_RESULT_MARKER'
            reports[name]={'status':status,'marker':str(marker)};continue
        gate=json.loads(marker.read_text())
        if not gate.get(case['gate_key'],False):
            reports[name]={'status':'SKIPPED_UNSUCCESSFUL_OR_INCOMPLETE_SOURCE','marker':str(marker)};continue
        # Source completion marker must predate any H5 opening here.
        rows,sums,ranks,steps=report_discrete.load(ledger,ledger/'energy-ledger.h5')
        index,q,summary=report_discrete.calculate(rows,sums,steps)
        values,pct,meta=construct(rows,sums,index,q)
        target=root/name
        fingerprint={str(path):sha(path) for path in [marker,ledger/'energy-ledger.h5',*sorted(ledger.glob('discrete-field-rank-*'))]}
        old=target/'SOURCE_SHA256.json'
        if old.exists():
            assert json.loads(old.read_text())==fingerprint,'Completed source changed: '+name
        meta['source_ledger']=str(ledger);meta['source_closed_marker']=str(marker);meta['mpi_ranks']=ranks
        write_tables(target,index,values,pct,meta);plot(target,index,pct,meta,case.get('plot_title',case['title']))
        (target/'SOURCE_SHA256.json').write_text(json.dumps(fingerprint,indent=2,sort_keys=True)+'\n')
        # Invert the percentages to verify no change in the reported energies.
        den=np.longdouble(meta['fixed_denominator_decimal_J'])
        error=max(float(np.max(abs(pct[k.removesuffix('_J')+'_pct']/100*den-v))) for k,v in values.items())
        assert all(sha(path)==digest for path,digest in fingerprint.items())
        reports[name]={'status':'RELATIVE_PLOT_READY','fixed_denominator_J':meta['fixed_denominator_J'],
          'final_total_difference_percent':meta['final_total_residual_pct_fixed_scale'],
          'maximum_total_difference_percent':meta['maximum_absolute_total_residual_pct_fixed_scale'],
          'normalization_roundtrip_maximum_abs_error_J':error,'raw_inputs_unchanged':True}
    root.mkdir(parents=True,exist_ok=True)
    (root/'STATUS.json').write_text(json.dumps(reports,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    lines=['# 相对能量收支图与数据','',
      '图和主CSV使用百分数。含粒子工况统一以整段已完成记录中 `D=max|ΔK_particle|` 为固定分母；所有曲线共用同一个D。',
      'ΔK为boost系同一粒子集合（含累计移除事件）的稳定动能差。D不是束流初始总动能，也不是每条曲线各自的最大值。',
      '没有粒子动能变化时，依次使用初始、峰值原生场能；全零工况的相对量未定义。',
      '分母为整个显示区间的固定数值，过零处不作逐时刻相除。符号保持；百分数与J均可追溯。',
      '每组 `energy_budget.csv` 是百分数表，`energy_budget_absolute.csv` 保存对应J值，`normalization.json`说明分母。',
      '正常化只改变表达，不改变物理结论；时间离散与Cowan项仍保留为数值通道。','']
    for case in plan['cases']:
        name=case['name'];report=reports[name]
        lines+=['## '+case['title'],'']
        if report['status']=='RELATIVE_PLOT_READY':
            lines += [f"终点总差额为 **{report['final_total_difference_percent']:+.6g}%**，全程最大绝对差额为 **{report['maximum_total_difference_percent']:.6g}%**，均相对于该工况的D。",'',
              f'![相对能量收支]({name}/energy_budget.png)','',
              f'[百分数数据]({name}/energy_budget.csv) · [分母说明]({name}/normalization.json) · [PDF]({name}/energy_budget.pdf)','']
        elif report['status']=='WAITING_FOR_CLOSED_RESULT_MARKER':
            lines += ['尚未出现可用的已关闭结果标记；未读取该工况运行中的H5。','']
        else:lines += ['原作业已结束或来源检查未成功；没有可用的完整关闭结果，本项未绘图、未读取其H5。','']
    (root/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
