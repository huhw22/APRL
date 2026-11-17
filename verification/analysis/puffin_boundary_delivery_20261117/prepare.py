"""Reuse closed native outputs and independently test the wrap-around interpretation."""
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from common import SOURCE,OLD,sha,save,table

def background(root):
    root=Path(root);assert not (root/'PACKAGE.json').exists()
    bg=root/'background';bg.mkdir(exist_ok=True)
    old=json.loads((OLD/'RESULT.json').read_text())
    save(bg/'SOURCE_OPERATOR_RESULT.json',{k:old[k] for k in ['grid','modal_factors','tail_mask','validations','sources']})
    shutil.copytree(OLD/'evidence',bg/'evidence',dirs_exist_ok=True)
    dest=bg/'data';dest.mkdir(exist_ok=True)
    for name in ['native_transverse_mask.csv','native_tail_masks.csv','native_tail_actual_rank.csv','tail_mask_comparison.csv','modal_factors.csv']:
        shutil.copyfile(OLD/'data'/name,dest/name)
    (dest/'README.md').write_text('''# 源码背景数据

原生掩膜函数的横向/纵向数值、只更改几何长度括号的参考，以及边缘支路谱系数；它们不是本次传播的反射率。实际19834传播表在根目录data中。旧规定位置扫描未作为新版素材重复打包。
''')
    stems=['01_transverse_mask','02_mask_profile','03_modal_factors','04_longitudinal_length_issue','05_longitudinal_grid_sensitivity']
    for lang in ['zh-CN','en']:
        target=bg/'figures'/lang;target.mkdir(parents=True,exist_ok=True)
        for stem in stems:
            for suffix in ['.png','.pdf']:shutil.copyfile(OLD/'figures'/lang/(stem+suffix),target/(stem+suffix))
        for stem in ['06_field_vs_energy','07_phase_example','08_prescribed_input','09_single_boundary_output']:
            for suffix in ['.png','.pdf']:
                extra=target/(stem+suffix)
                if extra.exists():extra.unlink()
    for name in ['previous_source_README.md','previous_source_RESULT.json']:
        extra=bg/name
        if extra.exists():extra.unlink()
    extra=dest/'prescribed_position_scan.csv'
    if extra.exists():extra.unlink()
    (bg/'README.md').write_text('''# 源码背景与本次传播测试的关系

PPT主素材是[实际传播说明](../README.md)。这里仅保留理解边界实现所需的源码资料，不再重复旧版规定输入位置扫描或单次算子演示。

## 横向边界算法

Puffin将边缘场取出，进行依赖波数和频率的谱滤波，再与保留部分相加。掩膜M表示进入滤波分支的场权重，不是反射率。输入项sBeta是吸收参数，不是Lorentz速度。相关[实际源码行](evidence/diffraction.f90.txt)与[掩膜函数行](evidence/masks.f90.txt)保留来源SHA。

可用[横向掩膜](figures/zh-CN/01_transverse_mask.png)、[掩膜剖面](figures/zh-CN/02_mask_profile.png)及[支路系数](figures/zh-CN/03_modal_factors.png)作附录。支路D或D²不能直接解释为整束损失或反射率。

## 纵向长度表达式是另一项问题

既有原生库函数查询确认纵向mask总长度写作nx−1×dx，与坐标几何(nx−1)×dx不同。对应[曲线](figures/zh-CN/04_longitudinal_length_issue.png)和[总网格长度敏感性](figures/zh-CN/05_longitudinal_grid_sensitivity.png)仅作独立源码背景。只改括号的几何参考不是已验收的边界修复。

**19834采用原生单波长周期纵向分支，纵向末端mask被置零。因此这处纵向表达式没有参与本次0.457%中央返回结果。** 本次测试隔离了横向边界及FFT有限数组绕回，不能用它验证纵向mask修复收益。

[源码数值摘要](SOURCE_OPERATOR_RESULT.json) · [背景数据](data/README.md)。旧版完整说明和位置扫描在原r1交付中保留，新版主讲不使用它们替代实际传播。
''')

def main(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False)
    data=root/'data';data.mkdir();meta=root/'metadata';meta.mkdir()
    result=json.loads((SOURCE/'analysis/RESULT.json').read_text())
    assert result['complete']
    plan=json.loads((SOURCE/'PLAN.json').read_text())
    closed=json.loads((SOURCE/'PHYSICS_COMPLETE.json').read_text());assert closed['complete']
    identity=json.loads((SOURCE/'analysis/PROTECTED_OUTPUTS.json').read_text())
    protected={}
    for path,s in identity.items():
        p=Path(path);st=p.stat()
        assert [st.st_ino,st.st_size,st.st_mtime_ns]==[s['inode'],s['size'],s['mtime_ns']],path
        protected[path]=dict(**s,sha256=sha(p))
    for path,h in json.loads((SOURCE/'FROZEN.json').read_text()).items():assert sha(path)==h,path
    for path,h in plan['source_hashes'].items():assert sha(path)==h,path
    for f in ['comparison.csv','history.csv','RESULT.json','visual_arrays.npz']:
        shutil.copyfile(SOURCE/'analysis'/f,data/f)
    for f in ['PLAN.json','PHYSICS_COMPLETE.json','ANALYSIS_COMPLETE.json','FROZEN.json','SOURCE_COMMIT.json']:
        shutil.copyfile(SOURCE/f,meta/('SOURCE_'+f))
    for c in plan['cases']:
        d=SOURCE/'cases'/c['name'];m=json.loads((d/'META.json').read_text());q=json.loads((d/'COMPLETE.json').read_text())
        assert q['steps']==c['steps'] and json.loads((d/'LAUNCH.json').read_text())['returncode']==0
        save(meta/(c['name']+'.json'),dict(configuration=c,meta=m,completion=q,
             elapsed_s=json.loads((d/'LAUNCH.json').read_text())['elapsed_s']))
    a=np.genfromtxt(data/'comparison.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    arrays=np.load(data/'visual_arrays.npz');x=arrays['x_m'];zs=arrays['snapshot_z_m'];X,Y=np.meshgrid(x,x)
    ref=arrays['reference'];noabs=arrays['no_absorption'];native=arrays['native']
    dx=40e-6;u0=json.loads((SOURCE/'cases/wide-reference/COMPLETE.json').read_text())['initial_norm2_m2']
    qmeta=json.loads((SOURCE/'cases/wide-reference/META.json').read_text());theta=qmeta['theta_x_rad']
    sig=200e-6;kx=30000.;k=kx/theta;length=81*dx
    # Independent periodized analytic Gaussian. This is a reference calculation,
    # not another field solver and not a fitted/rephased numerical field.
    periodized=[]
    for z in zs:
        q=1+1j*z/(2*k*sig**2);sx=np.zeros(81,complex);sy=np.zeros(81,complex)
        for m in range(-4,5):
            xx=x+m*length
            sx+=np.exp(1j*kx*xx-1j*kx*kx*z/(2*k))*np.exp(-(xx+300e-6-theta*z)**2/(4*sig**2*q))
            yy=x+m*length;sy+=np.exp(-yy*yy/(4*sig**2*q))
        periodized.append(sy[:,None]*sx[None,:]/q)
    periodized=np.array(periodized)
    err=np.sqrt(np.sum(abs(noabs-periodized)**2,axis=(1,2))*dx**2/u0)
    # Finite initial tails outside the stored small grid account for the small
    # initial discrepancy; the total discrepancy remains constant under propagation.
    assert float(np.max(err))<1e-5
    drift=float(np.max(abs(err-err[0])));assert drift<1e-10
    table(data/'periodic_image_check.csv',[dict(z_m=z,field_difference_over_input=v) for z,v in zip(zs,err)])
    summary={}
    for name in ['small-native','small-no-absorption']:
        q=a[a['case']==name];late=q['z_m']>=1.3
        i=int(np.argmax(np.where(late,q['roi_norm2_over_initial'],-1)))
        thresholds=[]
        for threshold in [.001,.01]:
            inds=np.flatnonzero(q['roi_field_L2_difference_over_initial']>=threshold)
            if len(inds):
                j=int(inds[0]);thresholds.append(dict(field_fraction=threshold,
                    first_saved_z_m=float(q['z_m'][j]),previous_saved_z_m=float(q['z_m'][max(0,j-1)])))
        summary[name]=dict(**result['metrics'][name],central_return_peak_z_m=float(q['z_m'][i]),
              reference_central_at_return_peak=float(q['roi_reference_norm2_over_initial'][i]),
              central_field_difference_at_return_peak=float(q['roi_field_L2_difference_over_initial'][i]),
              first_saved_threshold_crossings=thresholds)
    summary['reference']=dict(**result['validation'])
    summary['periodic_image_validation']=dict(max_field_difference_over_input=float(err.max()),
       initial_field_difference_over_input=float(err[0]),drift_from_initial=float(drift),image_orders=list(range(-4,5)),
       interpretation='Closed no-absorption fields agree with the independently periodized analytic beam to the finite-input-tail discrepancy; this supports periodic wrap-around.')
    summary['geometry']=dict(dx_um=40,small_nodes=[81,81],fft_period_um=3240,
       cell_boundary_um=1620,mask_inner_nodes_um=1000,roi_half_width_um=520,
       nominal_center_enter_layer_z_m=(.001+.0003)/theta,
       nominal_center_cross_period_boundary_z_m=(.00162+.0003)/theta,
       nominal_center_enter_ROI_after_wrap_z_m=(.00324-.00052+.0003)/theta,
       note='Geometrical center markers, not the first arrival of the finite-width Gaussian tails.')
    summary['scope']=dict(real_native_propagation=True,transverse_only=True,longitudinal_mask_tested=False,
        particles=False,FEL_gain_tested=False,statistic_denominator='Initial full-domain integrated intensity of the common wide-domain beam.')
    save(root/'SUMMARY.json',summary)
    save(meta/'PROTECTED_SOURCES.json',protected)
    save(meta/'SOURCE_SHA256.json',{str(p):sha(p) for p in [SOURCE/'analysis/RESULT.json',SOURCE/'analysis/visual_arrays.npz',OLD.with_suffix('.zip')]})
    background(root)
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main(sys.argv[1])
