"""Standalone two-method curves from closed six-segment ledger records."""
import csv
import hashlib
import json
import sys
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker

LD = np.longdouble
CASES = [('F3','energy-binomial3'), ('MM','energy-momentum')]


def read(p):
    return json.loads(Path(p).read_text())


def save(p, value):
    Path(p).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''):
            h.update(block)
    return h.hexdigest()


def table(path, data):
    with path.open('w') as f:
        writer=csv.writer(f);writer.writerow(list(data))
        writer.writerows(zip(*data.values()))


def get_case(campaign, root, tag, case):
    source = root/'three-point' if tag=='F3' else campaign/'momentum-ledger-20261116-r1'
    result=read(source/'RESULT.json')
    assert result['complete'] and result['case']==case
    with np.load(source/'boost/cache.npz') as f:
        a=f['records']
    assert len(a)==28200 and np.array_equal(a['step'],np.arange(28200))
    with np.load(source/'whole'/case/'whole_cache.npz') as z:
        cum=lambda k:np.r_[LD(0),np.cumsum(z[k],dtype=LD)]
        WJ,PJ,GE=cum('fixed_end_WJ_particle_J'),cum('fixed_end_PJ_particle_J'),cum('fixed_end_GaussE_J')
        stored=z['particle_vs_deposited_current_Lorentz_gap_J']
        complete_sum=z['full_window_Lorentz_accounting_difference_J']
    G,beta=LD(a['boost_gamma'][0]),LD(a['boost_beta'][0])
    ext=a['prescribed_boris_work_J'].astype(LD)
    energy=G*(a['delta_particle_kinetic_box_J'].astype(LD)-ext-WJ)
    current=G*beta*(a['delta_particle_cPz_box_J'].astype(LD)+ext/beta+PJ)
    force=current-G*beta*GE
    reconstructed=energy+current
    assert np.max(abs(reconstructed-stored))<LD('2e-18')
    assert np.max(abs(reconstructed-complete_sum))<LD('2e-18')
    residual=complete_sum.copy()
    independent_energy=G*(a['particle_grid_work_J'].astype(LD)-WJ)
    assert np.max(abs(independent_energy-energy))<LD('1e-15')
    boost=np.genfromtxt(source/'boost/energy_budget_absolute.csv',delimiter=',',names=True,dtype=LD)
    lab=np.genfromtxt(source/'lab'/case/'lab_energy_absolute.csv',delimiter=',',names=True,dtype=LD)
    data=dict(step=a['step'],time_box_s=a['time_box_s'],
        energy_pairing_J=energy,longitudinal_current_pairing_J=current,
        longitudinal_force_pairing_J=force,lab_on_boost_slice_sum_J=residual,
        boost_discrete_sum_J=boost['total_residual_J'])
    table(root/'data'/f'{tag}_time_curves.csv',data)
    shutil.copy2(source/'lab'/case/'lab_energy_absolute.csv',root/'data'/f'{tag}_lab_position_curves.csv')
    shutil.copy2(source/'boost/energy_budget_absolute.csv',root/'data'/f'{tag}_boost_components.csv')
    shutil.copy2(source/'whole'/case/'whole_lab_slice.csv',root/'data'/f'{tag}_lab_slice_components.csv')
    complete=np.isfinite(lab['physical_lab_sum_J'])
    denominator=max(abs(lab['particle_net_transport_J'][complete]))
    summary=dict(physics_job=result['physics_job'],case=case,source_directory=str(source),
        final_energy_pairing_J=float(energy[-1]),
        final_longitudinal_current_pairing_J=float(current[-1]),
        final_longitudinal_force_pairing_J=float(force[-1]),
        final_lab_on_boost_slice_sum_J=float(residual[-1]),
        final_boost_discrete_sum_J=float(boost['total_residual_J'][-1]),
        max_absolute_boost_discrete_sum_J=float(max(abs(boost['total_residual_J']))),
        final_fixed_lab_sum_J=float(lab['physical_lab_sum_J'][-1]) if complete[-1] else None,
        final_fixed_lab_particle_net_J=float(lab['particle_net_transport_J'][-1]),
        fixed_lab_fixed_denominator_J=float(denominator),
        final_fixed_lab_sum_percent_fixed_D=float(100*lab['physical_lab_sum_J'][-1]/denominator) if complete[-1] else None,
        all_six_fixed_lab_budgets_complete=bool(complete.all()),
        maximum_pairing_split_vs_stored_J=float(np.max(abs(reconstructed-stored))),
        independent_energy_check_max_J=float(np.max(abs(energy-independent_energy))))
    protected={str(p):sha(p) for p in [source/'RESULT.json',source/'boost/cache.npz',
        source/'whole'/case/'whole_cache.npz',source/'lab'/case/'lab_energy_absolute.csv']}
    return data,boost,lab,summary,protected


def main(root,font):
    root,font=Path(root).resolve(),Path(font).resolve()
    campaign=root.parent
    assert (root/'three-point/RESULT.json').exists()
    (root/'data').mkdir(exist_ok=False)
    cases={};summaries={};sources={}
    for tag,case in CASES:
        data,boost,lab,summary,protected=get_case(campaign,root,tag,case)
        cases[tag]=(data,boost,lab);summaries[tag]=summary;sources.update(protected)
    assert np.array_equal(cases['F3'][0]['time_box_s'],cases['MM'][0]['time_box_s'])
    font_manager.fontManager.addfont(str(font))
    family=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.family':[family,'DejaVu Sans'],'font.size':12,
        'axes.titlesize':16,'axes.labelsize':13,'axes.spines.top':False,'axes.spines.right':False,
        'text.parse_math':False,'axes.unicode_minus':False,'pdf.fonttype':42,
        'figure.facecolor':'white','savefig.facecolor':'white'})
    records=[]
    def plot(language,name,series,title,xlabel,ylabel):
        target=root/'figures'/language;target.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
        lines=[]
        for i,(label,x,y,color,style) in enumerate(series):
            xx,yy=np.asarray(x,float),np.asarray(y,float);good=np.isfinite(xx)&np.isfinite(yy)
            line,=ax.plot(xx[good],yy[good],label=label,color=color,linestyle=style,lw=1.8,
                         marker='o' if len(xx)<20 else None,markersize=5)
            assert np.array_equal(line.get_ydata(),yy[good])
            lines.append(dict(label=label,x=xx[good].tolist(),y=yy[good].tolist()))
        ax.axhline(0,color='#888888',lw=.7);ax.grid(alpha=.2)
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos:f'{y:g}'))
        ax.legend(frameon=False,ncol=2 if len(series)>4 else 1,fontsize=10)
        fig.canvas.draw();renderer=fig.canvas.get_renderer();hidden=set()
        for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            low,high=sorted(limits)
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not low<=tick.get_loc()<=high:hidden.update([id(tick.label1),id(tick.label2)])
        for text in fig.findobj(matplotlib.text.Text):
            if not text.get_visible() or not text.get_text() or id(text) in hidden:continue
            assert '$' not in text.get_text() and '_' not in text.get_text()
            box=text.get_window_extent(renderer)
            assert box.x0>=-2 and box.y0>=-2 and box.x1<=fig.bbox.width+2 and box.y1<=fig.bbox.height+2,text.get_text()
        for suffix in ['png','pdf']:fig.savefig(target/(name+'.'+suffix),dpi=180)
        plt.close(fig);records.append(dict(language=language,name=name,axes=1,series=lines))
    for language in ['zh-CN','en']:
        zh=language=='zh-CN'
        labels={'F3':'三点核' if zh else 'EC + three-point kernel',
                'MM':'动量匹配（五次补偿）' if zh else 'EC + compensated quintic'}
        boxx='计算系模拟时间（ps）' if zh else 'Computational-frame time (ps)'
        labx='实验室系探测面位置（m）' if zh else 'Laboratory detector position (m)'
        micro='带符号能量（µJ）' if zh else 'Signed energy (µJ)'
        specs=[
            ('energy_pairing','energy_pairing_J',1e9,'电功配对差的变化' if zh else 'Electric-work pairing gap',
             '投影后的能量（nJ）' if zh else 'Projected energy (nJ)'),
            ('longitudinal_current_pairing','longitudinal_current_pairing_J',1e6,
             '粒子与电流的纵向配对差' if zh else 'Particle-current longitudinal pairing gap',micro),
            ('longitudinal_force_pairing','longitudinal_force_pairing_J',1e6,
             '计入原生电荷电力后的纵向配对差' if zh else 'Longitudinal pairing gap including native charge force',micro),
            ('lab_slice_sum','lab_on_boost_slice_sum_J',1e6,
             '计算系等时面上的实验室系账本加和' if zh else 'Laboratory ledger sum on computational-time slices',micro),
            ('boost_discrete_sum','boost_discrete_sum_J',1e18,
             '计算系完整离散能量账本加和' if zh else 'Complete computational discrete energy sum',
             '带符号能量（aJ）' if zh else 'Signed energy (aJ)')]
        for name,key,scale,title,ylabel in specs:
            curves=[(labels[tag],cases[tag][0]['time_box_s']*1e12,cases[tag][0][key]*scale,
                     '#d37b28' if tag=='F3' else '#57916f','-' if tag=='F3' else '--') for tag,_ in CASES]
            plot(language,name,curves,title,boxx,ylabel)
        plot(language,'fixed_lab_sum',[(labels[tag],cases[tag][2]['z_lab_m'],
            cases[tag][2]['physical_lab_sum_J']*1e6,'#d37b28' if tag=='F3' else '#57916f',
            '-' if tag=='F3' else '--') for tag,_ in CASES],
            '固定实验室观测管能量账本加和' if zh else 'Fixed laboratory worldtube energy sum',labx,micro)
        for tag,_ in CASES:
            data,boost,lab=cases[tag]
            curves=[]
            for key,cn,en,color,style in [
                ('particle_change_J','粒子能量变化','Particle energy change','#2679ae','-'),
                ('field_boundary_source_prediction_J','场与边界独立预测','Independent field and boundary prediction','#d37b28','--'),
                ('total_residual_J','全部离散分项加和','Sum of all discrete terms','#57916f','-.')]:
                curves.append((cn if zh else en,data['time_box_s']*1e12,boost[key]*1e6,color,style))
            plot(language,tag+'_boost_balance',curves,
                 labels[tag]+('：计算系能量收支' if zh else ': computational energy balance'),boxx,micro)
            curves=[]
            for key,cn,en,color,style in [
                ('particle_net_transport_J','粒子净输运（含侧向）','Net particle transport including sides','#2679ae','-'),
                ('field_ports_net_J','场出口减入口','Field exit minus entrance','#d37b28','--'),
                ('field_transverse_net_J','场四侧净能流','Net field energy through four sides','#57916f','-.'),
                ('field_time_caps_net_J','场晚端面减早端面','Late minus early field time cap','#a45d88',':'),
                ('physical_lab_sum_J','全部物理分项加和','Sum of all physical terms','#666666','-')]:
                curves.append((cn if zh else en,lab['z_lab_m'],lab[key]*1e6,color,style))
            plot(language,tag+'_fixed_lab_balance',curves,
                 labels[tag]+('：固定实验室管能量收支' if zh else ': fixed laboratory energy balance'),labx,micro)
    for path,expected in sources.items():assert sha(Path(path))==expected,path
    save(root/'RESULT.json',dict(complete=True,methods=summaries,no_new_physics=True,
        static_external_impulse_algebraically_canceled_not_independently_measured=True,
        actual_particle_force_substep_trace_available=False,
        fixed_lab_and_boost_slice_residuals_are_different_geometries=True,
        original_cached_momentum_statistics_reused=True))
    save(root/'PLOTS.json',dict(single_axes_only=True,figures=records))
    save(root/'SOURCE_DATA_SHA256.json',sources)
    print(json.dumps(summaries,indent=2),flush=True)


if __name__=='__main__':main(*sys.argv[1:])
