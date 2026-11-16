"""Bilingual standalone figures; no physics evolution or fitted normalization."""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,ticker,colors

def read(p):return json.loads(Path(p).read_text())
def data(p):return np.genfromtxt(p,delimiter=',',names=True,dtype=np.longdouble)

def main(root,font,only_boost=False):
    only_sums=only_boost=='sums-only';only_boost=only_boost=='boost-only' or only_boost is True
    root=Path(root).resolve();result=read(root/'RESULT.json')
    font_manager.fontManager.addfont(font);family=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.family':[family,'DejaVu Sans'],'font.size':12,'axes.titlesize':16,
        'axes.labelsize':13,'axes.spines.top':False,'axes.spines.right':False,
        'text.parse_math':False,'axes.unicode_minus':False,'pdf.fonttype':42,
        'figure.facecolor':'white','savefig.facecolor':'white'})
    plots=[]
    def plot(lang,name,series,title,xlabel,ylabel):
        out=root/'figures'/lang;out.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(10.8,6.6),layout='constrained');records=[]
        for i,(label,x,y) in enumerate(series):
            x,y=np.asarray(x,float),np.asarray(y,float);good=np.isfinite(x)&np.isfinite(y)
            color='#111111' if '加和' in label or label=='Sum' else None
            ax.plot(x[good],y[good],label=label,lw=1.8,linestyle=['-','--','-.',':'][i%4],
                    marker='o' if len(x)<20 else None,markersize=4,color=color)
            records.append(dict(label=label,x=x[good].tolist(),y=y[good].tolist()))
        ax.axhline(0,color='#888888',lw=.7);ax.grid(alpha=.18)
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:g}'))
        ax.legend(frameon=False,ncol=2 if len(series)>4 else 1,fontsize=10)
        fig.canvas.draw()
        for obj in [ax.title,ax.xaxis.label,ax.yaxis.label]+ax.get_legend().texts:
            assert '$' not in obj.get_text() and '_' not in obj.get_text()
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=170)
        plots.append(dict(name=name,language=lang,data_axes=1,series=records));plt.close(fig)

    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';time='计算系模拟时间（ps）' if zh else 'Computational time (ps)'
        position='实验室系探测面位置（m）' if zh else 'Laboratory detector position (m)'
        unit='带符号能量（fJ）' if zh else 'Signed energy (fJ)'
        for case in result['cases']:
            if only_sums:continue
            tag=case['tag'];o=root/'cases'/tag
            label=('三点核' if case['method']=='F3' else '动量匹配（五次补偿）') if zh else case['method']
            kind=('弯铁' if case['kind']=='bend' else '直线') if zh else case['kind']
            suffix=('样本一' if case['seed']==20261111 else '样本二') if zh else str(case['seed'])
            p=data(o/'boost/physical_components.csv')
            keys=['particle_change_J','field_change_J','front_outward_minus_injection_J','back_outward_plus_cropped_J',
                  'side_outward_J','moving_interface_J','external_and_constraint_work_signed_J','physical_component_sum_J']
            labels=['粒子能量变化','内场能量变化','前端向外减注入','后端向外加裁出','四侧净传出','移窗接口项','外场及约束源功（负号）','物理分项加和'] if zh else [
                'Particle change','Field change','Front outward minus injection','Back outward plus cropped',
                'Transverse outward','Moving interface','Negative external and constraint work','Sum']
            plot(lang,tag+'_boost_balance',[(lb,p['time_box_s']*1e12,p[k]*1e15) for lb,k in zip(labels,keys)],
                 f'{kind} {label} {suffix}'+('：计算系分项收支' if zh else ': computational components'),time,unit)
            if only_boost:continue
            q=data(o/'lab/complete_components.csv')
            keys=['particle_net_transport_J','field_ports_net_J','field_transverse_net_J','field_time_caps_net_J','physical_lab_sum_J']
            labels=['粒子净输运','场出口减入口','场四侧净流','晚减早时间端面','完整加和'] if zh else [
                'Net particle transport','Field exit minus entrance','Transverse field flux','Late minus early time cap','Sum']
            plot(lang,tag+'_lab_balance',[(lb,q['z_lab_m'],q[k]*1e15) for lb,k in zip(labels,keys)],
                 f'{kind} {label} {suffix}'+('：实验室系分项收支' if zh else ': laboratory components'),position,unit)
        if only_boost:continue
        for seed in [20261111,20261112]:
            for kind in ['bend','straight']:
                name=f'{seed}_{kind}'
                for key,filename,xkey,scale,xlab,title,ylab in [
                    ('physical_component_sum_J','boost/physical_components.csv','time_box_s',1e15,time,
                     '计算系物理分项加和' if zh else 'Computational physical component sum',unit),
                    ('total_residual_J','boost/complete_components.csv','time_box_s',1e24,time,
                     '计算系离散加和' if zh else 'Computational discrete sum','带符号能量（1e-24 J）' if zh else 'Signed energy (1e-24 J)'),
                    ('physical_lab_sum_J','lab/complete_components.csv','z_lab_m',1e15,position,
                     '实验室系完整加和' if zh else 'Complete laboratory sum',unit)]:
                    series=[]
                    for method,label in [('F3','三点核' if zh else 'F3'),('MM','五次补偿' if zh else 'Compensated quintic')]:
                        a=data(root/'cases'/f'{method}_{name}'/filename)
                        series.append((label,a[xkey]*(1e12 if xkey=='time_box_s' else 1),a[key]*scale))
                    plot(lang,name+'_'+key.replace('_J',''),series,title+f' {seed} {kind}',xlab,ylab)
            if only_sums:continue
            spec=[]
            for method,label in [('F3','三点核' if zh else 'F3'),('MM','五次补偿' if zh else 'Compensated quintic')]:
                a=data(root/'radiation'/f'{method}_{seed}_axis_spectrum.csv')
                spec.append((label,a['photon_eV'],a['PIC_paired_J_per_eV_m2']*1e6))
            spec.append(('规定轨迹 LW' if zh else 'Prescribed-orbit LW',a['photon_eV'],a['prescribed_orbit_LW_J_per_eV_m2']*1e6))
            plot(lang,f'{seed}_axis_spectrum',spec,
                '轴上弯转减直线频谱' if zh else 'On-axis bend-minus-straight spectrum',
                '光子能量（eV）' if zh else 'Photon energy (eV)',
                '谱通量密度（µJ / eV / m2）' if zh else 'Spectral flux density (µJ / eV / m2)')
            for mapkey,title in [('LW_full','规定轨迹 LW 总场' if zh else 'Prescribed-orbit full LW'),
                                 ('LW_acc','规定轨迹 LW 加速度项' if zh else 'Prescribed-orbit acceleration LW'),
                                 ('F3','自洽三点核总场' if zh else 'Self-consistent F3 total field'),
                                 ('MM','自洽五次补偿总场' if zh else 'Self-consistent quintic total field')]:
                with np.load(root/'radiation'/f'F3_{seed}_maps.npz') as z:
                    x,y=z['x'],z['y'];reference=z['LW_full'];image=z[mapkey] if mapkey.startswith('LW') else None
                with np.load(root/'radiation'/f'MM_{seed}_maps.npz') as z:mm=z['PIC']
                with np.load(root/'radiation'/f'F3_{seed}_maps.npz') as z:ff=z['PIC']
                if mapkey=='F3':image=ff
                if mapkey=='MM':image=mm
                vmax=max(reference[10:-10,10:-10].max(),ff[10:-10,10:-10].max(),mm[10:-10,10:-10].max())*1e6
                fig,ax=plt.subplots(figsize=(9.5,6.6))
                fig.subplots_adjust(left=.12,right=.83,bottom=.14,top=.89)
                im=ax.imshow(image[10:-10,10:-10]*1e6,origin='lower',aspect='equal',cmap='magma',
                    extent=[(x[10]-2.5e-6)*1e6,(x[-11]+2.5e-6)*1e6,(y[10]-2.5e-6)*1e6,(y[-11]+2.5e-6)*1e6],
                    norm=colors.PowerNorm(gamma=.45,vmin=0,vmax=vmax))
                ax.set(title=title+f' 8 eV {seed}',xlabel='原实验室横坐标（µm）' if zh else 'Original laboratory x (µm)',
                       ylabel='横向位置 y（µm）' if zh else 'Transverse y (µm)')
                color_axis=fig.add_axes([.86,.15,.025,.70])
                fig.colorbar(im,cax=color_axis,label='谱场强度（µJ / eV / m2）' if zh else 'Spectral field intensity (µJ / eV / m2)')
                for ext in ['png','pdf']:fig.savefig(root/'figures'/lang/f'{seed}_{mapkey}_8eV.png'.replace('.png','.'+ext),dpi=170)
                plots.append(dict(name=f'{seed}_{mapkey}_8eV',language=lang,data_axes=1,color_scale_max=vmax));plt.close(fig)
            line=[]
            for method,label in [('F3','三点核' if zh else 'F3'),('MM','五次补偿' if zh else 'Compensated quintic')]:
                with np.load(root/'radiation'/f'{method}_{seed}_maps.npz') as z:
                    iy=np.sort(np.argsort(abs(z['y']))[:2]);xx=z['x'];yy=z['PIC'][iy].mean(axis=0)
                    ref=z['LW_full'][iy].mean(axis=0)
                    line.append((label,xx[10:-10]*1e6,yy[10:-10]*1e6))
            line.append(('完整 LW 参照' if zh else 'Full LW reference',xx[10:-10]*1e6,ref[10:-10]*1e6))
            plot(lang,f'{seed}_8eV_horizontal_profile',line,'8 eV 横向强度剖面' if zh else '8 eV transverse intensity profile',
                '原实验室横坐标（µm）' if zh else 'Original laboratory x (µm)',
                '谱场强度（µJ / eV / m2）' if zh else 'Spectral field intensity (µJ / eV / m2)')
    if only_boost or only_sums:
        replaced={(p['name'],p['language']) for p in plots}
        plots+=[p for p in read(root/'PLOTS.json')['figures'] if (p['name'],p['language']) not in replaced]
    (root/'PLOTS.json').write_text(json.dumps(dict(standalone_only=True,figures=plots),ensure_ascii=False,indent=2)+'\n')
    print('FIGURES',len(plots),flush=True)

if __name__=='__main__':main(*sys.argv[1:])
