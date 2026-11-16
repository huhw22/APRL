"""Export figures and explicitly labelled prescribed-input operator video."""
import json
import subprocess
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FFMpegWriter
import numpy as np
from prepare import sha, save

FONT=Path('/public/home/wfollowm/APRL_Workspace/results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf')
BIN=Path('/public/home/wfollowm/micromamba/envs/physics/bin')
BLUE='#1769a0';ORANGE='#cf661f';GREEN='#308367'

def main(root,mode='all'):
    root=Path(root);r=json.loads((root/'RESULT.json').read_text())
    font_manager.fontManager.addfont(str(FONT))
    plt.rcParams.update({'font.family':[font_manager.FontProperties(fname=str(FONT)).get_name(),'DejaVu Sans'],
       'font.size':13,'axes.titlesize':17,'axes.labelsize':13,'legend.fontsize':11,
       'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.unicode_minus':False,
       'figure.facecolor':'white','savefig.facecolor':'white','text.parse_math':False})
    plt.rcParams['animation.ffmpeg_path']=str(BIN/'ffmpeg')
    n=np.load(root/'data/operator_illustration.npz');x=n['x_m']*1e6;y=n['y_m']*1e6
    tails=np.genfromtxt(root/'data/tail_mask_comparison.csv',names=True,delimiter=',')
    rows=np.genfromtxt(root/'data/prescribed_position_scan.csv',names=True,delimiter=',')
    records=[]
    videos=[] if mode=='all' else json.loads((root/'RENDER_VALIDATION.json').read_text())['videos']
    def figure():
        f,a=plt.subplots(figsize=(12.8,7.2));f.subplots_adjust(left=.12,right=.93,bottom=.18,top=.84)
        return f,a
    def footer(f,text): f.text(.5,.035,text,ha='center',fontsize=10,color='#444444')
    def finish(f,name,folder):
        f.canvas.draw();ren=f.canvas.get_renderer()
        for a in f.axes:
            objs=[a.title,a.xaxis.label,a.yaxis.label]
            if a.get_legend():objs+=list(a.get_legend().texts)
            for t in objs:
                b=t.get_window_extent(ren)
                assert b.x0>=-1 and b.y0>=-1 and b.x1<=f.bbox.width+1 and b.y1<=f.bbox.height+1,(name,t.get_text())
        f.savefig(folder/(name+'.png'),dpi=150);f.savefig(folder/(name+'.pdf'))
        records.append(dict(path=str((folder/(name+'.png')).relative_to(root)),sha256=sha(folder/(name+'.png'))))
        plt.close(f)
    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';out=root/'figures'/lang;out.mkdir(parents=True,exist_ok=True)
        f,a=figure();im=a.imshow(n['mask'],origin='lower',extent=[-1600,1600,-1600,1600],vmin=0,vmax=1,cmap='viridis',aspect='equal',interpolation='nearest')
        a.set(xlabel='x (µm)',ylabel='y (µm)',title='Puffin 横向边界的实际掩膜权重' if zh else 'Actual transverse mask from the Puffin library')
        a.axvline(-1000,color='white',ls='--',lw=.8);a.axvline(1000,color='white',ls='--',lw=.8)
        f.colorbar(im,ax=a,label='进入滤波分支的权重 M' if zh else 'Weight M sent to the filter branch')
        footer(f,'81 × 81 节点，40 µm 间隔；每侧 16 节点。M 不是反射率。' if zh else
            '81 × 81 nodes, 40 µm spacing, 16 boundary nodes per side. M is not reflectivity.')
        finish(f,'01_transverse_mask',out)
        f,a=figure();a.plot(x,n['mask'][40],color=BLUE,lw=2,label='y = 0')
        a.set(xlabel='x (µm)',ylabel='M',ylim=(-.05,1.1),title='吸收权重在边缘渐增，中央区域为零' if zh else 'The boundary weight rises near the transverse edges')
        a.axvspan(-520,520,color=GREEN,alpha=.1,label='既有比较孔径 ±520 µm' if zh else 'Existing comparison aperture ±520 µm')
        a.legend();a.grid(alpha=.2);footer(f,'这里显示的是场的分支权重，不是每步直接删除的能量份额。' if zh else
            'This is a field-splitting weight, not the fraction of energy deleted in one step.')
        finish(f,'02_mask_profile',out)
        f,a=figure();k=np.linspace(0,60000,301);g=r['grid'];scale=json.loads((root/'evidence/scale.json').read_text())
        kz=g['carrier_kz2_bar']
        d=np.exp(-g['h_bar']*k*scale['transverse_scale_m']/np.sqrt(2*kz))
        a.plot(k,d,color=BLUE,lw=2,label='边缘支路场系数 D' if zh else 'Boundary-branch amplitude factor D')
        a.plot(k,d*d,color=ORANGE,lw=2,ls='--',label='同一傅里叶分量的平方范数系数 D²' if zh else 'Squared-norm factor D² of the same Fourier component')
        a.set(xlabel='kx (1/m), ky = 0',ylabel='保留比例' if zh else 'Retained fraction',ylim=(.7,1.025),
          title='边界滤波依赖横向波数与频率' if zh else 'Boundary filtering depends on wave number and frequency')
        a.legend();a.grid(alpha=.2)
        footer(f,'β = 1；单次衍射步约 7 mm；仅为边缘分支系数，不是整束能量或反射率。' if zh else
            'β = 1, one nominal 7 mm step. These branch factors are neither whole-beam loss nor reflectivity.')
        finish(f,'03_modal_factors',out)
        s=tails[tails['nz']==96000]
        f,a=figure();a.plot(s['distance_from_tail_um'],s['mask_native'],color=ORANGE,lw=2,marker='o',ms=3,label='原生库实际值' if zh else 'Values from the existing native library')
        a.plot(s['distance_from_tail_um'],s['geometric_parentheses_reference'],color=BLUE,ls='--',lw=2,
            label='只改长度括号的几何参考' if zh else 'Geometric reference: length parentheses only')
        a.set(xlabel='距纵向末节点的距离（µm）' if zh else 'Distance from the final longitudinal node (µm)',
          ylabel='纵向掩膜权重' if zh else 'Longitudinal mask weight',ylim=(-.04,1.08),
          title='纵向掩膜的长度表达式使曲线发生错位' if zh else 'The longitudinal length expression shifts the mask profile')
        a.legend(loc='upper right');a.grid(alpha=.2)
        footer(f,'原式 nx − 1 × dx；几何参考 (nx − 1) × dx。参考曲线尚未作为修复验收。' if zh else
            'Original: nx − 1 × dx. Reference: (nx − 1) × dx. The reference is not an accepted fix.')
        finish(f,'04_longitudinal_length_issue',out)
        f,a=figure()
        for nz,col,sty in [(95980,GREEN,':'),(96000,ORANGE,'-'),(96020,BLUE,'--')]:
            s=tails[tails['nz']==nz];a.plot(s['distance_from_tail_um'],s['mask_native'],color=col,ls=sty,lw=2,label=f'Nz2 = {nz}')
        a.set(xlabel='距纵向末节点的距离（µm）' if zh else 'Distance from the final longitudinal node (µm)',ylabel='M',ylim=(-.04,1.08),
          title='仅增减一个波长的网格长度也会改变末端掩膜' if zh else 'Adding or removing one wavelength changes the end mask')
        a.legend(loc='lower center',bbox_to_anchor=(.5,1.025),ncol=3,frameon=False)
        a.set_title(a.get_title(),pad=48)
        a.grid(alpha=.2);footer(f,'同一 dx；每波长 20 间隔。这里只重新评价函数，没有重跑 FEL。' if zh else
            'Fixed dx, 20 intervals per wavelength. Function evaluations only; no FEL rerun.')
        finish(f,'05_longitudinal_grid_sensitivity',out)
        f,a=figure();a.plot(rows['x0_um'],rows['complex_field_L2_difference_percent'],color=BLUE,lw=2,label='复场 L2 差 / 输入 L2' if zh else 'Complex-field L2 change / input L2')
        a.plot(rows['x0_um'],rows['energy_proxy_change_percent'],color=ORANGE,ls='--',lw=2,label='平方范数净变化 / 输入平方范数' if zh else 'Net squared-norm change / input squared norm')
        a.axvline(1000,color='#777777',ls=':');a.axhline(0,color='#777777',lw=.7)
        a.set(xlabel='规定输入的中心位置 x0（µm）' if zh else 'Prescribed input center x0 (µm)',ylabel='%',
           title='场变化与能量变化必须分别计算' if zh else 'Field changes and energy changes are different measurements')
        a.legend(loc='lower left');a.grid(alpha=.2);footer(f,'独立规定场的一次边界操作；不是反射率、沿程增长或生产工况误差。' if zh else
            'One boundary operation on independent prescribed fields; not reflectivity or a production-run error.')
        finish(f,'06_field_vs_energy',out)
        f,a=figure();phi=np.linspace(0,np.pi,301)
        a.plot(phi*180/np.pi,2*np.sin(phi/2)*100,color=BLUE,lw=2,label='复场相对差' if zh else 'Relative complex-field difference')
        a.plot(phi*180/np.pi,np.zeros_like(phi),color=ORANGE,lw=2,ls='--',label='能量相对变化' if zh else 'Relative energy change')
        a.scatter([60],[100],color=BLUE,zorder=3)
        a.set(xlabel='整体相位差（度）' if zh else 'Uniform phase difference (degrees)',ylabel='%',
           title='相差 60°：复场差可达 100%，能量完全相同' if zh else 'A 60° phase shift gives 100% field difference and unchanged energy')
        a.legend();a.grid(alpha=.2);footer(f,'严格数学示例，用于解释指标；不是 Puffin 已测得 60° 相位误差。' if zh else
            'Exact mathematical example; it does not claim a measured 60° Puffin phase error.')
        finish(f,'07_phase_example',out)
        for name,array,title in [('08_prescribed_input',n['input_A'][128],('规定输入场强度' if zh else 'Prescribed input intensity')),
                                 ('09_single_boundary_output',n['output_A'][128],('一次边界操作后的场强度' if zh else 'Intensity after one boundary operation'))]:
            f,a=figure();im=a.imshow(abs(array)**2,origin='lower',extent=[-1600,1600,-1600,1600],vmin=0,vmax=1,cmap='magma',aspect='equal',interpolation='nearest')
            a.set(xlabel='x (µm)',ylabel='y (µm)',title=title+' | x0 = 1280 µm')
            f.colorbar(im,ax=a,label='|A|²（同一绝对色标，任意单位）' if zh else '|A|² (same fixed scale, arbitrary units)')
            footer(f,'按源码公式的单次响应；没有衍射传播、粒子源或纵向掩膜。' if zh else
                'Single response from the source formula; no diffraction march, particles, or longitudinal mask.')
            finish(f,name,out)
        if mode=='static':continue
        # Each frame is an independent prescribed input, never a time/propagation step.
        f,axes=plt.subplots(1,2,figsize=(12.8,7.2));f.subplots_adjust(left=.09,right=.96,bottom=.18,top=.77,wspace=.30)
        lines=[]
        for ax in axes:
            ax.set(xlabel='x (µm)',ylabel='Re A (a.u.)',xlim=(-1600,1600),ylim=(-1.1,1.1))
            ax.axvspan(-1600,-1000,color='#aaaaaa',alpha=.2);ax.axvspan(1000,1600,color='#aaaaaa',alpha=.2);ax.grid(alpha=.2)
        axes[0].set_title('规定输入复场的实部' if zh else 'Real part of the prescribed input')
        axes[1].set_title('一次边界操作后的实部' if zh else 'Real part after one boundary operation')
        l0,=axes[0].plot(x,n['input_A'][0,40].real,color=BLUE,lw=2)
        l1,=axes[1].plot(x,n['output_A'][0,40].real,color=ORANGE,lw=2)
        sup=f.suptitle('',y=.94,fontsize=17);caption=f.text(.5,.085,'',ha='center',fontsize=12)
        footer(f,'规定位置扫描；每帧独立计算；不是 Puffin 沿程仿真，也不是反射率测量。' if zh else
           'Prescribed-position scan: independent frames, not Puffin propagation or a reflectivity measurement.')
        def update(i):
            l0.set_ydata(n['input_A'][i,40].real);l1.set_ydata(n['output_A'][i,40].real)
            sup.set_text(('边缘场的一次滤波响应' if zh else 'One boundary-filter response')+f" | x0 = {rows['x0_um'][i]:.0f} µm")
            caption.set_text((f"复场 L2 差 = {rows['complex_field_L2_difference_percent'][i]:.3f}%    平方范数净变化 = {rows['energy_proxy_change_percent'][i]:+.3f}%" if zh else
                f"Complex-field L2 change = {rows['complex_field_L2_difference_percent'][i]:.3f}%    Net squared-norm change = {rows['energy_proxy_change_percent'][i]:+.3f}%"))
        path=out/'prescribed_boundary_response.mp4';writer=FFMpegWriter(fps=15,codec='libx264',bitrate=1800,extra_args=['-pix_fmt','yuv420p','-threads','2'])
        with writer.saving(f,str(path),100):
            for i in range(len(rows)):
                update(i)
                if i in [0,96,128,140]:f.savefig(out/f'video_frame_{i:03d}.png',dpi=120)
                writer.grab_frame()
        plt.close(f)
        stream=json.loads(subprocess.check_output([str(BIN/'ffprobe'),'-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,duration','-of','json',str(path)],text=True))['streams'][0]
        assert int(stream['nb_frames'])==141 and abs(float(stream['duration'])-9.4)<1e-5
        videos.append(dict(path=str(path.relative_to(root)),sha256=sha(path),**stream,independent_prescribed_frames=True,not_propagation=True))
        print('RENDERED',lang,flush=True)
    source_hashes={str(p.relative_to(root)):sha(p) for p in [root/'RESULT.json',root/'data/operator_illustration.npz',root/'data/tail_mask_comparison.csv',root/'data/prescribed_position_scan.csv']}
    save(root/'RENDER_VALIDATION.json',dict(complete=True,figures=records,videos=videos,sources=source_hashes,
         all_static_title_axis_legend_bounds_checked=True,pdfs_same_source_not_separately_rendered=True))

if __name__=='__main__':main(*sys.argv[1:])
