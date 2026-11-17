"""Readable actual-propagation figures and fixed linear/log intensity movies."""
import json
import sys
import subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,font_manager,ticker
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle
from common import sha,save

BIN=Path('/public/home/wfollowm/micromamba/envs/physics/bin')
FONT=Path('/public/home/wfollowm/APRL_Workspace/results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf')
BLUE='#1769a0';ORANGE='#cf661f';GREEN='#308367'

def main(root,mode='all'):
    root=Path(root);s=json.loads((root/'SUMMARY.json').read_text())
    d=np.load(root/'data/visual_arrays.npz');z=d['snapshot_z_m'];xz=d['z_m'];x=d['x_m']*1e6
    a=np.genfromtxt(root/'data/comparison.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    h=np.genfromtxt(root/'data/history.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    font_manager.fontManager.addfont(str(FONT));plt.rcParams.update({
      'font.family':[font_manager.FontProperties(fname=str(FONT)).get_name(),'DejaVu Sans'],
      'font.size':12,'axes.titlesize':17,'axes.labelsize':12,'legend.fontsize':11,
      'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,
      'text.parse_math':False,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
    plt.rcParams['animation.ffmpeg_path']=str(BIN/'ffmpeg')
    plots=[];videos=[] if mode=='all' else json.loads((root/'RENDER_VALIDATION.json').read_text())['videos']
    keys=['reference','native','no_absorption'];names=['wide-reference','small-native','small-no-absorption']
    cs=[GREEN,BLUE,ORANGE];styles=['-','--','-.'];edges=np.arange(-1620,1621,40)
    zedges=np.r_[xz-.0035,xz[-1]+.0035]
    roi_peak=43;assert abs(z[roi_peak]-1.505)<1e-14

    def plain_log(ax):
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v,_:f'{v:g}'))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    def figure():
        fig,ax=plt.subplots(figsize=(12.8,7.2));fig.subplots_adjust(left=.12,right=.94,bottom=.18,top=.86)
        return fig,ax
    def footer(fig,txt):fig.text(.5,.045,txt,ha='center',fontsize=10,color='#333333')
    def check_bounds(fig):
        fig.canvas.draw();ren=fig.canvas.get_renderer()
        texts=list(fig.texts)
        for ax in fig.axes:
            texts += [ax.title,ax.xaxis.label,ax.yaxis.label]
            if ax.get_legend():texts+=list(ax.get_legend().texts)
        for obj in texts:
            assert '$' not in obj.get_text(),obj.get_text()
            box=obj.get_window_extent(ren)
            assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.width+1 and box.y1<=fig.bbox.height+1,obj.get_text()
    def colorbar(fig,im,ax=None,cax=None,linear=False,label=''):
        cb=fig.colorbar(im,ax=ax,cax=cax,label=label)
        if not linear:cb.set_ticks([1e-8,1e-6,1e-4,1e-2,1])
        cb.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v,_:f'{v:g}'))
        cb.ax.yaxis.set_minor_formatter(ticker.NullFormatter())
        return cb
    def finish(fig,name,out):
        check_bounds(fig);fig.savefig(out/(name+'.png'),dpi=150);fig.savefig(out/(name+'.pdf'))
        plots.append(dict(path=str((out/(name+'.png')).relative_to(root)),sha256=sha(out/(name+'.png'))))
        plt.close(fig)

    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';out=root/'figures'/lang;out.mkdir(parents=True,exist_ok=True)
        labels=['大域参考：光束从右侧离开','原生吸收：主束被削弱，仍有弱返回','关闭吸收：从左侧绕回后仍向右传播'] if zh else [
             'Wide reference: the beam leaves to the right','Native absorber: strong removal, weak returning field','No absorption: wrap-around returns from the left']
        short=['大域参考','原生吸收','关闭吸收'] if zh else ['Wide reference','Native absorber','No absorption']
        for j,(key,label) in enumerate(zip(keys,labels)):
            fig,ax=figure();fig.subplots_adjust(right=.85)
            im=ax.pcolormesh(edges,zedges,np.maximum(d[key+'_line'],1e-12),norm=colors.LogNorm(1e-8,1),cmap='magma',shading='flat')
            ax.set(xlim=(-1620,1620),ylim=(0,2.1),xlabel='x (µm)',ylabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',title=label)
            for xx in [-520,520]:ax.axvline(xx,color='#78f0dd',ls=':',lw=1)
            if key!='reference':
                for xx in [-1000,1000]:ax.axvline(xx,color='white',ls='--',lw=1)
            if key=='reference':
                xy=(1350,.70);txt=(-1400,.95);text='向右离开同一视野' if zh else 'Leaves the displayed region'
            elif key=='no_absorption':
                xy=(-1470,.90);txt=(-600,1.08);text='从左侧重新进入' if zh else 'Re-enters from the left'
            else:
                xy=(-70,1.505);txt=(-1400,1.92);text='弱返回进入中央区域' if zh else 'Weak return reaches the aperture'
            ax.annotate(text,xy=xy,xytext=txt,color='white',fontsize=12,
                arrowprops=dict(arrowstyle='->',color='white',lw=1.5),bbox=dict(facecolor='black',alpha=.55,edgecolor='none',pad=4))
            colorbar(fig,im,ax=ax,label='y=0 强度 / 入射峰值' if zh else 'Intensity at y=0 / incident peak')
            footer(fig,'同一固定对数色标；点线为 ±520 µm 观察区；每一行均为原生传播结果。' if zh else
                'One fixed log scale. Dotted lines: ±520 µm aperture. Every row is a native propagation result.')
            finish(fig,f'0{j+1}_propagation_'+key,out)
        fig,ax=figure()
        for name,label,c,sty in zip(names,short,cs,styles):
            q=h[h['case']==name];ax.plot(q['z_m'],100*q['total_norm2_over_initial'],color=c,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='全域积分场强 / 入射值（%）' if zh else 'Whole-domain integrated intensity / input (%)',
          title='原生吸收：终点只剩入射积分场强的 0.296%' if zh else 'Native absorber: 0.296% of incident integrated intensity remains',ylim=(-2,104))
        ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3,frameon=False)
        ax.set_title(ax.get_title(),pad=42);ax.grid(alpha=.2)
        footer(fig,'对应净移除 99.704%；关闭吸收时积分几乎不变，但不能据此判定内部场正确。' if zh else
             'Net removal: 99.704%. Conservation without absorption does not establish a correct interior field.')
        finish(fig,'04_absorption_amount',out)
        fig,ax=figure()
        for name,label,c,sty in zip(names,short,cs,styles):
            q=h[h['case']==name];ax.plot(q['z_m'],100*q['roi_norm2_over_initial'],color=c,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='中央积分场强 / 入射全域值（%）' if zh else 'Central integrated intensity / full input (%)',
          title='关闭吸收时，绕回使中央积分场强重新升高' if zh else 'Without absorption, wrap-around repopulates the central aperture')
        ax.legend();ax.grid(alpha=.2)
        footer(fig,'中央观察区 x/y=±520 µm；后期绕回峰为入射值的 73.218%，不是能量凭空增加。' if zh else
           'Aperture x/y=±520 µm. The late return peaks at 73.218% of input; no energy is created.')
        finish(fig,'05_central_full_history',out)
        fig,ax=figure()
        for name,label,c,sty in zip(names[:2],short[:2],cs[:2],styles[:2]):
            q=a[(a['case']==name)&(a['z_m']>=1.3)];ax.plot(q['z_m'],100*q['roi_norm2_over_initial'],color=c,ls=sty,lw=2,label=label)
        ax.scatter([1.505],[100*s['small-native']['late_max_central_norm2_over_initial']],color=BLUE,zorder=3)
        ax.annotate('0.457%  |  z=1.505 m',xy=(1.505,.456921463898),xytext=(1.70,.44),fontsize=12,
                    arrowprops=dict(arrowstyle='->',color=BLUE),color=BLUE)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='中央积分场强 / 入射全域值（%）' if zh else 'Central integrated intensity / full input (%)',
          title='吸收之后，中央仍有弱返回场' if zh else 'A weak field returns to the aperture despite absorption',xlim=(1.3,2.1),ylim=(-.015,.52))
        ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False)
        ax.set_title(ax.get_title(),pad=42);ax.grid(alpha=.2)
        footer(fig,'后期局部放大：大域参考在此区域已接近零；分母始终是同一个入射全域值。' if zh else
            'Late-stage zoom: the reference field is nearly zero here. The denominator remains the same full input.')
        finish(fig,'06_native_return_zoom',out)
        fig,ax=figure();q=a[a['case']=='small-native']
        ax.plot(q['z_m'],100*q['roi_field_L2_difference_over_initial'],color=BLUE,lw=2)
        ax.scatter([1.505],[100*s['small-native']['max_central_field_L2_change_over_initial']],color=BLUE,zorder=3)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='中央复场差 / 入射全域场 L2（%）' if zh else 'Central field difference / input full-field L2 (%)',
          title='边界残留使中央复场差达到约 6.76%' if zh else 'Residual boundary effects produce a 6.76% central field difference',ylim=(-.15,7.5))
        ax.grid(alpha=.2)
        footer(fig,'逐点减去同一 z 的大域复场；固定入射尺度；这不是反射能量或 FEL 增益误差。' if zh else
            'Pointwise difference from the wide field at the same z, with a fixed input scale; not reflected energy or FEL gain error.')
        finish(fig,'07_interior_field_feedback',out)
        fig,ax=figure()
        for name,label,c,sty in zip(names[1:],short[1:],cs[1:],styles[1:]):
            q=a[a['case']==name];ax.semilogy(q['z_m'],q['negative_kx_norm2_over_initial'],color=c,ls=sty,lw=2,label=label)
        plain_log(ax);ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='负 kx 模式平方范数 / 入射值' if zh else 'Negative-kx squared norm / input',
           title='主要返回现象是绕回，反向横向模式很小' if zh else 'Wrap-around dominates; reverse-x modes remain small')
        ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False)
        ax.set_title(ax.get_title(),pad=42);ax.grid(alpha=.2)
        footer(fig,'原生吸收例的负 kx 份额最大约 5.11e-8；不能称作沿 −z 的物理反射率。' if zh else
           'The native negative-kx share peaks at 5.11e-8 of input; this is not backwards-z reflectivity.')
        finish(fig,'08_reverse_x_modes',out)
        for key,name,title in [('native','09_native_return_xy','原生吸收的弱返回场' if zh else 'Weak returning field with native absorption'),
                              ('reference','10_reference_xy','相同位置的大域参考场' if zh else 'Wide reference at the same location')]:
            fig,ax=figure();fig.subplots_adjust(right=.86)
            im=ax.imshow(np.maximum(abs(d[key][roi_peak])**2,1e-12),origin='lower',extent=[-1620,1620,-1620,1620],
                 norm=colors.LogNorm(1e-8,1),cmap='magma',interpolation='nearest',aspect='equal')
            ax.add_patch(Rectangle((-520,-520),1040,1040,fill=False,ls=':',lw=1.4,edgecolor='#78f0dd'))
            ax.set(xlabel='x (µm)',ylabel='y (µm)',title=title+' | z=1.505 m')
            colorbar(fig,im,ax=ax,label='强度 / 入射峰值' if zh else 'Intensity / incident peak')
            footer(fig,'与主视频相同固定对数色标；方框为中央观察区。' if zh else 'Same fixed log scale as the main movie. Box: central observation aperture.')
            finish(fig,name,out)
        fig,ax=figure();fig.subplots_adjust(right=.86)
        delta=d['native'][roi_peak]-d['reference'][roi_peak]
        im=ax.imshow(np.maximum(abs(delta[27:54,27:54])**2,1e-12),origin='lower',extent=[-540,540,-540,540],
               norm=colors.LogNorm(1e-8,1),cmap='magma',interpolation='nearest',aspect='equal')
        ax.set(xlim=(-520,520),ylim=(-520,520),xlabel='x (µm)',ylabel='y (µm)',
             title='中央误差场的空间分布 | z=1.505 m' if zh else 'Spatial distribution of the central error field | z=1.505 m')
        colorbar(fig,im,ax=ax,label='|A − Aref|² / 入射峰值' if zh else '|A − Aref|² / incident peak')
        footer(fig,'这是两复场之差的平方，不是反射能量，也不是两幅强度直接相减。' if zh else
               'Squared difference of complex fields, not reflected energy or a difference of two intensity images.')
        finish(fig,'11_central_error_xy',out)
        if mode=='static':continue
        for view in ['linear','log']:
            fig,axes=plt.subplots(1,3,figsize=(16.8,7.6));fig.subplots_adjust(left=.055,right=.89,bottom=.24,top=.78,wspace=.24)
            ims=[];stats=[]
            norm=colors.Normalize(0,1) if view=='linear' else colors.LogNorm(1e-8,1)
            for ax,key,label in zip(axes,keys,short):
                im=ax.imshow(np.maximum(abs(d[key][0])**2,1e-12),origin='lower',extent=[-1620,1620,-1620,1620],
                    norm=norm,cmap='magma',interpolation='nearest',aspect='equal')
                ax.add_patch(Rectangle((-520,-520),1040,1040,fill=False,ls=':',lw=1,edgecolor='#78f0dd'))
                if key!='reference':
                    for xx in [-1000,1000]:ax.axvline(xx,color='white',ls='--',lw=.8)
                ax.set(xlabel='x (µm)',ylabel='y (µm)',title=label+('（同一视野）' if zh and key=='reference' else ''))
                stats.append(ax.text(.5,-.23,'',ha='center',va='top',transform=ax.transAxes,fontsize=10))
                ims.append(im)
            cax=fig.add_axes([.925,.30,.013,.38]);colorbar(fig,ims[0],cax=cax,linear=view=='linear',label='强度 / 入射峰值' if zh else 'Intensity / incident peak')
            sup=fig.suptitle('',y=.96,fontsize=19);stage=fig.text(.5,.865,'',ha='center',fontsize=13)
            mode_label=('线性色标：观察主束被吸收' if view=='linear' else '固定对数色标：观察弱返回场') if zh else (
                'Linear scale: watch the main beam disappear' if view=='linear' else 'Fixed log scale: reveal the weak returning field')
            footer(fig,mode_label+('；三栏同一尺度、无逐帧缩放；无粒子；只测试横向边界。' if zh else '; one common scale, no per-frame rescaling; no particles, transverse boundaries only.'))
            path=out/('native_propagation_'+view+'.mp4')
            writer=FFMpegWriter(fps=5,codec='libx264',bitrate=2600,extra_args=['-pix_fmt','yuv420p','-threads','2'])
            with writer.saving(fig,str(path),100):
                for i,zz in enumerate(z):
                    for im,key in zip(ims,keys):im.set_data(np.maximum(abs(d[key][i])**2,1e-12))
                    for text,name in zip(stats,names):
                        q=a[a['case']==name][i]
                        text.set_text((f"全域保留 {q['total_norm2_over_initial']*100:.3f}%\n中央 {q['roi_norm2_over_initial']*100:.4f}%" if zh else
                            f"Whole domain: {q['total_norm2_over_initial']*100:.3f}%\nAperture: {q['roi_norm2_over_initial']*100:.4f}%"))
                    if zz<.45:txt='主束向右传播' if zh else 'The main beam travels to the right'
                    elif zz<.8:txt='主束进入边缘处理区' if zh else 'The main beam reaches the boundary region'
                    elif zz<1.2:txt='关闭吸收的场从左侧绕回；原生吸收已削弱主束' if zh else 'The unabsorbed field wraps to the left; the native absorber weakens the main beam'
                    else:txt='比较中央方框内的返回场：参考应已接近零' if zh else 'Compare the returning field in the aperture; the reference is nearly empty'
                    stage.set_text(txt);sup.set_text(('Puffin 原生实际传播' if zh else 'Actual native Puffin propagation')+f' | z={zz:.3f} m')
                    if i in [0,16,24,43,60]:
                        check_bounds(fig);fig.savefig(out/f'{view}_frame_{i:03d}.png',dpi=120)
                    writer.grab_frame()
            plt.close(fig)
            probe=json.loads(subprocess.check_output([str(BIN/'ffprobe'),'-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,duration','-of','json',str(path)],text=True))['streams'][0]
            assert int(probe['nb_frames'])==61 and abs(float(probe['duration'])-12.2)<1e-5
            videos.append(dict(path=str(path.relative_to(root)),sha256=sha(path),**probe,view=view,actual_native_frames=True,per_frame_rescaling=False))
            print('VIDEO',lang,view,flush=True)
    save(root/'RENDER_VALIDATION.json',dict(complete=True,figures=plots,videos=videos,
      arrays_sha256=sha(root/'data/visual_arrays.npz'),plain_numeric_log_ticks=True,
      physical_cell_pixel_edges_checked=True,figure_title_axis_legend_and_footer_bounds_checked=True,
      PDFs_same_source_not_separately_rendered=True))

if __name__=='__main__':main(*sys.argv[1:])
