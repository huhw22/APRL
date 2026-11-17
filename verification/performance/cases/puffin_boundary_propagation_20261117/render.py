"""Show actual native propagation, absorption, wrap-around and interior changes."""
import json
import subprocess
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors
from matplotlib.animation import FFMpegWriter
import numpy as np
from prepare import sha,save

BIN=Path('/public/home/wfollowm/micromamba/envs/physics/bin')
FONT=Path('/public/home/wfollowm/APRL_Workspace/results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf')

def main(root):
    root=Path(root);out=root/'analysis';r=json.loads((out/'RESULT.json').read_text());d=np.load(out/'visual_arrays.npz')
    font_manager.fontManager.addfont(str(FONT));plt.rcParams.update({'font.family':[font_manager.FontProperties(fname=str(FONT)).get_name(),'DejaVu Sans'],
        'font.size':12,'axes.titlesize':15,'axes.labelsize':12,'legend.fontsize':10,'axes.unicode_minus':False,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'text.parse_math':False})
    plt.rcParams['animation.ffmpeg_path']=str(BIN/'ffmpeg')
    a=np.genfromtxt(out/'comparison.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    h=np.genfromtxt(out/'history.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    keys=['reference','native','no_absorption'];z=d['snapshot_z_m'];xx=d['x_m']*1e6
    files=[];videos=[]
    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';target=out/'figures'/lang;target.mkdir(parents=True,exist_ok=True)
        labels=['大计算域参考','原生吸收边界','关闭吸收：周期绕回控制'] if zh else ['Wide-domain reference','Native absorbing boundary','No absorption: periodic wrap-around']
        def savefig(fig,name):
            fig.savefig(target/(name+'.png'),dpi=150);fig.savefig(target/(name+'.pdf'))
            files.append(str((target/(name+'.png')).relative_to(root)));plt.close(fig)
        for key,label in zip(keys,labels):
            fig,ax=plt.subplots(figsize=(10.8,7.2));fig.subplots_adjust(left=.13,right=.86,bottom=.16,top=.86)
            image=ax.imshow(np.maximum(d[key+'_line'],1e-10),origin='lower',extent=[-1600,1600,0,d['z_m'][-1]],
                      aspect='auto',norm=colors.LogNorm(1e-8,1),cmap='magma',interpolation='nearest')
            if key!='reference':
                ax.axvline(-1000,color='white',ls='--',lw=1);ax.axvline(1000,color='white',ls='--',lw=1)
            ax.axvline(-520,color='#7de2de',ls=':',lw=1);ax.axvline(520,color='#7de2de',ls=':',lw=1)
            ax.set(xlabel='x (µm)',ylabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',title=label)
            fig.colorbar(image,ax=ax,label='y=0 场强 / 入射峰值（固定对数色标）' if zh else 'Intensity at y=0 / incident peak (fixed log scale)')
            fig.text(.5,.035,'小域虚线：边缘处理区内缘；点线：中央观测区。每行来自实际原生推进。' if zh else
                         'Small-domain dashed lines: boundary-layer inner edges. Dotted: observation aperture. Native propagation.',ha='center',fontsize=10)
            savefig(fig,'spacetime_'+key)
        case_names=['wide-reference','small-native','small-no-absorption'];cols=['#308367','#1769a0','#cf661f'];styles=['-', '--', '-.']
        fig,ax=plt.subplots(figsize=(11,7));fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.86)
        for name,label,col,sty in zip(case_names,labels,cols,styles):
            q=h[h['case']==name];ax.plot(q['z_m'],q['total_norm2_over_initial'],color=col,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='全域积分场强 / 相同入射值' if zh else 'Whole-domain integrated intensity / common input',
               title='能量代理实际减少多少' if zh else 'How much integrated field intensity is removed');ax.grid(alpha=.2);ax.legend()
        fig.text(.5,.035,'同一入射，无逐帧缩放；这是前向模型的场强积分，不是完整 Maxwell 能量账本。' if zh else
                  'Same incident field, no frame-by-frame rescaling. A forward-model intensity integral, not a full Maxwell ledger.',ha='center',fontsize=10)
        savefig(fig,'total_intensity')
        fig,ax=plt.subplots(figsize=(11,7));fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.86)
        for name,label,col,sty in zip(case_names,labels,cols,styles):
            q=h[h['case']==name];ax.plot(q['z_m'],q['roi_norm2_over_initial'],color=col,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='中央积分场强 / 相同入射值' if zh else 'Central integrated intensity / common input',
               title='离开后的场是否重新进入中央区域' if zh else 'Does the field re-enter the central aperture?');ax.grid(alpha=.2);ax.legend()
        fig.text(.5,.035,'共同观测区：x/y = ±520 µm；横向绕回与传播方向反转分开统计。' if zh else
                   'Common aperture: x/y = ±520 µm. Wrap-around and reversal of transverse direction are distinct.',ha='center',fontsize=10)
        savefig(fig,'central_intensity')
        fig,ax=plt.subplots(figsize=(11,7));fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.86)
        for name,label,col,sty in zip(case_names[1:],labels[1:],cols[1:],styles[1:]):
            q=a[a['case']==name];ax.plot(q['z_m'],100*q['roi_field_L2_difference_over_initial'],color=col,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='中央复场差 / 入射全域 L2（%）' if zh else 'Central field difference / input full-domain L2 (%)',
               title='边界对内部计算域的影响' if zh else 'Boundary effects inside the observation aperture');ax.grid(alpha=.2);ax.legend()
        fig.text(.5,.035,'逐点减去同一时刻的大域复场；这不是反射能量，未包含粒子反馈。' if zh else
                  'Pointwise difference from the wide-domain field at the same z. Not reflected energy; no electron feedback.',ha='center',fontsize=10)
        savefig(fig,'interior_field_difference')
        fig,ax=plt.subplots(figsize=(11,7));fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.86)
        for name,label,col,sty in zip(case_names[1:],labels[1:],cols[1:],styles[1:]):
            q=a[a['case']==name];ax.semilogy(q['z_m'],np.maximum(q['negative_kx_norm2_over_initial'],1e-16),color=col,ls=sty,lw=2,label=label)
        ax.set(xlabel='传播距离 z（m）' if zh else 'Propagation distance z (m)',ylabel='负 kx 模式范数 / 入射范数' if zh else 'Negative-kx mode norm / incident norm',
               title='检查是否生成横向反向模式' if zh else 'Check for newly generated reverse-x modes');ax.grid(alpha=.2);ax.legend()
        fig.text(.5,.035,'完整小域傅里叶基底中的方向诊断；不是沿 −z 的物理反射率。' if zh else
                  'Directional diagnostic in the full small-domain Fourier basis; not a backwards-z reflectivity.',ha='center',fontsize=10)
        savefig(fig,'negative_kx_modes')
        # Three-panel movie: every saved field is from native along-z propagation.
        fig,axes=plt.subplots(1,3,figsize=(16.8,7.2));fig.subplots_adjust(left=.055,right=.9,bottom=.20,top=.82,wspace=.24)
        ims=[]
        for ax,key,label in zip(axes,keys,labels):
            im=ax.imshow(np.maximum(abs(d[key][0])**2,1e-10),origin='lower',extent=[-1600,1600,-1600,1600],
                norm=colors.LogNorm(1e-8,1),cmap='magma',interpolation='nearest',aspect='equal')
            ax.set(title=label,xlabel='x (µm)',ylabel='y (µm)')
            edges=[(520,'#7de2de',':')]
            if key!='reference':edges.append((1000,'white','--'))
            for e,col,sty in edges:
                ax.axvline(-e,color=col,ls=sty,lw=.7);ax.axvline(e,color=col,ls=sty,lw=.7)
            ims.append(im)
        cax=fig.add_axes([.925,.27,.012,.43]);fig.colorbar(ims[0],cax=cax,label='|A|² / 入射峰值' if zh else '|A|² / incident peak')
        sup=fig.suptitle('',y=.95,fontsize=19)
        fig.text(.5,.08,'同一输入的实际 Puffin 原生传播；固定色标；左图仅显示大域中的同一块区域。' if zh else
                'Actual native Puffin propagation of the same input; fixed scale; reference shows the same subregion.',ha='center',fontsize=11)
        fig.text(.5,.035,'逐帧增加的是传播距离 z，不是人为移动输入位置；本控制只测试横向边界，不含粒子。' if zh else
                'Frames advance in propagation distance z, not prescribed input position. Transverse boundaries only; no particles.',ha='center',fontsize=10)
        movie=target/'native_propagation.mp4'
        writer=FFMpegWriter(fps=6,codec='libx264',bitrate=2600,extra_args=['-pix_fmt','yuv420p','-threads','2'])
        with writer.saving(fig,str(movie),100):
            for i,zz in enumerate(z):
                for im,key in zip(ims,keys):im.set_data(np.maximum(abs(d[key][i])**2,1e-10))
                sup.set_text(('传播、吸收与计算域内的返回场' if zh else 'Propagation, absorption and returning field')+f' | z = {zz:.3f} m')
                if i in [0,16,24,40,60]:fig.savefig(target/f'frame_{i:03d}.png',dpi=120)
                writer.grab_frame()
        plt.close(fig)
        stream=json.loads(subprocess.check_output([str(BIN/'ffprobe'),'-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,duration','-of','json',str(movie)],text=True))['streams'][0]
        assert int(stream['nb_frames'])==len(z)
        videos.append(dict(path=str(movie.relative_to(root)),sha256=sha(movie),**stream))
        print('RENDERED',lang,flush=True)
    save(out/'RENDER_VALIDATION.json',dict(complete=True,figures=files,videos=videos,actual_native_propagation=True,
          same_input_common_color_scale=True,awaiting_visual_review=True))

if __name__=='__main__':main(sys.argv[1])
