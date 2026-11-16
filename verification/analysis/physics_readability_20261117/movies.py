"""Re-export saved vacuum frames with plane/component titles and fixed scales."""
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle
import numpy as np
from render import configure, read_table, save, sha

C = 299792458.
BIN = Path('/public/home/wfollowm/micromamba/envs/physics/bin')
plt.rcParams['animation.ffmpeg_path'] = str(BIN/'ffmpeg')


def analytic_y(x,y,z,t,m):
    r=np.sqrt(x*x+y*y+z*z);rr=np.where(r>0,r,1.)
    sig=m['gaussian_sigma_m'];a=m['shell_radius_m']+C*t
    u=(r-a)/sig;v=(r+a)/sig
    gu=m['hertz_amplitude']*np.exp(-u*u/2);gv=m['hertz_amplitude']*np.exp(-v*v/2)
    f=gu-gv;fr=(-u*gu+v*gv)/sig;frr=((u*u-1)*gu-(v*v-1)*gv)/sig**2
    pr=fr/rr-f/rr**2;prr=frr/rr-2*fr/rr**2+2*f/rr**3
    ey=(prr-pr/rr)*(y/rr)**2-(prr+pr/rr)
    origin=4*m['hertz_amplitude']*np.exp(-a*a/(2*sig*sig))*(a/sig**4-a**3/(3*sig**6))
    return np.where(r>0,ey,origin)


def main(root,font,metadata,mode='all'):
    root=Path(root);out=root/'07_readable_figures/00_vacuum';out.mkdir(parents=True,exist_ok=True)
    configure(font)
    src=root/'01_reviewed_physics/01_vacuum_and_short_drift'
    metadata=Path(metadata);m=json.loads(metadata.read_text())
    inputs={str(p):sha(p) for p in [metadata,src/'vacuum_time_fields.npz',src/'vacuum_energy_and_error.csv',src/'cpml_time_fields.npz']}
    save(out/'VACUUM_METADATA.json',m)
    d=np.load(src/'vacuum_time_fields.npz');h=read_table(src/'vacuum_energy_and_error.csv')
    assert np.array_equal(h['step'],d['step'])
    x,z,y=d['x_m']*1e6,d['z_m']*1e6,d['y_m']*1e6;t=d['time_s']*1e12
    xz,xy=d['Ey_xz_V_per_m']/1e6,d['Ey_xy_V_per_m']/1e6
    X,Z=np.meshgrid(d['x_m'],d['z_m'])
    exact=np.array([.5*(analytic_y(X,-m['dy_m']/2,Z,tt,m)+analytic_y(X,m['dy_m']/2,Z,tt,m))/1e6 for tt in d['time_s']])
    assert np.max(abs(exact[0]-xz[0])) < np.max(abs(xz[0]))*1e-12
    ufull=np.asarray(h['U_full_J']/h['U_full_J'][0],float)
    uinside=np.asarray(h['U_interior_J']/h['U_interior_J'][0],float)
    # Interior curve retains its own documented initial interior-energy denominator.
    vmax=float(max(np.max(abs(xz[0])),np.max(abs(xy[0]))))
    extents=[[x[0],x[-1],z[0],z[-1]],[x[0],x[-1],z[0],z[-1]],
             [x[0],x[-1],y[0],y[-1]]]
    arrays=[xz,exact,xy]
    records=[];energy_only=mode=='energy-only'
    def export(fig, path, frames, update):
        if energy_only and path.name not in ['vacuum_shell_absorption.mp4','04_field_energy.mp4']:
            plt.close(fig)
            return
        update(frames[0]);fig.canvas.draw()
        # LogFormatterMathtext otherwise prints literal TeX when parse_math is disabled.
        for ax in fig.axes:
            objs=[ax.title,ax.xaxis.label,ax.yaxis.label]+ax.get_xticklabels()+ax.get_yticklabels()
            if ax.get_legend() is not None:objs+=list(ax.get_legend().texts)
            for text in objs:
                assert '$' not in text.get_text() and '_' not in text.get_text(),text.get_text()
        writer=FFMpegWriter(fps=20,codec='libx264',bitrate=1800,
                           extra_args=['-pix_fmt','yuv420p','-threads','4'])
        with writer.saving(fig,str(path),100):
            for k in frames:
                update(k)
                if k in [frames[0],frames[len(frames)//2],frames[-1]]:
                    fig.savefig(path.with_name(path.stem+f'_frame{k:04d}.png'),dpi=110)
                writer.grab_frame()
        plt.close(fig)
        probe=json.loads(subprocess.check_output([str(BIN/'ffprobe'),'-v','error',
              '-select_streams','v:0','-show_entries','stream=width,height,nb_frames,duration',
              '-of','json',str(path)],text=True))['streams'][0]
        assert int(probe['nb_frames'])==len(frames)
        assert abs(float(probe['duration'])-len(frames)/20)<1e-5
        assert int(probe['width'])%2==0 and int(probe['height'])%2==0
        records.append(dict(path=str(path.relative_to(root)),**probe,
                            actual_saved_frames=True,sha256=sha(path)))
        print('VIDEO',path.name,flush=True)
    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';target=out/lang;target.mkdir(exist_ok=True)
        titles=['xz 平面：仿真电场 Ey','xz 平面：理论电场 Ey','xy 平面：仿真电场 Ey','场能量随时间的变化'] if zh else [
                'xz plane: simulated Ey','xz plane: theoretical Ey','xy plane: simulated Ey','Field energy versus time']
        def image(ax,i):
            im=ax.imshow(arrays[i][0],origin='lower',extent=extents[i],cmap='RdBu_r',
                         vmin=-vmax,vmax=vmax,aspect='equal',interpolation='nearest')
            ax.set(title=titles[i],xlabel='x (µm)',ylabel=('z' if i<2 else 'y')+' (µm)')
            if i<2:
                zi=(256-40)*m['dz_m']*1e6
                ax.add_patch(Rectangle((-270,-zi),540,2*zi,fill=False,ls='--',lw=1,edgecolor='#555555'))
            else:
                ax.add_patch(Rectangle((-270,-270),540,540,fill=False,ls='--',lw=1,edgecolor='#555555'))
            return im
        def energy(ax):
            ax.semilogy(t,np.maximum(ufull,1e-12),color='#1769a0',label='全网格' if zh else 'Whole grid')
            ax.semilogy(t,np.maximum(uinside,1e-12),'--',color='#cf661f',label='吸收层内侧' if zh else 'Physical interior')
            cursor=ax.axvline(0,color='#151b23')
            ax.set(xlabel='模拟时间（ps）' if zh else 'Simulation time (ps)',
                   ylabel='场能量 / 各自初值' if zh else 'Field energy / respective initial value',
                   title=titles[3],ylim=(1e-10,2));ax.legend(frameon=False);ax.grid(alpha=.2)
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value,_:f'{value:g}'))
            ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            return cursor
        fig,axes=plt.subplots(2,2,figsize=(13.6,10.8))
        fig.subplots_adjust(left=.07,right=.94,bottom=.09,top=.88,hspace=.35,wspace=.37)
        ims=[image(axes.flat[i],i) for i in range(3)]
        for ax,im in zip(list(axes.flat)[:3],ims):fig.colorbar(im,ax=ax,shrink=.87,label='Ey (MV/m)')
        cursor=energy(axes[1,1]);suptitle=fig.suptitle('',y=.965,fontsize=18)
        fig.text(.5,.025,'虚线框内为物理区域；各时刻与两个切面使用同一色标' if zh else
                 'Dashed box: physical interior. One fixed color scale for both planes and all times.',ha='center',fontsize=11)
        def update_all(k):
            for i,im in enumerate(ims):im.set_data(arrays[i][k])
            cursor.set_xdata([t[k]]*2)
            suptitle.set_text(('真空传播与边界吸收' if zh else 'Vacuum propagation and boundary absorption')+f'  |  t = {t[k]:.3f} ps')
        export(fig,target/'vacuum_shell_absorption.mp4',list(range(len(t))),update_all)
        for i,name in enumerate(['01_xz_simulation','02_xz_theory','03_xy_simulation','04_field_energy']):
            fig,ax=plt.subplots(figsize=(9.6,7.2));fig.subplots_adjust(left=.14,right=.87,bottom=.14,top=.85)
            sup=fig.suptitle('',y=.96,fontsize=17)
            if i<3:
                im=image(ax,i);fig.colorbar(im,ax=ax,label='Ey (MV/m)',shrink=.90)
                def update(k,i=i,im=im,sup=sup):
                    im.set_data(arrays[i][k]);sup.set_text(f't = {t[k]:.3f} ps')
            else:
                cur=energy(ax)
                def update(k,cur=cur,sup=sup):
                    cur.set_xdata([t[k]]*2);sup.set_text(f't = {t[k]:.3f} ps')
            export(fig,target/(name+'.mp4'),list(range(len(t))),update)
        q=np.load(src/'cpml_time_fields.npz');tx=q['time_s']*1e12;xx=q['x_m']*1e6
        fig,ax=plt.subplots(figsize=(11,6.6));fig.subplots_adjust(left=.12,right=.97,bottom=.18,top=.80)
        l0,=ax.plot(xx,q['PEC_Ey_V_per_m'][0]/1e6,color='#cf661f',label='反射壁' if zh else 'Reflecting wall')
        l1,=ax.plot(xx,q['CPML_Ey_V_per_m'][0]/1e6,'--',color='#1769a0',label='吸收层' if zh else 'Absorbing layer')
        ax.axvspan(0,50,color='#dddddd');ax.axvspan(750,800,color='#dddddd')
        ax.set(xlabel='x (µm)',ylabel='Ey (MV/m)',ylim=(-1.15,1.15),
               title='沿 x 方向的电场剖面：吸收层与反射壁' if zh else 'Electric-field profile along x: absorber versus reflecting wall')
        ax.legend(frameon=False);ax.grid(alpha=.18);sup=fig.suptitle('',y=.95,fontsize=17)
        def update_cpml(k):
            l0.set_ydata(q['PEC_Ey_V_per_m'][k]/1e6);l1.set_ydata(q['CPML_Ey_V_per_m'][k]/1e6)
            sup.set_text(f't = {tx[k]:.3f} ps')
        export(fig,target/'cpml_vs_reflecting_wall.mp4',list(range(0,len(tx),5)),update_cpml)
    for p,expected in inputs.items():assert sha(p)==expected,p
    if energy_only:
        old=json.loads((out/'MOVIE_VALIDATION.json').read_text())
        replaced={v['path'] for v in records}
        records += [v for v in old['videos'] if v['path'] not in replaced]
        assert len(records)==len(old['videos'])==12
        for v in records:assert sha(root/v['path'])==v['sha256']
    save(out/'MOVIE_VALIDATION.json',dict(complete=True,videos=records,inputs=inputs,
         initial_theory_slice_match=True,field_color_scale_fixed_MV_per_m=vmax,
         wavefront_spherical_not_isotropic=True,source_tables_unchanged=True,
         corrected_plain_log_ticks=True,no_math_markup_in_titles_axes_or_ticks=True))


if __name__=='__main__':main(*sys.argv[1:])
