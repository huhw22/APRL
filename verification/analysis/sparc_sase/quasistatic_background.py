"""Offline source preparation and coherent subtraction; never changes a PIC run."""
import gc
import hashlib
import json
import math
import struct
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml
from scipy import fft
from scipy.constants import c, e, m_e, mu_0, epsilon_0

from analyze_explicit import lagrange


def read(p):
    return json.loads(Path(p).read_text())


def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()) + '\n')


EVENT = np.dtype([(n, '<f8') for n in
    ['t', 'x', 'y', 'ux', 'uy', 'uz', 'q', 'mass', 'weight']] + [('id', '<u8')])
REQUEST = np.dtype([(n, '<f8') for n in ['t_box', 't_lab', 'z_box']])


def particles(path):
    with h5py.File(path, 'r') as f:
        g = f['particle_plane']
        assert int(g['complete'][()]) == 1
        assert int(g['committed_records'][()]) == len(g['records'])
        a = g['records'][:]
    a = a[np.argsort(a['source_id'])]
    assert np.all(np.diff(a['source_id']) > 0)
    assert np.all(a['crossing_direction'] == 1)
    return a


def box_beta(u, gamma_boost, beta_boost):
    gamma = np.sqrt(1 + np.sum(u*u, axis=1))
    plus = gamma + u[:, 2]
    minus = (1 + u[:, 0]**2 + u[:, 1]**2) / plus
    rapidity = gamma_boost * (1 + beta_boost)
    bp, bm = plus / rapidity, minus * rapidity
    return (bp-bm)/(bp+bm)


def prepare(top, out):
    top, out = Path(top).resolve(), Path(out).resolve()
    root = top / 'compact-production-r1'
    out.mkdir(exist_ok=False)
    plan = read(root/'plan.json')
    gamma = plan['physics']['boost_gamma']
    cfg = dict(root=str(out), nx=80, ny=80, dx=plan['geometry']['dx_m'],
      dy=plan['geometry']['dy_m'], dz=plan['geometry']['dz_m'],
      dt=plan['geometry']['dt_s'], boost_gamma=gamma, threads=8,
      physical_guard=10, guards=[10, 10, 0], block_layers=32,
      frame_z_origin_lab_m=0., source_rewind='common', plane_z_lab_m=15.246,
      cases={})
    manifest = {}
    for case, cache in [
        ('boris', top/'explicit_review_20261104/radiation_phase_r1'),
        ('energy', top/'energy_review_20261104/radiation_comparison_r1')]:
        assert read(root/case/'PIC_COMPLETE.json')['run']['exit_code'] == 0
        directory = out/case; directory.mkdir()
        source = particles(root/case/'detectors/far-particles.h5')
        w = source['weight'];q = source['charge_C']
        assert np.max(abs(q/(-e*w)-1)) < 2e-14
        assert np.max(abs(source['mass_kg']/(m_e*w)-1)) < 2e-14
        with h5py.File(root/case/'detectors/exit-transport-energy.h5', 'r') as f:
            g = f['energy_plane'];assert int(g['complete'][()]) == 1
            times = g['records'][:]
            t0 = float(g.attrs['time_lab_origin_s'])
            beta = float(g.attrs['boost_beta'])
        z = 15.246
        tb = times['time_box_s'];tl = times['time_lab_s']
        assert len(tb) == 23985 and np.all(np.diff(tb)>0)
        assert np.max(abs(tl - (t0+tb/gamma+beta*z/c))) < 2e-23
        # Independently recover the final frame translation from the original
        # exact safe-exit diagnostic; solver translation_bits is beta^2, not this.
        import re
        log = (root/case/'aprl.log').read_text()
        exit_lab = float(re.search(r'fixed_lab_particle_exit_m=([^ ]+)', log)[1])
        offset = int(re.search(r'initial_offset=(\d+)', log)[1])
        meta = (root/case/'reference/metadata.txt').read_text()
        bits = re.search(r'^origin_bits (.+)$', meta, re.M)[1].split()
        old_origin = [struct.unpack('>d',bytes.fromhex(x))[0] for x in bits]
        # Use long-double FMA equivalent, rounded once to binary64.
        fma = lambda a,b,d: float(np.longdouble(a)*np.longdouble(b)+np.longdouble(d))
        active_origin = fma(offset,cfg['dz'],old_origin[2])
        safe = fma(plan['geometry']['active_nz']-5,cfg['dz'],active_origin)
        box_origin = safe-exit_lab/gamma
        assert abs(box_origin - .08709711907273354) < 1e-15
        uzbeta = box_beta(source['proper_velocity'],gamma,beta)
        mass = source['mass_kg'].astype(np.longdouble)
        common = float(np.sum(mass*uzbeta)/np.sum(mass))
        event_tb = gamma*((source['time_s']-t0)-beta*z/c)
        event_zb = box_origin+gamma*(z-beta*c*(source['time_s']-t0))
        reference_time = float(np.sum(mass*event_tb)/np.sum(mass))
        cloud = event_zb+common*c*(reference_time-event_tb)
        zb = box_origin+z/gamma-beta*c*tb
        query = zb-common*c*(tb-reference_time)
        lower = min(float(cloud.min()),float(query.min()))
        upper = max(float(cloud.max()),float(query.max()))
        # Original longitudinal lattice phase, with the existing 64/128-cell
        # spatial bases. No new PIC window or simulation-time padding.
        low_index = math.floor((lower-old_origin[2])/cfg['dz'])-64
        origin_z = fma(low_index,cfg['dz'],old_origin[2])
        nz = int(math.ceil((upper-origin_z)/cfg['dz']))+128
        nz = int(math.ceil(nz/32))*32
        request = np.empty(len(tb),REQUEST)
        request['t_box']=tb;request['t_lab']=tl;request['z_box']=zb
        request.tofile(directory/'requests.bin')
        events = np.empty(len(source),EVENT)
        for dst, field in [('t','time_s'),('q','charge_C'),('mass','mass_kg'),('weight','weight'),('id','source_id')]:
            events[dst] = source[field]
        events['x']=source['position_m'][:,0];events['y']=source['position_m'][:,1]
        for axis,key in enumerate(['ux','uy','uz']):events[key]=source['proper_velocity'][:,axis]
        events.tofile(directory/'events.bin')
        with h5py.File(root/case/'detectors/far.h5','r') as f:
            g=f['field_plane'];assert int(g['complete'][()])==1
            raw_t=g['time_s'][:]
            assert np.array_equal(raw_t,tl[np.isin(tl,raw_t)])
            assert len(raw_t)==23984
            assert float(g.attrs['plane_z_m'])==z
        with np.load(cache/'far-spectra.npz') as a:
            assert np.max(abs((a['retarded_time_s']+z/c)-tl))<3e-23
        assert np.max(abs(events['x']))<.0012 and np.max(abs(events['y']))<.0012
        row = dict(nz=nz,origin=[*old_origin[:2],origin_z],
          frame_time_origin_lab_s=t0,frame_box_origin_m=box_origin,
          solver_beta=common,translation_beta=common,ref_time_box_s=reference_time,
          n_particles=len(source),n_samples=len(tb))
        cfg['cases'][case]=row
        coverage=float(abs(np.sum(q,dtype=np.longdouble))/plan['physics']['charge_C'])
        result=dict(case=case,cache=str(cache),source_particles=len(source),
          input_particles=plan['input_particles'],charge_coverage_fraction=coverage,
          source_charge_C=float(np.sum(q,dtype=np.longdouble)),
          beta_box_range=[float(uzbeta.min()),float(uzbeta.max())],
          weighted_beta_box_std=float(np.sqrt(np.sum(mass*(uzbeta-common)**2)/np.sum(mass))),
          grid_origin=cfg['cases'][case]['origin'],grid_nz=nz,
          cloud_z_range_m=[float(cloud.min()),float(cloud.max())],
          query_z_range_m=[float(query.min()),float(query.max())],
          frame_box_origin_m=box_origin,original_grid_phase_preserved=True,
          source_rewind='common; preserves every supplied crossing event',
          missing_tail_not_assumed_zero=True,
          scope='Observed far-crossing charge only; transverse event positions frozen; one common longitudinal velocity')
        save(directory/'INPUT.json',result)
        for p in [directory/'events.bin',directory/'requests.bin']:
            manifest[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
        del source,events
    (out/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    manifest[str(out/'config.yaml')]=hashlib.sha256((out/'config.yaml').read_bytes()).hexdigest()
    save(out/'INPUT_FROZEN.json',manifest)
    prepare_controls(out,cfg)
    print(json.dumps({k:dict(v) for k,v in cfg['cases'].items()},indent=2))


def prepare_controls(out,cfg):
    test = dict(cfg);test['root']=str(out/'controls')
    test.update(nx=24,ny=24,nz=128,physical_guard=4,guards=[4,4,0],
      origin=[-12*cfg['dx'],-12*cfg['dy'],-64*cfg['dz']],
      frame_time_origin_lab_s=0.,frame_box_origin_m=0.,plane_z_lab_m=0.,
      ref_time_box_s=0.,n_particles=3,n_samples=17,cases={})
    directory=out/'controls';directory.mkdir()
    gb=cfg['boost_gamma'];bb=np.sqrt((1-1/gb)*(1+1/gb))
    u=np.array([[0.,0.,np.sqrt(298.2608843802298**2-1)]]*3)
    beta=float(box_beta(u,gb,bb)[0])
    test.update(solver_beta=beta,translation_beta=beta)
    for case,scale in [('base',1),('double',2)]:
        d=directory/case;d.mkdir()
        events=np.zeros(3,EVENT);events['id']=[1,2,3]
        events['t']=np.array([-2,0,2])*cfg['dt']/gb
        events['x']=[-cfg['dx']*.5,cfg['dx']*.5,0]
        events['y']=[0,cfg['dy']*.5,-cfg['dy']*.5]
        events['uz']=u[:,2];events['weight']=np.array([2,3,5])*scale
        events['q']=-e*events['weight'];events['mass']=m_e*events['weight']
        events.tofile(d/'events.bin')
        req=np.empty(17,REQUEST);req['t_box']=np.arange(-8,9)*cfg['dt']
        req['t_lab']=req['t_box']/gb;req['z_box']=-bb*c*req['t_box']
        req.tofile(d/'requests.bin')
        test['cases'][case]={}
    (out/'controls.yaml').write_text(yaml.safe_dump(test,sort_keys=False))


def controls(out):
    out=Path(out)
    a=np.fromfile(out/'controls/base/fields.bin',float)
    b=np.fromfile(out/'controls/double/fields.bin',float)
    assert len(a)==len(b)==17*16*16*6
    assert np.isfinite(a).all() and np.max(abs(a))>0
    error=float(np.max(abs(b-2*a))/np.max(abs(2*a)))
    assert error<5e-12
    for case in ['base','double']:
        r=read(out/'controls'/case/'COMPLETE.json')
        assert r['complete'] and r['deposition_relative_charge_error']<1e-10
        assert r['original_event_consistency_zmax_m']<1e-13
    save(out/'CONTROLS.json',dict(linearity_relative_error=error,passed=True,
      scope='Offline reader and source scaling control, not physical radiation acceptance'))


def reduce_case(top,out,case):
    top,out=Path(top),Path(out);directory=out/case
    complete=read(directory/'COMPLETE.json');assert complete['complete']
    inp=read(directory/'INPUT.json');cache=Path(inp['cache'])
    original=top/'compact-production-r1'/case/'detectors/far.h5'
    p=read(cache/'pic.json')[-1]
    if case=='energy':p=read(cache/'pic.json')[-1]
    N=p['fft_samples'];dt=p['dt_s'];dx=4e-5;area=dx*dx
    assert list(complete['shape'])==[N,60,60,6]
    assert (directory/'fields.bin').stat().st_size==N*3600*6*8
    background=np.memmap(directory/'fields.bin',dtype=float,mode='r',shape=(N,3600,6))
    with h5py.File(original,'r') as f:
        g=f['field_plane'];assert int(g['complete'][()])==1
        t=g['time_s'][:];ticks=np.rint((t-t[0])/dt).astype(int)
        total=np.empty((N,3600,6))
        for a in range(0,len(t),32):
            b=min(len(t),a+32);v=g['fields'][a:b]
            total[ticks[a:b],:,:3]=v['electric_V_per_m'][:,10:70,10:70,:].reshape(b-a,-1,3)
            total[ticks[a:b],:,3:]=v['magnetic_T'][:,10:70,10:70,:].reshape(b-a,-1,3)
        missing=np.setdiff1d(np.arange(N),ticks)
        for k in missing:
            neighbors=ticks[np.argsort(abs(ticks-k))[:8]]
            total[k]=np.tensordot(lagrange(neighbors,k),total[neighbors],axes=(0,0))
    assert np.isfinite(total).all()
    energy=read(cache/'far-summary.json')
    assert energy['fft_samples']==N
    with np.load(cache/'far-spectra.npz') as a:
        frequency=a['energy_eV'].copy()
    de=frequency[1]-frequency[0];one=np.full(len(frequency),2.);one[0]=1.
    if N%2==0:one[-1]=1.
    band=(frequency>=2.397330687668754)&(frequency<=2.649681286370728)
    arrays={'energy_eV':frequency};summary={}
    for window,win in [('rect',np.ones(N)),('hann',np.hanning(N))]:
        win=win/np.sqrt(np.mean(win*win))
        sums={k:np.zeros(len(frequency)) for k in ['total','background','residual','cross','total_E','background_E','residual_E']}
        direct={k:np.longdouble(0.) for k in ['total','background','residual']}
        for start in range(0,3600,64):
            stop=min(start+64,3600)
            a=np.array(total[:,start:stop,:].transpose(1,2,0),order='C',copy=True)
            b=np.array(background[:,start:stop,:].transpose(1,2,0),order='C',copy=True)
            a*=win[None,None,:];b*=win[None,None,:]
            for key,v in [('total',a),('background',b),('residual',a-b)]:
                direct[key]+=np.sum(v[:,0,:]*v[:,4,:]-v[:,1,:]*v[:,3,:],dtype=np.longdouble)*dt*area/mu_0
            A=fft.rfft(a,axis=-1,workers=8,overwrite_x=True)*dt
            B=fft.rfft(b,axis=-1,workers=8,overwrite_x=True)*dt
            R=A-B
            for key,F in [('total',A),('background',B),('residual',R)]:
                sums[key]+=np.sum((F[:,0,:]*F[:,4,:].conj()-F[:,1,:]*F[:,3,:].conj()).real,axis=0)
                sums[key+'_E']+=np.sum(abs(F[:,:2,:])**2,axis=(0,1))
            sums['cross']+=np.sum((R[:,0,:]*B[:,4,:].conj()+B[:,0,:]*R[:,4,:].conj()
                -R[:,1,:]*B[:,3,:].conj()-B[:,1,:]*R[:,3,:].conj()).real,axis=0)
        for key in sums:
            scale=one*area/(N*dt*de)
            sums[key]*=scale*(epsilon_0*c if key.endswith('_E') else 1/mu_0)
            arrays[window+'_'+key+'_J_per_eV']=sums[key]
        scale=max(np.max(abs(sums['total'])),np.max(abs(sums['background'])),1e-300)
        closure=float(np.max(abs(sums['total']-sums['background']-sums['residual']-sums['cross']))/scale)
        assert closure<2e-12
        parseval={}
        for key in ['total','background','residual']:
            diff=abs(float(sums[key].sum()*de-direct[key]))
            parseval[key]=float(diff/max(abs(direct[key]),np.sum(abs(sums[key]))*de,1e-300))
            assert parseval[key]<2e-10
        with np.load(cache/'far-spectra.npz') as previous:
            for key,oldkey in [('total','native_flux'),('total_E','E_proxy')]:
                ref=previous[window+'_'+oldkey+'_J_per_eV']
                assert np.max(abs(sums[key][band]-ref[band]))/np.max(abs(ref[band]))<5e-10
        values={k:float(v[band].sum()*de) for k,v in sums.items()}
        assert values['total']>0 and values['total_E']>0 and values['background_E']>=0
        q=values['background_E']/values['total_E']
        summary[window]=dict(band_energy_J=values,
          record_mean_W={k:v/(N*dt) for k,v in values.items()},
          background_E_proxy_fraction=q,E_proxy_worst_phase_relative_change_bound=2*np.sqrt(q)+q,
          residual_over_total=values['residual']/values['total'],
          relative_subtraction_change=values['residual']/values['total']-1,
          spectral_decomposition_relative_error=closure,parseval_relative=parseval)
        print(case,window,json.dumps(summary[window]),flush=True)
    np.savez_compressed(directory/'background_spectra.npz',**arrays)
    save(directory/'SPECTRUM.json',dict(case=case,results=summary,record_ps=N*dt*1e12,
      input=inp,near_field_estimate=True,unique_radiation_separation=False,
      scope='Raw total native E/B minus a conditional finite-domain Cowan quasi-static reference'))
    del total,background;gc.collect()


def summarize(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);rows=[]
    plt.rcParams.update({'font.size':12,'axes.titlesize':13,'legend.fontsize':10,
      'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    for case,title in [('boris','Explicit'),('energy','Energy coupled')]:
        directory=out/case;result=read(directory/'SPECTRUM.json')
        T=result['record_ps']*1e-12
        with np.load(directory/'background_spectra.npz') as a:
            x=a['energy_eV']
            fig,axes=plt.subplots(2,2,figsize=(12.4,8.2),layout='constrained')
            for i,window in enumerate(['rect','hann']):
                v=result['results'][window]
                for j in range(2):
                    ax=axes[i,j];lim=(0,5) if j==0 else (2.397330687668754,2.649681286370728)
                    m=(x>=lim[0])&(x<=lim[1])
                    hi=max(float(np.max(a[window+'_'+k+'_J_per_eV'][m]/T)) for k in ['total','background','residual'])
                    for key,label,color,style in [
                        ('total','PIC total E x B','#24648c','-'),
                        ('background','Quasi-static estimate','#777777',':'),
                        ('residual','PIC minus reference','#d47928','--')]:
                        y=a[window+'_'+key+'_J_per_eV']/T
                        ax.semilogy(x,np.where(y>0,y,np.nan),style,color=color,label=label)
                    ax.set_xlim(lim);ax.set_ylim(max(hi*1e-12,1e-25),hi*2)
                    ax.axvline(2.523505987019741,color='#777777',ls=':',lw=.8)
                    ax.set_xlabel('Photon energy (eV)');ax.set_ylabel('Record-mean spectrum (W/eV)')
                    ax.set_title(f'({chr(97+2*i+j)}) '+('Rectangular' if i==0 else 'Energy-normalized Hann')+
                        (': wide view' if j==0 else ': resonant band'))
                    ax.grid(alpha=.18)
                    if i==0 and j==0:ax.legend()
                rows.append(dict(case=case,window=window,
                  raw_total_W=v['record_mean_W']['total'],
                  background_W=v['record_mean_W']['background'],
                  residual_W=v['record_mean_W']['residual'],
                  cross_W=v['record_mean_W']['cross'],
                  subtraction_change_pct=100*v['relative_subtraction_change'],
                  background_E_proxy_pct=100*v['background_E_proxy_fraction'],
                  E_proxy_worst_phase_change_pct=100*v['E_proxy_worst_phase_relative_change_bound'],
                  charge_coverage_pct=100*result['input']['charge_coverage_fraction']))
            fig.suptitle(title+': conditional quasi-static background at 15.246 m',fontsize=15)
            fig.savefig(out/(case+'_background_spectra.png'),dpi=210,bbox_inches='tight',pad_inches=.12);plt.close(fig)
    import csv
    with (out/'BACKGROUND_COMPARISON.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    save(out/'RESULT.json',dict(rows=rows,physical_acceptance=False,
      unique_radiation_separation=False,scope='Conditional finite-domain common-beta Cowan reference only'))
    hann={row['case']:row for row in rows if row['window']=='hann'}
    before=hann['energy']['raw_total_W']/hann['boris']['raw_total_W']
    after=hann['energy']['residual_W']/hann['boris']['residual_W']
    ratios={}
    for window in ['rect','hann']:
        a=read(out/'boris/SPECTRUM.json')['results'][window]
        b=read(out/'energy/SPECTRUM.json')['results'][window]
        q_a=a['background_E_proxy_fraction'];q_b=b['background_E_proxy_fraction']
        base=b['band_energy_J']['total_E']/a['band_energy_J']['total_E']
        ratios[window]=dict(E_proxy_before_ratio=base,
          unknown_phase_E_proxy_ratio_lower_bound=base*((1-np.sqrt(q_b))/(1+np.sqrt(q_a)))**2)
    save(out/'RATIO_BOUNDS.json',ratios)
    text=r'''# 出口准静态背景能谱估算

以两版各自的真实出口粒子过面记录重建刚性准静态参考。没有重新推进PIC，也没有用平滑高斯替代实际微聚束和横向分布。
'''
    text+=f"\n本次参考下，目标频带Hann背景约为显式版{hann['boris']['background_W']:.3f} W、能量版{hann['energy']['background_W']:.3f} W；场级扣除使各自功率仅下降{-hann['boris']['subtraction_change_pct']:.6f}%、{-hann['energy']['subtraction_change_pct']:.6f}%。两版倍率{before:.6f}→{after:.6f}，几乎不变。因此这种出口准静态束缚场不足以解释约6.5倍的推进差异。此结论限于本次条件性参考，不代表所有非辐射背景已经排除。\n"
    text+=r'''
## 同一探测面上的比较

下表都是2.397331–2.649681 eV频带、2.4×2.4 mm²孔径及同一3.930766 ps记录的平均功率；背景和余场均为有符号原生E×B口径。交叉项可以为负，不能把两条功率谱直接相减当作余场。对数图显示正值，带符号原始频谱完整保存在NPZ，积分没有截去负值。

|版本|窗口|PIC总场|准静态背景|场级相减后|交叉项|相减变化|
|---|---|---:|---:|---:|---:|---:|
'''
    for row in rows:
        text+=f"|{row['case']}|{row['window']}|{row['raw_total_W']/1000:.6g} kW|{row['background_W']/1000:.6g} kW|{row['residual_W']/1000:.6g} kW|{row['cross_W']/1000:.6g} kW|{row['subtraction_change_pct']:+.6g}%|\n"
    text+=r'''
![显式版四宫格](boris_background_spectra.png)

![能量耦合版四宫格](energy_background_spectra.png)

## 为什么在复场层面相减

对同一像素和同一时间权重做傅里叶变换：

$$
\widetilde{\mathbf E}_{r}=\widetilde{\mathbf E}_{t}-\widetilde{\mathbf E}_{q},
\qquad
\widetilde{\mathbf B}_{r}=\widetilde{\mathbf B}_{t}-\widetilde{\mathbf B}_{q}.
$$

$$
S_t=S_r+S_q+S_{\rm cross},\qquad
S_{\rm cross}\propto\Re[(\widetilde{\mathbf E}_r\times\widetilde{\mathbf B}_q^*
 +\widetilde{\mathbf E}_q\times\widetilde{\mathbf B}_r^*)_z].
$$

下标t、q、r分别表示总场、准静态参考、余场。逐像素计算后才对全平面积分；保留矩形窗和能量归一Hann窗，DC也保持。四项由各自复场直接算出，并独立核对恒等式，不用剩余差额定义交叉项。

另报横向E²正定代理的背景比例q及最坏相位界限：

$$\left|\frac{W_{E,r}}{W_{E,t}}-1\right|\le 2\sqrt q+q,\qquad q=\frac{W_{E,q}}{W_{E,t}}.$$

这个界限只针对正定横向E²代理，不冒充一般有符号Poynting能流的界限。实际Poynting变化已经通过保留相位的场级相减计算。

## 参考模型与采样

固定lab面给出的是不同时刻的过面事件，不能把所有粒子直接当成同时快照。将事件变到原boost坐标，选共同参考时刻t'_0，按共同纵向速度beta_q回推：

$$z'_{i,0}=z'_i+\beta_q c(t'_0-t'_i).$$

再按同一beta_q平移。这使每个输入事件的过面位置和相位原样保留。若先按各自速度回推、再按平均速度传播，会额外涂抹微聚束，本次没有采用那种做法。实际宏粒子电荷已含权重，不再次乘weight。

使用既有BlockedCowanReference一次求解，不通过长时间传播等待稳定。保持生产横向/纵向网格分辨率、线性形函数、有限Dirichlet横向边界；纵向参考范围仅按源云及实际查询范围推导。保留原boost坐标、原格点相位、原生E整步/B存储半步、参考平移插值及Yee探测插值，随后变回lab。actual time_box_s直接来自完整exit-transport记录，未用近似时间原点。

该参考的连续物理对应是共同速度的束缚场，而不是加速度辐射。匀速速度场及洛伦兹变换参见[Feynman II第26章](https://www.feynmanlectures.caltech.edu/II_26.html)。不能仅用E≈cB判断辐射，因为近光速束缚场也可满足此关系。高频背景需保留横向结构，不能仅由全局BF决定；参见[Venturini 2008](https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.11.034401)。

## 解释范围

- 两版far出口事件覆盖约99.949%和99.960%的输入电荷；缺的是低权重尾部，未把它当零电荷，也未声称完整重建全部束团。
- 参考保持各粒子的过面位置、到达时间及电荷；采用共同纵向速度，横向位置冻结。它不包含真实加速度历史或横向电流的全部效应。
- 这是与现有Cowan离散及有限横向边界匹配的背景估计；没有在本轮宣称独立开边界连续解或网格收敛。491 nm速度场的特征横向尺度约23 μm，小于40 μm网格，必须保留这一限制。
- 余场不能直接改名为“纯物理辐射”。即使该准静态背景不足以解释差异，启动数值波、边界残留和其它推进误差仍需另查。

程序控制和数值证据保存在CONTROLS.json、各版INPUT/SOURCE/REFERENCE/COMPLETE/SPECTRUM.json。源电荷、事件回代、完整时间记录、原谱复算、Parseval及四项分解分别检查，程序检查不冒充物理模型验收。
'''
    (out/'README.md').write_text(text)
    save(out/'REDUCTION_COMPLETE.json',dict(complete=True,physical_acceptance=False,
      rows=rows,images_visually_inspected=False))
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare(*sys.argv[2:4])
    elif action=='controls':controls(sys.argv[2])
    elif action=='reduce':reduce_case(*sys.argv[2:5])
    elif action=='summarize':summarize(sys.argv[2])
