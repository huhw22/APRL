"""Phase-resolved origin of early bunching; closed energy-coupled data only."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c, e, m_e

from entrance_ballistic import save, table, digest, read_plane, phases, csum, lag


def mean(x, w):
    return np.sum(w*x)/np.sum(w)


def rms(x, w):
    return float(np.sqrt(mean(abs(x)**2, w)))


def correlation(x, y, w):
    x, y = x-mean(x,w), y-mean(y,w)
    return float(mean(x*y,w)/np.sqrt(mean(x*x,w)*mean(y*y,w)))


def complex_value(z):
    return dict(real=float(z.real), imag=float(z.imag), abs=float(abs(z)),
                phase_rad=float(np.angle(z)))


def reduce(values, weights, group, ng=None):
    if ng is None: ng = int(group.max()+1)
    W = np.bincount(group, weights=weights, minlength=ng)
    if np.iscomplexobj(values):
        total = csum(values*weights, group, ng)
    else:
        total = np.bincount(group, weights=values*weights, minlength=ng)
    return total/W, W


def summarize(p, w, groups):
    local, W = reduce(p,w,groups)
    b = mean(p,w)
    return dict(global_b=float(abs(b)), local_rms=rms(local,W),
        phase_alignment=float(abs(b)/mean(abs(local),W))), local


def stable_dgamma(u0,u1,g0,g1):
    return np.sum((u1-u0)*(u1+u0),axis=1)/(g1+g0)


def compare_energy_planes(directory):
    rows=[]
    for rawname, transport in [('entry','entry-transport'),('und-entry','und-entry-transport')]:
        a=read_plane(directory/(rawname+'-particles.h5'))
        u=a['proper_velocity'].astype(np.longdouble)
        gamma=np.sqrt(1+np.sum(u*u,axis=1))
        rest=a['mass_kg'].astype(np.longdouble)*np.longdouble(c)**2
        rawK=np.sum(rest*(gamma-1)); Ne=np.sum(a['weight'].astype(np.longdouble))
        rawQ=np.sum(a['charge_C'].astype(np.longdouble))
        with h5py.File(directory/(transport+'-energy.h5'),'r') as f:
            g=f['energy_plane'];assert int(g['complete'][()])==1
            r=g['records'][-1]
            assert int(r['particle_downstream_count'])==len(a)
            assert int(r['particle_upstream_count'])==0
            assert abs(np.longdouble(r['particle_charge_C'])-rawQ)/abs(rawQ)<1e-12
            K=np.longdouble(r['particle_kinetic_lab_J'])+np.longdouble(r['particle_kinetic_lab_J_low'])
        rows.append(dict(plane=rawname,count=len(a),
            raw_mean_kinetic_eV=float(rawK/Ne/e),
            scalar_mean_kinetic_eV=float(K/Ne/e),
            raw_minus_scalar_eV=float((rawK-K)/Ne/e)))
    return dict(rows=rows,
        raw_mean_change_eV=rows[1]['raw_mean_kinetic_eV']-rows[0]['raw_mean_kinetic_eV'],
        scalar_mean_change_eV=rows[1]['scalar_mean_kinetic_eV']-rows[0]['scalar_mean_kinetic_eV'],
        interpretation='Equal counts and charge. Raw whole-step three-momentum reconstruction and scalar substep four-momentum interpolation disagree. This establishes diagnostic inconsistency; it does not yet assign the harmonic spectrum or BF time error.')


def run(top, out):
    top, out = Path(top).resolve(), Path(out).resolve()
    out.mkdir(exist_ok=False)
    root = top/'compact-production-r1'
    directory = root/'energy/detectors'
    plan = json.loads((root/'plan.json').read_text())
    assert json.loads((root/'energy/PIC_COMPLETE.json').read_text())['run']['exit_code']==0
    source = top/'compact-preflight-r1/particles.h5'
    paths = [source,root/'plan.json',root/'energy/run.yaml',root/'energy/PIC_COMPLETE.json',
        directory/'entry-particles.h5',directory/'und-entry-particles.h5',
        directory/'p1-particles.h5',directory/'entry-transport-energy.h5',
        directory/'und-entry-transport-energy.h5']
    before = [digest(p) for p in paths]
    with h5py.File(source,'r') as f: a = f['particles/records'][:]
    n, clones = len(a), plan['NBINS']
    assert clones==16 and n==plan['input_particles']
    ids=a['source_id']; assert np.array_equal(ids,np.arange(1,n+1))
    w=a['macro_weight'].copy(); Q=w.sum()
    u0=a['proper_velocity'].copy(); g0=np.sqrt(1+np.sum(u0*u0,axis=1))
    assert np.all(u0.reshape(-1,clones,3)==u0.reshape(-1,clones,3)[:,:1])
    t0=a['arrival_time_offset_s'].copy(); anchor=float(np.median(t0))
    wl=plan['physics']['wavelength_m']; k=2*np.pi/wl
    K, period=plan['physics']['K'],plan['physics']['period_m']; ku=2*np.pi/period
    prep=.01; length=77*period; end=prep+length
    q0=1+u0[:,0]**2+u0[:,1]**2
    yentry=a['plane_position_m'][:,1]+prep*u0[:,1]/u0[:,2]
    I2=length/2-5*period/8
    gref=plan['physics']['gamma']
    index=np.arange(n); family=index//clones; ngfamily=n//clones
    group=index//(plan['NPART']*10); ngroup=group[-1]+1
    optical_slice=index//plan['NPART']
    basephase=phases(t0,anchor,wl)
    stage_rows=[]; moment_rows=[]; arrays={}; validation={}
    entry=read_plane(directory/'entry-particles.h5')
    assert np.array_equal(entry['source_id'],ids)
    offset=mean(entry['time_s'].astype(np.longdouble)-t0.astype(np.longdouble),w)
    validation['entry_time_nonconstant_max_s']=float(np.max(abs(
        entry['time_s'].astype(np.longdouble)-t0.astype(np.longdouble)-offset)))
    validation['entry_velocity_max_error']=float(np.max(abs(entry['proper_velocity']-u0)))
    del entry
    previous_gamma=g0
    phase_in=None; p1_data=None
    for name,z in [('und-entry',prep),('p1',end)]:
        p=read_plane(directory/(name+'-particles.h5'))
        assert np.array_equal(p['source_id'],ids) and np.array_equal(p['weight'],w)
        u=p['proper_velocity']; gamma=np.sqrt(1+np.sum(u*u,axis=1))
        dg=stable_dgamma(u0,u,g0,gamma)
        dq=(u[:,0]-u0[:,0])*(u[:,0]+u0[:,0])+(u[:,1]-u0[:,1])*(u[:,1]+u0[:,1])
        dangle=(u[:,0]**2+u[:,1]**2)/u[:,2]**2-(u0[:,0]**2+u0[:,1]**2)/u0[:,2]**2
        d,l=min(z,prep),max(0.,z-prep)
        if l:
            extra=K*K*I2/l*(1+ku*ku*yentry*yentry)
            flight=d*lag(g0,q0)+l*lag(g0,q0+extra)
        else: flight=z*lag(g0,q0)
        reference=d*lag(gref,1.)+l*lag(gref,1.+K*K/2)
        pred=basephase*np.exp(1j*k*(flight-reference))
        ds=np.asarray(np.longdouble(c)*(p['time_s'].astype(np.longdouble)
            -offset-t0.astype(np.longdouble))-np.longdouble(z)-flight.astype(np.longdouble),float)
        common=float(mean(ds,w))
        dp=k*(ds-common)
        actual=pred*np.exp(1j*dp)
        # Check against raw time phases with exactly the same common reference.
        raw=phases(p['time_s'],np.longdouble(anchor)+offset+np.longdouble((z+reference+common)/c),wl)
        validation[name+'_phase_reconstruction_max']=float(np.max(abs(raw-actual)))
        assert validation[name+'_phase_reconstruction_max'] < 1e-7
        phase_raw=phases(p['time_s'],np.longdouble(anchor)+offset+np.longdouble((z+reference)/c),wl)
        if phase_in is None: phase_in=phase_raw.copy()
        pf,FW=reduce(dp,w,family,ngfamily); within=dp-pf[family]
        gf,_=reduce(dg,w,family,ngfamily); gwithin=dg-gf[family]
        af,_=reduce(dangle,w,family,ngfamily); awithin=dangle-af[family]
        pcommon=pred*np.exp(1j*pf[family]); pwithin=pred*np.exp(1j*within)
        summaries={}
        for label,ph in [('reference',pred),('family_common_only',pcommon),
                         ('within_family_only',pwithin),('actual',actual)]:
            ss,local=summarize(ph,w,group);summaries[label]=ss
            stage_rows.append(dict(plane=name,z_m=z,phase_model=label,**ss))
            arrays[name+'_'+label+'_local']=local
        stage=dict(plane=name,z_m=z,common_delay_s=common/c,
            residual_phase_rms_rad=rms(dp,w),family_phase_rms_rad=rms(pf,FW),
            within_phase_rms_rad=rms(within,w),
            delta_gamma_mean=float(mean(dg,w)),delta_gamma_std=rms(dg-mean(dg,w),w),
            delta_energy_mean_eV=float(mean(dg,w))*m_e*c*c/e,
            delta_energy_std_eV=rms(dg-mean(dg,w),w)*m_e*c*c/e,
            within_delta_gamma_rms=rms(gwithin,w),
            corr_phase_energy=correlation(dp,dg,w),
            corr_within_phase_energy=correlation(within,gwithin,w),
            corr_within_phase_angle2=correlation(within,awithin,w),
            bunching=summaries)
        if name=='und-entry':
            columns=[]
            for harmonic in [10,20,30]:
                carrier=phase_in**harmonic
                for column in [carrier.real,carrier.imag]:
                    fm,_=reduce(column,w,family,ngfamily)
                    columns.append(column-fm[family])
            X=np.column_stack(columns)
            coefficients=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w*gwithin))
            fit=X@coefficients
            stage['harmonic_10_20_30_fit']=dict(
                within_energy_variance_fraction=float(1-mean((gwithin-fit)**2,w)/mean(gwithin**2,w)),
                residual_energy_rms_eV=rms(gwithin-fit,w)*m_e*c*c/e,
                limitation='Descriptive fit to raw-plane reconstructed energy, not independent physical work.')
            del X,columns
        for harmonic in range(1,33):
            carrier=phase_in**harmonic
            for label,value in [('delta_gamma_within',gwithin),('phase_residual_within',within)]:
                component=mean(value*carrier,w)
                local,_=reduce(value*carrier,w,group)
                moment_rows.append(dict(plane=name,harmonic=harmonic,quantity=label,
                    global_abs=float(abs(component)),global_real=float(component.real),
                    global_imag=float(component.imag),local_rms=rms(local,np.bincount(group,weights=w))))
        arrays[name+'_group_mean_delta_gamma']=reduce(dg,w,group)[0]
        arrays[name+'_group_phase_rms']=np.sqrt(reduce(dp*dp,w,group)[0])
        if name=='p1':
            p1_data=dict(stage=stage,pred=pred,actual=actual,dp=dp,within=within,
                dgamma=dg,gwithin=gwithin,dangle=dangle,awithin=awithin,
                delta_uz=u[:,2]-u0[:,2],dq=dq,gamma=gamma,
                u=u.copy(),common_phase=pf[family],pcommon=pcommon,pwithin=pwithin)
        save(out/(name+'_stage.json'),stage)
        del p
    d=p1_data;pred,actual,dp=d['pred'],d['actual'],d['dp']
    bref,bpic=mean(pred,w),mean(actual,w)
    diff=actual-pred
    b0local,W=reduce(pred,w,group);b1local,_=reduce(actual,w,group)
    r0,r1=rms(b0local,W),rms(b1local,W)
    # Signed exact magnitude changes, rather than sums of absolute contributions.
    global_gain=(w/Q)*np.real(diff*np.conj(bpic+bref))/(abs(bpic)+abs(bref))
    local_gain=(w/Q)*np.real(diff*np.conj((b1local+b0local)[group]))/(r1+r0)
    validation['global_gain_closure']=float(global_gain.sum()-(abs(bpic)-abs(bref)))
    validation['local_gain_closure']=float(local_gain.sum()-(r1-r0))
    assert abs(validation['global_gain_closure'])<1e-14
    assert abs(validation['local_gain_closure'])<1e-14
    family_common=mean(d['pcommon'],w);within_only=mean(d['pwithin'],w)
    # Symmetric ordering average avoids assigning a mixed contribution arbitrarily.
    common_shapley=.5*((abs(family_common)-abs(bref))+(abs(bpic)-abs(within_only)))
    within_shapley=.5*((abs(within_only)-abs(bref))+(abs(bpic)-abs(family_common)))
    attribution=dict(reference_global=complex_value(bref),actual_global=complex_value(bpic),
        global_BF_gain=float(global_gain.sum()),local_RMS_gain=float(local_gain.sum()),
        family_common_symmetric_amplitude_change=float(common_shapley),
        within_family_symmetric_amplitude_change=float(within_shapley),
        interpretation='Exact phase-rearrangement attribution, not force-source causality.')
    # Region and quantile contributions sum to the exact signed gain.
    contributions=[]
    def bins(label,bins,binlabels=None):
        for g in np.unique(bins):
            m=bins==g
            contributions.append(dict(grouping=label,bin=int(g),
                label=str(binlabels[int(g)]) if binlabels else str(int(g)),
                charge_fraction=float(w[m].sum()/Q),
                global_BF_gain=float(global_gain[m].sum()),
                global_gain_fraction=float(global_gain[m].sum()/global_gain.sum()),
                local_RMS_gain=float(local_gain[m].sum()),
                local_gain_fraction=float(local_gain[m].sum()/local_gain.sum()),
                delta_phase_mean_rad=float(mean(dp[m],w[m])),
                delta_phase_rms_rad=rms(dp[m],w[m]),
                delta_energy_mean_eV=float(mean(d['dgamma'][m],w[m]))*m_e*c*c/e))
    bins('longitudinal_10slice',group)
    bins('longitudinal_decile',np.minimum(optical_slice*10//plan['active_slices'],9))
    for axis,j in [('x',0),('y',1)]:
        val=a['plane_position_m'][:,j]
        order=np.argsort(val,kind='stable');cum=np.cumsum(w[order])/Q
        qb=np.empty(n,int);qb[order]=np.minimum((cum*10).astype(int),9)
        bins(axis+'_charge_decile',qb)
    # Top individual contributions: diagnostic concentration, never a new cohort.
    top_share={}
    for fraction in [.001,.01,.1]:
        count=int(np.ceil(n*fraction));sel=np.argpartition(abs(global_gain),-count)[-count:]
        top_share[str(fraction)]=dict(count=count,charge_fraction=float(w[sel].sum()/Q),
            signed_gain_fraction=float(global_gain[sel].sum()/global_gain.sum()),
            absolute_contribution_fraction=float(abs(global_gain[sel]).sum()/abs(global_gain).sum()))
    phase_bins=np.minimum((np.mod(np.angle(phase_in),2*np.pi)/(2*np.pi)*48).astype(int),47)
    phase_rows=[]
    for j in range(48):
        mask=phase_bins==j
        phase_rows.append(dict(phase_bin=j,phase_rad=(j+.5)*2*np.pi/48,
            charge_fraction=float(w[mask].sum()/Q),
            phase_within_mean_rad=float(mean(d['within'][mask],w[mask])),
            energy_within_mean_eV=float(mean(d['gwithin'][mask],w[mask]))*m_e*c*c/e,
            global_BF_gain=float(global_gain[mask].sum())))
    # Endpoint regressions are descriptive: no unobserved work history is inferred.
    x=d['gwithin']; y=d['within'];xx=x-mean(x,w);yy=y-mean(y,w)
    slope=float(mean(xx*yy,w)/mean(xx*xx,w))
    correlations=dict(within_phase_vs_energy=correlation(y,x,w),
        within_phase_vs_angle2=correlation(y,d['awithin'],w),
        phase_vs_delta_uz=correlation(dp,d['delta_uz'],w),
        endpoint_regression_rad_per_delta_gamma=slope,
        endpoint_regression_explained_variance=correlation(y,x,w)**2,
        limitation='Endpoint correlations and descriptive regression; no force or time-history causal fraction.')
    # Optical/release-step alias is expected by construction; do not call it grid locking.
    with h5py.File(directory/'entry-transport-energy.h5','r') as f:
        assert int(f['energy_plane/complete'][()])==1
        time_origin=float(f['energy_plane'].attrs['time_lab_origin_s'])
    entry_time=t0.astype(np.longdouble)+offset
    boost=plan['physics']['boost_gamma'];dt=plan['geometry']['dt_s']
    tau=np.asarray(np.longdouble(boost)*(entry_time-np.longdouble(time_origin))/np.longdouble(dt),float)
    optical=np.asarray((entry_time-np.longdouble(anchor)-offset)*np.longdouble(c/wl),float)
    rate=boost*wl/(c*dt)
    alias=dict(time_lab_origin_s=time_origin,steps_per_optical_period=float(rate),
        optical10_release1_coherence=float(abs(mean(np.exp(2j*np.pi*(tau-10*optical)),w))),
        interpretation='Release-step phase and optical phase are nearly harmonically aliased; correlation alone cannot identify numerical grid locking.')
    # Figure data are small; do not persist per-particle trajectories or phase arrays.
    hx=np.asarray(x*m_e*c*c/e/1e3,float)
    bounds=[np.quantile(hx,[.002,.998]),np.quantile(y,[.002,.998])]
    H,xe,ye=np.histogram2d(hx,y,bins=110,range=bounds,weights=w/Q)
    arrays.update(histogram_energy_phase=H,histogram_energy_edges_keV=xe,
        histogram_phase_edges_rad=ye,group_charge_e=W,
        group_global_BF_gain=np.bincount(group,weights=global_gain),
        group_local_RMS_gain=np.bincount(group,weights=local_gain))
    table(out/'stage_bunching.csv',stage_rows)
    table(out/'harmonics.csv',moment_rows)
    table(out/'contributions.csv',contributions)
    table(out/'phase_bins.csv',phase_rows)
    np.savez_compressed(out/'figure_data.npz',**arrays)
    detector_comparison=compare_energy_planes(directory)
    save(out/'detector_energy_comparison.json',detector_comparison)
    summary=dict(case='energy',all_input_particles=n,represented_electrons=float(Q),
        wavelength_m=wl,first_segment_end_m=end,attribution=attribution,
        p1=d['stage'],correlations=correlations,top_individual_contributions=top_share,
        release_phase_alias=alias,detector_energy_comparison=detector_comparison,no_new_simulation=True,
        reference_limit='Same averaged undulator and paraxial bulk focusing plus short-end approximation as previous entrance analysis; not exact 3D trajectories.',
        energy_limit='Raw whole-Maxwell-step momentum interpolation disagrees with independent substep scalar energy at the undulator entrance by about 645 eV per electron. Raw gamma correlations are diagnostic associations, not established physical energy modulation. Raw arrival time uses the same whole-step chord; its BF error is not quantified here.',
        missing_history='No trajectories and no particle plane between z=.01 and 2.166 m. Force-source causality and onset location remain unavailable.')
    save(out/'RESULT.json',summary)
    after=[digest(p) for p in paths];assert before==after
    validation['source_manifest']=before;validation['all_sources_unchanged']=True
    validation['only_energy_case_read']=True
    save(out/'VALIDATION.json',validation)
    render(out,summary,stage_rows,moment_rows,phase_rows,arrays,contributions)
    document(out,summary,stage_rows,moment_rows,contributions)
    save(out/'COMPLETE.json',dict(complete=True,case='energy',new_simulation=False))
    print(json.dumps(summary,indent=2),flush=True)


def render(out,result,stages,harmonics,phase_rows,data,contributions):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.7),layout='constrained')
    p=[r for r in stages if r['plane']=='p1'];labels=['Reference','Family common','Within family','Actual PIC']
    axes[0,0].semilogy(range(4),[r['local_rms'] for r in p],'o-',color='#24648c',label='Local RMS')
    axes[0,0].semilogy(range(4),[r['global_b'] for r in p],'s--',color='#b18428',label='Global')
    axes[0,0].set_xticks(range(4),labels,rotation=15)
    axes[0,0].set_ylabel('Bunching factor');axes[0,0].legend()
    axes[0,0].set_title('(a) Separate common and within-family phase changes')
    W=data['group_charge_e'];x=np.arange(len(W))
    axes[0,1].plot(x,np.cumsum(W)/W.sum(),':',color='#777777',label='Charge')
    for key,label,color,line in [('group_global_BF_gain','Global BF gain','#24648c','-'),
                               ('group_local_RMS_gain','Local RMS gain','#8e4b75','--')]:
        v=data[key];axes[0,1].plot(x,np.cumsum(v)/v.sum(),line,color=color,label=label)
    axes[0,1].set_xlabel('Initial longitudinal group (10 optical slices)')
    axes[0,1].set_ylabel('Cumulative signed fraction');axes[0,1].legend(fontsize=9)
    axes[0,1].set_title('(b) Where the new bunching is accumulated')
    H=data['histogram_energy_phase'];xe=data['histogram_energy_edges_keV'];ye=data['histogram_phase_edges_rad']
    im=axes[1,0].pcolormesh(xe,ye,np.ma.masked_less_equal(H.T,0),
        norm=LogNorm(vmin=max(H[H>0].min(),1e-8),vmax=H.max()),cmap='magma',shading='auto')
    xline=np.array([xe[0],xe[-1]])
    slope=result['correlations']['endpoint_regression_rad_per_delta_gamma']/(m_e*c*c/e/1e3)
    axes[1,0].plot(xline,slope*xline,color='#59b9c7',lw=1.4)
    axes[1,0].set_xlim(xe[0],xe[-1]);axes[1,0].set_ylim(ye[0],ye[-1])
    axes[1,0].set_xlabel('Within-family raw-plane energy change (keV)')
    axes[1,0].set_ylabel('Within-family arrival-phase residual (rad)')
    axes[1,0].set_title('(c) Endpoint association: r = '+f"{result['correlations']['within_phase_vs_energy']:.3f}")
    fig.colorbar(im,ax=axes[1,0],label='Charge fraction / bin')
    for plane,label,color,line in [('und-entry','Before undulator','#777777',':'),
                                 ('p1','After 77 periods','#24648c','-')]:
        rr=[r for r in harmonics if r['plane']==plane and r['quantity']=='delta_gamma_within']
        axes[1,1].semilogy([r['harmonic'] for r in rr],
            [r['global_abs']*m_e*c*c/e for r in rr],line+'o',ms=3,color=color,label=label)
    axes[1,1].set_xlabel('Harmonic of entrance optical phase')
    axes[1,1].set_ylabel('Raw-plane energy harmonic magnitude (eV)')
    axes[1,1].set_title('(d) Raw-plane energy structure at the two faces')
    axes[1,1].legend(fontsize=9)
    for ax in axes.flat:ax.grid(alpha=.18)
    fig.suptitle('Energy-coupled PIC: phase origin of first-segment bunching\n'
        'Full input; 16-clone families preserved; raw energy and arrival-time diagnostics require consistency checks',fontsize=12)
    fig.savefig(out/'energy_startup_origin.png',dpi=220,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(10.8,6.6),layout='constrained',sharex=True)
    angles=np.array([r['phase_rad'] for r in phase_rows])
    axes[0].plot(angles,[r['energy_within_mean_eV']/1e3 for r in phase_rows],'o-',color='#24648c')
    axes[1].plot(angles,[r['phase_within_mean_rad'] for r in phase_rows],'o-',color='#8e4b75')
    axes[0].set_ylabel('Mean raw energy change\nwithin family (keV)')
    axes[1].set_ylabel('Mean phase residual\nwithin family (rad)')
    axes[1].set_xlabel('Actual undulator-entrance optical phase (rad)')
    for ax in axes:ax.grid(alpha=.18);ax.axhline(0,color='#777777',lw=.8)
    fig.suptitle('Phase-resolved endpoint changes; energy-coupled case only',fontsize=13)
    fig.savefig(out/'phase_resolved_changes.png',dpi=220,bbox_inches='tight');plt.close(fig)


def document(out,result,stages,harmonics,contributions):
    text='# 能量守恒版：首段新增BF的相位来源\n\n'
    text+='仅使用原始输入和能量版已闭合的入口、磁入口、首段出口；本轮不读取显式版，不启动新仿真。\n\n'
    text+='现有过面记录中的新增BF，主要需要原本同动量的16个相位克隆之间产生差别位移；整组共同迟滞无法解释。束团内部贡献占主导。端点重建能量与到达相位存在较强负相关，但本轮同时查实粒子面能量与独立能量面的统计不一致，不能将这种相关直接解释为真实场做功。\n\n'
    dc=result['detector_energy_comparison']
    text+=f"**优先核对探测口径：入口到磁入口，原粒子面给出平均能量变化{dc['raw_mean_change_eV']:.3f} eV，同面独立标量能量面给出{dc['scalar_mean_change_eV']:+.3f} eV。磁入口两口径相差{dc['rows'][1]['raw_minus_scalar_eV']:.3f} eV/电子，粒子数与电荷一致。** 这已证明能量重建口径存在差异，但尚未测出原过面时间插值对BF的误差。\n\n"
    text+='![新增BF来源](energy_startup_origin.png)\n\n'
    text+='|首段末采用的相位变化|全局BF|10切片局部RMS|\n|---|---:|---:|\n'
    names={'reference':'无反馈参考','family_common_only':'仅每16粒子组共同位移',
        'within_family_only':'仅组内差别位移','actual':'实际能量版'}
    for r in stages:
        if r['plane']=='p1':text+=f"|{names[r['phase_model']]}|{r['global_b']:.9g}|{r['local_rms']:.9g}|\n"
    text+='\n## 精确相位分解\n\n'
    text+=r'对每个ID计算实际过面时刻减无反馈飞行预测，扣除一次全束公共迟滞。令参考相位为$p_i$，额外相位为$\Delta\phi_i$，则$\Delta b=\sum_i w_i p_i(e^{i\Delta\phi_i}-1)/Q$。所有分组共用同一相位锚点。'+'\n\n'
    text+=r'16粒子组内原始动量与横向位置逐位相同。把$\Delta\phi_i$分成组平均与组内偏差，分别构造相位反事实。对全局BF幅值的带符号粒子贡献使用$h_i=(w_i/Q)\operatorname{Re}[(p_i^{PIC}-p_i^{ref})(b_{PIC}+b_{ref})^*]/(|b_{PIC}|+|b_{ref}|)$，精确加和为幅值差；局部RMS采用对应分组二次恒等式。不能把粒子绝对贡献相加，也不能把相位反事实命名为受力因果份额。'+'\n\n'
    middle=[r for r in contributions if r['grouping']=='longitudinal_decile' and 1<=r['bin']<=8]
    text+=f"中间80%原始纵向索引范围含{100*sum(r['charge_fraction'] for r in middle):.3f}%电荷，贡献{100*sum(r['global_gain_fraction'] for r in middle):.3f}%全局新增幅值、{100*sum(r['local_gain_fraction'] for r in middle):.3f}%局部RMS新增幅值。\n\n"
    p1=result['p1'];co=result['correlations']
    text+=f"去公共迟滞后额外相位rms为{p1['residual_phase_rms_rad']:.6g} rad，其中组内部分{p1['within_phase_rms_rad']:.6g} rad。组内相位与组内能量变化的相关系数为{co['within_phase_vs_energy']:.6f}，与组内角度平方变化为{co['within_phase_vs_angle2']:.6f}。这些是端点关联，不是能量历史或做功分解。\n\n"
    text+='## 作用区间与谐波\n\n'
    text+='|已有面|额外相位rms/rad|组内能量变化rms/eV|一阶组内能量谐波模/eV|\n|---|---:|---:|---:|\n'
    for plane in ['und-entry','p1']:
        st=json.loads((out/(plane+'_stage.json')).read_text())
        h=next(r for r in harmonics if r['plane']==plane and r['harmonic']==1 and r['quantity']=='delta_gamma_within')
        text+=f"|{plane}|{st['residual_phase_rms_rad']:.8g}|{st['within_delta_gamma_rms']*m_e*c*c/e:.8g}|{h['global_abs']*m_e*c*c/e:.8g}|\n"
    early=json.loads((out/'und-entry_stage.json').read_text())
    text+='\n一阶谐波按实际磁入口相位计算。谐波模不是余弦峰值，纯正弦时峰值为其2倍。目标一阶结构在现有记录中主要于0.01–2.166 m之间发展；目前这段没有中间粒子面，不能确定更精确起始位置。磁入口重建能量主要包含10、20、30次谐波，三者描述性拟合解释了'+f"{100*early['harmonic_10_20_30_fit']['within_energy_variance_fraction']:.3f}%"+'的组内方差；结合独立能量面对照，不应将它直接当作真实前端能损。\n\n'
    text+='## 查实的探测口径差异\n\n'
    text+='|同一固定面|raw粒子面平均动能/eV|独立标量能量面平均动能/eV|raw减标量/eV|\n|---|---:|---:|---:|\n'
    for row in dc['rows']:
        text+=f"|{row['plane']}|{row['raw_mean_kinetic_eV']:.6f}|{row['scalar_mean_kinetic_eV']:.6f}|{row['raw_minus_scalar_eV']:.6f}|\n"
    text+='\n源码核对：`particle_energy_audit.h`逐粒子子步通知能量面，`simulation.cpp`提交这些子步过面事件；原粒子面则在整个Maxwell活动段两端之间，对位置、到达时刻和boost三动量作弦插值，再重建lab动量。日志记录每Maxwell步含7个粒子子步。\n\n'
    text+=r'同端点同插值比例$a$下，由凸性有$\gamma_{raw}-\gamma_{endpoint}=\gamma_b\{\sqrt{1+|(1-a)\mathbf u_0^\prime+a\mathbf u_1^\prime|^2}-(1-a)\gamma_0^\prime-a\gamma_1^\prime\}\le0$。实际两个探测器的端点跨度也不同，因此不能把645 eV全归给这个单独公式。'+'\n\n'
    text+='下一步最小核验应保持推进不变，让原粒子面与独立能量面采用同一子步事件，并在首段出口增加同口径对照，分别量化到达时间/BF与能量的诊断变化；本轮尚未修改探测器或启动该测试。\n\n'
    text+='![相位分箱](phase_resolved_changes.png)\n\n'
    text+='## 解释限制\n\n'
    text+='- 无反馈参考为上一轮的平均波荡器及近轴常幅作用量+短端部近似，不是原三维外场精确轨迹。\n'
    text+='- 粒子面动量经过过面插值，再变换回lab；由它重建的能量变化用于关联，不等于独立测得的严格逐粒子做功。\n'
    text+='- 到达相位依赖沿途能量与角度历史；端点回归不能单独证明哪项力造成聚束。\n'
    text+=f"- 每光学周期对应{result['release_phase_alias']['steps_per_optical_period']:.10f}个Maxwell步，释放步相位与光学第10谐波的相关幅值为{result['release_phase_alias']['optical10_release1_coherence']:.6f}。这属于设计性混淆，不能由相关性直接判网格锁定。\n"
    text+='- 相位克隆分裂本身既可能来自真实场反馈，也可能来自数值作用，本轮没有将其判为错误。\n\n'
    text+='数据：stage_bunching.csv、harmonics.csv、contributions.csv、phase_bins.csv、figure_data.npz、RESULT.json、VALIDATION.json。源文件读取前后SHA256保持。\n'
    (out/'README.md').write_text(text)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('top');p.add_argument('out')
    a=p.parse_args();run(a.top,a.out)
