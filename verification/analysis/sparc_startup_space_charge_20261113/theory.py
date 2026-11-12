"""Parameter-based LSC estimates for the closed random-startup campaign.

No PIC/GENESIS tracking, fitted force coefficient, power renormalization or
large field-file read. Cold coherent-mode transfers are not SASE RMS predictions.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import c,e,m_e,epsilon_0
from scipy.integrate import quad,solve_ivp
from scipy.linalg import expm
from scipy.special import exp1,roots_laguerre

W=Path('/public/home/wfollowm/APRL_Workspace')


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,
        default=lambda x:x.item() if isinstance(x,np.generic) else x)+'\n')


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def form_factor(k,gamma_z,sx,sy,order=48):
    u,w=roots_laguerre(order);a=k/gamma_z
    return float(np.sum(w*sx*sy/np.sqrt((sx*sx+u/a**2)*(sy*sy+u/a**2))))


def cold_fel(L,g,kp,H=1.):
    # A'=g H b; b'=-i g P; P'=-g H A-i kp² b/g.
    matrix=np.array([[0,g*H,0],[0,0,-1j*g],[-g*H,-1j*kp**2/g,0]],complex)
    return expm(matrix*L)@np.array([0,1,0],complex)


def main(out):
    out=Path(out);out.mkdir(exist_ok=False)
    campaign=W/'results/sparc_random_startup_20261112_r2'
    source=campaign/'PLAN.json';p=json.loads(source.read_text());s=p['source_config']['phase_space']
    I=p['source_config']['profile']['peak_current_A'];G=s['gamma_mean'];K=p['K'];lam=p['wavelength_m']
    sx,sy=s['sigma_x_m'],s['sigma_y_m'];Lu=p['period_m'];ku=2*np.pi/Lu;k=2*np.pi/lam;L=77*Lu
    IA=4*np.pi*epsilon_0*m_e*c**3/e;A=1+K*K/2;Gz=G/np.sqrt(A)
    kap=k*A/G**3;R=form_factor(k,Gz,sx,sy)
    S=I/IA/(k*sx*sy)*R;kp=np.sqrt(kap*S)
    beta=np.sqrt(1-1/Gz**2);dz=Gz*(1+beta)*lam/20
    theta_sc=Gz*(k/beta-beta*k)*dz;theta_rad=Gz*k*(1-beta)*dz
    Fsc=np.cos(theta_sc/2)**2;Frad=np.cos(theta_rad/2)**2
    Hrad=Frad*np.sinc(theta_rad/(2*np.pi))**2
    estimates=[]
    for label,gkernel in [('straight-beam kernel',G),('period-average kernel',Gz)]:
        r=form_factor(k,gkernel,sx,sy);omega=np.sqrt(I/IA*A/(G**3*sx*sy)*r)
        estimates.append(dict(model=label,transverse_R=r,kp_per_m=omega,plasma_phase_rad=omega*L,
            cold_density_transfer=float(np.cos(omega*L)),
            cold_integrated_radiation_energy_ratio=float(np.sinc(omega*L/np.pi)**2)))
    ex=s['emittance_n_x_m']/G;ey=s['emittance_n_y_m']/G
    def covariance(sig,eps,alpha):
        return np.array([[sig**2,-eps*alpha],[-eps*alpha,eps**2/sig**2*(1+alpha**2)]])
    Cx=covariance(sx,ex,s['alpha_x']);Cy=covariance(sy,ey,s['alpha_y'])
    assert abs(np.linalg.det(Cx)/ex**2-1)<1e-12 and abs(np.linalg.det(Cy)/ey**2-1)<1e-12
    def kval(z):
        a=min(z,L-z)
        return K*(.5*(1-np.cos(np.pi*a/Lu)) if a<Lu else 1.)
    def coefficients(z,M):
        kk=kval(z);aa=1+kk*kk/2
        row=np.array([1.,z]);sigx=np.sqrt(row@Cx@row);sigy=np.sqrt((M@Cy@M.T)[0,0])
        rr=form_factor(k,G/np.sqrt(aa),sigx,sigy)
        return kk,k*aa/G**3,I/IA/(k*sigx*sigy)*rr,sigx,sigy
    def rhs(z,state):
        M=state[:4].reshape(2,2);kk,kappa,force,_,_=coefficients(z,M)
        J=np.array([[0,1],[-kk**2*ku**2/(2*G**2),0]])
        # b,h are two cold coherent-mode density/energy quadratures, g=-i h.
        result=[*(J@M).ravel()]
        for i,factor in [(4,1.),(7,Fsc**2)]:
            b,h,integral=state[i:i+3]
            result.extend([-kappa*h,force*factor*b,b])
        return result
    trajectory=solve_ivp(rhs,[0,L],[1,0,0,1,1,0,0,1,0,0],rtol=2e-10,atol=1e-12,max_step=.01,dense_output=True)
    assert trajectory.success
    z=np.linspace(0,L,301);v=trajectory.sol(z);curves=[]
    for j,x in enumerate(z):
        _,_,_,xx,yy=coefficients(x,v[:4,j].reshape(2,2))
        curves.append(dict(z_m=x,sigma_x_m=xx,sigma_y_m=yy,
            constant_density_transfer=np.cos(kp*x),envelope_density_transfer=v[4,j],
            envelope_F3_only_density_transfer=v[7,j],
            constant_coherent_energy_ratio=1. if x==0 else np.sinc(kp*x/np.pi)**2,
            envelope_coherent_energy_ratio=1. if x==0 else (v[6,j]/x)**2,
            envelope_F3_only_coherent_energy_ratio=1. if x==0 else (v[9,j]/x)**2))
    cold=[]
    for Lg in [.4466785032352482,.67]:
        gain=1/(np.sqrt(3)*Lg);base=cold_fel(L,gain,0)
        for name,omega,H in [('LSC only',kp,1.),('source/gather shape only',0,Hrad),
            ('LSC plus approximate F3/source/gather',kp*Fsc,Hrad)]:
            result=cold_fel(L,gain,omega,H)
            roots=np.roots([1,0,omega**2,-1j*gain**3*H**2])
            cold.append(dict(model=name,reference_power_gain_length_m=Lg,gain_scale_per_m=gain,
                density_amplitude_ratio=float(abs(result[1]/base[1])),radiation_energy_ratio=float(abs(result[0]/base[0])**2),
                growing_eigenvalue_real_per_m=float(max(roots.real)),
                assumptions='Cold, on-resonance, single coherent density mode, zero initial field and energy harmonic. Effective0.67m case is illustrative; not a 3D derivation.'))
    # Independent numerical checks of analytic limits and an alternative quadrature.
    alternate=quad(lambda u:np.exp(-u)*sx*sy/np.sqrt((sx*sx+u/(k/Gz)**2)*(sy*sy+u/(k/Gz)**2)),0,np.inf,epsabs=1e-13)[0]
    sig=np.sqrt(sx*sy);X=(k*sig/Gz)**2
    round_error=abs(form_factor(k,Gz,sig,sig)-X*np.exp(X)*exp1(X))
    solution=solve_ivp(lambda x,y:[-kap*y[1],S*y[0]],[0,L],[1.,0.],rtol=1e-11,atol=1e-13)
    cos_error=abs(solution.y[0,-1]-np.cos(kp*L))
    assert abs(R-alternate)<1e-11 and round_error<1e-11 and cos_error<1e-10
    assert abs(theta_rad-2*np.pi/20)<1e-10
    b0=.01
    angular_x=(s['emittance_n_x_m']/sx)**2*(1+s['alpha_x']**2)
    angular_y=(s['emittance_n_y_m']/sy)**2*(1+s['alpha_y']**2)
    angle_a=k*L/(2*G**2)
    result=dict(physical_acceptance=False,calculation_only=True,no_new_PIC_or_GENESIS=True,
        inputs=dict(current_A=I,gamma=G,K=K,sigma_x_m=sx,sigma_y_m=sy,wavelength_m=lam,undulator_length_m=L),
        constants=dict(Alfven_current_A=IA,gamma_z=Gz,k_per_m=k,kappa_per_m_per_gamma=kap),
        fixed_size_estimates=estimates,
        varying_envelope=dict(density_ratio=float(v[4,-1]),coherent_energy_ratio=float((v[6,-1]/L)**2),
            F3_only_density_ratio=float(v[7,-1]),F3_only_coherent_energy_ratio=float((v[9,-1]/L)**2)),
        field_scale=dict(reference_complex_b=b0,Ez_complex_coefficient_V_per_m_per_b=S*m_e*c*c/e,
            Ez_real_peak_for_reference_b_V_per_m=2*S*m_e*c*c/e*b0,
            perturbative_energy_harmonic_eV=S*L*b0*m_e*c*c/e),
        shape_symbols=dict(static_density_cells_per_wave=2*np.pi/theta_sc,radiation_cells_per_wave=2*np.pi/theta_rad,
            F3_source_and_gather_SC_factor=Fsc**2,F3_radiation_source_power_factor=Frad**2,
            F3_plus_linear_tangential_source_power_factor=Hrad**2,
            scope='Exact F3 symbol; linear tangential source main band only. No aliases, full discrete Maxwell Green tensor, transverse grid or detector transfer are included. No power correction is applied.'),
        common_dephasing=dict(energy_spread_free_density_ratio=float(np.exp(-.5*(kap*s['gamma_sigma']*L)**2)),
            straight_constant_angle_density_ratio=float(((1+(2*angle_a*angular_x)**2)*(1+(2*angle_a*angular_y)**2))**(-.25)),
            scope='Already represented in both tracking models; ballistic constant-angle number is a scale, not an additional missing PIC-only effect.'),
        cold_FEL_examples=cold,checks=dict(screened_integral_error=abs(R-alternate),round_Gaussian_identity_error=round_error,constant_ODE_vs_cos_error=cos_error),
        limits=['No fit to the measured BF or energy ratios.',
                'Single coherent density mode with fixed Gaussian shape is not the measured multi-mode SASE RMS BF.',
                'Period-averaged screened potential at resonance is an estimate; the long-wavelength ESASE impedance is not asserted valid at lambda_r.',
                'Cold coherent line-integral radiation model is not a bound for a finite-aperture stochastic field with cancellation.',
                'Macro-particle noise, physical space charge and numerical self-force are not identified by these estimates alone.'])
    save(out/'RESULT.json',result);table(out/'curves.csv',curves);table(out/'cold_FEL_models.csv',cold)
    fixed=W/'results/binomial3_fel_20261106/resonance_origin_20261108/genesis_reassessment_r1/adi_only_control_r3/source/fixed/src'
    beam=fixed/'Core/BeamSolver.cpp';text=beam.read_text()
    code_lines=text.splitlines();calls=[(j+1,x.strip()) for j,x in enumerate(code_lines) if 'shortRange(' in x]
    assert calls and all(x.startswith('//') for _,x in calls)
    paths=[source,Path(__file__),beam,fixed/'Core/EFieldSolver.cpp',W/'source/aprl/src/longitudinal_shape.h',W/'source/aprl/src/source_coupling.cpp']
    save(out/'PROVENANCE.json',dict(source_hashes={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},
        genesis_effective_optical_short_range=False,genesis_short_range_occurrences=calls,
        literature=[dict(title='Huang and Kim2007, Eqs5-6 and electron phase equations',url='https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.10.034801'),
        dict(title='Wang, Litvinenko and Webb2012, space-charge parameter and dispersion, p1',url='https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.15.120701'),
        dict(title='Geloni et al2007, explicitly long-wavelength asymptote; validity warning',url='https://arxiv.org/abs/0706.2280')]))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    colors=['#626975','#167887','#b66324'];labels=['Fixed entrance size','Natural-focusing envelope','Envelope + F3 only (estimate)']
    for key,other,color,label in zip(['constant_density_transfer','envelope_density_transfer','envelope_F3_only_density_transfer'],
            ['constant_coherent_energy_ratio','envelope_coherent_energy_ratio','envelope_F3_only_coherent_energy_ratio'],colors,labels):
        axes[0].plot(z,[x[key] for x in curves],color=color,label=label)
        axes[1].plot(z,[x[other] for x in curves],color=color,label=label)
    for ax in axes:ax.set(xlabel='Undulator distance (m)',ylim=(.9,1.005));ax.grid(alpha=.18);ax.spines[['right','top']].set_visible(False)
    axes[0].set_ylabel('Cold density-mode amplitude / initial amplitude')
    axes[1].set_ylabel('Coherent line-integral energy / no-SC energy')
    axes[0].legend(fontsize=8,loc='lower left');axes[0].set_title('Space charge gives a restoring response')
    axes[1].set_title('Integrated radiation has a smaller suppression')
    fig.suptitle('Theory estimates only: no fitted coefficients, no new particle simulation',fontsize=12)
    fig.savefig(out/'theory_transfer.png',dpi=170);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    main(parser.parse_args().output)
