import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.constants import c,e,m_e

WS=Path('/public/home/wfollowm/APRL_Workspace')
GAMMA=1+151.9e6*e/(m_e*c*c)
LU=.028
K=e*.78796*LU/(2*np.pi*m_e*c)
LAMBDA=LU*(1+K*K/2)/(2*GAMMA*GAMMA)
T=12.67e-12
Q=450e-12
PREP=.01
IP=1.5*Q/T
PUBLISHED_T=T
PUBLISHED_Q=Q
PUBLISHED_IP=IP
DX=DY=20e-6
NX=NY=128
EXEC_NODES=2
TOTAL_SLIP=6*77*LAMBDA
NBINS=16
NPART=1024
GBOOST=GAMMA/np.sqrt(1+K*K/2)
GENESIS=Path('/public/software/Genesis-1.3-Version4-4.6.4/build/genesis4')
def save(p,v):
    Path(p).write_text(json.dumps(v,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(x,msg):
    if not x:raise RuntimeError(msg)
def native(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
    if isinstance(x,list):return [native(v) for v in x]
    return x

def configure(settings):
    global T,Q,IP,DX,DY,NX,NY,EXEC_NODES
    multiple=settings.get('slippage_multiple')
    IP=PUBLISHED_IP
    T=PUBLISHED_T if multiple is None else float(multiple)*TOTAL_SLIP/c
    Q=2*IP*T/3
    DX=DY=float(settings.get('transverse_cell_um',20))*1e-6
    NX=NY=int(settings.get('transverse_cells',128))
    EXEC_NODES=int(settings.get('nodes',2))
    need(T>0 and DX>0 and NX>20 and EXEC_NODES in [1,2],'Invalid benchmark settings')

def stage_scripts(root,settings):
    import re,shutil
    root=Path(root)
    save(root/'benchmark.json',settings)
    shutil.copytree(Path(__file__).parent,root/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
    nodes=int(settings.get('nodes',2))
    for name in ['preflight.sbatch','main.sbatch']:
        p=root/'scripts'/name;text=p.read_text()
        text=re.sub(r'#SBATCH --nodes=\d+',f'#SBATCH --nodes={nodes}',text)
        text=re.sub(r'#SBATCH --ntasks=\d+',f'#SBATCH --ntasks={16*nodes}',text)
        p.write_text(text)

_settings_path=Path(__file__).resolve().parents[1]/'benchmark.json'
if _settings_path.exists():configure(load(_settings_path))

def slice_profile():
    ds=LAMBDA
    n=int(np.ceil(c*T/ds-1e-10))
    a=np.arange(n)*ds;b=np.minimum(a+ds,c*T)
    # Integral of I_p [1-(2s/L-1)^2], preserving Q exactly.
    def primitive(s):return IP*(2*s*s/(c*T)-4*s**3/(3*(c*T)**2))
    current=(primitive(b)-primitive(a))/ds
    current=np.maximum(current,0.)
    need(abs(current.sum()*ds/c/Q-1)<1e-11,'Charge/profile mismatch')
    return current

def main_lattice(directory):
    """One-period internal ramps; independently precomputed phase matching."""
    from scipy.integrate import quad
    d=Path(directory);d.mkdir(parents=True,exist_ok=True)
    period=LU;length=77*LU;ku=2*np.pi/LU;ks=2*np.pi/LAMBDA
    beta=np.sqrt(1-1/GAMMA**2)
    gap_lag=.4*(1/beta-1)
    def envelope(z):
        a=min(z,length-z)
        return .5*(1-np.cos(np.pi*a/LU)) if a<LU else 1.
    def delay(phase):
        def f(z):
            ux=-K*envelope(z)*np.cos(ku*z+phase)
            uz=np.sqrt(GAMMA**2-1-ux*ux)
            return (1+ux*ux)/(uz*(GAMMA+uz))
        return sum(quad(f,j*LU,(j+1)*LU,epsabs=1e-17,epsrel=1e-11)[0] for j in range(77))
    phases=[0.]
    kicks=[];closure=[]
    aw=K/np.sqrt(2)
    ramp_aw=[aw*.5*(1-np.cos(np.pi*(j+.5)/4)) for j in range(4)]
    def averaged_lag(a):
        u=np.sqrt(GAMMA*GAMMA-1-a*a)
        return (1+a*a)/(u*(GAMMA+u))
    g4_lag=75*LU*averaged_lag(aw)+2*sum(LU/4*averaged_lag(a) for a in ramp_aw)
    for i in range(5):
        lag=delay(phases[-1])
        delta=ks*(lag+gap_lag)-2*np.pi*77
        # Local G4 4.6.4 uses autophase=-phi once per non-undulator
        # interval. Its fractional free-drift auto-phasing is commented out.
        kick=-delta
        phases.append(float(np.angle(np.exp(1j*(phases[-1]+delta)))))
        kicks.append(float(np.angle(np.exp(1j*kick))))
        closure.append(float(ks*(lag-g4_lag)))
    magnets=[];defs=[];line=[]
    aw=K/np.sqrt(2)
    # G4 uses a period-averaged ramp sampled at quarter-period midpoints.
    ramp=[]
    for j in range(4):
        value=.5*(1-np.cos(np.pi*(j+.5)/4))
        name=f'R{j}'
        defs.append(f'{name}: UNDULATOR={{lambdau={LU:.17g},nwig=0.25,aw={aw*value:.17g}}};')
        ramp.append(name)
    defs += [f'PREP: DRIFT={{l={PREP}}};',f'U: UNDULATOR={{lambdau={LU},nwig=75,aw={aw:.17g}}};',
      'D1: DRIFT={l=.105};',
      f'Q: QUADRUPOLE={{l=.053,k1={e*9/(m_e*c*np.sqrt(GAMMA**2-1)):.17g}}};',
      'FAR: DRIFT={l=.3};']
    line.append('PREP')
    for i in range(6):
        z=PREP+i*(length+.4)
        magnets.append(dict(type='planar-undulator',entrance_z=z,period=LU,periods=77,
          strength_parameter=K,polarization_angle_rad=float(np.pi/2),gaussian_fringe=False,
          end_ramp_periods=1.,phase_rad=phases[i]))
        line.extend([*ramp,'U',*ramp[::-1]])
        if i<5:
            magnets.append(dict(type='quadrupole',entrance_z=z+length+.105,length=.053,
                                gradient_T_per_m=-9.))
            # Replace the remaining 0.242 m drift by one finite interval.
            # Zero-length phase shifters are not found by G4 findElement.
            defs.append(f'P{i}: PHASESHIFTER={{l=.242,phi={kicks[i]:.17g}}};')
            line.extend(['D1','Q',f'P{i}'])
    line.append('FAR')
    defs.append('FEL: LINE={'+','.join(line)+'};')
    (d/'reference.lat').write_text('\n'.join(defs)+'\n')
    save(d/'lattice.json',dict(magnets=magnets,segment_phases_rad=phases,
      genesis_gap_phase_kicks_rad=kicks,physical_line_length_m=6*length+5*.4,
      genesis_phase_implementation='Finite final 0.242 m gap interval; one autophase=-phi kick, no extra drift length; local G4 fractional auto-phasing is disabled',
      nominal_phase_closure_G4_minus_exact_rad=closure,
      field_far_z_m=PREP+6*length+5*.4+.3,preparation_drift_m=PREP,
      quadrupole_location_convention='0.105 m leading edge downstream of preceding undulator',
      phase_matching='Nominal magnetic trajectory and free-drift flight delay, no radiation fitting',
      model_scope='Published SPARC beam parameters with explicit reconstructed ends, quadrupole placement and phase-matched lattice; not an exact original input deck'))
    return magnets
