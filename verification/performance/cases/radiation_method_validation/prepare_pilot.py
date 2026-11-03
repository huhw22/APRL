"""Freeze existing closed dipole evidence and inexpensive GENESIS design pilots."""
import argparse
import math
import shutil
from pathlib import Path
import yaml
from prepare import require,save,sha

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('workspace',type=Path)
    ap.add_argument('output',type=Path)
    a=ap.parse_args()
    w=a.workspace.resolve();root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    scripts=root/'scripts';shutil.copytree(Path(__file__).parent,scripts)
    old=w/'results/dipole_gaussian_radiation_20261102/r3'
    require((old/'PIC_COMPLETE.json').exists(),'Old PIC did not close')
    shutil.copy2(old/'PLAN.json',root/'PLAN.json')
    g=1174.171072271572;ku=1.092714657005;period=.0234
    wavelength=period*(1+ku*ku/2)/(2*g*g)
    configs=[]
    for sigma in [20e-6,50e-6]:
        for current in [10.,100.,1000.]:
            name=f'sigma{sigma*1e6:g}-I{current:g}'
            d=root/'genesis-design'/name;d.mkdir(parents=True)
            # Small nonzero emittance only defines a Gaussian size in GENESIS.
            emit=1e-10;beta=sigma*sigma*g/emit
            deck=f"""&setup
rootname=design
lattice=design.lat
beamline=FEL
lambda0={wavelength:.17g}
gamma0={g:.17g}
delz={period/4:.17g}
shotnoise=false
npart=2048
nbins=16
&end
&field
power=0
dgrid=0.00075
ngrid=101
waist_size={2*sigma:.17g}
&end
&beam
current={current}
delgam=0
ex={emit}
ey={emit}
betax={beta}
betay={beta}
bunch=0.01
&end
&track
output_step=4
&end
"""
            (d/'design.in').write_text(deck)
            (d/'design.lat').write_text(
                f'UND: UNDULATOR = {{ lambdau={period}, nwig=400, aw={ku/math.sqrt(2):.17g}, helical=false }};\nFEL: LINE={{UND}};\n')
            configs.append({'directory':str(d),'name':name,'sigma_m':sigma,'current_A':current})
    program=Path('/public/software/Genesis-1.3-Version4-4.6.4/build/genesis4')
    save(root/'pilot-plan.json',{'workspace':str(w),'old_closed_dipole':str(old),
      'genesis_binary':str(program),'genesis_sha256':sha(program),
      'genesis_design':configs,'resonant_wavelength_m':wavelength,
      'scope':'Existing dipole absolute-reference repair and FEL parameter planning; not three-method acceptance'})
    paths=[root/'pilot-plan.json',root/'PLAN.json',program,*scripts.glob('*.py')]
    paths.extend(d/p for d in [Path(q['directory']) for q in configs] for p in ['design.in','design.lat'])
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths})
    print(root)

if __name__=='__main__':main()
