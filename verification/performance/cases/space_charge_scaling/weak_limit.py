#!/usr/bin/env python3
"""Independent leading weak-charge expansion at fixed laboratory planes."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.constants import epsilon_0,m_e,e,c

def main(root):
    root=Path(root).resolve();sys.path.insert(0,str(root/'scripts'))
    from model import Sphere,G_BEAM,B_BEAM
    model=Sphere('uniform');R=np.sqrt(5)*model.sigma;T0=model.entry_center_time
    theory=json.loads((root/'THEORY.json').read_text());z=theory[0]['last_plane_z_m']
    T1=T0+z/(G_BEAM*B_BEAM*c);rows=[]
    for item in theory:
        omega2=item['electrons']*e**2/(4*np.pi*epsilon_0*m_e*R**3)
        # Keep the fixed-plane event-time change through second order as well.
        # Rest kinetic energy, event time, and changing force give 3/10-1/10+1/5.
        predicted=2*G_BEAM*m_e*omega2**2*R**2*(T1*T1-T0*T0)/(5*e)
        rows.append(dict(charge_factor=item['charge_factor'],full_sphere_gain_eV=item['fixed_plane_lab_mean_gain_eV'],
                         cold_weak_limit_gain_eV=predicted,
                         weak_limit_minus_full_pct=100*(predicted/item['fixed_plane_lab_mean_gain_eV']-1)))
    result=dict(rows=rows,T0_s=T0,T1_s=T1,
                formula='2*Gamma*m*omega^4*R^2*(T1^2-T0^2)/5; omega^2=N e^2/(4 pi eps0 m R^3).',
                derivation='Keep the first-order plane-crossing time shift in the second-order mean; coefficient 3/10-1/10+1/5=2/5. Replaces the former ultrarelativistic (3+beta^2)/10 auxiliary estimate.',
                scope='Cold uniform sphere, fixed lab planes, leading weak expansion. Not exact at original strong charge.')
    (root/'WEAK_LIMIT_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(sys.argv[1])
