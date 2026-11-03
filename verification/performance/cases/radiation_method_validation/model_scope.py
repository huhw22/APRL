"""A priori cold rms-envelope estimate; never fitted to PIC or GENESIS."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import c,e,m_e,epsilon_0
from scipy.integrate import solve_ivp
from prepare import save

def main():
    ws=Path('/public/home/wfollowm/APRL_Workspace')
    base=ws/'results/radiation_method_validation_20261103';root=base/'fel-r4'
    plan=json.loads((root/'campaign.json').read_text());p=plan['parameters']
    IA=4*np.pi*epsilon_0*m_e*c**3/e
    gamma=p['gamma'];I=p['peak_current_A'];sigma=p['sigma_x_m']
    ky=p['undulator_K']**2*(2*np.pi/p['period_m'])**2/(2*gamma**2)
    L=p['periods']*p['period_m']
    def rhs(z,q):
        sx,sy,ax,ay=q
        common=I/(IA*gamma**3*(sx+sy))
        return [ax,ay,common,common-ky*sy if z<L else common]
    initial=[sigma,sigma,0.,0.]
    sol=solve_ivp(rhs,[0,L+.3],initial,method='DOP853',rtol=1e-11,atol=1e-15,dense_output=True,max_step=.005)
    result=dict(input_current_A=I,input_gamma=gamma,
       initial_round_estimate_fractional_expansion=I*L*L/(4*IA*gamma**3*sigma*sigma),
       assumptions=['Long beam in mean rest frame; cold rms closure',
                    'Neglect generated emittance and slice/end variations',
                    'Period-averaged natural undulator focusing; no fitted PIC parameters'],
       at_undulator_exit=dict(sigma_x_m=float(sol.sol(L)[0]),sigma_y_m=float(sol.sol(L)[1])),
       at_far_plane=dict(sigma_x_m=float(sol.sol(L+.3)[0]),sigma_y_m=float(sol.sol(L+.3)[1])))
    save(base/'TRANSVERSE_SC_MODEL_SCOPE.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
