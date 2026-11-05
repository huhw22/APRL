"""Pure-coordinate regression for terminal GENESIS field-dump positions."""
from pathlib import Path
import json
import sys

import numpy as np

helper = Path(__file__).resolve().parents[3]/'analysis/sparc_sase'
sys.path.insert(0,str(helper))
from analyze_explicit import lattice_dump_positions


def run():
    dz=np.full(309,.007)
    dz[0]=.01
    z=np.r_[0.,np.cumsum(dz[:-1])]
    original_steps=np.array([0,1,77,154,308])
    original=z[original_steps]
    interior=lattice_dump_positions(z,dz,original_steps)
    assert np.array_equal(interior,original)
    terminal=lattice_dump_positions(z,dz,np.array([309]))[0]
    assert terminal==z[-1]+dz[-1]
    assert abs(terminal-2.166)<1e-12
    # The common old full lattice already contains the first-segment exit.
    # Its existing coordinate and the short-run terminal boundary agree.
    longer_z=np.r_[z,terminal,terminal+.007]
    longer_dz=np.r_[dz,.007,.007]
    assert lattice_dump_positions(longer_z,longer_dz,np.array([309]))[0]==terminal
    for bad in [np.array([-1]),np.array([310]),np.array([1.5])]:
        try:lattice_dump_positions(z,dz,bad)
        except RuntimeError:pass
        else:raise AssertionError('Invalid dump coordinate accepted')
    return dict(passed=True,synthetic_lattice_steps=309,
        terminal_step=309,terminal_z_m=float(terminal),
        old_interior_coordinates_bitwise_preserved=True,
        short_terminal_matches_long_interior=True,
        invalid_negative_beyond_terminal_fractional_steps_rejected=True,
        scope='Coordinate lookup only; no HDF5 read, simulation or physical acceptance')


if __name__=='__main__':
    result=run()
    if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
