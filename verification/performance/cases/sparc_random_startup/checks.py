"""Small offline sampling/mapping checks, no PIC or GENESIS tracking."""
import copy
import tempfile
from pathlib import Path

import h5py
import numpy as np

from common import frozen, load
from make_inputs import generate
from analyze import moments


def main(root):
    original=load(Path(root)/'PLAN.json')
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp);(tmp/'inputs').mkdir()
        plan=copy.deepcopy(original);plan['source_config']['count_mode']='diagnostic_fixed_n'
        plan['source_config']['electrons']=65536
        for entry,weight in zip(plan['cases'],(8,8,2)):
            entry['directory']=str(tmp/entry['name']);Path(entry['directory']).mkdir()
            for name in ('genesis-off','genesis-longrange'):(Path(entry['directory'])/name).mkdir()
            entry['target_weight']=weight
        generate(tmp,plan)
        arrays={}
        for entry in plan['cases']:
            d=Path(entry['directory']);info=load(d/'INPUT_VALIDATED.json');frozen(d/'FROZEN.json')
            assert info['zero_transverse_duplicates']
            with h5py.File(info['input'],'r') as f:a=f['particles/records'][:]
            arrays[entry['name']]=a
            with h5py.File(d/'common.par.h5','r') as f:
                ids=[];ne=0
                for index in load(d/'NATIVE_MAPPING.json')['active_groups']:
                    g=f[f'slice{index:06d}'];s=g['source_id'][:]
                    rows=np.searchsorted(a['source_id'],s)
                    assert np.array_equal(a['source_id'][rows],s)
                    assert np.array_equal(a['plane_position_m'][rows,0],g['x'][:])
                    assert np.array_equal(a['plane_position_m'][rows,1],g['y'][:])
                    assert np.array_equal(a['proper_velocity'][rows,:2],np.c_[g['px'][:],g['py'][:]])
                    ids.extend(s.tolist());ne+=g['current'][0]*plan['wavelength_m']/299792458/1.602176634e-19
                assert len(set(ids))==len(a) and abs(ne/info['physical_electrons']-1)<1e-12
        small=arrays['w512-a'];large=arrays['w64-a'];rows=np.searchsorted(large['source_id'],small['source_id'])
        for key in ('source_id','arrival_time_offset_s','plane_position_m','proper_velocity'):
            assert np.array_equal(small[key],large[key][rows])
        assert not np.array_equal(small['plane_position_m'],arrays['w512-b']['plane_position_m'])
        x=moments(np.arange(1,5,dtype='u8'),np.ones(4),np.array([1,1,-1,-1],complex),np.array([0,2,4]),0,2)
        assert x['local_rms_b1']==1 and x['global_b1']==0 and x['common_phase_bf_fraction']==0
    print('PASS: independent sampling, nested IDs/6D, separate seeds, weighted current, native mapping, BF projection',flush=True)


if __name__=='__main__':
    import sys
    main(sys.argv[1])
