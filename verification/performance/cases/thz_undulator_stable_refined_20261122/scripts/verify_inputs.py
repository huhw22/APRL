"""Check frozen inputs and the refined mesh and local inlet stability constraint before starting MPI."""
import hashlib,json,math
from pathlib import Path
import h5py,numpy as np,yaml
from scipy.fft import next_fast_len
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(1<<20),b''):h.update(x)
    return h.hexdigest()
def main():
    frozen=json.loads((ROOT/'INPUT_SHA256.json').read_text())
    for name,digest in frozen['files'].items():assert sha(ROOT/name)==digest,name
    for name,digest in frozen['external_files'].items():assert sha(Path(name))==digest,name
    plan=json.loads((ROOT/'PLAN.json').read_text())
    old=json.loads((Path(plan['old_campaign'])/'PLAN.json').read_text())
    assert plan['source_revision']=='d279179a895a999ceb75412f98d6b86146e8c9b3'
    assert plan['particle_input_sha256']==old['particle_input_sha256']==sha(ROOT/'beam.h5')
    assert plan['particles']==80000 and plan['charge_C']==old['charge_C']
    with h5py.File(ROOT/'beam.h5','r') as f:
        p=f['particles/records'][:]
        assert len(p)==80000 and len(np.unique(p['source_id']))==80000
    for row in plan['rows']:
        before=next(r for r in old['rows'] if r['divisor']==row['divisor'])
        i=row['divisor'];case=ROOT/f'f0-div-{i}'
        cfg=yaml.safe_load((case/'run.yaml').read_text());pre=yaml.safe_load((case/'preflight/run.yaml').read_text())
        for key in ('K','B0_T','boost_gamma','target_frequency_Hz'):assert row[key]==before[key],key
        nx,ny,nz=cfg['mesh']['cells'];dx,dy,dz=cfg['mesh']['cell_size'];G=cfg['mesh']['boost_gamma']
        wavelength=(G+math.sqrt(G*G-1))*299792458./8.19e12
        assert wavelength/dz>=40-1e-10 and dz<=min(dx,dy)
        assert [nx,ny]==[180,120] and dx==dy==32e-6
        assert math.isclose(dz/dx,.6,rel_tol=1e-14)
        assert dz/min(dx,dy)<math.sqrt(math.sqrt(2)-1)
        r=(dz/dx)**2;assert 0<2*r*(6-r*r)/(r+2)<2
        assert math.isclose(nx*dx,.00576,rel_tol=1e-14) and math.isclose(ny*dy,.00384,rel_tol=1e-14)
        assert row['actual_forward_points_per_cycle']==wavelength/dz
        assert cfg['boundary']['cells']==[20,20,0]
        assert math.isclose(20*dx,.00064,rel_tol=1e-14)
        mw=cfg['mesh']['moving_window'];assert mw['cpml_startup']=='osiris' and mw['open_transverse_reference']
        assert mw['lower_guard_cells']>=4 and mw['upper_guard_cells']>=8 and mw['cells']//16>=8
        first=row['initial_offset_estimate']-256;N=next_fast_len(mw['cells']+514)
        assert first>=1 and first+N<nz
        assert cfg['mesh']['particle_steps_per_undulator_period']==512
        assert cfg['mesh']['longitudinal_shape']==cfg['mesh']['transverse_x_shape']=='linear-binomial3'
        assert cfg['mesh']['self_field_feedback'] and cfg['stop']['mode']=='after-last-element'
        screen=cfg['radiation_screen'];assert screen['nx']==25 and screen['ny']==19 and screen['z_m']==1.01
        assert abs(1/(screen['time_step_s']*8.19e12)-40)<1e-10
        assert pre['mesh']['cells']==cfg['mesh']['cells'] and pre['mesh']['cell_size']==cfg['mesh']['cell_size']
        assert pre['mesh']['moving_window']['cpml_startup']=='standard' and pre['stop']==cfg['stop']
        gate=yaml.safe_load((case/'stability-gate/run.yaml').read_text())
        assert gate['mesh']['self_field_feedback'] and gate['mesh']['cell_size']==cfg['mesh']['cell_size']
        assert gate['mesh']['moving_window']['inlet_charge_match'] and gate['stop']==cfg['stop']
        assert gate['mesh']['moving_window']['cpml_startup_directory']==str(case/'startup/arrays')
        assert 'radiation_screen' not in gate and not gate['detectors']['enabled']
    print('INPUTS VERIFIED: unchanged 80000 particles/K; >=40 forward samples at 8.19 THz; dz/H=0.6; 32 um transverse mesh; fresh per-grid OSIRIS startup and 256-step feedback gate',flush=True)
if __name__=='__main__':main()
