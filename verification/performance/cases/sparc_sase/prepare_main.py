"""Derive immutable long cases from a completed, bounded full-geometry preflight."""
import argparse,copy,shutil
from pathlib import Path
import yaml
from common import *
def main():
    ap=argparse.ArgumentParser();ap.add_argument('preflight',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
    pre=a.preflight.resolve();root=a.output.resolve()
    marker=pre/'PREFLIGHT_VALIDATED.json'
    if not marker.exists():marker=pre/'PREFLIGHT_COMPLETE.json'
    need(load(marker)['pass_preflight'],'Preflight has not passed')
    need(all(sha(p)==h for p,h in load(pre/'FROZEN_CARDS.json').items()),'Preflight cards or input changed')
    root.mkdir(exist_ok=False)
    settings=load(pre/'benchmark.json') if (pre/'benchmark.json').exists() else {}
    stage_scripts(root,settings)
    plan=copy.deepcopy(load(pre/'plan.json'));plan['preflight']=str(pre);plan['root']=str(root)
    cases=[];far=plan['geometry']['field_planes'][-1][1]
    for case in plan['cases']:
        old=Path(case['directory']);d=root/case['name'];d.mkdir()
        cfg=yaml.safe_load((old/'run.yaml').read_text())
        def replace(x):
            if isinstance(x,dict):return {k:replace(v) for k,v in x.items()}
            if isinstance(x,list):return [replace(v) for v in x]
            return x.replace(str(old),str(d)) if isinstance(x,str) else x
        cfg=replace(cfg)
        cfg['detectors']['field_planes'].insert(-1,dict(name='free-015',z=far-.15,
          rhythm=plan['geometry']['dt_s']/GBOOST,buffer_samples=8,compression=1,
          particle_background={'enabled':False}))
        cfg['detectors']['particle_planes'].insert(-1,dict(name='free-015-particles',z=far-.15,
          buffer_records=32768,compression=0))
        cfg['detectors']['particle_planes'].insert(1,dict(name='und-entry-particles',z=PREP,
          buffer_records=32768,compression=0))
        cfg['detectors']['energy_planes']=[dict(name=n,z=z,
          rhythm=plan['geometry']['dt_s']/GBOOST,buffer_samples=64,compression=0)
          for n,z in [('entry-transport',0.),('und-entry-transport',PREP),('exit-transport',far)]]
        (d/'run.yaml').write_text(yaml.safe_dump(native(cfg),sort_keys=False))
        cases.append(dict(name=case['name'],directory=str(d),timeout_seconds=136*3600))
    plan['cases']=cases;plan['launch_stage']='One physical noise seed, full six-segment SASE comparison'
    nodes=int(settings.get('nodes',2))
    plan['execution']=dict(nodes=nodes,ranks=16*nodes,threads_per_rank=8,case_order=['boris','energy','momentum'],
      array_max_concurrent=1,wall_limit_hours_per_case=144,
      reason='User-approved reduced bunch and transverse grid; fixed layout across methods, no fixed node selection.')
    plan['geometry']['field_planes'].insert(-1,['free-015',far-.15])
    plan['geometry']['uncompressed_EB_upper_GiB_per_case']=8*(plan['geometry']['steps']+1)*plan['geometry']['nx']*plan['geometry']['ny']*48/2**30
    gen=root/'genesis';gen.mkdir()
    # Only zero-length diagnostic markers are added to the frozen reference.
    magnets=main_lattice(gen)
    need(magnets==cfg['sources']['magnetic_elements'],'Reconstructed physical lattice changed')
    text=(gen/'reference.lat').read_text()
    text=text.replace('FEL: LINE={','M: MARKER={dumpfield=1};\nFEL: LINE={')
    text=text.replace('R3,R2,R1,R0,','R3,R2,R1,R0,M,')
    text=text.replace('R3,R2,R1,R0,FAR','R3,R2,R1,R0,M,FAR')
    text=text.replace('FAR: DRIFT={l=.3};','F15: DRIFT={l=.15};\nFAR: DRIFT={l=.15};')
    text=text.replace(',FAR};',',F15,M,FAR};')
    (gen/'reference.lat').write_text(text)
    (gen/'common.par.h5').symlink_to(pre/'genesis/common.par.h5')
    (gen/'reference.in').write_text(f'''&setup
rootname=reference
lattice=reference.lat
beamline=FEL
lambda0={LAMBDA:.17g}
gamma0={GAMMA:.17g}
delz={LU/4:.17g}
shotnoise=false
npart={NPART}
nbins={NBINS}
beam_global_stat=true
field_global_stat=true
&end
&importbeam
file=common.par.h5
time=true
&end
&field
power=0
dgrid={plan['geometry'].get('physical_half_width_m',.00108):.17g}
ngrid=161
&end
&efield
longrange=true
rmax=2
nz=1
nphi=0
ngrid=100
&end
&track
output_step=16
&end
&write
field=far
&end
''')
    save(root/'plan.json',plan)
    files=[plan['binary'],plan['genesis'],root/'plan.json',root/'benchmark.json',gen/'reference.in',gen/'reference.lat',gen/'lattice.json',
      gen/'common.par.h5',pre/'particles.h5',*list((root/'scripts').iterdir()),
      *[Path(q['directory'])/'run.yaml' for q in cases]]
    save(root/'FROZEN.json',{str(Path(p)):sha(p) for p in files if Path(p).is_file()})
    print(root)
if __name__=='__main__':main()
