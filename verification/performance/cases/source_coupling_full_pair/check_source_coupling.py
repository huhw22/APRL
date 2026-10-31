import csv,json,math,sys
from pathlib import Path
import yaml

def check_sources(directory,steps,ranks):
    directory=Path(directory);files=sorted(directory.glob('source-coupling-rank-*.csv'))
    assert [p.name for p in files]==[f'source-coupling-rank-{i:05d}.csv' for i in range(ranks)]
    maximum_force=0.;maximum_work=0.;maximum_iterations=0
    for p in files:
        assert p.with_suffix('.complete').read_text().split()==['1',str(steps)]
        count=0
        with p.open() as f:
            for row in csv.DictReader(f):
                count+=1;assert int(row['step'])==count
                force=float(row['self_force_defect']);work=float(row['relative_fixed_point_work_defect']);n=int(row['iterations'])
                assert math.isfinite(force) and 0<=force<1e-13 and math.isfinite(work) and 0<=work<1e-8 and 1<=n<=8
                maximum_force=max(maximum_force,force);maximum_work=max(maximum_work,work);maximum_iterations=max(maximum_iterations,n)
        assert count==steps
    return dict(passed=True,ranks=ranks,steps=steps,maximum_force_defect=maximum_force,maximum_relative_work_defect=maximum_work,maximum_iterations=maximum_iterations)

if __name__=='__main__':
    r=Path(sys.argv[1]);p=json.loads((r/'run-plan.json').read_text());closed=json.loads((r/'completion-check.json').read_text());assert closed['passed']
    cfg=yaml.safe_load((r/'run.yaml').read_text());assert cfg['mesh']['particle_coupling']=='conservative-source'
    assert '[source-coupling] conservative-source' in (r/'aprl.log').read_text()
    result=check_sources(cfg['mesh']['moving_window']['reference_directory'],closed['final_step'],p['ranks'])
    result['physics_accepted']=False
    (r/'source-coupling-completion.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
