"""Replay the original three space-charge physical references without retuning."""
import argparse,json,shutil
from pathlib import Path
import yaml
from prepare import sha,save,require
def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--build-record',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();built=a.build_record.resolve();r=a.output.resolve();r.mkdir(parents=True,exist_ok=False);proof=json.loads((built/'RESULT.json').read_text());require(len(proof['units'])==3 and all(x['weighted_work_relative_max']<5e-12 and x['weighted_magnetic_impulse_relative_max']<5e-12 for x in proof['units']),'MPI interpolation gate');require(all(x['complete'] for x in proof['periodic']),'Closed exchange execution gate');binary=Path(json.loads((built/'PLAN.json').read_text())['programs']['aprl']);require(sha(binary)==json.loads((built/'SOURCE_BUILD.json').read_text())['program_sha256']['aprl'],'Build identity')
    baseline=w/'results/source_coupling_exit_20261031/physics-r1';scripts=r/'scripts';shutil.copytree(baseline/'scripts',scripts);here=Path(__file__).resolve().parent;shutil.copy2(here/'physics.sbatch',scripts/'physics.sbatch');cases={}
    for name in ['uniform','converging','gaussian']:
        old=baseline/name;d=r/name;d.mkdir();(d/'momentum').mkdir();cfg=yaml.safe_load((old/'run.yaml').read_text());beam=cfg['beam']
        def change(v):
            if isinstance(v,dict):return {k:change(x) for k,x in v.items()}
            if isinstance(v,list):return [change(x) for x in v]
            return v.replace(str(old),str(d)) if isinstance(v,str) else v
        cfg=change(cfg);cfg['beam']=beam;(d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));shutil.copy2(old/'case.json',d/'case.json');(d/'reference_input').symlink_to((old/'reference_input').resolve(),target_is_directory=True)
        cases[name]=dict(directory=str(d),original=str(old),steps=json.loads((d/'case.json').read_text())['geometry']['expected_stop_step'],timeout_s=2700,first_source_reference=str(old/'FIRST_SOURCE_STATES.npy'))
    plan=dict(workspace=str(w),build=str(binary.parent),binary=str(binary),binary_sha256=sha(binary),binding_probe=str(w/'build/cluster-tools/binding_probe'),preflight=str(built),cases=cases,jobs={n:[n] for n in cases},nominal_dt_s=1.2543262415359034e-6/299792458.,selected_cases=list(cases),scope='Matched self-field gather candidate: original charge/input/grid/initial reference/device-free stop/theory unchanged; exploratory physical comparison, not acceptance',periodic_energy_and_magnetic_exchange_verified=True,all_periodic_momentum_cases_accepted=False,production_replaced=False);save(r/'plan.json',plan)
    files=[r/'plan.json',binary,*[p for p in scripts.iterdir() if p.is_file()]]
    for case in cases.values():
        d=Path(case['directory']);files.extend([d/'run.yaml',d/'case.json',d/'reference_input/expected_events.npz',Path(yaml.safe_load((d/'run.yaml').read_text())['beam']['input']['file'])])
    save(r/'FROZEN.json',{str(p):sha(p) for p in files});print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
