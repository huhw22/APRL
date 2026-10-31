import json,os,sys
from pathlib import Path
from prepare import sha
r=Path(sys.argv[1]).resolve();p=json.loads((r/'run-plan.json').read_text())
for name,expected in json.loads((r/'frozen.json').read_text()).items():assert sha(name)==expected,name
pre=json.loads(Path(p['preflight_pass']).read_text());assert pre['pass'] and pre['binary_sha256']==p['binary_sha256']
assert pre['card_sha256'][str(r/'run.yaml')]==sha(r/'run.yaml')
st=Path(p['input_file']).stat();assert (st.st_size,st.st_mtime_ns)==(p['input_size_bytes'],p['input_mtime_ns'])
assert (int(os.environ['SLURM_JOB_NUM_NODES']),int(os.environ['SLURM_NTASKS']),int(os.environ['SLURM_CPUS_PER_TASK']))==(p['nodes'],p['ranks'],8)
assert not any(os.environ.get(k) for k in ['APRL_DIAGNOSTIC_STEPS','APRL_FORCE_SPLIT_DIRECTORY','APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_WORK_GEOMETRY'])
print(json.dumps(dict(passed=True,binary_sha256=p['binary_sha256'],input_sha256=p['input_sha256'],input_stat_unchanged=True,job_id=os.environ['SLURM_JOB_ID'],nodes=p['nodes'],ranks=p['ranks'],threads=8,diagnostic_environment_disabled=True),indent=2))
