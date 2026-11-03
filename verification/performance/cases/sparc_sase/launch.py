import os,signal,subprocess,time
from pathlib import Path
from common import save,need,EXEC_NODES
def run(cmd,directory,limit,log,env=None,expected_duration_steps=None):
    path=Path(directory);env=os.environ.copy() if env is None else env
    rec=dict(command=cmd,cwd=str(path),limit_seconds=limit,start_unix=time.time(),
             environment={k:env.get(k) for k in ['OMP_NUM_THREADS','OMPI_MCA_btl','HDF5_USE_FILE_LOCKING']})
    with (path/log).open('x') as stream:
        start=time.monotonic()
        p=subprocess.Popen(['timeout','--signal=TERM','--kill-after=20',str(limit),*cmd],
           cwd=path,env=env,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:code=p.wait(timeout=limit+30)
        except BaseException:
            os.killpg(p.pid,signal.SIGTERM)
            try:p.wait(timeout=20)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
            raise
    rec.update(exit_code=code,seconds=time.monotonic()-start)
    save(path/(log+'.launch.json'),rec)
    if expected_duration_steps is None:
        need(code==0,'Nonzero bounded process: '+str(path/log))
    else:
        import duration_status
        status=duration_status.bounded_duration_status((path/log).read_text(),
          dict(exit_code=code,seconds=rec['seconds'],limit_s=limit),
          expected_duration_steps,duration_status)
        save(path/'BOUNDED_DURATION.json',status)
    return rec
def srun(nodes=None,ranks=None):
    nodes=EXEC_NODES if nodes is None else nodes
    ranks=16*nodes if ranks is None else ranks
    return ['srun','--mpi='+os.environ.get('SRUN_MPI','pmix_v3'),'--exact',
      '--nodes='+str(nodes),'--ntasks='+str(ranks),'--ntasks-per-node=16','--cpus-per-task=8',
      '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
