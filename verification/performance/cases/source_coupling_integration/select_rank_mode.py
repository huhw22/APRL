import os,sys
rank=int(os.environ['SLURM_PROCID'])
os.execv(sys.argv[1],[sys.argv[1],sys.argv[2] if rank==0 else sys.argv[3]])
