#!/usr/bin/env python3
"""Independently audit a completed 256-rank two-node checkpoint campaign.

Reads existing HDF5 files in <=8 MiB chunks. Never launches computations or
modifies source evidence. Creates only a new JSON review in the result root.
This checks schema v1 on the recorded x86_64 ABI: 16-byte x87 long doubles
contain 80 value bits and 48 unspecified padding bits. The exact 19 stored
long doubles have only those padding bytes masked; all effective bits and
all other scalar bytes must match. Changing the native schema must update
and revalidate this decoder, rather than weakening its size assertions.

Example: python3 review.py --root /absolute/results/parallel-final
To repeat without replacing a prior review, pass --output-name REVIEW_2.json.
"""
import argparse,datetime,hashlib,json,re,struct
from pathlib import Path
import h5py,numpy as np,yaml
CHUNK=8<<20

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(CHUNK),b''):h.update(b)
 return h.hexdigest()
def require(x,m):
 if not x:raise RuntimeError(m)
def finite(a,name):
 if a.dtype.names:
  for k in a.dtype.names:finite(a[k],name+'.'+k)
 elif a.dtype.kind in 'fc':require(np.isfinite(a).all(),'nonfinite '+name)
def bits(a,b,name):
 require(a.shape==b.shape and a.dtype==b.dtype,'shape/type '+name)
 if a.dtype.names:
  for k in a.dtype.names:bits(a[k],b[k],name+'.'+k)
 else:require(a.tobytes()==b.tobytes(),'bitwise '+name)
def rank_blocks(f,name):
 g=f['blocks/'+name];counts=[int(v) for v in g['rank_counts'][:]];data=g['data'][:];pos=0;out=[]
 for n in counts:out.append(data[pos:pos+n].tobytes());pos+=n
 require(pos==len(data),'rank blocks '+name);return out

def inspect_snapshot(path):
 result={'path':str(path),'file_bytes':path.stat().st_size,'datasets':{}}
 with h5py.File(path) as f:
  header=f['header'][:];require(header.tolist()==[0x4150524c434b5031,1,1,256,8,8,16,64,16384,1,1,1],'header '+str(path))
  names=[];f.visititems(lambda n,x:names.append(n) if isinstance(x,h5py.Dataset) else None)
  for name in names:
   ds=f[name];h=hashlib.sha256();nonzero=0
   if ds.shape:
    width=max(1,ds.dtype.itemsize*int(np.prod(ds.shape[1:],dtype=np.int64)))
    chunk=max(1,CHUNK//width)
    for start in range(0,ds.shape[0],chunk):
     a=ds[start:start+chunk];finite(a,name);h.update(a.tobytes())
     if a.dtype.kind in 'fiu':nonzero+=int(np.count_nonzero(a))
   else:
    a=np.asarray(ds[()]);finite(a,name);h.update(a.tobytes())
   result['datasets'][name]={'shape':list(ds.shape),'dtype':repr(ds.dtype.descr if ds.dtype.names else ds.dtype.str),'bytes':int(ds.size*ds.dtype.itemsize),'sha256':h.hexdigest(),'nonzero_elements':nonzero}
  payload=0
  for name in names:
   if not name.startswith('blocks/') or not name.endswith('/data'):continue
   g=f[name.rsplit('/',1)[0]];cs=g['rank_counts'][:]
   require(cs.shape==(256,),'256 rank counts '+name)
   require(sum(map(int,cs))==len(g['data']),'counts total '+name)
   require(int(g['kind'][0]) in (1,2,3,4),'kind '+name)
   payload+=g['data'].size*g['data'].dtype.itemsize
  result['payload_bytes']=int(payload)
  state=rank_blocks(f,'state/scalars');geo=rank_blocks(f,'state/geometry')
  require(all(len(b)==1114 for b in state),'scalar schema size')
  require(all(len(b)==337 for b in geo),'geometry schema size')
  clocks=[struct.unpack_from('=dQ',b,4) for b in state]
  require(all(c==clocks[0] for c in clocks),'rank clocks differ')
  result['time_box_s'],result['step']=clocks[0]
  require(all(struct.unpack_from('=I',b)[0]==1 for b in state+geo),'state schemas')
  # Native x87 long double occupies 16 bytes: only the initial 10 are value
  # bits. There are 19 consecutive long doubles from scalar offset 32.
  mask=np.ones(1114,dtype=bool)
  for i in range(19):mask[32+i*16+10:32+(i+1)*16]=False
  meaningful=hashlib.sha256()
  for b in state:meaningful.update(np.frombuffer(b,dtype='u1')[mask].tobytes())
  result['scalar_meaningful_sha256']=meaningful.hexdigest()
  result['scalar_padding_excluded_bytes_per_rank']=int((~mask).sum())
  result['particle_count']=len(f['blocks/particles/data'])
  result['frozen_particles']=int(np.count_nonzero(f['blocks/particles/data']['frozen']))
  result['carrier_counts']={n:len(f['blocks/carriers/'+n+'/data']) for n in ['cpml','retirement']}
  prefixes=[yaml.safe_load(b.decode()) for b in rank_blocks(f,'metadata/diagnostic_prefixes')]
  require(all(v['step']==result['step'] and v['time_box_s']==result['time_box_s'] for v in prefixes),'prefix rank clock')
  result['prefixes']=prefixes[0]['outputs']
 return result

def review_directory(P, output_name):
 require(P.is_dir(), 'Result directory is missing')
 require(Path(output_name).name == output_name and output_name.endswith('.json'), 'Output must be a JSON basename')
 require(not (P/output_name).exists(), 'Review exists; choose another --output-name')
 start=datetime.datetime.fromisoformat((P/'JOB_STARTED_AT').read_text().strip())
 end=datetime.datetime.fromisoformat((P/'JOB_ENDED_AT').read_text().strip())
 allocation=(P/'allocation.txt').read_text()
 job_match=re.search(r'JobId=(\d+)',allocation)
 require(job_match is not None, 'Allocation job identity is missing')
 job_evidence={'job':int(job_match.group(1)), 'batch_exit_code':0,
               'start':start.isoformat(),'end':end.isoformat(),
               'elapsed_seconds':(end-start).total_seconds(),
               'source':'Local JOB_STARTED_AT/JOB_ENDED_AT/JOB_EXIT_CODE, allocation and binding evidence; no scheduler query'}
 require(job_evidence['elapsed_seconds']>0,'Job timestamp order')
 plan=json.loads((P/'plan.json').read_text());frozen=json.loads((P/'frozen.json').read_text())
 for path,value in frozen.items():require(sha(path)==value,'frozen identity '+path)
 require(sha(plan['binary'])==plan['binary_sha256'],'parallel binary changed')
 require((P/'JOB_EXIT_CODE').read_text().strip()=='0','batch exit')
 require((P/'PASS').is_file(),'parent PASS missing')
 bind=(P/'binding.log').read_text();rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)',bind)
 require(len(rows)==256 and {int(r[0]) for r in rows}==set(range(256)),'binding rank coverage')
 require(all(r[1]=='0' and r[7]=='PASS' for r in rows),'binding status')
 hosts=sorted({r[2] for r in rows})
 require(len(hosts)==2 and all(sum(r[2]==h for r in rows)==128 for h in hosts),'binding nodes')
 require(len({(r[2],r[4],r[5]) for r in rows})==256,'unique physical core binding')
 paths=[P/'continuous/checkpoints/step-000000000003.h5',P/'continuous/checkpoints/step-000000000009.h5',P/'restart/checkpoints/step-000000000009.h5']
 snaps=[]
 for path in paths:
  print('Reading snapshot',path,flush=True);snaps.append(inspect_snapshot(path))
 a,b=snaps[1:];same=[]
 for name,record in a['datasets'].items():
  if name.startswith(('blocks/fields/','blocks/cpml/','blocks/particles/','blocks/carriers/','blocks/state/geometry/','blocks/state/detectors/')) or name=='header':
   require(record==b['datasets'][name],'snapshot bitwise '+name);same.append(name)
 require(a['scalar_meaningful_sha256']==b['scalar_meaningful_sha256'],'meaningful native scalars')
 require(a['step']==b['step']==9 and snaps[0]['step']==3,'checkpoint steps')
 require(abs(snaps[0]['time_box_s']-3*.5e-6/299792458)<1e-28,'clock consistency')
 require(snaps[1]['time_box_s']==snaps[2]['time_box_s'],'restart preserved accumulated box time')
 # Exact diagnostic suffixes start at the persisted parent prefix lengths.
 identity={'run_id','run_manifest','configuration_digest_fnv1a64','source_revision','restart_parent','restart_step','restart_time_box_s','output_scope'}
 output_reviews=[]
 for relative,group,commit,data in [('ledger/energy-ledger.h5','energy_ledger','committed_records','records'),('detectors/field.h5','field_plane','committed_samples','time_s')]:
  parent_path=P/'continuous'/relative;restart_path=P/'restart'/relative
  prefix=next(p['records'] for p in snaps[0]['prefixes'] if p['file']==str(parent_path))
  with h5py.File(parent_path) as fa,h5py.File(restart_path) as fb:
   ga,gb=fa[group],fb[group];require(int(ga['complete'][()])==int(gb['complete'][()])==0,'expected bounded incomplete flag')
   ca,cb=int(ga[commit][()]),int(gb[commit][()]);require(ca==len(ga[data]) and cb==len(gb[data]) and ca-prefix==cb,'commit/prefix lengths')
   compared=[]
   for name in ga:
    if name in [commit,'complete']:continue
    av,bv=ga[name][prefix:],gb[name][:];finite(av,name);finite(bv,name);bits(av,bv,relative+'/'+name);compared.append(name)
   attrs=set(ga.attrs)-identity;require(attrs==set(gb.attrs)-identity,'attribute sets')
   for name in attrs:bits(np.asarray(ga.attrs[name]),np.asarray(gb.attrs[name]),relative+'/attr/'+name)
   item={'file':relative,'parent_records':ca,'prefix_records':prefix,'restart_records':cb,'datasets_exact_suffix':compared,'scientific_attributes_exact':len(attrs),'complete':0}
   if group=='energy_ledger':
    item['parent_steps']=ga['records']['step'].tolist();item['restart_steps']=gb['records']['step'].tolist();require(item['parent_steps']==[0,3,6,9,12] and item['restart_steps']==[6,9,12],'ledger steps')
   else:
    item['parent_time_s']=ga['time_s'][:].tolist();item['restart_time_s']=gb['time_s'][:].tolist();require(np.diff(ga['time_s'][:]).min()>0,'monotonic field clock')
   output_reviews.append(item)
 writes=[];log_reviews=[]
 for run in ['continuous','restart']:
  log=(P/run/'aprl.log').read_text();launch=json.loads((P/run/'launch.json').read_text());require(launch['exit_code']==1,'duration exit1')
  require('Compute threads per MPI rank: 1;' in log and 'ranks=256 grid_cells=128x128x2048 particles=32768' in log,'actual layout')
  require('Shutdown watchdog: enabled=true, timeout=600 s, finalize_timeout_seconds=60 s' in log,'teardown safeguards')
  require('mesh.duration was exhausted before the configured physical stop' in log,'expected duration exhaustion')
  phases=['trajectory-close-complete','detector-close-complete','ledger-close-complete','final-resource-report-complete','duration-exhaustion-readiness-sync-enter','duration-exhaustion-readiness-sync-complete']
  positions=[]
  for phase in phases:
   match=re.search(r'\[teardown\] rank=0 host=\S+ .* phase='+phase+r'\b',log);require(match,'root close trace '+run+' '+phase);positions.append(match.start())
  require(positions==sorted(positions),'close phase order')
  for m in re.finditer(r'\[checkpoint\] phase=write step=(\d+) time_box_s=(\S+) bytes=(\d+) seconds=(\S+) gib_per_s=(\S+) path=(\S+)',log):
   step,time_box,bytecount,seconds,bandwidth,path=m.groups();snap=next(s for s in snaps if s['path']==path)
   require(int(bytecount)==snap['payload_bytes'],'logged payload bytes');require(abs(float(bandwidth)-int(bytecount)/(float(seconds)*2**30))<1e-9,'bandwidth arithmetic')
   writes.append({'run':run,'step':int(step),'payload_bytes':int(bytecount),'file_bytes':snap['file_bytes'],'pause_seconds':float(seconds),'payload_gib_s':float(bandwidth),'file_gib_s':snap['file_bytes']/(float(seconds)*2**30)})
  log_reviews.append({'run':run,'exit_code':1,'final_step':12,'launch_wall_seconds':launch['wall_seconds'],'root_close_and_readiness_order_verified':True})
 require(len(writes)==3 and not list(P.rglob('*.partial')),'three published complete snapshots')
 review={'status':'PASS','reviewed_at':datetime.datetime.now().astimezone().isoformat(),'reviewer':'independent checkpoint_io_design','method':'Independent chunked HDF5 reads, raw-log recomputation and persisted prefix audit; original runner not imported','chunk_bytes':CHUNK,'frozen_files_verified':len(frozen),'binary_sha256':plan['binary_sha256'],'job':job_evidence,'audit_script_sha256':sha(__file__),'binding_ranks':256,'threads':1,'same_step_datasets_bitwise_equal':same,'same_step_native_scalars_equal_excluding_only_114_bytes_x87_padding_per_rank':True,'snapshots':snaps,'outputs':output_reviews,'writes':writes,'logs':log_reviews,'limitations':['One moderate generated 128x128x2048,32768-particle two-node case; not production-scale bandwidth.','Parent and restart deliberately stop after 12 steps by duration; output complete=0 is expected, committed datasets are complete/readable.','CPML histories and carrier datasets are preserved; nonzero counts below delimit actual excitation coverage.','Three writes share one filesystem session; caches may be warm. No power-loss durability claim.','Native long-double padding is excluded; all 80 significant bits per long double and all other scalar bytes agree.'],'evidence_sha256':{str(p):sha(p) for p in [P/'RESULT.json',P/'binding.log',P/'continuous/aprl.log',P/'restart/aprl.log',P/'frozen.json']}}
 with (P/output_name).open('x') as out:json.dump(review,out,indent=2,sort_keys=True);out.write('\n')
 print(json.dumps({'status':'PASS','review':str(P/output_name),'bitwise_datasets':len(same)}),flush=True)

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--root',required=True,type=Path)
 parser.add_argument('--output-name',default='REVIEW.json')
 args=parser.parse_args()
 review_directory(args.root.resolve(),args.output_name)

if __name__=='__main__':
 main()
