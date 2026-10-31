import json,os,sys
cores={};record={}
for line in open('/proc/cpuinfo'):
 if ':' in line:
  k,v=line.split(':',1);record[k.strip()]=v.strip()
 elif record:
  cores[int(record['processor'])]=(record.get('physical id','0'),record.get('core id',record['processor']));record={}
if record:cores[int(record['processor'])]=(record.get('physical id','0'),record.get('core id',record['processor']))
affinity=sorted(os.sched_getaffinity(0))
open(sys.argv[1],'w').write(json.dumps({'logical_cpu_affinity':affinity,'distinct_physical_cores':len({cores[i] for i in affinity})},indent=2)+'\n')
