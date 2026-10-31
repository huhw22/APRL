import re
def binding(path,nodes,ranks=16,threads=8):
 text=path.read_text();rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)',text)
 assert '[binding] result=PASS failures=0' in text and len(rows)==ranks*threads,'Incomplete binding check'
 assert all(x[-1]=='PASS' for x in rows)
 assert {(int(x[0]),int(x[1])) for x in rows}=={(r,t) for r in range(ranks) for t in range(threads)}
 assert len({x[2] for x in rows})==nodes
 assert len({(x[2],x[4],x[5]) for x in rows})==ranks*threads,'Physical core overlap'
 return {'passed':True,'nodes':sorted({x[2] for x in rows}),'ranks':ranks,'threads':threads,'physical_cores':ranks*threads}
