"""Audit only closed historical input H5; no simulation and no output H5 reads."""
import json
from pathlib import Path
import sys
import h5py
import numpy as np


def main(workspace,output):
    workspace=Path(workspace).resolve();output=Path(output).resolve()
    sys.path.insert(0,str(workspace/'source/aprl/verification/performance/cases/sparc_sase'))
    from loading_guard import inspect_cell
    files=[workspace/'results/sparc_sase_comparison_20261103/compact-preflight-r1/particles.h5',
        workspace/'results/sparc_sase_comparison_20261103/campaign-r4/particles.h5',
        workspace/'results/sparc_sase_11segment_2slip_20261111_r1/particles.h5']
    results=[]
    for path in files:
        item=dict(input=str(path),checked_records=0,checked_cells=0,column_particles=0,min_unique_xy=1024,
            max_unique_xy=0,closed_historical_input_only=True,no_running_output_read=True)
        with h5py.File(path,'r') as f:
            ds=f['particles/records'];assert len(ds)%1024==0
            for start in range(0,len(ds),65536):
                block=ds[start:start+65536]
                for lo in range(0,len(block),1024):
                    a=block[lo:lo+1024];xy=a['plane_position_m']
                    v=inspect_cell(xy,np.column_stack([xy,a['proper_velocity']]),a['arrival_time_offset_s'])
                    item['checked_records']+=len(a);item['checked_cells']+=1;item['column_particles']+=v['column_particles']
                    item['min_unique_xy']=min(item['min_unique_xy'],v['unique_xy']);item['max_unique_xy']=max(item['max_unique_xy'],v['unique_xy'])
        item['column_particle_fraction']=item['column_particles']/item['checked_records']
        item['independent_loading_admission_passed']=item['column_particles']==0
        results.append(item);print(json.dumps(item),flush=True)
    (output/'HISTORICAL_LOADING_AUDIT.json').write_text(json.dumps(dict(complete=True,results=results),indent=2)+'\n')


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
