#!/usr/bin/env python3
"""Sparse-cadence and corrupted-input checks using isolated derived fixtures."""
import argparse,csv,json
from pathlib import Path
import h5py
import numpy as np
from report import load,calculate

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    a=p.parse_args();source=a.source.resolve();work=a.work.resolve();work.mkdir()
    with h5py.File(source/'energy-ledger.h5') as f:
        g=f['energy_ledger'];rows=g['records'][:];attrs=dict(g.attrs)
    indices=np.unique(np.r_[np.arange(0,len(rows),7),len(rows)-1]);fixture=work/'energy-ledger.h5'
    with h5py.File(fixture,'w') as f:
        g=f.create_group('energy_ledger')
        for k,v in attrs.items():g.attrs[k]=v
        g.attrs['derived_test_fixture']='7-step cadence; copied numerical records, not a new PIC run'
        g.create_dataset('records',data=rows[indices]);g.create_dataset('complete',data=np.uint8(1));g.create_dataset('committed_records',data=np.uint64(len(indices)))
    for path in source.glob('discrete-field-rank-*'): (work/path.name).symlink_to(path)
    full,sums,_,steps=load(source,source/'energy-ledger.h5');full_i,full_q,_=calculate(full,sums,steps)
    sparse,sparse_sums,_,_=load(work,fixture);sparse_i,sparse_q,_=calculate(sparse,sparse_sums,steps)
    assert np.array_equal(sparse_i,full_i[indices])
    assert all(np.array_equal(v,full_q[k][indices]) for k,v in sparse_q.items())
    checks={'sample_cadence_7_matches_every_selected_full_record':True}
    def rejected(label):
        try:load(work,fixture)
        except (RuntimeError,FileNotFoundError):checks[label]=True
        else:raise AssertionError('Accepted '+label)
    last=work/'discrete-field-rank-00015.csv';target=last.resolve();last.unlink();rejected('missing_rank_rejected');last.symlink_to(target)
    marker=work/'discrete-field-rank-00000.complete';target=marker.resolve();marker.unlink();marker.write_text('0 257\n');rejected('incomplete_rank_rejected');marker.unlink();marker.symlink_to(target)
    path=work/'discrete-field-rank-00000.csv';target=path.resolve()
    with target.open() as f:reader=csv.DictReader(f);records=list(reader);names=reader.fieldnames
    path.unlink();records[100]['time_electric_J']='1e-16'
    with path.open('w') as f:writer=csv.DictWriter(f,fieldnames=names);writer.writeheader();writer.writerows(records)
    rejected('corrupted_numerical_term_rejected');path.unlink();path.symlink_to(target)
    with h5py.File(fixture,'r+') as f:f['energy_ledger'].attrs['frame']='lab frame'
    rejected('wrong_frame_rejected')
    result={'checks_pass':True,'checks':checks,'physical_input_unchanged':True,'fixture_scope':'Derived records and symlinks; no new simulation'}
    (work/'READER_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
