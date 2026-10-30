#!/usr/bin/env python3
"""Review closed PIC output separately from a failed postprocessing assertion.

Original FAILED.json and Slurm state are evidence and are never changed here.
"""
import argparse,json,sys
from pathlib import Path
import h5py,numpy as np

def main(root,audit):
    sys.path.insert(0,str(root/'scripts'))
    from prepare import sha,save,require
    import report_discrete
    frozen=json.loads((root/'FROZEN.json').read_text())
    require(all(sha(p)==h for p,h in frozen.items()),'Original frozen file changed')
    for kind in ('uniform','gaussian'):
        source=json.loads((audit/(kind+'_SOURCE_CHECK.json')).read_text());require(source['pass_all'],'Raw source check failed')
        d=root/(kind+'_full');original_failed=(d/'FAILED.json').read_text()
        require(json.loads((d/'aprl.log.launch.json').read_text())['exit_code']==0,'PIC did not exit normally')
        result=json.loads((d/'analysis/RESULT.json').read_text());require(result['execution_integrity_pass'],'Original physical output integrity check failed')
        files={}
        for sub in ('ledger','trajectory','detectors'):
            for p in sorted((d/sub).glob('*.h5')):
                with h5py.File(p) as h:
                    group=h[next(iter(h))];require(int(group['complete'][()])==1,'Unclosed H5: '+str(p))
                    a=group['records'];require(int(group['committed_records'][()])==len(a),'Commit/extent mismatch')
                    for lo in range(0,len(a),65536):
                        data=a[lo:lo+65536]
                        require(all(np.isfinite(data[n]).all() for n in data.dtype.names),'Non-finite output')
                    files[str(p.relative_to(d))]=len(a)
        previous=root.parent/'reference_surface_current_20261030'/(kind+'_full')
        with h5py.File(d/'detectors/entry.h5') as h,h5py.File(previous/'detectors/entry.h5') as g:
            a=h['particle_plane/records'][:];b=g['particle_plane/records'][:]
            a=a[np.argsort(a['source_id'])];b=b[np.argsort(b['source_id'])];require(a.dtype==b.dtype and len(a)==len(b)==32768,'Wrong entry schema/count')
            differences={}
            for n in a.dtype.names:
                if a[n].tobytes()==b[n].tobytes():continue
                require(n=='proper_velocity','Unexpected entry diagnostic difference')
                aa,bb=a[n],b[n];changed=aa!=bb
                require((np.signbit(aa[changed])==np.signbit(bb[changed])).all(),'Sign mismatch')
                ulps=np.abs(np.abs(aa[changed]).view(np.int64)-np.abs(bb[changed]).view(np.int64))
                differences[n]=dict(changed_components=int(changed.sum()),max_abs_difference=float(np.max(np.abs(aa-bb))),max_ulp_difference=int(ulps.max()),changed_axes=np.where(changed)[1].tolist())
            require(set(differences)=={'proper_velocity'},'Unexpected comparison result')
            require(max(differences['proper_velocity']['changed_axes'])<2,'Longitudinal entry momentum changed')
            # This reports diagnostic differences, without declaring entry records bitwise.
            save(d/'REVIEWED_ENTRY_DIAGNOSTIC.json',dict(all_record_fields_bitwise=False,other_fields_bitwise=True,records=32768,differences=differences,actual_first_push_source_check=str(audit/(kind+'_SOURCE_CHECK.json')),
                mechanism='Input-plane endpoint rounding can cause the plane sampler to interpolate a tiny active segment; that sampled momentum depends on the first post-entry force. Source states are checked directly and bitwise in separate short replays.'))
        r,s,_,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');_,_,summary=report_discrete.calculate(r,s,steps)
        save(d/'REVIEWED_LEDGER_CHECK.json',summary)
        require((d/'FAILED.json').read_text()==original_failed,'Original failure modified')
        save(d/'REVIEWED_CLOSED_OUTPUT.json',dict(status='CLOSED_PIC_OUTPUT_INDEPENDENTLY_REVIEWED',original_wrapper_failed=True,original_failure_preserved=True,pic_exit_code=0,steps=result['steps'],files_complete_committed_finite=files,actual_source_states_bitwise=True,source_audit_job=source['job_id'],entry_sampler_records_bitwise=False,physics_accepted=False))
    require(all(sha(p)==h for p,h in frozen.items()),'Original frozen file changed')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('audit',type=Path);a=p.parse_args();main(a.root.resolve(),a.audit.resolve())
