#!/usr/bin/env python3
"""Focused validation of particle identity and the post-check-only amendment.

All manufactured records live in a temporary directory and are not physics
evidence. Production scans remain a separate compute-node action.
"""
import argparse,importlib.util,json,sys,tempfile
from pathlib import Path
import h5py,numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--checker',type=Path,required=True)
    p.add_argument('--helpers',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    sys.path.insert(0,str(a.helpers))
    spec=importlib.util.spec_from_file_location('repaired_checker',a.checker)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    outcomes={}
    def expect(name,accepted,call):
        try:call();actual=True
        except (ValueError,KeyError,OSError,AssertionError):actual=False
        assert actual==accepted,(name,actual,accepted)
        outcomes[name]={'accepted':actual,'expected_acceptance':accepted}
    with tempfile.TemporaryDirectory(prefix='aprl-id-repair-') as d:
        base=Path(d);case=base/'case';case.mkdir();scripts=base/'scripts';scripts.mkdir()
        identity={'run_id':'manufactured-ID-rule-only','source_revision':'test','configuration_digest_fnv1a64':'test'}
        source=np.array([104097203,7,100003,9999,2,59000000],dtype='u8')
        dtype=np.dtype([('particle_id','u8'),('source_id','u8'),('crossing_direction','i4'),('value','f8')])
        records=np.zeros(6,dtype=dtype);records['particle_id']=[4,1,6,2,5,3]
        records['source_id']=source[records['particle_id'].astype(int)-1];records['crossing_direction']=1
        for name in ['sparse_source_valid','within_chunk_duplicate','cross_chunk_duplicate','source_mapping_corrupt',
                     'particle_id_zero','particle_id_above_N','missing_record','backward_crossing','nonfinite_value','incomplete_transaction']:
            values=records.copy();complete=1
            if name=='within_chunk_duplicate':values[1]=values[0]
            elif name=='cross_chunk_duplicate':values[3]=values[0]
            elif name=='source_mapping_corrupt':values['source_id'][1]+=1
            elif name=='particle_id_zero':values['particle_id'][0]=0
            elif name=='particle_id_above_N':values['particle_id'][0]=7
            elif name=='missing_record':values=values[:-1]
            elif name=='backward_crossing':values['crossing_direction'][0]=-1
            elif name=='nonfinite_value':values['value'][0]=np.nan
            elif name=='incomplete_transaction':complete=0
            path=base/(name+'.h5')
            with h5py.File(path,'w') as f:
                g=f.create_group('particle_plane');g.attrs.update(identity);g['records']=values
                g['committed_records']=len(values);g['complete']=complete
            expect(name,name=='sparse_source_valid',lambda:m.particle_file(path,6,identity,source,chunk_rows=2))
        checker=scripts/'check_completion.py';checker.write_text('original checker\n')
        card=case/'run.yaml';card.write_text('original physics card\n')
        helper=scripts/'helper.py';helper.write_text('original helper\n')
        original={str(f):m.digest(f) for f in (checker,card,helper)}
        frozen=case/'frozen.json';frozen.write_text(json.dumps(original))
        expect('strict_original_freeze',True,lambda:m.validate_frozen(case))
        amended=case/'postcheck-amendments';amended.mkdir()
        archive=amended/'check_completion.original.py';archive.write_bytes(checker.read_bytes());checker.write_text('fixed ID checker\n')
        expect('checker_change_without_record',False,lambda:m.validate_frozen(case))
        record={'schema':'aprl-postcheck-only-amendment-v1','role':'post-close-validator-only','case_directory':str(case),
          'checker_path':str(checker),'original_frozen_manifest_sha256':m.digest(frozen),
          'original_checker_sha256':original[str(checker)],'original_checker_archive':str(archive),'replacement_checker_sha256':m.digest(checker)}
        record_path=amended/'source-id-rule.json'
        def restore():
            archive.write_text('original checker\n');checker.write_text('fixed ID checker\n');card.write_text('original physics card\n');helper.write_text('original helper\n')
            record_path.write_text(json.dumps(record))
        restore();expect('documented_checker_only_change',True,lambda:m.validate_frozen(case))
        card.write_text('altered physical card\n');expect('card_change_rejected',False,lambda:m.validate_frozen(case));restore()
        helper.write_text('altered helper\n');expect('helper_change_rejected',False,lambda:m.validate_frozen(case));restore()
        archive.write_text('altered original archive\n');expect('archive_change_rejected',False,lambda:m.validate_frozen(case));restore()
        checker.write_text('second unrecorded change\n');expect('replacement_change_rejected',False,lambda:m.validate_frozen(case));restore()
        bad={**record,'original_frozen_manifest_sha256':'0'*64};record_path.write_text(json.dumps(bad))
        expect('wrong_original_manifest_rejected',False,lambda:m.validate_frozen(case));restore()
        bad={**record,'case_directory':str(base/'different-case')};record_path.write_text(json.dumps(bad))
        expect('different_case_rejected',False,lambda:m.validate_frozen(case));restore()
        bad={**record,'original_checker_archive':str(helper)};record_path.write_text(json.dumps(bad))
        expect('unrelated_archive_rejected',False,lambda:m.validate_frozen(case))
    result={'checks_pass':True,'checks':outcomes,'check_count':len(outcomes),'new_physics_runs':0,'scope':'Post-close identity and immutable-input guard only; not physical accuracy'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
