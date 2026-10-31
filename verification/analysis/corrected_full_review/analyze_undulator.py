#!/usr/bin/env python3
"""Review closed 19099 without hiding its incomplete full-input exit coverage.

Field and scalar/strict ledgers are complete. Particle statistics explicitly
condition both entry and exit on IDs that reached the exit; absent IDs retain
their own list and entry four-momentum. No event or lost-particle energy is
invented, and the original postcheck failure remains unchanged.
"""
import argparse,json,os,re,sys
from pathlib import Path
import h5py,numpy as np,yaml
import analyze as shared


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();assert os.environ.get('SLURM_JOB_ID'),'Large scans require a compute allocation'
    root=a.run.resolve();out=a.out.resolve();assert not out.exists();out.mkdir()
    plan=json.loads((root/'run-plan.json').read_text());assert plan['case']=='undulator_4p5m' and plan['ranks']==32
    prior=json.loads((root/'completion-check.json').read_text())
    assert prior.get('error')=='Particle plane incomplete: '+str(root/'detectors/exit-particles.h5')
    assert prior['field']['complete'] and prior['field']['all_values_finite']
    log=(root/'aprl.log').read_text()
    assert all(t in log for t in ['phase=solve-complete','phase=mpi-finalize-complete','phase=detector-close-complete'])
    assert 'Simulation failure' not in log
    sys.path.insert(0,str(Path(__file__).resolve().parent));import completion_validator
    frozen=completion_validator.validate_frozen(root)
    manifest=yaml.safe_load((root/'run-manifest.yaml').read_text());identities={}
    for path in sorted((root/'detectors').glob('*.h5')):
        with h5py.File(path) as f:
            group=f[next(iter(f))];completion_validator.check_identity(group,manifest)
            assert int(group['complete'][()])==1
            if 'committed_records' in group:
                assert len(group['records'])==int(group['committed_records'][()])
                identities[path.name]={'complete':True,'records':len(group['records'])}
            else:
                assert len(group['fields'])==int(group['committed_samples'][()])
                identities[path.name]={'complete':True,'samples':len(group['fields'])}
    n=plan['records'];exit_n=identities['exit-particles.h5']['records'];assert identities['entry-particles.h5']['records']==n
    assert 0<exit_n<n
    boundary_line=next(t for t in log.splitlines() if t.startswith('CPML physical-domain exits:'))
    counts={face:int(value) for face,value in re.findall(r'([xyz][+-])=(\d+)',boundary_line)}
    transverse=sum(v for face,v in counts.items() if not face.startswith('z'))
    assert n-exit_n<=transverse
    paths=[root/'completion-check.json',root/'aprl.log',root/'run.yaml',root/'run-manifest.yaml',*(root/'detectors').glob('*.h5'),root/'ledger/energy-ledger.h5']
    before={str(q):[q.stat().st_size,q.stat().st_mtime_ns] for q in paths}
    print(f'Closed run: entry={n}, exit={exit_n}, absent={n-exit_n}, total transverse exits={transverse}',flush=True)
    particles=shared.scan_particles(root,out,plan,allow_exit_subset=True)
    coverage=json.loads((out/'particle_statistics.json').read_text())
    assert coverage['available_counts']=={'entry':n,'exit':exit_n,'missing_exit':n-exit_n}
    print('Exact matched-cohort particle statistics complete',flush=True)
    energy=shared.read_energy(root,out,n,coverage=coverage)
    ledger=shared.read_ledger(root,out,'19099 | corrected full-input undulator, 4.5 m')
    with h5py.File(root/'ledger/energy-ledger.h5') as f:
        rows=f['energy_ledger/records'][:]
        assert np.all(rows['active_macroparticles']+rows['removed_macroparticles']==n)
    boundary={'original_input':n,'particle_exit_records':exit_n,'missing_exit_count':n-exit_n,
       'missing_exit_fraction':(n-exit_n)/n,'transverse_exit_counts':counts,'transverse_exit_total':transverse,
       'missing_IDs_proven_to_match_boundary_events':False,
       'limitation':'Only aggregate boundary counters were saved. The difference of counts is compatible with some transverse exits occurring after the plane, but no per-ID exit ledger proves this partition.'}
    shared.save(out/'CORE_RESULT.json',{'particles':particles,'coverage':coverage,'energy_planes':energy,'ledger':ledger,'boundary':boundary})
    shared.save(out/'COVERAGE.json',boundary)
    print('Complete ledgers and matched-cohort results saved; starting FFT',flush=True)
    spectrum=shared.spectrum(root,out,simulation_job=19099,
      completion_scope='Complete closed raw field; matched arriving-particle cohort checked. Full input-to-exit particle coverage remains short and is not silently accepted.')
    assert all([Path(q).stat().st_size,Path(q).stat().st_mtime_ns]==value for q,value in before.items())
    shared.save(out/'COMPLETE.json',{'passed':True,'scope':'Integrity of analysis of available closed outputs and exact matched arriving cohort',
       'full_input_exit_coverage_passed':False,'source_simulation_job':19099,'original_batch_exit_code':1,
       'original_postcheck_error':prior['error'],'entry_particles':n,'matched_particles':exit_n,'missing_exit_count':n-exit_n,
       'field_samples':spectrum['samples'],'diagnostic_identities':identities,'frozen_validation':frozen,
       'source_stats_unchanged':before,'new_pic_simulations':0,'analysis_job':os.environ['SLURM_JOB_ID']})
    print('COMPLETE: 19099 field/ledger and conditional particle analysis; full exit coverage remains flagged',flush=True)


if __name__=='__main__':main()
