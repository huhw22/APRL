#!/usr/bin/env python3
"""Check the closed full pair, including two particle and energy planes."""
import argparse,json,os,re,sys
from pathlib import Path
import h5py,numpy as np,yaml
from completion_helpers import digest,finite,require,string,check_identity,check_binding,check_time_axis,check_window_csv
import report_discrete

def particle_file(path,expected,manifest,input_source_ids,chunk_rows=65536):
    require(len(input_source_ids)==expected,'Input source-ID count differs')
    require(chunk_rows>0,'Particle scan chunk size must be positive')
    with h5py.File(path) as h:
        g=h['particle_plane'];check_identity(g,manifest)
        data=g['records'];count=int(g['committed_records'][()])
        require(int(g['complete'][()])==1 and count==len(data)==expected,'Particle plane incomplete: '+str(path))
        seen=np.zeros(expected,dtype=bool)
        for start in range(0,count,chunk_rows):
            a=data[start:start+chunk_rows];finite(a,str(path))
            # particle_id indexes input rows; source_id preserves the arbitrary
            # SDDS particleID, including sparse IDs left by input cropping.
            ids=a['particle_id'];require(np.all((ids>=1)&(ids<=expected)),'Particle row IDs out of range')
            index=ids.astype(np.int64)-1
            require(np.unique(index).size==len(index) and not np.any(seen[index]),'Duplicate particle row IDs')
            require(np.array_equal(a['source_id'],input_source_ids[index]),'Source IDs differ from input rows')
            seen[index]=True
            require(np.all(a['crossing_direction']==1),'Backward particle crossing')
        require(np.all(seen),'Missing/duplicate particle row IDs')
    return dict(file=str(path),records=count,unique_input_IDs=expected,complete=True,
                identity_rule='Unique particle_id indexes input rows; source_id matches the original arbitrary SDDS IDs')

def energy_file(path,expected,manifest):
    with h5py.File(path) as h:
        g=h['energy_plane'];check_identity(g,manifest);a=g['records'][:]
        require(int(g['complete'][()])==1 and int(g['committed_records'][()])==len(a)>0,'Energy plane incomplete')
        finite(a,str(path));require(np.all(np.diff(a['time_box_s'])>0),'Energy plane time not increasing')
        last=a[-1];require(int(last['particle_downstream_count'])==expected and int(last['particle_upstream_count'])==0,'Energy-plane particle coverage differs')
        values={n:str(np.longdouble(last[n])+np.longdouble(last[n+'_low'])) for n in ['particle_kinetic_lab_J','particle_total_lab_J','particle_cpz_lab_J','particle_kinetic_box_J','particle_cpz_box_J']}
        return dict(file=str(path),complete=True,records=len(a),particles=expected,
                    field_samples=int(last['field_sample_count']),outside_field_window_steps=int(last['outside_field_window_steps']),
                    field_gap_count=int(last['field_gap_count']),extended_particle_totals=values,
                    interpretation='High+low totals retained; raw H5 remains authoritative. Missing field coverage is reported, not silently integrated.')

def validate_frozen(root):
    """Retain the original freeze, with at most one documented checker repair.

    Only the post-close checker may change. The original checker is archived;
    every binary, card, helper and other original dependency remains strict.
    """
    root=Path(root).resolve();manifest_path=root/'frozen.json'
    original=json.loads(manifest_path.read_text());checker=(root.parent/'scripts/check_completion.py').resolve()
    changes=[name for name,expected in original.items() if digest(name)!=expected]
    result={'original_files':len(original),'original_frozen_manifest_sha256':digest(manifest_path),
            'unchanged_files':len(original)-len(changes),'amendment_applied':False}
    if not changes:return result
    require(changes==[str(checker)],'Frozen file changed outside the post-close checker: '+str(changes))
    amendment_path=root/'postcheck-amendments/source-id-rule.json'
    a=json.loads(amendment_path.read_text())
    require(a['schema']=='aprl-postcheck-only-amendment-v1','Unknown postcheck amendment')
    require(a['role']=='post-close-validator-only','Amendment changes simulation role')
    require(Path(a['case_directory']).resolve()==root,'Amendment case differs')
    require(Path(a['checker_path']).resolve()==checker,'Amendment checker differs')
    require(a['original_frozen_manifest_sha256']==digest(manifest_path),'Original freeze manifest changed')
    require(a['original_checker_sha256']==original[str(checker)],'Amendment original checker differs')
    archive=Path(a['original_checker_archive']).resolve()
    require(archive.is_relative_to(root/'postcheck-amendments'),'Original checker archive outside amendment directory')
    require(digest(archive)==original[str(checker)],'Archived original checker changed')
    require(digest(checker)==a['replacement_checker_sha256'],'Replacement checker changed')
    result.update(amendment_applied=True,amendment_file=str(amendment_path),amendment_sha256=digest(amendment_path),
                  original_checker_archive=str(archive),replacement_checker_sha256=digest(checker))
    return result

def check(root,output=None):
    root=Path(root).resolve();out=dict(pass_=False,passed=False,physics_accepted=False)
    try:
        require(bool(os.environ.get('SLURM_JOB_ID')),'Full output scan must run on compute allocation')
        p=json.loads((root/'run-plan.json').read_text());cfg=yaml.safe_load((root/'run.yaml').read_text());geo=p['geometry']
        out['frozen_input_validation']=validate_frozen(root)
        require(p['format']=='aprl-corrected-full-pair-v1','Wrong plan schema')
        st=Path(p['input_file']).stat();require((st.st_size,st.st_mtime_ns)==(p['input_size_bytes'],p['input_mtime_ns']),'Input stat changed')
        identity=json.loads((root/'identity-check.json').read_text());require(identity['passed'],'Startup identity failed')
        manifest=yaml.safe_load((root/'run-manifest.yaml').read_text())
        require(manifest['source_revision']==p['compiled_revision'] and manifest['mpi_size']==p['ranks'],'Binary/layout differs')
        require(yaml.safe_load(manifest['configuration_yaml'])==cfg,'Manifest card differs')
        require(manifest['resolved_mesh']['cells']==[160,160,geo['active_cells_z']],'Active mesh changed')
        out['binding']=check_binding(root/'binding.log',p)
        mw=cfg['mesh']['moving_window'];require(mw['reference_surface_current'] and mw['inlet_charge_match'],'Expected corrections disabled')
        require(cfg['energy_ledger']['discrete_operator_audit'],'Discrete audit disabled')
        log=(root/'aprl.log').read_text();require('Simulation failure' not in log and 'mesh.duration was exhausted' not in log,'Simulation failed/duration exhausted')
        require(re.search(r'Pending admission queue: waiting=\d+ resident=\d+;',log) is not None,'Pending admission fallback used')
        require('Reference transverse surface current:' in log and 'Inlet charge completion:' in log,'Corrections not active')
        final=re.findall(r'\[resource\] phase=final step=(\d+)\b',log)
        require(len(final)==1 and 'phase=solve-complete' in log and 'phase=mpi-finalize-complete' in log,'Normal completion missing')
        step=int(final[0]);require(0<step<geo['reference_query_steps'],'Step outside planned duration')
        stop=re.findall(r'Configured stop reached: all (\d+) valid particles passed the last interaction boundary at z=([^ ]+) m',log)
        require(len(stop)==1 and int(stop[0][0])>0 and np.isclose(float(stop[0][1]),p['field_plane_z_m'],rtol=2e-6),'Physical stop differs')
        reference=Path(mw['reference_directory']);prep=json.loads((reference/'preparation.json').read_text())
        require(prep['complete'] and prep['particles']==p['records'],'Reference preparation incomplete')
        require((reference/'transverse-current.bin').stat().st_size>0,'Surface-current reference missing')
        out['moving_window']=check_window_csv(reference,p['ranks'],step,cfg['energy_ledger']['sample_interval_steps'],geo)
        field=cfg['detectors']['field_planes'][0];path=root/'detectors'/(field['name']+'.h5')
        with h5py.File(path) as h:
            g=h['field_plane'];check_identity(g,manifest);data=g['fields'];ts=g['time_s'][:]
            require(int(g['complete'][()])==1 and int(g['committed_samples'][()])==len(data)==len(ts)==step+1,'Field transaction incomplete')
            require(data.shape==(step+1,160,160),'Field aperture differs')
            require(string(g.attrs['external_background_subtracted'])=='none','Expected raw total field')
            require('adjacent MPI half-cell planes included' in string(g.attrs['spatial_sampling']),'MPI sampling correction absent')
            require('stored leapfrog half time' in string(g.attrs['magnetic_time_stagger']),'B time staggering differs')
            cadence=geo['dt_box_s']/geo['gamma_boost'];check_time_axis(ts,cadence)
            for start in range(0,len(data),24):finite(data[start:start+24],str(path))
            out['field']=dict(file=str(path),complete=True,samples=len(ts),time_start_s=float(ts[0]),time_end_s=float(ts[-1]),all_values_finite=True)
        expected_files={path}
        out['particle_planes']=[]
        with h5py.File(p['input_file']) as h:input_source_ids=h['particles/records'].fields('source_id')[:]
        require(len(input_source_ids)==p['records'],'Input ID count differs')
        for plane in cfg['detectors']['particle_planes']:
            path=root/'detectors'/(plane['name']+'.h5');expected_files.add(path)
            out['particle_planes'].append(particle_file(path,p['records'],manifest,input_source_ids))
        out['energy_planes']=[]
        for plane in cfg['detectors']['energy_planes']:
            path=root/'detectors'/(plane['name']+'-energy.h5');expected_files.add(path)
            out['energy_planes'].append(energy_file(path,p['records'],manifest))
        require(set((root/'detectors').glob('*.h5'))==expected_files,'Unexpected/missing detector file')
        ledger=root/'ledger';rows,sums,ranks,steps=report_discrete.load(ledger,ledger/'energy-ledger.h5')
        require(ranks==p['ranks'] and steps==step,'Strict ledger dimensions differ')
        _,_,summary=report_discrete.calculate(rows,sums,steps);out['discrete_ledger']=summary
        finite(rows,'physical ledger')
        require(np.all(rows['active_macroparticles']+rows['removed_macroparticles']==p['records']),'Ledger particle inventory differs')
        if p['case']=='free_drift_0p5m':require(np.max(np.abs(rows['prescribed_boris_work_J']))==0,'Free drift external work nonzero')
        out.update(pass_=True,passed=True,final_step=step,binary_sha256=p['binary_sha256'],output_transactions_complete=True)
        out['limitations']=['Complete finite data and discrete identity do not prove pure radiation, physical accuracy, or a closed lab control volume.','Energy-plane scalar four-momentum interpolation differs from raw particle-plane momentum interpolation.']
    except (ValueError,KeyError,OSError,TypeError,IndexError,AssertionError) as error:out['error']=str(error) or type(error).__name__
    destination=Path(output) if output is not None else root/'completion-check.json'
    destination.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    return out
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path,help='Fresh sidecar result; preserves an earlier completion record')
    args=parser.parse_args()
    if args.output is not None and args.output.exists():parser.error('Use a new output path; keep previous completion records')
    result=check(args.root,args.output);print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
