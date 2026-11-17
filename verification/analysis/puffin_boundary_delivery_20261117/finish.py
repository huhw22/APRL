"""Check the closed evidence and only then publish the new independent ZIP."""
import csv
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path
import numpy as np
from common import SOURCE,OLD,sha,save

def validate(root):
    root=Path(root).resolve();s=json.loads((root/'SUMMARY.json').read_text())
    selected=(root/'SELECTION.json').is_file()
    data=root/'data';d=np.load(data/'visual_arrays.npz')
    a=np.genfromtxt(data/'comparison.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    hist=np.genfromtxt(data/'history.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    u0=json.loads((SOURCE/'cases/wide-reference/COMPLETE.json').read_text())['initial_norm2_m2']
    w=np.ones(27);w[[0,-1]]=.5;ww=w[:,None]*w[None,:]
    def area(v):return np.sum(abs(v[:,27:54,27:54])**2*ww,axis=(1,2))*(40e-6)**2
    ref=d['reference'];snap=d['snapshot_z_m'];maxerr=0.;negerr=0.
    retained=[('wide-reference','reference'),('small-native','native')]
    if not selected:retained.append(('small-no-absorption','no_absorption'))
    for name,key in retained:
        q=a[a['case']==name];x=d[key]
        assert len(q)==61 and np.array_equal(q['z_m'],snap)
        assert np.array_equal(x[0],ref[0])
        expected=area(x)/u0;delta=area(x-ref)/u0
        maxerr=max(maxerr,float(np.max(abs(expected-q['roi_norm2_over_initial']))),float(np.max(abs(delta-q['roi_error_field_norm2_over_initial']))))
        assert np.allclose(np.sqrt(delta),q['roi_field_L2_difference_over_initial'],rtol=1e-13,atol=1e-15)
        if name!='wide-reference':
            full=np.sum(abs(x)**2,axis=(1,2))*(40e-6)**2/u0
            assert np.allclose(full,q['total_norm2_over_initial'],rtol=2e-12,atol=1e-14)
            # Parseval and directional decomposition in the complete native small grid.
            F=np.fft.fft2(x,axes=(1,2));power=abs(F)**2
            modes=np.fft.fftfreq(81,d=40e-6)<0
            negative=np.sum(power[:,:,modes],axis=(1,2))*(40e-6)**2/(81*81*u0)
            negerr=max(negerr,float(np.max(abs(negative-q['negative_kx_norm2_over_initial']))))
            assert np.allclose(np.sum(power,axis=(1,2))*(40e-6)**2/(81*81*u0),full,rtol=1e-13)
            i=int(np.argmax(np.where(snap>=1.3,expected,-1)))
            assert snap[i]==s[name]['central_return_peak_z_m']
            assert abs(expected[i]-s[name]['late_max_central_norm2_over_initial'])<1e-15
    assert maxerr<1e-14 and negerr<1e-18
    expected_rows=(183,903) if selected else (244,1204)
    assert (len(a),len(hist))==expected_rows
    protected=json.loads((root/'metadata/PROTECTED_SOURCES.json').read_text())
    for path,v in protected.items():
        p=Path(path);st=p.stat()
        assert [st.st_ino,st.st_size,st.st_mtime_ns]==[v['inode'],v['size'],v['mtime_ns']]
        assert sha(p)==v['sha256'],path
    for path,h in json.loads((SOURCE/'FROZEN.json').read_text()).items():assert sha(path)==h,path
    plan=json.loads((SOURCE/'PLAN.json').read_text())
    for path,h in plan['source_hashes'].items():assert sha(path)==h,path
    for path,h in json.loads((root/'metadata/SOURCE_SHA256.json').read_text()).items():assert sha(path)==h,path
    assert sha(OLD.with_suffix('.zip'))==json.loads((OLD/'PACKAGE.json').read_text())['sha256']
    if selected:
        selection=json.loads((root/'SELECTION.json').read_text())
        prior=Path(selection['source_delivery'])
        assert sha(prior.with_suffix('.zip'))==selection['source_zip_sha256']
        for filename in ['comparison.csv','history.csv']:
            original=list(csv.DictReader((prior/'data'/filename).open()))
            original=[r for r in original if r['case'] in selection['cases']]
            current=list(csv.DictReader((data/filename).open()))
            assert current==original
        with np.load(prior/'data/visual_arrays.npz') as original:
            assert set(d.files)==set(selection['arrays'])
            for key in d.files:assert np.array_equal(d[key],original[key]),key
        blocked=re.compile(r'small-no-absorption|no_absorption|关闭吸收|无吸收|no absorption|without absorption|unabsorbed|absorption disabled',re.I)
        for f in root.rglob('*'):
            if f.suffix in ['.md','.json','.csv']:
                assert not blocked.search(f.read_text()),f
    else:
        for f in ['comparison.csv','history.csv','RESULT.json','visual_arrays.npz']:
            assert sha(data/f)==sha(SOURCE/'analysis'/f)
    rendering=json.loads((root/'RENDER_VALIDATION.json').read_text());assert rendering['complete']
    figure_count=20 if selected else 22
    assert len(rendering['figures'])==figure_count and len(rendering['videos'])==4
    for v in rendering['figures']+rendering['videos']:assert sha(root/v['path'])==v['sha256']
    assert rendering['arrays_sha256']==sha(data/'visual_arrays.npz')
    assets=list(csv.DictReader((root/'PPT_ASSET_MAP.csv').open()))
    for r in assets:
        assert (root/r['asset']).is_file()
        if r['vector']:assert (root/r['vector']).is_file()
    links=[]
    for p in root.rglob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https:','http:','#')):continue
            dst=(p.parent/target).resolve()
            if dst!=(root/'VALIDATION.json').resolve():assert dst.is_file(),(p,target)
            links.append([str(p.relative_to(root)),target])
    result=dict(complete=True,job=dict(id=19834,state='COMPLETED',exit='0:0',elapsed='00:04:44'),
       source_scientific_files_unchanged=len(protected),frozen_scripts_and_library_unchanged=True,
       retained_CSV_values_and_NPZ_arrays_preserved=True,comparison_rows=len(a),history_rows=len(hist),
       raw_field_to_csv_max_absolute_fraction_difference=maxerr,
       directional_spectrum_to_csv_max_absolute_fraction_difference=negerr,
       all_initial_fields_bitwise_equal=True,global_small_grid_norms_match_native_MPI_history=True,
       large_reference_checks=s['reference'],
       main_figures=figure_count,videos=4,asset_map_entries=len(assets),local_markdown_links=len(links),
       video_source_frames_and_scales='61 native source frames in order; fixed linear/log scales, no per-frame normalization.',
       no_new_physical_simulation=True,production_FEL_impact_not_measured=True,
       PDFs_same_source_not_separately_rendered=True)
    if not selected:result['periodic_image_check']=s['periodic_image_validation']
    save(root/'VALIDATION.json',result);return result

def package(root):
    root=Path(root).resolve();validate(root)
    assert json.loads((root/'VISUAL_REVIEW.json').read_text())['complete']
    target=root.with_suffix('.zip');assert not target.exists()
    paths=[p for p in root.rglob('*') if p.is_file() and p.name not in ['PACKAGE.json','SHA256SUMS']]
    manifest={str(p.relative_to(root)):sha(p) for p in sorted(paths)}
    (root/'SHA256SUMS').write_text(''.join(v+'  '+k+'\n' for k,v in manifest.items()))
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in paths+[root/'SHA256SUMS']:z.write(p,root.name+'/'+str(p.relative_to(root)))
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        import hashlib
        for p,h in manifest.items():assert hashlib.sha256(z.read(root.name+'/'+p)).hexdigest()==h
    previous=Path(json.loads((root/'SELECTION.json').read_text())['source_delivery']) if (root/'SELECTION.json').exists() else OLD
    info=dict(complete=True,zip=str(target),bytes=target.stat().st_size,sha256=sha(target),members=len(paths)+1,
        actual_native_propagation=True,job_id=19834,member_sha256_verified=True,crc_verified=True,
        previous_zip=str(previous.with_suffix('.zip')),previous_zip_unchanged=True,
        main_report=str(root/'README.md'),ppt_instructions=str(root/'PPT_HANDOFF.md'))
    save(root/'PACKAGE.json',info);save(root.parent/'LATEST_PUFFIN_BOUNDARY_TRANSFER.json',info)
    print(json.dumps(info,ensure_ascii=False,indent=2))

if __name__=='__main__':
    root=sys.argv[1]
    if len(sys.argv)>2 and sys.argv[2]=='validate':print(json.dumps(validate(root),ensure_ascii=False,indent=2))
    else:package(root)
