"""Make a delivery containing only native absorption and its large references."""
import csv
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from common import SOURCE,sha,save,table

PREVIOUS=SOURCE.parent/'puffin_boundary_explanation_20261117_r2'
CASES=['small-native','wide-reference','wider-reference']

def main(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False)
    (root/'data').mkdir();(root/'metadata').mkdir()
    old_package=json.loads((PREVIOUS/'PACKAGE.json').read_text())
    assert sha(PREVIOUS.with_suffix('.zip'))==old_package['sha256']
    kept={}
    for name in ['comparison.csv','history.csv']:
        with (PREVIOUS/'data'/name).open() as f:rows=list(csv.DictReader(f))
        rows=[r for r in rows if r['case'] in CASES]
        table(root/'data'/name,rows);kept[name]=len(rows)
    with np.load(PREVIOUS/'data/visual_arrays.npz') as a:
        selected={k:a[k] for k in ['x_m','z_m','snapshot_z_m','reference','native','reference_line','native_line']}
        np.savez_compressed(root/'data/visual_arrays.npz',**selected)
    result=json.loads((PREVIOUS/'data/RESULT.json').read_text())
    result['metrics']={k:v for k,v in result['metrics'].items() if k in CASES}
    result['validation'].pop('no_absorption_global_norm_drift',None)
    result['delivery_projection']=True
    save(root/'data/RESULT.json',result)
    summary=json.loads((PREVIOUS/'SUMMARY.json').read_text())
    summary.pop('small-no-absorption');summary.pop('periodic_image_validation')
    summary['reference'].pop('no_absorption_global_norm_drift',None)
    save(root/'SUMMARY.json',summary)
    save(root/'SELECTION.json',dict(cases=CASES,displayed_cases=['wide-reference','small-native'],arrays=list(selected),
        retained_source_values_unchanged=True,source_delivery=str(PREVIOUS),source_zip_sha256=old_package['sha256']))
    for f in (PREVIOUS/'metadata').iterdir():
        if not f.is_file() or f.name.startswith('small-no-absorption'):continue
        if f.name=='PROTECTED_SOURCES.json':
            raw=json.loads(f.read_text());raw={k:v for k,v in raw.items() if '/small-no-absorption/' not in k}
            save(root/'metadata'/f.name,raw)
        elif f.name in ['SOURCE_PLAN.json','SOURCE_PHYSICS_COMPLETE.json']:
            raw=json.loads(f.read_text());raw['cases']=[c for c in raw['cases'] if c['name'] in CASES+['plane-preflight']]
            raw['delivery_projection']=True;save(root/'metadata'/f.name,raw)
        else:shutil.copyfile(f,root/'metadata'/f.name)
    shutil.copytree(PREVIOUS/'background',root/'background')
    save(root/'metadata/SELECTION_PROVENANCE.json',dict(previous_zip_sha256=old_package['sha256'],
        source_files={str(PREVIOUS/'data'/f):sha(PREVIOUS/'data'/f) for f in ['comparison.csv','history.csv','visual_arrays.npz','RESULT.json']},
        rows=kept,arrays=list(selected),case_values_preserved=True,
        json_metadata_is_selected_projection=True))
    print(json.dumps(dict(rows=kept,arrays=list(selected)),indent=2))

if __name__=='__main__':main(sys.argv[1])
