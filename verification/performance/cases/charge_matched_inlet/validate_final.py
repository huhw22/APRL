#!/usr/bin/env python3
"""Independent final identity/card/report checks, without rereading running H5."""
import argparse,csv,hashlib,json,re
from pathlib import Path
import numpy as np,yaml

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def differences(a,b,p=''):
    if isinstance(a,dict) and isinstance(b,dict):
        return [x for k in a.keys()|b.keys() for x in differences(a.get(k),b.get(k),p+'.'+k)]
    return [] if a==b else [dict(path=p,old=a,new=b)]
def main(root,aux,audit):
    frozen={}
    for d in (root,aux,audit):
        original=read(d/'FROZEN.json');assert all(sha(p)==h for p,h in original.items()),d
        frozen[str(d)]=len(original)
    cards={};allowed={'.trajectory.directory','.detectors.directory','.mesh.moving_window.inlet_charge_match','.mesh.moving_window.reference_directory','.energy_ledger.directory','.output.manifest'}
    for k in ('uniform','converging','gaussian'):
        old=yaml.safe_load((root.parent/'reference_surface_current_20261030'/(k+'_full/run.yaml')).read_text());new=yaml.safe_load((root/(k+'_full/run.yaml')).read_text());ds=differences(old,new)
        assert {x['path'] for x in ds}==allowed
        assert new['mesh']['moving_window']['lower_guard_cells']==64
        cards[k]=ds
    (root/'CARD_AUDIT.json').write_text(json.dumps(cards,indent=2)+'\n')
    result=read(root/'FINAL.json');assert not result['production_spectrum_accepted']
    for row in csv.DictReader((root/'comparison.csv').open()):
        values=result['cases'][row['case']]
        for n,v in row.items():
            if n!='case':assert np.isfinite(float(v)) and float(v)==values[n]
    links=[]
    for p in (root/'README.md',root/'PROGRAM_VALIDATION.md'):
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            assert (p.parent/target).exists(),target
            links.append(target)
    artifacts={str(p.relative_to(root)):sha(p) for p in [root/'README.md',root/'PROGRAM_VALIDATION.md',root/'FINAL.json',root/'comparison.csv',root/'CARD_AUDIT.json',*(root/'figures').glob('*')] if p.is_file()}
    (root/'ARTIFACTS.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    check=dict(pass_all=True,original_frozen_hashes=frozen,full_cards_only_inlet_flag_and_output_paths_changed=True,report_links_checked=len(links),csv_json_exact=True,figure_files=list(artifacts),small_case_physics_convergence_accepted=False)
    (root/'VALIDATION.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('aux',type=Path);p.add_argument('audit',type=Path);a=p.parse_args();main(a.root.resolve(),a.aux.resolve(),a.audit.resolve())
