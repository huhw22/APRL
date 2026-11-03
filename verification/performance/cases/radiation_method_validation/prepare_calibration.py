"""Only physical cards change: weaken bend, refine transverse x, halve dz/dt."""
import argparse,json,sys
from pathlib import Path
import yaml
from prepare import save,sha,require

def main():
    ap=argparse.ArgumentParser();ap.add_argument('workspace',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args();w=a.workspace.resolve();root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).parent/'prepare_point.py';base=source.read_text()
    cases=[]
    for name in ['weak-bend','x-refined','zdt-refined']:
        d=root/name
        code=base.replace('WS=ROOT.parents[2]',f'WS=Path({str(w)!r})')
        code=code.replace("'rhythm':DT/G","'rhythm':.999*DT/G")
        if name=='weak-bend':
            code=code.replace("constants(length=.0024,screen=.3024)",
                              "constants(rho=12.,length=.0024,screen=.3024)")
            code=code.replace("'registered_band_eV':[50.,100.]","'registered_band_eV':[8.,12.]")
        if name=='x-refined':
            code=code.replace('cells=[96,80,full]','cells=[240,80,full]')
            code=code.replace("'nx':96,'ny':80","'nx':240,'ny':80")
        if name=='zdt-refined':
            code=code.replace('DZ=1.2543262415359034e-6','DZ=6.271631207679517e-7')
        argv=sys.argv
        try:
            sys.argv=[str(source),str(d)]
            exec(compile(code,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
        finally:sys.argv=argv
        (d/'scripts/generated_prepare.py').write_text(code)
        for method in ['boris','energy','momentum']:
            p=d/method/'run.yaml';cfg=yaml.safe_load(p.read_text())
            if name=='x-refined':
                cfg['mesh']['cell_size'][0]=2e-6
                cfg['boundary']['cells'][0]=25 # Preserve the physical 50 um CPML.
            if name=='zdt-refined':cfg['mesh']['cell_size'][2]=6.271631207679517e-7
            p.write_text(yaml.safe_dump(cfg,sort_keys=False))
        campaign=json.loads((d/'campaign.json').read_text())
        campaign['cases']=[q for q in campaign['cases'] if q['name']=='energy']
        campaign['cases'][0]['timeout_seconds']=900
        campaign['scope']='Single-electron calibration, same coupling: distinguish transverse phase-velocity effects from longitudinal/time refinement'
        save(d/'campaign.json',campaign)
        frozen=json.loads((d/'FROZEN.json').read_text())
        frozen[str(d/'scripts/generated_prepare.py')]=sha(d/'scripts/generated_prepare.py')
        for path in list(frozen):frozen[path]=sha(path)
        save(d/'FROZEN.json',frozen)
        cases.append(str(d))
    scripts=root/'scripts';scripts.mkdir()
    (scripts/'run_sync.py').write_text(
       "import subprocess,sys\nfrom pathlib import Path\n"
       "root=Path(__file__).resolve().parents[1]\n"
       "for name in ['weak-bend','x-refined','zdt-refined']:\n"
       " print('CALIBRATION_START',name,flush=True)\n"
       " subprocess.run([sys.executable,'-B',str(root/name/'scripts/run_sync.py')],check=True)\n")
    save(root/'PLAN.json',dict(cases=cases,scope='Additional radiation calibration; no C++ changes'))
    print(root)

if __name__=='__main__':main()
