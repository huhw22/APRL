"""Repair postprocessing only; never touch old cards, outputs, or failed records."""
import argparse,json
from pathlib import Path
from prepare import require,save
from analyze_fel import analyze
from genesis_aperture import compare

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('--only',nargs='+')
    a=ap.parse_args();root=a.root.resolve()
    plan=json.loads((root/'campaign.json').read_text())
    records=[]
    for case in plan['cases']:
        if a.only and case['name'] not in a.only:continue
        d=Path(case['directory'])
        require((d/'PIC_COMPLETE.json').exists(),'PIC has not closed: '+case['name'])
        launch=json.loads((d/'PIC_COMPLETE.json').read_text())
        require(launch['exit_code']==0,'Nonzero PIC exit')
        out=d/'analysis-regularized-r1'
        if not (out/'RESULT.json').exists():
            require(not out.exists(),'Incomplete repair output; inspect before retry')
            analyze(root,case,out.name)
        compare(root,case,out.name)
        records.append(dict(case=case['name'],analysis=str(out),PIC_closed=True,physical_acceptance=False))
        save(root/'REPAIRED_ANALYSIS_STATUS.json',records)
    print(json.dumps(records,indent=2))

if __name__=='__main__':main()
