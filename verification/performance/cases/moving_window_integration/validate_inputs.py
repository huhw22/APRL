#!/usr/bin/env python3
"""Bounded rejection checks; never enters Maxwell stepping or modifies input files."""
import argparse, copy, json, os, pathlib, signal, subprocess, time
import yaml
p=argparse.ArgumentParser();p.add_argument('--card',type=pathlib.Path,required=True);p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--binary',type=pathlib.Path,required=True);p.add_argument('--disabled-binary',type=pathlib.Path,required=True);a=p.parse_args()
root=a.root.resolve();root.mkdir();base=yaml.safe_load(a.card.read_text())
cases=[
 ('unknown-key',('mesh','moving_window','unknown'),1,'unknown'),
 ('negative-cells',('mesh','moving_window','cells'),-1,'positive'),
 ('zero-block',('mesh','moving_window','block_layers'),0,'positive'),
 ('missing-directory',('mesh','moving_window','reference_directory'),'','reference_directory'),
 ('filter',('mesh','current_filter'),'compensated-bilinear','current_filter'),
 ('z-cpml',('boundary','cells'),[10,10,2],'zero z CPML'),
 ('thin-xy-cpml',('boundary','cells'),[1,10,0],'at least two cells'),
 ('not-held',('beam','reference','frozen_lab_plane_injection'),False,'frozen_lab_plane_injection'),
 ('not-held-transverse',('beam','reference','hold_transverse_until_activation'),False,'hold_transverse'),
 ('checkpoint',('checkpoint',),{'steps':[1]},'checkpoint'),
 ('restart',('checkpoint',),{'restart_from':'unused'},'checkpoint'),
 ('wrong-initial-model',('initial_self_field','model'),'relativistic-poisson','cowan'),
 ('ballistic',('mesh','ballistic_outside_transverse'),True,'ballistic'),
 ('too-short-window',('mesh','moving_window','cells'),100,'cells must cover'),
 ('reference-too-short',('mesh','duration'),base['mesh']['duration']*3,'reference domain'),
 ('build-option-off',(),None,'without APRL_ENABLE_MOVING_WINDOW')]
results=[]
for name,path,value,needle in cases:
    directory=root/name;directory.mkdir();card=copy.deepcopy(base)
    card['mesh']['moving_window']['reference_directory']=str(directory/'unused-reference')
    node=card
    for key in path[:-1]:node=node[key]
    if path:node[path[-1]]=value
    (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    binary=a.disabled_binary if name=='build-option-off' else a.binary
    started=time.monotonic()
    with (directory/'run.log').open('w') as log:
        proc=subprocess.Popen(['mpiexec','-n','1','--bind-to','core',str(binary.resolve()),str(directory/'run.yaml')],cwd=directory,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:code=proc.wait(timeout=30)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGKILL);proc.wait();raise RuntimeError('Timed out: '+name)
    log=(directory/'run.log').read_text()
    ok=code!=0 and needle.lower() in log.lower() and not (directory/'unused-reference').exists()
    results.append({'case':name,'exit_code':code,'expected_message':needle,'pass':ok,'seconds':time.monotonic()-started})
    (root/'summary.json').write_text(json.dumps(results,indent=2))
    if not ok:raise RuntimeError('Unexpected validation result: '+name)
(root/'PASS.json').write_text(json.dumps({'status':'PASS','cases':len(results)},indent=2))
