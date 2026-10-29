#!/usr/bin/env python3
"""Organize completed successful-test evidence, retain old paths, build Markdown.

Only the completed discrete_energy_ledger_20261029 result tree is relocated.
The original tree is replaced by a compatibility symlink. No simulation runs.
"""
import argparse,csv,hashlib,json,os,re,shutil
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-test-record-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
NAME='discrete_energy_ledger_20261029'

def save(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+'\n')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def inventory(root):
    rows=[]
    for parent,dirs,files in os.walk(root,followlinks=False):
        for name in dirs+files:
            p=Path(parent)/name;s=p.lstat()
            item={'relative':str(p.relative_to(root)),'device':s.st_dev,'inode':s.st_ino,'size':s.st_size,'mtime_ns':s.st_mtime_ns}
            item['kind']='link' if p.is_symlink() else ('directory' if p.is_dir() else 'file')
            if p.is_symlink():
                item['link_target']=os.readlink(p);item['target_existed']=p.exists()
                if not os.path.isabs(item['link_target']):raise RuntimeError('Inspect relative link before moving: '+str(p))
            elif p.is_file() and p.suffix in ('.json','.yaml','.yml','.txt','.log','.py','.cpp','.png','.pdf','.csv','.complete'):
                item['sha256']=sha(p)
            rows.append(item)
    return rows

def summary(root,stage,mode):
    directory=root/('campaign-r2' if stage!='two' else 'mpi-r1')/'execution'
    name=mode+('-1rank' if stage=='one' else ('-4rank' if stage=='four' else ''))
    p=directory/name/'summary.json';d=json.loads(p.read_text());assert d['pass'];return d

def log_bars(ax,series,labels,title,ylabel):
    colors=['#4b78a4','#e5a04a','#639879'];x=np.arange(3);v=np.array(series)
    assert np.all(v>0)
    ax.bar(x,v,color=colors);ax.set_yscale('log');ax.set_xticks(x);ax.set_xticklabels(labels,fontsize=8)
    ax.set_title(title,fontsize=11);ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.2)
    ax.set_ylim(v.min()/4,v.max()*8)
    for i,y in enumerate(v):ax.annotate(f'{y:.2e}',(i,y),xytext=(0,6),textcoords='offset points',ha='center',fontsize=9)

def make_plots(root,output):
    labels=['1 rank\n1 node','4 ranks\n1 node','4 ranks\n2 nodes'];stages=['one','four','two'];data={}
    pulse=[summary(root,s,'pulse') for s in stages];data['pulse']=pulse
    fig,axes=plt.subplots(1,3,figsize=(12,4.3),layout='constrained')
    for ax,key,title in zip(axes,['maximum_translation_relative_error','maximum_analytic_flux_relative_error','max_normalized_operator_residual'],['Known field translation','Native face exchange','Field operator accounting']):
        log_bars(ax,[d[key] for d in pulse],labels,title,'Relative / normalized difference')
    fig.savefig(output/'01_vacuum_pulse.png',dpi=180);fig.savefig(output/'01_vacuum_pulse.pdf');plt.close(fig)
    manufactured=[summary(root,s,'manufactured') for s in stages];data['manufactured']=manufactured
    fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    log_bars(axes[0],[d['max_normalized_operator_residual'] for d in manufactured],labels,'3-D field operator accounting','Normalized difference')
    log_bars(axes[1],[d['max_electric_prediction_relative_error'] for d in manufactured],labels,'Electric increment prediction','Relative difference')
    fig.savefig(output/'02_manufactured_fields.png',dpi=180);fig.savefig(output/'02_manufactured_fields.pdf');plt.close(fig)
    injected=[summary(root,s,'injected') for s in stages];data['injected']=injected
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained');x=np.arange(3)
    for ax,key,title in zip(axes,['max_normalized_operator_residual','max_electric_prediction_relative_error'],['Accounting detects the perturbation','Pointwise prediction detects the perturbation']):
        for dx,ds,label,color in [(-.18,manufactured,'Expected numerical state','#4b78a4'),(.18,injected,'Known injected E perturbation','#d88639')]:
            ax.bar(x+dx,[d[key] for d in ds],width=.36,label=label,color=color)
        ax.set_yscale('log');ax.set_xticks(x);ax.set_xticklabels(labels);ax.set_title(title,fontsize=10)
        ax.set_ylabel('Relative / normalized difference');ax.grid(axis='y',alpha=.2)
        ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,1.24),frameon=False)
    fig.savefig(output/'03_perturbation_detection.png',dpi=180);fig.savefig(output/'03_perturbation_detection.pdf');plt.close(fig)
    data['particle_series']={};fig,axes=plt.subplots(2,1,figsize=(10,7.5),sharex=True,layout='constrained')
    for name,style,color in [('single','-','#4b78a4'),('mpi','--','#d88639')]:
        path=root/f'report-{name}-final/discrete_energy.csv'
        with path.open() as f:rows=list(csv.DictReader(f))
        fields=['field_operator_residual_J','discrete_total_residual_J','particle_grid_coupling_gap_J']
        selected=[{'step':int(d['step']),**{n:d[n] for n in fields}} for d in rows];data['particle_series'][name]=selected
        steps=[d['step'] for d in selected]
        axes[0].plot(steps,[float(d['field_operator_residual_J'])/1e-31 for d in selected],style,label='Single node' if name=='single' else 'Two nodes',color=color)
        axes[1].plot(steps,[float(d['discrete_total_residual_J'])/1e-20 for d in selected],style,label='Total: single node' if name=='single' else 'Total: two nodes',color=color)
        if name=='single':axes[1].plot(steps,[float(d['particle_grid_coupling_gap_J'])/1e-20 for d in selected],':',color='#639879',label='Particle-grid exchange: single node')
    axes[0].set_title('Particle free drift: completed discrete-accounting checks')
    axes[0].set_ylabel(r'Cumulative field check ($10^{-31}$ J)');axes[1].set_ylabel(r'Cumulative total / exchange check ($10^{-20}$ J)');axes[1].set_xlabel('Maxwell step')
    axes[1].set_xlim(0,257);axes[1].set_xticks([0,50,100,150,200,257])
    for ax in axes:ax.axhline(0,color='#999999',linewidth=.7);ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.savefig(output/'04_particle_free_drift.png',dpi=180);fig.savefig(output/'04_particle_free_drift.pdf');plt.close(fig)
    final=json.loads((root/'FINAL.json').read_text());perf=final['performance'];memory={}
    for run in ['campaign-r2','mpi-r1']:
        memory[run]={}
        for flag in ['off','on']:
            text=(root/run/flag/'aprl.log').read_text();values=[float(x) for x in re.findall(r'current_rss_total_mib=([0-9.]+)',text)]
            memory[run][flag]=max(values)
    data['performance']=perf;data['observed_total_RSS_MiB']=memory
    fig,axes=plt.subplots(1,2,figsize=(10,4.8),layout='constrained');x=np.arange(2)
    for dx,flag,color,label in [(-.18,'off','#4b78a4','Audit disabled'),(.18,'on','#d88639','Audit enabled')]:
        t=[perf[n][flag+'_seconds_per_step'] for n in ['single_node','two_node']];m=[memory[n][flag] for n in ['campaign-r2','mpi-r1']]
        for ax,v,fmt in [(axes[0],t,'.3f'),(axes[1],m,'.0f')]:
            ax.bar(x+dx,v,width=.36,label=label,color=color)
            for i,y in enumerate(v):ax.annotate(format(y,fmt),(i+dx,y),xytext=(0,5),textcoords='offset points',ha='center',fontsize=9)
    axes[0].set_title('Measured step time (40-240)');axes[0].set_ylabel('Seconds / step');axes[1].set_title('Maximum observed concurrent RSS');axes[1].set_ylabel('Total RSS (MiB)')
    for ax in axes:ax.set_xticks(x);ax.set_xticklabels(['Single node','Two nodes']);ax.set_ylim(0,ax.get_ylim()[1]*1.2);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=9)
    fig.savefig(output/'05_audit_cost.png',dpi=180);fig.savefig(output/'05_audit_cost.pdf');plt.close(fig)
    save(output/'PLOT_DATA.json',data)

def render_documents(root,target):
    """Refresh authored documents without modifying archived test evidence."""
    final=json.loads((root/'FINAL.json').read_text())
    report=json.loads((root/'report-single-final/summary.json').read_text());perf=final['performance']
    values={'FIELD_MAX':f"{final['max_cumulative_field_operator_residual_J']:.3g}",
      'TOTAL_ABS':f"{abs(report['final']['discrete_total_residual_J']):.3g}",
      'TOTAL_PERCENT':f"{100*abs(final['new_final_over_abs_net_particle_change']):.5f}"}
    for label,key in [('SINGLE','single_node'),('MPI','two_node')]:
        values[label+'_OFF']=f"{perf[key]['off_seconds_per_step']:.5f}";values[label+'_ON']=f"{perf[key]['on_seconds_per_step']:.5f}";values[label+'_COST']=f"{perf[key]['increment_percent']:.2f}"
    for source,name in [('report_template.md','README.md'),('software_validation_template.md','PROGRAM_VALIDATION.md')]:
        template=(HERE/source).read_text()
        for key,value in values.items():template=template.replace('{{'+key+'}}',value)
        assert '{{' not in template and '旧差额' not in template
        (target/name).write_text(template)
    save(target/'metadata/SUCCESS_METRICS.json',{'source':str(root/'FINAL.json'),'source_sha256':sha(root/'FINAL.json'),'values':values,'test_program_commit':'c7fcd20','new_simulations':0})
    save(target/'metadata/RENDER_SHA256.json',{
      str(q.relative_to(target)):sha(q) for q in [target/'README.md',target/'PROGRAM_VALIDATION.md',*(target/'figures').iterdir()] if q.is_file()})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',type=Path,required=True);p.add_argument('--target',type=Path,required=True);p.add_argument('--confirmed-idle',action='store_true',required=True)
    a=p.parse_args();ws=a.workspace.resolve();target=a.target.absolute();old=ws/'results'/NAME
    assert old.is_dir() and not old.is_symlink() and not target.exists()
    final=json.loads((old/'FINAL.json').read_text());assert final['checks_pass']
    for stage in ['campaign-r2','mpi-r1']:assert json.loads((old/stage/'COMPLETION.json').read_text())['checks_pass']
    assert old.stat().st_dev==target.parent.stat().st_dev
    rows=inventory(old);identity=old.stat();target.mkdir()
    for name in ['evidence','figures','tests','scripts','programs','metadata']:(target/name).mkdir()
    new=target/'evidence'/NAME
    save(target/'metadata/ORIGINAL_INVENTORY.json',{'original_root':str(old),'files':rows})
    save(target/'metadata/PLAN.json',{'operation':'same-filesystem rename + old-path compatibility symlink','original':str(old),'destination':str(new),'scope':'Completed discrete-ledger validation; scientific raw outputs retained byte-for-byte','idle_confirmed_before_operation':True})
    old.rename(new)
    try:old.symlink_to(os.path.relpath(new,old.parent),target_is_directory=True)
    except BaseException:
        if not os.path.lexists(old):new.rename(old)
        raise
    assert new.stat().st_ino==identity.st_ino and old.resolve()==new.resolve()
    save(target/'metadata/MIGRATION.json',{'original':str(old),'destination':str(new),'compatibility_link':os.readlink(old),'root_device':identity.st_dev,'root_inode':identity.st_ino,'same_filesystem':True,'no_scientific_files_deleted':True})
    for label,child in [('01_single_node','campaign-r2'),('02_two_nodes','mpi-r1'),('03_reader_checks','reader-fixtures')]:
        (target/'tests'/label).symlink_to(os.path.relpath(new/child,target/'tests'),target_is_directory=True)
    source=ws/'source/aprl/verification/performance/cases/discrete_energy_ledger'
    shutil.copytree(source,target/'scripts/discrete_energy_ledger',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    programs={'aprl':Path(final['binary']),'operator_checks':ws/'build/discrete-energy-ledger-20261029/operator_checks'}
    index={}
    for name,path in programs.items():
        (target/'programs'/name).symlink_to(os.path.relpath(path,target/'programs'))
        index[name]={'target':str(path),'sha256':sha(path)}
    save(target/'programs/INDEX.json',index)
    make_plots(new,target/'figures')
    render_documents(new,target)
    save(target/'metadata/SNAPSHOT_SHA256.json',{str(q.relative_to(target)):sha(q) for q in (target/'scripts').rglob('*') if q.is_file()})
    shutil.copy2(HERE/'validate.py',target/'metadata/validate.py')
    import importlib.util
    spec=importlib.util.spec_from_file_location('validate_bundle',HERE/'validate.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.validate(target)
    print(target/'README.md')
if __name__=='__main__':main()
