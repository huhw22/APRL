#!/usr/bin/env python3
"""Plot completed moving-window energy and diagnostic-cost acceptance outputs.

Only reads the frozen campaign and writes a new analysis directory. No FFT,
windowing, near-field subtraction, residual offset removal, or PIC execution.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-energy-mpl-'+str(os.getuid()))

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter

LD = np.longdouble
MAGNETS = ('undulator', 'free-drift')
TITLES = {'undulator': 'One-period undulator', 'free-drift': 'Free drift; no undulator'}
VARIANTS = ('old-off', 'new-off', 'old-on', 'new-on')
VARIANT_LABELS = {'old-off': 'Old\nledger off', 'new-off': 'New\nledger off',
                  'old-on': 'Old\nledger on', 'new-on': 'New\nledger on'}
COLOURS = {'old-off': '#697386', 'new-off': '#2A788E', 'old-on': '#8B6A9C', 'new-on': '#C06A24'}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def unit_for(values):
    peak = max((float(np.max(np.abs(value))) for value in values if len(value)), default=0.)
    if peak >= 1.e-9:
        return 1.e-9, 'nJ'
    if peak >= 1.e-12:
        return 1.e-12, 'pJ'
    return 1.e-15, 'fJ'


def ledger(path):
    with h5py.File(path, 'r') as stream:
        group = stream['energy_ledger']
        require(int(group['complete'][()]) == 1, 'Ledger did not close at physical completion: '+str(path))
        count = int(group['committed_records'][()])
        require(count == len(group['records']), 'Ledger has uncommitted or truncated records')
        rows = group['records'][:count]
        frame = group.attrs['frame']
        if isinstance(frame, bytes):
            frame = frame.decode()
        require(frame == 'boosted computational frame', 'Unexpected energy-frame definition')
        require(int(group.attrs['format_version']) >= 4, 'Moving-volume schema required')
    require(len(rows) > 1 and int(rows['step'][0]) == 0 and np.all(np.diff(rows['step']) > 0), 'Invalid ledger step axis')
    def a(name):
        require(name in rows.dtype.names, 'Required ledger field absent: '+name)
        result = np.asarray(rows[name], dtype=LD)
        require(np.isfinite(result).all(), 'Nonfinite ledger field: '+name)
        return result
    kinetic = a('particle_kinetic_active_J')+a('particle_kinetic_removed_J')
    u = a('field_energy_interior_J')
    faces = a('outward_field_energy_face_J')
    raw = a('balance_residual_raw_J')
    boundary_work = a('boundary_source_work_J')
    remap = a('moving_remap_field_change_J')
    overwrite = a('inlet_overwrite_work_J')
    frozen = a('frozen_constraint_work_J')
    carrier = a('carrier_constraint_work_J')
    boundary_residual = raw-boundary_work-remap-overwrite
    full = a('balance_residual_J')
    longitudinal = a('longitudinal_boundary_outward_J')
    transverse = np.sum(faces[:,:4], axis=1)
    prescribed = a('prescribed_work_J')
    delta_k, delta_u = kinetic-kinetic[0], u-u[0]
    reconstructed = delta_k+delta_u+transverse+longitudinal-prescribed-boundary_work-frozen-carrier
    source_scale = abs(kinetic)+abs(kinetic[0])+abs(u)+abs(u[0])+np.sum(abs(faces),axis=1)
    source_scale += abs(prescribed)+abs(boundary_work)+abs(remap)+abs(overwrite)+abs(frozen)+abs(carrier)
    bound = 128*LD(np.finfo(np.float64).eps)*np.maximum(source_scale,LD(np.finfo(np.float64).tiny))
    require(np.all(abs(full-reconstructed)<=bound), 'Plotted energy terms do not reconstruct the recorded residual')
    require(np.all(abs(full-(boundary_residual-frozen-carrier))<=bound), 'Residual decomposition differs from ledger')
    exited = np.flatnonzero(rows['removed_macroparticles']>0)
    series = {'delta_particle_J':delta_k, 'delta_field_J':delta_u, 'transverse_outward_J':transverse,
        'longitudinal_outward_J':longitudinal, 'prescribed_work_J':prescribed, 'boundary_source_work_J':boundary_work,
        'frozen_constraint_work_J':frozen, 'carrier_constraint_work_J':carrier,
        'raw_residual_J':raw, 'boundary_corrected_residual_J':boundary_residual, 'fully_corrected_residual_J':full,
        'longitudinal_propagation_J':np.sum(faces[:,4:],axis=1), 'remap_outward_J':-remap,
        'overwrite_outward_J':-overwrite, 'remap_interface_correction_J':a('moving_remap_interface_correction_J')}
    return {'step':rows['step'], 'series':series, 'frame':frame,
        'first_exit_step':int(rows['step'][exited[0]]) if len(exited) else None,
        'closure_valid_values':np.unique(rows['closure_valid']).tolist(),
        'records':len(rows), 'final_step':int(rows['step'][-1]), 'source':str(path),
        'maximum_reconstruction_difference_J':float(np.max(abs(full-reconstructed))),
        'initial_particle_J':float(kinetic[0]), 'initial_field_J':float(u[0])}


def style_axis(axis):
    axis.axhline(0, color='0.45', linewidth=.65, zorder=0)
    axis.grid(axis='both', alpha=.2, linewidth=.6)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    axis.set_xlabel('Complete boost-frame step')
    axis.spines[['top','right']].set_visible(False)
    fmt = ScalarFormatter(useOffset=False)
    fmt.set_powerlimits((-3,4))
    axis.yaxis.set_major_formatter(fmt)


def save_figure(figure, output, stem):
    files=[]
    for suffix in ('png','pdf'):
        path=output/(stem+'.'+suffix)
        figure.savefig(path, dpi=180 if suffix=='png' else None, facecolor='white')
        files.append(path)
    plt.close(figure)
    return files


def energy_plot(data, root, output, late=False):
    figure,axes=plt.subplots(2,4,figsize=(19,9.5),squeeze=False)
    panels = [
        [('delta_particle_J',r'$\Delta K_{active+removed}$','#276BA4','-'),
         ('prescribed_work_J',r'$W_{external}$','#C06A24','--')],
        [('delta_field_J',r'$\Delta U_{field}$','#C06A24','-'),
         ('transverse_outward_J',r'$Q_{\perp,out}$','#32836C','--'),
         ('longitudinal_outward_J',r'$Q_{z,out}$ (combined)','#8756A0','-.')],
        [('frozen_constraint_work_J',r'$W_{frozen}$ (estimate)','#C06A24','--'),
         ('carrier_constraint_work_J',r'$W_{carrier}$ (estimate)','#8756A0','-.'),
         ('boundary_source_work_J',r'$W_{other\ boundary}$','#32836C',':')],
        [('raw_residual_J','Raw residual','#7D8491',':'),
         ('boundary_corrected_residual_J','After moving-boundary terms','#C06A24','--'),
         ('fully_corrected_residual_J','After constraint-supply estimates','#276BA4','-')]]
    titles=['Particle change and external-device work','Field change and boundary exchange',
        'Constrained-source supply estimates','Residual: no offset removal']
    final_notes=[]
    for r,magnet in enumerate(MAGNETS):
        item=data[magnet]
        mask=item['step']>=150 if late and magnet=='undulator' else np.ones(len(item['step']),dtype=bool)
        require(np.any(mask),'Requested late step region is empty')
        x=item['step'][mask]
        row_title=TITLES[magnet]
        if late:
            row_title+='; steps '+str(int(x[0]))+'–'+str(int(x[-1]))
        for c, specs in enumerate(panels):
            axis=axes[r,c]
            factor,unit=(1.e-15,'fJ') if late else unit_for([item['series'][key][mask] for key,_,_,_ in specs])
            for key,label,colour,style in specs:
                axis.plot(x,np.asarray(item['series'][key][mask]/factor,dtype=float),label=label,color=colour,linestyle=style,linewidth=1.35)
            if item['first_exit_step'] is not None and x[0]<=item['first_exit_step']<=x[-1]:
                axis.axvline(item['first_exit_step'],color='0.65',linestyle=':',linewidth=.7,zorder=0)
            style_axis(axis)
            axis.set_ylabel('Boost-frame energy ('+unit+')')
            axis.set_title(('('+chr(97+r*4+c)+') '+row_title+'\n'+titles[c]),fontsize=10)
            axis.legend(loc='best',fontsize=7.5,framealpha=.9)
        report=item['series']
        peak_index=int(np.argmax(np.abs(report['fully_corrected_residual_J'])))
        peak_factor,peak_unit=unit_for([report['fully_corrected_residual_J']])
        final_notes.append(TITLES[magnet]+': final R = '+format(float(report['fully_corrected_residual_J'][-1]/1.e-15),'+.5g')+
            ' fJ; full-course max |R| = '+format(float(abs(report['fully_corrected_residual_J'][peak_index])/peak_factor),'.5g')+
            ' '+peak_unit+' at step '+str(int(item['step'][peak_index]))+
            '; closure_valid = '+','.join(map(str,item['closure_valid_values'])))
    title='Late undulator stage and full free drift — unchanged cumulative energy reference' if late else 'Moving-window online energy budget — full history of both short cases'
    figure.suptitle(title,fontsize=15,y=.985)
    figure.text(.5,.103,'\n'.join(final_notes),ha='center',fontsize=8.5)
    figure.text(.5,.064,r'$R=\Delta K+\Delta U+Q_{\perp,out}+Q_{z,out}-W_{external}-W_{other\ boundary}-W_{frozen}-W_{carrier}$',
        ha='center',fontsize=11)
    region_note='Undulator view is restricted to steps >=150; the step-92 peak remains in energy_budget_comparison. No re-zeroing at step 150.' if late else 'Panels select units independently (pJ/fJ); no curve is normalized by its own peak. Full history includes the undulator step-92 residual peak.'
    figure.text(.5,.021,'Source: '+root.name+'/SCIENCE_REVIEW.json and science/*-new-on/ledger/energy-ledger.h5\n'+
        region_note+'\nBoost frame; no FFT/window; signed linear axes. Dotted vertical line: first particle exit. Constraint work remains approximate.',
        ha='center',fontsize=8,color='0.3')
    figure.subplots_adjust(left=.055,right=.985,bottom=.19,top=.91,hspace=.43,wspace=.31)
    return save_figure(figure,output,'late_energy_balance' if late else 'energy_budget_comparison')


def boundary_plot(data, root, output):
    figure,axes=plt.subplots(1,2,figsize=(12,4.7))
    specs=[('longitudinal_propagation_J','Propagation through z faces','#276BA4','-'),
        ('remap_outward_J','Moving-volume removal minus supply','#C06A24','--'),
        ('overwrite_outward_J','Negative inlet overwrite supply','#8756A0','-.'),
        ('longitudinal_outward_J','Combined longitudinal outward exchange','#222222','-')]
    for axis,magnet in zip(axes,MAGNETS):
        item=data[magnet];factor,unit=unit_for([item['series'][key] for key,_,_,_ in specs])
        for key,label,colour,style in specs:
            axis.plot(item['step'],np.asarray(item['series'][key]/factor,dtype=float),label=label,color=colour,
                linestyle=style,linewidth=1.7 if key=='longitudinal_outward_J' else 1.1)
        style_axis(axis);axis.set(title=TITLES[magnet],ylabel='Signed outward energy ('+unit+')')
        axis.legend(fontsize=8,loc='best')
    figure.suptitle('Longitudinal exchange: propagation + remap + inlet overwrite',fontsize=13,y=.98)
    figure.text(.5,.035,'Source: '+root.name+'/science/*-new-on/ledger/energy-ledger.h5\n'+
        'Boost frame, no FFT/window. Remap includes inserted/cropped cells and staggered-grid interface correction; negative values mean net supply.',
        ha='center',fontsize=8,color='0.3')
    figure.subplots_adjust(left=.075,right=.985,bottom=.19,top=.85,wspace=.27)
    return save_figure(figure,output,'longitudinal_boundary_exchange')


def performance_plot(review, root, output):
    figure,axes=plt.subplots(2,2,figsize=(13,8.5))
    comparisons=[('increment_on_new_vs_old','New / old\nledger on'),
        ('disabled_new_vs_old','New / old\nledger off'),
        ('new_total_diagnostic_on_vs_off','New on / off\nall diagnostics'),
        ('old_total_diagnostic_on_vs_off','Old on / off\nall diagnostics')]
    validation={}
    for r,magnet in enumerate(MAGNETS):
        data=review['performance'][magnet]; axis=axes[r,0]
        for i,variant in enumerate(VARIANTS):
            entry=data['configurations'][variant]
            values=np.asarray(entry['samples_s_per_step'],dtype=float)*1000
            require(len(values)>=3 and np.isfinite(values).all() and np.all(values>0),'Need at least three finite positive timing samples')
            middle=statistics.median(values.tolist())
            require(np.isclose(middle,entry['median_s_per_step']*1000,rtol=1e-14),'Stored timing median differs')
            spread=np.linspace(-.12,.12,len(values))
            axis.vlines(i,float(values.min()),float(values.max()),color=COLOURS[variant],linewidth=1.8,zorder=2)
            axis.scatter(i+spread,values,color=COLOURS[variant],s=25,alpha=.85,zorder=3)
            axis.plot([i-.18,i+.18],[middle,middle],color=COLOURS[variant],linewidth=2.5,zorder=4)
            axis.annotate(format(middle,'.4g'),(i,float(values.max())),xytext=(0,6),textcoords='offset points',ha='center',fontsize=8)
        axis.set_xticks(range(4),[VARIANT_LABELS[v] for v in VARIANTS])
        axis.set_ylabel('Step time (ms/step)');axis.set_title(TITLES[magnet]+'\nDots: repeats; bar: median; line: min–max',fontsize=10.5)
        # Explicit headroom: setting the lower limit disables autoscaling,
        # so a later margins() call cannot protect the top value labels.
        maximum=max(max(entry['samples_s_per_step']) for entry in data['configurations'].values())*1000
        axis.set_ylim(0,maximum*1.25)
        axis=axes[r,1]
        validation[magnet]={}
        for i,(name,label) in enumerate(comparisons):
            entry=data['comparisons'][name]
            values=np.asarray(entry['paired_percent_changes'],dtype=float)
            require(np.isfinite(values).all() and len(values)>=3,'Invalid paired changes')
            median=statistics.median(values.tolist())
            require(np.isclose(median,entry['median_paired_percent_change'],rtol=1e-13,atol=1e-12),'Stored paired median differs')
            axis.vlines(i,float(values.min()),float(values.max()),color='#276BA4',linewidth=1.8)
            axis.scatter(i+np.linspace(-.12,.12,len(values)),values,color='#276BA4',s=25,alpha=.85,zorder=3)
            axis.plot([i-.18,i+.18],[median,median],color='#16374F',linewidth=2.5,zorder=4)
            axis.annotate(format(median,'+.3g')+'%',(i,float(values.max())),xytext=(0,6),textcoords='offset points',ha='center',fontsize=8)
            validation[magnet][name]={'paired_percent_changes':values.tolist(),'median_paired_percent_change':median}
        axis.axhline(0,color='0.4',linewidth=.75)
        axis.set_xticks(range(4),[label for _,label in comparisons]);axis.set_ylabel('Paired step-time change (%)')
        axis.set_title('Same-repeat comparisons\nPositive = slower numerator configuration',fontsize=10.5)
        axis.margins(y=.2)
    for axis in axes.flat:
        axis.grid(axis='y',alpha=.2,linewidth=.6);axis.spines[['top','right']].set_visible(False)
        axis.tick_params(axis='x',labelsize=8)
    figure.suptitle('Energy-diagnostic cost — fixed 16 MPI × 8 threads, one-node comparisons',fontsize=14,y=.985)
    figure.text(.5,.035,'Source: '+root.name+'/PERFORMANCE_REVIEW.json; all recorded repeats retained, no profile.\n'+
        '40–240 step interval excludes preparation; initial reference preparation is reported separately in CSV. Small-case measurements are not production-overhead estimates.',
        ha='center',fontsize=8,color='0.3')
    figure.subplots_adjust(left=.075,right=.985,bottom=.14,top=.90,hspace=.48,wspace=.24)
    return save_figure(figure,output,'diagnostic_performance'),validation


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True,help='Frozen campaign containing completed SCIENCE_REVIEW.json')
    parser.add_argument('--output',type=Path,help='New analysis directory; defaults to ROOT/figures')
    parser.add_argument('--science-only',action='store_true',help='Create only completed science figures before performance is available')
    args=parser.parse_args();root=args.root.resolve();output=args.output.resolve() if args.output else root/'figures'
    require(root.is_dir(),'Campaign does not exist')
    require(not output.exists(),'Use a fresh figure directory; existing artifacts are not overwritten')
    science_path=root/'SCIENCE_REVIEW.json';plan_path=root/'plan.json'
    science=read_json(science_path);plan=read_json(plan_path)
    require(science['checks_pass'],'Science data/invariance/accounting checks have not passed')
    require(science['binary_sha256']==plan['binary_sha256'],'Science binary identity differs from frozen plan')
    sources=[science_path,plan_path,Path(__file__).resolve()]
    data={}
    for magnet in MAGNETS:
        entries=[x for x in plan['cases']['science'] if x['magnet']==magnet and x['variant']=='new-on']
        require(len(entries)==1,'Expected exactly one new-on science case per physical configuration')
        item=entries[0]
        require(science['cases'][item['name']]['checks_pass'],'Selected case has not passed structural checks')
        path=Path(item['directory'])/'ledger/energy-ledger.h5'
        require(path.resolve().is_relative_to(root),'Selected science ledger lies outside campaign')
        data[magnet]=ledger(path);sources.append(path)
    performance=None
    if not args.science_only:
        performance_path=root/'PERFORMANCE_REVIEW.json';performance=read_json(performance_path)
        require(performance['checks_pass'] and performance['binary_sha256']==plan['binary_sha256'],'Performance science/identity check not passed')
        sources.append(performance_path)
    output.mkdir(parents=True)
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,
        'axes.unicode_minus':True,'axes.formatter.useoffset':False})
    figures=energy_plot(data,root,output)+energy_plot(data,root,output,late=True)+boundary_plot(data,root,output)
    summary={'format':'aprl-moving-energy-figures-v1','campaign':str(root),'binary_sha256':plan['binary_sha256'],
        'sources_sha256':{str(p):sha(p) for p in sources},'cases':{},'science_only':args.science_only,
        'method':'Native boost-frame online energies; signed linear scales; delta K includes removed particles; no FFT, taper, subtraction, smoothing or forced residual offset.',
        'scope':'Plotting does not add physical-conservation acceptance or alter closure_valid.'}
    summary['late_view']={'undulator_first_step':150,'free_drift':'complete history',
        'energy_reference':'Original step-zero reference retained; no offset reset at the start of displayed region.',
        'full_history_figure':'energy_budget_comparison.png','late_figure':'late_energy_balance.png',
        'unit_selection':'Main panels select nJ/pJ/fJ independently; late view uses fJ throughout.'}
    for magnet,item in data.items():
        names=list(item['series'])
        with (output/(magnet+'_energy_series.csv')).open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['step']+names)
            writer.writerows([[int(step)]+[float(item['series'][name][i]) for name in names] for i,step in enumerate(item['step'])])
        summary['cases'][magnet]={k:v for k,v in item.items() if k not in ('step','series')}
        summary['cases'][magnet]['final_values_J']={name:float(values[-1]) for name,values in item['series'].items()}
    if performance:
        new_figures,metrics=performance_plot(performance,root,output);figures+=new_figures
        summary['performance']=metrics
        with (output/'performance_cases.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['case','magnet','repeat','variant','step_s_40_240','reference_preparation_s','launcher_s'])
            for item in plan['cases']['performance']:
                record=performance['cases'][item['name']]
                writer.writerow([item['name'],item['magnet'],item['repeat'],item['variant'],record['steps_40_240_seconds_per_step'],
                    record['reference_preparation_seconds'],record['launcher_seconds_including_preparation_and_close']])
    summary['figures_sha256']={str(path):sha(path) for path in figures}
    (output/'PLOT_REVIEW.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(output/'PLOT_REVIEW.json')


if __name__=='__main__':
    main()
