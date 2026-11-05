"""Small independent statistical and conversion checks; no PIC tracking.

Outputs two retained roundtrip fixtures and a compact CHECKS.json. Ensemble
masters are temporary and removed after their statistics have been inspected.
"""
import argparse
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile

import h5py
import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2

from generate import generate
from source_math import random_by_id,validate_config
from genesis import convert
from verify import C,QE,DTYPE,LD,need,recover_parameter_plane,verify


def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def configuration(n=4096,seed=20261106):
    wave=4.913172351123507e-7;T=12*wave/float(C);lo=2*wave/float(C);hi=10*wave/float(C)
    # Direct integral of 4 I_peak (t/T)(1-t/T), independent of sampling code.
    integrated_unit=2*T*((hi/T)**2-(lo/T)**2)-4*T/3*((hi/T)**3-(lo/T)**3)
    peak=n*float(QE)/integrated_unit
    return dict(schema='APRL_physical_electron_source_v1',electrons=n,count_mode='physical_fixed_n',
        seed=seed,wavelength_m=wave,chunk_records=257,compression=0,
        profile=dict(type='parabolic_time',full_duration_s=T,peak_current_A=peak,
            window_start_s=lo,window_stop_s=hi,bin_width_s=wave/float(C)),
        phase_space=dict(gamma_mean=298.2608843802298,gamma_sigma=298.2608843802298*.0002,
            sigma_x_m=132e-6,sigma_y_m=75e-6,emittance_n_x_m=2.5e-6,emittance_n_y_m=2.9e-6,
            alpha_x=.938,alpha_y=-.705,parameter_plane_z_m=.01,output_plane_z_m=0.))


def quiet_generate(config,path,chunk=None):
    with contextlib.redirect_stdout(io.StringIO()):return generate(config,path,chunk_records=chunk)


def expected_time_mode(config,harmonic):
    p=config['profile'];T=p['full_duration_s'];lo=p['window_start_s']/T;hi=p['window_stop_s']/T
    norm=(3*hi**2-2*hi**3)-(3*lo**2-2*lo**3)
    frequency=2*np.pi*float(C)/config['wavelength_m']*harmonic*T
    # Fourier quadrature of an analytic polynomial, independent of inverse CDF.
    real=quad(lambda x:6*x*(1-x),lo,hi,weight='cos',wvar=frequency,epsabs=1e-13)[0]/norm
    imag=quad(lambda x:6*x*(1-x),lo,hi,weight='sin',wvar=frequency,epsabs=1e-13)[0]/norm
    return complex(real,imag)


def mode_covariance(config,width=None):
    p=config['phase_space'];a1=a2=1.
    if width is not None:
        sx=p['sigma_x_m'];sy=p['sigma_y_m']
        a1=1/np.sqrt((1+sx*sx/(2*width*width))*(1+sy*sy/(2*width*width)))
        a2=1/np.sqrt((1+sx*sx/(width*width))*(1+sy*sy/(width*width)))
    mean=a1*expected_time_mode(config,1);second=a2*expected_time_mode(config,2);power=a2
    covariance=np.array([[(power+second.real)/2-mean.real**2,second.imag/2-mean.real*mean.imag],
        [second.imag/2-mean.real*mean.imag,(power-second.real)/2-mean.imag**2]])
    return mean,covariance


def ensemble_check(output,draws=128,n=2048):
    config=configuration(n);modes={'plain':None,'gaussian_100um':100e-6}
    values={key:[] for key in modes};counts=[];probabilities=None;all_counts={}
    with tempfile.TemporaryDirectory(prefix='source-ensemble-',dir=output) as temp:
        for index in range(draws):
            cfg=copy.deepcopy(config);cfg['seed']=20270000+index
            path=Path(temp)/'sample.h5';quiet_generate(cfg,path)
            with h5py.File(path,'r') as f:
                a=f['particles/records'][:];counts.append(f['generation/bin_counts'][:])
                if probabilities is None:
                    edges=f['generation/bin_edges_time_s'][:];T=cfg['profile']['full_duration_s'];r=edges/T
                    integrals=2*r*r-(4/3)*r*r*r
                    probabilities=np.diff(integrals)/np.diff(integrals).sum()
            variables,time,_=recover_parameter_plane(a,cfg['phase_space'])
            exponential=np.exp(2j*np.pi*np.asarray(time,dtype=float)*float(C)/cfg['wavelength_m'])
            for key,width in modes.items():
                weight=1. if width is None else np.exp(-np.asarray(variables[:,0]**2+variables[:,1]**2,float)/(4*width*width))
                values[key].append(np.mean(weight*exponential))
            path.unlink()
    results={}
    for key,width in modes.items():
        mean,covariance=mode_covariance(config,width);b=np.asarray(values[key])
        centered=np.column_stack([(b.real-mean.real),(b.imag-mean.imag)])*np.sqrt(n)
        statistic=float(np.einsum('ni,ij,nj->',centered,np.linalg.inv(covariance),centered))
        interval=chi2.ppf([1e-5,1-1e-5],2*draws)
        need(interval[0]<statistic<interval[1],'Physical mode-noise variance rejected: '+key)
        whitened=np.linalg.solve(np.linalg.cholesky(covariance),centered.T).T
        mean_z=np.mean(whitened,axis=0)*np.sqrt(draws)
        need(np.max(abs(mean_z))<5.5,'Mode ensemble mean differs from analytic finite envelope: '+key)
        results[key]=dict(samples=draws,particles_each=n,analytic_mean=[mean.real,mean.imag],
            analytic_single_particle_covariance=covariance.tolist(),
            aggregate_mahalanobis_chi_square=statistic,chi_square_degrees_of_freedom=2*draws,
            acceptance_interval_approximate_CLT=interval.tolist(),whitened_mean_standard_errors=mean_z.tolist(),
            mean_centered_power_times_N=float(n*np.mean(abs(b-mean)**2)),
            theoretical_centered_power_times_N=float(np.trace(covariance)))
    counts=np.asarray(counts,dtype=float);expect=n*probabilities
    cov_expected=n*(np.diag(probabilities)-np.outer(probabilities,probabilities))
    observed=np.cov(counts,rowvar=False,ddof=1)
    z=(counts.mean(axis=0)-expect)/np.sqrt(np.diag(cov_expected)/draws)
    variance_ratio=np.diag(observed)/np.diag(cov_expected)
    bounds=chi2.ppf([1e-5,1-1e-5],draws-1)/(draws-1)
    need(np.max(abs(z))<5.5,'Multinomial mean failed')
    need(np.all((variance_ratio>bounds[0])&(variance_ratio<bounds[1])),'Multinomial diagonal covariance failed')
    covariance_se=np.sqrt((cov_expected**2+np.outer(np.diag(cov_expected),np.diag(cov_expected)))/(draws-1))
    covariance_z=(observed-cov_expected)/covariance_se
    need(np.max(abs(covariance_z))<6,'Multinomial cross-bin covariance failed')
    need(np.all(counts.sum(axis=1)==n),'Total electron count fluctuated in fixed-N mode')
    results['fixed_N_counts']=dict(analytic_mean=expect.tolist(),analytic_covariance=cov_expected.tolist(),
        measured_covariance=observed.tolist(),mean_standard_errors=z.tolist(),
        variance_ratios=variance_ratio.tolist(),approximate_variance_ratio_acceptance=bounds.tolist(),
        maximum_covariance_approximate_standard_errors=float(np.max(abs(covariance_z))),
        off_diagonal_covariance_max_abs=float(np.max(abs(observed-cov_expected-np.diag(np.diag(observed-cov_expected))))),
        total_N_exact=True)
    save(output/'ENSEMBLE.json',results)
    return results


def fixture(path):
    period=4.913172351123507e-7/float(C)
    times=np.array([-3.2,-3.1,-1.,np.nextafter(0.,-1.),0.,np.nextafter(1.,0.),1.,3.])*period
    a=np.zeros(len(times),dtype=DTYPE);a['source_id']=np.arange(1,len(times)+1);a['macro_weight']=1.
    a['arrival_time_offset_s']=times;a['plane_position_m'][:,0]=np.arange(len(times))*1e-6
    a['plane_position_m'][:,1]=-np.arange(len(times))*2e-6
    a['proper_velocity'][:,0]=np.linspace(-.01,.02,len(times));a['proper_velocity'][:,1]=np.linspace(.03,-.01,len(times));a['proper_velocity'][:,2]=10+np.arange(len(times))*.1
    with h5py.File(path,'x') as f:
        g=f.create_group('particles');g.attrs['format_version']=4;g.attrs['input_plane_z_m']=0.;g.attrs['electron_count']=len(a)
        g.attrs['wavelength_m']=4.913172351123507e-7;g.attrs['charge_C']=len(a)*float(QE)
        g.create_dataset('records',data=a);m=f.create_group('generation');m['complete']=np.uint8(1);m['committed_records']=np.uint64(len(a))
        m.attrs['record_payload_sha256']=hashlib.sha256(a.tobytes()).hexdigest()
    return a


def expect_failure(name,callback):
    try:callback()
    except (ValueError,RuntimeError,OSError,KeyError):return dict(passed=True,rejected=name)
    raise AssertionError('Negative control accepted: '+name)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    parser.add_argument('--ensemble-draws',type=int,default=128);args=parser.parse_args()
    out=args.output.resolve();need(not out.exists(),'Refusing existing check directory');out.mkdir(parents=True)
    config=configuration();save(out/'roundtrip-config.json',config)
    quiet_generate(config,out/'roundtrip-master.h5',37)
    quiet_generate(config,out/'chunk-comparison.h5',501)
    with h5py.File(out/'roundtrip-master.h5','r') as f,h5py.File(out/'chunk-comparison.h5','r') as g:
        a=f['particles/records'][:];b=g['particles/records'][:]
        need(np.array_equal(a,b),'Source records changed with chunk size')
        need(np.array_equal(f['generation/bin_counts'][:],g['generation/bin_counts'][:]),'Counts changed with chunk size')
    need(len(np.unique(a['plane_position_m'],axis=0))==len(a),'A transverse position was copied into a rod')
    whole_time,whole_normal=random_by_id(config['seed'],0,128)
    recovered_time=np.empty_like(whole_time);recovered_normal=np.empty_like(whole_normal)
    for start,stop in [(77,128),(11,77),(0,11)]:
        t,z=random_by_id(config['seed'],start,stop-start);recovered_time[start:stop]=t;recovered_normal[start:stop]=z
    need(np.array_equal(whole_time,recovered_time) and np.array_equal(whole_normal,recovered_normal),
         'Random access by source ID changed with processing order')
    # Statistical independence across long source-ID lags, rather than family labels.
    variables,_,_=recover_parameter_plane(a,config['phase_space']);cross=[]
    for lag in [1,4,16,64]:
        q=np.asarray(variables[:,0],float);value=float(np.corrcoef(q[:-lag],q[lag:])[0,1]);cross.append(value)
        need(abs(value)<.12,'Source-ID correlation exceeds the independent-sample check')
    convert(out/'roundtrip-master.h5',out/'roundtrip-genesis.h5',config['wavelength_m'],pad_slices=4,ranks_multiple=16,chunk_records=257)
    normal=verify(out/'roundtrip-master.h5',out/'roundtrip-genesis.h5',chunk_records=193)
    save(out/'ROUNDTRIP.json',normal)
    cold=configuration(512);cold['phase_space']['gamma_sigma']=0.;cold['phase_space']['emittance_n_y_m']=0.
    quiet_generate(cold,out/'cold-master.h5',29)
    cold_report=verify(out/'cold-master.h5',chunk_records=31)
    save(out/'COLD.json',cold_report)
    fixture(out/'edge-master.h5')
    convert(out/'edge-master.h5',out/'edge-genesis.h5',config['wavelength_m'],pad_slices=4,ranks_multiple=16,chunk_records=3)
    edge=verify(out/'edge-master.h5',out/'edge-genesis.h5',chunk_records=3,check_distribution=False)
    need(edge['genesis']['empty_slices']>0,'Fixture did not exercise empty slices')
    save(out/'EDGE_ROUNDTRIP.json',edge)
    negatives={}
    for name,edit in [('incomplete',lambda f:f['generation/complete'].__setitem__((),0)),
                      ('weight_not_one',lambda f:f['particles/records'].__setitem__(0,mutated_record(f,'macro_weight',2.))),
                      ('nonfinite_time',lambda f:f['particles/records'].__setitem__(0,mutated_record(f,'arrival_time_offset_s',np.nan))),
                      ('noncontiguous_ID',lambda f:f['particles/records'].__setitem__(0,mutated_record(f,'source_id',4)))]:
        path=out/(name+'.h5');shutil.copy2(out/'edge-master.h5',path)
        with h5py.File(path,'r+') as f:edit(f)
        negatives[name]=expect_failure(name,lambda p=path:convert(p,out/(name+'-unexpected.h5'),config['wavelength_m']))
    negatives['overwrite']=expect_failure('overwrite',lambda:quiet_generate(config,out/'roundtrip-master.h5'))
    invalid=copy.deepcopy(config);invalid['electrons']+=1
    negatives['physical_count_mismatch']=expect_failure('physical_count_mismatch',lambda:validate_config(invalid))
    path=out/'finite-time-outside-profile.h5';shutil.copy2(out/'roundtrip-master.h5',path)
    with h5py.File(path,'r+') as f:
        a=f['particles/records'][0];a['arrival_time_offset_s']+=1e-3;f['particles/records'][0]=a
    negatives['finite_time_outside_profile']=expect_failure('finite_time_outside_profile',lambda:verify(path))
    invalid=copy.deepcopy(config);invalid['phase_space']['gamma_mean']=1.00001;invalid['phase_space']['gamma_sigma']=0.;invalid['phase_space']['emittance_n_x_m']=.001
    negatives['unphysical_shell']=expect_failure('unphysical_shell',lambda:quiet_generate(invalid,out/'unphysical-shell.h5'))
    need(not (out/'unphysical-shell.h5').exists(),'Failed shell source was published')
    ensemble=ensemble_check(out,args.ensemble_draws)
    result=dict(passed=True,normal_master=str(out/'roundtrip-master.h5'),normal_genesis=str(out/'roundtrip-genesis.h5'),N=config['electrons'],
        fixtures=dict(edge_master=str(out/'edge-master.h5'),edge_genesis=str(out/'edge-genesis.h5')),
        tests=dict(chunk_size_bitwise=True,bin_count_chunk_size_bitwise=True,out_of_order_ID_access_bitwise=True,unique_transverse_positions=True,
            source_ID_lag_correlations=cross,analytic_distribution=normal['master']['distribution'],
            deterministic_master_genesis=normal['genesis'],edge_empty_variableN=edge['genesis'],
            degenerate_cold_variables=True,negative_controls=negatives,independent_multiseed_statistics=True),
        limitations=['Small analytic tests are not FEL physics acceptance','Native PIC reader and actual Genesis import/write are tested separately by the parent runner'])
    save(out/'CHECKS.json',result);print(json.dumps(dict(passed=True,N=config['electrons'],output=str(out)),indent=2))


def mutated_record(f,key,value):
    a=f['particles/records'][0];a[key]=value;return a


if __name__=='__main__':main()
