#!/usr/bin/env python3
"""Predeclared source-scale roundoff gates for the complete cohort energy ledger.

No gate is fitted to candidate results. Axial chirp/uncorrelated spread use an
explicit affine-z diagnostic approximation and are reported without an equality gate.
"""
import h5py
import numpy as np

LD=np.longdouble
EPS64=LD(np.finfo(np.float64).eps)
ULD=LD(np.finfo(np.longdouble).eps)/2


def thresholds(count):
    g=LD(count)*ULD/(1-LD(count)*ULD)
    return {'source_records':count,'eps64':float(EPS64),'longdouble_unit_roundoff':float(ULD),
        'longdouble_mantissa_bits':int(np.finfo(np.longdouble).nmant+1),
        'gamma_N_longdouble':float(g),'source_sum_relative_bound':float(128*EPS64+2*g),
        'source_sum_derivation':'128eps64 plus two gamma_N(u_longdouble) naive positive-sum bounds; N is fixed initial cohort, not measured discrepancies.',
        'integer_time_contract':'All integer fields and time_box_s exact; active_macroparticles includes pending; active+removed==initial cohort every sample.',
        'primitive_scale':'Positive energy/weight on own source scale; gamma means on positive gamma scale; signed work on initial/current beam+field energy; face flux on maximum field/face energy; boundary work on absolute boundary work.',
        'residual_bound':'Sum of primitive absolute bounds for K(t)+Kremoved(t)+F(t)-K(0)-F(0)+out-work[-boundary work], plus export/arithmetic roundoff.',
        'ratio_bound':'Numerator/denominator absolute-bound propagation when denominator resolved; otherwise cross-multiplied residual gate with denominator uncertainty, explicitly marked ill-conditioned.',
        'sigma_bound':'Compare variance, bound=(5e-15+2r)*reference_variance+mean_gamma_bound^2+(64eps64*max(1,mean_gamma))^2. Last term covers existing zero-resolution clamp.',
        'axial_report_only':['linear_chirp_gamma_per_m','uncorrelated_sigma_gamma_lab'],
        'axial_scope':'Pending z statistics use long-double affine diagnostic trajectories; exact repeated-addition coordinates are used at admission. No claim of bitwise axial diagnostic equivalence.',
        'field_relative_l2_unchanged':2e-12,'particle_x_u_t_eps64_factor_unchanged':128}


def require(test,message):
    if not test: raise RuntimeError(message)


def checked(a,b,bound,label,scale=None):
    x,y=np.asarray(a,dtype=LD),np.asarray(b,dtype=LD)
    d=abs(y-x); t=np.broadcast_to(np.asarray(bound,dtype=LD),d.shape)
    bad=d>t
    result={'pass':not bool(np.any(bad)),'gate':label,'max_abs_difference':float(np.max(d)),
        'rms_difference':float(np.sqrt(np.mean(d*d))),'max_absolute_bound':float(np.max(t)),
        'min_absolute_bound':float(np.min(t)),'elements_over_bound':int(np.count_nonzero(bad)),
        'bitwise_equal':a.dtype==b.dtype and a.tobytes()==b.tobytes()}
    if scale is not None:
        s=np.broadcast_to(np.asarray(scale,dtype=LD),d.shape)
        ratio=np.divide(d,s,out=np.zeros_like(d),where=s>0)
        result['max_error_over_source_scale']=float(np.max(ratio))
        result['zero_source_nonzero_error_elements']=int(np.count_nonzero((s==0)&(d!=0)))
    return result


def compare(reference,candidate,count):
    with h5py.File(reference/'ledger/energy-ledger.h5') as f: a=f['energy_ledger/records'][:]
    with h5py.File(candidate/'ledger/energy-ledger.h5') as f: b=f['energy_ledger/records'][:]
    require(a.dtype==b.dtype and a.shape==b.shape and len(a)>0,'Ledger schemas/sample counts differ')
    require(all(np.isfinite(a[n]).all() and np.isfinite(b[n]).all() for n in a.dtype.names),'Nonfinite ledger value')
    contract=thresholds(count);r=LD(contract['source_sum_relative_bound'])
    report={'pass':True,'thresholds':contract,'fields':{},'count_includes_pending':True}
    bound={};handled=set()
    def gate(name,t,label,scale=None):
        report['fields'][name]=checked(a[name],b[name],t,label,scale)
        bound[name]=np.asarray(t,dtype=LD);handled.add(name)
        report['pass'] &= report['fields'][name]['pass']
    for name in a.dtype.names:
        if a[name].dtype.kind not in 'fc' or name=='time_box_s':
            passed=np.array_equal(a[name],b[name]);report['fields'][name]={'pass':bool(passed),'gate':'exact','bitwise_equal':bool(passed)}
            handled.add(name);report['pass'] &= bool(passed)
    for records,label in ((a,'reference'),(b,'candidate')):
        require(int(records['active_macroparticles'][0])==count,'Initial ledger omits source/pending particles: '+label)
        require(np.all(records['active_macroparticles']+records['removed_macroparticles']==count),'Ledger cohort count not conserved: '+label)
    positive=('active_represented_electrons','particle_kinetic_active_J','particle_total_active_J',
              'particle_kinetic_removed_J','particle_total_removed_J','field_energy_interior_J')
    for name in positive:
        require(np.all(a[name]>=0) and np.all(b[name]>=0),'Negative positive ledger primitive: '+name)
        x=np.asarray(a[name],dtype=LD);gate(name,r*abs(x),'positive source sum',abs(x))
    ka=np.asarray(a['particle_kinetic_active_J'],dtype=LD);kr=np.asarray(a['particle_kinetic_removed_J'],dtype=LD)
    fe=np.asarray(a['field_energy_interior_J'],dtype=LD)
    dynamic=abs(ka)+abs(kr)+abs(fe)+abs(ka[0])+abs(fe[0])
    face='outward_field_energy_face_J';faces=np.asarray(a[face],dtype=LD)
    field_scale=np.maximum(np.max(abs(fe)),np.max(abs(faces),axis=1))
    gate(face,r*field_scale[:,None],'signed face flux uses field/face source scale',field_scale[:,None])
    outward_scale=np.maximum(np.sum(abs(faces),axis=1),field_scale)
    gate('outward_field_energy_J',6*r*outward_scale,'sum of six signed face flux bounds',outward_scale)
    gate('prescribed_work_J',r*np.maximum(abs(a['prescribed_work_J']),dynamic),'signed work uses beam/field source energy',dynamic)
    absolute=np.asarray(a['boundary_source_abs_work_J'],dtype=LD)
    gate('boundary_source_abs_work_J',r*abs(absolute),'absolute source work',abs(absolute))
    gate('boundary_source_work_J',r*abs(absolute),'signed source work on absolute work scale',abs(absolute))
    mean_bounds={};sigma_bounds={}
    for frame in ('lab','box'):
        mean='mean_gamma_'+frame;mu=np.asarray(a[mean],dtype=LD)
        source=np.maximum(abs(mu),np.where(a['active_macroparticles']>0,1,0))
        mean_bounds[frame]=2*r*source
        gate(mean,mean_bounds[frame],'positive weighted gamma numerator/weight propagation',source)
        sigma='sigma_gamma_'+frame
        require(np.all(a[sigma]>=0) and np.all(b[sigma]>=0),'Negative gamma spread')
        s=np.asarray(a[sigma],dtype=LD);variance=s*s
        resolution=64*EPS64*source
        bv=(LD('5e-15')+2*r)*variance+mean_bounds[frame]**2+resolution**2
        sigma_bounds[frame]=np.sqrt(variance+bv)-s
        item=checked(np.asarray(a[sigma],dtype=LD)**2,np.asarray(b[sigma],dtype=LD)**2,bv,
                     'central variance + mean uncertainty + existing zero clamp',variance+resolution**2)
        item['max_sigma_abs_difference']=float(np.max(abs(np.asarray(b[sigma],dtype=LD)-s)))
        item['max_sigma_absolute_bound']=float(np.max(sigma_bounds[frame]))
        report['fields'][sigma]=item;handled.add(sigma);report['pass'] &= item['pass']
    for name in ('minimum_gamma_lab','maximum_gamma_lab'):
        source=np.maximum(abs(np.asarray(a[name],dtype=LD)),abs(np.asarray(a['mean_gamma_lab'],dtype=LD)))
        gate(name,128*EPS64*source,'per-particle gamma extremum source scale',source)
    raw_names=('particle_kinetic_active_J','particle_kinetic_removed_J','field_energy_interior_J','outward_field_energy_J','prescribed_work_J')
    residual_scale=dynamic+np.sum(abs(faces),axis=1)+abs(a['prescribed_work_J'])
    br=sum(bound[n] for n in raw_names)+bound['particle_kinetic_active_J'][0]+bound['field_energy_interior_J'][0]+128*EPS64*residual_scale
    gate('balance_residual_raw_J',br,'propagated subtractand absolute roundoff',residual_scale)
    corrected_scale=residual_scale+abs(absolute);brc=br+bound['boundary_source_work_J']+128*EPS64*corrected_scale
    gate('balance_residual_J',brc,'raw residual plus source work absolute bound',corrected_scale)
    loss=abs((ka+kr)-ka[0])
    for name in ('balance_residual_raw_J','balance_residual_J'):
        d=abs(np.asarray(b[name],dtype=LD)-np.asarray(a[name],dtype=LD))
        resolved=loss>brc
        report['fields'][name]['particle_energy_change_resolved_samples']=int(np.count_nonzero(resolved))
        report['fields'][name]['max_error_over_resolved_particle_energy_change']=float(np.max(np.divide(d,loss,out=np.zeros_like(d),where=resolved)))
    def ratio_gate(name,num,den,bnum,bden,label):
        q=np.asarray(a[name],dtype=LD);qc=np.asarray(b[name],dtype=LD)
        den=np.broadcast_to(np.asarray(den,dtype=LD),q.shape);bden=np.broadcast_to(np.asarray(bden,dtype=LD),q.shape)
        num=np.asarray(num,dtype=LD);bnum=np.broadcast_to(np.asarray(bnum,dtype=LD),q.shape)
        resolved=abs(den)>bden
        t=np.zeros_like(q)
        t[resolved]=(bnum[resolved]+abs(q[resolved])*bden[resolved])/(abs(den[resolved])-bden[resolved])+128*EPS64*abs(q[resolved])
        diff=abs(qc-q)
        # If denominator cancellation is below its inherited source-energy
        # error, quotient equality is not numerically meaningful. Test its
        # cross-multiplied effect and record this resolution limit explicitly.
        product_t=bnum+(abs(q)+abs(qc))*bden+128*EPS64*(abs(num)+abs(q*den))
        near_ok=abs((qc-q)*den)<=product_t
        passed=np.where(resolved,diff<=t,near_ok)
        known_zero=(den==0)&(bden==0)&(num==0)&(bnum==0)
        passed[known_zero]=(q[known_zero]==0)&(qc[known_zero]==0)
        report['fields'][name]={'pass':bool(np.all(passed)),'gate':label,'max_abs_difference':float(np.max(diff)),
            'resolved_denominator_samples':int(np.count_nonzero(resolved)),
            'ill_conditioned_denominator_samples':int(np.count_nonzero(~resolved)),
            'max_resolved_absolute_bound':float(np.max(t)),'max_unresolved_cross_product_bound':float(np.max(product_t[~resolved])) if np.any(~resolved) else 0.,
            'elements_over_bound':int(np.count_nonzero(~passed)),
            'scope':'Unresolved denominator uses propagated cross-product gate; does not certify precise relative residual.',
            'bitwise_equal':a[name].tobytes()==b[name].tobytes()}
        handled.add(name);report['pass'] &= bool(np.all(passed))
    exchange=loss+abs(fe-fe[0])+abs(a['outward_field_energy_J'])+abs(a['prescribed_work_J'])+abs(a['boundary_source_work_J'])
    ratio_gate('relative_balance_exchange',a['balance_residual_J'],exchange,brc,brc,'residual/exchange propagated absolute bounds')
    initial=abs(ka[0]+fe[0]);bi=bound['particle_kinetic_active_J'][0]+bound['field_energy_interior_J'][0]+128*EPS64*initial
    ratio_gate('relative_balance_initial',a['balance_residual_J'],initial,brc,bi,'residual/initial dynamic energy bounds')
    ratio_gate('field_fraction_kinetic_plus_field',fe,ka+fe,bound['field_energy_interior_J'],bound['particle_kinetic_active_J']+bound['field_energy_interior_J'],'field/(kinetic+field) source bounds')
    total=np.asarray(a['particle_total_active_J'],dtype=LD)
    ratio_gate('field_fraction_including_rest',fe,total+fe,bound['field_energy_interior_J'],bound['particle_total_active_J']+bound['field_energy_interior_J'],'field/(total+field) source bounds')
    ratio_gate('relative_energy_spread_lab',a['sigma_gamma_lab'],a['mean_gamma_lab'],sigma_bounds['lab'],mean_bounds['lab'],'sigma/mean source moment propagation')
    for name in contract['axial_report_only']:
        x=np.asarray(a[name],dtype=LD);d=np.asarray(b[name],dtype=LD)-x
        report['fields'][name]={'reported_only':True,'max_abs_difference':float(np.max(abs(d))),
            'rms_difference':float(np.sqrt(np.mean(d*d))),'bitwise_equal':a[name].tobytes()==b[name].tobytes(),
            'scope':contract['axial_scope']};handled.add(name)
    require(handled==set(a.dtype.names),'Unclassified ledger fields: '+str(set(a.dtype.names)-handled))
    return report
