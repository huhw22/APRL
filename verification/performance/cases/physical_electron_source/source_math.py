"""Physical fixed-N independent-electron source: pure sampling operations.

An ID consumes eight uint64 Philox values, independent of chunking. Six values
form three Box--Muller pairs, one draws arrival time, and one is reserved. Count
draws use a separate key. No beamlets, empirical moment matching or BF fitting.
"""
import copy
import hashlib
import json
import math

import numpy as np

C = 299792458.0
E = 1.602176634e-19
SCHEMA = 'APRL_physical_electron_source_v1'
ALGORITHM = 'fixed-n-parabolic-conditional-philox4x64-boxmuller-v1'
DTYPE = np.dtype([('plane_position_m','<f8',(2,)),
    ('arrival_time_offset_s','<f8'),('proper_velocity','<f8',(3,)),
    ('source_id','<u8'),('macro_weight','<f8')])
assert DTYPE.itemsize == 64


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha_json(value):
    return hashlib.sha256(canonical_json(value).encode('utf-8')).hexdigest()


def finite_float(value, name):
    require(not isinstance(value, bool), 'Boolean is not a number: '+name)
    result = float(value)
    require(math.isfinite(result), 'Nonfinite '+name)
    return result


def cdf_fraction(t, duration):
    """CDF of I(t)=4 I_peak (t/T)(1-t/T), 0<=t<=T."""
    y = np.asarray(t, dtype=np.longdouble)/np.longdouble(duration)
    return y*y*(3-2*y)


def inverse_cdf_fraction(probability):
    p = np.asarray(probability, dtype=np.float64)
    require(np.all((p >= 0) & (p <= 1)), 'CDF argument outside [0,1]')
    return .5-np.sin(np.arcsin(1-2*p)/3)


def nominal_electrons(config):
    p = config['profile']
    fraction = cdf_fraction(p['window_stop_s'], p['full_duration_s'])-cdf_fraction(p['window_start_s'], p['full_duration_s'])
    charge = np.longdouble(2)/3*np.longdouble(p['peak_current_A'])*np.longdouble(p['full_duration_s'])*fraction
    return float(charge/np.longdouble(str(E)))


def validate_config(config):
    """Return a normalized copy; never change the caller's configuration."""
    require(isinstance(config, dict), 'Configuration must be a JSON object')
    cfg = copy.deepcopy(config)
    allowed = {'schema','electrons','count_mode','seed','profile','phase_space',
               'wavelength_m','chunk_records','compression'}
    require(not set(cfg)-allowed, 'Unknown source option(s): '+', '.join(sorted(set(cfg)-allowed)))
    require(cfg.get('schema') == SCHEMA, 'Wrong source schema')
    for key in ['electrons','seed']:
        require(isinstance(cfg.get(key), int) and not isinstance(cfg[key], bool), key+' must be an integer')
    require(0 < cfg['electrons'] <= 2**53, 'electrons must be in [1,2^53] for exact PIC charge normalization')
    require(0 <= cfg['seed'] < 2**64, 'seed must be uint64')
    cfg.setdefault('count_mode','physical_fixed_n')
    require(cfg['count_mode'] in ['physical_fixed_n','diagnostic_fixed_n'], 'Unknown count_mode')
    p = cfg['profile']
    required_profile = {'type','full_duration_s','peak_current_A','window_start_s','window_stop_s','bin_width_s'}
    require(set(p) == required_profile, 'Profile fields differ from parabolic-time contract')
    require(p['type'] == 'parabolic_time', 'Only declared parabolic_time is supported')
    for key in required_profile-{'type'}:p[key] = finite_float(p[key], 'profile.'+key)
    require(p['full_duration_s'] > 0 and p['peak_current_A'] > 0 and p['bin_width_s'] > 0, 'Profile scales must be positive')
    require(0 <= p['window_start_s'] < p['window_stop_s'] <= p['full_duration_s'], 'Profile crop lies outside full support')
    required_phase = {'gamma_mean','gamma_sigma','sigma_x_m','sigma_y_m','emittance_n_x_m',
        'emittance_n_y_m','alpha_x','alpha_y','parameter_plane_z_m','output_plane_z_m'}
    s = cfg['phase_space']
    require(set(s) == required_phase, 'Phase-space fields differ from analytic Gaussian contract')
    for key in required_phase:s[key] = finite_float(s[key], 'phase_space.'+key)
    require(s['gamma_mean'] > 1 and s['gamma_sigma'] >= 0, 'Invalid Lorentz-factor distribution')
    require(s['sigma_x_m'] > 0 and s['sigma_y_m'] > 0, 'Transverse widths must be positive')
    require(s['emittance_n_x_m'] >= 0 and s['emittance_n_y_m'] >= 0, 'Normalized emittances must be nonnegative')
    cfg['wavelength_m'] = finite_float(cfg['wavelength_m'], 'wavelength_m')
    require(cfg['wavelength_m'] > 0, 'wavelength_m must be positive')
    cfg.setdefault('chunk_records', 262144)
    cfg.setdefault('compression', 0)
    require(isinstance(cfg['chunk_records'], int) and not isinstance(cfg['chunk_records'], bool) and cfg['chunk_records'] > 0, 'chunk_records must be a positive integer')
    require(isinstance(cfg['compression'], int) and not isinstance(cfg['compression'], bool) and 0 <= cfg['compression'] <= 9, 'compression must be integer0..9')
    expected = nominal_electrons(cfg)
    require(expected > 0, 'Nominal crop charge underflow')
    if cfg['count_mode'] == 'physical_fixed_n':
        require(cfg['electrons'] == int(np.rint(expected)), 'Physical fixed-N must equal round(integrated profile charge/e); use explicit diagnostic_fixed_n for low-charge tests')
    return cfg


def philox_key(seed, namespace):
    # Separate namespaces are algorithm constants, not Python hash values.
    sequence = np.random.SeedSequence([seed & 0xffffffff, seed >> 32, namespace])
    return sequence.generate_state(2, dtype=np.uint64)


def plan_bins(config):
    """Multinomial counts at the parameter plane, not after the drift map."""
    cfg = validate_config(config);p = cfg['profile']
    ratio = (p['window_stop_s']-p['window_start_s'])/p['bin_width_s']
    nearest = round(ratio)
    snapped = abs(ratio-nearest) <= 64*np.finfo(float).eps*max(1.,abs(ratio))
    bins = int(nearest if snapped else math.ceil(ratio))
    require(1 <= bins <= 1000000, 'Generation bin count outside [1,1000000]')
    edges = p['window_start_s']+np.arange(bins+1,dtype=float)*p['bin_width_s']
    edges[-1] = p['window_stop_s']
    require(np.all(np.diff(edges) > 0), 'Nonpositive-width generation bin')
    cumulative = cdf_fraction(edges, p['full_duration_s'])
    probabilities = np.asarray(np.diff(cumulative)/(cumulative[-1]-cumulative[0]),float)
    probabilities /= probabilities.sum()
    require(np.all(probabilities > 0), 'Nonpositive generation-bin probability')
    rng = np.random.Generator(np.random.Philox(key=philox_key(cfg['seed'],0x434f554e)))
    counts = rng.multinomial(cfg['electrons'], probabilities)
    require(int(counts.sum()) == cfg['electrons'], 'Multinomial total differs')
    return dict(edges=edges, cdf_values=np.asarray(cumulative,float), probabilities=probabilities,
        counts=counts.astype(np.uint64), ratio_before_snap=ratio, ratio_snapped_to_integer=snapped)


def random_by_id(seed, start_index, count):
    """Zero-based IDs -> one open uniform plus six independent N(0,1)."""
    require(isinstance(start_index,(int,np.integer)) and start_index >= 0, 'Invalid first random index')
    require(isinstance(count,(int,np.integer)) and count >= 0, 'Invalid random count')
    raw = np.random.Philox(counter=2*int(start_index), key=philox_key(seed,0x50415254)).random_raw(8*int(count)).reshape((-1,8))
    # 52 bits makes both endpoints strictly excluded after float rounding.
    values = ((raw >> np.uint64(12)).astype(float)+.5)*2.**-52
    normal = np.empty((count,6),dtype=float)
    for j in range(3):
        radius = np.sqrt(-2*np.log(values[:,2*j]));angle = 2*np.pi*values[:,2*j+1]
        normal[:,2*j] = radius*np.cos(angle);normal[:,2*j+1] = radius*np.sin(angle)
    return values[:,6],normal


def sample_records(config, start_index, count, cdf_low, cdf_high):
    """Independent particle records conditional on one generation-time bin.

    config must already have been validated. No statistic is forced to its
    sample target, and impossible on-shell realizations cause a hard failure.
    """
    u,z = random_by_id(config['seed'],start_index,count)
    tmag = config['profile']['full_duration_s']*inverse_cdf_fraction(cdf_low+(cdf_high-cdf_low)*u)
    s = config['phase_space'];x=s['sigma_x_m']*z[:,0];y=s['sigma_y_m']*z[:,2]
    ux=s['emittance_n_x_m']/s['sigma_x_m']*(-s['alpha_x']*z[:,0]+z[:,1])
    uy=s['emittance_n_y_m']/s['sigma_y_m']*(-s['alpha_y']*z[:,2]+z[:,3])
    gamma=s['gamma_mean']+s['gamma_sigma']*z[:,4]
    uz2=(gamma-1)*(gamma+1)-ux*ux-uy*uy
    require(np.all(np.isfinite(gamma)) and np.all(gamma > 1) and np.all(uz2 > 0), 'Unphysical Gaussian draw: gamma/forward longitudinal momentum; no clipping or hidden resampling')
    uz=np.sqrt(uz2);distance=s['parameter_plane_z_m']-s['output_plane_z_m']
    result=np.empty(count,dtype=DTYPE)
    result['plane_position_m'][:,0]=x-distance*ux/uz
    result['plane_position_m'][:,1]=y-distance*uy/uz
    result['proper_velocity'][:,0]=ux;result['proper_velocity'][:,1]=uy;result['proper_velocity'][:,2]=uz
    # The v4 master stores u, not the sampled gamma scalar. Use the on-shell
    # gamma of those stored float64 momenta for the event worldline itself;
    # its <=few-ULP difference from the sampled gamma is reported separately.
    reconstructed=np.sqrt(1+np.sum(result['proper_velocity']**2,axis=1))
    result['arrival_time_offset_s']=tmag-distance*reconstructed/(C*uz)
    result['source_id']=np.arange(start_index+1,start_index+count+1,dtype=np.uint64)
    result['macro_weight']=1.
    require(all(np.isfinite(result[n]).all() for n in ['plane_position_m','arrival_time_offset_s','proper_velocity']), 'Nonfinite generated event')
    return result,float(np.max(abs(reconstructed-gamma))) if count else 0.
