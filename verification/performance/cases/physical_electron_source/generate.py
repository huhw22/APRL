"""Generate one unit-electron fixed-lab-plane source, without PIC tracking."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import h5py
import numpy as np

from source_math import ALGORITHM,C,DTYPE,E,canonical_json,nominal_electrons,plan_bins,require,sample_records,sha_json,validate_config


def generate(config,output,chunk_records=None):
    """Publish a new complete master; refuse existing final or partial files.

    A closed complete .partial is hard-linked into the final name atomically
    without overwriting a racing writer. Failed generation leaves only the
    unpublished partial with its diagnostic state.
    """
    cfg=validate_config(config)
    chunk=cfg['chunk_records'] if chunk_records is None else chunk_records
    require(isinstance(chunk,int) and not isinstance(chunk,bool) and chunk>0, 'Runtime chunk_records must be positive integer')
    output=Path(output).resolve();partial=output.with_name(output.name+'.partial')
    require(not output.exists() and not partial.exists(), 'Refusing existing final or partial output')
    require(output.parent.is_dir(), 'Output parent must already exist')
    start=time.monotonic();bins=plan_bins(cfg);count=cfg['electrons'];nominal=nominal_electrons(cfg)
    spec={k:v for k,v in cfg.items() if k not in ['chunk_records','compression']}
    spec['algorithm']=ALGORITHM
    payload=hashlib.sha256();cursor=0;max_shell=0.
    kwargs={} if not cfg['compression'] else dict(compression='gzip',compression_opts=cfg['compression'])
    # Bin boundaries need not align with dataset chunks. Keep two ordinary
    # record chunks cached to avoid repeated disk read/modify/write, with an
    # explicit cap for unusually large user chunk requests.
    cache_bytes=max(1024**2,min(256*1024**2,2*min(chunk,count)*DTYPE.itemsize))
    with h5py.File(partial,'x',rdcc_nbytes=cache_bytes) as f:
        group=f.create_group('particles');group.attrs['format_version']=4
        attrs=dict(input_plane_z_m=cfg['phase_space']['output_plane_z_m'],electron_count=np.uint64(count),
            charge_C=count*E,wavelength_m=cfg['wavelength_m'],seed=np.uint64(cfg['seed']),
            config_sha256=sha_json(cfg),source_spec_sha256=sha_json(spec),
            scope='Fixed integer N independent unit-electron6D entrance events; no clones or shot-noise rescaling',
            frame='lab fixed-plane crossing events',count_mode=cfg['count_mode'])
        for k,v in attrs.items():group.attrs[k]=v
        records=group.create_dataset('records',shape=(count,),dtype=DTYPE,chunks=(min(chunk,count),),**kwargs)
        meta=f.create_group('generation');complete=meta.create_dataset('complete',data=np.uint8(0))
        meta.create_dataset('committed_records',data=np.uint64(0))
        meta.create_dataset('config_json',data=canonical_json(cfg),dtype=h5py.string_dtype('utf-8'))
        meta.create_dataset('bin_edges_time_s',data=bins['edges'])
        meta.create_dataset('bin_probabilities',data=bins['probabilities'])
        meta.create_dataset('bin_counts',data=bins['counts'])
        meta.attrs['bin_plane_z_m']=cfg['phase_space']['parameter_plane_z_m']
        meta.attrs['bin_time_scope']='Before exact ballistic map to input plane; actual input arrival bins may differ'
        meta.attrs['algorithm']=ALGORITHM;meta.attrs['numpy_version']=np.__version__;meta.attrs['h5py_version']=h5py.__version__
        meta.attrs['nominal_integrated_electrons']=nominal
        meta.attrs['actual_electrons']=np.uint64(count)
        meta.attrs['actual_charge_over_nominal_charge']=count/nominal
        meta.attrs['bin_ratio_before_snap']=bins['ratio_before_snap']
        meta.attrs['bin_ratio_snapped_to_integer']=bins['ratio_snapped_to_integer']
        meta.attrs['runtime_chunk_records']=chunk
        meta.attrs['hdf5_raw_chunk_cache_bytes']=cache_bytes
        meta.attrs['random_contract']='Separate count key; each zero-based particle ID consumes exactly8 uint64 values from counter2*ID;52bitopen uniform,3 Box-Muller pairs'
        meta.attrs['reproducibility_scope']='Pinned NumPy and platform/libm; record values independent of chunk size or evaluation order. Cross-platform libm bitwise identity is not promised.'
        meta.attrs['time_origin_definition']='Parameter-plane full parabolic profile starts at t=0. Exact ballistic flight time to output plane retained; no common reference flight delay subtracted.'
        for j,n in enumerate(bins['counts']):
            remaining=int(n)
            while remaining:
                size=min(chunk,remaining)
                block,shell=sample_records(cfg,cursor,size,bins['cdf_values'][j],bins['cdf_values'][j+1])
                records[cursor:cursor+size]=block
                payload.update(block.tobytes(order='C'))
                cursor+=size;remaining-=size;max_shell=max(max_shell,shell)
            meta['committed_records'][()]=np.uint64(cursor)
            if j%16==0 or j+1==len(bins['counts']):
                f.flush()
                print(json.dumps(dict(stage='generate',bin=j+1,bins=len(bins['counts']),records=cursor,total=count,elapsed_seconds=time.monotonic()-start)),flush=True)
        require(cursor==count,'Generated total differs from fixedN')
        meta.attrs['record_payload_sha256']=payload.hexdigest()
        meta.attrs['max_abs_reconstructed_gamma_minus_sampled']=max_shell
        meta.attrs['elapsed_seconds_before_publication']=time.monotonic()-start
        complete[()]=np.uint8(1);f.flush()
        handle=f.id.get_vfd_handle()
        if isinstance(handle,int):os.fsync(handle)
    # Atomic no-clobber publication on the same filesystem. Neither the
    # caller nor a concurrent worker can overwrite an existing final source.
    os.link(partial,output)
    partial.unlink()
    result=dict(complete=True,output=str(output),electron_count=count,charge_C=count*E,
        nominal_integrated_electrons=nominal,actual_charge_over_nominal_charge=count/nominal,
        records_bytes=count*DTYPE.itemsize,file_bytes=output.stat().st_size,
        record_payload_sha256=payload.hexdigest(),config_sha256=sha_json(cfg),source_spec_sha256=sha_json(spec),
        algorithm=ALGORITHM,count_mode=cfg['count_mode'],bins=len(bins['counts']),
        max_abs_reconstructed_gamma_minus_sampled=max_shell,elapsed_seconds=time.monotonic()-start)
    print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--chunk-records',type=int)
    args=parser.parse_args()
    generate(json.loads(args.config.read_text()),args.output,args.chunk_records)


if __name__=='__main__':main()
