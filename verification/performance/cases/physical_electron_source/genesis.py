"""Deterministic, bounded-memory unit-electron source to Genesis native import.

The input is a closed APRL v4 fixed-lab-plane event master, not a distribution
to resample. Only native ``&importbeam`` is suitable: no random number generator,
beamlet replication or extra shot noise appears in this conversion.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

import h5py
import numpy as np

C = 299792458.0
E = 1.602176634e-19
GENESIS_CE = 4.8032045e-11
SCHEMA = 'APRL_genesis_unit_electron_v1'
LD = np.longdouble
TWO_PI_LD = 2*np.arccos(LD(-1))
TWO_PI = float(TWO_PI_LD)
FIELDS = ('x', 'y', 'px', 'py', 'gamma', 'theta')


def require(value, message):
    if not value:
        raise ValueError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x))+'\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def stat(path):
    s = Path(path).stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns,
                inode=s.st_ino, device=s.st_dev)


def scalar_int(value, name, minimum=0):
    require(not isinstance(value, (bool, np.bool_)), name+' must be an integer')
    converted = int(value)
    require(converted == value and converted >= minimum,
            name+' must be an integer >= '+str(minimum))
    return converted


def gamma_from_u(proper):
    """The stable scaled mass shell used by the PIC frame implementation."""
    u = np.asarray(proper, dtype=np.float64)
    scale = np.maximum(1., np.max(np.abs(u), axis=1))
    inverse = 1./scale
    a, b, d = (u[:, j]*inverse for j in range(3))
    return scale*np.sqrt(inverse*inverse+a*a+b*b+d*d)


def phase_coordinates(arrival_s, wavelength_m):
    """Absolute head-directed slice integer and optical fractional phase.

    The common arrival-time offset zero is the anchor. This permits counting
    and writing in two event scans, without a full-N sort or data-dependent
    optical origin. Native metadata additionally records the integer origin.
    """
    coordinate = -LD(C)*np.asarray(arrival_s, dtype=LD)/LD(wavelength_m)
    require(np.isfinite(coordinate).all(), 'Nonfinite optical coordinate')
    require(np.all(abs(coordinate) < LD(2)**62),
            'Arrival offsets exceed representable native slice range')
    whole = np.floor(coordinate)
    integer = np.asarray(whole, dtype=np.int64)
    phase = np.asarray(TWO_PI_LD*(coordinate-whole), dtype=np.float64)
    # A fraction strictly below one can round to float64 2pi. Preserve its
    # half-open slice with a one-ULP representational correction, reported below.
    rounded = phase >= TWO_PI
    phase[rounded] = np.nextafter(TWO_PI, 0.)
    require(np.all(phase >= 0.) and np.all(phase < TWO_PI),
            'Optical phase outside its half-open native slice')
    return integer, phase, int(np.count_nonzero(rounded))


def master_group(source):
    require('generation/complete' in source and
            int(source['generation/complete'][()]) == 1,
            'Master is not normally completed; no event records were read')
    require('particles/records' in source, 'Master has no particle records')
    group = source['particles']
    require(int(group.attrs.get('format_version', -1)) == 4,
            'Expected APRL fixed-lab-plane event format_version=4')
    require('input_plane_z_m' in group.attrs and
            np.isfinite(group.attrs['input_plane_z_m']),
            'Master must name one finite common lab plane')
    data = group['records']
    require(data.ndim == 1 and len(data) > 0, 'Master contains no event records')
    if 'generation/committed_records' in source:
        require(int(source['generation/committed_records'][()]) == len(data),
                'Master committed record count differs')
    expected = {'plane_position_m', 'arrival_time_offset_s', 'proper_velocity',
                'source_id', 'macro_weight'}
    require(data.dtype.names is not None and expected <= set(data.dtype.names),
            'Master does not contain the v4 crossing-event fields')
    require(data.dtype.itemsize == 64 and data.dtype['source_id'].kind == 'u' and
            data.dtype['source_id'].itemsize == 8,
            'Master must use the native 64-byte event/uint64 source-ID contract')
    require(data.dtype['plane_position_m'].shape == (2,) and
            data.dtype['proper_velocity'].shape == (3,), 'Master vector layout differs')
    require(scalar_int(group.attrs.get('electron_count', -1), 'electron_count', 1)
            == len(data), 'Unit-electron count differs from master record count')
    return group, data


def checked_chunk(data, first, last):
    a = data[first:last]
    require(np.array_equal(a['source_id'],
            np.arange(first+1, last+1, dtype=np.uint64)),
            'Master source IDs must be contiguous 1..N in record order')
    require(np.all(a['macro_weight'] == 1.),
            'Only one physical electron per record is accepted')
    for name in ['plane_position_m', 'arrival_time_offset_s', 'proper_velocity']:
        require(np.isfinite(a[name]).all(), 'Nonfinite master field: '+name)
    require(np.all(a['proper_velocity'][:, 2] > 0.),
            'Master contains a nonforward crossing event')
    gamma = gamma_from_u(a['proper_velocity'])
    require(np.isfinite(gamma).all() and np.all(gamma >= 1.),
            'Master contains an invalid mass shell')
    return a, gamma


def convert(master, output, wavelength_m, pad_slices=16, ranks_multiple=16,
            chunk_records=65536):
    """Convert one closed physical-electron master without changing its sample.

    Memory is O(chunk_records + represented native slices). The first pass
    validates the master and counts true arrival-time bins; the second pass
    performs stable within-chunk grouping and contiguous HDF5 writes. Every
    native slice also stores source_id for input-only independent auditing.
    """
    started = time.monotonic()
    master, output = Path(master).resolve(), Path(output).resolve()
    partial = output.with_name(output.name+'.partial')
    manifest = output.with_name(output.name+'.conversion.json')
    require(master.is_file(), 'Master does not exist')
    require(master != output and master != partial, 'Output would replace master')
    require(not output.exists() and not partial.exists() and not manifest.exists(),
            'Refusing to overwrite conversion output or partial output')
    require(np.isfinite(wavelength_m) and wavelength_m > 0.,
            'wavelength_m must be finite and positive')
    wavelength_m = float(wavelength_m)
    pad_slices = scalar_int(pad_slices, 'pad_slices')
    ranks_multiple = scalar_int(ranks_multiple, 'ranks_multiple', 1)
    chunk_records = scalar_int(chunk_records, 'chunk_records', 1)
    before = stat(master)
    counts_map = {}
    record_hash = hashlib.sha256()
    gamma_sum = LD(0)
    minimum_time, maximum_time = math.inf, -math.inf
    with h5py.File(master, 'r') as source:
        group, data = master_group(source)
        count = len(data)
        input_plane_z_m = float(group.attrs['input_plane_z_m'])
        for first in range(0, count, chunk_records):
            last = min(count, first+chunk_records)
            a, gamma = checked_chunk(data, first, last)
            record_hash.update(a.tobytes(order='C'))
            gamma_sum += np.sum(gamma.astype(LD), dtype=LD)
            minimum_time = min(minimum_time, float(a['arrival_time_offset_s'].min()))
            maximum_time = max(maximum_time, float(a['arrival_time_offset_s'].max()))
            indices, _, _ = phase_coordinates(a['arrival_time_offset_s'], wavelength_m)
            low, high = int(indices.min()), int(indices.max())
            require(high-low < 1_000_000, 'Chunk spans too many native slices')
            block_counts = np.bincount(indices-low, minlength=high-low+1)
            for j in np.flatnonzero(block_counts):
                index = low+int(j)
                counts_map[index] = counts_map.get(index, 0)+int(block_counts[j])
        first_bin, last_bin = min(counts_map), max(counts_map)
        logical_origin = first_bin-pad_slices
        physical_span = last_bin-first_bin+1
        require(physical_span+2*pad_slices <= 1_000_000,
                'More than one million native slices are not supported')
        slices = ((physical_span+2*pad_slices+ranks_multiple-1)//ranks_multiple)*ranks_multiple
        require(slices < np.iinfo(np.int32).max, 'Native slice count exceeds int32')
        counts = np.zeros(slices, dtype=np.int64)
        for index, number in counts_map.items():counts[index-logical_origin] = number
        require(int(counts.sum()) == count, 'Slice counts lost source records')
        require(np.all(counts <= np.iinfo(np.int32).max),
                'A native slice exceeds Genesis int32 particle indexing')
        refposition = float(LD(logical_origin)*LD(wavelength_m))
        source_hash = sha(master)
        require(stat(master) == before, 'Master changed during counting')
        output.parent.mkdir(parents=True, exist_ok=True)
        upper_rounding = 0
        with h5py.File(partial, 'x') as native:
            native.attrs['converter_schema'] = SCHEMA
            native.attrs['frame'] = 'lab-fixed-plane'
            native.attrs['input_plane_z_m'] = input_plane_z_m
            native.attrs['time_anchor_s'] = 0.
            native.attrs['logical_slice_origin_integer'] = logical_origin
            native.attrs['electron_count'] = count
            native.attrs['source_master_sha256'] = source_hash
            native.attrs['source_record_payload_sha256'] = record_hash.hexdigest()
            native.attrs['source_id_contract'] = 'Input audit only: source_id equals master row index plus one; Genesis does not track this extra field'
            native.attrs['phase_definition'] = 's=-c*arrival_time_offset_s; absolute bin=floor(s/lambda); theta=2pi*fraction(s/lambda)'
            native.attrs['momentum_units'] = 'px and py are p/(m_e*c), not angles'
            native.attrs['gamma_definition'] = 'Scaled mass shell from master proper_velocity; no resampling of sampled gamma'
            for name, value in [('refposition', refposition), ('slicelength', wavelength_m),
                                ('slicespacing', wavelength_m)]:native[name] = [value]
            for name, value in [('slicecount', slices), ('beamletsize', 1), ('one4one', 1)]:
                native.create_dataset(name, data=np.array([value], dtype='<i4'))
            audit = native.create_group('conversion')
            audit['complete'] = np.uint8(0)
            audit['slice_counts'] = counts
            audit.attrs['time_anchor_s'] = 0.
            audit.attrs['logical_slice_origin_integer'] = logical_origin
            audit.attrs['padding_requested_each_side'] = pad_slices
            audit.attrs['padding_added_to_rank_multiple'] = slices-(physical_span+2*pad_slices)
            audit.attrs['rank_multiple'] = ranks_multiple
            audit.attrs['sample'] = 1
            audit.attrs['source_master_path'] = str(master)
            audit.attrs['source_id_contiguous_first'] = 1
            audit.attrs['source_id_contiguous_last'] = count
            audit.attrs['charge_C_exact_SI'] = count*E
            audit.attrs['genesis_ce_relative_difference'] = GENESIS_CE/(E*C)-1.
            audit.attrs['random_draws_or_particle_replicas_added'] = 0
            groups = []
            for j, number in enumerate(counts):
                g = native.create_group(f'slice{j+1:06d}')
                g['current'] = [int(number)*E*C/wavelength_m]
                g.attrs['absolute_slice_index'] = logical_origin+j
                g.attrs['input_electrons'] = int(number)
                for name in FIELDS:
                    g.create_dataset(name, shape=(int(number),), dtype='<f8',
                        chunks=(min(chunk_records, int(number)),) if number else None)
                g.create_dataset('source_id', shape=(int(number),), dtype='<u8',
                    chunks=(min(chunk_records, int(number)),) if number else None)
                groups.append(g)
            cursor = np.zeros(slices, dtype=np.int64)
            for first in range(0, count, chunk_records):
                last = min(count, first+chunk_records)
                a, gamma = checked_chunk(data, first, last)
                integer, phase, rounded = phase_coordinates(a['arrival_time_offset_s'], wavelength_m)
                index = integer-logical_origin
                require(np.all(index >= 0) and np.all(index < slices),
                        'Source moved outside counted native window')
                upper_rounding += rounded
                order = np.argsort(index, kind='stable')
                ordered_index = index[order]
                cuts = np.r_[0, np.flatnonzero(np.diff(ordered_index))+1, len(order)]
                for left, right in zip(cuts[:-1], cuts[1:]):
                    j = int(ordered_index[left]);take = order[left:right]
                    start, stop = int(cursor[j]), int(cursor[j])+len(take)
                    require(stop <= counts[j], 'Source slice changed between scans')
                    g = groups[j]
                    for name, values in [('x', a['plane_position_m'][:, 0]),
                                         ('y', a['plane_position_m'][:, 1]),
                                         ('px', a['proper_velocity'][:, 0]),
                                         ('py', a['proper_velocity'][:, 1]),
                                         ('gamma', gamma), ('theta', phase),
                                         ('source_id', a['source_id'])]:
                        g[name][start:stop] = values[take]
                    cursor[j] = stop
            require(np.array_equal(cursor, counts), 'Incomplete native slice writes')
            audit.attrs['theta_upper_boundary_roundoff_count'] = upper_rounding
            require(stat(master) == before, 'Master changed during conversion')
            audit['complete'][()] = np.uint8(1)
            native.flush()
            handle = native.id.get_vfd_handle()
            if isinstance(handle, int):os.fsync(handle)
    require(stat(master) == before, 'Master changed before conversion publication')
    # Same-filesystem no-clobber publication, including against a racing writer.
    os.link(partial, output)
    partial.unlink()
    result = dict(schema=SCHEMA, complete=True, source_master=str(master),
        source_master_sha256=source_hash, source_master_identity=before,
        source_record_payload_sha256=record_hash.hexdigest(), output=str(output),
        output_sha256=sha(output), output_bytes=output.stat().st_size,
        electron_count=count, charge_C=count*E, mean_gamma=float(gamma_sum/LD(count)),
        input_plane_z_m=input_plane_z_m,
        wavelength_m=wavelength_m, time_anchor_s=0., refposition_m=refposition,
        logical_slice_origin_integer=logical_origin, native_slices=slices,
        nonempty_slices=int(np.count_nonzero(counts)), empty_slices=int(np.count_nonzero(counts == 0)),
        slice_counts=counts.tolist(), pad_slices=pad_slices, ranks_multiple=ranks_multiple,
        chunk_records=chunk_records, arrival_time_range_s=[minimum_time, maximum_time],
        particle_payload_bytes=56*count, source_id_audit_bytes=8*count,
        native_six_coordinate_bytes=48*count, random_draws_added=0,
        theta_upper_boundary_roundoff_count=upper_rounding,
        seconds=time.monotonic()-started,
        limitations=['Native input stores source_id only for initial auditing; Genesis does not propagate that ID.',
            'Gamma is reconstructed from the mass shell of master proper velocity; equality to the pre-mass-shell random draw is not asserted.',
            'Native slice refposition and phase are floating-point; use the stored integer slice origin and common time anchor for independent event reconstruction.',
            'Genesis one4one sorting uses its historical e*c constant; the relative difference is recorded, not hidden by changing the sample.'])
    with manifest.open('x') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def make_import_deck(nativefile, run_dir):
    """Prepare a real Genesis import/write-only check, with no tracking."""
    nativefile, run_dir = Path(nativefile).resolve(), Path(run_dir).resolve()
    require(not run_dir.exists(), 'Refusing to overwrite a native import test')
    with h5py.File(nativefile, 'r') as f:
        require('conversion/complete' in f and int(f['conversion/complete'][()]) == 1,
                'Native converter output is not complete')
        wavelength = float(f['slicelength'][0])
        require(int(f['one4one'][0]) == 1 and int(f['beamletsize'][0]) == 1,
                'Import test expects unit-electron native input')
    run_dir.mkdir(parents=True)
    (run_dir/'roundtrip.lat').write_text('D: DRIFT={l=0.001};\nTEST: LINE={D};\n')
    deck = f'''&setup
rootname=import-check
lattice=roundtrip.lat
beamline=TEST
lambda0={wavelength:.17g}
gamma0=298.2608843802298
delz=0.001
one4one=true
shotnoise=false
npart=1
nbins=1
&end
&importbeam
file={nativefile}
time=true
&end
&write
beam=roundtrip
&end
'''
    (run_dir/'roundtrip.in').write_text(deck)
    result = dict(native_input=str(nativefile), native_input_sha256=sha(nativefile),
        deck=str(run_dir/'roundtrip.in'), expected_dump=str(run_dir/'roundtrip.par.h5'),
        method='Setup, native importbeam and write only; no field generation, noise loading, sorting or track block')
    save(run_dir/'PLAN.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('master', nargs='?');p.add_argument('output', nargs='?')
    p.add_argument('--wavelength-m', type=float)
    p.add_argument('--pad-slices', type=int, default=16)
    p.add_argument('--ranks-multiple', type=int, default=16)
    p.add_argument('--chunk-records', type=int, default=65536)
    p.add_argument('--write-test-deck', nargs=2, metavar=('NATIVE_FILE','NEW_DIRECTORY'))
    a = p.parse_args()
    if a.write_test_deck:
        require(not a.master and not a.output, 'Deck preparation takes no master/output arguments')
        result = make_import_deck(*a.write_test_deck)
    else:
        require(a.master and a.output and a.wavelength_m is not None,
                'MASTER OUTPUT and --wavelength-m are required')
        result = convert(a.master,a.output,a.wavelength_m,a.pad_slices,
                         a.ranks_multiple,a.chunk_records)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':main()
