#!/usr/bin/env python3
"""List immutable diagnostic prefixes in a checkpoint ancestry (read-only).

Only checkpoint metadata, dataset shapes, scalar commit markers and identity
attributes are read. Scientific arrays are never copied or concatenated.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys
from unittest import mock

import h5py
import yaml

MAGIC = 0x4150524C434B5031
FORMAT_VERSION = 1
MAX_METADATA_BYTES = 64 * 1024 * 1024
GROUPS = {'/field_plane', '/particle_plane', '/ballistic_reference',
          '/trajectory', '/energy_ledger'}


class PrefixError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise PrefixError(message)


def nonnegative_integer(value, label):
    require(isinstance(value, int) and not isinstance(value, bool) and value >= 0,
            label + ' must be a nonnegative integer')
    return value


def absolute_file(value, label):
    require(isinstance(value, str) and value and '\0' not in value,
            label + ' must name a file')
    path = Path(value)
    require(path.is_absolute(), label + ' must be absolute; relative ancestry is ambiguous: ' + value)
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise PrefixError(label + ' is missing or inaccessible: ' + value) from error
    require(resolved.is_file(), label + ' is not a regular file: ' + value)
    return resolved


def identity(path):
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def text_attribute(value):
    return value.decode('utf-8') if isinstance(value, bytes) else str(value)


def rank_texts(file, name, ranks):
    group_path = '/blocks/metadata/' + name
    require(group_path in file, 'Missing checkpoint metadata: ' + name)
    group = file[group_path]
    require(group['kind'].shape == (1,) and int(group['kind'][0]) == 1,
            'Checkpoint metadata is not a byte block: ' + name)
    counts_data, data = group['rank_counts'], group['data']
    require(counts_data.shape == (ranks,) and counts_data.dtype.kind == 'u',
            'Invalid checkpoint rank counts: ' + name)
    require(data.ndim == 1 and data.dtype.kind == 'u' and data.dtype.itemsize == 1,
            'Invalid checkpoint byte dataset: ' + name)
    counts = [int(count) for count in counts_data[:]]
    require(sum(counts) == data.shape[0] and sum(counts) <= MAX_METADATA_BYTES,
            'Invalid or oversized checkpoint metadata length: ' + name)
    offset, texts = 0, []
    for count in counts:
        texts.append(data[offset:offset + count].tobytes().decode('utf-8'))
        offset += count
    return texts


def common_text(file, name, ranks):
    texts = rank_texts(file, name, ranks)
    require(texts and all(value == texts[0] for value in texts),
            'Rank metadata differs: ' + name)
    return texts[0]


def read_checkpoint(path):
    require(not path.name.endswith('.partial'), 'Unpublished .partial checkpoint is not readable history')
    with h5py.File(path, 'r') as file:
        require('/header' in file and file['/header'].shape == (12,), 'Invalid checkpoint header')
        header = [int(value) for value in file['/header'][:]]
        require(header[0] == MAGIC and header[1] == FORMAT_VERSION and header[2] == 1,
                'Checkpoint has an unsupported schema or is not committed: ' + str(path))
        ranks = header[3]
        require(0 < ranks <= MAX_METADATA_BYTES // 8, 'Invalid checkpoint MPI size')
        parent = common_text(file, 'parent_checkpoint', ranks)
        source_run_id = common_text(file, 'source_run_id', ranks)
        config = yaml.safe_load(common_text(file, 'configuration', ranks))
        require(isinstance(config, dict), 'Checkpoint configuration is not a mapping')
        checkpoint_config = config.get('checkpoint') or {}
        require(isinstance(checkpoint_config, dict), 'Checkpoint configuration block is not a mapping')
        reset = checkpoint_config.get('reset_diagnostics', False)
        require(isinstance(reset, bool), 'Invalid reset_diagnostics in checkpoint configuration')
        prefixes, step, time_box = [], None, None
        for rank, text in enumerate(rank_texts(file, 'diagnostic_prefixes', ranks)):
            record = yaml.safe_load(text)
            require(isinstance(record, dict) and isinstance(record.get('outputs'), list),
                    'Invalid diagnostic prefix metadata on rank ' + str(rank))
            this_step = nonnegative_integer(record.get('step'), 'Checkpoint step')
            this_time = record.get('time_box_s')
            require(isinstance(this_time, (int, float)) and not isinstance(this_time, bool) and
                    math.isfinite(this_time) and this_time >= 0, 'Invalid checkpoint boost time')
            require(record.get('parent_checkpoint') == parent, 'Prefix parent differs from checkpoint parent')
            if step is None:
                step, time_box = this_step, this_time
            require(this_step == step and this_time == time_box, 'Checkpoint rank step/time mismatch')
            for output in record['outputs']:
                require(isinstance(output, dict), 'Invalid diagnostic output prefix')
                group = output.get('group')
                require(isinstance(group, str), 'Missing diagnostic group')
                group = '/' + group.lstrip('/')
                require(group in GROUPS, 'Unknown diagnostic group: ' + group)
                records = nonnegative_integer(output.get('records'), 'Prefix records')
                validation = nonnegative_integer(output.get('validation_records', 0), 'Validation prefix records')
                require(group == '/ballistic_reference' or validation == 0,
                        'Validation prefix on a non-reference output')
                prefixes.append(dict(rank=rank, file=output.get('file'), group=group,
                                     records=records, validation_records=validation))
    return dict(checkpoint=str(path), parent_checkpoint=parent, step=step,
                time_box_s=time_box, ranks=ranks, source_run_id=source_run_id,
                reset_diagnostics=reset, outputs=prefixes)


def scalar_count(group, name):
    require(name in group and group[name].shape == () and group[name].dtype.kind == 'u',
            'Missing or invalid scalar commit marker: ' + group.name + '/' + name)
    return int(group[name][()])


def verify_prefix(prefix, generation):
    """Validate sufficient committed source length, without reading array data."""
    with h5py.File(prefix['file'], 'r') as file:
        require(prefix['group'] in file, 'Missing source group: ' + prefix['file'] + prefix['group'])
        group = file[prefix['group']]
        require(text_attribute(group.attrs.get('run_id', '')) == generation['source_run_id'],
                'Source output run_id differs from its checkpoint: ' + prefix['file'])
        parent = generation['parent_checkpoint']
        output_parent = text_attribute(group.attrs.get('restart_parent', ''))
        if parent:
            require(output_parent and absolute_file(output_parent, 'Source restart_parent') == Path(parent),
                    'Source output points to a different restart parent: ' + prefix['file'])
            require(int(group.attrs.get('restart_step', -1)) == generation['start_step'] and
                    float(group.attrs.get('restart_time_box_s', -1)) == generation['start_time_box_s'],
                    'Source output restart step/time differs: ' + prefix['file'])
        else:
            require(not output_parent, 'Original output unexpectedly has a restart parent')
        expected = prefix['records']
        if prefix['group'] == '/field_plane':
            committed = scalar_count(group, 'committed_samples')
            dimensions = {'time_s': 1, 'fields': 3}
        else:
            committed = scalar_count(group, 'committed_records')
            dimensions = {'records': 1}
        require(committed >= expected, 'Source committed length is below the frozen prefix: ' + prefix['file'])
        datasets = {}
        for name, ndim in dimensions.items():
            require(name in group and group[name].ndim == ndim and group[name].shape[0] >= committed,
                    'Source dataset is shorter than its commit marker: ' + prefix['file'] + '/' + name)
            datasets[name] = dict(start=0, stop=expected, source_shape=list(group[name].shape))
        if prefix['group'] == '/ballistic_reference':
            count = prefix['validation_records']
            if 'committed_validation_records' in group or 'validation_records' in group or count:
                validation_committed = scalar_count(group, 'committed_validation_records')
                require(validation_committed >= count and 'validation_records' in group and
                        group['validation_records'].ndim == 1 and
                        group['validation_records'].shape[0] >= validation_committed,
                        'Invalid source validation prefix: ' + prefix['file'])
                datasets['validation_records'] = dict(start=0, stop=count,
                    source_shape=list(group['validation_records'].shape))
        return dict(status='verified', committed_records=committed, datasets=datasets,
                    plane_name=text_attribute(group.attrs.get('name', '')))


def build_manifest(checkpoint, verify_sources=False, max_depth=256):
    """Return generations oldest-first and distinct slices; never merge them."""
    require(isinstance(max_depth, int) and max_depth > 0, 'max_depth must be positive')
    path = Path(checkpoint).expanduser().resolve(strict=True)
    require(path.is_file(), 'Checkpoint is not a regular file')
    requested_checkpoint = str(path)
    reverse, seen_checkpoints = [], set()
    while True:
        key = identity(path)
        require(key not in seen_checkpoints, 'Checkpoint ancestry cycle or repeated file: ' + str(path))
        require(len(reverse) < max_depth, 'Checkpoint ancestry exceeds max_depth')
        seen_checkpoints.add(key)
        generation = read_checkpoint(path)
        reverse.append(generation)
        if not generation['parent_checkpoint']:
            break
        path = absolute_file(generation['parent_checkpoint'], 'Parent checkpoint')
        generation['parent_checkpoint'] = str(path)
    generations, slices, seen_outputs = list(reversed(reverse)), [], set()
    previous = None
    for index, generation in enumerate(generations):
        generation['index'] = index
        generation['start_step'] = previous['step'] if previous else 0
        generation['start_time_box_s'] = previous['time_box_s'] if previous else 0.0
        if previous:
            require(generation['step'] >= previous['step'] and
                    generation['time_box_s'] >= previous['time_box_s'], 'Ancestry runs backwards in boost time')
            require(generation['ranks'] == previous['ranks'], 'Checkpoint ancestry changes MPI rank count')
        for raw in generation.pop('outputs'):
            path = absolute_file(raw['file'], 'Diagnostic source file')
            key = (identity(path), raw['group'])
            require(key not in seen_outputs, 'Duplicate diagnostic file/group prefix in ancestry: ' + str(path))
            seen_outputs.add(key)
            prefix = dict(raw, file=str(path), generation=index,
                          checkpoint=generation['checkpoint'], source_run_id=generation['source_run_id'],
                          reset_diagnostics=generation['reset_diagnostics'])
            if verify_sources:
                prefix['source_validation'] = verify_prefix(prefix, generation)
            slices.append(prefix)
        previous = generation
    return dict(format='aprl-checkpoint-diagnostic-segments-v1', checkpoint=requested_checkpoint,
                sources_verified=verify_sources, generations=generations, slices=slices,
                interpretation='Oldest generation first; each slice is [0:records], with a separate validation prefix. '
                               'Segments are listed, not merged; reset_diagnostics may change plane definitions or cadence. '
                               'Only history through the requested snapshot is included; later resumed output is not discovered.')


def self_test(root):
    """Read existing small acceptance output; failures are injected only in memory."""
    root = Path(root).resolve(strict=True)
    two_generation = root / 'gaussian-cpml/restart-7/checkpoints/step-000000000017.h5'
    validation = root / 'ballistic-validation/continuous/checkpoints/step-000000000007.h5'
    chain = build_manifest(two_generation, True)
    require([g['step'] for g in chain['generations']] == [7, 17], 'Expected the saved two-generation 7 -> 17 chain')
    require({p['rank'] for p in chain['slices']} == {0, 1, 2, 3}, 'Expected all four rank prefix lists')
    reference = build_manifest(validation, True)
    references = [p for p in reference['slices'] if p['group'] == '/ballistic_reference']
    require(len(references) == 1 and references[0]['records'] > references[0]['validation_records'],
            'Validation fixture must straddle reference entrance and exit')
    original = read_checkpoint
    path = two_generation.resolve()
    fixture = original(path)
    def rejected(mutator, expected):
        altered = copy.deepcopy(fixture)
        mutator(altered)
        def injected(target):
            return copy.deepcopy(altered) if target == path else original(target)
        with mock.patch.object(sys.modules[__name__], 'read_checkpoint', injected):
            try:
                build_manifest(path, True)
            except (PrefixError, FileNotFoundError) as error:
                require(expected in str(error), 'Unexpected rejection: ' + str(error))
            else:
                raise AssertionError('Malformed ancestry was accepted: ' + expected)
    rejected(lambda g: g.update(parent_checkpoint=str(path)), 'cycle')
    rejected(lambda g: g['outputs'].append(copy.deepcopy(g['outputs'][0])), 'Duplicate')
    missing = root / '__checkpoint_segments_self_test_missing__' / 'absent.h5'
    require(not missing.exists(), 'Missing-file test path unexpectedly exists')
    rejected(lambda g: g['outputs'][0].update(file=str(missing)), 'missing')
    rejected(lambda g: g['outputs'][0].update(records=2**63), 'below the frozen prefix')
    return dict(status='PASS', tests=6, two_generation_checkpoint=str(two_generation),
                validation_checkpoint=str(validation), generations=len(chain['generations']),
                slices=len(chain['slices']), scope='Metadata/commit/shape reads plus in-memory rejection fixtures; no writes or simulation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path, nargs='?')
    parser.add_argument('--verify-sources', action='store_true',
                        help='Also verify source run identity, commit markers and dataset shapes')
    parser.add_argument('--max-depth', type=int, default=256)
    parser.add_argument('--self-test', type=Path, metavar='SCIENCE_RESULT_ROOT',
                        help='Read existing science-r1 two-generation and ballistic-validation fixtures')
    args = parser.parse_args()
    if args.self_test and args.checkpoint:
        parser.error('Use checkpoint or --self-test, not both')
    if not args.self_test and args.checkpoint is None:
        parser.error('Supply a checkpoint or --self-test')
    try:
        result = self_test(args.self_test) if args.self_test else build_manifest(
            args.checkpoint, args.verify_sources, args.max_depth)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, UnicodeError, yaml.YAMLError) as error:
        print('Checkpoint prefix error: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
