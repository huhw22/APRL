#!/usr/bin/env python3
"""Read only HDF5 metadata and small energy ledgers; never read particle/field arrays."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np


def native(value):
    if isinstance(value, bytes):
        return value.decode(errors='replace')
    if isinstance(value, np.ndarray):
        return [native(x) for x in value.tolist()]
    if isinstance(value, np.generic):
        return native(value.item())
    if isinstance(value, (list, tuple)):
        return [native(x) for x in value]
    return value


def inspect(root):
    report = {'source': str(root.resolve()), 'files': []}
    for path in sorted(root.rglob('*.h5')):
        item = {'path': str(path.relative_to(root)), 'size_bytes': path.stat().st_size}
        try:
            with h5py.File(path, 'r') as f:
                groups, datasets = {}, {}
                def visit(name, obj):
                    if isinstance(obj, h5py.Group):
                        groups[name] = {k: native(v) for k, v in obj.attrs.items()}
                    elif isinstance(obj, h5py.Dataset):
                        d = {'shape': obj.shape, 'dtype': str(obj.dtype), 'chunks': obj.chunks}
                        if obj.shape == ():
                            d['value'] = native(obj[()])
                        datasets[name] = d
                f.visititems(visit)
                item.update(root_attributes={k: native(v) for k, v in f.attrs.items()},
                            groups=groups, datasets=datasets, readable_metadata=True)
                if 'energy_ledger' in f:
                    group = f['energy_ledger']
                    n = int(group['committed_records'][()])
                    ds = group['records']
                    item['ledger_endpoints'] = [
                        {name: native(ds[index][name]) for name in ds.dtype.names}
                        for index in (0, n-1)] if n else []
        except (OSError, ValueError) as exc:
            item.update(readable_metadata=False, error=str(exc))
        report['files'].append(item)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(inspect(args.source), indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
