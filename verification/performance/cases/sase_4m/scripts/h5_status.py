#!/usr/bin/env python3
"""Return success only when an APRL HDF5 product carries complete=1."""

import sys
from pathlib import Path

import h5py
import numpy as np


def scalar(value: object) -> int:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError("completion marker is not scalar")
    return int(array.reshape(-1)[0])


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: h5_status.py FILE GROUP", file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    if not path.is_file():
        print("HDF5 product is missing: %s" % path, file=sys.stderr)
        return 1
    try:
        with h5py.File(path, "r") as handle:
            group_name = sys.argv[2]
            if group_name not in handle:
                print(
                    "HDF5 product %s has no group %s" % (path, group_name),
                    file=sys.stderr,
                )
                return 1
            group = handle[group_name]
            if "complete" in group:
                value = scalar(group["complete"][()])
            elif "complete" in group.attrs:
                value = scalar(group.attrs["complete"])
            else:
                print(
                    "HDF5 group %s:%s has no completion marker"
                    % (path, group_name),
                    file=sys.stderr,
                )
                return 1
    except Exception as error:
        print("HDF5 status check failed: %s" % error, file=sys.stderr)
        return 1
    if value != 1:
        print(
            "HDF5 group %s:%s is incomplete (complete=%s)"
            % (path, sys.argv[2], value),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
