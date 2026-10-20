#!/usr/bin/env python3
"""Select a central longitudinal window from an Elegant-derived APRL HDF5 v4 file.

The source is read in bounded blocks.  Its arrival times, IDs, weights and all
phase-space coordinates are copied verbatim for retained particles.
"""

import argparse
import math
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np


LIGHT_SPEED_M_S = 299_792_458.0
ELEMENTARY_CHARGE_C = 1.602_176_634e-19


def _positive_finite(text):
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return value


def _finite(text):
    value = float(text)
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError("must be finite")
    return value


def _positive_int(text):
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def _check_record_type(records):
    if records.ndim != 1 or len(records) == 0:
        raise ValueError("/particles/records must be a nonempty 1-D dataset")
    dtype = records.dtype
    required = {
        "plane_position_m": ("f", (2,)),
        "arrival_time_offset_s": ("f", ()),
        "proper_velocity": ("f", (3,)),
        "source_id": ("u", ()),
        "macro_weight": ("f", ()),
    }
    if dtype.names is None or not required.keys() <= set(dtype.names):
        raise ValueError("/particles/records is missing required HDF5 v4 fields")
    for name, (kind, shape) in required.items():
        field = dtype.fields[name][0]
        base, actual_shape = field.subdtype or (field, ())
        if base.kind != kind or actual_shape != shape:
            raise ValueError(f"/particles/records field {name} has the wrong type or shape")


def crop(input_path, output_path, width_um=60.0, center_time_offset_fs=0.0,
         block_records=262_144):
    """Copy records with |beta_z*c*(tau-center)| <= width/2; return statistics."""
    input_path = Path(input_path).expanduser().resolve(strict=True)
    output_path = Path(output_path).expanduser().absolute()
    if output_path.exists() or output_path.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    if input_path == output_path.resolve():
        raise ValueError("Input and output must be different files")
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"Output directory does not exist: {output_path.parent}")
    if not math.isfinite(width_um) or width_um <= 0:
        raise ValueError("width_um must be finite and positive")
    if not math.isfinite(center_time_offset_fs):
        raise ValueError("center_time_offset_fs must be finite")
    if block_records <= 0:
        raise ValueError("block_records must be positive")

    half_width_m = width_um * 1e-6 / 2
    center_s = center_time_offset_fs * 1e-15
    temporary_path = None
    try:
        with h5py.File(input_path, "r") as source:
            if "particles" not in source or "records" not in source["particles"]:
                raise ValueError("Expected /particles/records in the input HDF5 file")
            source_group = source["particles"]
            if source_group.attrs.get("format_version") != 4:
                raise ValueError("Only Elegant-derived HDF5 format_version=4 is supported")
            if set(source.keys()) != {"particles"} or set(source_group.keys()) != {"records"}:
                raise ValueError("Unexpected HDF5 objects; refusing to silently discard them")
            records = source_group["records"]
            _check_record_type(records)

            handle, temporary_name = tempfile.mkstemp(
                prefix=output_path.name + ".", suffix=".partial", dir=output_path.parent)
            os.close(handle)
            temporary_path = Path(temporary_name)
            scanned_weight = 0.0
            kept_weight = 0.0
            all_unit_weights = True
            kept = 0
            total = len(records)
            next_progress = max(total // 10, 1)

            with h5py.File(temporary_path, "w") as target:
                for name, value in source.attrs.items():
                    target.attrs[name] = value
                target_group = target.create_group("particles")
                for name, value in source_group.attrs.items():
                    target_group.attrs[name] = value
                output = target_group.create_dataset(
                    "records", shape=(0,), maxshape=(None,), dtype=records.dtype,
                    chunks=(min(block_records, 65_536),))
                for start in range(0, total, block_records):
                    block = records[start:start + block_records]
                    tau = block["arrival_time_offset_s"]
                    velocity = block["proper_velocity"]
                    weights = block["macro_weight"]
                    if (not np.isfinite(tau).all() or not np.isfinite(velocity).all()
                            or not np.isfinite(block["plane_position_m"]).all()
                            or not np.isfinite(weights).all() or np.any(weights <= 0)
                            or np.any(velocity[:, 2] <= 0)):
                        raise ValueError(f"Invalid particle values near input row {start + 1}")
                    gamma = np.sqrt(1.0 + np.einsum("ij,ij->i", velocity, velocity))
                    if not np.isfinite(gamma).all():
                        raise ValueError(f"Invalid proper velocity near input row {start + 1}")
                    beta_z = velocity[:, 2] / gamma
                    longitudinal_m = -beta_z * LIGHT_SPEED_M_S * (tau - center_s)
                    # Admit values on the mathematical boundary despite the
                    # round-off in converting u -> beta and seconds -> metres.
                    tolerance_m = 8 * np.spacing(half_width_m)
                    mask = np.abs(longitudinal_m) <= half_width_m + tolerance_m
                    chosen = block[mask]
                    if len(chosen):
                        output.resize((kept + len(chosen),))
                        output[kept:kept + len(chosen)] = chosen
                        kept += len(chosen)
                        kept_weight += float(np.sum(weights[mask], dtype=np.longdouble))
                    scanned_weight += float(np.sum(weights, dtype=np.longdouble))
                    all_unit_weights = all_unit_weights and bool(np.all(weights == 1.0))
                    scanned = min(start + block_records, total)
                    if scanned >= next_progress or scanned == total:
                        print(f"  scanned {scanned:,}/{total:,}; kept {kept:,}", flush=True)
                        next_progress = scanned + max(total // 10, 1)

                if kept == 0:
                    raise ValueError("The requested window contains no particles")
                target_group.attrs["crop_method"] = "-beta_z*c*(tau-center), inclusive"
                target_group.attrs["crop_width_m"] = 2 * half_width_m
                target_group.attrs["crop_center_arrival_time_offset_s"] = center_s
                target_group.attrs["crop_source_records"] = total
                target_group.attrs["crop_selected_records"] = kept
                target_group.attrs["crop_source_weight_sum"] = scanned_weight
                target_group.attrs["crop_selected_weight_sum"] = kept_weight
                target_group.attrs["crop_source_file"] = str(input_path)
                target.flush()

        # A hard link publishes the completed file atomically without replacing
        # an output that another process may have created in the meantime.
        os.link(temporary_path, output_path)
        temporary_path.unlink()
        temporary_path = None
        return total, kept, scanned_weight, kept_weight, all_unit_weights
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="original Elegant HDF5 v4")
    parser.add_argument("--output", required=True, type=Path, help="new cropped HDF5 v4")
    parser.add_argument("--width-um", type=_positive_finite, default=60.0,
                        help="full longitudinal window in micrometres (default: 60)")
    parser.add_argument("--center-time-offset-fs", type=_finite, default=0.0,
                        help="window centre relative to the source arrival-time reference (default: 0)")
    parser.add_argument("--block-records", type=_positive_int, default=262_144,
                        help="maximum records held in each input block")
    args = parser.parse_args(argv)
    try:
        total, kept, source_weight, kept_weight, all_unit_weights = crop(
            args.input, args.output, args.width_um,
            args.center_time_offset_fs, args.block_records)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(f"Wrote {kept:,} / {total:,} particles to {args.output}")
    print(f"Selected weight fraction: {kept_weight / source_weight:.12g}")
    print("All retained records and the original arrival-time reference are unchanged.")
    if all_unit_weights:
        print("All relative macro weights are 1. If the upstream beam used "
              "one physical electron per macro-particle, this selection represents "
              f"{kept:,} electrons ({kept * ELEMENTARY_CHARGE_C:.12g} C).")
    else:
        print("For non-unit weights, set APRL beam.input.electrons to the "
              "cropped physical total; do not retain the full-bunch value.")


if __name__ == "__main__":
    main()
