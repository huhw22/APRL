#!/usr/bin/env python3
"""Compare two fixed-laboratory particle planes by source_id.

The two HDF5 files may list crossings in different orders.  Only the first
downstream crossing of each source_id is compared.  A disk-backed table keeps
memory bounded for the roughly ten million particles in the 1 m drift test.
"""

import argparse
import json
import math
import shutil
import tempfile
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    raise SystemExit("This comparison needs h5py and numpy") from exc


FIELDS = ("source_id", "position_m", "proper_velocity", "crossing_direction")
ENDPOINT_METRICS = ("x_um", "y_um", "r_um", "xp_urad", "yp_urad")
DIFFERENCE_METRICS = (
    "dx_um", "dy_um", "dr_um", "dxy_um", "dxp_urad", "dyp_urad"
)


class Metric:
    def __init__(self):
        self.count = 0
        self.minimum = math.inf
        self.maximum = -math.inf
        self.total = 0.0
        self.square_total = 0.0
        self.sample = []

    def add(self, values, sample_mask):
        if not len(values):
            return
        self.count += len(values)
        self.minimum = min(self.minimum, float(np.min(values)))
        self.maximum = max(self.maximum, float(np.max(values)))
        self.total += float(np.sum(values, dtype=np.float64))
        self.square_total += float(np.sum(np.square(values), dtype=np.float64))
        if np.any(sample_mask):
            self.sample.append(np.abs(values[sample_mask]))

    def report(self):
        if not self.count:
            return None
        sample = np.concatenate(self.sample) if self.sample else np.empty(0)
        quantiles = np.quantile(sample, [0.5, 0.9, 0.99, 0.999]) if len(
            sample) else None
        return {
            "count": self.count,
            "min": self.minimum,
            "max": self.maximum,
            "mean": self.total / self.count,
            "rms": math.sqrt(self.square_total / self.count),
            "abs_quantile_sample_count": int(len(sample)),
            "abs_p50": float(quantiles[0]) if quantiles is not None else None,
            "abs_p90": float(quantiles[1]) if quantiles is not None else None,
            "abs_p99": float(quantiles[2]) if quantiles is not None else None,
            "abs_p99_9": float(quantiles[3]) if quantiles is not None else None,
        }


def plane_header(path, h5file, allow_incomplete):
    group = h5file["/particle_plane"]
    records = group["records"]
    missing = sorted(set(FIELDS) - set(records.dtype.fields or {}))
    if missing:
        raise ValueError(f"{path}: missing particle-plane fields: {missing}")
    complete = bool(group["complete"][()])
    if not complete and not allow_incomplete:
        raise ValueError(f"{path}: incomplete detector; use --allow-incomplete "
                         "only for diagnostic comparison")
    committed = int(group["committed_records"][()])
    if not 0 <= committed <= len(records):
        raise ValueError(f"{path}: invalid committed_records={committed}")
    z = float(group.attrs["plane_z_m"])
    if not math.isfinite(z):
        raise ValueError(f"{path}: non-finite plane_z_m")
    return {
        "path": str(path.resolve()), "plane_z_m": z,
        "complete": complete, "committed_records": committed,
        "records": records,
    }


def blocks(records, committed, chunk):
    for start in range(0, committed, chunk):
        yield records[start:min(start + chunk, committed)]


def maximum_source_id(headers, chunk):
    maximum = 0
    for header in headers:
        for block in blocks(header["records"], header["committed_records"], chunk):
            ids = block["source_id"]
            if np.any(ids == 0):
                raise ValueError(f"{header['path']}: source_id zero is invalid")
            maximum = max(maximum, int(np.max(ids)))
    return maximum


def new_counts(header):
    return {
        "file": header["path"], "complete": header["complete"],
        "committed_records": header["committed_records"],
        "downstream_records": 0, "upstream_records": 0,
        "duplicate_downstream_records": 0, "unique_downstream_particles": 0,
    }


def first_downstream(block, seen, maximum, plane_z, counts, filename):
    direction = block["crossing_direction"]
    if np.any((direction != 1) & (direction != -1)):
        raise ValueError(f"{filename}: crossing_direction must be +1 or -1")
    counts["upstream_records"] += int(np.count_nonzero(direction == -1))
    forward = block[direction == 1]
    counts["downstream_records"] += len(forward)
    if not len(forward):
        return np.empty(0, dtype=np.int64), forward
    ids = np.asarray(forward["source_id"], dtype=np.uint64)
    if np.any(ids == 0) or np.any(ids > maximum):
        raise ValueError(f"{filename}: invalid source_id")
    position = forward["position_m"]
    velocity = forward["proper_velocity"]
    if position.shape != (len(forward), 3) or velocity.shape != (len(forward), 3):
        raise ValueError(f"{filename}: invalid position/velocity shape")
    if (not np.all(np.isfinite(position)) or
            not np.all(np.isfinite(velocity)) or
            not np.all(velocity[:, 2] > 0) or
            not np.all(np.abs(position[:, 2] - plane_z) <= 1e-9)):
        raise ValueError(f"{filename}: invalid downstream phase-space record")
    _, first = np.unique(ids, return_index=True)
    first.sort()  # Keep the first event, not sorted source-ID order.
    counts["duplicate_downstream_records"] += len(ids) - len(first)
    ids = ids[first]
    forward = forward[first]
    unobserved = seen[ids] == 0
    counts["duplicate_downstream_records"] += int(np.count_nonzero(~unobserved))
    ids = ids[unobserved].astype(np.int64, copy=False)
    forward = forward[unobserved]
    seen[ids] = 1
    counts["unique_downstream_particles"] += len(ids)
    return ids, forward


def endpoint_values(position, velocity):
    x = position[:, 0] * 1e6
    y = position[:, 1] * 1e6
    return {
        "x_um": x,
        "y_um": y,
        "r_um": np.hypot(x, y),
        "xp_urad": velocity[:, 0] / velocity[:, 2] * 1e6,
        "yp_urad": velocity[:, 1] / velocity[:, 2] * 1e6,
    }


def compare(reference, variant, scratch, chunk, quantile_sample,
            max_scratch_gib):
    maximum = maximum_source_id((reference, variant), chunk)
    if maximum == 0:
        raise ValueError("neither file contains particle crossings")
    # These detectors come from the one-based uniform-bunch generator.  This
    # guards against a corrupt, extremely sparse ID causing huge allocation.
    limit = max(50_000_000, 4 * max(reference["committed_records"],
                                   variant["committed_records"]))
    if maximum > limit:
        raise ValueError(f"maximum source_id {maximum} is implausibly sparse "
                         f"for these files (limit {limit})")
    scratch_bytes = (maximum + 1) * (4 * 8 + 2)
    if scratch_bytes > max_scratch_gib * 1024 ** 3:
        raise ValueError(f"disk-backed join needs {scratch_bytes / 1024**3:.3g} "
                         "GiB; increase --max-scratch-gib deliberately")
    if shutil.disk_usage(scratch).free < int(scratch_bytes * 1.1):
        raise ValueError("insufficient free space in --scratch-dir")
    stride = max(1, math.ceil(maximum / quantile_sample))
    with tempfile.TemporaryDirectory(prefix="aprl-paired-planes-", dir=scratch) as tmp:
        path = Path(tmp)
        ref_values = np.memmap(path / "reference.bin", dtype=np.float64,
                               mode="w+", shape=(maximum + 1, 4))
        ref_seen = np.memmap(path / "reference-seen.bin", dtype=np.uint8,
                             mode="w+", shape=(maximum + 1,))
        var_seen = np.memmap(path / "variant-seen.bin", dtype=np.uint8,
                             mode="w+", shape=(maximum + 1,))
        ref_seen[:] = 0
        var_seen[:] = 0
        ref_counts = new_counts(reference)
        var_counts = new_counts(variant)
        for block in blocks(reference["records"], reference["committed_records"],
                            chunk):
            ids, first = first_downstream(block, ref_seen, maximum,
                                          reference["plane_z_m"], ref_counts,
                                          reference["path"])
            if len(ids):
                position = first["position_m"]
                velocity = first["proper_velocity"]
                ref_values[ids, :2] = position[:, :2]
                ref_values[ids, 2:] = velocity[:, :2] / velocity[:, 2, None]

        metric_groups = {
            "reference": {name: Metric() for name in ENDPOINT_METRICS},
            "variant": {name: Metric() for name in ENDPOINT_METRICS},
            "difference": {name: Metric() for name in DIFFERENCE_METRICS},
        }
        matched = 0
        for block in blocks(variant["records"], variant["committed_records"],
                            chunk):
            ids, first = first_downstream(block, var_seen, maximum,
                                          variant["plane_z_m"], var_counts,
                                          variant["path"])
            if not len(ids):
                continue
            common = ref_seen[ids] != 0
            ids = ids[common]
            first = first[common]
            if not len(ids):
                continue
            matched += len(ids)
            baseline = ref_values[ids]
            variant_position = first["position_m"]
            variant_velocity = first["proper_velocity"]
            baseline_position = np.empty((len(ids), 3), dtype=np.float64)
            baseline_position[:, :2] = baseline[:, :2]
            baseline_position[:, 2] = reference["plane_z_m"]
            baseline_velocity = np.empty((len(ids), 3), dtype=np.float64)
            baseline_velocity[:, :2] = baseline[:, 2:]
            baseline_velocity[:, 2] = 1.0
            a = endpoint_values(baseline_position, baseline_velocity)
            b = endpoint_values(variant_position, variant_velocity)
            dx = b["x_um"] - a["x_um"]
            dy = b["y_um"] - a["y_um"]
            delta = {
                "dx_um": dx,
                "dy_um": dy,
                "dr_um": b["r_um"] - a["r_um"],
                "dxy_um": np.hypot(dx, dy),
                "dxp_urad": b["xp_urad"] - a["xp_urad"],
                "dyp_urad": b["yp_urad"] - a["yp_urad"],
            }
            sampled = ids % stride == 0
            for group_name, values in (("reference", a), ("variant", b),
                                       ("difference", delta)):
                for name, value in values.items():
                    metric_groups[group_name][name].add(value, sampled)
        reference_only = int(np.count_nonzero((ref_seen != 0) & (var_seen == 0)))
        variant_only = int(np.count_nonzero((var_seen != 0) & (ref_seen == 0)))
        del ref_values, ref_seen, var_seen
    return {
        "plane_z_m": reference["plane_z_m"],
        "maximum_source_id": maximum,
        "scratch_bytes_required": scratch_bytes,
        "quantile_method": "absolute-value quantiles for source_id modulo stride = 0",
        "quantile_sample_stride": stride,
        "reference": ref_counts,
        "variant": var_counts,
        "matching": {
            "matched_particles": matched,
            "reference_only_particles": reference_only,
            "variant_only_particles": variant_only,
            "same_unique_population": reference_only == 0 and variant_only == 0,
            "input_completeness_verified": False,
            "input_completeness_note": (
                "Detector-only comparison cannot detect IDs absent from both "
                "files; check each detector against its source particles.h5 "
                "with analyze_particle_plane.py. A common source_id also only "
                "identifies the same physical input particle if both runs "
                "used identical generated bunches."
            ),
        },
        "metrics": {
            group: {name: metric.report() for name, metric in entries.items()}
            for group, entries in metric_groups.items()
        },
        "metric_definition": {
            "reference_variant": "matched first downstream crossings at the same lab z",
            "dr_um": "variant radial coordinate minus reference radial coordinate",
            "dxy_um": "magnitude of transverse variant-reference displacement",
            "slope": "proper_velocity_x_or_y / proper_velocity_z",
            "p99": "approximately sampled P99 of absolute value",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--variant", required=True, type=Path)
    parser.add_argument("--scratch-dir", required=True, type=Path)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument("--chunk", type=int, default=65536)
    parser.add_argument("--quantile-sample", type=int, default=250000)
    parser.add_argument("--max-scratch-gib", type=float, default=5.0)
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="diagnostic only; never call a partial comparison final")
    args = parser.parse_args()
    if args.reference.resolve() == args.variant.resolve():
        parser.error("reference and variant must be distinct files")
    if args.chunk < 1 or args.quantile_sample < 1 or not math.isfinite(
            args.max_scratch_gib) or args.max_scratch_gib <= 0:
        parser.error("chunk, quantile-sample and max-scratch-gib must be positive")
    if not args.scratch_dir.is_dir():
        parser.error("--scratch-dir must already exist")
    if not args.reference.is_file() or not args.variant.is_file():
        parser.error("both particle-plane files must exist")
    with h5py.File(args.reference, "r") as reference_file, h5py.File(
            args.variant, "r") as variant_file:
        reference = plane_header(args.reference, reference_file,
                                 args.allow_incomplete)
        variant = plane_header(args.variant, variant_file,
                               args.allow_incomplete)
        if not math.isclose(reference["plane_z_m"], variant["plane_z_m"],
                            rel_tol=0, abs_tol=1e-9):
            raise ValueError("particle planes are at different lab z positions")
        report = compare(reference, variant, args.scratch_dir, args.chunk,
                         args.quantile_sample, args.max_scratch_gib)
    match = report["matching"]
    print("matched={:,}, reference-only={:,}, variant-only={:,}".format(
        match["matched_particles"], match["reference_only_particles"],
        match["variant_only_particles"]))
    for name in DIFFERENCE_METRICS:
        value = report["metrics"]["difference"][name]
        if value is None:
            continue
        print("  {}: mean={:.6g}, RMS={:.6g}, |value| P99≈{}".format(
            name, value["mean"], value["rms"],
            "n/a" if value["abs_p99"] is None else
            "{:.6g}".format(value["abs_p99"])))
    for group in ("reference", "variant"):
        print(f"{group} endpoint:")
        for name in ENDPOINT_METRICS:
            value = report["metrics"][group][name]
            if value is None:
                continue
            print("  {}: mean={:.6g}, RMS={:.6g}, |value| P99≈{}".format(
                name, value["mean"], value["rms"],
                "n/a" if value["abs_p99"] is None else
                "{:.6g}".format(value["abs_p99"])))
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                             encoding="utf-8")
    print(f"JSON: {args.json_out.resolve()}")


if __name__ == "__main__":
    main()
