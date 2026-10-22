#!/usr/bin/env python3
"""Compare APRL particle-plane crossings with inertial input-plane tracks.

This is intentionally specific to the independent uniform-bunch generator:
its HDF5-v4 source IDs must equal their one-based input row numbers.  That
contract permits a bounded-memory join of detector records, which arrive in
crossing-time order, without loading the full input bunch into RAM.  A
temporary disk-backed array holds four predicted values per source particle.
"""

import argparse
import json
import math
import tempfile
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "This analysis needs h5py and numpy in the selected Python environment"
    ) from exc


METRICS = (
    "dx_um", "dy_um", "dr_um", "dxp_urad", "dyp_urad"
)


class StreamingMetric:
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
        self.square_total += float(np.sum(values * values, dtype=np.float64))
        if np.any(sample_mask):
            self.sample.append(np.abs(values[sample_mask]))

    def report(self):
        if not self.count:
            return None
        sampled = np.concatenate(self.sample) if self.sample else np.empty(0)
        quantiles = np.quantile(sampled, [0.5, 0.9, 0.99, 0.999]) if len(
            sampled) else [None] * 4
        return {
            "count": self.count,
            "min": self.minimum,
            "max": self.maximum,
            "max_abs": max(abs(self.minimum), abs(self.maximum)),
            "mean": self.total / self.count,
            "rms": math.sqrt(self.square_total / self.count),
            "abs_quantile_sample_count": int(len(sampled)),
            "abs_p50": float(quantiles[0]) if len(sampled) else None,
            "abs_p90": float(quantiles[1]) if len(sampled) else None,
            "abs_p99": float(quantiles[2]) if len(sampled) else None,
            "abs_p99_9": float(quantiles[3]) if len(sampled) else None,
        }


def require_fields(dataset, required, name):
    fields = dataset.dtype.fields or {}
    absent = sorted(set(required) - set(fields))
    if absent:
        raise ValueError(f"{name} is missing fields: {', '.join(absent)}")


def build_prediction(input_path, predicted, input_plane_z, target_z,
                     offset_x, offset_y, chunk):
    with h5py.File(input_path, "r") as source:
        group = source["/particles"]
        version = int(group.attrs["format_version"])
        if version != 4:
            raise ValueError(f"expected input HDF5 format_version=4; got {version}")
        if int(group.attrs.get("complete", 0)) != 1:
            raise ValueError("input HDF5 particle file is incomplete")
        records = group["records"]
        require_fields(records, ("plane_position_m", "proper_velocity",
                                 "source_id"), "/particles/records")
        if len(records) != len(predicted):
            raise ValueError("input length changed while building predictions")
        distance = target_z - input_plane_z
        for start in range(0, len(records), chunk):
            stop = min(start + chunk, len(records))
            block = records[start:stop]
            source_ids = block["source_id"]
            expected = np.arange(start + 1, stop + 1, dtype=np.uint64)
            if not np.array_equal(source_ids, expected):
                raise ValueError(
                    "this low-memory analyzer requires generator source_id "
                    "= one-based input row; source IDs are not sequential "
                    f"in rows {start + 1}..{stop}"
                )
            position = block["plane_position_m"]
            velocity = block["proper_velocity"]
            if position.shape != (stop - start, 2) or velocity.shape != (
                    stop - start, 3):
                raise ValueError("unexpected input coordinate vector shape")
            if (not np.all(np.isfinite(position)) or
                    not np.all(np.isfinite(velocity)) or
                    not np.all(velocity[:, 2] > 0)):
                raise ValueError("input has non-finite values or nonpositive uz")
            slope_x = velocity[:, 0] / velocity[:, 2]
            slope_y = velocity[:, 1] / velocity[:, 2]
            predicted[start:stop, 0] = position[:, 0] + offset_x + distance * slope_x
            predicted[start:stop, 1] = position[:, 1] + offset_y + distance * slope_y
            predicted[start:stop, 2] = slope_x
            predicted[start:stop, 3] = slope_y
    predicted.flush()


def inspect_plane(plane_path, predicted, target_z, chunk, sample_stride,
                  scratch_path, allow_incomplete):
    seen = np.memmap(scratch_path, dtype=np.uint8, mode="w+",
                     shape=(len(predicted),))
    seen[:] = 0
    metrics = {name: StreamingMetric() for name in METRICS}
    upstream = duplicate_downstream = 0
    total_downstream = 0
    with h5py.File(plane_path, "r") as source:
        group = source["/particle_plane"]
        plane_z = float(group.attrs["plane_z_m"])
        if not math.isclose(plane_z, target_z, rel_tol=0, abs_tol=1e-9):
            raise ValueError(
                f"{plane_path}: plane is at {plane_z:.12g} m, "
                f"not requested {target_z:.12g} m"
            )
        complete = bool(group["complete"][()])
        if not complete and not allow_incomplete:
            raise ValueError(
                f"{plane_path}: detector file is incomplete; "
                "pass --allow-incomplete only for diagnostics"
            )
        committed = int(group["committed_records"][()])
        records = group["records"]
        if committed < 0 or committed > len(records):
            raise ValueError(f"{plane_path}: invalid committed_records")
        require_fields(records, ("source_id", "position_m",
                                 "proper_velocity", "crossing_direction"),
                       "/particle_plane/records")
        for start in range(0, committed, chunk):
            stop = min(start + chunk, committed)
            block = records[start:stop]
            direction = block["crossing_direction"]
            if np.any((direction != 1) & (direction != -1)):
                raise ValueError(f"{plane_path}: invalid crossing direction")
            upstream += int(np.count_nonzero(direction == -1))
            forward = block[direction == 1]
            if not len(forward):
                continue
            total_downstream += len(forward)
            ids = np.asarray(forward["source_id"], dtype=np.uint64)
            if np.any(ids < 1) or np.any(ids > len(predicted)):
                raise ValueError(f"{plane_path}: source_id outside input range")
            indices = ids - 1
            # Preserve the first downstream crossing of each source ID.
            _, first = np.unique(indices, return_index=True)
            first.sort()
            duplicate_downstream += len(indices) - len(first)
            indices = indices[first]
            forward = forward[first]
            not_seen = seen[indices] == 0
            duplicate_downstream += int(np.count_nonzero(~not_seen))
            indices = indices[not_seen]
            forward = forward[not_seen]
            if not len(indices):
                continue
            seen[indices] = 1
            expected = predicted[indices]
            position = forward["position_m"]
            velocity = forward["proper_velocity"]
            if position.shape != (len(indices), 3) or velocity.shape != (
                    len(indices), 3):
                raise ValueError("unexpected detector coordinate vector shape")
            if (not np.all(np.isfinite(position)) or
                    not np.all(np.isfinite(velocity)) or
                    not np.all(velocity[:, 2] > 0)):
                raise ValueError("detector has non-finite values or nonpositive uz")
            if not np.all(np.abs(position[:, 2] - target_z) <= 1e-9):
                raise ValueError("detector record position does not match its plane")
            dx = (position[:, 0] - expected[:, 0]) * 1e6
            dy = (position[:, 1] - expected[:, 1]) * 1e6
            values = {
                "dx_um": dx,
                "dy_um": dy,
                "dr_um": np.hypot(dx, dy),
                "dxp_urad": (velocity[:, 0] / velocity[:, 2] - expected[:, 2]) * 1e6,
                "dyp_urad": (velocity[:, 1] / velocity[:, 2] - expected[:, 3]) * 1e6,
            }
            sample_mask = ids[first][not_seen] % sample_stride == 0
            for name, value in values.items():
                if not np.all(np.isfinite(value)):
                    raise ValueError(f"{plane_path}: non-finite {name}")
                metrics[name].add(value, sample_mask)
    unique_forward = int(np.count_nonzero(seen))
    del seen
    return {
        "plane_file": str(plane_path.resolve()),
        "plane_z_m": plane_z,
        "complete": complete,
        "committed_records": committed,
        "downstream_records": total_downstream,
        "unique_downstream_particles": unique_forward,
        "duplicate_downstream_records": duplicate_downstream,
        "upstream_records": upstream,
        "missing_input_particles": len(predicted) - unique_forward,
        "metrics": {name: metric.report() for name, metric in metrics.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="generator HDF5-v4 particle input")
    parser.add_argument("--plane", required=True, action="append", type=Path,
                        help="APRL particle detector HDF5; repeat for cases")
    parser.add_argument("--target-z-m", type=float,
                        help="defaults to the first detector file's plane_z_m; "
                             "current 1 m-from-centre test uses 1.000010 m")
    parser.add_argument("--input-plane-z-m", type=float, default=0.0)
    parser.add_argument("--offset-x-m", type=float, default=0.0)
    parser.add_argument("--offset-y-m", type=float, default=0.0)
    parser.add_argument("--chunk", type=int, default=65536)
    parser.add_argument("--quantile-sample", type=int, default=250000,
                        help="deterministic source-ID sample; quantiles approximate")
    parser.add_argument("--scratch-dir", type=Path,
                        help="disk with space for about 33 bytes/input particle")
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    if args.chunk <= 0 or args.quantile_sample <= 0:
        parser.error("--chunk and --quantile-sample must be positive")
    if not args.input.is_file() or any(not path.is_file() for path in args.plane):
        parser.error("input and all particle-plane files must exist")
    target_z = args.target_z_m
    if target_z is None:
        with h5py.File(args.plane[0], "r") as first_plane:
            target_z = float(first_plane["/particle_plane"].attrs["plane_z_m"])
    if not all(map(math.isfinite, (target_z, args.input_plane_z_m,
                                 args.offset_x_m, args.offset_y_m))):
        parser.error("coordinates must be finite")
    if target_z <= args.input_plane_z_m:
        parser.error("target must be downstream of input plane")
    with h5py.File(args.input, "r") as source:
        group = source["/particles"]
        if int(group.attrs.get("complete", 0)) != 1:
            parser.error("input /particles is incomplete or lacks complete=1")
        count = len(group["records"])
    if count < 1:
        parser.error("input particle dataset is empty")
    stride = max(1, math.ceil(count / args.quantile_sample))
    scratch_root = str(args.scratch_dir) if args.scratch_dir else None
    with tempfile.TemporaryDirectory(prefix="aprl-plane-analysis-",
                                     dir=scratch_root) as temporary:
        path = Path(temporary)
        predicted = np.memmap(path / "predicted.bin", dtype=np.float64,
                              mode="w+", shape=(count, 4))
        build_prediction(args.input, predicted, args.input_plane_z_m,
                         target_z, args.offset_x_m, args.offset_y_m,
                         args.chunk)
        results = []
        for index, plane in enumerate(args.plane):
            results.append(inspect_plane(
                plane, predicted, target_z, args.chunk, stride,
                path / f"seen-{index}.bin", args.allow_incomplete))
        del predicted
    report = {
        "input_file": str(args.input.resolve()),
        "input_particles": count,
        "input_plane_z_m": args.input_plane_z_m,
        "target_z_m": target_z,
        "input_transverse_offset_m": [args.offset_x_m, args.offset_y_m],
        "quantile_method": "deterministic source_id modulo sample_stride; absolute errors",
        "quantile_sample_stride": stride,
        "units": {"position_error": "um", "slope_error": "urad"},
        "planes": results,
    }
    for result in results:
        print(Path(result["plane_file"]).name + ": "
              f"{result['unique_downstream_particles']:,}/{count:,} "
              "unique downstream particles; "
              f"duplicates={result['duplicate_downstream_records']:,}, "
              f"upstream={result['upstream_records']:,}")
        for name in METRICS:
            metric = result["metrics"][name]
            if metric:
                p99 = (f"{metric['abs_p99']:.6g}" if
                       metric["abs_p99"] is not None else "n/a")
                p999 = (f"{metric['abs_p99_9']:.6g}" if
                        metric["abs_p99_9"] is not None else "n/a")
                print(f"  {name}: mean={metric['mean']:.6g}, "
                      f"RMS={metric['rms']:.6g}, "
                      f"|error| P99≈{p99}, P99.9≈{p999}")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2,
                                            allow_nan=False) + "\n",
                                 encoding="utf-8")


if __name__ == "__main__":
    main()
