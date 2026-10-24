#!/usr/bin/env python3
"""Project observed laboratory crossings to the last direct aperture arrival.

Run the full reference scan on a Slurm node.  The reusable API only needs one
chunk at a time and never modifies its input.  Old front-anchored trajectories
are planning evidence, not a rigorous bound on a new held-start trajectory.
"""

import argparse
import json
import math
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

C = 299792458.0
REFERENCE_RELATIVE = (
    "results/elegant_4m_analysis_20261024/original/runs/"
    "20261021T070618Z-861230/detectors/exit-field-ballistic-reference.h5"
)


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def _finite_positive(name, value):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")


def analyze_reference(reference_file, *, input_plane_z_m,
                      retirement_exit_z_m, detector_z_m,
                      physical_half_x_m, physical_half_y_m,
                      dx_m=5e-6, dy_m=5e-6, chunk_rows=262144,
                      expected_records=None, group_name=None):
    """Return direct-arrival extrema from complete lab crossing records.

    ``proper_velocity`` is laboratory gamma*v/c.  Retarded time is ALWAYS
    t_lab - (z - input_plane_z_m)/c.  No common-time or boost-origin offset is
    added: both are already represented in the stored laboratory event.

    The projected velocity is constant from the reference plane to the taper
    exit.  For a pre-retirement plane, any remaining active push is therefore
    an explicit model approximation.  Transverse CPML entry is counted, not
    silently treated as a successfully completed longitudinal retirement.
    """
    reference_file = Path(reference_file).resolve()
    before = reference_file.stat()
    for name, value in (("physical_half_x_m", physical_half_x_m),
                        ("physical_half_y_m", physical_half_y_m),
                        ("dx_m", dx_m), ("dy_m", dy_m)):
        _finite_positive(name, value)
    if not all(math.isfinite(v) for v in
               (input_plane_z_m, retirement_exit_z_m, detector_z_m)):
        raise ValueError("Plane positions must be finite")
    if chunk_rows <= 0 or int(chunk_rows) != chunk_rows:
        raise ValueError("chunk_rows must be a positive integer")
    if expected_records is not None and expected_records <= 0:
        raise ValueError("expected_records must be positive")
    half = np.array([physical_half_x_m, physical_half_y_m], dtype=float)
    spacing = np.array([dx_m, dy_m], dtype=float)
    pixel_half = half - spacing / 2
    if np.any(pixel_half <= 0):
        raise ValueError("Physical aperture must contain pixel centres")
    pixel_counts = 2 * half / spacing
    if not np.allclose(pixel_counts, np.rint(pixel_counts), rtol=0, atol=1e-9):
        raise ValueError("Physical aperture must contain an integer pixel count")

    extrema = {}
    latest = {key: -np.longdouble(np.inf)
              for key in ("retirement", "boundary", "pixel")}
    witnesses = {}
    outside_reference = outside_retirement = 0
    records_seen = 0
    max_plane_error = 0.0

    def record_range(name, values):
        low, high = float(np.min(values)), float(np.max(values))
        if name not in extrema:
            extrema[name] = [low, high]
        else:
            extrema[name][0] = min(extrema[name][0], low)
            extrema[name][1] = max(extrema[name][1], high)

    with h5py.File(reference_file, "r") as handle:
        if group_name is None:
            groups = [name for name in ("ballistic_reference", "particle_plane")
                      if name in handle]
            if len(groups) != 1:
                raise ValueError("Specify one particle-plane/reference group")
            group_name = groups[0]
        group = handle[group_name]
        frame = _text(group.attrs.get("frame", ""))
        unit = _text(group.attrs.get("proper_velocity_unit", ""))
        if frame != "laboratory" or unit != "gamma*v/c":
            raise ValueError("Reference must store lab events and lab gamma*v/c")
        plane_z = float(group.attrs["plane_z_m"])
        if not (math.isfinite(plane_z) and
                input_plane_z_m < plane_z < retirement_exit_z_m < detector_z_m):
            raise ValueError("Require input < reference < retirement exit < detector")
        data = group["records"]
        records = len(data)
        committed = int(group["committed_records"][()])
        complete = int(group["complete"][()])
        if complete != 1 or committed != records or records == 0:
            raise ValueError("Reference must be complete with all records committed")
        if expected_records is not None and records != expected_records:
            raise ValueError("Reference particle count differs from expected_records")
        required = {"particle_id", "source_id", "time_s", "position_m",
                    "proper_velocity", "charge_C", "mass_kg", "weight",
                    "crossing_direction"}
        if not required.issubset(data.dtype.names or ()):
            raise ValueError("Unexpected reference record schema")
        metadata = {name: _text(group.attrs[name]) for name in
                    ("run_id", "source_revision", "configuration_digest_fnv1a64",
                     "run_manifest", "role") if name in group.attrs}
        flight_z = retirement_exit_z_m - plane_z
        photon_z = detector_z_m - retirement_exit_z_m
        wide_c = np.longdouble(C)
        reference_light_time = (np.longdouble(plane_z) -
                                np.longdouble(input_plane_z_m)) / wide_c
        retirement_light_time = (np.longdouble(retirement_exit_z_m) -
                                 np.longdouble(input_plane_z_m)) / wide_c
        detector_light_time = (np.longdouble(detector_z_m) -
                               np.longdouble(input_plane_z_m)) / wide_c
        plane_tolerance = 128 * np.finfo(float).eps * max(1.0, abs(plane_z))

        for begin in range(0, records, chunk_rows):
            rows = data[begin:begin + chunk_rows]
            for key in required:
                if rows[key].dtype.kind in "fc" and not np.isfinite(rows[key]).all():
                    raise ValueError(f"Nonfinite {key} at row chunk {begin}")
            if (np.any(rows["crossing_direction"] != 1) or
                    np.any(rows["proper_velocity"][:, 2] <= 0) or
                    np.any(rows["mass_kg"] <= 0) or np.any(rows["weight"] <= 0)):
                raise ValueError(f"Invalid direction/velocity/mass/weight at {begin}")
            plane_error = float(np.max(np.abs(rows["position_m"][:, 2] - plane_z)))
            max_plane_error = max(max_plane_error, plane_error)
            if plane_error > plane_tolerance:
                raise ValueError(f"Reference positions disagree with plane_z_m at {begin}")
            u = rows["proper_velocity"]
            gamma = np.sqrt(1 + np.sum(u * u, axis=1))
            beta_z = u[:, 2] / gamma
            if (not np.isfinite(gamma).all() or not np.isfinite(beta_z).all()
                    or np.any(beta_z <= 0) or np.any(beta_z > 1)):
                raise ValueError(f"Invalid derived gamma/velocity at {begin}")
            # Stable 1/beta_z - 1; direct subtraction loses significant digits.
            delay = (1 + u[:, 0] ** 2 + u[:, 1] ** 2) / (
                u[:, 2] * (gamma + u[:, 2]))
            xy = rows["position_m"][:, :2]
            xy_exit = xy + flight_z * u[:, :2] / u[:, 2, None]
            tau_reference = rows["time_s"].astype(np.longdouble) - reference_light_time
            tau_exit = tau_reference + np.longdouble(flight_z) / wide_c * delay
            outside_reference += int(np.count_nonzero(np.any(np.abs(xy) > half, axis=1)))
            outside_retirement += int(np.count_nonzero(np.any(np.abs(xy_exit) > half, axis=1)))
            record_range("reference_lab_time_s", rows["time_s"])
            record_range("reference_retarded_relative_input_s", tau_reference)
            record_range("gamma_lab", gamma)
            record_range("beta_z_lab", beta_z)
            record_range("one_over_beta_z_minus_one", delay)
            record_range("retirement_retarded_relative_input_s", tau_exit)
            record_range("projected_retirement_x_m", xy_exit[:, 0])
            record_range("projected_retirement_y_m", xy_exit[:, 1])
            for key, aperture in (("retirement", None), ("boundary", half),
                                  ("pixel", pixel_half)):
                if aperture is None:
                    geometrical_delay = np.zeros(len(rows))
                    tau_arrival = tau_exit
                    light_time = retirement_light_time
                else:
                    far_xy = np.abs(xy_exit) + aperture
                    rho2 = np.sum(far_xy * far_xy, axis=1)
                    distance = np.sqrt(photon_z * photon_z + rho2)
                    # sqrt(D^2+rho^2)-D evaluated without cancellation.
                    geometrical_delay = rho2 / (distance + photon_z) / C
                    tau_arrival = tau_exit + geometrical_delay.astype(np.longdouble)
                    light_time = detector_light_time
                    record_range(key + "_geometric_delay_s", geometrical_delay)
                    record_range(key + "_retarded_relative_input_s", tau_arrival)
                if not (np.isfinite(tau_arrival).all() and np.isfinite(xy_exit).all()):
                    raise ValueError(f"Nonfinite projection at {begin}")
                index = int(np.argmax(tau_arrival))
                if tau_arrival[index] > latest[key]:
                    latest[key] = tau_arrival[index]
                    row = rows[index]
                    witness = {
                        "row_index": begin + index,
                        "particle_id": int(row["particle_id"]),
                        "source_id": int(row["source_id"]),
                        "reference_lab_time_s": float(row["time_s"]),
                        "reference_position_m": row["position_m"].tolist(),
                        "reference_proper_velocity": row["proper_velocity"].tolist(),
                        "gamma_lab": float(gamma[index]),
                        "one_over_beta_z_minus_one": float(delay[index]),
                        "retirement_position_m": [*xy_exit[index].tolist(), retirement_exit_z_m],
                        "retirement_event_lab_s": float(tau_exit[index] + retirement_light_time),
                        "retirement_retarded_relative_input_s": float(tau_exit[index]),
                        "geometric_delay_s": float(geometrical_delay[index]),
                        "arrival_lab_s": float(tau_arrival[index] + light_time),
                        "retarded_relative_input_s": float(tau_arrival[index]),
                    }
                    if aperture is not None:
                        corner = np.where(xy_exit[index] >= 0, -aperture, aperture)
                        witness["farthest_detector_point_m"] = [*corner.tolist(), detector_z_m]
                    witnesses[key] = witness
            records_seen += len(rows)

    after = reference_file.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Reference changed during scan")
    if records_seen != records:
        raise ValueError("Reference scan did not cover every record")
    aperture_rho2 = float(np.sum((2 * half) ** 2))
    aperture_delay = aperture_rho2 / (
        math.sqrt(photon_z * photon_z + aperture_rho2) + photon_z) / C
    independent_aperture_target = latest["retirement"] + np.longdouble(aperture_delay)
    aperture_target = max(latest["boundary"], independent_aperture_target)
    return {
        "format": "aprl-observed-direct-arrival-v1", "pass": True,
        "source_file": str(reference_file), "source_group": group_name,
        "source_size_bytes": before.st_size, "source_mtime_ns": before.st_mtime_ns,
        "source_metadata": metadata, "records": records,
        "input_plane_z_m": input_plane_z_m, "reference_plane_z_m": plane_z,
        "retirement_exit_z_m": retirement_exit_z_m, "detector_z_m": detector_z_m,
        "latest_retirement_event_lab_s": witnesses["retirement"]["arrival_lab_s"],
        "latest_retirement_retarded_relative_input_s": float(latest["retirement"]),
        "latest_boundary_arrival_lab_s": witnesses["boundary"]["arrival_lab_s"],
        "latest_boundary_retarded_relative_input_s": float(latest["boundary"]),
        "latest_pixel_arrival_lab_s": witnesses["pixel"]["arrival_lab_s"],
        "latest_pixel_retarded_relative_input_s": float(latest["pixel"]),
        "latest_aperture_envelope_arrival_lab_s": float(aperture_target + detector_light_time),
        "latest_aperture_envelope_retarded_relative_input_s": float(aperture_target),
        "aperture_envelope": {
            "latest_observed_retirement_retarded_relative_input_s": float(latest["retirement"]),
            "full_physical_aperture_diagonal_m": math.sqrt(aperture_rho2),
            "independent_maximum_geometric_delay_s": aperture_delay,
            "independent_time_plus_geometry_bound_s": float(independent_aperture_target),
            "observed_per_particle_boundary_bound_s": float(latest["boundary"]),
            "selected_bound": ("independent_full_aperture" if independent_aperture_target >= latest["boundary"]
                               else "observed_per_particle_boundary"),
            "time_witness_particle_id": witnesses["retirement"]["particle_id"],
            "definition": "max(observed per-particle boundary arrival, latest observed retirement retarded time + longest path delay between opposite physical-aperture corners)",
            "limitation": "The independent spatial bound covers new in-aperture transverse positions at the observed latest retirement time; it does not prove that a new held trajectory has no later retirement time. The observed per-particle term retains observed out-of-aperture projections.",
        },
        "quality": {
            "complete": complete, "committed_records": committed,
            "all_records_scanned": True, "finite_required_fields": True,
            "all_downstream": True, "positive_longitudinal_velocity": True,
            "positive_mass_and_weight": True, "source_unchanged": True,
            "maximum_reference_plane_error_m": max_plane_error,
            "outside_physical_aperture_at_reference": outside_reference,
            "projected_outside_physical_aperture_at_retirement_exit": outside_retirement,
        },
        "aperture": {
            "physical_half_extents_m": half.tolist(),
            "outermost_physical_pixel_centres_m": pixel_half.tolist(),
            "pixel_spacing_m": spacing.tolist(),
            "physical_pixel_counts": np.rint(pixel_counts).astype(int).tolist(),
            "centre_xy_m": [0.0, 0.0],
        },
        "ranges": extrema, "witnesses": witnesses,
        "projection": {
            "retarded_time_definition": "t_lab - (z - input_plane_z_m)/c",
            "applied_extra_time_offset_s": 0.0,
            "common_time_note": "Input common-time and boost origins are already encoded in stored lab events; do not subtract them again.",
            "motion": "constant observed lab velocity from reference plane to retirement exit",
            "signal": "last direct vacuum arrival at the farthest point of the physical aperture",
            "latest_endpoint_reason": "Along any subluminal trajectory d(t+distance/c)/dt=1-n dot beta > 0, so the endpoint is latest for each detector point.",
            "cpml_policy": "All records are retained in the projection; transverse exits are counted. Actual CPML carrier disappearance is not modelled.",
            "limitations": [
                "Historical front-anchored trajectories do not rigorously bound a new held-start trajectory.",
                "A pre-retirement reference plane leaves an active segment before the taper; this projection does not evolve its self fields.",
                "Observed crossing momentum uses the solver's crossing interpolation, not a continuous exact orbit.",
                "This is a direct-arrival model, not proof of pure radiation, negligible reflections, or zero residual numerical field.",
            ],
        },
    }


def make_tiny_fixture(path, group_name="ballistic_reference", plane_z_m=4.03):
    """Write three small analytic records, for tests only; return their count."""
    dtype = np.dtype([("particle_id", "<u8"), ("source_id", "<u8"),
                      ("time_s", "<f8"), ("position_m", "<f8", (3,)),
                      ("proper_velocity", "<f8", (3,)), ("charge_C", "<f8"),
                      ("mass_kg", "<f8"), ("weight", "<f8"),
                      ("crossing_direction", "<i4")])
    rows = np.zeros(3, dtype=dtype)
    rows["particle_id"] = [1, 2, 3]
    rows["source_id"] = [11, 12, 13]
    rows["position_m"] = [[0, 0, plane_z_m], [.0003, -.0002, plane_z_m],
                          [-.0001, .00005, plane_z_m]]
    rows["proper_velocity"] = [[0, 0, 1200], [.02, -.01, 1150], [-.03, .02, 1180]]
    rows["time_s"] = (plane_z_m + .05) / C + np.array([0, 1e-15, 2e-15])
    rows["charge_C"] = -1.602176634e-19
    rows["mass_kg"] = 9.1093837139e-31
    rows["weight"] = 1
    rows["crossing_direction"] = 1
    with h5py.File(path, "w") as h:
        g = h.create_group(group_name)
        g.attrs.update(frame="laboratory", proper_velocity_unit="gamma*v/c", plane_z_m=plane_z_m)
        g.create_dataset("records", data=rows)
        g.create_dataset("committed_records", data=len(rows))
        g.create_dataset("complete", data=1)
    return len(rows)


def self_test():
    with tempfile.TemporaryDirectory(prefix="aprl-arrival-fixture-") as directory:
        file = Path(directory) / "reference.h5"
        make_tiny_fixture(file)
        kwargs = dict(input_plane_z_m=-.05, retirement_exit_z_m=4.33,
                      detector_z_m=5.673099411756358, physical_half_x_m=.00035,
                      physical_half_y_m=.00035, expected_records=3)
        result = analyze_reference(file, chunk_rows=1, **kwargs)
        with h5py.File(file) as h:
            rows = h["ballistic_reference/records"][:]
        # Independent direct calculation: enumerate every aperture corner,
        # compute physical flight times and distances, then maximize.
        for key, half in (("boundary", .00035), ("pixel", .0003475)):
            arrivals = []
            for row in rows:
                u = row["proper_velocity"].astype(np.longdouble)
                gamma = np.sqrt(1 + np.dot(u, u))
                beta = u / gamma
                travel = (np.longdouble(4.33) - np.longdouble(4.03)) / (C * beta[2])
                position = row["position_m"].astype(np.longdouble) + C * beta * travel
                for x in (-half, half):
                    for y in (-half, half):
                        delta = np.array([x, y, kwargs["detector_z_m"]], dtype=np.longdouble) - position
                        arrivals.append(np.longdouble(row["time_s"]) + travel + np.sqrt(np.dot(delta, delta)) / C)
            expected = max(arrivals)
            actual = result["latest_" + key + "_arrival_lab_s"]
            if abs(np.longdouble(actual) - expected) > 8e-24:
                raise AssertionError("Corner enumeration differs from streaming extrema")
            tau = expected - (np.longdouble(kwargs["detector_z_m"]) - np.longdouble(-.05)) / C
            if abs(np.longdouble(result["latest_" + key + "_retarded_relative_input_s"]) - tau) > 8e-24:
                raise AssertionError("Retarded time origin mismatch")
        # Independent spatial envelope: allow the latest source event to be
        # at any aperture corner, independently of its observed x/y.
        aperture_diagonal = np.longdouble(.0007) * np.sqrt(np.longdouble(2))
        axial = np.longdouble(kwargs["detector_z_m"]) - np.longdouble(4.33)
        spatial_delay = (np.sqrt(axial * axial + aperture_diagonal ** 2) - axial) / C
        expected_envelope = max(
            np.longdouble(result["latest_boundary_retarded_relative_input_s"]),
            np.longdouble(result["latest_retirement_retarded_relative_input_s"]) + spatial_delay)
        if abs(np.longdouble(result["latest_aperture_envelope_retarded_relative_input_s"]) - expected_envelope) > 8e-24:
            raise AssertionError("Independent full-aperture bound differs")
        if result["latest_aperture_envelope_retarded_relative_input_s"] < result["latest_boundary_retarded_relative_input_s"]:
            raise AssertionError("Envelope does not retain observed boundary arrival")
        shifted_by = 123e-15
        with h5py.File(file, "r+") as h:
            rows["time_s"] += shifted_by
            h["ballistic_reference/records"][:] = rows
        shifted = analyze_reference(file, chunk_rows=2, **kwargs)
        observed_shift = shifted["latest_boundary_retarded_relative_input_s"] - result["latest_boundary_retarded_relative_input_s"]
        if abs(observed_shift - shifted_by) > 8e-24:
            raise AssertionError("Lab clock shift was not preserved")
        with h5py.File(file, "r+") as h:
            rows["position_m"][0, 0] = .002
            h["ballistic_reference/records"][:] = rows
        outside = analyze_reference(file, **kwargs)
        if (outside["aperture_envelope"]["selected_bound"] != "observed_per_particle_boundary"
                or outside["quality"]["projected_outside_physical_aperture_at_retirement_exit"] == 0):
            raise AssertionError("Out-of-aperture source bound was discarded")
        with h5py.File(file, "r+") as h:
            rows["proper_velocity"][0, 2] = 0
            h["ballistic_reference/records"][:] = rows
        try:
            analyze_reference(file, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("Nonpositive velocity was accepted")
        return {"pass": True, "records": 3,
                "checks": ["streaming_vs_direct_corner_enumeration", "retarded_time_origin",
                           "independent_full_aperture_envelope", "lab_clock_shift",
                           "out_of_aperture_bound_retained", "nonpositive_velocity_rejected"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--output", "--out", type=Path)
    parser.add_argument("--reference-file", type=Path)
    parser.add_argument("--input-plane-z", type=float)
    parser.add_argument("--retirement-exit-z", type=float, default=4.33)
    parser.add_argument("--detector-z", type=float, default=5.673099411756358)
    parser.add_argument("--physical-half-x", type=float, default=.00035)
    parser.add_argument("--physical-half-y", type=float, default=.00035)
    parser.add_argument("--dx", type=float, default=5e-6)
    parser.add_argument("--dy", type=float, default=5e-6)
    parser.add_argument("--chunk-rows", type=int, default=262144)
    parser.add_argument("--expected-records", type=int)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    if not os.environ.get("SLURM_JOB_ID", "").isdigit():
        raise RuntimeError("Full reference scans require a real Slurm job; use --self-test for the tiny fixture")
    if args.output is None or (args.workspace is None and args.reference_file is None):
        parser.error("Need --output and --workspace or --reference-file")
    reference = args.reference_file or args.workspace.resolve() / REFERENCE_RELATIVE
    if args.input_plane_z is None:
        if args.workspace is None:
            parser.error("Explicit --reference-file requires --input-plane-z or --workspace")
        historical = args.workspace.resolve() / REFERENCE_RELATIVE
        plan = json.loads((historical.parent.parent / "input-plan.json").read_text())
        input_z = float(plan["input_plane_z_m"])
        expected_records = args.expected_records or int(plan["records"])
    else:
        input_z = args.input_plane_z
        expected_records = args.expected_records
    if args.output.exists():
        raise ValueError("Refusing to overwrite an observed-arrival report")
    result = analyze_reference(reference, input_plane_z_m=input_z,
        retirement_exit_z_m=args.retirement_exit_z, detector_z_m=args.detector_z,
        physical_half_x_m=args.physical_half_x, physical_half_y_m=args.physical_half_y,
        dx_m=args.dx, dy_m=args.dy, chunk_rows=args.chunk_rows,
        expected_records=expected_records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: result[key] for key in
        ("pass", "records", "latest_boundary_retarded_relative_input_s",
         "latest_pixel_retarded_relative_input_s",
         "latest_aperture_envelope_retarded_relative_input_s")}, indent=2))


if __name__ == "__main__":
    main()
