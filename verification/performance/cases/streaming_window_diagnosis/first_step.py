#!/usr/bin/env python3
"""Audit the first CPML error using recorded fields, without advancing a solver.

Independent of analyze.py and the boundary adapter. This proves a local first
step mechanism, not that it explains the complete detector error or supplies a
usable replacement boundary. The roundoff checks below are algebra checks,
not acceptance tolerances for the moving-window physics.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct

import numpy as np


NAMES = ("Ex", "Ey", "Ez", "Bx", "By", "Bz")
C = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
DZ = 1e-6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def exact(stream, size):
    block = stream.read(size)
    require(len(block) == size, "Truncated native section file")
    return block


def read_sections(directory, metrics):
    manifest = json.loads((directory / "field-sections.json").read_text())
    require(manifest["format"] == "aprl-native-field-xz-sections-v1" and
            manifest["completed"] and manifest["byte_order"] == "little",
            "Unsupported native section contract")
    steps = manifest["maximum_state_step"]
    require(2 <= steps <= metrics["steps"] and manifest["records"] == 1 + 3 * steps,
            "Invalid section coverage")
    nx, ny, nz = metrics["nx"], metrics["ny"], metrics["window_nz"]
    shapes = ((nx, ny+1, nz+1), (nx+1, ny, nz+1), (nx+1, ny+1, nz),
              (nx+1, ny, nz), (nx, ny+1, nz), (nx, ny, nz+1))
    expected = [(0, 0)] + [(n, s) for n in range(1, steps+1) for s in (1, 2, 3)]
    records = {}
    with (directory / "field-sections.bin").open("rb") as stream:
        require(exact(stream, 8) == b"APRLXZ01", "Incorrect section magic")
        for key in expected:
            n, stage, offset, count = struct.unpack("<4Q", exact(stream, 32))
            require((n, stage) == key and count == 6, "Incorrect record order")
            lattice_step = n-1 if stage in (1, 2) else n
            expected_offset = (0 if metrics["mode"] == "full" else
                               math.floor(218-.8*lattice_step)+84-24-20)
            require(offset == expected_offset, "Unexpected native lattice offset")
            values = []
            for component, shape in enumerate(shapes):
                a, cx, cy, cz, j = struct.unpack("<5Q", exact(stream, 40))
                require(a == component and (cx, cy, cz) == shape and
                        j == manifest["section_j"] == 14, "Native section shape mismatch")
                value = np.frombuffer(exact(stream, 8*cx*cz), "<f8").reshape(cz, cx)
                require(np.isfinite(value).all(), "Nonfinite native field")
                values.append(value)
            records[key] = (offset, values)
        require(not stream.read(1), "Trailing native section bytes")
    return records, steps


def max_abs(value):
    return float(np.max(np.abs(value)))


def rms(value):
    return float(np.sqrt(np.mean(np.square(value, dtype=np.longdouble))))


def roundoff_check(error, operands):
    scale = max(max_abs(value) for value in operands)
    limit = 1024*np.finfo(float).eps*max(scale, np.finfo(float).tiny)
    return {"max_abs_error": max_abs(error), "roundoff_limit": float(limit),
            "passed": bool(max_abs(error) <= limit)}


def cached_translation(initial, indices, alpha):
    # Reconstruct the frozen cache's actual coordinate/floor arithmetic.
    coordinate = np.asarray(indices, dtype=float)-alpha
    lower = np.floor(coordinate).astype(int)
    fraction = coordinate-lower
    require(lower.min() >= 0 and lower.max()+1 < initial.shape[0],
            "Reference reconstruction outside cache")
    return ((1-fraction[:, None])*initial[lower] +
            fraction[:, None]*initial[lower+1])


def bit_failures(left, right):
    require(left.shape == right.shape, "Bitwise comparison shape mismatch")
    return int(np.count_nonzero(left.view(np.uint64) != right.view(np.uint64)))


def audit(root, full_name, candidate_name):
    full_dir, candidate_dir = root/full_name, root/candidate_name
    fm = json.loads((full_dir/"metrics.json").read_text())
    cm = json.loads((candidate_dir/"metrics.json").read_text())
    require(fm["completed"] and cm["completed"] and fm["mode"] == "full" and
            cm["mode"] == "residual", "Expected completed original full/residual cases")
    for key in ("beta_box", "guard_cells", "dt_s", "steps", "nx", "ny", "seed",
                "particles", "cpml_cells", "history_policy", "full_nz"):
        require(fm[key] == cm[key], "Mismatched case parameter: "+key)
    require(cm["beta_box"] == .23 and cm["guard_cells"] == 24 and
            cm["nx"] == cm["ny"] == 24 and cm["window_nz"] == 104 and
            fm["window_nz"] == 424 and cm["cpml_cells"] == [0, 0, 20],
            "This algebra audit is restricted to the unchanged beta .23/g24 fixture")
    full, full_steps = read_sections(full_dir, fm)
    short, short_steps = read_sections(candidate_dir, cm)
    require(full_steps == short_steps, "Mismatched snapshot coverage")
    beta, dt = cm["beta_box"], cm["dt_s"]
    alpha, scale = beta*C*dt/DZ, beta*C*EPS0*MU0
    require(0 < alpha < 1, "First-step interpolation formula requires 0 < alpha < 1")
    f0, f1 = full[0, 0][1], full[1, 1][1]
    offset, local1 = short[1, 1]
    first_magnetic = {}
    for a, name in enumerate(NAMES):
        actual = local1[a]
        baseline = f1[a][offset:offset+actual.shape[0]]
        first_magnetic[name] = {"bitwise_failures": bit_failures(actual, baseline),
                                "max_abs_difference": max_abs(actual-baseline)}

    # Stay four or more vertices away from the full domain's CPML interfaces;
    # the smoothing/difference commutation proof applies in this uniform region.
    k = np.arange(24, 398)
    identity = {}
    algebra_checks = []
    for b, e, sign in ((3, 1, -1.), (4, 0, 1.)):
        columns = np.arange(2, f0[b].shape[1]-2)
        ix = np.ix_(k, columns)
        e_here, e_next = f0[e][ix], f0[e][np.ix_(k+1, columns)]
        e_previous = f0[e][np.ix_(k-1, columns)]
        initial_formula = sign*scale*((1-alpha)/2*e_here+(1+alpha)/2*e_next)
        # Rbar is the actual first Cowan B step, not a proposed boundary model.
        rbar = sign*scale*((1+alpha)/2*e_here+(1-alpha)/2*e_next)
        translated = cached_translation(f0[b], k, alpha)[:, columns]
        actual = f1[b][ix]
        predicted_defect = sign*scale*alpha*(1-alpha)/2*(e_previous-2*e_here+e_next)
        observed_defect = translated-actual
        checks = {
            "initial_half_step_formula": roundoff_check(f0[b][ix]-initial_formula,
                                                         [f0[b][ix], initial_formula]),
            "Rbar_equals_first_Cowan_step": roundoff_check(actual-rbar, [actual, rbar]),
            "frozen_translation_defect_identity": roundoff_check(
                observed_defect-predicted_defect, [actual, translated, predicted_defect])}
        algebra_checks.extend(value["passed"] for value in checks.values())
        identity[NAMES[b]] = {"source": NAMES[e], "sign": sign, "checks": checks,
            "sample_count": int(actual.size), "global_k_range_inclusive": [24, 397],
            "x_index_range_inclusive": [int(columns[0]), int(columns[-1])],
            "max_abs_defect_T": max_abs(observed_defect),
            "max_abs_c_defect_V_m": C*max_abs(observed_defect),
            "rms_c_defect_V_m": C*rms(observed_defect)}

    # Recompute only the first electric CPML correction (psi initially zero).
    # All B values here come from recorded full evolution, not another solver.
    # Compare native matching components; no field colocation is introduced.
    cpml_prediction = {}
    width, pml = cm["window_nz"], cm["cpml_z_cells"]
    indices = np.arange(offset-1, offset+width)
    for e, b, curl_sign in ((0, 4, -1.), (1, 3, 1.)):
        residual = f1[b][indices]-cached_translation(f0[b], indices, alpha)
        predicted, observed, total_targets = [], [], []
        for local_k in list(range(1, pml))+list(range(width-pml+1, width)):
            depth = max((pml-local_k)/pml, (local_k-(width-pml))/pml, 0.)
            graded = depth**3
            sigma = -(4*C*math.log(1e-8)/(2*pml*DZ))*graded
            kappa = 1+7*graded
            decay = math.exp(-sigma/kappa*dt)
            drive = (decay-1)/kappa  # alphaFraction=0, initial psi=0.
            correction = (1/kappa-1+drive)*(residual[local_k+1]-residual[local_k])
            prediction = curl_sign*dt/(MU0*EPS0*DZ)*correction
            baseline = full[1, 2][1][e][offset+local_k]
            actual = short[1, 2][1][e][local_k]
            valid = slice(1, -1) if e == 1 else slice(None)  # Existing transverse PEC.
            predicted.extend(prediction[valid])
            observed.extend((actual-baseline)[valid])
            total_targets.extend(actual[valid])
        predicted, observed, total_targets = map(np.asarray, (predicted, observed, total_targets))
        check = roundoff_check(predicted-observed, [total_targets, predicted, observed])
        algebra_checks.append(check["passed"])
        cpml_prediction[NAMES[e]] = {"sample_count": int(predicted.size),
            "predicted_max_abs_V_m": max_abs(predicted),
            "observed_max_abs_V_m": max_abs(observed), "check": check}

    history = {}
    with (candidate_dir/"history.csv").open() as stream:
        for row in csv.DictReader(stream):
            history[int(row["step"])] = row
    promotions, total_bit_failures, retained_values = [], 0, 0
    for n in range(1, short_steps+1):
        old_offset, before = short[n, 2]
        new_offset, after = short[n, 3]
        if new_offset == old_offset:
            require(all(bit_failures(x, y) == 0 for x, y in zip(before, after)),
                    "No-move stage changed recorded fields")
            continue
        require(new_offset == old_offset-1, "Unexpected move direction/size")
        failures = [bit_failures(x[:-1], y[1:]) for x, y in zip(before, after)]
        values = sum(x[:-1].size for x in before)
        total_bit_failures += sum(failures)
        retained_values += values
        ordinary_entry_errors = {}
        for a, name in enumerate(NAMES):
            baseline = full[n, 2][1][a][old_offset+pml]
            ordinary_entry_errors[name] = max_abs(before[a][pml]-baseline)
        row = history[n]
        promotions.append({"state_step": n, "old_offset": old_offset, "new_offset": new_offset,
            "old_native_k": pml, "new_native_k": pml+1,
            "old_E_vertex_world_z_m": fm["full_origin_z_m"]+(old_offset+pml)*DZ,
            "retained_native_values": int(values), "retained_bitwise_failures": sum(failures),
            "promoted_six_native_rows_bitwise_failures": sum(
                bit_failures(x[pml], y[pml+1]) for x, y in zip(before, after)),
            "promoted_rows_full_difference_by_component": ordinary_entry_errors,
            "ordinary_gauss_before_shift": float(row["total_gauss_before_shift"]),
            "ordinary_gauss_after_shift": float(row["total_gauss_after_shift"])})

    events = {}
    for name in ("release-particles.csv", "detector-crossings.csv"):
        with (candidate_dir/name).open() as stream:
            rows = list(csv.DictReader(stream))
        events[name] = min(float(row["step"])+float(row["fraction"]) for row in rows)

    files = [Path(__file__).resolve(), root/"source-files.json"]
    for directory in (full_dir, candidate_dir):
        files.extend(directory/name for name in ("metrics.json", "field-sections.json",
                                                "field-sections.bin", "history.csv"))
    files.extend(candidate_dir/name for name in ("release-particles.csv", "detector-crossings.csv"))
    source_hashes = json.loads((root/"source-files.json").read_text())
    return {"format": "aprl-moving-window-first-step-proof-v1", "completed": True,
        "scope": "Recorded j=14 native x-z sections, first-step algebra and retained field mapping; no solver was run",
        "operator_evidence_checks_passed": bool(all(algebra_checks) and total_bit_failures == 0),
        "physical_acceptance": "NOT_ASSESSED; this local mechanism is not an attribution of the entire 20.67% B discrepancy",
        "root": str(root), "full_case": full_name, "candidate_case": candidate_name,
        "input_sha256": {str(path): sha256(path) for path in files},
        "recorded_core_source_sha256": {path: digest for path, digest in source_hashes.items()
            if Path(path).name in ("eb_field.cpp", "eb_initial_field.cpp", "eb_cpml.cpp")},
        "explicit_coefficients": {"c_m_s": C, "epsilon0": EPS0, "mu0": MU0,
            "dz_m": DZ, "dt_s": dt, "beta": beta, "alpha_beta_c_dt_over_dz": alpha,
            "magnetic_scale_beta_c_epsilon0_mu0": scale,
            "initial_B_weights_Ek_Ekplus1": [(1-alpha)/2, (1+alpha)/2],
            "Rbar_first_B_weights_Ek_Ekplus1": [(1+alpha)/2, (1-alpha)/2],
            "CPML": {"z_cells": pml, "order": 3, "target_reflection": 1e-8,
                     "kappa_max": 8, "alpha_fraction": 0, "initial_psi": 0}},
        "algebra": [
            "s=beta*c*epsilon0*mu0; alpha=beta*c*dt/dz; sign=-1 for Bx/Ey and +1 for By/Ex",
            "B(-1/2,k)=sign*s*((1-alpha)*E(0,k)+(1+alpha)*E(0,k+1))/2",
            "Uniform-interior Cowan commutation gives Rbar_B(+1/2,k)=sign*s*((1+alpha)*E(0,k)+(1-alpha)*E(0,k+1))/2",
            "Rbar is the recorded original first magnetic step, not an imposed propagation approximation",
            "Frozen cache R_B(1,k)=(1-alpha)*B(-1/2,k)+alpha*B(-1/2,k-1), modulo coordinate floating rounding",
            "R_B(1,k)-Rbar_B(+1/2,k)=sign*s*alpha*(1-alpha)/2*(E(0,k-1)-2*E(0,k)+E(0,k+1))",
            "First electric CPML correction: curl_sign*dt/(mu0*epsilon0*dz)*(1/kappa-1+drive)*(deltaB[k]-deltaB[k-1]); deltaB=recorded_B(+1/2)-R_B(1), psi starts zero"],
        "first_magnetic_candidate_vs_full": first_magnetic,
        "reference_first_step_defect": identity, "first_electric_CPML_prediction": cpml_prediction,
        "remap": {"retained_native_values_checked": int(retained_values),
                  "retained_bitwise_failures": total_bit_failures, "promotions": promotions},
        "first_event_fractional_steps": events,
        "limitations": [
            "Native components retain different Yee positions and times; cB is only a unit conversion, not a lab transform or power",
            "The first CPML prediction excludes the explicitly clamped exterior E faces and transverse PEC edges",
            "Only j=14 is recorded; local algebra is not a complete three-dimensional energy balance",
            "Post-remap ordinary Gauss changes its sampled world region; a jump is not proof that copying corrupts stored fields",
            "At electric vertex old k=20 to new k=21, a row previously at the CPML interface enters the ordinary Gauss mask; magnetic component staggering is retained",
            "Incoming new-layer and halo/face discrepancies are separate possible later contributors; factorial replay is needed to separate their complete detector impact",
            "Replacing cached B by Rbar at one time would close this first-step identity but does not establish a valid long-time reference or boundary method"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--full", default="b023-full-g24")
    parser.add_argument("--candidate", default="b023-residual-g24")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root/"diagnosis"/"first_step.json"
    require(not output.exists(), "Refusing to overwrite "+str(output))
    result = audit(root, args.full, args.candidate)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(str(output))
    raise SystemExit(0 if result["operator_evidence_checks_passed"] else 1)


if __name__ == "__main__":
    main()
