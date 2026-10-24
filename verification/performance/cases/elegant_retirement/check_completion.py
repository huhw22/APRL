#!/usr/bin/env python3
"""Check normal retirement completion and observed-plane projected arrival coverage."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tempfile

import h5py
import numpy as np
import yaml

C = 299792458.0
E = 1.602176634e-19


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def check(root):
    from observed_arrival import analyze_reference
    report = dict(format="aprl-retirement-completion-v2", pass_=False, complete=False,
        output_transactions_complete=False, scientific_validation=False, pure_radiation_validated=False,
        limitations=[
            "Coverage projects observed pre-retirement particles at constant velocity to the taper exit, then light to the furthest physical detector boundary.",
            "The short drift before retirement entrance still includes self-field evolution; the exact post-record orbit is not recorded here.",
            "Acceptance does not establish pure radiation, absence of retirement transients, or physical convergence."])
    try:
        card = yaml.safe_load((root / "run.yaml").read_text())
        plan = json.loads((root / "run-plan.json").read_text())
        geometry = json.loads((root / "geometry-preflight.json").read_text())
        manifest = yaml.safe_load((root / "run-manifest.yaml").read_text())
        require(card["units"] == {"length": "m", "time": "s"}, "Checker requires the production SI card")
        core = plan.get("compiled_core_revision", plan.get("core_revision"))
        require(core and manifest["source_revision"] == core, "Manifest compute revision differs from frozen plan")
        require(digest(plan["binary"]) == plan["binary_sha256"], "Executable hash differs from frozen plan")
        require(manifest["mpi_size"] == plan["ranks"], "Manifest MPI layout differs from frozen plan")
        require(yaml.safe_load(manifest["configuration_yaml"]) == card, "Manifest configuration differs from frozen run card")
        cells = [*card["mesh"]["cells"][:2], geometry["estimated_cells_z"]]
        mesh = manifest["resolved_mesh"]
        require(mesh["cells"] == cells and mesh["cell_size_m"] == card["mesh"]["cell_size"], "Resolved grid differs from frozen geometry/card")
        require(mesh["center_m"] == card["mesh"]["center"] and not mesh["ballistic_outside_transverse"], "Unexpected mesh center or transverse ballistic mode")
        held = manifest["resolved_beam_reference"]
        require(held["startup_mode"] == "held" and held["frozen_lab_plane_injection"] and held["hold_transverse_until_activation"], "Manifest does not preserve held startup")
        for key in ("input_plane_z", "activation_plane_z"):
            require(held[key + "_m"] == card["beam"]["reference"][key], "Held lab plane differs from card")
        report["identity"] = dict(binary=plan["binary"], binary_sha256=plan["binary_sha256"], core_revision=core, mpi_size=manifest["mpi_size"], cells=cells)
        log = (root / "aprl.log").read_text(errors="replace")
        require("Configured stop reached: boost reference centre" in log, "No configured reference stop in APRL log")
        require("Simulation failure" not in log and "mesh.duration was exhausted" not in log, "APRL ended with a simulation failure")
        final = re.findall(r"\[resource\] phase=final step=(\d+)\b", log)
        require(len(final) == 1, "Missing or ambiguous final complete-step record")
        final_step = int(final[0]); expected = int(geometry["records"])
        activation = re.findall(r"Frozen lab-plane injection: activations=(\d+).*?still frozen=(\d+)", log)
        retirement = re.findall(r"Experimental particle retirement: entries=(\d+).*?completed=(\d+), residual Q at taper exit=([^ ]+)", log)
        require(len(activation) == len(retirement) == 1, "Missing or ambiguous held/retirement report")
        activated, frozen = map(int, activation[0]); entries, completed = map(int, retirement[0][:2])
        residual = float(retirement[0][2]); bound = 64 * np.finfo(float).eps * float(geometry["input_weight_sum"]) * E
        report["particles"] = dict(expected=expected, activated=activated, still_frozen=frozen,
            retirement_entries=entries, completed_carriers=completed, residual_exit_charge_C=residual if np.isfinite(residual) else None,
            exit_charge_roundoff_bound_C=bound, boundary_reports=[s for s in log.splitlines() if s.startswith(("CPML physical-domain exits:", "Direct outer-domain exits:", "Residual CPML-carrier cleanup:"))])
        require((activated, frozen) == (expected, 0), "Held injection did not activate the complete cohort")
        require((entries, completed) == (expected, expected), "Full-cohort retirement requirement not met; inspect physical/CPML exit counts. This does not by itself mean saved data are corrupt")
        require(np.isfinite(residual) and abs(residual) <= bound, "Retirement residual exceeds arithmetic-roundoff bound")
        fields = card["detectors"]["field_planes"]
        require(len(fields) == 1 and not fields[0].get("particle_background", {}).get("enabled", False), "Expected one retirement field plane without analytic background")
        products = [(root / "detectors" / (p["name"] + ".h5"), group)
            for key, group in (("field_planes", "field_plane"), ("particle_planes", "particle_plane")) for p in card["detectors"][key]]
        products.append((root / "ledger" / card["energy_ledger"]["filename"], "energy_ledger"))
        report["products"] = []; last_field_time = None; dt_box = card["mesh"]["cell_size"][2] / C
        for path, name in products:
            with h5py.File(path, "r") as file:
                g = file[name]; field = name == "field_plane"
                for key in ("run_id", "source_revision", "configuration_digest_fnv1a64"):
                    if key in manifest:
                        value = g.attrs.get(key, "")
                        value = value.decode() if isinstance(value, bytes) else str(value)
                        require(value == str(manifest[key]), "Diagnostic identity differs: " + str(path) + ":" + key)
                n = int(g["committed_samples" if field else "committed_records"][()]); data = g["fields" if field else "records"]
                require(int(g["complete"][()]) == 1 and n > 0 and len(data) == n, "Output incomplete: " + str(path))
                item = dict(file=str(path), group=name, records=n, complete=True)
                if field:
                    t = g["time_s"][:]
                    require(t.shape == (n,) and data.shape == (n, *cells[:2]), "Field time/grid/data shape mismatch")
                    require(np.isfinite(t).all() and np.all(np.diff(t) > 0), "Invalid field lab time axis")
                    cadence = dt_box / card["mesh"]["boost_gamma"]
                    # Subtracting ~19 ns binary64 stamps to recover ~4.5 as
                    # increments is limited by timestamp ULPs, not dt's ULP.
                    cadence_roundoff = 4 * abs(np.spacing(float(np.max(np.abs(t))))) + 64 * np.finfo(float).eps * abs(cadence)
                    require(n > 1 and np.all(np.abs(np.diff(t)-cadence) <= cadence_roundoff), "Field sampling cadence exceeds timestamp-roundoff bound")
                    tail = data[-1]
                    require(all(np.isfinite(tail[k]).all() for k in ("electric_V_per_m", "magnetic_T")), "Nonfinite last field sample")
                    last_field_time = float(t[-1]); item.update(first_lab_time_s=float(t[0]), last_lab_time_s=last_field_time,
                        field_finite_check_scope="last sample", cadence_expected_s=cadence, cadence_timestamp_roundoff_bound_s=cadence_roundoff)
                if name == "energy_ledger":
                    a = data[:]
                    require(all(np.isfinite(a[k]).all() for k in a.dtype.names if a[k].dtype.kind in "fc"), "Nonfinite energy ledger")
                    require(int(a[-1]["step"]) == final_step and np.all(np.diff(a["step"].astype(np.int64)) > 0), "Ledger final step or ordering differs from completed run")
                    require(np.isclose(float(a[-1]["time_box_s"]), final_step * dt_box, rtol=1e-12, atol=0), "Ledger final time differs from completed step")
                    active, removed = int(a[-1]["active_macroparticles"]), int(a[-1]["removed_macroparticles"])
                    report["particles"].update(ledger_final_active=active, ledger_final_removed=removed, ledger_final_step=final_step)
                    require(active == 0 and removed == expected, "Ledger does not show complete particle retirement")
                report["products"].append(item)
        report["output_transactions_complete"] = True
        pre = next(p for p in card["detectors"]["particle_planes"] if p["name"] == "pre-retirement-particles")
        entrance = card["particle_retirement"]["entrance_z"]
        require(pre["z"] < entrance and card["mesh"]["center"][:2] == [0, 0], "Expected upstream reference plane and centered transverse aperture")
        dx, dy = card["mesh"]["cell_size"][:2]; cpml = card["boundary"]["cells"]
        arrival = analyze_reference(root / "detectors" / (pre["name"] + ".h5"),
            input_plane_z_m=card["beam"]["reference"]["input_plane_z"], retirement_exit_z_m=entrance + card["particle_retirement"]["length"],
            detector_z_m=fields[0]["z"], physical_half_x_m=(cells[0] / 2 - cpml[0]) * dx,
            physical_half_y_m=(cells[1] / 2 - cpml[1]) * dy, dx_m=dx, dy_m=dy, expected_records=expected)
        latest = float(arrival["latest_boundary_arrival_lab_s"]); require(np.isfinite(latest), "Nonfinite projected latest arrival")
        covered = last_field_time >= latest
        report["observed_reference_arrival"] = arrival
        report["arrival_coverage"] = dict(pass_=bool(covered), last_field_lab_time_s=last_field_time,
            latest_projected_boundary_arrival_lab_s=latest, remaining_window_fs=(last_field_time-latest)*1e15,
            comparison_tolerance_s=0, reference_to_retirement_entrance_m=entrance-pre["z"],
            criterion="Last saved field reaches the latest projected retirement-exit photon at the furthest physical detector boundary; no arbitrary quiet tail.")
        require(covered, "Saved field window ends before the observed-reference projected latest arrival")
        report.update(pass_=True, complete=True)
    except (ValueError, KeyError, OSError, StopIteration, TypeError) as error:
        report["error"] = str(error) or type(error).__name__
    save(root / "completion-check.json", report)
    return report


def self_test():
    from observed_arrival import analyze_reference, make_tiny_fixture
    with tempfile.TemporaryDirectory(prefix="aprl-completion-test-") as temporary:
        root = Path(temporary); (root / "detectors").mkdir(); (root / "ledger").mkdir()
        particle = root / "detectors/pre-retirement-particles.h5"; plane, entrance, exit_z, detector = 4.027, 4.03, 4.33, 5.673
        n = make_tiny_fixture(particle, group_name="particle_plane", plane_z_m=plane)
        observed = analyze_reference(particle, input_plane_z_m=-.05, retirement_exit_z_m=exit_z,
            detector_z_m=detector, physical_half_x_m=.00035, physical_half_y_m=.00035, dx_m=5e-6, dy_m=5e-6, expected_records=n)
        binary = root / "binary"; binary.write_bytes(b"synthetic identity only; never executed")
        card = dict(units=dict(length="m", time="s"), mesh=dict(cells=[160,160,"auto"], cell_size=[5e-6,5e-6,1.2543262415359034e-6], center=[0,0,0], boost_gamma=929.1305492858543),
            beam=dict(reference=dict(input_plane_z=-.05, activation_plane_z=-.05)), boundary=dict(cells=[10,10,20]),
            particle_retirement=dict(entrance_z=entrance, length=exit_z-entrance),
            detectors=dict(field_planes=[dict(name="field", z=detector)], particle_planes=[dict(name="pre-retirement-particles", z=plane)]),
            energy_ledger=dict(filename="energy-ledger.h5"))
        plan = dict(binary=str(binary), binary_sha256=digest(binary), compiled_core_revision="synthetic-core", ranks=256)
        manifest = dict(source_revision="synthetic-core", mpi_size=256, configuration_yaml=yaml.safe_dump(card),
            resolved_mesh=dict(cells=[160,160,4096], cell_size_m=card["mesh"]["cell_size"], center_m=[0,0,0], ballistic_outside_transverse=False),
            resolved_beam_reference=dict(startup_mode="held", frozen_lab_plane_injection=True, hold_transverse_until_activation=True, input_plane_z_m=-.05, activation_plane_z_m=-.05))
        (root / "run.yaml").write_text(yaml.safe_dump(card)); (root / "run-manifest.yaml").write_text(yaml.safe_dump(manifest))
        save(root / "run-plan.json", plan); save(root / "geometry-preflight.json", dict(records=n, input_weight_sum=n, estimated_cells_z=4096))
        log = (f"Frozen lab-plane injection: activations={n}, still frozen=0.\n"
               f"Experimental particle retirement: entries={n} (Q=0 C), completed={n}, residual Q at taper exit=0 C\n"
               "[resource] phase=final step=2\nConfigured stop reached: boost reference centre reached target\n")
        (root / "aprl.log").write_text(log); dt = card["mesh"]["cell_size"][2]/C
        cadence = dt/card["mesh"]["boost_gamma"]
        times = observed["latest_boundary_arrival_lab_s"] + 1e-14 - np.arange(2,-1,-1)*cadence
        with h5py.File(root / "detectors/field.h5", "w") as f:
            g=f.create_group("field_plane"); g["complete"]=1; g["committed_samples"]=3; g["time_s"]=times
            g.create_dataset("fields", shape=(3,160,160), chunks=(1,160,160), dtype=[("electric_V_per_m","f8",(3,)),("magnetic_T","f8",(3,))])
        with h5py.File(root / "ledger/energy-ledger.h5", "w") as f:
            g=f.create_group("energy_ledger"); g["complete"]=1; g["committed_records"]=3
            a=np.zeros(3,dtype=[("step","u8"),("time_box_s","f8"),("active_macroparticles","u8"),("removed_macroparticles","u8")])
            a["step"]=[0,1,2]; a["time_box_s"]=np.arange(3)*dt; a["removed_macroparticles"]=n; g["records"]=a
        for path, group in ((particle,"particle_plane"),(root/"detectors/field.h5","field_plane"),(root/"ledger/energy-ledger.h5","energy_ledger")):
            with h5py.File(path,"r+") as f: f[group].attrs["source_revision"]="synthetic-core"
        accepted=check(root); require(accepted["pass_"], "Synthetic success rejected: "+accepted.get("error",""))
        require(1e-8 < abs(times[-1]) < 3e-8 and 4e-18 < cadence < 5e-18, "Fixture no longer exercises high-gamma timestamp cancellation")
        with h5py.File(root / "detectors/field.h5", "r+") as f: f["field_plane/time_s"][1] = times[1]+cadence*.1
        require(not check(root)["pass_"], "Real cadence change accepted")
        with h5py.File(root / "detectors/field.h5", "r+") as f: f["field_plane/time_s"][:]=times
        with h5py.File(root / "detectors/field.h5", "r+") as f: f["field_plane/time_s"][:]=times-2e-14
        failed=check(root)
        require(not failed["pass_"] and failed["output_transactions_complete"] and not failed["arrival_coverage"]["pass_"], "Insufficient arrival window accepted")
        with h5py.File(root / "detectors/field.h5", "r+") as f: f["field_plane/time_s"][:]=times
        (root / "aprl.log").write_text(log.replace(f"completed={n}",f"completed={n-1}"))
        require(not check(root)["pass_"], "Incomplete retirement accepted"); (root / "aprl.log").write_text(log)
        manifest["source_revision"]="wrong-core"; (root / "run-manifest.yaml").write_text(yaml.safe_dump(manifest))
        require(not check(root)["pass_"], "Wrong compute revision accepted")
    return dict(pass_=True, tests=["high-gamma complete observed-reference coverage within timestamp ULPs", "real cadence change rejected", "insufficient window rejected", "incomplete retirement rejected", "wrong compute revision rejected"])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("root",type=Path,nargs="?"); parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if not args.self_test and args.root is None: parser.error("root is required")
    result=self_test() if args.self_test else check(args.root.resolve())
    print(json.dumps(result,indent=2,allow_nan=False))
    raise SystemExit(0 if result["pass_"] else 1)
