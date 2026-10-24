#!/usr/bin/env python3
"""Freeze one production ELEGANT retirement run from a checked historical card."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import yaml

BINARY_SHA = "56cc1e0931103687dcc440c3ad390a923b2e74f5da0d68ef1806e3992ccbfeac"
CORE_REVISION = "ac58e67e33d6b4ecfd12c3d00f4a9947d0e8500d"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace", type=Path, required=True)
    p.add_argument("--geometry", type=Path, required=True)
    p.add_argument("--preflight-result", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    w = args.workspace.resolve()
    root = args.output_dir.resolve()
    if root.exists():
        raise ValueError("Production output directory must be new")
    g = json.loads(args.geometry.read_text())
    if not g["pass"] or g["initial_outside_physical_count"] != 0:
        raise ValueError("Full-input geometry preflight must pass with all initial particles inside")
    policy = g.get("arrival_policy", {})
    if not policy.get("observed_arrival_report") or policy.get("quiet_tail_fs") != 0:
        raise ValueError("Compact production requires observed arrival coverage with no arbitrary quiet tail")
    if not (args.preflight_result / "PASS").is_file():
        raise ValueError("Completed MPI/retirement functional preflight is required")
    binary = w / "build/candidate/aprl"
    if sha(binary) != BINARY_SHA:
        raise ValueError("Accepted binary hash differs")
    original = w / "results/elegant_4m_analysis_20261024/original"
    historical = original / "runs/20261021T070618Z-861230/run.yaml"
    base = yaml.safe_load(historical.read_text())
    if base["mesh"]["cells"] != [160, 160, "auto"]:
        raise ValueError("Unexpected historical transverse geometry")
    particle = original / "particles.h5"
    stat = particle.stat()
    if (str(particle), stat.st_size, stat.st_mtime_ns) != (
            g["input_file"], g["input_size_bytes"], g["input_mtime_ns"]):
        raise ValueError("Input metadata differs from full-input geometry preflight")
    card = copy.deepcopy(base)
    card["mesh"]["ballistic_outside_transverse"] = False
    card["mesh"]["longitudinal"]["lower_guard_cells"] = g["recommended_lower_guard_cells"]
    card["mesh"]["longitudinal"]["upper_guard_cells"] = g["recommended_upper_guard_cells"]
    card["initial_self_field"]["include_transverse_cpml_particles"] = False
    reference = card["beam"]["reference"]
    reference.update(frozen_lab_plane_injection=True, hold_transverse_until_activation=True,
                     activation_plane_z=reference["input_plane_z"])
    card["beam"]["input"]["file"] = str(particle)
    card["particle_retirement"] = {
        "enabled": True, "entrance_z": g["retirement_entrance_z_m"],
        "length": g["retirement_length_m"]}
    field = card["detectors"]["field_planes"][0]
    field["name"] = "post-retirement-field"
    field["z"] = g["detector_z_m"]
    field["particle_background"] = {"enabled": False}
    field["retirement_frequency_protection"] = {
        "enabled": True, "minimum_photon_energy_eV": 50.0, "cycles": 2}
    card["detectors"]["particle_planes"].append({
        "name": "pre-retirement-particles",
        "z": g["retirement_entrance_z_m"] - 2 * base["mesh"]["boost_gamma"] * base["mesh"]["cell_size"][2],
        "buffer_records": 65536, "compression": 0})
    card["stop"] = {"mode": "reference-center-z", "z": g["stop_z_m"]}
    if g["expected_stop_box_s"] >= card["mesh"]["duration"]:
        raise ValueError("Physical stop exceeds existing duration guard")
    card["output"]["manifest"] = str(root / "run-manifest.yaml")
    card["detectors"]["directory"] = str(root / "detectors")
    card["energy_ledger"]["directory"] = str(root / "ledger")
    root.mkdir(parents=True)
    (root / "scripts").mkdir()
    here = Path(__file__).resolve().parent
    for name in ("prepare_case.py", "run.sbatch", "check_completion.py", "geometry_preflight.py",
                 "compact_geometry.py", "observed_arrival.py"):
        shutil.copy2(here / name, root / "scripts" / name)
    shutil.copy2(args.geometry, root / "geometry-preflight.json")
    shutil.copy2(policy["observed_arrival_report"], root / "arrival-projection.json")
    shutil.copy2(historical, root / "historical-run.yaml")
    (root / "run.yaml").write_text(yaml.safe_dump(card, sort_keys=False))
    copy_manifest = json.loads((original.parent / "copy-manifest.json").read_text())
    input_sha = next(f["sha256"] for f in copy_manifest["files"] if f["path"] == "particles.h5")
    plan = dict(
        binary=str(binary), binary_sha256=BINARY_SHA, compiled_core_revision=CORE_REVISION,
        preparation_revision=subprocess.check_output(
            ["git", "-C", str(w / "source/aprl"), "rev-parse", "HEAD"], text=True).strip(),
        input_file=str(particle), input_size_bytes=stat.st_size,
        input_mtime_ns=stat.st_mtime_ns, input_sha256=input_sha,
        preflight_result=str(args.preflight_result.resolve()),
        geometry_preflight=str(args.geometry.resolve()),
        nodes=2, ranks=256, threads_per_rank=1, exclusive=True, time_limit="144:00:00",
        signal_only=True, zero_field_control=False, transient_validation=False,
        user_interpretation_assumption="Treat the downstream field after retirement and propagation as radiation; residual transients are not tested in this run.",
        changes=[
            "Enable the existing C2 retirement layer after the undulator.",
            "Disable analytic particle-background recording.",
            "Move the fixed lab field plane downstream of retirement and its causal clearance.",
            "Use observed exit-event arrivals plus retirement and finite-aperture propagation to set the physical stop; no arbitrary quiet tail or dynamics allowance.",
            "Restore the original 64/128 guard floor and left anchoring; align only the total mesh to 256 cells.",
            "Use existing held frozen-plane startup; disable incompatible transverse ballistic injection and CPML initial sources.",
            "Replace the old ballistic-reference particle record with an ordinary pre-retirement plane two lab-equivalent cells before the retirement entrance.",
            "Use the accepted optimized executable; keep the historical beam, charge, transverse mesh, dz, boost, magnetic device, subcycling and sampling cadence."])
    (root / "run-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    files = [root / "run.yaml", root / "run-plan.json", root / "geometry-preflight.json", root / "arrival-projection.json",
             *sorted((root / "scripts").iterdir()), binary, w / "build/cluster-tools/binding_probe"]
    (root / "frozen.sha256").write_text("".join(f"{sha(f)}  {f}\n" for f in files))
    print(root)


if __name__ == "__main__":
    main()
