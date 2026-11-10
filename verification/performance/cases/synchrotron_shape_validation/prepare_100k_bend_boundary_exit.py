"""Prepare bend-only controls whose particles leave through the right boundary."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import yaml


NX = 80
NY = 128
DX = 5e-6
CENTER_X = -20e-6
CPML_CELLS = 10
APERTURE_RADIUS = 75e-6
DETECTOR_Z = .3524
REFERENCE_NPZ = Path(
    "/public/home/wfollowm/APRL_Workspace/results/"
    "synchrotron_sampling_grid_controls_20261110_r1/"
    "transverse_single_frequency_8eV_separate_r1/transverse_maps_8eV.npz")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def rewrite(value, old, new):
    if isinstance(value, dict):
        return {key: rewrite(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, old, new) for item in value]
    if isinstance(value, str):
        return value.replace(str(old), str(new))
    return value


def main(output, specifications):
    output = Path(output).resolve()
    require(not output.exists(), "Output already exists")
    require(REFERENCE_NPZ.is_file(), "Reference transverse map is missing")
    roots = []
    for specification in specifications:
        label, text = specification.split("=", 1)
        root = Path(text).resolve()
        require((root / "PIC_COMPLETE.json").is_file(),
                f"Incomplete source campaign: {root}")
        roots.append((label, root))
    require(len(roots) == 2, "Exactly two independent 100k roots are required")

    first_campaign = json.loads((roots[0][1] / "campaign.json").read_text())
    binary = Path(first_campaign["binary"])
    require(sha(binary) == first_campaign["binary_sha256"], "Binary changed")
    require(all(json.loads((root / "campaign.json").read_text())["binary_sha256"]
                == first_campaign["binary_sha256"] for _, root in roots),
            "Source campaigns use different binaries")

    output.mkdir(parents=True)
    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ("prepare_100k_bend_boundary_exit.py",
                 "run_100k_bend_boundary_exit.py",
                 "analyze_100k_bend_boundary_exit.py",
                 "bend_boundary_exit.sbatch"):
        shutil.copy2(here / name, scripts / name)
    shutil.copy2(roots[0][1] / "scripts/launcher.py", scripts / "launcher.py")
    shutil.copy2(roots[0][1] / "scripts/prepare_paired.py",
                 scripts / "prepare.py")

    full_left = CENTER_X - NX * DX / 2
    full_right = CENTER_X + NX * DX / 2
    physical_right = full_right - CPML_CELLS * DX
    theta = 0.0024 / 3.0
    bend_exit_x = 3.0 * (1.0 - math.cos(theta))
    beam_center = bend_exit_x + (DETECTOR_Z - .0524) * math.sin(theta)
    boundary_crossing_z = .0524 + (full_right - bend_exit_x) / math.sin(theta)
    sigma_x = 20e-6
    expected_beyond = .5 * math.erfc(
        (full_right - beam_center) / (math.sqrt(2) * sigma_x))

    cases = []
    for label, root in roots:
        source = root / "bend"
        target = output / label / "bend"
        target.mkdir(parents=True)
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source / "run.yaml").read_text())), source, target)
        config["mesh"]["cells"][0] = NX
        config["mesh"]["cells"][1] = NY
        config["mesh"]["center"][0] = CENTER_X
        config["output"]["manifest"] = str(target / "run-manifest.yaml")
        config["mesh"]["moving_window"]["reference_directory"] = str(
            target / "reference")
        config["detectors"]["directory"] = str(target / "detectors")
        config["energy_ledger"]["directory"] = str(target / "ledger")
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

        plan = rewrite(json.loads((source / "PLAN.json").read_text()),
                       source, target)
        plan["geometry"].update(nx=NX, ny=NY)
        plan.update(
            case_directory=str(target), label=label,
            purpose=("Test whether naturally absorbing the bent particles at the "
                     "right transverse boundary removes their velocity near field "
                     "from the central radiation aperture"),
            sole_geometry_change=(
                "nx 160 -> 80 and center_x 120 um -> -20 um; all longitudinal "
                "geometry, particles, push, deposition, solver, CPML thickness and "
                "detector timing are unchanged"),
            transverse_domain_m=[full_left, full_right],
            right_cpml_start_m=physical_right,
            aperture_radius_m=APERTURE_RADIUS,
            predicted_beam_center_at_detector_m=beam_center,
            predicted_reference_boundary_crossing_z_m=boundary_crossing_z,
            predicted_fraction_beyond_right_boundary_at_detector=expected_beyond,
            physical_acceptance=False)
        save(target / "PLAN.json", plan)
        cases.append({"label": label, "directory": str(target),
                      "timeout_seconds": 1800,
                      "expected_particles": int(plan["expected_particles"])})

    campaign = {
        "binary": str(binary),
        "binary_sha256": sha(binary),
        "binding_probe": first_campaign["binding_probe"],
        "cases": cases,
        "source_roots": {label: str(root) for label, root in roots},
        "reference_npz": str(REFERENCE_NPZ),
        "reference_npz_sha256": sha(REFERENCE_NPZ),
        "aperture_radius_m": APERTURE_RADIUS,
        "transverse_domain_m": [full_left, full_right],
        "right_cpml_start_m": physical_right,
        "predicted_beam_center_at_detector_m": beam_center,
        "predicted_reference_boundary_crossing_z_m": boundary_crossing_z,
        "predicted_propagation_after_boundary_exit_m":
            DETECTOR_Z - boundary_crossing_z,
        "scope": (
            "Bend-only boundary-exit trial. It is not presumed to be cleaner: "
            "transition radiation and CPML/particle-removal transients are retained."),
    }
    save(output / "campaign.json", campaign)
    frozen_paths = [output / "campaign.json", binary, REFERENCE_NPZ]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
        input_path = Path(yaml.safe_load((directory / "run.yaml").read_text())
                          ["beam"]["input"]["file"])
        frozen_paths.append(input_path)
    save(output / "FROZEN.json",
         {str(path): sha(path) for path in dict.fromkeys(frozen_paths)})
    save(output / "PREPARED.json", {
        "complete": True,
        "cases": len(cases),
        "new_PIC_submitted_by_prepare": False,
        "right_boundary_m": full_right,
        "right_cpml_start_m": physical_right,
        "aperture_to_cpml_clearance_m": physical_right - APERTURE_RADIUS,
        "predicted_beam_center_at_detector_m": beam_center,
        "predicted_fraction_beyond_boundary_at_detector": expected_beyond,
        "predicted_propagation_after_boundary_exit_m":
            DETECTOR_Z - boundary_crossing_z,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: prepare_100k_bend_boundary_exit.py OUTPUT LABEL1=ROOT1 LABEL2=ROOT2")
    main(sys.argv[1], sys.argv[2:])
