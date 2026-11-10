"""Prepare the same-geometry one-electron bend/straight SRW control."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np
import yaml


SOURCE = Path("/public/home/wfollowm/APRL_Workspace/results/"
              "synchrotron_independent_electron_20261110_r1")
GAMMA = 1174.171072271572
REFERENCE_PADDING_CELLS = 64


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
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def rewrite(value, old, new):
    if isinstance(value, dict):
        return {key: rewrite(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, old, new) for item in value]
    if isinstance(value, str):
        return value.replace(str(old), str(new))
    return value


def generate_input(path):
    dtype = np.dtype([
        ("plane_position_m", "<f8", (2,)),
        ("arrival_time_offset_s", "<f8"),
        ("proper_velocity", "<f8", (3,)),
        ("source_id", "<u8"),
        ("macro_weight", "<f8"),
    ])
    records = np.zeros(1, dtype=dtype)
    records["proper_velocity"][0, 2] = math.sqrt((GAMMA - 1) * (GAMMA + 1))
    records["source_id"] = 1
    records["macro_weight"] = 1.0
    with h5py.File(path, "x") as handle:
        group = handle.create_group("particles")
        group.attrs["format_version"] = 4
        group.attrs["scope"] = (
            "One physical electron on the reference axis; cold longitudinal "
            "momentum; fixed laboratory input-plane event")
        group.attrs["electron_count"] = 1
        group.attrs["charge_C"] = 1.602176634e-19
        group.create_dataset("records", data=records)
        generation = handle.create_group("generation")
        generation["complete"] = np.uint8(1)
        generation["committed_records"] = np.uint64(1)
        generation.attrs["schema"] = "APRL_single_physical_electron_v1"
        generation.attrs["record_payload_sha256"] = hashlib.sha256(
            records.tobytes()).hexdigest()
    return records


def main(output):
    output = Path(output).resolve()
    require(not output.exists(), "Output already exists")
    require((SOURCE / "PIC_COMPLETE.json").is_file(), "Source campaign incomplete")
    output.mkdir(parents=True)
    inputs = output / "inputs"
    inputs.mkdir()
    input_path = inputs / "single_axis_electron.h5"
    generate_input(input_path)

    source_campaign = json.loads((SOURCE / "campaign.json").read_text())
    binary = Path(source_campaign["binary"])
    require(sha(binary) == source_campaign["binary_sha256"], "Binary changed")
    here = Path(__file__).resolve().parent
    scripts = output / "scripts"
    scripts.mkdir()
    for name in ["prepare_single_electron_synchrotron.py",
                 "run_single_electron_case.py",
                 "single_electron_case.sbatch",
                 "analyze_single_electron_synchrotron.py"]:
        shutil.copy2(here / name, scripts / name)
    for name in ["launcher.py", "prepare_paired.py"]:
        shutil.copy2(SOURCE / "scripts" / name, scripts / name)
    shutil.copy2(SOURCE / "scripts" / "prepare.py", scripts / "prepare.py")

    cases = []
    for kind in ("straight", "bend"):
        source_case = SOURCE / kind
        target = output / kind
        target.mkdir()
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source_case / "run.yaml").read_text())), SOURCE, output)
        config["beam"]["input"]["file"] = str(input_path)
        config["beam"]["input"]["electrons"] = 1.0
        dz = float(config["mesh"]["cell_size"][2])
        config["mesh"]["cells"][2] = (
            int(config["mesh"]["cells"][2]) + REFERENCE_PADDING_CELLS)
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        plan = rewrite(json.loads((source_case / "PLAN.json").read_text()),
                       SOURCE, output)
        plan.update(
            binary=str(binary), binary_sha256=sha(binary),
            macroparticles=1, represented_electrons=1.0,
            input=str(input_path), input_sha256=sha(input_path),
            case_directory=str(target), expected_particles=1,
            source_sampling="one physical electron at x=y=t_offset=0",
            sole_physics_input_change=(
                "Replace the 10000-electron Gaussian realization by one on-axis "
                "physical electron; grid, boost, solver, magnet, boundaries, "
                "detectors and no-feedback mode are unchanged"),
            continuum_reference="native single-electron SRW; no bunch form factor",
            purpose=(
                "Provide a no-form-factor APRL/SRW control and isolate particle "
                "sampling from the field/detector chain"),
            reference_source_mesh_padding_cells=REFERENCE_PADDING_CELLS,
            reference_source_mesh_padding_m=REFERENCE_PADDING_CELLS * dz,
            reference_translation_note=(
                "The explicit full reference/source mesh gains 32 native z cells "
                "on each side to cover the single-event synchronization. Because "
                "the compact moving window retains its explicit 497 cells and its "
                "integer start offset increases by the 32-cell lower padding, the "
                "active E/B window physical location, resolution, magnet and "
                "detectors are unchanged."))
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 1800})

    campaign = rewrite(copy.deepcopy(source_campaign), SOURCE, output)
    campaign.update(
        cases=cases, source_campaign=str(SOURCE), input=str(input_path),
        input_sha256=sha(input_path), particle_count=1,
        represented_electrons=1, seed=None,
        sampling="one on-axis unit electron",
        reference_source_mesh_padding_cells=REFERENCE_PADDING_CELLS,
        scope=(
            "Same four-plane no-feedback linear-notch bend/straight control with "
            "one physical electron. This isolates the APRL field/detector chain "
            "from Gaussian sampling and all bunch form factors. The finite full "
            "reference/source mesh is padded by 32 native cells on each side; the "
            "compact active E/B window remains 497 cells and is invariant under "
            "the paired lower-origin/integer-offset change."))
    save(output / "campaign.json", campaign)
    frozen_paths = [output / "campaign.json", input_path, binary]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json", {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True,
        "cases": 2,
        "particles": 1,
        "macro_weight": 1.0,
        "same_grid_solver_orbit_boundaries_detectors": True,
        "reference_source_mesh_padding_cells": REFERENCE_PADDING_CELLS,
        "active_moving_window_expected_unchanged": True,
        "bunch_form_factor_applied": False,
        "new_PIC_submitted_by_prepare": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 2,
            "usage: prepare_single_electron_synchrotron.py OUTPUT")
    main(sys.argv[1])

