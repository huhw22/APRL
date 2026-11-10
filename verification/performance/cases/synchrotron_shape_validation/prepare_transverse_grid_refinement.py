"""Prepare a paired 2x transverse-grid refinement of the closed 10k run."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import yaml


SOURCE = Path("/public/home/wfollowm/APRL_Workspace/results/"
              "synchrotron_independent_electron_20261110_r1")


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


def main(output):
    output = Path(output).resolve()
    require(not output.exists(), "Output already exists")
    require((SOURCE / "PIC_COMPLETE.json").is_file(), "Source campaign incomplete")
    output.mkdir(parents=True)
    inputs = output / "inputs"
    inputs.mkdir()
    source_campaign = json.loads((SOURCE / "campaign.json").read_text())
    source_input = Path(source_campaign["input"])
    input_path = inputs / source_input.name
    shutil.copy2(source_input, input_path)
    shutil.copy2(SOURCE / "SOURCE_AUDIT.json", output / "SOURCE_AUDIT.json")

    scripts = output / "scripts"
    scripts.mkdir()
    for source in (SOURCE / "scripts").iterdir():
        if source.is_file():
            shutil.copy2(source, scripts / source.name)
    sbatch = scripts / "notch_propagation.sbatch"
    sbatch_text = sbatch.read_text()
    require("#SBATCH --time=01:00:00" in sbatch_text,
            "Unexpected source batch time limit")
    sbatch.write_text(sbatch_text.replace("#SBATCH --time=01:00:00",
                                          "#SBATCH --time=03:00:00"))

    binary = Path(source_campaign["binary"])
    require(sha(binary) == source_campaign["binary_sha256"], "Binary changed")
    cases = []
    for kind in ("straight", "bend"):
        source_case = SOURCE / kind
        target = output / kind
        target.mkdir()
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source_case / "run.yaml").read_text())), SOURCE, output)
        config["mesh"]["cells"][0:2] = [320, 256]
        config["mesh"]["cell_size"][0:2] = [2.5e-6, 2.5e-6]
        config["boundary"]["cells"][0:2] = [20, 20]
        config["beam"]["input"]["file"] = str(input_path)
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

        plan = rewrite(json.loads((source_case / "PLAN.json").read_text()),
                       SOURCE, output)
        plan["geometry"].update(nx=320, ny=256)
        plan.update(
            input=str(input_path), input_sha256=sha(input_path),
            case_directory=str(target), timeout_seconds=5400,
            transverse_grid_refinement={
                "factor": 2,
                "source_dx_dy_m": [5e-6, 5e-6],
                "refined_dx_dy_m": [2.5e-6, 2.5e-6],
                "source_nx_ny": [160, 128],
                "refined_nx_ny": [320, 256],
                "physical_extent_invariant_m": [8e-4, 6.4e-4],
                "source_cpml_cells": [10, 10],
                "refined_cpml_cells": [20, 20],
                "cpml_physical_thickness_invariant": True,
            },
            sole_numerical_change=(
                "Halve transverse dx/dy while doubling nx/ny and transverse CPML "
                "cells; preserve physical extent, CPML thickness, particle input, "
                "dz/dt, active longitudinal window, solver and detectors"),
            purpose=(
                "Separate transverse-grid response from fixed 10k particle sampling"))
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 5400})

    campaign = rewrite(copy.deepcopy(source_campaign), SOURCE, output)
    campaign.update(
        cases=cases, source_campaign=str(SOURCE), input=str(input_path),
        input_sha256=sha(input_path),
        scope=(
            "Exact closed 10k independent-electron input on a 2x refined "
            "transverse grid. Physical transverse extent and CPML thickness, dz/dt, "
            "active longitudinal window, orbit, source coupling and detectors are "
            "held invariant."),
        transverse_grid_refinement_factor=2)
    save(output / "campaign.json", campaign)

    frozen_paths = [output / "campaign.json", output / "SOURCE_AUDIT.json",
                    input_path, binary]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json", {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True,
        "source_campaign": str(SOURCE),
        "same_particle_input_sha256": sha(input_path) == sha(source_input),
        "same_particle_count": True,
        "same_dz_dt_longitudinal_window": True,
        "same_physical_transverse_extent": True,
        "same_physical_cpml_thickness": True,
        "transverse_grid_refinement_factor": 2,
        "new_PIC_submitted_by_prepare": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 2,
            "usage: prepare_transverse_grid_refinement.py OUTPUT")
    main(sys.argv[1])
