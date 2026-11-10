"""Prepare a deposited-current-filter control for the synchrotron axis notch."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import yaml


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rewrite(value, old, new):
    if isinstance(value, dict):
        return {key: rewrite(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, old, new) for item in value]
    if isinstance(value, str):
        return value.replace(str(old), str(new))
    return value


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def main(unfiltered_root, output):
    unfiltered_root, output = map(lambda value: Path(value).resolve(),
                                  (unfiltered_root, output))
    require(not output.exists(), "Output already exists")
    require((unfiltered_root / "PIC_COMPLETE.json").is_file(),
            "Unfiltered propagation campaign is incomplete")
    campaign = json.loads((unfiltered_root / "campaign.json").read_text())
    frozen = json.loads((unfiltered_root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Unfiltered frozen inputs changed")
    binary = Path(campaign["binary"])
    require(sha(binary) == campaign["binary_sha256"], "Binary changed")

    output.mkdir(parents=True)
    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ["prepare_notch_current_filter.py",
                 "run_notch_propagation.py", "notch_propagation.sbatch",
                 "prepare_paired.py"]:
        shutil.copy2(here / name, scripts / name)
    shutil.copy2(here / "prepare_paired.py", scripts / "prepare.py")
    shutil.copy2(unfiltered_root / "scripts/launcher.py", scripts / "launcher.py")

    cases = []
    for kind in ["straight", "bend"]:
        source = unfiltered_root / kind
        target = output / kind
        target.mkdir()
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source / "run.yaml").read_text())), source, target)
        require(config["mesh"]["current_filter"] == "none",
                "Source current filter is not none")
        config["mesh"]["current_filter"] = "compensated-bilinear"
        config["detectors"]["directory"] = str(target / "detectors")
        config["energy_ledger"]["directory"] = str(target / "ledger")
        config["mesh"]["moving_window"]["reference_directory"] = str(
            target / "reference")
        config["output"]["manifest"] = str(target / "run-manifest.yaml")
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        plan = rewrite(json.loads((source / "PLAN.json").read_text()),
                       source, target)
        plan.update(
            case_directory=str(target),
            current_filter="compensated-bilinear",
            sole_physics_card_change="mesh.current_filter: none -> compensated-bilinear",
            purpose="Test whether deposited transverse Nyquist current creates the axis notch")
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 1800})

    control = copy.deepcopy(campaign)
    control.update(
        cases=cases,
        source_unfiltered_root=str(unfiltered_root),
        current_filter="compensated-bilinear",
        scope=("Exact clone of the closed four-plane no-feedback 7-point pair; "
               "the sole physics-card change is mesh.current_filter from none "
               "to compensated-bilinear. Particle input, external orbit, grid, "
               "Cowan field solver, boundaries and detectors are unchanged."))
    save(output / "campaign.json", control)
    frozen_paths = [output / "campaign.json", binary]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json", {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True,
        "cases": len(cases),
        "sole_physics_card_change": "mesh.current_filter",
        "source_push_grid_solver_boundaries_unchanged": True,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 3,
            "usage: prepare_notch_current_filter.py UNFILTERED_ROOT OUTPUT")
    main(*sys.argv[1:])
