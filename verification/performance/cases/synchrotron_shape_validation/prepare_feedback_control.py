"""Prepare strict no-grid-feedback controls from the closed synchrotron campaign."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import yaml


NAMES = ["linear", "quintic-compensated", "linear-notch", "linear-binomial3"]


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


def main(workspace, full_root, output):
    workspace, full_root, output = map(lambda x: Path(x).resolve(),
                                       (workspace, full_root, output))
    require(not output.exists(), "Output directory already exists")
    campaign = json.loads((full_root / "campaign.json").read_text())
    frozen = json.loads((full_root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Closed full-feedback campaign changed")
    require(campaign["shapes"] == NAMES, "Unexpected shape list")
    require(all((full_root / f"ARRAY_COMPLETE_{index}.json").is_file()
                for index in range(4)), "Full-feedback campaign is incomplete")
    require((full_root / "axis_analysis_r2/FINAL.json").is_file(),
            "Validated full-feedback axis analysis is absent")
    require((full_root / "srw_analysis/FINAL.json").is_file(),
            "SRW comparison is absent")
    binary = Path(campaign["binary"])
    require(sha(binary) == campaign["common"]["binary_sha256"],
            "Simulation binary changed")

    output.mkdir(parents=True)
    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ["prepare_feedback_control.py", "run_feedback_control.py",
                 "analyze_feedback_control.py", "feedback_control.sbatch",
                 "feedback_post.sbatch", "prepare_paired.py"]:
        shutil.copy2(here / name, scripts / name)
    # The historic launcher imports helpers under this compatibility name.
    shutil.copy2(here / "prepare_paired.py", scripts / "prepare.py")
    # The launcher is a preserved campaign helper rather than a source template.
    shutil.copy2(full_root / "scripts/launcher.py", scripts / "launcher.py")

    cases = []
    for shape in NAMES:
        for kind in ["straight", "bend"]:
            source = full_root / shape / kind
            target = output / shape / kind
            target.mkdir(parents=True)
            config = yaml.safe_load((source / "run.yaml").read_text())
            config = rewrite(copy.deepcopy(config), source, target)
            config["mesh"]["self_field_feedback"] = False
            config["output"]["manifest"] = str(target / "run-manifest.yaml")
            config["mesh"]["moving_window"]["reference_directory"] = str(
                target / "reference")
            config["detectors"]["directory"] = str(target / "detectors")
            config["energy_ledger"]["directory"] = str(target / "ledger")
            (target / "run.yaml").write_text(
                yaml.safe_dump(config, sort_keys=False))
            original_plan = json.loads((source / "PLAN.json").read_text())
            plan = rewrite(original_plan, source, target)
            plan.update(case_directory=str(target), self_field_feedback=False,
                        particle_force="prescribed external magnetic field only",
                        maxwell_sources="unchanged particle deposition and Maxwell evolution")
            save(target / "PLAN.json", plan)
            cases.append(dict(shape=shape, kind=kind, directory=str(target),
                              source_full_feedback_directory=str(source),
                              timeout_seconds=1800))

    control = copy.deepcopy(campaign)
    control["cases"] = cases
    control["full_feedback_root"] = str(full_root)
    control["mode"] = "self_field_feedback_false"
    control["scope"] = (
        "Same input, grid, external magnet, deposition and Maxwell solver as the closed "
        "four-shape campaign; only grid E/B force on particles is disabled.")
    save(output / "campaign.json", control)

    frozen_paths = [output / "campaign.json", binary,
                    Path(campaign["common"]["input"]),
                    full_root / "campaign.json",
                    full_root / "axis_analysis_r2/axis_spectra.csv",
                    full_root / "srw_analysis/gaussian_comparison.csv"]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json", {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True,
        "cases": len(cases),
        "only_physics_card_change": "mesh.self_field_feedback: true(default) -> false",
        "initial_self_field_retained": True,
        "deposition_and_Maxwell_retained": True,
        "external_magnetic_field_retained": True,
        "binary_rebuilt": False,
        "source_core_changed": False,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: prepare_feedback_control.py WORKSPACE FULL_ROOT OUTPUT")
    main(*sys.argv[1:])
