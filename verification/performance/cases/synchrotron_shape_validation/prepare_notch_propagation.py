"""Prepare a four-plane no-feedback 7-point notch propagation control."""
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


def main(feedback_root, output):
    feedback_root, output = map(lambda value: Path(value).resolve(),
                                (feedback_root, output))
    require(not output.exists(), "Output already exists")
    campaign = json.loads((feedback_root / "campaign.json").read_text())
    frozen = json.loads((feedback_root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Feedback campaign inputs changed")
    binary = Path(campaign["binary"])
    require(sha(binary) == campaign["common"]["binary_sha256"],
            "Binary changed")
    output.mkdir(parents=True)
    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ["prepare_notch_propagation.py", "run_notch_propagation.py",
                 "notch_propagation.sbatch", "prepare_paired.py"]:
        shutil.copy2(here / name, scripts / name)
    shutil.copy2(here / "prepare_paired.py", scripts / "prepare.py")
    shutil.copy2(feedback_root / "scripts/launcher.py", scripts / "launcher.py")

    planes = [("z0060", .0600), ("z0100", .1000),
              ("z0200", .2000), ("z03524", .3524)]
    cases = []
    for kind in ["straight", "bend"]:
        source = feedback_root / "linear-notch/no-feedback" / kind
        target = output / kind
        target.mkdir()
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source / "run.yaml").read_text())), source, target)
        original = config["detectors"]["field_planes"][0]
        config["detectors"]["field_planes"] = [
            dict(copy.deepcopy(original), name=name, z=z) for name, z in planes]
        config["detectors"]["directory"] = str(target / "detectors")
        config["energy_ledger"]["directory"] = str(target / "ledger")
        config["mesh"]["moving_window"]["reference_directory"] = str(
            target / "reference")
        config["output"]["manifest"] = str(target / "run-manifest.yaml")
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        plan = rewrite(json.loads((source / "PLAN.json").read_text()),
                       source, target)
        plan.update(case_directory=str(target), field_planes=[
            {"name": name, "z_m": z} for name, z in planes],
            purpose="Track axis notches and transverse high-k content with distance")
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 1800})
    control = {
        "binary": str(binary),
        "binary_sha256": sha(binary),
        "binding_probe": campaign["binding_probe"],
        "source_feedback_root": str(feedback_root),
        "shape": "linear-notch",
        "mode": "no-feedback",
        "planes": [{"name": name, "z_m": z} for name, z in planes],
        "cases": cases,
        "scope": ("Exact clone of the closed no-feedback 7-point bend/straight pair; "
                  "only three read-only field planes are added before the original screen."),
    }
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
        "added_read_only_field_planes": len(planes) - 1,
        "source_push_deposition_grid_boundaries_unchanged": True,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 3,
            "usage: prepare_notch_propagation.py FEEDBACK_ROOT OUTPUT")
    main(*sys.argv[1:])
