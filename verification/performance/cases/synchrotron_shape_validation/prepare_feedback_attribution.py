"""Prepare same-binary full/no-feedback synchrotron attribution controls."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import yaml


SHAPES = ["linear", "quintic-compensated", "linear-notch", "linear-binomial3"]
MODES = ["full-feedback", "no-feedback"]


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


def main(workspace, accepted_root, output, binary, build_record):
    workspace, accepted_root, output, binary, build_record = map(
        lambda value: Path(value).resolve(),
        (workspace, accepted_root, output, binary, build_record))
    require(not output.exists(), "Output directory already exists")
    accepted = json.loads((accepted_root / "campaign.json").read_text())
    accepted_frozen = json.loads((accepted_root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in accepted_frozen.items()),
            "Accepted campaign changed")
    require(accepted["shapes"] == SHAPES, "Unexpected accepted shape list")
    require(all((accepted_root / f"ARRAY_COMPLETE_{index}.json").is_file()
                for index in range(4)), "Accepted campaign is incomplete")
    require((accepted_root / "axis_analysis_r2/FINAL.json").is_file(),
            "Accepted axis analysis is absent")
    require((accepted_root / "srw_analysis/FINAL.json").is_file(),
            "SRW comparison is absent")
    identity = json.loads(build_record.read_text())
    require(identity.get("passed") is True and
            Path(identity["binary"]).resolve() == binary,
            "Unaccepted attribution binary identity")
    require(sha(binary) == identity["binary_sha256"],
            "Attribution binary changed")

    output.mkdir(parents=True)
    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ["prepare_feedback_attribution.py",
                 "run_feedback_attribution.py",
                 "analyze_feedback_control.py",
                 "feedback_attribution.sbatch",
                 "feedback_attribution_post.sbatch", "prepare_paired.py"]:
        shutil.copy2(here / name, scripts / name)
    shutil.copy2(here / "prepare_paired.py", scripts / "prepare.py")
    shutil.copy2(accepted_root / "scripts/launcher.py", scripts / "launcher.py")

    cases = []
    for shape in SHAPES:
        for mode in MODES:
            for kind in ["straight", "bend"]:
                source = accepted_root / shape / kind
                target = output / shape / mode / kind
                target.mkdir(parents=True)
                config = yaml.safe_load((source / "run.yaml").read_text())
                config = rewrite(copy.deepcopy(config), source, target)
                if mode == "full-feedback":
                    config["mesh"]["self_field_feedback"] = True
                    config["mesh"].pop("one_way_source_diagnostic", None)
                else:
                    config["mesh"]["self_field_feedback"] = False
                    config["mesh"]["one_way_source_diagnostic"] = True
                config["output"]["manifest"] = str(target / "run-manifest.yaml")
                config["mesh"]["moving_window"]["reference_directory"] = str(
                    target / "reference")
                config["detectors"]["directory"] = str(target / "detectors")
                config["energy_ledger"]["directory"] = str(target / "ledger")
                (target / "run.yaml").write_text(
                    yaml.safe_dump(config, sort_keys=False))
                plan = rewrite(json.loads((source / "PLAN.json").read_text()),
                               source, target)
                plan.update(case_directory=str(target), mode=mode,
                            self_field_feedback=(mode == "full-feedback"),
                            one_way_source_diagnostic=(mode == "no-feedback"),
                            particle_force=("external magnetic plus grid E/B" if
                                            mode == "full-feedback" else
                                            "prescribed external magnetic field only"),
                            maxwell_sources=("unchanged particle deposition and "
                                             "Maxwell evolution"))
                save(target / "PLAN.json", plan)
                cases.append(dict(shape=shape, mode=mode, kind=kind,
                                  directory=str(target),
                                  accepted_directory=str(source),
                                  timeout_seconds=1800))

    common = copy.deepcopy(accepted["common"])
    common.update(binary=str(binary), binary_sha256=sha(binary),
                  build_record=str(build_record))
    campaign = dict(common=common, shapes=SHAPES, modes=MODES, cases=cases,
                    binary=str(binary), binding_probe=accepted["binding_probe"],
                    accepted_full_feedback_root=str(accepted_root),
                    scope=("Same new executable for feedback-on and feedback-off; "
                           "identical input, grid, external magnet, deposition, "
                           "Maxwell solver and detectors. The feedback-off mode only "
                           "omits grid E/B from the particle force."))
    save(output / "campaign.json", campaign)
    frozen_paths = [output / "campaign.json", binary, build_record,
                    Path(common["input"]), accepted_root / "campaign.json",
                    accepted_root / "axis_analysis_r2/axis_spectra.csv",
                    accepted_root / "srw_analysis/gaussian_comparison.csv"]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json",
         {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True,
        "paired_cases": len(cases),
        "same_binary_for_both_modes": True,
        "accepted_full_feedback_reproduction_included": True,
        "only_physics_card_difference_between_modes": [
            "mesh.self_field_feedback", "mesh.one_way_source_diagnostic"],
        "initial_self_field_retained": True,
        "deposition_and_Maxwell_retained": True,
        "external_magnetic_field_retained": True,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 6,
            "usage: prepare_feedback_attribution.py WORKSPACE ACCEPTED_ROOT "
            "OUTPUT BINARY BUILD_RECORD")
    main(*sys.argv[1:])
