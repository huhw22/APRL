"""Prepare a paired synchrotron run with independent unit-electron sampling."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np
from scipy.constants import c
from scipy.special import ndtri
import yaml


SOURCE = Path("/public/home/wfollowm/APRL_Workspace/results/"
              "synchrotron_axis_notch_propagation_20261110_r1")
COUNT = 10_000
SEED = 20261110
SIGMA_XY = 20e-6
SIGMA_Z = 20e-9
GAMMA = 1174.171072271572


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


def load_physical_sampler(here):
    source_dir = here.parent / "physical_electron_source"
    sys.path.insert(0, str(source_dir))
    try:
        from source_math import random_by_id
    finally:
        sys.path.pop(0)
    return random_by_id, source_dir / "source_math.py"


def generate_input(path, here, count, seed):
    random_by_id, sampler_path = load_physical_sampler(here)
    uniform, normal = random_by_id(seed, 0, count)
    beta = math.sqrt((GAMMA - 1) * (GAMMA + 1)) / GAMMA
    dtype = np.dtype([
        ("plane_position_m", "<f8", (2,)),
        ("arrival_time_offset_s", "<f8"),
        ("proper_velocity", "<f8", (3,)),
        ("source_id", "<u8"),
        ("macro_weight", "<f8"),
    ])
    records = np.zeros(count, dtype=dtype)
    records["plane_position_m"][:, 0] = SIGMA_XY * normal[:, 0]
    records["plane_position_m"][:, 1] = SIGMA_XY * normal[:, 2]
    records["arrival_time_offset_s"] = (SIGMA_Z / (beta * c)) * ndtri(uniform)
    records["proper_velocity"][:, 2] = math.sqrt((GAMMA - 1) * (GAMMA + 1))
    records["source_id"] = np.arange(1, count + 1, dtype=np.uint64)
    records["macro_weight"] = 1.0
    require(np.isfinite(records["plane_position_m"]).all() and
            np.isfinite(records["arrival_time_offset_s"]).all(),
            "Nonfinite independent sample")
    require(len(np.unique(records["plane_position_m"], axis=0)) == count,
            "A transverse position was duplicated")
    with h5py.File(path, "x") as handle:
        group = handle.create_group("particles")
        group.attrs["format_version"] = 4
        group.attrs["scope"] = (
            "Independent physical unit-electron cold Gaussian sample; no tensor-product "
            "Gauss-Hermite rods and no forced sample moments")
        group.attrs["electron_count"] = count
        group.attrs["charge_C"] = count * 1.602176634e-19
        group.create_dataset("records", data=records)
        generation = handle.create_group("generation")
        generation["complete"] = np.uint8(1)
        generation["committed_records"] = np.uint64(count)
        generation.attrs["schema"] = "APRL_physical_electron_source_v1"
        generation.attrs["algorithm"] = (
            "fixed-n-independent-gaussian-philox4x64-boxmuller-v1")
        generation.attrs["seed"] = np.uint64(seed)
        generation.attrs["sampler_source_sha256"] = sha(sampler_path)
        generation.attrs["record_payload_sha256"] = hashlib.sha256(
            records.tobytes()).hexdigest()
    return records, sampler_path


def source_audit(records):
    xy = records["plane_position_m"]
    arrival = records["arrival_time_offset_s"]
    beta = math.sqrt((GAMMA - 1) * (GAMMA + 1)) / GAMMA
    z = arrival * beta * c
    correlation = np.corrcoef(np.c_[xy, z], rowvar=False)
    energies = np.linspace(4.0, 12.0, 321)
    omega = energies * 1.602176634e-19 / 1.054571817e-34
    sample_form = np.mean(np.exp(1j * omega[:, None] * arrival[None, :]), axis=1)
    analytic_form = np.exp(-0.5 * (omega * SIGMA_Z / (beta * c)) ** 2)
    sample_power = abs(sample_form) ** 2
    analytic_power = analytic_form ** 2
    return {
        "records": len(records),
        "unique_transverse_positions": len(np.unique(xy, axis=0)),
        "all_macro_weights_one": bool(np.all(records["macro_weight"] == 1.0)),
        "sample_mean_x_y_z_m": [float(value) for value in
                                  (xy[:, 0].mean(), xy[:, 1].mean(), z.mean())],
        "sample_sigma_x_y_z_m": [float(value) for value in
                                   (xy[:, 0].std(), xy[:, 1].std(), z.std())],
        "maximum_absolute_xyz_correlation": float(np.max(abs(
            correlation - np.eye(3)))),
        "longitudinal_form_factor_power_relative_L2_to_gaussian_4_12eV":
            float(np.linalg.norm(sample_power - analytic_power) /
                  np.linalg.norm(analytic_power)),
        "longitudinal_form_factor_power_ratio_min_4_12eV":
            float(np.min(sample_power / analytic_power)),
        "longitudinal_form_factor_power_ratio_max_4_12eV":
            float(np.max(sample_power / analytic_power)),
        "energy_grid_eV": energies.tolist(),
        "sample_longitudinal_form_factor_power": sample_power.tolist(),
        "analytic_longitudinal_form_factor_power": analytic_power.tolist(),
    }


def main(output, count=COUNT, seed=SEED):
    output = Path(output).resolve()
    require(not output.exists(), "Output already exists")
    require((SOURCE / "PIC_COMPLETE.json").is_file(), "Source campaign incomplete")
    output.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    inputs = output / "inputs"
    inputs.mkdir()
    input_path = inputs / f"independent_gaussian_{count}_electrons.h5"
    records, sampler_path = generate_input(input_path, here, count, seed)
    audit = source_audit(records)
    save(output / "SOURCE_AUDIT.json", audit)

    source_campaign = json.loads((SOURCE / "campaign.json").read_text())
    binary = Path(source_campaign["binary"])
    require(sha(binary) == source_campaign["binary_sha256"], "Binary changed")
    scripts = output / "scripts"
    scripts.mkdir()
    for name in ["prepare_independent_electron_synchrotron.py",
                 "analyze_independent_electron_synchrotron.py"]:
        shutil.copy2(here / name, scripts / name)
    for name in ["notch_propagation.sbatch", "launcher.py", "prepare_paired.py"]:
        shutil.copy2(SOURCE / "scripts" / name, scripts / name)
    shutil.copy2(here / "run_independent_electron_synchrotron.py",
                 scripts / "run_notch_propagation.py")
    # The preserved launcher imports this historical module name.
    shutil.copy2(SOURCE / "scripts" / "prepare_paired.py",
                 scripts / "prepare.py")
    shutil.copy2(sampler_path, scripts / "physical_source_math.py")

    cases = []
    for kind in ("straight", "bend"):
        source_case = SOURCE / kind
        target = output / kind
        target.mkdir()
        config = rewrite(copy.deepcopy(yaml.safe_load(
            (source_case / "run.yaml").read_text())), SOURCE, output)
        config["beam"]["input"]["file"] = str(input_path)
        config["beam"]["input"]["electrons"] = float(count)
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        plan = rewrite(json.loads((source_case / "PLAN.json").read_text()),
                       SOURCE, output)
        plan.update(
            binary=str(binary), binary_sha256=sha(binary),
            macroparticles=count, represented_electrons=float(count),
            input=str(input_path), input_sha256=sha(input_path),
            case_directory=str(target), expected_particles=count,
            source_sampling=(
                f"{count} independent unit-electron Gaussian draws; unique x/y per record; "
                "same cold momentum; no tensor-product Gauss-Hermite replication"),
            sole_physics_input_change=(
                "Replace weighted 8x8x16 GH source by an independent unit-electron "
                f"sample of the same Gaussian distribution; N={count}, with absolute "
                "continuum comparison using the corresponding explicit N^2 factor"),
            continuum_reference=(
                "N^2 Gaussian SRW continuum plus separately reported finite-N sample "
                "longitudinal form-factor deviation; no fitted normalization"),
            purpose=(
                "Test synchrotron/SRW agreement with the current independent physical-"
                "electron sampling method"))
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 1800})

    campaign = rewrite(copy.deepcopy(source_campaign), SOURCE, output)
    campaign.update(
        cases=cases, source_campaign=str(SOURCE),
        input=str(input_path), input_sha256=sha(input_path),
        particle_count=count, represented_electrons=count, seed=seed,
        sampling="independent unit-electron Gaussian Philox/Box-Muller",
        source_audit=str(output / "SOURCE_AUDIT.json"),
        scope=(
            "Exact four-plane no-feedback linear-notch bend/straight control; the only "
            "physics input change is weighted tensor-product GH -> independent physical "
            f"unit-electron sampling of the same Gaussian bunch with N={count}; all "
            "reported continuum spectra use that run's explicit N^2 factor."))
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
        "cases": 2,
        "particles": count,
        "unique_transverse_positions": audit["unique_transverse_positions"],
        "all_macro_weights_one": audit["all_macro_weights_one"],
        "same_total_electrons_as_GH_control": count == 10_000,
        "same_grid_solver_orbit_boundaries_detectors": True,
        "sole_physics_input_change": "particle sampling realization",
        "new_PIC_submitted_by_prepare": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) in (2, 4),
            "usage: prepare_independent_electron_synchrotron.py OUTPUT [COUNT SEED]")
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
