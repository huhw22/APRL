"""Prepare a 12x12 transverse Gauss-Hermite axis-notch source control."""
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
from scipy.special import roots_hermitenorm
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


def input_records(transverse_order, longitudinal_order, sigma_xy, sigma_z,
                  gamma):
    nodes, weights = [], []
    for order in [transverse_order, transverse_order, longitudinal_order]:
        node, weight = roots_hermitenorm(order)
        nodes.append(node)
        weights.append(weight / math.sqrt(2 * math.pi))
    xx, yy, zz = np.meshgrid(*nodes, indexing="ij")
    wx, wy, wz = np.meshgrid(*weights, indexing="ij")
    weight = (wx * wy * wz).ravel()
    weight /= weight.sum()
    count = len(weight)
    dtype = np.dtype([("plane_position_m", "<f8", (2,)),
                      ("arrival_time_offset_s", "<f8"),
                      ("proper_velocity", "<f8", (3,)),
                      ("source_id", "<u8"), ("macro_weight", "<f8")])
    particles = np.zeros(count, dtype=dtype)
    particles["plane_position_m"] = np.c_[xx.ravel() * sigma_xy,
                                            yy.ravel() * sigma_xy]
    beta = math.sqrt((gamma - 1) * (gamma + 1)) / gamma
    particles["arrival_time_offset_s"] = zz.ravel() * sigma_z / (beta * c)
    particles["proper_velocity"][:, 2] = math.sqrt((gamma - 1) * (gamma + 1))
    particles["source_id"] = np.arange(1, count + 1)
    particles["macro_weight"] = weight * count
    return particles


def source_alias_fraction(order, sigma, nx, ny, dx, dy):
    node, weight = roots_hermitenorm(order)
    weight /= math.sqrt(2 * math.pi)
    kx = 2 * math.pi * np.fft.fftfreq(nx, dx)
    ky = 2 * math.pi * np.fft.fftfreq(ny, dy)
    fx = np.sum(weight[:, None] * np.exp(
        -1j * node[:, None] * sigma * kx[None, :]), axis=0)
    fy = np.sum(weight[:, None] * np.exp(
        -1j * node[:, None] * sigma * ky[None, :]), axis=0)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    # First-order transverse deposition has the static tent response below.
    response = (np.sinc(kx_grid * dx / (2 * math.pi))**2 *
                np.sinc(ky_grid * dy / (2 * math.pi))**2)
    power = abs(fy[:, None] * fx[None, :] * response)**2
    near = ((abs(kx_grid) > .7 * math.pi / dx) |
            (abs(ky_grid) > .7 * math.pi / dy))
    return float(power[near].sum() / power.sum())


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

    source_card = yaml.safe_load((unfiltered_root / "bend/run.yaml").read_text())
    gamma = 1174.171072271572
    sigma_xy, sigma_z = 20e-6, 20e-9
    old_order, new_order, longitudinal_order = 8, 12, 16
    particles = input_records(new_order, longitudinal_order, sigma_xy,
                              sigma_z, gamma)

    output.mkdir(parents=True)
    inputs = output / "inputs"
    inputs.mkdir()
    input_path = inputs / "gaussian_plane_events_12x12x16.h5"
    with h5py.File(input_path, "x") as handle:
        group = handle.create_group("particles")
        group.attrs["format_version"] = 4
        group.attrs["scope"] = (
            "Product Gauss-Hermite weighted Gaussian continuum; 12x12 "
            "transverse and 16 longitudinal nodes; cold parallel fixed-lab-plane events")
        group.create_dataset("records", data=particles)

    # Longitudinal support and momentum must be unchanged from the 8x8x16 input.
    old_input = Path(source_card["beam"]["input"]["file"])
    with h5py.File(old_input, "r") as handle:
        old = handle["particles/records"][:]
    require(np.array_equal(np.unique(old["arrival_time_offset_s"]),
                           np.unique(particles["arrival_time_offset_s"])),
            "Longitudinal nodes changed")
    require(np.array_equal(np.unique(old["proper_velocity"][:, 2]),
                           np.unique(particles["proper_velocity"][:, 2])),
            "Beam momentum changed")

    scripts = output / "scripts"
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ["prepare_notch_transverse_sampling.py",
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
        config["beam"]["input"]["file"] = str(input_path)
        config["detectors"]["directory"] = str(target / "detectors")
        config["energy_ledger"]["directory"] = str(target / "ledger")
        config["mesh"]["moving_window"]["reference_directory"] = str(
            target / "reference")
        config["output"]["manifest"] = str(target / "run-manifest.yaml")
        (target / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        plan = rewrite(json.loads((source / "PLAN.json").read_text()),
                       source, target)
        plan.update(
            case_directory=str(target), expected_particles=len(particles),
            macroparticles=len(particles), input=str(input_path),
            input_sha256=sha(input_path), transverse_order=new_order,
            longitudinal_order=longitudinal_order,
            sole_physics_input_change="transverse Gauss-Hermite order 8x8 -> 12x12",
            purpose="Test whether sparse transverse macro sampling injects the axis-notch high-k mode")
        save(target / "PLAN.json", plan)
        cases.append({"kind": kind, "directory": str(target),
                      "timeout_seconds": 1800})

    nx, ny = source_card["mesh"]["cells"][:2]
    dx, dy = source_card["mesh"]["cell_size"][:2]
    alias = {
        "definition": ("Fraction of discrete transverse input form-factor power "
                       "above 0.7 transverse Nyquist after the static first-order "
                       "tent-deposition response; diagnostic, not a propagated-field prediction."),
        "order_8": source_alias_fraction(old_order, sigma_xy, nx, ny, dx, dy),
        "order_12": source_alias_fraction(new_order, sigma_xy, nx, ny, dx, dy),
    }
    control = copy.deepcopy(campaign)
    control.update(
        cases=cases, source_unfiltered_root=str(unfiltered_root),
        input=str(input_path), input_sha256=sha(input_path),
        transverse_order=new_order, longitudinal_order=longitudinal_order,
        macroparticles=len(particles), static_source_alias=alias,
        scope=("Exact clone of the closed four-plane no-feedback 7-point pair; "
               "only transverse Gauss-Hermite sampling changes from 8x8 to 12x12. "
               "Longitudinal nodes, total represented charge, orbit, deposition, "
               "Cowan solver, grid, boundaries and detectors are unchanged."))
    save(output / "campaign.json", control)
    frozen_paths = [output / "campaign.json", input_path, binary]
    frozen_paths += sorted(path for path in scripts.iterdir() if path.is_file())
    for case in cases:
        directory = Path(case["directory"])
        frozen_paths += [directory / "run.yaml", directory / "PLAN.json"]
    save(output / "FROZEN.json", {str(path): sha(path) for path in frozen_paths})
    save(output / "PREPARED.json", {
        "complete": True, "cases": len(cases),
        "old_macroparticles": len(old), "new_macroparticles": len(particles),
        "sole_physics_input_change": "transverse Gauss-Hermite order",
        "longitudinal_nodes_and_momentum_bitwise_unchanged": True,
        "deposition_grid_solver_boundaries_unchanged": True,
        "static_source_alias": alias,
        "production_default_changed": False,
    })
    print(output)


if __name__ == "__main__":
    require(len(sys.argv) == 3,
            "usage: prepare_notch_transverse_sampling.py UNFILTERED_ROOT OUTPUT")
    main(*sys.argv[1:])
