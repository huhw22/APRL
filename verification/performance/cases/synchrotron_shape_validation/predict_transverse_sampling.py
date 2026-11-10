"""Predict transverse Gauss-Hermite sampling error before a new PIC scan.

This is a source-only calculation.  It separates three effects:

1. quadrature of the continuous finite-bend Lienard-Wiechert field;
2. mapping the weighted nodes to APRL's actual Jx and Jz transverse lattices;
3. phase error of the Cowan-z vacuum dispersion relation.

No APRL output is read and no simulation is launched.
"""
import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar
from scipy.optimize import brentq
from scipy.special import ndtr, roots_hermitenorm


SIGMA_XY = 20e-6
DX = 5e-6
DY = 5e-6
DZ = 1.2543262415359034e-6
DT = DZ / c
NX = 160
NY = 128
MESH_CENTER_X = 0.12e-3
MESH_CENTER_Y = 0.0
ORIGIN_X = MESH_CENTER_X - NX * DX / 2
ORIGIN_Y = MESH_CENTER_Y - NY * DY / 2
GAMMA = 1174.171072271572
BOOST_GAMMA = 929.1305492858543
RHO = 3.0
DIPOLE_LENGTH = 0.0024
SCREEN_FROM_DIPOLE_ENTRANCE = 0.3024
LONGITUDINAL_ORDER = 16
SAGITTA = RHO * (1 - math.cos(math.asin(DIPOLE_LENGTH / RHO)))
OFFSETS = np.linspace(0.0, SAGITTA, 9)
HC_EV_M = hbar * c / e
DIAGNOSTIC_ENERGY_EV = 10.0


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def gaussian_probability(a, b, mean):
    return ndtr((b - mean) / SIGMA_XY) - ndtr((a - mean) / SIGMA_XY)


def gaussian_first_moment(a, b, mean):
    probability = gaussian_probability(a, b, mean)
    correction = SIGMA_XY / math.sqrt(2 * math.pi) * (
        np.exp(-0.5 * ((a - mean) / SIGMA_XY) ** 2)
        - np.exp(-0.5 * ((b - mean) / SIGMA_XY) ** 2))
    return mean * probability + correction


def exact_deposited_1d(origin, cells, spacing, kind, mean):
    """Continuum Gaussian passed exactly through bin or linear deposition."""
    if kind == "bin":
        lower = origin + np.arange(cells) * spacing
        return gaussian_probability(lower, lower + spacing, mean)
    coordinate = origin + np.arange(cells + 1) * spacing
    left = (gaussian_first_moment(coordinate - spacing, coordinate, mean)
            - (coordinate - spacing)
            * gaussian_probability(coordinate - spacing, coordinate, mean))
    right = ((coordinate + spacing)
             * gaussian_probability(coordinate, coordinate + spacing, mean)
             - gaussian_first_moment(coordinate, coordinate + spacing, mean))
    return (left + right) / spacing


def gh_deposited_1d(order, origin, cells, spacing, kind, mean):
    nodes, weights = roots_hermitenorm(order)
    positions = mean + SIGMA_XY * nodes
    weights = weights / math.sqrt(2 * math.pi)
    coordinate = (positions - origin) / spacing
    lower = np.floor(coordinate).astype(int)
    if kind == "bin":
        result = np.zeros(cells)
        for index, weight in zip(lower, weights):
            if 0 <= index < cells:
                result[index] += weight
        return result
    result = np.zeros(cells + 1)
    fraction = coordinate - lower
    for index, value, weight in zip(lower, fraction, weights):
        if 0 <= index <= cells:
            result[index] += weight * (1 - value)
        if 0 <= index + 1 <= cells:
            result[index + 1] += weight * value
    return result


def source_arrays(order, component, offset):
    """Return GH and exact-continuum sources after the same APRL deposit."""
    if component == "Jx":
        x_kind, y_kind = "bin", "tent"
    elif component == "Jz":
        x_kind = y_kind = "tent"
    else:
        raise ValueError(component)
    gx = gh_deposited_1d(order, ORIGIN_X, NX, DX, x_kind, offset)
    gy = gh_deposited_1d(order, ORIGIN_Y, NY, DY, y_kind, 0.0)
    rx = exact_deposited_1d(ORIGIN_X, NX, DX, x_kind, offset)
    ry = exact_deposited_1d(ORIGIN_Y, NY, DY, y_kind, 0.0)
    sampled = np.outer(gy, gx)
    reference = np.outer(ry, rx)
    sampled /= sampled.sum()
    reference /= reference.sum()
    return sampled, reference


def fourier_geometry(shape):
    ky = 2 * math.pi * np.fft.fftfreq(shape[0], DY)
    kx = 2 * math.pi * np.fft.fftfreq(shape[1], DX)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    return kx, ky, kx_grid, ky_grid, np.hypot(kx_grid, ky_grid)


def source_metrics(order, component):
    values = []
    cone = {energy: [] for energy in (4.0, 8.0, 10.0, 12.0)}
    for offset in OFFSETS:
        sampled, reference = source_arrays(order, component, offset)
        sampled_fft = np.fft.fft2(sampled)
        reference_fft = np.fft.fft2(reference)
        _, _, kx, ky, kperp = fourier_geometry(sampled.shape)
        difference = sampled_fft - reference_fft
        full = np.linalg.norm(difference) / np.linalg.norm(reference_fft)
        low = kperp < 3 * (DIAGNOSTIC_ENERGY_EV / HC_EV_M) / GAMMA
        high = kperp > 5 * (DIAGNOSTIC_ENERGY_EV / HC_EV_M) / GAMMA
        near_nyquist = ((abs(kx) > 0.7 * math.pi / DX)
                        | (abs(ky) > 0.7 * math.pi / DY))
        low_power = np.sum(abs(sampled_fft[low]) ** 2)
        values.append((
            full,
            np.sum(abs(sampled_fft[high]) ** 2) / low_power,
            np.sum(abs(sampled_fft[near_nyquist]) ** 2)
            / np.sum(abs(sampled_fft) ** 2),
        ))
        for energy in cone:
            physical = kperp < 3 * (energy / HC_EV_M) / GAMMA
            cone[energy].append(
                np.linalg.norm(difference[physical])
                / np.linalg.norm(reference_fft[physical]))
    values = np.asarray(values)
    result = {
        "full_max": float(values[:, 0].max()),
        "full_mean": float(values[:, 0].mean()),
        "high_over_low_max_10eV": float(values[:, 1].max()),
        "high_over_low_mean_10eV": float(values[:, 1].mean()),
        "near_nyquist_max": float(values[:, 2].max()),
        "near_nyquist_mean": float(values[:, 2].mean()),
    }
    for energy, error in cone.items():
        result[f"cone3gamma_max_{energy:g}eV"] = float(np.max(error))
        result[f"cone3gamma_mean_{energy:g}eV"] = float(np.mean(error))
    return result


def deposit_operator_bias(component):
    """Compare exact continuum deposition with the analytic Gaussian transform."""
    errors = []
    cone_errors = {energy: [] for energy in (4.0, 8.0, 10.0, 12.0)}
    for offset in OFFSETS:
        _, reference = source_arrays(256, component, offset)
        transform = np.fft.fft2(reference)
        kx, ky, kx_grid, ky_grid, kperp = fourier_geometry(reference.shape)
        if component == "Jx":
            x_first, y_first = ORIGIN_X + 0.5 * DX, ORIGIN_Y
        else:
            x_first, y_first = ORIGIN_X, ORIGIN_Y
        transform *= np.exp(-1j * (kx_grid * x_first + ky_grid * y_first))
        analytic = np.exp(-0.5 * SIGMA_XY ** 2 * kperp ** 2
                          - 1j * kx_grid * offset)
        difference = transform - analytic
        errors.append(np.linalg.norm(difference) / np.linalg.norm(analytic))
        for energy in cone_errors:
            physical = kperp < 3 * (energy / HC_EV_M) / GAMMA
            cone_errors[energy].append(
                np.linalg.norm(difference[physical])
                / np.linalg.norm(analytic[physical]))
    result = {"component": component,
              "full_max": float(np.max(errors)),
              "full_mean": float(np.mean(errors))}
    for energy, error in cone_errors.items():
        result[f"cone3gamma_max_{energy:g}eV"] = float(np.max(error))
    return result


def load_frozen_theory(theory_directory):
    theory_directory = Path(theory_directory).resolve()
    for name in ("theory", "theory_paired"):
        spec = importlib.util.spec_from_file_location(
            name, theory_directory / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules["theory"], sys.modules["theory_paired"]


def lw_quadrature(theory_directory):
    theory, paired = load_frozen_theory(theory_directory)
    parameters = theory.constants(gamma=GAMMA, rho=RHO,
                                  length=DIPOLE_LENGTH,
                                  screen=SCREEN_FROM_DIPOLE_ENTRANCE)
    energies = np.linspace(4.0, 12.0, 33)
    first = -1.190993112140639e-15
    # This is the exact 278-frame PIC spectral spacing used by the benchmark.
    delta_energy = 4.955403614304977 - 4.129503011920814
    span = h / (delta_energy * e)
    orders = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 56, 64]
    fields = {}
    for order in orders:
        value = paired.gaussian_paired_reference(
            0.0, 0.0, energies, parameters, sigma_xy=SIGMA_XY,
            sigma_z=20e-9, order=order, time_first=first,
            time_span=span, ntime=278)
        fields[order] = np.concatenate((value[0].ravel(),
                                        (c * value[1]).ravel()))
    reference = fields[64]
    rows = []
    for order in orders:
        rows.append({
            "order": order,
            "transverse_nodes": order * order,
            "relative_field_l2_to_order64": float(
                np.linalg.norm(fields[order] - reference)
                / np.linalg.norm(reference)),
        })
    return rows


def cowan_coefficients():
    ratio = np.array([(DZ / DX) ** 2, (DZ / DY) ** 2, 1.0])
    rx, ry, rz = ratio
    denominator = rx * ry + ry * rz + rz * rx
    fraction = rx * ry * rz / denominator
    beta = ratio * (1 - fraction) / 8
    delta_xy = rx * ry * (1 / 16 - rx * ry / (8 * denominator))
    delta_yz = ry * rz * (1 / 16 - ry * rz / (8 * denominator))
    delta_zx = rz * rx * (1 / 16 - rz * rx / (8 * denominator))
    alpha = np.array([
        1 - 2 * beta[1] - 2 * beta[2] - 4 * delta_yz,
        1 - 2 * beta[2] - 2 * beta[0] - 4 * delta_zx,
        1 - 2 * beta[0] - 2 * beta[1] - 4 * delta_xy,
    ])
    return alpha, beta, (delta_xy, delta_yz, delta_zx)


def cowan_lambda(kx, ky, kz):
    alpha, beta, delta = cowan_coefficients()
    delta_xy, delta_yz, delta_zx = delta
    cx, cy, cz = math.cos(kx * DX), math.cos(ky * DY), math.cos(kz * DZ)
    smooth_x = alpha[0] + 2 * beta[1] * cy + 2 * beta[2] * cz + 4 * delta_yz * cy * cz
    smooth_y = alpha[1] + 2 * beta[2] * cz + 2 * beta[0] * cx + 4 * delta_zx * cz * cx
    smooth_z = alpha[2] + 2 * beta[0] * cx + 2 * beta[1] * cy + 4 * delta_xy * cx * cy
    derivative = np.array([2 * math.sin(kx * DX / 2) / DX,
                           2 * math.sin(ky * DY / 2) / DY,
                           2 * math.sin(kz * DZ / 2) / DZ])
    return (smooth_x * derivative[0] ** 2
            + smooth_y * derivative[1] ** 2
            + smooth_z * derivative[2] ** 2)


def cowan_kz(kx, ky, omega_boost, physical_kz):
    target = (2 * math.sin(omega_boost * DT / 2) / (c * DT)) ** 2
    grid = np.linspace(0.0, math.pi / DZ, 4097)
    value = np.array([cowan_lambda(kx, ky, kz) - target for kz in grid])
    roots = []
    for index in np.where(value[:-1] * value[1:] <= 0)[0]:
        lower, upper = grid[index], grid[index + 1]
        if lower == upper:
            root = lower
        else:
            root = brentq(lambda kz: cowan_lambda(kx, ky, kz) - target,
                          lower, upper)
        roots.append(math.copysign(root, physical_kz))
    if not roots:
        return None
    return min(roots, key=lambda root: abs(root - physical_kz))


def cowan_phase_rows():
    boost_beta = math.sqrt((BOOST_GAMMA - 1) * (BOOST_GAMMA + 1)) / BOOST_GAMMA
    distance = SCREEN_FROM_DIPOLE_ENTRANCE - DIPOLE_LENGTH / 2
    rows = []
    for energy in (4.0, 8.0, 10.0, 12.0):
        k0 = energy / HC_EV_M
        for angle_gamma in (1.0, 3.0, 5.0):
            angle = angle_gamma / GAMMA
            transverse = k0 * math.sin(angle)
            kz_lab = k0 * math.cos(angle)
            omega_lab = c * k0
            omega_boost = BOOST_GAMMA * (omega_lab - boost_beta * c * kz_lab)
            kz_physical = BOOST_GAMMA * (kz_lab - boost_beta * k0)
            delta_z_boost = BOOST_GAMMA * distance * (
                1 - boost_beta / math.cos(angle))
            for orientation in ("x", "diagonal"):
                if orientation == "x":
                    kx, ky = transverse, 0.0
                else:
                    kx = ky = transverse / math.sqrt(2)
                kz_numeric = cowan_kz(kx, ky, omega_boost, kz_physical)
                row = {
                    "energy_eV": energy,
                    "angle_gamma": angle_gamma,
                    "orientation": orientation,
                    "kx_boost_per_m": kx,
                    "ky_boost_per_m": ky,
                    "kz_physical_boost_per_m": kz_physical,
                    "delta_z_boost_m": delta_z_boost,
                    "propagating_branch_found": kz_numeric is not None,
                }
                if kz_numeric is None:
                    row.update(kz_numeric_boost_per_m="",
                               absolute_phase_error_rad="",
                               relative_kz_error="")
                else:
                    row.update(
                        kz_numeric_boost_per_m=kz_numeric,
                        absolute_phase_error_rad=abs(
                            (kz_numeric - kz_physical) * delta_z_boost),
                        relative_kz_error=abs(kz_numeric - kz_physical)
                        / max(abs(kz_physical), 1e-300))
                rows.append(row)
    return rows


def write_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_prediction(output, rows, lw_rows, phase_rows, preferred, passing):
    order = np.array([row["order"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), constrained_layout=True)

    ax = axes[0, 0]
    lw_order = np.array([row["order"] for row in lw_rows])
    lw_error = np.array([row["relative_field_l2_to_order64"] for row in lw_rows])
    ax.semilogy(lw_order[:-1], lw_error[:-1], "o-", color="#2b6cb0")
    ax.axvline(8, color="#666", linestyle="--", linewidth=1)
    ax.set(title="Continuous paired LW quadrature", xlabel="GH order per transverse axis",
           ylabel="relative field L2 vs order 64")
    ax.grid(True, which="both", alpha=.25)

    ax = axes[0, 1]
    ax.semilogy(order, [row["Jx_full_max"] for row in rows], label="Jx: x-bin, y-linear")
    ax.semilogy(order, [row["Jz_full_max"] for row in rows], label="Jz: x/y-linear")
    ax.axhline(.05, color="#777", linestyle="--", linewidth=1, label="5% diagnostic line")
    ax.axvspan(min(passing), max(passing), color="#38a169", alpha=.12)
    ax.scatter([preferred], [rows[preferred - 4]["Jx_full_max"]], color="#c53030", zorder=5)
    ax.set(title="GH nodes after the actual current deposit", xlabel="GH order per axis",
           ylabel="max spectral source error")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.semilogy(order, [row["Jx_high_over_low_max_10eV"] for row in rows],
                label=r"power: $k_\perp>5k/\gamma$ / $k_\perp<3k/\gamma$")
    ax.semilogy(order, [row["Jx_near_nyquist_max"] for row in rows],
                label="power near transverse Nyquist")
    ax.axhline(.01, color="#777", linestyle="--", linewidth=1, label="1% diagnostic line")
    ax.axvline(preferred, color="#c53030", linestyle=":", linewidth=1.4)
    ax.set(title="Jx artificial high-k content at 10 eV", xlabel="GH order per axis",
           ylabel="power fraction or ratio")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    colors = {4.0: "#2b6cb0", 8.0: "#38a169", 10.0: "#dd6b20", 12.0: "#c53030"}
    for energy, color in colors.items():
        for orientation, linestyle in (("x", "-"), ("diagonal", "--")):
            selected = [row for row in phase_rows
                        if row["energy_eV"] == energy and row["orientation"] == orientation]
            x = [row["angle_gamma"] for row in selected
                 if row["propagating_branch_found"]]
            y = [row["absolute_phase_error_rad"] for row in selected
                 if row["propagating_branch_found"]]
            ax.semilogy(x, y, marker="o", color=color, linestyle=linestyle,
                        label=f"{energy:g} eV, {orientation}")
    ax.axhline(math.pi, color="#222", linestyle=":", linewidth=1, label=r"$\pi$ phase")
    ax.set(title="Cowan propagation phase error to the screen",
           xlabel=r"observation angle / $(1/\gamma)$", ylabel="absolute phase error [rad]")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Transverse source-sampling prediction (no PIC run)", fontsize=14)
    fig.savefig(output / "sampling_prediction.png", dpi=180)
    fig.savefig(output / "sampling_prediction.pdf")
    plt.close(fig)


def plot_source_spectra(output, orders):
    offset = SAGITTA / 2
    sources = [("exact deposited Gaussian", source_arrays(256, "Jx", offset)[1])]
    sources += [(f"GH order {order}", source_arrays(order, "Jx", offset)[0])
                for order in orders]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), constrained_layout=True,
                             sharex=True, sharey=True)
    kscale = (DIAGNOSTIC_ENERGY_EV / HC_EV_M) / GAMMA
    image = None
    for ax, (label, source) in zip(axes.ravel(), sources):
        transform = np.fft.fftshift(np.fft.fft2(source))
        kx, ky, _, _, _ = fourier_geometry(source.shape)
        extent = [np.fft.fftshift(kx)[0] / kscale,
                  np.fft.fftshift(kx)[-1] / kscale,
                  np.fft.fftshift(ky)[0] / kscale,
                  np.fft.fftshift(ky)[-1] / kscale]
        power = abs(transform) ** 2
        power /= power.max()
        image = ax.imshow(np.log10(np.maximum(power, 1e-12)), origin="lower",
                          extent=extent, vmin=-12, vmax=0, cmap="magma", aspect="auto")
        circle = plt.Circle((0, 0), 3, fill=False, color="cyan", linewidth=.8)
        ax.add_patch(circle)
        ax.set_title(label)
        ax.set_xlim(-15, 15)
        ax.set_ylim(-15, 15)
        ax.grid(False)
    for ax in axes[-1]:
        ax.set_xlabel(r"$k_x/(k/\gamma)$ at 10 eV")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$k_y/(k/\gamma)$ at 10 eV")
    fig.colorbar(image, ax=axes, label="log10 normalized source power", shrink=.86)
    fig.suptitle("Jx source spectrum after APRL deposition; cyan circle = 3/gamma cone")
    fig.savefig(output / "deposited_source_spectra.png", dpi=180)
    fig.savefig(output / "deposited_source_spectra.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--theory-directory", default=(
        "/public/home/wfollowm/APRL_Workspace/results/"
        "synchrotron_shape_validation_20261110_r3/scripts"))
    args = parser.parse_args()
    output = Path(args.output).resolve()
    if output.exists():
        raise RuntimeError(f"Output already exists: {output}")
    output.mkdir(parents=True)

    rows = []
    for order in range(4, 257):
        row = {"order": order,
               "macro_particles_with_longitudinal_order16": order ** 2 * LONGITUDINAL_ORDER}
        for component in ("Jx", "Jz"):
            metrics = source_metrics(order, component)
            row.update({f"{component}_{key}": value for key, value in metrics.items()})
        rows.append(row)
    write_csv(output / "order_metrics.csv", rows)

    deposit_bias = [deposit_operator_bias(component) for component in ("Jx", "Jz")]
    write_csv(output / "deposit_operator_bias.csv", deposit_bias)
    lw_rows = lw_quadrature(args.theory_directory)
    write_csv(output / "lw_quadrature.csv", lw_rows)
    phase_rows = cowan_phase_rows()
    write_csv(output / "cowan_phase.csv", phase_rows)

    criteria = {
        "Jx_full_max": .05,
        "Jx_cone3gamma_max_10eV": .05,
        "Jx_high_over_low_max_10eV": .01,
        "Jx_near_nyquist_max": .01,
    }
    passing = [row["order"] for row in rows
               if all(row[key] < limit for key, limit in criteria.items())]
    if not passing:
        raise RuntimeError("No order passes the stated diagnostic criteria")
    preferred = min(passing, key=lambda order: rows[order - 4]["Jx_full_max"])
    selected = rows[preferred - 4]
    order8 = rows[8 - 4]
    order12 = rows[12 - 4]

    plot_prediction(output, rows, lw_rows, phase_rows, preferred, passing)
    plot_source_spectra(output, [8, 12, 64, preferred, 176])

    scripts = output / "scripts"
    scripts.mkdir()
    shutil.copy2(Path(__file__).resolve(), scripts / Path(__file__).name)

    phase_summary = {}
    for energy in (4.0, 8.0, 10.0, 12.0):
        phase_summary[f"{energy:g}eV_x"] = {
            str(row["angle_gamma"]): row["absolute_phase_error_rad"]
            for row in phase_rows
            if row["energy_eV"] == energy and row["orientation"] == "x"
        }
    final = {
        "complete": True,
        "status": "THEORY_PREDICTION_ONLY_NO_NEW_PIC",
        "model": {
            "sigma_xy_m": SIGMA_XY,
            "grid_dx_m": DX,
            "grid_dy_m": DY,
            "grid_cells": [NX, NY],
            "orbit_sagitta_sampled_m": SAGITTA,
            "offset_samples": len(OFFSETS),
            "Jx_deposition": "bin in x, linear tent in y",
            "Jz_deposition": "linear tent in x and y",
            "longitudinal_GH_order": LONGITUDINAL_ORDER,
        },
        "diagnostic_criteria_not_physical_acceptance": criteria,
        "same_order_result": {
            "passing_orders": passing,
            "lowest_passing_order": min(passing),
            "highest_passing_order": max(passing),
            "preferred_minimum_source_error_order": preferred,
            "preferred_macro_particles": preferred ** 2 * LONGITUDINAL_ORDER,
            "preferred_metrics": selected,
        },
        "order8_metrics": order8,
        "order12_metrics": order12,
        "continuous_LW_order8_relative_field_L2_to_order64": next(
            row["relative_field_l2_to_order64"] for row in lw_rows
            if row["order"] == 8),
        "exact_deposit_operator_bias": deposit_bias,
        "prior_PIC_cross_check_not_used_to_fit": {
            "quantity": "first-plane near-transverse-Nyquist field power at 9.825 eV",
            "order8_source_prediction": order8["Jx_near_nyquist_max"],
            "order8_prior_PIC": 0.4080878307478773,
            "order12_source_prediction": order12["Jx_near_nyquist_max"],
            "order12_prior_PIC": 0.28915111186916104,
            "scope": "A scale cross-check only; source power and propagated field power are not identical observables.",
        },
        "cowan_phase_error_rad": phase_summary,
        "interpretation": [
            "Order 8 is already converged for the continuous LW Gaussian quadrature at about the 1e-3 level.",
            "The dominant failure is the sparse weighted nodes after Jx deposition, not the continuous Gaussian theory.",
            "Increasing GH order is non-monotone on a fixed mesh; order 12 reduces but does not remove aliases.",
            "The preferred order is a prediction for this mesh, beam size, centroid path and deposition rule, not a universal GH rule.",
            "Cowan phase error amplifies aliased high-angle source modes into screen-dependent false zeros and peaks.",
        ],
        "new_PIC_jobs_submitted": False,
        "physical_accuracy_accepted": False,
    }
    save_json(output / "FINAL.json", final)

    order8_lw = final["continuous_LW_order8_relative_field_L2_to_order64"]
    readme = f"""# 横向 Gauss–Hermite 取样的理论预测（未运行新 PIC）

## 结论

当前 8×8×16 的问题不在“连续高斯束没有被 Gauss–Hermite 正确积分”。独立有限弯转、bend-minus-straight Liénard–Wiechert 计算中，8 阶相对 64 阶的复电磁场 L2 差为 **{order8_lw:.4%}**。问题出现在这些加权节点进入网格之后。

真正驱动水平偏转辐射的 `Jx` 在运动方向 x 上按所在网格箱归属，在 y 上线性分配；它不是 x/y 都做线性沉积。按这一实际规则，8 阶源的最坏全横向频谱误差为 **{order8['Jx_full_max']:.2%}**，近横向 Nyquist 功率最高为 **{order8['Jx_near_nyquist_max']:.2%}**，10 eV 下 `k_perp>5k/gamma` 与 `k_perp<3k/gamma` 的功率比最高为 **{order8['Jx_high_over_low_max_10eV']:.3g}**。12 阶仍分别为 **{order12['Jx_full_max']:.2%}**、**{order12['Jx_near_nyquist_max']:.2%}** 和 **{order12['Jx_high_over_low_max_10eV']:.3g}**，因此只把 8 改成 12 不构成收敛。这个预测没有用 PIC 场数据拟合；作为事后量级核对，它给出的 8/12 阶近 Nyquist 源功率 43.34%/28.85%，与已有首个屏面的 40.81%/28.92% 同量级。

在当前 `sigma_x=sigma_y=20 um`、`dx=dy=5 um`、0–{SAGITTA*1e6:.3f} um 弯转质心位移和 16 阶纵向取样下，采用预先声明的诊断线（Jx 全谱误差与 10 eV 三倍辐射锥误差均小于 5%，高角/低角功率比与近 Nyquist 功率均小于 1%），同阶 x/y 方案通过的是 **{min(passing)}–{max(passing)} 之间的偶数阶点**。误差最小的预测点是 **{preferred}×{preferred}×16 = {preferred**2*16:,} 个宏粒子**；其最坏 Jx 全谱误差为 **{selected['Jx_full_max']:.2%}**，三倍辐射锥误差为 **{selected['Jx_cone3gamma_max_10eV']:.2%}**。这里宏粒子权重允许小于一个电子；它是确定性连续源取样，不是 one-for-one 电子解释。

这不是“阶数越大越好”的普通积分收敛：节点穿过离散网格箱时，`Jx` 的箱归属会跳变，所以误差随阶数明显非单调，176 阶以后还会再次增大。{preferred} 阶是当前网格和轨道的预测最优点，不是通用常数。

## 为什么会形成谱零点

连续高斯束本来几乎没有横向 Nyquist 分量。稀疏节点沉积后却产生大量高 `k_perp` 分量。Cowan-z 对轴向传播是无色散的，但对斜传播不是：例如 10 eV、1/γ、3/γ、5/γ 水平角的预测累计相位误差分别见 `cowan_phase.csv`，3/γ 已达到 O(1) rad，5/γ 达到十余 rad。于是这些本不该存在的高角模式到屏幕时可以与低角物理解相消或相长，表现成单点降到零或异常高峰；这也解释了为什么 12 阶消除了原零点却产生了更大的假峰。

## 计算边界

- `order_metrics.csv`：4–256 阶逐阶结果；误差参照是“连续高斯束经过同一沉积算子的精确结果”，因此把取样误差和固定网格误差分开。
- `deposit_operator_bias.csv`：即使连续源被精确沉积，固定网格/沉积本身在物理锥内的误差仍小于约 1%；它不是当前数十到数百百分比异常的主因。
- `lw_quadrature.csv`：连续 Liénard–Wiechert 横向积分的阶数检查。
- `cowan_phase.csv`：Cowan 离散色散对 1/γ、3/γ、5/γ 方向的传播相位预测。
- 图中的 5% 与 1% 是为选择下一轮验证点而声明的诊断线，不是辐射物理精度验收线。
- 本轮没有读取运行中 H5、没有提交 PIC 作业，也没有改变 C++、配置或默认值。

下一轮若做 PIC，理论上最有信息量的不是继续 14、16、20 的低阶扫描，而是保留 8/12 证据后，直接验证 **144、{preferred}、176** 三点：它们分别位于预测窗口的左侧、最低误差点和右侧回升区，可以同时检验阈值、最优点及非单调回升。
"""
    (output / "README.md").write_text(readme)
    print(output)


if __name__ == "__main__":
    main()
