"""Scan angular acceptance and project bend-only fields onto forward vacuum modes."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar, mu_0


ENERGY_EV = 8.0
SEEDS = ("seed20261111", "seed20261112")
THETA_MRAD = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0, 16.0)
DISPLAY_THETA_MRAD = 1.0
TAPER_CELLS = 10


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def spectral_intensity(ef, bf):
    transverse = (abs(ef[:, :, 0])**2 + abs(ef[:, :, 1])**2 +
                  c**2 * (abs(bf[:, :, 0])**2 + abs(bf[:, :, 1])**2))
    return e / (h * mu_0 * c) * transverse


def edge_taper(size, cells):
    result = np.ones(size)
    phase = (np.arange(cells) + .5) / cells
    ramp = np.sin(.5 * np.pi * phase)**2
    result[:cells] = ramp
    result[-cells:] = ramp[::-1]
    return result


def mode_geometry(x, y):
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    kx = 2 * np.pi * np.fft.fftfreq(len(x), d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(len(y), d=dy)
    kxx, kyy = np.meshgrid(kx, ky, indexing="xy")
    k0 = ENERGY_EV * e / (hbar * c)
    kperp2 = kxx**2 + kyy**2
    propagating = kperp2 < k0**2
    kz = np.sqrt(np.maximum(k0**2 - kperp2, 0.))
    khat = np.stack((kxx / k0, kyy / k0, kz / k0), axis=-1)
    theta = np.arcsin(np.minimum(np.sqrt(kperp2) / k0, 1.))
    taper = edge_taper(len(y), TAPER_CELLS)[:, None] * \
        edge_taper(len(x), TAPER_CELLS)[None, :]
    return dx, dy, khat, theta, propagating, taper


def spatial_modes(ef, bf_stored, dt, taper):
    omega = ENERGY_EV * e / hbar
    bf = bf_stored * np.exp(-.5j * omega * dt)
    e_k = np.fft.fft2(ef * taper[:, :, None], axes=(0, 1), norm="ortho")
    b_k = np.fft.fft2(bf * taper[:, :, None], axes=(0, 1), norm="ortho")
    return e_k, b_k


def forward_vacuum_modes(e_k, b_k, khat, propagating):
    e_t = e_k - khat * np.sum(khat * e_k, axis=-1)[:, :, None]
    b_t = b_k - khat * np.sum(khat * b_k, axis=-1)[:, :, None]
    e_forward = .5 * (e_t - c * np.cross(khat, b_t))
    e_forward[~propagating] = 0.
    return e_forward


def reconstruct(e_modes, khat, mask):
    selected = e_modes * mask[:, :, None]
    b_selected = np.cross(khat, selected) / c
    ef = np.fft.ifft2(selected, axes=(0, 1), norm="ortho")
    bf = np.fft.ifft2(b_selected, axes=(0, 1), norm="ortho")
    return ef, bf, spectral_intensity(ef, bf)


def reconstruct_raw(e_modes, b_modes, mask):
    ef = np.fft.ifft2(e_modes * mask[:, :, None],
                      axes=(0, 1), norm="ortho")
    bf = np.fft.ifft2(b_modes * mask[:, :, None],
                      axes=(0, 1), norm="ortho")
    return spectral_intensity(ef, bf)


def projection_self_check():
    direction = np.array([0., 0., 1.])
    electric = np.array([1. + .2j, .3 - .1j, 0.])
    magnetic_forward = np.cross(direction, electric) / c
    magnetic_backward = -magnetic_forward
    projected_forward = .5 * (
        electric - c * np.cross(direction, magnetic_forward))
    projected_backward = .5 * (
        electric - c * np.cross(direction, magnetic_backward))
    forward_error = float(np.max(abs(projected_forward - electric)))
    backward_error = float(np.max(abs(projected_backward)))
    if not (forward_error < 1e-15 and backward_error < 1e-15):
        raise RuntimeError("Forward/backward vacuum-mode projection self-check failed")
    return {"forward_mode_max_abs_error": forward_error,
            "backward_mode_max_abs_leakage": backward_error}


def cosine(first, second):
    return float(np.vdot(first.ravel(), second.ravel()).real /
                 (np.linalg.norm(first) * np.linalg.norm(second)))


def shape_l2(first, second):
    return float(np.linalg.norm(first / first.sum() -
                                second / second.sum()) /
                 np.linalg.norm(second / second.sum()))


def main(npz_path, output):
    npz_path, output = Path(npz_path).resolve(), Path(output).resolve()
    data = np.load(npz_path)
    self_check = projection_self_check()
    x, y = data["x_m"], data["y_m"]
    dt = float(data["dt_s"])
    dx, dy, khat, theta, propagating, taper = mode_geometry(x, y)
    area = dx * dy
    rows, display, projection_rows = [], [], []
    for seed in SEEDS:
        pic_e = data[f"{seed}_bend_PIC_E"]
        pic_b = data[f"{seed}_bend_PIC_B"]
        lw_e = data[f"{seed}_acceleration_LW_E"]
        lw_b = data[f"{seed}_acceleration_LW_B"]
        pic_e_modes, pic_b_modes = spatial_modes(pic_e, pic_b, dt, taper)
        lw_e_modes, lw_b_modes = spatial_modes(lw_e, lw_b, dt, taper)
        pic_modes = forward_vacuum_modes(
            pic_e_modes, pic_b_modes, khat, propagating)
        lw_modes = forward_vacuum_modes(
            lw_e_modes, lw_b_modes, khat, propagating)
        _, _, pic_all = reconstruct(pic_modes, khat, propagating)
        _, _, lw_all = reconstruct(lw_modes, khat, propagating)
        raw_pic = reconstruct_raw(pic_e_modes, pic_b_modes, propagating)
        raw_lw = reconstruct_raw(lw_e_modes, lw_b_modes, propagating)
        shared_metrics = {
            "raw_tapered_PIC_J_per_eV": float(raw_pic.sum() * area),
            "raw_tapered_LW_J_per_eV": float(raw_lw.sum() * area),
            "all_forward_PIC_fraction_of_raw": float(
                pic_all.sum() / raw_pic.sum()),
            "all_forward_LW_fraction_of_raw": float(
                lw_all.sum() / raw_lw.sum()),
        }
        projection_row = {
            "label": seed,
            "PIC_over_acceleration_LW": float(pic_all.sum() / lw_all.sum()),
            "intensity_cosine": cosine(pic_all, lw_all),
            "unit_integral_shape_L2": shape_l2(pic_all, lw_all),
            **shared_metrics,
        }
        projection_rows.append(projection_row)
        display.append({"label": seed, "PIC": pic_all, "LW": lw_all,
                        "metrics": projection_row})
        for theta_mrad in THETA_MRAD:
            mask = propagating & (theta <= theta_mrad * 1e-3)
            _, _, pic_i = reconstruct(pic_modes, khat, mask)
            _, _, lw_i = reconstruct(lw_modes, khat, mask)
            raw_pic_i = reconstruct_raw(pic_e_modes, pic_b_modes, mask)
            raw_lw_i = reconstruct_raw(lw_e_modes, lw_b_modes, mask)
            row = {
                "label": seed,
                "theta_max_mrad": theta_mrad,
                "retained_mode_count": int(mask.sum()),
                "PIC_J_per_eV": float(pic_i.sum() * area),
                "acceleration_LW_J_per_eV": float(lw_i.sum() * area),
                "PIC_over_acceleration_LW": float(pic_i.sum() / lw_i.sum()),
                "intensity_cosine": cosine(pic_i, lw_i),
                "unit_integral_shape_L2": shape_l2(pic_i, lw_i),
                "PIC_fraction_of_all_forward_modes":
                    float(pic_i.sum() / pic_all.sum()),
                "LW_fraction_of_all_forward_modes":
                    float(lw_i.sum() / lw_all.sum()),
                "raw_angle_PIC_J_per_eV": float(raw_pic_i.sum() * area),
                "raw_angle_acceleration_LW_J_per_eV":
                    float(raw_lw_i.sum() * area),
                "raw_angle_PIC_over_acceleration_LW":
                    float(raw_pic_i.sum() / raw_lw_i.sum()),
                "raw_angle_intensity_cosine": cosine(raw_pic_i, raw_lw_i),
                "raw_angle_unit_integral_shape_L2":
                    shape_l2(raw_pic_i, raw_lw_i),
                "raw_angle_PIC_fraction_of_all":
                    float(raw_pic_i.sum() / raw_pic.sum()),
                "raw_angle_LW_fraction_of_all":
                    float(raw_lw_i.sum() / raw_lw.sum()),
                **shared_metrics,
            }
            rows.append(row)
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.5),
                             layout="constrained")
    for seed, marker in zip(SEEDS, ("o", "s")):
        selected = [row for row in rows if row["label"] == seed]
        angle = [row["theta_max_mrad"] for row in selected]
        axes[0].plot(angle,
                     [row["raw_angle_PIC_over_acceleration_LW"]
                      for row in selected],
                     marker=marker, linestyle="--", label=seed + " angle only")
        axes[0].plot(angle,
                     [row["PIC_over_acceleration_LW"] for row in selected],
                     marker=marker, label=seed + " projected")
        axes[1].plot(angle,
                     [row["raw_angle_intensity_cosine"] for row in selected],
                     marker=marker, linestyle="--", label=seed + " angle only")
        axes[1].plot(angle,
                     [row["intensity_cosine"] for row in selected],
                     marker=marker, label=seed + " projected")
        axes[2].plot(angle,
                     [row["raw_angle_LW_fraction_of_all"] for row in selected],
                     marker=marker, linestyle="--", label=seed + " angle only")
        axes[2].plot(angle,
                     [row["LW_fraction_of_all_forward_modes"] for row in selected],
                     marker=marker, label=seed + " projected")
    axes[0].axhline(1., color="black", linewidth=.7)
    axes[0].set_ylabel("PIC / acceleration LW")
    axes[1].set_ylabel("2D intensity cosine")
    axes[2].set_ylabel("retained acceleration-LW fraction")
    for axis in axes:
        axis.set_xscale("log")
        axis.set_xlabel(r"$\theta_{max}$ [mrad]")
        axis.grid(alpha=.22)
        axis.legend(fontsize=7)
    fig.suptitle(
        "8 eV bend-only forward-vacuum-mode angular acceptance scan")
    scan_png = output / "theta_scan_8eV.png"
    scan_pdf = output / "theta_scan_8eV.pdf"
    fig.savefig(scan_png, dpi=210)
    fig.savefig(scan_pdf)
    plt.close(fig)

    vmax = max(float(values[key].max()) for values in display
               for key in ("PIC", "LW"))
    differences = [(values["PIC"] - values["LW"]) /
                   max(float(values["LW"].max()), 1e-300)
                   for values in display]
    difference_limit = max(float(abs(value).max()) for value in differences)
    norm = colors.PowerNorm(gamma=.45, vmin=0., vmax=vmax * 1e6)
    extent = [x[0] * 1e6, x[-1] * 1e6,
              y[0] * 1e6, y[-1] * 1e6]
    fig, axes = plt.subplots(2, 3, figsize=(14.7, 7.7),
                             layout="constrained")
    for index, (values, difference) in enumerate(zip(display, differences)):
        image = axes[index, 0].imshow(
            values["LW"] * 1e6, origin="lower", extent=extent,
            aspect="equal", cmap="magma", norm=norm)
        axes[index, 1].imshow(
            values["PIC"] * 1e6, origin="lower", extent=extent,
            aspect="equal", cmap="magma", norm=norm)
        difference_image = axes[index, 2].imshow(
            difference, origin="lower", extent=extent, aspect="equal",
            cmap="coolwarm", vmin=-difference_limit,
            vmax=difference_limit)
        for column in (0, 1):
            fig.colorbar(image, ax=axes[index, column], shrink=.82,
                         label=r"$\mu$J/(eV m$^2$)")
        fig.colorbar(difference_image, ax=axes[index, 2], shrink=.82,
                     label="(PIC - LW) / LW peak")
        row = values["metrics"]
        axes[index, 2].text(
            .02, .98,
            f"PIC/LW = {row['PIC_over_acceleration_LW']:.5f}\n"
            f"cosine = {row['intensity_cosine']:.5f}\n"
            f"forward/raw PIC = {100 * row['all_forward_PIC_fraction_of_raw']:.2f}%\n"
            f"forward/raw LW = {100 * row['all_forward_LW_fraction_of_raw']:.2f}%",
            transform=axes[index, 2].transAxes, va="top", fontsize=8.5,
            bbox={"facecolor": "white", "alpha": .84, "edgecolor": "none"})
        axes[index, 0].set_ylabel(values["label"] + "\ny [µm]")
        for axis in axes[index]:
            axis.set_xlabel("x [µm]")
    axes[0, 0].set_title("acceleration LW")
    axes[0, 1].set_title("raw PIC projected to forward modes")
    axes[0, 2].set_title("projected intensity difference")
    fig.suptitle(
        "8 eV bend-only projection onto all representable forward vacuum modes")
    modes_png = output / "propagating_mode_projection_8eV.png"
    modes_pdf = output / "propagating_mode_projection_8eV.pdf"
    fig.savefig(modes_png, dpi=210)
    fig.savefig(modes_pdf)
    plt.close(fig)

    csv_path = output / "theta_scan_8eV.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    record = {
        "complete": True,
        "photon_energy_eV": ENERGY_EV,
        "method": (
            "2D complex E/B Fourier decomposition on one detector plane; "
            "Yee half-step B phase corrected; projected onto forward vacuum "
            "Maxwell modes and then limited by theta_max"),
        "not_an_exact_closed_Huygens_transform": True,
        "projection_self_check": self_check,
        "transverse_edge_taper_cells": TAPER_CELLS,
        "taper_amplitude_normalization_applied": False,
        "theta_scan_mrad": list(THETA_MRAD),
        "projection_display_angular_cutoff": None,
        "angular_resolution_mrad": {
            "x": float(2 * np.pi / (len(x) * dx) /
                       (ENERGY_EV * e / (hbar * c)) * 1e3),
            "y": float(2 * np.pi / (len(y) * dy) /
                       (ENERGY_EV * e / (hbar * c)) * 1e3),
        },
        "rows": rows,
        "projection_only_rows": projection_rows,
        "theta_scan_png": scan_png.name,
        "theta_scan_png_sha256": sha256(scan_png),
        "theta_scan_pdf": scan_pdf.name,
        "theta_scan_pdf_sha256": sha256(scan_pdf),
        "projection_png": modes_png.name,
        "projection_png_sha256": sha256(modes_png),
        "projection_pdf": modes_pdf.name,
        "projection_pdf_sha256": sha256(modes_pdf),
        "csv": csv_path.name,
        "csv_sha256": sha256(csv_path),
        "source_npz": str(npz_path),
        "source_npz_sha256": sha256(npz_path),
        "acceptance_pending_review": True,
    }
    (output / "MODE_PROJECTION.json").write_text(
        json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: review_100k_bend_modes.py NPZ OUTPUT")
    main(sys.argv[1], sys.argv[2])
