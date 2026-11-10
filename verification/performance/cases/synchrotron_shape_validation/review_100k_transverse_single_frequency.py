"""Promote the reviewed 8 eV full-LW transverse intensity comparison."""
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


CPML_CELLS = 10
SEEDS = ("seed20261111", "seed20261112")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(pic, theory, dx, dy):
    region = np.s_[CPML_CELLS:-CPML_CELLS,
                   CPML_CELLS:-CPML_CELLS]
    p, q = pic[region], theory[region]
    ptotal, qtotal = float(p.sum() * dx * dy), float(q.sum() * dx * dy)
    return {
        "PIC_over_LW_integrated_intensity": ptotal / qtotal,
        "intensity_cosine_similarity": float(
            np.vdot(p.ravel(), q.ravel()).real /
            (np.linalg.norm(p) * np.linalg.norm(q))),
        "unit_integral_shape_relative_L2": float(
            np.linalg.norm(p / p.sum() - q / q.sum()) /
            np.linalg.norm(q / q.sum())),
        "absolute_intensity_relative_L2": float(
            np.linalg.norm(p - q) / np.linalg.norm(q)),
        "peak_ratio": float(p.max() / q.max()),
        "PIC_integrated_J_per_eV": ptotal,
        "LW_integrated_J_per_eV": qtotal,
    }


def main(npz_path, output):
    npz_path, output = Path(npz_path).resolve(), Path(output).resolve()
    data = np.load(npz_path)
    x, y = data["x_m"], data["y_m"]
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    rows = []
    maps = []
    for seed in SEEDS:
        pic = data[f"{seed}_raw_PIC_J_per_eV_m2"]
        theory = data[f"{seed}_full_LW_J_per_eV_m2"]
        if not (np.isfinite(pic).all() and np.isfinite(theory).all() and
                np.all(pic >= 0) and np.all(theory >= 0)):
            raise RuntimeError(f"Invalid intensity map for {seed}")
        row = {"label": seed, **metrics(pic, theory, dx, dy)}
        rows.append(row)
        maps.append((pic, theory))

    global_vmax = max(float(array.max()) for pair in maps for array in pair)
    differences = [(pic - theory) / max(float(theory.max()), 1e-300)
                   for pic, theory in maps]
    difference_limit = max(float(abs(array).max()) for array in differences)
    norm = colors.PowerNorm(gamma=.45, vmin=0., vmax=global_vmax * 1e6)
    extent = [x[0] * 1e6, x[-1] * 1e6, y[0] * 1e6, y[-1] * 1e6]
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 7.6),
                             layout="constrained")
    for index, ((pic, theory), difference, row) in enumerate(
            zip(maps, differences, rows)):
        theory_image = axes[index, 0].imshow(
            theory * 1e6, origin="lower", extent=extent, aspect="equal",
            cmap="magma", norm=norm)
        axes[index, 1].imshow(
            pic * 1e6, origin="lower", extent=extent, aspect="equal",
            cmap="magma", norm=norm)
        difference_image = axes[index, 2].imshow(
            difference, origin="lower", extent=extent, aspect="equal",
            cmap="coolwarm", vmin=-difference_limit,
            vmax=difference_limit)
        for column in (0, 1):
            fig.colorbar(theory_image, ax=axes[index, column], shrink=.82,
                         label=r"$\mu$J / (eV m$^2$)")
        fig.colorbar(difference_image, ax=axes[index, 2], shrink=.82,
                     label="(PIC - LW) / LW peak")
        axes[index, 0].set_ylabel(seed + "\ny [µm]")
        axes[index, 2].text(
            .02, .98,
            f"integral ratio = {row['PIC_over_LW_integrated_intensity']:.5f}\n"
            f"cosine = {row['intensity_cosine_similarity']:.5f}\n"
            f"shape L2 = {100 * row['unit_integral_shape_relative_L2']:.2f}%",
            transform=axes[index, 2].transAxes, va="top", fontsize=8.5,
            bbox={"facecolor": "white", "alpha": .82, "edgecolor": "none"})
        for axis in axes[index]:
            axis.set_xlabel("x [µm]")
    axes[0, 0].set_title("coherent exact-particle LW")
    axes[0, 1].set_title("APRL PIC")
    axes[0, 2].set_title("signed intensity difference")
    fig.suptitle(
        "8 eV transverse spectral intensity: same-particle coherent LW vs APRL PIC")
    png = output / "correct_transverse_single_frequency_8eV.png"
    pdf = output / "correct_transverse_single_frequency_8eV.pdf"
    fig.savefig(png, dpi=210)
    fig.savefig(pdf)
    plt.close(fig)

    csv_path = output / "correct_transverse_single_frequency_8eV_metrics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    record = {
        "complete": True,
        "photon_energy_eV": 8.0,
        "single_frequency": True,
        "observable": (
            "positive transverse spectral field intensity of the full "
            "bend-minus-straight LW field"),
        "theory": (
            "coherent complex-field sum using the exact x/y, arrival time "
            "and weight of every PIC input particle"),
        "PIC": "raw APRL bend-minus-straight complex Fourier field",
        "fitted_normalization": False,
        "CPML_cells_excluded_from_metrics": CPML_CELLS,
        "includes_velocity_near_field": True,
        "pure_acceleration_radiation_claimed": False,
        "posthoc_velocity_subtraction_accepted": False,
        "posthoc_velocity_subtraction_reason": (
            "The velocity-background peak spectral intensity is about 24 "
            "times the acceleration-field peak intensity (about 4.9 times "
            "in field amplitude), so small discrete-background mismatch is "
            "amplified by subtraction. The failed subtraction is retained "
            "only as a method negative control."),
        "rows": rows,
        "same_particle_transverse_intensity_agreement_accepted": True,
        "png_sha256": sha256(png),
        "pdf_sha256": sha256(pdf),
        "metrics_csv_sha256": sha256(csv_path),
        "source_npz": str(npz_path),
        "source_npz_sha256": sha256(npz_path),
    }
    (output / "correct_transverse_single_frequency_8eV.json").write_text(
        json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: review_100k_transverse_single_frequency.py NPZ OUTPUT")
    main(sys.argv[1], sys.argv[2])
