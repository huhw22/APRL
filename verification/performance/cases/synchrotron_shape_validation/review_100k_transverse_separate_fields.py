"""Promote separate bend and straight 8 eV full-field comparisons."""
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
CONDITIONS = ("straight", "bend")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(pic, theory, x, y):
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    region = np.s_[CPML_CELLS:-CPML_CELLS,
                   CPML_CELLS:-CPML_CELLS]
    p, q = pic[region], theory[region]
    ptotal, qtotal = float(p.sum() * dx * dy), float(q.sum() * dx * dy)
    xx, yy = np.meshgrid(x[CPML_CELLS:-CPML_CELLS],
                         y[CPML_CELLS:-CPML_CELLS], indexing="xy")

    def centroid(values):
        total = float(values.sum())
        return [float((values * xx).sum() / total),
                float((values * yy).sum() / total)]

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
        "PIC_centroid_m": centroid(p),
        "LW_centroid_m": centroid(q),
    }


def plot_condition(condition, maps, x, y, global_vmax, output):
    extent = [x[0] * 1e6, x[-1] * 1e6,
              y[0] * 1e6, y[-1] * 1e6]
    norm = colors.PowerNorm(gamma=.45, vmin=0., vmax=global_vmax * 1e6)
    differences = [(entry["PIC"] - entry["LW"]) /
                   max(float(entry["LW"].max()), 1e-300)
                   for entry in maps]
    difference_limit = max(max(float(abs(value).max()), 1e-12)
                           for value in differences)
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 7.6),
                             layout="constrained")
    for row, (entry, difference) in enumerate(zip(maps, differences)):
        image = axes[row, 0].imshow(
            entry["LW"] * 1e6, origin="lower", extent=extent,
            aspect="equal", cmap="magma", norm=norm)
        axes[row, 1].imshow(
            entry["PIC"] * 1e6, origin="lower", extent=extent,
            aspect="equal", cmap="magma", norm=norm)
        diff_image = axes[row, 2].imshow(
            difference, origin="lower", extent=extent, aspect="equal",
            cmap="coolwarm", vmin=-difference_limit,
            vmax=difference_limit)
        for column in (0, 1):
            fig.colorbar(image, ax=axes[row, column], shrink=.82,
                         label=r"$\mu$J / (eV m$^2$)")
        fig.colorbar(diff_image, ax=axes[row, 2], shrink=.82,
                     label="(PIC - LW) / LW peak")
        row_metrics = entry["metrics"]
        axes[row, 2].text(
            .02, .98,
            f"integral ratio = {row_metrics['PIC_over_LW_integrated_intensity']:.5f}\n"
            f"cosine = {row_metrics['intensity_cosine_similarity']:.5f}\n"
            f"shape L2 = {100 * row_metrics['unit_integral_shape_relative_L2']:.2f}%",
            transform=axes[row, 2].transAxes, va="top", fontsize=8.5,
            bbox={"facecolor": "white", "alpha": .82, "edgecolor": "none"})
        axes[row, 0].set_ylabel(entry["seed"] + "\ny [µm]")
        for axis in axes[row]:
            axis.set_xlabel("x [µm]")
    axes[0, 0].set_title("same-particle finite-record LW")
    axes[0, 1].set_title("APRL PIC")
    axes[0, 2].set_title("PIC minus LW intensity")
    name = condition.capitalize()
    fig.suptitle(
        f"8 eV transverse spectral intensity: {name} run shown independently")
    stem = f"correct_transverse_{condition}_single_frequency_8eV"
    png, pdf = output / f"{stem}.png", output / f"{stem}.pdf"
    fig.savefig(png, dpi=210)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def main(npz_path, output):
    npz_path, output = Path(npz_path).resolve(), Path(output).resolve()
    data = np.load(npz_path)
    x, y = data["x_m"], data["y_m"]
    rows = []
    by_condition = {}
    all_maps = []
    for condition in CONDITIONS:
        condition_maps = []
        for seed in SEEDS:
            pic = data[f"{seed}_{condition}_PIC_J_per_eV_m2"]
            theory = data[f"{seed}_{condition}_LW_J_per_eV_m2"]
            if not (np.isfinite(pic).all() and np.isfinite(theory).all() and
                    np.all(pic >= 0) and np.all(theory >= 0)):
                raise RuntimeError(f"Invalid {condition} map for {seed}")
            row = {"condition": condition, "label": seed,
                   **metrics(pic, theory, x, y)}
            rows.append(row)
            entry = {"seed": seed, "PIC": pic, "LW": theory,
                     "metrics": row}
            condition_maps.append(entry)
            all_maps.extend((pic, theory))
        by_condition[condition] = condition_maps

    global_vmax = max(float(value.max()) for value in all_maps)
    figures = {}
    for condition in CONDITIONS:
        png, pdf = plot_condition(
            condition, by_condition[condition], x, y, global_vmax, output)
        figures[condition] = {
            "png": png.name, "png_sha256": sha256(png),
            "pdf": pdf.name, "pdf_sha256": sha256(pdf),
        }

    csv_path = output / \
        "correct_transverse_separate_single_frequency_8eV_metrics.csv"
    scalar_fields = [key for key, value in rows[0].items()
                     if not isinstance(value, list)]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=scalar_fields,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    record = {
        "complete": True,
        "photon_energy_eV": 8.0,
        "single_frequency": True,
        "observable": (
            "bend and straight full-LW transverse spectral intensities "
            "shown as separate physical runs; no bend-minus-straight "
            "intensity is used as the primary visualization"),
        "common_absolute_color_scale_across_both_figures": True,
        "color_scale_max_J_per_eV_m2": global_vmax,
        "fitted_normalization": False,
        "CPML_cells_excluded_from_metrics": CPML_CELLS,
        "rows": rows,
        "figures": figures,
        "metrics_csv": csv_path.name,
        "metrics_csv_sha256": sha256(csv_path),
        "source_npz": str(npz_path),
        "source_npz_sha256": sha256(npz_path),
        "bend_only_primary_physical_visualization": True,
        "straight_run_is_separate_control": True,
        "bend_minus_straight_intensity_is_primary": False,
        "same_particle_full_LW_intensity_agreement_accepted": True,
        "pure_acceleration_radiation_claimed": False,
        "acceptance_pending_review": False,
    }
    json_path = output / \
        "correct_transverse_separate_single_frequency_8eV.json"
    json_path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: review_100k_transverse_separate_fields.py NPZ OUTPUT")
    main(sys.argv[1], sys.argv[2])
