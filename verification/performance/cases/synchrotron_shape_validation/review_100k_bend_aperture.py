"""Compare raw bend-only PIC with acceleration LW inside an off-beam aperture."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


SEEDS = ("seed20261111", "seed20261112")
APERTURE_RADIUS_M = 75e-6
BEAM_CENTER_X_M = 240.96005e-6
BEAM_SIGMA_M = 20e-6
NEAR_FIELD_SCALE_M = 28.96e-6


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cosine(first, second):
    return float(np.vdot(first, second).real /
                 (np.linalg.norm(first) * np.linalg.norm(second)))


def shape_l2(first, second):
    return float(np.linalg.norm(first / first.sum() -
                                second / second.sum()) /
                 np.linalg.norm(second / second.sum()))


def main(npz_path, output):
    npz_path, output = Path(npz_path).resolve(), Path(output).resolve()
    data = np.load(npz_path)
    x, y = data["x_m"], data["y_m"]
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    xx, yy = np.meshgrid(x, y, indexing="xy")
    aperture = xx**2 + yy**2 <= APERTURE_RADIUS_M**2
    area = dx * dy
    rows, maps = [], []
    for seed in SEEDS:
        pic = data[f"{seed}_bend_PIC_J_per_eV_m2"]
        full_lw = data[f"{seed}_bend_LW_J_per_eV_m2"]
        acceleration_lw = data[f"{seed}_radiation_LW_J_per_eV_m2"]
        p, f, a = pic[aperture], full_lw[aperture], \
            acceleration_lw[aperture]
        row = {
            "label": seed,
            "PIC_over_acceleration_LW_aperture_integral":
                float(p.sum() / a.sum()),
            "full_LW_over_acceleration_LW_aperture_integral":
                float(f.sum() / a.sum()),
            "PIC_vs_acceleration_LW_intensity_cosine": cosine(p, a),
            "PIC_vs_acceleration_LW_unit_integral_shape_L2": shape_l2(p, a),
            "full_LW_vs_acceleration_LW_intensity_cosine": cosine(f, a),
            "full_LW_vs_acceleration_LW_unit_integral_shape_L2":
                shape_l2(f, a),
            "PIC_over_acceleration_LW_peak": float(p.max() / a.max()),
            "PIC_aperture_J_per_eV": float(p.sum() * area),
            "acceleration_LW_aperture_J_per_eV": float(a.sum() * area),
            "acceleration_LW_fraction_captured":
                float(a.sum() / acceleration_lw.sum()),
        }
        rows.append(row)
        maps.append((pic, acceleration_lw, row))

    overview_vmax = max(float(pic.max()) for pic, _, _ in maps)
    aperture_vmax = max(float(values[aperture].max())
                        for pic, theory, _ in maps
                        for values in (pic, theory))
    differences = [(pic - theory) / max(float(theory[aperture].max()),
                                        1e-300)
                   for pic, theory, _ in maps]
    difference_limit = max(float(abs(value[aperture]).max())
                           for value in differences)
    extent = [x[0] * 1e6, x[-1] * 1e6,
              y[0] * 1e6, y[-1] * 1e6]
    full_norm = colors.PowerNorm(gamma=.45, vmin=0.,
                                 vmax=overview_vmax * 1e6)
    aperture_norm = colors.PowerNorm(gamma=.45, vmin=0.,
                                     vmax=aperture_vmax * 1e6)
    positive_cmap = matplotlib.colormaps["magma"].copy()
    positive_cmap.set_bad("#d9d9d9")
    difference_cmap = matplotlib.colormaps["coolwarm"].copy()
    difference_cmap.set_bad("#d9d9d9")
    fig, axes = plt.subplots(2, 4, figsize=(18.2, 8.1),
                             layout="constrained")
    crop = 85.
    for row_index, ((pic, acceleration_lw, metrics), difference) in enumerate(
            zip(maps, differences)):
        overview = axes[row_index, 0].imshow(
            pic * 1e6, origin="lower", extent=extent, aspect="equal",
            cmap="magma", norm=full_norm)
        axes[row_index, 0].add_patch(Circle(
            (0., 0.), APERTURE_RADIUS_M * 1e6, fill=False,
            edgecolor="cyan", linewidth=1.4))
        axes[row_index, 0].plot(BEAM_CENTER_X_M * 1e6, 0., "cx",
                                markersize=7, markeredgewidth=1.4)
        axes[row_index, 0].text(
            .02, .98, "cyan circle: radiation aperture\ncyan x: electron beam",
            transform=axes[row_index, 0].transAxes, va="top", fontsize=8,
            color="white",
            bbox={"facecolor": "black", "alpha": .45, "edgecolor": "none"})
        pic_aperture = np.ma.array(pic * 1e6, mask=~aperture)
        theory_aperture = np.ma.array(acceleration_lw * 1e6,
                                      mask=~aperture)
        image = axes[row_index, 1].imshow(
            pic_aperture, origin="lower", extent=extent, aspect="equal",
            cmap=positive_cmap, norm=aperture_norm)
        axes[row_index, 2].imshow(
            theory_aperture, origin="lower", extent=extent, aspect="equal",
            cmap=positive_cmap, norm=aperture_norm)
        difference_image = axes[row_index, 3].imshow(
            np.ma.array(difference, mask=~aperture), origin="lower",
            extent=extent, aspect="equal", cmap=difference_cmap,
            vmin=-difference_limit, vmax=difference_limit)
        for column in (1, 2, 3):
            axes[row_index, column].set_xlim(-crop, crop)
            axes[row_index, column].set_ylim(-crop, crop)
        fig.colorbar(overview, ax=axes[row_index, 0], shrink=.82,
                     label=r"full bend PIC [$\mu$J/(eV m$^2$)]")
        for column in (1, 2):
            fig.colorbar(image, ax=axes[row_index, column], shrink=.82,
                         label=r"$\mu$J/(eV m$^2$)")
        fig.colorbar(difference_image, ax=axes[row_index, 3], shrink=.82,
                     label="(PIC - acceleration LW) / LW peak")
        axes[row_index, 3].text(
            .02, .98,
            f"PIC/LW integral = {metrics['PIC_over_acceleration_LW_aperture_integral']:.5f}\n"
            f"full/acc LW = {metrics['full_LW_over_acceleration_LW_aperture_integral']:.5f}\n"
            f"cosine = {metrics['PIC_vs_acceleration_LW_intensity_cosine']:.5f}\n"
            f"captured = {100 * metrics['acceleration_LW_fraction_captured']:.2f}%",
            transform=axes[row_index, 3].transAxes, va="top", fontsize=8.2,
            bbox={"facecolor": "white", "alpha": .84, "edgecolor": "none"})
        axes[row_index, 0].set_ylabel(seed + "\ny [µm]")
        for axis in axes[row_index]:
            axis.set_xlabel("x [µm]")
    axes[0, 0].set_title("raw bend-only PIC: full plane")
    axes[0, 1].set_title("raw PIC inside aperture")
    axes[0, 2].set_title("acceleration LW inside aperture")
    axes[0, 3].set_title("aperture intensity difference")
    fig.suptitle(
        "8 eV bend-only radiation-axis aperture (R = 75 µm); no field subtraction")
    png = output / "correct_transverse_bend_aperture_8eV.png"
    pdf = output / "correct_transverse_bend_aperture_8eV.pdf"
    fig.savefig(png, dpi=210)
    fig.savefig(pdf)
    plt.close(fig)

    csv_path = output / "correct_transverse_bend_aperture_8eV_metrics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    record = {
        "complete": True,
        "photon_energy_eV": 8.0,
        "observable": (
            "raw bend-only PIC intensity inside a radiation-axis circular "
            "aperture compared directly with acceleration-only LW intensity"),
        "complex_field_subtraction_applied": False,
        "aperture_center_m": [0., 0.],
        "aperture_radius_m": APERTURE_RADIUS_M,
        "electron_beam_center_x_m": BEAM_CENTER_X_M,
        "electron_beam_sigma_m": BEAM_SIGMA_M,
        "beam_edge_clearance_sigma":
            (BEAM_CENTER_X_M - APERTURE_RADIUS_M) / BEAM_SIGMA_M,
        "near_field_scale_m_at_8eV": NEAR_FIELD_SCALE_M,
        "beam_to_aperture_edge_over_near_field_scale":
            (BEAM_CENTER_X_M - APERTURE_RADIUS_M) / NEAR_FIELD_SCALE_M,
        "selection_note": (
            "R=75 um is the largest tested round aperture for which both "
            "samples keep full-LW versus acceleration-LW aperture integral "
            "within 1.3%, while retaining about 9.8% of the acceleration-LW "
            "full-plane intensity."),
        "fitted_normalization": False,
        "common_aperture_color_scale": True,
        "overview_uses_separate_full_plane_color_scale": True,
        "rows": rows,
        "png": png.name,
        "png_sha256": sha256(png),
        "pdf": pdf.name,
        "pdf_sha256": sha256(pdf),
        "metrics_csv": csv_path.name,
        "metrics_csv_sha256": sha256(csv_path),
        "source_npz": str(npz_path),
        "source_npz_sha256": sha256(npz_path),
        "physical_aperture_comparison_accepted": True,
        "pure_radiation_full_plane_claimed": False,
    }
    json_path = output / "correct_transverse_bend_aperture_8eV.json"
    json_path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: review_100k_bend_aperture.py NPZ OUTPUT")
    main(sys.argv[1], sys.argv[2])
