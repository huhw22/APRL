"""Compare bend fields after particle boundary exit with the no-exit control."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar, mu_0


ENERGY_EV = 8.0
DISPLAY_RADIUS = 80e-6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


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


def load_field(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    with h5py.File(path, "r") as handle:
        group = handle["field_plane"]
        require(int(group["complete"][()]) == 1, f"Incomplete {path}")
        time = group["time_s"][:]
        nx, ny = int(group.attrs["nx"]), int(group.attrs["ny"])
        dx, dy = float(group.attrs["dx_m"]), float(group.attrs["dy_m"])
        x = float(group.attrs["x_first_m"]) + np.arange(nx) * dx
        y = float(group.attrs["y_first_m"]) + np.arange(ny) * dy
        dt = float((time[-1] - time[0]) / (len(time) - 1))
        phase = np.exp(1j * ENERGY_EV * e / hbar * (time - time[0])) * dt
        fields = group["fields"]
        transformed = {}
        for key, name in (("E", "electric_V_per_m"), ("B", "magnetic_T")):
            result = np.zeros((ny, nx, 3), complex)
            for first in range(0, len(time), 64):
                last = min(first + 64, len(time))
                result += np.einsum(
                    "t,tyxd->yxd", phase[first:last],
                    fields[name][first:last], optimize=True)
            transformed[key] = result
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            f"Closed detector changed while reading: {path}")
    return {"x": x, "y": y, "dx": dx, "dy": dy, "dt": dt,
            "intensity": spectral_intensity(transformed["E"], transformed["B"])}


def index_subset(full, subset):
    indices = []
    for value in subset:
        matches = np.flatnonzero(np.isclose(full, value, rtol=0, atol=1e-15))
        require(len(matches) == 1, "New grid is not an exact old-grid subset")
        indices.append(int(matches[0]))
    return np.asarray(indices)


def metrics(pic, theory, control, mask, area):
    p, q, r = pic[mask], theory[mask], control[mask]
    def cosine(first, second):
        return float(np.vdot(first, second).real /
                     (np.linalg.norm(first) * np.linalg.norm(second)))
    def shape_l2(first, second):
        return float(np.linalg.norm(first / first.sum() - second / second.sum()) /
                     np.linalg.norm(second / second.sum()))
    return {
        "boundary_PIC_J_per_eV": float(p.sum() * area),
        "acceleration_LW_J_per_eV": float(q.sum() * area),
        "collision_free_PIC_J_per_eV": float(r.sum() * area),
        "boundary_PIC_over_acceleration_LW": float(p.sum() / q.sum()),
        "boundary_PIC_vs_acceleration_LW_cosine": cosine(p, q),
        "boundary_PIC_vs_acceleration_LW_shape_L2": shape_l2(p, q),
        "boundary_PIC_over_collision_free_PIC": float(p.sum() / r.sum()),
        "boundary_PIC_vs_collision_free_PIC_cosine": cosine(p, r),
        "boundary_PIC_vs_collision_free_PIC_shape_L2": shape_l2(p, r),
    }


def main(root):
    root = Path(root).resolve()
    output = root / "analysis"
    require(not output.exists(), "Analysis output already exists")
    output.mkdir()
    campaign = json.loads((root / "campaign.json").read_text())
    complete = json.loads((root / "PIC_COMPLETE.json").read_text())
    reference_path = Path(campaign["reference_npz"])
    require(sha256(reference_path) == campaign["reference_npz_sha256"],
            "Reference map changed")
    reference = np.load(reference_path)
    x_old, y_old = reference["x_m"], reference["y_m"]
    rows, plots = [], []
    for case, closed in zip(campaign["cases"], complete["cases"]):
        label = case["label"]
        field = load_field(Path(case["directory"]) / "detectors/z03524.h5")
        ix, iy = index_subset(x_old, field["x"]), index_subset(y_old, field["y"])
        control = reference[f"{label}_bend_PIC_J_per_eV_m2"][np.ix_(iy, ix)]
        theory = reference[f"{label}_radiation_LW_J_per_eV_m2"][np.ix_(iy, ix)]
        xx, yy = np.meshgrid(field["x"], field["y"], indexing="xy")
        aperture = xx**2 + yy**2 <= campaign["aperture_radius_m"]**2
        row = {"label": label,
               **metrics(field["intensity"], theory, control, aperture,
                         field["dx"] * field["dy"]),
               "entry_records": closed["particle_plane_records"]["entry"],
               "bend_exit_records": closed["particle_plane_records"]["bend-exit"],
               "screen_exit_records": closed["particle_plane_records"]["exit"],
               "screen_exit_survival_fraction": closed["exit_survival_fraction"]}
        rows.append(row)
        plots.append({"label": label, "field": field, "control": control,
                      "theory": theory, "aperture": aperture, "metrics": row})

    csv_path = output / "bend_boundary_exit_8eV_metrics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    crop = ((abs(plots[0]["field"]["y"]) <= DISPLAY_RADIUS)[:, None] &
            (abs(plots[0]["field"]["x"]) <= DISPLAY_RADIUS)[None, :])
    crop_values = [values[crop] for entry in plots for values in
                   (entry["control"], entry["field"]["intensity"],
                    entry["theory"])]
    crop_vmax = max(float(value.max()) for value in crop_values) * 1e6
    crop_norm = colors.PowerNorm(gamma=.45, vmin=0., vmax=crop_vmax)
    differences = [(entry["field"]["intensity"] - entry["theory"]) /
                   max(float(entry["theory"][entry["aperture"]].max()), 1e-300)
                   for entry in plots]
    difference_limit = max(float(abs(value[entry["aperture"]]).max())
                           for entry, value in zip(plots, differences))
    fig, axes = plt.subplots(2, 5, figsize=(20, 7.7),
                             layout="constrained", squeeze=False)
    for row_index, (entry, difference) in enumerate(zip(plots, differences)):
        field = entry["field"]
        extent = [field["x"][0] * 1e6, field["x"][-1] * 1e6,
                  field["y"][0] * 1e6, field["y"][-1] * 1e6]
        full_norm = colors.PowerNorm(
            gamma=.35, vmin=0., vmax=float(field["intensity"].max() * 1e6))
        full_image = axes[row_index, 0].imshow(
            field["intensity"] * 1e6, origin="lower", extent=extent,
            aspect="equal", cmap="magma", norm=full_norm)
        axes[row_index, 0].axvline(campaign["right_cpml_start_m"] * 1e6,
                                   color="cyan", linestyle="--", linewidth=1)
        axes[row_index, 0].axvline(campaign["transverse_domain_m"][1] * 1e6,
                                   color="white", linestyle=":", linewidth=1)
        axes[row_index, 0].plot(
            campaign["predicted_beam_center_at_detector_m"] * 1e6, 0,
            marker="x", color="lime", clip_on=False)
        axes[row_index, 0].set_xlim(extent[0],
                                   campaign["predicted_beam_center_at_detector_m"] * 1e6 + 15)
        fig.colorbar(full_image, ax=axes[row_index, 0], shrink=.78,
                     label=r"$\mu$J/(eV m$^2$)")
        for column, values in enumerate(
                (entry["control"], field["intensity"], entry["theory"]), 1):
            image = axes[row_index, column].imshow(
                values * 1e6, origin="lower", extent=extent, aspect="equal",
                cmap="magma", norm=crop_norm)
            axes[row_index, column].set_xlim(-DISPLAY_RADIUS * 1e6,
                                             DISPLAY_RADIUS * 1e6)
            axes[row_index, column].set_ylim(-DISPLAY_RADIUS * 1e6,
                                             DISPLAY_RADIUS * 1e6)
            fig.colorbar(image, ax=axes[row_index, column], shrink=.78,
                         label=r"$\mu$J/(eV m$^2$)")
        diff_image = axes[row_index, 4].imshow(
            difference, origin="lower", extent=extent, aspect="equal",
            cmap="coolwarm", vmin=-difference_limit, vmax=difference_limit)
        axes[row_index, 4].set_xlim(-DISPLAY_RADIUS * 1e6,
                                    DISPLAY_RADIUS * 1e6)
        axes[row_index, 4].set_ylim(-DISPLAY_RADIUS * 1e6,
                                    DISPLAY_RADIUS * 1e6)
        fig.colorbar(diff_image, ax=axes[row_index, 4], shrink=.78,
                     label="(boundary PIC - LW) / LW aperture peak")
        circle = plt.Circle((0, 0), campaign["aperture_radius_m"] * 1e6,
                            fill=False, color="cyan", linewidth=.8)
        axes[row_index, 0].add_patch(circle)
        for column in range(1, 5):
            axes[row_index, column].add_patch(plt.Circle(
                (0, 0), campaign["aperture_radius_m"] * 1e6,
                fill=False, color="cyan", linewidth=.8))
        row = entry["metrics"]
        axes[row_index, 4].text(
            .02, .98,
            f"survives = {100 * row['screen_exit_survival_fraction']:.3f}%\n"
            f"PIC/LW = {row['boundary_PIC_over_acceleration_LW']:.5f}\n"
            f"cosine = {row['boundary_PIC_vs_acceleration_LW_cosine']:.5f}",
            transform=axes[row_index, 4].transAxes, va="top", fontsize=8.2,
            bbox={"facecolor": "white", "alpha": .84, "edgecolor": "none"})
        axes[row_index, 0].set_ylabel(entry["label"] + "\ny [µm]")
        for axis in axes[row_index]:
            axis.set_xlabel("x [µm]")
    titles = ("boundary-exit PIC: full plane", "collision-free PIC: aperture",
              "boundary-exit PIC: aperture", "acceleration LW: aperture",
              "boundary-exit PIC - LW")
    for axis, title in zip(axes[0], titles):
        axis.set_title(title)
    fig.suptitle(
        "8 eV bend-only control after particles hit the right transverse boundary")
    png_path = output / "bend_boundary_exit_8eV.png"
    pdf_path = output / "bend_boundary_exit_8eV.pdf"
    fig.savefig(png_path, dpi=210)
    fig.savefig(pdf_path)
    plt.close(fig)

    record = {
        "complete": True,
        "photon_energy_eV": ENERGY_EV,
        "aperture_radius_m": campaign["aperture_radius_m"],
        "rows": rows,
        "interpretation_guard": (
            "This trial retains transition radiation and numerical/CPML transients "
            "from particle removal; improvement is not assumed a priori."),
        "physical_accuracy_accepted": False,
        "figure": png_path.name,
        "figure_sha256": sha256(png_path),
        "pdf": pdf_path.name,
        "pdf_sha256": sha256(pdf_path),
        "csv": csv_path.name,
        "csv_sha256": sha256(csv_path),
        "reference_npz": str(reference_path),
        "reference_npz_sha256": sha256(reference_path),
    }
    (output / "FINAL.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 2,
            "usage: analyze_100k_bend_boundary_exit.py ROOT")
    main(sys.argv[1])
