"""Render the accepted two-seed 100k same-particle LW comparison."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(source_csv, output):
    source_csv = Path(source_csv).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = np.genfromtxt(source_csv, delimiter=",", names=True)
    energy = np.asarray(data["photon_eV"])
    if not (len(energy) == 321 and np.array_equal(
            energy, np.linspace(4., 12., 321))):
        raise RuntimeError("Unexpected photon-energy grid")
    names = ["100k_seed11", "100k_seed12"]
    for name in names:
        for suffix in ("PIC_J_per_eV_m2", "sample_LW_J_per_eV_m2",
                       "continuum_J_per_eV_m2"):
            values = np.asarray(data[f"{name}_{suffix}"])
            if not (np.isfinite(values).all() and np.all(values > 0)):
                raise RuntimeError(f"Invalid spectrum: {name}_{suffix}")

    scale = 1e16 / 100000.**2
    colors = {"100k_seed11": "#1f77b4", "100k_seed12": "#d97706"}
    labels = {"100k_seed11": "100k seed 20261111",
              "100k_seed12": "100k seed 20261112"}
    fig, axis = plt.subplots(figsize=(10.5, 6.4), layout="constrained")
    axis.plot(energy, data["100k_seed11_continuum_J_per_eV_m2"] * scale,
              color="#222222", lw=2.4,
              label="ideal continuous Gaussian SRW")
    for name in names:
        color = colors[name]
        axis.plot(energy, data[f"{name}_sample_LW_J_per_eV_m2"] * scale,
                  color=color, lw=2.5,
                  label=f"{labels[name]}: same-particle LW")
        axis.plot(energy, data[f"{name}_PIC_J_per_eV_m2"] * scale,
                  color=color, lw=1.8, ls="--",
                  label=f"{labels[name]}: APRL PIC")
    axis.set(
        title="100k physical-electron samples: APRL follows the matched-particle LW field",
        xlabel="photon energy [eV]",
        ylabel=r"$10^{-16}$ J / (eV m$^2$ electron$^2$)")
    axis.grid(alpha=.22)
    axis.legend(frameon=True, fontsize=9, ncol=1)
    axis.text(
        .015, .02,
        "Solid: coherent LW sum of the exact PIC particles   "
        "Dashed: raw APRL bend-minus-straight spectrum",
        transform=axis.transAxes, fontsize=8.5, color="#444444")
    png = output / "correct_100k_sample_matched.png"
    pdf = output / "correct_100k_sample_matched.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)

    csv_path = output / "correct_100k_sample_matched.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "photon_eV", "ideal_continuous_Gaussian_scaled",
            "seed20261111_same_particle_LW_scaled", "seed20261111_APRL_scaled",
            "seed20261112_same_particle_LW_scaled", "seed20261112_APRL_scaled"])
        writer.writerows(zip(
            energy,
            data["100k_seed11_continuum_J_per_eV_m2"] * scale,
            data["100k_seed11_sample_LW_J_per_eV_m2"] * scale,
            data["100k_seed11_PIC_J_per_eV_m2"] * scale,
            data["100k_seed12_sample_LW_J_per_eV_m2"] * scale,
            data["100k_seed12_PIC_J_per_eV_m2"] * scale))
    record = {
        "complete": True,
        "scope": "two independent 100k physical-electron samples only",
        "source_csv": str(source_csv),
        "source_csv_sha256": sha256(source_csv),
        "normalization": "display only: spectrum * 1e16 / 100000^2",
        "fitted_normalization": False,
        "ideal_continuous_curve_is_acceptance_reference": False,
        "same_particle_LW_is_registered_reference": True,
        "png_sha256": sha256(png),
        "pdf_sha256": sha256(pdf),
        "csv_sha256": sha256(csv_path),
    }
    (output / "correct_100k_sample_matched.json").write_text(
        json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: plot_100k_sample_matched.py SOURCE_CSV OUTPUT")
    main(sys.argv[1], sys.argv[2])
