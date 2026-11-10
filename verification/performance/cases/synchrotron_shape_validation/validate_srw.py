"""Independent arithmetic validation for the SRW comparison products."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.array([float(row[name]) for row in rows])
            for name in rows[0]}


def main(root, result):
    root, result = Path(root).resolve(), Path(result).resolve()
    single = load(result / "single_electron_comparison.csv")
    gaussian = load(result / "gaussian_comparison.csv")
    final = json.loads((result / "FINAL.json").read_text())
    campaign = json.loads((root / "campaign.json").read_text())
    original = load(root / "axis_analysis_r2/axis_spectra.csv")
    energies = gaussian["photon_eV"]
    require(len(energies) == 10 and np.array_equal(
            energies, single["photon_eV"]), "Energy grids differ")
    require(np.array_equal(energies, original["photon_eV"]),
            "SRW and original PIC grids differ")
    require(all(np.isfinite(values).all() for values in single.values()),
            "Nonfinite single-electron value")
    require(all(np.isfinite(values).all() for values in gaussian.values()),
            "Nonfinite Gaussian value")

    number = final["coherent_gaussian"]["represented_electrons"]
    rebuilt = (single["SRW_J_per_eV_m2_per_electron"] * number**2 *
               gaussian["longitudinal_form_factor_squared"] *
               gaussian["LW_transverse_correction"])
    reconstruction = float(np.max(abs(
        rebuilt - gaussian["SRW_hybrid_transverse_corrected_J_per_eV_m2"]) /
        np.maximum(abs(rebuilt), 1e-300)))
    require(reconstruction < 3e-15, "Hybrid SRW reconstruction failed")
    require(np.array_equal(gaussian["finite_arc_LW_J_per_eV_m2"],
                           original["theory_J_per_eV_m2"]),
            "Original LW curve changed")
    for name in ["linear", "quintic-compensated", "linear-notch",
                 "linear-binomial3"]:
        require(np.array_equal(gaussian[name], original[name]),
                f"Original PIC curve changed for {name}")

    band_limits = campaign["common"]["registered_band_eV"]
    band = (energies >= band_limits[0]) & (energies <= band_limits[1])
    spacing = energies[1] - energies[0]
    band_value = float(np.sum(
        gaussian["SRW_hybrid_transverse_corrected_J_per_eV_m2"][band]) * spacing)
    stored_value = final["coherent_gaussian"]["SRW_hybrid_band_J_per_m2"]
    band_error = abs(band_value - stored_value) / stored_value
    require(band_error < 3e-15, "Stored SRW band integral differs")

    summary_rows = list(csv.DictReader((result / "comparison_summary.csv").open()))
    require(len(summary_rows) == 5, "Unexpected summary row count")
    stored = {row["curve"]: row for row in
              final["coherent_gaussian"]["curves"]}
    summary_error = 0.
    for row in summary_rows:
        name = row["curve"]
        for key in ["band_J_per_m2", "band_over_SRW",
                    "registered_spectral_relative_L2_to_SRW"]:
            lhs, rhs = float(row[key]), float(stored[name][key])
            summary_error = max(summary_error,
                                abs(lhs - rhs) / max(abs(rhs), 1e-300))
    require(summary_error < 3e-15, "CSV and JSON summaries differ")
    require((result / "srw_comparison.png").stat().st_size > 100000,
            "PNG is unexpectedly small")
    require((result / "srw_comparison.pdf").stat().st_size > 10000,
            "PDF is unexpectedly small")

    validation = {
        "passed": True,
        "energy_bins": len(energies),
        "identical_original_frequency_grid": True,
        "finite_all_numeric_columns": True,
        "hybrid_reconstruction_maximum_relative_error": reconstruction,
        "band_integral_relative_error": band_error,
        "summary_csv_vs_json_maximum_relative_error": summary_error,
        "original_LW_and_four_PIC_curves_bitwise_preserved": True,
        "files_sha256": {
            name: sha256(result / name) for name in [
                "single_electron_comparison.csv", "gaussian_comparison.csv",
                "comparison_summary.csv", "FINAL.json", "srw_comparison.png",
                "srw_comparison.pdf", "RESULT_REVIEW.md"]
        },
    }
    (result / "VALIDATION.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: validate_srw.py CAMPAIGN_ROOT RESULT")
    main(sys.argv[1], sys.argv[2])
