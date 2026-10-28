#!/usr/bin/env python3
"""Independently audit small saved spectra, without rereading fields or FFTs.

Use extended-precision integration of the saved per-eV arrays and compare it
with the producer's independently accumulated time-domain sums. This verifies
the reduction and units, not complete pulse coverage or radiation separation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

H_EV_S = np.longdouble("4.135667696e-15")
MU0 = np.longdouble("1.25663706212e-6")
C = np.longdouble(299792458)
TINY = np.longdouble("1e-300")
CASES = ("undulator", "free_drift")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def review_case(directory):
    summary_file = directory / "summary.json"
    spectra_file = directory / "spectra.npz"
    summary = json.loads(summary_file.read_text())
    n = int(summary["samples"])
    dt = np.longdouble(summary["sample_interval_s"])
    duration = n * dt
    de = H_EV_S / duration
    require(n >= 3 and dt > 0, "Invalid sample count or cadence")
    with np.load(spectra_file) as saved:
        energy = saved["photon_energy_eV"].astype(np.longdouble)
        frequency = saved["frequency_Hz"].astype(np.longdouble)
        require(len(energy) == n // 2 + 1, "One-sided frequency extent differs")
        require(frequency.shape == energy.shape, "Frequency/energy shapes differ")
        require(np.isfinite(energy).all() and np.isfinite(frequency).all(),
                "Nonfinite frequency axis")
        require(np.all(np.diff(energy) > 0), "Energy bins are not increasing")
        expected = np.arange(len(energy), dtype=np.longdouble) * de
        axis_error_bins = np.max(abs(energy - expected)) / de
        require(axis_error_bins < 1e-8, "Frequency grid does not match N*dt")
        hf_error = np.max(abs(energy - H_EV_S * frequency)) / max(energy)
        require(hf_error < 1e-13, "Energy axis differs from h*f")
        require(abs(de - summary["frequency_bin_eV"]) / de < 1e-13,
                "Summary frequency bin differs")
        area = np.longdouble(summary["physical_area_m2"])
        pixels = np.prod(summary["physical_shape_yx"])
        attrs = summary["source_attributes"]
        expected_area = pixels * np.longdouble(attrs["dx_m"]) * np.longdouble(attrs["dy_m"])
        require(area > 0 and abs(area - expected_area) / area < 1e-13,
                "Physical area differs from pixel count and spacing")
        windows = {}
        for window in ("rect", "hann"):
            electric = saved[window + "_E_psd_V2_m2_per_eV"].astype(np.longdouble)
            magnetic = saved[window + "_cB_psd_V2_m2_per_eV"].astype(np.longdouble)
            flux = saved[window + "_native_flux_J_per_eV"].astype(np.longdouble)
            require(electric.shape == magnetic.shape == (len(energy), 3),
                    "Field spectrum dimensions differ")
            require(flux.shape == (len(energy),), "Flux spectrum dimensions differ")
            require(all(np.isfinite(a).all() for a in (electric, magnetic, flux)),
                    "Nonfinite saved spectral value")
            require((electric >= 0).all() and (magnetic >= 0).all(),
                    "Negative squared-field spectral density")
            time_summary = summary["windows"][window]
            direct = np.array(time_summary["mean_square_components_E_then_cB"],
                              dtype=np.longdouble)
            integrated = np.r_[electric.sum(axis=0) * de,
                               magnetic.sum(axis=0) * de]
            error = abs(integrated - direct) / np.maximum(abs(direct), TINY)
            time_flux = np.longdouble(time_summary["native_flux_time_sum_J"])
            total_flux = flux.sum() * de
            flux_error = abs(total_flux - time_flux) / max(
                np.sum(abs(flux)) * de, abs(time_flux), TINY)
            require(max(error) < 2e-13 and flux_error < 2e-13,
                    "Independent Parseval check failed")
            bound = area * duration / (MU0 * C) * np.sqrt(
                electric[:, :2].sum(axis=1) * magnetic[:, :2].sum(axis=1))
            require(np.all(flux[bound == 0] == 0), "Nonzero flux with zero Cauchy bound")
            ratio = np.divide(abs(flux), bound, out=np.zeros_like(flux), where=bound > 0)
            require(ratio.max() < 1 + 1e-12, "Per-bin flux Cauchy bound failed")
            target = (energy >= 50) & (energy <= 115)
            indices = np.flatnonzero(target)
            require(len(indices) > 0, "No samples in requested peak interval")
            peak = int(indices[np.argmax(electric[target].sum(axis=1))])
            flux_peak = int(indices[np.argmax(flux[target])])
            bands = {}
            for label, lo, hi in (("91p3_to_91p7", 91.3, 91.7),
                                  ("50_to_115", 50, 115), ("below_1", 0, 1)):
                selected = (energy >= lo) & (energy <= hi)
                require(selected.any(), "Empty requested integration interval")
                band_energy = flux[selected].sum() * de
                bands[label] = {
                    "request_eV": [lo, hi],
                    "bin_selection": "Full bins whose centres lie in inclusive interval; no partial-edge interpolation",
                    "selected_bins": int(selected.sum()),
                    "first_last_bin_centres_eV": [float(energy[selected][0]), float(energy[selected][-1])],
                    "integrated_bin_width_eV": float(selected.sum() * de),
                    "native_flux_energy_J": float(band_energy),
                    "record_average_native_flux_W": float(band_energy / duration),
                    "E_mean_square_V2_m2": float(electric[selected].sum() * de),
                    "cB_mean_square_V2_m2": float(magnetic[selected].sum() * de),
                }
            windows[window] = {
                "peak_50_to_115_eV": float(energy[peak]),
                "peak_E_sum_density_V2_m2_per_eV": float(electric[peak].sum()),
                "native_flux_peak_50_to_115_eV": float(energy[flux_peak]),
                "native_flux_at_E_peak_J_per_eV": float(flux[peak]),
                "parseval_relative_error_components": error.astype(float).tolist(),
                "native_flux_parseval_relative_error": float(flux_error),
                "native_flux_total_J": float(total_flux),
                "native_flux_total_record_average_W": float(total_flux / duration),
                "cauchy_bound_max_ratio": float(ratio.max()),
                "cauchy_bound_formula": "abs(dW/de) <= area*T/(mu0*c)*sqrt((D_Ex+D_Ey)*(D_cBx+D_cBy))",
                "bands": bands,
            }
        array_keys = list(saved.keys())
    return {
        "samples": n, "duration_Ndt_s": float(duration), "frequency_bin_eV": float(de),
        "frequency_array_max_deviation_bins": float(axis_error_bins),
        "photon_energy_h_times_f_max_relative_error": float(hf_error),
        "array_keys": array_keys, "summary_sha256": digest(summary_file),
        "spectra_sha256": digest(spectra_file), "windows": windows,
        "hann_to_rect_91p3_to_91p7_energy_ratio":
            windows["hann"]["bands"]["91p3_to_91p7"]["native_flux_energy_J"] /
            windows["rect"]["bands"]["91p3_to_91p7"]["native_flux_energy_J"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = {}
    for case in CASES:
        directory = args.root / case
        if (directory / "summary.json").exists() and (directory / "spectra.npz").exists():
            cases[case] = review_case(directory)
    require(cases, "Neither complete saved spectral pair exists")
    pending = [case for case in CASES if case not in cases]
    report = {
        "status": "PARTIAL" if pending else "PASS",
        "status_scope": "Available derived spectral arrays only; not production-job, pure-radiation or quiet-tail acceptance",
        "pending_cases": pending,
        "read_scope": "Small summaries and NPZ only; no HDF5 or FFT. Independent longdouble integration versus stored time-domain sums.",
        "longdouble_precision_bits": int(np.finfo(np.longdouble).nmant + 1),
        "script_sha256": digest(Path(__file__)), "cases": cases,
        "units_review": {
            "E_and_cB": "Pixel- and record-time-averaged squared-field density per eV: (V/m)^2/eV. cB includes c scaling. Integrate with dE once.",
            "native_flux": "Signed whole-physical-plane E x B spectrum, already area integrated: J/eV. Multiply dE once for J; divide Ndt once for record-average W. Never multiply area again.",
            "window": "Rect retains raw data. Hann divides samples by sqrt(mean(window^2)); it changes weighting of a nonstationary finite record, not pulse recovery.",
            "limitations": "Native time stagger remains. Total fields are not radiation-only; the 91.3-91.7 eV band includes background with no pedestal subtraction.",
        },
    }
    if not pending:
        ratios = {}
        for window in ("rect", "hann"):
            a = cases["undulator"]["windows"][window]["bands"]["91p3_to_91p7"]
            b = cases["free_drift"]["windows"][window]["bands"]["91p3_to_91p7"]
            ratios[window] = {
                "undulator_over_free_drift_native_band_energy": a["native_flux_energy_J"] / b["native_flux_energy_J"],
                "undulator_over_free_drift_record_average_native_band_power": a["record_average_native_flux_W"] / b["record_average_native_flux_W"],
                "comparison_scope": "Different Ndt and bin grids retained; not radiation enhancement or exact background subtraction",
            }
        report["band_comparison"] = ratios
    metrics_file = args.root / "figures" / "METRICS.json"
    if metrics_file.exists():
        metrics = json.loads(metrics_file.read_text())
        checks = 0
        for case, reviewed in cases.items():
            plotted = metrics["runs"][case]
            for name in ("summary", "spectra"):
                require(plotted["input_hashes"][name + (".json" if name == "summary" else ".npz")] ==
                        reviewed[name + "_sha256"], "Plot input hash differs")
                checks += 1
            for window, values in reviewed["windows"].items():
                other = plotted["windows"][window]
                band = other["bands"]["91p3_to_91p7_eV"]
                ours = values["bands"]["91p3_to_91p7"]
                for a, b in (
                    (ours["native_flux_energy_J"], band["native_energy_J"]),
                    (ours["record_average_native_flux_W"], band["record_average_power_W"]),
                    (values["native_flux_total_J"], other["native_total_energy_J"]),
                    (values["native_flux_total_record_average_W"], other["record_average_total_power_W"]),
                    (values["peak_50_to_115_eV"], other["peaks"]["target_50_115_eV_E_density"]["energy_eV"]),
                    (values["native_flux_peak_50_to_115_eV"], other["peaks"]["target_50_115_eV_native_flux"]["energy_eV"]),
                ):
                    require(abs(a - b) <= 1e-13 * max(abs(a), abs(b), 1e-300),
                            "Plot metric differs from independent computation")
                    checks += 1
        report["plot_metrics_review"] = {
            "status": "PASS", "checks": checks, "metrics_sha256": digest(metrics_file),
            "scope": "Input hashes, target-band energy/power, full-band energy/power and both field/flux peak positions",
        }
    output = args.output or args.root / "INDEPENDENT_REVIEW.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": report["status"], "pending_cases": pending, "output": str(output)}))


if __name__ == "__main__":
    main()
