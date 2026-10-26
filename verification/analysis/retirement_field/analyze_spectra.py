#!/usr/bin/env python3
"""Time-only spectra of raw lab E/B; no near-field or angular filtering.

Read each HDF5 once in time order, keep the physical aperture in memory, then
FFT spatial batches. Production reads belong on an allocated compute node.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import os
from pathlib import Path
import time

import h5py
import numpy as np
from scipy import fft
import yaml

C = 299792458.0
H_EV_S = 4.135667696e-15
MU0 = 1.25663706212e-6


def plain(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    raise TypeError(type(v).__name__)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, default=plain, allow_nan=False) + "\n")


def reduce_cube(cube, dt, area_cell, workers=1, pixel_batch=256):
    """cube has shape [time,pixel,6], with E in V/m and B in tesla."""
    n, pixels, components = cube.shape
    assert components == 6 and n >= 3 and pixels > 0
    duration = n * dt
    freq = fft.rfftfreq(n, dt)
    energy = H_EV_S * freq
    de = H_EV_S / duration
    one = np.full(freq.size, 2.0)
    one[0] = 1.0
    if n % 2 == 0:
        one[-1] = 1.0
    arrays = dict(frequency_Hz=freq, photon_energy_eV=energy)
    reviews = {}
    for name, window in (("rect", np.ones(n)), ("hann", np.hanning(n))):
        u = float(np.mean(window ** 2))
        window = window / np.sqrt(u)
        spectra = np.zeros((6, freq.size))
        flux = np.zeros(freq.size)
        time_ms = np.zeros(6)
        time_flux = 0.0
        for start in range(0, pixels, pixel_batch):
            stop = min(start + pixel_batch, pixels)
            values = np.array(cube[:, start:stop, :].transpose(1, 2, 0),
                              order="C", copy=True)
            values[:, 3:, :] *= C
            values *= window[None, None, :]
            time_ms += np.sum(values * values, axis=(0, 2)) / (pixels * n)
            time_flux += float(np.sum(values[:, 0] * values[:, 4] -
                                      values[:, 1] * values[:, 3])) * dt * area_cell / (MU0 * C)
            transformed = fft.rfft(values, axis=-1, workers=workers, overwrite_x=True)
            del values
            transformed *= dt
            spectra += np.sum(transformed.real ** 2 + transformed.imag ** 2, axis=0)
            flux += np.sum((transformed[:, 0] * np.conjugate(transformed[:, 4]) -
                            transformed[:, 1] * np.conjugate(transformed[:, 3])).real, axis=0)
            del transformed
            if start % (pixel_batch * 8) == 0:
                print(f"FFT {name}: {stop}/{pixels} pixels", flush=True)
        spectra *= one[None, :] / (pixels * duration * H_EV_S)
        flux *= one * area_cell / (MU0 * C * H_EV_S)
        spectral_ms = spectra.sum(axis=1) * de
        spectral_flux = float(flux.sum() * de)
        ms_error = np.abs(spectral_ms - time_ms) / np.maximum(time_ms, 1e-300)
        flux_scale = max(float(np.sum(np.abs(flux)) * de), abs(time_flux), 1e-300)
        flux_error = abs(spectral_flux - time_flux) / flux_scale
        if np.any(ms_error > 2e-10) or flux_error > 2e-10:
            raise ValueError(f"Parseval failed: {ms_error}, flux {flux_error}")
        arrays[f"{name}_E_psd_V2_m2_per_eV"] = spectra[:3].T
        arrays[f"{name}_cB_psd_V2_m2_per_eV"] = spectra[3:].T
        arrays[f"{name}_native_flux_J_per_eV"] = flux
        etotal = spectra[:3].sum(axis=0)
        btotal = spectra[3:].sum(axis=0)
        bands = {}
        for label, lo, hi in (("below_1_eV", 0., 1.), ("below_10_eV", 0., 10.),
                              ("10_to_50_eV", 10., 50.), ("50_to_114p352_eV", 50., 114.3520538),
                              ("114p352_to_200_eV", 114.3520538, 200.),
                              ("above_200_eV", 200., np.inf)):
            mask = (energy >= lo) & (energy < hi)
            es = float(etotal[mask].sum() * de)
            bs = float(btotal[mask].sum() * de)
            bands[label] = dict(E_mean_square_V2_m2=es, cB_mean_square_V2_m2=bs,
                               E_fraction=es / float(time_ms[:3].sum()),
                               cB_fraction=bs / float(time_ms[3:].sum()),
                               native_flux_J=float(flux[mask].sum() * de))
        band = (energy >= 50.) & (energy <= 114.3520538)
        peak = np.flatnonzero(band)[np.argmax(etotal[band])]
        reviews[name] = dict(window=name, original_window_mean_square=u,
                            amplitude_normalization="divide by sqrt(mean(window**2))",
                            mean_square_components_E_then_cB=time_ms,
                            spectral_mean_square_components=spectral_ms,
                            parseval_relative_error_components=ms_error,
                            native_flux_parseval_relative_error=flux_error,
                            native_flux_time_sum_J=time_flux,
                            native_flux_spectral_sum_J=spectral_flux,
                            dc_E_mean_square_V2_m2=float(etotal[0] * de),
                            dc_E_fraction=float(etotal[0] * de / time_ms[:3].sum()),
                            peak_50_to_114p352_eV=float(energy[peak]), bands=bands)
    return arrays, reviews


def analyze(run_dir, field_file, output_dir, slab_cells, workers, pixel_batch):
    started = time.monotonic()
    config = yaml.safe_load((run_dir / "run.yaml").read_text())
    boundary = config.get("boundary", {})
    cx, cy = boundary.get("cells", [0, 0, 0])[:2] if boundary.get("type") == "cpml" else (0, 0)
    with h5py.File(field_file, "r") as handle:
        group = handle["field_plane"]
        assert int(group["complete"][()]) == 1
        n = int(group["committed_samples"][()])
        data = group["fields"]
        assert 3 <= n <= data.shape[0]
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in group.attrs.items()}
        assert attrs["frame"] == "laboratory"
        assert attrs["external_background_subtracted"] == "none"
        times = group["time_s"][:n]
        assert np.isfinite(times).all() and np.all(np.diff(times) > 0)
        dt = float((times[-1] - times[0]) / (n - 1))
        linear = times[0] + np.arange(n) * dt
        jitter = float(np.max(np.abs(times - linear)) / dt)
        if jitter > 1e-4:
            raise ValueError(f"Time stamps are not uniform enough for a direct FFT: {jitter} samples")
        ny, nx = data.shape[1:]
        assert 0 <= 2 * cx < nx and 0 <= 2 * cy < ny
        px, py = nx - 2 * cx, ny - 2 * cy
        expected_bytes = n * px * py * 6 * 8
        if expected_bytes > 80 * 1024 ** 3:
            raise ValueError("Physical field cube exceeds the 80 GiB memory budget")
        cube = np.empty((n, py * px, 6), dtype=np.float64)
        for start in range(0, n, 32):
            stop = min(start + 32, n)
            records = data[start:stop]
            electric, magnetic = records["electric_V_per_m"], records["magnetic_T"]
            if not np.isfinite(electric).all() or not np.isfinite(magnetic).all():
                raise ValueError(f"Nonfinite field at block {start}")
            cube[start:stop, :, :3] = electric[:, cy:ny-cy, cx:nx-cx].reshape(stop-start, py*px, 3)
            cube[start:stop, :, 3:] = magnetic[:, cy:ny-cy, cx:nx-cx].reshape(stop-start, py*px, 3)
            if start % 4096 == 0:
                print(f"read {field_file.name}: {stop}/{n} frames", flush=True)
        del records, electric, magnetic
    read_seconds = time.monotonic() - started
    arrays, reviews = reduce_cube(cube, dt, float(attrs["dx_m"] * attrs["dy_m"]), workers, pixel_batch)
    del cube
    gc.collect()
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "spectra.npz", **arrays)
    columns = dict(photon_energy_eV=arrays["photon_energy_eV"], frequency_Hz=arrays["frequency_Hz"])
    for name in ("rect", "hann"):
        for kind in ("E", "cB"):
            vals = arrays[f"{name}_{kind}_psd_V2_m2_per_eV"]
            for k, letter in enumerate("xyz"):
                columns[f"{name}_{kind}{letter}_V2_m2_per_eV"] = vals[:, k]
        columns[f"{name}_native_flux_J_per_eV"] = arrays[f"{name}_native_flux_J_per_eV"]
    with (output_dir / "spectra.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    summary = dict(source_file=str(field_file.resolve()), run_directory=str(run_dir.resolve()),
                   source_size_bytes=field_file.stat().st_size, source_attributes=attrs,
                   samples=n, sample_interval_s=dt, time_start_s=float(times[0]), time_end_s=float(times[-1]),
                   recorded_span_fs=float((times[-1]-times[0])*1e15), fft_period_fs=n*dt*1e15,
                   maximum_timestamp_deviation_samples=jitter,
                   time_handling="Direct FFT of native samples; rounded timestamps fit uniform indices; no resampling",
                   frequency_bin_eV=H_EV_S/(n*dt), theoretical_nyquist_eV=H_EV_S/(2*dt),
                   last_frequency_eV=float(arrays["photon_energy_eV"][-1]),
                   cpml_excluded_cells_xy=[int(cx), int(cy)], physical_shape_yx=[py, px],
                   physical_area_m2=float(px*py*attrs["dx_m"]*attrs["dy_m"]),
                   slab_cells=slab_cells, sampling_comb_spacing_eV=H_EV_S/(slab_cells*dt),
                   cube_bytes=expected_bytes, fft_workers=workers, pixel_batch=pixel_batch,
                   read_seconds=read_seconds, total_seconds=time.monotonic()-started,
                   windows=reviews,
                   notes=["Raw fields: no near-field subtraction, mean removal, detrending, spike removal, angular projection or band filtering.",
                          "Each pixel is Fourier transformed before squared amplitudes are spatially averaged. DC remains in band accounting.",
                          "E and cB spectra are field mean-square densities, not separated radiation-energy spectra.",
                          "The Hann view uses identical energy-normalized windows on all six components; it tests sensitivity to finite-record edges, not complete pulse recovery.",
                          "Native spectral flux inherits Yee time staggering and the known MPI-interface diagnostic defect; retain its sign.",
                          "The FFT period is N*dt. Parseval uses discrete time sums, not the previous trapezoidal energy integrals.",
                          "Nyquist is a sampling bound, not evidence of physical/numerical resolution at high photon energies."])
    dump(output_dir / "summary.json", summary)
    print(f"COMPLETE {output_dir}: {summary['total_seconds']:.2f} s", flush=True)


def self_test():
    # Odd/even grids exercise the last one-sided bin; offset checks retained DC.
    for n in (99, 100):
        dt, area = 2e-17, 2e-12
        phase = 2*np.pi*7*np.arange(n)/n
        e = 3 + 2*np.cos(phase)
        cube = np.zeros((n, 3, 6))
        for pixel, scale in enumerate((1., -1., 2.)):
            cube[:, pixel, 0] = scale * e
            cube[:, pixel, 4] = scale * e / C
        arrays, reviews = reduce_cube(cube, dt, area, pixel_batch=2)
        ms = float(np.mean(cube[:, :, 0]**2))
        de = H_EV_S/(n*dt)
        assert np.isclose(arrays['rect_E_psd_V2_m2_per_eV'][:, 0].sum()*de, ms, rtol=1e-13)
        assert np.isclose(reviews['rect']['dc_E_mean_square_V2_m2'], 18., rtol=1e-13)
        assert np.isclose(reviews['rect']['native_flux_time_sum_J'], ms*3*area*n*dt/(MU0*C), rtol=1e-13)
        assert np.argmax(arrays['rect_E_psd_V2_m2_per_eV'][1:, 0])+1 == 7
    print("PASS: odd/even one-sided factors, retained DC, opposite-phase pixels, forward-wave flux, both-window Parseval", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--field-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--slab-cells", type=int)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--pixel-batch", type=int, default=256)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if any(v is None for v in (args.run_dir, args.field_file, args.output_dir, args.slab_cells)):
        parser.error("run-dir, field-file, output-dir and slab-cells are required")
    if min(args.workers, args.pixel_batch, args.slab_cells) < 1:
        parser.error("workers, pixel-batch and slab-cells must be positive")
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Full-field scans must run through Slurm")
    analyze(args.run_dir, args.field_file, args.output_dir, args.slab_cells, args.workers, args.pixel_batch)


if __name__ == "__main__":
    main()
