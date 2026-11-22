#!/usr/bin/env python3
"""Compute signed complex E/B spectra from time-bin means of an APRL radiation screen."""
import argparse
import os
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np


# Tile all three potentially large sets of indices. A field tile is at most
# 128 x 2048 x 3 doubles and each accumulator holds at most 16 frequency planes.
TIME_BLOCK = 128
PIXEL_BLOCK = 2048
FREQUENCY_BLOCK = 16
UNITS = {"E": ("V/m", "V s/m"), "B": ("T", "T s")}


def real_scalar(attrs, key):
    if key not in attrs:
        raise ValueError(f"missing scalar attribute {key}")
    value = np.asarray(attrs[key])
    if value.shape or value.dtype.kind not in "fiu" or not np.isfinite(value):
        raise ValueError(f"{key} must be a finite scalar")
    return float(value)


def validate_screen(src):
    if real_scalar(src.attrs, "complete") != 1:
        raise ValueError("input radiation screen did not complete")
    start = real_scalar(src.attrs, "time_start_s")
    dt = real_scalar(src.attrs, "time_step_s")
    if dt <= 0:
        raise ValueError("time_step_s must be positive")
    real_scalar(src.attrs, "z_m")
    for name in ("time_s", "x_m", "y_m"):
        axis = src.get(name)
        if not isinstance(axis, h5py.Dataset) or axis.ndim != 1 or not len(axis) or axis.dtype.kind != "f":
            raise ValueError(f"{name} must be a nonempty floating-point coordinate vector")
    nt, ny, nx = (len(src[name]) for name in ("time_s", "y_m", "x_m"))
    for name in ("x_m", "y_m"):
        axis = src[name][:]
        if not np.all(np.isfinite(axis)) or np.any(np.diff(axis) <= 0):
            raise ValueError(f"invalid pixel centres in {name}")
    for field, (unit, _) in UNITS.items():
        ds = src.get(field + "_time")
        if not isinstance(ds, h5py.Dataset) or ds.shape != (nt, ny, nx, 3) or ds.dtype.kind != "f":
            raise ValueError(f"{field}_time must have shape (time, y, x, 3) and floating-point values")
        stored_unit = ds.attrs.get("units")
        if isinstance(stored_unit, bytes):
            stored_unit = stored_unit.decode("ascii")
        if stored_unit != unit:
            raise ValueError(f"{field}_time must have units {unit}")
    # The native writer stores the rounded first centre then uses fma for
    # each subsequent centre. At large absolute lab times the displayed time
    # coordinate can lose precision compared with the explicit bin width.
    first = start + 0.5 * dt
    for i in range(0, nt, 16384):
        stop = min(i + 16384, nt)
        expected = (np.longdouble(first) + np.arange(i, stop, dtype=np.longdouble) * dt).astype("f8")
        observed = src["time_s"][i:stop]
        tol = 2 * np.spacing(np.maximum(np.abs(expected), abs(dt)))
        if not np.all(np.isfinite(observed)) or not np.all(np.abs(observed - expected) <= tol):
            raise ValueError("time_s centres do not agree with time_start_s and time_step_s")
    return start, dt, nt, ny, nx


def requested_frequencies(src, explicit):
    if explicit is None:
        for name in ("requested_frequency_hz", "frequency_hz"):
            if name in src:
                ds = src[name]
                if not isinstance(ds, h5py.Dataset) or ds.ndim != 1 or ds.dtype.kind not in "fiu":
                    raise ValueError(f"{name} must be a numeric frequency vector")
                explicit = ds[:]
                break
    if explicit is None or len(explicit) == 0:
        raise ValueError("no frequency list in input; give --frequencies-hz 0 1e12 or other frequencies in Hz")
    freq = np.asarray(explicit, dtype=np.float64)
    if not np.all(np.isfinite(freq)) or np.any(freq < 0) or np.unique(freq).size != freq.size:
        raise ValueError("frequencies in Hz must be finite, nonnegative and distinct")
    return freq


def bin_weights(freq, start, dt, origin, t0, t1):
    # Integrate exp(+i*2*pi*f*(T-origin)) over each constant-field bin.
    # Reduce the phase of the reference offset before adding small bin delays.
    frequency = freq.astype(np.longdouble)
    cycles = frequency * (np.longdouble(start) - np.longdouble(origin))
    phase0 = cycles - np.rint(cycles)
    fdt = frequency * np.longdouble(dt)
    centres = np.arange(t0, t1, dtype=np.longdouble) + 0.5
    cycles = fdt[:, None] * centres
    phase = (cycles - np.rint(cycles) + phase0[:, None]).astype(np.float64)
    return (dt * np.sinc(fdt)).astype(np.float64)[:, None] * np.exp((2j * np.pi) * phase)


def transform(src, dst, freq, origin, start, dt, nt, ny, nx):
    shape = (len(freq), ny, nx, 3, 2)
    out = {}
    for field, (_, unit) in UNITS.items():
        out[field] = dst.create_dataset(field + "_spectrum", shape=shape, dtype="<f8")
        out[field].attrs["units"] = np.bytes_(unit)
    tile_x = min(nx, PIXEL_BLOCK)
    tile_y = max(1, PIXEL_BLOCK // tile_x)
    for i in range(0, len(freq), FREQUENCY_BLOCK):
        selected = freq[i:i + FREQUENCY_BLOCK]
        for y in range(0, ny, tile_y):
            for x in range(0, nx, tile_x):
                sy = slice(y, min(y + tile_y, ny))
                sx = slice(x, min(x + tile_x, nx))
                tile = (len(selected), sy.stop - y, sx.stop - x, 3)
                acc = {field: np.zeros(tile, dtype=np.complex128) for field in UNITS}
                for t0 in range(0, nt, TIME_BLOCK):
                    t1 = min(t0 + TIME_BLOCK, nt)
                    weights = bin_weights(selected, start, dt, origin, t0, t1)
                    for field in UNITS:
                        values = src[field + "_time"][t0:t1, sy, sx, :]
                        if not np.all(np.isfinite(values)):
                            raise ValueError(f"nonfinite value in {field}_time")
                        acc[field] += np.einsum("ft,tyxc->fyxc", weights, values)
                for field in UNITS:
                    if not np.all(np.isfinite(acc[field])):
                        raise ValueError(f"nonfinite Fourier integral in {field}_spectrum")
                    out[field][i:i + len(selected), sy, sx, :, :] = acc[field].view(np.float64).reshape(tile + (2,))


def write_spectra(input_path, output_path, frequencies=None, origin=None):
    if os.path.lexists(output_path):
        raise FileExistsError(f"output already exists: {output_path}; supply a new output filename")
    with h5py.File(input_path, "r") as src:
        start, dt, nt, ny, nx = validate_screen(src)
        freq = requested_frequencies(src, frequencies)
        if origin is None:
            origin = real_scalar(src.attrs, "time_origin_s") if "time_origin_s" in src.attrs else start
        if not np.isfinite(origin):
            raise ValueError("time origin must be finite")
        # Publish only a completed file; concurrent invocations cannot replace
        # a previous result or truncate the input via an alternate path.
        with tempfile.TemporaryDirectory(prefix=".screen-spectrum-", dir=Path(output_path).parent) as tmpdir:
            tmp = Path(tmpdir) / "spectrum.h5"
            with h5py.File(tmp, "x") as dst:
                for key in ("frame", "model", "z_m"):
                    if key in src.attrs:
                        dst.attrs[key] = src.attrs[key]
                dst.attrs["time_origin_s"] = origin
                dst.attrs["time_start_s"] = start
                dst.attrs["time_step_s"] = dt
                dst.attrs["source_time_bins"] = nt
                dst.attrs["source_screen"] = str(Path(input_path).resolve())
                dst.attrs["spectrum_order"] = np.bytes_("frequency,y,x,component,real_imag")
                dst.attrs["fourier_transform"] = np.bytes_("integral F(T) exp(+i*2*pi*f*(T-time_origin_s)) dT")
                dst.attrs["reconstruction"] = np.bytes_("piecewise constant bin means; integral over saved arrival window")
                src.copy("x_m", dst)
                src.copy("y_m", dst)
                dst.create_dataset("frequency_hz", data=freq, dtype="<f8")
                transform(src, dst, freq, origin, start, dt, nt, ny, nx)
                dst.attrs["complete"] = 1
            os.link(tmp, output_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="completed time-domain radiation screen HDF5")
    parser.add_argument("output", type=Path, help="new spectral HDF5 file in an existing directory")
    parser.add_argument("--frequencies-hz", nargs="+", type=float, metavar="HZ",
                        help="nonnegative frequencies, overriding requested_frequency_hz or legacy frequency_hz")
    parser.add_argument("--time-origin-s", type=float, metavar="SECONDS",
                        help="Fourier phase origin in lab seconds (default: input time_origin_s or first bin edge)")
    args = parser.parse_args()
    try:
        write_spectra(args.input, args.output, args.frequencies_hz, args.time_origin_s)
    except (ValueError, OSError, TypeError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
