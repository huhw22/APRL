#!/usr/bin/env python3
"""Check the off-grid bin integral, time origin, fallback and file-safety contract."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True  # Keep the project source tree clean when CTest loads the tool.

try:
    import h5py
    import numpy as np
except ModuleNotFoundError as exc:
    if exc.name in ("numpy", "h5py"):
        print(f"Skipping optional screen spectrum checks: {exc}")
        sys.exit(77)
    raise


def make_screen(path, start, dt):
    ntime, ny, nx = 7, 3, 5
    # Distinct signed fields at all pixels and field components.
    base = np.arange(1., 1. + ny * nx * 3).reshape(ny, nx, 3)
    electric = np.stack([factor * base for factor in (0, 2, -1, 0, 0, -2, 0)])
    magnetic = np.stack([factor * base[..., ::-1] for factor in (1, 0, 0, -1, 0, 2, 0)])
    with h5py.File(path, "x") as f:
        f.attrs.update(complete=1, z_m=0.42, time_start_s=start, time_step_s=dt,
                       time_origin_s=start-dt/8, model="finite-distance vacuum Lienard-Wiechert acceleration term")
        f.create_dataset("time_s", data=start+(np.arange(ntime)+.5)*dt)
        f.create_dataset("x_m", data=np.linspace(-2.e-5,2.e-5,nx))
        f.create_dataset("y_m", data=np.linspace(-1.e-5,1.e-5,ny))
        f.create_dataset("E_time", data=electric).attrs["units"] = np.bytes_("V/m")
        f.create_dataset("B_time", data=magnetic).attrs["units"] = np.bytes_("T")
    return electric, magnetic


def analytic(bins, freq, start, dt, origin):
    """Independent antiderivative on each interval, using the bin-edge clocks."""
    offset = start-origin
    output = np.zeros((len(freq),)+bins.shape[1:], dtype=complex)
    for k, value in enumerate(freq):
        if value==0:
            output[k] = dt * bins.sum(axis=0)
            continue
        iw = 2j*np.pi*value
        for t, sample in enumerate(bins):
            left, right = offset+dt*t, offset+dt*(t+1)
            output[k] += sample * (np.exp(iw*right)-np.exp(iw*left))/iw
    return output


def run(script, test_dir):
    spec = importlib.util.spec_from_file_location("aprl_finite_screen_spectrum", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Exercise x, y, time and frequency tiling with a tiny manufactured screen.
    module.TIME_BLOCK, module.PIXEL_BLOCK, module.FREQUENCY_BLOCK = 3, 4, 2
    path, output = test_dir/"radiation.h5", test_dir/"spectra.h5"
    start, dt = 2.**-13, 2.**-45
    arrays = make_screen(path, start, dt)
    frequencies = np.array([0, 0.23/dt, 0.5/dt, 1.0/dt, 0.76/dt])
    with h5py.File(path, "a") as src:
        src.create_dataset("requested_frequency_hz", data=frequencies)
        src.create_dataset("frequency_hz", data=[2.e12])
    module.write_spectra(path, output)
    with h5py.File(output) as result, h5py.File(path) as src:
        assert result.attrs["complete"] == 1
        assert result.attrs["z_m"] == src.attrs["z_m"]
        assert result.attrs["time_origin_s"] == src.attrs["time_origin_s"]
        for axis in ("x_m", "y_m"):
            np.testing.assert_array_equal(result[axis][:], src[axis][:])
        np.testing.assert_array_equal(result["frequency_hz"][:], frequencies)
        for key, field, unit in zip(("E_spectrum","B_spectrum"), arrays, (b"V s/m", b"T s")):
            ds = result[key]
            assert ds.shape == (5,3,5,3,2) and ds.attrs["units"] == unit
            data = ds[...,0] + 1j*ds[...,1]
            expect = analytic(field, frequencies, start, dt, result.attrs["time_origin_s"])
            np.testing.assert_allclose(data, expect, rtol=5e-14, atol=dt*1e-12)
            np.testing.assert_allclose(data[0].imag, 0.0, atol=0.0)
            np.testing.assert_allclose(data[3], 0.0, atol=dt*1e-13)  # full cycle per bin

    def invalid(destination, frequencies=None):
        try:
            module.write_spectra(path, destination, frequencies)
        except (ValueError, OSError):
            return
        raise AssertionError("invalid conversion accepted")

    original_size = output.stat().st_size
    invalid(output)
    assert output.stat().st_size == original_size
    invalid(path)
    linked = test_dir/"same-input.h5"
    os.link(path, linked)
    invalid(linked)

    # Existing archives may list frequencies as frequency_hz only.
    with h5py.File(path, "a") as src:
        del src["requested_frequency_hz"]
    archived = test_dir/"archived.h5"
    module.write_spectra(path, archived)
    with h5py.File(archived) as f:
        np.testing.assert_array_equal(f["frequency_hz"][:], [2.e12])

    # Explicit CLI frequencies/phase origin override metadata.
    override = test_dir/"override.h5"
    new_origin = start+dt/4
    subprocess.run([sys.executable, "-B", str(script), str(path), str(override),
                    "--time-origin-s", str(new_origin), "--frequencies-hz", "0", str(frequencies[1])], check=True)
    with h5py.File(override) as f:
        np.testing.assert_array_equal(f["frequency_hz"][:], frequencies[:2])
        data = f["B_spectrum"][...,0]+1j*f["B_spectrum"][...,1]
        expect = analytic(arrays[1], frequencies[:2], start, dt, new_origin)
        np.testing.assert_allclose(data, expect, rtol=5e-14, atol=dt*1e-12)

    never_created = test_dir/"invalid.h5"
    with h5py.File(path, "a") as src:
        del src["frequency_hz"]
    invalid(never_created)
    for bad in ([-1.], [np.nan], [1.,1.]):
        invalid(never_created, bad)
    with h5py.File(path, "a") as src:
        src["E_time"][-1,-1,-1,-1] = np.nan
    invalid(never_created, [0.])
    assert not never_created.exists() and not list(test_dir.glob(".screen-spectrum-*"))
    with h5py.File(path, "a") as src:
        src["E_time"][-1,-1,-1,-1] = 0.
        src["time_s"][3] += dt
    invalid(never_created, [0.])
    with h5py.File(path, "a") as src:
        src["time_s"][3] -= dt
        src.attrs["complete"] = 0
    invalid(never_created, [0.])
    with h5py.File(path, "a") as src:
        src.attrs["complete"] = 1
        src["B_time"].attrs["units"] = "V/m"
    invalid(never_created, [0.])
    with h5py.File(path, "a") as src:
        src["B_time"].attrs["units"] = "T"
        del src.attrs["time_origin_s"]
    no_origin = test_dir/"default-origin.h5"
    module.write_spectra(path, no_origin, [frequencies[1]])
    with h5py.File(no_origin) as f:
        assert f.attrs["time_origin_s"] == start
        complex_e = f["E_spectrum"][...,0] + 1j*f["E_spectrum"][...,1]
        expect = analytic(arrays[0], frequencies[1:2], start, dt, start)
        np.testing.assert_allclose(complex_e, expect, rtol=5e-14, atol=dt*1e-12)


if __name__ == "__main__":
    script = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/"postprocess/finite_screen_spectrum.py"
    with tempfile.TemporaryDirectory(prefix="aprl-spectrum-test-") as test_dir:
        run(script, Path(test_dir))
    print("Passed finite-screen waveform spectrum checks.")
