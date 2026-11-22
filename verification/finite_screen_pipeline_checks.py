#!/usr/bin/env python3
"""Check that a real committed-step screen leaves both release solvers unchanged."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
try:
    import h5py
    import numpy as np
except ModuleNotFoundError as exc:
    if exc.name in ("numpy", "h5py"):
        print(f"Skipping optional screen pipeline checks: {exc}")
        sys.exit(77)
    raise


def run(command, work):
    result = subprocess.run([str(x) for x in command], cwd=work, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
    assert result.returncode == 0, result.stdout
    return result.stdout


def compare(left, right):
    ignored = {"moving_energy_diagnostic_rank_seconds", "constraint_energy_diagnostic_rank_seconds"}
    with h5py.File(left) as a, h5py.File(right) as b:
        def check(name, obj):
            if not isinstance(obj, h5py.Dataset):
                return
            assert name in b and obj.shape == b[name].shape, name
            x, y = obj[...], b[name][...]
            for field in obj.dtype.names or [None]:
                if field in ignored:
                    continue
                xx, yy = (x, y) if field is None else (x[field], y[field])
                np.testing.assert_array_equal(xx, yy, err_msg=f"{name}/{field}")
            if name.endswith("/complete"):
                assert int(x) == 1, name
        a.visititems(check)


def main(aprl, generator, card, spectrum_script, work):
    work.mkdir(parents=True, exist_ok=True)
    module_spec = importlib.util.spec_from_file_location("screen_spectrum", spectrum_script)
    spectrum = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(spectrum)
    with tempfile.TemporaryDirectory(prefix="finite-screen-", dir=work) as temporary:
        root = Path(temporary)
        beam = root/"particles.h5"
        run([generator, "--output", beam, "--length-m", "1e-6", "--current-A", "0.001",
             "--central-momentum-MeV-c", "2", "--emit-x-m", "1e-12", "--emit-y-m", "1e-12",
             "--beta-x-m", "1", "--beta-y-m", "1", "--sigma-dp", "0"], root)
        for shape in ("linear-binomial3", "quintic-compensated"):
            cases = []
            for enabled in (False, True):
                case = root/(shape+("-on" if enabled else "-off"))
                case.mkdir()
                text = card.read_text().replace("@ALGORITHM@", "  longitudinal_shape: "+shape)
                text = text.replace("@PARTICLE_FILE@", str(beam)).replace("@REFERENCE_DIR@", str(case/"reference"))
                if enabled:
                    text += """
radiation_screen:
  output_file: radiation.h5
  nx: 5
  ny: 3
  x0_m: -5.0e-6
  y0_m: -5.0e-6
  dx_m: 2.5e-6
  dy_m: 5.0e-6
  z_m: 4.0e-5
  time_start_s: -1.0e-13
  time_step_s: 2.0e-16
  time_bins: 2048
"""
                path = case/"run.yaml"
                path.write_text(text)
                log = run([aprl, path], case)
                (case/"run.log").write_text(log)
                assert "Configured stop reached:" in log, log
                cases.append(case)
            off, on = cases
            for relative in ("detectors/field.h5", "detectors/particles.h5", "ledger/energy-ledger.h5"):
                compare(off/relative, on/relative)
            assert not (off/"radiation.h5").exists()
            source = on/"radiation.h5"
            with h5py.File(source) as f:
                assert int(f.attrs["complete"]) == 1
                assert int(f.attrs["accelerated_segments"]) > 0
                for key, unit in (("E_time", b"V/m"), ("B_time", b"T")):
                    assert f[key].shape == (2048, 3, 5, 3)
                    assert f[key].attrs["units"] == unit
                    assert np.isfinite(f[key][:]).all()
                    assert np.max(abs(f[key][:])) > 0
                assert all(key not in f for key in ("frequency_hz", "E_spectrum", "B_spectrum"))
            output = on/"offline.h5"
            spectrum.write_spectra(source, output, np.array([0., 2.75e12]))
            with h5py.File(source) as f, h5py.File(output) as s:
                assert int(s.attrs["complete"]) == 1
                for field in ("E", "B"):
                    dc = np.sum(f[field+"_time"][:], axis=0)*f.attrs["time_step_s"]
                    np.testing.assert_allclose(s[field+"_spectrum"][0, ..., 0], dc, rtol=2e-13, atol=1e-30)
                    np.testing.assert_array_equal(s[field+"_spectrum"][0, ..., 1], 0.)
            print(shape+": identical particle/field/physical ledger data; finite six-component screen; offline spectrum passed")


if __name__ == "__main__":
    main(*(Path(arg).resolve() for arg in sys.argv[1:]))
