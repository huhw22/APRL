"""Compare 8 eV transverse radiation fields for two exact 100k samples.

The PIC bend-minus-straight Fourier field contains the outgoing-bunch velocity
near field.  The registered optical comparison subtracts a same-particle,
finite-record, first-order transverse-deposition proxy for that velocity term
at the complex-field level.  It then compares the result with the coherent
acceleration LW field of the exact particle realization.  The uncorrected full
paired maps are retained separately as a diagnostic, not called optical light.
"""
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
from scipy.signal import fftconvolve
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "radiation_method_validation"))
from theory import arc_spectrum, velocity_spectrum  # noqa: E402
from theory_paired import straight_velocity_spectrum  # noqa: E402
from analyze_sample_matched_lw import load_particles  # noqa: E402


ENERGY_EV = 8.0
FINE_STEP = 1.25e-6
COARSE_STEP = 2.5e-6
NTIME_MULTIPLIER = 8
NTHETA = 512
CPML_CELLS = 10


def progress(message):
    print(f"[8 eV transverse] {message}", flush=True)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def interpolation_matrix(x, y, xy):
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    ix = np.floor((xy[:, 0] - x[0]) / dx).astype(int)
    iy = np.floor((xy[:, 1] - y[0]) / dy).astype(int)
    ix = np.where(ix == len(x) - 1, ix - 1, ix)
    iy = np.where(iy == len(y) - 1, iy - 1, iy)
    require(np.all((ix >= 0) & (ix < len(x) - 1) &
                   (iy >= 0) & (iy < len(y) - 1)),
            "Particle lies outside transverse interpolation grid")
    fx = (xy[:, 0] - x[ix]) / dx
    fy = (xy[:, 1] - y[iy]) / dy
    rows = np.repeat(np.arange(len(xy)), 4)
    columns = np.stack([
        iy * len(x) + ix, iy * len(x) + ix + 1,
        (iy + 1) * len(x) + ix, (iy + 1) * len(x) + ix + 1,
    ], axis=1).ravel()
    values = np.stack([
        (1 - fx) * (1 - fy), fx * (1 - fy),
        (1 - fx) * fy, fx * fy,
    ], axis=1).ravel()
    matrix = coo_matrix((values, (rows, columns)),
                        shape=(len(xy), len(x) * len(y))).tocsr()
    require(np.max(abs(np.asarray(matrix.sum(axis=1)).ravel() - 1)) < 4e-15,
            "Transverse interpolation weights do not sum to one")
    return matrix


def complex_source(x, y, sample):
    matrix = interpolation_matrix(x, y, sample["xy"])
    omega = ENERGY_EV * e / hbar
    amplitude = sample["weights"] * np.exp(1j * omega * sample["arrival"])
    return np.asarray(matrix.T @ amplitude).reshape(len(y), len(x))


def convolve_same_nodes(source, kernel):
    """Return sum_s source[s] kernel[o-s] on an identical regular grid."""
    ny, nx = source.shape
    require(kernel.shape[:2] == (2 * ny - 1, 2 * nx - 1),
            "Invalid same-node convolution kernel")
    result = np.empty((ny, nx, kernel.shape[2]), complex)
    for component in range(kernel.shape[2]):
        full = fftconvolve(source, kernel[:, :, component], mode="full")
        result[:, :, component] = full[ny - 1:2 * ny - 1,
                                                nx - 1:2 * nx - 1]
    return result


def convolve_vertices_to_centres(source, kernel, ny, nx):
    """Map (ny+1,nx+1) vertex sources to ny,nx cell-centre samples."""
    require(source.shape == (ny + 1, nx + 1) and
            kernel.shape[:2] == (2 * ny, 2 * nx),
            "Invalid vertex-to-centre convolution geometry")
    result = np.empty((ny, nx, kernel.shape[2]), complex)
    for component in range(kernel.shape[2]):
        full = fftconvolve(source, kernel[:, :, component], mode="full")
        result[:, :, component] = full[ny:2 * ny, nx:2 * nx]
    return result


def convolution_self_check():
    rng = np.random.default_rng(731)
    source = rng.normal(size=(4, 5)) + 1j * rng.normal(size=(4, 5))
    kernel = rng.normal(size=(7, 9, 1)) + 1j * rng.normal(size=(7, 9, 1))
    got = convolve_same_nodes(source, kernel)[:, :, 0]
    expected = np.zeros_like(source)
    for oy in range(4):
        for ox in range(5):
            for sy in range(4):
                for sx in range(5):
                    expected[oy, ox] += source[sy, sx] * kernel[
                        oy - sy + 3, ox - sx + 4, 0]
    first = np.max(abs(got - expected))
    vertices = rng.normal(size=(5, 6)) + 1j * rng.normal(size=(5, 6))
    half_kernel = rng.normal(size=(8, 10, 1)) + \
        1j * rng.normal(size=(8, 10, 1))
    got_half = convolve_vertices_to_centres(vertices, half_kernel, 4, 5)[:, :, 0]
    expected_half = np.zeros((4, 5), complex)
    for oy in range(4):
        for ox in range(5):
            for sy in range(5):
                for sx in range(6):
                    expected_half[oy, ox] += vertices[sy, sx] * half_kernel[
                        oy - sy + 4, ox - sx + 5, 0]
    second = np.max(abs(got_half - expected_half))
    require(first < 2e-13 and second < 2e-13,
            "FFT convolution indexing self-check failed")
    return {"same_node_max_abs_error": float(first),
            "vertex_to_centre_max_abs_error": float(second)}


def load_pic_field(root):
    result = None
    shared = None
    for kind, sign in (("bend", 1.), ("straight", -1.)):
        path = root / kind / "detectors/z03524.h5"
        before = (path.stat().st_size, path.stat().st_mtime_ns)
        with h5py.File(path, "r") as handle:
            group = handle["field_plane"]
            require(int(group["complete"][()]) == 1, f"Incomplete {path}")
            time = group["time_s"][:]
            nx, ny = int(group.attrs["nx"]), int(group.attrs["ny"])
            dx, dy = float(group.attrs["dx_m"]), float(group.attrs["dy_m"])
            x = float(group.attrs["x_first_m"]) + np.arange(nx) * dx
            y = float(group.attrs["y_first_m"]) + np.arange(ny) * dy
            if shared is None:
                shared = (time, x, y, dx, dy)
                dt = float((time[-1] - time[0]) / (len(time) - 1))
                phase = np.exp(1j * ENERGY_EV * e / hbar *
                               (time - time[0])) * dt
                result = {"E": np.zeros((ny, nx, 3), complex),
                          "B": np.zeros((ny, nx, 3), complex)}
            else:
                require(np.array_equal(shared[0], time) and
                        np.array_equal(shared[1], x) and
                        np.array_equal(shared[2], y),
                        "Bend and straight detector grids differ")
            fields = group["fields"]
            for key, name in (("E", "electric_V_per_m"),
                              ("B", "magnetic_T")):
                values = fields[name][:]
                result[key] += sign * np.einsum(
                    "t,tyxd->yxd", phase, values, optimize=True)
                del values
        require(before == (path.stat().st_size, path.stat().st_mtime_ns),
                f"Closed detector changed while reading: {path}")
    time, x, y, dx, dy = shared
    result.update({"time": time, "x": x, "y": y, "dx": dx, "dy": dy,
                   "dt": dt, "span": len(time) * dt,
                   "samples": len(time)})
    return result


def acceleration_kernel(xfine, yfine, parameters):
    rx = np.arange(-(len(xfine) - 1), len(xfine)) * FINE_STEP
    ry = np.arange(-(len(yfine) - 1), len(yfine)) * FINE_STEP
    xx, yy = np.meshgrid(rx, ry, indexing="xy")
    points = np.c_[xx.ravel(), yy.ravel(),
                   np.full(xx.size, parameters["screen"])]
    ef, bf = arc_spectrum(points, np.array([ENERGY_EV]), parameters,
                          ntheta=NTHETA)
    shape = (len(ry), len(rx), 3)
    return ef[0].reshape(shape), bf[0].reshape(shape)


def velocity_kernel(nx, ny, dx, dy, parameters, first, span, ntime):
    rx = (np.arange(-nx, nx) + .5) * dx
    ry = (np.arange(-ny, ny) + .5) * dy
    xx, yy = np.meshgrid(rx, ry, indexing="xy")
    points = np.c_[xx.ravel(), yy.ravel(),
                   np.full(xx.size, parameters["screen"])]
    ve, vb = velocity_spectrum(
        points, np.array([ENERGY_EV]), parameters, time_first=first,
        time_span=span, ntime=ntime)
    se, sb = straight_velocity_spectrum(
        points, np.array([ENERGY_EV]), parameters, time_first=first,
        time_span=span, ntime=ntime)
    shape = (2 * ny, 2 * nx, 3)
    return ((ve - se)[0].reshape(shape),
            (vb - sb)[0].reshape(shape))


def spectral_intensity(ef, bf):
    transverse = (abs(ef[:, :, 0])**2 + abs(ef[:, :, 1])**2 +
                  c**2 * (abs(bf[:, :, 0])**2 + abs(bf[:, :, 1])**2))
    return e / (h * mu_0 * c) * transverse


def spectral_poynting(ef, bf):
    return 2 * e / (h * mu_0) * np.real(
        ef[:, :, 0] * np.conj(bf[:, :, 1]) -
        ef[:, :, 1] * np.conj(bf[:, :, 0]))


def field_l2(first_e, first_b, second_e, second_b, region):
    first_e, first_b = first_e[region], first_b[region]
    second_e, second_b = second_e[region], second_b[region]
    numerator = np.linalg.norm(first_e - second_e)**2 + \
        c**2 * np.linalg.norm(first_b - second_b)**2
    denominator = np.linalg.norm(second_e)**2 + \
        c**2 * np.linalg.norm(second_b)**2
    return float(np.sqrt(numerator / max(denominator, 1e-300)))


def map_metrics(pic_i, theory_i, pic_e, pic_b, theory_e, theory_b,
                x, y, dx, dy):
    region = np.s_[CPML_CELLS:-CPML_CELLS, CPML_CELLS:-CPML_CELLS]
    p, q = pic_i[region], theory_i[region]
    area = dx * dy
    ptotal, qtotal = float(p.sum() * area), float(q.sum() * area)
    cosine = float(np.vdot(p.ravel(), q.ravel()).real /
                   (np.linalg.norm(p) * np.linalg.norm(q)))
    pshape = p / p.sum()
    qshape = q / q.sum()
    shape_l2 = float(np.linalg.norm(pshape - qshape) /
                     np.linalg.norm(qshape))
    xx, yy = np.meshgrid(x[CPML_CELLS:-CPML_CELLS],
                         y[CPML_CELLS:-CPML_CELLS], indexing="xy")
    def moments(values):
        total = values.sum()
        cx, cy = float((values * xx).sum() / total), \
            float((values * yy).sum() / total)
        sx = float(np.sqrt((values * (xx - cx)**2).sum() / total))
        sy = float(np.sqrt((values * (yy - cy)**2).sum() / total))
        return [cx, cy], [sx, sy]
    pc, ps = moments(p)
    qc, qs = moments(q)
    return {
        "PIC_over_LW_integrated_intensity": ptotal / qtotal,
        "intensity_cosine_similarity": cosine,
        "unit_integral_shape_relative_L2": shape_l2,
        "complex_field_relative_L2": field_l2(
            pic_e, pic_b, theory_e, theory_b, region),
        "PIC_integrated_J_per_eV": ptotal,
        "LW_integrated_J_per_eV": qtotal,
        "PIC_centroid_m": pc, "LW_centroid_m": qc,
        "PIC_rms_size_m": ps, "LW_rms_size_m": qs,
    }


def plot_maps(cases, output, suffix, title, pic_key, theory_key):
    fig, axes = plt.subplots(len(cases), 3, figsize=(14.5, 7.5),
                             layout="constrained", squeeze=False)
    for row, case in enumerate(cases):
        pic = case[pic_key] * 1e6
        theory = case[theory_key] * 1e6
        difference = (pic - theory) / max(float(theory.max()), 1e-300)
        vmax = max(float(pic.max()), float(theory.max()))
        norm = colors.PowerNorm(gamma=.45, vmin=0., vmax=vmax)
        extent = [case["x"][0] * 1e6, case["x"][-1] * 1e6,
                  case["y"][0] * 1e6, case["y"][-1] * 1e6]
        first = axes[row, 0].imshow(theory, origin="lower", extent=extent,
                                    aspect="equal", norm=norm, cmap="magma")
        axes[row, 1].imshow(pic, origin="lower", extent=extent,
                            aspect="equal", norm=norm, cmap="magma")
        limit = max(float(abs(difference).max()), 1e-12)
        third = axes[row, 2].imshow(
            difference, origin="lower", extent=extent, aspect="equal",
            cmap="coolwarm", vmin=-limit, vmax=limit)
        for column in (0, 1):
            fig.colorbar(first, ax=axes[row, column], shrink=.82,
                         label=r"$\mu$J / (eV m$^2$)")
        fig.colorbar(third, ax=axes[row, 2], shrink=.82,
                     label="(PIC - LW) / LW peak")
        axes[row, 0].set_ylabel(case["label"] + "\ny [µm]")
        for axis in axes[row]:
            axis.set_xlabel("x [µm]")
    axes[0, 0].set_title("coherent same-particle LW")
    axes[0, 1].set_title("APRL PIC")
    axes[0, 2].set_title("signed difference")
    fig.suptitle(title)
    fig.savefig(output / f"{suffix}.png", dpi=210)
    fig.savefig(output / f"{suffix}.pdf")
    plt.close(fig)


def main(output, specifications):
    output = Path(output).resolve()
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    checks = convolution_self_check()
    cases = []
    parameters = None
    progress("reading exact particles and completed full transverse PIC planes")
    for specification in specifications:
        label, text = specification.split("=", 1)
        root = Path(text).resolve()
        sample = load_particles(root)
        pic = load_pic_field(root)
        plan = json.loads((root / "bend/PLAN.json").read_text())
        current = plan["parameters"]
        if parameters is None:
            parameters = current
        else:
            require(all(np.isclose(parameters[key], current[key], rtol=0,
                                   atol=1e-15)
                        for key in ("gamma", "beta", "rho", "length",
                                    "screen")),
                    "Theory geometry differs among cases")
        cases.append({"label": label, "root": root, "sample": sample,
                      "pic": pic, "first": float(
                          plan["screen_first_relative_u_s"])})
    require(len(cases) == 2, "Exactly two 100k cases are required")
    reference = cases[0]["pic"]
    for case in cases:
        pic = case["pic"]
        require(len(case["sample"]["xy"]) == 100000,
                "This registered comparison requires 100k physical electrons")
        require(np.array_equal(pic["x"], reference["x"]) and
                np.array_equal(pic["y"], reference["y"]) and
                np.isclose(pic["dx"], 5e-6, rtol=0, atol=1e-18) and
                np.isclose(pic["dy"], 5e-6, rtol=0, atol=1e-18),
                "Unexpected transverse detector geometry")

    x, y = reference["x"], reference["y"]
    dx, dy = reference["dx"], reference["dy"]
    factor = int(round(dx / FINE_STEP))
    require(factor == 4 and np.isclose(dy / FINE_STEP, factor),
            "Fine theory grid must nest the PIC detector grid")
    xfine = x[0] + np.arange((len(x) - 1) * factor + 1) * FINE_STEP
    yfine = y[0] + np.arange((len(y) - 1) * factor + 1) * FINE_STEP
    progress(f"building one-electron acceleration kernel on "
             f"{2 * len(xfine) - 1}x{2 * len(yfine) - 1} displacements")
    acceleration_e, acceleration_b = acceleration_kernel(
        xfine, yfine, parameters)
    progress("finished acceleration kernel")

    results = []
    omega = ENERGY_EV * e / hbar
    for case in cases:
        label, sample, pic = case["label"], case["sample"], case["pic"]
        progress(f"coherently summing acceleration field for {label}")
        source_fine = complex_source(xfine, yfine, sample)
        radiation_e_fine = convolve_same_nodes(source_fine, acceleration_e)
        radiation_b_fine = convolve_same_nodes(source_fine, acceleration_b)
        radiation_e = radiation_e_fine[::factor, ::factor]
        radiation_b = radiation_b_fine[::factor, ::factor]
        require(radiation_e.shape[:2] == (len(y), len(x)),
                "Fine-grid radiation sampling failed")

        # Nested 2.5 um convergence check for the smooth acceleration field.
        coarse_factor = int(round(COARSE_STEP / FINE_STEP))
        xc, yc = xfine[::coarse_factor], yfine[::coarse_factor]
        source_coarse = complex_source(xc, yc, sample)
        coarse_e = convolve_same_nodes(
            source_coarse, acceleration_e[::coarse_factor,
                                          ::coarse_factor])
        coarse_b = convolve_same_nodes(
            source_coarse, acceleration_b[::coarse_factor,
                                          ::coarse_factor])
        sample_factor = int(round(dx / COARSE_STEP))
        coarse_e = coarse_e[::sample_factor, ::sample_factor]
        coarse_b = coarse_b[::sample_factor, ::sample_factor]

        progress(f"building finite-record deposited velocity background for {label}")
        xv = x[0] - dx / 2 + np.arange(len(x) + 1) * dx
        yv = y[0] - dy / 2 + np.arange(len(y) + 1) * dy
        source_vertex = complex_source(xv, yv, sample)
        velocity_e_kernel, velocity_b_kernel = velocity_kernel(
            len(x), len(y), dx, dy, parameters, case["first"], pic["span"],
            NTIME_MULTIPLIER * pic["samples"])
        velocity_e = convolve_vertices_to_centres(
            source_vertex, velocity_e_kernel, len(y), len(x))
        velocity_b = convolve_vertices_to_centres(
            source_vertex, velocity_b_kernel, len(y), len(x))

        # Convert the theory's absolute u=t-z/c Fourier phase to the detector
        # convention exp(i omega (t-t_first)).  B is stored at t-dt/2 but
        # transformed with the E timestamp, hence the additional half-step.
        electric_phase = np.exp(-1j * omega * case["first"])
        magnetic_phase = electric_phase * np.exp(1j * omega * pic["dt"] / 2)
        radiation_e *= electric_phase
        radiation_b *= magnetic_phase
        coarse_e *= electric_phase
        coarse_b *= magnetic_phase
        velocity_e *= electric_phase
        velocity_b *= magnetic_phase

        pic_radiation_e = pic["E"] - velocity_e
        pic_radiation_b = pic["B"] - velocity_b
        theory_full_e = radiation_e + velocity_e
        theory_full_b = radiation_b + velocity_b
        radiation_pic_i = spectral_intensity(pic_radiation_e,
                                             pic_radiation_b)
        radiation_lw_i = spectral_intensity(radiation_e, radiation_b)
        raw_pic_i = spectral_intensity(pic["E"], pic["B"])
        full_lw_i = spectral_intensity(theory_full_e, theory_full_b)
        metrics = map_metrics(
            radiation_pic_i, radiation_lw_i, pic_radiation_e,
            pic_radiation_b, radiation_e, radiation_b, x, y, dx, dy)
        region = np.s_[CPML_CELLS:-CPML_CELLS,
                       CPML_CELLS:-CPML_CELLS]
        metrics.update({
            "fine_vs_2p5um_acceleration_field_relative_L2": field_l2(
                coarse_e, coarse_b, radiation_e, radiation_b, region),
            "raw_full_complex_field_relative_L2": field_l2(
                pic["E"], pic["B"], theory_full_e, theory_full_b, region),
            "raw_velocity_background_peak_over_radiation_peak": float(
                spectral_intensity(velocity_e, velocity_b)[region].max() /
                radiation_lw_i[region].max()),
            "radiation_PIC_negative_poynting_pixel_fraction": float(np.mean(
                spectral_poynting(pic_radiation_e, pic_radiation_b)[region] < 0)),
            "radiation_LW_negative_poynting_pixel_fraction": float(np.mean(
                spectral_poynting(radiation_e, radiation_b)[region] < 0)),
        })
        results.append({
            "label": label, "root": str(case["root"]), "x": x, "y": y,
            "radiation_PIC": radiation_pic_i,
            "radiation_LW": radiation_lw_i,
            "raw_PIC": raw_pic_i, "full_LW": full_lw_i,
            "metrics": metrics,
        })
        progress(f"finished {label}")

    plot_maps(
        results, output, "transverse_radiation_8eV",
        "8 eV transverse synchrotron radiation; velocity background removed before intensity",
        "radiation_PIC", "radiation_LW")
    plot_maps(
        results, output, "transverse_raw_paired_8eV",
        "8 eV raw bend-minus-straight field; velocity near field retained (diagnostic only)",
        "raw_PIC", "full_LW")
    np.savez_compressed(
        output / "transverse_maps_8eV.npz", x_m=x, y_m=y,
        **{f"{case['label']}_{key}_J_per_eV_m2": case[key]
           for case in results for key in
           ("radiation_PIC", "radiation_LW", "raw_PIC", "full_LW")})
    with (output / "transverse_metrics.csv").open("w", newline="") as stream:
        rows = [{"label": case["label"], **case["metrics"]}
                for case in results]
        fieldnames = ["label"] + [key for key, value in rows[0].items()
                                  if key != "label" and
                                  not isinstance(value, list)]
        writer = csv.DictWriter(stream, fieldnames=fieldnames,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    final = {
        "complete": True,
        "photon_energy_eV": ENERGY_EV,
        "single_frequency": True,
        "frequency_band_integrated": False,
        "low_frequency_included": False,
        "registered_observable": (
            "positive transverse spectral field intensity at 8 eV after "
            "complex-field subtraction of the same-particle finite-record "
            "deposited velocity background"),
        "theory": (
            "coherent acceleration LW field of the exact particle x/y, "
            "arrival time and weight; 1.25 um transverse convolution grid"),
        "raw_full_map_is_optical_acceptance": False,
        "CPML_cells_excluded_from_metrics": CPML_CELLS,
        "convolution_self_check": checks,
        "rows": [{"label": case["label"], "root": case["root"],
                  **case["metrics"]} for case in results],
        "physical_accuracy_accepted": False,
        "acceptance_note": (
            "Acceptance requires reviewing the generated maps and the nested "
            "acceleration-grid convergence. The vertex-deposited velocity "
            "background is a transverse first-order proxy, not an exact "
            "component-staggered Cowan discrete Green function."),
    }
    (output / "FINAL.json").write_text(json.dumps(final, indent=2) + "\n")
    hashes = {path.name: sha256(path) for path in output.iterdir()
              if path.is_file() and path.name != "SOURCE_AUDIT.json"}
    (output / "SOURCE_AUDIT.json").write_text(json.dumps({
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "outputs_sha256": hashes,
    }, indent=2) + "\n")
    progress(f"complete: {output}")
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: analyze_100k_transverse_single_frequency.py OUTPUT "
            "LABEL1=ROOT1 LABEL2=ROOT2")
    main(sys.argv[1], sys.argv[2:])
