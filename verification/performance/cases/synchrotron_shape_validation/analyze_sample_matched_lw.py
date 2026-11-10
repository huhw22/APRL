"""Sample-matched finite-record LW reference for independent-electron PIC runs.

The reference reads the exact x/y/arrival-time/weight records used by PIC.  A
converged regular library of single-electron complex bend-minus-straight fields
is interpolated to every particle and summed coherently before the Poynting
spectrum is formed.  No continuum bunch form factor or fitted normalization is
used for the registered sample-matched comparison.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar, mu_0
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "radiation_method_validation"))
from theory import arc_spectrum, velocity_spectrum  # noqa: E402
from theory_paired import straight_velocity_spectrum  # noqa: E402


ENERGIES = np.linspace(4., 12., 321)
FINE_STEP = 1.25e-6
COARSE_STEP = 2.5e-6
NTIME_MULTIPLIER = 8
NTHETA = 512


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def progress(message):
    print(f"[sample-matched LW] {message}", flush=True)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def poynting(ef, bf):
    return 2 * e / (h * mu_0) * np.real(
        ef[..., 0] * np.conj(bf[..., 1]) -
        ef[..., 1] * np.conj(bf[..., 0]))


def field_relative_l2(first_e, first_b, second_e, second_b):
    numerator = np.linalg.norm(first_e - second_e)**2 + \
        c**2 * np.linalg.norm(first_b - second_b)**2
    denominator = np.linalg.norm(second_e)**2 + \
        c**2 * np.linalg.norm(second_b)**2
    return float(np.sqrt(numerator / max(denominator, 1e-300)))


def load_particles(root):
    campaign = json.loads((root / "campaign.json").read_text())
    path = Path(campaign["input"])
    require(sha(path) == campaign["input_sha256"], f"Input hash changed: {path}")
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    with h5py.File(path, "r") as handle:
        records = handle["particles/records"][:]
        require(int(handle["generation/complete"][()]) == 1,
                f"Incomplete particle input: {path}")
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            f"Particle input changed while reading: {path}")
    xy = np.asarray(records["plane_position_m"], dtype=float)
    arrival = np.asarray(records["arrival_time_offset_s"], dtype=float)
    weights = np.asarray(records["macro_weight"], dtype=float)
    velocity = np.asarray(records["proper_velocity"], dtype=float)
    require(np.isfinite(xy).all() and np.isfinite(arrival).all() and
            np.isfinite(weights).all() and np.all(weights > 0),
            "Invalid particle input")
    require(np.max(np.ptp(velocity, axis=0)) < 1e-13,
            "The sample-matched library requires the common cold momentum")
    require(np.isclose(weights.sum(), campaign["represented_electrons"],
                       rtol=2e-15, atol=0.), "Represented electron mismatch")
    return {
        "path": path,
        "sha256": campaign["input_sha256"],
        "xy": xy,
        "arrival": arrival,
        "weights": weights,
        "represented": float(weights.sum()),
    }


def read_pic_axis(root):
    values = {}
    shared_time = None
    shared_detector_xy = None
    for kind in ("bend", "straight"):
        path = root / kind / "detectors/z03524.h5"
        before = (path.stat().st_size, path.stat().st_mtime_ns)
        with h5py.File(path, "r") as handle:
            group = handle["field_plane"]
            require(int(group["complete"][()]) == 1, f"Incomplete {path}")
            time = group["time_s"][:]
            if shared_time is None:
                shared_time = time
            else:
                require(np.array_equal(shared_time, time),
                        "Bend/straight detector times differ")
            nx, ny = int(group.attrs["nx"]), int(group.attrs["ny"])
            dx, dy = float(group.attrs["dx_m"]), float(group.attrs["dy_m"])
            x = float(group.attrs["x_first_m"]) + np.arange(nx) * dx
            y = float(group.attrs["y_first_m"]) + np.arange(ny) * dy
            ix = np.sort(np.argsort(abs(x))[:2])
            iy = np.sort(np.argsort(abs(y))[:2])
            require(np.max(abs(x[ix])) <= dx / 2 + 1e-18 and
                    np.max(abs(y[iy])) <= dy / 2 + 1e-18,
                    "Axis is not bracketed by the central samples")
            detector_xy = np.array([[x_value, y_value]
                                    for y_value in y[iy]
                                    for x_value in x[ix]])
            if shared_detector_xy is None:
                shared_detector_xy = detector_xy
            else:
                require(np.array_equal(shared_detector_xy, detector_xy),
                        "Bend/straight central detector coordinates differ")
            patch = group["fields"][:, iy[0]:iy[-1] + 1,
                                    ix[0]:ix[-1] + 1]
            values[kind] = {
                "E": patch["electric_V_per_m"].mean(axis=(1, 2)),
                "B": patch["magnetic_T"].mean(axis=(1, 2)),
            }
        require(before == (path.stat().st_size, path.stat().st_mtime_ns),
                f"Closed detector changed while reading: {path}")
    dt = float((shared_time[-1] - shared_time[0]) /
               (len(shared_time) - 1))
    jitter = np.max(abs(shared_time - shared_time[0] -
                        np.arange(len(shared_time)) * dt)) / dt
    require(jitter < 1e-4, "Nonuniform detector times")
    phase = np.exp(1j * ENERGIES[:, None] * e / hbar *
                   (shared_time - shared_time[0])[None, :]) * dt
    ef = phase @ (values["bend"]["E"] - values["straight"]["E"])
    bf = phase @ (values["bend"]["B"] - values["straight"]["B"])
    return {
        "E": ef,
        "B": bf,
        "spectrum": poynting(ef, bf),
        "samples": len(shared_time),
        "dt_s": dt,
        "time_span_s": len(shared_time) * dt,
        "detector_xy": shared_detector_xy,
    }


def interpolation_matrix(nodes, xy):
    step = float(nodes[1] - nodes[0])
    ix = np.floor((xy[:, 0] - nodes[0]) / step).astype(int)
    iy = np.floor((xy[:, 1] - nodes[0]) / step).astype(int)
    ix = np.where(ix == len(nodes) - 1, ix - 1, ix)
    iy = np.where(iy == len(nodes) - 1, iy - 1, iy)
    require(np.all((ix >= 0) & (ix < len(nodes) - 1) &
                   (iy >= 0) & (iy < len(nodes) - 1)),
            "A particle lies outside the LW interpolation library")
    fx = (xy[:, 0] - nodes[ix]) / step
    fy = (xy[:, 1] - nodes[iy]) / step
    rows = np.repeat(np.arange(len(xy)), 4)
    columns = np.stack([
        iy * len(nodes) + ix,
        iy * len(nodes) + ix + 1,
        (iy + 1) * len(nodes) + ix,
        (iy + 1) * len(nodes) + ix + 1,
    ], axis=1).ravel()
    data = np.stack([
        (1 - fx) * (1 - fy), fx * (1 - fy),
        (1 - fx) * fy, fx * fy,
    ], axis=1).ravel()
    matrix = coo_matrix((data, (rows, columns)),
                        shape=(len(xy), len(nodes)**2)).tocsr()
    require(np.max(abs(np.asarray(matrix.sum(axis=1)).ravel() - 1)) < 3e-15,
            "Interpolation weights do not sum to one")
    return matrix


def coherent_sum(library_e, library_b, nodes, sample, detector_xy,
                 energy_batch=16):
    matrices = [interpolation_matrix(nodes, sample["xy"] - detector)
                for detector in detector_xy]
    flat_e = library_e.reshape(len(ENERGIES), len(nodes)**2, 3)
    flat_b = library_b.reshape(len(ENERGIES), len(nodes)**2, 3)
    total_e = np.empty((len(ENERGIES), 3), complex)
    total_b = np.empty_like(total_e)
    for start in range(0, len(ENERGIES), energy_batch):
        stop = min(start + energy_batch, len(ENERGIES))
        omega = ENERGIES[start:stop] * e / hbar
        particle_phase = (np.exp(1j * omega[:, None] *
                                 sample["arrival"][None, :]) *
                          sample["weights"][None, :])
        node_coefficients = sum(
            np.asarray(matrix.T @ particle_phase.T).T
            for matrix in matrices) / len(matrices)
        total_e[start:stop] = np.einsum(
            "bn,bnd->bd", node_coefficients, flat_e[start:stop],
            optimize=True)
        total_b[start:stop] = np.einsum(
            "bn,bnd->bd", node_coefficients, flat_b[start:stop],
            optimize=True)
    return total_e, total_b


def interpolate_selected(library, nodes, xy, detector_xy):
    matrices = [interpolation_matrix(nodes, xy - detector)
                for detector in detector_xy]
    flat = library.reshape(len(ENERGIES), len(nodes)**2, 3)
    result = np.empty((len(ENERGIES), len(xy), 3), complex)
    for index in range(len(ENERGIES)):
        result[index] = sum(matrix @ flat[index]
                            for matrix in matrices) / len(matrices)
    return result


def direct_fields(xy, parameters, first, span, ntime):
    points = np.c_[-xy[:, 0], -xy[:, 1],
                   np.full(len(xy), parameters["screen"])]
    ae, ab = arc_spectrum(points, ENERGIES, parameters, ntheta=NTHETA)
    ve, vb = velocity_spectrum(points, ENERGIES, parameters,
                               time_first=first, time_span=span,
                               ntime=ntime)
    se, sb = straight_velocity_spectrum(points, ENERGIES, parameters,
                                        time_first=first, time_span=span,
                                        ntime=ntime)
    return ae + ve - se, ab + vb - sb


def direct_patch_fields(xy, detector_xy, parameters, first, span, ntime):
    fields = [direct_fields(xy - detector, parameters, first, span, ntime)
              for detector in detector_xy]
    return (sum(value[0] for value in fields) / len(fields),
            sum(value[1] for value in fields) / len(fields))


def band_metrics(pic, reference, low, high):
    selected = (ENERGIES >= low) & (ENERGIES <= high)
    return {
        "relative_L2": float(np.linalg.norm(pic[selected] -
                                             reference[selected]) /
                             np.linalg.norm(reference[selected])),
        "integral_ratio": float(np.trapz(pic[selected], ENERGIES[selected]) /
                                np.trapz(reference[selected],
                                         ENERGIES[selected])),
        "mean_pointwise_ratio": float(np.mean(pic[selected] /
                                                   reference[selected])),
    }


def main(srw_csv, output, specifications):
    srw_csv, output = Path(srw_csv).resolve(), Path(output).resolve()
    require(not output.exists(), "Output already exists")
    roots = []
    progress("reading immutable particle inputs and completed PIC detectors")
    for specification in specifications:
        require("=" in specification, "Each case must be LABEL=ROOT")
        label, text = specification.split("=", 1)
        root = Path(text).resolve()
        require((root / "PIC_COMPLETE.json").is_file(),
                f"PIC campaign is incomplete: {root}")
        roots.append({"label": label, "root": root,
                      "sample": load_particles(root),
                      "pic": read_pic_axis(root)})
    parameters = json.loads((roots[0]["root"] / "bend/PLAN.json").read_text())[
        "parameters"]
    for item in roots[1:]:
        other = json.loads((item["root"] / "bend/PLAN.json").read_text())[
            "parameters"]
        require(all(np.isclose(parameters[key], other[key], rtol=0, atol=1e-15)
                    for key in ("gamma", "beta", "rho", "length", "screen")),
                "Theory geometry differs among cases")

    maximum = max(float(np.max(abs(item["sample"]["xy"][:, None, :] -
                                   item["pic"]["detector_xy"][None, :, :])))
                  for item in roots)
    limit = max(100e-6, np.ceil((maximum + 2 * FINE_STEP) / FINE_STEP) *
                FINE_STEP)
    count = int(round(2 * limit / FINE_STEP)) + 1
    nodes = np.linspace(-limit, limit, count)
    require(len(nodes) % 2 == 1 and
            np.isclose(nodes[1] - nodes[0], FINE_STEP, rtol=0, atol=1e-18),
            "Invalid fine library nodes")
    coarse_nodes = nodes[::2]
    require(np.isclose(coarse_nodes[1] - coarse_nodes[0], COARSE_STEP,
                       rtol=0, atol=2e-18), "Invalid coarse library nodes")
    xx, yy = np.meshgrid(nodes, nodes, indexing="xy")
    points = np.c_[-xx.ravel(), -yy.ravel(),
                   np.full(xx.size, parameters["screen"])]
    acceleration_e, acceleration_b = arc_spectrum(
        points, ENERGIES, parameters, ntheta=NTHETA)
    progress(f"built acceleration library on {len(nodes)}x{len(nodes)} nodes")
    shape = (len(ENERGIES), len(nodes), len(nodes), 3)
    acceleration_e = acceleration_e.reshape(shape)
    acceleration_b = acceleration_b.reshape(shape)

    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    require(np.array_equal(ENERGIES, np.asarray(srw["photon_eV"])),
            "SRW and sample-matched energy grids differ")
    library_cache = {}
    reference_cache = {}
    rows = []
    spectra = {}
    validation = {}
    for item in roots:
        root, sample, pic = item["root"], item["sample"], item["pic"]
        plan = json.loads((root / "bend/PLAN.json").read_text())
        first = float(plan["screen_first_relative_u_s"])
        span = float(pic["time_span_s"])
        ntime = NTIME_MULTIPLIER * int(pic["samples"])
        window_key = (first, span, ntime)
        if window_key not in library_cache:
            progress(
                f"building finite-record velocity library: {ntime} time samples")
            velocity_e, velocity_b = velocity_spectrum(
                points, ENERGIES, parameters, time_first=first,
                time_span=span, ntime=ntime)
            straight_e, straight_b = straight_velocity_spectrum(
                points, ENERGIES, parameters, time_first=first,
                time_span=span, ntime=ntime)
            library_cache[window_key] = (
                acceleration_e + velocity_e.reshape(shape) -
                straight_e.reshape(shape),
                acceleration_b + velocity_b.reshape(shape) -
                straight_b.reshape(shape))
            progress("finished one finite-record velocity library")
        library_e, library_b = library_cache[window_key]
        detector_key = tuple(np.asarray(pic["detector_xy"]).ravel())
        reference_key = (sample["sha256"], window_key, detector_key)
        if reference_key not in reference_cache:
            progress(f"coherently summing exact particles for {item['label']}")
            fine_e, fine_b = coherent_sum(
                library_e, library_b, nodes, sample, pic["detector_xy"])
            coarse_e, coarse_b = coherent_sum(
                library_e[:, ::2, ::2], library_b[:, ::2, ::2],
                coarse_nodes, sample, pic["detector_xy"])
            # Direct single-particle checks use deterministic quantiles spanning
            # the actual sample, never a random validation subset.
            order = np.argsort(sample["xy"][:, 0] +
                               np.sqrt(2) * sample["xy"][:, 1])
            indices = order[np.linspace(0, len(order) - 1, 24,
                                        dtype=int)]
            direct_e, direct_b = direct_patch_fields(
                sample["xy"][indices], pic["detector_xy"], parameters,
                first, span, ntime)
            interpolated_e = interpolate_selected(
                library_e, nodes, sample["xy"][indices],
                pic["detector_xy"])
            interpolated_b = interpolate_selected(
                library_b, nodes, sample["xy"][indices],
                pic["detector_xy"])
            interpolation_error = field_relative_l2(
                interpolated_e, interpolated_b, direct_e, direct_b)

            # Test the common-window phase approximation at the actual earliest
            # and latest particle events, including their actual x/y positions.
            endpoint_errors = []
            for index in (int(np.argmin(sample["arrival"])),
                          int(np.argmax(sample["arrival"]))):
                xy_one = sample["xy"][index:index + 1]
                tau = float(sample["arrival"][index])
                base_e, base_b = direct_patch_fields(
                    xy_one, pic["detector_xy"], parameters, first, span,
                    ntime)
                shifted_e, shifted_b = direct_patch_fields(
                    xy_one, pic["detector_xy"], parameters, first - tau,
                    span, ntime)
                phase = np.exp(1j * ENERGIES * e / hbar * tau)[:, None]
                approximate_e = phase * base_e[:, 0]
                approximate_b = phase * base_b[:, 0]
                exact_e = phase * shifted_e[:, 0]
                exact_b = phase * shifted_b[:, 0]
                endpoint_errors.append(field_relative_l2(
                    approximate_e, approximate_b, exact_e, exact_b))
            reference_cache[reference_key] = {
                "E": fine_e, "B": fine_b,
                "spectrum": poynting(fine_e, fine_b),
                "coarse_spectrum": poynting(coarse_e, coarse_b),
                "coarse_field_relative_L2": field_relative_l2(
                    coarse_e, coarse_b, fine_e, fine_b),
                "direct_interpolation_field_relative_L2":
                    interpolation_error,
                "maximum_extreme_arrival_window_field_relative_L2":
                    float(max(endpoint_errors)),
            }
            progress(f"finished exact-particle reference for {item['label']}")
        reference = reference_cache[reference_key]
        continuum = (np.asarray(srw[
            "SRW_point_per_electron_J_per_eV_m2"]) *
            sample["represented"]**2 *
            np.asarray(srw["longitudinal_form_factor_squared"]) *
            np.asarray(srw["transverse_correction_interpolated"]))
        row = {
            "label": item["label"], "root": str(root),
            "particle_count": len(sample["xy"]),
            "represented_electrons": sample["represented"],
            "input_sha256": sample["sha256"],
            "detector_xy_m": pic["detector_xy"].tolist(),
            "library_time_samples": ntime,
            "coarse_vs_fine_field_relative_L2":
                reference["coarse_field_relative_L2"],
            "coarse_vs_fine_spectrum_4_12eV_relative_L2":
                band_metrics(reference["coarse_spectrum"],
                             reference["spectrum"], 4., 12.)["relative_L2"],
            "coarse_over_fine_4_12eV_integral_ratio":
                band_metrics(reference["coarse_spectrum"],
                             reference["spectrum"], 4., 12.)["integral_ratio"],
            "direct_vs_interpolated_field_relative_L2":
                reference["direct_interpolation_field_relative_L2"],
            "extreme_arrival_window_field_relative_L2":
                reference["maximum_extreme_arrival_window_field_relative_L2"],
        }
        for low, high in ((4., 8.), (8., 12.), (4., 12.)):
            prefix = f"{low:g}_{high:g}eV"
            for key, value in band_metrics(
                    pic["spectrum"], reference["spectrum"], low, high).items():
                row[f"PIC_over_sample_LW_{prefix}_{key}"] = value
            for key, value in band_metrics(
                    reference["spectrum"], continuum, low, high).items():
                row[f"sample_LW_over_continuum_{prefix}_{key}"] = value
        rows.append(row)
        spectra[item["label"]] = {
            "represented": sample["represented"],
            "pic": pic["spectrum"],
            "sample_lw": reference["spectrum"],
            "continuum": continuum,
        }
        validation[item["label"]] = {
            key: row[key] for key in row if
            key.endswith("field_relative_L2")
        }

    output.mkdir(parents=True)
    progress("writing tables and comparison figure")
    with (output / "sample_matched_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output / "sample_matched_spectra.csv").open("w", newline="") as stream:
        names = ["photon_eV"]
        columns = [ENERGIES]
        for label in spectra:
            names += [f"{label}_PIC_J_per_eV_m2",
                      f"{label}_sample_LW_J_per_eV_m2",
                      f"{label}_continuum_J_per_eV_m2"]
            columns += [spectra[label]["pic"], spectra[label]["sample_lw"],
                        spectra[label]["continuum"]]
        writer = csv.writer(stream)
        writer.writerow(names)
        writer.writerows(zip(*columns))

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5),
                             layout="constrained")
    colors = plt.cm.tab10(np.linspace(0, .8, len(rows)))
    first_curve = next(iter(spectra.values()))
    axes[0, 0].plot(ENERGIES,
                    first_curve["continuum"] /
                    first_curve["represented"]**2 * 1e16,
                    "k-", lw=2.2, label=r"continuous SRW / $N^2$")
    for row, color in zip(rows, colors):
        curve = spectra[row["label"]]
        scale = curve["represented"]**2
        axes[0, 0].plot(ENERGIES, curve["sample_lw"] / scale * 1e16,
                        color=color, lw=2, label=row["label"] + " sample LW")
        axes[0, 0].plot(ENERGIES, curve["pic"] / scale * 1e16,
                        color=color, lw=1.2, ls="--",
                        label=row["label"] + " PIC")
        axes[0, 1].plot(ENERGIES, curve["pic"] / curve["sample_lw"],
                        color=color, label=row["label"])
        axes[1, 0].plot(ENERGIES, curve["sample_lw"] / curve["continuum"],
                        color=color, label=row["label"])
    axes[0, 0].set(title="Exact finite sample: coherent scale",
                   ylabel=r"$10^{-16}$ J / (eV m$^2$ electron$^2$)")
    axes[0, 1].axhline(1, color="black", lw=1)
    axes[0, 1].set(title="Raw APRL / same-particle LW", ylabel="ratio")
    axes[1, 0].axhline(1, color="black", lw=1)
    axes[1, 0].set(title="Same-particle LW / continuous Gaussian SRW",
                   ylabel="ratio")
    indices = np.arange(len(rows))
    width = .36
    axes[1, 1].bar(indices - width / 2,
                   [row["PIC_over_sample_LW_4_8eV_integral_ratio"]
                    for row in rows], width, label="4–8 eV")
    axes[1, 1].bar(indices + width / 2,
                   [row["PIC_over_sample_LW_8_12eV_integral_ratio"]
                    for row in rows], width, label="8–12 eV")
    axes[1, 1].axhline(1, color="black", lw=1)
    axes[1, 1].set(title="Signed band integral: APRL / same-particle LW",
                   ylabel="ratio", xticks=indices,
                   xticklabels=[row["label"] for row in rows])
    axes[1, 1].tick_params(axis="x", labelrotation=20)
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xlabel("photon energy [eV]")
        axis.legend(fontsize=6 if axis is axes[0, 0] else 8)
    fig.suptitle("Sample-matched finite-record LW control; exact PIC particle realization")
    fig.savefig(output / "sample_matched_lw_comparison.png", dpi=200)
    fig.savefig(output / "sample_matched_lw_comparison.pdf")
    plt.close(fig)

    final = {
        "complete": True,
        "sample_matched": True,
        "particle_fields_summed_coherently_before_poynting": True,
        "continuum_form_factor_used_for_registered_comparison": False,
        "angular_filter_applied": False,
        "aperture_average_applied": False,
        "fitted_normalization_applied": False,
        "library": {
            "fine_step_m": FINE_STEP,
            "coarse_step_m": COARSE_STEP,
            "transverse_limit_m": limit,
            "fine_nodes_per_axis": len(nodes),
            "ntheta": NTHETA,
            "time_quadrature_multiplier_over_PIC_samples": NTIME_MULTIPLIER,
        },
        "validation": validation,
        "rows": rows,
        "physical_accuracy_accepted": False,
        "acceptance_note": (
            "This removes finite-particle realization bias from the reference. "
            "Physical acceptance still requires judging the resulting PIC/LW error "
            "and its coupled particle/grid convergence."),
    }
    (output / "FINAL.json").write_text(json.dumps(final, indent=2) + "\n")
    progress(f"complete: {output}")
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) >= 4,
            "usage: analyze_sample_matched_lw.py SRW_CSV OUTPUT LABEL=ROOT ...")
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
