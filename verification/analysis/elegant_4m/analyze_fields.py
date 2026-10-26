#!/usr/bin/env python3
"""Stream APRL laboratory field planes and audit the separate boosted ledger.

This reads existing results only; it does not start an APRL simulation. Full
field scans belong on an allocated compute node. Only committed, complete
datasets are accepted. Raw fields are never relabelled as pure radiation.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import tempfile

import h5py
import numpy as np

C = 299792458.0
MU0 = 1.25663706212e-6
EPS0 = 1.0 / (MU0 * C * C)
QE = 1.602176634e-19
FACE_NAMES = ("x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus")
NOTES = [
    "Raw laboratory Maxwell fields include particle near/bound fields; no background was subtracted.",
    "Poynting flux uses stored E/B. Yee B has a half-step time stagger; this is not an exact simultaneous-field invariant.",
    "Lorentz-transformed laboratory E and B both mix staggered boosted components. Shifting only laboratory B cannot exactly undo the stagger.",
    "epsilon0*c*(Ex^2+Ey^2) is an instantaneous forward plane-wave proxy, not a measurement of separated radiation.",
    "Positive/negative flux split is pointwise max/min(Sz,0) before area integration; it is not a modal forward/backward wave decomposition.",
    "All integrals cover only the recorded time window and finite aperture. A complete HDF5 file does not prove that the pulse tail was recorded.",
    "The full recorded aperture includes CPML cells; the physical aperture excludes configured transverse CPML layers.",
    "The central waveform averages the nearest one or four cell centres; it is not a directly sampled geometric-axis value on an even grid.",
]


def plain(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def dump_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               default=plain, allow_nan=False) + "\n")


def trapezoid(y, x, axis=0):
    fn = getattr(np, "trapezoid", np.trapz)
    return fn(y, x=x, axis=axis)


def committed(group, dataset, marker):
    if int(group["complete"][()]) != 1:
        raise ValueError(f"Incomplete input: {group.file.filename}:{group.name}")
    n = int(group[marker][()])
    if n <= 0 or n > group[dataset].shape[0]:
        raise ValueError("Invalid committed extent")
    return n


def save_csv(path, columns):
    names = list(columns)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(zip(*(columns[n] for n in names)))


def save_figure(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=170, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.grid": True,
                         "grid.alpha": .22, "figure.dpi": 120})
    return plt


def analyze_ledger(run_dir, output_dir, plots=True):
    path = run_dir / "ledger/energy-ledger.h5"
    with h5py.File(path, "r") as handle:
        group = handle["energy_ledger"]
        n = committed(group, "records", "committed_records")
        records = group["records"][:n]
        attrs = {k: plain(v) for k, v in group.attrs.items()}
    columns = {}
    for key in records.dtype.names:
        value = records[key]
        if value.ndim == 1:
            columns[key] = value
        else:
            for index in range(value.shape[1]):
                suffix = FACE_NAMES[index] if key == "outward_field_energy_face_J" else str(index)
                columns[f"{key}_{suffix}"] = value[:, index]
    if not all(np.isfinite(v).all() for v in columns.values()):
        raise ValueError("Non-finite energy ledger values")
    t = records["time_box_s"]
    if n > 1 and not np.all(np.diff(t) > 0):
        raise ValueError("Ledger time must be strictly increasing")
    accounted = records["particle_kinetic_active_J"] + records["particle_kinetic_removed_J"]
    delta_k = accounted - accounted[0]
    delta_u = records["field_energy_interior_J"] - records["field_energy_interior_J"][0]
    outward = records["outward_field_energy_J"]
    prescribed = records["prescribed_work_J"]
    boundary = records["boundary_source_work_J"] if "boundary_source_work_J" in records.dtype.names else np.zeros(n)
    residual_raw = delta_k + delta_u + outward - prescribed
    residual = residual_raw - boundary
    error = residual - records["balance_residual_J"]
    face_error = np.sum(records["outward_field_energy_face_J"], axis=1) - outward
    scale = max(float(np.max(np.abs(accounted))), float(np.max(np.abs(prescribed))), 1e-300)
    atol = 256 * np.finfo(float).eps * scale
    if np.max(np.abs(error)) > atol or np.max(np.abs(face_error)) > atol:
        raise ValueError("Ledger algebraic identity disagrees with stored columns")
    columns.update(delta_accounted_kinetic_J=delta_k, delta_interior_field_J=delta_u,
                   recomputed_raw_residual_J=residual_raw, recomputed_residual_J=residual,
                   residual_identity_difference_J=error)
    save_csv(output_dir / "boost_energy_ledger.csv", columns)
    np.savez_compressed(output_dir / "boost_energy_ledger.npz", **columns)
    electrons = float(records["active_represented_electrons"][0])
    final = {k: plain(records[k][-1]) for k in records.dtype.names}
    summary = {
        "source_file": str(path.resolve()), "frame": "boosted computational frame",
        "committed_records": n, "attributes": attrs,
        "closure_valid_all": bool(np.all(records["closure_valid"] == 1)),
        "initial": {k: plain(records[k][0]) for k in records.dtype.names},
        "final": final,
        "delta_accounted_kinetic_J": float(delta_k[-1]),
        "particle_kinetic_loss_J": float(-delta_k[-1]),
        "particle_kinetic_loss_eV_per_initial_represented_electron": float(-delta_k[-1] / electrons / QE),
        "delta_interior_field_J": float(delta_u[-1]),
        "sum_field_change_and_outflow_J": float(delta_u[-1] + outward[-1]),
        "external_work_minus_kinetic_change_J": float(prescribed[-1] + boundary[-1] - delta_k[-1]),
        "residual_identity_max_abs_error_J": float(np.max(np.abs(error))),
        "face_sum_max_abs_error_J": float(np.max(np.abs(face_error))),
        "identity_tolerance_J": atol,
        "transverse_fraction_of_net_outward_energy": float(np.sum(records["outward_field_energy_face_J"][-1, :4]) / outward[-1]) if outward[-1] else None,
        "notes": [
            "Delta(K_active+K_removed)+Delta(U_interior)+E_out-W_prescribed-W_boundary=R.",
            "U_interior excludes CPML and includes all grid field energy, not only particle-bound fields.",
            "E_out is cumulative signed flux through six inner CPML surfaces, accumulated every Maxwell step.",
            "Prescribed magnetic devices have an electric component and can do work in the boosted frame.",
            "closure_valid=0 forbids interpreting a small algebraic residual as a whole-bunch conservation proof.",
            "Lab momentum moments in this ledger are on equal boosted-time slices, not the fixed-location particle diagnostics.",
            "Do not combine boosted ledger energy with laboratory detector-plane energy in a scalar conservation equation.",
        ],
    }
    dump_json(output_dir / "boost_energy_summary.json", summary)
    if plots:
        plt = plotting()
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        tx = t * 1e12
        axes[0, 0].plot(tx, delta_k * 1e9, label="Particle kinetic change")
        axes[0, 0].plot(tx, prescribed * 1e9, label="Prescribed-device work")
        axes[0, 0].set_ylabel("Energy [nJ]")
        axes[0, 1].plot(tx, delta_u * 1e9, label="Interior field change")
        axes[0, 1].plot(tx, outward * 1e9, label="Cumulative net field outflow")
        axes[0, 1].plot(tx, (prescribed + boundary - delta_k) * 1e9,
                        label="Device + boundary work - kinetic change", lw=1)
        axes[0, 1].set_ylabel("Energy [nJ]")
        for i, name in enumerate(FACE_NAMES):
            axes[1, 0].plot(tx, records["outward_field_energy_face_J"][:, i] * 1e9,
                            label=name.replace("minus", "-").replace("plus", "+"))
        axes[1, 0].set_ylabel("Cumulative signed field outflow [nJ]")
        axes[1, 1].plot(tx, residual_raw * 1e9, label="Raw residual")
        axes[1, 1].plot(tx, residual * 1e9, label="Boundary-source-corrected residual", ls="--")
        axes[1, 1].plot(tx, boundary * 1e9, label="Boundary source work")
        axes[1, 1].set_ylabel("Energy [nJ]")
        for ax in axes.flat:
            ax.set_xlabel("Boosted computational time [ps]")
            ax.legend(fontsize=8)
        fig.suptitle("Boosted-frame energy accounting | approximate boundary injection; closure_valid=0"
                     if not summary["closure_valid_all"] else "Boosted-frame energy accounting")
        save_figure(fig, output_dir / "boost_energy_accounting")
    return summary


def cpml_layers(run_dir):
    import yaml
    path = run_dir / "run.yaml"
    if not path.exists():
        return 0, 0, "run.yaml absent; physical aperture assumed equal to recorded aperture"
    config = yaml.safe_load(path.read_text())
    boundary = config.get("boundary", {})
    if boundary.get("type") != "cpml":
        return 0, 0, "no CPML exclusion"
    cells = boundary.get("cells", [0, 0, 0])
    return int(cells[0]), int(cells[1]), "transverse CPML excluded using run.yaml"


def analyze_field(run_dir, output_dir, block_frames=8, plots=True, field_file=None):
    path = field_file if field_file is not None else run_dir / "detectors/exit-field.h5"
    cx, cy, aperture_note = cpml_layers(run_dir)
    with h5py.File(path, "r") as handle:
        group = handle["field_plane"]
        n = committed(group, "fields", "committed_samples")
        attrs = {k: plain(v) for k, v in group.attrs.items()}
        data = group["fields"]
        times = group["time_s"][:n]
        if n < 2 or not np.isfinite(times).all() or not np.all(np.diff(times) > 0):
            raise ValueError("Field plane needs at least two finite increasing times")
        _, ny, nx = data.shape
        if nx <= 2 * cx or ny <= 2 * cy:
            raise ValueError("CPML exclusion leaves no physical aperture")
        dx, dy = float(attrs["dx_m"]), float(attrs["dy_m"])
        x = float(attrs["x_first_m"]) + np.arange(nx) * dx
        y = float(attrs["y_first_m"]) + np.arange(ny) * dy
        physical = (slice(cy, ny - cy), slice(cx, nx - cx))
        centre_x = np.flatnonzero(np.isclose(np.abs(x), np.min(np.abs(x)), rtol=1e-10, atol=dx*1e-10))
        centre_y = np.flatnonzero(np.isclose(np.abs(y), np.min(np.abs(y)), rtol=1e-10, atol=dy*1e-10))
        dt = np.diff(times)
        weights = np.empty(n)
        weights[0], weights[-1] = dt[0] / 2, dt[-1] / 2
        weights[1:-1] = (dt[:-1] + dt[1:]) / 2
        peak_e = np.zeros((ny, nx, 3))
        peak_b = np.zeros_like(peak_e)
        integrated_e2 = np.zeros_like(peak_e)
        integrated_b2 = np.zeros_like(peak_e)
        peak_e_magnitude = np.zeros((ny, nx))
        peak_b_magnitude = np.zeros((ny, nx))
        fluence = np.zeros((ny, nx))
        proxy_fluence = np.zeros((ny, nx))
        columns = {"time_lab_s": times, "time_from_first_fs": (times - times[0]) * 1e15}
        for prefix in ("recorded", "physical"):
            for key in ("net_power_W", "positive_power_W", "negative_power_W", "forward_E_proxy_power_W", "rms_E_V_per_m", "rms_B_T", "peak_E_V_per_m"):
                columns[f"{prefix}_{key}"] = np.empty(n)
        for letter in "xyz":
            columns[f"physical_rms_E{letter}_V_per_m"] = np.empty(n)
        for key in ("peak_E_V_per_m", "peak_B_T", "axis_Ex_V_per_m", "axis_Ey_V_per_m", "axis_Ez_V_per_m", "axis_Bx_T", "axis_By_T", "axis_Bz_T"):
            columns[key] = np.empty(n)
        global_peak = -1.0
        peak_index = 0
        peak_frame = None
        for start in range(0, n, block_frames):
            block = data[start:min(n, start + block_frames)]
            for offset, frame in enumerate(block):
                index = start + offset
                electric = frame["electric_V_per_m"]
                magnetic = frame["magnetic_T"]
                if not np.isfinite(electric).all() or not np.isfinite(magnetic).all():
                    raise ValueError(f"Non-finite field at sample {index}")
                e2, b2 = electric * electric, magnetic * magnetic
                em2, bm2 = e2.sum(axis=2), b2.sum(axis=2)
                np.maximum(peak_e, np.abs(electric), out=peak_e)
                np.maximum(peak_b, np.abs(magnetic), out=peak_b)
                np.maximum(peak_e_magnitude, np.sqrt(em2), out=peak_e_magnitude)
                np.maximum(peak_b_magnitude, np.sqrt(bm2), out=peak_b_magnitude)
                integrated_e2 += weights[index] * e2
                integrated_b2 += weights[index] * b2
                sz = (electric[:, :, 0] * magnetic[:, :, 1] - electric[:, :, 1] * magnetic[:, :, 0]) / MU0
                proxy = EPS0 * C * (e2[:, :, 0] + e2[:, :, 1])
                fluence += weights[index] * sz
                proxy_fluence += weights[index] * proxy
                for prefix, view in (("recorded", (slice(None), slice(None))), ("physical", physical)):
                    s, p = sz[view], proxy[view]
                    columns[f"{prefix}_net_power_W"][index] = s.sum() * dx * dy
                    columns[f"{prefix}_positive_power_W"][index] = np.maximum(s, 0).sum() * dx * dy
                    columns[f"{prefix}_negative_power_W"][index] = np.minimum(s, 0).sum() * dx * dy
                    columns[f"{prefix}_forward_E_proxy_power_W"][index] = p.sum() * dx * dy
                    columns[f"{prefix}_rms_E_V_per_m"][index] = np.sqrt(em2[view].mean())
                    columns[f"{prefix}_rms_B_T"][index] = np.sqrt(bm2[view].mean())
                    columns[f"{prefix}_peak_E_V_per_m"][index] = np.sqrt(em2[view].max())
                for component, letter in enumerate("xyz"):
                    columns[f"physical_rms_E{letter}_V_per_m"][index] = np.sqrt(e2[physical][:, :, component].mean())
                columns["peak_E_V_per_m"][index] = np.sqrt(em2.max())
                columns["peak_B_T"][index] = np.sqrt(bm2.max())
                central_e = electric[np.ix_(centre_y, centre_x)].mean(axis=(0, 1))
                central_b = magnetic[np.ix_(centre_y, centre_x)].mean(axis=(0, 1))
                for component, letter in enumerate("xyz"):
                    columns[f"axis_E{letter}_V_per_m"][index] = central_e[component]
                    columns[f"axis_B{letter}_T"][index] = central_b[component]
                if columns["peak_E_V_per_m"][index] > global_peak:
                    global_peak = columns["peak_E_V_per_m"][index]
                    peak_index, peak_frame = index, frame.copy()
            if start % (block_frames * 128) == 0:
                print(f"field: {min(n, start + block_frames)}/{n} frames", flush=True)
    duration = times[-1] - times[0]
    rms_e, rms_b = np.sqrt(integrated_e2 / duration), np.sqrt(integrated_b2 / duration)
    save_csv(output_dir / "lab_field_timeseries.csv", columns)
    np.savez_compressed(output_dir / "lab_field_timeseries.npz", **columns)
    np.savez_compressed(output_dir / "lab_field_maps.npz", x_m=x, y_m=y,
                        peak_abs_E_components_V_per_m=peak_e, peak_abs_B_components_T=peak_b,
                        peak_E_magnitude_V_per_m=peak_e_magnitude, peak_B_magnitude_T=peak_b_magnitude,
                        time_rms_E_components_V_per_m=rms_e, time_rms_B_components_T=rms_b,
                        native_signed_fluence_J_per_m2=fluence, forward_E_proxy_fluence_J_per_m2=proxy_fluence,
                        peak_frame_E_V_per_m=peak_frame["electric_V_per_m"], peak_frame_B_T=peak_frame["magnetic_T"])
    peak_location = np.unravel_index(np.argmax(np.sum(peak_frame["electric_V_per_m"]**2, axis=2)), (ny, nx))
    summary = {"source_file": str(path.resolve()), "attributes": attrs, "committed_samples": n,
               "frame": "laboratory", "field_is_raw_total": True, "background_subtraction": "none",
               "time_start_s": float(times[0]), "time_end_s": float(times[-1]), "window_duration_fs": float(duration * 1e15),
               "dt_mean_s": float(dt.mean()), "dt_min_s": float(dt.min()), "dt_max_s": float(dt.max()),
               "nyquist_photon_energy_eV": float(4.135667696e-15 / (2 * dt.mean())),
               "nx": nx, "ny": ny, "dx_m": dx, "dy_m": dy,
               "transverse_cpml_cells": [cx, cy], "aperture_note": aperture_note,
               "recorded_aperture_area_m2": dx * dy * nx * ny,
               "physical_aperture_area_m2": dx * dy * (nx - 2*cx) * (ny - 2*cy),
               "peak_E_magnitude_V_per_m": float(global_peak), "peak_B_magnitude_T": float(peak_b_magnitude.max()),
               "peak_abs_E_components_V_per_m": peak_e.max(axis=(0, 1)).tolist(),
               "peak_abs_B_components_T": peak_b.max(axis=(0, 1)).tolist(),
               "peak_E_time_lab_s": float(times[peak_index]), "peak_E_sample_index": peak_index,
               "peak_E_xy_m": [float(x[peak_location[1]]), float(y[peak_location[0]])],
               "central_sample_x_m": x[centre_x].tolist(), "central_sample_y_m": y[centre_y].tolist(),
               "last_frame_peak_E_fraction_of_global_peak": float(columns["peak_E_V_per_m"][-1] / global_peak) if global_peak else 0.0,
               "notes": NOTES}
    for prefix, view in (("recorded", (slice(None), slice(None))), ("physical", physical)):
        sub = {}
        for key in ("net", "positive", "negative"):
            power = columns[f"{prefix}_{key}_power_W"]
            sub[f"{key}_energy_J"] = float(trapezoid(power, times))
        power = columns[f"{prefix}_net_power_W"]
        sub.update(peak_net_power_W=float(power.max()), minimum_net_power_W=float(power.min()),
                   peak_E_magnitude_V_per_m=float(peak_e_magnitude[view].max()),
                   peak_abs_E_components_V_per_m=peak_e[view].max(axis=(0, 1)).tolist(),
                   initial_net_power_W=float(power[0]), final_net_power_W=float(power[-1]),
                   forward_E_proxy_energy_J=float(trapezoid(columns[f"{prefix}_forward_E_proxy_power_W"], times)),
                   space_time_rms_E_V_per_m=float(np.sqrt(np.sum(rms_e[view]**2, axis=2).mean())),
                   space_time_rms_B_T=float(np.sqrt(np.sum(rms_b[view]**2, axis=2).mean())))
        flux_atol = 1e-12 * (abs(sub["positive_energy_J"]) + abs(sub["negative_energy_J"])) + 1e-300
        if not np.isclose(sub["net_energy_J"], fluence[view].sum()*dx*dy, rtol=1e-11, atol=flux_atol):
            raise ValueError("Time/space integration order consistency failed")
        if not np.isclose(sub["net_energy_J"], sub["positive_energy_J"]+sub["negative_energy_J"], rtol=1e-11, atol=flux_atol):
            raise ValueError("Positive/negative flux partition consistency failed")
        summary[prefix + "_aperture"] = sub
    dump_json(output_dir / "lab_field_summary.json", summary)
    if plots:
        plot_fields(output_dir, columns, summary, x, y, peak_e, peak_b,
                    peak_e_magnitude, peak_b_magnitude, rms_e, fluence, proxy_fluence)
    return summary


def plot_fields(out, columns, summary, x, y, peak_e, peak_b, peak_em, peak_bm, rms_e, fluence, proxy):
    plt = plotting()
    from matplotlib.patches import Rectangle
    extent = [(x[0]-summary["dx_m"]/2)*1e6, (x[-1]+summary["dx_m"]/2)*1e6,
              (y[0]-summary["dy_m"]/2)*1e6, (y[-1]+summary["dy_m"]/2)*1e6]
    cx, cy = summary["transverse_cpml_cells"]
    def heat(ax, values, title, unit, signed=False):
        options = {"cmap": "coolwarm" if signed else "magma"}
        if signed:
            limit = max(float(np.max(np.abs(values))), 1e-300)
            options.update(vmin=-limit, vmax=limit)
        img = ax.imshow(values, origin="lower", extent=extent, aspect="equal", **options)
        ax.set(title=title, xlabel="x [um]", ylabel="y [um]")
        ax.grid(False)
        if cx or cy:
            ax.add_patch(Rectangle((extent[0]+cx*summary["dx_m"]*1e6, extent[2]+cy*summary["dy_m"]*1e6),
                         (summary["nx"]-2*cx)*summary["dx_m"]*1e6,
                         (summary["ny"]-2*cy)*summary["dy_m"]*1e6,
                         fill=False, edgecolor="cyan", linewidth=.7))
        plt.colorbar(img, ax=ax, label=unit, shrink=.85)
    for kind, component, magnitude, divisor, unit in (("electric", peak_e, peak_em, 1e6, "MV/m"), ("magnetic", peak_b, peak_bm, 1, "T")):
        fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
        heat(axes[0, 0], magnitude/divisor, "Peak vector magnitude over recorded time", unit)
        for index, ax in enumerate(list(axes.flat)[1:]):
            symbol = "E" if kind == "electric" else "B"
            heat(ax, component[:, :, index]/divisor, "Peak absolute " + symbol + "xyz"[index], unit)
        fig.suptitle("Raw lab " + kind + " field | cyan: physical aperture")
        save_figure(fig, out / ("lab_" + kind + "_peak_maps"))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), constrained_layout=True)
    heat(axes[0], np.sqrt((rms_e*rms_e).sum(axis=2))/1e6, "Time RMS electric magnitude", "MV/m")
    heat(axes[1], fluence, "Native signed Poynting fluence", "J/m2", signed=True)
    heat(axes[2], proxy, "Forward electric-field proxy fluence", "J/m2")
    fig.suptitle("Raw total field; recorded window only; no near-field subtraction")
    save_figure(fig, out / "lab_field_rms_fluence")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    tx = columns["time_from_first_fs"]
    for letter in "xyz":
        axes[0, 0].plot(tx, columns[f"axis_E{letter}_V_per_m"]/1e6, label="E"+letter, lw=.65)
    axes[0, 0].set_ylabel("Central-cell electric field [MV/m]")
    for prefix, label in (("recorded", "Full aperture"), ("physical", "Physical aperture")):
        axes[0, 1].plot(tx, columns[f"{prefix}_rms_E_V_per_m"]/1e6, label=label)
        axes[1, 0].plot(tx, columns[f"{prefix}_net_power_W"]/1e6, label=label+" native net")
    axes[0, 1].set_ylabel("Spatial RMS electric magnitude [MV/m]")
    axes[1, 0].plot(tx, columns["physical_forward_E_proxy_power_W"]/1e6, label="Physical aperture E proxy", ls="--", lw=1)
    axes[1, 0].set_ylabel("Power [MW]")
    axes[1, 1].plot(tx, columns["physical_positive_power_W"]/1e6, label="Physical positive Sz")
    axes[1, 1].plot(tx, columns["physical_negative_power_W"]/1e6, label="Physical negative Sz")
    axes[1, 1].set_ylabel("Pointwise signed power split [MW]")
    for ax in axes.flat:
        ax.set_xlabel("Laboratory time from first sample [fs]")
        ax.legend(fontsize=8)
    fig.suptitle("Fixed lab field plane z = %.6f m | raw fields; native Yee staggering" % summary["attributes"]["plane_z_m"])
    save_figure(fig, out / "lab_field_time_traces")
    central_amplitude2 = sum(columns[f"axis_E{letter}_V_per_m"]**2 for letter in "xyz")
    centre_time = tx[int(np.argmax(central_amplitude2))]
    selected = np.abs(tx - centre_time) <= .5
    if selected.sum() < 3:
        selected = np.ones(tx.size, dtype=bool)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for letter in "xyz":
        axes[0].plot(tx[selected], columns[f"axis_E{letter}_V_per_m"][selected]/1e6, label="E"+letter)
        axes[1].plot(tx[selected], columns[f"axis_B{letter}_T"][selected], label="B"+letter)
    axes[0].set_ylabel("Central-cell electric field [MV/m]")
    axes[1].set_ylabel("Central-cell magnetic field [T]")
    for ax in axes:
        ax.set_xlabel("Laboratory time from first sample [fs]")
        ax.legend()
    fig.suptitle("Central-cell waveform near its peak | native E/B; raw field")
    save_figure(fig, out / "lab_central_waveform_zoom")


def self_test():
    with tempfile.TemporaryDirectory(prefix="aprl-field-analysis-test-") as root:
        run = Path(root) / "run"
        (run / "detectors").mkdir(parents=True)
        out = Path(root) / "out"
        out.mkdir()
        dtype = np.dtype([("electric_V_per_m", "<f8", (3,)), ("magnetic_T", "<f8", (3,))])
        a = np.zeros((5, 4, 6), dtype=dtype)
        amp = np.arange(1., 6.)
        a["electric_V_per_m"][:, :, :, 1] = amp[:, None, None]
        a["magnetic_T"][:, :, :, 0] = -amp[:, None, None]/C
        with h5py.File(run / "detectors/exit-field.h5", "w") as handle:
            g = handle.create_group("field_plane")
            g["fields"], g["time_s"] = a, np.arange(5.)*1e-9
            g["complete"], g["committed_samples"] = 1, 5
            g.attrs.update(dx_m=2e-6, dy_m=3e-6, x_first_m=-5e-6, y_first_m=-4.5e-6, plane_z_m=1.)
        result = analyze_field(run, out, block_frames=2, plots=False)
        exact = EPS0*C*trapezoid(amp**2, np.arange(5.)*1e-9)*24*6e-12
        assert np.isclose(result["recorded_aperture"]["net_energy_J"], exact, rtol=1e-13, atol=0)
        assert np.isclose(result["recorded_aperture"]["forward_E_proxy_energy_J"], exact, rtol=1e-13, atol=0)
        assert result["recorded_aperture"]["negative_energy_J"] == 0
        assert result["peak_E_magnitude_V_per_m"] == 5
        assert np.isclose(result["recorded_aperture"]["space_time_rms_E_V_per_m"], np.sqrt(10.5))
        print("PASS: analytic forward plane wave, signed flux, trapezoid integration, RMS and peaks")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--ledger-only", action="store_true")
    parser.add_argument("--field-only", action="store_true")
    parser.add_argument("--field-file", type=Path, help="Explicit field-plane HDF5 file; defaults to detectors/exit-field.h5")
    parser.add_argument("--no-plots", action="store_true", help="Write numerical summaries for separate plotting")
    parser.add_argument("--block-frames", type=int, default=8)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.run_dir is None or args.output_dir is None:
        parser.error("--run-dir and --output-dir are required")
    if args.block_frames < 1:
        parser.error("--block-frames must be positive")
    if args.field_only and args.ledger_only:
        parser.error("--field-only and --ledger-only are mutually exclusive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {}
    if not args.field_only:
        summary["ledger"] = analyze_ledger(args.run_dir, args.output_dir, plots=not args.no_plots)
    if not args.ledger_only:
        summary["field"] = analyze_field(args.run_dir, args.output_dir, args.block_frames,
                                        plots=not args.no_plots, field_file=args.field_file)
    dump_json(args.output_dir / "field_and_ledger_summary.json", summary)


if __name__ == "__main__":
    main()
