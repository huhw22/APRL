#!/usr/bin/env python3
"""Time-resolved narrow-band power in the settled half of the short SPARC bunch.

The original growth analysis Fourier transformed the complete detector record.  That
is useful for pulse energy, but it is not the right estimator for the long-bunch
central saturation power: this short bunch is only twice the total eleven-module
slippage length.  Here source IDs recover the longitudinal slice number, the later
half of the bunch is selected, and configurable guards remove the field turn-on and
turn-off portions before power statistics are formed.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy import fft
from scipy.constants import c, epsilon_0, mu_0


RHO = 2.88e-3
PUBLISHED_DURATION_S = 12.67e-12
PUBLISHED_FINAL_ENERGY_J = (96e-6, 111e-6)
GUARD_FRACTIONS = (0.0, 0.05, 0.10, 0.15, 0.20)
MAIN_GUARD = 0.10


def need(value, message):
    if not value:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(
        value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()) + "\n")


def lagrange(nodes, target):
    nodes = np.asarray(nodes, dtype=float)
    weights = np.ones(len(nodes))
    for j in range(len(nodes)):
        for k in range(len(nodes)):
            if j != k:
                weights[j] *= (target - nodes[k]) / (nodes[j] - nodes[k])
    return weights


def particle_interval_and_bunching(path, module, z_m, n_part, active_slices,
                                   wavelength_m, guard_fraction):
    with h5py.File(path, "r") as handle:
        group = handle["particle_plane"]
        need(int(group["complete"][()]) == 1, "unclosed particle plane")
        records = group["records"][:]
    source_id = records["source_id"].astype(np.int64)
    slice_id = (source_id - 1) // n_part
    need(np.all((slice_id >= 0) & (slice_id < active_slices)), "bad source slice")
    time_s = records["time_s"]
    weight = records["weight"]
    half = active_slices // 2
    need(active_slices % 2 == 0, "analysis expects an even active-slice count")
    first_last = half - 1
    second_first = half
    first_times = time_s[slice_id == first_last]
    second_times = time_s[slice_id == second_first]
    need(len(first_times) and len(second_times), "missing central source slices")
    # The midpoint between the two adjacent slice-centre crossing times is the
    # physical support boundary.  The nominal duration is preferable to extrema,
    # which are broadened by energy spread and late-stage FEL dynamics.
    t_start = 0.5 * (float(np.median(first_times)) + float(np.median(second_times)))
    second_slices = active_slices - half
    half_duration = second_slices * wavelength_m / c
    t_end = t_start + half_duration
    core_start = t_start + guard_fraction * half_duration
    core_end = t_end - guard_fraction * half_duration

    guard_slices = int(math.ceil(guard_fraction * second_slices))
    selected = ((slice_id >= half + guard_slices) &
                (slice_id < active_slices - guard_slices))
    need(np.count_nonzero(selected) > 0, "empty guarded particle selection")
    t = time_s[selected]
    w = weight[selected]
    phase_cycles = c * (t - z_m / c) / wavelength_m
    phasor = np.exp(2j * np.pi * phase_cycles)
    total = np.sum(w, dtype=np.longdouble)
    global_b1 = abs(np.sum(w * phasor, dtype=np.complex128) / total)
    bins = np.floor(phase_cycles - np.min(phase_cycles)).astype(np.int64)
    charge = np.bincount(bins, weights=w)
    real = np.bincount(bins, weights=w * phasor.real)
    imag = np.bincount(bins, weights=w * phasor.imag)
    use = charge > 0
    local_squared = (real[use] ** 2 + imag[use] ** 2) / charge[use] ** 2
    local_rms = math.sqrt(float(np.sum(charge[use] * local_squared) /
                                np.sum(charge[use])))
    return dict(
        module=module, t_start_s=t_start, t_end_s=t_end,
        core_start_s=core_start, core_end_s=core_end,
        half_duration_s=half_duration, half_first_slice=half,
        half_last_slice=active_slices - 1, guard_slices=guard_slices,
        core_first_slice=half + guard_slices,
        core_last_slice=active_slices - guard_slices - 1,
        selected_particle_count=int(np.count_nonzero(selected)),
        selected_weight=float(total), global_b1=float(global_b1),
        local_rms_b1=float(local_rms),
        particle_time_min_s=float(np.min(t)), particle_time_max_s=float(np.max(t)))


def detector_cube(path, cpml):
    with h5py.File(path, "r") as handle:
        group = handle["field_plane"]
        need(int(group["complete"][()]) == 1, "unclosed field plane")
        count = int(group["committed_samples"][()])
        times = group["time_s"][:count]
        dt0 = float(np.median(np.diff(times)))
        ticks = np.rint((times - times[0]) / dt0).astype(np.int64)
        tick_ld = ticks.astype(np.longdouble)
        time_ld = times.astype(np.longdouble)
        tick_center = tick_ld - np.mean(tick_ld)
        time_center = time_ld - np.mean(time_ld)
        dt = float(np.sum(tick_center * time_center) /
                   np.sum(tick_center * tick_center))
        offset = float(np.mean(time_ld) - np.longdouble(dt) * np.mean(tick_ld))
        residual_ticks = np.max(abs((times - (offset + dt * ticks)) / dt))
        need(np.all(np.diff(ticks) > 0) and residual_ticks < 1e-5,
             "bad detector time grid")
        n = int(ticks[-1] + 1)
        missing = np.setdiff1d(np.arange(n), ticks)
        need(len(missing) <= 2, "too many detector gaps")
        data = group["fields"]
        ny, nx = data.shape[1:3]
        py, px = ny - 2 * cpml, nx - 2 * cpml
        need(py > 0 and px > 0, "empty physical aperture")
        area = float(group.attrs["dx_m"]) * float(group.attrs["dy_m"])
        cube = np.empty((n, py * px, 6), dtype=np.float64)
        for start in range(0, count, 32):
            stop = min(count, start + 32)
            block = data[start:stop]
            electric = block["electric_V_per_m"][:, cpml:ny-cpml, cpml:nx-cpml, :]
            magnetic = block["magnetic_T"][:, cpml:ny-cpml, cpml:nx-cpml, :]
            need(np.isfinite(electric).all() and np.isfinite(magnetic).all(),
                 "nonfinite field")
            cube[ticks[start:stop], :, :3] = electric.reshape(stop-start, -1, 3)
            cube[ticks[start:stop], :, 3:] = magnetic.reshape(stop-start, -1, 3)
    for gap in missing:
        nodes = ticks[np.argsort(abs(ticks-gap), kind="stable")[:8]]
        cube[gap] = np.tensordot(lagrange(nodes, gap), cube[nodes], axes=(0, 0))
    uniform_time = offset + dt * np.arange(n)
    return uniform_time, dt, area, cube, count, len(missing)


def band_envelope(path, cpml, photon_frequency_hz):
    time_s, dt, area, cube, native_samples, missing_samples = detector_cube(path, cpml)
    n, pixels, _ = cube.shape
    frequency = fft.rfftfreq(n, dt)
    band = abs(frequency / photon_frequency_hz - 1.0) <= 0.05
    band_indices = np.flatnonzero(band)
    need(len(band_indices) >= 2, "target band is unresolved")
    power_e = np.zeros(n, dtype=np.float64)
    power_s = np.zeros(n, dtype=np.float64)
    # Build the analytic signal from positive-frequency target-band bins.  Its
    # squared modulus gives the optical-cycle-averaged narrow-band envelope.
    for start in range(0, pixels, 64):
        stop = min(pixels, start + 64)
        values = np.array(cube[:, start:stop, :][:, :, [0, 1, 3, 4]].transpose(1, 2, 0),
                          order="C", copy=True)
        transform = fft.rfft(values, axis=-1, workers=8, overwrite_x=True)
        analytic_spectrum = np.zeros((stop-start, 4, n), dtype=np.complex128)
        analytic_spectrum[:, :, band_indices] = 2.0 * transform[:, :, band]
        envelope = fft.ifft(analytic_spectrum, axis=-1, workers=8,
                            overwrite_x=True)
        power_e += epsilon_0 * c * area * 0.5 * np.sum(
            abs(envelope[:, 0, :]) ** 2 + abs(envelope[:, 1, :]) ** 2, axis=0)
        power_s += area / (2.0 * mu_0) * np.sum(
            (envelope[:, 0, :] * envelope[:, 3, :].conj() -
             envelope[:, 1, :] * envelope[:, 2, :].conj()).real, axis=0)
        del values, transform, analytic_spectrum, envelope
    del cube
    need(np.isfinite(power_e).all() and np.isfinite(power_s).all(),
         "nonfinite band envelope")
    return dict(time_s=time_s, dt_s=dt, power_e_W=power_e, power_s_W=power_s,
                native_samples=native_samples, samples=n,
                missing_samples=missing_samples,
                band_bins=int(np.count_nonzero(band)))


def power_statistics(envelope, interval):
    time_s = envelope["time_s"]
    e_power = envelope["power_e_W"]
    s_power = envelope["power_s_W"]
    result = {}
    for guard in GUARD_FRACTIONS:
        duration = interval["half_duration_s"]
        start = interval["t_start_s"] + guard * duration
        end = interval["t_end_s"] - guard * duration
        use = (time_s >= start) & (time_s < end)
        need(np.count_nonzero(use) > 100, "guarded field interval is too short")
        e = e_power[use]
        s = s_power[use]
        result[f"guard_{int(round(100*guard)):02d}"] = dict(
            guard_fraction=guard, samples=int(np.count_nonzero(use)),
            start_s=start, end_s=end, duration_s=float(np.count_nonzero(use) * envelope["dt_s"]),
            mean_E_proxy_W=float(np.mean(e)), median_E_proxy_W=float(np.median(e)),
            p90_E_proxy_W=float(np.quantile(e, 0.90)), peak_E_proxy_W=float(np.max(e)),
            mean_signed_flux_W=float(np.mean(s)), median_signed_flux_W=float(np.median(s)),
            signed_over_E_mean=float(np.mean(s) / np.mean(e)))
    return result


def ledger_analysis(case_dir, output):
    path = case_dir / "ledger" / "energy-ledger.h5"
    with h5py.File(path, "r") as handle:
        group = handle["energy_ledger"]
        committed = int(group["committed_records"][()])
        need(int(group["complete"][()]) == 1, "unclosed energy ledger")
        rows = group["records"][:committed]
        attrs = {key: (value.decode() if isinstance(value, bytes) else value.item()
                       if isinstance(value, np.generic) else value)
                 for key, value in group.attrs.items()}
    a = lambda key: np.asarray(rows[key], dtype=np.longdouble)
    delta_u = a("field_energy_interior_J") - a("field_energy_interior_J")[0]
    source_terms = (a("prescribed_boris_work_J") + a("boundary_source_work_J") +
                    a("moving_remap_field_change_J") + a("inlet_overwrite_work_J") +
                    a("frozen_constraint_work_J") + a("carrier_constraint_work_J"))
    reconstructed = (a("delta_particle_kinetic_box_J") + delta_u +
                     a("outward_field_energy_J") - source_terms)
    residual = a("balance_residual_boris_J")
    reconstruction_error = np.max(abs(reconstructed - residual))
    need(reconstruction_error < np.longdouble("1e-18"), "ledger equation mismatch")

    columns = {
        "step": rows["step"], "time_box_s": rows["time_box_s"],
        "delta_particle_kinetic_box_J": a("delta_particle_kinetic_box_J"),
        "delta_field_energy_J": delta_u,
        "outward_field_energy_J": a("outward_field_energy_J"),
        "negative_source_terms_J": -source_terms,
        "balance_residual_boris_J": residual,
        "field_exchange_residual_J": a("field_exchange_residual_J"),
        "particle_work_balance_residual_J": a("particle_work_balance_residual_J"),
        "relative_balance_exchange": a("relative_balance_exchange"),
        "relative_balance_initial": a("relative_balance_initial"),
        "delta_particle_kinetic_lab_J": a("delta_particle_kinetic_lab_J")}
    with (output / "energy_ledger.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for index in range(len(rows)):
            writer.writerow([columns[key][index] for key in columns])

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "lines.linewidth": 1.7})
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 9.0), sharex=True,
                             layout="constrained")
    step = np.asarray(rows["step"])
    axes[0].plot(step, np.asarray(a("delta_particle_kinetic_box_J")*1e6, float),
                 label=r"$\Delta K'_{particles}$", color="#24648c")
    axes[0].plot(step, np.asarray(delta_u*1e6, float),
                 label=r"$\Delta U'_{field}$", color="#cf8a24")
    axes[0].plot(step, np.asarray(a("outward_field_energy_J")*1e6, float),
                 label="outward field", color="#75884b")
    axes[0].plot(step, np.asarray(-source_terms*1e6, float),
                 label="negative source/remap terms", color="#8e4b75")
    axes[0].set_ylabel("signed ledger term [µJ]")
    axes[0].set_title("(a) Boost-frame moving-window energy terms")
    axes[0].legend(ncol=2)

    axes[1].plot(step, np.asarray(residual*1e9, float), label="total Boris ledger residual",
                 color="#8e4b75")
    axes[1].plot(step, np.asarray(a("field_exchange_residual_J")*1e9, float),
                 label="field-exchange residual", color="#cf8a24", linestyle="--")
    axes[1].plot(step, np.asarray(a("particle_work_balance_residual_J")*1e9, float),
                 label="particle-work residual", color="#24648c", linestyle=":")
    axes[1].set_ylabel("residual [nJ]")
    axes[1].set_title("(b) Accounting residuals (legacy field ledger)")
    axes[1].legend(ncol=2)

    axes[2].plot(step, np.asarray(100*a("relative_balance_exchange"), float),
                 label="residual / exchange scale", color="#8e4b75")
    axes[2].plot(step, np.asarray(100*a("relative_balance_initial"), float),
                 label="residual / initial kinetic", color="#24648c", linestyle="--")
    axes[2].axhline(0, color="#777777", linewidth=0.8)
    axes[2].set_ylabel("relative residual [%]")
    axes[2].set_xlabel("Maxwell step")
    axes[2].set_title("(c) Relative residual; closure_valid remains false")
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.22)
    fig.savefig(output / "energy_ledger.png", dpi=210, bbox_inches="tight")
    fig.savefig(output / "energy_ledger.pdf", bbox_inches="tight")
    plt.close(fig)

    result = dict(
        complete=True, records=committed, first_step=int(rows["step"][0]),
        final_step=int(rows["step"][-1]), format_version=int(attrs["format_version"]),
        frame=attrs["frame"], discrete_operator_audit=False,
        closure_valid_all=bool(np.all(rows["closure_valid"] != 0)),
        final_residual_J=float(residual[-1]),
        maximum_absolute_residual_J=float(np.max(abs(residual))),
        final_relative_exchange_percent=float(100*a("relative_balance_exchange")[-1]),
        maximum_absolute_relative_exchange_percent=float(
            100*np.max(abs(a("relative_balance_exchange")))),
        final_particle_work_residual_J=float(a("particle_work_balance_residual_J")[-1]),
        final_field_exchange_residual_J=float(a("field_exchange_residual_J")[-1]),
        final_lab_particle_kinetic_change_J=float(a("delta_particle_kinetic_lab_J")[-1]),
        equation_reconstruction_max_error_J=float(reconstruction_error),
        interpretation=("Complete format-v5 boosted-frame moving-window diagnostic, but "
                        "discrete_operator_audit was disabled and closure_valid is false; "
                        "this is not a strict discrete or closed laboratory-frame ledger."))
    save(output / "ENERGY_LEDGER_RESULT.json", result)
    return result


def render_power(output, rows, profiles, theory):
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "lines.linewidth": 1.8})
    z = np.array([row["z_m"] for row in rows])
    mean = np.array([row["main"]["mean_E_proxy_W"] for row in rows])
    p90 = np.array([row["main"]["p90_E_proxy_W"] for row in rows])
    signed = np.array([row["main"]["mean_signed_flux_W"] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), layout="constrained")
    axes[0].semilogy(z, mean, "o-", color="#24648c", label="core mean E proxy")
    axes[0].semilogy(z, signed, "s--", color="#75884b", label="core mean signed flux")
    axes[0].semilogy(z, p90, ":", color="#cf8a24", label="core 90th percentile")
    axes[0].axhspan(theory["published_mean_power_W"][0],
                    theory["published_mean_power_W"][1], color="#8e4b75", alpha=0.14,
                    label="published 28 m pulse-mean scale")
    axes[0].axhline(theory["rho_peak_beam_power_W"], color="#333333", linestyle="-.",
                    label=r"$\rho P_{beam,peak}$")
    axes[0].set_xlabel("undulator-exit z [m]")
    axes[0].set_ylabel("±5% narrow-band power [W]")
    axes[0].set_title("(a) Later-half core after 10% guards at both ends")
    axes[0].grid(True, which="both", alpha=0.22)
    axes[0].legend(fontsize=9)

    peak_index = int(np.argmax(mean))
    profile = profiles[peak_index]
    interval = rows[peak_index]["interval"]
    coordinate = ((profile["time_s"] - interval["t_start_s"]) /
                  interval["half_duration_s"])
    within = (coordinate >= -0.08) & (coordinate <= 1.08)
    axes[1].plot(coordinate[within], profile["power_e_W"][within] * 1e-6,
                 color="#24648c", label="E-proxy envelope")
    axes[1].plot(coordinate[within], profile["power_s_W"][within] * 1e-6,
                 color="#75884b", linestyle="--", label="signed-flux envelope")
    axes[1].axvspan(0, MAIN_GUARD, color="#999999", alpha=0.18)
    axes[1].axvspan(1-MAIN_GUARD, 1, color="#999999", alpha=0.18,
                    label="discarded edge transients")
    axes[1].axvline(MAIN_GUARD, color="#555555", linewidth=0.8)
    axes[1].axvline(1-MAIN_GUARD, color="#555555", linewidth=0.8)
    axes[1].set_xlabel("position within later half (0=center boundary, 1=bunch end)")
    axes[1].set_ylabel("band-power envelope [MW]")
    axes[1].set_title(f"(b) Module {rows[peak_index]['module']} envelope and retained core")
    axes[1].grid(alpha=0.22)
    axes[1].legend(fontsize=9)
    fig.savefig(output / "stable_half_power.png", dpi=210, bbox_inches="tight")
    fig.savefig(output / "stable_half_power.pdf", bbox_inches="tight")
    plt.close(fig)


def main(root, output):
    root = Path(root).resolve()
    output = Path(output).resolve()
    need(not output.exists(), f"output already exists: {output}")
    marker = json.loads((root / "energy-binomial3" / "PIC_COMPLETE.json").read_text())
    need(marker["complete"] is True and marker["run"]["exit_code"] == 0,
         "energy-binomial3 is not complete")
    plan = json.loads((root / "plan.json").read_text())
    run_card = yaml.safe_load((root / "energy-binomial3" / "run.yaml").read_text())
    output.mkdir(parents=True)

    physics = plan["physics"]
    wavelength = float(physics["wavelength_m"])
    f0 = c / wavelength
    n_part = int(plan["NPART"])
    active_slices = int(plan["active_slices"])
    detector_dir = root / "energy-binomial3" / "detectors"
    rows = []
    profiles = []
    plane_positions = dict(plan["geometry"]["field_planes"])
    for module in range(1, 12):
        plane = f"p{module}"
        z_m = float(plane_positions[plane])
        interval = particle_interval_and_bunching(
            detector_dir / f"{plane}-particles.h5", module, z_m, n_part,
            active_slices, wavelength, MAIN_GUARD)
        envelope = band_envelope(detector_dir / f"{plane}.h5", 6, f0)
        stats = power_statistics(envelope, interval)
        main_stats = stats[f"guard_{int(round(100*MAIN_GUARD)):02d}"]
        row = dict(module=module, z_m=z_m, interval=interval,
                   detector=dict(samples=envelope["samples"],
                                 native_samples=envelope["native_samples"],
                                 missing_samples=envelope["missing_samples"],
                                 dt_s=envelope["dt_s"], band_bins=envelope["band_bins"]),
                   guard_statistics=stats, main=main_stats)
        rows.append(row)
        profiles.append(dict(time_s=envelope["time_s"],
                             power_e_W=envelope["power_e_W"],
                             power_s_W=envelope["power_s_W"]))
        print(json.dumps({"module": module, "z_m": z_m, "main": main_stats,
                          "local_rms_b1": interval["local_rms_b1"]}), flush=True)

    kinetic_voltage = float(physics["kinetic_energy_MeV"]) * 1e6
    peak_beam_power = float(physics["peak_current_A"]) * kinetic_voltage
    rho_power = RHO * peak_beam_power
    published_power = [value / PUBLISHED_DURATION_S
                       for value in PUBLISHED_FINAL_ENERGY_J]
    theory = dict(
        pierce_parameter=RHO, peak_current_A=float(physics["peak_current_A"]),
        kinetic_voltage_V=kinetic_voltage, peak_beam_power_W=peak_beam_power,
        rho_peak_beam_power_W=rho_power,
        published_duration_s=PUBLISHED_DURATION_S,
        published_final_energy_J=list(PUBLISHED_FINAL_ENERGY_J),
        published_mean_power_W=published_power,
        comparison_note=("Published 28 m values are whole-pulse energies from a longer "
                         "450 pC bunch; their energy/duration is only a pulse-mean scale, "
                         "not an exact time-local estimator."))
    mean_power = np.array([row["main"]["mean_E_proxy_W"] for row in rows])
    peak_index = int(np.argmax(mean_power))
    peak_row = rows[peak_index]
    summary = dict(
        complete=True, selection=("source-ID slices 847..1693 define the later bunch half; "
                                  "the primary statistic discards 10% of that interval at "
                                  "both ends before averaging the time-resolved ±5% envelope"),
        main_guard_fraction=MAIN_GUARD, rows=rows, theory=theory,
        peak_module=int(peak_row["module"]), peak_z_m=float(peak_row["z_m"]),
        peak_core_mean_power_W=float(peak_row["main"]["mean_E_proxy_W"]),
        peak_core_mean_signed_flux_W=float(peak_row["main"]["mean_signed_flux_W"]),
        peak_power_over_rho_peak_beam=float(
            peak_row["main"]["mean_E_proxy_W"] / rho_power),
        peak_power_over_published_mean_range=[
            float(peak_row["main"]["mean_E_proxy_W"] / max(published_power)),
            float(peak_row["main"]["mean_E_proxy_W"] / min(published_power))],
        limitations=["one physical SASE noise seed",
                     "raw total field; no near-field subtraction",
                     "published whole-pulse energy and current time-local power are not identical observables",
                     "10% edge guard is reported together with 0/5/15/20% sensitivity"])

    flat_fields = ["module", "z_m", "local_rms_b1", "global_b1",
                   "mean_E_proxy_W", "mean_signed_flux_W", "median_E_proxy_W",
                   "p90_E_proxy_W", "peak_E_proxy_W", "signed_over_E_mean"]
    with (output / "stable_half_growth.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=flat_fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(module=row["module"], z_m=row["z_m"],
                local_rms_b1=row["interval"]["local_rms_b1"],
                global_b1=row["interval"]["global_b1"], **{
                    key: row["main"][key] for key in flat_fields[4:]}))
    with (output / "guard_sensitivity.csv").open("w", newline="") as stream:
        fields = ["module", "z_m", "guard_fraction", "samples", "duration_s",
                  "mean_E_proxy_W", "mean_signed_flux_W", "median_E_proxy_W",
                  "p90_E_proxy_W", "peak_E_proxy_W", "signed_over_E_mean"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for stats in row["guard_statistics"].values():
                writer.writerow({"module": row["module"], "z_m": row["z_m"],
                                 **{key: stats[key] for key in fields[2:]}})
    save(output / "RESULT.json", summary)
    render_power(output, rows, profiles, theory)
    ledger = ledger_analysis(root / "energy-binomial3", output)

    readme = f"""# SPARC 11-segment settled-half power and energy ledger

The primary radiation statistic uses source-ID slices 847..1693 (the later half of
the 1694-slice bunch) and removes 10% of that half at each end.  The field is first
restricted to ±5% around the nominal resonance, then converted to a time-resolved
analytic-signal power envelope; no full-record energy divided by bunch duration is
used.  Guard sensitivities of 0, 5, 10, 15 and 20% are retained in
`guard_sensitivity.csv`.

The peak retained-core mean is {summary['peak_core_mean_power_W']/1e6:.6g} MW at
module {summary['peak_module']} (z={summary['peak_z_m']:.6g} m), or
{summary['peak_power_over_rho_peak_beam']:.6g} times rho times peak beam power.
The comparison with published 28 m whole-pulse energy is only a scale comparison,
because the observables and bunch durations differ.

`energy_ledger.png` is a separate boosted-frame moving-window diagnostic.  Its HDF5
is complete, but this run had `discrete_operator_audit: false` and every stored
`closure_valid` flag is false.  The final residual is
{ledger['final_residual_J']:.9g} J ({ledger['final_relative_exchange_percent']:.6g}%
of the stored exchange scale).  It must not be presented as a strict discrete
closure proof or as a closed laboratory-frame energy balance.
"""
    (output / "README.md").write_text(readme)
    save(output / "VALIDATION.json", dict(
        complete=True, case_complete=True, field_planes=11,
        main_guard_fraction=MAIN_GUARD,
        finite_main_rows=bool(np.isfinite(mean_power).all()),
        signed_flux_agreement_max_relative=float(max(
            abs(row["main"]["signed_over_E_mean"] - 1.0) for row in rows)),
        energy_ledger_equation_reconstruction_max_error_J=
            ledger["equation_reconstruction_max_error_J"],
        discrete_operator_audit=bool(run_card["energy_ledger"]["discrete_operator_audit"]),
        figures=["stable_half_power.png", "energy_ledger.png"]))


if __name__ == "__main__":
    main(*sys.argv[1:3])
