#!/usr/bin/env python3
"""Central high-current power after the short-bunch start-up half.

The window is defined from the known parabolic beam profile, not from observed
radiation amplitude.  This avoids selecting only bright SASE spikes while removing
the low-current field turn-off portion that contaminated a fixed 10%-guard average.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c

import analyze_stable_half as base


CURRENT_FLOORS = (0.6, 0.7, 0.8, 0.9)
MAIN_CURRENT_FLOOR = 0.8
INNER_GUARD_FRACTION = 0.10
RHO = 2.88e-3
PUBLISHED_GAIN_LENGTH_M = 0.67


def need(value, message):
    if not value:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(
        value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()) + "\n")


def stats(envelope, interval, current_floor):
    # For the later half of I(u)=Ipk*(1-u^2), u=0 is the bunch centre and
    # u=1 the tail.  Select I/Ipk >= floor, then remove 10% of that selected
    # interval at both ends so its boundary transitions do not set the mean.
    u_max = math.sqrt(1.0 - current_floor)
    low = INNER_GUARD_FRACTION * u_max
    high = (1.0 - INNER_GUARD_FRACTION) * u_max
    u = ((envelope["time_s"] - interval["t_start_s"]) /
         interval["half_duration_s"])
    use = (u >= low) & (u < high)
    need(np.count_nonzero(use) > 100, "central current window is too short")
    pe = envelope["power_e_W"][use]
    ps = envelope["power_s_W"][use]
    return dict(
        current_floor_fraction=current_floor, u_low=low, u_high=high,
        minimum_actual_current_fraction=float(1.0-high*high),
        maximum_actual_current_fraction=float(1.0-low*low),
        samples=int(np.count_nonzero(use)),
        duration_s=float(np.count_nonzero(use) * envelope["dt_s"]),
        mean_E_proxy_W=float(np.mean(pe)), median_E_proxy_W=float(np.median(pe)),
        p10_E_proxy_W=float(np.quantile(pe, 0.10)),
        p90_E_proxy_W=float(np.quantile(pe, 0.90)), peak_E_proxy_W=float(np.max(pe)),
        mean_signed_flux_W=float(np.mean(ps)),
        signed_over_E_mean=float(np.mean(ps)/np.mean(pe)))


def main(root, output):
    root = Path(root).resolve()
    output = Path(output).resolve()
    need(not output.exists(), f"output already exists: {output}")
    plan = json.loads((root / "plan.json").read_text())
    marker = json.loads((root / "energy-binomial3" / "PIC_COMPLETE.json").read_text())
    need(marker["complete"] is True and marker["run"]["exit_code"] == 0,
         "energy-binomial3 is incomplete")
    output.mkdir(parents=True)

    physics = plan["physics"]
    wavelength = float(physics["wavelength_m"])
    detector_dir = root / "energy-binomial3" / "detectors"
    positions = dict(plan["geometry"]["field_planes"])
    rows = []
    profiles = []
    for module in range(1, 12):
        plane = f"p{module}"
        z_m = float(positions[plane])
        interval = base.particle_interval_and_bunching(
            detector_dir / f"{plane}-particles.h5", module, z_m,
            int(plan["NPART"]), int(plan["active_slices"]), wavelength, 0.0)
        envelope = base.band_envelope(detector_dir / f"{plane}.h5", 6, c/wavelength)
        sensitivity = {f"floor_{int(round(100*floor)):02d}":
                       stats(envelope, interval, floor) for floor in CURRENT_FLOORS}
        main = sensitivity[f"floor_{int(round(100*MAIN_CURRENT_FLOOR)):02d}"]
        rows.append(dict(module=module, z_m=z_m, interval=interval,
                         sensitivity=sensitivity, main=main))
        profiles.append(dict(time_s=envelope["time_s"],
                             power_e_W=envelope["power_e_W"],
                             power_s_W=envelope["power_s_W"]))
        print(json.dumps({"module": module, "z_m": z_m, "main": main}), flush=True)

    peak_current = float(physics["peak_current_A"])
    kinetic_voltage = float(physics["kinetic_energy_MeV"])*1e6
    beam_power = peak_current*kinetic_voltage
    l1d = float(physics["period_m"])/(4*math.pi*math.sqrt(3)*RHO)
    xie_power = (1.6*RHO*(l1d/PUBLISHED_GAIN_LENGTH_M)**2*beam_power)
    rho_power = RHO*beam_power
    mean_power = np.array([row["main"]["mean_E_proxy_W"] for row in rows])
    peak_index = int(np.argmax(mean_power))

    flat = ["module", "z_m", "mean_E_proxy_W", "mean_signed_flux_W",
            "median_E_proxy_W", "p10_E_proxy_W", "p90_E_proxy_W",
            "peak_E_proxy_W", "signed_over_E_mean"]
    with (output / "central_power.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=flat)
        writer.writeheader()
        for row in rows:
            writer.writerow({"module": row["module"], "z_m": row["z_m"],
                             **{key: row["main"][key] for key in flat[2:]}})
    sensitivity_fields = ["module", "z_m", "current_floor_fraction", "u_low", "u_high",
                          "minimum_actual_current_fraction", "maximum_actual_current_fraction",
                          "samples", "duration_s", "mean_E_proxy_W", "mean_signed_flux_W",
                          "median_E_proxy_W", "p10_E_proxy_W", "p90_E_proxy_W",
                          "peak_E_proxy_W", "signed_over_E_mean"]
    with (output / "current_floor_sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=sensitivity_fields)
        writer.writeheader()
        for row in rows:
            for value in row["sensitivity"].values():
                writer.writerow({"module": row["module"], "z_m": row["z_m"],
                                 **{key: value[key] for key in sensitivity_fields[2:]}})

    theory = dict(
        pierce_parameter=RHO, peak_beam_power_W=beam_power,
        rho_peak_beam_power_W=rho_power, gain_length_1d_m=l1d,
        published_ming_xie_gain_length_m=PUBLISHED_GAIN_LENGTH_M,
        xie_3d_saturation_power_W=xie_power,
        xie_equation="1.6*rho*(Lg1D/Lg)^2*Pbeam")
    summary = dict(
        complete=True, selection=("later half after source-ID midpoint; central parabolic-profile "
                                  "region with I/Ipk>=0.8; remove 10% of that selected interval "
                                  "at each end; no radiation-amplitude threshold"),
        main_current_floor_fraction=MAIN_CURRENT_FLOOR,
        inner_guard_fraction=INNER_GUARD_FRACTION,
        rows=rows, theory=theory,
        peak_module=int(rows[peak_index]["module"]),
        peak_z_m=float(rows[peak_index]["z_m"]),
        peak_central_mean_power_W=float(mean_power[peak_index]),
        peak_central_median_power_W=float(rows[peak_index]["main"]["median_E_proxy_W"]),
        maximum_observed_central_peak_power_W=float(max(
            row["main"]["peak_E_proxy_W"] for row in rows)),
        peak_mean_over_xie=float(mean_power[peak_index]/xie_power),
        peak_mean_over_rho_beam=float(mean_power[peak_index]/rho_power),
        limitations=["one physical SASE noise seed", "raw total field; no near-field subtraction",
                     "current floor is scanned from 0.6 to 0.9",
                     "Xie is a steady-state central-slice estimate, not a pulse-integrated prediction"])
    save(output / "RESULT.json", summary)

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "lines.linewidth": 1.8})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), layout="constrained")
    z = np.array([row["z_m"] for row in rows])
    median = np.array([row["main"]["median_E_proxy_W"] for row in rows])
    p10 = np.array([row["main"]["p10_E_proxy_W"] for row in rows])
    p90 = np.array([row["main"]["p90_E_proxy_W"] for row in rows])
    axes[0].semilogy(z, mean_power, "o-", color="#24648c", label="central mean")
    axes[0].semilogy(z, median, "s--", color="#75884b", label="central median")
    axes[0].fill_between(z, p10, p90, color="#cf8a24", alpha=0.18,
                         label="central 10th–90th percentile")
    axes[0].axhline(xie_power, color="#8e4b75", linestyle="-.",
                    label=f"Xie 3D saturation: {xie_power/1e6:.2f} MW")
    axes[0].axvspan(18, 20, color="#777777", alpha=0.10,
                    label="published saturation z")
    axes[0].set_xlabel("undulator-exit z [m]")
    axes[0].set_ylabel("±5% narrow-band power [W]")
    axes[0].set_title("(a) Central high-current region after short-bunch start-up")
    axes[0].grid(True, which="both", alpha=0.22)
    axes[0].legend(fontsize=9)

    row = rows[peak_index]
    profile = profiles[peak_index]
    u = ((profile["time_s"] - row["interval"]["t_start_s"]) /
         row["interval"]["half_duration_s"])
    use_plot = (u >= -0.08) & (u <= 1.03)
    axes[1].plot(u[use_plot], profile["power_e_W"][use_plot]*1e-6,
                 color="#24648c", label="E-proxy envelope")
    axes[1].plot(u[use_plot], profile["power_s_W"][use_plot]*1e-6,
                 color="#75884b", linestyle="--", label="signed-flux envelope")
    axes[1].axvspan(row["main"]["u_low"], row["main"]["u_high"],
                    color="#cf8a24", alpha=0.16, label="retained central window")
    axes[1].set_xlabel("position within later half (0=bunch centre, 1=tail)")
    axes[1].set_ylabel("band-power envelope [MW]")
    axes[1].set_title(f"(b) Module {row['module']} and current-defined window")
    axes[1].grid(alpha=0.22)
    axes[1].legend(fontsize=9)
    fig.savefig(output / "central_power.png", dpi=210, bbox_inches="tight")
    fig.savefig(output / "central_power.pdf", bbox_inches="tight")
    plt.close(fig)

    max_flux_error = max(abs(row["main"]["signed_over_E_mean"]-1) for row in rows)
    save(output / "VALIDATION.json", dict(
        complete=True, field_planes=11, finite=bool(np.isfinite(mean_power).all()),
        main_current_floor_fraction=MAIN_CURRENT_FLOOR,
        inner_guard_fraction=INNER_GUARD_FRACTION,
        signed_flux_agreement_max_relative=float(max_flux_error),
        amplitude_based_window_selection=False))
    (output / "README.md").write_text(
        f"""# Central high-current saturation-power comparison

The main window lies in the later half of the short bunch, requires the known
parabolic current profile to satisfy I/Ipk >= {MAIN_CURRENT_FLOOR}, and removes
{INNER_GUARD_FRACTION:.0%} of that selected interval at each end.  It is defined
without looking at the radiation amplitude.  Current floors 0.6, 0.7, 0.8 and 0.9
are retained for sensitivity.

The maximum central-window mean is {summary['peak_central_mean_power_W']/1e6:.6g} MW
at module {summary['peak_module']} (z={summary['peak_z_m']:.6g} m), compared with
the Ming-Xie 3D central-slice estimate {xie_power/1e6:.6g} MW.  This remains a
single-noise-seed raw-field comparison.
""")


if __name__ == "__main__":
    main(*sys.argv[1:3])
