#!/usr/bin/env python3
"""Stream fixed-laboratory-plane APRL particle diagnostics; never launch APRL.

The input directory must contain run.yaml and detectors/. Production work is
intended for an allocated compute node. --self-test uses only tiny synthetic
files. Each plane is read twice; entry references are indexed by particle ID in
a temporary disk-backed array, not held as complete particle tables in RAM.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

C = 299_792_458.0
QE = 1.602176634e-19
ME = 9.1093837139e-31
H_EV_S = 4.135667696923859e-15
MEC2_EV = ME * C * C / QE
PLANES = (
    ("entry", "entry-particles.h5", "particle_plane"),
    ("1m", "after-1m-particles.h5", "particle_plane"),
    ("2m", "after-2m-particles.h5", "particle_plane"),
    ("3m", "after-3m-particles.h5", "particle_plane"),
    ("exit", "exit-field-ballistic-reference.h5", "ballistic_reference"),
)
FIELDS = ("particle_id", "source_id", "time_s", "position_m", "proper_velocity",
          "charge_C", "mass_kg", "weight", "crossing_direction")


def plain(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    return value


def dump_json(path, value):
    path.write_text(json.dumps(plain(value), ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n")


def write_csv(path, columns):
    names = list(columns)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(zip(*(columns[name] for name in names)))


def checked_group(path, group_name):
    file = h5py.File(path, "r")
    try:
        group = file[group_name]
        ds = group["records"]
        if int(group["complete"][()]) != 1:
            raise ValueError(f"Incomplete particle file: {path}")
        n = int(group["committed_records"][()])
        if n <= 0 or ds.shape != (n,):
            raise ValueError(f"Commit/shape mismatch or empty file: {path}")
        if not set(FIELDS).issubset(ds.dtype.names or ()):
            raise ValueError(f"Missing particle fields: {path}")
        if plain(group.attrs.get("frame")) != "laboratory":
            raise ValueError(f"Expected laboratory frame: {path}")
        if plain(group.attrs.get("proper_velocity_unit")) != "gamma*v/c":
            raise ValueError(f"Unexpected momentum units: {path}")
        return file, group, ds, n
    except Exception:
        file.close()
        raise


def chunks(ds, n, size):
    for start in range(0, n, size):
        yield ds[start:min(start + size, n)]


def accumulate(sums, key, value):
    sums[key] = sums.get(key, np.longdouble(0)) + np.sum(value, dtype=np.longdouble)


def weighted_variance(second, first, weight):
    result = float(second / weight - (first / weight) ** 2)
    tolerance = 2e-12 * max(1.0, abs(float(second / weight)))
    if result < -tolerance:
        raise ValueError(f"Negative weighted variance beyond rounding: {result}")
    return max(result, 0.0)


def scan_plane(path, group_name, reference, entry, n_expected, retarded_origin,
               gamma_anchor, physical_half_width, grid_half_width, chunk_rows):
    file, group, ds, n = checked_group(path, group_name)
    try:
        if n != n_expected:
            raise ValueError(f"Different population count in {path}: {n} vs {n_expected}")
        z = float(group.attrs["plane_z_m"])
        seen = np.zeros(n_expected, dtype=np.bool_)
        sums = {}
        lower = {name: math.inf for name in ("time_s", "gamma", "loss_keV")}
        upper = {name: -math.inf for name in lower}
        raw_weight_min, raw_weight_max = math.inf, -math.inf
        raw_q_min, raw_q_max = math.inf, -math.inf
        mass_per_electron = None
        for batch in chunks(ds, n, chunk_rows):
            ids_u = batch["particle_id"]
            if np.any(ids_u < 1) or np.any(ids_u > n_expected):
                raise ValueError(f"Particle ID outside 1..N in {path}")
            ids = ids_u.astype(np.int64) - 1
            if np.unique(ids).size != ids.size or np.any(seen[ids]):
                raise ValueError(f"Duplicate particle ID in {path}")
            if np.any(batch["crossing_direction"] != 1):
                raise ValueError(f"Upstream or invalid crossing in {path}; selection needs review")
            for name in ("time_s", "position_m", "proper_velocity", "charge_C", "mass_kg", "weight"):
                if not np.all(np.isfinite(batch[name])):
                    raise ValueError(f"Nonfinite {name} in {path}")
            if (np.any(batch["charge_C"] >= 0) or np.any(batch["mass_kg"] <= 0)
                    or np.any(batch["weight"] <= 0)):
                raise ValueError(f"Expected positive-mass, positive-weight electron records: {path}")
            if not np.allclose(batch["position_m"][:, 2], z, rtol=0, atol=1e-12):
                raise ValueError(f"Records do not lie on the configured plane: {path}")
            # charge and mass ALREADY include represented particle count.
            w = -batch["charge_C"] / QE
            this_mass = batch["mass_kg"] / w
            if mass_per_electron is None:
                mass_per_electron = float(this_mass[0])
            if not np.allclose(this_mass, mass_per_electron, rtol=1e-12, atol=0):
                raise ValueError(f"Mixed mass/charge species in {path}")
            if not math.isclose(mass_per_electron, ME, rel_tol=1e-7):
                raise ValueError(f"Unexpected electron rest mass in {path}")
            gamma = np.sqrt(1.0 + np.einsum("ij,ij->i", batch["proper_velocity"], batch["proper_velocity"]))
            if entry:
                reference["gamma"][ids] = gamma
                reference["charge"][ids] = batch["charge_C"]
                reference["mass"][ids] = batch["mass_kg"]
                reference["source_id"][ids] = batch["source_id"]
                loss = np.zeros(gamma.shape)
            else:
                if (not np.array_equal(reference["source_id"][ids], batch["source_id"])
                        or not np.array_equal(reference["charge"][ids], batch["charge_C"])
                        or not np.array_equal(reference["mass"][ids], batch["mass_kg"])):
                    raise ValueError(f"ID-matched source identity or represented weights changed: {path}")
                loss = (reference["gamma"][ids] - gamma) * mass_per_electron * C * C / QE / 1e3
            seen[ids] = True
            tau = (batch["time_s"] - z / C) - retarded_origin
            s_um = C * tau * 1e6  # positive toward later arrivals (tail)
            dg = gamma - gamma_anchor
            x, y = batch["position_m"][:, 0], batch["position_m"][:, 1]
            for key, value in {"w": w, "w2": w*w, "g": w*dg, "g2": w*dg*dg,
                               "s": w*s_um, "s2": w*s_um*s_um, "sg": w*s_um*dg,
                               "x": w*x, "x2": w*x*x, "y": w*y, "y2": w*y*y,
                               "loss": w*loss, "loss2": w*loss*loss,
                               "losing": w*(loss > 0), "gaining": w*(loss < 0),
                               "physical_aperture": w*((np.abs(x) <= physical_half_width[0]) & (np.abs(y) <= physical_half_width[1])),
                               "grid_aperture": w*((np.abs(x) <= grid_half_width[0]) & (np.abs(y) <= grid_half_width[1]))}.items():
                accumulate(sums, key, value)
            for key, value in {"time_s": tau, "gamma": gamma, "loss_keV": loss}.items():
                lower[key] = min(lower[key], float(np.min(value)))
                upper[key] = max(upper[key], float(np.max(value)))
            raw_weight_min = min(raw_weight_min, float(batch["weight"].min()))
            raw_weight_max = max(raw_weight_max, float(batch["weight"].max()))
            raw_q_min = min(raw_q_min, float(batch["charge_C"].min()))
            raw_q_max = max(raw_q_max, float(batch["charge_C"].max()))
        if not np.all(seen):
            raise ValueError(f"Missing particle IDs: {path}")
        w = sums["w"]
        varg = weighted_variance(sums["g2"], sums["g"], w)
        vars_ = weighted_variance(sums["s2"], sums["s"], w)
        cov = float(sums["sg"] / w - sums["s"] * sums["g"] / w**2)
        slope = cov / vars_ if vars_ > 0 else 0.0
        residual_variance = max(0.0, varg - cov*cov/vars_) if vars_ > 0 else varg
        mean_gamma = gamma_anchor + float(sums["g"] / w)
        rest_eV = mass_per_electron * C * C / QE
        result = {
            "source_file": str(path.resolve()), "group": group_name,
            "attrs": plain(dict(group.attrs)), "z_m": z, "records": n,
            "quality": {"complete": True, "commit_equals_shape": True,
                        "unique_ids_1_through_N": True, "all_downstream": True,
                        "finite_all_required_numeric_fields": True,
                        "same_ID_cohort_as_entry": True, "same_source_ids_and_represented_weights_as_entry": True,
                        "raw_relative_weight_min": raw_weight_min, "raw_relative_weight_max": raw_weight_max,
                        "charge_C_min": raw_q_min, "charge_C_max": raw_q_max},
            "represented_electrons": float(w), "absolute_charge_pC": float(w)*QE*1e12,
            "electron_rest_mass_kg": mass_per_electron,
            "effective_particle_count": float(w*w/sums["w2"]),
            "shot_noise_b_rms_reference": math.sqrt(float(sums["w2"]))/float(w),
            "mean_gamma": mean_gamma, "mean_total_energy_MeV": mean_gamma*rest_eV/1e6,
            "mean_kinetic_energy_MeV": (mean_gamma-1)*rest_eV/1e6,
            "sigma_gamma": math.sqrt(varg), "sigma_energy_keV": math.sqrt(varg)*rest_eV/1e3,
            "relative_rms_energy_spread": math.sqrt(varg)/mean_gamma,
            "linear_detrended_sigma_energy_keV": math.sqrt(residual_variance)*rest_eV/1e3,
            "linear_chirp_keV_per_um_tail_positive": slope*rest_eV/1e3,
            "mean_ID_matched_loss_keV": float(sums["loss"] / w),
            "rms_ID_matched_loss_about_mean_keV": math.sqrt(weighted_variance(sums["loss2"], sums["loss"], w)),
            "total_ID_matched_loss_J": float(sums["loss"])*1e3*QE,
            "fraction_losing_energy": float(sums["losing"] / w),
            "fraction_gaining_energy": float(sums["gaining"] / w),
            "sigma_x_um": math.sqrt(weighted_variance(sums["x2"], sums["x"], w))*1e6,
            "sigma_y_um": math.sqrt(weighted_variance(sums["y2"], sums["y"], w))*1e6,
            "physical_xy_aperture_charge_fraction": float(sums["physical_aperture"] / w),
            "full_grid_xy_charge_fraction": float(sums["grid_aperture"] / w),
            "retarded_time_min_s": lower["time_s"], "retarded_time_max_s": upper["time_s"],
            "gamma_min": lower["gamma"], "gamma_max": upper["gamma"],
            "ID_matched_loss_min_keV": lower["loss_keV"], "ID_matched_loss_max_keV": upper["loss_keV"],
        }
        return result
    finally:
        file.close()


def profiles_plane(path, group_name, summary, reference, origin, gamma_anchor,
                   wavelength_entry, wavelength_local, args):
    file, group, ds, n = checked_group(path, group_name)
    try:
        z = summary["z_m"]
        tmin, tmax = summary["retarded_time_min_s"], summary["retarded_time_max_s"]
        slice_dt = args.slice_fs * 1e-15
        left = math.floor(tmin / slice_dt) * slice_dt
        nslice = max(1, math.floor((tmax-left)/slice_dt)+1)
        dt = H_EV_S / (args.spectrum_max_energy_ev * args.spectrum_samples_per_cycle)
        nbin = max(1, math.floor((tmax-tmin)/dt)+1)
        nfft = args.fft_zero_padding * (1 << max(0, (nbin-1).bit_length()))
        if nfft > args.maximum_fft_points:
            raise ValueError(f"FFT needs {nfft} points; adjust --maximum-fft-points explicitly")
        histogram = np.zeros(nbin)
        sw, sg, sg2, sc, sw2 = (np.zeros(nslice) for _ in range(5))
        slice_b_real, slice_b_imag = np.zeros((3,nslice)), np.zeros((3,nslice))
        harmonics = np.arange(1, 4)
        exact_entry = np.zeros(3, np.clongdouble)
        exact_local = np.zeros(3, np.clongdouble)
        elo, ehi = summary["ID_matched_loss_min_keV"], summary["ID_matched_loss_max_keV"]
        if ehi-elo < 1e-10:
            elo, ehi = elo-1e-6, ehi+1e-6
        loss_edges = np.linspace(elo, ehi, 202)
        loss_hist = np.zeros(201)
        phase_t_edges = np.linspace(tmin, tmax if tmax>tmin else tmin+slice_dt, 161)
        glo, ghi = summary["gamma_min"], summary["gamma_max"]
        if ghi-glo < 1e-12:
            glo, ghi = glo-1e-6, ghi+1e-6
        phase_g_edges = np.linspace(glo, ghi, 161)
        phase_hist = np.zeros((160, 160))
        rest_eV = summary["electron_rest_mass_kg"]*C*C/QE
        for batch in chunks(ds, n, args.chunk_rows):
            ids = batch["particle_id"].astype(np.int64)-1
            w = -batch["charge_C"]/QE
            gamma = np.sqrt(1+np.einsum("ij,ij->i", batch["proper_velocity"], batch["proper_velocity"]))
            dg = gamma-gamma_anchor
            tau = (batch["time_s"]-z/C)-origin
            index = np.floor((tau-left)/slice_dt).astype(np.int64)
            sw += np.bincount(index, weights=w, minlength=nslice)
            sg += np.bincount(index, weights=w*dg, minlength=nslice)
            sg2 += np.bincount(index, weights=w*dg*dg, minlength=nslice)
            sc += np.bincount(index, minlength=nslice)
            sw2 += np.bincount(index, weights=w*w, minlength=nslice)
            bins = np.floor((tau-tmin)/dt).astype(np.int64)
            histogram += np.bincount(bins, weights=w, minlength=nbin)
            phase = 2*math.pi*C*(tau-tmin)/wavelength_entry
            local_phase = phase*(wavelength_entry/wavelength_local)
            for h in harmonics:
                real, imag = w*np.cos(h*phase), -w*np.sin(h*phase)
                exact_entry[h-1] += np.sum(real+1j*imag, dtype=np.clongdouble)
                slice_b_real[h-1] += np.bincount(index, weights=real, minlength=nslice)
                slice_b_imag[h-1] += np.bincount(index, weights=imag, minlength=nslice)
                if wavelength_local != wavelength_entry:
                    exact_local[h-1] += np.sum(w*(np.cos(h*local_phase)-1j*np.sin(h*local_phase)), dtype=np.clongdouble)
            loss = (reference["gamma"][ids]-gamma)*rest_eV/1e3
            loss_hist += np.histogram(loss, bins=loss_edges, weights=w)[0]
            phase_hist += np.histogram2d(tau, gamma, bins=(phase_t_edges, phase_g_edges), weights=w)[0]
        if wavelength_local == wavelength_entry:
            exact_local[:] = exact_entry
        occupied = sw>0
        mean = np.full(nslice, np.nan)
        sigma = np.full(nslice, np.nan)
        mean[occupied] = gamma_anchor + sg[occupied]/sw[occupied]
        sigma[occupied] = np.sqrt(np.maximum(0, sg2[occupied]/sw[occupied]-(sg[occupied]/sw[occupied])**2))
        slice_b = np.full((3,nslice),np.nan)
        slice_noise = np.full(nslice,np.nan)
        slice_b[:,occupied] = np.hypot(slice_b_real[:,occupied],slice_b_imag[:,occupied])/sw[occupied]
        slice_noise[occupied] = np.sqrt(sw2[occupied])/sw[occupied]
        energy = np.fft.rfftfreq(nfft, dt)*H_EV_S
        select = energy <= args.spectrum_max_energy_ev
        correction = np.sinc(energy[select]*dt/H_EV_S)
        bf = np.abs(np.fft.rfft(histogram, nfft)[select])/histogram.sum()/correction
        center = left+(np.arange(nslice)+.5)*slice_dt
        summary["bunching"] = {
            "fixed_entry_resonant_wavelength_nm": wavelength_entry*1e9,
            "local_resonant_wavelength_nm": wavelength_local*1e9,
            "exact_b_h1_h2_h3_fixed_entry_resonance": np.abs(exact_entry/summary["represented_electrons"]).astype(float).tolist(),
            "exact_b_h1_h2_h3_local_resonance": np.abs(exact_local/summary["represented_electrons"]).astype(float).tolist(),
            "charge_weighted_slice_rms_b_h1_h2_h3_fixed_entry_resonance": np.sqrt(np.sum(slice_b[:,occupied]**2*sw[occupied],axis=1)/sw.sum()).tolist(),
            "uniform_full_slice_window_b_h1_h2_h3": np.abs(np.sinc(C*slice_dt*harmonics/wavelength_entry)).tolist(),
            "histogram_bin_width_s": dt, "histogram_bins": nbin, "fft_points": nfft,
            "nyquist_energy_eV": H_EV_S/(2*dt), "display_energy_range_eV": [0, args.spectrum_max_energy_ev],
            "fft_energy_spacing_eV": H_EV_S/(nfft*dt),
            "finite_window_resolution_eV": H_EV_S/max(tmax-tmin, dt),
            "window": "rectangular complete observed bunch; no taper and no mean subtraction",
            "bin_response_correction": "divide magnitude by sinc(f * histogram_dt); approximate sub-bin deconvolution",
            "zero_padding_note": "interpolates Fourier samples; does not improve finite-bunch spectral resolution",
            "time_coordinate": "arrival time at fixed z; s=c*(t-z/c-retarded_origin), positive toward tail",
        }
        summary["slicing"] = {"width_fs": args.slice_fs, "bins": nslice,
                              "occupied_bins": int(np.count_nonzero(occupied)),
                              "slice_particle_count_min": int(sc[occupied].min()),
                              "definition": "Eulerian arrival-time slices; membership can change between planes"}
        return {
            "slice_time_fs": center*1e15, "slice_s_um": center*C*1e6,
            "slice_records": sc.astype(np.int64), "slice_electrons": sw,
            "slice_current_A": sw*QE/slice_dt,
            "slice_mean_energy_MeV": mean*rest_eV/1e6,
            "slice_sigma_energy_keV": sigma*rest_eV/1e3,
            "slice_b_h1": slice_b[0], "slice_b_h2": slice_b[1], "slice_b_h3": slice_b[2],
            "slice_shot_noise_b_rms_reference": slice_noise,
            "spectrum_energy_eV": energy[select], "spectrum_b": bf,
            "loss_edges_keV": loss_edges, "loss_weighted_electrons": loss_hist,
            "phase_time_edges_fs": phase_t_edges*1e15,
            "phase_energy_edges_MeV": phase_g_edges*rest_eV/1e6,
            "phase_weighted_electrons": phase_hist,
        }
    finally:
        file.close()


def plots(output, summaries, profiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .23,
                         "savefig.dpi": 180, "pdf.fonttype": 42})
    labels = [f"z = {s['z_m']:g} m" for s in summaries]
    colors = plt.cm.viridis(np.linspace(.05,.9,len(labels)))
    z = [s["z_m"] for s in summaries]
    def save(fig, name):
        fig.savefig(output / f"{name}.png", bbox_inches="tight")
        fig.savefig(output / f"{name}.pdf", bbox_inches="tight")
        plt.close(fig)
    fig, axes = plt.subplots(2,2,figsize=(11,7.2),constrained_layout=True)
    axes[0,0].plot(z,[s["mean_ID_matched_loss_keV"] for s in summaries],"o-")
    axes[0,0].axhline(0,color="k",lw=.7)
    axes[0,0].set_ylabel("Mean energy loss from entry [keV/electron]")
    axes[0,1].plot(z,[s["sigma_energy_keV"] for s in summaries],"o-",label="Projected RMS")
    axes[0,1].plot(z,[s["linear_detrended_sigma_energy_keV"] for s in summaries],"s--",label="Linear trend removed")
    axes[0,1].set_ylabel("Energy spread [keV]"); axes[0,1].legend()
    for h in range(3):
        axes[1,0].semilogy(z,[s["bunching"]["exact_b_h1_h2_h3_fixed_entry_resonance"][h] for s in summaries],"o-",label=f"h = {h+1}")
    axes[1,0].axhline(summaries[0]["shot_noise_b_rms_reference"],color="gray",ls=":",label="Independent-electron RMS reference")
    axes[1,0].set_ylabel("Exact bunching factor at entry resonance"); axes[1,0].legend(fontsize=8)
    for key,label in (("sigma_x_um",r"$\sigma_x$"),("sigma_y_um",r"$\sigma_y$")):
        axes[1,1].plot(z,[s[key] for s in summaries],"o-",label=label)
    axes[1,1].set_ylabel("Transverse RMS [µm]"); axes[1,1].legend()
    for ax in axes.flat: ax.set_xlabel("Fixed laboratory plane z [m]")
    fig.suptitle("Fixed-plane particle evolution — identical particle IDs at every plane")
    save(fig,"particle-evolution")
    fig,axes = plt.subplots(4,1,figsize=(11,10.5),sharex=True,constrained_layout=True)
    for p,label,color in zip(profiles,labels,colors):
        for ax,key in zip(axes,("slice_current_A","slice_mean_energy_MeV","slice_sigma_energy_keV","slice_b_h1")):
            ax.plot(p["slice_s_um"],p[key],color=color,label=label)
    for ax,title in zip(axes,("Current [A]","Slice mean total energy [MeV]","Slice RMS energy spread [keV]","Slice bunching at fixed entry λ")):
        ax.set_ylabel(title)
    axes[0].legend(ncol=3); axes[-1].set_xlabel("Retarded arrival coordinate c(t - z/c - t₀) [µm]; tail positive")
    axes[-1].axhline(summaries[0]["bunching"]["uniform_full_slice_window_b_h1_h2_h3"][0],color="gray",ls=":",label="Uniform full-window form factor")
    axes[-1].set_ylim(bottom=0); axes[-1].legend(fontsize=8)
    fig.suptitle(f"Arrival-time slices ({summaries[0]['slicing']['width_fs']:g} fs); slices are not fixed particle cohorts")
    save(fig,"particle-slices")
    fig,axes=plt.subplots(2,1,figsize=(11,7),constrained_layout=True)
    for s,p,label,color in zip(summaries,profiles,labels,colors):
        axes[0].semilogy(p["spectrum_energy_eV"][1:],np.maximum(p["spectrum_b"][1:],1e-16),color=color,lw=.8,label=label)
        if s["z_m"]!=summaries[0]["z_m"]:
            edges=p["loss_edges_keV"]
            density=p["loss_weighted_electrons"]/s["represented_electrons"]/np.diff(edges)
            axes[1].plot((edges[1:]+edges[:-1])/2,density,color=color,label=label)
    axes[0].axhline(summaries[0]["shot_noise_b_rms_reference"],ls=":",color="gray")
    axes[0].set_xlabel("Photon-equivalent frequency hf [eV]"); axes[0].set_ylabel("Bunching factor |b(f)|")
    axes[0].legend(ncol=3); axes[0].set_ylim(bottom=1e-7,top=1.1)
    axes[1].set_xlabel("ID-matched energy loss from entry [keV]; negative = gain")
    axes[1].set_ylabel("Charge-weighted probability density [1/keV]"); axes[1].legend()
    fig.suptitle("Full-bunch arrival-time form factor and particle energy-loss distribution")
    save(fig,"particle-bunching-and-loss")
    fig,axes=plt.subplots(1,len(profiles),figsize=(16,3.8),constrained_layout=True)
    positives=[p["phase_weighted_electrons"][p["phase_weighted_electrons"]>0] for p in profiles]
    norm=LogNorm(vmin=min(float(p.min()) for p in positives),vmax=max(2.,max(float(p.max()) for p in positives)))
    shared_energy_limits=(min(p["phase_energy_edges_MeV"][0] for p in profiles),max(p["phase_energy_edges_MeV"][-1] for p in profiles))
    for ax,p,label in zip(np.atleast_1d(axes),profiles,labels):
        h=p["phase_weighted_electrons"].T
        im=ax.pcolormesh(p["phase_time_edges_fs"],p["phase_energy_edges_MeV"],np.ma.masked_equal(h,0),norm=norm,cmap="magma",rasterized=True)
        ax.set_title(label); ax.set_xlabel("Retarded arrival time [fs]"); ax.grid(False)
        ax.set_ylim(*shared_energy_limits)
    np.atleast_1d(axes)[0].set_ylabel("Total energy [MeV]")
    fig.colorbar(im,ax=np.atleast_1d(axes).tolist(),label="Represented electrons / bin (shared scale)",shrink=.8)
    save(fig,"particle-longitudinal-phase-space")


def analyze(args):
    run=Path(args.run_dir).resolve(); output=Path(args.output_dir).resolve()
    output.mkdir(parents=True,exist_ok=True)
    config=yaml.safe_load((run/"run.yaml").read_text())
    magnetic=config["sources"]["magnetic_elements"]
    undulator=next(e for e in magnetic if e["type"]=="planar-undulator")
    K=float(undulator["strength_parameter"]); period=float(undulator["period"])
    cells=np.asarray(config["mesh"]["cells"][:2],dtype=float)
    spacing=np.asarray(config["mesh"]["cell_size"][:2],dtype=float)
    pml=np.asarray(config["boundary"]["cells"][:2],dtype=float)
    center=np.asarray(config["mesh"].get("center",[0,0,0])[:2],dtype=float)
    if np.any(center!=0):
        raise ValueError("This aperture diagnostic requires mesh transverse center=[0,0]")
    grid_half=cells*spacing/2; physical_half=grid_half-pml*spacing
    if np.any(physical_half<=0):
        raise ValueError("Invalid physical aperture")
    paths=[(tag,run/"detectors"/name,group) for tag,name,group in PLANES]
    metadata=[]
    for tag,path,group_name in paths:
        file,group,ds,n=checked_group(path,group_name)
        metadata.append((n,plain(group.attrs.get("run_id")),plain(group.attrs.get("configuration_digest_fnv1a64")),plain(group.attrs.get("source_revision"))))
        file.close()
    if any(item!=metadata[0] for item in metadata):
        raise ValueError("Plane files differ in count, run identity, or configuration digest")
    n=metadata[0][0]
    if (run/"input-plan.json").exists():
        plan=json.loads((run/"input-plan.json").read_text())
        if int(plan["records"])!=n:
            raise ValueError("Expected input particle population differs from detector count")
    file,group,ds,_=checked_group(paths[0][1],paths[0][2])
    first=ds[:min(args.chunk_rows,n)]
    origin=float(np.min(first["time_s"])-float(group.attrs["plane_z_m"])/C)
    gamma_anchor=float(np.mean(np.sqrt(1+np.einsum("ij,ij->i",first["proper_velocity"],first["proper_velocity"]))))
    file.close(); del first
    summaries=[]; all_profiles=[]
    # No input is modified. Temporary per-ID references are removed on exit.
    with tempfile.TemporaryDirectory(prefix="particle-reference-",dir=output) as temporary:
        dtype=np.dtype([("gamma","f8"),("charge","f8"),("mass","f8"),("source_id","u8")])
        reference=np.memmap(Path(temporary)/"entry-by-id.bin",mode="w+",dtype=dtype,shape=(n,))
        entry_lambda=None
        for i,(tag,path,group_name) in enumerate(paths):
            print(f"Scanning {tag}: {n:,} records",flush=True)
            result=scan_plane(path,group_name,reference,i==0,n,origin,gamma_anchor,physical_half,grid_half,args.chunk_rows)
            result["label"]=tag
            local_lambda=period*(1+K*K/2)/(2*result["mean_gamma"]**2)
            if entry_lambda is None:
                entry_lambda=local_lambda
            print(f"Profiles and BF: {tag}",flush=True)
            profile=profiles_plane(path,group_name,result,reference,origin,gamma_anchor,entry_lambda,local_lambda,args)
            summaries.append(result); all_profiles.append(profile)
            np.savez_compressed(output/f"particles-{tag}.npz",**profile)
            write_csv(output/f"slices-{tag}.csv",{k:v for k,v in profile.items() if k.startswith("slice_")})
            write_csv(output/f"bunching-{tag}.csv",{k:v for k,v in profile.items() if k.startswith("spectrum_")})
            write_csv(output/f"loss-{tag}.csv",{"bin_left_keV":profile["loss_edges_keV"][:-1],"bin_right_keV":profile["loss_edges_keV"][1:],"represented_electrons":profile["loss_weighted_electrons"]})
        del reference
    result={
        "format":"aprl-fixed-lab-particle-analysis-v1", "run_dir":str(run),
        "parameters":vars(args), "retarded_time_origin_s":origin,
        "undulator_K":K,"undulator_period_m":period,
        "physical_xy_half_width_m":physical_half.tolist(),"grid_xy_half_width_m":grid_half.tolist(),
        "planes":summaries,
        "definitions":{
            "gamma":"sqrt(1+ux^2+uy^2+uz^2), proper velocity u=gamma*v/c",
            "weights":"abs(charge_C)/e; charge and mass already include represented electrons; do not multiply by relative weight again",
            "energy":"per-electron total gamma*m_e*c^2; kinetic=(gamma-1)*m_e*c^2; energy loss is entry minus plane for identical particle_id",
            "projected_energy_spread":"charge-weighted population RMS about weighted mean; relative spread sigma_gamma/mean_gamma",
            "linear_detrending":"weighted least-squares gamma against tail-positive s=c*(t-z/c-t0); residual is not an irreversible-heating or true slice-spread estimator",
            "bunching":"|sum(w*exp(-i*2*pi*c*(t-tref)/lambda))/sum(w)|; all z coordinates on a plane are equal",
            "resonance":"planar on-axis estimate lambda_u*(1+K^2/2)/(2*mean_gamma^2); fixed entry resonance compares like frequencies; local resonance is also tabulated",
            "aperture":"charge fraction whose lab x,y at the crossing lie inside the simulation transverse inner-CPML/full-grid rectangle; not a trajectory-integrated participation fraction",
        },
        "limitations":[
            "Fixed-z crossings occur at different laboratory times; these are transported beam diagnostics, not simultaneous laboratory or boost snapshots.",
            "Energy loss is descriptive. No control run is used to separate radiation, space-charge, boundary, or other contributions.",
            "Global full-bunch BF may cancel local microbunching phases; this script also saves slices but does not equate low global BF with no local modulation.",
            "The finite-window, histogram BF spectrum is approximate; harmonics are summed exactly from all unbinned arrival times.",
            "ID-matched losses use direct differences of double-precision gamma, then extended-precision reduction; changes close to machine precision times gamma are not resolved.",
            "Slice BF uses the same entry resonance but finite 1 fs windows by default; compare its local noise and finite-window leakage, not only full-bunch noise.",
            "A uniform complete rectangular slice has nonzero BF=abs(sinc(c*slice_dt*h/lambda)); its h1-h3 values are recorded. Partially filled edge slices and current-envelope variation can increase this. Raw slice BF is not all microbunching; amplitudes must not be subtracted as a correction.",
            "Independent-electron 1/sqrt(N_eff) is a reference, not an uncertainty estimate or proof of independent phases in these data.",
            "High-frequency particle form factors are not proof of resolved electromagnetic radiation or high-harmonic convergence.",
            "Read run.yaml, input metadata, and the field/ledger analysis for crop, aperture-injection, initial-field, and closure limitations.",
        ],
    }
    dump_json(output/"particle-summary.json",result)
    numeric=[k for k,v in summaries[0].items() if isinstance(v,(int,float))]
    write_csv(output/"particle-summary.csv",{"plane":[s["label"] for s in summaries],**{k:[s[k] for s in summaries] for k in numeric}})
    if not args.no_plots:
        plots(output,summaries,all_profiles)
    print(f"Particle analysis complete: {output}",flush=True)
    return result


def self_test(keep_output=None):
    """Check cohort matching, moments, time-phase BF and rejection of duplicates."""
    with tempfile.TemporaryDirectory(prefix="aprl-particle-self-test-") as tmp:
        root=Path(tmp); run=root/"run"; (run/"detectors").mkdir(parents=True)
        cfg={"mesh":{"cells":[160,160,256],"cell_size":[5e-6,5e-6,1e-6]},"boundary":{"cells":[10,10,20]},"sources":{"magnetic_elements":[{"type":"planar-undulator","period":.0234,"strength_parameter":1.092714657005}]}}
        (run/"run.yaml").write_text(yaml.safe_dump(cfg))
        dtype=np.dtype([("particle_id","u8"),("source_id","u8"),("time_s","f8"),("position_m","f8",3),("proper_velocity","f8",3),("charge_C","f8"),("mass_kg","f8"),("weight","f8"),("crossing_direction","i4")])
        n=96; ids=np.arange(1,n+1); g0=1174+np.tile([-1.,0.,1.],n//3)
        wavelength=.0234*(1+1.092714657005**2/2)/(2*1174**2)
        arrival=np.arange(n)*wavelength/C
        for i,((tag,name,group_name),z) in enumerate(zip(PLANES,[-.03,1,2,3,4.03])):
            a=np.zeros(n,dtype=dtype); a["particle_id"]=ids; a["source_id"]=ids+1000
            a["time_s"]=arrival+z/C; a["position_m"][:,2]=z
            a["proper_velocity"][:,2]=np.sqrt((g0-i*.1)**2-1)
            a["charge_C"]=-QE; a["mass_kg"]=ME; a["weight"]=1; a["crossing_direction"]=1
            a=a[np.random.default_rng(i).permutation(n)]
            with h5py.File(run/"detectors"/name,"w") as f:
                g=f.create_group(group_name); g.attrs.update({"frame":"laboratory","proper_velocity_unit":"gamma*v/c","plane_z_m":z,"run_id":"synthetic","configuration_digest_fnv1a64":"synthetic"})
                g.create_dataset("records",data=a); g["committed_records"]=n; g["complete"]=1
        args=parser().parse_args(["--run-dir",str(run),"--output-dir",str(keep_output or root/"out"),"--chunk-rows","31"] + ([] if keep_output else ["--no-plots"]))
        result=analyze(args)
        for i,s in enumerate(result["planes"]):
            assert abs(s["mean_gamma"]-(1174-i*.1))<1e-10
            assert abs(s["sigma_gamma"]-math.sqrt(2/3))<1e-10
            assert abs(s["mean_ID_matched_loss_keV"]-i*.1*MEC2_EV/1e3)<1e-7
            assert abs(s["bunching"]["exact_b_h1_h2_h3_fixed_entry_resonance"][0]-1)<1e-8
        # A duplicate must fail rather than silently alter the population.
        target=run/"detectors"/PLANES[1][1]
        with h5py.File(target,"r+") as f:
            d=f["particle_plane/records"]; a=d[:]; a["particle_id"][1]=a["particle_id"][0]; d[:]=a
        try:
            args.output_dir=str(root/"bad")
            args.no_plots=True
            analyze(args)
        except ValueError as exc:
            assert "Duplicate particle ID" in str(exc)
        else:
            raise AssertionError("Duplicate ID was not rejected")
    print("SELF-TEST PASS: moments, matched loss, exact time-phase BF, duplicate rejection")


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir"); p.add_argument("--output-dir")
    p.add_argument("--chunk-rows",type=int,default=262144)
    p.add_argument("--slice-fs",type=float,default=1.)
    p.add_argument("--spectrum-max-energy-ev",type=float,default=500.)
    p.add_argument("--spectrum-samples-per-cycle",type=int,default=16)
    p.add_argument("--fft-zero-padding",type=int,default=4)
    p.add_argument("--maximum-fft-points",type=int,default=16777216)
    p.add_argument("--no-plots",action="store_true")
    p.add_argument("--self-test",action="store_true")
    return p


def main():
    p=parser(); args=p.parse_args()
    if args.self_test:
        self_test(args.output_dir); return
    if not args.run_dir or not args.output_dir:
        p.error("--run-dir and --output-dir are required")
    if min(args.chunk_rows,args.slice_fs,args.spectrum_max_energy_ev,args.spectrum_samples_per_cycle,args.fft_zero_padding,args.maximum_fft_points)<=0:
        p.error("All numeric analysis options must be positive")
    analyze(args)


if __name__=="__main__":
    main()
