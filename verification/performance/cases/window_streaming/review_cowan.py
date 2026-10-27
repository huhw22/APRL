#!/usr/bin/env python3
"""Independent same-release-state comparison of Cowan and controlled held.

The frozen PLAN defines the research screens. Held is a comparison, not truth;
ordinary field energy is a diagnostic, not a Cowan exact discrete invariant.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from review_startup import (C, particle_moments, read_csv, relative,
                            sampled_gaussian_force)

QE = 1.602176634e-19
ME = 9.1093837139e-31
REST_EV = np.longdouble(ME)*C*C/QE
AXES = ("x", "y", "z")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def velocities(p):
    return np.array([p["u"+a] for a in AXES], dtype=np.longdouble).T


def delta_electron_energy(p, reference):
    u, v = velocities(p), velocities(reference)
    gamma = np.sqrt(1+np.sum(u*u, axis=1))
    gamma0 = np.sqrt(1+np.sum(v*v, axis=1))
    return REST_EV*np.sum((u-v)*(u+v), axis=1)/(gamma+gamma0)


def weighted_stats(value, particles):
    weight = np.abs(particles["q_C"].astype(np.longdouble))/QE
    mean = np.average(value, weights=weight)
    return dict(mean=float(mean), sigma=float(np.sqrt(np.average((value-mean)**2, weights=weight))),
                rms=float(np.sqrt(np.average(value*value, weights=weight))),
                minimum=float(np.min(value)), maximum=float(np.max(value)))


def energy_stats(p, initial):
    u = velocities(p)
    u2 = np.sum(u*u, axis=1)
    kinetic = REST_EV*u2/(np.sqrt(1+u2)+1)
    return dict(kinetic_eV_per_physical_electron=weighted_stats(kinetic, p),
                delta_kinetic_eV_per_physical_electron=weighted_stats(delta_electron_energy(p, initial), p))


def orbit_arrays(history):
    sigma = np.array([np.sqrt(history["var_"+a+"_m2"]/history["var_"+a+"_m2"][0])-1 for a in AXES]).T
    cov = np.array([history["cov_"+a+"_u"+a+"_m"] for a in AXES]).T
    return sigma, cov


def per_axis_difference(a, b):
    amplitude = np.max(np.abs(b), axis=0)
    difference = np.max(np.abs(a-b), axis=0)
    return dict(maximum_absolute_by_axis=difference.tolist(),
                reference_amplitude_by_axis=amplitude.tolist(),
                relative_by_axis=(difference/np.maximum(amplitude, 1e-300)).tolist())


def analyze(root, require_complete=False):
    plan = json.loads((root/"PLAN.json").read_text())
    limits = plan["functional_limits"]
    cases, raw, functional, incomplete = {}, {}, [], []
    def check(case, name, passed):
        functional.append(dict(case=case, name=name, passed=bool(passed)))

    for folder in sorted((root/"r1").iterdir()):
        if not folder.is_dir() or folder.name == "frozen": continue
        needed = [folder/f for f in ("metadata.json", "history.csv", "initial-particles.csv", "release-particles.csv", "final-particles.csv", "EXIT_CODE")]
        if not all(p.is_file() for p in needed):
            incomplete.append(folder.name)
            continue
        name = folder.name
        meta = json.loads((folder/"metadata.json").read_text())
        h = read_csv(folder/"history.csv")
        initial, release, final = [read_csv(folder/f"{tag}-particles.csv") for tag in ("initial", "release", "final")]
        check(name, "runner_exit_zero", (folder/"EXIT_CODE").read_text().strip() == "0")
        check(name, "physical_core_binding", "[binding] result=PASS failures=0" in (folder/"binding.log").read_text())
        check(name, "finite_complete_history", len(h) == meta["steps"]+1 and h["step"][-1] == meta["steps"] and all(np.isfinite(h[n]).all() for n in h.dtype.names))
        check(name, "release_state_bitwise_unchanged", (folder/"initial-particles.csv").read_bytes() == (folder/"release-particles.csv").read_bytes())
        check(name, "particle_identity_count_charge_mass", len(initial) == meta["particles"] and len(final) == len(initial) and all(np.array_equal(initial[k], final[k]) for k in ("id", "q_C", "m_kg")))
        pm = particle_moments(final, initial)
        check(name, "independent_final_delta_K", relative(pm["delta_K_J"], h["delta_K_J"][-1]) < 1e-10)
        check(name, "independent_final_variances", relative(pm["var_x"], np.array([h["var_"+a+"_m2"][-1] for a in AXES])) < 1e-10)
        check(name, "independent_final_covariance", relative(pm["cov_xu"], np.array([h["cov_"+a+"_u"+a+"_m"][-1] for a in AXES])) < 1e-10)
        check(name, "independent_analytic_Gaussian_moment", relative(sampled_gaussian_force(initial, meta), meta["analytic_sampled_force_covariance_m_per_s"]) < 1e-12)
        check(name, "initial_Gauss", h["Gauss_relative_max"][0] <= limits["initial_gauss_max_relative"])
        check(name, "continuity", np.max(h["continuity_relative_rho_over_dt"]) <= limits["continuity_max_relative"])
        check(name, "release_energy_starts_at_zero", h["delta_K_J"][0] == 0 and h["JdotE_cumulative_J"][0] == 0)
        check(name, "fixed_requested_physics", meta["seed"] == plan["physics"]["seed"] and meta["particles"] == plan["physics"]["macroparticles"] and meta["charge_C"] == plan["physics"]["charge_C"] and meta["rest_sigma_m"] == plan["physics"]["rest_sigma_m"] and abs(h["time_s"][-1]/plan["physics"]["release_duration_s"]-1) < 1e-14)
        held = None
        if meta["held_prehistory_steps"]:
            hh = read_csv(folder/"held-history.csv")
            start = read_csv(folder/"held-start-particles.csv")
            hsteps = meta["held_prehistory_steps"]
            check(name, "held_complete_finite_history", len(hh) == hsteps+1 and hh["time_s"][-1] == 0 and all(np.isfinite(hh[k]).all() for k in hh.dtype.names))
            check(name, "held_duration", abs(hh["time_s"][0]+meta["held_prehistory_duration_s"]) <= 1e-28)
            check(name, "held_current_work_bookkeeping", relative(hh["JdotE_cumulative_J"][-1], meta["held_grid_work_J"]) < 1e-12)
            check(name, "held_continuity", max(np.max(hh["continuity_relative_rho_over_dt"]), meta["held_max_continuity_relative"]) <= limits["continuity_max_relative"])
            u = velocities(initial)
            velocity = np.asarray(C*u/np.sqrt(1+np.sum(u*u, axis=1))[:, None], float)
            predicted = np.array([initial[a+"_m"] for a in AXES]).T - meta["held_prehistory_duration_s"]*velocity
            observed = np.array([start[a+"_m"] for a in AXES]).T
            check(name, "held_start_is_prescribed_backtrace", np.max(np.abs(predicted-observed)) <= 64*np.finfo(float).eps*meta["rest_sigma_m"] and all(np.array_equal(start["u"+a], initial["u"+a]) for a in AXES))
            held = dict(duration_s=meta["held_prehistory_duration_s"], steps=hsteps,
                        grid_work_J=meta["held_grid_work_J"], field_energy_change_native_J=float(hh["UE_J"][-1]+hh["UB_native_minus_half_J"][-1]-hh["UE_J"][0]-hh["UB_native_minus_half_J"][0]),
                        max_gauss_relative=float(np.max(hh["Gauss_relative_max"])),
                        max_continuity_relative=float(np.max(hh["continuity_relative_rho_over_dt"])))
        if meta["static_source_only"]:
            check(name, "static_sources_do_no_work", np.array_equal(initial, final) and np.all(h["delta_K_J"] == 0) and np.all(h["JdotE_cumulative_J"] == 0))
        field = h["UE_J"]+h["UB_centered_J"]
        du = field-field[0]
        if meta["mode"] == "cowan-fieldonly":
            check(name, "Cowan_static_field_energy_drift", np.max(np.abs(du))/field[0] <= limits["static_field_energy_drift_over_initial"])
        previous_reference = None
        if meta["mode"] == "qs":
            tag = "0" if meta["beta"] == 0 else "023"
            old = root.parent/"startup_self_field_20261027"/"r1"/f"qs-beta{tag}-{meta['grid']}"
            check(name, "old_qs_history_bitwise_regression", old.is_dir() and (old/"history.csv").read_bytes() == (folder/"history.csv").read_bytes())
            check(name, "old_qs_particles_bitwise_regression", old.is_dir() and (old/"final-particles.csv").read_bytes() == (folder/"final-particles.csv").read_bytes())
            previous_reference = dict(directory=str(old), history_sha256=digest(old/"history.csv") if old.is_dir() else None,
                                      frozen_source=str(root.parent/"startup_self_field_20261027"/"r1"/"frozen"/"startup_self_field_probe.cpp"))
        sigma, cov = orbit_arrays(h)
        cases[name] = dict(meta=meta, directory=str(folder), previous_qs_regression_reference=previous_reference,
            initial_particle_sha256=digest(folder/"initial-particles.csv"),
            final_delta_K_J=float(h["delta_K_J"][-1]), minimum_delta_K_J=float(np.min(h["delta_K_J"])),
            maximum_delta_K_J=float(np.max(h["delta_K_J"])),
            final_sigma_growth_fraction=sigma[-1].tolist(), final_covariance_xu_m=cov[-1].tolist(),
            energy_statistics=energy_stats(final, initial),
            initial_field_energy_J=float(field[0]), final_field_energy_change_J=float(du[-1]),
            max_field_energy_change_over_initial=float(np.max(np.abs(du))/field[0]),
            final_grid_work_J=float(h["JdotE_cumulative_J"][-1]),
            final_particle_minus_grid_work_J=float(h["delta_K_J"][-1]-h["JdotE_cumulative_J"][-1]),
            final_field_plus_grid_work_J=float(du[-1]+h["JdotE_cumulative_J"][-1]),
            final_particle_plus_field_change_J=float(h["delta_K_J"][-1]+du[-1]),
            max_particle_plus_field_change_over_Uinf=float(np.max(np.abs(h["delta_K_J"]+du))/meta["gaussian_rest_self_energy_J"]),
            initial_UE_J=float(h["UE_J"][0]), final_UE_J=float(h["UE_J"][-1]),
            initial_UB_centered_J=float(h["UB_centered_J"][0]), final_UB_centered_J=float(h["UB_centered_J"][-1]),
            initial_UB_native_J=float(h["UB_native_plus_half_J"][0]), final_UB_native_J=float(h["UB_native_plus_half_J"][-1]),
            max_gauss_relative=float(np.max(h["Gauss_relative_max"])),
            max_continuity_relative=float(np.max(h["continuity_relative_rho_over_dt"])),
            magnetic_divergence={k:v for k,v in meta.items() if "divB" in k},
            preparation_timing_s=dict(initialization=meta["initialization_wall_s"], held_prehistory=meta["held_prehistory_wall_s"], total=meta["initialization_wall_s"]+meta["held_prehistory_wall_s"]),
            held_history=held)
        raw[(meta["beta"],meta["grid"],meta["mode"])] = dict(name=name, h=h, particles=final, sigma=sigma, cov=cov)

    for beta in plan["physics"]["beta"]:
        matching = [c for c in cases.values() if c["meta"]["beta"] == beta]
        if matching: check(str(beta), "all_cases_have_same_release_state", len({c["initial_particle_sha256"] for c in matching}) == 1)
    comparisons = {}
    for beta in plan["physics"]["beta"]:
        for grid in plan["grids"]:
            group = {mode:raw.get((beta,grid,mode)) for mode in ("qs", "cowan", "held50", "held100")}
            if not all(group.values()): continue
            ref, short = group["held100"], group["held50"]
            ref_k = ref["h"]["delta_K_J"]
            held_envelope = float(np.max(np.abs(short["h"]["delta_K_J"]-ref_k)))
            amplitude = float(np.max(np.abs(ref_k)))
            energy_limit = held_envelope+.05*amplitude
            out = dict(held_reference_name=ref["name"], held_history_envelope_J=held_envelope,
                       held100_exchange_amplitude_J=amplitude, provisional_energy_difference_limit_J=energy_limit,
                       held_history_envelope_over_exchange=held_envelope/max(amplitude,1e-300),
                       provisional_orbit_relative_limit=plan["comparison_policy"]["provisional_orbit_screen_relative"], pairs={})
            for mode in ("cowan", "qs", "held50"):
                item = group[mode]
                same_t = np.array_equal(item["h"]["time_s"], ref["h"]["time_s"])
                check(item["name"], "same_grid_comparison_time_samples", same_t)
                if not same_t: continue
                diff = item["h"]["delta_K_J"]-ref_k
                sigma_diff = per_axis_difference(item["sigma"], ref["sigma"])
                cov_diff = per_axis_difference(item["cov"], ref["cov"])
                orbit_error = max(sigma_diff["relative_by_axis"]+cov_diff["relative_by_axis"])
                electron_diff = delta_electron_energy(item["particles"], ref["particles"])
                energy_difference = float(np.max(np.abs(diff)))
                out["pairs"][mode+"_minus_held100"] = dict(
                    maximum_absolute_delta_K_difference_J=energy_difference,
                    final_delta_K_difference_J=float(diff[-1]),
                    maximum_difference_over_Uinf=energy_difference/cases[item["name"]]["meta"]["gaussian_rest_self_energy_J"],
                    max_additional_loss_or_reduced_gain_J=float(max(0,-np.min(diff))),
                    max_additional_actual_negative_excursion_J=float(max(0,np.max(np.maximum(-item["h"]["delta_K_J"],0)-np.maximum(-ref_k,0)))),
                    sigma_growth_difference=sigma_diff, covariance_xu_difference=cov_diff,
                    final_particle_energy_difference_eV_per_physical_electron=weighted_stats(electron_diff,item["particles"]),
                    final_energy_spread_difference_eV=cases[item["name"]]["energy_statistics"]["kinetic_eV_per_physical_electron"]["sigma"]-cases[ref["name"]]["energy_statistics"]["kinetic_eV_per_physical_electron"]["sigma"],
                    initial_field_energy_difference_J=cases[item["name"]]["initial_field_energy_J"]-cases[ref["name"]]["initial_field_energy_J"],
                    final_field_energy_change_difference_J=cases[item["name"]]["final_field_energy_change_J"]-cases[ref["name"]]["final_field_energy_change_J"],
                    final_grid_work_difference_J=cases[item["name"]]["final_grid_work_J"]-cases[ref["name"]]["final_grid_work_J"],
                    energy_screen=energy_difference <= energy_limit,
                    orbit_screen=orbit_error <= plan["comparison_policy"]["provisional_orbit_screen_relative"])
            candidate = out["pairs"].get("cowan_minus_held100")
            out["candidate_research_screen_pass"] = bool(candidate and candidate["energy_screen"] and candidate["orbit_screen"])
            comparisons[f"beta{beta}-{grid}"] = out

    changed_grid = {}
    for beta in plan["physics"]["beta"]:
        for mode in ("qs", "cowan", "held50", "held100"):
            base = raw.get((beta,"base",mode))
            if base is None: continue
            for grid in ("wide", "fine"):
                item = raw.get((beta,grid,mode))
                if item is None: continue
                stride, remainder = divmod(len(item["h"])-1,len(base["h"])-1)
                if remainder or not stride:
                    check(item["name"], "grid_comparison_integer_time_stride", False)
                    continue
                same_time = np.allclose(item["h"]["time_s"][::stride],base["h"]["time_s"],rtol=1e-14,atol=0)
                check(item["name"], "grid_comparison_same_physical_times", same_time)
                if not same_time: continue
                diff = item["h"]["delta_K_J"][::stride]-base["h"]["delta_K_J"]
                energy_amplitude = float(np.max(np.abs(base["h"]["delta_K_J"])))
                changed_grid[f"{item['name']}_minus_{base['name']}"] = dict(
                    matching_time_samples=len(base["h"]), finer_time_stride=stride,
                    maximum_absolute_delta_K_difference_J=float(np.max(np.abs(diff))),
                    maximum_difference_over_base_exchange=float(np.max(np.abs(diff)))/max(energy_amplitude,1e-300),
                    final_delta_K_difference_J=float(diff[-1]),
                    sigma_growth_difference=per_axis_difference(item["sigma"][::stride],base["sigma"]),
                    covariance_xu_difference=per_axis_difference(item["cov"][::stride],base["cov"]),
                    final_particle_energy_difference_eV_per_physical_electron=weighted_stats(delta_electron_energy(item["particles"],base["particles"]),item["particles"]),
                    interpretation="Diagnostic only; finite boundary or mesh/time and fixed-N shape resolution changed. This does not replace the previous convergence gate.")

    expected = {(b,g,m) for b in plan["physics"]["beta"] for g in plan["grids"] for m in plan["modes"] if "fieldonly" not in m or b == 0}
    missing = sorted(expected-set(raw))
    complete = not missing and not incomplete
    functional_pass = all(c["passed"] for c in functional)
    result = dict(status=("COMPLETE" if complete else "PARTIAL")+("_FUNCTIONAL_OK" if functional_pass else "_FUNCTIONAL_FAILURE"),
                  complete=complete, completed_cases=len(cases), expected_cases=len(expected),
                  missing_cases=[dict(beta=b,grid=g,mode=m) for b,g,m in missing], incomplete_directories=incomplete,
                  all_available_functional_checks_pass=functional_pass, functional_checks=functional,
                  cases=cases, comparisons=comparisons, changed_grid_diagnostics=changed_grid, production_ready=False,
                  all_comparative_research_screens_pass=bool(complete and comparisons and all(c["candidate_research_screen_pass"] for c in comparisons.values())),
                  interpretation="Same-release-state comparisons to controlled held histories, not physical truth. Passing a wide held-history envelope is not high-precision energy acceptance and does not overturn the previous mesh/ordinary-energy gate. Initial field energy belongs to initial state; held source work is separately recorded. Ordinary Cowan field energy is not an exact discrete invariant. Single-run initialization costs are not production speedups; no moving-window injection has been accepted.")
    (root/"REVIEW.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(status=result["status"], completed_cases=len(cases), expected_cases=len(expected),
        functional_checks=len(functional), failures=[c for c in functional if not c["passed"]],
        candidate_screens={name:dict(passed=c["candidate_research_screen_pass"],
            energy_difference_pJ=c["pairs"]["cowan_minus_held100"]["maximum_absolute_delta_K_difference_J"]*1e12,
            energy_limit_pJ=c["provisional_energy_difference_limit_J"]*1e12) for name,c in comparisons.items()}),indent=2))
    if not functional_pass or (require_complete and not complete): raise SystemExit(1)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("root",type=Path)
    parser.add_argument("--require-complete",action="store_true")
    args=parser.parse_args()
    analyze(args.root.resolve(),args.require_complete)
