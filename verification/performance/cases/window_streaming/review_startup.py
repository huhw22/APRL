#!/usr/bin/env python3
"""Independent release-state, particle-moment and energy diagnostics.

Research thresholds come from the pre-run PLAN, never from the observed values.
Zero/scaled fields are deliberately invalid initial constraints, not references.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import erf

C = 299792458.0
EPS0 = 8.8541878128e-12


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def particle_moments(p, initial):
    m = p["m_kg"].astype(np.longdouble)
    x = np.array([p[a] for a in ("x_m", "y_m", "z_m")], dtype=np.longdouble).T
    u = np.array([p[a] for a in ("ux", "uy", "uz")], dtype=np.longdouble).T
    u0 = np.array([initial[a] for a in ("ux", "uy", "uz")], dtype=np.longdouble).T
    g = np.sqrt(1 + np.sum(u*u, axis=1))
    g0 = np.sqrt(1 + np.sum(u0*u0, axis=1))
    mx = np.average(x, weights=m, axis=0)
    mu = np.average(u, weights=m, axis=0)
    dx, du = x-mx, u-mu
    return dict(delta_K_J=float(np.sum(m*C*C*np.sum((u-u0)*(u+u0), axis=1)/(g+g0))),
                var_x=np.asarray(np.average(dx*dx, weights=m, axis=0), float),
                cov_xu=np.asarray(np.average(dx*du, weights=m, axis=0), float))


def sampled_gaussian_force(p, meta):
    x = np.array([p[a] for a in ("x_m", "y_m", "z_m")]).T
    rest = x.copy()
    rest[:, 2] *= meta["gamma"]
    radius = np.linalg.norm(rest, axis=1)
    a = radius/meta["rest_sigma_m"]
    enclosed = erf(a/np.sqrt(2)) - np.sqrt(2/np.pi)*a*np.exp(-a*a/2)
    e = rest*(meta["charge_C"]*enclosed/(4*np.pi*EPS0*radius**3))[:, None]
    # Transverse E transforms by gamma; transverse Lorentz force by 1/gamma.
    e[:, :2] /= meta["gamma"]
    center = np.average(x, weights=p["m_kg"], axis=0)
    return np.sum((x-center)*p["q_C"][:, None]*e, axis=0)/(C*np.sum(p["m_kg"]))


def relative(a, b):
    return float(np.max(np.abs(np.asarray(a)-b))/max(float(np.max(np.abs(b))), 1e-300))


def analyze(root):
    plan = json.loads((root/"PLAN.json").read_text())
    limits = plan["predeclared_screening"]
    cases, functional = {}, []
    for path in sorted(root.glob("r[12]/*/metadata.json")):
        folder = path.parent
        meta = json.loads(path.read_text())
        h = read_csv(folder/"history.csv")
        p0, p1 = [read_csv(folder/f"{tag}-particles.csv") for tag in ("initial", "final")]
        checks = []
        def check(name, ok):
            checks.append(dict(name=name, passed=bool(ok)))
        check("runner_exit_zero", (folder/"EXIT_CODE").read_text().strip() == "0")
        check("bound_one_core", "[binding] result=PASS failures=0" in (folder/"binding.log").read_text())
        check("finite_history", all(np.isfinite(h[n]).all() for n in h.dtype.names))
        check("complete_history", len(h) == meta["steps"]+1 and h["step"][-1] == meta["steps"])
        check("particle_count_and_identity", len(p0) == meta["particles"] and np.array_equal(p0["id"], p1["id"]))
        pm = particle_moments(p1, p0)
        check("independent_delta_K", relative(pm["delta_K_J"], h["delta_K_J"][-1]) < 1e-10)
        for key, cols in (("var_x", ("var_x_m2", "var_y_m2", "var_z_m2")),
                          ("cov_xu", ("cov_x_ux_m", "cov_y_uy_m", "cov_z_uz_m"))):
            check("independent_"+key, relative(pm[key], np.array([h[n][-1] for n in cols])) < 1e-10)
        analytic = sampled_gaussian_force(p0, meta)
        check("independent_gaussian_force", relative(analytic, meta["analytic_sampled_force_covariance_m_per_s"]) < 1e-12)
        scale = meta["initial_field_scale"]
        check("expected_initial_Gauss", abs(h["Gauss_relative_L2"][0]-abs(1-scale)) < 1e-7)
        check("continuity", h["continuity_relative_rho_over_dt"].max() <= limits["continuity_relative_max"])
        if meta["mode"] == "zero" and meta["beta"] == 0:
            check("rest_zero_no_dynamics", np.array_equal(p0, p1) and np.all(h["delta_K_J"] == 0))
        if meta["mode"] == "fieldonly":
            check("fixed_source_no_particle_work", np.array_equal(p0, p1) and np.all(h["delta_K_J"] == 0) and np.all(h["JdotE_cumulative_J"] == 0))
        field = h["UE_J"]+h["UB_centered_J"]
        delta_u = field-field[0]
        closure = h["delta_K_J"]+delta_u
        exchange = max(float(np.max(np.abs(h["delta_K_J"]))), float(np.max(np.abs(delta_u))), 1e-300)
        gridforce = np.array(meta["grid_sampled_force_covariance_m_per_s"])
        continuum = meta["continuum_cov_xu_slope_m_per_s"]
        screening = {}
        if meta["mode"] == "qs":
            screening = dict(initial_gauss=bool(h["Gauss_relative_max"][0] <= limits["initial_qs_gauss_relative_max"]),
                gaussian_force=bool(np.max(np.abs(gridforce/continuum-1)) <= limits["rest_gaussian_initial_force_moment_relative"]),
                outward_expansion=bool(np.all(pm["cov_xu"] > 0)),
                energy_exchange_diagnostic=bool(np.max(np.abs(closure))/exchange <= limits["energy_exchange_residual_fraction_diagnostic"]))
        cases[folder.name] = dict(meta=meta, functional=checks,
            initial_particle_sha256=hashlib.sha256((folder/"initial-particles.csv").read_bytes()).hexdigest(),
            force_relative_to_continuum=(gridforce/continuum-1).tolist(),
            force_relative_to_sampled=(gridforce/analytic-1).tolist(),
            final_delta_K_J=float(h["delta_K_J"][-1]),
            final_delta_K_over_Uinf=float(h["delta_K_J"][-1]/meta["gaussian_rest_self_energy_J"]),
            sigma_growth_fraction=[float(np.sqrt(h[n][-1]/h[n][0])-1) for n in ("var_x_m2", "var_y_m2", "var_z_m2")],
            final_cov_xu=pm["cov_xu"].tolist(),
            initial_field_energy_J=float(field[0]), final_delta_field_J=float(delta_u[-1]),
            final_energy_balance_J=float(closure[-1]),
            max_energy_balance_over_exchange=float(np.max(np.abs(closure))/exchange),
            max_energy_balance_over_Uinf=float(np.max(np.abs(closure))/meta["gaussian_rest_self_energy_J"]),
            final_grid_work_J=float(h["JdotE_cumulative_J"][-1]),
            final_particle_minus_grid_work_J=float(h["delta_K_J"][-1]-h["JdotE_cumulative_J"][-1]),
            max_gauss_relative=float(h["Gauss_relative_max"].max()),
            max_continuity_relative=float(h["continuity_relative_rho_over_dt"].max()),
            screening=screening)
        functional.extend(dict(case=folder.name, **c) for c in checks)
    comparisons = {}
    for name, case in cases.items():
        if case["meta"]["mode"] != "fieldonly": continue
        counterpart = cases.get("qs-beta0-"+case["meta"]["grid"])
        functional.append(dict(case=name, name="same_initial_field_as_coupled_case", passed=bool(counterpart and case["initial_field_energy_J"] == counterpart["initial_field_energy_J"] and case["initial_particle_sha256"] == counterpart["initial_particle_sha256"])))
    for beta in (0, .23):
        chosen = {n:c for n,c in cases.items() if c["meta"]["beta"] == beta}
        hashes = {c["initial_particle_sha256"] for c in chosen.values()}
        functional.append(dict(case=str(beta), name="same_release_state_all_modes_grids", passed=len(hashes) == 1))
        base = next((c for c in chosen.values() if c["meta"]["mode"] == "qs" and c["meta"]["grid"] == "base"), None)
        if base is None: continue
        for name, other in chosen.items():
            if other["meta"]["mode"] != "qs" or other["meta"]["grid"] == "base": continue
            diffs = {key:relative(other[key], base[key]) for key in ("final_delta_K_J", "final_cov_xu", "sigma_growth_fraction")}
            comparisons[name] = dict(relative_to_base=diffs,
                delta_K_difference_over_Uinf=(other["final_delta_K_J"]-base["final_delta_K_J"])/other["meta"]["gaussian_rest_self_energy_J"],
                passed=all(v <= limits["qs_domain_and_resolution_observable_relative_change"] for v in diffs.values()))
    result = dict(scope=plan["scope"], functional_pass=all(c["passed"] for c in functional),
        functional_checks=functional, cases=cases, comparisons=comparisons,
        limits=limits, production_ready=False,
        interpretation="Energy closure is a convergent diagnostic, not an exact Cowan invariant. Fixed-N refinement is not particle-number convergence. No moving injection has been implemented.")
    (root/"REVIEW.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(functional_pass=result["functional_pass"], cases={n:dict(dK_pJ=c["final_delta_K_J"]*1e12, expansion=c["sigma_growth_fraction"][0], energy_residual_fraction=c["max_energy_balance_over_exchange"], screening=c["screening"]) for n,c in cases.items()}, comparisons=comparisons), indent=2))
    if not result["functional_pass"]: raise SystemExit(1)


if __name__ == "__main__":
    analyze(Path(sys.argv[1]).resolve())
