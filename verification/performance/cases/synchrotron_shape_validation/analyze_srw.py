"""Independent SRW check of the finite-dipole on-axis spectrum.

Run this with the Python interpreter from the micromamba ``srw`` environment.
SRW's bending-magnet single-electron result is kept separate from the small
finite transverse-size correction imported from the converged LW quadrature.
"""
from array import array
import csv
import json
from pathlib import Path
import platform
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar

import srwpy
from srwpy.srwlib import (SRWLPartBeam, SRWLPrtTrj, SRWLWfr,
                          SRWLMagFldC, SRWLMagFldM, srwl)


SHAPES = ["linear", "quintic-compensated", "linear-notch",
          "linear-binomial3"]
LABELS = {"linear": "linear", "quintic-compensated": "quintic",
          "linear-notch": "7-point", "linear-binomial3": "3-point"}
COLORS = {"linear": "#4c78a8", "quintic-compensated": "#f58518",
          "linear-notch": "#54a24b", "linear-binomial3": "#b279a2"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def relative_l2(value, reference):
    return float(np.linalg.norm(value - reference) /
                 max(np.linalg.norm(reference), 1e-300))


def make_beam(gamma, current):
    beam = SRWLPartBeam()
    beam.Iavg = current
    beam.partStatMom1.x = 0.
    beam.partStatMom1.y = 0.
    beam.partStatMom1.z = 0.
    beam.partStatMom1.xp = 0.
    beam.partStatMom1.yp = 0.
    beam.partStatMom1.gamma = gamma
    return beam


def make_magnet(field, length):
    # Centre the SRW hard-edge magnet so that it occupies 0 <= z <= length.
    element = SRWLMagFldM(field, 1, "n", length, 0.)
    return SRWLMagFldC([element], [0.], [0.], [length / 2])


def srw_point_spectrum(energies, p, margin, ntrajectory, current=1.):
    """Return raw SRW intensity and per-electron spectral fluence.

    CalcIntFromElecField type 0 is reported by SRW as
    ph/s/0.1%bw/mm^2.  Dividing by I/e, by the 0.001*E bandwidth,
    converting mm^-2 to m^-2, and multiplying by E*e gives the energy-
    independent factor e^2*1e9/I for J/(eV m^2) per electron passage.
    """
    beam = make_beam(p["gamma"], current)
    magnet = make_magnet(p["field_T"], p["length"])
    wavefront = SRWLWfr()
    wavefront.allocate(len(energies), 1, 1)
    wavefront.mesh.zStart = p["screen"]
    wavefront.mesh.eStart = float(energies[0])
    wavefront.mesh.eFin = float(energies[-1])
    wavefront.mesh.xStart = wavefront.mesh.xFin = 0.
    wavefront.mesh.yStart = wavefront.mesh.yFin = 0.
    wavefront.partBeam = beam
    precision = [0, 1e-6, -margin, p["length"] + margin,
                 int(ntrajectory), 1, -1]
    srwl.CalcElecFieldSR(wavefront, 0, magnet, precision)
    raw = array("f", [0.] * len(energies))
    srwl.CalcIntFromElecField(raw, wavefront, 6, 0, 0,
                             float(energies[0]), 0., 0.)
    raw = np.asarray(raw, dtype=float)
    require(np.all(np.isfinite(raw)) and np.all(raw > 0),
            "SRW returned non-positive or non-finite intensity")
    return raw, raw * e**2 * 1e9 / current


def srw_exit_angle(p, margin):
    beam = make_beam(p["gamma"], 1.)
    trajectory = SRWLPrtTrj(_np=40001, _ctStart=0.,
                            _ctEnd=p["length"] + margin,
                            _partInitCond=beam.partStatMom1)
    srwl.CalcPartTraj(trajectory, make_magnet(p["field_T"], p["length"]), [1])
    return float(np.arctan2(trajectory.arXp[-1], trajectory.arZp[-1]))


def write_csv(path, columns):
    names = list(columns)
    rows = zip(*(columns[name] for name in names))
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        writer.writerows(rows)


def main(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    require(not output.exists(), "SRW output already exists")
    output.mkdir(parents=True)

    scripts = root / "scripts"
    sys.path.insert(0, str(scripts))
    from theory import arc_spectrum, constants, energy_density_spectrum
    from theory_paired import gaussian_paired_reference

    campaign = json.loads((root / "campaign.json").read_text())
    validation = json.loads((root / "axis_analysis_r2/AXIS_VALIDATION.json").read_text())
    require(validation["passed"] and validation["frames_per_screen"] == 278,
            "Closed axis result is not the validated 278-frame result")
    common = campaign["common"]
    supplied = common["parameters"]
    p = constants(gamma=supplied["gamma"], rho=supplied["rho"],
                  length=supplied["length"], screen=supplied["screen"])
    for key in ["gamma", "beta", "rho", "length", "screen", "field_T"]:
        require(np.isclose(p[key], supplied[key], rtol=3e-15, atol=0.),
                f"Parameter mismatch for {key}")

    table = np.genfromtxt(root / "axis_analysis_r2/axis_spectra.csv",
                          delimiter=",", names=True)
    energies = np.asarray(table["photon_eV"])
    require(len(energies) == 10 and np.allclose(np.diff(energies),
            np.diff(energies)[0], rtol=2e-14, atol=0.),
            "SRW requires the existing uniform 10-bin energy grid")
    delta_energy = float(energies[1] - energies[0])
    time_span = h / (delta_energy * e)
    first = common["screen_first_relative_u_s"]
    frames = validation["frames_per_screen"]
    window = dict(time_first=first, time_span=time_span, ntime=8 * frames)

    # Independent SRW point-electron result and three numerical checks.
    margin = 1.5e-3
    raw, srw_point = srw_point_spectrum(energies, p, margin, 100000, 1.)
    raw_half, srw_half = srw_point_spectrum(energies, p, margin, 100000, .5)
    _, srw_margin_low = srw_point_spectrum(energies, p, 1e-3, 100000, 1.)
    _, srw_margin_high = srw_point_spectrum(energies, p, 2e-3, 100000, 1.)
    _, srw_ntrajectory_low = srw_point_spectrum(energies, p, margin, 20000, 1.)
    current_scaling_error = float(np.max(abs(raw_half / raw - .5)))
    converted_current_error = relative_l2(srw_half, srw_point)
    margin_envelope = float(max(relative_l2(srw_margin_low, srw_point),
                                relative_l2(srw_margin_high, srw_point)))
    ntrajectory_error = relative_l2(srw_ntrajectory_low, srw_point)
    require(current_scaling_error < 2e-6 and converted_current_error < 2e-6,
            "SRW current/unit scaling check failed")
    require(margin_envelope < 2e-3, "SRW endpoint margin is not converged")
    require(ntrajectory_error < 1e-2, "SRW trajectory sampling is not converged")

    # Two independent LW point-electron references: finite-record paired and
    # acceleration-field-only.  Neither calls SRW or APRL.
    point_values = gaussian_paired_reference(
        0., 0., energies, p, sigma_xy=0., sigma_z=0., order=1, **window)
    lw_point_paired = energy_density_spectrum(point_values[0], point_values[1])
    arc_e, arc_b = arc_spectrum(np.array([[0., 0., p["screen"]]]),
                                energies, p, ntheta=1024)
    lw_point_arc = energy_density_spectrum(arc_e[:, 0], arc_b[:, 0])

    # Map the single-electron result to the deterministic coherent Gaussian
    # continuum.  The longitudinal factor is analytic.  The small transverse
    # correction is measured from the converged LW quadrature and kept as a
    # separately reported hybrid factor, never presented as native SRW ME output.
    sigma_z = common["sigma_z_m"]
    sigma_xy = common["sigma_x_m"]
    represented = common["represented_electrons"]
    omega = energies * e / hbar
    form2 = np.exp(-(omega * sigma_z / (p["beta"] * c))**2)
    longitudinal_values = gaussian_paired_reference(
        0., 0., energies, p, sigma_xy=0., sigma_z=sigma_z,
        order=1, **window)
    lw_point_longitudinal = energy_density_spectrum(
        longitudinal_values[0], longitudinal_values[1])
    lw_full = np.asarray(table["theory_J_per_eV_m2"])
    transverse_correction = lw_full / (represented**2 * lw_point_longitudinal)
    numeric_longitudinal_factor = lw_point_longitudinal / lw_point_paired
    longitudinal_factor_error = relative_l2(numeric_longitudinal_factor, form2)
    require(np.all((transverse_correction > .95) &
                   (transverse_correction < 1.05)),
            "LW transverse correction is too large for the hybrid mapping")
    srw_gaussian_point = srw_point * represented**2 * form2
    srw_gaussian_hybrid = srw_gaussian_point * transverse_correction

    exit_angle = srw_exit_angle(p, margin)
    angle_error = abs(exit_angle) - p["theta"]
    require(abs(angle_error) < 2e-6, "SRW orbit does not match the target bend")

    single_columns = {
        "photon_eV": energies,
        "SRW_J_per_eV_m2_per_electron": srw_point,
        "LW_paired_J_per_eV_m2_per_electron": lw_point_paired,
        "LW_acceleration_J_per_eV_m2_per_electron": lw_point_arc,
        "SRW_margin_1mm": srw_margin_low,
        "SRW_margin_2mm": srw_margin_high,
        "SRW_ntrajectory_20000": srw_ntrajectory_low,
    }
    write_csv(output / "single_electron_comparison.csv", single_columns)

    gaussian_columns = {
        "photon_eV": energies,
        "SRW_point_scaled_J_per_eV_m2": srw_gaussian_point,
        "SRW_hybrid_transverse_corrected_J_per_eV_m2": srw_gaussian_hybrid,
        "finite_arc_LW_J_per_eV_m2": lw_full,
        "LW_transverse_correction": transverse_correction,
        "longitudinal_form_factor_squared": form2,
    }
    for shape in SHAPES:
        # numpy.genfromtxt deletes punctuation from CSV field names.
        gaussian_columns[shape] = np.asarray(table[shape.replace("-", "")])
    write_csv(output / "gaussian_comparison.csv", gaussian_columns)

    band = ((energies >= common["registered_band_eV"][0]) &
            (energies <= common["registered_band_eV"][1]))
    def integrate(values):
        return float(np.sum(values[band]) * delta_energy)

    srw_band = integrate(srw_gaussian_hybrid)
    summary = [{
        "curve": "finite-arc LW",
        "band_J_per_m2": integrate(lw_full),
        "band_over_SRW": integrate(lw_full) / srw_band,
        "registered_spectral_relative_L2_to_SRW":
            relative_l2(lw_full[band], srw_gaussian_hybrid[band]),
    }]
    for shape in SHAPES:
        spectrum = gaussian_columns[shape]
        summary.append({
            "curve": shape,
            "band_J_per_m2": integrate(spectrum),
            "band_over_SRW": integrate(spectrum) / srw_band,
            "registered_spectral_relative_L2_to_SRW":
                relative_l2(spectrum[band], srw_gaussian_hybrid[band]),
        })
    with (output / "comparison_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    # A compact four-panel view keeps the independent and hybrid comparisons
    # visually separate.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    axes[0, 0].plot(energies, srw_point * 1e16, "o-", lw=2,
                    label="SRW 4.2.1")
    axes[0, 0].plot(energies, lw_point_paired * 1e16, "s--",
                    label="finite-record paired LW")
    axes[0, 0].plot(energies, lw_point_arc * 1e16, ":", lw=2,
                    label="acceleration-only LW")
    axes[0, 0].set(title="Independent point-electron check",
                   ylabel=r"$10^{-16}$ J / (eV m$^2$)")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(energies, srw_point / lw_point_paired, "o-",
                    label="SRW / paired LW")
    axes[0, 1].plot(energies, srw_point / lw_point_arc, "s--",
                    label="SRW / acceleration LW")
    axes[0, 1].axhline(1., color="black", lw=1)
    axes[0, 1].set(title="Point-electron absolute ratios", ylabel="ratio")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(energies, lw_full * 1e12, "k--", lw=2,
                    label="finite-arc LW")
    axes[1, 0].plot(energies, srw_gaussian_point * 1e12, color="#e45756",
                    ls=":", lw=2, label="SRW + longitudinal factor")
    axes[1, 0].plot(energies, srw_gaussian_hybrid * 1e12,
                    color="#e45756", marker="o", lw=2,
                    label="SRW + transverse correction")
    for shape in SHAPES:
        axes[1, 0].plot(energies, gaussian_columns[shape] * 1e12,
                        color=COLORS[shape], alpha=.68, marker=".",
                        label=LABELS[shape])
    axes[1, 0].set(title="Coherent Gaussian bunch comparison",
                   ylabel=r"pJ / (eV m$^2$)")
    axes[1, 0].legend(fontsize=7, ncol=2)

    x = np.arange(len(summary))
    axes[1, 1].bar(x, [row["band_over_SRW"] for row in summary],
                   color=["#777777"] + [COLORS[s] for s in SHAPES])
    axes[1, 1].axhline(1., color="black", lw=1)
    axes[1, 1].set_xticks(x, ["LW"] + [LABELS[s] for s in SHAPES],
                          rotation=15)
    axes[1, 1].set(title="Integrated 8–12 eV absolute fluence",
                   ylabel="curve / SRW")
    for axis in axes.flat:
        axis.set_xlabel("Photon energy (eV)" if axis is not axes[1, 1] else "")
        axis.grid(alpha=.2)
    for axis in [axes[1, 0], axes[1, 1]]:
        axis.axvspan(*common["registered_band_eV"], color="grey", alpha=.10) \
            if axis is axes[1, 0] else None
    fig.suptitle("SRW cross-check of the finite 2.4 mm dipole; no fitted normalization")
    fig.savefig(output / "srw_comparison.png", dpi=200)
    fig.savefig(output / "srw_comparison.pdf")
    plt.close(fig)

    final = {
        "complete": True,
        "status": "INDEPENDENT_SRW_CROSS_CHECK_SINGLE_PIC_GRID",
        "srw": {
            "python": platform.python_version(),
            "module_path": str(Path(srwpy.__file__).resolve()),
            "environment_python": str(Path(sys.executable).resolve()),
            "calculation": "SRWL CalcElecFieldSR, manual method, terminating terms enabled",
            "raw_intensity_unit": "ph/s/0.1%bw/mm^2",
            "per_electron_conversion": "raw * elementary_charge^2 * 1e9 / Iavg",
            "integration_margin_m": margin,
            "trajectory_points": 100000,
            "exit_angle_rad": exit_angle,
            "absolute_exit_angle_minus_target_rad": angle_error,
        },
        "checks": {
            "raw_half_ampere_over_one_ampere_max_error_from_half": current_scaling_error,
            "converted_half_vs_one_ampere_relative_L2": converted_current_error,
            "one_vs_one_point_five_vs_two_mm_margin_relative_L2_envelope": margin_envelope,
            "ntrajectory_20000_vs_100000_relative_L2": ntrajectory_error,
            "numeric_vs_analytic_longitudinal_factor_relative_L2": longitudinal_factor_error,
            "LW_transverse_correction_min": float(np.min(transverse_correction)),
            "LW_transverse_correction_max": float(np.max(transverse_correction)),
        },
        "single_electron": {
            "SRW_vs_finite_record_paired_LW_relative_L2":
                relative_l2(srw_point, lw_point_paired),
            "SRW_vs_acceleration_only_LW_relative_L2":
                relative_l2(srw_point, lw_point_arc),
            "pointwise_SRW_over_paired_LW_min":
                float(np.min(srw_point / lw_point_paired)),
            "pointwise_SRW_over_paired_LW_max":
                float(np.max(srw_point / lw_point_paired)),
        },
        "coherent_gaussian": {
            "represented_electrons": represented,
            "sigma_xy_m": sigma_xy,
            "sigma_z_m": sigma_z,
            "SRW_hybrid_band_J_per_m2": srw_band,
            "finite_arc_LW_band_J_per_m2": integrate(lw_full),
            "finite_arc_LW_over_SRW": integrate(lw_full) / srw_band,
            "curves": summary,
            "mapping_note": ("Native SRW point-electron spectrum times analytic longitudinal "
                             "form factor and N^2; the separately tabulated 0.98-scale "
                             "transverse factor comes from converged finite-arc LW quadrature."),
        },
        "normalization": "none; all spectra and band integrals are absolute",
        "physical_accuracy_accepted": False,
        "acceptance_note": ("SRW independently checks the radiation reference, but the APRL "
                            "comparison still uses one PIC mesh and time step."),
    }
    dump(output / "FINAL.json", final)

    curve_rows = {row["curve"]: row for row in summary}
    review = f"""# SRW 独立同步辐射对照

使用已安装的 SRW Python 环境，对同一个 2.4 mm 硬边偶极磁铁、同一个入口切线
和同一组 10 个频率点进行了独立计算。所有结果都是绝对量，没有幅度归一化。

## 结果

- 单电子 SRW 与有限记录配对 LW 谱的相对 L2 差为
  `{final['single_electron']['SRW_vs_finite_record_paired_LW_relative_L2']:.6g}`；逐点
  SRW/LW 范围为
  `{final['single_electron']['pointwise_SRW_over_paired_LW_min']:.6g}`--
  `{final['single_electron']['pointwise_SRW_over_paired_LW_max']:.6g}`。
- 将 SRW 单电子谱乘以解析纵向高斯形状因子和 N^2 后，再施加单独列出的横向
  修正 `{np.min(transverse_correction):.6g}`--`{np.max(transverse_correction):.6g}`，
  得到 8--12 eV 轴上通量 `{srw_band:.9e} J/m2`。
- 原有限轨道 LW 为 `{integrate(lw_full):.9e} J/m2`，即 SRW 的
  `{integrate(lw_full) / srw_band:.6g}` 倍。二者在该带内给出相近但不完全相同的
  独立理论基准。

| 曲线 | 8--12 eV (J/m2) | 曲线/SRW | 带内谱形 L2（相对 SRW） |
|---|---:|---:|---:|
"""
    for row in summary:
        review += (f"| {row['curve']} | {row['band_J_per_m2']:.9e} | "
                   f"{row['band_over_SRW']:.6g} | "
                   f"{row['registered_spectral_relative_L2_to_SRW']:.6g} |\n")
    review += f"""

## 数值检查与限制

SRW 原始强度对 0.5 A/1 A 的线性缩放误差为 `{current_scaling_error:.3e}`；转换为
单电子通量后两者相对 L2 为 `{converted_current_error:.3e}`。积分端点从 1.0/1.5/
2.0 mm 变化的最大相对 L2 为 `{margin_envelope:.3e}`，轨迹点从 20,000 提高到
100,000 的相对变化为 `{ntrajectory_error:.3e}`。SRW 出口偏转角绝对值与目标角
相差 `{angle_error:.3e} rad`。

SRW 原生多电子强度卷积是非相干处理，不能直接代表这里的确定性相干高斯束团；
因此主记录同时保留纯单电子 SRW、只乘解析纵向因子的曲线，以及带有明确 LW
横向修正的混合曲线。横向修正不足 2%，但它仍不是 SRW 原生结果。这个检查支持
“辐射理论基准本身大体正确”，不把单网格 APRL 结果升级为物理精度验收；后者仍
需要 dz/dt 和横向网格扫描。
"""
    (output / "RESULT_REVIEW.md").write_text(review)
    print(json.dumps(final, indent=2), flush=True)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: analyze_srw.py CAMPAIGN_ROOT OUTPUT")
    main(sys.argv[1], sys.argv[2])
