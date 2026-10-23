#!/usr/bin/env python3
"""Prepare a 4 m APRL run from a streamed Elegant HDF5 plane file."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    raise SystemExit("需要 h5py 和 numpy；请在含这两个包的 Python 环境运行 submit.sh") from exc


C = 299_792_458.0
ROOT = Path(__file__).resolve().parent
GAMMA = 1174.169835810956
UNDULATOR_K = 1.092714657005
PERIOD_M = 0.0234
PERIODS = 170
FRINGE_CUTOFF = 1.0e-9
DZ_LAB_M = 1.35e-9
MESH_TRANSVERSE_CELLS = 104
MESH_TRANSVERSE_CELL_M = 5e-6
EXIT_PARTICLE_Z_M = 4.03
INTERMEDIATE_PARTICLE_Z_M = (1.0, 2.0, 3.0)
FIELD_RHYTHM_S = DZ_LAB_M / C
BOOST_GAMMA = GAMMA / math.sqrt(1.0 + 0.5 * UNDULATOR_K**2)
FRINGE_M = math.sqrt(-2.0 * math.log(FRINGE_CUTOFF)) / (2.0 * math.pi / PERIOD_M)
INTERACTION_START_M = -FRINGE_M
INTERACTION_END_M = PERIODS * PERIOD_M + FRINGE_M
ONE_CELL_MARGIN_M = BOOST_GAMMA**2 * DZ_LAB_M
INPUT_PLANE_Z_M = -0.05  # Elegant WATCH: 50 mm before nominal RAD entrance.
PROJECTION_CLEARANCE_M = 20e-6
CHUNK_ROWS = 262_144


def quoted(path: Path) -> str:
    return json.dumps(str(path), ensure_ascii=False)


def read_particle_plan(path: Path) -> dict[str, float | int | bool]:
    if not path.is_file():
        raise ValueError(f"未找到输入文件: {path}")
    before = path.stat()
    with h5py.File(path, "r") as source:
        group = source.get("/particles")
        if group is None or int(group.attrs.get("format_version", -1)) != 4:
            raise ValueError("particles.h5 必须是 Elegant 转换器生成的 HDF5 v4 过面事件文件")
        records = group.get("records")
        required = {"plane_position_m", "arrival_time_offset_s",
                    "proper_velocity", "source_id", "macro_weight"}
        if records is None or records.ndim != 1 or not required.issubset(records.dtype.names or ()):
            raise ValueError("/particles/records 缺少 v4 必需字段")
        count = len(records)
        if count == 0:
            raise ValueError("/particles/records 是空的")
        if records.dtype.fields["proper_velocity"][0].shape != (3,):
            raise ValueError("proper_velocity 必须是长度为 3 的向量")
        if records.dtype.fields["plane_position_m"][0].shape != (2,):
            raise ValueError("plane_position_m 必须是长度为 2 的向量")

        weight_sum = 0.0
        beta_weight_sum = 0.0
        beta_time_weight_sum = 0.0
        earliest = math.inf
        latest = -math.inf
        maximum_gamma = 1.0
        all_unit_weights = True
        selected = records.fields(["arrival_time_offset_s", "proper_velocity", "macro_weight"])
        for begin in range(0, count, CHUNK_ROWS):
            block = selected[begin:min(begin + CHUNK_ROWS, count)]
            tau = block["arrival_time_offset_s"]
            proper_velocity = block["proper_velocity"]
            weight = block["macro_weight"]
            uz = proper_velocity[:, 2]
            if (not np.isfinite(tau).all() or not np.isfinite(proper_velocity).all()
                    or not np.isfinite(weight).all() or np.any(weight <= 0)
                    or np.any(uz <= 0)):
                raise ValueError(f"粒子记录 {begin + 1} 起存在无效时间、速度或权重")
            gamma = np.sqrt(1.0 + np.sum(proper_velocity**2, axis=1))
            if not np.isfinite(gamma).all():
                raise ValueError(f"粒子记录 {begin + 1} 起存在溢出的 Lorentz γ")
            beta_z = uz / gamma
            maximum_gamma = max(maximum_gamma, float(np.max(gamma)))
            weight_sum += float(np.sum(weight, dtype=np.float64))
            beta_weight_sum += float(np.sum(weight * beta_z, dtype=np.float64))
            beta_time_weight_sum += float(np.sum(weight * beta_z * tau, dtype=np.float64))
            earliest = min(earliest, float(np.min(tau)))
            latest = max(latest, float(np.max(tau)))
            all_unit_weights &= bool(np.all(weight == 1.0))

    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("particles.h5 在检查期间发生变化；请等转换完成后再提交")

    mean_beta = beta_weight_sum / weight_sum
    mean_beta_time = beta_time_weight_sum / weight_sum
    # APRL requires T >= every plane-crossing time. Advance just enough to
    # include the latest record, plus a small, explicit numerical clearance.
    minimum_advance = C * (mean_beta * latest - mean_beta_time)
    advance = max(0.0, minimum_advance) + PROJECTION_CLEARANCE_M
    center = INPUT_PLANE_Z_M + advance
    common_time = (advance / C + mean_beta_time) / mean_beta
    if common_time <= latest:
        raise ValueError("计算出的共同时间仍早于最晚过面粒子")

    # z_i(T) = z_plane + beta_z*c*(T-tau_i). Since beta_z < 1, this bound
    # safely overestimates the head without a second full-file read.
    head_upper_bound = INPUT_PLANE_Z_M + C * (common_time - earliest)
    head_limit = INTERACTION_START_M - ONE_CELL_MARGIN_M
    if head_upper_bound >= head_limit:
        raise ValueError(
            f"束团头部上界 {head_upper_bound:.9g} m 已接近边缘场；"
            f"必须小于 {head_limit:.9g} m。请增加 WATCH 到波荡器的距离。"
        )

    # APRL puts the particle-background reference at field_z - gamma_max*rho,
    # where rho is the diagonal of the *full* transverse mesh width. Put that
    # reference exactly at the former exit-particle detector instead of writing
    # a second, redundant full-bunch particle file at the same location.
    reference_rho = math.hypot(
        MESH_TRANSVERSE_CELLS * MESH_TRANSVERSE_CELL_M,
        MESH_TRANSVERSE_CELLS * MESH_TRANSVERSE_CELL_M,
    )
    field_z = EXIT_PARTICLE_Z_M + maximum_gamma * reference_rho
    if EXIT_PARTICLE_Z_M <= INTERACTION_END_M + ONE_CELL_MARGIN_M:
        raise ValueError("出口参考粒子面过于接近波荡器边缘场")

    requested_electrons = os.environ.get("ELECTRONS", "").strip()
    if requested_electrons:
        electrons = float(requested_electrons)
        if not math.isfinite(electrons) or electrons <= 0:
            raise ValueError("ELECTRONS 必须是正的物理电子总数")
    elif all_unit_weights:
        electrons = float(count)
    else:
        raise ValueError("输入含非单位相对权重；请显式设置 ELECTRONS=物理电子总数")

    return {
        "records": count,
        "electrons": electrons,
        "all_unit_weights": all_unit_weights,
        "earliest_arrival_s": earliest,
        "latest_arrival_s": latest,
        "input_plane_z_m": INPUT_PLANE_Z_M,
        "initial_center_z_m": center,
        "common_time_after_reference_s": common_time,
        "head_upper_bound_z_m": head_upper_bound,
        "head_limit_z_m": head_limit,
        "interaction_start_z_m": INTERACTION_START_M,
        "one_cell_margin_m": ONE_CELL_MARGIN_M,
        "maximum_particle_gamma": maximum_gamma,
        "ballistic_reference_z_m": EXIT_PARTICLE_Z_M,
        "field_plane_z_m": field_z,
        "field_reference_distance_m": maximum_gamma * reference_rho,
        "field_rhythm_s": FIELD_RHYTHM_S,
    }


def make_card(input_file: Path, run_dir: Path, plan: dict[str, float | int | bool]) -> str:
    photon_energy_ev = 1.239841984e-6 / (
        PERIOD_M * (1.0 + 0.5 * UNDULATOR_K**2) / (2.0 * GAMMA**2)
    )
    return f"""# Generated from the original 4 m dz1p35-B0768 card.
# Elegant WATCH is 50 mm upstream of the nominal undulator entrance z=0.
units:
  length: m
  time: s

runtime:
  mode: throughput
  resource_monitor:
    enabled: true
    progress_interval_steps: 50
    calibration_steps: 5
    memory_safety_factor: 1.25
    time_safety_factor: 1.25
  shutdown_watchdog:
    enabled: true
    timeout_seconds: 600
    finalize_timeout_seconds: 60

output:
  overwrite: false
  manifest: {quoted(run_dir / 'run-manifest.yaml')}

mesh:
  field_solver: cowan-z
  current_filter: none
  cells: [{MESH_TRANSVERSE_CELLS}, {MESH_TRANSVERSE_CELLS}, auto]
  cell_size: [{MESH_TRANSVERSE_CELL_M:.17g}, {MESH_TRANSVERSE_CELL_M:.17g}, {BOOST_GAMMA * DZ_LAB_M:.17g}]
  center: [0.0, 0.0, 0.0]
  duration: 5.0e-10
  boost_gamma: {BOOST_GAMMA:.17g}
  particle_steps_per_undulator_period: 80
  maximum_particle_substeps: 4096
  longitudinal:
    mode: observation-bounded
    lower_guard_cells: 64
    upper_guard_cells: 128
    self_field_padding_rms: 0.0
    ballistic_time_safety_factor: 1.1
    alignment_cells: 256

radiation_resolution:
  enabled: true
  maximum_photon_energy_eV: {1.25 * photon_energy_ev:.17g}
  warning_grid_points_per_wavelength: 8.0
  warning_maxwell_samples_per_cycle: 8.0
  warning_detector_samples_per_cycle: 8.0

boundary:
  type: cpml
  cells: [10, 10, 20]
  polynomial_order: 3.0
  target_reflection: 1.0e-8
  kappa_max: 8.0
  alpha_fraction: 0.0

initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-8
  maximum_mean_transverse_beta: 1.0e-7
  maximum_iterations: 10000

beam:
  reference:
    input_plane_z: {plan['input_plane_z_m']:.17g}
    initial_center_z: {plan['initial_center_z_m']:.17g}
  input:
    type: hdf5
    file: {quoted(input_file)}
    electrons: {plan['electrons']:.17g}
    position_offset: [0.0, 0.0, 0.0]

sources:
  incident_waves: []
  magnetic_elements:
    - type: planar-undulator
      characteristic: true
      strength_parameter: {UNDULATOR_K:.17g}
      period: {PERIOD_M:.17g}
      periods: {PERIODS}
      entrance_z: 0.0
      polarization_angle_rad: 1.5707963267948966
      gaussian_fringe: true
      fringe_relative_cutoff: {FRINGE_CUTOFF:.17g}

detectors:
  enabled: true
  directory: {quoted(run_dir / 'detectors')}
  field_planes:
    - name: exit-field
      z: {plan['field_plane_z_m']:.17g}
      rhythm: {FIELD_RHYTHM_S:.17g}
      buffer_samples: 8
      compression: 0
      particle_background:
        enabled: true
        buffer_records: 65536
        compression: 0
        validation:
          enabled: false
  particle_planes:
    - name: entry-particles
      z: -0.03
      buffer_records: 65536
      compression: 0
    - name: after-1m-particles
      z: {INTERMEDIATE_PARTICLE_Z_M[0]:.17g}
      buffer_records: 65536
      compression: 0
    - name: after-2m-particles
      z: {INTERMEDIATE_PARTICLE_Z_M[1]:.17g}
      buffer_records: 65536
      compression: 0
    - name: after-3m-particles
      z: {INTERMEDIATE_PARTICLE_Z_M[2]:.17g}
      buffer_records: 65536
      compression: 0

stop:
  mode: after-last-element

trajectory:
  enabled: false

energy_ledger:
  enabled: true
  directory: {quoted(run_dir / 'ledger')}
  filename: energy-ledger.h5
  sample_interval_steps: 10
  buffer_records: 128
  compression: 0
  warning_relative_tolerance: 0.05
"""


def make_reconstruction_card(run_dir: Path) -> str:
    return f"""# Run later with field_reconstruction; the APRL job keeps the raw field.
input:
  field_file: {quoted(run_dir / 'detectors' / 'exit-field.h5')}
  particle_file: {quoted(run_dir / 'detectors' / 'exit-field-ballistic-reference.h5')}
  require_complete: true
  particle_read_chunk: 65536
model:
  fft_threads: 1
  padding_factor: [1, 1, 1]
  transverse_smoothing_m: 0.0
  maximum_relative_gamma_spread: 0.05
  maximum_rms_transverse_beta: 0.05
  maximum_outside_charge_fraction: 0.01
output:
  file: {quoted(run_dir / 'detectors' / 'exit-field-reconstructed.h5')}
  compression: 0
  overwrite: false
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    input_file = ROOT / "particles.h5"
    run_dir = args.run_dir.expanduser().resolve()
    if run_dir.exists() and any(run_dir.iterdir()):
        raise SystemExit(f"运行目录非空，拒绝覆盖: {run_dir}")
    try:
        plan = read_particle_plan(input_file)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"输入检查失败: {exc}") from exc
    run_dir.mkdir(parents=True, exist_ok=True)
    stat = input_file.stat()
    plan.update({
        "input_file": str(input_file),
        "input_size_bytes": stat.st_size,
        "input_mtime_ns": stat.st_mtime_ns,
        "run_dir": str(run_dir),
    })
    (run_dir / "run.yaml").write_text(make_card(input_file, run_dir, plan), encoding="utf-8")
    (run_dir / "reconstruction.yaml").write_text(
        make_reconstruction_card(run_dir), encoding="utf-8"
    )
    (run_dir / "input-plan.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"粒子记录: {plan['records']:,}; 物理电子数: {plan['electrons']:.17g}")
    print(f"Elegant 输出面 z={INPUT_PLANE_Z_M:.9g} m; 波荡器名义入口 z=0 m")
    print(f"APRL 初始中心 z={plan['initial_center_z_m']:.9g} m")
    print(f"边缘场起点 z={INTERACTION_START_M:.9g} m; 一格安全距离={ONE_CELL_MARGIN_M:.9g} m")
    print(f"中途粒子面 z={', '.join(f'{z:g}' for z in INTERMEDIATE_PARTICLE_Z_M)} m")
    print(f"出口参考粒子面 z={EXIT_PARTICLE_Z_M:.9g} m; E/B 场面 z={plan['field_plane_z_m']:.9g} m")
    print(f"E/B 采样间隔={FIELD_RHYTHM_S:.9g} s; 原始场与参考粒子保留，重构留待后处理")
    print(f"输入卡: {run_dir / 'run.yaml'}")
    print(f"后处理卡: {run_dir / 'reconstruction.yaml'}")


if __name__ == "__main__":
    main()
