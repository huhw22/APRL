#!/usr/bin/env python3
"""Bounded algebra check for a transversely modal Cowan reference solve.

This is NOT PIC, a production field supplier, or a performance benchmark.
The dense 3-D DST oracle is confined to a fixed small test box. The candidate
accepts source blocks, stores forward boundary summaries and reconstructs
solution blocks in reverse order; it never owns a full longitudinal field.
Summaries currently remain in RAM and grow with length/block size. Moving
them to disk, particle-to-source streaming, and a physical entrance boundary
are deliberately not implemented here.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

# Set before importing numerical libraries; no process/system configuration is
# changed. scipy.fft calls additionally request one worker explicitly.
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
import scipy
from scipy.fft import dstn


REPO = Path(__file__).resolve().parents[4]
WORK = REPO.parent.parent
DEFAULT_METRICS = WORK / 'results/held_quasisteady_20261028/cowan/cowan0/metrics.json'
NX, NY, NZ = 12, 14, 80  # Interior vertices, not cell counts.
BLOCK_SIZES = (1, 7, 16, 33)
FIELD_TOLERANCE = 5.0e-11
RESIDUAL_TOLERANCE = 5.0e-11
SYMBOL_TOLERANCE = 5.0e-14
NEGATIVE_CONTROL_MINIMUM = 1.0e-7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def squared_norm(values):
    values = np.asarray(values, dtype=np.longdouble)
    return np.sum(values * values, dtype=np.longdouble)


def relative_norm(difference, reference):
    denominator = squared_norm(reference)
    require(denominator > 0, 'Zero reference norm in bounded nonzero fixture')
    return float(np.sqrt(squared_norm(difference) / denominator))


def read_constants(path):
    text = Path(path).read_text()
    result = {}
    for name in ('c', 'epsilon0', 'mu0'):
        match = re.search(r'const Double ' + name + r'\s*=\s*([0-9.eE+\-]+)\s*;', text)
        require(match is not None, 'Cannot identify original SI constant ' + name)
        result[name] = float(match.group(1))
    return result


def cowan_coefficients(spacing):
    """Literal scalar formulas from EBCowanCoefficients::forZDispersion."""
    dx, dy, dz = spacing
    require(min(spacing) > 0 and dx >= dz and dy >= dz, 'Invalid production spacing')
    rx, ry, rz = (dz / dx) ** 2, (dz / dy) ** 2, 1.0
    denominator = rx * ry + ry * rz + rz * rx
    fraction = rx * ry * rz / denominator
    beta = [r * (1.0 - fraction) / 8.0 for r in (rx, ry, rz)]
    delta_xy = rx * ry * (1.0 / 16.0 - rx * ry / (8.0 * denominator))
    delta_yz = ry * rz * (1.0 / 16.0 - ry * rz / (8.0 * denominator))
    delta_zx = rz * rx * (1.0 / 16.0 - rz * rx / (8.0 * denominator))
    alpha = [1.0 - 2.0 * beta[1] - 2.0 * beta[2] - 4.0 * delta_yz,
             1.0 - 2.0 * beta[2] - 2.0 * beta[0] - 4.0 * delta_zx,
             1.0 - 2.0 * beta[0] - 2.0 * beta[1] - 4.0 * delta_xy]
    return {'alpha': alpha, 'beta': beta,
            'delta_by_derivative_axis': [delta_yz, delta_zx, delta_xy],
            'ratio': [rx, ry, rz]}


def stencil27(spacing, coefficient, translation):
    """Literal 27-point assembly from the production CowanPoissonOperator."""
    result = np.zeros((3, 3, 3), dtype=np.float64)  # z,y,x
    for axis in range(3):
        ta, tb = (axis + 1) % 3, (axis + 2) % 3
        inverse_spacing2 = 1.0 / spacing[axis] ** 2
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                if a == 0:
                    weight = coefficient['alpha'][axis] if b == 0 else coefficient['beta'][tb]
                else:
                    weight = coefficient['beta'][ta] if b == 0 else coefficient['delta_by_derivative_axis'][axis]
                point = [1, 1, 1]
                point[ta] += a
                point[tb] += b
                result[tuple(reversed(point))] += 2.0 * weight * inverse_spacing2
                point[axis] -= 1
                result[tuple(reversed(point))] -= weight * inverse_spacing2
                point[axis] += 2
                result[tuple(reversed(point))] -= weight * inverse_spacing2
    amount = translation / spacing[2] ** 2
    result[1, 1, 1] -= 2.0 * amount
    result[0, 1, 1] += amount
    result[2, 1, 1] += amount
    return result


def modal_tridiagonal(stencil, nx, ny):
    """Transverse sine modes leave diagonal a and equal z couplings -b."""
    theta_x = np.pi * np.arange(1, nx + 1) / (nx + 1)
    theta_y = np.pi * np.arange(1, ny + 1) / (ny + 1)
    planes = np.zeros((3, ny, nx))
    for z in range(3):
        for y in range(3):
            for x in range(3):
                planes[z] += stencil[z, y, x] * np.cos((y - 1) * theta_y)[:, None] * np.cos((x - 1) * theta_x)[None, :]
    scale = float(np.max(np.abs(planes)))
    require(np.max(np.abs(planes[0] - planes[2])) <= SYMBOL_TOLERANCE * scale,
            'Lower/upper modal couplings are not symmetric')
    return planes[1], -0.5 * (planes[0] + planes[2])


def analytic_spectrum(spacing, coefficient, translation, nx, ny, nz):
    """Independent symbol formula used by the existing 3-D DST oracle."""
    angles = [np.pi * np.arange(1, n + 1) / (n + 1) for n in (nx, ny, nz)]
    shapes = [(1, 1, nx), (1, ny, 1), (nz, 1, 1)]
    cosines = [np.cos(t).reshape(shape) for t, shape in zip(angles, shapes)]
    waves2 = [(2.0 * np.sin(0.5 * t) / h) ** 2 for t, h in zip(angles, spacing)]
    result = np.zeros((nz, ny, nx))
    for axis in range(3):
        p, q = (axis + 1) % 3, (axis + 2) % 3
        symbol = (coefficient['alpha'][axis] + 2.0 * coefficient['beta'][p] * cosines[p]
                  + 2.0 * coefficient['beta'][q] * cosines[q]
                  + 4.0 * coefficient['delta_by_derivative_axis'][axis] * cosines[p] * cosines[q])
        if axis == 2:
            symbol = symbol - translation
        result += symbol * waves2[axis].reshape(shapes[axis])
    require(np.all(result > 0), 'Nonpositive finite-box spectrum')
    return result


class Source:
    """Reproducible block source; stores transverse arrays only, never Nz planes."""
    def __init__(self, kind, omit_far=False):
        require(kind in ('asymmetric', 'lowest_mode'), 'Unknown source case')
        self.kind, self.omit_far = kind, omit_far
        x, y = np.arange(1, NX + 1)[None, :], np.arange(1, NY + 1)[:, None]
        self.local_xy = np.exp(-0.5 * ((x - 3.4) / 1.4) ** 2 - 0.5 * ((y - 10.3) / 1.7) ** 2)
        self.remote_xy = np.exp(-0.5 * ((x - 9.2) / 1.2) ** 2 - 0.5 * ((y - 4.1) / 1.6) ** 2)
        self.mode_xy = np.sin(np.pi * x / (NX + 1)) * np.sin(np.pi * y / (NY + 1))

    def block(self, first, stop):
        require(0 <= first < stop <= NZ, 'Out-of-bounds source request')
        z = np.arange(first + 1, stop + 1)[:, None, None]
        local_z = np.exp(-0.5 * ((z - 16.3) / 2.1) ** 2)
        far_z = np.exp(-0.5 * ((z - 67.8) / 2.6) ** 2)
        if self.kind == 'lowest_mode':
            return (local_z + (0.0 if self.omit_far else 0.8) * far_z) * self.mode_xy
        rhs = local_z * self.local_xy
        if not self.omit_far:
            rhs = rhs + 0.7 * far_z * self.remote_xy
            if first <= 60 < stop:
                rhs[60 - first, 10, 4] += 0.35
        return rhs


def dst2(block):
    return dstn(block, type=1, axes=(1, 2), norm='ortho', workers=1)


class BlockwiseModalSolver:
    """Exact block-replayed Thomas elimination for each transverse mode.

    Retains 2*A*ceil(Nz/B) boundary-summary scalars in RAM, plus O(A*B)
    reconstruction state. Source is read twice. Future disk-backed summaries
    could remove their longitudinal RAM growth but are not implemented.
    """
    def __init__(self, diagonal, coupling_b, block_size):
        require(1 <= block_size <= 33, 'Outside bounded block size')
        require(diagonal.shape == (NY, NX), 'Unexpected modal shape')
        self.diagonal, self.off = diagonal, -coupling_b
        self.block_size = block_size
        self.checkpoints = []
        self.source_planes_read = 0
        self.minimum_pivot = math.inf

    def solve(self, source, consume):
        previous_c = np.zeros((NY, NX))
        previous_d = np.zeros((NY, NX))
        for first in range(0, NZ, self.block_size):
            stop = min(NZ, first + self.block_size)
            self.checkpoints.append((first, stop, previous_c.copy(), previous_d.copy()))
            rhs = dst2(source.block(first, stop))
            self.source_planes_read += stop - first
            for plane in rhs:
                pivot = self.diagonal - self.off * previous_c
                require(np.all(pivot > 0), 'Thomas elimination lost a positive pivot')
                self.minimum_pivot = min(self.minimum_pivot, float(np.min(pivot)))
                previous_c = self.off / pivot
                previous_d = (plane - self.off * previous_d) / pivot

        next_phi = np.zeros((NY, NX))  # Upper zero Dirichlet endpoint.
        for first, stop, start_c, start_d in reversed(self.checkpoints):
            rhs = dst2(source.block(first, stop))
            self.source_planes_read += stop - first
            cprime, dprime = np.empty_like(rhs), np.empty_like(rhs)
            previous_c, previous_d = start_c, start_d
            for index, plane in enumerate(rhs):
                pivot = self.diagonal - self.off * previous_c
                require(np.all(pivot > 0), 'Block replay lost a positive pivot')
                previous_c = self.off / pivot
                previous_d = (plane - self.off * previous_d) / pivot
                cprime[index], dprime[index] = previous_c, previous_d
            for index in range(stop - first - 1, -1, -1):
                dprime[index] -= cprime[index] * next_phi
                next_phi = dprime[index].copy()
            consume(first, stop, dst2(dprime))
        require(self.source_planes_read == 2 * NZ, 'Unexpected source pass count')
        area = NX * NY
        return {
            'block_size': self.block_size, 'blocks': len(self.checkpoints),
            'checkpoint_numpy_bytes': sum(a.nbytes + b.nbytes for _, _, a, b in self.checkpoints),
            'four_largest_block_arrays_bytes': 4 * self.block_size * area * 8,
            'source_planes_read': self.source_planes_read,
            'transverse_transform_plane_count': 3 * NZ,
            'minimum_elimination_pivot': self.minimum_pivot,
            'candidate_retains_complete_phi_or_rhs': False,
            'memory_scope': 'Recorded numeric arrays, not measured RSS/peak; excludes FFT workspace, NumPy expression temporaries, source transverse profiles and Python metadata. Checkpoints currently grow as O(A*Nz/B) in RAM; block workspace O(A*B). No constant-memory production claim.',
        }


class StreamingAudit:
    """Small-test oracle is dense; candidate residual needs only two z planes."""
    def __init__(self, reference, source, stencil):
        self.reference, self.source, self.stencil = reference, source, stencil
        self.next_stop = NZ
        self.pending_z = None
        self.pending = None
        self.upper = np.zeros((NY, NX))
        self.error2 = np.longdouble(0)
        self.reference2 = np.longdouble(0)
        self.residual2 = np.longdouble(0)
        self.source2 = np.longdouble(0)
        self.near_error2 = np.longdouble(0)
        self.near_reference2 = np.longdouble(0)
        self.maximum_difference = 0.0
        self.rows = 0

    def residual_plane(self, lower, center, upper, z):
        result = np.zeros_like(center)
        for iz, plane in enumerate((lower, center, upper)):
            padded = np.pad(plane, ((1, 1), (1, 1)))
            for iy in range(3):
                for ix in range(3):
                    result += self.stencil[iz, iy, ix] * padded[iy:iy + NY, ix:ix + NX]
        rhs = self.source.block(z, z + 1)[0]
        self.residual2 += squared_norm(result - rhs)
        self.source2 += squared_norm(rhs)

    def consume(self, first, stop, values):
        require(stop == self.next_stop and values.shape == (stop - first, NY, NX),
                'Reverse reconstruction omitted, duplicated or reordered a block')
        require(np.all(np.isfinite(values)), 'Nonfinite reconstructed field')
        difference = values - self.reference[first:stop]
        self.error2 += squared_norm(difference)
        self.reference2 += squared_norm(self.reference[first:stop])
        self.maximum_difference = max(self.maximum_difference, float(np.max(np.abs(difference))))
        for index in range(stop - first - 1, -1, -1):
            z, plane = first + index, values[index]
            if 8 <= z < 25:
                self.near_error2 += squared_norm(plane - self.reference[z])
                self.near_reference2 += squared_norm(self.reference[z])
            if self.pending_z is not None:
                require(self.pending_z == z + 1, 'Streaming residual lost an adjacent plane')
                self.residual_plane(plane, self.pending, self.upper, self.pending_z)
                self.upper = self.pending
            self.pending = plane.copy()
            self.pending_z = z
            self.rows += 1
        self.next_stop = first

    def finish(self):
        require(self.next_stop == 0 and self.pending_z == 0 and self.rows == NZ,
                'Incomplete reconstructed field')
        self.residual_plane(np.zeros((NY, NX)), self.pending, self.upper, 0)
        return {
            'phi_relative_L2': float(np.sqrt(self.error2 / self.reference2)),
            'maximum_abs_phi_difference': self.maximum_difference,
            'spatial_27point_relative_residual': float(np.sqrt(self.residual2 / self.source2)),
            'near_region_z_indices_half_open': [8, 25],
            'near_region_phi_relative_L2': float(np.sqrt(self.near_error2 / self.near_reference2)),
            'reconstructed_planes': self.rows,
        }


def production_modes(stencil, metrics, dz):
    guards = metrics['cpml_cells']
    nx = int(metrics['nx']) - 2 * int(guards[0]) - 1
    ny = int(metrics['ny']) - 2 * int(guards[1]) - 1
    require((nx, ny) == (139, 139), 'Unexpected production transverse physical domain')
    diagonal, b = modal_tridiagonal(stencil, nx, ny)
    require(np.all(diagonal > 2.0 * np.abs(b)), 'No exponentially decaying infinite-z branch for some mode')
    require(np.all(b > 0), 'This bounded report expects positive production modal coupling')
    exponent = np.arccosh(diagonal / (2.0 * b))
    lengths = 1.0 / exponent
    slow = np.unravel_index(int(np.argmax(lengths)), lengths.shape)
    slow_length = float(lengths[slow])
    return {
        'physical_cell_counts_xy': [nx + 1, ny + 1], 'interior_mode_counts_xy': [nx, ny],
        'array_order': 'row transverse y mode 1..139, column x mode 1..139',
        'equation': '-b*phi[k-1] + a*phi[k] - b*phi[k+1] = rho_mode[k]/epsilon0',
        'evaluated_parameter_regime': 'All modes for these production spacings and common beta have b>0 and a>2*b. The reported positive-root decay formula is restricted to this checked regime, not asserted for arbitrary Cowan parameters.',
        'a_diagonal_per_m2': diagonal.tolist(), 'b_coupling_per_m2': b.tolist(),
        'b_over_a': (b / diagonal).tolist(), 'e_fold_length_cells': lengths.tolist(),
        'b_over_a_minimum': float(np.min(b / diagonal)),
        'b_over_a_maximum': float(np.max(b / diagonal)),
        'slowest_mode_xy': [int(slow[1] + 1), int(slow[0] + 1)],
        'slowest_e_fold_cells': slow_length, 'slowest_e_fold_m': slow_length * dz,
        'slowest_mode_amplitude_at_distance_cells': {
            str(distance): math.exp(-distance / slow_length) for distance in (64, 128, 256, 512)},
        'scope': 'Homogeneous infinite-z decay of the same finite transverse zero-Dirichlet Cowan elliptic operator with a common beta. Not an error bound for omitted beam tails, a moving Maxwell boundary, CPML history, free-space radiation, or actual source modal amplitudes.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--metrics', type=Path, default=DEFAULT_METRICS)
    args = parser.parse_args()
    require(not args.output.exists(), 'Refusing to overwrite output: ' + str(args.output))
    metrics = json.loads(args.metrics.read_text())
    spacing = [float(metrics[k]) for k in ('dx_m', 'dy_m', 'dz_m')]
    beta = float(metrics['particle_beta_mean'])
    require(0 < beta < 1, 'Unexpected common source beta')
    constants_path = REPO / 'src/eb_field.cpp'
    constants = read_constants(constants_path)
    translation = (beta * constants['c']) ** 2 * constants['epsilon0'] * constants['mu0']
    coefficient = cowan_coefficients(spacing)
    stencil = stencil27(spacing, coefficient, translation)
    diagonal, b = modal_tridiagonal(stencil, NX, NY)
    eigenvalues = analytic_spectrum(spacing, coefficient, translation, NX, NY, NZ)
    theta_z = np.pi * np.arange(1, NZ + 1) / (NZ + 1)
    modal_eigenvalues = diagonal[None] - 2.0 * b[None] * np.cos(theta_z)[:, None, None]
    symbol_difference = float(np.max(np.abs(modal_eigenvalues - eigenvalues)) / np.max(np.abs(eigenvalues)))
    require(symbol_difference <= SYMBOL_TOLERANCE, 'Independent spectral and stencil symbols differ')
    require(np.all(diagonal > 2.0 * np.abs(b)), 'Unstable modal elimination fixture')

    cases, negative_controls = [], []
    for name in ('asymmetric', 'lowest_mode'):
        source = Source(name)
        rhs = source.block(0, NZ)  # Dense SMALL-BOX oracle only.
        reference = dstn(dstn(rhs, type=1, norm='ortho', workers=1) / eigenvalues,
                         type=1, norm='ortho', workers=1)
        for block_size in BLOCK_SIZES:
            audit = StreamingAudit(reference, source, stencil)
            memory = BlockwiseModalSolver(diagonal, b, block_size).solve(source, audit.consume)
            result = audit.finish()
            require(result['phi_relative_L2'] < FIELD_TOLERANCE, 'Blocked solve differs from 3-D DST')
            require(result['spatial_27point_relative_residual'] < RESIDUAL_TOLERANCE,
                    'Blocked solve fails independent spatial 27-point residual')
            cases.append({'source': name, 'result': result, 'candidate_state': memory})
        truncated = Source(name, omit_far=True)
        audit = StreamingAudit(reference, truncated, stencil)
        memory = BlockwiseModalSolver(diagonal, b, 16).solve(truncated, audit.consume)
        result = audit.finish()
        require(result['spatial_27point_relative_residual'] < RESIDUAL_TOLERANCE,
                'Negative control must still solve its own truncated-source equation')
        require(result['near_region_phi_relative_L2'] > NEGATIVE_CONTROL_MINIMUM,
                'Missing-far-source negative control did not affect the near region')
        negative_controls.append({'source': name, 'omitted': 'remote Gaussian and asymmetric remote point when present',
                                  'expected_difference_detected': True, 'result': result, 'candidate_state': memory})

    report = {
        'format': 'aprl-modal-reference-algebra-v1', 'status': 'PASS',
        'classification': 'Bounded algebra verification only; not PIC or a measured production optimization',
        'small_box_interior_nodes_xyz': [NX, NY, NZ],
        'small_box_cell_counts_xyz': [NX + 1, NY + 1, NZ + 1],
        'boundaries': 'Homogeneous zero potential on all six finite-box faces',
        'production_spacing_m': spacing, 'common_source_beta': beta,
        'source_beta_is_not_window_beta': True,
        'original_double_SI_constants': constants, 'translation_coefficient': translation,
        'cowan_coefficients': coefficient, 'stencil27_zyx_per_m2': stencil.tolist(),
        'independent_symbol_max_difference_over_max_eigenvalue': symbol_difference,
        'predeclared_tolerances': {'phi_relative_L2': FIELD_TOLERANCE,
                                 'spatial_residual_relative_L2': RESIDUAL_TOLERANCE,
                                 'independent_symbol': SYMBOL_TOLERANCE,
                                 'missing_far_source_near_difference_minimum': NEGATIVE_CONTROL_MINIMUM},
        'cases': cases, 'negative_controls': negative_controls,
        'production_transverse_modes': production_modes(stencil, metrics, spacing[2]),
        'input_sha256': {str(args.metrics.resolve()): digest(args.metrics),
                         str(constants_path): digest(constants_path),
                         str(REPO / 'src/eb_initial_field.cpp'): digest(REPO / 'src/eb_initial_field.cpp'),
                         str(Path(__file__).resolve()): digest(__file__)},
        'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__},
        'thread_environment': {name: os.environ[name] for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        'limitations': [
            'Dense 3-D oracle and its source exist only in the fixed small algebra fixture.',
            'Candidate summaries are currently in RAM and grow as A*ceil(Nz/block); disk-backed summary I/O is not implemented.',
            'Candidate reads each source plane twice and performs three transverse transforms per plane; source generation, sorting, disk I/O and MPI costs are not measured.',
            'A changing source or velocity distribution requires additional work; this script solves one fixed common-beta source state only.',
            'No E/B reconstruction, leapfrog phase update, moving inlet coupling, particle feedback or CPML history is validated.',
            'Near-only source truncation is deliberately a negative control; far-source contributions cannot be discarded merely because only a local field is requested.',
            'No production speedup, constant total memory, or full streaming complexity claim follows from this check.',
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'status': report['status'], 'classification': report['classification'],
                      'cases': len(cases), 'negative_controls': len(negative_controls),
                      'slowest_production_e_fold_cells': report['production_transverse_modes']['slowest_e_fold_cells'],
                      'output': str(args.output.resolve())}, sort_keys=True))


if __name__ == '__main__':
    main()
