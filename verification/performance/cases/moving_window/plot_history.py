#!/usr/bin/env python3
"""Plots of measured moving-CPML history-policy differences and costs."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


COLORS = {'world-position': '#1565c0', 'boundary-relative-z': '#d55e00'}
LABELS = {'world-position': 'History at physical position',
          'boundary-relative-z': 'History at CPML depth'}


def save(fig, directory, stem):
    fig.savefig(directory / (stem + '.png'))
    fig.savefig(directory / (stem + '.pdf'))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review', type=Path)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--long-review', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    review = json.loads(args.review.read_text())
    root = args.root or Path(review['root'])
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'axes.axisbelow': True,
                         'grid.alpha': .23, 'savefig.dpi': 180})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), layout='constrained')
    for row, (depth, world, boundary) in enumerate([(20, 'wave20', 'wave-boundary20'),
                                                   (8, 'wave', 'wave-boundary')]):
        data = {policy: np.genfromtxt(root / case / 'wave-plane.csv', delimiter=',', names=True)
                for policy, case in [('world-position', world), ('boundary-relative-z', boundary)]}
        reference = data['world-position']
        norm = np.max(reference['combined_reference_rms_V_m'])
        axes[row, 0].plot(reference['step'], reference['Ey_reference_center_V_m'],
                          color='#333333', linewidth=2.5, label='Common large domain')
        for policy, d in data.items():
            axes[row, 0].plot(d['step'], d['Ey_candidate_center_V_m'], color=COLORS[policy],
                              linestyle='--' if policy == 'world-position' else ':', label=LABELS[policy])
            axes[row, 1].semilogy(d['step'], np.maximum(d['combined_difference_rms_V_m'] / norm, 1e-16),
                                  color=COLORS[policy], label=LABELS[policy])
            late = d['step'] >= 100
            axes[row, 2].plot(d['step'][late], 100 * d['combined_difference_rms_V_m'][late] / norm,
                              color=COLORS[policy], label=LABELS[policy])
        axes[row, 0].set_xlim(0, 99)
        axes[row, 0].set_title(f'z CPML = {depth} cells: main pulse')
        axes[row, 0].set_ylabel('Center native Ey (V/m)')
        axes[row, 1].set_title(f'z CPML = {depth}: full observation')
        axes[row, 1].set_ylabel('Plane error RMS / reference peak')
        axes[row, 1].axvspan(100, 176, color='#999999', alpha=.12)
        axes[row, 1].text(.02, .03, 'Zero values shown at 1e-16', transform=axes[row, 1].transAxes, fontsize=8)
        axes[row, 2].set_title(f'z CPML = {depth}: late response')
        axes[row, 2].set_ylabel('Plane error RMS / reference peak (%)')
        for ax in axes[row]:
            ax.set_xlabel('Step')
            ax.legend(fontsize=8)
    fig.suptitle('Same moving boundary and unchanged Cowan + CPML recurrence; only longitudinal history mapping differs\n'
                 'Finite-interval field differences against a common large domain; not reflection coefficients', fontsize=12)
    save(fig, args.output, 'cpml_history_wave_comparison')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout='constrained')
    widths = [64, 128, 256]
    x = np.arange(len(widths))
    for policy, offset in [('world-position', -.18), ('boundary-relative-z', .18)]:
        costs = [next(c for c in review['costs'] if c['width'] == width and c['policy'] == policy)
                 for width in widths]
        # Each point is the ratio within one repetition of its own case.
        ratios = [[100 * value for value in c['paired_total_overhead_fractions']] for c in costs]
        med = np.array([np.median(values) for values in ratios])
        lo, hi = np.array([min(v) for v in ratios]), np.array([max(v) for v in ratios])
        axes[0].bar(x + offset, med, width=.34, color=COLORS[policy], label=LABELS[policy],
                     yerr=np.array([med - lo, hi - med]), capsize=3)
        step_ms = np.array([c['modes']['compact']['move_s'] / 64 * 1000 for c in costs])
        axes[1].bar(x + offset, step_ms, width=.34, color=COLORS[policy], label=LABELS[policy])
    axes[0].set_ylabel('Extra total time over own fixed-window run (%)')
    axes[0].set_title('Median of 3 paired repeat ratios; bars show range')
    axes[1].set_ylabel('Measured migration time per step (ms)')
    axes[1].set_title('Median of 3 runs; includes new-layer supply')
    for ax in axes:
        ax.set_xticks(x, widths)
        ax.set_xlabel('Window z cells')
        ax.legend(fontsize=8)
    fig.suptitle('1 MPI rank / 1 physical core; same width and original update kernels\n'
                 'Policies ran as separate cases: cross-policy timing differences are descriptive', fontsize=12)
    save(fig, args.output, 'cpml_history_migration_cost')

    if args.long_review:
        extended = json.loads(args.long_review.read_text())
        long_root = Path(extended['root'])
        fig, axes = plt.subplots(2, 2, figsize=(11.8, 8), layout='constrained')
        for row, (depth, world, boundary) in enumerate([
                (20, 'wave-long20', 'wave-boundary-long20'),
                (8, 'wave-long', 'wave-boundary-long')]):
            for policy, case in [('world-position', world), ('boundary-relative-z', boundary)]:
                d = np.genfromtxt(long_root / case / 'wave-plane.csv', delimiter=',', names=True)
                norm = np.max(d['combined_reference_rms_V_m'])
                error = d['combined_difference_rms_V_m'] / norm
                ref_square = np.sum(d['combined_reference_rms_V_m'] ** 2)
                accumulated = np.cumsum(d['combined_difference_rms_V_m'] ** 2) / ref_square
                axes[row, 0].plot(d['step'], error * 100, color=COLORS[policy], label=LABELS[policy])
                axes[row, 1].plot(d['step'], accumulated * 100, color=COLORS[policy], label=LABELS[policy])
                peak = int(np.argmax(error))
                axes[row, 0].plot(d['step'][peak], error[peak] * 100, 'o', color=COLORS[policy], markersize=4)
            for ax in axes[row]:
                ax.set_xlim(100, 384)
                ax.axvline(176, color='#555555', linestyle='--', linewidth=1, label='Original observation end')
                ax.set_xlabel('Step')
                ax.legend(fontsize=8)
            axes[row, 0].set_title(f'z CPML = {depth}: extended tail; dots mark peaks')
            axes[row, 0].set_ylabel('Plane error RMS / reference peak (%)')
            axes[row, 1].set_title(f'z CPML = {depth}: accumulated squared amplitude')
            axes[row, 1].set_ylabel('Cumulative squared difference / reference (%)')
        fig.suptitle('Same candidate windows observed for 384 steps; common reference enlarged to 1024 cells\n'
                     'Separate delay from amplitude reduction; squared-amplitude sums are not energy or reflection coefficients', fontsize=12)
        save(fig, args.output, 'cpml_history_extended_tail')

        fig, ax = plt.subplots(figsize=(8, 4.2), layout='constrained')
        for policy, case in [('world-position', 'wave-long20'),
                             ('boundary-relative-z', 'wave-boundary-long20')]:
            d = np.genfromtxt(long_root / case / 'wave-plane.csv', delimiter=',', names=True)
            error = 100 * d['combined_difference_rms_V_m'] / np.max(d['combined_reference_rms_V_m'])
            ax.plot(d['step'], error, color=COLORS[policy], label=LABELS[policy])
            peak = int(np.argmax(error))
            ax.plot(d['step'][peak], error[peak], 'o', color=COLORS[policy])
            ax.annotate(f'{error[peak]:.2f}%', (d['step'][peak], error[peak]),
                        xytext=(8, 5), textcoords='offset points', color=COLORS[policy])
        ax.axvline(176, color='#777777', linestyle='--', label='Original observation end')
        ax.set(xlim=(100, 384), xlabel='Step', ylabel='Plane difference RMS / main pulse peak (%)',
               title='20-cell longitudinal CPML: delayed response is larger')
        ax.legend(fontsize=9)
        save(fig, args.output, 'cpml_history_20_tail')


if __name__ == '__main__':
    main()
