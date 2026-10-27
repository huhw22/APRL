#!/usr/bin/env python3
"""Plot measured probe errors and counted work, never production speedup."""
import argparse
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-window-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('campaign', type=Path)
args = p.parse_args()
root = args.campaign.resolve()
review = json.loads((root/'SCIENCE_REVIEW.json').read_text())
oracle = json.loads((root/'probe-r1/window/oracle-results.json').read_text())
assert review['pass_']
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
cases = [c for c in oracle['cases'] if c['beta'] > .99]
t = np.array([c['steps'] for c in cases])
reference = np.array([c['reference_cell_updates_per_grid'] for c in cases])
window = np.array([c['window_cell_updates_per_grid'] for c in cases])
axes[0].plot(t, reference/reference[0], 'o-', color='#64748b', label='Full reference: T(T + 80)')
axes[0].plot(t, window/window[0], 's-', color='#087f8c', label='Window interior: 48T')
axes[0].set(xlabel='Number of field steps', ylabel='Cell updates / value at 48 steps',
            title='Counted work per grid')
axes[0].set_xticks(t)
axes[0].legend(loc='upper left', frameon=False, fontsize=9)
axes[0].set_title('Counted work per grid\nWindow-only work excludes reference cost', fontsize=11)
for high, color, label in [(False, '#64748b', 'Window speed = 0.9c'),
                           (True, '#087f8c', 'Window speed = production boost')]:
    rows = [c for c in review['analytic'] if (c['beta'] > .99) == high]
    axes[1].semilogy([c['steps'] for c in rows], [c['weak_error_over_weak'] for c in rows],
                     'o-', color=color, label=label)
axes[1].axhline(1e-5, color='#c65d21', linestyle='--', label='Declared limit')
axes[1].set(xlabel='Number of field steps', ylabel='Max error / weak-field peak',
            title='Independent analytic vacuum boundary', ylim=(1e-9, 1e-4))
axes[1].set_xticks(t)
axes[1].legend(loc='upper right', frameon=False, fontsize=8)
for ax in axes:
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=.18)
fig.suptitle('Fixed-window feasibility: conditional replay and known-wave tests', fontsize=13)
fig.savefig(root/'feasibility.png', dpi=180)
fig.savefig(root/'feasibility.pdf')
plt.close(fig)
