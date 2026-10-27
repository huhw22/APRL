#!/usr/bin/env python3
"""Plot completed self-field startup probes as a standalone scientific figure.

Usage: plot_startup.py RESULT_ROOT [--require-complete]
RESULT_ROOT contains r1/<case>/{history.csv,metadata.json,EXIT_CODE}.
The figure is descriptive; it does not turn numerical completion into physics
acceptance or compare implementation throughput.
"""
import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aprl-startup-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


STYLES = {
    ("zero", "base"): dict(label="Zero field (negative control)", color="#4B5563",
                            linestyle="--", marker="x"),
    ("qs", "base"): dict(label="Quasi-steady: base", color="#1565A8",
                          linestyle="-", marker="o"),
    ("qs", "wide"): dict(label="Quasi-steady: wider domain", color="#C47418",
                          linestyle="-.", marker="s"),
    ("qs", "fine"): dict(label="Quasi-steady: finer mesh", color="#7955A1",
                          linestyle=":", marker="^"),
}
REQUIRED = (
    "step", "time_s", "var_x_m2", "delta_K_J", "UE_J", "UB_centered_J"
)


def load(root):
    cases = {}
    skipped = []
    for path in sorted((root / "r1").glob("*/metadata.json")):
        meta = json.loads(path.read_text())
        kind = (meta["mode"], meta["grid"])
        if kind not in STYLES:
            continue
        beta = float(meta["beta"])
        if beta not in (0.0, 0.23):
            raise ValueError(f"Unexpected beta in {path}: {beta}")
        exit_path = path.parent / "EXIT_CODE"
        if not meta.get("completed") or not exit_path.exists() or exit_path.read_text().strip() != "0":
            skipped.append(str(path.parent))
            continue
        rows = np.atleast_1d(np.genfromtxt(path.parent / "history.csv",
                                          delimiter=",", names=True))
        if not set(REQUIRED).issubset(rows.dtype.names or ()):
            raise ValueError(f"Missing required history columns: {path.parent}")
        for name in REQUIRED:
            if not np.isfinite(rows[name]).all():
                raise ValueError(f"Non-finite {name}: {path.parent}")
        steps = int(meta["steps"])
        if len(rows) != steps + 1 or not np.array_equal(rows["step"], np.arange(steps + 1)):
            raise ValueError(f"Incomplete or discontinuous history: {path.parent}")
        if not np.allclose(rows["time_s"], np.arange(steps + 1) * meta["dt_s"],
                           rtol=4e-15, atol=0):
            raise ValueError(f"Unexpected physical time grid: {path.parent}")
        if np.any(rows["var_x_m2"] <= 0):
            raise ValueError(f"Nonpositive transverse variance: {path.parent}")
        if int(meta["particles"]) != 4096 or int(meta["seed"]) != 710271:
            raise ValueError(f"Unexpected particle population: {path.parent}")
        key = (beta,) + kind
        if key in cases:
            raise ValueError(f"Duplicate completed case: {key}")
        cases[key] = (rows, meta, path.parent)
    return cases, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--require-complete", action="store_true",
                        help="Require all eight plotted beta/mode/grid cases")
    args = parser.parse_args()
    root = args.root.resolve()
    cases, skipped = load(root)
    expected = [(beta,) + kind for beta in (0.0, 0.23) for kind in STYLES]
    missing = [key for key in expected if key not in cases]
    if not cases:
        raise ValueError("No completed startup cases found")
    if args.require_complete and missing:
        raise ValueError(f"Missing complete plotted cases: {missing}")
    end_times = [rows["time_s"][-1] for rows, _, _ in cases.values()]
    if not np.allclose(end_times, end_times[0], rtol=4e-15, atol=0):
        raise ValueError("Cases end at different physical times")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "legend.fontsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0), sharex=True)
    fig.subplots_adjust(left=.076, right=.984, bottom=.19, top=.79,
                        hspace=.36, wspace=.29)
    for row, beta in enumerate((0.0, 0.23)):
        for kind, style in STYLES.items():
            key = (beta,) + kind
            if key not in cases:
                continue
            data, _, _ = cases[key]
            time_fs = data["time_s"] * 1e15
            growth = 100 * (np.sqrt(data["var_x_m2"] / data["var_x_m2"][0]) - 1)
            kinetic = data["delta_K_J"] * 1e12
            field_energy = data["UE_J"] + data["UB_centered_J"]
            total_change = (data["delta_K_J"] + field_energy - field_energy[0]) * 1e12
            for col, values in enumerate((growth, kinetic, total_change)):
                axes[row, col].plot(time_fs, values, **style, linewidth=1.7,
                                    markersize=4.0, markerfacecolor="white",
                                    markevery=max(1, (len(time_fs)-1)//5))
        for ax in axes[row]:
            ax.axhline(0, color="#B8BEC7", linewidth=.7, zorder=0)
            ax.grid(axis="y", color="#DCE0E5", linewidth=.6, alpha=.8)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
            ax.set_xlim(0, max(end_times) * 1e15)
            ax.set_xticks([0, 10, 20, 30, 40, 50])
        axes[row, 0].text(.025, .96, rf"$\beta = {beta:g}$",
                          transform=axes[row, 0].transAxes, va="top", fontsize=11,
                          bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=2))
        axes[row, 0].set_ylabel(r"$100[\sigma_x(t)/\sigma_x(0)-1]$ (%)")
        axes[row, 1].set_ylabel(r"$\Delta K$ (pJ)")
        axes[row, 2].set_ylabel(r"$\Delta K+\Delta U_{\rm field}$ (pJ)")
    for ax, title in zip(axes[0], ("Transverse RMS size", "Particle energy change",
                                   "Particle + field energy change")):
        ax.set_title(title, pad=10)
    for ax in axes[1]:
        ax.set_xlabel("Time after release (fs)")

    # Preserve the full zero-field energy deficit on the main linear axis,
    # while showing the small positive quasi-steady changes at readable scale.
    moving_qs = [(kind, cases[(.23,) + kind]) for kind in STYLES
                 if kind[0] == "qs" and (.23,) + kind in cases]
    if moving_qs:
        inset = axes[1, 1].inset_axes([.52, .38, .44, .47])
        maximum = 0.
        for kind, (data, _, _) in moving_qs:
            values = data["delta_K_J"] * 1e12
            inset.plot(data["time_s"] * 1e15, values, **STYLES[kind],
                       linewidth=1.4, markersize=3, markerfacecolor="white",
                       markevery=max(1, (len(values)-1)//4))
            maximum = max(maximum, float(np.max(values)))
        inset.set_xlim(0, max(end_times) * 1e15)
        inset.set_ylim(0, max(1e-12, maximum) * 1.12)
        inset.set_xticks([0, 25, 50])
        inset.tick_params(labelsize=7, length=2)
        inset.set_title("QS detail (pJ)", fontsize=8, pad=3)
        inset.grid(axis="y", alpha=.2)

    legend = [Line2D([], [], **style, linewidth=1.7, markersize=4,
                     markerfacecolor="white") for style in STYLES.values()]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.53, .89),
               ncol=4, frameon=False, columnspacing=1.8)
    fig.suptitle("Self-field startup: fixed 4096-particle Gaussian beam",
                 fontsize=15, y=.975)
    fig.text(.5, .932,
             "Identical release positions and momenta in each comparison; no held interval or moving window",
             ha="center", fontsize=10, color="#374151")
    fig.text(.076, .118,
             r"$\Delta U_{\rm field}$ uses Yee electric energy and time-centered magnetic energy. "
             "PEC boundaries; wider-domain and finer-mesh comparisons are separate.",
             ha="left", fontsize=9, color="#374151")
    fig.text(.076, .084,
             "Zero initial field violates Gauss's law. The ordinary Cowan field-energy sum is not an exact discrete invariant. "
             "Y-axis scales vary by row.",
             ha="left", fontsize=9, color="#374151")
    last_note = "Source: r1/*/history.csv and metadata.json. Fixed particle count; this is not a particle-number convergence study."
    if missing:
        last_note = f"PARTIAL FIGURE: {len(cases)}/8 plotted cases complete; missing cases are omitted, not zero-filled."
    fig.text(.076, .050, last_note, ha="left", fontsize=9, color="#374151")
    for suffix in ("png", "pdf"):
        fig.savefig(root / f"startup.{suffix}", dpi=200)
    plt.close(fig)
    print(json.dumps({"figure": str(root / "startup.png"),
                      "completed_cases": [str(value[2]) for value in cases.values()],
                      "missing_cases": missing, "skipped_incomplete": skipped}, indent=2))


if __name__ == "__main__":
    main()
