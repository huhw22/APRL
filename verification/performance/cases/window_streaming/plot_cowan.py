#!/usr/bin/env python3
"""Standalone Cowan-reference / synchronized-held comparison figures.

Usage: plot_cowan.py RESULT_ROOT [--require-complete] [--require-all-grids]
Only completed r1/<case> history.csv + metadata.json are plotted. The first
option requires the eight base cases; the second requires all24 primary cases.
Missing cases are explicitly reported, never filled with zeros. This script
does not grant scientific acceptance or publish to an external service.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aprl-cowan-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


STYLES = {
    "qs": dict(label="Ordinary Poisson", color="#1565A8", linestyle="--", marker="o"),
    "cowan": dict(label="Cowan reference", color="#007A67", linestyle="-", marker="^"),
    "held50": dict(label="Held: 50 fs history", color="#C47418", linestyle="-.", marker="s"),
    "held100": dict(label="Held: 100 fs history", color="#454B54", linestyle=":", marker="D"),
}
BETAS = (0., .23)
GRIDS = ("base", "wide", "fine")
REQUIRED = ("step", "time_s", "var_x_m2", "delta_K_J", "UE_J", "UB_centered_J")


def load(root):
    cases, skipped = {}, []
    for path in sorted((root/"r1").glob("*/metadata.json")):
        meta = json.loads(path.read_text())
        if meta["mode"] not in STYLES or meta["grid"] not in GRIDS:
            continue
        beta = float(meta["beta"])
        if beta not in BETAS:
            raise ValueError(f"Unexpected beta in {path}: {beta}")
        exit_path = path.parent/"EXIT_CODE"
        if not meta.get("completed") or not exit_path.exists() or exit_path.read_text().strip() != "0":
            skipped.append(str(path.parent))
            continue
        data = np.atleast_1d(np.genfromtxt(path.parent/"history.csv", delimiter=",", names=True))
        if not set(REQUIRED).issubset(data.dtype.names or ()):
            raise ValueError(f"Missing plot columns: {path.parent}")
        for name in REQUIRED:
            if not np.isfinite(data[name]).all():
                raise ValueError(f"Nonfinite {name}: {path.parent}")
        steps = int(meta["steps"])
        if len(data) != steps+1 or not np.array_equal(data["step"], np.arange(steps+1)):
            raise ValueError(f"Discontinuous/incomplete history: {path.parent}")
        if not np.allclose(data["time_s"], np.arange(steps+1)*meta["dt_s"], rtol=4e-15, atol=0):
            raise ValueError(f"Unexpected time grid: {path.parent}")
        if np.any(data["var_x_m2"] <= 0):
            raise ValueError(f"Nonpositive transverse variance: {path.parent}")
        if int(meta["particles"]) != 4096 or int(meta["seed"]) != 710271:
            raise ValueError(f"Unexpected particle population: {path.parent}")
        key = beta, meta["grid"], meta["mode"]
        if key in cases:
            raise ValueError(f"Duplicate completed case: {key}")
        state_hash = hashlib.sha256((path.parent/"initial-particles.csv").read_bytes()).hexdigest()
        cases[key] = dict(data=data, meta=meta, path=path.parent, state_hash=state_hash)
    for beta in BETAS:
        hashes = {case["state_hash"] for key, case in cases.items() if key[0] == beta}
        if len(hashes) > 1:
            raise ValueError(f"Release positions/momenta differ for beta={beta}")
    if cases:
        ends = [case["data"]["time_s"][-1] for case in cases.values()]
        if not np.allclose(ends, ends[0], rtol=4e-15, atol=0):
            raise ValueError("Cases end at unequal physical times")
    return cases, skipped


def configure():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "savefig.facecolor": "white"})


def save(fig, root, stem):
    for extension in ("png", "pdf"):
        fig.savefig(root/f"{stem}.{extension}", dpi=200)
    plt.close(fig)


def primary(root, cases):
    expected = [(beta, "base", mode) for beta in BETAS for mode in STYLES]
    missing = [key for key in expected if key not in cases]
    present = [key for key in expected if key in cases]
    if not present:
        return None
    end_fs = max(cases[key]["data"]["time_s"][-1] for key in present)*1e15
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.), sharex=True)
    fig.subplots_adjust(left=.076, right=.984, bottom=.19, top=.79, hspace=.36, wspace=.30)
    for row, beta in enumerate(BETAS):
        for mode, style in STYLES.items():
            key = beta, "base", mode
            if key not in cases:
                continue
            data = cases[key]["data"]
            time = data["time_s"]*1e15
            growth = 100*(np.sqrt(data["var_x_m2"]/data["var_x_m2"][0])-1)
            kinetic = data["delta_K_J"]*1e12
            field = data["UE_J"]+data["UB_centered_J"]
            balance = (data["delta_K_J"]+field-field[0])*1e12
            for col, values in enumerate((growth, kinetic, balance)):
                axes[row,col].plot(time, values, **style, linewidth=1.7,
                    markersize=4, markerfacecolor="white", markevery=max(1,(len(time)-1)//5))
        for ax in axes[row]:
            ax.axhline(0, color="#B8BEC7", linewidth=.7, zorder=0)
            ax.grid(axis="y", color="#DCE0E5", linewidth=.6, alpha=.8)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
            ax.set_xlim(0,end_fs)
            ax.set_xticks([0,10,20,30,40,50])
        axes[row,0].text(.025,.96,rf"$\beta={beta:g}$", transform=axes[row,0].transAxes,
            va="top", fontsize=11, bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=2))
        axes[row,0].set_ylabel(r"$100[\sigma_x(t)/\sigma_x(0)-1]$ (%)")
        axes[row,1].set_ylabel(r"$\Delta K$ (pJ)")
        axes[row,2].set_ylabel(r"$\Delta K+\Delta U_{\rm field}$ (pJ)")
    for ax,title in zip(axes[0],("Transverse RMS size", "Particle energy change", "Particle + field energy diagnostic")):
        ax.set_title(title,pad=10)
    for ax in axes[1]:
        ax.set_xlabel("Time after release (fs)")
    handles = [Line2D([],[],**style,linewidth=1.7,markersize=4,markerfacecolor="white") for style in STYLES.values()]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.53,.89),ncol=4,frameon=False,columnspacing=1.9)
    fig.suptitle("Cowan reference compared with prescribed held histories",fontsize=15,y=.975)
    fig.text(.5,.932,"Same release positions and momenta; 4096-particle Gaussian, 1 pC; base mesh",ha="center",fontsize=10,color="#374151")
    fig.text(.076,.118,"Held is a prescribed straight-trajectory history followed by synchronized release, not production lab-plane release.",fontsize=9,color="#374151")
    fig.text(.076,.084,"Ordinary field energy uses Yee E and time-centered B; it is not an exact Cowan invariant. Each branch starts its energy account at release.",fontsize=9,color="#374151")
    note = "Source: r1/*/{history.csv, metadata.json}. PEC boundaries; fixed particle count. Y-axis scales vary by row; no production acceptance claim."
    if missing:
        note = f"PARTIAL FIGURE: {len(present)}/8 base cases complete; absent cases are omitted, not zero-filled."
    fig.text(.076,.050,note,fontsize=9,color="#374151")
    save(fig,root,"cowan-reference")
    return dict(figure=str(root/"cowan-reference.png"),missing=missing,cases=[str(cases[k]["path"]) for k in present])


def comparison(root, cases):
    complete = [(beta,grid) for beta in BETAS for grid in GRIDS
                if all((beta,grid,mode) in cases for mode in STYLES)]
    if not complete:
        return None
    fig, axes = plt.subplots(2,1,figsize=(10.8,7.8),sharex=True)
    fig.subplots_adjust(left=.105,right=.97,top=.79,bottom=.20,hspace=.32)
    points = []
    for row,beta in enumerate(BETAS):
        ax=axes[row]
        for x,grid in enumerate(GRIDS):
            if (beta,grid) not in complete:
                continue
            ref = cases[(beta,grid,"held100")]["data"]["delta_K_J"][-1]
            held50 = cases[(beta,grid,"held50")]["data"]["delta_K_J"][-1]
            sensitivity=abs(held50-ref)*1e12
            ax.fill_between([x-.35,x+.35],-sensitivity,sensitivity,color="#D5D9DE",alpha=.8,zorder=0)
            for mode,offset in (("qs",-.12),("cowan",.12)):
                data=cases[(beta,grid,mode)]["data"]
                delta=(data["delta_K_J"][-1]-ref)*1e12
                style=STYLES[mode]
                ax.scatter(x+offset,delta,color=style["color"],marker=style["marker"],s=44,edgecolors="white",linewidths=.6,zorder=3)
                points.append(dict(beta=beta,grid=grid,mode=mode,deltaK_difference_pJ=float(delta),held_history_sensitivity_pJ=float(sensitivity)))
        ax.axhline(0,color="#454B54",linewidth=.9)
        ax.grid(axis="y",color="#DCE0E5",linewidth=.6)
        ax.set_xlim(-.6,2.6)
        ax.ticklabel_format(axis="y",style="plain",useOffset=False)
        ax.set_ylabel(r"$\Delta K_{\rm method}-\Delta K_{\rm held100}$ (pJ)")
        ax.text(.025,.94,rf"$\beta={beta:g}$",transform=ax.transAxes,va="top",fontsize=11,
            bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=2))
    axes[1].set_xticks(range(3),["Base: 48³, 1 µm", "Wide: 64³, 1 µm", "Fine: 96³, 0.5 µm"])
    axes[1].set_xlabel("Separate same-grid comparisons; fine mesh also halves the timestep")
    handles=[Line2D([],[],color=STYLES[m]["color"],marker=STYLES[m]["marker"],linestyle="none",label=STYLES[m]["label"],markersize=6) for m in ("qs","cowan")]
    handles.append(Patch(facecolor="#D5D9DE",label="± |held50 − held100| endpoint difference"))
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.52,.884),ncol=3,frameon=False)
    fig.suptitle("Particle-energy difference from the held100 reference",fontsize=15,y=.972)
    fig.text(.5,.924,"End of the same 50 fs release interval; energy changes start at each branch's release state",ha="center",fontsize=10,color="#374151")
    fig.text(.105,.123,"Gray intervals show measured sensitivity to held-history length, not confidence intervals or a physical truth bound.",fontsize=9,color="#374151")
    fig.text(.105,.087,"Synchronized held is a controlled reference. Endpoint differences do not replace the full-history acceptance checks.",fontsize=9,color="#374151")
    note="Source: r1 histories. Fixed 4096-particle samples; independent y-axis scales. No production acceptance is implied."
    if len(complete)<6:
        note=f"PARTIAL COMPARISON: {len(complete)}/6 beta/grid groups complete; incomplete groups are omitted."
    fig.text(.105,.051,note,fontsize=9,color="#374151")
    save(fig,root,"held-comparison")
    return dict(figure=str(root/"held-comparison.png"),complete_groups=complete,points=points)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("--require-complete",action="store_true")
    parser.add_argument("--require-all-grids",action="store_true")
    args=parser.parse_args()
    root=args.root.resolve()
    cases,skipped=load(root)
    expected=[(b,g,m) for b in BETAS for g in (GRIDS if args.require_all_grids else ("base",)) for m in STYLES]
    missing=[k for k in expected if k not in cases]
    if not cases:
        raise ValueError("No completed physical-comparison cases")
    if (args.require_complete or args.require_all_grids) and missing:
        raise ValueError(f"Missing required completed cases: {missing}")
    configure()
    output=dict(primary=primary(root,cases),comparison=comparison(root,cases),skipped_incomplete=skipped,
                input_scope="r1 completed cases only; fixed release-state hashes checked",physical_acceptance="not_assessed_by_plot")
    (root/"plot-cowan-manifest.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    main()
