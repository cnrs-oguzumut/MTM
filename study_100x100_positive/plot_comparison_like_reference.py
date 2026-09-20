#!/usr/bin/env python3
"""
Generate Pre-Yield vs Post-Yield comparison plot matching user's reference figure:
- Two subplots side-by-side
- Panel 1: Flip energy jump |Eafter - Ebefore| in logarithmic bins, Probability per logarithmic bin
- Panel 2: Flip stress jump |\\sigma_{12,after} - \\sigma_{12,before}| in logarithmic bins, Probability per logarithmic bin
- Step histograms for pre-yield and post-yield
- Highlighted Mean, Min, and Max values for both regimes
"""

import csv
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def setup_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 11.0,
        "axes.titlesize": 12.0,
        "axes.labelsize": 11.5,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 9.2,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    })

def find_yield_alpha(sim_dir):
    log_file = sim_dir / "energy_stress_log.csv"
    with open(log_file, "r") as f:
        rows = list(csv.DictReader(f))
    max_row = max(rows, key=lambda r: float(r["PostStress"]))
    return float(max_row["Alpha"]), float(max_row["PostStress"])

def load_nodal_data(vtk_path):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(vtk_path))
    r.Update()
    grid = r.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    ne = vtk_to_numpy(grid.GetPointData().GetArray("NodalEnergy"))
    stress_arr = grid.GetPointData().GetArray("NodalCauchyStress")
    if stress_arr is not None:
        ns = vtk_to_numpy(stress_arr)[:, 1] # sigma_12 shear component
    else:
        ns = np.zeros_like(ne)
    return pts, ne, ns

def collect_data(sim_dir, alpha_yield):
    log_file = sim_dir / "energy_stress_log.csv"
    step_to_alpha = {}
    with open(log_file, "r") as f:
        for r in csv.DictReader(f):
            step_to_alpha[int(r["Iteration"])] = float(r["Alpha"])

    surgery_dir = sim_dir / "vtk_surgery"
    before_files = sorted(surgery_dir.glob("*_1_before_remesh.vtk"))
    print(f"Processing {len(before_files)} surgery passes...")

    pre_dE, post_dE = [], []
    pre_dS, post_dS = [], []
    zeros_E_pre, zeros_E_post = 0, 0
    zeros_S_pre, zeros_S_post = 0, 0

    for bf in before_files:
        af = Path(str(bf).replace("_1_before_remesh.vtk", "_2_after_remesh.vtk"))
        if not af.exists():
            continue

        try:
            step_idx = int(bf.name.split("_")[1])
        except:
            step_idx = 0
        alpha_val = step_to_alpha.get(step_idx, 0.14)

        pts1, ne1, ns1 = load_nodal_data(bf)
        pts2, ne2, ns2 = load_nodal_data(af)

        h2 = { (round(p[0], 3), round(p[1], 3)): (ne2[i], ns2[i]) for i, p in enumerate(pts2) }

        for i, p in enumerate(pts1):
            k = (round(p[0], 3), round(p[1], 3))
            if k in h2:
                e2, s2 = h2[k]
                de = abs(e2 - ne1[i])
                ds = abs(s2 - ns1[i])

                if alpha_val < alpha_yield:
                    if de > 1e-9:
                        pre_dE.append(de)
                    else:
                        zeros_E_pre += 1
                    
                    if ds > 1e-6:
                        pre_dS.append(ds)
                    else:
                        zeros_S_pre += 1
                else:
                    if de > 1e-9:
                        post_dE.append(de)
                    else:
                        zeros_E_post += 1
                    
                    if ds > 1e-6:
                        post_dS.append(ds)
                    else:
                        zeros_S_post += 1

    return {
        "pre_dE": np.array(pre_dE),
        "post_dE": np.array(post_dE),
        "pre_dS": np.array(pre_dS),
        "post_dS": np.array(post_dS),
        "zeros_E_pre": zeros_E_pre,
        "zeros_E_post": zeros_E_post,
        "zeros_S_pre": zeros_S_pre,
        "zeros_S_post": zeros_S_post,
    }

def format_n(n):
    exp = int(np.floor(np.log10(n)))
    coef = n / (10**exp)
    return rf"{coef:.1f} \times 10^{{{exp}}}"

def main():
    setup_style()
    curr_dir = Path(__file__).resolve().parent
    if (curr_dir / "energy_stress_log.csv").exists():
        sim_dir = curr_dir
        repo_root = curr_dir.parent
    else:
        sim_dir = curr_dir / "study_100x100_positive"
        repo_root = curr_dir

    alpha_yield, peak_stress = find_yield_alpha(sim_dir)
    print(f"Alpha yield = {alpha_yield:.5f} (Peak stress = {peak_stress:.5f})")

    data = collect_data(sim_dir, alpha_yield)

    pre_dE = data["pre_dE"]
    post_dE = data["post_dE"]
    pre_dS = data["pre_dS"]
    post_dS = data["post_dS"]

    # Calculate statistics
    stats = {
        "pre_E": {
            "n": len(pre_dE),
            "mean": np.mean(pre_dE),
            "min": np.min(pre_dE),
            "max": np.max(pre_dE),
            "median": np.median(pre_dE)
        },
        "post_E": {
            "n": len(post_dE),
            "mean": np.mean(post_dE),
            "min": np.min(post_dE),
            "max": np.max(post_dE),
            "median": np.median(post_dE)
        },
        "pre_S": {
            "n": len(pre_dS),
            "mean": np.mean(pre_dS),
            "min": np.min(pre_dS),
            "max": np.max(pre_dS),
            "median": np.median(pre_dS)
        },
        "post_S": {
            "n": len(post_dS),
            "mean": np.mean(post_dS),
            "min": np.min(post_dS),
            "max": np.max(post_dS),
            "median": np.median(post_dS)
        }
    }

    print("\n--- Summary Statistics ---")
    for k, v in stats.items():
        print(f"{k}: n={v['n']:,}, Min={v['min']:.6e}, Mean={v['mean']:.6e}, Max={v['max']:.6e}")

    # Create figure matching reference
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.97, top=0.92, bottom=0.14)

    color_pre = "#1f77b4"   # blue
    color_post = "#d95f02"  # orange/vermilion

    # =========================================================================
    # LEFT PANEL: Flip energy jump |Eafter - Ebefore|
    # =========================================================================
    bins_E = np.logspace(-8, 0, 29)
    h_pre_E, _ = np.histogram(pre_dE, bins=bins_E)
    h_post_E, _ = np.histogram(post_dE, bins=bins_E)

    prob_pre_E = h_pre_E / len(pre_dE)
    prob_post_E = h_post_E / len(post_dE)

    ax1.stairs(prob_pre_E, bins_E, color=color_pre, lw=2.2,
               label=r"pre-yield")
    ax1.stairs(prob_post_E, bins_E, color=color_post, lw=2.2,
               label=r"post-yield")

    # Reference threshold (Zanzotto single-element maximum barrier at gamma = 0.5)
    E_max_ref = 0.0481  # Single-element E_max at gamma = 0.5
    ax1.axvline(E_max_ref, color="#222222", ls="--", lw=1.3,
                label=r"single-element $E_{\mathrm{max}}$ ($\gamma = 0.5$)")

    ax1.set_xscale("log")
    ax1.set_xlim(1e-9, 1e0)
    ax1.set_ylim(0, max(np.max(prob_pre_E), np.max(prob_post_E)) * 1.25)
    ax1.set_xlabel(r"Nodal energy jump upon remeshing $|\Delta e_a^{(k)}|$")
    ax1.set_ylabel(r"Probability per logarithmic bin")
    ax1.legend(frameon=True, facecolor="white", edgecolor="#d0d0d0",
               framealpha=0.92, loc="upper left", fontsize=9.2)

    # Panel label (a) at bottom left
    ax1.text(0.04, 0.21, r"$\mathbf{(a)}$", transform=ax1.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left")

    # Subtitle / Stats callout box for Min, Mean, Max
    stats_text_E = (
        f"Pre : Min = {stats['pre_E']['min']:.1e}   Mean = {stats['pre_E']['mean']:.2e}   Max = {stats['pre_E']['max']:.3f}\n"
        f"Post: Min = {stats['post_E']['min']:.1e}   Mean = {stats['post_E']['mean']:.2e}   Max = {stats['post_E']['max']:.3f}"
    )
    ax1.text(0.04, 0.05, stats_text_E, transform=ax1.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

    # =========================================================================
    # RIGHT PANEL: Nodal stress jump upon remeshing |\Delta \sigma_{xy, a}^{(k)}|
    # =========================================================================
    bins_S = np.logspace(-6, 0, 25)
    h_pre_S, _ = np.histogram(pre_dS, bins=bins_S)
    h_post_S, _ = np.histogram(post_dS, bins=bins_S)

    prob_pre_S = h_pre_S / len(pre_dS)
    prob_post_S = h_post_S / len(post_dS)

    ax2.stairs(prob_pre_S, bins_S, color=color_pre, lw=2.2,
               label=r"pre-yield")
    ax2.stairs(prob_post_S, bins_S, color=color_post, lw=2.2,
               label=r"post-yield")

    # Reference threshold (Loss of ellipticity at gamma = 0.1322)
    S_ellipticity_ref = 0.334
    ax2.axvline(S_ellipticity_ref, color="#222222", ls="--", lw=1.3,
                label=r"$\sigma_{xy}$ at loss of ellipticity ($\gamma = 0.1322$)")

    ax2.set_xscale("log")
    ax2.set_xlim(5e-7, 1.2e0)
    ax2.set_ylim(0, max(np.max(prob_pre_S), np.max(prob_post_S)) * 1.25)
    ax2.set_xlabel(r"Nodal stress jump upon remeshing $|\Delta \sigma_{xy, a}^{(k)}|$")
    ax2.set_ylabel(r"Probability per logarithmic bin")
    ax2.legend(frameon=True, facecolor="white", edgecolor="#d0d0d0",
               framealpha=0.92, loc="upper left", fontsize=9.2)

    # Panel label (b) at bottom left
    ax2.text(0.04, 0.21, r"$\mathbf{(b)}$", transform=ax2.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left")

    # Subtitle / Stats callout box for Min, Mean, Max
    stats_text_S = (
        f"Pre : Min = {stats['pre_S']['min']:.1e}   Mean = {stats['pre_S']['mean']:.2e}   Max = {stats['pre_S']['max']:.3f}\n"
        f"Post: Min = {stats['post_S']['min']:.1e}   Mean = {stats['post_S']['mean']:.2e}   Max = {stats['post_S']['max']:.3f}"
    )
    ax2.text(0.04, 0.05, stats_text_S, transform=ax2.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

    out_dir = repo_root / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "flip_jump_distributions_comparison.png"
    out_pdf = out_dir / "flip_jump_distributions_comparison.pdf"

    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"\nFigure saved to:\n  {out_png}\n  {out_pdf}")

if __name__ == "__main__":
    main()
