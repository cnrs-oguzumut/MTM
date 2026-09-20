#!/usr/bin/env python3
"""
Publication Figures: Local Nodal Energy Change Distributions (Pre-Yield vs Post-Yield)
in Positive Continuous Loading.

Yield threshold is determined dynamically from the simulation data as the strain alpha
corresponding to the maximum post-optimization shear stress:
    alpha_yield = argmax_alpha( PostStress(alpha) )

Generates two separate publication-quality figures:
1. figures/nodal_distribution_pre_yield.{png,pdf}
2. figures/nodal_distribution_post_yield.{png,pdf}
Both figures prominently highlight:
  - Max (\\Delta e_max)
  - Min (\\Delta e_min)
  - Mean (\\langle \\Delta e_a \\rangle)
together with median, sample size, and energetic breakdown.
"""

import sys
import csv
import glob
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable

def setup_publication_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "figure.titlesize": 11.5,
        "lines.linewidth": 1.2,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })

import vtk
from vtk.util.numpy_support import vtk_to_numpy

def load_pts_tri_energy(vtk_path):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(vtk_path))
    r.Update()
    grid = r.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    ne = vtk_to_numpy(grid.GetPointData().GetArray("NodalEnergy"))
    
    # Get triangles connectivity
    try:
        cells_data = grid.GetCells()
        if hasattr(cells_data, "GetConnectivityArray"):
            tri = vtk_to_numpy(cells_data.GetConnectivityArray()).reshape(-1, 3)
        else:
            raw = vtk_to_numpy(cells_data.GetData())
            tri = raw.reshape(-1, 4)[:, 1:4]
    except Exception:
        tri = np.empty((0, 3), dtype=int)
        
    return pts, tri, ne

def find_yield_alpha(sim_dir):
    log_file = sim_dir / "energy_stress_log.csv"
    if not log_file.exists():
        raise FileNotFoundError(f"Cannot find {log_file}")
    with open(log_file, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Calculate argmax of PostStress
    max_row = max(rows, key=lambda r: float(r["PostStress"]))
    alpha_yield = float(max_row["Alpha"])
    peak_stress = float(max_row["PostStress"])
    return alpha_yield, peak_stress

def extract_nodal_data(dirs, alpha_yield):
    pre_yield_diffs = []
    post_yield_diffs = []
    spatial_pre = None
    spatial_post = None

    seen_passes = set()

    for d in dirs:
        d_path = Path(d)
        surgery_dir = d_path / "vtk_surgery"
        log_file = d_path / "energy_stress_log.csv"

        if not surgery_dir.exists():
            continue

        step_to_alpha = {}
        if log_file.exists():
            with open(log_file) as f:
                for r in csv.DictReader(f):
                    step_to_alpha[int(r["Iteration"])] = abs(float(r["Alpha"]))

        before_files = sorted(surgery_dir.glob("*_1_before_remesh.vtk"))
        print(f"Scanning {d_path.name}: found {len(before_files)} surgery passes...")

        for bf in before_files:
            af = Path(str(bf).replace("_1_before_remesh.vtk", "_2_after_remesh.vtk"))
            if not af.exists():
                continue

            fname = bf.name
            if fname in seen_passes:
                continue
            seen_passes.add(fname)

            try:
                step_idx = int(fname.split("_")[1])
            except:
                step_idx = 0

            alpha_val = step_to_alpha.get(step_idx, 0.14)

            pts1, tri1, ne1 = load_pts_tri_energy(bf)
            pts2, tri2, ne2 = load_pts_tri_energy(af)

            if len(pts1) == 0 or len(ne1) == 0 or len(pts2) == 0 or len(ne2) == 0:
                continue

            h2 = { (round(p[0], 3), round(p[1], 3)): ne2[i] for i, p in enumerate(pts2) }

            diffs = np.zeros(len(pts1))
            matched_mask = np.zeros(len(pts1), dtype=bool)
            for i, p in enumerate(pts1):
                k = (round(p[0], 3), round(p[1], 3))
                if k in h2:
                    diffs[i] = h2[k] - ne1[i]
                    matched_mask[i] = True

            active = diffs[matched_mask & (np.abs(diffs) > 1e-6)]
            if len(active) == 0:
                continue

            if alpha_val < alpha_yield:
                pre_yield_diffs.extend(active)
                if spatial_pre is None or len(active) > spatial_pre.get("count", 0):
                    spatial_pre = {
                        "pts": pts1, "tri": tri1, "diff": diffs,
                        "alpha": alpha_val, "step": step_idx, "count": len(active)
                    }
            else:
                post_yield_diffs.extend(active)
                if spatial_post is None or len(active) > spatial_post.get("count", 0):
                    spatial_post = {
                        "pts": pts1, "tri": tri1, "diff": diffs,
                        "alpha": alpha_val, "step": step_idx, "count": len(active)
                    }

    return np.array(pre_yield_diffs), np.array(post_yield_diffs), spatial_pre, spatial_post

def plot_single_regime_figure(diffs, spatial, regime_title, alpha_condition_tex, regime_tag, out_path_base, col_primary="#1f77b4", col_dark="#08519c"):
    setup_publication_style()
    out_path_base.parent.mkdir(parents=True, exist_ok=True)

    if len(diffs) == 0:
        print(f"Warning: No data for {regime_tag}")
        return

    # Compute statistics
    val_max = float(np.max(diffs))
    val_min = float(np.min(diffs))
    val_mean = float(np.mean(diffs))
    val_median = float(np.median(diffs))
    val_std = float(np.std(diffs))
    val_pos_pct = float(np.mean(diffs > 0) * 100)
    val_neg_pct = float(np.mean(diffs < 0) * 100)

    print(f"\n=======================================================")
    print(f"STATISTICS FOR REGIME: {regime_title} ({alpha_condition_tex})")
    print(f"=======================================================")
    print(f"  Active node events: {len(diffs):,}")
    print(f"  Max value:          {val_max:+.6f}")
    print(f"  Min value:          {val_min:+.6f}")
    print(f"  Mean value:         {val_mean:+.6f}")
    print(f"  Median value:       {val_median:+.6f}")
    print(f"  Std deviation:      {val_std:.6f}")
    print(f"  Barrier jump (>0):  {val_pos_pct:.2f}%")
    print(f"  Local release (<0): {val_neg_pct:.2f}%")
    print(f"=======================================================\n")

    fig, axs = plt.subplots(2, 2, figsize=(7.8, 6.4))
    fig.subplots_adjust(hspace=0.40, wspace=0.34, left=0.10, right=0.96, top=0.92, bottom=0.10)

    # -------------------------------------------------------------
    # Panel (a): Semi-Log Distribution with Max, Min, Mean Markers
    # -------------------------------------------------------------
    ax = axs[0, 0]
    hist_max_edge = min(max(val_max * 1.05, 0.05), 0.25)
    bins = np.linspace(val_min * 1.05, hist_max_edge, 70)
    ax.hist(diffs, bins=bins, density=True, color=col_primary, alpha=0.65, edgecolor=col_dark, lw=0.6)
    
    ax.set_yscale("log")
    ax.axvline(0.0, color="gray", ls="--", lw=0.9, alpha=0.6)
    
    # Highlight Mean, Max, Min with explicit labels
    ax.axvline(val_mean, color="#d95f02", ls="-", lw=1.8, label=rf"$\mathbf{{Mean}} = {val_mean:+.4f}$")
    ax.axvline(val_max, color="#b30000", ls="--", lw=1.4, label=rf"$\mathbf{{Max}} = {val_max:+.3f}$")
    ax.axvline(val_min, color="#7570b3", ls=":", lw=1.6, label=rf"$\mathbf{{Min}} = {val_min:+.3f}$")

    ax.set_xlabel(r"Local nodal energy change $\Delta e_a$")
    ax.set_ylabel(r"Probability density $p(\Delta e_a)$ (log)")
    ax.set_title(rf"(a) $\Delta e_a$ Distribution [{alpha_condition_tex}]", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.92, loc="upper right", fontsize=8.0)
    ax.set_ylim(1e-2, 1.5e3)
    ax.grid(True, which="both", ls=":", alpha=0.35)

    # -------------------------------------------------------------
    # Panel (b): Log-Log Magnitude Spectrum |\Delta e_a|
    # -------------------------------------------------------------
    ax = axs[0, 1]
    abs_vals = np.abs(diffs)
    log_bins = np.logspace(-6, np.log10(np.max(abs_vals)), 50)
    ax.hist(abs_vals, bins=log_bins, density=True, color=col_primary, alpha=0.7, edgecolor=col_dark, lw=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Magnitude $|\Delta e_a|$")
    ax.set_ylabel(r"Probability density $p(|\Delta e_a|)$")
    ax.set_title(rf"(b) Magnitude Spectrum [{alpha_condition_tex}]", fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.4)

    # Annotate prominent stats callout box in panel (b)
    stats_box = (
        f"Active Events: {len(diffs):,}\n"
        f"MAX:    {val_max:+.5f}\n"
        f"MIN:    {val_min:+.5f}\n"
        f"MEAN:   {val_mean:+.5f}\n"
        f"Median: {val_median:+.5f}\n"
        f"Std:    {val_std:.5f}"
    )
    ax.text(
        0.05, 0.08, stats_box, transform=ax.transAxes,
        fontsize=8.0, family="monospace", va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f7f7f7", edgecolor="#cccccc", lw=0.8)
    )

    # -------------------------------------------------------------
    # Panel (c): Representative Spatial Map \Delta e_a(x, y)
    # -------------------------------------------------------------
    ax = axs[1, 0]
    if spatial is not None:
        pts = spatial["pts"]
        tri = spatial["tri"]
        diff_field = spatial["diff"]
        alpha_s = spatial["alpha"]
        
        tri_obj = mtri.Triangulation(pts[:, 0], pts[:, 1], tri)
        vlim = 0.03
        im = ax.tripcolor(tri_obj, diff_field, shading="gouraud", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
        ax.set_aspect("equal")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_title(rf"(c) Spatial Surgery Zone ($\alpha = {alpha_s:.3f}$)", fontweight="bold", loc="left")

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(r"$\Delta e_a$", fontsize=9.0)
        cb.ax.tick_params(labelsize=8)
    else:
        ax.text(0.5, 0.5, "No spatial map available", ha="center", va="center")

    # -------------------------------------------------------------
    # Panel (d): Energy Character Breakdown
    # -------------------------------------------------------------
    ax = axs[1, 1]
    bars = ax.bar([0, 1], [val_pos_pct, val_neg_pct], width=0.5,
                  color=["#e6550d", "#31a354"], edgecolor="black", lw=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"Barrier Jump ($\Delta e_a > 0$)", r"Local Release ($\Delta e_a < 0$)"])
    ax.set_ylabel(r"Active node fraction (%)")
    ax.set_ylim(0, 75)
    ax.set_title(rf"(d) Energy Balance [{alpha_condition_tex}]", fontweight="bold", loc="left")

    for b, val in zip(bars, [val_pos_pct, val_neg_pct]):
        h = b.get_height()
        ax.text(b.get_x() + b.get_width()/2.0, h + 2.0, f"{val:.1f}%",
                ha="center", va="bottom", fontsize=9.0, fontweight="bold")

    png_path = out_path_base.with_suffix(".png")
    pdf_path = out_path_base.with_suffix(".pdf")
    plt.savefig(pdf_path, dpi=300)
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Saved {regime_title} figure to:\n  {png_path}\n  {pdf_path}")
    return png_path, pdf_path, {
        "count": len(diffs),
        "max": val_max,
        "min": val_min,
        "mean": val_mean,
        "median": val_median,
        "std": val_std,
        "pos_pct": val_pos_pct,
        "neg_pct": val_neg_pct
    }

if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parent
    primary_dir = repo_root / "study_100x100_positive"
    
    # 1. Calculate yield point dynamically from PostStress maximum
    alpha_yield, peak_stress = find_yield_alpha(primary_dir)
    print("=================================================================")
    print(f"DYNAMIC YIELD DETECTION (Max PostStress):")
    print(f"  alpha_yield = {alpha_yield:.5f}")
    print(f"  sigma_yield = {peak_stress:.5f}")
    print("=================================================================")

    # 2. Extract surgery passes from the positive loading run
    dirs = [primary_dir]
    diffs_pre, diffs_post, sp_pre, sp_post = extract_nodal_data(dirs, alpha_yield=alpha_yield)
    
    out_dir = repo_root / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 3. Figure 1: Pre-Yield Figure
    _, _, stats_pre = plot_single_regime_figure(
        diffs_pre, sp_pre,
        regime_title="Pre-Yield",
        alpha_condition_tex=rf"$\alpha < {alpha_yield:.3f}$",
        regime_tag="pre_yield",
        out_path_base=out_dir / "nodal_distribution_pre_yield",
        col_primary="#1f77b4", col_dark="#08519c"
    )

    # 4. Figure 2: Post-Yield Figure
    _, _, stats_post = plot_single_regime_figure(
        diffs_post, sp_post,
        regime_title="Post-Yield",
        alpha_condition_tex=rf"$\alpha \geq {alpha_yield:.3f}$",
        regime_tag="post_yield",
        out_path_base=out_dir / "nodal_distribution_post_yield",
        col_primary="#d95f02", col_dark="#a63603"
    )
