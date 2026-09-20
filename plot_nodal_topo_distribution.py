#!/usr/bin/env python3
r"""
plot_nodal_topo_distribution.py

Calculates and visualizes the distribution of LOCAL NODAL ENERGY CHANGES
\Delta e_a = e_a^{\text{after remesh}} - e_a^{\text{before remesh}}
across all avalanches and remeshing passes in 100x100 crystal lattices.

Features:
- Extracts exact nodal energy differences at frozen coordinates.
- Plots publication-grade 4-panel visualization:
  (a) Linear distribution PDF of \Delta e_a (centered around 0 with positive net mean)
  (b) Log-scale distribution of magnitude |\Delta e_a| showing heavy-tailed defect behavior
  (c) Spatial map of local nodal energy change \Delta e_a(x, y) across a representative avalanche
  (d) Breakdown of local barrier increases (\Delta e_a > 0) vs local releases (\Delta e_a < 0)
"""

import os
import sys
import glob
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable

def setup_publication_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "axes.titlesize": 10.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "mathtext.fontset": "cm",
        "lines.antialiased": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

def load_pts_tri_energy(vtk_path):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(vtk_path))
    r.Update()
    grid = r.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    ne = vtk_to_numpy(grid.GetPointData().GetArray('NodalEnergy'))
    cells = grid.GetCells()
    try:
        tri = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 3)
    except:
        tri = vtk_to_numpy(cells.GetData()).reshape(-1, 4)[:, 1:]
    return pts, tri, ne

def collect_nodal_data(surgery_dir):
    before_files = sorted(surgery_dir.glob("*_1_before_remesh.vtk"))
    all_diffs = []
    sample_spatial_data = None

    print(f"Processing {len(before_files)} surgery passes in {surgery_dir}...")
    for bf in before_files:
        af = Path(str(bf).replace('_1_before_remesh.vtk', '_2_after_remesh.vtk'))
        if not af.exists():
            continue
        pts1, tri1, ne1 = load_pts_tri_energy(bf)
        pts2, tri2, ne2 = load_pts_tri_energy(af)

        # Coordinate matching
        h2 = { (round(p[0], 3), round(p[1], 3)): ne2[i] for i, p in enumerate(pts2) }
        
        diffs = np.zeros(len(pts1))
        matched_mask = np.zeros(len(pts1), dtype=bool)
        for i, p in enumerate(pts1):
            k = (round(p[0], 3), round(p[1], 3))
            if k in h2:
                diffs[i] = h2[k] - ne1[i]
                matched_mask[i] = True

        active = diffs[matched_mask & (np.abs(diffs) > 1e-6)]
        all_diffs.extend(active)

        # Save step 165 pass 1 as the representative spatial example
        if "step_00165_pass_01" in str(bf):
            sample_spatial_data = {
                "pts": pts1,
                "tri": tri1,
                "diff": diffs,
                "name": "Step 165, Pass 1"
            }

    all_diffs = np.array(all_diffs)
    return all_diffs, sample_spatial_data

def plot_nodal_figure(all_diffs, sample_spatial, out_dir):
    setup_publication_style()
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(2, 2, figsize=(7.8, 6.4))
    fig.subplots_adjust(hspace=0.38, wspace=0.32, left=0.10, right=0.96, top=0.93, bottom=0.10)

    # -------------------------------------------------------------
    # Panel (a): Linear Distribution of Nodal Delta e
    # -------------------------------------------------------------
    ax = axs[0, 0]
    # Filter to [-0.05, 0.10] for clear visualization of the central peak
    bins = np.linspace(-0.04, 0.08, 61)
    ax.hist(all_diffs, bins=bins, density=True, color="#1f77b4", alpha=0.7, edgecolor="black", lw=0.6)

    # Pure numpy Gaussian KDE
    x_grid = np.linspace(-0.04, 0.08, 250)
    bw = 0.003
    dens = np.mean(
        np.exp(-0.5 * ((x_grid[:, None] - all_diffs[None, :]) / bw) ** 2)
        / (bw * np.sqrt(2 * np.pi)),
        axis=1
    )
    ax.plot(x_grid, dens, color="#08519c", lw=1.8, label="KDE density")

    mean_val = np.mean(all_diffs)
    median_val = np.median(all_diffs)
    ax.axvline(0.0, color="gray", ls="--", lw=1.0, alpha=0.7)
    ax.axvline(mean_val, color="#d95f02", ls="-", lw=1.6, label=f"Mean: +{mean_val:.4f}")
    ax.axvline(median_val, color="#2ca02c", ls=":", lw=1.6, label=f"Median: {median_val:.5f}")

    ax.set_xlabel(r"Local nodal energy change $\Delta e_a$")
    ax.set_ylabel(r"Probability density $p(\Delta e_a)$")
    ax.set_title(r"(a) Distribution of Local Nodal $\Delta e_a$", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")
    ax.set_xlim(-0.04, 0.08)

    # -------------------------------------------------------------
    # Panel (b): Log-scale magnitude |\Delta e_a|
    # -------------------------------------------------------------
    ax = axs[0, 1]
    abs_diffs = np.abs(all_diffs)
    log_bins = np.logspace(-6, np.log10(np.max(abs_diffs)), 50)
    ax.hist(abs_diffs, bins=log_bins, density=True, color="#3182bd", alpha=0.7, edgecolor="black", lw=0.6)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Magnitude $|\Delta e_a|$")
    ax.set_ylabel(r"Probability density $p(|\Delta e_a|)$")
    ax.set_title(r"(b) Magnitude Spectrum (Log-Log)", fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.5)

    # Annotate max value
    max_val = np.max(all_diffs)
    ax.axvline(max_val, color="#b30000", ls="--", lw=1.4, label=rf"Max $\Delta e_a$: {max_val:.3f}")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")

    # -------------------------------------------------------------
    # Panel (c): Spatial Map \Delta e_a(x, y)
    # -------------------------------------------------------------
    ax = axs[1, 0]
    if sample_spatial is not None:
        pts = sample_spatial["pts"]
        tri = sample_spatial["tri"]
        diff = sample_spatial["diff"]
        
        tri_obj = mtri.Triangulation(pts[:, 0], pts[:, 1], tri)
        vlim = 0.03
        im = ax.tripcolor(tri_obj, diff, shading="gouraud", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
        ax.set_aspect("equal")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_title(r"(c) Spatial Map $\Delta e_a(x,y)$ (Step 165, Pass 1)", fontweight="bold", loc="left")

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(r"$\Delta e_a$", fontsize=9.0)
        cb.ax.tick_params(labelsize=8)

    # -------------------------------------------------------------
    # Panel (d): Positive vs. Negative Fractions & Summary
    # -------------------------------------------------------------
    ax = axs[1, 1]
    pos_diffs = all_diffs[all_diffs > 0]
    neg_diffs = all_diffs[all_diffs < 0]

    f_pos = len(pos_diffs) / len(all_diffs) * 100
    f_neg = len(neg_diffs) / len(all_diffs) * 100

    bars = ax.bar(
        [0, 1], [f_pos, f_neg], width=0.5,
        color=["#e6550d", "#31a354"], edgecolor="black", lw=0.8
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"Barrier Increase ($\Delta e_a > 0$)", r"Local Relaxation ($\Delta e_a < 0$)"])
    ax.set_ylabel(r"Percentage of active nodes (%)")
    ax.set_ylim(0, 75)
    ax.set_title(r"(d) Local Reconnection Character", fontweight="bold", loc="left")

    for b, val in zip(bars, [f_pos, f_neg]):
        h = b.get_height()
        ax.text(b.get_x() + b.get_width()/2.0, h + 2.0, f"{val:.1f}%",
                ha="center", va="bottom", fontsize=9.0, fontweight="bold")

    # Inset stats box in panel d
    stats_text = (
        f"Active nodes: {len(all_diffs):,}\n"
        f"Max $\\Delta e_a$: +{np.max(all_diffs):.3f}\n"
        f"Min $\\Delta e_a$: {np.min(all_diffs):.3f}\n"
        f"Mean $\\Delta e_a$: +{np.mean(all_diffs):.5f}"
    )
    ax.text(
        0.52, 0.40, stats_text, transform=ax.transAxes,
        fontsize=8.0, family="monospace", va="center", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f7f7f7", edgecolor="#cccccc", lw=0.8)
    )

    pdf_path = out_dir / "nodal_energy_distribution_100x100.pdf"
    png_path = out_dir / "nodal_energy_distribution_100x100.png"
    plt.savefig(pdf_path, dpi=300)
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Saved publication figures to:\n  {pdf_path}\n  {png_path}")
    return png_path, pdf_path

if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parent
    surgery_dir = repo_root / "study_100x100_10avalanches" / "vtk_surgery"
    all_diffs, sample_spatial = collect_nodal_data(surgery_dir)
    out_dir = repo_root / "figures"
    plot_nodal_figure(all_diffs, sample_spatial, out_dir)
