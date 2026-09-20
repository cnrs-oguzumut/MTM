#!/usr/bin/env python3
r"""
plot_avalanche_morphology.py

Generates publication-quality figures illustrating the spatial evolution and
topological mechanics of plastic avalanches during minimization:

Option 1: Filmstrip sequence of 4 milestone states (A, B, C, D) showing \sigma_{xy}
Option 2: Difference maps showing net stress relaxation \Delta\sigma_{xy} and energy dissipation
Option 3: Microscopic zoom-in on the topological surgery zone showing edge reconnection
Option 4: Composite figure linking energy/stress minimization curves directly with spatial states
"""

import sys
import os
from pathlib import Path
import csv
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Style Setup
# -----------------------------------------------------------------------------
def setup_publication_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "axes.titlesize": 9.5,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 7.8,
        "mathtext.fontset": "cm",
        "lines.antialiased": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

# -----------------------------------------------------------------------------
# VTK Loading Helper
# -----------------------------------------------------------------------------
def load_vtk_data(filepath):
    path = str(filepath)
    if not os.path.exists(path):
        raise FileNotFoundError(f"VTK file not found: {path}")
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    grid = reader.GetOutput()
    
    # Points
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    
    # Cells (triangles)
    cells = grid.GetCells()
    try:
        conn = vtk_to_numpy(cells.GetConnectivityArray())
        triangles = conn.reshape(-1, 3)
    except:
        conn = vtk_to_numpy(cells.GetData())
        triangles = conn.reshape(-1, 4)[:, 1:]
        
    # Point data
    pd = grid.GetPointData()
    s_tensor = vtk_to_numpy(pd.GetArray("NodalCauchyStress"))
    s_xy = s_tensor[:, 1]
    s_xx = s_tensor[:, 0]
    s_yy = s_tensor[:, 4]
    
    nodal_e = vtk_to_numpy(pd.GetArray("NodalEnergy"))
    
    return {
        "points": pts,
        "triangles": triangles,
        "s_xy": s_xy,
        "s_xx": s_xx,
        "s_yy": s_yy,
        "energy": nodal_e,
        "bounds": grid.GetBounds()
    }

# -----------------------------------------------------------------------------
# Option 1: 4-Panel Filmstrip Sequence
# -----------------------------------------------------------------------------
def plot_option1_sequence(surgery_dir, output_base):
    setup_publication_style()
    files = [
        ("step_00165_pass_01_1_before_remesh.vtk", r"$\mathbf{A}$ ($k=395$): Pre-avalanche"),
        ("step_00165_pass_01_3_relaxed_accepted.vtk", r"$\mathbf{B}$ ($k=770$): Pass 1"),
        ("step_00165_pass_02_3_relaxed_accepted.vtk", r"$\mathbf{C}$ ($k=1015$): Pass 2"),
        ("step_00165_pass_04_3_relaxed_accepted.vtk", r"$\mathbf{D}$ ($k=1381$): Final")
    ]
    labels = ["(a)", "(b)", "(c)", "(d)"]
    
    datasets = [load_vtk_data(surgery_dir / f[0]) for f in files]
    
    vmin = -0.28
    vmax = 0.28
    
    fig, axes = plt.subplots(1, 4, figsize=(6.8, 2.25), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.08, "left": 0.07, "right": 0.89, "top": 0.86, "bottom": 0.18})
    
    im_ref = None
    for i, (ax, data, (fname, title), lab) in enumerate(zip(axes, datasets, files, labels)):
        triang = mtri.Triangulation(data["points"][:, 0], data["points"][:, 1], data["triangles"])
        im = ax.tripcolor(triang, data["s_xy"], shading="gouraud", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        im_ref = im
        ax.set_aspect("equal")
        ax.set_xlim(-18, 105)
        ax.set_ylim(-3, 103)
        ax.set_title(title, fontsize=7.8, pad=3)
        
        # Subplot letter aligned at bottom-left
        ax.text(0.05, 0.06, lab, transform=ax.transAxes, fontsize=8.5, fontweight="bold",
                bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"))
        
        ax.set_xlabel(r"$x$", labelpad=1)
        if i == 0:
            ax.set_ylabel(r"$y$", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
            
    # Shared vertical colorbar
    cbar_ax = fig.add_axes([0.905, 0.20, 0.015, 0.64])
    cb = fig.colorbar(im_ref, cax=cbar_ax)
    cb.set_label(r"Shear stress $\sigma_{xy}$", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_ticks([-0.2, -0.1, 0.0, 0.1, 0.2])
    
    pdf_path = output_base.with_name("step_00165_avalanche_sequence.pdf")
    png_path = output_base.with_name("step_00165_avalanche_sequence.png")
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"Saved Option 1: {pdf_path} and {png_path}")

# -----------------------------------------------------------------------------
# Option 2: Avalanche Difference Maps
# -----------------------------------------------------------------------------
def plot_option2_diffmap(surgery_dir, output_base):
    setup_publication_style()
    data_init = load_vtk_data(surgery_dir / "step_00165_pass_01_1_before_remesh.vtk")
    data_final = load_vtk_data(surgery_dir / "step_00165_pass_04_3_relaxed_accepted.vtk")
    
    p1 = data_init["points"]
    t1 = data_init["triangles"]
    s1 = data_init["s_xy"]
    e1 = data_init["energy"]
    
    p4 = data_final["points"]
    t4 = data_final["triangles"]
    s4 = data_final["s_xy"]
    e4 = data_final["energy"]
    
    # Interpolate initial fields at final positions
    tri1 = mtri.Triangulation(p1[:, 0], p1[:, 1], t1)
    tri4 = mtri.Triangulation(p4[:, 0], p4[:, 1], t4)
    
    interp_s1 = mtri.LinearTriInterpolator(tri1, s1)
    interp_e1 = mtri.LinearTriInterpolator(tri1, e1)
    
    s1_at_p4 = interp_s1(p4[:, 0], p4[:, 1])
    e1_at_p4 = interp_e1(p4[:, 0], p4[:, 1])
    
    diff_s = np.where(np.isnan(s4 - s1_at_p4), 0.0, s4 - s1_at_p4)
    diff_e = np.where(np.isnan(e4 - e1_at_p4), 0.0, e4 - e1_at_p4)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 3.2), sharey=True,
                                   gridspec_kw={"wspace": 0.35, "left": 0.08, "right": 0.88, "top": 0.88, "bottom": 0.15})
    
    # (a) Net stress change
    ax1.set_aspect("equal")
    im1 = ax1.tripcolor(tri4, diff_s, shading="gouraud", cmap="RdBu_r", vmin=-0.22, vmax=0.22)
    ax1.set_xlim(-18, 105)
    ax1.set_ylim(-3, 103)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax1.set_title(r"Net stress change $\Delta \sigma_{xy}$", fontsize=9.5, pad=4)
    ax1.text(0.05, 0.06, "(a)", transform=ax1.transAxes, fontsize=10, fontweight="bold",
             bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"))
    
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.08)
    cb1 = fig.colorbar(im1, cax=cax1)
    cb1.set_label(r"$\Delta \sigma_{xy}$", fontsize=9.0)
    cb1.ax.tick_params(labelsize=8)
    cb1.set_ticks([-0.2, -0.1, 0.0, 0.1, 0.2])
    
    # (b) Net energy change
    ax2.set_aspect("equal")
    e_bound = 0.025
    im2 = ax2.tripcolor(tri4, diff_e, shading="gouraud", cmap="PRGn_r", vmin=-e_bound, vmax=e_bound)
    ax2.set_xlim(-18, 105)
    ax2.set_ylim(-3, 103)
    ax2.set_xlabel(r"$x$")
    ax2.set_title(r"Energy change $\Delta e = e^{\rm final} - e^{\rm initial}$", fontsize=9.5, pad=4)
    ax2.text(0.05, 0.06, "(b)", transform=ax2.transAxes, fontsize=10, fontweight="bold",
             bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"))
    ax2.tick_params(labelleft=False)
    
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="5%", pad=0.08)
    cb2 = fig.colorbar(im2, cax=cax2)
    cb2.set_label(r"$\Delta e$ (nodal energy)", fontsize=9.0)
    cb2.ax.tick_params(labelsize=8)
    
    pdf_path = output_base.with_name("step_00165_avalanche_diffmap.pdf")
    png_path = output_base.with_name("step_00165_avalanche_diffmap.png")
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"Saved Option 2: {pdf_path} and {png_path}")

# -----------------------------------------------------------------------------
# Option 3: Mesh Topology Surgery Zoom
# -----------------------------------------------------------------------------
def plot_option3_zoom_surgery(surgery_dir, output_base):
    setup_publication_style()
    data1 = load_vtk_data(surgery_dir / "step_00165_pass_01_1_before_remesh.vtk")
    data2 = load_vtk_data(surgery_dir / "step_00165_pass_01_2_after_remesh.vtk")
    data4 = load_vtk_data(surgery_dir / "step_00165_pass_04_3_relaxed_accepted.vtk")
    
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.7), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.10, "left": 0.08, "right": 0.89, "top": 0.87, "bottom": 0.16})
    
    xlim = (41, 51)
    ylim = (47, 57)
    vmin = -0.28
    vmax = 0.28
    
    panels = [
        (data1, r"$\mathbf{A}$: Pre-reconnection ($k=395$)", "(a)"),
        (data2, r"$\mathbf{A}^*$: Reconnected ($k=396$)", "(b)"),
        (data4, r"$\mathbf{D}$: Final relaxed ($k=1381$)", "(c)"),
    ]
    
    im_last = None
    for i, (ax, (data, title, lab)) in enumerate(zip(axes, panels)):
        triang = mtri.Triangulation(data["points"][:, 0], data["points"][:, 1], data["triangles"])
        im = ax.tripcolor(triang, data["s_xy"], shading="gouraud", cmap="RdBu_r", vmin=vmin, vmax=vmax, alpha=0.85)
        im_last = im
        ax.triplot(triang, color="black", lw=0.55, alpha=0.75)
        ax.scatter(data["points"][:, 0], data["points"][:, 1], color="black", s=6, zorder=5)
        
        # Highlight reconnected edge in panel (b)
        if i == 1:
            ax.plot([46.14, 46.95], [50.30, 50.81], color="crimson", lw=2.4, zorder=6,
                    label=r"Flipped edge ($\Delta E_{\rm topo}$)")
            ax.scatter([46.14, 46.95], [50.30, 50.81], color="crimson", s=18, zorder=7)
            ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none",
                      framealpha=0.9, fontsize=7.2, handlelength=1.4)
            
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=8.5, pad=3)
        ax.set_xlabel(r"$x$", labelpad=1)
        ax.text(0.05, 0.06, lab, transform=ax.transAxes, fontsize=9.5, fontweight="bold",
                bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"))
        if i == 0:
            ax.set_ylabel(r"$y$", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
            
    cbar_ax = fig.add_axes([0.905, 0.18, 0.015, 0.67])
    cb = fig.colorbar(im_last, cax=cbar_ax)
    cb.set_label(r"$\sigma_{xy}$", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    cb.set_ticks([-0.2, -0.1, 0.0, 0.1, 0.2])
    
    pdf_path = output_base.with_name("step_00165_avalanche_zoom_surgery.pdf")
    png_path = output_base.with_name("step_00165_avalanche_zoom_surgery.png")
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"Saved Option 3: {pdf_path} and {png_path}")

# -----------------------------------------------------------------------------
# Option 4: Composite Figure (Minimization Curve + 4 Spatial Snapshots)
# -----------------------------------------------------------------------------
def plot_option4_composite(csv_path, surgery_dir, output_base):
    setup_publication_style()
    csv_file = Path(csv_path)
    rows = []
    with open(csv_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
            
    x = np.array([int(r["global_micro_step"]) for r in rows])
    energy = np.array([float(r["energy"]) for r in rows])
    stress = np.array([float(r["stress"]) for r in rows])
    grad_norm = np.array([float(r["grad_norm"]) for r in rows])
    event_type = [r["event_type"] for r in rows]
    phase = [r["phase"] for r in rows]
    
    # 4 VTK milestone datasets
    files = [
        ("step_00165_pass_01_1_before_remesh.vtk", r"$\mathbf{A}$: $k=395$"),
        ("step_00165_pass_01_3_relaxed_accepted.vtk", r"$\mathbf{B}$: $k=770$"),
        ("step_00165_pass_02_3_relaxed_accepted.vtk", r"$\mathbf{C}$: $k=1015$"),
        ("step_00165_pass_04_3_relaxed_accepted.vtk", r"$\mathbf{D}$: $k=1381$")
    ]
    datasets = [load_vtk_data(surgery_dir / f[0]) for f in files]
    milestone_k = [395, 770, 1015, 1381]
    letters = ["A", "B", "C", "D"]
    
    mask_s = stress != 0.0
    x_s = x[mask_s]
    s_vals = stress[mask_s]
    
    fig = plt.figure(figsize=(6.8, 4.8))
    subfigs = fig.subfigures(2, 1, height_ratios=[1.15, 1.0], hspace=0.05)
    
    # =========================================================================
    # Subfig 0: Curves (2 columns: Energy on left, Stress on right)
    # =========================================================================
    axes_top = subfigs[0].subplots(1, 2, gridspec_kw={"wspace": 0.28, "left": 0.08, "right": 0.96, "top": 0.86, "bottom": 0.18})
    ax_e, ax_s = axes_top[0], axes_top[1]
    
    # --- Top Left: Energy ---
    col_e = "#1b4d3e"
    ax_e.plot(x, energy, color=col_e, lw=1.2, label=r"$E(k)$", zorder=3)
    
    # Mark topological jumps (\Delta E_topo)
    reconn_idx = [i for i, ev in enumerate(event_type) if ev == "AFTER_REMESH"]
    for idx in reconn_idx:
        k_step = x[idx]
        e_pre = energy[idx-1]
        e_post = energy[idx]
        ax_e.plot([k_step, k_step], [e_pre, e_post], color="darkorange", lw=1.4, zorder=4)
        ax_e.scatter([k_step], [e_post], color="darkorange", s=14, zorder=5)
        
    # Mark final accepted state with single solid red dot
    final_idx = [i for i, ev in enumerate(event_type) if ev == "REMESH_ACCEPTED"][-1]
    final_k = x[final_idx]
    ax_e.scatter([final_k], [energy[final_idx]], color="red", s=28, zorder=6)
    
    # Annotate milestones A, B, C, D on Energy curve
    text_offsets_e = {
        "A": (0, 7),
        "B": (0, -13),
        "C": (0, -13),
        "D": (-12, 7)
    }
    for k_val, let in zip(milestone_k, letters):
        e_val = energy[np.where(x == k_val)[0][0]]
        ax_e.scatter([k_val], [e_val], color="blue", marker="o", s=20, zorder=7)
        dx, dy = text_offsets_e[let]
        ax_e.annotate(r"$\mathbf{" + let + "}$", (k_val, e_val),
                      textcoords="offset points", xytext=(dx, dy),
                      ha="center", fontsize=8.5, fontweight="bold", color="navy")
        
    ax_e.set_ylabel(r"Internal energy $E$", labelpad=2)
    ax_e.set_xlabel(r"Inner minimization step $k$", labelpad=2)
    ax_e.grid(True, linestyle=":", alpha=0.5, color="gray")
    ax_e.text(-0.06, -0.155, "(a)", transform=ax_e.transAxes, fontsize=10, fontweight="bold", va="center")
    
    # Top Left Inset: Zoom on reconnection jump
    ax_ins_e = ax_e.inset_axes([0.48, 0.46, 0.48, 0.48])
    k_first_jump = x[reconn_idx[0]]
    ins_mask = (x >= k_first_jump - 15) & (x <= k_first_jump + 80)
    ax_ins_e.plot(x[ins_mask], energy[ins_mask], color=col_e, lw=1.2)
    ax_ins_e.plot([k_first_jump, k_first_jump], [energy[reconn_idx[0]-1], energy[reconn_idx[0]]],
                  color="darkorange", lw=1.5)
    ax_ins_e.scatter([k_first_jump], [energy[reconn_idx[0]]], color="darkorange", s=14, zorder=5)
    ax_ins_e.annotate(r"$\Delta E_{\rm topo}$",
                      xy=(k_first_jump, 0.5 * (energy[reconn_idx[0]-1] + energy[reconn_idx[0]])),
                      xytext=(k_first_jump + 15, 0.5 * (energy[reconn_idx[0]-1] + energy[reconn_idx[0]])),
                      arrowprops=dict(arrowstyle="->", color="darkorange", lw=0.8),
                      fontsize=7.2, color="darkorange", va="center")
    ax_ins_e.tick_params(labelsize=6.5, direction="in")
    ax_ins_e.grid(True, linestyle=":", alpha=0.4)
    
    # --- Top Right: Stress ---
    col_s = "#003366"
    ax_s.plot(x_s, s_vals, color=col_s, lw=1.2, label=r"$\sigma(k)$", zorder=3)
    ax_s.scatter([final_k], [s_vals[-1]], color="red", s=28, zorder=6)
    
    text_offsets_s = {
        "A": (0, 7),
        "B": (10, 5),
        "C": (-8, 7),
        "D": (10, 5)
    }
    for k_val, let in zip(milestone_k, letters):
        idx_near = np.argmin(np.abs(x_s - k_val))
        s_val = s_vals[idx_near]
        ax_s.scatter([x_s[idx_near]], [s_val], color="blue", marker="o", s=20, zorder=7)
        dx, dy = text_offsets_s[let]
        ax_s.annotate(r"$\mathbf{" + let + "}$", (x_s[idx_near], s_val),
                      textcoords="offset points", xytext=(dx, dy),
                      ha="center", fontsize=8.5, fontweight="bold", color="navy")
        
    ax_s.set_ylabel(r"Shear stress $\sigma_{xy}$", labelpad=2)
    ax_s.set_xlabel(r"Inner minimization step $k$", labelpad=2)
    ax_s.grid(True, linestyle=":", alpha=0.5, color="gray")
    ax_s.text(-0.06, -0.155, "(b)", transform=ax_s.transAxes, fontsize=10, fontweight="bold", va="center")
    
    # Top Right Inset: Gradient norm (raised to upper-left clear of the curve)
    ax_ins_s = ax_s.inset_axes([0.16, 0.52, 0.42, 0.38])
    grad_mask = (grad_norm > 0) & np.isfinite(grad_norm)
    if np.any(grad_mask):
        ax_ins_s.semilogy(x[grad_mask], grad_norm[grad_mask], color="#8b0000", lw=0.85)
        ax_ins_s.set_ylabel(r"$\|\nabla E\|_\infty$", fontsize=7.2, labelpad=1)
        ax_ins_s.set_title(r"Residual $\|\nabla E\|_\infty$", fontsize=7.2, pad=2)
        ax_ins_s.tick_params(labelsize=6.5, direction="in")
        ax_ins_s.grid(True, linestyle=":", alpha=0.4)
        
    # Top header legend
    leg_elems = [
        Line2D([0], [0], color=col_e, lw=1.2, label="Minimization"),
        Line2D([0], [0], color="darkorange", lw=1.4, marker="o", markersize=3.5, label=r"$\Delta E_{\rm topo}$ reconnection"),
        Line2D([0], [0], color="blue", marker="o", ls="none", markersize=4.5, label=r"Surgery states $\mathbf{A}$--$\mathbf{D}$"),
        Line2D([0], [0], color="red", marker="o", ls="none", markersize=4.5, label="Final accepted state")
    ]
    subfigs[0].legend(handles=leg_elems, loc="upper center", bbox_to_anchor=(0.5, 1.01),
                      ncol=4, frameon=False, fontsize=7.5)

    # =========================================================================
    # Subfig 1: 4 Spatial Snapshots (Bottom Row)
    # =========================================================================
    axes_bot = subfigs[1].subplots(1, 4, sharex=True, sharey=True,
                                   gridspec_kw={"wspace": 0.08, "left": 0.07, "right": 0.89, "top": 0.86, "bottom": 0.18})
    sub_labels = ["(c)", "(d)", "(e)", "(f)"]
    vmin = -0.28
    vmax = 0.28
    im_spatial = None
    
    for i, (ax_m, (fname, title), lab) in enumerate(zip(axes_bot, files, sub_labels)):
        data = datasets[i]
        triang = mtri.Triangulation(data["points"][:, 0], data["points"][:, 1], data["triangles"])
        im_spatial = ax_m.tripcolor(triang, data["s_xy"], shading="gouraud", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax_m.set_aspect("equal")
        ax_m.set_xlim(-18, 105)
        ax_m.set_ylim(-3, 103)
        ax_m.set_title(title, fontsize=8.0, pad=3)
        ax_m.set_xlabel(r"$x$", labelpad=1)
        
        ax_m.text(0.05, 0.06, lab, transform=ax_m.transAxes, fontsize=8.5, fontweight="bold",
                  bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"))
        
        if i == 0:
            ax_m.set_ylabel(r"$y$", labelpad=1)
        else:
            ax_m.tick_params(labelleft=False)
            
    # Shared colorbar for bottom spatial snapshots
    cbar_ax = subfigs[1].add_axes([0.905, 0.20, 0.015, 0.64])
    cb = subfigs[1].colorbar(im_spatial, cax=cbar_ax)
    cb.set_label(r"$\sigma_{xy}$", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_ticks([-0.2, -0.1, 0.0, 0.1, 0.2])
    
    pdf_path = output_base.with_name("step_00165_avalanche_composite.pdf")
    png_path = output_base.with_name("step_00165_avalanche_composite.png")
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"Saved Option 4: {pdf_path} and {png_path}")

# -----------------------------------------------------------------------------
# Main Runner
# -----------------------------------------------------------------------------
def main():
    base_dir = Path("study_100x100_10avalanches/vtk_surgery")
    csv_path = base_dir / "step_00165_avalanche_surgery.csv"
    output_base = base_dir / "step_00165"
    
    print("Generating Option 1: 4-Panel Filmstrip Sequence...")
    plot_option1_sequence(base_dir, output_base)
    
    print("Generating Option 2: Avalanche Difference Maps...")
    plot_option2_diffmap(base_dir, output_base)
    
    print("Generating Option 3: Mesh Topology Surgery Zoom...")
    plot_option3_zoom_surgery(base_dir, output_base)
    
    print("Generating Option 4: Complete Composite Figure...")
    plot_option4_composite(csv_path, base_dir, output_base)
    
    print("All 4 options successfully generated!")

if __name__ == "__main__":
    main()
