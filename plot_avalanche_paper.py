#!/usr/bin/env python3
r"""
plot_avalanche_paper.py

Publication-ready visualization of avalanche micro-surgery traces tailored for
academic journals (Physical Review, Acta Materialia, JMPS).

Features:
- Sized for standard journal column widths (single-column: 3.37 in, double-column: 6.8 in)
- Computer Modern math & serif typography (mathtext.fontset = 'cm')
- Boxed axes with inward ticks on all four sides (direction='in', top=True, right=True)
- Minimalist editorial aesthetics (no pastel rainbow background stripes)
- Discrete topological reconnection jumps marked clearly (\Delta E_topo)
- Vector PDF export (for LaTeX \includegraphics) and 300+ DPI PNG preview
"""

import sys
import argparse
from pathlib import Path
import csv
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def setup_publication_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "mathtext.fontset": "cm",
        "lines.antialiased": True,
    })

def plot_paper_figure(csv_path, output_base=None, layout="stacked", col_width="double", forces="inset"):
    """
    layout: 'stacked' (2 rows: Energy on top, Stress on bottom)
            'sidebyside' (1 row, 2 cols: Energy on left, Stress on right)
    col_width: 'double' (~6.8 inches / 17.3 cm) or 'single' (~3.37 inches / 8.6 cm)
    forces: 'inset' (inner inset plot inside panel b, default)
            'panel' (dedicated 3rd panel for residual force)
            'none'  (no force plot)
    """
    setup_publication_style()
    csv_file = Path(csv_path)
    if not csv_file.exists():
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)

    rows = []
    with open(csv_file, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    if not rows:
        print(f"Error: Empty CSV: {csv_file}")
        sys.exit(1)

    load_step = rows[0]['load_step']
    x = np.array([int(r['global_micro_step']) for r in rows])
    energy = np.array([float(r['energy']) for r in rows])
    stress = np.array([float(r['stress']) for r in rows])
    grad_norm = np.array([float(r['grad_norm']) for r in rows])
    event_type = [r['event_type'] for r in rows]
    phase = [r['phase'] for r in rows]

    # Detect phases
    phase_boundaries = []
    curr_phase = phase[0]
    start_step = x[0]
    for i in range(1, len(x)):
        if phase[i] != curr_phase:
            phase_boundaries.append((curr_phase, start_step, x[i-1]))
            curr_phase = phase[i]
            start_step = x[i]
    phase_boundaries.append((curr_phase, start_step, x[-1]))

    # Key events
    reconn_idx = [i for i, ev in enumerate(event_type) if ev == 'AFTER_REMESH']
    accepted_idx = [i for i, ev in enumerate(event_type) if ev == 'REMESH_ACCEPTED']
    rejected_idx = [i for i, ev in enumerate(event_type) if ev == 'REMESH_REJECTED']

    # Non-zero stress checkpoints
    mask_s = stress != 0.0
    x_s = x[mask_s]
    s_vals = stress[mask_s]

    mask_g = grad_norm > 0.0

    # Colors & Markers
    col_energy = "#154360"    # Deep navy for minimization path
    col_reconn = "#2e86c1"    # Cyan/blue dotted stem for topological jump
    col_accept = "#c0392b"    # Solid red dot for accepted state
    col_reject = "#7f8c8d"    # Gray cross for rejected state
    col_stress = "#78281f"    # Wine red for stress
    col_force  = "#5b2c6f"    # Deep purple for inner residual force

    # Dimensions & Axis Setup
    if forces == "panel":
        if layout == "sidebyside":
            fig_w = 7.2
            fig_h = 2.7
            fig, (ax_e, ax_s, ax_f) = plt.subplots(1, 3, figsize=(fig_w, fig_h), dpi=300)
        else:
            fig_w = 6.8 if col_width == "double" else 3.37
            fig_h = 5.2 if col_width == "double" else 5.8
            fig, (ax_e, ax_s, ax_f) = plt.subplots(3, 1, figsize=(fig_w, fig_h), sharex=True, dpi=300)
        all_axes = (ax_e, ax_s, ax_f)
    else:
        if layout == "sidebyside":
            fig_w = 6.8 if col_width == "double" else 5.5
            fig_h = 2.7
            fig, (ax_e, ax_s) = plt.subplots(1, 2, figsize=(fig_w, fig_h), dpi=300)
        else:
            fig_w = 6.8 if col_width == "double" else 3.37
            fig_h = 3.8 if col_width == "double" else 4.2
            fig, (ax_e, ax_s) = plt.subplots(2, 1, figsize=(fig_w, fig_h), sharex=True, dpi=300)
        all_axes = (ax_e, ax_s)

    for ax in all_axes:
        ax.tick_params(direction="in", top=True, right=True, which="both")
        ax.tick_params(which="major", length=4.5, width=0.7)
        ax.tick_params(which="minor", length=2.2, width=0.5)

    # Headroom for labels and legend
    y_min_e, y_max_e = np.min(energy), np.max(energy)
    y_range_e = y_max_e - y_min_e
    ax_e.set_ylim(y_min_e - 0.04 * y_range_e, y_max_e + 0.18 * y_range_e)

    # --- Panel (a): Energy Trajectory ---
    ax_e.plot(x, energy, color=col_energy, lw=1.25, zorder=3, label="Minimization path")

    # Topological reconnection jumps (Delta E_topo)
    for idx in reconn_idx:
        ax_e.plot([x[idx-1], x[idx]], [energy[idx-1], energy[idx]],
                  color=col_reconn, lw=1.0, ls=":", zorder=4)
        ax_e.plot(x[idx], energy[idx], marker="^", color=col_reconn,
                  markerfacecolor="white", markeredgewidth=1.1, markersize=4.5, zorder=5)

    # Accepted states (RED DOTS)
    for idx in accepted_idx:
        ax_e.plot(x[idx], energy[idx], marker="o", color=col_accept,
                  markersize=4.5, zorder=6)

    # Rejected states
    for idx in rejected_idx:
        ax_e.plot(x[idx], energy[idx], marker="x", color=col_reject,
                  markersize=5.0, mew=1.5, zorder=6)

    ax_e.set_ylabel(r"Internal energy $E$")
    ax_e.text(0.025, 0.92, "(a)", transform=ax_e.transAxes,
              fontsize=11, fontweight="bold", va="top")

    # Phase dividers and top labels
    has_short_phase = any((pb[2] - pb[1]) < 250 for pb in phase_boundaries)
    for p_idx, (p_name, p_start, p_end) in enumerate(phase_boundaries):
        mid_x = 0.5 * (p_start + p_end)
        if layout == "sidebyside" or has_short_phase:
            p_label = f"P{p_idx}" if p_idx > 0 else "Init"
        else:
            p_label = f"Pass {p_idx}" if p_idx > 0 else "Initial"
        ax_e.text(mid_x, y_max_e + 0.05 * y_range_e, p_label,
                  ha="center", va="bottom", fontsize=7.5, color="#555555")
        if p_idx < len(phase_boundaries) - 1:
            for ax in all_axes:
                ax.axvline(p_end, color="#cccccc", ls="--", lw=0.6, zorder=1)

    # Clean, horizontal legend placed OUTSIDE at the top
    legend_elements = [
        Line2D([0], [0], color=col_energy, lw=1.25, label="Minimization"),
        Line2D([0], [0], marker="^", color=col_reconn, markerfacecolor="white", ls=":", lw=1.0, markersize=4.5, label=r"$\Delta E_{\rm topo}$ reconnection"),
        Line2D([0], [0], marker="o", color=col_accept, ls="none", markersize=4.5, label="Accepted state"),
    ]
    if rejected_idx:
        legend_elements.append(
            Line2D([0], [0], marker="x", color=col_reject, ls="none", markersize=5.0, mew=1.5, label="Rejected state")
        )
    if forces in ("inset", "panel"):
        legend_elements.append(
            Line2D([0], [0], color=col_force, lw=1.0, label=r"Gradient $\|\nabla E\|_\infty$")
        )

    fig.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 0.99),
               ncol=len(legend_elements), frameon=False, fontsize=7.5)

    # --- Panel (b): Stress Evolution ---
    if len(x_s) > 0:
        ax_s.plot(x_s, s_vals, color=col_stress, lw=1.2, marker="s",
                  markersize=3.5, markerfacecolor="white", markeredgewidth=1.0, zorder=3)
    else:
        ax_s.plot(x, stress, color=col_stress, lw=1.2, zorder=3)

    ax_s.set_ylabel(r"Shear stress $\sigma_{xy}$")
    ax_s.text(0.025, 0.92, "(b)", transform=ax_s.transAxes,
              fontsize=11, fontweight="bold", va="top")

    # --- Inner Force / Gradient Inset inside Panel (b) ---
    if forces == "inset" and mask_g.any():
        if layout == "sidebyside":
            # Upper-left of panel (b) avoids the stress curve with proper headroom
            ax_ins = ax_s.inset_axes([0.18, 0.42, 0.42, 0.40])
        else:
            # Lower-left of stacked panel (b)
            ax_ins = ax_s.inset_axes([0.18, 0.42, 0.35, 0.40])
        ax_ins.semilogy(x[mask_g], grad_norm[mask_g], color=col_force, lw=0.85, zorder=3)
        ax_ins.tick_params(direction="in", top=True, right=True, which="both", labelsize=6.5)
        ax_ins.tick_params(which="major", length=2.5, width=0.5)
        ax_ins.tick_params(which="minor", length=1.2, width=0.4)
        ax_ins.set_ylabel(r"$\|\nabla E\|_\infty$", fontsize=7.5, labelpad=1)
        ax_ins.set_xlabel(r"$k$", fontsize=7.5, labelpad=1)
        ax_ins.set_title(r"Residual $\|\nabla E\|_\infty$", fontsize=7.5, pad=2)
        ax_ins.grid(True, ls=":", lw=0.4, color="#bbbbbb")
    elif forces == "panel":
        if mask_g.any():
            ax_f.semilogy(x[mask_g], grad_norm[mask_g], color=col_force, lw=1.0, zorder=3)
        ax_f.set_ylabel(r"Residual $\|\nabla E\|_\infty$")
        ax_f.text(0.025, 0.92, "(c)", transform=ax_f.transAxes,
                  fontsize=11, fontweight="bold", va="top")
        ax_f.grid(True, ls=":", lw=0.5, color="#dddddd")

    xlabel_text = r"Inner minimization step $k$"
    if forces == "panel":
        if layout == "sidebyside":
            for ax in all_axes:
                ax.set_xlabel(xlabel_text)
            fig.tight_layout(rect=[0.01, 0, 0.98, 0.92], pad=0.6)
        else:
            ax_f.set_xlabel(xlabel_text)
            fig.tight_layout(rect=[0, 0, 1, 0.94], pad=0.5)
    else:
        if layout == "sidebyside":
            ax_e.set_xlabel(xlabel_text)
            ax_s.set_xlabel(xlabel_text)
            fig.tight_layout(rect=[0, 0, 1, 0.92], pad=0.6)
        else:
            ax_s.set_xlabel(xlabel_text)
            fig.tight_layout(rect=[0, 0, 1, 0.94], pad=0.5)

    # Output naming
    if output_base is None:
        out_stem = csv_file.parent / f"{csv_file.stem}_paper_{layout}"
    else:
        out_stem = Path(output_base)

    out_png = out_stem.with_suffix(".png")
    out_pdf = out_stem.with_suffix(".pdf")

    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"✓ Saved publication figure: {out_png}")
    print(f"✓ Saved vector PDF:         {out_pdf}")
    return out_png, out_pdf

def main():
    parser = argparse.ArgumentParser(description="Generate publication-style avalanche surgery figures for papers.")
    parser.add_argument("csv_file", type=str, help="Path to avalanche_trace CSV file or directory")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output base path (without extension)")
    parser.add_argument("--layout", choices=["stacked", "sidebyside"], default="stacked",
                        help="Figure layout: 'stacked' (default) or 'sidebyside'")
    parser.add_argument("--col", choices=["double", "single"], default="double",
                        help="Column sizing: 'double' (~6.8 in, default) or 'single' (~3.37 in)")
    parser.add_argument("--forces", choices=["inset", "panel", "none"], default="inset",
                        help="Force representation: 'inset' (inner inset plot, default), 'panel' (dedicated panel), or 'none'")
    args = parser.parse_args()

    input_path = Path(args.csv_file)
    if input_path.is_dir():
        csv_files = sorted(list(input_path.glob("*.csv")))
        if not csv_files:
            print(f"No CSV files found in {input_path}")
            return
        for f in csv_files:
            out_base = f.parent / f"{f.stem}_paper_{args.layout}" if args.output is None else Path(args.output) / f"{f.stem}_paper_{args.layout}"
            plot_paper_figure(f, output_base=out_base, layout=args.layout, col_width=args.col, forces=args.forces)
    else:
        plot_paper_figure(input_path, output_base=args.output, layout=args.layout, col_width=args.col, forces=args.forces)

if __name__ == "__main__":
    main()
