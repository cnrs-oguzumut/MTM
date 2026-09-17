#!/usr/bin/env python3
"""
plot_avalanche_surgery.py

Visualizes the detailed micro-event surgery of an avalanche from an
avalanche_trace/step_XXXXX_surgery.csv file.

Panels:
1. Total lattice energy E vs micro-step (continuous L-BFGS slopes + discrete topological jumps)
2. Macroscopic shear stress σ_xy vs micro-step
3. Force equilibrium metric ||g||_inf on a logarithmic scale
"""

import sys
import argparse
from pathlib import Path
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_surgery(csv_path, output_png=None):
    csv_file = Path(csv_path)
    if not csv_file.exists():
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)

    rows = []
    with open(csv_file, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    if not rows:
        print(f"Error: Empty CSV: {csv_file}")
        sys.exit(1)

    load_step = rows[0]['load_step']

    # Extract columns
    global_micro_step = np.array([int(r['global_micro_step']) for r in rows])
    energy = np.array([float(r['energy']) for r in rows])
    stress = np.array([float(r['stress']) for r in rows])
    grad_norm = np.array([float(r['grad_norm']) for r in rows])
    energy_change = np.array([float(r['energy_change']) for r in rows])
    event_type = [r['event_type'] for r in rows]

    phase = [r['phase'] for r in rows]

    # Style
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax_energy, ax_stress, ax_grad) = plt.subplots(3, 1, figsize=(13, 10), sharex=True)

    x = global_micro_step

    # Colors
    color_lbfgs = '#1f77b4'
    color_remesh_line = '#d62728'
    color_accept_line = '#2ca02c'

    # Shaded phase bands
    phase_colors = ['#f8f9fa', '#e8f4f8', '#fef9e7', '#f4ecf7', '#e8f8f5', '#fdf2e9']
    phase_boundaries = []
    curr_phase = phase[0]
    start_step = x[0]
    for i in range(1, len(x)):
        if phase[i] != curr_phase:
            phase_boundaries.append((curr_phase, start_step, x[i-1]))
            curr_phase = phase[i]
            start_step = x[i]
    phase_boundaries.append((curr_phase, start_step, x[-1]))

    for p_idx, (p_name, p_start, p_end) in enumerate(phase_boundaries):
        bg_col = phase_colors[p_idx % len(phase_colors)]
        for ax in (ax_energy, ax_stress, ax_grad):
            ax.axvspan(p_start, p_end, facecolor=bg_col, alpha=0.6, zorder=0)

    # 1. Energy panel
    ax_energy.plot(x, energy, color=color_lbfgs, lw=2.2, label="Total Energy $E$", zorder=3)
    ax_energy.set_ylabel("Energy $E$", fontsize=12, fontweight='bold')
    ax_energy.set_title(f"Avalanche Micro-Surgery Trace — Load Step {load_step}", fontsize=14, fontweight='bold', pad=12)

    # Annotate phase labels at top of energy panel
    y_min, y_max = np.min(energy), np.max(energy)
    y_label_pos = y_max + (y_max - y_min) * 0.04
    for p_idx, (p_name, p_start, p_end) in enumerate(phase_boundaries):
        clean_name = p_name.replace("_", " ").title()
        mid_x = 0.5 * (p_start + p_end)
        ax_energy.text(mid_x, y_label_pos, clean_name, ha='center', va='bottom',
                       fontsize=9, fontweight='bold', color='#495057',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ced4da', alpha=0.85))
    ax_energy.set_ylim(y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.18)

    # 2. Stress panel
    mask_stress = stress != 0.0
    if mask_stress.any():
        ax_stress.plot(x[mask_stress], stress[mask_stress], color='#ff7f0e', lw=2.0, marker='o', markersize=3, label="Shear Stress $\\sigma_{xy}$", zorder=3)
    else:
        ax_stress.plot(x, stress, color='#ff7f0e', lw=2.0, label="Shear Stress $\\sigma_{xy}$", zorder=3)
    ax_stress.set_ylabel("Stress $\\sigma_{xy}$", fontsize=12, fontweight='bold')

    # 3. Gradient Norm panel (log scale)
    mask_grad = grad_norm > 0.0
    if mask_grad.any():
        ax_grad.semilogy(x[mask_grad], grad_norm[mask_grad], color='#9467bd', lw=1.8, marker='.', markersize=4, label="$\\|\\nabla E\\|_\\infty$", zorder=3)
    ax_grad.set_ylabel("$\\|\\nabla E\\|_\\infty$ (log)", fontsize=12, fontweight='bold')
    ax_grad.set_xlabel("Global Micro-Step", fontsize=12, fontweight='bold')

    # Add vertical event markers for remeshing events
    for i, ev in enumerate(event_type):
        step_val = x[i]
        if ev == 'BEFORE_REMESH':
            for ax in (ax_energy, ax_stress, ax_grad):
                ax.axvline(step_val, color=color_remesh_line, linestyle='--', alpha=0.7, lw=1.2, zorder=2)
        elif ev == 'AFTER_REMESH':
            # Annotate topological jump on energy plot
            de_topo = energy_change[i]
            ax_energy.annotate(f"$\\Delta E_{{topo}} = {de_topo:+.2e}$",
                               xy=(step_val, energy[i]),
                               xytext=(step_val + max(1, len(x)*0.015), energy[i]),
                               arrowprops=dict(arrowstyle="->", color=color_remesh_line, lw=1.5),
                               fontsize=9, fontweight='bold', color=color_remesh_line,
                               bbox=dict(boxstyle="round,pad=0.3", fc="#ffebee", ec=color_remesh_line, alpha=0.9),
                               zorder=5)
        elif ev == 'REMESH_ACCEPTED':
            for ax in (ax_energy, ax_stress, ax_grad):
                ax.axvline(step_val, color=color_accept_line, linestyle=':', alpha=0.8, lw=1.5, zorder=2)

    # Identify the final accepted state
    final_accepted_idx = None
    for i in range(len(event_type) - 1, -1, -1):
        if event_type[i] == 'REMESH_ACCEPTED':
            final_accepted_idx = i
            break
    if final_accepted_idx is None:
        for i in range(len(event_type) - 1, -1, -1):
            if event_type[i] == 'LBFGS_CONVERGED':
                final_accepted_idx = i
                break
    if final_accepted_idx is None:
        final_accepted_idx = len(x) - 1

    final_x = x[final_accepted_idx]
    final_E = energy[final_accepted_idx]
    final_stress = stress[final_accepted_idx]

    # Annotate Final Accepted State on Energy panel
    ax_energy.scatter([final_x], [final_E], color='#d4ac0d', edgecolor='#7d6608', s=240, marker='*', zorder=10)
    # Position text cleanly
    text_offset_x = -max(1, len(x) * 0.14) if final_x > len(x) * 0.6 else max(1, len(x) * 0.05)
    ax_energy.annotate(f"Final Accepted State\n$E = {final_E:.6f}$",
                       xy=(final_x, final_E),
                       xytext=(final_x + text_offset_x, final_E + (y_max - y_min) * 0.12),
                       arrowprops=dict(arrowstyle="->", color='#b7950b', lw=2),
                       fontsize=9, fontweight='bold', color='#7d6608',
                       bbox=dict(boxstyle="round,pad=0.35", fc="#fef9e7", ec="#d4ac0d", lw=1.5, alpha=0.95),
                       zorder=11)

    if final_stress != 0.0:
        ax_stress.scatter([final_x], [final_stress], color='#d4ac0d', edgecolor='#7d6608', s=180, marker='*', zorder=10)

    # Custom legend for events
    custom_lines = [
        Line2D([0], [0], color=color_lbfgs, lw=2, label="L-BFGS Relaxation"),
        Line2D([0], [0], color=color_remesh_line, linestyle='--', lw=1.5, label="Topology Reconnected"),
        Line2D([0], [0], color=color_accept_line, linestyle=':', lw=1.5, label="Remesh Accepted"),
        Line2D([0], [0], marker='*', color='#fef9e7', markerfacecolor='#d4ac0d', markeredgecolor='#7d6608', markersize=14, label="Final Accepted State")
    ]
    ax_energy.legend(handles=custom_lines, loc="lower left", frameon=True, framealpha=0.9)

    for ax in (ax_energy, ax_stress, ax_grad):
        ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    if output_png is None:
        output_png = csv_file.with_suffix('.png')
    else:
        output_png = Path(output_png)

    plt.savefig(output_png, dpi=200)
    plt.close()
    print(f"✓ Avalanche surgery plot saved to: {output_png}")

def main():
    parser = argparse.ArgumentParser(description="Plot avalanche surgery micro-trace")
    parser.add_argument("csv_file", help="Path to step_XXXXX_surgery.csv")
    parser.add_argument("-o", "--output", help="Output PNG path", default=None)
    args = parser.parse_args()

    plot_surgery(args.csv_file, args.output)

if __name__ == "__main__":
    main()
