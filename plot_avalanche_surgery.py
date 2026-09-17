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

    import re
    match = re.search(r'config_(\d+)_to_(\d+)_step_(\d+)', csv_file.stem)
    if match:
        pre_cfg, post_cfg, step_str = match.groups()
        title_str = f"Avalanche Micro-Surgery Trace — Load Step {int(step_str)} (configuration_{pre_cfg}.vtk → configuration_{post_cfg}.vtk)"
    else:
        title_str = f"Avalanche Micro-Surgery Trace — Load Step {load_step}"

    # 1. Energy panel
    ax_energy.plot(x, energy, color=color_lbfgs, lw=2.2, label="L-BFGS Relaxation ($E$)", zorder=3)
    ax_energy.set_ylabel("Energy $E$", fontsize=12, fontweight='bold')
    ax_energy.set_title(title_str, fontsize=14, fontweight='bold', pad=12)

    # Collect Delta E_topo for each phase
    phase_delta_e = {}
    for i, ev in enumerate(event_type):
        if ev == 'AFTER_REMESH':
            phase_delta_e[phase[i]] = energy_change[i]

    y_min, y_max = np.min(energy), np.max(energy)
    y_range = y_max - y_min

    # Clean, unboxed phase labels along the top margin with Delta E (no boxes, no arrows)
    has_narrow_phase = any((pb[2] - pb[1]) < 220 for pb in phase_boundaries)
    for p_idx, (p_name, p_start, p_end) in enumerate(phase_boundaries):
        clean_name = p_name.replace("_", " ").title().replace("Remesh Pass", "Re-construct")
        if p_name in phase_delta_e:
            de_val = phase_delta_e[p_name]
            label_text = f"{clean_name}\n({de_val:+.2f})"
        else:
            label_text = clean_name

        mid_x = 0.5 * (p_start + p_end)
        # Alternate tiers if any phase is narrow so wider labels like "Re-construct" never collide
        y_pos = y_max + (y_range * 0.12 if (has_narrow_phase and p_idx % 2 == 1) else y_range * 0.04)
        ax_energy.text(mid_x, y_pos, label_text, ha='center', va='bottom',
                       fontsize=8.0, fontweight='bold', color='#495057', zorder=8)

    # Vertical event lines (clean dashed/dotted lines, no arrows)
    for i, ev in enumerate(event_type):
        step_val = x[i]
        if ev == 'BEFORE_REMESH':
            for ax in (ax_energy, ax_stress, ax_grad):
                ax.axvline(step_val, color=color_remesh_line, linestyle='--', alpha=0.75, lw=1.2, zorder=2)
        elif ev == 'REMESH_ACCEPTED':
            for ax in (ax_energy, ax_stress, ax_grad):
                ax.axvline(step_val, color=color_accept_line, linestyle=':', alpha=0.85, lw=1.4, zorder=2)

    ax_energy.set_ylim(y_min - y_range * 0.05, y_max + y_range * 0.22)

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

    # Mark Rejected remesh pass (e.g. Pass 6) with red X (no box, no arrow)
    for i, ev in enumerate(event_type):
        if ev == 'REMESH_REJECTED':
            rej_x = x[i]
            rej_E = energy[i]
            ax_energy.scatter([rej_x], [rej_E], color='#dc3545', marker='X', s=90, edgecolor='#721c24', lw=1.2, zorder=9)
            ax_energy.text(rej_x, rej_E + y_range * 0.025, "Rejected", ha='center', va='bottom',
                           fontsize=8.0, fontweight='bold', color='#dc3545', zorder=9)

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

    # Annotate Final Accepted State on Energy panel with golden star & vertical guideline (no arrows, no boxes)
    ax_energy.scatter([final_x], [final_E], color='#d4ac0d', edgecolor='#7d6608', s=240, marker='*', zorder=10)
    for ax in (ax_energy, ax_stress, ax_grad):
        ax.axvline(final_x, color='#d4ac0d', linestyle=':', lw=1.6, alpha=0.85, zorder=2)

    if final_stress != 0.0:
        ax_stress.scatter([final_x], [final_stress], color='#d4ac0d', edgecolor='#7d6608', s=180, marker='*', zorder=10)

    # Custom legend for events with exact values in the legend (saving plot space)
    custom_lines = [
        Line2D([0], [0], color=color_lbfgs, lw=2.2, label="L-BFGS Relaxation ($E$)"),
        Line2D([0], [0], color=color_remesh_line, linestyle='--', lw=1.5, label="Topology Reconnected"),
        Line2D([0], [0], color=color_accept_line, linestyle=':', lw=1.5, label="Remesh Accepted"),
        Line2D([0], [0], marker='*', color='#fef9e7', markerfacecolor='#d4ac0d', markeredgecolor='#7d6608',
               markersize=14, label=f"Final Accepted State ($E = {final_E:.6f}$)")
    ]
    ax_energy.legend(handles=custom_lines, loc="lower left", frameon=True, framealpha=0.92, fontsize=8.5)

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
