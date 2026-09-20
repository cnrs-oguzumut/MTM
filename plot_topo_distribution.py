#!/usr/bin/env python3
r"""
plot_topo_distribution.py

Collects topological reconnection energy jumps \Delta E_{\text{topo}}^{(k)} across multiple
avalanches and remeshing passes in 100x100 crystal lattices.

Generates publication-quality figures showing:
  (a) Empirical distribution (histogram + KDE) of reconnection jumps \Delta E_{\text{topo}}
  (b) Topological jump size vs. pass number inside avalanches (Pass 1 vs subsequent passes)
  (c) Complementary cumulative distribution function (survival function) P(\Delta E > x)
  (d) Net topological energy cost per avalanche along the loading trajectory (\alpha)
"""

import os
import sys
import glob
import csv
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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

def parse_surgery_csv(filepath):
    """Extracts all AFTER_REMESH events and their load steps, phases, and energy changes."""
    records = []
    with open(filepath, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r.get('event_type') == 'AFTER_REMESH':
                phase_str = r.get('phase', '')
                try:
                    pass_num = int(phase_str.split('_')[-1])
                except:
                    pass_num = 1
                records.append({
                    'file': str(filepath),
                    'load_step': int(r.get('load_step', 0)),
                    'phase': phase_str,
                    'pass_num': pass_num,
                    'delta_E_topo': float(r['energy_change']),
                    'energy': float(r.get('energy', 0.0)),
                    'stress': float(r.get('stress', 0.0)),
                })
    return records

def collect_all_data(repo_root):
    """Collects all surgery records from 100x100 simulation folders."""
    # 1. study_100x100_10avalanches (alpha in [0.14, 0.20])
    dir1 = repo_root / "study_100x100_10avalanches" / "vtk_surgery"
    files1 = sorted(dir1.glob("*_avalanche_surgery.csv"))
    
    # 2. test_100x100 (alpha in [0.14, 1.0])
    dir2 = repo_root / "test_100x100" / "avalanche_trace"
    files2 = sorted(dir2.glob("*_surgery.csv"))
    
    # Load step to alpha mapping for dir1
    step_to_alpha_1 = {}
    csv1_log = repo_root / "study_100x100_10avalanches" / "energy_stress_log.csv"
    if csv1_log.exists():
        with open(csv1_log) as f:
            for row in csv.DictReader(f):
                step_to_alpha_1[int(row['Iteration'])] = abs(float(row['Alpha']))

    # Load step to alpha mapping for dir2
    step_to_alpha_2 = {}
    csv2_log = repo_root / "test_100x100" / "energy_stress_log.csv"
    if csv2_log.exists():
        with open(csv2_log) as f:
            for row in csv.DictReader(f):
                step_to_alpha_2[int(row['Iteration'])] = abs(float(row['Alpha']))

    dataset_1 = []
    for f in files1:
        recs = parse_surgery_csv(f)
        for r in recs:
            r['alpha'] = step_to_alpha_1.get(r['load_step'], 0.0)
            r['dataset'] = "study_100x100_10avalanches"
        dataset_1.extend(recs)

    dataset_2 = []
    for f in files2:
        recs = parse_surgery_csv(f)
        for r in recs:
            r['alpha'] = step_to_alpha_2.get(r['load_step'], 0.0)
            r['dataset'] = "test_100x100"
        dataset_2.extend(recs)

    return dataset_1, dataset_2

def plot_topological_distribution(dataset_1, dataset_2, out_dir):
    setup_publication_style()
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract values
    jumps_1 = np.array([r['delta_E_topo'] for r in dataset_1])
    passes_1 = np.array([r['pass_num'] for r in dataset_1])
    steps_1 = np.array([r['load_step'] for r in dataset_1])
    alphas_1 = np.array([r['alpha'] for r in dataset_1])

    jumps_2 = np.array([r['delta_E_topo'] for r in dataset_2])
    alphas_2 = np.array([r['alpha'] for r in dataset_2])
    # Subset of test_100x100 up to alpha=0.30
    sub03_mask = alphas_2 <= 0.30 + 1e-4
    jumps_2_sub03 = jumps_2[sub03_mask]

    # Combined dataset (all 100x100 passes)
    all_jumps = np.concatenate([jumps_1, jumps_2])

    print("==================================================")
    print("TOPOLOGICAL RECONNECTION ENERGY JUMP STATISTICS")
    print("==================================================")
    print(f"Dataset 1 (study_100x100, 11 avalanches, alpha in [0.14, 0.20]):")
    print(f"  Total passes: {len(jumps_1)}")
    print(f"  Mean:   {np.mean(jumps_1):.4f}")
    print(f"  Std:    {np.std(jumps_1):.4f}")
    print(f"  Median: {np.median(jumps_1):.4f}")
    print(f"  IQR:    [{np.percentile(jumps_1, 25):.4f}, {np.percentile(jumps_1, 75):.4f}]")
    print(f"  Range:  [{np.min(jumps_1):.4f}, {np.max(jumps_1):.4f}]")
    print(f"\nDataset 2 (test_100x100, 19 avalanches, alpha in [0.14, 1.00]):")
    print(f"  Total passes: {len(jumps_2)} (up to alpha=0.30: {len(jumps_2_sub03)})")
    print(f"  Mean:   {np.mean(jumps_2):.4f}")
    print(f"  Median: {np.median(jumps_2):.4f}")
    print(f"  Range:  [{np.min(jumps_2):.4f}, {np.max(jumps_2):.4f}]")
    print(f"\nCombined (121 reconnection passes across 30 avalanches):")
    print(f"  Mean:   {np.mean(all_jumps):.4f}")
    print(f"  Std:    {np.std(all_jumps):.4f}")
    print(f"  Median: {np.median(all_jumps):.4f}")
    print(f"  Range:  [{np.min(all_jumps):.4f}, {np.max(all_jumps):.4f}]")
    print("==================================================")

    # 4-Panel Figure
    fig, axs = plt.subplots(2, 2, figsize=(7.8, 6.2))
    fig.subplots_adjust(hspace=0.40, wspace=0.34, left=0.10, right=0.96, top=0.93, bottom=0.11)

    col_primary = "#1f77b4"
    col_accent = "#d95f02"
    col_green = "#2ca02c"
    col_purple = "#7570b3"

    # -------------------------------------------------------------
    # Panel (a): Histogram & KDE of \Delta E_topo
    # -------------------------------------------------------------
    ax = axs[0, 0]
    bins = np.linspace(0, 13, 27)
    counts, edges, patches = ax.hist(
        jumps_1, bins=bins, density=True, color=col_primary, alpha=0.65,
        edgecolor="black", linewidth=0.8, label=r"Study 100$\times$100 ($\alpha \leq 0.20$)"
    )
    
    # Pure numpy Kernel Density Estimate
    bandwidth = 0.45
    x_kde = np.linspace(0, 13, 300)
    # Gaussian kernel density:
    dens = np.mean(
        np.exp(-0.5 * ((x_kde[:, None] - jumps_1[None, :]) / bandwidth) ** 2)
        / (bandwidth * np.sqrt(2 * np.pi)),
        axis=1
    )
    ax.plot(x_kde, dens, color="#08519c", lw=2.0, label="KDE density")

    # Add vertical lines for mean and median
    mean_val = np.mean(jumps_1)
    median_val = np.median(jumps_1)
    ax.axvline(mean_val, color="#b30000", ls="--", lw=1.6, label=f"Mean: {mean_val:.2f}")
    ax.axvline(median_val, color="#006d2c", ls=":", lw=1.8, label=f"Median: {median_val:.2f}")

    ax.set_xlabel(r"Topological jump $\Delta E_{\rm topo}^{(k)}$")
    ax.set_ylabel(r"Probability density $p(\Delta E)$")
    ax.set_title(r"(a) Distribution of Reconnection Jumps", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")
    ax.set_xlim(0, 13)

    # -------------------------------------------------------------
    # Panel (b): Jump size vs. Pass Index k inside Avalanches
    # -------------------------------------------------------------
    ax = axs[0, 1]
    unique_passes = sorted(list(set(passes_1)))
    pass_groups = [jumps_1[passes_1 == p] for p in unique_passes]

    # Boxplot
    bp = ax.boxplot(
        pass_groups, positions=unique_passes, widths=0.55, patch_artist=True,
        boxprops=dict(facecolor="#c6dbef", edgecolor="#2171b5", lw=1.2),
        medianprops=dict(color="#b30000", lw=2.0),
        whiskerprops=dict(color="#2171b5", lw=1.2),
        capprops=dict(color="#2171b5", lw=1.2),
        flierprops=dict(marker="o", markerfacecolor=col_accent, markeredgecolor="none", markersize=4.5)
    )

    # Overlay individual points with slight jitter
    np.random.seed(42)
    for p in unique_passes:
        vals = jumps_1[passes_1 == p]
        jitter = np.random.normal(0, 0.06, size=len(vals))
        ax.scatter(p + jitter, vals, color="#08519c", s=18, alpha=0.75, zorder=4)

    # Annotate pass 1 vs pass > 1
    p1_mean = np.mean(jumps_1[passes_1 == 1])
    sub_mean = np.mean(jumps_1[passes_1 > 1])
    ax.axhline(p1_mean, color="#d95f02", ls="--", lw=1.2, alpha=0.7)
    ax.axhline(sub_mean, color="#7570b3", ls=":", lw=1.2, alpha=0.7)

    ax.set_xlabel(r"Remeshing pass index $k$ inside avalanche")
    ax.set_ylabel(r"$\Delta E_{\rm topo}^{(k)}$")
    ax.set_title(r"(b) Jump Size Across Surgery Passes", fontweight="bold", loc="left")
    ax.set_xticks(unique_passes)
    ax.set_xticklabels([f"$k={p}$" for p in unique_passes])
    ax.set_ylim(0, 13)

    # -------------------------------------------------------------
    # Panel (c): Complementary CDF (Survival Function) P(\Delta E > x)
    # -------------------------------------------------------------
    ax = axs[1, 0]
    # Sorted jumps
    sorted_jumps1 = np.sort(jumps_1)
    ccdf_1 = 1.0 - np.arange(len(sorted_jumps1)) / float(len(sorted_jumps1))
    
    sorted_jumps_comb = np.sort(all_jumps)
    ccdf_comb = 1.0 - np.arange(len(sorted_jumps_comb)) / float(len(sorted_jumps_comb))

    ax.step(sorted_jumps1, ccdf_1, where="post", color=col_primary, lw=2.2, label=r"Study 100$\times$100 ($N=51$)")
    ax.step(sorted_jumps_comb, ccdf_comb, where="post", color=col_accent, ls="--", lw=1.8, label=r"All 100$\times$100 runs ($N=121$)")

    ax.set_yscale("log")
    ax.set_xlabel(r"Threshold $x$")
    ax.set_ylabel(r"Survival probability $P(\Delta E > x)$")
    ax.set_title(r"(c) Exceedance Probability (CCDF)", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")
    ax.set_xlim(0, 13)
    ax.set_ylim(8e-3, 1.2)
    ax.grid(True, which="both", ls=":", alpha=0.5)

    # -------------------------------------------------------------
    # Panel (d): Cumulative topological cost per avalanche vs alpha
    # -------------------------------------------------------------
    ax = axs[1, 1]
    # Group by avalanche load step
    avalanche_dict = {}
    for r in dataset_1:
        st = r['load_step']
        if st not in avalanche_dict:
            avalanche_dict[st] = {
                'alpha': r['alpha'],
                'jumps': [],
                'passes': 0
            }
        avalanche_dict[st]['jumps'].append(r['delta_E_topo'])
        avalanche_dict[st]['passes'] += 1

    av_steps = sorted(avalanche_dict.keys())
    av_alphas = [avalanche_dict[s]['alpha'] for s in av_steps]
    av_total_jumps = [sum(avalanche_dict[s]['jumps']) for s in av_steps]
    av_num_passes = [avalanche_dict[s]['passes'] for s in av_steps]

    # Bar plot with color representing number of surgery passes
    bars = ax.bar(
        range(len(av_steps)), av_total_jumps, width=0.6,
        color="#3182bd", edgecolor="black", linewidth=0.8, alpha=0.85
    )

    # Annotate pass count on top of each bar
    for i, (b, n_p) in enumerate(zip(bars, av_num_passes)):
        h = b.get_height()
        ax.text(b.get_x() + b.get_width()/2.0, h + 0.5, f"{n_p}p",
                ha="center", va="bottom", fontsize=7.5, color="#08519c", fontweight="bold")

    ax.set_xlabel(r"Avalanche event index (ordered by strain $\alpha$)")
    ax.set_ylabel(r"Total avalanche cost $\sum_k \Delta E_{\rm topo}^{(k)}$")
    ax.set_title(r"(d) Net Topo Cost per Avalanche", fontweight="bold", loc="left")
    ax.set_xticks(range(len(av_steps)))
    ax.set_xticklabels([f"{a:.3f}" for a in av_alphas], rotation=45, ha="right", fontsize=7.5)
    ax.set_ylim(0, max(av_total_jumps) * 1.18)

    pdf_path = out_dir / "topo_energy_distribution_100x100.pdf"
    png_path = out_dir / "topo_energy_distribution_100x100.png"
    plt.savefig(pdf_path, dpi=300)
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Saved publication figures to:\n  {pdf_path}\n  {png_path}")
    return png_path, pdf_path

if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parent
    d1, d2 = collect_all_data(repo_root)
    out_dir = repo_root / "figures"
    plot_topological_distribution(d1, d2, out_dir)
