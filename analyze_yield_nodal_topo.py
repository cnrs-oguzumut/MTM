#!/usr/bin/env python3
r"""
analyze_yield_nodal_topo.py

Performs comprehensive comparative analysis of LOCAL NODAL ENERGY CHANGES
\Delta e_a = e_a^{\text{after remesh}} - e_a^{\text{before remesh}}
for plastic avalanches occurring BEFORE / AT YIELD (\alpha <= 0.20)
versus AFTER YIELD (\alpha > 0.20) in 100x100 crystal lattices.

Generates publication-quality 4-panel comparative visualization:
  (a) Stress-Strain curve \sigma_{xy}(\alpha) with yield point and avalanche markers
  (b) Linear probability density p(\Delta e_a) before vs. after yield
  (c) Log-Log magnitude spectrum p(|\Delta e_a|) comparing barrier tails
  (d) Local energy balance comparison (% barrier increase vs % local relaxation)
"""

import os
import sys
import glob
import csv
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def setup_publication_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "axes.titlesize": 10.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.2,
        "mathtext.fontset": "cm",
        "lines.antialiased": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

def load_pts_energy(vtk_path):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(vtk_path))
    r.Update()
    grid = r.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    ne = vtk_to_numpy(grid.GetPointData().GetArray('NodalEnergy'))
    return pts, ne

def load_stress_curve(csv_path):
    alphas, stresses, steps = [], [], []
    if not os.path.exists(csv_path):
        return np.array([]), np.array([]), np.array([])
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            alphas.append(abs(float(r['Alpha'])))
            stresses.append(float(r['PostStress']))
            steps.append(int(r['Iteration']))
    return np.array(steps), np.array(alphas), np.array(stresses)

def collect_nodal_by_regime(directories, yield_alpha_cutoff=0.20):
    diffs_before_yield = []
    diffs_after_yield = []
    
    avalanche_events = []
    seen_passes = set()

    for d in directories:
        surgery_dir = Path(d) / "vtk_surgery"
        log_file = Path(d) / "energy_stress_log.csv"
        
        step_to_alpha = {}
        if log_file.exists():
            with open(log_file) as f:
                for r in csv.DictReader(f):
                    step_to_alpha[int(r['Iteration'])] = abs(float(r['Alpha']))

        before_files = sorted(surgery_dir.glob("*_1_before_remesh.vtk"))
        print(f"Directory {d}: found {len(before_files)} surgery passes.")
        
        for bf in before_files:
            af = Path(str(bf).replace('_1_before_remesh.vtk', '_2_after_remesh.vtk'))
            if not af.exists():
                continue
            
            # Deduplicate by filename
            fname = bf.name
            if fname in seen_passes:
                continue
            seen_passes.add(fname)

            # Determine load step from filename (step_XXXXX)
            try:
                parts = fname.split('_')
                step_idx = int(parts[1])
            except:
                step_idx = 0
            
            alpha_val = step_to_alpha.get(step_idx, 0.140)

            pts1, ne1 = load_pts_energy(bf)
            pts2, ne2 = load_pts_energy(af)

            # Coordinate matching
            h2 = { (round(p[0], 3), round(p[1], 3)): ne2[i] for i, p in enumerate(pts2) }
            
            diffs = []
            for i, p in enumerate(pts1):
                k = (round(p[0], 3), round(p[1], 3))
                if k in h2:
                    diffs.append(h2[k] - ne1[i])
            
            diffs = np.array(diffs)
            active = diffs[np.abs(diffs) > 1e-6]
            
            if len(active) == 0:
                continue

            is_pre_yield = (alpha_val <= yield_alpha_cutoff)
            if is_pre_yield:
                diffs_before_yield.extend(active)
            else:
                diffs_after_yield.extend(active)

            avalanche_events.append({
                'step': step_idx,
                'alpha': alpha_val,
                'is_pre_yield': is_pre_yield,
                'active_count': len(active),
                'max_jump': np.max(active),
                'mean_jump': np.mean(active),
                'sum_jump': np.sum(diffs)
            })

    diffs_before_yield = np.array(diffs_before_yield)
    diffs_after_yield = np.array(diffs_after_yield)
    return diffs_before_yield, diffs_after_yield, avalanche_events

def plot_yield_comparison(diffs_pre, diffs_post, steps, alphas, stresses, out_dir, cutoff_alpha=0.20):
    setup_publication_style()
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(2, 2, figsize=(7.8, 6.4))
    fig.subplots_adjust(hspace=0.38, wspace=0.34, left=0.10, right=0.96, top=0.93, bottom=0.10)

    col_pre = "#1f77b4"     # Blue: Pre/at-yield
    col_post = "#d95f02"    # Orange: Post-yield flow

    # -------------------------------------------------------------
    # Panel (a): Stress-Strain Curve \sigma_{xy}(\alpha)
    # -------------------------------------------------------------
    ax = axs[0, 0]
    if len(alphas) > 0:
        ax.plot(alphas, stresses, color="#252525", lw=1.2, label=r"$\sigma_{xy}(\alpha)$ trajectory")
        
        # Shade pre-yield vs post-yield regions
        ax.axvspan(alphas[0], cutoff_alpha, color="#c6dbef", alpha=0.45, label=r"Yield regime ($\alpha \leq 0.20$)")
        ax.axvspan(cutoff_alpha, alphas[-1], color="#fdd0a2", alpha=0.35, label=r"Flow regime ($\alpha > 0.20$)")
        
        # Mark yield threshold line
        ax.axvline(cutoff_alpha, color="#b30000", ls="--", lw=1.4)
        
        ax.set_xlabel(r"Applied shear strain $\alpha$")
        ax.set_ylabel(r"Shear stress $\sigma_{xy}$")
        ax.set_title(r"(a) Shear Stress Trajectory & Yield", fontweight="bold", loc="left")
        ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="lower left")
        ax.set_xlim(alphas[0], min(alphas[-1], 1.0))

    # -------------------------------------------------------------
    # Panel (b): Linear Probability Density p(\Delta e_a)
    # -------------------------------------------------------------
    ax = axs[0, 1]
    bins = np.linspace(-0.03, 0.06, 51)
    
    if len(diffs_pre) > 0:
        ax.hist(diffs_pre, bins=bins, density=True, color=col_pre, alpha=0.55, edgecolor="#08519c", lw=0.6,
                label=rf"Pre/At Yield ($N={len(diffs_pre):,}$)")
    if len(diffs_post) > 0:
        ax.hist(diffs_post, bins=bins, density=True, color=col_post, alpha=0.55, edgecolor="#a63603", lw=0.6,
                label=rf"Post-Yield ($N={len(diffs_post):,}$)")
    
    ax.axvline(0.0, color="gray", ls="--", lw=1.0, alpha=0.7)
    if len(diffs_pre) > 0:
        ax.axvline(np.mean(diffs_pre), color="#08519c", ls="-", lw=1.5, label=rf"Mean pre: +{np.mean(diffs_pre):.4f}")
    if len(diffs_post) > 0:
        ax.axvline(np.mean(diffs_post), color="#a63603", ls=":", lw=1.8, label=rf"Mean post: +{np.mean(diffs_post):.4f}")

    ax.set_yscale("log")
    ax.set_xlabel(r"Local nodal energy change $\Delta e_a$")
    ax.set_ylabel(r"Probability density $p(\Delta e_a)$ (log)")
    ax.set_title(r"(b) Local $\Delta e_a$ Distribution (Semi-Log)", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")
    ax.set_xlim(-0.04, 0.08)
    ax.set_ylim(5e-2, 1.5e3)
    ax.grid(True, which="both", ls=":", alpha=0.35)

    # -------------------------------------------------------------
    # Panel (c): Log-Log Magnitude Spectrum p(|\Delta e_a|)
    # -------------------------------------------------------------
    ax = axs[1, 0]
    log_bins = np.logspace(-6, 0, 45)
    if len(diffs_pre) > 0:
        ax.hist(np.abs(diffs_pre), bins=log_bins, density=True, color=col_pre, alpha=0.55, edgecolor="#08519c", lw=0.6,
                label=r"Pre/At Yield ($\alpha \leq 0.20$)")
    if len(diffs_post) > 0:
        ax.hist(np.abs(diffs_post), bins=log_bins, density=True, color=col_post, alpha=0.55, edgecolor="#a63603", lw=0.6,
                label=r"Post-Yield ($\alpha > 0.20$)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Magnitude $|\Delta e_a|$")
    ax.set_ylabel(r"Probability density $p(|\Delta e_a|)$")
    ax.set_title(r"(c) Barrier Magnitude Spectrum (Log-Log)", fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")
    ax.set_xlim(1e-6, 1.0)

    # -------------------------------------------------------------
    # Panel (d): Local Energy Balance (Bar Comparison)
    # -------------------------------------------------------------
    ax = axs[1, 1]
    
    def get_fractions(data):
        if len(data) == 0:
            return 50.0, 50.0
        pos = np.mean(data > 0) * 100
        neg = np.mean(data < 0) * 100
        return pos, neg

    pre_pos, pre_neg = get_fractions(diffs_pre)
    post_pos, post_neg = get_fractions(diffs_post)

    x_idx = np.array([0, 1])
    width = 0.35

    b1 = ax.bar(x_idx - width/2, [pre_pos, pre_neg], width=width, color=col_pre, alpha=0.8, edgecolor="black", lw=0.8,
                label=r"Pre/At Yield ($\alpha \leq 0.20$)")
    b2 = ax.bar(x_idx + width/2, [post_pos, post_neg], width=width, color=col_post, alpha=0.8, edgecolor="black", lw=0.8,
                label=r"Post-Yield ($\alpha > 0.20$)")

    ax.set_xticks(x_idx)
    ax.set_xticklabels([r"Barrier Increase ($\Delta e_a > 0$)", r"Local Relaxation ($\Delta e_a < 0$)"])
    ax.set_ylabel(r"Active node fraction (%)")
    ax.set_ylim(0, 75)
    ax.set_title(r"(d) Energetic Character Comparison", fontweight="bold", loc="left")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, loc="upper right")

    for bars in [b1, b2]:
        for b in bars:
            h = b.get_height()
            if h > 0:
                ax.text(b.get_x() + b.get_width()/2.0, h + 1.5, f"{h:.1f}%",
                        ha="center", va="bottom", fontsize=8.0, fontweight="bold")

    pdf_path = out_dir / "nodal_topo_before_after_yield.pdf"
    png_path = out_dir / "nodal_topo_before_after_yield.png"
    plt.savefig(pdf_path, dpi=300)
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Saved comparative figure to:\n  {pdf_path}\n  {png_path}")
    return png_path, pdf_path

if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parent
    dirs_to_check = [
        repo_root / "study_100x100_resumed_to_1",
        repo_root / "study_100x100_seed43",
        repo_root / "study_100x100_10avalanches"
    ]
    
    diffs_pre, diffs_post, av_events = collect_nodal_by_regime(dirs_to_check, yield_alpha_cutoff=0.20)
    
    # Load long stress trajectory from test_100x100
    steps, alphas, stresses = load_stress_curve(repo_root / "test_100x100" / "energy_stress_log.csv")
    
    out_dir = repo_root / "figures"
    plot_yield_comparison(diffs_pre, diffs_post, steps, alphas, stresses, out_dir, cutoff_alpha=0.20)
