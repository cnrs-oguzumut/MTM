#!/usr/bin/env python3
"""
Plot Alpha vs PostEnergy and Alpha vs PostStress comparison for 300x300 system:
- Reconnecting: /Users/usalman/runs_300X300/run_01_positive_seed42_alpha3
- Non-reconnecting: /Users/usalman/runs_300X300/run_01_positive_seed42_noremesh_alpha3
- Truncated at alpha = 1.25
- Labels (a) and (b) at bottom-left
"""

import csv
from pathlib import Path
import numpy as np
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
        "legend.fontsize": 9.5,
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

def load_data(csv_path, alpha_max=1.25):
    with open(csv_path, "r") as f:
        reader = list(csv.DictReader(f))
    alpha = np.array([float(r["Alpha"]) for r in reader])
    energy = np.array([float(r["PostEnergy"]) for r in reader])
    stress = np.array([float(r["PostStress"]) for r in reader])
    area = np.array([float(r["PostArea"]) for r in reader])
    
    mask = alpha <= alpha_max
    return alpha[mask], energy[mask], stress[mask], area[mask]

def main():
    setup_style()

    base_dir = Path("/Users/usalman/runs_300X300")
    dir_remesh = base_dir / "run_01_positive_seed42_alpha3"
    dir_noremesh = base_dir / "run_01_positive_seed42_noremesh_alpha3"

    file_remesh = dir_remesh / "energy_stress_log.csv"
    file_noremesh = dir_noremesh / "energy_stress_log.csv"

    alpha_max = 1.25
    print("Loading data up to alpha =", alpha_max)
    a_r, e_r, s_r, _ = load_data(file_remesh, alpha_max)
    a_nr, e_nr, s_nr, _ = load_data(file_noremesh, alpha_max)

    print(f"Reconnecting: {len(a_r)} points, alpha [{a_r.min():.3f}, {a_r.max():.3f}]")
    print(f"Non-reconnecting: {len(a_nr)} points, alpha [{a_nr.min():.3f}, {a_nr.max():.3f}]")

    # Dynamic yield detection
    idx_peak_r = np.argmax(s_r)
    alpha_yield_r = a_r[idx_peak_r]
    stress_yield_r = s_r[idx_peak_r]

    idx_peak_nr = np.argmax(s_nr)
    alpha_yield_nr = a_nr[idx_peak_nr]
    stress_yield_nr = s_nr[idx_peak_nr]

    print(f"Yield (Reconnecting):     alpha = {alpha_yield_r:.5f}, sigma = {stress_yield_r:.5f}")
    print(f"Yield (Non-reconnecting): alpha = {alpha_yield_nr:.5f}, sigma = {stress_yield_nr:.5f}")

    # Color palette
    color_remesh = "#1f77b4"     # blue
    color_noremesh = "#d95f02"   # vermilion / orange

    # =========================================================================
    # 2-Panel Side-by-Side Figure
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=300)
    fig.subplots_adjust(wspace=0.25, left=0.08, right=0.96, top=0.92, bottom=0.14)

    # -------------------------------------------------------------------------
    # Panel (a): Alpha vs PostEnergy
    # -------------------------------------------------------------------------
    ax1.plot(a_r, e_r, color=color_remesh, lw=1.6, label="Reconnecting")
    ax1.plot(a_nr, e_nr, color=color_noremesh, lw=1.6, ls="--", label="Non-reconnecting")

    ax1.set_xlabel(r"Shear strain $\alpha$")
    ax1.set_ylabel(r"Post-relaxation energy $E_{\mathrm{post}}$")
    ax1.set_xlim(0.14, alpha_max)
    ax1.set_ylim(-100, 4800)
    ax1.grid(True, ls=":", color="#cccccc", alpha=0.7)
    ax1.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.92)

    # Bottom left label (a)
    ax1.text(0.04, 0.08, r"$\mathbf{(a)}$", transform=ax1.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left",
             bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.9))

    # -------------------------------------------------------------------------
    # Panel (b): Alpha vs PostStress
    # -------------------------------------------------------------------------
    ax2.plot(a_r, s_r, color=color_remesh, lw=1.4, label="Reconnecting")
    ax2.plot(a_nr, s_nr, color=color_noremesh, lw=1.4, ls="--", label="Non-reconnecting")

    ax2.set_xlabel(r"Shear strain $\alpha$")
    ax2.set_ylabel(r"Post-relaxation shear stress $\sigma_{xy, \mathrm{post}}$")
    ax2.set_xlim(0.14, alpha_max)
    ax2.set_ylim(-0.012, 0.235)
    ax2.grid(True, ls=":", color="#cccccc", alpha=0.7)
    ax2.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.92)

    # Bottom left label (b)
    ax2.text(0.04, 0.08, r"$\mathbf{(b)}$", transform=ax2.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left",
             bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.9))

    # Save destinations
    repo_fig_dir = Path(__file__).resolve().parent / "figures"
    repo_fig_dir.mkdir(parents=True, exist_ok=True)

    out_png1 = repo_fig_dir / "energy_stress_comparison_300x300.png"
    out_pdf1 = repo_fig_dir / "energy_stress_comparison_300x300.pdf"

    plt.savefig(out_pdf1)
    plt.savefig(out_png1)
    plt.close()

    # Also save inside the runs folder for direct access
    out_png2 = base_dir / "energy_stress_comparison_300x300.png"
    out_pdf2 = base_dir / "energy_stress_comparison_300x300.pdf"
    import shutil
    shutil.copyfile(out_png1, out_png2)
    shutil.copyfile(out_pdf1, out_pdf2)

    print(f"\nSaved combined comparison plot to:\n  {out_png1}\n  {out_pdf1}")
    print(f"Copied to runs folder:\n  {out_png2}\n  {out_pdf2}")

if __name__ == "__main__":
    main()
