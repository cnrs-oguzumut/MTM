#!/usr/bin/env python3
"""
Plot Alpha vs PostEnergy and Alpha vs PostStress comparison for 300x300 system:
- Remeshing: /Users/usalman/runs_300X300/run_01_positive_seed42_alpha3
- No-Remeshing: /Users/usalman/runs_300X300/run_01_positive_seed42_noremesh_alpha3
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

def load_data(csv_path):
    with open(csv_path, "r") as f:
        reader = list(csv.DictReader(f))
    alpha = np.array([float(r["Alpha"]) for r in reader])
    energy = np.array([float(r["PostEnergy"]) for r in reader])
    stress = np.array([float(r["PostStress"]) for r in reader])
    area = np.array([float(r["PostArea"]) for r in reader])
    return alpha, energy, stress, area

def main():
    setup_style()

    base_dir = Path("/Users/usalman/runs_300X300")
    dir_remesh = base_dir / "run_01_positive_seed42_alpha3"
    dir_noremesh = base_dir / "run_01_positive_seed42_noremesh_alpha3"

    file_remesh = dir_remesh / "energy_stress_log.csv"
    file_noremesh = dir_noremesh / "energy_stress_log.csv"

    print("Loading data...")
    a_r, e_r, s_r, area_r = load_data(file_remesh)
    a_nr, e_nr, s_nr, area_nr = load_data(file_noremesh)

    print(f"Remeshing: {len(a_r)} points, alpha [{a_r.min():.3f}, {a_r.max():.3f}]")
    print(f"No-Remeshing: {len(a_nr)} points, alpha [{a_nr.min():.3f}, {a_nr.max():.3f}]")

    # Dynamic yield detection
    idx_peak_r = np.argmax(s_r)
    alpha_yield_r = a_r[idx_peak_r]
    stress_yield_r = s_r[idx_peak_r]

    idx_peak_nr = np.argmax(s_nr)
    alpha_yield_nr = a_nr[idx_peak_nr]
    stress_yield_nr = s_nr[idx_peak_nr]

    print(f"Yield (Remesh):    alpha = {alpha_yield_r:.5f}, sigma = {stress_yield_r:.5f}")
    print(f"Yield (No-Remesh): alpha = {alpha_yield_nr:.5f}, sigma = {stress_yield_nr:.5f}")

    # Color palette
    color_remesh = "#1f77b4"     # blue
    color_noremesh = "#d95f02"   # vermilion / orange

    # =========================================================================
    # 2-Panel Side-by-Side Figure
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8), dpi=300)
    fig.subplots_adjust(wspace=0.25, left=0.08, right=0.96, top=0.92, bottom=0.13)

    # -------------------------------------------------------------------------
    # Panel (a): Alpha vs PostEnergy
    # -------------------------------------------------------------------------
    ax1.plot(a_r, e_r, color=color_remesh, lw=1.6, label="With remeshing")
    ax1.plot(a_nr, e_nr, color=color_noremesh, lw=1.6, ls="--", label="Without remeshing")

    ax1.set_xlabel(r"Shear strain $\alpha$")
    ax1.set_ylabel(r"Post-relaxation energy $E_{\mathrm{post}}$")
    ax1.set_xlim(0.14, max(a_r.max(), a_nr.max()) * 1.01)
    ax1.set_ylim(0, max(e_nr.max(), e_r.max()) * 1.05)
    ax1.grid(True, ls=":", color="#cccccc", alpha=0.7)
    ax1.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.92)

    ax1.text(0.04, 0.88, r"$\mathbf{(a)}$", transform=ax1.transAxes,
             fontsize=13.0, fontweight="bold", va="top", ha="left")

    # Inset or callout on final energy difference
    callout_E = (
        f"End energy (remesh):    {e_r[-1]:.1f}\n"
        f"End energy (no-remesh): {e_nr[-1]:.1f}\n"
        f"Energy reduction:       {(1.0 - e_r[-1]/e_nr[-1])*100:.1f}%"
    )
    ax1.text(0.42, 0.25, callout_E, transform=ax1.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

    # -------------------------------------------------------------------------
    # Panel (b): Alpha vs PostStress
    # -------------------------------------------------------------------------
    ax2.plot(a_r, s_r, color=color_remesh, lw=1.4, label="With remeshing")
    ax2.plot(a_nr, s_nr, color=color_noremesh, lw=1.4, ls="--", label="Without remeshing")

    # Mark yield peaks
    ax2.plot(alpha_yield_r, stress_yield_r, marker="o", markersize=5.0, color=color_remesh, zorder=5)
    ax2.plot(alpha_yield_nr, stress_yield_nr, marker="s", markersize=5.0, color=color_noremesh, zorder=5)

    ax2.set_xlabel(r"Shear strain $\alpha$")
    ax2.set_ylabel(r"Post-relaxation shear stress $\sigma_{xy, \mathrm{post}}$")
    ax2.set_xlim(0.14, max(a_r.max(), a_nr.max()) * 1.01)
    ax2.set_ylim(-0.01, max(s_r.max(), s_nr.max()) * 1.15)
    ax2.grid(True, ls=":", color="#cccccc", alpha=0.7)
    ax2.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.92)

    ax2.text(0.04, 0.88, r"$\mathbf{(b)}$", transform=ax2.transAxes,
             fontsize=13.0, fontweight="bold", va="top", ha="left")

    callout_S = (
        f"Peak $\\sigma_{{xy}}$ (remesh):    {stress_yield_r:.4f} ($\\alpha={alpha_yield_r:.3f}$)\n"
        f"Peak $\\sigma_{{xy}}$ (no-remesh): {stress_yield_nr:.4f} ($\\alpha={alpha_yield_nr:.3f}$)"
    )
    ax2.text(0.04, 0.08, callout_S, transform=ax2.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

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
